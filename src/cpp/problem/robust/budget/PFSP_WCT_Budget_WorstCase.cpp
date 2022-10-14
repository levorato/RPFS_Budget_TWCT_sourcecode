//
// Created by mlevorato on 5/1/20.
//

#include "PFSP_WCT_Budget_WorstCase.h"
#include "../../PFSP_Parameters.h"
#include "../../../util/precision.h"
#include "../../deterministic/PFSProblem.h"
#include "../../../util/include/TimeDateUtil.h"
#include <ilcplex/ilocplex.h>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o


// cplex time limit set to 4 hours (for large instances)
#define WORSTCASE_MIP_TIMELIMIT 14400 // 7200

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    typedef IloArray<IloNumVarArray> NumVarMatrix;
    namespace ublas = boost::numeric::ublas;
    using boost::timer::cpu_timer;

    // static variables initilization for class PFSP_WCT_Budget_WorstCase
    double PFSP_WCT_Budget_WorstCase::time_spent_worstcase_mip = 0.0;

    template<typename T>
    std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
    {
        os << "[";
        int neg_count = 0;
        for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
        {
            if((*ii) < 0) {
                neg_count++;
            } else {
                if(neg_count > 0) {
                    os << " [" << neg_count << " blank(s) ]";
                    neg_count = 0;
                }
                os << " " << *ii;
            }
        }
        os << " ]";
        return os;
    }

    /**
     * Calculates the completion time matrix C, given a problem instance (m x n), processing time matrices P_bar and
     * P_hat, a permutation pi and budget scenario m_dev.
     * @param m number of machines
     * @param n number of jobs
     * @param P_bar matrix of nominal processing times
     * @param P_hat matrix of processing time variations
     * @param pi the job permutation (job number starts at 1)
     * @param m_dev processing time scenario (i.e., for each machine r, which jobs will have their proc. time deviated )
     * @return completion time matrix C (zero-based)
     */
    ublas::matrix<double> calculate_cmax_given_scenario(int m, int n, const ublas::matrix<double> &P_bar,
                                                        const ublas::matrix<double> &P_hat, const std::vector<Job> &pi,
                                                        const std::vector<std::vector<int>> &m_dev) {
        ublas::matrix<double> P_sum(P_bar);
        for (int i = 0; i < m; i++) {
            for (int j : m_dev[i]) {
                if (j > 0) {
                    P_sum(i + 1, j) += P_hat(i + 1, j);
                }
            }
        }
        std::vector<int> perm;
        for(int x = 0; x < pi.size(); x++) {
            perm.push_back(pi[x].id);
        }
        ublas::matrix<double> C = PFSProblem::makespan(perm, m, P_sum);
        return C;
    }

    /**
     * Calculates the BigM value to be used in one of the Rob-PFSP-WCT worstcase MILP models, according to the type:
     * - type == 0 : standard upper bound for BigM (naive);
     * - type == 2 : Calculate the Makespan (Cmax) for the worst-case scenario (on all machines r, all jobs j have their
     *                  processing times deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
     * - type == 3 : Given all machines r: calculate the maximum difference between job completion times (C_{r, j}) of
     *                  a job j and any other job k, on the same machine r.
     * @param m number of machines
     * @param n number of jobs
     * @param P_bar matrix of nominal processing times
     * @param P_hat matrix of processing time variations
     * @return bigM value.
     */
    ublas::matrix<double> calculate_bigM(int m, int n, const ublas::matrix<double> &P_bar,
                                            const ublas::matrix<double> &P_hat,
                                            const std::vector<Job> &pi, const unsigned& bigM_type = 7) {
        std::vector<std::vector<int>> worst_Gamma_scenario(m, std::vector<int>());
        std::vector<std::vector<int>> best_Gamma_scenario(m, std::vector<int>());
        for (int r = 0; r < m; r++) {
            for (int i = 0; i < n; i++) {
                worst_Gamma_scenario[r].push_back(i + 1);
            }
        }
        ublas::matrix<double> max_C = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario);
        ublas::matrix<double> min_C = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, best_Gamma_scenario);
        double max_value = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario)(m,n) + 1;
        // std::cout << "max_value = " << max_value << "\n";
        // std::cout << "bigM_type = " << bigM_type << "\n";
        ublas::matrix<double> bigM = ublas::matrix<double>(m, n, max_value);
        if(bigM_type <= 1) {
            max_value = 0.0;
            for (int r = 0; r < m; r++) {
                for (int i = 0; i < n; i++) {
                    max_value += P_bar(r + 1, i + 1) + P_hat(r + 1, i + 1);
                }
            }
            bigM = ublas::matrix<double>(m, n, max_value);
        } else if(bigM_type == 2) {
            max_value = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario)(m,n) + 1;
            bigM = ublas::matrix<double>(m, n, max_value);
        } else if(bigM_type == 3) {
            double max_value = 0.0;
            for (int r = 1; r <= m; r++) {
                for (int j = 1; j <= n; j++) {
                    for(int k = 1; k <= n; k++) {
                        if(j != k) {
                            max_value = std::max(max_value, (max_C(r, j) - max_C(r, k)));
                        }
                    }
                }
            }
            bigM = ublas::matrix<double>(m, n, max_value);
        } else if(bigM_type == 4) {
            for(int i = 2; i <= m; i++) {  // [i, k-1] x [i-1, k]
                for (int k = 2; k <= n; k++) {
                    bigM(i - 1, k - 1) = max_C(i, k);
                }
            }
        } else if(bigM_type == 5) {
            for(int i = 2; i <= m; i++) {  // [i, k-1] x [i-1, k]
                for (int k = 2; k <= n; k++) {
                    if((i == 2) || (k == 2))
                        bigM(i - 1, k - 1) = max_C(i, k);
                    else
                        bigM(i - 1, k - 1) = 2 * (max_C(i, k) - std::min(min_C(i-1, k-2), min_C(i-2, k-1))) + 1;
                }
            }
        } else if(bigM_type == 6) {
            for(int i = 2; i <= m; i++) {  // [i, k-1] x [i-1, k]
                for (int k = 2; k <= n; k++) {
                    bigM(i - 1, k - 1) = 2 * (max_C(i, k) - std::min(min_C(i, k-1), min_C(i-1, k))) + 1;
                }
            }
        } else if(bigM_type == 7) {
            for(int i = 2; i <= m; i++) {  // [i, k-1] x [i-1, k]
                for (int k = 2; k <= n; k++) {
                    bigM(i - 1, k - 1) = 2 * (std::max(max_C(i, k-1) - min_C(i-1, k), max_C(i-1, k) - min_C(i, k-1))) + 1;
                }
            }
        }
        return bigM;
    }

    /**
     * Solve the Wilson-based MILP model for the RobPFSP-WCT Budget worstcase scenario.
     * For more information, see the file `src/julia/wct/robust_pfsp_wilson_worstcase_wct_mip.jl`.
     */
    NumVarMatrix populate_wilson_model(IloModel model, IloRangeArray c, IloNumVar wct,
                                         int m, int n, const ublas::matrix<double> &P_bar,
                                         const ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                         const std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w,
                                         const std::pair<double,double>& cutoff_value) {
        IloEnv env = model.getEnv();
        ublas::matrix<double> z = ublas::zero_matrix<double>(n, n);
        for (IloInt i = 0, j = 1; i < n; i++, j++) {
            z(pi[i].id - 1, j - 1) = 1.0;
        }
        ublas::matrix<double> bigM = calculate_bigM(m, n, P_bar, P_hat, pi, 7);
        // create the variables
        NumVarMatrix B(env, m);  // B >= 0
        NumVarMatrix d(env, m);  // d >= 0
        NumVarMatrix e(env, m);  // e >= 0
        NumVarMatrix min_d_e(env, m);  // min_d_e >= 0
        NumVarMatrix dev(env, m);  // 0 <= dev <= 1
        NumVarMatrix disj1(env, m);  // disk1 is bin
        NumVarMatrix disj2(env, m);  // disj2 is bin
        NumVarMatrix abs_dif_d_e(env, m);  // abs_dif_d_e >= 0
        std::stringstream ss;
        for (IloInt r = 0; r < m; r++) {
            B[r] = IloNumVarArray(env, n, 0, IloInfinity);  // B >= 0
            ss << "B(" << r << ")";
            B[r].setNames(ss.str().c_str());
            ss.str(std::string());
            d[r] = IloNumVarArray(env, n, 0, IloInfinity);  // d >= 0
            ss << "d(" << r << ")";
            d[r].setNames(ss.str().c_str());
            ss.str(std::string());
            e[r] = IloNumVarArray(env, n, 0, IloInfinity);  // e >= 0
            ss << "e(" << r << ")";
            e[r].setNames(ss.str().c_str());
            ss.str(std::string());
            min_d_e[r] = IloNumVarArray(env, n, 0, IloInfinity);  // min_d_e >= 0
            ss << "min_d_e(" << r << ")";
            min_d_e[r].setNames(ss.str().c_str());
            ss.str(std::string());
            dev[r] = IloNumVarArray(env, n, 0, 1);  // 0 <= dev <= 1
            ss << "dev(" << r << ")";
            dev[r].setNames(ss.str().c_str());
            ss.str(std::string());
            disj1[r] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // disj1 is bin
            ss << "disj1(" << r << ")";
            disj1[r].setNames(ss.str().c_str());
            ss.str(std::string());
            disj2[r] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // disj2 is bin
            ss << "disj2(" << r << ")";
            disj2[r].setNames(ss.str().c_str());
            ss.str(std::string());
            abs_dif_d_e[r] = IloNumVarArray(env, n, 0, IloInfinity);  // abs_dif_d_e >= 0
            ss << "abs_dif_d_e(" << r << ")";
            abs_dif_d_e[r].setNames(ss.str().c_str());
            ss.str(std::string());
        }
        // c1 : B[1, j] + sum_((P_hat[1, i] * dev[1, i]) * Z[i, j] for i in Jobs) + sum(P_bar[1, i] * Z[i, j] for i in Jobs) == B[1, j+1]
        // B[1][j] + sum_i ( T[1][i] * z(i,j) ) = B[1][j+1]            for all j in 1..n-1
        for (IloInt j = 0; j < n - 1; j++) {
            IloExpr expr1(env);
            expr1 += B[0][j];
            for (IloInt i = 0; i < n; i++) {
                expr1 += (P_hat(1, i+1) * z(i,j) * dev[0][i]) + (P_bar(1, i+1) * z(i,j));
            }
            model.add(expr1 == B[0][j+1]);
            expr1.end();
        }
        c.add(B[0][0] == 0);
        // c2 : B[r, 1] + sum((P_bar[r, i] + P_hat[r, i] * dev[r, i]) * Z[i, 1] for i in Jobs) == B[r+1, 1]
        // B[r][1] + sum_i ( T[r][i] * z[i][1] ) = B[r+1][1]            for r in 1..m-1
        for (IloInt r = 0; r < m - 1; r++) {
            IloExpr expr1(env);
            expr1 += B[r][0];
            for (IloInt i = 0; i < n; i++) {
                expr1 += (P_bar(r+1, i+1) + P_hat(r+1, i+1) * dev[r][i]) * z(i,0);
            }
            model.add(expr1 == B[r + 1][0]);
            expr1.end();
        }
        // For r in 2..m ; For j in 2:n :
        // c3 : B[r-1, j] + sum((P_bar[r-1, i] + P_hat[r-1, i] * dev[r-1, i]) * Z[i, j] for i in Jobs) == d[r, j]
        // c4 : B[r, j-1] + sum((P_bar[r, i] + P_hat[r, i] * dev[r, i]) * Z[i, j-1] for i in Jobs) == e[r, j]
        // B[r][j] + sum_i ( T[r][i] * z(i,j) ) <= B[r+1][j]           for r in 1..m-1, j in 2..n
        // B[r][j] + sum_i ( T[r][i] * z(i,j) ) <= B[r][j+1]           for r in 2..m, j in 1..n-1
        for (IloInt r = 1; r < m; r++) {
            for (IloInt j = 1; j < n; j++) {
                IloExpr expr1(env);
                expr1 += B[r-1][j];
                for (IloInt i = 0; i < n; i++) {
                    expr1 += (P_bar(r-1 + 1, i + 1) + P_hat(r-1 + 1, i + 1) * dev[r-1][i]) * z(i,j);
                }
                model.add(expr1 == d[r][j]);
                expr1.end();
                IloExpr expr2(env);
                expr2 += B[r][j-1];
                for (IloInt i = 0; i < n; i++) {
                    expr2 += (P_bar(r + 1, i + 1) + P_hat(r + 1, i + 1) * dev[r][i]) * z(i,j-1);
                }
                model.add(expr2 == e[r][j]);
                expr2.end();
                // c5 : min_d_e[r, j] <= d[r, j]
                model.add(min_d_e[r][j] <= d[r][j]);
                // c6 : min_d_e[r, j] <= e[r, j]
                model.add(min_d_e[r][j] <= e[r][j]);
                // c7 : abs_dif_d_e[r, j] <= (d[r, j] - e[r, j]) + bigM * disj1[r, j]
                model.add(abs_dif_d_e[r][j] <= (d[r][j] - e[r][j]) + bigM(r, j) * disj1[r][j]);
                // c8 : abs_dif_d_e[r, j] <= -(d[r, j] - e[r, j]) + bigM * disj2[r, j]
                model.add(abs_dif_d_e[r][j] <= -(d[r][j] - e[r][j]) + bigM(r, j) * disj2[r][j]);
                // c9 : disj1[r, j] + disj2[r, j] == 1
                model.add(disj1[r][j] + disj2[r][j] == 1);
                // c10 : B[r, j] <= min_d_e[r, j] + abs_dif_d_e[r, j]
                model.add(B[r][j] <= min_d_e[r][j] + abs_dif_d_e[r][j]);
            }
        }
        if(rob_params.budget_type == "machine") {
            // c11 : for r in M : sum(dev[r, i] for i in Jobs) <= Γ[r]
            ublas::vector<double> Gamma = ublas::zero_vector<double>(n);
            for(int r = 0; r < m; r++) {  // Convert the budget Gamma parameters from percentage to absolute values
                Gamma(r) = ((n * rob_params.budget_gamma[r]) / 100.0);
            }
            for (IloInt r = 0; r < m; r++) {
                IloExpr expr1(env);
                for (IloInt i = 0; i < n; i++) {
                    expr1 += dev[r][i];
                }
                c.add(expr1 == Gamma(r));
                expr1.end();
            }
        } else if (rob_params.budget_type == "global") {
            // c11 : sum(dev[r, i] for r in M, i in Jobs) <= Γ
            double Gamma = ((n * m * rob_params.budget_gamma[0]) / 100.0);
            IloExpr expr_gamma(env);
            for (IloInt r = 0; r < m; r++) {
                for (IloInt i = 0; i < n; i++) {
                    expr_gamma += dev[r][i];
                }
            }
            c.add(expr_gamma == Gamma);
            expr_gamma.end();
        } else {
            throw std::invalid_argument( "Invalid budget type!");
        }
        // set the objective function
        // c12 : (Max) sum(w[i] * sum((B[m, j] + P_bar[m, i] + P_hat[m, i] * dev[m, i]) * Z[i, j] for j in Jobs) for i in Jobs)
        // Cmax == B[m][n] + sum_i ( T[m][i] * z[i][n] )
        IloExpr expr0(env);
        for (IloInt i = 0; i < n; i++) {
            for (IloInt j = 0; j < n; j++) {
                expr0 += w(i+1) * ((B[m-1][j] + P_bar(m-1 + 1, i + 1) + P_hat(m-1 + 1, i + 1) * dev[m-1][i]) * z(i,j));
            }
        }
        model.add(c);  // add all the constraints defined above
        model.add(wct == expr0);
        expr0.end();
        if(cutoff_value.first > 0) {
            model.add(wct >= cutoff_value.first);   // lb
        }
        if(cutoff_value.second > 0) {
            // model.add(wct <= cutoff_value.second);  // ub
        }
        model.add(IloMaximize(env, wct));

        return dev;
    }

    /**
     * Solve the alternative MILP model for the RobPFSP-WCT Budget worstcase scenario (faster than Wilson model).
     * For more information, consult the file `src/julia/wct/robust_pfsp_alt_worstcase_wct_mip.jl`.
     */
    NumVarMatrix populate_alt_model(IloModel model, IloRangeArray c, IloNumVar wct,
                                int m, int n, const ublas::matrix<double> &P_bar,
                                const ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                const std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w,
                                const std::pair<double,double>& cutoff_value) {
        bool fractional_budget = false;
        IloEnv env = model.getEnv();
        ublas::matrix<double> z = ublas::zero_matrix<double>(n, n);
        for (IloInt i = 0, j = 1; i < n; i++, j++) {
            z(pi[i].id - 1, j - 1) = 1.0;
        }
        ublas::matrix<double> bigM = calculate_bigM(m, n, P_bar, P_hat, pi, 7);
        // create the variables
        NumVarMatrix C(env, m);  // C >= 0
        NumVarMatrix min_d_e(env, m);  // min_d_e >= 0
        NumVarMatrix dev(env, m);  // 0 <= dev <= 1
        NumVarMatrix disj(env, m);  // disj is bin
        NumVarMatrix abs_dif_d_e(env, m);  // abs_dif_d_e >= 0
        std::stringstream ss;
        for (IloInt r = 0; r < m; r++) {
            C[r] = IloNumVarArray(env, n, 0, IloInfinity);  // C >= 0
            ss << "C(" << r << ")";
            C[r].setNames(ss.str().c_str());
            ss.str(std::string());
            min_d_e[r] = IloNumVarArray(env, n, 0, IloInfinity);  // min_d_e >= 0
            ss << "min_d_e(" << r << ")";
            min_d_e[r].setNames(ss.str().c_str());
            ss.str(std::string());
            if(fractional_budget) {
                dev[r] = IloNumVarArray(env, n, 0, 1);  // 0 <= dev <= 1
            } else {
                dev[r] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // dev is bin
            }
            ss << "dev(" << r << ")";
            dev[r].setNames(ss.str().c_str());
            ss.str(std::string());
            disj[r] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // disj is bin
            ss << "disj(" << r << ")";
            disj[r].setNames(ss.str().c_str());
            ss.str(std::string());
            abs_dif_d_e[r] = IloNumVarArray(env, n, 0, IloInfinity);  // abs_dif_d_e >= 0
            ss << "abs_dif_d_e(" << r << ")";
            abs_dif_d_e[r].setNames(ss.str().c_str());
            ss.str(std::string());
        }
        // (C1) Special case (1): when i = 1 (machine 1)
        IloExpr expr_c1(env);
        // C[1][1] == sum(P_bar[1][j] * Z[j][1] for j in Jobs) + sum(P_hat[1][j] * dev[1][j] * Z[j][1] for j in Jobs)
        for (IloInt j = 0; j < n; j++) {
            expr_c1 += P_bar(1,j+1) * z(j,1-1) + P_hat(1,j+1) * dev[1-1][j] * z(j,1-1);
        }
        model.add(expr_c1 == C[1-1][1-1]);
        expr_c1.end();
        // (C2) for k in 2:n => C[1, k] == C[1, k - 1] + sum(P_bar[1][j] * Z[j][k] for j in Jobs) + sum(P_hat[1][j] * dev[1][j] * Z[j, k] for j in Jobs)
        for (IloInt k = 2 - 1; k <= n - 1; k++) {
            IloExpr expr1(env);
            expr1 += C[1-1][k - 1];
            for (IloInt j = 0; j < n; j++) {
                expr1 += P_bar(1,j+1) * z(j,k) + P_hat(1,j+1) * dev[1-1][j] * z(j,k);
            }
            model.add(expr1 == C[1-1][k]);
            expr1.end();
        }
        // (C3) Special case (2): when k == 1 (first position of each machine)
        // for i in 2:m => C[i, 1] <= C[i - 1, 1] + sum(P_bar[i, j] * Z[j, 1] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, 1] for j in Jobs)
        for (IloInt i = 2 - 1; i <= m - 1; i++) {
            IloExpr expr1(env);
            expr1 += C[i - 1][1-1];
            for (IloInt j = 0; j < n; j++) {
                expr1 += P_bar(i+1,j+1) * z(j,1-1) + P_hat(i+1,j+1) * dev[i][j] * z(j,1-1);
            }
            model.add(expr1 >= C[i][1-1]);
            expr1.end();
            ss.str(std::string());
        }

        // General case: (k > 2, i > 2)
        // For k in 2..n ; For i in 2:m :
        // (d) : d[i, k] == C[i, k - 1] + sum(P_bar[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs)
        // (e) : e[i, k] == C[i - 1, k] + sum(P_bar[i, j] * Z[j, k] for j in Jobs) + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs)
        for (IloInt k = 2 - 1; k <= n - 1; k++) {
            for (IloInt i = 2 - 1; i <= m - 1; i++) {
                // c5 : min_d_e[i, k] <= C[i, k - 1]
                model.add(min_d_e[i][k] <= C[i][k - 1]);
                // c6 : min_d_e[i, k] <= C[i - 1, k]
                model.add(min_d_e[i][k] <= C[i - 1][k]);
                // c7 : abs_dif_d_e[i, k] >= (C[i, k - 1] - C[i - 1, k])
                model.add(abs_dif_d_e[i][k] >= (C[i][k - 1] - C[i - 1][k]));
                // c8 : abs_dif_d_e[i, k] >= -(C[i, k - 1] - C[i - 1, k])
                model.add(abs_dif_d_e[i][k] >= -(C[i][k - 1] - C[i - 1][k]));
                // c9 : abs_dif_d_e[i, k] <= (C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * disj[i, k]
                model.add(abs_dif_d_e[i][k] <= (C[i][k - 1] - C[i - 1][k]) + bigM(i, k) * disj[i][k]);
                // c10 : abs_dif_d_e[i, k] <= -(C[i, k - 1] - C[i - 1, k]) + bigM[i, k] * (1 - disj[i, k])
                model.add(abs_dif_d_e[i][k] <= -(C[i][k - 1] - C[i - 1][k]) + bigM(i, k) * (1 - disj[i][k]));
                // c11 : C[i, k] == min_d_e[i, k] + abs_dif_d_e[i, k] + sum(T[i, j] * Z[j, k] for j in Jobs)
                //                                    + sum(P_hat[i, j] * dev[i, j] * Z[j, k] for j in Jobs);
                IloExpr expr1(env);  // sum((P_bar[i, j] + P_hat[i, j] * dev[i, j]) * Z[j, k] for j in Jobs);
                expr1 += min_d_e[i][k] + abs_dif_d_e[i][k];
                for (IloInt j = 0; j < n; j++) {
                    expr1 += (P_bar(i+1,j+1) + P_hat(i+1,j+1) * dev[i][j]) * z(j,k);
                }
                model.add(C[i][k] == expr1);
                expr1.end();
            }
        }
        if(rob_params.budget_type == "machine") {
            // c11 : for r in M : sum(dev[r, i] for i in Jobs) <= Γ[r]
            ublas::vector<double> Gamma = ublas::zero_vector<double>(n);
            for(int r = 0; r < m; r++) {  // Convert the budget Gamma parameters from percentage to absolute values
                if(fractional_budget) {
                    Gamma(r) = ((n * rob_params.budget_gamma[r]) / 100.0);
                } else {
                    Gamma(r) = ceil(((n * rob_params.budget_gamma[r]) / 100.0));
                }
            }
            for (IloInt r = 0; r < m; r++) {
                IloExpr expr1(env);
                for (IloInt i = 0; i < n; i++) {
                    expr1 += dev[r][i];
                }
                c.add(expr1 == Gamma(r));
                expr1.end();
            }
        } else if (rob_params.budget_type == "global") {
            // c11 : sum(dev[r, i] for r in M, i in Jobs) <= Γ
            double Gamma = ((n * m * rob_params.budget_gamma[0]) / 100.0);
            if(! fractional_budget) {
                Gamma = ceil(Gamma);
            }
            IloExpr expr_gamma(env);
            for (IloInt r = 0; r < m; r++) {
                for (IloInt i = 0; i < n; i++) {
                    expr_gamma += dev[r][i];
                }
            }
            c.add(expr_gamma == Gamma);
            expr_gamma.end();
        } else {
            throw std::invalid_argument( "Invalid budget type!");
        }
        // set the objective function
        // c12 : (Max) sum(w[j] * sum(C[m, k] * Z[j, k] for k in Jobs) for j in Jobs)
        IloExpr expr0(env);
        for (IloInt j = 0; j < n; j++) {
            for (IloInt k = 0; k < n; k++) {
                expr0 += w(j+1) * (C[m-1][k] * z(j,k));
            }
        }
        model.add(c);  // add all the constraints defined above
        model.add(wct == expr0);
        ss.str(std::string());
        expr0.end();
        if(cutoff_value.first > 0) {
            model.add(wct >= cutoff_value.first);   // lb
        }
        if(cutoff_value.second > 0) {
            //model.add(wct <= cutoff_value.second);  // ub
        }
        model.add(IloMaximize(env, wct));
        return dev;
    }

    void PFSP_WCT_Budget_WorstCase::release_cplex_model() {
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);
        model.end();
        // std::cout << "[PFSP_WCT_Budget_WorstCase] CPLEX model memory released successfully." << std::endl;
    }

    /**
     * Available MILP models:
     * Model == 0 : Alt MILP model.
     * Model == 1 : Wilson MILP model.
     */
    PFSPBudgetScenario PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(const unsigned &n,
                             const boost::numeric::ublas::matrix<double> &P_bar,
                             const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                             const std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w,
                             const bool& first_call, const bool& release_cplex, const unsigned& mip_model,
                             const bool& verbose, const std::pair<double,double>& cutoff_value) {
        // Create an IloCplex instance and load the model.
        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);
        IloNum start;
        boost::timer::cpu_timer timer;
        timer.start();  // Time measurement
        try {
            if(! verbose)  cplex.setOut(env.getNullStream());  // mute cplex output
            /* Sets the rule for selecting the next node to process when the search is backtracking.
              The depth-first search strategy chooses the most recently created node.
              The best-bound strategy chooses the node with the best objective function for the associated LP relaxation.
              The best-estimate strategy selects the node with the best estimate of the integer objective value that would be
              obtained from a node once all integer infeasibilities are removed. An alternative best-estimate search is also available. */
            /*  0   CPX_NODESEL_DFS 	Depth-first search
                1 	CPX_NODESEL_BESTBOUND 	Best-bound search; default
                2 	CPX_NODESEL_BESTEST 	Best-estimate search
                3 	CPX_NODESEL_BESTEST_ALT */
            // cplex.setParam(IloCplex::NodeSel, CPX_NODESEL_BESTEST);  // IloCplex.NodeSelect.BestBound , BestEst, BestEstAlt, DFS

            // IMPORTANT: set a tolerance so that we don't split hairs
            cplex.setParam( IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.001 );  // default: 1e-06
            cplex.setParam( IloCplex::Param::TimeLimit, WORSTCASE_MIP_TIMELIMIT );
            cplex.setParam( IloCplex::Param::ClockType, 2 );
            // apply cutoff_value (maximum worst-case value based on upper bound)
            // https://www.ibm.com/support/knowledgecenter/SSSA5P_20.1.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/CutLo.html
            // if(cutoff_value > 0) {
                // cplex.setParam( IloCplex::Param::MIP::Tolerances::LowerCutoff, cutoff_value );
            // }

            IloRangeArray con(env);
            IloNumVar wct(env);
            unsigned m = P_bar.size1() - 1;
            NumVarMatrix dev;
            if(first_call || true) {
                if(mip_model == 0) {  // Alt model
                    dev = populate_alt_model(model, con, wct, m, n, P_bar, P_hat, rob_params, pi, w, cutoff_value);
                } else {  // model == 1, Wilson model
                    dev = populate_wilson_model(model, con, wct, m, n, P_bar, P_hat, rob_params, pi, w, cutoff_value);
                }
            } else {
                // TODO implement model reoptimization (update constraints after change of z matrix
            }
            //cplex.exportModel("wct_worstcase.lp");

            // Solve the model.
            if(verbose) {
                std::cout << "Solving the worst-case model through CPLEX..." << std::endl;
                std::cout << "Budget type: " << rob_params.budget_type << std::endl;
                std::cout << "Permutation: " << pi << std::endl;
                std::cout << "Gamma: " << rob_params.budget_gamma << std::endl;
                if(cutoff_value.first > 0) {
                    std::cerr << "cutoff_value (lb): " << cutoff_value.first << std::endl;
                }
                if(cutoff_value.second > 0) {
                    std::cerr << "cutoff_value (ub): " << cutoff_value.second << std::endl;
                }
            }
            std::vector< std::vector<int> > dev_scenario;
            std::vector< std::vector<double> > dev_multiplier;
            double wct_value = 0.0;
            if ( cplex.solve() ) {
                // fetch the solution
                wct_value = cplex.getObjValue();
                for (int r = 0; r < m; r++) {
                    dev_scenario.push_back(std::vector<int>());
                    dev_multiplier.push_back(std::vector<double>());
                    for (int i = 0; i < n; i++) {
                        // dev[r][i] = 1, if job i deviates to its worst-case value on Machine r
                        double dev_ri = cplex.getValue(dev[r][i]);
                        if(dev_ri >= EPS) {
                            dev_scenario[r].push_back(i + 1);
                            dev_multiplier[r].push_back(dev_ri);
                        }
                    }
                    if(verbose) {
                        std::cout << "dev[" << r << "] = (" << dev_scenario[r] << " ; " << dev_multiplier[r] << ")\n";
                    }
                }
                double time_spent = TimeDateUtil::calculate_time_spent(timer);
                time_spent_worstcase_mip = time_spent;  // (double) cplex.getTime() - start;
                if(verbose) {
                    std::cout << "Objective: " << wct_value << std::endl;
                    std::cout << "CPLEX explored " << cplex.getNnodes() << " nodes.\n";
                    std::cout << "Time spent (s): " << time_spent_worstcase_mip << std::endl;
                }
            } else {
                if((cutoff_value.first == 0) && (cutoff_value.second == 0)) {
                    std::cerr << "[PFSP_WCT_Budget_WorstCase] ERROR : MILP INFEASIBLE !!!" << std::endl;
                    // Export model to file (DEBUG)
                    cplex.extract(model);
                    cplex.exportModel("debug_wct_worstcase.lp");
                    std::cerr << "[PFSP_WCT_Budget_WorstCase] MILP model file exported to debug_wct_worstcase.lp.\n";
                } else {
                    if(verbose) {
                        std::cout << "Leaf node worst-case infeasible due to cutoff value." << std::endl;
                    }
                    if(release_cplex) {  // release cplex model, variables and memory after solve
                        model.end();
                    }
                    PFSPBudgetScenario scenario(rob_params.budget_gamma, dev_scenario, dev_multiplier);
                    scenario.value = std::numeric_limits<double>::max();
                    return scenario;
                }
            }
            if(release_cplex) {  // release cplex model, variables and memory after solve
                model.end();
            }
            PFSPBudgetScenario scenario(rob_params.budget_gamma, dev_scenario, dev_multiplier);
            scenario.value = wct_value;
            return scenario;
        } catch (IloException &e) {
            std::cerr << "[PFSP_WCT_Budget_WorstCase] IloException: " << e << std::endl;
            cplex.end();
            model.end();
            env.end();
            throw;
        } catch (...) {
            std::cerr << "[PFSP_WCT_Budget_WorstCase] Unknown exception in CPLEX !\n";
            cplex.end();
            model.end();
            env.end();
            // ex.printStackTrace();
            throw;
        }
    }
}
