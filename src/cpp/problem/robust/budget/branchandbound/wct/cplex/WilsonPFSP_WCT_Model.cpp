#include "WilsonPFSP_WCT_Model.h"
#include "../../../../../deterministic/PFSProblem.h"

#include <glog/logging.h>
#include <numeric>

namespace robust {
    namespace hybrid {
        using namespace std;
        using namespace boost;
        // boost::numeric::ublas::matrix

        /*
        * PFSP Wilson Model
        * The MIP model will look like the following, where
        *   z[i][j] = 1, if job i is assigned to sequence position j, 0 otherwise
        *   b[r][j] = start (begin) time of job in sequence position j on machine r
        *
        * minimize B[m][n] + sum_i ( T[m][i] * z[i][n] )
        * s.t.
        *   sum_j z[i][j] = 1 for all i (assign every sequence position i a job j)
        *   sum_i z[i][j] = 1 for all j (assign every job j a sequence position i)
        *
        *   B[1][j] + sum_i ( T[1][i] * z[i][j] ) = B[1][j+1]            for all j in 1..n-1
        *   B[1][1] = 0
        *   B[r][1] + sum_i ( T[r][i] * z[i][1] ) = B[r+1][1]            for r in 1..m-1
        *   B[r][j] + sum_i ( T[r][i] * z[i][j] ) <= B[r+1][j]           for r in 1..m-1, j in 2..n
        *   B[r][j] + sum_i ( T[r][i] * z[i][j] ) <= B[r][j+1]           for r in 2..m, j in 1..n-1
        */
        NumVarMatrix WilsonPFSP_WCT_Model::populate_model(IloCplex &cplex, IloModel &model, IloRangeArray &c,
                                                IloNumVar &obj,
                                                const int &m, const int &n, const ublas::matrix<double> &P_bar,
                                                const ublas::matrix<double> &P_hat,
                                                const ublas::vector<double> &w,
                                                const std::vector<int> &initial_permutation,
                                                const std::vector< ublas::matrix<double> > &cut_list,
                                                const bool &relax) {
            IloEnv env = model.getEnv();
            std::vector<int> pi(n);
            std::iota (std::begin(pi), std::end(pi), 1); // Fill with 0, 1, ..., n.
            double bigM = calculate_bigM(m, n, P_bar, P_hat, pi, 3);
            // 1. First stage variables and constraints
            NumVarMatrix z(env, n);  // z is bin => here and now variable, first stage
            for (IloInt i = 0; i < n; i++) {
                if(relax) {
                    z[i] = IloNumVarArray(env, n, 0, 1);  // Float
                } else {
                    z[i] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // z is bin
                }
                std::stringstream sz;
                sz << "z(" << i << ")";
                z[i].setNames(sz.str().c_str());
            }
            for (IloInt i = 0; i < n; i++) {  // Assignment 1
                c.add(IloSum(z[i]) == 1); // single assignment, job i
            }
            for (IloInt j = 0; j < n; j++) {  // Assignment 2
                IloExpr expr1(env);
                for (IloInt i = 0; i < n; i++) {
                    expr1 += z[i][j];
                }
                c.add(expr1 == 1);
                expr1.end();
            }
            model.add(c);  // add all the constraints defined above
            // 2. Second stage variables and constraints
            for(int cut_count = 0; cut_count < cut_list.size(); cut_count++) {
                ublas::matrix<double> sep = cut_list[cut_count];  // the matrix containing current cut
                LOG(INFO) << "Processing cut: " << sep;
                cout << "Processing cut: " << sep << "\n";
                NumVarMatrix b(env, m);  // b >= 0
                // E[i] : start time of job i on the last machine (M_m)
                IloNumVarArray E(env, n, 0, IloInfinity); // E[N] >= 0
                for (IloInt r = 0; r < m; r++) {
                    b[r] = IloNumVarArray(env, n, 0, IloInfinity);
                    std::stringstream sb;
                    sb << "b_" << cut_count << "(" << r << ")";
                    b[r].setNames(sb.str().c_str());
                }
                for (IloInt i = 0; i < n; i++) {
                    std::stringstream sb;
                    sb << "E_" << cut_count << "(" << i << ")";
                    E[i].setName(sb.str().c_str());
                }
                // add the constraints
                // B[1][j] + sum_i ( (P_bar[1][i] + sep[1][i] * P_hat[1][i]) * z[i][j] ) = B[1][j+1]            for all j in 1..n-1
                for (IloInt j = 0; j < n - 1; j++) {
                    IloExpr expr1(env);
                    expr1 += b[0][j];
                    for (IloInt i = 0; i < n; i++) {
                        expr1 += (P_bar(1, i+1) + P_hat(1, i+1) * sep(1, i+1)) * z[i][j];
                    }
                    model.add(expr1 == b[0][j+1]);
                    expr1.end();
                }
                model.add(b[0][0] == 0);
                // B[r][1] + sum_i ( (P_bar[r][i] + P_hat[r][i] * sep[r][i]) * z[i][1] ) = B[r+1][1]            for r in 1..m-1
                for (IloInt r = 0; r < m - 1; r++) {
                    IloExpr expr1(env);
                    expr1 += b[r][0];
                    for (IloInt i = 0; i < n; i++) {
                        expr1 += (P_bar(r+1, i+1) + P_hat(r+1, i+1) * sep(r+1, i+1) ) * z[i][0];
                    }
                    model.add(expr1 == b[r + 1][0]);
                    expr1.end();
                }
                // B[r][j] + sum_i ( (P_bar[r][i] + P_hat[r][i] * sep[r][i] ) * z[i][j] ) <= B[r+1][j]           for r in 1..m-1, j in 2..n
                for (IloInt r = 0; r < m - 1; r++) {
                    for (IloInt j = 1; j < n; j++) {
                        IloExpr expr1(env);
                        expr1 += b[r][j];
                        for (IloInt i = 0; i < n; i++) {
                            expr1 += (P_bar(r+1, i+1) + P_hat(r+1, i+1) * sep(r+1, i+1)) * z[i][j];
                        }
                        model.add(expr1 <= b[r + 1][j]);
                        expr1.end();
                    }
                }
                // B[r][j] + sum_i ( (P_bar[r][i] + P_hat[r][i] * sep[r][i]) * z[i][j] ) <= B[r][j+1]           for r in 2..m, j in 1..n-1
                for (IloInt r = 1; r < m; r++) {
                    for (IloInt j = 0; j < n - 1; j++) {
                        IloExpr expr1(env);
                        expr1 += b[r][j];
                        for (IloInt i = 0; i < n; i++) {
                            expr1 += (P_bar(r+1, i+1) + P_hat(r+1, i+1) * sep(r+1, i+1)) * z[i][j];
                        }
                        model.add(expr1 <= b[r][j + 1]);
                        expr1.end();
                    }
                }
                // Disjuctive constraints used to determine the start time of job i on the last machine m
                for (IloInt i = 0; i < n; i++) {
                    for (IloInt j = 0; j < n; j++) {
                        // @constraint(model, E[i] >= B[m, j] - bigM * (1 - b1[i, j]) )
                        model.add(E[i] >= b[m - 1][j] - bigM * (1 - z[i][j]));
                        // @constraint(model, b1[i, j] <= Z[i, j])
                        //// model.add(b1[i][j] <= z[i][j]);
                    }
                }
                // add the objective function
                // wct == sum(w[i] * E[i] for i in Jobs) + sum(w[i] * sum( ( P_bar[m, i] + P_hat[m, i] * sep[m, i] ) * Z[i, j] for j in Jobs) for i in Jobs)
                IloExpr expr0(env);
                for (IloInt i = 0; i < n; i++) {
                    // sum(w[i] * E[i] for i in Jobs)
                    expr0 += w(i + 1) * E[i];
                    // sum(w[i] * sum( ( P_bar[m, i] + P_hat[m, i] * sep[m, i] ) * Z[i, j] for j in Jobs) for i in Jobs)
                    for (IloInt j = 0; j < n; j++) {
                        expr0 += w(i + 1) * (P_bar(m, i+1) + P_hat(m, i+1) * sep(m, i+1) ) * z[i][j];
                    }
                }
                model.add(obj >= expr0);
                expr0.end();
            }
            model.add(IloMinimize(env, obj));
            // set mip start using the given permutation
            set_warm_start_solution(cplex, model, n, initial_permutation, z);
            return z;
        }

        void WilsonPFSP_WCT_Model::fix_partial_permutation_in_model (IloModel &model, NumVarMatrix &z, IloRangeArray &c,
                                                                const std::vector<int> &seq_1,
                                                                const int &n) {
            IloEnv env = model.getEnv();
            for(int position_to_fix = 0; position_to_fix < seq_1.size(); position_to_fix++) {  // seq_1
                int job_to_fix = seq_1[position_to_fix] - 1;  // job number must be zero-based in the MIP model !
                for (int job = 0; job < n; job++) {  // for each job 'job'
                    if (job == job_to_fix) {  // 'job_to_fix' is in 'position_to_fix': force variable z[i][j] to 1
                        c.add(z[job][position_to_fix] == 1);
                    } else {  // no other job can occupy position 'position_to_fix': force variables z[i][j] to 0
                        c.add(z[job][position_to_fix] == 0);
                    }
                }
                for (int p = 0; p < n; p++) {  // for each position p
                    if (p != position_to_fix) {  // job 'job_to_fix' cannot occupy any other position other than position_to_fix: force variables to be 0
                        c.add(z[job_to_fix][p] == 0);
                    }
                }
            }
            model.add(c);
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
        ublas::matrix<double> WilsonPFSP_WCT_Model::calculate_cmax_given_scenario(int m, int n, const ublas::matrix<double> &P_bar,
                                                            const ublas::matrix<double> &P_hat, const std::vector<int> &pi,
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
                perm.push_back(pi[x]);
            }
            ublas::matrix<double> C = PFSProblem::makespan(perm, m, P_sum);
            return C;
        }

        /**
         * Calculates the BigM value to be used in one of the Rob-PFSP-WCT MILP models, according to the type:
         * - type == 0 : standard upper bound for BigM (naive);
         * - type == 2 : Calculate the Makespan (Cmax) for the worst-case scenario (on all machines r, all jobs j have their
         *                  processing times deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
         * - type == 3: Calculate the upper bound of the makespan considering the (m + n - 1) largest operations.
         * @param m number of machines
         * @param n number of jobs
         * @param P_bar matrix of nominal processing times
         * @param P_hat matrix of processing time variations
         * @return bigM value.
         */
        double WilsonPFSP_WCT_Model::calculate_bigM(int m, int n, const ublas::matrix<double> &P_bar,
                                                const ublas::matrix<double> &P_hat,
                                                const std::vector<int> &pi, const unsigned& bigM_type = 3) {
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
                return max_value;
            } else if(bigM_type == 2) {
                return calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario)(m,n) + 1;
            } else {  // bigM_type == 3
                std::vector<double> op_list;
                for (int r = 0; r < m; r++) {
                    for (int i = 0; i < n; i++) {
                        op_list.push_back(P_bar(r + 1, i + 1) + P_hat(r + 1, i + 1));
                    }
                }
                // Get the (m + n - 1) largest elements from op_list
                std::vector<double> largest(m + n - 1);
                std::partial_sort_copy(
                    std::begin(op_list), std::end(op_list),
                    std::begin(largest), std::end(largest),
                    std::greater<double>()
                );
                // Sum the (m + n - 1) largest elements from largest
                max_value = std::accumulate(largest.begin(), largest.end(), 0.0);
                return max_value + 1;
            }
        }

        /**
         * Set MIP start solution using the Upper Bound Solution given by permutation perm.
         */
        void WilsonPFSP_WCT_Model::set_warm_start_solution(IloCplex &cplex, IloModel &model, const int &n,
                const std::vector<int>& perm, NumVarMatrix &z) {
            IloEnv env = model.getEnv();
            IloNumVarArray startVar(env);
            IloNumArray startVal(env);
            for (int j = 0; j < n; ++j) {  // z[i][j] = 1, if job i is assigned to sequence position j, 0 otherwise
                int job_i = perm[j] - 1;  // all indices in PFSP solutions are 1-based !
                startVar.add(z[job_i][j]);
                startVal.add(1);
                for (int i = 0; i < n; ++i) {
                    if(i != job_i) {
                        startVar.add(z[i][j]);
                        startVal.add(0);
                    }
                }
            }
            cplex.addMIPStart(startVar, startVal);
            startVal.end();
            startVar.end();
            LOG(INFO) << "Wilson PFSP-WCT MP warm start solution set.";
        }

        std::vector<int> WilsonPFSP_WCT_Model::getSolutionAsPermutation(IloCplex &cplex, IloModel &model,
                const int &n, NumVarMatrix &z) {
            Permutation permutation;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    // z[i][j] = 1, if job i is assigned to sequence position j, 0 otherwise
                    if(cplex.getValue(z[i][j]) >= 0.95)
                        permutation.push_back(i + 1);
                }
            }
            return permutation;
        }

    } /* namespace hybrid */
} /* namespace robust */
