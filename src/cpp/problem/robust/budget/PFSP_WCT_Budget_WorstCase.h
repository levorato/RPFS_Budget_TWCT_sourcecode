//
// Created by mlevorato on 5/1/20.
//

#ifndef PFSP_PFSP_WCT_BUDGET_WORSTCASE_H
#define PFSP_PFSP_WCT_BUDGET_WORSTCASE_H

#include "../../deterministic/PFSProblem.h"
#include "../../PFSSolution.h"
#include "../../PFSP_Parameters.h"
#include "../../metaheuristic/common/Job.h"

#include "PFSPBudgetScenario.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define MAX_MACHINES 502

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;
    using namespace problem::pfsp;

    class JobDeviationList {
    public:
        std::vector<std::vector<int>> list;

        JobDeviationList() : list(MAX_MACHINES) {
        }
        JobDeviationList(const int &m) : list(m) {
        }
        JobDeviationList(const JobDeviationList& o) : list(o.list) {
        }
        std::vector<int> &operator[](int i) {
            if( i >= list.size() ) {
                std::cerr << "Index out of bounds" << std::endl;
                // return first element.
                return list[0];
            }
            return list[i];
        }
        void add(int r, int job_number) {
            list[r].push_back(job_number);
        }
        unsigned size(int r) {
            if(list.size() == 0) {
                return 0;
            }
            return list[r].size();
        }
        std::string str() {
            stringstream ss;
            ss << "[ ";
            for(std::vector<int> l : list) {
                ss << "[";
                for(int x : l) {
                    ss << " " << x;
                }
                ss << " ], ";
            }
            ss << " ]";
            return ss.str();
        }
        virtual ~JobDeviationList() {}
    };

    /**
     * Portions of the class code were ported from Julia, file src/julia/wct/robust_pfsp_budget_worstcase_wct_v3.jl
     */
    class PFSP_WCT_Budget_WorstCase {
    public:
        // stores statistics about the branch-and-bound process
        static double time_spent_worstcase_mip;

        /**
         * Procedure Worst Case - Continuous time intervals.
         * Calculates the worst-case makespan given a sequence / permutation pi and
         * budget parameter Gamma. Select between the machine budget type and global budget type.
         * @param nUsedJobs is the number of jobs in the partial permutation pi, where nUsedJobs <= n
         */
        static double worst_case_wct_dp(const unsigned &nUsedJobs,
                                         const boost::numeric::ublas::matrix<double> &P_bar,
                                         const boost::numeric::ublas::matrix<double> &P_hat,
                                         const Robust_Parameters &rob_params,
                                         const std::vector<Job> &pi,
                                         const boost::numeric::ublas::vector<double> &w) {
            if(rob_params.budget_type == "machine") {
                return worst_case_wct_dp_machine_budget(nUsedJobs, P_bar, P_hat, rob_params, pi, w);
            } else if (rob_params.budget_type == "global") {
                return worst_case_wct_dp_global_budget(nUsedJobs, P_bar, P_hat, rob_params, pi, w);
            } else {
                throw std::invalid_argument( "Invalid budget type!");
            }
        }

        /**
         * Procedure Worst Case Weighted Sum of Completion Times (WCT) - Continuous time intervals.
         * Calculates the worst-case WCT given a sequence / permutation pi and
         * budget parameters T1 and T2. Uses the EXACT approach based on a MILP model.
         * @param n number of jobs
         * @param P_bar matrix of nominal processing times
         * @param P_hat matrix of processing time variations
         * @param rob_params Robust PFSP problem parameters
         * @param pi the job permutation
         * @param w vector of job weights
         * @param first_call
         * @param release_cplex
         * @param verbose
         * @return
         */
        static PFSPBudgetScenario worst_case_wct_mip(const unsigned &n,
                                const boost::numeric::ublas::matrix<double> &P_bar,
                                const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                const std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w,
                                const bool& first_call, const bool& release_cplex = false, const unsigned& mip_model = 0,
                                const bool& verbose = false,
                                const std::pair<double,double>& cutoff_value = std::pair<double,double>());

        static void release_cplex_model();

        /**
         * Procedure Worst Case Weighted Sum of Completion Times (WCT) - Continuous time intervals.
         * Approximate DP for worst-case WCT.
         * Calculates the worst-case WCT given a sequence / permutation pi and
         * budget parameters T1 and T2. Uses a dynamic programming based on a recursion.
         * WARNING: For the RobPFSP-WCT problem, this dynamic programming approach DOES NOT yield the exact value for the worst case.
         * Attention: original formulas (in Julia) were made for 1-based arrays (e.g. pi[k] with k = 1, ..., n).
         * Complexity: O(n ^ 3)
         */
        static double worst_case_wct_dp_machine_budget(const unsigned &n,
                                         const boost::numeric::ublas::matrix<double> &P_bar,
                                         const boost::numeric::ublas::matrix<double> &P_hat,
                                         const Robust_Parameters &rob_params,
                                         const std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w) {
            unsigned k = 0, y1 = 0, y2 = 0, q = 0;
            unsigned m = P_bar.size1() - 1;
            //std::cout << "*** The number of machines is " << m << "\n";
            // TODO FIXME In case this procedure is invoked from a metaheuristic to find the worst case of a partial
            // solution pi where size(pi) < n, should we still apply T1 and T2 as a percentage of the sequence size
            // or should we make T1 = min(T1, instance.n) ?
            std::vector<unsigned> T(n, 0);
            for(y1 = 0; y1 < n; y1++) {
                T[y1] = (unsigned)ceil((n * rob_params.budget_gamma[y1]) / 100.0);
            }
            boost::multi_array<double, 3> C(boost::extents[m + 1][n + 2][n + 2]);
            std::fill_n(C.data(), C.num_elements(), 0);
            // LOG(INFO) << "[DP] Calculating worst-case Cmax given a permutation pi = $(pi), m = $(m)"
            //   << "and n = $(n) and budget parameters T = $((T1, T2))...";
            // Filling matrix C[1]
            for(y1 = 0; y1 <= n; y1++) {  // T1
                C[1][0][y1] = 0;
            }
            for(k = 1; k <= n; k++) {
                for (y1 = 0; y1 <= T[(1) - 1]; y1++) {  // for y1 in 0:T[1]
                    C[1][k][y1] = std::max(C[1][k][y1], C[1][k - 1][y1] + P_bar(1, pi[(k) - 1].id));
                    C[1][k][y1 + 1] = std::max(C[1][k][y1 + 1], C[1][k - 1][y1] + P_hat(1, pi[(k) - 1].id) + P_bar(1, pi[(k) - 1].id));
                }
            }
            // Filling matrix C[q]
            for(q = 2; q <= m; q++) {  // For each machine [q] in 2:m
                y1 = T[(q - 1) - 1];  // y1 = T[q - 1]
                for (k = 1; k <= n; k++) {  // for k in 1:n
                    for (y2 = 0; y2 <= T[(q) - 1]; y2++) {  // for y2 in 0:T[q]
                        C[q][k][y2] = std::max(C[q][k][y2], std::max(C[q - 1][k][y1], C[q][k - 1][y2])
                                                            + P_bar(q, pi[(k) - 1].id));
                        C[q][k][y2 + 1] = std::max(C[q][k][y2 + 1], std::max(C[q - 1][k][y1], C[q][k - 1][y2])
                                                            + P_bar(q, pi[(k) - 1].id) + P_hat(q, pi[(k)- 1].id));
                    }
                }
            }
            double worst_cmax = C[m][n][T[(m) - 1]];
            double wct = 0.0;
            for(int j = 1; j <= n; j++) {  // C[m, j] is the final completion time of job j
                k = pi[j - 1].id;
                wct += C[m][j][T[(m) - 1]] * w[k];  // C[m, j] * w[k]
            }
            // LOG(INFO) << "[DP] The worst-case Cmax is " << worst_cmax << ".\n";
            return wct;
        }

        /**
         * Procedure Worst Case - Continuous time intervals - Global Budget.
         * Approximate DP for worst-case WCT.
         * Based on global-budget type, calculates the worst-case WCT given a sequence / permutation pi and
         * budget parameter Gamma. Uses a dynamic programming based on a recursion.
         * Attention: original formulas (in Julia) were made for 1-based arrays (e.g. pi[k] with k = 1, ..., n).
         * Complexity: O(n ^ 3)
         * @param nUsedJobs is the number of jobs in the partial permutation pi, where nUsedJobs <= n
         */
        static double worst_case_wct_dp_global_budget(const unsigned &nUsedJobs,
                                                      const boost::numeric::ublas::matrix<double> &P_bar,
                                                      const boost::numeric::ublas::matrix<double> &P_hat,
                                                      const Robust_Parameters &rob_params,
                                                      const std::vector<Job> &pi,
                                                      const boost::numeric::ublas::vector<double> &w) {
            unsigned k = 0, y = 0, q = 0;
            unsigned m = P_bar.size1() - 1;
            unsigned n = P_bar.size2() - 1;  // n is the total number of jobs of the problem instance
            // IMPORTANT: In case this procedure is invoked from a metaheuristic to find the worst case of a partial
            // solution pi where size(pi) < n, then we should apply T[i] as T[i] = min(T[i], nUsedJobs)
            unsigned Gamma = (unsigned)ceil((n * m * rob_params.budget_gamma[0]) / 100.0);
            int max_gamma = n * m;
            boost::multi_array<double, 3> C(boost::extents[m + 1][n + 2][max_gamma + 2]);
            std::fill_n(C.data(), C.num_elements(), 0);
            // Filling matrix C[1]
            for(y = 0; y <= n; y++) {  // Gamma
                C[1][0][y] = 0;
            }
            for(k = 1; k <= n; k++) {
                for(y = 0; y <= max_gamma; y++) {  // for y in 0:(T+1)
                    // Ca[1, k, y, 0] = max(Ca[1, k, y, 0], Ca[1, k - 1, y, 0] + P_bar[1, pi[k]])
                    C[1][k][y] = std::max(C[1][k][y], C[1][k - 1][y] + P_bar(1, pi[(k) - 1].id));
                    // Ca[1, k, y + 1, 0] = max(Ca[1, k, y + 1, 0], Ca[1, k - 1, y, 0] + P_hat[1, pi[k]] + P_bar[1, pi[k]])
                    C[1][k][y + 1] = std::max(C[1][k][y + 1], C[1][k - 1][y] + P_hat(1, pi[(k) - 1].id) + P_bar(1, pi[(k) - 1].id));
                }
            }
            // Filling matrix C[q]
            for(q = 2; q <= m; q++) {  // For each machine [q] in 2:m
                for(k = 1; k <= n; k++) {  // for k in 1:nUsedJobs
                    for(y = 0; y <= max_gamma; y++) {  // for y in 0:Gamma
                        // Ca[q, k, y, 0] = max(Ca[q, k, y, 0], max(Ca[q - 1, k, y, 0], Ca[q, k - 1, y, 0]) + P_bar[q, pi[k]])
                        C[q][k][y] = std::max(C[q][k][y], std::max(C[q - 1][k][y], C[q][k - 1][y])
                                                          + P_bar(q, pi[(k) - 1].id));
                        // Ca[q, k, y + 1, 0] = max(Ca[q, k, y + 1, 0], max(Ca[q - 1, k, y, 0], Ca[q, k - 1, y, 0]) + P_bar[q, pi[k]] + P_hat[q, pi[k]])
                        C[q][k][y + 1] = std::max(C[q][k][y + 1], std::max(C[q - 1][k][y], C[q][k - 1][y])
                                                                  + P_bar(q, pi[(k) - 1].id) + P_hat(q, pi[(k)- 1].id));
                    }
                }
            }
            double worst_wct = C[m][nUsedJobs][Gamma];
            // LOG(INFO) << "[DP] The worst-case WCT is " << worst_wct << ".\nUsedJobs";
            return worst_wct;
        }

        // WARN : THE CODE BELOW IS DEPRECATED !

        static int calculate_max_time_deviation(const unsigned &m, const unsigned &n,
                                                const boost::numeric::ublas::matrix<double> &P_bar,
                                                const boost::numeric::ublas::matrix<double> &P_hat,
                                                const std::vector<int> &T) {
            int max_time_deviation = 0;
            int max_gamma = 0;
            for(int r = 0; r < m; r++) {
                int gamma = int(ceil(T[r]));
                max_gamma = std::max(max_gamma, gamma);
            }
            for(int r = 1; r <= m; r++) {
                int max_dev_r = 0;
                for(int i = 1; i <= n; i++) {  // min(max_gamma + 1, m)
                    max_dev_r += int(ceil(P_hat(r, i)));
                }
                max_dev_r = int(ceil(max_dev_r)) + 1;
                max_time_deviation = std::max(max_time_deviation, max_dev_r) + 1;
            }
            return max_time_deviation;
        }

        /**
            # f(pi, r, k, C, deviated_jobs, m, n, P_bar, P_hat, w)
            # pi                : the sequence / permutation pi
            # r                 : the current machine index r in 1:m
            # k                 : until which sequence position should the objective function be calculated ?
            # C                 : current completion time matrix
            # deviated_jobs     : the list of deviated jobs in each machine i in 1:m
            # m                 : the number of machines m
            # n                 : the number of jobs n
            # P_har             : the nominal processing times matrix
            # P_hat             : the processing time deviations matrix
            # w                 : the weight array for the weighted completion time calculation
         */
        static double fp(const std::vector<Job> &pi, int r, int k, boost::numeric::ublas::matrix<double> &C,
                JobDeviationList &deviated_jobs, int m, int n,
                const boost::numeric::ublas::matrix<double> &P_bar,
                const boost::numeric::ublas::matrix<double> &P_hat, const boost::numeric::ublas::vector<double> &w,
                bool verbose = false) {
            // println("Deviated jobs : $(deviated_jobs)")
            // indices r and k (machine and sequence position) are already iterated inside the DP procedure
            unsigned x = pi[k - 1].id;
            if(r == 1) {
                C(1, k) = (k >= 2 ? C(1, k - 1) : 0) + P_bar(1, x);
                if (std::find(deviated_jobs[1].begin(), deviated_jobs[1].end(), x) != deviated_jobs[1].end()) {  // (x in deviated_jobs[1]) {
                    if(verbose)
                        std::cout << "Deviating job " << x << ": " << (C(1, k)) << " + " << (P_hat(1, x)) << "\n";
                    C(1, k) += P_hat(1, x);
                }
            } else {
                C(r, k) = std::max((k >= 2 ? C(r, k - 1) : 0), C(r - 1, k)) + P_bar(r, x);
                if (std::find(deviated_jobs[r].begin(), deviated_jobs[r].end(), x) != deviated_jobs[r].end()) {  // (x in deviated_jobs[r]) {
                    if(verbose)
                        std::cout << "Deviating job " << x << ": " << (C(r, k)) << " + " << (P_hat(r, x)) << "\n";
                    C(r, k) += P_hat(r, x);
                }
            }
            // New: propagate the makespan calculation
            for(int j = (k+1); j <= n; j++) {
                unsigned x = pi[j - 1].id;
                C(r, j) = std::max(C(r, j - 1), (r >= 2 ? C(r - 1, j) : 0)) + P_bar(r, x);
            }
            for(int j = 1; j <= n; j++) {
                unsigned x = pi[j - 1].id;
                for (int i = (r + 1); i <= m; i++) {
                    C(i, j) = std::max((j >= 2 ? C(i, j - 1) : 0), C(i - 1, j)) + P_bar(i, x);
                }
            }
            double sum_wct = 0;  // Calculate the partial completion time for the current machine, until sequence index k
            for(int j = 1; j <= n; j++) {  // sequence_pos  # C[r, j] is the final completion time of job on position j on machine r
                unsigned k = pi[j - 1].id;   // Job k is in position j of the permutation
                sum_wct += C(m, j) * w(k);
            }
            // if(verbose)
            //     std::cout << "[f!] wct = " << sum_wct << "; C = " << C << "\n";
            return sum_wct;
        }

        /**
            # ff(1, ptk, pi, 1, n, P_bar, P_hat, w)
            # pi                : the sequence / permutation pi
            # r                 : the machine index r in 1:m
            # sequence_pos      : until which sequence position should the objective function be calculated ?
            # current_deviation : current processing time deviation for machine r
            # deviated_jobs     : the list of deviated jobs in each machine i in 1:m
            # n                 : the number of jobs n
            # P_har             : the nominal processing times matrix
            # P_hat             : the processing time deviations matrix
            # w                 : the weight array for the weighted completion time calculation
        */
        static std::pair<double,boost::numeric::ublas::matrix<double>> ff(const std::vector<Job> &pi, int r, int sequence_pos,
                  int current_deviation,
                  JobDeviationList &deviated_jobs, int m, int n,
                  const boost::numeric::ublas::matrix<double> &P_bar,
                  const boost::numeric::ublas::matrix<double> &P_hat, const boost::numeric::ublas::vector<double> &w,
                  bool verbose = false) {
            // println("Deviated jobs : $(deviated_jobs)")
            r = m;
            boost::numeric::ublas::matrix<double> C = ublas::zero_matrix<double>(m + 1, n + 1);
            unsigned k = pi[0].id;
            C(1, 1) = P_bar(1, k);
            if(std::find(deviated_jobs[1].begin(), deviated_jobs[1].end(), k) != deviated_jobs[1].end()) {  // if(k in deviated_jobs[1]) {
                C(1, 1) += P_hat(1, k);
            }
            for(int i = 2; i <= r; i++) {
                C(i, 1) = C(i - 1, 1) + P_bar(i, k);
                if(std::find(deviated_jobs[i].begin(), deviated_jobs[i].end(), k) != deviated_jobs[i].end()) {  // if(k in deviated_jobs[i]) {
                   C(i, 1) += P_hat(i, k);
                }
            }
            for(int j = 2; j <= n; j++) {
                k = pi[j - 1].id;
                C(1, j) = C(1, j - 1) + P_bar(1, k);
                if(std::find(deviated_jobs[1].begin(), deviated_jobs[1].end(), k) != deviated_jobs[1].end()) {  // if(k in deviated_jobs[1]) {
                   C(1, j) += P_hat(1, k);
                }
                for(int i = 2; i <= r; i++) {
                    C(i, j) = std::max(C(i, j - 1), C(i - 1, j)) + P_bar(i, k);
                    if(std::find(deviated_jobs[i].begin(), deviated_jobs[i].end(), k) != deviated_jobs[i].end()) {  // if(k in deviated_jobs[i]) {
                        C(i, j) += P_hat(i, k);
                    }
                }
            }
            double sum_wct = 0;
            for(int j = 1; j <= n; j++) {  // sequence_pos  # C[r, j] is the final completion time of job on position j on machine r
                unsigned k = pi[j - 1].id;   // Job k is in position j of the permutation
                sum_wct += C(m, j) * w(k);
            }
            if(sum_wct < 0)
                sum_wct = 0;
            if(verbose)
                std::cout << "[ff] wct = " << sum_wct << "; C = " << C << "\n";
            return std::make_pair(sum_wct, C);
        }

        static double worst_case_wct_dp_v2(const unsigned &n,
                                         const boost::numeric::ublas::matrix<double> &P_bar,
                                         const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                         std::vector<Job> &pi, const boost::numeric::ublas::vector<double> &w) {

            // function worst_case_wct_dp_m_machines_2(pi, m, n, Γ, P_bar, P_hat, w, verbose = false)
            unsigned m = P_bar.size1() - 1;
            bool verbose = false;
            std::vector<int> T(n, 0);
            stringstream ss_T, ss_pi;
            // std::vector<int> pi_int{ 1, 5, 10, 2, 3, 7, 4, 6, 9, 8 };
            if(rob_params.budget_gamma.size() > 1) {
                for(int y1 = 0; y1 < m; y1++) {
                    T[y1] = (unsigned)ceil((n * rob_params.budget_gamma[y1]) / 100.0);
                    //ss_T << T[y1] << " ";
                }
            } else {
                for(int y1 = 0; y1 < m; y1++) {
                    T[y1] = (unsigned)ceil((n * rob_params.budget_gamma[0]) / 100.0);
                    //ss_T << T[y1] << " ";
                }
            }
            /*
            for(int y1 = 0; y1 < n; y1++) {
                // pi[y1].setId(pi_int[y1]);
                ss_pi << pi[y1].id << " ";
            } */
            // calculate max_time_deviation for this solution
            int max_time_deviation = calculate_max_time_deviation(m, n, P_bar, P_hat, T);
            int max_gamma = 0;
            // sort process_time_deviation and pick the gamma bigger ones
            for(int r = 0; r < m; r++) {
                max_gamma = std::max(max_gamma, T[r]);
            }
            if(verbose)
                std::cout << "[WCT_DP] Gamma = [ " << ss_T.str() << "]; pi = [ " << ss_pi.str() << "]; max_gamma = "
                            << max_gamma << "; max_time_deviation = " << max_time_deviation << "\n";
            // Create dynamic programming matrix to store results
            // First index : which machine r in 1:m
            // Second index : which position k in 1:n
            // Third index : number of deviated processing times g in 1:max_gamma+1
            // Forth index : accumulated processing time deviation (discretized) d in 1:max_time_deviation+1
            boost::multi_array<double, 4> dynmatrix(boost::extents[m + 1][n + 2][max_gamma+2][max_time_deviation+2]);
            boost::multi_array<JobDeviationList, 4> deviated_jobs(boost::extents[m + 1][n + 2][max_gamma+2][max_time_deviation+2]);
            boost::multi_array<boost::numeric::ublas::matrix<double>, 4> C(boost::extents[m + 1][n + 2][max_gamma+2][max_time_deviation+2]);
            std::fill_n(dynmatrix.data(), dynmatrix.num_elements(), 0);
            JobDeviationList emptyJobDeviationList(m + 1);
            std::fill_n(deviated_jobs.data(), deviated_jobs.num_elements(), emptyJobDeviationList);
            ublas::zero_matrix<double> my_zero_matrix(m + 1, n + 1);
            std::fill_n(C.data(), C.num_elements(), my_zero_matrix);

            int gamma = max_gamma;
            boost::numeric::ublas::vector<double> max_wct(ublas::zero_vector<double>(m + 1));
            std::vector<JobDeviationList> max_delayed_jobs(m + 1);
            std::vector<std::vector<JobDeviationList>> max_delayed_jobs_all(m + 1);
            std::vector<JobDeviationList> delayed_jobs_all(m + 1);

            // Initialize the matrices for the first position (k = 1)
            std::pair<double,boost::numeric::ublas::matrix<double>> ret_0 = ff(pi, m, 1, 0, emptyJobDeviationList, m, n, P_bar, P_hat, w, verbose);
            double wct_0 = ret_0.first;
            boost::numeric::ublas::matrix<double> C_mtx(ret_0.second);
            for(int r = 1; r <= m; r++) {    // For each machine r
                for(int g = 1; g <= (gamma + 1); g++) {  // For each gamma value g
                    for(int d = 1; d <= max_time_deviation + 1; d++) {     // For each accumulated proc time deviation
                        dynmatrix[r][1][g][d] = wct_0;
                        deviated_jobs[r][1][g][d] = JobDeviationList(emptyJobDeviationList);
                        C[r][1][g][d] = boost::numeric::ublas::matrix<double>(C_mtx);
                        dynmatrix[r][1][g][d] = fp(pi, r, 1, C[r][1][g][d], deviated_jobs[r][1][g][d], m, n, P_bar, P_hat, w, verbose);
                    }
                }
            }
            // **************  RULES FOR MACHINE r = 1 *******************************
            // calculate processing time of job in position 1
            unsigned jobk = pi[1 - 1].id;
            double ptk = P_hat(1, jobk);
            deviated_jobs[1][1][1][int(ceil(ptk)) + 1] = JobDeviationList(emptyJobDeviationList);
            dynmatrix[1][1][1][int(ceil(ptk) + 1)] = fp(pi, 1, 1, C[1][1][1][int(ceil(ptk) + 1)], deviated_jobs[1][1][1][int(ceil(ptk) + 1)], m, n, P_bar, P_hat, w, verbose);
            for(int g = 2; g <= (gamma + 1); g++) {
                deviated_jobs[1][1][g][int(ceil(ptk)) + 1] = JobDeviationList(emptyJobDeviationList);
                deviated_jobs[1][1][g][int(ceil(ptk)) + 1].add(1, jobk);
                //std::cout << "temp(r = 1, k = 1) = " << deviated_jobs[1][1][g][int(ceil(ptk)) + 1].str() << "\n";
                std::pair<double,boost::numeric::ublas::matrix<double>> ret_t  = ff(pi, 1, 1, ptk, deviated_jobs[1][1][g][int(ceil(ptk)) + 1], m, n, P_bar, P_hat, w, verbose);
                dynmatrix[1][1][g][int(ceil(ptk)) + 1] = ret_t.first;
                C[1][1][g][int(ceil(ptk) + 1)] = ret_t.second;
            }

            deviated_jobs[1][1][1][1] = JobDeviationList(emptyJobDeviationList);
            C[1][1][1][1] = ublas::zero_matrix<double>(m + 1, n + 1);
            std::pair<double,boost::numeric::ublas::matrix<double>> ret_t = ff(pi, 1, 1, 0, deviated_jobs[1][1][1][1], m, n, P_bar, P_hat, w, verbose);
            dynmatrix[1][1][1][1] = ret_t.first;
            C[1][1][1][1] = ret_t.second;
            for(int k = 2; k <= n; k++) {
                deviated_jobs[1][k][1][1] = JobDeviationList(deviated_jobs[1][k - 1][1][1]);
                C[1][k][1][1] = boost::numeric::ublas::matrix<double>(C[1][k - 1][1][1]);
                dynmatrix[1][k][1][1] = fp(pi, 1, k, C[1][k][1][1], deviated_jobs[1][k][1][1], m, n, P_bar, P_hat, w, verbose);
                for(int d = 2; d <= max_time_deviation + 1; d++) {
                    deviated_jobs[1][k][1][d] = JobDeviationList(emptyJobDeviationList);
                    C[1][k][1][d] = boost::numeric::ublas::matrix<double>(C[1][k][1][d - 1]);
                    dynmatrix[1][k][1][d] = fp(pi, 1, k, C[1][k][1][d], deviated_jobs[1][k][1][d], m, n, P_bar, P_hat, w, verbose);
                }
                // EXTRA...
                unsigned jobk = pi[k - 1].id;
                double ptk = P_hat(1, jobk);
                for(int g = 2; g <= (gamma + 1); g++) {
                    deviated_jobs[1][k][g][int(ceil(ptk)) + 1] = JobDeviationList(emptyJobDeviationList);
                    deviated_jobs[1][k][g][int(ceil(ptk)) + 1].add(1, jobk);
                    //std::cout << "temp(r = 1, k) = " << deviated_jobs[1][k][g][int(ceil(ptk)) + 1].str() << "\n";
                    std::pair<double,boost::numeric::ublas::matrix<double>> ret_t  = ff(pi, 1, 1, ptk, deviated_jobs[1][k][g][int(ceil(ptk)) + 1], m, n, P_bar, P_hat, w, verbose);
                    dynmatrix[1][k][g][int(ceil(ptk)) + 1] = ret_t.first;
                    C[1][k][g][int(ceil(ptk) + 1)] = ret_t.second;
                }
            }
            boost::numeric::ublas::matrix<double> max_C_prev_machine = ublas::zero_matrix<double>(m + 1, n + 1);
            for(int r = 1; r <= m; r++) {
                // std::cout << "*** r = " << r << "\n";
                std::vector<JobDeviationList> max_delayed_jobs_all_combinations;
                if (r == 1) {
                    max_delayed_jobs_all_combinations.push_back(JobDeviationList(emptyJobDeviationList));
                } else {  // if(r >= 2) {
                    max_delayed_jobs_all_combinations = max_delayed_jobs_all[r - 1];
                }
                max_wct(r) = 0;
                max_delayed_jobs[r] = JobDeviationList(emptyJobDeviationList);
                max_delayed_jobs_all[r] = std::vector<JobDeviationList>();  // JobDeviationList(emptyJobDeviationList);  // = [ ];

                for(JobDeviationList max_delayed_jobs_prev_machine : max_delayed_jobs_all_combinations) {
                    if(r >= 2) {
                        // Machine r >= 2, k = 1 : no deviation (base case)
                        deviated_jobs[r][1][1][1] = JobDeviationList(max_delayed_jobs_prev_machine);
                        C[r][1][1][1] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);
                        dynmatrix[r][1][1][1] = fp(pi, r, 1, C[r][1][1][1], deviated_jobs[r][1][1][1], m, n, P_bar, P_hat, w, verbose);
                        // calculate processing time of job in position 1
                        jobk = pi[1 - 1].id;
                        ptk = P_hat(r, jobk);
                        deviated_jobs[r][1][1][int(ceil(ptk)) + 1] = JobDeviationList(max_delayed_jobs_prev_machine);
                        C[r][1][1][int(ceil(ptk)) + 1] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);
                        dynmatrix[r][1][1][int(ceil(ptk)) + 1] = fp(pi, r, 1, C[r][1][1][int(ceil(ptk)) + 1], deviated_jobs[r][1][1][int(ceil(ptk)) + 1], m, n, P_bar, P_hat, w, verbose);
                        for(int g = 1; g <= (gamma + 1); g++) {
                            for(int d = 1; d <= max_time_deviation + 1; d++) {
                                C[r][1][g][d] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);
                                deviated_jobs[r][1][g][d] = JobDeviationList(max_delayed_jobs_prev_machine);
                            }
                        }
                        for(int g = 2; g <= (gamma + 1); g++) {
                            deviated_jobs[r][1][g][int(ceil(ptk)) + 1] = JobDeviationList(max_delayed_jobs_prev_machine);
                            JobDeviationList temp = deviated_jobs[r][1][g][int(ceil(ptk)) + 1];
                            temp.add(r, jobk);
                            deviated_jobs[r][1][g][int(ceil(ptk)) + 1] = temp;
                            // std::cout << "temp(k = 1) = " << deviated_jobs[r][1][g][int(ceil(ptk)) + 1].str() << "\n";
                            C[r][1][g][int(ceil(ptk)) + 1] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);
                            dynmatrix[r][1][g][int(ceil(ptk)) + 1] = fp(pi, r, 1, C[r][1][g][int(ceil(ptk)) + 1], deviated_jobs[r][1][g][int(ceil(ptk)) + 1], m, n, P_bar, P_hat, w, verbose);
                        }

                        deviated_jobs[r][1][1][1] = JobDeviationList(max_delayed_jobs_prev_machine);
                        C[r][1][1][1] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);
                        dynmatrix[r][1][1][1] = fp(pi, r, 1, C[r][1][1][1], deviated_jobs[r][1][1][1], m, n, P_bar, P_hat, w, verbose);
                        for(int k = 2; k <= n; k++) {
                            // Machine r >= 2, position k : no deviation (base case)
                            deviated_jobs[r][k][1][1] = JobDeviationList(deviated_jobs[r][k - 1][1][1]);
                            C[r][k][1][1] = boost::numeric::ublas::matrix<double>(C[r][k - 1][1][1]);
                            dynmatrix[r][k][1][1] = fp(pi, r, k, C[r][k][1][1], deviated_jobs[r][k][1][1], m, n, P_bar, P_hat, w, verbose);
                            for(int d = 2; d <= max_time_deviation + 1; d++) {
                                // FIXME verificar linha abaixo!!!
                                deviated_jobs[r][k][1][d] = JobDeviationList(max_delayed_jobs_prev_machine);   //  copy(deviated_jobs[r - 1, k, 1, d])
                                C[r][k][1][d] = boost::numeric::ublas::matrix<double>(max_C_prev_machine);  // #  deepcopy(C[r - 1, k, 1, d])
                                dynmatrix[r][k][1][d] = fp(pi, r, k, C[r][k][1][d], deviated_jobs[r][k][1][d], m, n, P_bar, P_hat, w, verbose);
                            }
                        }
                    }
                    for(int k = 2; k <= n; k++) {
                        // std::cout << "*** k = " << k << "\n";
                        for(int g = 2; g <= (gamma + 1); g++) {
                            for(int d = 1; d <= max_time_deviation + 1; d++) {
                                // calculate processing time of job in position k
                                jobk = pi[k - 1].id;
                                ptk = int(ceil(P_hat(r, jobk)));
                                if(d - ptk - 1 < 0) {
                                    deviated_jobs[r][k][g][d] = JobDeviationList(deviated_jobs[r][k - 1][g][d]);
                                    C[r][k][g][d] = boost::numeric::ublas::matrix<double>(C[r][k - 1][g][d]);
                                    dynmatrix[r][k][g][d] = fp(pi, r, k, C[r][k][g][d], deviated_jobs[r][k][g][d], m, n, P_bar, P_hat, w, verbose);
                                } else {
                                    deviated_jobs[r][k][g][d] = JobDeviationList(deviated_jobs[r][k - 1][g - 1][d - ptk]);
                                    JobDeviationList temp = deviated_jobs[r][k][g][d];
                                    temp.add(r, jobk);
                                    deviated_jobs[r][k][g][d] = temp;
                                    // std::cout << "temp(k) = " << deviated_jobs[r][k][g][d].str() << "\n";
                                    // deviated_jobs[r][k][g][d][r].push_back(jobk);
                                    C[r][k][g][d] = boost::numeric::ublas::matrix<double>(C[r][k - 1][g - 1][d - ptk]);
                                    dynmatrix[r][k][g][d] = fp(pi, r, k, C[r][k][g][d], deviated_jobs[r][k][g][d], m, n, P_bar, P_hat, w, verbose);
                                    if (dynmatrix[r][k - 1][g][d] > dynmatrix[r][k][g][d]) {
                                        deviated_jobs[r][k][g][d] = JobDeviationList(deviated_jobs[r][k - 1][g][d]);
                                        C[r][k][g][d] = boost::numeric::ublas::matrix<double>(C[r][k - 1][g][d]);
                                        dynmatrix[r][k][g][d] = fp(pi, r, k, C[r][k][g][d], deviated_jobs[r][k][g][d], m, n, P_bar, P_hat, w, verbose);
                                    }
                                }
                            }
                        }
                    }
                }
                for(int d = max_time_deviation + 1; d >= 1; d--) {
                    if(deviated_jobs[r][n][int(ceil(T[(r) - 1])) + 1][d].size(r) >= int(ceil(T[(r) - 1]))) {
                        // std::cout << "r = " << r << "; n = " << n << "; T[(r) - 1] = " << T[(r) - 1] << "; d = " << d << "; ";
                        // std::cout << "dynmatrix[r][n][int(ceil(T[(r) - 1])) + 1][d] = " << dynmatrix[r][n][int(ceil(T[(r) - 1])) + 1][d] << "\n";
                        if(max_wct(r) < dynmatrix[r][n][int(ceil(T[(r) - 1])) + 1][d]) {
                            max_wct(r) = dynmatrix[r][n][int(ceil(T[(r) - 1])) + 1][d];
                            max_delayed_jobs[r] = JobDeviationList(deviated_jobs[r][n][int(ceil(T[(r) - 1])) + 1][d]);
                            max_delayed_jobs_all[r] = std::vector<JobDeviationList>();
                            max_delayed_jobs_all[r].push_back(JobDeviationList(deviated_jobs[r][n][int(ceil(T[(r) - 1])) + 1][d]));
                            max_C_prev_machine = boost::numeric::ublas::matrix<double>(C[r][n][int(ceil(T[(r) - 1])) + 1][d]);
                        }
                        /*else if(max_wct(r) <= dynmatrix[r][n][int(ceil(T[(r) - 1])) + 1][d]) {  // tie
                            JobDeviationList who = deviated_jobs[r][n][int(ceil(T[(r) - 1])) + 1][d];
                            if(std::find(max_delayed_jobs_all[r].begin(), max_delayed_jobs_all[r].end(), who) == max_delayed_jobs_all[r].end()) {  // if(k in deviated_jobs[i]) {
                            // if(!(deviated_jobs[r][n][int(ceil(Γ[r])) + 1][d] in max_delayed_jobs_all[r])) {
                                max_delayed_jobs_all[r].push_back(JobDeviationList(deviated_jobs[r][n][int(ceil(T[(r) - 1])) + 1][d]));
                                max_C_prev_machine = boost::numeric::ublas::matrix<double>(C[r][n][int(ceil(T[(r) - 1])) + 1][d]);
                            }
                        } */
                    }
                }
                // std::cout << "max_wct(r) = " << max_wct(r);
                // std::cout << "; max_delayed_jobs_all[r] = " << max_delayed_jobs_all[r][0].str() << "\n";
            } /*
            std::vector<JobDeviationList> Γ_scenario_list;
            for(JobDeviationList worst_Γ_scenario : max_delayed_jobs_all[m]) {
                bool fits = true;
                for(int i = 1; i <= m; i++) {
                    if(worst_Γ_scenario[i].size() != T[i]) {
                        fits = false;
                    }
                }
                if(fits)
                    Γ_scenario_list.push_back(worst_Γ_scenario);
            }
            if(verbose)
                std::cout << "[WCT_DP] The m-machine worst-case WCT is " << std::setprecision(10) << max_wct(m) << "\n"; */
            return max_wct(m);  // , Γ_scenario_list;
            // return worst_case_wct_dp(n, P_bar, P_hat, rob_params, pi, w);
        }
    };
}

#endif //PFSP_PFSP_WCT_BUDGET_WORSTCASE_H
