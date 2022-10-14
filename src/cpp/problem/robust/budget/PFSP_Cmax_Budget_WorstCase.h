//
// Created by mlevorato on 1/13/20.
//

#ifndef PFSP_PFSP_CMAX_BUDGET_WORSTCASE_H
#define PFSP_PFSP_CMAX_BUDGET_WORSTCASE_H

#include "../../deterministic/PFSProblem.h"
#include "../../PFSSolution.h"
#include "../../PFSP_Parameters.h"
#include "../../metaheuristic/common/Job.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;
    using namespace problem::pfsp;

    class PFSP_Cmax_Budget_WorstCase {
    public:
        /**
         * Procedure Worst Case - Continuous time intervals.
         * Calculates the worst-case makespan given a sequence / permutation pi and
         * budget parameter Gamma. Select between the machine budget type and global budget type.
         * @param nUsedJobs is the number of jobs in the partial permutation pi, where nUsedJobs <= n
         */
        static double worst_case_cmax_dp(const unsigned &nUsedJobs,
                                         const boost::numeric::ublas::matrix<double> &P_bar,
                                         const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                         const std::vector<Job> &pi) {
            if(rob_params.budget_type == "machine") {
                return worst_case_cmax_dp_machine_budget(nUsedJobs, P_bar, P_hat, rob_params, pi);
            } else if (rob_params.budget_type == "global") {
                return worst_case_cmax_dp_global_budget(nUsedJobs, P_bar, P_hat, rob_params, pi);
            } else {
                throw std::invalid_argument( "Invalid budget type!");
            }
        }

        /**
         * Procedure Worst Case - Continuous time intervals - Machine Budget.
         * Based on machine-budget type, calculates the worst-case makespan given a sequence / permutation pi and
         * budget parameter Gamma. Uses a dynamic programming based on a recursion.
         * Attention: original formulas (in Julia) were made for 1-based arrays (e.g. pi[k] with k = 1, ..., n).
         * Complexity: O(n ^ 3)
         * @param nUsedJobs is the number of jobs in the partial permutation pi, where nUsedJobs <= n
         */
        static double worst_case_cmax_dp_machine_budget(const unsigned &nUsedJobs,
                const boost::numeric::ublas::matrix<double> &P_bar,
                const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                const std::vector<Job> &pi) {
            unsigned k = 0, y1 = 0, y2 = 0, q = 0;
            unsigned m = P_bar.size1() - 1;
            unsigned n = P_bar.size2() - 1;  // n is the total number of jobs of the problem instance
            //std::cout << "*** The number of machines is " << m << "\nUsedJobs";
            // IMPORTANT: In case this procedure is invoked from a metaheuristic to find the worst case of a partial
            // solution pi where size(pi) < n, then we should apply T[i] as T[i] = min(T[i], nUsedJobs)
            std::vector<unsigned> T(m, 0);
            int max_gamma = n;
            for(y1 = 0; y1 < m; y1++) {
                T[y1] = (unsigned)ceil((n * rob_params.budget_gamma[y1]) / 100.0);
            }
            boost::multi_array<double, 3> C(boost::extents[m + 1][n + 2][max_gamma + 2]);
            std::fill_n(C.data(), C.num_elements(), 0);
            // LOG(INFO) << "[DP] Calculating worst-case Cmax given a permutation pi = $(pi), m = $(m)"
            //   << "and nUsedJobs = $(nUsedJobs) and budget parameters T = $((T1, T2))...";
            // Filling matrix C[1]
            for(y1 = 0; y1 <= max_gamma; y1++) {  // T1
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
                for (k = 1; k <= n; k++) {  // for k in 1:nUsedJobs
                    for (y2 = 0; y2 <= T[(q) - 1]; y2++) {  // for y2 in 0:T[q]
                        C[q][k][y2] = std::max(C[q][k][y2], std::max(C[q - 1][k][y1], C[q][k - 1][y2])
                                             + P_bar(q, pi[(k) - 1].id));
                        C[q][k][y2 + 1] = std::max(C[q][k][y2 + 1], std::max(C[q - 1][k][y1], C[q][k - 1][y2])
                                             + P_bar(q, pi[(k) - 1].id) + P_hat(q, pi[(k)- 1].id));
                    }
                }
            }
            double worst_cmax = C[m][nUsedJobs][T[(m) - 1]];
            // LOG(INFO) << "[DP] The worst-case Cmax is " << worst_cmax << ".\nUsedJobs";
            return worst_cmax;
        }

        /**
         * Procedure Worst Case - Continuous time intervals - Global Budget.
         * Based on global-budget type, calculates the worst-case makespan given a sequence / permutation pi and
         * budget parameter Gamma. Uses a dynamic programming based on a recursion.
         * Attention: original formulas (in Julia) were made for 1-based arrays (e.g. pi[k] with k = 1, ..., n).
         * Complexity: O(n ^ 3)
         * @param nUsedJobs is the number of jobs in the partial permutation pi, where nUsedJobs <= n
         */
        static double worst_case_cmax_dp_global_budget(const unsigned &nUsedJobs,
                                                        const boost::numeric::ublas::matrix<double> &P_bar,
                                                        const boost::numeric::ublas::matrix<double> &P_hat,
                                                        const Robust_Parameters &rob_params,
                                                        const std::vector<Job> &pi) {
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
            double worst_cmax = C[m][nUsedJobs][Gamma];
            // LOG(INFO) << "[DP] The worst-case Cmax is " << worst_cmax << ".\nUsedJobs";
            return worst_cmax;
        }

        static double worst_case_cmax_dp(const unsigned &nUsedJobs,
                                         const boost::numeric::ublas::matrix<double> &P_bar,
                                         const boost::numeric::ublas::matrix<double> &P_hat, const Robust_Parameters &rob_params,
                                         const std::vector<int> &pi) {
            std::vector<Job> pi2;
            for(int j : pi) {
                pi2.push_back(Job(j));
            }
            return worst_case_cmax_dp(nUsedJobs, P_bar, P_hat, rob_params, pi2);
        }

        /**
         * Procedure Worst Case - Continuous time intervals.
         *   Worst-case scenario calculation based on a heuristic.
         *   Use a simple greedy algorithm to obtain a lower bound on the worst-case solution value.
         *   Complexity: O(max(T1, T2) * n^2)
         * @return Worst-case makespan scenario (set of processing times) and corresponding value of makespan.
         */
        static double worst_case_heuristic(const unsigned &n,
                const boost::numeric::ublas::matrix<double> &p_bar, const boost::numeric::ublas::matrix<double> &p_hat,
                const Robust_Parameters &rob_params, const std::vector<Job> &s) {
            // LOG(INFO) << "[PFSP_Cmax_Budget_WorstCase] Calculating worst-case scenario (heuristic)...";
            unsigned m = 2;
            double T1 = (n * rob_params.budget_gamma[0]) / 100;
            double T2 = (n * rob_params.budget_gamma[1]) / 100;
            int floor_T1 = (int)floor(T1);  int ceil_T1 = (int)ceil(T1);
            int floor_T2 = (int)floor(T2);  int ceil_T2 = (int)ceil(T2);
            double max_makespan = 0.0;
            // {Set of jobs with deviated processing times in Machine M1} => |setT1| <= T1
            std::vector<int> deviated_M1;  // empty set
            std::vector<int> not_deviated_M1(n);
            std::iota (std::begin(not_deviated_M1), std::end(not_deviated_M1), 1); // Fill with 1, ..., n.
            // {Set of jobs with deviated processing times in Machine M2} => |setT2| <= T2
            std::vector<int> deviated_M2;  // empty set
            std::vector<int> not_deviated_M2(n);
            std::iota (std::begin(not_deviated_M2), std::end(not_deviated_M2), 1); // Fill with 1, ..., n.
            boost::numeric::ublas::matrix<double> p_time(p_bar);
            // Choose a pair of jobs (j1, j2) to have their processing time deviated,
            // in machines M1 and M2, respectively (i.e. j1 in M1 and j2 in M2).
            // (pick the pair of jobs that provides the maximum increase in the makespan)
            unsigned max_t = min(ceil_T1, ceil_T2);
            unsigned iter = 0;
            for(unsigned t = 1; t <= max_t; t++) {
                double max_deviation = 0.0, dev_1 = 0.0, dev_2 = 0.0;
                int which_j1 = 0, which_j2 = 0;
                double which_dev_1 = 0.0, which_dev_2 = 0.0;
                for (int j1 : not_deviated_M1) {  // For all pairs of deviated jobs (j1, j2)
                    for (int j2 : not_deviated_M2) {
                        // Calculate Cmax value when we consider the maximum deviation of processing time
                        if (iter == ceil_T1 && floor_T1 < ceil_T1) {
                            dev_1 = (T1 - floor_T1) * p_hat(1, j1);
                        } else {  // deviate processing time of job j1 to its maximum
                            dev_1 = p_hat(1, j1);
                        }
                        if (iter == ceil_T2 && floor_T2 < ceil_T2) {
                            dev_2 = (T2 - floor_T2) * p_hat(2, j2);
                        } else {  // deviate processing time of job j1 to its maximum
                            dev_2 = p_hat(2, j2);
                        }
                        p_time(1, j1) += dev_1;
                        p_time(2, j2) += dev_2;
                        double cmax_s = PFSProblem::calculateCmax(s, m, n, p_time);
                        if (cmax_s > max_deviation) {
                            max_deviation = cmax_s;
                            which_j1 = j1;
                            which_j2 = j2;
                            which_dev_1 = dev_1;
                            which_dev_2 = dev_2;
                        }
                        // undo changes to the matrix p_time
                        p_time(1, j1) -= dev_1;
                        p_time(2, j2) -= dev_2;
                    }
                }
                if (max_deviation >= max_makespan) {
                    // Include in deviated_T1 and deviated_T2 the jobs which caused the maximum increase in Cmax
                    deviated_M1.push_back(which_j1);
                    deviated_M2.push_back(which_j2);
                    p_time(1, which_j1) += which_dev_1;
                    p_time(2, which_j2) += which_dev_2;
                    not_deviated_M1.erase(std::remove(not_deviated_M1.begin(), not_deviated_M1.end(), which_j1), not_deviated_M1.end());
                    not_deviated_M2.erase(std::remove(not_deviated_M2.begin(), not_deviated_M2.end(), which_j2), not_deviated_M2.end());
                    iter++;
                } else {
                    // TODO if makespan is not augmented to the previous iteration, set the proc time of the
                    // first job pair (j1, j2) in the sequennce (not already deviated) to its maximum deviation
                }
            }
            for(unsigned t1 = iter + 1; t1 <= ceil_T1; t1++) {
                double max_deviation = 0.0, dev_1 = 0.0, which_dev_1 = 0.0;
                int which_j1 = 0;
                for (int j1 : not_deviated_M1) {  // For all possible deviated job j1
                    if (t1 == ceil_T1 && floor_T1 < ceil_T1) {
                        dev_1 = (T1 - floor_T1) * p_hat(1, j1);
                    } else {  // deviate processing time of job j1 to its maximum
                        dev_1 = p_hat(1, j1);
                    }
                    p_time(1, j1) += dev_1;
                    double cmax_s = PFSProblem::calculateCmax(s, m, n, p_time);
                    if (cmax_s > max_deviation) {
                        max_deviation = cmax_s;
                        which_j1 = j1;
                        which_dev_1 = dev_1;
                    }
                    // undo changes to the matrix p_time
                    p_time(1, j1) -= dev_1;
                }
                if (max_deviation >= max_makespan) {
                    // Include in deviated_T1 the job which caused the maximum increase in Cmax
                    deviated_M1.push_back(which_j1);
                    p_time(1, which_j1) += which_dev_1;
                    not_deviated_M1.erase(std::remove(not_deviated_M1.begin(), not_deviated_M1.end(), which_j1), not_deviated_M1.end());
                } else {
                    // if makespan is not augmented to the previous iteration, set the proc time of the
                    // first job in the sequennce (not already deviated) to its maximum deviation
                    int which_j1 = not_deviated_M1[0];
                    if (t1 == ceil_T1 && floor_T1 < ceil_T1) {
                        dev_1 = (T1 - floor_T1) * p_hat(1, which_j1);
                    } else {  // deviate processing time of job j1 to its maximum
                        dev_1 = p_hat(1, which_j1);
                    }
                    deviated_M1.push_back(which_j1);
                    p_time(1, which_j1) += dev_1;
                    not_deviated_M1.erase(std::remove(not_deviated_M1.begin(), not_deviated_M1.end(), which_j1), not_deviated_M1.end());
                }
            }
            for(unsigned t2 = iter + 1; t2 <= ceil_T2; t2++) {
                double max_deviation = 0.0, dev_2 = 0.0, which_dev_2 = 0.0;
                int which_j2 = 0;
                for (int j2 : not_deviated_M2) {  // For all possible deviated job j2
                    if (t2 == ceil_T2 && floor_T2 < ceil_T2) {
                        dev_2 = (T2 - floor_T2) * p_hat(2, j2);
                    } else {  // deviate processing time of job j1 to its maximum
                        dev_2 = p_hat(2, j2);
                    }
                    p_time(2, j2) += dev_2;
                    double cmax_s = PFSProblem::calculateCmax(s, m, n, p_time);
                    if (cmax_s > max_deviation) {
                        max_deviation = cmax_s;
                        which_j2 = j2;
                        which_dev_2 = dev_2;
                    }
                    // undo changes to the matrix p_time
                    p_time(2, j2) -= dev_2;
                }
                if (max_deviation >= max_makespan) {
                    // Include in deviated_T2 the job which caused the maximum increase in Cmax
                    deviated_M2.push_back(which_j2);
                    p_time(2, which_j2) += which_dev_2;
                    not_deviated_M2.erase(std::remove(not_deviated_M2.begin(), not_deviated_M2.end(), which_j2), not_deviated_M2.end());
                } else {
                    // if makespan is not augmented to the previous iteration, set the proc time of the
                    // first job in the sequennce (not already deviated) to its maximum deviation
                    int which_j2 = not_deviated_M2[0];
                    if (t2 == ceil_T2 && floor_T2 < ceil_T2) {
                        dev_2 = (T2 - floor_T2) * p_hat(2, which_j2);
                    } else {  // deviate processing time of job j1 to its maximum
                        dev_2 = p_hat(2, which_j2);
                    }
                    deviated_M2.push_back(which_j2);
                    p_time(2, which_j2) += dev_2;
                    not_deviated_M2.erase(std::remove(not_deviated_M2.begin(), not_deviated_M2.end(), which_j2), not_deviated_M2.end());
                }
            }
            assert(deviated_M1.size() == ceil_T1);
            assert(deviated_M2.size() == ceil_T2);
            // Assemble the scenario obtained
            /*
            deviated_M1.insert(deviated_M1.end(), not_deviated_M1.begin(), not_deviated_M1.end());
            deviated_M2.insert(deviated_M2.end(), not_deviated_M2.begin(), not_deviated_M2.end());
            std::vector<Permutation> pi_scen;
            pi_scen.push_back(deviated_M1);  // Machine M_1
            pi_scen.push_back(deviated_M2);  // Machine M_2
            PFSPBudgetScenario scenario(rob_params.budget_gamma, pi_scen);
            scenario.value = PFSProblem::calculateCmax(s, m, n, p_time); */
            return PFSProblem::calculateCmax(s, m, n, p_time);
        }

        static double worst_case_heuristic(const unsigned &n,
                                           const boost::numeric::ublas::matrix<double> &p_bar, const boost::numeric::ublas::matrix<double> &p_hat,
                                           const Robust_Parameters &rob_params, const std::vector<int> &s) {
            std::vector<Job> s2;
            for(int j : s) {
                s2.push_back(Job(j));
            }
            return worst_case_heuristic(n, p_bar, p_hat, rob_params, s2);
        }
    };
}

#endif //PFSP_PFSP_CMAX_BUDGET_WORSTCASE_H
