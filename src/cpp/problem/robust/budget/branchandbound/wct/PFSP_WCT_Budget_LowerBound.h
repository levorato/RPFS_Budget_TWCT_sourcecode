//
// Created by mlevorato on 3/1/21.
//

#ifndef PFSP_PFSP_WCT_BUDGET_LOWERBOUND_H
#define PFSP_PFSP_WCT_BUDGET_LOWERBOUND_H

#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o
#include <glog/logging.h>

#include "../../../../deterministic/PFSProblem.h"
#include "../../../RobPFSInstance_WCT.h"
#include "PFSP_WCT_Budget_NodeData.h"
#include "../../PFSP_WCT_Budget_WorstCase.h"
#include "../../../../../util/random.h"
#include "../../../../../util/include/TimeDateUtil.h"
#include "../../../../../util/NumericUtil.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using boost::timer::cpu_timer;

    class PFSP_WCT_Budget_LowerBound {
    public:
        typedef boost::random::beta_distribution<double> BETA_DISTRIBUTION;  // Beta distribution
        // stores statistics about the branch-and-bound process
        long num_leaf_nodes, num_leaf_nodes_infeasible, num_leaf_nodes_dp;
        double total_time_worstcase_mip, total_time_worstcase_dp;

        PFSP_WCT_Budget_LowerBound() : num_leaf_nodes(0), num_leaf_nodes_infeasible(0), 
            total_time_worstcase_mip(0.0), total_time_worstcase_dp(0.0), num_leaf_nodes_dp(0) {

        }

        /**
         * Calculate the objective function value of a given leaf node, induced by permutation partial_seq.
         */
        double calculate_obj_leaf_node(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                         const Permutation &partial_seq, Time current_lb) {
            double wct = 0.0;
            num_leaf_nodes++;
            boost::timer::cpu_timer timer;
            timer.start();  // Time measurement
            // Use deterministic objective validation if gamma == 0 or gamma == 100 (saves time)
            if((rob_params.budget_gamma[0] == 0) || (rob_params.budget_gamma[0] == 100)) {
                wct = PFSProblem::calculate_wct(partial_seq, I.m, I.n, I.get_oscillated_p_matrix(), I.w);
            } else {
                std::vector<Job> pi;
                for (int i = 0; i < I.n; i++)
                    pi.push_back(Job(partial_seq[i] - 1, I.m));
                // First, let's try to calculate the objective with the lazy obj function (relaxation)
                // TODO FIXME Check if using the relaxation DP is valid
                wct = I.calculate_worst_case_wct(pi, I.n, 4, false);
                total_time_worstcase_dp += TimeDateUtil::calculate_time_spent(timer);
                if(NumericUtil::TestsetIsGE(wct, ub)) {  // wct >= ub
                    num_leaf_nodes_dp++;
                    return wct;
                }
                // If the relaxation fails, use the full obj function calculation
                wct = I.calculate_worst_case_wct(pi, I.n, 0, false);  //, std::make_pair(current_lb, ub));
                total_time_worstcase_mip += PFSP_WCT_Budget_WorstCase::time_spent_worstcase_mip;
            }
            if(wct == std::numeric_limits<double>::max()) {
                num_leaf_nodes_infeasible++;
            }
            // std::cout << "Leaf node wlb = " << wlb << " vs. wct = " << wct << "\n";
            return wct;
        }

        double lb(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                         const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list, 
                         PFSP_WCT_Budget_NodeData &node_data) {
            double lb = 0.0;
            if((rob_params.budget_gamma[0] == 0) || (rob_params.budget_gamma[0] == 100)) {
                lb = lb_1(I, ub, rob_params, partial_seq, unscheduled_job_list);
                if(partial_seq.size() == I.n) {  // leaf node
                    lb = calculate_obj_leaf_node(I, ub, rob_params, partial_seq, lb);
                }
                return lb;
            }
            Robust_Parameters rob_params2(rob_params);
            // rob_params2.budget_gamma[0] = 0;
            lb = lb_2(I, ub, rob_params2, partial_seq, unscheduled_job_list);
            int src[] = { 6, 2, 3, 5, 10, 7, 8 };
            int n = sizeof(src) / sizeof(src[0]);
            std::vector<int> v(src, src + n);
            if(v == partial_seq) {
                std::cerr << "****** Node: lb = " << lb << "\n";
            }
            /* if(lb < ub) {
                double newlb = lb_2(I, ub, rob_params, partial_seq, unscheduled_job_list);
                // node_data.mem_lb_sum, node_data.mem_lb_count);
                if(newlb > lb)
                    lb = newlb;
            } */
            if(partial_seq.size() == I.n) {  // leaf node
                lb = calculate_obj_leaf_node(I, ub, rob_params, partial_seq, lb);
            }
            return lb;
        }

        /**
         * Sort the unscheduled jobs from unscheduled_job_list according to their processing time on
         * machine k, in ascending order.
         * @param unscheduled_job_list
         * @param p_time
         * @return
         */
        std::vector<int> sort_unscheduled_jobs_by_ptime(const std::vector<int> &unscheduled_job_list,
                                                      const ublas::matrix<double>& p_time,
                                                      const unsigned& k) {
            // Let p([s+1], k), p([s+2], k), ..., p([n, k) denote a permutation of the processing times of the jobs in
            // unscheduled_job_list U on machine k, which satisfies p([s+1], k) <= p([s+2], k) <= ... <= p([n], k):
            std::vector< std::tuple<int, int> > v;
            for(int j : unscheduled_job_list) {
                v.emplace_back(p_time(k, j), j);
            }
            // Using sort() function to sort by 1st element of tuple, in non-descending order
            //std::sort(v.begin(), v.end());
            std::sort(v.begin(), v.end(), [](const std::tuple<int, int> &i, const std::tuple<int, int> &j) {
                return i < j;
            });
            std::vector<int> permutation;
            stringstream ss;
            for (auto i : v){
                permutation.push_back(std::get<1>(i));
                ss << std::get<0>(i) << " ";
            }
            //std::cout << "sort_unscheduled_jobs_by_ptime: " << ss.str() << "\n";
            return permutation;
        }

        /**
         * Sort the unscheduled jobs from unscheduled_job_list according to their weight w, in descending order.
         * @param unscheduled_job_list
         * @param w
         * @return
         */
        std::vector<int> sort_unscheduled_jobs_by_job_weight(const std::vector<int> &unscheduled_job_list,
                                                                    const ublas::vector<double>& w) {
            // Let w[s+1], w[s+2], ..., w[n] denote a permutation of the weights of the jobs in U
            // which satisfies w[s+1] >= w[s+2] >= ... >= w[n] :
            std::vector< std::tuple<int, int> > v;
            for(int j : unscheduled_job_list) {
                v.emplace_back(w(j), j);
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(v.begin(), v.end(), [](const std::tuple<int, int> &i, const std::tuple<int, int> &j) {
                return i > j;
            });
            std::vector<int> permutation;
            stringstream ss;
            for (auto i : v){
                permutation.push_back(std::get<1>(i));
                ss << std::get<0>(i) << " ";
            }
            //std::cout << "sort_unscheduled_jobs_by_job_weight: " << ss.str() << "\n";
            return permutation;
        }

        /**
         * Calculate which operations (r, j) will have their processing times oscillated on LB_2.
         * First includes as many operations related to unscheduled jobs from unscheduled_job_list on machine k.
         * Then includes as many operations as possible, related to scheduled jobs given all machines r. 
         * Total number of oscillated operations follows available gamma budget parameter. 
         * Returns a list of the oscillated operations (r, j) regarding jobs j on machines r.
         * @param unscheduled_job_list
         * @param p_time
         * @return
         */
        std::vector< std::tuple<int,int> > get_oscillated_operations_by_ptime_budget_lb2(const Permutation &partial_seq, 
                                                      const std::vector<int> &unscheduled_job_list,
                                                      const ublas::matrix<double>& p_bar,
                                                      const ublas::matrix<double>& p_hat,
                                                      const unsigned& gamma_value,
                                                      const unsigned& k,
                                                      const unsigned& m, 
                                                      const bool& ascending = false) {
            // 1. Try to oscillate the maximum number of operations on machine k, according to available budget
            //    Simple greedy algorithm based on largest oscillation values p_hat(k, j)
            std::vector< std::tuple<int, int> > op_k;
            op_k.reserve(unscheduled_job_list.size());
            for(int j : unscheduled_job_list) {
                op_k.emplace_back(p_bar(k, j) + p_hat(k, j), j);
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(op_k.begin(), op_k.end(), [ascending](const std::tuple<int, int> &i, const std::tuple<int, int> &j) {
                if(ascending)  return i < j;
                else  return i > j;
            });
            std::vector< std::tuple<int, int> > oscillated_operations;
            stringstream ss;
            int used_budget = 0;
            for (auto i : op_k){  // select at most gamma operations (k, j) from list op_k
                if(NumericUtil::TestsetIsGE(used_budget, gamma_value))  break;
                oscillated_operations.emplace_back(k, std::get<1>(i));
                ss << "operation (" << k << ", " << std::get<1>(i) << ") -> p_time = " << std::get<0>(i) << ", ";
                used_budget++;
            }
            ss << "TOTAL = " << used_budget << "\n";
            //std::cout << "oscillared_operations: " << ss.str();
            // 2. Try to oscillate operations regarding scheduled jobs (partial permutation), if remaining budget allows
            std::vector< std::tuple<int, int> > op_s;
            op_s.reserve(partial_seq.size() * m);
            for(int j : partial_seq) {
                for(int r = 1; r <= m; r++) {
                    op_s.emplace_back(r, j);
                }
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(op_s.begin(), op_s.end(), [p_bar, p_hat, ascending](const std::tuple<int, int> &x, const std::tuple<int, int> &y) {
                if(ascending)  return NumericUtil::TestsetIsLT(p_bar(std::get<0>(x), std::get<1>(x)) + p_hat(std::get<0>(x), std::get<1>(x)), 
                    p_bar(std::get<0>(y), std::get<1>(y)) + p_hat(std::get<0>(y), std::get<1>(y)));  // a < b
                else  return NumericUtil::TestsetIsGT(p_bar(std::get<0>(x), std::get<1>(x)) + p_hat(std::get<0>(x), std::get<1>(x)), 
                    p_bar(std::get<0>(y), std::get<1>(y)) + p_hat(std::get<0>(y), std::get<1>(y)));  // a > b
            });
            for (auto i : op_s){  // select at most (gamma - used_budget) operations (r, j) from list op_s
                if(NumericUtil::TestsetIsGE(used_budget, gamma_value))  break;
                oscillated_operations.emplace_back(std::get<0>(i), std::get<1>(i));
                ss << "operation (" << std::get<0>(i) << ", " << std::get<1>(i) << ") -> p_time = " 
                    << (p_bar(std::get<0>(i), std::get<1>(i)) + p_hat(std::get<0>(i), std::get<1>(i))) << ", ";
                used_budget++;
            }
            ss << "TOTAL = " << used_budget << "\n";
            //std::cout << "oscillared_operations: " << ss.str();
            return oscillated_operations;
        }

        /**
        * Procedure LB_1 - Lower Bound.
        * Based on the paper:
        * Chung, C. S., Flynn, J., & Kirca, O. (2002). A branch and bound algorithm to minimize the total flow time
        * for m-machine permutation flowshop problems. International Journal of Production Economics, 79(3), 185-196.
        * Complexity: O(max(T1, T2) * n^2)
        * @return Lower Bound on the worst-case makespan.
        */
        double lb_1(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                           const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list) {
            // LOG(INFO) << "[PFSPBudgetScenario] Calculating lower bound...";
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            // If Gamma == 0, the following 2 lines will use the nominal P_bar matrix as reference
            // If Gamma == 100%, the 2 calculations above will be based on (P_bar + P_hat) matrix as P
            const ublas::matrix<double>& p_time = (rob_params.budget_gamma[0] > 0) ? I.get_oscillated_p_matrix() : I.get_p_bar();
            const boost::multi_array<double, 3>& q = I.get_l_kl_matrix();
            // Completion time matrix of partial schedule \sigma
            ublas::matrix<double> C = PFSProblem::makespan(partial_seq, I.m, p_time);
            // std::cout << "C_matrix: " << C << "\n";
            // std::cout << "partial_seq: " << partial_seq << "\n";
            // Total weighted flow time for jobs in the partial schedule \sigma
            // TODO Maybe include worst-case calculation in this part of the formula
            double g_w_sigma = 0.0;
            for(int t = 1; t <= s; t++) {
                g_w_sigma += I.w(partial_seq[t - 1]) * C(I.m, t);
            }
            // double g_w_sigma = PFSProblem::calculate_wct(partial_seq, I.m, partial_seq.size(), p_time, I.w);

            std::vector<int> unscheduled_permutation_w_sorted = sort_unscheduled_jobs_by_job_weight(unscheduled_job_list, I.w);

            // Etk is an underestimate of the earliest start time of the t-th job on machine k
            // R(t, k) is an underestimate of the earliest completion time of the t-th job on machine k
            // Rules for machine 1
            I.E(s+1, 1) = C(1, s);  // L_1(\sigma) : Completion time of the last job of \sigma on machine 1
            I.R(s, 1) = C(1, s);
            std::vector<int> unscheduled_permutation_m1 = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, 1);
            for(int t = s + 1; t <= I.n; t++) {
                double sum_p_r1 = 0.0;
                for (int r = s + 1; r <= t; r++) {
                    sum_p_r1 += p_time(1, unscheduled_permutation_m1[r - (s + 1)]);
                }
                I.R(t, 1) = I.R(t-1, 1) + p_time(1, unscheduled_permutation_m1[t - (s + 1)]);  // I.R(s, 1) + sum_p_r1;
                double aux = I.R(s, 1) + sum_p_r1;
                if(fabs(aux - I.R(t, 1)) > 1e-5) {
                    std::cerr << "Error on I.R(t, 1) calculation! " << aux << " x " << I.R(t, 1) << "\n";
                }
            }
            for(int t = s + 2; t <= I.n; t++) {
                I.E(t, 1) = I.R(t-1, 1);
            }
            for(int k = 2; k <= I.m; k++) {  // for each machine k
                I.E(s+1, k) = C(k, s);  // L_k(\sigma) : Completion time of the last job of \sigma on machine k;
                for(int r = 1; r <= k - 1; r++) {
                    double smallest_q = std::numeric_limits<double>::max();
                    if(unscheduled_job_list.size() == 0) {
                        //std::cerr << "LB WARN: unscheduled_job_list.size() == 0\n";
                        smallest_q = std::numeric_limits<double>::min();
                    }
                    for(int j : unscheduled_job_list) {
                        smallest_q = std::min(smallest_q, q[j][r][k - 1]);
                    }
                    I.E(s+1, k) = std::max(I.E(s+1, k), I.E(s+1, r) + smallest_q, less_eps);  // C(r, s) + smallest_q);  //
                }
                std::vector<int> unscheduled_permutation_mk = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, k);
                I.R(s, k) = I.E(s + 1, k);
                for(int t = s + 1; t <= I.n; t++) {
                    double sum_p_rk = 0.0;
                    for(int r = s + 1; r <= t; r++) {
                        sum_p_rk += p_time(k, unscheduled_permutation_mk[r - (s + 1)]);
                    }
                    I.R(t, k) = I.R(t-1, k) + p_time(k, unscheduled_permutation_mk[t - (s + 1)]);  // I.R(s, k) + sum_p_rk;
                    double aux = I.R(s, k) + sum_p_rk;
                    if(fabs(aux - I.R(t, k)) > 1e-5) {
                        std::cerr << "Error on I.R(t, k) calculation! " << aux << " x " << I.R(t, k) << "\n";
                    }
                }
                for(int t = s + 2; t <= I.n; t++) {
                    I.E(t, k) = std::max(I.R(t-1, k), I.R(t, k-1), less_eps);
                }
            }
            // WLB below is a lower bound on the weighted flow time at node \sigma:
            double max_wlb_k = 0.0;
            for(int k = 1; k <= I.m; k++) {
                double wlb_k = 0.0;
                for(int t = s + 1; t <= I.n; t++) {
                    wlb_k += I.w(unscheduled_permutation_w_sorted[t - (s + 1)]) * I.E(t, k);
                }
                for(int j : unscheduled_job_list) {
                    wlb_k += I.w(j) * q[j][k][I.m];
                }
                max_wlb_k = std::max(max_wlb_k, wlb_k, less_eps);
            }
            double wlb = max_wlb_k + g_w_sigma;
            // std::cout << "Leaf node wlb = " << wlb << "\n";
            return wlb;
            // PFSPBudgetScenario scenario = PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(I.n, I.p_bar, I.p_hat,
            //                                                                                rob_params, pi, I.w, true,
            //                                                                                false, 0, true);
            // return PFSP_WCT_Budget_WorstCase::worst_case_wct_dp(I.n, I.p_bar, I.p_hat, rob_params, pi, I.w);
        }

        /**
        * Procedure LB_2 - Lower Bound. LB_1 with special modification on oscillated jobs in critical path relaxation.
        * Based on the paper:
        * Chung, C. S., Flynn, J., & Kirca, O. (2002). A branch and bound algorithm to minimize the total flow time
        * for m-machine permutation flowshop problems. International Journal of Production Economics, 79(3), 185-196.
        * Complexity: O(max(T1, T2) * n^2)
        * @return Lower Bound on the worst-case makespan.
        */
        double lb_2(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                           const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list) {
            // LOG(INFO) << "[PFSPBudgetScenario] Calculating lower bound...";
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            ublas::matrix<double> p_time = I.get_p_bar();  // I.get_oscillated_p_matrix();
            const ublas::matrix<double>& p_bar = I.get_p_bar();
            const ublas::matrix<double>& p_hat = I.get_p_hat();
            double gamma = ((I.m * I.n * rob_params.budget_gamma[0]) / 100.0);
            int max_gamma = (int)floor(gamma);
            // Retrieve the q matrix, whose calculation is based on non-oscillated processing times p_bar
            boost::multi_array<double, 3> q = I.get_l_kl_matrix();
            // std::cout << "partial_seq: " << partial_seq << "\n";
            // Total weighted flow time for jobs in the partial schedule \sigma
            // TODO Maybe include worst-case calculation in this part of the formula
            ublas::vector<double> g_w_sigma = ublas::zero_vector<double>(I.m + 1);
            ublas::vector<double> w_q_unscheduled = ublas::zero_vector<double>(I.m + 1);
            std::vector<int> unscheduled_permutation_w_sorted = sort_unscheduled_jobs_by_job_weight(unscheduled_job_list, I.w);

            // Etk is an underestimate of the earliest start time of the t-th job on machine k
            // R(t, k) is an underestimate of the earliest completion time of the t-th job on machine k
            // Rules for machine 1
            std::vector< std::tuple<int,int> > oscillated_operations = get_oscillated_operations_by_ptime_budget_lb2(
                    partial_seq, unscheduled_job_list, p_bar, p_hat, max_gamma, 1, I.m, false);
            // FIXME TEST 1
            ////// ublas::matrix<double> C = PFSProblem::makespan(partial_seq, I.m, p_time);
            // k = 1 : update p_time and q matrices with oscillated operations
            for(std::tuple<int,int> op : oscillated_operations) {
                p_time(std::get<0>(op), std::get<1>(op)) += p_hat(std::get<0>(op), std::get<1>(op));
                int k1 = 1;
                for(int k2 = k1; k2 <= I.m; k2++) {
                    if(NumericUtil::TestsetIsGE(std::get<0>(op), k1) && NumericUtil::TestsetIsLE(std::get<0>(op), k2)) {  // if r >= k1 && r <= k2
                        q[std::get<1>(op)][k1][k2] += p_hat(std::get<0>(op), std::get<1>(op));
                    }
                }
            }
            // FIXME TEST 1
            ublas::matrix<double> C = PFSProblem::makespan(partial_seq, I.m, p_time);
            for(int t = 1; t <= s; t++) {
                g_w_sigma[1] += I.w(partial_seq[t - 1]) * C(I.m, t);
            }
            I.E(s+1, 1) = C(1, s);  // L_1(\sigma) : Completion time of the last job of \sigma on machine 1
            I.R(s, 1) = C(1, s);
            // FIXME TEST 2
            std::vector<int> unscheduled_permutation_m1 = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, 1);
            for(int t = s + 1; t <= I.n; t++) {
                double sum_p_r1 = 0.0;
                for (int r = s + 1; r <= t; r++) {
                    sum_p_r1 += p_time(1, unscheduled_permutation_m1[r - (s + 1)]);
                }
                //// I.R(t, 1) = I.R(t-1, 1) + p_time(1, unscheduled_permutation_m1[t - (s + 1)]);  // I.R(s, 1) + sum_p_r1;
                I.R(t, 1) = I.R(t-1, 1) + p_time(1, unscheduled_permutation_m1[t - (s + 1)]);  // I.R(s, 1) + sum_p_r1;
                double aux = I.R(s, 1) + sum_p_r1;
                if(fabs(aux - I.R(t, 1)) > 1e-5) {
                    std::cerr << "Error on I.R(t, 1) calculation! " << aux << " x " << I.R(t, 1) << "\n";
                }
            }
            for(int t = s + 2; t <= I.n; t++) {
                I.E(t, 1) = I.R(t-1, 1);
            }
            for(int j : unscheduled_job_list) {
                w_q_unscheduled(1) += I.w(j) * q[j][1][I.m];
            }
            // k = 1 : re-update p_time and q matrices with oscillated operations
            for(std::tuple<int,int> op : oscillated_operations) {
                p_time(std::get<0>(op), std::get<1>(op)) -= p_hat(std::get<0>(op), std::get<1>(op));
                int k1 = 1;
                for(int k2 = k1; k2 <= I.m; k2++) {
                    if(NumericUtil::TestsetIsGE(std::get<0>(op), k1) && NumericUtil::TestsetIsLE(std::get<0>(op), k2)) {  // if r >= k1 && r <= k2
                        q[std::get<1>(op)][k1][k2] -= p_hat(std::get<0>(op), std::get<1>(op));
                    }
                }
            }
            for(int k = 2; k <= I.m; k++) {  // for each machine k
                oscillated_operations = get_oscillated_operations_by_ptime_budget_lb2(
                    partial_seq, unscheduled_job_list, p_bar, p_hat, max_gamma, k, I.m, false);
                // FIXME TEST 1
                //// C = PFSProblem::makespan(partial_seq, I.m, p_time);
                // update p_time matrix with oscillated operations O_(r, j)
                for(std::tuple<int,int> op : oscillated_operations) {
                    p_time(std::get<0>(op), std::get<1>(op)) += p_hat(std::get<0>(op), std::get<1>(op));
                    for(int k1 = 1; k1 <= k; k1++) {
                        for(int k2 = k1; k2 <= I.m; k2++) {
                            if(NumericUtil::TestsetIsGE(std::get<0>(op), k1) && NumericUtil::TestsetIsLE(std::get<0>(op), k2)) {  // if r >= k1 && r <= k2
                                q[std::get<1>(op)][k1][k2] += p_hat(std::get<0>(op), std::get<1>(op));
                            }
                        }
                    }
                }
                // TODO Deviate the processing times of the remaining (gamma - |oscillated_mk|) jobs, on the other machines != k
                //     ( i.e., for current C matrix calculation )
                // FIXME TEST 1: Completion time matrix of partial schedule \sigma
                C = PFSProblem::makespan(partial_seq, I.m, p_time);
                // std::cout << "C_matrix: " << C << "\n";
                for(int t = 1; t <= s; t++) {
                    g_w_sigma[k] += I.w(partial_seq[t - 1]) * C(I.m, t);
                }
                // double g_w_sigma = PFSProblem::calculate_wct(partial_seq, I.m, partial_seq.size(), p_time, I.w);

                I.E(s+1, k) = C(k, s);  // L_k(\sigma) : Completion time of the last job of \sigma on machine k;
                for(int r = 1; r <= k - 1; r++) {
                    double smallest_q = std::numeric_limits<double>::max();
                    if(unscheduled_job_list.size() == 0) {
                        //std::cerr << "LB WARN: unscheduled_job_list.size() == 0\n";
                        smallest_q = 0.0;
                    }
                    for(int j : unscheduled_job_list) {
                        smallest_q = std::min(smallest_q, q[j][r][k - 1]);
                    }
                    I.E(s+1, k) = std::max(I.E(s+1, k), I.E(s+1, r) + smallest_q, less_eps);  // C(r, s) + smallest_q);  //
                    ////////I.E(s+1, k) = std::max(I.E(s+1, k), C(r, s) + smallest_q);
                }
                // iterate over remaining unscheduled jobs t
                // FIXME TEST 2
                std::vector<int> unscheduled_permutation_mk = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, k);
                I.R(s, k) = I.E(s + 1, k);
                for(int t = s + 1; t <= I.n; t++) {  
                    double sum_p_rk = 0.0;
                    for(int r = s + 1; r <= t; r++) {
                        sum_p_rk += p_time(k, unscheduled_permutation_mk[r - (s + 1)]);
                    }
                    //// I.R(t, k) = I.R(t-1, k) + p_time(k, unscheduled_permutation_mk[t - (s + 1)]);  // I.R(s, k) + sum_p_rk;
                    I.R(t, k) = I.R(t-1, k) + p_time(k, unscheduled_permutation_mk[t - (s + 1)]);  // I.R(s, k) + sum_p_rk;
                    double aux = I.R(s, k) + sum_p_rk;
                    if(fabs(aux - I.R(t, k)) > 1e-5) {
                        std::cerr << "Error on I.R(t, k) calculation! " << aux << " x " << I.R(t, k) << "\n";
                    }
                }
                for(int t = s + 2; t <= I.n; t++) {
                    I.E(t, k) = std::max(I.R(t-1, k), I.R(t, k-1), less_eps);
                }
                for(int j : unscheduled_job_list) {
                    w_q_unscheduled(k) += I.w(j) * q[j][k][I.m];
                }
                // update p_time matrix with oscillated operations O_(r, j)
                for(std::tuple<int,int> op : oscillated_operations) {
                    p_time(std::get<0>(op), std::get<1>(op)) -= p_hat(std::get<0>(op), std::get<1>(op));
                    for(int k1 = 1; k1 <= k; k1++) {
                        for(int k2 = k1; k2 <= I.m; k2++) {
                            if(NumericUtil::TestsetIsGE(std::get<0>(op), k1) && NumericUtil::TestsetIsLE(std::get<0>(op), k2)) {  // if r >= k1 && r <= k2
                                q[std::get<1>(op)][k1][k2] -= p_hat(std::get<0>(op), std::get<1>(op));
                            }
                        }
                    }
                }
            }
            // WLB below is a lower bound on the weighted flow time at node \sigma:
            double max_wlb_k = 0.0;
            for(int k = 1; k <= I.m; k++) {
                double wlb_k = g_w_sigma[k];
                for(int t = s + 1; t <= I.n; t++) {
                    wlb_k += I.w(unscheduled_permutation_w_sorted[t - (s + 1)]) * I.E(t, k);
                }
                wlb_k += w_q_unscheduled(k);  // replaces: wlb_k += I.w(j) * q[j][k][I.m];
                max_wlb_k = std::max(max_wlb_k, wlb_k, less_eps);
            }
            double wlb = max_wlb_k;
            // std::cout << "Leaf node wlb = " << wlb << "\n";
            return wlb;
            // PFSPBudgetScenario scenario = PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(I.n, I.p_bar, I.p_hat,
            //                                                                                rob_params, pi, I.w, true,
            //                                                                                false, 0, true);
            // return PFSP_WCT_Budget_WorstCase::worst_case_wct_dp(I.n, I.p_bar, I.p_hat, rob_params, pi, I.w);
        }

        /**
         * Multi-scenario lower bound for the PFSP-TWCT Master Problem (C&CG algorithm).
         * Returns the highest lower bound, given all scenarios / processing time matrices in P_time_scenario_list.
         * Besides, P_time_scenario_list, this method also uses pre-calculated l_kl matrices, stored in l_kl_scenario_list.
         */
        double lb_multi_scenario(RobPFSInstance_WCT &I, Time &ub, 
                           const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list, 
                           std::vector< ublas::matrix<double> > &P_time_scenario_list,
                           std::vector< boost::multi_array<double, 3> > &l_kl_scenario_list) {
            // LOG(INFO) << "[PFSPBudgetScenario] Calculating lower bound...";
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            std::vector<double> wlb_list;
            for(int scenario = 0; scenario < P_time_scenario_list.size(); scenario++) {
                const ublas::matrix<double>& p_time = P_time_scenario_list[scenario];
                const boost::multi_array<double, 3>& q = l_kl_scenario_list[scenario];
                
                // Completion time matrix of partial schedule \sigma
                ublas::matrix<double> C = PFSProblem::makespan(partial_seq, I.m, p_time);
                // std::cout << "C_matrix: " << C << "\n";
                // std::cout << "partial_seq: " << partial_seq << "\n";
                // Total weighted flow time for jobs in the partial schedule \sigma
                // TODO Maybe include worst-case calculation in this part of the formula
                double g_w_sigma = 0.0;
                for(int t = 1; t <= s; t++) {
                    g_w_sigma += I.w(partial_seq[t - 1]) * C(I.m, t);
                }
                // double g_w_sigma = PFSProblem::calculate_wct(partial_seq, I.m, partial_seq.size(), p_time, I.w);

                std::vector<int> unscheduled_permutation_w_sorted = sort_unscheduled_jobs_by_job_weight(unscheduled_job_list, I.w);

                // Etk is an underestimate of the earliest start time of the t-th job on machine k
                // R(t, k) is an underestimate of the earliest completion time of the t-th job on machine k
                // Rules for machine 1
                I.E(s+1, 1) = C(1, s);  // L_1(\sigma) : Completion time of the last job of \sigma on machine 1
                I.R(s, 1) = C(1, s);
                std::vector<int> unscheduled_permutation_m1 = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, 1);
                for(int t = s + 1; t <= I.n; t++) {
                    double sum_p_r1 = 0.0;
                    for (int r = s + 1; r <= t; r++) {
                        sum_p_r1 += p_time(1, unscheduled_permutation_m1[r - (s + 1)]);
                    }
                    I.R(t, 1) = I.R(t-1, 1) + p_time(1, unscheduled_permutation_m1[t - (s + 1)]);  // I.R(s, 1) + sum_p_r1;
                    double aux = I.R(s, 1) + sum_p_r1;
                    if(fabs(aux - I.R(t, 1)) > 1e-5) {
                        LOG(ERROR) << "Error on I.R(t, 1) calculation! " << aux << " x " << I.R(t, 1) << "\n";
                    }
                }
                for(int t = s + 2; t <= I.n; t++) {
                    I.E(t, 1) = I.R(t-1, 1);
                }
                for(int k = 2; k <= I.m; k++) {  // for each machine k
                    I.E(s+1, k) = C(k, s);  // L_k(\sigma) : Completion time of the last job of \sigma on machine k;
                    for(int r = 1; r <= k - 1; r++) {
                        double smallest_q = std::numeric_limits<double>::max();
                        if(unscheduled_job_list.size() == 0) {
                            //std::cerr << "LB WARN: unscheduled_job_list.size() == 0\n";
                            smallest_q = std::numeric_limits<double>::min();
                        }
                        for(int j : unscheduled_job_list) {
                            smallest_q = std::min(smallest_q, q[j][r][k - 1]);
                        }
                        I.E(s+1, k) = std::max(I.E(s+1, k), I.E(s+1, r) + smallest_q, less_eps);  // C(r, s) + smallest_q);  //
                    }
                    std::vector<int> unscheduled_permutation_mk = sort_unscheduled_jobs_by_ptime(unscheduled_job_list, p_time, k);
                    I.R(s, k) = I.E(s + 1, k);
                    for(int t = s + 1; t <= I.n; t++) {
                        double sum_p_rk = 0.0;
                        for(int r = s + 1; r <= t; r++) {
                            sum_p_rk += p_time(k, unscheduled_permutation_mk[r - (s + 1)]);
                        }
                        I.R(t, k) = I.R(t-1, k) + p_time(k, unscheduled_permutation_mk[t - (s + 1)]);  // I.R(s, k) + sum_p_rk;
                        double aux = I.R(s, k) + sum_p_rk;
                        if(fabs(aux - I.R(t, k)) > 1e-5) {
                            LOG(ERROR) << "Error on I.R(t, k) calculation! " << aux << " x " << I.R(t, k) << "\n";
                        }
                    }
                    for(int t = s + 2; t <= I.n; t++) {
                        I.E(t, k) = std::max(I.R(t-1, k), I.R(t, k-1), less_eps);
                    }
                }
                // WLB below is a lower bound on the weighted flow time at node \sigma in the current scenario
                double max_wlb_k = 0.0;
                for(int k = 1; k <= I.m; k++) {
                    double wlb_k = 0.0;
                    for(int t = s + 1; t <= I.n; t++) {
                        wlb_k += I.w(unscheduled_permutation_w_sorted[t - (s + 1)]) * I.E(t, k);
                    }
                    for(int j : unscheduled_job_list) {
                        wlb_k += I.w(j) * q[j][k][I.m];
                    }
                    max_wlb_k = std::max(max_wlb_k, wlb_k, less_eps);
                }
                double wlb_scenario = max_wlb_k + g_w_sigma;
                // std::cout << "Leaf node wlb = " << wlb << "\n";
                wlb_list.push_back(wlb_scenario);
            }
            // Return the highest lower bound, considering all provided scenarios
            double wlb = *max_element(std::begin(wlb_list), std::end(wlb_list), less_eps); // c++11;
            return wlb;
        }

        bool dominance_rules(RobPFSInstance_WCT &instance, Time &ub, const Robust_Parameters &rob_params,
                                    Permutation &partial_seq, const std::vector<int> &unscheduled_job_list) {
            
            if(rob_params.dominance_rules) {
                if((rob_params.budget_gamma[0] == 0) || (rob_params.budget_gamma[0] == 100)) {
                    return dominance_1(instance, ub, rob_params, partial_seq, unscheduled_job_list);
                }
                bool phatom = dominance_2(instance, ub, rob_params, partial_seq, unscheduled_job_list);
                /**
                int src[] = { 6, 2, 3, 5, 10, 7, 8 };
                int n = sizeof(src) / sizeof(src[0]);
                std::vector<int> v(src, src + n);
                if(v == partial_seq) {
                    std::cerr << "****** Node: dominance = " << phatom << "\n";
                    std::cerr << result.first << " >= " << result.second << "\n";
                    phatom = false;
                    std::cerr << "****** Node: dominance = " << phatom << "\n";
                } */
                if(phatom) {
                    return true;
                }
            }
            return false;
        }

        /**
        * Procedure Dominance_1 - Dominance Rules - Returns true if sequence s violates dominance rules.
        * Based on the paper:
        * Chung, C. S., Flynn, J., & Kirca, O. (2002). A branch and bound algorithm to minimize the total flow time
        * for m-machine permutation flowshop problems. International Journal of Production Economics, 79(3), 185-196.
        * WARNING: This dominance rule 1 DOES NOT WORK for the robust case, except when Gamma == 100%.
        * @return true if sequence s violates dominance rules.
        */
        bool dominance_1(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                        Permutation partial_seq, const std::vector<int> &unscheduled_job_list) {
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            // If Gamma == 0, the following calculation will use the nominal P_bar matrix as reference
            // If Gamma == 100%, the calculation below will be based on (P_bar + P_hat) matrix as P
            const ublas::matrix<double>& p_time = (rob_params.budget_gamma[0] > 0) ? I.get_oscillated_p_matrix() : I.get_p_bar();
            double g_w_sigma_ij = 0.0;
            double g_w_sigma_ji = 0.0;
            double sum_w_U = 0.0;
            // Completion time matrix of partial schedule \sigma
            boost::numeric::ublas::matrix<double> C_sigma_ij = boost::numeric::ublas::zero_matrix<double>(I.m + 1, I.n + 2);
            PFSProblem::makespan(C_sigma_ij, partial_seq, I.m, p_time);
            boost::numeric::ublas::matrix<double> C_sigma_ji(C_sigma_ij);
            // std::cout << "C_matrix: " << C << "\n";
            // std::cout << "partial_seq: " << partial_seq << "\n";
            // Total weighted flow time for jobs in the partial schedule \sigma
            // TODO Maybe include worst-case calculation in this part of the formula
            // Tests only the special case of consective jobs
            // Let \sigma be a partial schedule for the set of jobs S, let s = |S|, and let i and j be distinct jobs
            // in N \ S. If Gw(\sigma.ij) - Gw(\sigma.ji) >=
            //              sum(r \in U - {i, j}, w_r) x max(1<=k<=m, Lk(\sigma.ji) - Lk(\sigma.ij))
            // then \sigma.ji dominates \sigma.ij.

            Permutation sigma_ij(partial_seq);
            sigma_ij.push_back(0);
            //sigma_ij.push_back(0);
            Permutation sigma_ji(sigma_ij);
            for(int i : unscheduled_job_list) {
                sum_w_U += I.w(i);
            }
            // If the permutation is empty, return false
            if(s == 0)  return false;

            // Property 1 : For two jobs i and j, where i precedes j
            // for (unsigned pi = 0; pi < partial_seq.size() - 1; pi++) {
                // int i = partial_seq[pi];
                int i = partial_seq[s - 1];
                sum_w_U += I.w(i);
                //for (unsigned pj = pi + 1; pj < partial_seq.size(); pj++) {
                for (int j : unscheduled_job_list) {
                    // int j = partial_seq[pj];
                    // sigma_ij[s] = i;
                    // sigma_ij[s+1] = j;
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+2);
                    //sigma_ji[s] = j;
                    //sigma_ji[s+1] = i;
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+2);
                    //g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+1) + I.w(j) * C_sigma_ij(I.m, s+2);
                    //g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+1) + I.w(i) * C_sigma_ji(I.m, s+2);
                    sigma_ij[s-1] = i;
                    sigma_ij[s+0] = j;
                    //// PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+1);
                    PFSProblem::makespan(C_sigma_ij, sigma_ij, I.m, p_time);
                    sigma_ji[s-1] = j;
                    sigma_ji[s+0] = i;
                    //// PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s);
                    //// PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+1);
                    PFSProblem::makespan(C_sigma_ji, sigma_ji, I.m, p_time);
                    g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+0) + I.w(j) * C_sigma_ij(I.m, s+1);
                    g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+0) + I.w(i) * C_sigma_ji(I.m, s+1);
                    double max_lk_diff = 0.0;
                    for(int k = 1; k <= I.m; k++) {
                        //max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+2) - C_sigma_ij(k, s+2));
                        max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+1) - C_sigma_ij(k, s+1), less_eps);
                    }
                    if(NumericUtil::TestsetIsGE(g_w_sigma_ij - g_w_sigma_ji, (sum_w_U - I.w(i) - I.w(j)) * max_lk_diff)) {
                        return true;  // \sigma.ji dominates \sigma.ij -> violation
                    }
                }
            //}
            return false;  // No violation
        }

        /**
         * Calculate which operations (r, j) will have their processing times oscillated on Dominance_2.
         * First includes as many operations related to last scheduled job partial_seq[end].
         * Then includes as many operations as possible, related to unscheduled job 'ui'. 
         * Total number of oscillated operations follows available gamma budget parameter. 
         * Returns a list of the oscillated operations (r, j) regarding jobs j on machines r.
         * @param unscheduled_job_list
         * @param p_time
         * @return
         */
        std::vector< std::tuple<int,int> > get_oscillated_operations_by_ptime_budget_dom2(const Permutation &partial_seq, 
                                                      const std::vector<int> &unscheduled_job_list,
                                                      const ublas::matrix<double>& p_bar,
                                                      const ublas::matrix<double>& p_hat,
                                                      const unsigned& gamma_value,
                                                      const unsigned& ui,
                                                      const unsigned& m, 
                                                      const bool& ascending = false) {
            // 1. Try to oscillate the maximum number of operations related to job 'ui', according to available budget
            //    Simple greedy algorithm based on largest oscillation values p_hat(r, ui)
            std::vector< std::tuple<int, int> > op_ui;
            op_ui.reserve(m);
            for(int r = 1; r <= m; r++) {
                op_ui.emplace_back(p_bar(r, ui) + p_hat(r, ui), r);
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(op_ui.begin(), op_ui.end(), [ascending](const std::tuple<int, int> &i, const std::tuple<int, int> &j) {
                if(ascending)  return i < j;
                else  return i > j;
            });
            std::vector< std::tuple<int, int> > oscillated_operations;
            stringstream ss;
            int used_budget = 0;
            for (auto i : op_ui){  // select at most gamma operations (k, j) from list op_ui
                if(NumericUtil::TestsetIsGE(used_budget, gamma_value))  break;
                oscillated_operations.emplace_back(std::get<1>(i), ui);
                ss << "operation (" << std::get<1>(i) << ", " << ui << ") -> p_time = " << std::get<0>(i) << ", ";
                used_budget++;
            }
            ss << "TOTAL = " << used_budget << "\n";
            //std::cout << "oscillared_operations: " << ss.str();
            // 2. Try to oscillate operations regarding the last job in partial_seq, if remaining budget allows
            std::vector< std::tuple<int, int> > op_s;
            op_s.reserve(m);
            int j = partial_seq[partial_seq.size() - 1];
            for(int r = 1; r <= m; r++) {
                op_s.emplace_back(r, j);
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(op_s.begin(), op_s.end(), [p_bar, p_hat, ascending](const std::tuple<int, int> &x, const std::tuple<int, int> &y) {
                if(ascending)  return NumericUtil::TestsetIsLT(p_bar(std::get<0>(x), std::get<1>(x)) + p_hat(std::get<0>(x), std::get<1>(x)), 
                    p_bar(std::get<0>(y), std::get<1>(y)) + p_hat(std::get<0>(y), std::get<1>(y)));  // a < b
                else  return NumericUtil::TestsetIsGT(p_bar(std::get<0>(x), std::get<1>(x)) + p_hat(std::get<0>(x), std::get<1>(x)), 
                    p_bar(std::get<0>(y), std::get<1>(y)) + p_hat(std::get<0>(y), std::get<1>(y)));  // a > b
            });
            for (auto i : op_s){  // select at most (gamma - used_budget) operations (r, j) from list op_s
                if(NumericUtil::TestsetIsGE(used_budget, gamma_value))  break;
                oscillated_operations.emplace_back(std::get<0>(i), std::get<1>(i));
                ss << "operation (" << std::get<0>(i) << ", " << std::get<1>(i) << ") -> p_time = " 
                    << (p_bar(std::get<0>(i), std::get<1>(i)) + p_hat(std::get<0>(i), std::get<1>(i))) << ", ";
                used_budget++;
            }
            ss << "TOTAL = " << used_budget << "\n";
            //std::cout << "oscillared_operations: " << ss.str();
            return oscillated_operations;
        }

        /**
        * Procedure Dominance_2 - Dominance_1 with special modification on oscillated jobs in processing time matrix,
        * to be used in the robust budget case, when 0 < Gamma < 100%.
        * Returns true if sequence s violates dominance rules.
        * Based on the paper:
        * Chung, C. S., Flynn, J., & Kirca, O. (2002). A branch and bound algorithm to minimize the total flow time
        * for m-machine permutation flowshop problems. International Journal of Production Economics, 79(3), 185-196.
        * @return true if sequence s violates dominance rules.
        */
        bool dominance_2(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                           Permutation partial_seq, const std::vector<int> &unscheduled_job_list) {
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            ublas::matrix<double> p_time_low = I.get_p_bar();  // I.get_oscillated_p_matrix();
            ublas::matrix<double> p_time_high = I.get_p_bar();  // I.get_oscillated_p_matrix();
            const ublas::matrix<double>& p_bar = I.get_p_bar();
            const ublas::matrix<double>& p_hat = I.get_p_hat();
            double gamma = ((I.m * I.n * rob_params.budget_gamma[0]) / 100.0);
            int max_gamma = (int)floor(gamma);
            double g_w_sigma_ij = 0.0;
            double g_w_sigma_ji = 0.0;
            double sum_w_U = 0.0;
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };

            // FIXME Remove this test code            
            int src[] = { 6, 2, 3, 5, 10, 7, 8 };
            int n = sizeof(src) / sizeof(src[0]);
            std::vector<int> v(src, src + n);
            // partial_seq = v;
            // s = partial_seq.size();
            // =============================


            // Completion time matrix of partial schedule \sigma AT NOMINAL PROCESSING TIMES P_bar
            boost::numeric::ublas::matrix<double> C_sigma_ij_low = boost::numeric::ublas::zero_matrix<double>(I.m + 1, I.n + 2);
            PFSProblem::makespan(C_sigma_ij_low, partial_seq, I.m, p_bar);
            boost::numeric::ublas::matrix<double> C_sigma_ji_low(C_sigma_ij_low);
            // std::cout << "C_matrix: " << C << "\n";
            // std::cout << "partial_seq: " << partial_seq << "\n";
            // Total weighted flow time for jobs in the partial schedule \sigma
            // TODO Maybe include worst-case calculation in this part of the formula
            // Tests only the special case of consective jobs
            // Let \sigma be a partial schedule for the set of jobs S, let s = |S|, and let i and j be distinct jobs
            // in N \ 􏰀S. If Gw(\sigma.ij) - Gw(\sigma.ji) >=
            //              sum(r \in U - {i, j}, w_r) x max(1<=k<=m, Lk(\sigma.ji) - Lk(\sigma.ij))
            // then \sigma.ji dominates \sigma.ij.

            Permutation sigma_ij(partial_seq);
            sigma_ij.push_back(0);
            //sigma_ij.push_back(0);
            Permutation sigma_ji(sigma_ij);
            for(int i : unscheduled_job_list) {
                sum_w_U += I.w(i);
            }

            // Property 1 : For two jobs i and j, where i precedes j
            // for (unsigned pi = 0; pi < partial_seq.size() - 1; pi++) {
                // int i = partial_seq[pi];
                int i = partial_seq[s - 1];
                sum_w_U += I.w(i);
                // update p_time matrix with oscillated operations O_(r, j)
                std::vector< std::tuple<int,int> > oscillated_operations_high = get_oscillated_operations_by_ptime_budget_dom2(
                     partial_seq, unscheduled_job_list, p_bar, p_hat, max_gamma, i, I.m, false);
                for(std::tuple<int,int> op : oscillated_operations_high) {
                    p_time_high(std::get<0>(op), std::get<1>(op)) += p_hat(std::get<0>(op), std::get<1>(op));
                }
                std::vector< std::tuple<int,int> > oscillated_operations_low = get_oscillated_operations_by_ptime_budget_dom2(
                     partial_seq, unscheduled_job_list, p_bar, p_hat, max_gamma, i, I.m, true);
                //for(std::tuple<int,int> op : oscillated_operations_low) {
                //    p_time_low(std::get<0>(op), std::get<1>(op)) += p_hat(std::get<0>(op), std::get<1>(op));
                //}
                // Completion time matrix of partial schedule \sigma WITH oscillated processing times
                boost::numeric::ublas::matrix<double> C_sigma_ij = boost::numeric::ublas::zero_matrix<double>(I.m + 1, I.n + 2);
                PFSProblem::makespan(C_sigma_ij, partial_seq, I.m, p_time_high);
                boost::numeric::ublas::matrix<double> C_sigma_ji(C_sigma_ij);
                //for (unsigned pj = pi + 1; pj < partial_seq.size(); pj++) {
                for (int j : unscheduled_job_list) {
                    // int j = partial_seq[pj];
                    // sigma_ij[s] = i;
                    // sigma_ij[s+1] = j;
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+2);
                    //sigma_ji[s] = j;
                    //sigma_ji[s+1] = i;
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+2);
                    //g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+1) + I.w(j) * C_sigma_ij(I.m, s+2);
                    //g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+1) + I.w(i) * C_sigma_ji(I.m, s+2);
                    sigma_ij[s-1] = i;
                    sigma_ij[s+0] = j;
                    PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time_high, s+1);
                    sigma_ji[s-1] = j;
                    sigma_ji[s+0] = i;
                    PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time_high, s);
                    PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time_high, s+1);

                    /*
                    g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+0) + I.w(j) * C_sigma_ij(I.m, s+1);
                    g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+0) + I.w(i) * C_sigma_ji(I.m, s+1);
                    */
                    PFSProblem::makespan(C_sigma_ij_low, sigma_ij, I.m, p_time_low);
                    PFSProblem::makespan(C_sigma_ji_low, sigma_ji, I.m, p_time_low);
                    // TODO Calculate the 2 values below by nominal on (1) AND oscillating max on (2)
                    g_w_sigma_ij = I.w(i) * C_sigma_ij_low(I.m, s+0) + I.w(j) * C_sigma_ij_low(I.m, s+1);  // (1)
                    g_w_sigma_ji = I.w(j) * C_sigma_ji_low(I.m, s+0) + I.w(i) * C_sigma_ji_low(I.m, s+1);  // (2)
                    //// g_w_sigma_ij = I.w(i) * C_sigma_ij_nominal(I.m, s+0) + I.w(j) * C_sigma_ij_nominal(I.m, s+1);  // (1)
                    //// g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+0) + I.w(i) * C_sigma_ji(I.m, s+1);  // (2)


                    double max_lk_diff = std::numeric_limits<double>::min();
                    for(int k = 1; k <= I.m; k++) {
                        //max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+2) - C_sigma_ij(k, s+2));
                        //// max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+1) - C_sigma_ij(k, s+1));
                        max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+1) - C_sigma_ij_low(k, s+1), less_eps);
                    }
                    // Try to make (g_w_sigma_ij - g_w_sigma_ji) as small as possible, by assuming nominal processing times
                    if(NumericUtil::TestsetIsGE(g_w_sigma_ij - g_w_sigma_ji, (sum_w_U - I.w(i) - I.w(j)) * max_lk_diff)) {
                        // return true;  // \sigma.ji dominates \sigma.ij -> violation
                        /**
                        if(v == partial_seq) {
                            std::cerr << "****** Dominance_2: + g_w_sigma_ij = " << g_w_sigma_ij << "\n";
                            std::cerr << "****** Dominance_2: - g_w_sigma_ji = " << g_w_sigma_ji << "\n";

                            std::cerr << "****** Dominance_2: > max_lk_diff = " << max_lk_diff << "\n";
                            std::cerr << "****** Dominance_2: * (sum_w_U = " << sum_w_U << "\n";
                            std::cerr << "****** Dominance_2: - I.w(i) = " << I.w(i) << "\n";
                            std::cerr << "****** Dominance_2: - I.w(j) ) = " << I.w(j) << "\n";
                        } */

                        return true;
                    }
                }
            //}
            return false;  // No violation
        }

        /**
         * Multi-scenario dominance rule for the PFSP-TWCT Master Problem (C&CG algorithm).
         * Returns true (prune current node), if dominance rule applies to given all scenarios / processing time matrices 
         * in P_time_scenario_list.
         * Besides, P_time_scenario_list, this method also uses pre-calculated l_kl matrices, stored in l_kl_scenario_list.
         */
        bool dominance_multi_scenario(RobPFSInstance_WCT &I, 
                                    Permutation &partial_seq, const std::vector<int> &unscheduled_job_list,
                                    std::vector< ublas::matrix<double> > &P_time_scenario_list,
                                    std::vector< boost::multi_array<double, 3> > &l_kl_scenario_list) {
            auto less_eps = [&](const double& a, const double& b) {
                return NumericUtil::TestsetIsLT(a, b);
            };
            // If the partial permutation is empty, return false
            if(partial_seq.size() == 0) {
                return false;
            }
            unsigned node_level = partial_seq.size();
            unsigned s = node_level;
            std::vector<bool> dominance_list;
            for(int scenario = 0; scenario < P_time_scenario_list.size(); scenario++) {
                const ublas::matrix<double>& p_time = P_time_scenario_list[scenario];
                const boost::multi_array<double, 3>& q = l_kl_scenario_list[scenario];
                double g_w_sigma_ij = 0.0;
                double g_w_sigma_ji = 0.0;
                double sum_w_U = 0.0;
                // Completion time matrix of partial schedule \sigma
                boost::numeric::ublas::matrix<double> C_sigma_ij = boost::numeric::ublas::zero_matrix<double>(I.m + 1, I.n + 2);
                PFSProblem::makespan(C_sigma_ij, partial_seq, I.m, p_time);
                boost::numeric::ublas::matrix<double> C_sigma_ji(C_sigma_ij);
                // Tests only the special case of consective jobs
                // Let \sigma be a partial schedule for the set of jobs S, let s = |S|, and let i and j be distinct jobs
                // in N \ S. If Gw(\sigma.ij) - Gw(\sigma.ji) >=
                //              sum(r \in U - {i, j}, w_r) x max(1<=k<=m, Lk(\sigma.ji) - Lk(\sigma.ij))
                // then \sigma.ji dominates \sigma.ij.
                Permutation sigma_ij(partial_seq);
                sigma_ij.push_back(0);
                Permutation sigma_ji(sigma_ij);
                for(int i : unscheduled_job_list) {
                    sum_w_U += I.w(i);
                }  
                // Property 1 : For two jobs i and j, where i precedes j
                int i = partial_seq[s - 1];
                sum_w_U += I.w(i);
                bool scenario_dominated = false;
                for (int j : unscheduled_job_list) {
                    // int j = partial_seq[pj];
                    // sigma_ij[s] = i;
                    // sigma_ij[s+1] = j;
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+2);
                    //sigma_ji[s] = j;
                    //sigma_ji[s+1] = i;
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+1);
                    //PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+2);
                    //g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+1) + I.w(j) * C_sigma_ij(I.m, s+2);
                    //g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+1) + I.w(i) * C_sigma_ji(I.m, s+2);
                    sigma_ij[s-1] = i;
                    sigma_ij[s+0] = j;
                    //// PFSProblem::incremental_makespan(C_sigma_ij, sigma_ij, I.m, p_time, s+1);
                    PFSProblem::makespan(C_sigma_ij, sigma_ij, I.m, p_time);
                    sigma_ji[s-1] = j;
                    sigma_ji[s+0] = i;
                    //// PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s);
                    //// PFSProblem::incremental_makespan(C_sigma_ji, sigma_ji, I.m, p_time, s+1);
                    PFSProblem::makespan(C_sigma_ji, sigma_ji, I.m, p_time);
                    g_w_sigma_ij = I.w(i) * C_sigma_ij(I.m, s+0) + I.w(j) * C_sigma_ij(I.m, s+1);
                    g_w_sigma_ji = I.w(j) * C_sigma_ji(I.m, s+0) + I.w(i) * C_sigma_ji(I.m, s+1);
                    double max_lk_diff = 0.0;
                    for(int k = 1; k <= I.m; k++) {
                        //max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+2) - C_sigma_ij(k, s+2));
                        max_lk_diff = std::max(max_lk_diff, C_sigma_ji(k, s+1) - C_sigma_ij(k, s+1), less_eps);
                    }
                    if(NumericUtil::TestsetIsGE(g_w_sigma_ij - g_w_sigma_ji, (sum_w_U - I.w(i) - I.w(j)) * max_lk_diff)) {
                        scenario_dominated = true;  // \sigma.ji dominates \sigma.ij -> violation
                        break;
                    }
                }
                dominance_list.push_back(scenario_dominated);  // No violation
            }
            bool all_scenarios_dominated = true;
            for(bool dom : dominance_list) {  // dominance rule must apply to all scenarios to be valid
                all_scenarios_dominated &= dom;
            }
            return all_scenarios_dominated;
        }

        /**
         * LB_3 : Pseudo-random lower bound, based on a spinning wheel.
         * TODO usar a roleta e memoria nas etapas construtivas do algoritmo GRASP / IGA
         * Uses a memory (function parameter) that, for each job j, stores the sum of lower bound values and the
         * total number of times the job j was chosen to assume its worst-case proc. time value.
         * For each job j ==> [ sum(LB) ; sum(count) ]  ( LB memory )
         * The memory in inherited from the parent B & B node.
         * @return
         */
        double randomized_lb(RobPFSInstance_WCT &I, Time &ub, const Robust_Parameters &rob_params,
                                    const Permutation &partial_seq,
                                    const std::vector<int> &unscheduled_job_list,
                                    boost::numeric::ublas::vector<double> &mem_lb_sum,
                                    boost::numeric::ublas::vector<unsigned> &mem_lb_count) {
            // Beta dist (5,2) - left-skewed distribution (most values on the right)
            BETA_DISTRIBUTION dist(5, 2);
            double x = dist(rng);  // spinning wheel: sample a value in the [0, 1] range using beta dist.

            // For each unassigned job j, calculate its absolute grade value (count x LB)
            double denom = 0;
            boost::numeric::ublas::vector<double> grade(I.n + 1, 0);
            for(int j : unscheduled_job_list) {
                grade(j) = mem_lb_count(j) * mem_lb_sum(j);
                denom += grade(j);
            }
            // For each unassigned job j, normalize its grade value (count x LB) between 0 and 1
            for(int j : unscheduled_job_list) {
                grade(j) /= denom;
            }
            // select the job whose normalized grade is the closest to the sampled value x
            int which_job = unscheduled_job_list[0];
            double min_distance = abs(x - grade(which_job));
            for(int j : unscheduled_job_list) {
                double distance = abs(x - grade(j));
                if(distance < min_distance) {
                    min_distance = distance;
                    which_job = j;
                }
            }
            // select the job 'which_job' to oscillate to its worst value
            double lb = 0;


            // save the lb results for the future (memory)
            mem_lb_count(which_job)++;
            mem_lb_sum(which_job) += lb;
            return lb;
        }
    };
}

#endif //PFSP_PFSP_WCT_BUDGET_LOWERBOUND_H
