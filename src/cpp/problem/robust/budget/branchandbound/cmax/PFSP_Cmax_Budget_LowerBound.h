//
// Created by mlevorato on 1/14/20.
//

#ifndef PFSP_PFSP_CMAX_BUDGET_LOWERBOUND_H
#define PFSP_PFSP_CMAX_BUDGET_LOWERBOUND_H

#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>

#include "../../../RobPFSInstance_Cmax.h"
#include "PFSP_Cmax_Budget_NodeData.h"
#include "../../PFSP_Cmax_Budget_WorstCase.h"
#include "../../../../deterministic/heuristic/include/JohnsonSolver.h"
#include "../../../../../util/random.h"
#include "../../../../../util/include/TimeDateUtil.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace problem::pfsp::heuristic;
    using boost::timer::cpu_timer;

    class PFSP_Cmax_Budget_LowerBound {
    public:
        typedef boost::random::beta_distribution<double> BETA_DISTRIBUTION;  // Beta distribution
        // stores statistics about the branch-and-bound process
        long num_leaf_nodes, num_leaf_nodes_infeasible, num_leaf_nodes_dp;
        double total_time_worstcase_mip, total_time_worstcase_dp;

        PFSP_Cmax_Budget_LowerBound() : num_leaf_nodes(0), num_leaf_nodes_infeasible(0), 
            total_time_worstcase_mip(0.0), total_time_worstcase_dp(0.0), num_leaf_nodes_dp(0) {

        }

        double lb(RobPFSInstance_Cmax &I, Time &ub, const Robust_Parameters &rob_params, const Permutation &partial_seq,
                         const std::vector<int> &unscheduled_job_list, PFSP_Cmax_Budget_NodeData &node_data) {
            double lb = lb_1(I, ub, rob_params, partial_seq, unscheduled_job_list);
            /*
            if(lb < ub) {
                double lb_2 = randomized_lb(I, ub, rob_params, partial_seq, unscheduled_job_list,
                        node_data.mem_lb_sum, node_data.mem_lb_count);
                if(lb_2 > lb)
                    lb = lb_2;
            } */
            return lb;
        }

        /**
        * Procedure LB_1 - Lower Bound.
        * Given a partial sequence [ 1 - 2 - ... - l ] (e.g. B & B node):
        * 1. Assign jobs with processing times of the nominal scenario: p_{ij} = p_bar{ij};
        * 2. Arrange the (n - l) positions of the sequence s, by manipulating the unassigned jobs
        *    in Johnson order according to p_bar{ij} and appending them to current partial schedule;
        * 3. Run procedure worst-case based on sequence s.
        * Complexity: O(max(T1, T2) * n^2)
        * @return Lower Bound on the worst-case makespan.
        */
        double lb_1(RobPFSInstance_Cmax &I, Time &ub, const Robust_Parameters &rob_params,
                           const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list) {
            LOG(INFO) << "[PFSPBudgetScenario] Calculating lower bound...";
            unsigned node_level = partial_seq.size();
            boost::timer::cpu_timer timer;
            timer.start();  // Time measurement
            // Arrange the (n - l) unassigned jobs in Johnson order according to p_bar{ij}
            Permutation s = JohnsonSolver::remaining_johnson_order(I.n, partial_seq, I.p_bar);
            double lb = 0.0;
            if(node_level == 0) {
                lb = PFSP_Cmax_Budget_WorstCase::worst_case_cmax_dp(I.n, I.p_bar, I.p_hat, rob_params, s);
            } else {
                lb = PFSP_Cmax_Budget_WorstCase::worst_case_heuristic(I.n, I.p_bar, I.p_hat, rob_params, s);
            }
            total_time_worstcase_dp += TimeDateUtil::calculate_time_spent(timer);
            if(partial_seq.size() == I.n) {  // leaf node
                num_leaf_nodes++;
                if(lb >= ub) {
                    num_leaf_nodes_dp++;
                }
            }
            return lb;
        }

        /**
         * LB_2 : Pseudo-random lower bound, based on a spinning wheel.
         * TODO usar a roleta e memoria nas etapas construtivas do algoritmo grasp / IGA
         * Uses a memory (function parameter) that, for each job j, stores the sum of lower bound values and the
         * total number of times the job j was chosen to assume its worst-case proc. time value.
         * For each job j ==> [ sum(LB) ; sum(count) ]  ( LB memory )
         * The memory in inherited from the parent B & B node.
         * @return
         */
        double randomized_lb(RobPFSInstance_Cmax &I, Time &ub, const Robust_Parameters &rob_params,
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

#endif //PFSP_PFSP_CMAX_BUDGET_LOWERBOUND_H
