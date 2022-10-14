//
// Created by mlevorato on 7/30/19.
//

#ifndef FLOWSHOP_SOLVER_PFSP_CMAX_BUDGET_CONTINUOUS_H
#define FLOWSHOP_SOLVER_PFSP_CMAX_BUDGET_CONTINUOUS_H

#include <cassert>
#include <algorithm>

#include "../../../../deterministic/PFSProblem.h"
#include "../../../../PFSSolution.h"
#include "../../../RobPFSInstance_Cmax.h"
#include "../../../../PFSP_Parameters.h"
#include "../../PFSPBudgetScenario.h"
#include "PFSP_Cmax_Budget_NodeData.h"
#include "PFSP_Cmax_Budget_LowerBound.h"
#include "../../metaheuristic/cmax/PFSP_Cmax_Budget_UpperBound_Alt.h"
#include "../../PFSP_Cmax_Budget_WorstCase.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    /**
     * Combinatorial branch-and-bound algorithm for the 2-machine robust permutation flowshop problem
     * with budget constraints, to minimize the restricted worst-case makespan (min max Cmax).
     */
    class PFSP_Budget_Cmax_BranchBound {
    public:
        PFSP_Budget_Cmax_BranchBound() : lower_bound() {  }

        void initialize(RobPFSInstance_Cmax &instance, const Robust_Parameters &rob_params) {

        }

        /**
         * Dominance Rules - Continuous time intervals. Returns true if sequence s violates dominance rules.
         */
        static bool dominance_rules(RobPFSInstance_Cmax &instance, Time &ub, const Robust_Parameters &rob_params,
                                    Permutation &partial_seq, const std::vector<int> &unscheduled_job_list,
                                    PFSP_Cmax_Budget_NodeData &node_data) {
            return false;
        }

        Time updatelb(RobPFSInstance_Cmax &instance, Time &ub, const Robust_Parameters &rob_params,
                      const std::vector<int> &unscheduled_job_list, PFSP_Cmax_Budget_NodeData &node_data) {
            Permutation partial_seq;
            return updatelb(instance, ub, rob_params, partial_seq, unscheduled_job_list, node_data);
        }

        Time updatelb(RobPFSInstance_Cmax &instance, Time &ub, const Robust_Parameters &rob_params,
                      const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list,
                      PFSP_Cmax_Budget_NodeData &node_data) {
            return lower_bound.lb(instance, ub, rob_params, partial_seq, unscheduled_job_list, node_data);
        }

        PFSSolution obtainInitialSolution(RobPFSInstance_Cmax &instance, const Time &lb, const PFSP_Parameters &pfsp_params,
                                          const Robust_Parameters &rob_params, const UpperBound_Parameters &upper_bound_params) {
            PFSP_Cmax_Budget_UpperBound_Alt budget_alt_ub;
            return budget_alt_ub.calculate(instance, pfsp_params, rob_params, upper_bound_params);
        }

        PFSPScenario calculate_objective(RobPFSInstance_Cmax &I, Permutation &s) {
            std::vector<Job> pi;
            for (int i = 0; i < I.n; i++)
                pi.push_back(Job(s[i] - 1, I.m));
            double obj = I.calcTotalCosts(pi, I.n, 0);
            PFSPScenario scenario;
            scenario.value = obj;
            return scenario;
        }

        string get_problem_name() {
            return string("PFSP_Cmax_Budget");
        }

        void finalize() {
            LOG(INFO) << "[Robust PFSP Budget] Cmax branch-and-bound DONE.";
        }

        PFSP_Cmax_Budget_LowerBound get_lower_bound() {
            return lower_bound;
        }

    private:
        PFSP_Cmax_Budget_LowerBound lower_bound;

    };

    typedef RobustBranchBound<RobPFSInstance_Cmax, PFSP_Budget_Cmax_BranchBound, PFSP_Cmax_Budget_NodeData> PFSP_Cmax_Budget_Continuous_BB;
}

#endif //FLOWSHOP_SOLVER_PFSP_CMAX_BUDGET_CONTINUOUS_H
