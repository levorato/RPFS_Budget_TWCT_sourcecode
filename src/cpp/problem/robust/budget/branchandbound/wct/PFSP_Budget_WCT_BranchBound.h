//
// Created by mlevorato on 3/1/21.
//

#ifndef FLOWSHOP_SOLVER_PFSP_WCT_BUDGET_CONTINUOUS_H
#define FLOWSHOP_SOLVER_PFSP_WCT_BUDGET_CONTINUOUS_H

#include <cassert>
#include <algorithm>

#include "../../../../deterministic/PFSProblem.h"
#include "../../../../PFSSolution.h"
#include "../../../RobPFSInstance_WCT.h"
#include "../../../../PFSP_Parameters.h"
#include "../../PFSPBudgetScenario.h"
#include "PFSP_WCT_Budget_NodeData.h"
#include "PFSP_WCT_Budget_LowerBound.h"
#include "../../metaheuristic/wct/GRASPSolver_RobPFSP_WCT.h"
#include "../../PFSP_WCT_Budget_WorstCase.h"
#include "../../../RobustBranchBound.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    /**
     * Combinatorial branch-and-bound algorithm for the m-machine robust permutation flowshop problem
     * with budget constraints, to minimize the restricted worst-case weighted completion time (min max TWCT).
     */
    class PFSP_Budget_WCT_BranchBound {
    public:
        PFSP_Budget_WCT_BranchBound() : lower_bound() {

        }

        void initialize(RobPFSInstance_WCT &instance, const Robust_Parameters &rob_params) {

        }

        /**
         * Dominance Rules - Continuous time intervals. Returns true if sequence s violates dominance rules.
         */
        bool dominance_rules(RobPFSInstance_WCT &instance, Time &ub, const Robust_Parameters &rob_params,
                                    Permutation &partial_seq, const std::vector<int> &unscheduled_job_list,
                                    PFSP_WCT_Budget_NodeData &node_data) {
            return lower_bound.dominance_rules(instance, ub, rob_params, partial_seq, unscheduled_job_list);
        }

        Time updatelb(RobPFSInstance_WCT &instance, Time &ub, const Robust_Parameters &rob_params,
                      const std::vector<int> &unscheduled_job_list, PFSP_WCT_Budget_NodeData &node_data) {
            Permutation partial_seq;
            return updatelb(instance, ub, rob_params, partial_seq, unscheduled_job_list, node_data);
        }

        Time updatelb(RobPFSInstance_WCT &instance, Time &ub, const Robust_Parameters &rob_params,
                      const Permutation &partial_seq, const std::vector<int> &unscheduled_job_list,
                      PFSP_WCT_Budget_NodeData &node_data) {
            return lower_bound.lb(instance, ub, rob_params, partial_seq, unscheduled_job_list, node_data);
        }

        PFSSolution obtainInitialSolution(RobPFSInstance_WCT &instance, const Time &lb, const PFSP_Parameters &pfsp_params,
                                          const Robust_Parameters &rob_params, const UpperBound_Parameters &upper_bound_params) {
            pfsp::GRASPSolver_RobPFSP_WCT grasp_wct;
            LOG(INFO) << "[Robust PFSP Budget] Solving Robust PFSP WCT instance with GRASP...";
            instance.robust_type = robust::RobustType::budget;
            UpperBound_Parameters new_ub_params(upper_bound_params);
            // disable GRASP objective validation if gamma == 0 or gamma == 100 (saves time)
            if((rob_params.budget_gamma[0] == 0) || (rob_params.budget_gamma[0] == 100)) {
                new_ub_params.validate_obj = false;
                LOG(INFO) << "[Robust PFSP Budget] Disabling Upper Bound obj validation.";
                return grasp_wct.solve(instance, 0, pfsp_params, new_ub_params);
            } else {
                // Solve GRASP for Deterministic PFSP with nominal processing times (Gamma == 0)
                new_ub_params.validate_obj = false;
                RobPFSInstance_WCT modified_instance(instance);
                modified_instance.rob_params.budget_gamma[0] = 0;
                LOG(INFO) << "[Robust PFSP Budget] Solving Deterministic GRASP with Gamma=0...";
                PFSSolution grasp_sol_gamma_0 = grasp_wct.solve(modified_instance, 0, pfsp_params, new_ub_params);
                // Solve GRASP for Deterministic PFSP with maximum processing times (Gamma == 100)
                modified_instance.rob_params.budget_gamma[0] = 100;
                LOG(INFO) << "[Robust PFSP Budget] Solving Deterministic GRASP with Gamma=100...";
                PFSSolution grasp_sol_gamma_100 = grasp_wct.solve(modified_instance, 0, pfsp_params, new_ub_params);
                // Gamma_0: Calculate real objective value (worstcase WCT) based on budget rob_params
                PFSPScenario s_worst = robust::PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(instance.n,
                                                                                    instance.p_bar, instance.p_hat,
                                                                                    rob_params,
                                                                                    grasp_sol_gamma_0.getJobs(),
                                                                                    instance.w, true, true, 0, false);
                grasp_sol_gamma_0.value = s_worst.value;
                // if both permutations are equal, we can save time by calculating the worst-case only once
                if (grasp_sol_gamma_0.permutation == grasp_sol_gamma_100.permutation) {
                    LOG(INFO) << "[Robust PFSP Budget] Final GRASP UB = " << grasp_sol_gamma_0.value;
                    std::cout << "==> Final GRASP UB = " << grasp_sol_gamma_0.value;
                    return grasp_sol_gamma_0;
                }
                // Gamma_100: Calculate real objective value based on budget rob_params
                s_worst = robust::PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(instance.n,
                                                                                    instance.p_bar, instance.p_hat,
                                                                                    rob_params,
                                                                                    grasp_sol_gamma_100.getJobs(),
                                                                                    instance.w, true, true, 0, false);
                grasp_sol_gamma_100.value = s_worst.value;
                if(grasp_sol_gamma_0.value <= grasp_sol_gamma_100.value) {
                    LOG(INFO) << "[Robust PFSP Budget] Final GRASP UB = " << grasp_sol_gamma_0.value;
                    std::cout << "==> Final GRASP UB = " << grasp_sol_gamma_0.value;
                    return grasp_sol_gamma_0;
                } else {
                    LOG(INFO) << "[Robust PFSP Budget] Final GRASP UB = " << grasp_sol_gamma_0.value;
                    std::cout << "==> Final GRASP UB = " << grasp_sol_gamma_100.value;
                    return grasp_sol_gamma_100;
                }
            }
        }

        PFSPScenario calculate_objective(RobPFSInstance_WCT &I, Permutation &s) {
            std::vector<Job> pi;
            for (int i = 0; i < I.n; i++)
                pi.push_back(Job(s[i] - 1, I.m));
            double obj = I.calculate_worst_case_wct(pi, I.n, 0, true);
            PFSPScenario scenario;
            scenario.value = obj;
            return scenario;
        }

        string get_problem_name() {
            return string("PFSP_WCT_Budget");
        }

        void finalize() {
            // robust::PFSP_WCT_Budget_WorstCase::release_cplex_model();
            LOG(INFO) << "[Robust PFSP Budget] TWCT branch-and-bound DONE.";
        }

        PFSP_WCT_Budget_LowerBound get_lower_bound() {
            return lower_bound;
        }

    private:
        PFSP_WCT_Budget_LowerBound lower_bound;

    };

    typedef RobustBranchBound<RobPFSInstance_WCT, PFSP_Budget_WCT_BranchBound, PFSP_WCT_Budget_NodeData> PFSP_WCT_Budget_Continuous_BB;
}

#endif //FLOWSHOP_SOLVER_PFSP_WCT_BUDGET_CONTINUOUS_H
