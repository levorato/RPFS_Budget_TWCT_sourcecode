//
// Created by mlevorato on 1/13/20.
//

#ifndef PFSP_PFSP_CMAX_BUDGET_UPPERBOUND_ALT_H
#define PFSP_PFSP_CMAX_BUDGET_UPPERBOUND_ALT_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <glog/logging.h>
#include "../../../../deterministic/PFSProblem.h"
#include "../../../../PFSSolution.h"
#include "../../../RobPFSInstance_Cmax.h"
#include "../../PFSPBudgetScenario.h"

#include "../../PFSP_Cmax_Budget_WorstCase.h"
#include "GRASPSolver_RobPFSP_Cmax.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace problem::common;

    /**
     * Alternative Robust 2-machine PFSP Budget Upper Bounds, based on Worst-case procedure.
     */
    class PFSP_Cmax_Budget_UpperBound_Alt {
    public:
        /**
         * Calculate the Upper Bound for the robust PFS Problem, according to the method chosen.
         * @param instance
         * @param pfsp_params parameter object which includes the general programa parameters (e.g. output folder)
         * @param rob_params parameter object which includes the method to be used in UB calculation.
         * @return the solution corresponding to the upper bound.
         */
        PFSSolution calculate(RobPFSInstance_Cmax &instance, const PFSP_Parameters &pfsp_params,
                              const Robust_Parameters &rob_params, const UpperBound_Parameters &upper_bound_params) {
            LOG(INFO) << "UB version is: Robust PFSP MiniMax Budget Alt.";
            LOG(INFO) << "[Robust PFSP Budget] Solving Robust PFSP Cmax instance with GRASP...";
            pfsp::GRASPSolver_RobPFSP_Cmax grasp_solver;
            problem::common::PFSSolution grasp_sol = grasp_solver.solve(instance, 0, pfsp_params,
                                                                        upper_bound_params);
            double cmax = PFSP_Cmax_Budget_WorstCase::worst_case_cmax_dp(instance.n, instance.p_bar, instance.p_hat,
                                                                            rob_params, grasp_sol.permutation);

            if (fabs(cmax - grasp_sol.value) > 1e-4) {
                LOG(ERROR) << "[GRASP] Objective function error. Should be = " << cmax;
            }
            LOG(INFO) << "[PFSP Cmax Budget Continuous] Initial solution value  = " << grasp_sol.value;

            // quick test
            /*
            int pi_arr[] = {16, 11, 6, 5, 3, 2, 7, 14, 4, 18, 8, 1, 9, 15, 17, 10, 13, 12, 20, 19};
            std::vector<Job> pi;
            for(int i = 0; i < 20; i++) {
                pi.push_back(Job(pi_arr[i]));
            }
            Robust_Parameters rob_params_2(rob_params);
            rob_params_2.budget_gamma[0] = 12*100/20;
            double cmax_d = PFSP_Cmax_Budget_WorstCase::worst_case_cmax_dp(instance.n, instance.p_bar, instance.p_hat,
                                                                            rob_params_2, pi);
            LOG(INFO) << "[PFSP Cmax Budget Continuous] DEBUG solution value  = " << cmax_d; */

            return grasp_sol;
        }

    };
}

#endif //PFSP_PFSP_CMAX_BUDGET_UPPERBOUND_ALT_H
