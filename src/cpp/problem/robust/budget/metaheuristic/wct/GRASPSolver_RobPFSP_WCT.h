//
// Created by mlevorato on 5/1/20.
//

#ifndef PFSP_GRASPSOLVER_ROBPFSP_WCT_H
#define PFSP_GRASPSOLVER_ROBPFSP_WCT_H

#include <glog/logging.h>
#include <boost/timer/timer.hpp>
#include <string>
#include <chrono>

#include "../../../../../util/random.h"
#include "../../../RobPFSInstance_WCT.h"
#include "../../../../PFSSolution.h"
#include "FileUtil.h"
#include "../../../../PFSP_Parameters.h"
#include "../../PFSPBudgetScenario.h"
#include "../../../../metaheuristic/grasp/GRASP.h"
#include "../../../../metaheuristic/common/Inputs.h"
#include "../../../../metaheuristic/common/OutputUtil.h"

namespace pfsp {

    using boost::timer::cpu_timer;
    using boost::timer::cpu_times;
    using boost::timer::nanosecond_type;
    using namespace parameters;
    using namespace std;
    using namespace robust;
    using namespace problem::common;
    using namespace util;

    class GRASPSolver_RobPFSP_WCT {
    public:
        problem::common::PFSSolution solve(RobPFSInstance_WCT &instance, double lb, const PFSP_Parameters &pfsp_params,
                                           const UpperBound_Parameters &upper_bound_params) {
            Inputs<RobPFSInstance_WCT> aInputs(instance.n, instance.m, instance);
            aInputs.initialize();
            // Setup algorithm
            string distribution("uniform");
            unsigned timefactor = upper_bound_params.grasp_tl;  // time per operations (ms) (default: 30)
            unsigned time_limit = (instance.n*instance.m*timefactor+999)/1000;
            LOG(INFO) << "[GRASPSolver_RobPFSP_WCT] Time factor  = " << timefactor;
            LOG(INFO) << "[GRASPSolver_RobPFSP_WCT] Time limit   = " << time_limit;

            Test aTest(instance.name, time_limit, upper_bound_params.grasp_maxiter,
                       distribution, upper_bound_params.beta1, upper_bound_params.beta2, upper_bound_params.seed,
                       upper_bound_params.vnd_permutation,
                       upper_bound_params.vnd_size, upper_bound_params.first_improvement,
                       upper_bound_params.random_vnd, upper_bound_params.adaptive_construction,
                       upper_bound_params.obj_update_freq);
            try {
                setupRandom(upper_bound_params.seed);

                // Invoke algorithm - 100000 iterations
                GRASP<RobPFSInstance_WCT> grasp(aTest, aInputs, false);  // disabled incremental updates for the obj func.
                LOG(INFO) << "[Robust PFSP Budget] GRASP initialize...";
                grasp.init();
                LOG(INFO) << "[Robust PFSP Budget] GRASP solve...";
                Outputs<RobPFSInstance_WCT> output = grasp.run();
                problem::common::PFSSolution neh_solution = convert_from_grasp_solution(output.getNehSol());
                SchedulingSolution<RobPFSInstance_WCT> grasp_sol = output.getOurBestSol();
                problem::common::PFSSolution solution = convert_from_grasp_solution(grasp_sol);
                solution.job_list = grasp_sol.getJobs();
                solution.time_spent = output.get_time_spent();
                // TODO Obtain worst-case scenario from GRASP solution
                std::vector<Permutation> p_pi(2, solution.permutation);
                PFSPBudgetScenario s_worst(instance.rob_params.budget_gamma, p_pi);
                if(upper_bound_params.validate_obj) {
                    LOG(INFO) << "[Robust PFSP Budget GRASP] Validating objective function...";
                    s_worst = robust::PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(instance.n,
                                                                                    instance.p_bar, instance.p_hat,
                                                                                    instance.rob_params,
                                                                                    grasp_sol.getJobs(),
                                                                                    instance.w, true, true, 0, false);
                    double wct = s_worst.value;
                    if (fabs(wct - solution.value) > 1e-4) {
                        LOG(ERROR) << "[GRASP] Objective function error. Obj = " << solution.value
                                                 << ". Should be = " << wct;
                        cout << "[GRASP] Objective function error. Obj = " << solution.value << ". Should be = " << wct
                             << "\n";
                        solution.value = wct;
                    }
                } else {
                    LOG(INFO) << "[GRASP] GRASP obj validation skipped.";
                }

                OutputUtil::print_robust_solution_summary(&instance, neh_solution.getPermutation(),
                                  output.getNehSol().getCosts(), output.getNehSol().getTime(),
                                  s_worst, solution, output.get_time_spent(),
                                  output.getOurBestSol().getTime(), output.get_number_iterations(),
                                  output.get_num_visited_solutions(), output.get_number_improvements(),
                                  output.get_time_results(), pfsp_params.outputFolder, pfsp_params.executionId, "GRASP",
                                  upper_bound_params, output.getAlphaHistory(),
                                  output.getAlphaQuality(), output.getAlphaProbability());


                LOG(INFO) << "========================================================================";
                LOG(INFO) << "              ROBUST  GRASP  FINISHED";
                LOG(INFO) << "========================================================================";
                LOG(INFO) << "[GRASPSolver_RobPFSP_WCT] Best solution value  = " << grasp_sol.getCosts();
                LOG(INFO) << "[GRASPSolver_RobPFSP_WCT] Time spent: " << output.get_time_spent() << " s";
                std::cout << " GRASP=" << grasp_sol.getCosts() << " ; Time(s) = " << output.get_time_spent() << "\n";
                LOG(INFO) << "========================================================================";

                return solution;
            } catch (const std::exception& e) {
                std::cerr << "GRASP: an unexpected exception has been caught : " << e.what() << std::endl;
            }catch (...)
            {
                std::cerr << "GRASP: an unexpected exception has been caught : " << std::endl;
            }
            std::cerr << "[ERROR] GRASP is returning an invalid solution! " << std::endl;
            return problem::common::PFSSolution();
        }

        problem::common::PFSSolution convert_from_grasp_solution(const SchedulingSolution<RobPFSInstance_WCT> &grasp_solution) {
            std::vector<Job> jobs = grasp_solution.getJobs(); // array of jobs in this solution
            Permutation perm;
            // job id starts at 1
            for(int i = 0; i < jobs.size(); i++) {
                perm.push_back(jobs[i].getId());
            }
            problem::common::PFSSolution pfsp_solution(perm, grasp_solution.getCosts(), grasp_solution.getNJobs());
            return pfsp_solution;
        }
    };
}

#endif //PFSP_GRASPSOLVER_ROBPFSP_WCT_H
