//
// Created by mlevorato on 11/28/19.
//

#ifndef FLOWSHOP_SOLVER_GRASPSOLVER_H
#define FLOWSHOP_SOLVER_GRASPSOLVER_H

#include <glog/logging.h>
#include <boost/timer/timer.hpp>
#include <string>
#include <chrono>

#include "../PFSProblem.h"
#include "../../../util/random.h"
#include "../PFSInstance.h"
#include "../../PFSSolution.h"
#include "FileUtil.h"
#include "../../PFSP_Parameters.h"
#include "../../metaheuristic/grasp/GRASP.h"
#include "../../metaheuristic/common/OutputUtil.h"

namespace problem {
    namespace pfsp {

        using boost::timer::cpu_timer;
        using boost::timer::cpu_times;
        using boost::timer::nanosecond_type;
        using namespace parameters;
        using namespace std;
        using namespace util;

        class GRASPSolver_PFSP_Cmax {
        public:
            problem::common::PFSSolution solve(PFSInstance &instance, double lb, const PFSP_Parameters &pfsp_params,
                                               const UpperBound_Parameters &upper_bound_params) {
                Inputs<PFSInstance> aInputs(instance.n, instance.m, instance);
                aInputs.initialize();
                // Setup algorithm
                string distribution("uniform");
                unsigned timefactor = upper_bound_params.grasp_tl;  // time per operations (ms) for IGA (default: 30)
                unsigned time_limit = (instance.n * instance.m * timefactor + 999) / 1000;

                Test aTest(instance.name, time_limit, upper_bound_params.grasp_maxiter,
                           distribution, upper_bound_params.beta1, upper_bound_params.beta2, upper_bound_params.seed,
                           upper_bound_params.vnd_permutation,
                           upper_bound_params.vnd_size, upper_bound_params.first_improvement,
                           upper_bound_params.random_vnd, upper_bound_params.adaptive_construction,
                           upper_bound_params.obj_update_freq);
                try {
                    setupRandom(upper_bound_params.seed);

                    // Invoke algorithm - 100000 iterations
                    GRASP<PFSInstance> grasp(aTest, aInputs);
                    grasp.init();
                    Outputs<PFSInstance> output = grasp.run();

                    SchedulingSolution<PFSInstance> grasp_sol = output.getOurBestSol();
                    problem::common::PFSSolution neh_solution = convert_from_grasp_solution(output.getNehSol());
                    problem::common::PFSSolution solution = convert_from_grasp_solution(grasp_sol);
                    solution.time_spent = output.get_time_spent();
                    // validate the objective function value
                    double cmax = PFSProblem::calculateCmax(solution.permutation, instance);
                    if (fabs(cmax - grasp_sol.getCosts()) > 1e-4) {
                        LOG(ERROR) << "[GRASP] Objective function error. Should be = " << cmax;
                    }

                    OutputUtil::print_deterministic_solution_summary(instance, neh_solution.getPermutation(),
                                                                     output.getNehSol().getCosts(),
                                                                     output.getNehSol().getTime(),
                                                                     solution, output.get_time_spent(),
                                                                     output.getOurBestSol().getTime(),
                                                                     output.get_number_iterations(),
                                                                     output.get_num_visited_solutions(),
                                                                     output.get_number_improvements(),
                                                                     output.get_time_results(),
                                                                     pfsp_params.outputFolder, pfsp_params.executionId,
                                                                     "GRASP", upper_bound_params,
                                                                     output.getAlphaHistory(),
                                                                     output.getAlphaQuality(),
                                                                     output.getAlphaProbability());

                    LOG(INFO)
                        << "========================================================================";
                    LOG(INFO) << "              GRASP  FINISHED";
                    LOG(INFO)
                        << "========================================================================";
                    LOG(INFO) << "[GRASPSolver_PFSP_Cmax] Best solution value  = "
                                            << grasp_sol.getCosts();
                    LOG(INFO) << "[GRASPSolver_PFSP_Cmax] Time spent: " << output.get_time_spent() << " s";
                    std::cout << " GRASP=" << grasp_sol.getCosts() << " ; Time(s) = " << output.get_time_spent() << "\n";
                    LOG(INFO)
                        << "========================================================================";

                    return solution;
                } catch (const std::exception &e) {
                    std::cout << "GRASP: an unexpected exception has been caught : " << e.what() << std::endl;
                } catch (...) {
                    std::cout << "GRASP: an unexpected exception has been caught : " << std::endl;
                }
                return problem::common::PFSSolution();
            }

            problem::common::PFSSolution convert_from_grasp_solution(const SchedulingSolution<PFSInstance> &grasp_solution) {
                std::vector<Job> jobs = grasp_solution.getJobs(); // array of jobs in this solution
                Permutation perm;
                // job id starts at 1
                for (int i = 0; i < jobs.size(); i++) {
                    perm.push_back(jobs[i].getId());
                }
                problem::common::PFSSolution pfsp_solution(perm, grasp_solution.getCosts(), grasp_solution.getNJobs());
                return pfsp_solution;
            }
        };
    }
}

#endif //FLOWSHOP_SOLVER_GRASPSOLVER_H
