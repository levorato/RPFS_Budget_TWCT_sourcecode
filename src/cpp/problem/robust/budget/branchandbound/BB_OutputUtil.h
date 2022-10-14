#ifndef PFSP_BB_ROB_OUTPUTUTIL_H
#define PFSP_BB_ROB_OUTPUTUTIL_H

#include "../../../deterministic/PFSInstance.h"
#include "../../../PFSSolution.h"
#include "../../RobPFSInstance.h"
#include "../PFSPBudgetScenario.h"
#include "../../../../util/include/FileUtil.h"
#include "../../../PFSP_Parameters.h"

#include <string>

namespace util {

    using namespace std;
    using namespace problem::common;
    using namespace robust;
    using namespace parameters;

    class BB_OutputUtil {
    public:
        static void print_deterministic_solution_summary(const PFSInstance &instance,
                                                 PFSSolution &solution, 
                                                 const string &base_folder, const string &executionId,
                                                 const string &method_name,
                                                 const UpperBound_Parameters &upper_bound_params) {
            stringstream sol_csv;
            double time_spent = solution.time_spent;
            long iterations = solution.num_iterations;
            long num_improvements = solution.num_improvements;
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "            B R A N C H    A N D     B O U N D     D O N E ";
            LOG(INFO) << "===========================================================================";
            std::cout << "\n=========================================================================\n";
            std::cout << "            B R A N C H    A N D     B O U N D     D O N E \n";
            std::cout << "=========================================================================\n";
            std::cout << "Solution: " << std::fixed << std::setprecision(10) << solution << "\n";
            std::cout << "Objective function value: " << solution.value << "\n";
            std::cout << "Time Spent (s): " << time_spent << "\n";
            std::cout << "Number of improvements: " << num_improvements << "\n";
            std::cout << "gap: " << solution.gap << "\n";
            std::cout << "is_optimal: " << solution.is_optimal << "\n";
            std::cout << "validated: " << solution.validated << "\n";
            std::cout << "execution_id: " << executionId << "\n";
            std::cout << "seed: " << upper_bound_params.seed << "\n";
            std::cout << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            LOG(INFO) << "[" << method_name << "] Solution: " << std::fixed << std::setprecision(10) << solution;
            LOG(INFO) << "[" << method_name << "] Objective function value: " << solution.value;
            LOG(INFO) << "[" << method_name << "] Time Spent (s): " << std::fixed << std::setprecision(2) << time_spent;
            LOG(INFO) << "[" << method_name << "] Iterations: " << iterations << "\n";
            LOG(INFO) << "[" << method_name << "] Number of improvements: " << num_improvements << "\n";
            LOG(INFO) << "[" << method_name << "] execution_id: " << executionId << "\n";
            LOG(INFO) << "[" << method_name << "] seed: " << upper_bound_params.seed << "\n";
            std::cout << "-----------------------------------------------------------------------------" << "\n";

            string exp_folder = util::FileUtil::create_output_folder(base_folder, "results");
            exp_folder = util::FileUtil::create_output_folder(exp_folder, executionId);
            // append csv result to all experiments csv output file
            boost::filesystem::path all_results_file_path(base_folder);
            all_results_file_path /= method_name + string("-bb_results.csv");
            if (!boost::filesystem::exists(all_results_file_path.string())) {
                // print CSV header to file
                sol_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "n, m"
                        << ", ";
                sol_csv << "solution_value, permutation, time_spent, iterations, ";
                sol_csv
                        << "num_improvements, first_improvement, vnd_size, vnd_permutation, ";
                sol_csv << "random_vnd, adaptive_construction, const_beta1, const_beta2, time_factor\n";
            }
            sol_csv << executionId << ", " << upper_bound_params.seed << ", " << method_name << ", "
                    << instance.name << ", "
                    << instance.n << ", " << instance.m << ", ";
            std::ostringstream vts3;
            std::copy(solution.getPermutation().begin(), solution.getPermutation().end(),
                        std::ostream_iterator<unsigned>(vts3, " "));
            sol_csv << solution.value << ", " << vts3.str() << ", " << time_spent << ", ";
            sol_csv << iterations << ", " << num_improvements << ", ";
            std::ostringstream vts2;
            std::copy(upper_bound_params.vnd_permutation.begin(), upper_bound_params.vnd_permutation.end(),
                        std::ostream_iterator<unsigned>(vts2, " "));
            sol_csv << upper_bound_params.first_improvement << ", " << upper_bound_params.vnd_size << ", "
                    << vts2.str() << ", " << upper_bound_params.random_vnd;
            sol_csv << ", " << upper_bound_params.adaptive_construction << ", " << upper_bound_params.beta1
                    << ", " << upper_bound_params.beta2 << ", " << upper_bound_params.grasp_tl
                    << "\n";
            LOG(INFO) << "[Upper Bound] Appending to all results file in "
                                    << all_results_file_path.string();
            util::FileUtil::saveTextContentsToFile(all_results_file_path.string(), sol_csv.str(), true);
        }

        static void print_robust_solution_summary(const RobPFSInstance *instance,
                                                  PFSSolution &pi_s_worst_star, 
                                                  const string &base_folder, const string &executionId,
                                                  const string &ub_name,
                                                  const UpperBound_Parameters &upper_bound_params) {
            stringstream sol_csv;
            double time_spent = pi_s_worst_star.time_spent;
            long iterations = pi_s_worst_star.num_iterations;
            long num_improvements = pi_s_worst_star.num_improvements;
            stringstream ss_T;
            for (int x = 0; x < instance->rob_params.budget_gamma.size(); x++) {
                ss_T << instance->rob_params.budget_gamma[x] << " ";
            }
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "            B R A N C H    A N D     B O U N D     D O N E ";
            LOG(INFO) << "===========================================================================";
            std::cout << "=========================================================================\n";
            std::cout << "            B R A N C H    A N D     B O U N D     D O N E \n";
            std::cout << "=========================================================================\n";
            std::cout << "Solution (pi_s_worst_star): " << std::fixed << std::setprecision(10) << pi_s_worst_star << "\n";
            std::cout << "Objective function value: " << pi_s_worst_star.value << "\n";
            std::cout << "Time Spent (s): " << time_spent << "\n";
            std::cout << "Iterations: " << iterations << "\n";
            std::cout << "Number of improvements: " << num_improvements << "\n";
            std::cout << "gap: " << pi_s_worst_star.gap << "\n";
            std::cout << "is_optimal: " << pi_s_worst_star.is_optimal << "\n";
            std::cout << "validated: " << pi_s_worst_star.validated << "\n";
            std::cout << "execution_id: " << executionId << "\n";
            std::cout << "seed: " << upper_bound_params.seed << "\n";
            std::cout << "budget-T: " << ss_T.str() << "\n";
            std::cout << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] Solution (pi_s_worst_star): " 
                        << std::fixed << std::setprecision(10) << pi_s_worst_star;
            // LOG(INFO) << "[Robust BB " << ub_name << "] Scenario (s_worst): " << s_worst;
            LOG(INFO) << "[Robust BB " << ub_name << "] Calculating worst solution value...";
            LOG(INFO) << "[Robust BB " << ub_name << "] Objective function value: " << std::fixed 
                        << std::setprecision(10) << pi_s_worst_star.value;
            LOG(INFO) << "[Robust BB " << ub_name << "] Time Spent (s): " << std::fixed << std::setprecision(2) << time_spent;
            LOG(INFO) << "[Robust BB " << ub_name << "] Iterations: " << iterations << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] Number of improvements: " << num_improvements
                                    << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] execution_id: " << executionId << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] seed: " << upper_bound_params.seed << "\n";
            std::cout << "-----------------------------------------------------------------------------" << "\n";

            string exp_folder = util::FileUtil::create_output_folder(base_folder, "results");
            exp_folder = util::FileUtil::create_output_folder(exp_folder, executionId);
            // append csv result to all experiments csv output file
            boost::filesystem::path all_results_file_path(base_folder);
            all_results_file_path /= ub_name + string("-rob-bb_results.csv");
            if (!boost::filesystem::exists(all_results_file_path.string())) {
                // print CSV header to file
                sol_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "alpha, n, m"
                        << ", ";
                sol_csv << "budget_gamma" << ", ";
                sol_csv << "solution_value, permutation, time_spent, iterations, ";
                sol_csv
                        << "num_improvements, first_improvement, vnd_size, vnd_permutation, ";
                sol_csv << "random_vnd, adaptive_construction, const_beta1, const_beta2, time_factor\n";
            }
            // print CSV contents to file
            sol_csv << executionId << ", " << upper_bound_params.seed << ", " << ub_name << ", " << instance->name
                    << ", " << instance->alpha << ", "
                    << instance->n << ", " << instance->m << ", ";
            sol_csv << ss_T.str() << ", " << pi_s_worst_star.value << ", " << pi_s_worst_star.getPermutation()
                    << ", ";
            sol_csv << time_spent << ", " << iterations << ", " << num_improvements << ", ";
            std::ostringstream vts2;
            std::copy(upper_bound_params.vnd_permutation.begin(), upper_bound_params.vnd_permutation.end(),
                        std::ostream_iterator<unsigned>(vts2, " "));
            sol_csv << upper_bound_params.first_improvement << ", " << upper_bound_params.vnd_size << ", ["
                    << vts2.str() << "], " << upper_bound_params.random_vnd;
            sol_csv << ", " << upper_bound_params.adaptive_construction << ", " << upper_bound_params.beta1
                    << ", " << upper_bound_params.beta2 << ", " << upper_bound_params.grasp_tl
                    << "\n";
            LOG(INFO) << "[Robust BB Upper Bound " << ub_name << "] Appending to all results file in "
                                    << all_results_file_path.string();
            util::FileUtil::saveTextContentsToFile(all_results_file_path.string(), sol_csv.str(), true);
        }

    };
}

#endif //PFSP_BB_ROB_OUTPUTUTIL_H
