//
// Created by mlevorato on 1/30/20.
//

#ifndef PFSP_OUTPUTUTIL_H
#define PFSP_OUTPUTUTIL_H

#include "../../deterministic/PFSInstance.h"
#include "../../PFSSolution.h"
#include "../../robust/RobPFSInstance.h"
#include "../../robust/budget/PFSPBudgetScenario.h"
#include "../../../util/include/FileUtil.h"
#include "../../PFSP_Parameters.h"

#include <string>

namespace util {

    using namespace std;
    using namespace problem::common;
    using namespace robust;
    using namespace parameters;

    class OutputUtil {
    public:
        static void print_deterministic_solution_summary(const PFSInstance &instance,
                                                 const Permutation &neh_permutation,
                                                 const double &neh_value, const double &neh_time,
                                                 PFSSolution &solution, const double &time_spent,
                                                 const double &time_to_best_sol, const long &iterations,
                                                 const long &num_visited_solutions,
                                                 const long &num_improvements,
                                                 const string &timeResults,
                                                 const string &base_folder, const string &executionId,
                                                 const string &method_name,
                                                 const UpperBound_Parameters &upper_bound_params,
                                                 const std::vector<double> &alpha_hist = std::vector<double>(),
                                                 const std::vector< std::vector<double> > &alpha_qual = std::vector< std::vector<double> >(),
                                                 const std::vector< std::vector<double> > &alpha_prob = std::vector< std::vector<double> >()) {
            stringstream solInfo, sol_csv, neh_csv;
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "            U P P E R   B O U N D     D O N E ";
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "[Upper Bound] " << method_name << " finished";
            solInfo << "permutation_neh: " << neh_permutation << "\n";
            solInfo << "cmax_neh: " << neh_value << "\n";
            solInfo << "time_neh: " << neh_time << "\n";
            solInfo << "permutation: " << solution.getPermutation() << "\n";
            solInfo << "cmax: " << solution.value << "\n";
            solInfo << "time: " << time_spent << "\n";
            solInfo << "time to best solution: " << time_to_best_sol << "\n";
            solInfo << "iterations: " << iterations << "\n";
            solInfo << "num_visited_solutions: " << num_visited_solutions << "\n";
            solInfo << "num_improvements: " << num_improvements << "\n";
            solInfo << "execution_id: " << executionId << "\n";
            solInfo << "seed: " << executionId << "\n";
            solInfo << "adaptive_construction: " << upper_bound_params.adaptive_construction << "\n";
            solInfo << "const_beta1: " << upper_bound_params.beta1 << "\n";
            solInfo << "const_beta2: " << upper_bound_params.beta2 << "\n";
            solInfo << "vnd_size: " << upper_bound_params.vnd_size << "\n";
            solInfo << "first_improvement: " << upper_bound_params.first_improvement << "\n";
            solInfo << "random_vnd: " << upper_bound_params.random_vnd << "\n";
            solInfo << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            std::cout << "\n=========================================================================\n";
            std::cout << "            U P P E R   B O U N D     D O N E \n";
            std::cout << "=========================================================================\n";
            std::cout << "[Upper Bound] " << method_name << " finished\n";
            std::cout << "NEH: c=" << neh_value << "; t=" << neh_time << "\n";
            std::cout << "Solution: " << solution << "\n";
            std::cout << "Objective function value: " << solution.value << "\n";
            std::cout << "Time Spent (s): " << time_spent << "\n";
            std::cout << "Time to best solution (s): " << time_to_best_sol << "\n";
            std::cout << "Iterations: " << iterations << "\n";
            std::cout << "Number of visited solutions: " << num_visited_solutions << "\n";
            std::cout << "Number of improvements: " << num_improvements << "\n";
            std::cout << "execution_id: " << executionId << "\n";
            std::cout << "seed: " << upper_bound_params.seed << "\n";
            std::cout << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            LOG(INFO) << "[" << method_name << "] NEH: c=" << neh_value << "; t=" << neh_time;
            LOG(INFO) << "[" << method_name << "] Solution: " << solution;
            LOG(INFO) << "[" << method_name << "] Objective function value: " << solution.value;
            LOG(INFO) << "[" << method_name << "] Time Spent (s): " << time_spent;
            LOG(INFO) << "[" << method_name << "] Time to best solution (s): " << time_to_best_sol
                                    << "\n";
            LOG(INFO) << "[" << method_name << "] Iterations: " << iterations << "\n";
            LOG(INFO) << "[" << method_name << "] Number of visited solutions: " << num_visited_solutions
                                    << "\n";
            LOG(INFO) << "[" << method_name << "] Number of improvements: " << num_improvements << "\n";
            LOG(INFO) << "[" << method_name << "] execution_id: " << executionId << "\n";
            LOG(INFO) << "[" << method_name << "] seed: " << upper_bound_params.seed << "\n";
            std::cout << "-----------------------------------------------------------------------------" << "\n";

            string exp_folder = util::FileUtil::create_output_folder(base_folder, "results");
            exp_folder = util::FileUtil::create_output_folder(exp_folder, executionId);
            // append csv result to all experiments csv output file
            if(! upper_bound_params.parametrize) {
                boost::filesystem::path all_results_file_path(base_folder);
                all_results_file_path /= method_name + string("-all_results.csv");
                if (!boost::filesystem::exists(all_results_file_path.string())) {
                    // print CSV header to file
                    sol_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "n, m"
                            << ", ";
                    sol_csv << "solution_value, permutation, time_spent, time_to_best_sol, iterations, ";
                    sol_csv
                            << "num_visited_solutions, num_improvements, first_improvement, vnd_size, vnd_permutation, ";
                    sol_csv << "random_vnd, adaptive_construction, const_beta1, const_beta2, time_factor\n";
                }
                boost::filesystem::path neh_results_file_path(base_folder);
                neh_results_file_path /= string("NEH-all_results.csv");
                if (!boost::filesystem::exists(neh_results_file_path.string())) {
                    // print CSV header to file
                    neh_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "n, m"
                            << ", solution_value, permutation, time_spent\n";
                }
                // print CSV contents to file
                if(neh_value > 0) {
                    std::ostringstream vts;
                    std::copy(neh_permutation.begin(), neh_permutation.end(),
                              std::ostream_iterator<unsigned>(vts, " "));
                    neh_csv << executionId << ", " << upper_bound_params.seed << ", " << "NEH" << ", "
                            << instance.name << ", "
                            << instance.n << ", " << instance.m << ", "
                            << neh_value << ", " << vts.str() << ", " << neh_time << "\n";
                }
                sol_csv << executionId << ", " << upper_bound_params.seed << ", " << method_name << ", "
                        << instance.name << ", "
                        << instance.n << ", " << instance.m << ", ";
                std::ostringstream vts3;
                std::copy(solution.getPermutation().begin(), solution.getPermutation().end(),
                          std::ostream_iterator<unsigned>(vts3, " "));
                sol_csv << solution.value << ", " << vts3.str() << ", " << time_spent << ", " << time_to_best_sol << ", ";
                sol_csv << iterations << ", " << num_visited_solutions << ", " << num_improvements << ", ";
                std::ostringstream vts2;
                std::copy(upper_bound_params.vnd_permutation.begin(), upper_bound_params.vnd_permutation.end(),
                          std::ostream_iterator<unsigned>(vts2, " "));
                sol_csv << upper_bound_params.first_improvement << ", " << upper_bound_params.vnd_size << ", "
                        << vts2.str() << ", " << upper_bound_params.random_vnd;
                sol_csv << ", " << upper_bound_params.adaptive_construction << ", " << upper_bound_params.beta1
                        << ", " << upper_bound_params.beta2 << ", " << upper_bound_params.grasp_tl
                        << "\n";
                LOG(INFO) << "[Upper Bound] Appending to NEH results file in "
                                        << neh_results_file_path.string();
                util::FileUtil::saveTextContentsToFile(neh_results_file_path.string(), neh_csv.str(), true);
                LOG(INFO) << "[Upper Bound] Appending to all results file in "
                                        << all_results_file_path.string();
                util::FileUtil::saveTextContentsToFile(all_results_file_path.string(), sol_csv.str(), true);
                // save time results output file
                boost::filesystem::path time_result_file_path(exp_folder);
                time_result_file_path /= instance.name + string("-") + method_name + string("-time_result.csv");
                LOG(INFO) << "[Upper Bound " << method_name << "] Saving time results file to "
                                        << time_result_file_path.string();
                util::FileUtil::saveTextContentsToFile(time_result_file_path.string(), timeResults, false);
                // save the other info results file
                boost::filesystem::path sol_info_file_path(exp_folder);
                sol_info_file_path /= instance.name + string("-") + method_name + string("-solution_info.csv");
                LOG(INFO) << "[Upper Bound " << method_name << "] Saving solution info file to "
                                        << sol_info_file_path.string();
                util::FileUtil::saveTextContentsToFile(sol_info_file_path.string(), solInfo.str(), false);
                // Save constructive phase statistics
                if (alpha_hist.size() > 0) {
                    stringstream constInfo;
                    int iter = alpha_hist.size();
                    constInfo << "iter, alpha, probability, quality\n";
                    for (int i = 0; i < iter; i++) {
                        constInfo << i << ", " << alpha_hist[i] << ", ";
                        for (double b : alpha_prob[i]) {
                            constInfo << b << " ";
                        }
                        constInfo << ", ";
                        for (double b : alpha_qual[i]) {
                            constInfo << b << " ";
                        }
                        constInfo << "\n";
                    }
                    boost::filesystem::path sol_const_file_path(exp_folder);
                    sol_const_file_path /= instance.name + string("-") + method_name + string("-constructive_info.csv");
                    LOG(INFO) << "[Robust BB " << method_name << "] Saving constructive info file to "
                                            << sol_const_file_path.string();
                    util::FileUtil::saveTextContentsToFile(sol_const_file_path.string(), constInfo.str(), false);
                }
            }
        }

        static void print_robust_solution_summary(const RobPFSInstance *instance,
                                                  const Permutation &neh_permutation,
                                                  const double &neh_value, const double &neh_time,
                                                  PFSPBudgetScenario &s_worst,
                                                  PFSSolution &pi_s_worst_star, const double &time_spent,
                                                  const double &time_to_best_sol, const long &iterations,
                                                  const long &num_visited_solutions, const long &num_improvements,
                                                  const string &timeResults,
                                                  const string &base_folder, const string &executionId,
                                                  const string &ub_name,
                                                  const UpperBound_Parameters &upper_bound_params,
                                                  const std::vector<double> &alpha_hist = std::vector<double>(),
                                                  const std::vector< std::vector<double> > &alpha_qual = std::vector< std::vector<double> >(),
                                                  const std::vector< std::vector<double> > &alpha_prob = std::vector< std::vector<double> >()) {
            stringstream solInfo, sol_csv, neh_csv;
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "            U P P E R    B O U N D     D O N E ";
            LOG(INFO) << "===========================================================================";
            LOG(INFO) << "[Robust Upper Bound] " << ub_name << " finished";
            stringstream ss_T;
            for (int x = 0; x < s_worst.T.size(); x++) {
                ss_T << s_worst.T[x] << " ";
            }
            // solInfo << "scenario: " << s_worst << "\n";
            solInfo << "permutation_neh: " << neh_permutation << "\n";
            solInfo << "cmax_neh: " << neh_value << "\n";
            solInfo << "time_neh: " << neh_time << "\n";
            solInfo << "permutation: " << pi_s_worst_star.getPermutation() << "\n";
            solInfo << "worst cmax: " << pi_s_worst_star.value << "\n";
            solInfo << "time: " << time_spent << "\n";
            solInfo << "time to best solution: " << time_to_best_sol << "\n";
            solInfo << "iterations: " << iterations << "\n";
            solInfo << "num_visited_solutions: " << num_visited_solutions << "\n";
            solInfo << "num_improvements: " << num_improvements << "\n";
            solInfo << "execution_id: " << executionId << "\n";
            solInfo << "seed: " << upper_bound_params.seed << "\n";
            solInfo << "adaptive_construction: " << upper_bound_params.adaptive_construction << "\n";
            solInfo << "const_beta1: " << upper_bound_params.beta1 << "\n";
            solInfo << "const_beta2: " << upper_bound_params.beta2 << "\n";
            solInfo << "vnd_size: " << upper_bound_params.vnd_size << "\n";
            solInfo << "first_improvement: " << upper_bound_params.first_improvement << "\n";
            solInfo << "random_vnd: " << upper_bound_params.random_vnd << "\n";
            solInfo << "budget-T: " << ss_T.str() << "\n";
            solInfo << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            std::cout << "=========================================================================\n";
            std::cout << "            U P P E R    B O U N D     D O N E \n";
            std::cout << "=========================================================================\n";
            std::cout << "[Robust Upper Bound] " << ub_name << " finished\n";
            // std::cout << "Scenario (s_worst): " << s_worst << "\n";
            std::cout << "NEH: c=" << neh_value << "; t=" << neh_time << "\n";
            std::cout << "Solution (pi_s_worst_star): " << pi_s_worst_star << "\n";
            std::cout << "Objective function value: " << pi_s_worst_star.value << "\n";
            std::cout << "Time Spent (s): " << time_spent << "\n";
            std::cout << "Time to best solution (s): " << time_to_best_sol << "\n";
            std::cout << "Iterations: " << iterations << "\n";
            std::cout << "Number of visited solutions: " << num_visited_solutions << "\n";
            std::cout << "Number of improvements: " << num_improvements << "\n";
            std::cout << "execution_id: " << executionId << "\n";
            std::cout << "seed: " << upper_bound_params.seed << "\n";
            std::cout << "budget-T: " << ss_T.str() << "\n";
            std::cout << "time-factor: " << upper_bound_params.grasp_tl << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] NEH: c=" << neh_value << "; t=" << neh_time;
            LOG(INFO) << "[Robust BB " << ub_name << "] Solution (pi_s_worst_star): " << pi_s_worst_star;
            // LOG(INFO) << "[Robust BB " << ub_name << "] Scenario (s_worst): " << s_worst;
            LOG(INFO) << "[Robust BB " << ub_name << "] Calculating worst solution value...";
            LOG(INFO) << "[Robust BB " << ub_name << "] Objective function value: " << pi_s_worst_star.value;
            LOG(INFO) << "[Robust BB " << ub_name << "] Time Spent (s): " << time_spent;
            LOG(INFO) << "[Robust BB " << ub_name << "] Time to best solution (s): " << time_to_best_sol
                                    << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] Iterations: " << iterations << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] Number of visited solutions: "
                                    << num_visited_solutions << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] Number of improvements: " << num_improvements
                                    << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] execution_id: " << executionId << "\n";
            LOG(INFO) << "[Robust BB " << ub_name << "] seed: " << upper_bound_params.seed << "\n";
            std::cout << "-----------------------------------------------------------------------------" << "\n";

            string exp_folder = util::FileUtil::create_output_folder(base_folder, "results");
            exp_folder = util::FileUtil::create_output_folder(exp_folder, executionId);
            if(! upper_bound_params.parametrize) {
                // append csv result to all experiments csv output file
                boost::filesystem::path all_results_file_path(base_folder);
                all_results_file_path /= ub_name + string("-all_results.csv");
                if (!boost::filesystem::exists(all_results_file_path.string())) {
                    // print CSV header to file
                    sol_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "alpha, n, m"
                            << ", ";
                    sol_csv << "budget_gamma" << ", ";
                    sol_csv << "solution_value, permutation, time_spent, time_to_best_sol, iterations, ";
                    sol_csv
                            << "num_visited_solutions, num_improvements, first_improvement, vnd_size, vnd_permutation, ";
                    sol_csv << "random_vnd, adaptive_construction, const_beta1, const_beta2, time_factor\n";
                }
                boost::filesystem::path neh_results_file_path(base_folder);
                neh_results_file_path /= string("NEH-all_results.csv");
                if (!boost::filesystem::exists(neh_results_file_path.string())) {
                    // print CSV header to file
                    neh_csv << "executionId, seed, ub_name" << ", " << "instance_name" << ", " << "alpha, n, m"
                            << ", budget_gamma" << ", solution_value, permutation, time_spent \n";
                }
                // print CSV contents to file
                if(neh_value > 0) {
                    std::ostringstream vts;
                    std::copy(neh_permutation.begin(), neh_permutation.end(),
                              std::ostream_iterator<unsigned>(vts, " "));
                    neh_csv << executionId << ", " << upper_bound_params.seed << ", " << "NEH" << ", " << instance->name
                            << ", " << instance->alpha << ", "
                            << instance->n << ", " << instance->m << ", "
                            << ss_T.str() << ", " << neh_value << ", " << vts.str() << ", " << neh_time << "\n";
                }
                sol_csv << executionId << ", " << upper_bound_params.seed << ", " << ub_name << ", " << instance->name
                        << ", " << instance->alpha << ", "
                        << instance->n << ", " << instance->m << ", ";
                sol_csv << ss_T.str() << ", " << pi_s_worst_star.value << ", " << pi_s_worst_star.getPermutation()
                        << ", ";
                sol_csv << time_spent << ", " << time_to_best_sol << ", ";
                sol_csv << iterations << ", " << num_visited_solutions << ", " << num_improvements << ", ";
                std::ostringstream vts2;
                std::copy(upper_bound_params.vnd_permutation.begin(), upper_bound_params.vnd_permutation.end(),
                          std::ostream_iterator<unsigned>(vts2, " "));
                sol_csv << upper_bound_params.first_improvement << ", " << upper_bound_params.vnd_size << ", ["
                        << vts2.str() << "], " << upper_bound_params.random_vnd;
                sol_csv << ", " << upper_bound_params.adaptive_construction << ", " << upper_bound_params.beta1
                        << ", " << upper_bound_params.beta2 << ", " << upper_bound_params.grasp_tl
                        << "\n";
                LOG(INFO) << "[Robust BB Upper Bound " << ub_name << "] Appending to NEH results file in "
                                        << neh_results_file_path.string();
                util::FileUtil::saveTextContentsToFile(neh_results_file_path.string(), neh_csv.str(), true);
                LOG(INFO) << "[Robust BB Upper Bound " << ub_name << "] Appending to all results file in "
                                        << all_results_file_path.string();
                util::FileUtil::saveTextContentsToFile(all_results_file_path.string(), sol_csv.str(), true);
                // save time results output file
                boost::filesystem::path time_result_file_path(exp_folder);
                time_result_file_path /= instance->name + string("-") + ub_name + string("-time_result.csv");
                LOG(INFO) << "[Robust BB " << ub_name << "] Saving time results file to "
                                        << time_result_file_path.string();
                util::FileUtil::saveTextContentsToFile(time_result_file_path.string(), timeResults, false);
                // save the other info results file
                boost::filesystem::path sol_info_file_path(exp_folder);
                sol_info_file_path /= instance->name + string("-") + ub_name + string("-solution_info.csv");
                LOG(INFO) << "[Robust BB " << ub_name << "] Saving solution info file to "
                                        << sol_info_file_path.string();
                util::FileUtil::saveTextContentsToFile(sol_info_file_path.string(), solInfo.str(), false);
                // Save constructive phase statistics
                if (alpha_hist.size() > 0) {
                    stringstream constInfo;
                    int iter = alpha_hist.size();
                    constInfo << "iter, alpha, probability, quality\n";
                    for (int i = 0; i < iter; i++) {
                        constInfo << i << ", " << alpha_hist[i] << ", ";
                        for (double b : alpha_prob[i]) {
                            constInfo << b << " ";
                        }
                        constInfo << ", ";
                        for (double b : alpha_qual[i]) {
                            constInfo << b << " ";
                        }
                        constInfo << "\n";
                    }
                    boost::filesystem::path sol_const_file_path(exp_folder);
                    sol_const_file_path /= instance->name + string("-") + ub_name + string("-constructive_info.csv");
                    LOG(INFO) << "[Robust BB " << ub_name << "] Saving constructive info file to "
                                            << sol_const_file_path.string();
                    util::FileUtil::saveTextContentsToFile(sol_const_file_path.string(), constInfo.str(), false);
                }
            }
        }

    };
}

#endif //PFSP_OUTPUTUTIL_H
