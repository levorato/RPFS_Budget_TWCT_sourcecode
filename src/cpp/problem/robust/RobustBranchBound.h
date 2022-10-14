//
// Created by mlevorato on 7/30/19.
//

#ifndef FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_H
#define FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_H

#include "../PFSP_Parameters.h"
#include "../PFSSolution.h"
#include "../../util/NumericUtil.h"
#include "../../util/include/FileUtil.h"
#include "../GlobalTypes.h"
#include "../deterministic/branchandbound/DFSSearchStrategy.h"
#include "../../util/include/TimeDateUtil.h"
#include "../../util/include/CollectionUtil.h"
#include "../../controller/include/PFSPFileController.h"
#include "PFSPScenario.h"
#include "../GlobalTypes.h"
#include "RobustBranchBound_Node.h"

#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o
#include <glog/logging.h>


namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;
    using namespace problem::common;
    using namespace problem::pfsp::bb;
    using namespace util;
    namespace fs = boost::filesystem;
    using boost::timer::cpu_timer;


    /**
     * Branch-and-bound algorithm for robust scheduling problems.
     *  Main elements:
     *       - Branching scheme
     *       - Lower bounds
     *       - Dominance conditions for the problem (if exist)
     * @tparam RobustProblem
     * @tparam Node
     */
    template<class Instance, class RobustProblem, class NodeData>
    class RobustBranchBound {
    public:
        typedef BaseNode<Instance, NodeData> Node;

        RobustBranchBound() : problem(), previous_time_stat(-1.0), last_node_statistics() {

        }

        PFSSolution solve(Instance &instance, const string &output_file,
                          const PFSP_Parameters &pfsp_params,
                          BranchBound_Parameters bb_params,
                          const Robust_Parameters &rob_params,
                          const UpperBound_Parameters &upper_bound_params) {
            stringstream logfile;
            LOG(INFO) << "[Robust BB] Processing instance file '" << instance.name << "'...";
            boost::timer::cpu_timer timer;
            boost::filesystem::path outfilePath(output_file);
            boost::system::error_code returnedError;
            PFSSolution best_sol;
            try {
                timer.start();  // Time measurement
                // ================= Solve ===================
                LOG(INFO) << "[Robust BB] Solving model with B & B...";
                best_sol = run(instance, output_file, pfsp_params, bb_params, rob_params, upper_bound_params);

                LOG(INFO) << "[Robust BB] Best solution value  = " << best_sol.value << endl;
                double time_spent = TimeDateUtil::calculate_time_spent(timer);
                LOG(INFO) << "[Robust BB] Global time spent: " << time_spent << " s";
                if(output_file.rfind("/") != string::npos) {  // output folder setup
                    boost::filesystem::path rootPath(output_file.substr(0, output_file.rfind("/")));
                    boost::system::error_code returnedError;
                    boost::filesystem::create_directories(rootPath, returnedError);
                    if (returnedError) {
                        LOG(FATAL) << "Cannot create output directory: " << rootPath.string();
                        return best_sol;
                    }
                }
                // Save result (solution) file
                outfilePath /= instance.name + string("-bb-solution.txt");
                PFSPFileController::writeSolutionToFile(outfilePath, best_sol, time_spent);
            } catch (...) {
                cerr << "Unknown Exception" << endl;
                LOG(FATAL) << "[Robust BB] Unknown Exception.";
            }
            return best_sol;
        }

        PFSSolution run(Instance &instance, const string &output_file,
                        const PFSP_Parameters &pfsp_params,
                        BranchBound_Parameters bb_params,
                        const Robust_Parameters &rob_params,
                        const UpperBound_Parameters &upper_bound_params) {
            double time_limit = pfsp_params.timeLimit;
            cpu_timer timer;  // Start time measure
            timer.start();  // Time measurement
            std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
            stringstream timeResults;  // Stringstream containing the best result found at each moment of time.
            LOG(INFO) << "[Robust BB]  Starting B & B (time_limit = "<< time_limit << "s)...\n\n";
            unsigned long node_uid_seq = 1, phatomed_count = 0, iteration = 0;
            unsigned int num_improvements = 0;
            ublas::vector<long> phatomed_by(8, 0);
            bool exact_solution = false;
            Time ub = std::numeric_limits<Time>::max();
            std::vector<int> J = util::CollectionUtil::create_initial_job_list(instance.n);
            // Initialize the required data structures of the problem
            problem.initialize(instance, rob_params);
            // Computation of the lower and upper bounds of the root node N_0
            NodeData node_data(instance);
            Time root_lb = 0;
            if(bb_params.lower_bound > 0) {  // lower bound given as input
                root_lb = bb_params.lower_bound;
                LOG(INFO) << "[Robust BB] Fixing root lower bound = " << root_lb;
            } else {
                root_lb = problem.updatelb(instance, ub, rob_params, J, node_data);
            }
            Node root(node_uid_seq++, J, instance.m, instance.n, root_lb, &instance, node_data);
            PFSSolution best_sol;
            if(bb_params.upper_bound > 0) {  // upper bound given as input
                ub = bb_params.upper_bound;
                LOG(INFO) << "[Robust BB] Fixing upper bound = " << root_lb;
                if(bb_params.permutation.size() > 0) {
                    bool read_error = false;
                    stringstream permutation_ss;
                    for (int x : bb_params.permutation) {
                        if (x <= 0 || x > instance.n) {
                            LOG(ERROR)
                                << "ERROR: Values in initial permutation are NOT in the required range of [1,n].";
                            cout << "ERROR: Values in initial permutation are NOT in the required range of [1,n].\n";
                            read_error = true;
                            break;
                        }
                        permutation_ss << x << " ";
                    }
                    if(! read_error) {
                        LOG(INFO) << "[Robust BB] Fixing initial permutation = " << permutation_ss.str();
                        best_sol.permutation = bb_params.permutation;
                    } else {
                        LOG(INFO) << "[Robust BB] Skipping initial permutation, due to read error.";
                    }
                }
            } else {
                best_sol = problem.obtainInitialSolution(instance, root.LB, pfsp_params, rob_params, upper_bound_params);
                ub = best_sol.value;
            }
            PFSSolution ub_sol(best_sol);
            notifyNewValue(ub, 0.0, iteration, timeResults);  // register new solution value
            Node best_node = root;
            Node N_p = root;
            print_root_node_stats(ub, root.LB);
            DFSSearchStrategy<Node> dfs(instance.n, start_time, false);
            unsigned max_level = 0;
            unsigned num_solutions = 1, num_solutions_so_far = 0;

            if (NumericUtil::TestisGE(root.LB, ub)) {  // root.lb >= ub : optimal schedule found on the root node
                PFSSolution solution = best_sol;
                //ub = PFSProblem::calculateCmax(solution.permutation, instance);
                solution.value = ub;
                if(NumericUtil::TestisGT(root.LB, solution.value)) {
                    LOG(ERROR) << "[Robust BB] Root LB is greater than ub ! LB = " << root.LB;
                    std:cerr << "[Robust BB] ERROR : Root LB is greater than ub ! LB = " << root.LB << "\n";
                }
                LOG(INFO) << "[Robust BB] Found solution in the B&B root node. Makespan: "
                                        << ub << ", Obtained solution: " << solution;
                phatomed_count++;  phatomed_by[1]++;
                best_sol = solution;
                exact_solution = true;
            } else {  // Step 2: Branching of node N_p
                bool bb_done = false;
                while (! bb_done) {  // *inside priority queue: node selection strategy
                    double time_spent = TimeDateUtil::get_time_spent_in_seconds(timer.elapsed());
                    if(time_spent >= time_limit) {
                        LOG(INFO)<< "[Robust BB] Time limit! t=" << time_spent << "s";
                        std::cout<< "******* Time limit! t=" << time_spent << "s";
                        break;
                    }
                    iteration++;
                    double gap = ub - N_p.LB;
                    if(iteration > 1)  print_statistics(instance, N_p, last_node_statistics, phatomed_by, max_level, ub, bb_params, dfs);
                    // (b) Branch this node completely and decide lower bounds of all nodes created
                    last_node_statistics = NodeStatistics(dfs.size(), phatomed_count,
                                                          N_p.get_level(), ub, iteration, time_spent, root.LB, gap,
                                                          max_level, num_solutions_so_far);
                    std::vector<Node> child_nodes;  // creation of Node N descendant of N_p
                    for (unsigned count = 0; count < N_p.unscheduled_job_list.size(); count++) {
                        int job = N_p.unscheduled_job_list[count];
                        Node N(node_uid_seq++, N_p, job, instance.m, instance.n);
                        N.seq_1.push_back(job);
                        // N.update_C_k_seq_1(instance, job, N.seq_1); // Update C_k,seq_1 of job
                        N.LB = max(N_p.LB, problem.updatelb(instance, ub, rob_params, N.seq_1, N.unscheduled_job_list, N.node_data));
                        if ( NumericUtil::TestisGE(N.LB, ub) ) {
                            phatomed_by[1]++;  phatomed_count++;
                        } else if(problem.dominance_rules(instance, ub, rob_params, N.seq_1, N.unscheduled_job_list, N.node_data)) {
                            phatomed_by[2]++;
                            phatomed_count++;
                        /* } else if(NumericUtil::TestisLE(fabs(ub - N.LB) / ub, BB_CONVERGENCE_EPS)) {
                            std::cerr << "BB WARN: Node phatomed by B&B convergence epsilon! \n";
                            phatomed_by[3]++;
                            phatomed_count++; // PRUNING BY OPTIMALITY DISABLED! */
                        } else {
                            child_nodes.push_back(N);
                        }
                    }
                    // default : normal next node selection strategy
                    dfs.add_nodes(child_nodes, N_p.get_level() + 1);
                    if (bb_done) break;
                    // Step 3:   N O D E     S E L E C T I O N    P O L I C Y
                    double previous_ub = ub;
                    std::list<Node> pruned_nodes;
                    // default : normal next node selection strategy
                    N_p = dfs.find_next_node(ub, phatomed_count, iteration, best_node, pruned_nodes, num_solutions, num_improvements);
                    if(ub < previous_ub) {
                        PFSSolution solution(best_node.seq_1, instance.n);
                        solution.value = ub;
                        best_sol = solution;  num_solutions_so_far++;
                        notifyNewValue(ub, time_spent, iteration, timeResults);
                    }
                    long queue_size =dfs.size();
                    if(N_p.get_id() < 0) {  //or queue_size == 0) {  //  ended processing all nodes : queue is empty!
                        std::cout << "Branch-and-bound done. \n";
                        exact_solution = true;  bb_done = true;  break;
                    }
                }
            }
            double time_spent = TimeDateUtil::calculate_time_spent(timer);
            if(! exact_solution) {
                // default : normal next node selection strategy
                best_node = dfs.getTopNode();
                // std::cout << " [" << best_node.LB << "]" << endl;
            }
            print_solution_summary(instance, exact_solution, best_sol, best_node, iteration, time_spent,
                                   phatomed_count, phatomed_by, timeResults, pfsp_params.outputFolder,
                                   max_level, pfsp_params.executionId, num_improvements, time_limit,
                                   rob_params.dominance_rules);
            best_sol.is_optimal = exact_solution;
            if (best_sol.value > ub_sol.value)  LOG(ERROR)<< "Exact solution worse than heuristic!!!";
            problem.finalize();
            return best_sol;
        }

        inline void notifyNewValue(const double& solution_value, const double& timeSpent,
                                              const int& iteration, stringstream& timeResults) {
            timeResults << fixed << setprecision(4) << timeSpent << "," << solution_value
                        << "," << (iteration+1) << "\n";
        }

        inline void print_root_node_stats(const double& UB, const double& LB) {
            double gap = (LB != 0) ? (100.0 * (UB - LB) / LB) : std::numeric_limits<double>::infinity();
            LOG(INFO) << "[Robust BB]  ***  Root node : UB = " << std::fixed << std::setprecision(2)
                                        << UB << " ; LB = " << LB
                                        << " ; gap% = " << gap;
            std::cout << "  *** Root node : UB = " << std::fixed << std::setprecision(2) << UB << " ; LB = " << LB
                                        << " ; gap% = " << gap << "\n";
        }

        inline void print_statistics(Instance &instance, Node &node, const NodeStatistics &stat,
                                     ublas::vector<long>& phatomed_by, unsigned &max_level, const double& ub,
                                     BranchBound_Parameters &bb_params, 
                                     DFSSearchStrategy<Node> &dfs) {
            if (node.get_level() > max_level) {
                max_level = node.get_level();
                LOG(INFO) << " BB Reached Level " << max_level << " / Objective* = " << ub;
            }
            if (previous_time_stat == -1.0 || stat.time_spent / 10 > previous_time_stat) {
                // Gap calculation : Obtain the lowest lower bound so far
                Time llb = 0;
                long queue_size = dfs.size();
                // default : normal next node selection strategy
                llb = dfs.getTopNode().LB;
                Time gap = stat.ub - llb;
                previous_time_stat = stat.time_spent / 10;
                LOG(INFO) << " Iteration " << stat.iteration;
                LOG(INFO) << " BB Queue size = " << queue_size << " / Phatomed so far: "
                                        << stat.phatomed_count;
                LOG(INFO) << "Best solution so far = " << stat.ub;
                LOG(INFO) << "Max Level so far = " << stat.max_level;
                LOG(INFO) << "Solution gap = " << gap;
                LOG(INFO) << "% gap = " << std::fixed << std::setprecision(2) << (100 * gap / stat.ub);
                LOG(INFO) << "[Robust BB] Time spent: " << stat.time_spent << " s";
            }
        }

        inline void print_solution_summary(Instance &instance, bool exact_solution, PFSSolution &best_sol, Node &best_node,
                                           int iteration, double time_spent, long phatomed_count, ublas::vector<long> &phatomed_by,
                                           stringstream &timeResults, const string &base_folder, unsigned max_level,
                                           string executionId, int num_improvements, double time_limit,
                                           bool dominance_rules) {
            stringstream solInfo, sol_csv;
            LOG(INFO) << "=========================================================================";
            LOG(INFO) << "            B R A N C H   A N D   B O U N D   D O N E ";
            LOG(INFO) << "=========================================================================";
            LOG(INFO) << "[Robust BB] Solution type: "
                                    << (exact_solution ? "Exact" : "Approximate");
            solInfo << "solution_type: " << (exact_solution ? "Exact" : "Approximate") << "\n";
            if(! exact_solution) {
                LOG(INFO) << "[Robust BB] Initial (heuristic) solution found to be optimal!";
                std::cerr << "Initial (heuristic) solution found to be optimal!\n";
                double gap = best_sol.value - best_node.LB;
                LOG(INFO) << "[Robust BB] Solution gap = " << gap;
                LOG(INFO) << "[Robust BB] % gap = " << std::fixed << std::setprecision(2) << (100 * gap / best_sol.value);
                solInfo << "solution_gap: " << gap << "\n";
                solInfo << "gap_perc: " << std::fixed << std::setprecision(2) << (100 * gap / best_node.LB) << "\n";
            }
            LOG(INFO) << "[Robust BB] Number of iterations: " << iteration;
            solInfo << "iterations: " << iteration << "\n";
            LOG(INFO) << "[Robust BB] Solution: " << best_sol;
            solInfo << "permutation: " << best_sol.getPermutationAsString() << "\n";
            LOG(INFO) << "[Robust BB] Time Spent (s): " << time_spent;
            solInfo << "time: " << time_spent << "\n";
            LOG(INFO) << "[Robust BB] Time Limit (s): " << time_limit;
            solInfo << "time_limit: " << time_limit << "\n";
            std::cout << "=========================================================================\n";
            std::cout << "            B R A N C H   A N D   B O U N D   D O N E \n";
            std::cout << "=========================================================================\n";
            std::cout << "Time Spent (s): " << time_spent << "\n";
            std::cout << "Solution: " << best_sol << "\n";
            LOG(INFO) << "[Robust BB] Calculating BB solution value...";
            LOG(INFO) << "[Robust BB] Worst case solution value: " << best_sol.value <<
                                    ", Obtained solution: " << best_sol;
            solInfo << "worst_case: " << best_sol.value << "\n";
            LOG(INFO) << "[Robust BB] Original objective: " << instance.objective_name();
            solInfo << "original_objective: " << instance.objective_name() << "\n";
            for(int i = 1; i <= 4; i++) {
                LOG(INFO) << "[BBSequentialPFSP] Num nodes phatomed by LB_" << i << ": " << phatomed_by(i);
                std::cout << "+ Num nodes phatomed by LB_" << i << ": " << phatomed_by(i) << "\n";
            }
            LOG(INFO) << "[Robust BB] Total number of phatomed nodes: " << phatomed_count;
            std::cout << "-----------------------------------------------------" << "\n";
            std::cout << "= Total number of phatomed nodes: " << phatomed_count << "\n";
            solInfo << "phatomed: " << phatomed_count << "\n";
            std::cout << "* Number of processed nodes (NN): " << iteration << "\n";
            solInfo << "NN: " << iteration << "\n";
            LOG(INFO) << "[Robust BB] * Number of processed nodes (NN): " << iteration;
            std::cout << "* Number of improvements: " << num_improvements << "\n";
            solInfo << "num_improvements: " << num_improvements << "\n";
            LOG(INFO) << "[Robust BB] * Number of improvements: " << num_improvements;

            if(! exact_solution) {
                double gap = (best_node.LB != 0) ? ((best_sol.value - best_node.LB) / best_node.LB) : std::numeric_limits<double>::infinity();
                std::cout << "gap: " << gap << "\n";
                std::cout << "is_optimal: false\n";
            } else {
                std::cout << "gap: 0\n";
                std::cout << "is_optimal: true\n";
            }
            // Objective function (worst-case) validation
            PFSPScenario wscenario = problem.calculate_objective(instance, best_sol.permutation);
            if(fabs(wscenario.value - best_sol.value) <= 1e-5) {
                std::cout << "validated: true\n";
            } else {
                std::cout << "validated: false\n";
                LOG(INFO) << "[Robust BB] Obj validation FAILED: best_sol.value = " << best_sol.value
                    << " x validation obj = " << wscenario.value;
            }

            // append csv result to all experiments csv output file
            boost::filesystem::path all_results_file_path(base_folder);
            all_results_file_path /= string("BB_") + problem.get_problem_name() + string("-all_results.csv");
            if (!boost::filesystem::exists(all_results_file_path.string())) {
                // print CSV header to file
                sol_csv << "executionId; problem_name; instance_name; n; m; ";
                sol_csv << "solution_value; permutation; time_spent; time_limit; NN; phatomed; ";
                sol_csv << "num_improvements; solution_type; dominance_rules\n";
            }
            sol_csv << executionId << "; " << problem.get_problem_name() << "; " << instance.name << "; "
                    << instance.n << "; " << instance.m << "; ";
            sol_csv << best_sol.value << "; " << best_sol.getPermutationAsString() << "; ";
            sol_csv << time_spent << "; " << time_limit << "; " << iteration << "; ";
            sol_csv << phatomed_count << "; " << num_improvements << "; ";
            sol_csv << (exact_solution ? "Exact" : "Approximate") << "; " << dominance_rules << "\n";
            LOG(INFO) << "[Robust BB] Appending to all results file in "
                                    << all_results_file_path.string();
            util::FileUtil::saveTextContentsToFile(all_results_file_path.string(), sol_csv.str(), true);
            
            // save time results output file
            boost::filesystem::path time_result_file_path(base_folder);
            time_result_file_path /= instance.name + string("-bb-time_result.csv");
            LOG(INFO) << "[Robust BB] Saving time results file to " << time_result_file_path.string();
            util::FileUtil::saveTextContentsToFile(time_result_file_path.string(), timeResults.str(), false);
            // save the other info results file
            boost::filesystem::path sol_info_file_path(base_folder);
            sol_info_file_path /= instance.name + string("-bb-solution_info.csv");
            LOG(INFO) << "[Robust BB] Saving solution info file to " << sol_info_file_path.string();
            util::FileUtil::saveTextContentsToFile(sol_info_file_path.string(), solInfo.str(), false);
        }

        RobustProblem getProblem() {
            return problem;
        }
        
    private:
        RobustProblem problem;
        double previous_time_stat;
        NodeStatistics last_node_statistics;
    };
}

#endif //FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_H
