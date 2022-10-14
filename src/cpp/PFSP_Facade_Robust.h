//
// Created by mlevorato on 11/25/19.
//

#ifndef FLOWSHOP_SOLVER_ROBUST_PFSP_FACADE_H
#define FLOWSHOP_SOLVER_ROBUST_PFSP_FACADE_H

#include <boost/filesystem.hpp>
#include "problem/PFSP_Parameters.h"
#include "problem/robust/budget/metaheuristic/cmax/PFSP_Cmax_Budget_UpperBound_Alt.h"
#include "problem/robust/budget/branchandbound/wct/PFSP_Budget_WCT_BranchBound.h"
#include "problem/robust/budget/branchandbound/wct/PFSP_WCT_Budget_Hybrid_BranchBound.h"
#include "problem/robust/budget/metaheuristic/wct/GRASPSolver_RobPFSP_WCT.h"
#include "ExecutionInfo.h"
#include "PFSPFileController.h"
#include "problem/robust/RobPFSInstance_Cmax.h"
#include "problem/robust/budget/PFSP_Cmax_Budget_WorstCase.h"
#include "problem/robust/budget/PFSP_WCT_Budget_WorstCase.h"

namespace facade {
    using namespace boost;
    using namespace parameters;
    using namespace robust;
    using namespace std;
    namespace fs = boost::filesystem;

    class PFSP_Facade_Robust {
    public:
        PFSP_Facade_Robust() {

        }

        void solve_problem(const PFSP_Parameters &pfsp_params, BranchBound_Parameters bb_params,
                Robust_Parameters &rob_params, const UpperBound_Parameters &upper_bound_params) {
            // output folder setup
            string base_dir = FileUtil::create_output_folder(pfsp_params.outputFolder);

            // Reads the instance from the specified text file
            // model versions: between 1 and 499 => Cmax objective; between 500 and 1000 => WCT objective;
            if(pfsp_params.model_version < 500) {  // makespan (Cmax)
                string file_type = string("kouvelis");
                if(pfsp_params.filePath.string().find("ying") !=std::string::npos)
                    file_type = string("ying");
                if(pfsp_params.filePath.string().find("taillard") !=std::string::npos)
                    file_type = string("ying");
                if(pfsp_params.filePath.string().find("inputs.txt") !=std::string::npos)
                    file_type = string("newformat");
                robust::RobPFSInstance_Cmax instance = PFSPFileController::readRobustInstance_Cmax_FromFile(pfsp_params,
                        rob_params, file_type);
                // Execution additional info
                string instance_name(instance.name);
                instance_name = instance_name.substr(0, instance_name.find('.'));
                ExecutionInfo info(pfsp_params.executionId, instance_name, pfsp_params.outputFolder);
                // Budgeted Uncertainty problems : repeat the budget value for the remaining (m - 1) machines
                if(instance.rob_params.budget_gamma.size() > 0) {
                    double gamma = instance.rob_params.budget_gamma[0];
                    if(instance.rob_params.budget_type == "machine"
                            && instance.rob_params.budget_gamma.size() < instance.m) {
                        instance.rob_params.budget_gamma.clear();
                        for (int x = 1; x <= instance.m; x++) {
                            instance.rob_params.budget_gamma.push_back(gamma);
                        }
                    }
                    stringstream budget_ss;
                    for(double x : instance.rob_params.budget_gamma) {
                        budget_ss << x << " ";
                    }
                    LOG(INFO) << "[Budget] The budget parameter is: " << budget_ss.str();
                }
                switch (pfsp_params.model_version) {
                    default:
                    case 104:  // Robust PFSP MiniMax Budget Alternative Upper Bound only -- GRASP
                    {
                        LOG(INFO)
                            << "Model version is: Robust PFSP Cmax Budget Alt Upper Bound only - GRASP.\n";
                        cout << "Model version is: Robust PFSP Cmax Budget Alt Upper Bound only - GRASP.\n";
                        PFSP_Cmax_Budget_UpperBound_Alt budget_ub;
                        instance.robust_type = robust::RobustType::budget;
                        PFSSolution sol = budget_ub.calculate(instance, pfsp_params, instance.rob_params, upper_bound_params);
                        break;
                    }
                }
            } else {  // weighted completion time (WCT)
                string file_type = string("ying");
                if(pfsp_params.filePath.string().find("inputs.txt") !=std::string::npos)
                    file_type = string("newformat");
                robust::RobPFSInstance_WCT instance = PFSPFileController::readRobustInstance_WCT_FromFile(pfsp_params,
                        rob_params, file_type);
                // Execution additional info
                string instance_name(instance.name);
                instance_name = instance_name.substr(0, instance_name.find('.'));
                ExecutionInfo info(pfsp_params.executionId, instance_name, pfsp_params.outputFolder);
                // Budgeted Uncertainty problems : repeat the budget value for the remaining (m - 1) machines
                if(instance.rob_params.budget_gamma.size() > 0) {
                    double gamma = instance.rob_params.budget_gamma[0];
                    if(instance.rob_params.budget_type == "machine"
                       && instance.rob_params.budget_gamma.size() < instance.m) {
                        instance.rob_params.budget_gamma.clear();
                        for (int x = 1; x <= instance.m; x++) {
                            instance.rob_params.budget_gamma.push_back(gamma);
                        }
                    }
                    stringstream budget_ss;
                    for(double x : instance.rob_params.budget_gamma) {
                        budget_ss << x << " ";
                    }
                    LOG(INFO) << "[Budget] The budget parameter is: " << budget_ss.str();
                }
                switch (pfsp_params.model_version) {
                    default:
                    case 604:  // Robust PFSP WCT Budget GRASP
                    {
                        LOG(INFO) << "Model version is: Robust PFSP WCT Budget GRASP only.\n";
                        cout << "Model version is: Robust PFSP WCT Budget GRASP only.\n";
                        pfsp::GRASPSolver_RobPFSP_WCT grasp_wct;
                        LOG(INFO) << "[Robust PFSP Budget] Solving Robust PFSP WCT instance with GRASP...";
                        instance.robust_type = robust::RobustType::budget;
                        PFSSolution grasp_sol = grasp_wct.solve(instance, 0, pfsp_params, upper_bound_params);
                        robust::PFSP_WCT_Budget_WorstCase::release_cplex_model();
                        break;
                    }
                    case 606:  // Robust PFSP PFSP_Budget_Continuous_BB WCT
                    {
                        LOG(INFO) << "Model version is: Robust PFSP WCT Budget Continuous B & B.\n";
                        cout << "Model version is: Robust PFSP WCT Budget Continuous B & B.\n";
                        instance.robust_type = robust::RobustType::budget;
                        robust::PFSP_WCT_Budget_Continuous_BB rob;
                        PFSSolution sol = rob.solve(instance, base_dir, pfsp_params, bb_params, instance.rob_params,
                                                    upper_bound_params);
                        LOG(INFO) << "total_time_worstcase_dp: " << rob.getProblem().get_lower_bound().total_time_worstcase_dp;
                        cout << "total_time_worstcase_dp: " << rob.getProblem().get_lower_bound().total_time_worstcase_dp << endl;
                        LOG(INFO) << "total_time_worstcase_mip: " << rob.getProblem().get_lower_bound().total_time_worstcase_mip;
                        cout << "total_time_worstcase_mip: " << rob.getProblem().get_lower_bound().total_time_worstcase_mip << endl;
                        LOG(INFO) << "num_leaf_nodes: " << rob.getProblem().get_lower_bound().num_leaf_nodes;
                        cout << "num_leaf_nodes: " << rob.getProblem().get_lower_bound().num_leaf_nodes << endl;
                        LOG(INFO) << "num_leaf_nodes_dp: " << rob.getProblem().get_lower_bound().num_leaf_nodes_dp;
                        cout << "num_leaf_nodes_dp: " << rob.getProblem().get_lower_bound().num_leaf_nodes_dp << endl;
                        LOG(INFO) << "num_leaf_nodes_infeasible: " << rob.getProblem().get_lower_bound().num_leaf_nodes_infeasible;
                        cout << "num_leaf_nodes_infeasible: " << rob.getProblem().get_lower_bound().num_leaf_nodes_infeasible << endl;
                        robust::PFSP_WCT_Budget_WorstCase::release_cplex_model();
                        break;
                    }
                    case 607:  // PFSP WCT Hybrid
                    {
                        LOG(INFO) << "\nInvoking CPLEX Hybrid Branch and bound for PFSP-WCT.\n";
                        std::cout << "Invoking CPLEX Hybrid Branch and bound for PFSP-WCT.\n";
                        robust::hybrid::PFSP_WCT_Budget_Hybrid_BranchBound rob_pfsp_hybrid;
                        rob_pfsp_hybrid.solve(instance, base_dir, pfsp_params, bb_params, instance.rob_params, upper_bound_params);
                        break;
                    }
                }
            }
        }
    };
}

#endif //FLOWSHOP_SOLVER_ROBUST_PFSP_FACADE_H
