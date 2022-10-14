//
// Created by mlevorato on 11/25/19.
//

#ifndef FLOWSHOP_SOLVER_DETERMINISTIC_PFSP_FACADE_H
#define FLOWSHOP_SOLVER_DETERMINISTIC_PFSP_FACADE_H

#include <boost/filesystem.hpp>

#include "problem/PFSP_Parameters.h"
#include "ExecutionInfo.h"
#include "PFSPFileController.h"
#include "problem/deterministic/PFSInstance.h"
#include "problem/deterministic/metaheuristic/GRASPSolver_PFSP_Cmax.h"


namespace facade {

    using namespace boost;
    using namespace parameters;
    namespace fs = boost::filesystem;

    class PFSP_Facade_Deterministic {
    public:
        PFSP_Facade_Deterministic() {

        }

        void solve_problem(const PFSP_Parameters &pfsp_params,
                           BranchBound_Parameters bb_params,
                           const UpperBound_Parameters &ub_params) {
            // Reads the instance from the specified text file
            PFSInstance instance = PFSPFileController::readFromFile(pfsp_params);
            if (pfsp_params.reverse) {  // generate reverse instance
                instance = instance.generate_reverse_instance();
            }
            // Execution additional info
            string instance_name(instance.name);
            instance_name = instance_name.substr(0, instance_name.find('.'));
            ExecutionInfo info(pfsp_params.executionId, instance_name, pfsp_params.outputFolder);
            if (bb_params.upper_bound < 0.0 and bb_params.ub_path != NULL) {
                // retrieve upper bound from file
                bb_params.upper_bound = retrieve_ub_from_file(instance_name, bb_params.ub_path);
                if (bb_params.upper_bound <= 0.0) {
                    LOG(INFO) << "UB not found in solution value file or equals zero.";
                    cerr << "WARN: UB not found in solution value file or equals zero.\n";
                }
            }
            // output folder setup
            string base_dir = FileUtil::create_output_folder(pfsp_params.outputFolder);
            switch (pfsp_params.model_version) {
                default:
                case 4:  // Deterministic PFSP_Cmax Upper Bound only
                {
                    LOG(INFO) << "Model version is: PFSP Upper Bound only.";
                    if(ub_params.grasp_maxiter > 0)
                    {
                        LOG(INFO) << "[PFSP] Solving PFSP instance with GRASP...";
                        problem::pfsp::GRASPSolver_PFSP_Cmax grasp_solver;
                        problem::common::PFSSolution grasp_sol = grasp_solver.solve(instance, 0, pfsp_params,
                                                                                    ub_params);
                        LOG(INFO) << "[PFSP Cmax GRASP] Solution value  = "
                                                << grasp_sol.value;
                    }
                    break;
                }
            }
        }

        double retrieve_ub_from_file(const string &which_instance, const boost::filesystem::path* ub_path) {
            std::ifstream infile(ub_path->string().c_str(), std::ios::in | std::ios::binary);
            LOG(INFO) << "Reading input UB file: " << ub_path->string();
            string fileId = ub_path->filename().string();
            double upper_bound = 0.0;  // if ub is not found in file

            if (infile) {
                // captures the first line of the file containing CSV header
                string line;
                std::getline(infile, line);
                char_separator<char> sep2(",");
                vector<string> vec;

                while (std::getline(infile, line)) {
                    trim(line);
                    tokenizer<char_separator<char> > tokens2(line, sep2);
                    vec.assign(tokens2.begin(), tokens2.end());
                    string instance_name = vec.at(0);
                    if(instance_name == which_instance) {
                        try {
                            sscanf(vec.at(1).c_str(), "%lf", &upper_bound);
                            break;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(1) = " << vec.at(1);
                        }
                    }
                }
                // LOG(INFO) << "Matrix read: \n" << matrix;
                infile.close();
                LOG(INFO) << "Successfully read UB input file.";
                LOG(INFO) << "Upper bound = " << upper_bound;

                if(upper_bound > 0.0)
                    return upper_bound;
                else
                    return -1.0;
            } else {
                LOG(ERROR) << "Failed to read input UB file!";
                throw "Failed to read input UB file!";
            }
        }
    };
}

#endif //FLOWSHOP_SOLVER_DETERMINISTIC_PFSP_FACADE_H
