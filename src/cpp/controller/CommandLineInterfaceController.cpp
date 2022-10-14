/*
 * CommandLineInterfaceController.cpp
 *
 *  Created on: Apr 30, 2013
 *      Author: mario
 */

#include "../problem/PFSP_Parameters.h"
#include "include/CommandLineInterfaceController.h"
#include "../util/include/TimeDateUtil.h"
#include "../util/include/EnumUtil.h"
#include "../util/include/ExecutionInfo.h"
#include "../util/include/RandomUtil.h"
#include "../PFSP_Facade_Deterministic.h"
#include "../PFSP_Facade_Robust.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/timer/timer.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/core.hpp>
#include <glog/logging.h>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/support/date_time.hpp>

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <execinfo.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <vector>
#include <cmath>

namespace controller {

    using namespace boost;
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;
    using namespace std;
    using namespace util;
    using namespace facade;

    // A helper function to simplify the main part.
    template<class T>
    ostream &operator<<(ostream &os, const std::vector<T> &v) {
        copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
        return os;
    }

    CommandLineInterfaceController::CommandLineInterfaceController() : logSeverity("info") {
        // TODO Auto-generated constructor stub

    }

    CommandLineInterfaceController::~CommandLineInterfaceController() {
        // TODO Auto-generated destructor stub
    }

    void CommandLineInterfaceController::PrintVariableMap(const boost::program_options::variables_map vm) {
        stringstream ss;
        for (po::variables_map::const_iterator it = vm.begin(); it != vm.end(); it++) {
            ss << "> " << it->first;
            if (((boost::any)it->second.value()).empty()) {
                ss << "(empty)";
            }
            if (vm[it->first].defaulted() || it->second.defaulted()) {
                ss << "(default)";
            }
            ss << "=";

            bool is_char = false;
            /*
            try {
                boost::any_cast<const char *>(it->second.value());
                is_char = true;
            } catch (const boost::bad_any_cast &) {
                is_char = false;
            } */
            bool is_str = false;
            /*
            try {
                boost::any_cast<std::string>(it->second.value());
                is_str = true;
            } catch (const boost::bad_any_cast &) {
                is_str = false;
            } */
            auto& value = it->second.value();
            if (auto v = boost::any_cast<uint32_t>(&value)) {
                std::cout << *v;
            } else if (auto v = boost::any_cast<char>(&value)) {
                std::cout << *v;
            } else if (auto v = boost::any_cast<std::string>(&value)) {
                //std::cout << *v;
                std::string temp = vm[it->first].as<std::string>();
                if (temp.size()) {
                    ss << temp << std::endl;
                } else {
                    ss << "true" << std::endl;
                }
            } else if (auto v = boost::any_cast<int>(&value)) {
                ss << vm[it->first].as<int>() << std::endl;
            } else if (auto v = boost::any_cast<unsigned>(&value)) {
                ss << vm[it->first].as<unsigned>() << std::endl;
            } else if (auto v = boost::any_cast<long>(&value)) {
                ss << vm[it->first].as<long>() << std::endl;
            } else if (auto v = boost::any_cast<bool>(&value)) {
                ss << vm[it->first].as<bool>() << std::endl;
            } else if (auto v = boost::any_cast<double>(&value)) {
                ss << vm[it->first].as<double>() << std::endl;
            } else if (auto v = boost::any_cast<const char *>(&value)) {
                ss << vm[it->first].as<const char * >() << std::endl;
            } else { // Assumes that the only remainder is vector<string>
                //ss << vm[it->first].as<const char * >() << std::endl;
                ss << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
            }
        }
        LOG(INFO) << "Parsed command-line args: \n" << ss.str();
    }

    void CommandLineInterfaceController::processInputFile(const PFSP_Parameters &pfsp_params,
                                                          const BranchBound_Parameters &bb_params,
                                                          Robust_Parameters &rob_params,
                                                          const UpperBound_Parameters &ub_params,
                                                          const bool &robust) {
        if (fs::exists(pfsp_params.filePath) && fs::is_regular_file(pfsp_params.filePath)) {
            if(not robust) {
                LOG(INFO) << "Facade: Processing deterministic problem.";
                cerr << "Facade: Processing deterministic problem.\n";
                //LOG(INFO) << "Model cuts are " << (enable_cuts ? "enabled." : "disabled.");
                // upper bound fixation
                if (bb_params.upper_bound < 0.0 and bb_params.ub_path != NULL) {
                    LOG(INFO) << "UB fixed by solution value file.";
                    cerr << "WARN: UB fixed by value file.\n";
                } else if (bb_params.upper_bound > 0.0) {
                    LOG(INFO) << "UB fixed by command-line argument.";
                    cerr << "WARN: UB fixed by command-line argument.\n";
                }
                // root lower bound fixation
                if(bb_params.lower_bound > 0.0) {
                    LOG(INFO) << "Root LB fixed by command-line argument.";
                    cerr << "WARN: Root LB fixed by command-line argument.\n";
                }
                PFSP_Facade_Deterministic det_facade;
                det_facade.solve_problem(pfsp_params, bb_params, ub_params);
            } else {  // Robust PFSP
                LOG(INFO) << "Facade: Processing robust problem.";
                cerr << "Facade: Processing robust problem.\n";
                if (bb_params.upper_bound > 0.0) {
                    LOG(INFO) << "WARN: UB fixed by command-line argument.";
                    cerr << "WARN: UB fixed by command-line argument.\n";
                }
                PFSP_Facade_Robust rob_facade;
                rob_facade.solve_problem(pfsp_params, bb_params, rob_params, ub_params);
            }
        } else {
            LOG(FATAL) << "Invalid file: '" << pfsp_params.filePath.string() << "'. Try again." << endl;
        }
    }

    bool CommandLineInterfaceController::testInputFile(fs::path filePath) {
        if (fs::exists(filePath) && fs::is_regular_file(filePath)) {
            return true;
        } else {
            LOG(FATAL) << "Invalid file: '" << filePath.string() << "'. Try again." << endl;
            return false;
        }
    }

// http://www.concentric.net/~Ttwang/tech/inthash.htm
    unsigned long mix(unsigned long a, unsigned long b, unsigned long c) {
        a = a - b;
        a = a - c;
        a = a ^ (c >> 13);
        b = b - c;
        b = b - a;
        b = b ^ (a << 8);
        c = c - a;
        c = c - b;
        c = c ^ (b >> 13);
        a = a - b;
        a = a - c;
        a = a ^ (c >> 12);
        b = b - c;
        b = b - a;
        b = b ^ (a << 16);
        c = c - a;
        c = c - b;
        c = c ^ (b >> 5);
        a = a - b;
        a = a - c;
        a = a ^ (c >> 3);
        b = b - c;
        b = b - a;
        b = b ^ (a << 10);
        c = c - a;
        c = c - b;
        c = c ^ (b >> 15);
        return c;
    }

    int CommandLineInterfaceController::processArgumentsAndExecute(int argc, char *argv[]) {

        // used for debugging purpose
        std::set_terminate(handler);

        cout << "\nPFSP solver" << endl;

        // reads the system properties from ini file
        this->readPropertiesFile();

        unsigned long the_seed = 0;
        bool debug = false;
        string inputFileDir;
        string outputFolder;
        double timeLimit = 3600.0;
        string jobid;
        int model_version = 1;
        bool enable_cuts = true;
        string solfname, outputsolfname;
        bool solrequired = false;
        double upper_bound = -1.0;
        double lower_bound = -1.0;
        int lb_strategy = 0;
        bool reverse = false;
        bool async_ub = false;
        bool robust = false;
        unsigned mip_depth = 0;
        bool mip_use_combinatorial_bounds = false;
        // ====== Robust PFSP Parameters
        Robust_UB_Type rob_ub_type = Robust_UB_Type::grasp;
        unsigned ub_iter_max = 0;  // maximum number of iterations (used by SA, IG robust upper bounds)
        // Budget parameters : T (one value per machine)
        std::string budget_gamma_str;
        // Budget parameters : budget_gamma[1 .. n]
        std::vector<double> budget_gamma;
        // Budget type (machine or global)
        string budget_type = "global";
        // Upper bound parameters
        unsigned grasp_tl = 30, grasp_maxiter = 100000;
        double const_beta1 = 0.6, const_beta2 = 1.0;
        // VND local search parameters
        std::vector<unsigned> vnd_permutation;
        unsigned vnd_size = 4;
        bool first_improvement = true;
        bool random_vnd = true;
        bool adaptive_construction = true;
        bool parametrize = false;
        bool dominance_rules = true;
        unsigned obj_update_freq = 0;
        std::vector<int> initial_permutation;
        // number of CPU cores to be used in parallel algorithms
        unsigned max_cores = 4;
        // Rob-PFSP Hybrid MP type
        string mp_model_name_str = string("hybrid-wilson");

        po::options_description desc("Available options:");
        desc.add_options()
                ("help", "show program options")
                ("time-limit,t", po::value<double>(&timeLimit)->default_value(3600.0), "maximum execution time (seconds)")
                ("input-file,f", po::value<std::vector<string> >(), "instance input file")
                ("debug", po::value<bool>(&debug)->default_value(false), "enable debug mode")
                ("input-file-dir", po::value<string>(&inputFileDir),
                 "input instance file directory (processes all files inside)")
                ("output-folder", po::value<string>(&outputFolder), "output folder for result files")
                ("jobid", po::value<string>(&jobid), "job identifier (string)")
                // model versions: between 1 and 499 => Cmax objective; between 500 and 1000 => WCT objective;
                ("model-version", po::value<int>(&model_version)->default_value(1),
                 "model version number (104 for Robust PFSP Cmax Budget GRASP; "
                 "604 for Robust PFSP WCT Budget GRASP; "
                 "606 for Robust PFSP PFSP_Budget_Continuous_BB WCT; 607 for PFSP WCT Hybrid MP Branch-and-Bound.")
                //("enable-cuts", po::value<bool>(&enable_cuts)->default_value(true), "enable TSP cuts in the model")
                ("solfname,o", po::value<string>(&solfname), "input precomputed solution for training (oracle)")
                ("sol", po::value<string>(&outputsolfname), "output file to write the solution (output solution filename)")
                ("ub", po::value<double>(&upper_bound)->default_value(-1.0), "fix Upper Bound value for instance")
                ("ub-file", po::value<std::vector<string> >(), "upper bound values input file")
                ("lb", po::value<double>(&lower_bound)->default_value(-1.0), "fix root Lower Bound value for instance")
                ("lb-strategy", po::value<int>(&lb_strategy)->default_value(0), "Lower Bound strategy (which LB is invoked according to node depth): "
                                                                                "0 for Ladhari default, 1 for Ladhari A1, 2 for Ritt.")
                ("sel-strategy", po::value<string>()->default_value("dfs"), "Node selection strategy: 'dfs' (Ladhari) or 'cplex' (CPLEX default).")
                ("branchingtype", po::value<string>()->default_value("alt"), "choose branching type from alt(ernating), sm(allest)d(ynamic)")
                ("rev", po::value<bool>(&reverse)->default_value(false), "run BB algorithm over the reserse instance (matrix P with reverse order of machines).")
                ("async-ub", po::value<bool>(&async_ub)->default_value(false), "run asynchronous Upper bound procedure (default = false).")
                ("mip-usage", po::value<string>()->default_value("no-mip"), "MIP Usage: 'no-mip', 'linear-relax' or 'full-mip'.")
                ("mip-depth", po::value<unsigned>(&mip_depth)->default_value(0), "Initial B & B tree depth to use auxiliary MIP model. Default = 0.")
                ("mip-use-c-bounds", po::value<bool>(&mip_use_combinatorial_bounds)->default_value(false), "true if combinatorial bounds (LB_1, etc) are to be "
                                                                                                           "used inside the PFSP MIP Model, false otherwise (default).")
                ("rob-ub-type", po::value<string>()->default_value("grasp"), "Robust PFSP Upper bound algorithm: grasp (Levorato).")
                ("ub-iter-max", po::value<unsigned>(&ub_iter_max)->default_value(50), "Robust upper bound max number of iterations. Default = 50 * n.")
                ("const-beta1", po::value<double>(&const_beta1)->default_value(0.6), "Metaheuristic constructive parameter beta1. Default = 0.6.")
                ("const-beta2", po::value<double>(&const_beta2)->default_value(1.0), "Metaheuristic constructive parameter beta2. Default = 1.0.")
                ("grasp-tf", po::value<unsigned>(&grasp_tl)->default_value(30), "GRASP upper bound Time Factor. Default = 30.")
                ("grasp-maxiter", po::value<unsigned>(&grasp_maxiter)->default_value(100000), "GRASP upper bound max number of iterations. Default = 100000.")
                ("budget-gamma", po::value<string>()->default_value(""), "Robust budget parameters T (in %), separated by spaces, eg: T[1] T[2] .. T[n]}.")
                ("budget-type", po::value<string>()->default_value("global"), "Robust budgeted uncertainty set type: 'global' (default) or 'machine'.")
                ("vnd-size", po::value<unsigned>(&vnd_size)->default_value(4), "Metaheuristic VND Local Search size. Default = 8.")
                ("first-improvement", po::value<bool>(&first_improvement)->default_value(true), "Metaheuristic VND First improvement ? Default = true.")
                ("vnd-permutation", po::value<string>()->default_value(""), "Metaheuristic VND permutation, separated by spaces, eg: 3 1 4 2 (for VND size 4).")
                ("random-vnd", po::value<bool>(&random_vnd)->default_value(true), "Enable randomized VND (RVND).")
                ("adaptive-construction", po::value<bool>(&adaptive_construction)->default_value(true), "Enable adaptive construction phase in metaheuristics.")
                ("seed", po::value<unsigned long>(&the_seed)->default_value(0), "Random seed, for usage in stochastic algorithms.")
                ("parametrize", po::value<bool>(&parametrize)->default_value(false), "Enable upper bound parametrization mode (disables solution output to file).")
                ("dominance", po::value<bool>(&dominance_rules)->default_value(true), "Enable dominance rules (if they exist) in combinatorial branch-and-bound.")
                ("initial-permutation", po::value<string>()->default_value(""), "Initial permutation (used by Branch-and-Bound procedure) - optional.")
                ("cuts-file", po::value<std::vector<string> >(), "Model cuts input file.")
                ("mp-model-name", po::value<string>()->default_value("hybrid-wilson"), "Robust Budget PFSP Hybrid Master Problem (MP) type: 'hybrid-wilson' (default) or 'hybrid-liao-you'.")
                ("max-cores", po::value<unsigned>(&max_cores)->default_value(4), "Maximum number of CPU cores to be used in parallel algorithms (e.g., MILP models).");
        po::positional_options_description p;
        p.add("input-file", -1);
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                options(desc).positional(p).run(), vm);
        po::notify(vm);

        // job id is obtained through command line parameter from PBS Scheduler
        //cout << "Job id is " << jobid << "\n";
        // initializes the logging subsystem
        if (jobid.length() == 0) {
            jobid = TimeDateUtil::generateRandomId();
        } /* else {
            jobid += "-" + TimeDateUtil::getDateAndTime();
        } */
        std::cout << "execution_id: " << jobid << "\n";
        try {
            this->initLogging(outputFolder, jobid, argv);
        }
        catch (std::exception &e) {
            cerr << "Fatal application error in log init.\n";
            cerr << "Abnormal program termination. Stracktrace: " << endl;
            cerr << e.what() << "\n";
            if (std::string const *stack = boost::get_error_info<stack_info>(e)) {
                cerr << stack << endl;
            }
            cerr << diagnostic_information(e);
            return 1;
        }
        // print passed arguments to log file
        PrintVariableMap(vm);

        // id used for output folders
        string executionId = jobid;
        try {
            if (vm.count("help")) {
                cout << "Usage: pfsp [options]\n";
                cout << desc;
                return 0;
            }

            if (not vm.count("output-folder")) {
                cout << "Please specify the output folder." << endl;
                return 1;
            }
            LOG(INFO) << "The maximum number of cores is " << max_cores;
            fs::path ub_file;
            bool ub_file_ok = false;
            if (vm.count("ub-file")) {
                LOG(INFO) << "Upper bound input file is: "
                                        << vm["ub-file"].as<std::vector<string> >().at(0) << "\n"
                                        << "Warning: best solution may be wrong, since explicit upper bound has been given.";
                std::cout << "Upper bound input file is: "
                          << vm["ub-file"].as<std::vector<string> >().at(0) << "\n";
                std::cerr << "Warning: best solution may be wrong, since explicit upper bound has been given.\n";
                fs::path filePath(vm["ub-file"].as<std::vector<string> >().at(0));
                ub_file = filePath;
                ub_file_ok = true;
            }
            fs::path cuts_file;
            if (vm.count("cuts-file")) {
                LOG(INFO) << "Model cuts input file is: "
                                        << vm["cuts-file"].as<std::vector<string> >().at(0);
                std::cout << "Model cuts input file is: "
                          << vm["cuts-file"].as<std::vector<string> >().at(0) << "\n";
                fs::path filePath(vm["cuts-file"].as<std::vector<string> >().at(0));
                cuts_file = filePath;
            }

            std::vector<fs::path> fileList;
            fs::path nodeseltrj_file, nodeprutrj_file, solfname_file, outputsolfname_file;
            if (vm.count("input-file")) {
                LOG(INFO) << "Input file is: "
                                        << vm["input-file"].as<std::vector<string> >().at(0) << "\n";
                fs::path filePath(vm["input-file"].as<std::vector<string> >().at(0));
                fileList.push_back(filePath);
            } else if (vm.count("f")) {
                LOG(INFO) << "Input file is: "
                                        << vm["f"].as<std::vector<string> >().at(0) << "\n";
                fs::path filePath(vm["f"].as<std::vector<string> >().at(0));
                fileList.push_back(filePath);
            } else if (vm.count("input-file-dir")) {
                LOG(INFO) << "Input file dir is: " << inputFileDir << endl;
                fs::path inputDir(inputFileDir);
                fs::directory_iterator end_iter;
                if (not(fs::is_directory(inputDir) or fs::is_symlink(inputDir)) and (not fs::exists(inputDir))) {
                    LOG(FATAL) << "Input file directory not found. Exiting." << endl;
                    return 1;
                }
                for (fs::directory_iterator dir_iter(inputDir); dir_iter != end_iter; ++dir_iter) {
                    string ext = dir_iter->path().extension().string();
                    boost::algorithm::to_lower(ext);
                    if (fs::is_regular_file(dir_iter->status())) {
                        fs::path filePath = *dir_iter;
                        fileList.push_back(filePath);
                    }
                }
                LOG(INFO) << "Will process " << fileList.size() << " files.";
            } else {
                cout << "Please specify at least one input file.\n";
                LOG(FATAL) << "Please specify at least one input file.";
                return 1;
            }
            if (solrequired and not(vm.count("solfname") or vm.count("o"))) {
                cout << "Missing optimal solution file\n";
                LOG(FATAL) << "Missing optimal solution file";
                return 1;
            }
            if (vm.count("solfname") or vm.count("o")) {
                LOG(INFO) << "Solution input file is: " << solfname << "\n";
                fs::path filePath(solfname);
                if (not testInputFile(filePath)) {
                    LOG(FATAL) << "Please specify a correct optimal solution file.";
                    cout << "Please specify a correct optimal solution file.";
                    return 1;
                }
                solfname_file = filePath;
            }
            // Random seed config
            if (the_seed == 0) {
                // default random seed used in the algorithms
                /*
                * Caveat: std::time(0) is not a very good truly-random seed.  When
                * called in rapid succession, it could return the same values, and
                * thus the same random number sequences could ensue.
                * Instead, we are using boost::random_device
                * http://stackoverflow.com/questions/4329284/c-boost-random-numeric-generation-problem
                * http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
                */
                the_seed = mix(clock(), time(NULL), getpid());
            }
            RandomUtil randomUtil(the_seed);
            randomUtil.setSeed(the_seed);
            LOG(INFO) << "Random seed = " << the_seed;
            cout << "Random seed = " << the_seed << "\n";
            if (parametrize) {
                LOG(INFO)
                    << "Upper bound parametrization mode ON. Solution will NOT be exported to file.";
                std::cout << "WARN: Upper bound parametrization mode ON. Solution will NOT be exported to file.\n";
            }

            if (vm.count("ub")) {
                if (upper_bound > 0.0) {
                    LOG(INFO) << "UB fixed to " << upper_bound;
                    cout << "UB fixed to " << upper_bound << "\n";
                    std::cerr << "Warning: best solution may be wrong, since explicit upper bound has been given.\n";
                    LOG(INFO)
                        << "Warning: best solution may be wrong, since explicit upper bound has been given.";
                }
            }
            if (vm.count("lb")) {
                if (lower_bound > 0.0) {
                    LOG(INFO) << "Root LB fixed to " << lower_bound;
                    cout << "Root LB fixed to " << lower_bound << "\n";
                    std::cerr
                            << "Warning: best solution may be wrong, since explicit root lower bound has been given.\n";
                    LOG(INFO)
                        << "Warning: best solution may be wrong, since explicit root lower bound has been given.";
                }
            }
            if (model_version >= 100) {
                robust = true;
            }
            NodeSelectionStrategy node_sel_strategy;
            BranchingType branching_type;
            if (vm["sel-strategy"].as<string>() == "dfs") {
                node_sel_strategy = NodeSelectionStrategy::dfs;
                LOG(INFO) << "Node selection strategy is DFS (Ladhari default).";
                cout << "Node selection strategy is DFS (Ladhari default).\n";
            } else if (vm["sel-strategy"].as<string>() == "cplex") {
                node_sel_strategy = NodeSelectionStrategy::cplex;
                LOG(INFO) << "Node selection strategy is CPLEX (default).";
                cout << "Node selection strategy is CPLEX (default).\n";
            } else {
                LOG(ERROR) << "Unknown Node selection strategy!";
                cerr << "Unknown Node selection strategy!";
                return -1;
            }
            if (not robust) {
                if (vm["branchingtype"].as<string>() == "alt") {
                    branching_type = BranchingType::alternating;
                    LOG(INFO) << "Branching type is alternating.";
                    cout << "Branching type is alternating.\n";
                } else if (vm["branchingtype"].as<string>() == "smd") {
                    branching_type = BranchingType::smallest_dynamic;
                    LOG(INFO) << "Branching type is smallest dynamic.";
                    cout << "Branching type is smallest dynamic.\n";
                } else {
                    LOG(ERROR) << "Unknown Branching type!";
                    cerr << "Unknown Branching type!";
                    return -1;
                }
            }
            // Initial permutation (used by BB algorithm)
            string initial_permutation_str = vm["initial-permutation"].as<string>();
            if (initial_permutation_str.size() > 0) {
                std::istringstream is(initial_permutation_str);
                initial_permutation = std::vector<int>(std::istream_iterator<int>(is), std::istream_iterator<int>());
                stringstream permutation_ss;
                for (int x : initial_permutation) {
                    if (x <= 0) {
                        LOG(ERROR)
                            << "ERROR: Values in initial permutation are NOT in the required range of j >= 1.";
                        cout << "ERROR: Values in initial permutation are NOT in the required range of j >= 1.\n";
                        return -1;
                    }
                    permutation_ss << x << " ";
                }
                LOG(INFO) << "[BB] Initial permutation is: " << permutation_ss.str();
                cout << "[BB] Initial permutation is: " << permutation_ss.str() << ".\n";
            }

            MIPUsage mip_usage = MIPUsage::no_mip;
            if (not robust) {
                LOG(INFO) << "Upper bound parameter values" << "===============================";
                LOG(INFO) << "GRASP time_factor=" << grasp_tl;
                cout << "GRASP time_factor=" << grasp_tl << "\n";
                LOG(INFO) << "GRASP max_iter=" << grasp_maxiter;
                cout << "GRASP max_iter=" << grasp_maxiter << "\n";

                LOG(INFO) << "LB strategy is "
                                        << (lb_strategy == 0 ? "Ladhari default." : (lb_strategy == 1 ? "Ladhari A1."
                                                                                                      : "Ritt."));
                cout << "LB strategy is "
                     << (lb_strategy == 0 ? "Ladhari default." : (lb_strategy == 1 ? "Ladhari A1." : "Ritt.")) << "\n";
                LOG(INFO) << "Use reverse instance ? " << reverse;
                cout << "Use reverse instance? " << reverse << "\n";
                LOG(INFO) << "Calculate Upper Bounds asynchronously ? " << async_ub;
                cout << "Calculate Upper Bounds asynchronously ? " << async_ub << "\n";
                LOG(INFO) << "Use combinatorial bounds inside MIP ? " << mip_use_combinatorial_bounds;
                cout << "Use combinatorial bounds inside MIP ? " << mip_use_combinatorial_bounds << "\n";
                LOG(INFO) << "Dominance rules are enabled ? " << dominance_rules;
                cout << "Dominance rules are enabled ? " << dominance_rules << "\n";

                if (vm["mip-usage"].as<string>() == "no-mip") {
                    mip_usage = MIPUsage::no_mip;
                    LOG(INFO) << "MIP Usage is NO MIP.";
                    cout << "MIP Usage is NO MIP.\n";
                } else if (vm["mip-usage"].as<string>() == "linear-relax") {
                    mip_usage = MIPUsage::linear_relax;
                    LOG(INFO) << "MIP Usage is Linear Relaxation. Depth = " << mip_depth;
                    cout << "MIP Usage is Linear Relaxation. Depth = " << mip_depth << "\n";
                } else if (vm["mip-usage"].as<string>() == "full-mip") {
                    mip_usage = MIPUsage::full_mip;
                    LOG(INFO) << "MIP Usage is Full MIP. Depth = " << mip_depth;
                    cout << "MIP Usage is Full MIP. Depth = " << mip_depth << "\n";
                } else {
                    LOG(ERROR) << "Unknown MIP Usage strategy! Using default no mip.";
                    cerr << "Unknown MIP Usage strategy! Using default no mip.\n";
                }
            } else {
                LOG(INFO) << "Requested the solution of a ROBUST PFS Problem.";
                cout << "Requested the solution of a ROBUST PFS Problem.\n";
                // Robust PFSP Parameters
                if (model_version >= 103
                    or model_version >= 604) {  // Robust Minimax Budget
                    string budget_gamma_str = vm["budget-gamma"].as<string>();
                    if (budget_gamma_str.size() == 0) {
                        LOG(ERROR) << "Budget parameter value is required !";
                        cerr << "Budget parameter value is required !\n";
                        return -1;
                    }
                    // parse values
                    std::istringstream is(budget_gamma_str);
                    budget_gamma = std::vector<double>(std::istream_iterator<double>(is),
                                                     std::istream_iterator<double>());
                    stringstream budget_ss;
                    for (double x : budget_gamma) {
                        if (x < 0 - 1e-5 or x > 100.0 + 1e-5) {
                            LOG(ERROR)
                                << "ERROR: Values in budget parameter Gamma are NOT in the required range of [0, 100].";
                            cout << "ERROR: Values in budget parameter Gamma NOT in the required range of [0, 100].\n";
                            return -1;
                        }
                        budget_ss << x << " ";
                    }
                    LOG(INFO) << "[Robust Minimax Budget] The budget parameter is: " << budget_ss.str();
                    cout << "[Robust Minimax Budget] The budget parameter Gamma is: " << budget_ss.str() << ".\n";

                    budget_type = vm["budget-type"].as<string>();
                    LOG(INFO) << "[Robust Minimax Budget] The budget type is: " << budget_type;
                    cout << "[Robust Minimax Budget] The budget type is: " << budget_type << ".\n";

                    LOG(INFO) << "[Robust Minimax Budget] Max number of Upper Bound iterations ub_iter_max="
                                            << ub_iter_max;
                    cout << "[Robust Minimax Budget] Max number of Upper Bound iterations ub_iter_max=" << ub_iter_max << "\n";
                    if (vm["rob-ub-type"].as<string>() == "grasp") {
                        rob_ub_type = Robust_UB_Type::grasp;
                        LOG(INFO) << "Robust upper bound type is GRASP.";
                        cout << "Robust upper bound type is GRASP.\n";
                    } else {
                        LOG(ERROR) << "Unknown Robust upper bound type!";
                        cerr << "Unknown Robust upper bound type!";
                        return -1;
                    }
                    // Rob-PFSP Hybrid MP type
                    mp_model_name_str = vm["mp-model-name"].as<string>();
                    LOG(INFO) << "[Robust Minimax Budget] The MP model type is: " << mp_model_name_str;
                    cout << "[Robust Minimax Budget] The MP model type is: " << mp_model_name_str << "\n";
                }
            }
            // Parse general metaheristic parameters (VND Neighborhood)
            if((not robust) or (vm["rob-ub-type"].as<string>() != "sa" and vm["rob-ub-type"].as<string>() != "ig")) {
                LOG(INFO) << "[VND] Variable Neighborhood Descent (VND) size is " << vnd_size;
                cout << "[VND] Variable Neighborhood Descent (VND) size is " << vnd_size << "\n";
                LOG(INFO) << "[VND] Variable Neighborhood Descent (VND) first improvement is enabled ? "
                                        << first_improvement;
                cout << "[VND] Variable Neighborhood Descent (VND) first improvement is enabled ? " << first_improvement
                     << "\n";
                string vnd_permutation_str = vm["vnd-permutation"].as<string>();
                vnd_permutation = std::vector<unsigned>(vnd_size, 0);
                if (vnd_permutation_str.size() == 0) {  // use default value for vnd_permutation : [1 2 3 ...]
                    std::iota(std::begin(vnd_permutation), std::end(vnd_permutation), 1);
                    std::ostringstream vts;
                    std::copy(vnd_permutation.begin(), vnd_permutation.end(),
                              std::ostream_iterator<unsigned>(vts, ", "));
                    vnd_permutation_str = vts.str();
                } else {  // parse values
                    std::istringstream is(vnd_permutation_str);
                    vnd_permutation = std::vector<unsigned>(std::istream_iterator<unsigned>(is),
                                                            std::istream_iterator<unsigned>());
                }
                LOG(INFO) << "[VND] The neighborhood permutation is : " << vnd_permutation_str;
                cout << "[VND] The neighborhood permutation is : " << vnd_permutation << ".\n";
                if(random_vnd) {
                    LOG(INFO) << "[VND] Random VND is enabled.";
                    cout << "[VND] Random VND is enabled." << "\n";
                } else {
                    LOG(INFO) << "[VND] Random VND is disabled.";
                    cout << "[VND] Random VND is disabled." << "\n";
                }
                if(adaptive_construction) {
                    LOG(INFO) << "Metahuristic Adaptive construction phase is enabled.";
                    cout << "Metahuristic Adaptive construction phase is enabled." << "\n";
                } else{
                    LOG(INFO) << "Metahuristic Adaptive construction phase is disabled. Construction alpha in [" << const_beta1 << ", " << const_beta2 << "].";
                    cout << "Metahuristic Adaptive construction phase is disabled. Construction alpha in [" << const_beta1 << ", " << const_beta2 << "]." << "\n";
                }
                LOG(INFO) << "Metahuristic objective function update frequency is: " << obj_update_freq;
                cout << "Metahuristic objective function update frequency is: " << obj_update_freq << "\n";
            }
            // -------------------  I N S T A N C E     F I L E S     P R O C E S S I N G -------------------------
            BranchBound_Parameters bb_params(lower_bound, upper_bound, lb_strategy,
                                             node_sel_strategy, branching_type, ub_file_ok ? &ub_file : NULL,
                                             mip_usage, mip_depth, mip_use_combinatorial_bounds, async_ub,
                                             initial_permutation, cuts_file);
            Robust_Parameters rob_params(dominance_rules, rob_ub_type, ub_iter_max, budget_gamma, budget_type, mp_model_name_str);
            UpperBound_Parameters ub_params(the_seed, const_beta1, const_beta2, grasp_tl, grasp_maxiter, 
                    vnd_permutation, vnd_size, first_improvement, random_vnd, adaptive_construction, parametrize, obj_update_freq);
            for (unsigned int i = 0; i < fileList.size(); i++) {
                fs::path filePath = fileList.at(i);
                cout << "\nProcessing file: " << filePath << "\n";
                PFSP_Parameters pfsp_params(filePath, outputFolder, executionId, timeLimit, model_version, reverse, max_cores);
                processInputFile(pfsp_params, bb_params, rob_params, ub_params, robust);
            }
            LOG(INFO) << "Done.";
            cout << "Done.\n";
        }
        catch (std::exception &e) {
            cerr << "Fatal application error.\n";
            LOG(FATAL) << "Abnormal program termination. Stracktrace: " << endl;
            LOG(FATAL) << e.what() << "\n";
            if (std::string const *stack = boost::get_error_info<stack_info>(e)) {
                LOG(FATAL) << stack << endl;
            }
            LOG(FATAL) << diagnostic_information(e);
            return 1;
        }
        return 0;
    }

    void CommandLineInterfaceController::readPropertiesFile() {
        namespace pt = boost::property_tree;
        pt::ptree propTree;

        if (!boost::filesystem::exists("config.ini")) {
            cout << "Can't find config.ini file! Assuming default properties." << endl;
        } else {
            read_ini("config.ini", propTree);

            boost::optional<string> severity = propTree.get_optional<std::string>("logging.severity");
            if (severity) {
                logSeverity = *severity;
            } else {
                cout << "WARNING: warn log level not specified, assuming warn level." << endl;
            }
        }
    }

    void CommandLineInterfaceController::initLogging(const string &outputFolder, const string &executionId, char *argv[]) {
        cout << "Init logging..." << endl;

        // Initialize Google’s logging library.
        FLAGS_logtostderr = 0;
        FLAGS_alsologtostderr = 0;
        FLAGS_minloglevel = 0;
        if(outputFolder.size() > 0) {
            FLAGS_log_dir = outputFolder + boost::filesystem::path::preferred_separator + string("log");
        } else {
            FLAGS_log_dir = "./log";
        }
        cout << "Writing log files to: " << FLAGS_log_dir << endl;
        string filename = FLAGS_log_dir + boost::filesystem::path::preferred_separator + executionId + string(".log");
        fs::create_directories(FLAGS_log_dir);
        google::SetLogDestination(0, filename.c_str());
        google::InitGoogleLogging("pfsp");
    }

    void CommandLineInterfaceController::handler() {
        void *trace_elems[20];
        int trace_elem_count(backtrace(trace_elems, 20));
        char **stack_syms(backtrace_symbols(trace_elems, trace_elem_count));
        for (int i = 0; i < trace_elem_count; ++i) {
            LOG(FATAL) << stack_syms[i] << "\n";
        }
        free(stack_syms);
    }


} // namespace controller
