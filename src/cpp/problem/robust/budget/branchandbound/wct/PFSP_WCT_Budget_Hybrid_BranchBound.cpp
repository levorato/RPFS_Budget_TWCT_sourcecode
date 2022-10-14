#include "PFSP_WCT_Budget_Hybrid_BranchBound.h"
#include "./cplex/PFSP_WCT_Brancher.h"
// #include "./cplex/PFSP_WCT_Brancher_v2.h"
#include "./cplex/PFSP_WCT_DFSNodeHandler.h"
#include "./cplex/PFSP_WCT_BranchInfo.h"

#include <ilcplex/ilocplex.h>
#include <glog/logging.h>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o
#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
#include <regex>
#include <thread>
#include <algorithm>

#include "cplex/PFSP_WCT_ModelSelector.h"
#include "cplex/PFSP_WCT_CPLEXEnvironment.h"
#include "../../../../deterministic/PFSProblem.h"
#include "../BB_OutputUtil.h"


namespace robust {
    namespace hybrid {

    using namespace boost;
    using namespace multichild;
    namespace ublas = boost::numeric::ublas;
    namespace fs = boost::filesystem;


    void PFSP_WCT_Budget_Hybrid_BranchBound::solve(RobPFSInstance_WCT &instance, 
                    const string &output_folder,
                    const PFSP_Parameters &pfsp_params,
                    const BranchBound_Parameters &bb_params,
                    const Robust_Parameters &rob_params, 
                    const UpperBound_Parameters &ub_params) {
        LOG(INFO) << "[Robust PFSP Hybrid BB] Initialize...";
        LOG(INFO) << "[Robust PFSP Hybrid BB] MP model type is " << rob_params.mp_model_name;
        PFSSolution solution = run(instance, output_folder, pfsp_params, bb_params, rob_params, ub_params);
        BB_OutputUtil::print_robust_solution_summary(&instance, solution, 
                            pfsp_params.outputFolder, pfsp_params.executionId, rob_params.mp_model_name, ub_params);

        LOG(INFO) << "========================================================================";
        LOG(INFO) << "      ROBUST  BB (" << rob_params.mp_model_name  << ") FINISHED";
        LOG(INFO) << "========================================================================";
        LOG(INFO) << std::fixed << std::setprecision(10) << " Best solution value  = " << solution.value;
        LOG(INFO) << " Time spent: " << solution.time_spent << " s";
        std::cout << std::fixed << std::setprecision(10) << " Cost=" 
                  << solution.value << " ; Time(s) = " << solution.time_spent << "\n";
        LOG(INFO) << "========================================================================";
    }

    std::vector<string> splitString(string str)
    {
        istringstream iss(str);
        std::vector<string> list;
        for (auto it = istream_iterator<string>(iss); it != istream_iterator<string>(); ++it)
            list.push_back(*it);
        return list;
    }

    std::string& replace(std::string& s, const std::string& from, const std::string& to)
    {
        if(!from.empty())
            for(size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
                s.replace(pos, from.size(), to);
        return s;
    }

    std::vector< ublas::matrix<double> > PFSP_WCT_Budget_Hybrid_BranchBound::read_model_cuts_from_file(const fs::path &cuts_file,
                const int &m, const int &n) {
        std::vector< ublas::matrix<double> > cut_list;
        // the first cut of the list is always the nil set (corresponding to nominal processing times)
        ublas::matrix<double> first_cut_matrix = ublas::zero_matrix<double>(m + 1, n + 1);
        cut_list.push_back(first_cut_matrix);
        if(fs::exists(cuts_file)) {
            std::ifstream file(cuts_file.string());
            std::string str_line; 
            double v;
            while (std::getline(file, str_line) && !str_line.empty()) {
                // remove brackets from line
                str_line = replace(str_line, "[", "");
                str_line = replace(str_line, "]", "");
                str_line = replace(str_line, ",", "");
                istringstream iss(str_line);
                ublas::matrix<double> cut_matrix = ublas::zero_matrix<double>(m + 1, n + 1);
                int i = 1, j = 1;
                std::vector<string> str_list = splitString(str_line);
                stringstream ss;
                for(string s : str_list) {
                    sscanf(s.c_str(), "%lf", &v);
                    if(fabs(v) <= 1e-5)  v = 0.0;
                    cut_matrix(i, j) = v;
                    ss << v << ", ";
                    j++;
                    if(j > n) {
                        j = 1;
                        i++;
                    }
                    if(i > m) {
                        break;
                    }
                }
                LOG(INFO) << "Line Split: " << ss.str();
                LOG(INFO) << "Line String: " << str_line;
                LOG(INFO) << "Line Matrix: " << cut_matrix;
                cut_list.push_back(cut_matrix);
            }
        }
        return cut_list;
    }

    PFSSolution PFSP_WCT_Budget_Hybrid_BranchBound::run(RobPFSInstance_WCT &instance, 
                    const string &output_folder,
                    const PFSP_Parameters &pfsp_params,
                    const BranchBound_Parameters &bb_params,
                    const Robust_Parameters &rob_params, 
                    const UpperBound_Parameters &ub_params) {
        // problem parameters
        int m = instance.m;
        int n = instance.n;
        PFSP_WCT_ModelSelector my_model(rob_params.mp_model_name);
        double time_limit = pfsp_params.timeLimit;
        LOG(INFO) << "Start PFSP-WCT Hybrid BB...";
        boost::timer::cpu_timer timer;

        try {  // Create an IloCplex instance and load the model.
            IloEnv env;
            IloModel model(env);
            IloCplex cplex(model);
            // cplex.setOut(env.getNullStream());
            /* Sets the rule for selecting the next node to process when the search is backtracking.
                The depth-first search strategy chooses the most recently created node.
                The best-bound strategy chooses the node with the best objective function for the associated LP relaxation.
                The best-estimate strategy selects the node with the best estimate of the integer objective value that would be
                obtained from a node once all integer infeasibilities are removed. An alternative best-estimate search is also available. */
            /*  0   CPX_NODESEL_DFS 	Depth-first search
                1 	CPX_NODESEL_BESTBOUND 	Best-bound search; default
                2 	CPX_NODESEL_BESTEST 	Best-estimate search
                3 	CPX_NODESEL_BESTEST_ALT */
            cplex.setParam( IloCplex::IntParam::NodeSel, CPX_NODESEL_BESTBOUND);  // CPX_NODESEL_BESTEST);  // IloCplex.NodeSelect.BestBound , BestEst, BestEstAlt, DFS

            // IMPORTANT: set a tolerance so that we don't split hairs
            cplex.setParam( IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1e-5);  /// 0.001 );
            // Setup model solve time limit
            cplex.setParam( IloCplex::NumParam::TiLim, time_limit );
            // IMPORTANT SET THE NUMBER OF THREADS CPLEX CAN USE: FORCES MULTITHREADING WHEN USING CALLBACKS
            unsigned processor_count = std::thread::hardware_concurrency();  // may return 0 when not able to detect
            if(processor_count == 0) {
                processor_count = pfsp_params.max_cores;
            }
            unsigned num_cores = std::min(pfsp_params.max_cores, processor_count);
            cplex.setParam( IloCplex::Param::Threads, num_cores );
            LOG(INFO) << "PFSP-WCT MILP: setting CPLEX thread count to " << num_cores;
            
            // Upper Bound calculation
            PFSSolution ub_sol;
            if (bb_params.upper_bound > 0.0) {
                ub_sol = PFSSolution();
                ub_sol.value = bb_params.upper_bound;
                ub_sol.permutation = bb_params.permutation;
            } else {
                ub_sol = obtainInitialSolution(instance, 0.0, pfsp_params, rob_params, ub_params);
            }
            LOG(INFO) << "[Hybrid-BB] PFSP-WCT MP upper bound calculated.";

            // Read model cut list from input cuts file
            std::vector< ublas::matrix<double> > cut_list = read_model_cuts_from_file(bb_params.cuts_file, m, n);
            IloRangeArray con(env);
            IloNumVar obj(env);
            const std::vector<int> &initial_permutation = ub_sol.getPermutation();
            NumVarMatrix bin_vars = my_model.populate_model (cplex, model, con, obj, m, n, instance.p_bar, instance.p_hat, 
                    instance.w, initial_permutation, cut_list);
            LOG(INFO) << "PFSP-WCT " << rob_params.mp_model_name << " MP MILP Model populated.";
            
            // Register an instance of the branch callback with IloCplex.
            bool verb = false;
            // Pre-calculate all scenario processing time matrices
            calculate_scenarios(instance.m, instance.n, cut_list, instance);
            // Now we get to setting up the callback.
            // We instantiate a PFSP_WCT_Brancher_v2 and set the wherefrom parameter.
            int numThreads = cplex.getNumCores();
            /*
            // NEW CPLEX BRANCH CALLBACK (GENERIC)
            PFSP_WCT_Brancher_v2 cb(cplex, numThreads, obj, z, instance, verb, 
                    bb_params.use_combinatorial_bounds, output_folder, 
                    true, ub_sol, bb_params.lower_bound, P_time_scenario_list, l_kl_scenario_list);
            CPXLONG contextmask = IloCplex::Callback::Context::Id::Branching
                            | IloCplex::Callback::Context::Id::ThreadUp
                            | IloCplex::Callback::Context::Id::ThreadDown;
            // We add the callback.
            cplex.use(&cb, contextmask);
            */

            //  ORIGINAL CODE CALLING CALLBACK v1 
            PFSP_WCT_Brancher *brancher = new PFSP_WCT_Brancher(env, obj, bin_vars, instance, verb, 
                    bb_params.use_combinatorial_bounds, output_folder, true, ub_sol, bb_params.lower_bound, 
                    P_time_scenario_list, l_kl_scenario_list, rob_params.mp_model_name, rob_params.dominance_rules);
            cplex.use(brancher);
            double UB = 0.0;
            cplex.use(new (env) PFSP_WCT_DFSNodeHandler(env, n, bb_params.node_sel_strategy, &UB, verb));
            cplex.use(brancher);
            
            // TODO Add PFSP_WCT_DFSNodeHandler to the new callback PFSP_WCT_Brancher_v2
            ///// cplex.use(new (env) PFSP_WCT_DFSNodeHandler(env, n, bb_params.node_sel_strategy, &PFSP_WCT_Brancher::UB, verb));

            // Solve the model.
            LOG(INFO) << "Solving the " << rob_params.mp_model_name << " model through CPLEX..." << std::endl;
            PFSSolution solution;
            if ( cplex.solve() ) {
                LOG(INFO) << "CPLEX solve completed.";
                LOG(INFO) << "[Hybrid PFSP model] Objective: " << std::fixed << std::setprecision(10) << cplex.getObjValue() << std::endl;
                // fetch the solution
                std::vector<int> permutation = my_model.getSolutionAsPermutation(cplex, model, n, bin_vars);
                LOG(INFO) << "\nPermutation: " << permutation << "\n";
                solution = PFSSolution(permutation, cplex.getObjValue(), n);
                LOG(INFO) << "CPLEX explored " << cplex.getNnodes() << " nodes.\n";
                double validation_wct = PFSProblem::calculate_wct(permutation, instance.m, instance.n, instance.p_bar, instance.w);
                LOG(INFO) << "Nominal WCT Cost = " << validation_wct;
                solution.validated = true;
                solution.time_spent = TimeDateUtil::calculate_time_spent(timer);
                solution.num_iterations = cplex.getNnodes();
                solution.num_improvements = brancher->getNumberOfImprovements();
                solution.gap = cplex.getMIPRelativeGap();
                solution.is_optimal = fabs(solution.gap) <= 1e-3;
            } else {
                LOG(INFO) << "MODEL INFEASIBLE" << std::endl;
            }
            cplex.end();
            model.end();
            env.end();
            return solution;
        } catch (IloException &e) {
            std::cerr << "IloException: " << e << std::endl;
            LOG(ERROR) << "IloException: " << e;
            throw;
        } catch (...) {
            std::cerr << "Unknown exception in CPLEX !\n";
            LOG(ERROR) << "Unknown exception in CPLEX !";
            // ex.printStackTrace();
            throw;
        }
    }

    /**
     * For each scenario in cut_list, calculate its corresponding processing time matrix.
     * Store the list of matrices in P_time_scenario_list.
     */
    void PFSP_WCT_Budget_Hybrid_BranchBound::calculate_scenarios(const int &m, const int &n,
            const std::vector< ublas::matrix<double> > &cut_list, const robust::RobPFSInstance_WCT &instance) {
        if(!scenarios_calculated) {
            for(int i = 0; i < cut_list.size(); i++) {
                ublas::matrix<double> sep = cut_list[i];
                // Calculation of processing time matrix for scenario i
                ublas::matrix<double> p_time(instance.p_bar);
                for(int r = 1; r <= instance.m; r++) {  // Machine r
                    for(int j = 1; j <= instance.n; j++) {  // Job j
                        p_time(r, j) += instance.p_hat(r, j) * sep(r, j);
                    }
                }
                P_time_scenario_list.push_back(p_time);
                // Calculation of l_kl matrix for scenario i
                // IMPORTANT: Initialization of l_kl auxiliary vectors used by lower bounds
                // For each machine pair (k, l)
                // Sum of processing times between machines k + 1 and l - 1 : Create a 3D array
                boost::multi_array<double, 3> l_kl(boost::extents[n + 2][m + 2][m + 2]);
                std::fill_n(l_kl.data(), l_kl.num_elements(), 0);
                for (int j = 1; j <= n; j++) {
                    for (int k = 1; k <= m; k++) {
                        for (int l = k; l <= m; l++) {
                            for (int u = k; u <= l; u++) {
                                l_kl[j][k][l] += p_time(u, j);
                            }
                        }
                    }
                }
                l_kl_scenario_list.push_back(l_kl);
            }
            scenarios_calculated = true;
        }
    }

    problem::common::PFSSolution PFSP_WCT_Budget_Hybrid_BranchBound::obtainInitialSolution(RobPFSInstance_WCT &instance, const double &lb,
            const parameters::PFSP_Parameters &pfsp_params, const Robust_Parameters &rob_params, const UpperBound_Parameters &ub_params) {
        LOG(INFO) << "[BB-Hybrid-PFSP-WCT] obtainInitialSolution: Current solution status = heuristic";
        pfsp::GRASPSolver_RobPFSP_WCT grasp_wct;
        LOG(INFO) << "[Robust PFSP Budget] Solving Robust PFSP WCT instance with GRASP...";
        instance.robust_type = robust::RobustType::budget;
        UpperBound_Parameters new_ub_params(ub_params);
        // disable GRASP objective validation if gamma == 0 or gamma == 100 (saves time)
        new_ub_params.validate_obj = false;
        problem::common::PFSSolution solution = grasp_wct.solve(instance, 0, pfsp_params, new_ub_params);
        LOG(INFO) << "[BB-Hybrid-PFSP-WCT] Best solution value  = " << solution.value << std::endl;
        return solution;
    }


    void PFSP_WCT_Budget_Hybrid_BranchBound::printResult() {
    /*
        IloNumArray vals(*env);
        env->out() << "Solution status = " << cplex->getStatus() << endl;
        env->out() << "Solution value  = " << cplex->getObjValue() << endl;
        env->out() << "GAP             = " << cplex->getMIPRelativeGap()*100 << "%"<< endl;

        cplex->getValues(vals, *variables);
        cout << "A = [ ";
        for( int i = 0; i < n; ++i )
            cout << vals[i] << " ";
        cout << "]" << endl;
        cout << "B = [ ";
        for( int i = 0; i < n; ++i )
            cout << vals[n+i] << " ";
        cout << "]" << endl;
        */
    }

    } /* namespace hybrid */
} /* namespace robust */
