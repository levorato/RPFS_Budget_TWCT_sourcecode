//
// Created by mlevorato on 11/25/19.
//

#ifndef FLOWSHOP_SOLVER_PFSP_PARAMETERS_H
#define FLOWSHOP_SOLVER_PFSP_PARAMETERS_H

#include <string>
#include <boost/filesystem.hpp>
#include "GlobalTypes.h"

namespace parameters {
    using namespace std;
    using namespace boost;
    namespace fs = boost::filesystem;

    class BranchBound_Parameters {
    public:
        BranchBound_Parameters(const double &p_lower_bound, const double &p_upper_bound,
                               const int &p_lb_strategy, NodeSelectionStrategy p_node_sel_strategy,
                               BranchingType p_branching_type, fs::path *p_ub_path,
                               MIPUsage p_mip_usage, const unsigned &p_mip_depth,
                               const bool &p_use_combinatorial_bounds, const bool &p_async_ub,
                               const std::vector<int>& p_permutation, const fs::path &p_cuts_file) :
                lower_bound(p_lower_bound), upper_bound(p_upper_bound), lb_strategy(p_lb_strategy),
                node_sel_strategy(p_node_sel_strategy), branching_type(p_branching_type),
                ub_path(p_ub_path), mip_usage(p_mip_usage),
                mip_depth(p_mip_depth), use_combinatorial_bounds(p_use_combinatorial_bounds),
                async_ub(p_async_ub), async_ub_poll_interval(200), async_ub_push_interval(200),
                permutation(p_permutation), cuts_file(p_cuts_file) {

        }

        double lower_bound;
        double upper_bound;
        std::vector<int> permutation;
        int lb_strategy;
        NodeSelectionStrategy node_sel_strategy;
        BranchingType branching_type;
        fs::path *ub_path;
        MIPUsage mip_usage;
        unsigned mip_depth;
        bool use_combinatorial_bounds;
        bool async_ub;
        unsigned async_ub_poll_interval;
        unsigned async_ub_push_interval;
        fs::path cuts_file;
    };

    class PFSP_Parameters {
    public:
        PFSP_Parameters(const fs::path &p_filePath, const string &p_outputFolder, const string &p_executionId,
                        const double &p_timeLimit, const int &p_model_version, 
                        const bool &p_reverse, const unsigned &p_max_cores)
                :
                filePath(p_filePath),
                outputFolder(p_outputFolder), executionId(p_executionId), timeLimit(p_timeLimit),
                model_version(p_model_version), reverse(p_reverse), max_cores(p_max_cores) {

        }

        fs::path filePath;
        string outputFolder;
        string executionId;
        double timeLimit;
        int model_version;
        bool reverse;
        unsigned max_cores;  // max number of CPU cores to use in model solution
    };

    class Robust_Parameters {
    public:
        Robust_Parameters(const bool &p_dominance_rules, const Robust_UB_Type p_ub_type, const unsigned &p_ub_iter_max,
                          const std::vector<double> &p_budget_gamma,
                          const string &p_budget_type, const string &p_mp_model_name):
                dominance_rules(p_dominance_rules), ub_type(p_ub_type), ub_iter_max(p_ub_iter_max), 
                budget_gamma(p_budget_gamma),
                budget_type(p_budget_type), mp_model_name(p_mp_model_name) {
        }

        // Enable dominance rules for combinatorial branch-and-bound
        bool dominance_rules;
        Robust_UB_Type ub_type;  // upper bound type : sa (simulated annealing) or ig (iterated greedy)
        unsigned ub_iter_max;  // maximum number of iterations (SA, IG)
        std::vector<double> budget_gamma;  // Budget parameter: budget_gamma[1 .. n](machine) or budget_gamma[1](global)
        string budget_type;  // Budget type: per 'machine' or 'global' (single budget parameter for all operations)
        string mp_model_name;  // Rob-PFSP Hybrid MP type: e.g. hybrid-wilson, hybrid-liao-you.
    };

    class UpperBound_Parameters {
    public:
        UpperBound_Parameters(const unsigned long &p_seed, 
                  const double &p_beta1, const double &p_beta2, const unsigned &p_grasp_tl,
                  const unsigned &p_grasp_maxiter, 
                  const std::vector<unsigned> &p_vnd_permutation, const unsigned &p_vnd_size,
                  const bool &p_first_improvement, const bool &p_random_vnd, const bool &p_adaptive_construction,
                  const bool& p_parametrize, const long &p_obj_update_freq) :
            seed(p_seed), grasp_tl(p_grasp_tl),
            beta1(p_beta1), beta2(p_beta2),
            grasp_maxiter(p_grasp_maxiter), 
            vnd_permutation(p_vnd_permutation), vnd_size(p_vnd_size), first_improvement(p_first_improvement),
            random_vnd(p_random_vnd), adaptive_construction(p_adaptive_construction), parametrize(p_parametrize),
            obj_update_freq(p_obj_update_freq), validate_obj(true) {

        }

        double beta1, beta2;  // constructive phase params, if not adaptive
        unsigned grasp_tl;
        unsigned grasp_maxiter;
        unsigned long seed;
        // Constructive phase parameters
        bool adaptive_construction;  // See RandNEHT class for further info
        // VND local search parameters
        std::vector<unsigned> vnd_permutation;
        unsigned vnd_size;
        bool first_improvement;
        bool random_vnd;
        // enable metaheuristic parametrization mode
        bool parametrize;
        long obj_update_freq;  // frequency to recalculate the objective function (used by metaheuristics)
        bool validate_obj;  // validate final objective function value by fully recalculating it
    };
}

#endif //FLOWSHOP_SOLVER_PFSP_PARAMETERS_H
