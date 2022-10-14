//
// Created by mlevorato on 7/12/19.
//

#ifndef FLOWSHOP_SOLVER_CPLEX_HYBRID_PFSP_WCT_H
#define FLOWSHOP_SOLVER_CPLEX_HYBRID_PFSP_WCT_H

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>

#include "../../../../PFSP_Parameters.h"
#include "../../../RobPFSInstance.h"
#include "../../../../PFSSolution.h"
#include "PFSP_Budget_WCT_BranchBound.h"

namespace robust {
    namespace hybrid {

        using namespace problem::common;
        namespace ublas = boost::numeric::ublas;
        using namespace parameters;
        using namespace robust;
        using namespace std;
        namespace fs = boost::filesystem;

        class PFSP_WCT_Budget_Hybrid_BranchBound {
        public:
            PFSP_WCT_Budget_Hybrid_BranchBound() : scenarios_calculated(false),
                P_time_scenario_list(), l_kl_scenario_list() { }

            void solve(RobPFSInstance_WCT &instance, 
                    const string &output_folder,
                    const PFSP_Parameters &pfsp_params,
                    const BranchBound_Parameters &bb_params,
                    const Robust_Parameters &rob_params, 
                    const UpperBound_Parameters &ub_params);

            PFSSolution run(RobPFSInstance_WCT &instance, const string &output_folder,
                            const PFSP_Parameters &pfsp_params,
                            const BranchBound_Parameters &bb_params,
                            const Robust_Parameters &rob_params,
                            const UpperBound_Parameters &ub_params);

            std::vector< ublas::matrix<double> > read_model_cuts_from_file(const fs::path &cuts_file,
                    const int &m, const int &n);

            problem::common::PFSSolution obtainInitialSolution(RobPFSInstance_WCT &instance, const double &lb,
                                                                const parameters::PFSP_Parameters &pfsp_params,
                                                                const Robust_Parameters &rob_params, 
                                                                const UpperBound_Parameters &ub_params);

            void printResult();

            void calculate_scenarios(const int &m, const int &n,
                const std::vector< ublas::matrix<double> > &cut_list, const robust::RobPFSInstance_WCT &instance);
                
            private:
                // The two members below store pre-calculations of values used in lower bound / dominance rules
                bool scenarios_calculated;
                std::vector< ublas::matrix<double> > P_time_scenario_list;
                std::vector< boost::multi_array<double, 3> > l_kl_scenario_list;
        };
    } /* namespace hybrid */
} /* namespace robust */


#endif //FLOWSHOP_SOLVER_CPLEX_HYBRID_PFSP_WCT_H
