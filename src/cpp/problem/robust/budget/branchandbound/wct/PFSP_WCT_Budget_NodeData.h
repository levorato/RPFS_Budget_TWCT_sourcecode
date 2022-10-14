//
// Created by mlevorato on 3/1/21.
//

#ifndef PFSP_PFSP_WCT_BUDGET_NODEDATA_H
#define PFSP_PFSP_WCT_BUDGET_NODEDATA_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "../../../RobPFSInstance_WCT.h"

namespace robust {

    using namespace std;
    using namespace boost;

    /**
     * Branch-and-bound node data for the robust PFSP budget.
     */
    class PFSP_WCT_Budget_NodeData {
    public:
        PFSP_WCT_Budget_NodeData() : mem_lb_sum(1, 0), mem_lb_count(1, 0),
            P_time_scenario_list(), l_kl_scenario_list() {

        }

        PFSP_WCT_Budget_NodeData(RobPFSInstance_WCT &instance) : mem_lb_sum(instance.n, 0), 
            mem_lb_count(instance.n, 0), P_time_scenario_list(), l_kl_scenario_list() {

        }

        // Copy constructor
        PFSP_WCT_Budget_NodeData(const PFSP_WCT_Budget_NodeData &nd) : mem_lb_sum(nd.mem_lb_sum), 
            mem_lb_count(nd.mem_lb_count),
            P_time_scenario_list(nd.P_time_scenario_list), l_kl_scenario_list(nd.l_kl_scenario_list) {

        }

        // Use a memory (function parameter) that, for each job j, stores the sum of lower bound values and the
        //  total number of times the job j was chosen to assume its worst-case proc. time value.
        //  For each job j ==> [ sum(LB) ; sum(count) ]  ( LB memory )
        //  The memory in inherited from the parent B & B node.
        boost::numeric::ublas::vector<double> mem_lb_sum;
        boost::numeric::ublas::vector<unsigned> mem_lb_count;
        // The two members below store pre-calculations of values used in lower bound / dominance rules
        std::vector< boost::numeric::ublas::matrix<double> > P_time_scenario_list;
        std::vector< boost::multi_array<double, 3> > l_kl_scenario_list;
    };
}

#endif //PFSP_PFSP_WCT_BUDGET_NODEDATA_H
