//
// Created by mlevorato on 1/14/20.
//

#ifndef PFSP_PFSP_CMAX_BUDGET_NODEDATA_H
#define PFSP_PFSP_CMAX_BUDGET_NODEDATA_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "../../../RobPFSInstance_Cmax.h"

namespace robust {

    using namespace std;
    using namespace boost;

    /**
     * Branch-and-bound node data for the robust PFSP budget.
     */
    class PFSP_Cmax_Budget_NodeData {
    public:
        PFSP_Cmax_Budget_NodeData() : mem_lb_sum(1, 0), mem_lb_count(1, 0) {

        }

        PFSP_Cmax_Budget_NodeData(RobPFSInstance_Cmax &instance) : mem_lb_sum(instance.n, 0), mem_lb_count(instance.n, 0) {

        }

        // Use a memory (function parameter) that, for each job j, stores the sum of lower bound values and the
        //  total number of times the job j was chosen to assume its worst-case proc. time value.
        //  For each job j ==> [ sum(LB) ; sum(count) ]  ( LB memory )
        //  The memory in inherited from the parent B & B node.
        boost::numeric::ublas::vector<double> mem_lb_sum;
        boost::numeric::ublas::vector<unsigned> mem_lb_count;
    };
}

#endif //PFSP_PFSP_CMAX_BUDGET_NODEDATA_H
