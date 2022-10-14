//
// Created by mlevorato on 8/15/19.
//

#ifndef FLOWSHOP_SOLVER_PFSPSCENARIO_H
#define FLOWSHOP_SOLVER_PFSPSCENARIO_H

#include <boost/numeric/ublas/matrix.hpp>

#include "../GlobalTypes.h"

namespace robust {

    class PFSPScenario {
    public:
        PFSPScenario(boost::numeric::ublas::matrix<Time> p, Time v) : p_time(p), value(v) {  }

        PFSPScenario() : p_time(), value(0) {  }


        boost::numeric::ublas::matrix<Time> p_time;
        /** objective value of the optimal solution for this scenario */
        Time value;  // can be worst cmax or worst deviation from optimum cmax (regret)
    };

}

#endif //FLOWSHOP_SOLVER_PFSPSCENARIO_H
