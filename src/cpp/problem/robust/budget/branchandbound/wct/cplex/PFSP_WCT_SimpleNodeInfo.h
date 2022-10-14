//
// Created by mlevorato on 7/23/19.
//

#ifndef FLOWSHOP_SOLVER_SIMPLENODEINFO_H
#define FLOWSHOP_SOLVER_SIMPLENODEINFO_H

#include <ilcplex/ilocplex.h>
#include "../../../../RobPFSInstance_WCT.h"

namespace multichild {

    class PFSP_WCT_SimpleNodeInfo {
    public:
        robust::RobPFSInstance_WCT *instance;
        IloCplex::NodeCallbackI::NodeId nodeNumber;
        int level;
        double LB;

        PFSP_WCT_SimpleNodeInfo() : nodeNumber(), LB(0), level(0), instance(NULL) {
            nodeNumber._id = -1;
        }

        PFSP_WCT_SimpleNodeInfo(IloCplex::NodeCallbackI::NodeId num, double lb, int lvl, robust::RobPFSInstance_WCT *p_instance) :
            nodeNumber(num), LB(lb),
            level(lvl), instance(p_instance) {

        }

        PFSP_WCT_SimpleNodeInfo(long node_id, const std::vector<int> &v, int n, int m, robust::RobPFSInstance_WCT *p_instance) :
                                    nodeNumber(),
                                    LB(0), level(0), instance(p_instance) {
            nodeNumber._id = node_id;
        }

        unsigned long get_level() const {
            return level;
        }

        long get_id() const {
            return nodeNumber._id;
        }

        // CBFS BB Priority Queue criteria used by M. Ritt.
        friend bool operator<(const PFSP_WCT_SimpleNodeInfo &a, const PFSP_WCT_SimpleNodeInfo &b) {
            // return a.LB > b.LB || (a.LB == b.LB && a.meanTime() > b.meanTime());
            // return a.lb > b.lb || (a.lb == b.lb && a.numFixed() < b.numFixed()
            return (a.LB > b.LB) || (a.LB == b.LB && a.level < b.level);
        }

        double getLB() const {
            return LB;
        }

        std::vector<int> getPartialSequence() {
            return std::vector<int>();  // FIXME
        }
    };

}

#endif //FLOWSHOP_SOLVER_SIMPLENODEINFO_H
