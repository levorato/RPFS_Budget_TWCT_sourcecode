//
// Created by mlevorato on 4/30/19.
//

#ifndef FLOWSHOP_SOLVER_SELECTIONSTRATEGY_H
#define FLOWSHOP_SOLVER_SELECTIONSTRATEGY_H

namespace problem {
    namespace pfsp {
        namespace bb {

/**
 * Base strategy: select a node solely based on the smallest lower bound.
 */
template<class Node>
class DescendantNodeSelectionStrategy {
    bool reverse;
public:
    DescendantNodeSelectionStrategy(const bool &revparam = false) { reverse = revparam; }

    // Order by the element which has the smallest lower bound value
    // In case of tie: select the element with the biggest level of the BB tree
    bool operator()(const Node &lhs, const Node &rhs) const {
        // return n.lb > c.lb || (n.lb == c.lb && n.numFixed() < c.numFixed()
        return (lhs.LB > rhs.LB) || (lhs.LB == rhs.LB && lhs.level < rhs.level);
    }
};


        } /* namespace bb */
    } //* namespace pfsp */
} /* namespace problem */

#endif //FLOWSHOP_SOLVER_SELECTIONSTRATEGY_H
