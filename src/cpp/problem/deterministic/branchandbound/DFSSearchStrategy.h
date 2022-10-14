//
// Created by mlevorato on 4/29/19.
//

#ifndef FLOWSHOP_SOLVER_DFSSEARCHSTRATEGY_H
#define FLOWSHOP_SOLVER_DFSSEARCHSTRATEGY_H

#include <boost/heap/pairing_heap.hpp>
#include <cassert>
#include <chrono>
#include <list>
#include <iostream>
#include <iomanip>
#include <glog/logging.h>

#include "../../GlobalTypes.h"
#include "SelectionStrategy.h"
#include "../../PFSSolution.h"
#include "../PFSInstance.h"
#include "../../../util/NumericUtil.h"

namespace problem {
    namespace pfsp {
        namespace bb {
            using namespace std;
            using namespace boost;
            using namespace problem::common;
            using namespace util;

            template<class Node>
            class DFSSearchStrategy {
            public:
                DFSSearchStrategy(const int& n, std::chrono::system_clock::time_point start_time,
                                  const bool& p_force_lb_update) : num_levels(n),
                                    elementCount(0), llb(), force_lb_update(p_force_lb_update),
                                    last(std::chrono::system_clock::now()), start(start_time),
                                    leafCount(0) {

                }
                ~DFSSearchStrategy() {   }

                void add(Node c) {
                    auto firstNonChild = llb.begin();
                    llb.insert(find_if(llb.begin(), firstNonChild, [&c](const Node& n) {
                        return NumericUtil::TestisGT(n.getLB(), c.getLB())
                            || (NumericUtil::TestsetIsEQ(n.getLB(), c.getLB()) && n.get_level() < c.get_level() );
                    }),c);
                    elementCount++;
                }

                void add_nodes(std::vector<Node> &node_list, const unsigned &level) {
                    for(const Node& c : node_list) {
                        auto firstNonChild = llb.begin();
                        llb.insert(find_if(llb.begin(), firstNonChild, [&c](const Node &n) {
                            return NumericUtil::TestisGT(n.getLB(), c.getLB())
                                || (NumericUtil::TestsetIsEQ(n.getLB(), c.getLB()) && n.get_level() < c.get_level());
                        }), c);
                        elementCount++;
                    }
                }

                template<typename TimeT>
                Node find_next_node(TimeT &UB, unsigned long &phatomed_count, long processed_nodes, Node &best,
                        std::list<Node> &phatomed_nodes, unsigned &num_solutions, unsigned &num_improvements) {

                    while (! llb.empty()) {
                        Node b = *llb.begin();
                        llb.erase(llb.begin());
                        elementCount--;
                        std::streamsize ss = std::cout.precision();
                        if (b.get_level() >= num_levels) {
                            leafCount++;
                        }

                        if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - last).count() > 3) {
                            last = std::chrono::system_clock::now();
                            double t_spent = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start).count();
                            cout << "Proc=" << setw(8) << processed_nodes << " Leaf=" << setw(8) << leafCount
                                 << " Q.size=" << setw(5) << llb.size() << " LB=" << b.getLB() << " "
                                 << "UB=" << UB << " t=" << setprecision(2) << t_spent << "s" << std::setprecision(ss) << endl;
                            // cout << "\t" << b << endl;
                        }
                        if (NumericUtil::TestisGE(b.getLB(), UB)) {
                            phatomed_count++;
                            phatomed_nodes.insert(phatomed_nodes.begin(), b);
                            continue;
                        } else {
                            if (b.get_level() == num_levels) {
                                num_solutions++;
                                if(force_lb_update) {  // force update LB
                                    b.LB = b.instance->calculate_objective(b.getPartialSequence());
                                }
                                if(NumericUtil::TestisLT(b.getLB(), UB)) {
                                    UB = b.getLB();  best = b;
                                    cout << "*  " << b.getLB() << " " << endl;
                                    LOG(INFO) << "[BB Leaf Node] Improved solution! UB=" << UB;
                                    num_improvements++;
                                }
                                continue;
                            }
                            return b;
                        }
                    }
                    return Node(-1, std::vector<int>(), 2, 2, nullptr);
                }

                Node getTopNode() {
                    if (llb.size() > 0) {
                        Node node_smallest_lb = *min_element(llb.begin(), llb.end(), [](const Node& a, const Node& b) {
                            return NumericUtil::TestisLT(a.getLB(), b.getLB());
                        });
                        return node_smallest_lb;
                    }
                    return Node(-1, std::vector<int>(), 2, 2, nullptr);
                }

                unsigned long size() {
                    return elementCount;
                }

                // force lower bound update on leaf node, using problem objective function
                bool force_lb_update;
                long elementCount;
                long leafCount;
                // bb queue based on smallest lower bound among recently created nodes (DFS)
                //heap::pairing_heap<Node> bb_queue;
                std::list<Node> llb;
                // time measurement
                std::chrono::system_clock::time_point start, last;
                int num_levels;
            };

        } /* namespace bb */
    } //* namespace pfsp */
} /* namespace problem */

#endif //FLOWSHOP_SOLVER_DFSSEARCHSTRATEGY_H
