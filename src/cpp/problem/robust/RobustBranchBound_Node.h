#ifndef FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_NODE_H
#define FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_NODE_H

#include "../GlobalTypes.h"
#include "../../util/NumericUtil.h"
#include <string>
#include <sstream>
#include <vector>
#include <chrono>


namespace robust {

    using namespace std;
    using namespace parameters;
    namespace fs = boost::filesystem;
    using boost::timer::cpu_timer;

    class NodeStatistics {
    public:
        NodeStatistics() { }

        NodeStatistics(long queue_s, long phat_cnt, long lvl, double ubound, long iter, long time, double root_LB,
                       double gap, int max_level, int num_solutions) :
                queue_size(queue_s), phatomed_count(phat_cnt), level(lvl), ub(ubound), iteration(iter),
                time_spent(time), rootlowerbound(root_LB), gap(gap),
                max_level(max_level), n_solutions(num_solutions)
        {
        }

        // the size of the BB queue when this node was created
        long queue_size;
        // number of phatomed nodes when this node was created
        long phatomed_count;
        // the level of the node in the tree
        long level;
        // the best upper bound when node was created
        double ub;
        // the iteration number of the BB loop when node was created
        long iteration;
        // the time spent so far (in seconds) when node was created
        long time_spent;
        // the lower bound value of the root node in the BB tree
        double rootlowerbound;
        // the gap in solution value (UB - LB)
        double gap;
        // maximum level explored so far in the BB tree
        int max_level;
        // number of feasible primal solutions found so far
        int n_solutions;
    };

    template<class Instance, class NodeData>
    class BaseNode {
    public:
        long nodeNumber;
        int level;
        Time LB;
        // set of unscheduled jobs J
        std::vector<int> unscheduled_job_list;
        // partial sequence of n1 jobs scheduled on top of the global sequence (job indices)
        // All indices start at 1, except for the permutation vector
        std::vector<int> seq_1;
        // Implementation-specific node data, inherited from parent node
        NodeData node_data;
        // Pointer to problem instance
        Instance *instance;

        // Root node c'tor
        BaseNode(long uid, std::vector<int> J, int m, int n, Time lb_value, Instance *p_instance, NodeData p_node_data):
                                                            nodeNumber(uid),
                                                            LB(lb_value), 
                                                            unscheduled_job_list(J), seq_1(),
                                                            node_data(p_node_data),
                                                            instance(p_instance) {
        }

        // Special c'tor to be used by the BB queue strategy
        BaseNode(long uid, std::vector<int> J, int m, int n, Instance *p_instance, Time lb_value = -50000) :
                nodeNumber(uid),
                LB(lb_value), 
                unscheduled_job_list(J), seq_1(),
                node_data(), instance(p_instance) {
        }

        /**
         * Build a B & B node from its parent, the job scheduled in the branch is job j.
         */
        BaseNode(long uid, BaseNode &parent, int j, int m, int n) : nodeNumber(uid), unscheduled_job_list(), seq_1(parent.seq_1),
                                                            LB(parent.LB), 
                                                            node_data(parent.node_data), instance(parent.instance)
        {
            // Remove job j from unscheduled job list of the descendant
            for(int job : parent.unscheduled_job_list) {
                if (job != j) {
                    unscheduled_job_list.push_back(job);
                }
            }
        }

        unsigned long get_level() const {
            return seq_1.size();
        }

        long get_id() const {
            return nodeNumber;
        }

        Time getLB() const {
            return LB;
        }

        std::vector<int> getPartialSequence() {
            return seq_1;
        }

        // CBFS BB Priority Queue criteria used by M. Ritt.
        friend bool operator<(const BaseNode &a, const BaseNode &b) {
            // return a.lb > b.lb || (a.lb == b.lb && a.meanTime() > b.meanTime());
            // return a.lb > b.lb || (a.lb == b.lb && a.numFixed() < b.numFixed()
            return NumericUtil::TestisGT(a.LB, b.LB) || (NumericUtil::TestsetIsEQ(a.LB, b.LB) && a.level < b.level);
        }
    };
}

#endif //FLOWSHOP_SOLVER_ROBUSTBRANCHBOUND_NODE_H
