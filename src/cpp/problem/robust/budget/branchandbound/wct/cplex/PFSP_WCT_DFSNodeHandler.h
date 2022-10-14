//
// Created by mlevorato on 7/12/19.
//

#ifndef FLOWSHOP_SOLVER_WCT_DFSNODEHANDLER_H
#define FLOWSHOP_SOLVER_WCT_DFSNODEHANDLER_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <glog/logging.h>

#include "PFSP_WCT_BranchInfo.h"
#include "PFSP_WCT_SimpleNodeInfo.h"
#include "../../../../../deterministic/branchandbound/DFSSearchStrategy.h"
#include "../../../../../GlobalTypes.h"


namespace multichild {

    namespace ublas = boost::numeric::ublas;
    using namespace std;
    using namespace problem::pfsp::bb;

    class PFSP_WCT_DFSNodeHandler : public IloCplex::NodeCallbackI {
    protected:
        int numLevels;
        bool verbose;
        NodeSelectionStrategy node_sel_strategy;
        double* UB;
        bool first_exec;
    public:
        static DFSSearchStrategy<PFSP_WCT_SimpleNodeInfo> dfs;
        
        static std::list<PFSP_WCT_SimpleNodeInfo> metanodes;
        static std::list<PFSP_WCT_SimpleNodeInfo> nodes_to_prune;
        static bool prune_next_node;

        PFSP_WCT_DFSNodeHandler(IloEnv env, int num_levels, NodeSelectionStrategy nsel_strategy, double *ub, bool verb) :
                IloCplex::NodeCallbackI(env), node_sel_strategy(nsel_strategy), UB(ub),
                numLevels(num_levels), verbose(verb), first_exec(true) {

        }

        /**
        * Function to duplicate this callback.
        * This function is required by the BranchCallbackI super class.
        */
        IloCplex::CallbackI *duplicateCallback() const {
            return new (getEnv()) PFSP_WCT_DFSNodeHandler(getEnv(), numLevels, node_sel_strategy, UB, verbose);
        }

        static void initialize_dfs(int num_levels) {
            std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
            PFSP_WCT_DFSNodeHandler::dfs = DFSSearchStrategy<PFSP_WCT_SimpleNodeInfo>(num_levels, start_time, false);
        }

        static void add_to_metanodes(const PFSP_WCT_SimpleNodeInfo &c) {
            PFSP_WCT_DFSNodeHandler::metanodes.insert(PFSP_WCT_DFSNodeHandler::metanodes.begin(), c);
        }

        static void add_to_prune(const PFSP_WCT_SimpleNodeInfo &c) {
            PFSP_WCT_DFSNodeHandler::nodes_to_prune.insert(PFSP_WCT_DFSNodeHandler::nodes_to_prune.begin(), c);
        }

        void main() {
            long remainingNodes = getNremainingNodes64();
            /*
            std::cout << "*** DFSNodeHandler ::: DFS Queue size=" << DFSNodeHandler::dfs.size() << " ; Metanode size="
                    << DFSNodeHandler::metanode_count << "; TOTAL = "
                    << (DFSNodeHandler::dfs.size() + DFSNodeHandler::metanode_count) << "\n";
            std::cout << "*** DFSNodeHandler ::: remainingNodes = " << remainingNodes << "\n";
             */
            if(first_exec) {
                if(node_sel_strategy == NodeSelectionStrategy::cplex) {
                    LOG(INFO) << "[CPLEXNodeHandler] Node selection strategy is CPLEX default.\n";
                } else if (node_sel_strategy==NodeSelectionStrategy::dfs) {
                    LOG(INFO) << "[CPLEXNodeHandler] Node selection strategy is DFS (Ladhari default).\n";
                }
                first_exec = false;
            }
            unsigned long phatomed_count;
            unsigned num_solutions = 1, num_improvements = 0;  // TODO FIXME make these variables static
            PFSP_WCT_SimpleNodeInfo best;
//            if(prune_nodes_if_exist())
//                return;

            
            if(node_sel_strategy == NodeSelectionStrategy::cplex) {
                // let CPLEX decide the next node to process
                return;
            } 
            if(process_metanodes_if_exist()) {
                //cerr << "Process metanode : ";
                return;         
            }
            // default : normal next node selection strategy
            // First, try to expand any normal nodes awaiting in the queue
            if (PFSP_WCT_DFSNodeHandler::dfs.size() > 0) {  // If no more metanodes are available for this level, expand the next normal node
                //assert(DFSNodeHandler::dfs.size() == remainingNodes);
                PFSP_WCT_SimpleNodeInfo which_node = PFSP_WCT_DFSNodeHandler::dfs.find_next_node(*UB, phatomed_count,
                                                                                getNnodes64(), best,
                                                                                PFSP_WCT_DFSNodeHandler::nodes_to_prune,
                                                                                num_solutions, num_improvements); // PFSPBrancher::UB
                PFSP_WCT_DFSNodeHandler::prune_next_node = false;
                if (which_node.nodeNumber._id > 0) {
                    selectNode(which_node.nodeNumber);
                    // std::cout << "*** DFSNodeHandler ::: Next Node Id is " << which_node.nodeNumber._id << " - normal node\n";
                    return;
                }
            }
//              if(prune_nodes_if_exist())
//                  return;
            if(process_metanodes_if_exist())
                return;
            /// LOG(ERROR) << "WARN 2: NO MORE NODES TO SELECT ! remainingNodes = " << remainingNodes << ";";
            /// LOG(ERROR) << "[DFS] SOLUTION VALUE = " << *UB << ";";
            ////// abort();
        }

        bool prune_nodes_if_exist() {
            if(!PFSP_WCT_DFSNodeHandler::nodes_to_prune.empty()) {
                PFSP_WCT_SimpleNodeInfo which_node = *PFSP_WCT_DFSNodeHandler::nodes_to_prune.begin();
                PFSP_WCT_DFSNodeHandler::nodes_to_prune.erase(PFSP_WCT_DFSNodeHandler::nodes_to_prune.begin());
                PFSP_WCT_DFSNodeHandler::prune_next_node = true;
                selectNode(which_node.nodeNumber);
                return true;
            }
            return false;
        }

        bool process_metanodes_if_exist() {
            if (!PFSP_WCT_DFSNodeHandler::metanodes.empty()) {
                PFSP_WCT_SimpleNodeInfo which_node = *PFSP_WCT_DFSNodeHandler::metanodes.begin();
                PFSP_WCT_DFSNodeHandler::metanodes.erase(PFSP_WCT_DFSNodeHandler::metanodes.begin());
                PFSP_WCT_DFSNodeHandler::prune_next_node = false;
                selectNode(which_node.nodeNumber);
                // std::cout << "*** DFSNodeHandler ::: Next Node Id is " << which_node.nodeNumber._id << " - meta node\n";
                return true;
            }
            return false;
        }
    };

}

#endif //FLOWSHOP_SOLVER_DFSNODEHANDLER_H
