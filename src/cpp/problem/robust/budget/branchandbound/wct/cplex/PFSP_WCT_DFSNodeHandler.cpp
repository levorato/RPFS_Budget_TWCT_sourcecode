//
// Created by mlevorato on 7/12/19.
//

#include "./PFSP_WCT_DFSNodeHandler.h"
#include <chrono>
#include <list>

namespace multichild {
    DFSSearchStrategy<PFSP_WCT_SimpleNodeInfo> PFSP_WCT_DFSNodeHandler::dfs = DFSSearchStrategy<PFSP_WCT_SimpleNodeInfo>(2, std::chrono::system_clock::now(), false);
    std::list<PFSP_WCT_SimpleNodeInfo> PFSP_WCT_DFSNodeHandler::metanodes = std::list<PFSP_WCT_SimpleNodeInfo>();
    std::list<PFSP_WCT_SimpleNodeInfo> PFSP_WCT_DFSNodeHandler::nodes_to_prune = std::list<PFSP_WCT_SimpleNodeInfo>();
    bool PFSP_WCT_DFSNodeHandler::prune_next_node = false;
}
