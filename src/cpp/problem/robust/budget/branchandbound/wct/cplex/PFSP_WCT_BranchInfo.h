//
// Created by mlevorato on 7/12/19.
//

#ifndef FLOWSHOP_SOLVER_PFSP_WCT_BRANCHINFO_H
#define FLOWSHOP_SOLVER_PFSP_WCT_BRANCHINFO_H

#include <ilcplex/ilocplex.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/foreach.hpp>

#include "../../../../RobPFSInstance_WCT.h"
#include "../PFSP_WCT_Budget_NodeData.h"
#include "PFSP_WCT_SimpleNodeInfo.h"
#include "../../../../RobustBranchBound_Node.h"

namespace multichild {

    namespace ublas = boost::numeric::ublas;
    // using namespace problem::pfsp::bb;

    // We use STL containers instead of subclasses of IloArray<> since
    // STL containers are more lightweight and easier to handle as long
    // as their size stays small.
    typedef std::vector<IloNumVar> VarVector;
    typedef std::vector<IloCplex::BranchDirection> DirVector;
    typedef std::vector<IloNum> BndVector;


/*
 * BranchInfo defines an object that can be attached to nodes and guides
 * branching decisions.
 * Description of one branch:
 * An instance of this class describes a single branch that should
 * be created. To keep the code short we assume that we only branch
 * on variables (and not on constraints).
 * Instances of this class are stored as user data in B&B nodes and
 * are used to emulate multi-way branching. Multiple instances of this
 * class are organized as singly linked lists, implemented by the prev
 * field.
 */
class PFSP_WCT_BranchInfo : public IloCplex::MIPCallbackI::NodeData {
    public:
        typedef robust::BaseNode<robust::RobPFSInstance_WCT, robust::PFSP_WCT_Budget_NodeData> PFSP_WCT_Node;
        PFSP_WCT_Node node;
        IloCplex::NodeCallbackI::NodeId parentNodeId;  // node id of parent node

        /**
         * Constructor for branch root node.
        */
        PFSP_WCT_BranchInfo(int iterCount, robust::RobPFSInstance_WCT *inst, 
                       IloCplex::NodeCallbackI::NodeId parentNodeId, std::vector<int> &J)
                : iterationCount(iterCount), parentNodeId(parentNodeId), processed(false), children(),
                  node(0, J, inst->m, inst->n, inst, 0)
        {
        }

        /**
         * Constructor for branch that branches on multiple variables (next multiple variable).
        */
        PFSP_WCT_BranchInfo(int iterCount, robust::RobPFSInstance_WCT *inst, 
                        IloCplex::NodeCallbackI::NodeId parentNodeId, PFSP_WCT_Node &parentNode, const int &fixed_job)
                : iterationCount(iterCount), parentNodeId(parentNodeId), processed(false), children(), 
                node(0, parentNode, fixed_job, inst->m, inst->n)
        {
        }

        /**
         * Constructor for branch that branches on multiple variables (copy of meta node == copy of previous node data).
        */
        PFSP_WCT_BranchInfo(int iterCount, robust::RobPFSInstance_WCT *inst, 
                       IloCplex::NodeCallbackI::NodeId parentNodeId, PFSP_WCT_Node &parentNode)
                : iterationCount(iterCount), parentNodeId(parentNodeId), processed(false), children(),
                  node(parentNode)
        {
        }
        /**
         * Destructor.
         * CPLEX will invoke this destructor whenever it deletes a node.
         */
        ~PFSP_WCT_BranchInfo() {
        }

        void setProcessed() {
            this->processed = true;
        }

        bool isProcessed() {
            return this->processed;
        }

        void addChildren(IloCplex::NodeCallbackI::NodeId childrenNodeId) {
            children.push_back(childrenNodeId);
        }

        std::vector<IloCplex::NodeCallbackI::NodeId> getChildren() {
            return children;
        }

        int getIterationCount() const {
            return iterationCount;
        }

        /**
         * Get the number of jobs.
         * @return number of jobs
         */
        int getNumberOfJobs() const {
            return node.instance->n;
        }

        int getNodeLevel() const {
            return node.get_level();
        }

        robust::RobPFSInstance_WCT* getProblemInstance() {
            return node.instance;
        }

        std::vector<int> getPartialSequence() {
            return node.getPartialSequence();
        }

        Time getLB() const {
            return node.getLB();
        }

        void setLB(const Time &LB) {
            node.LB = LB;
        }

        /**
         * Get the node ID from the data.
         * @return the node ID
         */
        IloCplex::NodeCallbackI::NodeId getId() const {
            return id;
        }

        /**
         * Set the node ID.
         * @param n the node ID to use
         */
        void setId(const IloCplex::NodeCallbackI::NodeId &n) {
            this->id = n;
        }

    private:
        // the number of iterations of the job fixing loop (in current level)
        int iterationCount;
        IloCplex::NodeCallbackI::NodeId id;     // ID of the node
        std::vector<IloCplex::NodeCallbackI::NodeId> children;  // node id of children nodes (2)
        bool processed;
    };

}
#endif //FLOWSHOP_SOLVER_PFSPBRANCHINFO_H
