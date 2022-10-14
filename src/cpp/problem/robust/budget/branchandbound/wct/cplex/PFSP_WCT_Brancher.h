//
// Created by mlevorato on 7/12/19.
//
/** \file nbranch.cpp Illustrate multi-way branching with CPLEX.
 * The code in this file illustrates how to perform multi-way branching in
 * CPLEX. CPLEX allows creation of at most two branches at a node. So in
 * order to create more than two branches we must emulate multi-way branching
 * by binary branching. To do that, instead of creating a branching like this
 * \verbatim
                N
        +---+---+---+---+
       C1  C2  C3  C4  C5
  \endverbatim
 * (create 5 children of node N) we create a branching like this
 * \verbatim
          N
        +-+-+
       C1   N
          +-+-+
         C2   N
            +-+-+
           C3   N
              +-+-+
             C4  C5
   \endverbatim
 * This branching produces exactly the same offsprings of node N but in
 * multiple levels. In this branching the left branch is always one of
 * the children to be created while the right branch is a duplicate of the
 * parent node N (expect the last level in which the right branch is also
 * one of the children to be created).
 *
 * To implement this strategy we proceed as follows:
 * Assume we are at a node N at which we want to create k > 2 children.
 * For each child to be created we specify the variables, directions and
 * bounds that define the child and store this information in an instance
 * of class BranchInfo. So one instance of class BranchInfo specifies one
 * child to be created. We put all those instances into a linked list and
 * store them as node data of N.
 * If we hit a node N that has node data stored then we know that we are
 * in the middle of emulating a multi-way branching. In that case we create
 * two branches:
 * 1. We extract the first BranchInfo instance from the linked list stored
 *    at N. From the information in that instance we create the left child.
 * 2. We create a right child that defines the same node as N (no additional
 *    bound changes). The node data of this child is the remaining list of
 *    BranchInfo instances.
 * In step 2, if there is exactly one BranchInfo left in the list then we
 * instead create the last child (which is defined by the remaining
 * BranchInfo instance).
 *
 * Source:
 * https://www.ibm.com/developerworks/community/forums/html/topic?id=40401984-2c29-4bea-88ce-1fd9f8b0cdc3
 */
#ifndef FLOWSHOP_SOLVER_PFSP_WCT_BRANCHER_H
#define FLOWSHOP_SOLVER_PFSP_WCT_BRANCHER_H

#include <iostream>
#include <sstream>
#include <vector>
#include <boost/foreach.hpp>
// #include <glog/logging.h>

#include <ilcplex/ilocplex.h>

#include "../../../../../../util/include/CollectionUtil.h"
#include "../../../../RobPFSInstance_WCT.h"
#include "../../../../../PFSSolution.h"
#include "PFSP_WCT_BranchInfo.h"
#include "PFSP_WCT_DFSNodeHandler.h"
#include "../../../../../PFSP_Parameters.h"
#include "../PFSP_WCT_Budget_LowerBound.h"
#include "../../../../../deterministic/PFSProblem.h"


namespace multichild {

    namespace ublas = boost::numeric::ublas;
    typedef IloArray<IloIntVarArray> IntVarMatrix;
    typedef IloArray<IloNumVarArray> NumVarMatrix;
    typedef IloArray<IloBoolVarArray> BoolVarMatrix;
    using namespace parameters;

/**
 * Branch callback that allows creation of multi-way branches.
 * Adaptation of the original Java multichild implementation by:
   Paul A. Rubin (http://about.me/paul.a.rubin)
 * https://orinanobworld.blogspot.com/2014/10/multiple-children-again.html
 * https://stackoverflow.com/questions/57244516/why-there-isnt-some-nodeid-in-nodecallback
 */
class PFSP_WCT_Brancher : public IloCplex::BranchCallbackI {
public:
    static IloCplex::MIPCallbackI::NodeId lastNodeId;
    std::vector< ublas::matrix<double> > P_time_scenario_list;
    std::vector< boost::multi_array<double, 3> > l_kl_scenario_list;
    
    /**
     * Constructor.
    */
    PFSP_WCT_Brancher(IloEnv env, const IloNumVar &cmax_v, NumVarMatrix v, const robust::RobPFSInstance_WCT &inst, 
                 bool verb, bool cbounds, 
                 const string &poutput_file, const bool &at_root, const problem::common::PFSSolution& initial_solution,
                 const double& root_lower_bound, const std::vector< ublas::matrix<double> > &p_P_time_scenario_list,
                 const std::vector< boost::multi_array<double, 3> > &p_l_kl_scenario_list,
                 const string &p_model_name, const bool &use_dominance) :
        IloCplex::BranchCallbackI(env), cmax_var(cmax_v), int_vars(v),
        numMachines(inst.m), numJobs(inst.n), atRoot(at_root), verbose(verb), instance(inst),
        use_combinatorial_bounds(cbounds), use_dominance_rules(use_dominance), 
        // UB(std::numeric_limits<double>::max()),
        output_file(poutput_file), ub_sol(initial_solution), rootLB(root_lower_bound),
        num_improvements(0), P_time_scenario_list(p_P_time_scenario_list), 
        l_kl_scenario_list(p_l_kl_scenario_list), model_name(p_model_name)
    {
            
    }

    /**
     * Function to duplicate this callback.
     * This function is required by the BranchCallbackI super class.
     */
    IloCplex::CallbackI *duplicateCallback() const {
        return new (getEnv()) PFSP_WCT_Brancher(getEnv(), cmax_var, int_vars, instance, verbose,
                use_combinatorial_bounds, output_file, atRoot, ub_sol, rootLB, P_time_scenario_list,
                l_kl_scenario_list, model_name, use_dominance_rules);

    }

    void check_for_CPLEX_bound_improvements(double &UB, double &LB) {
        double CPLEX_LB = getObjValue();  // Returns the obj value of the current solution of the LP relaxation
        double CPLEX_UB = getIncumbentObjValue(); // obj value of the current best integer solution found so far
        if (CPLEX_LB > LB) {
            // LOG(INFO)<<"*** CPLEX improved LB value! currentLB="<< LB << " vs CPLEX_LB=" << getObjValue() <<"\n";
            LB = CPLEX_LB;
        }
        if (CPLEX_UB < UB) {
            UB = CPLEX_UB;
            num_improvements++;
            // std::cout << "*** CPLEX improved UB value! UB=" << UB << "\n";
            // LOG(INFO) << "*** CPLEX improved UB value! UB=" << UB;
        }
    }

    long getNumberOfImprovements() {
        return num_improvements;
    }

    /**
     * Function that is invoked by CPLEX on each node.
     * It is responsible for creating new branches. If the function neither
     * explicitly prune()s the node nor explicitly creates at least one
     * branch then CPLEX will use the branches it would have created itself.
     * https://stackoverflow.com/questions/57244516/why-there-isnt-some-nodeid-in-nodecallback
     */
    void main() {
        // The objective function estimate for new nodes we create.
        double const estimate = getObjValue();

        // uncomment the next three lines for extra detail in the output
        //    std::cerr << "@@@ Entering the branch callback at node "
        //                       << getNodeId()._id << "\n";
        stringstream msg;
        msg << ">>> At node " << getNodeId()._id << ", LB = " << getObjValue() << ", UB = " << getIncumbentObjValue() << ", ";
        if (!use_combinatorial_bounds && !use_dominance_rules)
            return;
        // store last processed node id for use in NodeHandler (Node Selection Strategy)
        lastNodeId = getNodeId();
        IloCplex::NodeCallbackI::NodeId id;
        ublas::vector<long> lb_invocations(8, 0);
        robust::PFSP_WCT_Budget_LowerBound lower_bound;
        double UB = std::numeric_limits<double>::max();

        if (atRoot) {
            // create a new unscheduled job list
            std::vector<int> job_list = util::CollectionUtil::create_initial_job_list(numJobs);
            // create an attachment
            // https://www.ibm.com/developerworks/community/forums/html/topic?id=77777777-0000-0000-0000-000014856211
            PFSP_WCT_BranchInfo *root_info = new PFSP_WCT_BranchInfo(0, &instance, getNodeId(), job_list);
            // root_info->node.node_data.P_time_scenario_list = P_time_scenario_list;
            // root_info->node.node_data.l_kl_scenario_list = l_kl_scenario_list;
            double lb = lower_bound.lb_multi_scenario(instance, UB, root_info->node.seq_1, 
                            job_list, P_time_scenario_list, l_kl_scenario_list);
                            // root_info->node.node_data.P_time_scenario_list, 
                            // root_info->node.node_data.l_kl_scenario_list);
            root_info->setLB(lb);
            if(rootLB > 0.0) {
                root_info->setLB(max(root_info->getLB(), rootLB));
            }
            setNodeData(root_info);
            // NOTE: We CANNOT use ub_sol.value as UB, since this Obj Value does not consider existing MP cuts
            // Instead, we must ask for CPLEX best objective
            check_for_CPLEX_bound_improvements(UB, root_info->node.LB);
            atRoot = false;  // clear the root flag
            std::cout << "=> Initial Upper Bound is " << UB << "\n";
            std::cout << "=> Root Lower Bound is " << root_info->getLB() << "\n";

            PFSP_WCT_DFSNodeHandler::initialize_dfs(instance.getNumberOfJobs());
            // Note that we don't create the actual branches here.
            // Instead the code below will notice that info is now not
            // 0 and will start creating branches from that.
        }
        // Figure out if we are in a multi-way branching. The intermediate
        // nodes in a multi-way branching have node data attached to them that
        // describes the children that are yet to be created.
        PFSP_WCT_BranchInfo *currentNodeInfo = dynamic_cast<PFSP_WCT_BranchInfo *>(getNodeData());
        if ( !currentNodeInfo ) {
            // We don't have any node data. That is we need to decide what
            // to do: Either create the branches that CPLEX would create
            // check whether we are at the root node (which will have no attachment)
            cerr << "ERROR: No node info found, and we are not at root node !\n";
            LOG(ERROR) << "ERROR: No node info found, and we are not at root node !\n";
            abort();
        } else {
            try {
                // We have a node data object. That means that we are at an
                // intermediate node of a multi-level branch. We just create
                // the branch that is described in the node data object.
                currentNodeInfo->setProcessed();
                // TEST 2
                check_for_CPLEX_bound_improvements(UB, currentNodeInfo->node.LB);
                if (currentNodeInfo->node.LB >= UB || PFSP_WCT_DFSNodeHandler::prune_next_node) {
                    //cerr << " => pruned \n";
                    msg << " LB after CPLEX = " << currentNodeInfo->node.LB << ", ";
                    // cerr << msg.str() << " => CPLEX PRUNE! " << DFSNodeHandler::prune_next_node << "\n";
                    prune();
                    return;
                }
                //cerr << "\n";
                // If last tree level was reached, yield control to CPLEX
                if (currentNodeInfo->getNodeLevel() >= currentNodeInfo->getNumberOfJobs()) {
                    // TODO Calcular aqui o Upper Bound pelo Cmax e atualizar no CPLEX
                    cout << "Last B&B level reached !\n";
                    LOG(INFO) << "\n\n\n\n*** Last B&B level reached !\n\n\n\n\n\n";
                    if (currentNodeInfo->getLB() < UB) {
                        UB = currentNodeInfo->node.LB;
                        cout << "==>  " << currentNodeInfo->getLB() << " " << endl;
                        std::cout << "*** Improved solution! UB=" << UB;
                        LOG(INFO) << "*** Improved solution! UB=" << UB;
                    }
                    return;
                }
                int n = currentNodeInfo->getNumberOfJobs();
                int iterationCount = currentNodeInfo->getIterationCount();
                
                if (iterationCount >= currentNodeInfo->node.unscheduled_job_list.size()) {
                    cerr << "iterationCount >= J.size() !\n";
                }
                int currentJob = currentNodeInfo->node.unscheduled_job_list[iterationCount];
                if (currentJob > currentNodeInfo->getNumberOfJobs()) {
                    cerr << "currentJob > n !\n";
                    LOG(ERROR) << "currentJob > n !\n";
                }
                int currentLevel = currentNodeInfo->getNodeLevel();
                int children_count = 0;
                
                msg << "Level=" << currentNodeInfo->getNodeLevel() << ",Iteration="
                    << currentNodeInfo->getIterationCount()
                    << "/" << (currentNodeInfo->node.unscheduled_job_list.size() - 1)
                    << ", Unscheduled = " << currentNodeInfo->node.unscheduled_job_list.size();
                    // << ", Available sequence range is " << currentNodeInfo->node.seq_1.size() << " until "
                    // << (n - currentNodeInfo->node.seq_2.size() - 1);
                
                IloNum objValue = getObjValue();  // current node lower bound
                int nextPosition = currentNodeInfo->getNodeLevel();
                
                msg << ": Current job being positioned is " << currentJob << " on position " << nextPosition;
                if (currentLevel <= n) {
                    // CREATION OF THE FIRST NODE CHILD: First child will fix currentJob in the next available position
                    // Reset the iteration count, since remaining jobs will be fixed in a new sequence position
                    VarVector vars1;
                    BndVector bnds1;
                    DirVector dirs1;
                    cut_fix_job_in_position(vars1, bnds1, dirs1, currentJob - 1, nextPosition, 
                        currentNodeInfo->node.seq_1, currentNodeInfo->node.unscheduled_job_list);
                    PFSP_WCT_BranchInfo *nodeInfoChild1 = new PFSP_WCT_BranchInfo(0, &instance, getNodeId(),
                                                                        currentNodeInfo->node, currentJob);
                    nodeInfoChild1->node.seq_1.push_back(currentJob);
                    // Invoke Combinatorial Lower Bound calculations
                    bool prune_by_dominance = false;
                    if(use_dominance_rules) {
                        prune_by_dominance = lower_bound.dominance_multi_scenario(instance, currentNodeInfo->node.seq_1, 
                            currentNodeInfo->node.unscheduled_job_list, P_time_scenario_list, l_kl_scenario_list);
                    }
                    if (use_combinatorial_bounds) {
                        double lb = lower_bound.lb_multi_scenario(instance, UB, currentNodeInfo->node.seq_1, 
                            currentNodeInfo->node.unscheduled_job_list, P_time_scenario_list, l_kl_scenario_list);
                            // currentNodeInfo->node.node_data.P_time_scenario_list, 
                            // currentNodeInfo->node.node_data.l_kl_scenario_list);
                        nodeInfoChild1->node.LB = max(currentNodeInfo->node.LB, lb);
                    }
                    if (nodeInfoChild1->node.LB < UB && (!prune_by_dominance)) {
                        children_count++;
                        // job index used in the MIP model must start at zero!
                        ////// cut_lb_cmax(vars1, bnds1, dirs1, nodeInfoChild1->node.LB);   // FIXME Cut lb max
                        id = createBranch(vars1, bnds1, dirs1, nodeInfoChild1->node.LB, nodeInfoChild1);
                        nodeInfoChild1->setId(id);
                        currentNodeInfo->addChildren(id);
                        PFSP_WCT_DFSNodeHandler::dfs.add(
                                PFSP_WCT_SimpleNodeInfo(id, nodeInfoChild1->node.LB, nodeInfoChild1->node.get_level(), &instance));
                        // normal (non-leaf node), proceed to choose the next job to position or the next level
                        // next level of the BB tree will be processed
                        msg << "\n     -- adding one normal child (" << id._id << ") to fix job " << currentJob
                            << " on current tree level " << currentLevel;
                    } else {
                        //std::cerr << "PRUNE\n";
                        delete nodeInfoChild1;
                        msg << "\n     -- pruning normal child to fix job " << currentJob
                            << " on current tree level " << currentLevel;
                    }
                    vars1.end();
                    bnds1.end();
                    dirs1.end();
                    if (currentLevel == n) {  // Reached full solution (leaf) node, create a single child with exactly the job being fixed
                        // this case should only occur if the last tree level (n) has been reached
                        msg << "\n     -- (LEAF) adding SINGLE child (" << id._id << ") with exactly " << n
                            << " fixed jobs";
                        if (verbose) std::cout << msg.str();
                        // FIXME Trecho de debug ******************
                        // Calculate the objective function and validate it against CPLEX objective
                        // double wct_validation = PFSProblem::calculate_wct(pi, m, n, p_ij, w);
                        //// std::vector<int> pi;
                        //// getIncumbentObjValue();

                    }
                // Still positioning jobs in the current level ?
                } 
                if (iterationCount + 1 < currentNodeInfo->node.unscheduled_job_list.size()) {
                    children_count++;
                    // Second child clones the parent with adjusted node data that points to the next iteration on the current level
                    // job index used in the MIP model must start at zero!
                    VarVector vars2;
                    BndVector bnds2;
                    DirVector dirs2;
                    //cut_forbid_job_in_position(vars2, bnds2, dirs2, currentJob - 1, nextPosition);
                    PFSP_WCT_BranchInfo *nodeInfoChild2 = new PFSP_WCT_BranchInfo(iterationCount + 1, &instance, 
                                                                        getNodeId(), currentNodeInfo->node);
                    id = createBranch(vars2, bnds2, dirs2, nodeInfoChild2->node.LB, nodeInfoChild2);
                    nodeInfoChild2->setId(id);
                    currentNodeInfo->addChildren(id);
                    PFSP_WCT_DFSNodeHandler::add_to_metanodes(PFSP_WCT_SimpleNodeInfo(id, nodeInfoChild2->node.LB, nodeInfoChild2->node.get_level(), &instance));
                    msg << " and clone child (" << id._id << ") to process next iteration = " << (iterationCount + 1)
                        // << "/" << (currentNodeInfo->node.unscheduled_job_list.size() - 1)
                        << " in same tree level = " << currentLevel;
                    if (verbose) std::cout << msg.str() << "\n";
                    vars2.end();
                    bnds2.end();
                    dirs2.end();
                } else {  // Case B: Go to next level of the tree (unscheduled jobs will be put in next available position)
                    if (verbose)
                        std::cout << msg.str() << " (No Iteration Child, last iteration for the current level).\n";
                }
                // Update current node data with children info and processed flag
                setNodeData(currentNodeInfo);
                if(children_count == 0)
                    prune();
            } catch (IloException &e) {
                std::cerr << "IloException: " << e << std::endl;
                LOG(ERROR) << "IloException: " << e;
                throw;
            } catch (std::exception& e) {
                std::cerr << "Exception: " << e.what() << std::endl;
                LOG(ERROR) << "Exception: " << e.what();
                throw;
            } catch (...) {
                std::cerr << "UnknownException. " << std::endl;
                LOG(ERROR) << "UnknownException. ";
                // ex.printStackTrace();
                throw;
            }
        }
    }

private:
    /** 
     * Integer variables in the model.
     * This array is required to select the variables on which we branch.
     */
    IloNumVar cmax_var;
    NumVarMatrix int_vars;
    bool atRoot;  // are we at the root node?
    double rootLB;

    int numMachines;
    int numJobs;
    bool verbose;
    bool use_combinatorial_bounds;
    bool use_dominance_rules;
    robust::RobPFSInstance_WCT instance;
    string output_file;
    // Upper Bound
    double fixed_upper_bound;
    problem::common::PFSSolution best_sol, ub_sol;
    long num_improvements;
    // name of the underlying MILP model
    string model_name;

    /**
     * Create a branch.
     * The function translates the passed in STL containers into IloArray<>
     * instances and invokes BranchCallbackI::makeBranch() with those
     * arrays.
     */
    IloCplex::NodeCallbackI::NodeId createBranch(VarVector const &vars,
                      BndVector const &bnds,
                      DirVector const &dirs,
                      IloNum const &estimate,
                      PFSP_WCT_BranchInfo *info = 0) {
        // If we branch on a single variable then there is no need to
        // create IloArray<> instances since there is an explicit overload
        // of BranchCallbackI::makeBranch() for this special case.
        if ( vars.size() == 1 ) {
            return makeBranch(vars.front(), bnds.front(), dirs.front(), estimate, info);
        }
        else {
            IloNumVarArray v(getEnv());  // IloIntVarArray aVars;  // variables as array
            IloCplex::BranchDirectionArray d(getEnv());  // directions as array
            IloNumArray b(getEnv());  // bounds as doubles

            for (VarVector::size_type i = 0; i < vars.size(); ++i) {
                v.add(vars[i]);
                d.add(dirs[i]);
                b.add(bnds[i]);
            }
            IloCplex::NodeCallbackI::NodeId id = makeBranch(v, b, d, estimate, info);
            b.end();
            d.end();
            v.end();
            return id;
        }
    }

    bool is_dichotomy_model(const string &model_name) {
        return ((model_name == "hybrid-liao-you") || (model_name == "hybrid-manne"));
    }

    /**
    * Sets up the arguments for a cut that ensures job i
    * has been fixed in position j.
    *
    * https://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.3/ilog.odms.cplex.help/refdotnetcplex/html/T_ILOG_CPLEX_Cplex_BranchDirection.htm
    * @param job_to_fix the job being fixed
    * @param position_to_fix the target position in the sequence pi
    * @param partial_seq the partial job sequence (before fixing job 'job_to_fix')
    * @param unscheduled_job_list the list of unscheduled jobs, including 'job_to_fix'
    */
    void cut_fix_job_in_position(VarVector &vars,
                                 BndVector &bnds,
                                 DirVector &dirs, const int &job_to_fix, const int &position_to_fix,
                                 const std::vector<int> &partial_seq, const std::vector<int> &unscheduled_job_list) {
        if(is_dichotomy_model(model_name)) {
            // Cut generation for Manne and Liao-You models
            int pos_k = position_to_fix;
            int k = job_to_fix;
            // d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
            // Assign d[i][k] = 1 to all jobs i that come before job k in permutation seq_1
            for(int ii : partial_seq) {  // for each job i that comes before job k
                // job i is in position pos_i, and comes before job k
                int i = ii - 1;  // job number must be zero-based in the MIP model !
                // i = 1..n-1 , k = i+1..n
                if((i < numJobs - 1) && (k > i)) {
                    /// c.add(d[i][k] == 1);
                    vars.push_back(int_vars[i][k]);
                    bnds.push_back(1.0);
                    // Up: set the lower bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchUp);
                    // set upper and lower bounds of the variable to 1
                    vars.push_back(int_vars[i][k]);
                    bnds.push_back(1.0);
                    // BranchDirection Down: set the upper bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchDown);
                }
                if((k < numJobs - 1) && (i > k)) {
                    /// c.add(d[k][i] == 0);
                    vars.push_back(int_vars[k][i]);
                    bnds.push_back(0.0);
                    dirs.push_back(IloCplex::BranchDown);
                    vars.push_back(int_vars[k][i]);
                    bnds.push_back(0.0);
                    dirs.push_back(IloCplex::BranchUp);
                }
            }
            for(int ii : unscheduled_job_list) {  // for each job i that comes after job k
                int i = ii - 1;  // job number must be zero-based in the MIP model !
                if(i != job_to_fix) {  // job k is scheduled any time before job i
                    if((k < numJobs - 1) && (i > k)) {
                        /// c.add(d[k][i] == 1);
                        vars.push_back(int_vars[k][i]);
                        bnds.push_back(1.0);
                        dirs.push_back(IloCplex::BranchUp);
                        vars.push_back(int_vars[k][i]);
                        bnds.push_back(1.0);
                        dirs.push_back(IloCplex::BranchDown);
                    }
                    if((i < numJobs - 1) && (k > i)) {
                        /// c.add(d[i][k] == 0);
                        vars.push_back(int_vars[i][k]);
                        bnds.push_back(0.0);
                        dirs.push_back(IloCplex::BranchDown);
                        vars.push_back(int_vars[i][k]);
                        bnds.push_back(0.0);
                        dirs.push_back(IloCplex::BranchUp);
                    }
                }
            }
        } else {  // assignment model
            for (int job = 0; job < numJobs; job++) {  // for each job 'job'
                if (job == job_to_fix) {  // job 'job_to_fix' occupies position 'position_to_fix': force the position variable z[i][j] to be 1
                    vars.push_back(int_vars[job][position_to_fix]);
                    bnds.push_back(1.0);
                    // Up: set the lower bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchUp);
                    // set upper and lower bounds of the variable to 1
                    vars.push_back(int_vars[job][position_to_fix]);
                    bnds.push_back(1.0);
                    // BranchDirection Down: set the upper bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchDown);
                } else {  // no other job can occupy position 'position_to_fix': force variables z[i][j] to be 0
                    vars.push_back(int_vars[job][position_to_fix]);
                    bnds.push_back(0.0);
                    // BranchDirection Down: set the upper bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchDown);
                    vars.push_back(int_vars[job][position_to_fix]);
                    bnds.push_back(0.0);
                    dirs.push_back(IloCplex::BranchUp);
                }
            }
            for (int pos = 0; pos < numJobs; pos++) {  // for each position pos
                if (pos != position_to_fix) {  // job 'job_to_fix' cannot occupy any other position other than position_to_fix: force variables to be 0
                    vars.push_back(int_vars[job_to_fix][pos]);
                    bnds.push_back(0.0);
                    // BranchDirection Down: set the upper bound of the corresponding variable to the specified bound.
                    dirs.push_back(IloCplex::BranchDown);
                    vars.push_back(int_vars[job_to_fix][pos]);
                    bnds.push_back(0.0);
                    dirs.push_back(IloCplex::BranchUp);
                }
            }
        }
        // TODO Force the Objective Function (Makespan) to be greater than or equal to the Lower Bound (Cmax >= LB(Node))
    }

    /**
     * Sets up the arguments for a cut that forbids job i
     * of being fixed in position j.
     *
     * @param job_to_fix the job being fixed
     * @param position_to_fix the target position in the sequence pi
     */
    void cut_forbid_job_in_position(VarVector &vars,
                                    BndVector &bnds,
                                    DirVector &dirs, const int &job_to_fix, const int &position_to_fix) {
        if(is_dichotomy_model(model_name)) {
            // TODO Implement cuts for Manne and Liao-You models
        } else {  // assignment models
            // job 'job_to_fix' will NOT occupy position j in the sequence pi
            vars.push_back(int_vars[job_to_fix][position_to_fix]);
            bnds.push_back(0.0);
            // BranchDirection Down: set the upper bound of the corresponding variable to the specified bound.
            dirs.push_back(IloCplex::BranchDown);
        }
    }

    void cut_lb_cmax(VarVector &vars,
                     BndVector &bnds,
                     DirVector &dirs, const double &cmax_value) {

        vars.push_back(cmax_var);
        bnds.push_back(cmax_value);
        // BranchDirection Up: set the lower bound of the corresponding variable to the specified bound.
        dirs.push_back(IloCplex::BranchUp);
    }
};

}



#endif //FLOWSHOP_SOLVER_PFSP_WCT_BRANCHER_H
