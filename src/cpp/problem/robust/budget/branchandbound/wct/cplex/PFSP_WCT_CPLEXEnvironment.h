//
// Created by mlevorato on 7/26/19.
//

#ifndef FLOWSHOP_SOLVER_WCT_CPLEXENVIRONMENT_H
#define FLOWSHOP_SOLVER_WCT_CPLEXENVIRONMENT_H

#include "../../../../../../util/singleton.h"
#include <ilcplex/ilocplex.h>
#include <set>

/*
static void deepCopyIloModel (IloModel & oldModel_arg, IloModel & newModel ) {
    if ( oldModel_arg.getEnv().getImpl() == newModel.getEnv().getImpl() ) {
        std::cerr << "Models must be in different environments!" << std::endl;
        ::abort();
    }

    // Collect all IloExtractables that exist prior to cloning.
    std::set<IloExtractable,IloExtractableLess> exists;
    for (IloIterator<IloExtractable> it(oldModel_arg.getEnv()); it.ok(); ++it)
        exists.insert(*it);

    // Clone the model.
    IloEnv env = newModel.getEnv();
    newModel.end();
    newModel= IloGetClone(env, oldModel_arg);

    // Check that no new elements where created in the env of the original
    // model.
    for (IloIterator<IloExtractable> it(oldModel_arg.getEnv()); it.ok(); ++it) {
        if ( exists.find(*it) == exists.end() ) {
            std::cerr << "Extractable " << (*it).getId() << ", "
                      << (*it).getName() << ", "
                      << *it << " unexpectedly created"
                      << std::endl;
            ::abort();
        }
    }
}

class cpxnode {
public:

    // IloCplex should not be a
    // member of the object
    IloEnv env;
    IloModel model;

    //=======Constructor========
    cpxnode(): env(), model(env) {}

    //=====Copy constructor=======
    cpxnode(const cpxnode & original): env(), model(env){
        deepCopyIloModel(original.model, model);
    }
    //====Assignment operator=======
    cpxnode &operator=(const cpxnode &original) {
        env.end();
        env = IloEnv();
        model = IloModel(env);
        // Now the model is completely empty, and
        // we do not lose time, trying to remove
        // one-by-one its extractables
        deepCopyIloModel(original.model, model);
        return *this;
    }

    //Destructor
    ~cpxnode() {
        env.end();
    }
};
*/
/**
 * https://www.ibm.com/developerworks/community/forums/html/topic?id=a3d5edae-bfab-4cd5-bd8a-7bc90aa23f22&ps=25
 */
class PFSP_WCT_CPLEXEnvironment {
public:
    //cpxnode node;
    //int count;

    PFSP_WCT_CPLEXEnvironment() {  };  //: node(), count(0) {  };
    ~PFSP_WCT_CPLEXEnvironment() {  };

};

typedef Singleton<PFSP_WCT_CPLEXEnvironment> PFSP_WCT_CPLEXSingleton;   // Global declaration


#endif //FLOWSHOP_SOLVER_WCT_CPLEXENVIRONMENT_H
