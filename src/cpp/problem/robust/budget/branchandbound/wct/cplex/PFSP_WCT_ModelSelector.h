#ifndef FLOWSHOP_SOLVER_PFSP_WCT_MODELSELECTOR_H
#define FLOWSHOP_SOLVER_PFSP_WCT_MODELSELECTOR_H

#include <ilcplex/ilocplex.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <stdexcept>
#include <memory>

#include "WilsonPFSP_WCT_Model.h"
#include "LiaoYouPFSP_WCT_Model.h"
#include "MannePFSP_WCT_Model.h"
#include "TS2PFSP_WCT_Model.h"
#include "TS3PFSP_WCT_Model.h"
#include "TBAPFSP_WCT_Model.h"
#include "WST2PFSP_WCT_Model.h"
#include "PFSP_WCT_ModelStrategy.h"

namespace robust {
    namespace hybrid {

                typedef IloArray<IloNumVarArray> NumVarMatrix;
                namespace ublas = boost::numeric::ublas;

/**
 * This class represents the Context class of the strategy design pattern.
 */
class PFSP_WCT_ModelSelector {
public:
    PFSP_WCT_ModelSelector(const string &model_name) {
        set_model(model_name);
    }

    void set_model(const string &model_name) {
        if(model_name == "hybrid-wilson") {
            _strategy = std::make_unique<WilsonPFSP_WCT_Model>();
        } else if(model_name == "hybrid-liao-you") {
            _strategy = std::make_unique<LiaoYouPFSP_WCT_Model>();
        } else if(model_name == "hybrid-manne") {
            _strategy = std::make_unique<MannePFSP_WCT_Model>();
        } else if(model_name == "hybrid-ts2") {
            _strategy = std::make_unique<TS2PFSP_WCT_Model>();
        } else if(model_name == "hybrid-ts3") {
            _strategy = std::make_unique<TS3PFSP_WCT_Model>();
        } else if(model_name == "hybrid-tba") {
            _strategy = std::make_unique<TBAPFSP_WCT_Model>();
        } else if(model_name == "hybrid-wst2") {
            _strategy = std::make_unique<WST2PFSP_WCT_Model>();
        } else {
            throw std::runtime_error("Error! MP model not found: " + model_name);
        }
    }

    NumVarMatrix populate_model(IloCplex &cplex, IloModel &model, IloRangeArray &c, IloNumVar &cmax,
                                const int &m, const int &n,
                                const ublas::matrix<double> &P_bar,
                                const ublas::matrix<double> &P_hat,
                                const ublas::vector<double> &w,
                                const std::vector<int> &initial_permutation,
                                const std::vector< ublas::matrix<double> > &cut_list,
                                const bool &relax = false) {
        return _strategy->populate_model(cplex, model, c, cmax, m, n, P_bar, P_hat, w, initial_permutation,
            cut_list, relax);
    }

    void set_warm_start_solution(IloCplex &cplex, IloModel &model, const int &n,
                const std::vector<int>& perm, NumVarMatrix &bin_vars) {
        _strategy->set_warm_start_solution(cplex, model, n, perm, bin_vars);
    }

    std::vector<int> getSolutionAsPermutation(IloCplex &cplex, IloModel &model, const int &n,
            NumVarMatrix &bin_vars) {
        return _strategy->getSolutionAsPermutation(cplex, model, n, bin_vars);
    }

private:
    std::unique_ptr<PFSP_WCT_ModelStrategy> _strategy;
};

            } /* namespace hybrid */
} /* namespace robust */

#endif //FLOWSHOP_SOLVER_PFSP_WCT_MODELSELECTOR_H
