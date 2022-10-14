//
// Created by mlevorato on 8/14/19.
//

#ifndef FLOWSHOP_SOLVER_MANNE_PFSP_WCT_MODEL_H
#define FLOWSHOP_SOLVER_MANNE_PFSP_WCT_MODEL_H

#include <ilcplex/ilocplex.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

#include "PFSP_WCT_ModelStrategy.h"

namespace robust {
            namespace hybrid {

                namespace ublas = boost::numeric::ublas;

/**
 * This class consists of a model strategy (strategy design pattern).
 */
class MannePFSP_WCT_Model : public PFSP_WCT_ModelStrategy {
public:
    MannePFSP_WCT_Model() {  }
    ~MannePFSP_WCT_Model()  {  }

    NumVarMatrix populate_model(IloCplex &cplex, IloModel &model, IloRangeArray &c, IloNumVar &obj, 
                                            const int &m, const int &n,
                                            const ublas::matrix<double> &P_bar,
                                            const ublas::matrix<double> &P_hat, 
                                            const ublas::vector<double> &w, 
                                            const std::vector<int> &initial_permutation,
                                            const std::vector< ublas::matrix<double> > &cut_list,
                                            const bool &relax = false);

    void fix_partial_permutation_in_model(IloModel &model, NumVarMatrix &d, IloRangeArray &c,
                                                 const std::vector<int> &seq_1, const int &n);

    double calculate_bigM(int m, int n, const ublas::matrix<double> &P_bar,
                                                const ublas::matrix<double> &P_hat,
                                                const std::vector<int> &pi, const unsigned& bigM_type);

    ublas::matrix<double> calculate_cmax_given_scenario(int m, int n, const ublas::matrix<double> &P_bar,
                                                            const ublas::matrix<double> &P_hat, const std::vector<int> &pi,
                                                            const std::vector<std::vector<int>> &m_dev);

    void set_warm_start_solution(IloCplex &cplex, IloModel &model, const int &n, 
                const std::vector<int>& perm, NumVarMatrix &d);

    std::vector<int> getSolutionAsPermutation(IloCplex &cplex, IloModel &model, const int &n, NumVarMatrix &d);

    std::vector<int> topological_sort_by_dfs(const int &n, const ublas::matrix<int> &d_matrix);

    std::vector<int> topological_sort_by_dfs_v2(const int &n, const ublas::matrix<int> &d_matrix);
};

            } /* namespace hybrid */
} /* namespace robust */

#endif //FLOWSHOP_SOLVER_MANNE_PFSP_WCT_MODEL_H
