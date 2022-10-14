#ifndef FLOWSHOP_SOLVER_WCT_MODEL_STRATEGY_H
#define FLOWSHOP_SOLVER_WCT_MODEL_STRATEGY_H

#include <ilcplex/ilocplex.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

namespace robust {
            namespace hybrid {

                typedef IloArray<IloNumVarArray> NumVarMatrix;
                typedef IloArray<NumVarMatrix>   NumVar3Matrix;
                namespace ublas = boost::numeric::ublas;

                class PFSP_WCT_ModelStrategy {
                public:
                    virtual ~PFSP_WCT_ModelStrategy() = default;
                    
                    virtual NumVarMatrix populate_model(IloCplex &cplex, IloModel &model, IloRangeArray &c, IloNumVar &cmax, 
                                const int &m, const int &n,
                                const ublas::matrix<double> &P_bar,
                                const ublas::matrix<double> &P_hat, 
                                const ublas::vector<double> &w, 
                                const std::vector<int> &initial_permutation,
                                const std::vector< ublas::matrix<double> > &cut_list,
                                const bool &relax = false) = 0;

                    virtual void fix_partial_permutation_in_model(IloModel &model, NumVarMatrix &d, IloRangeArray &c,
                                                                const std::vector<int> &seq_1, const int &n) = 0;

                    virtual double calculate_bigM(int m, int n, const ublas::matrix<double> &P_bar,
                                                                const ublas::matrix<double> &P_hat,
                                                                const std::vector<int> &pi, const unsigned& bigM_type) = 0;

                    virtual ublas::matrix<double> calculate_cmax_given_scenario(int m, int n, const ublas::matrix<double> &P_bar,
                                                                            const ublas::matrix<double> &P_hat, const std::vector<int> &pi,
                                                                            const std::vector<std::vector<int>> &m_dev) = 0;

                    virtual void set_warm_start_solution(IloCplex &cplex, IloModel &model, const int &n, 
                                const std::vector<int>& perm, NumVarMatrix &d) = 0;

                    virtual std::vector<int> getSolutionAsPermutation(IloCplex &cplex, IloModel &model, const int &n, NumVarMatrix &d) = 0;
                };

            } /* namespace hybrid */
} /* namespace robust */

#endif //FLOWSHOP_SOLVER_WCT_MODEL_STRATEGY_H