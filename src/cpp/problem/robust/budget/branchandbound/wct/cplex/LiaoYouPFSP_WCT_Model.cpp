#include "LiaoYouPFSP_WCT_Model.h"
#include "../../../../../deterministic/PFSProblem.h"

#include <glog/logging.h>
#include <vector>
#include <stack>
#include <queue>
#include <algorithm>
#include <iostream> // std::cout
#include <utility> // std::pair
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>

namespace robust {
    namespace hybrid {
        using namespace std;
        using namespace boost;
        using namespace boost::graph;
        // boost::numeric::ublas::matrix

        /*
        * PFSP Liao-You Model
        * The deterministic MILP model will look like the following, where:
        *   d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
        *   s[r][i] = start (begin) time of job i on machine r
        *   q[r, i, k] = surplus variables related to machine r, job i and job k
        *   obj = objective function variable
        *
        * minimize obj
        * s.t.
        *   sum(w[i] * (S[m, i] + T[m, i]) for i in Jobs) <= obj
        *
        *   S[r, i] + T[r, i] <= S[r+1, i]                                    for r in 1..m-1, i in 1..n
        *   S[r, i] - S[r, k] + P * D[i, k] - T[r, k] == q[r, i, k]           for r in 1..m, i in 1..n-1, k in i+1..n
        *   P - T[r, i] - T[r, k] >= q[r, i, k]                               for r in 1..m, i in 1..n-1, k in i+1..n
        */
        NumVarMatrix LiaoYouPFSP_WCT_Model::populate_model(IloCplex &cplex, IloModel &model, IloRangeArray &c, 
                                                IloNumVar &obj,
                                                const int &m, const int &n, const ublas::matrix<double> &P_bar,
                                                const ublas::matrix<double> &P_hat, 
                                                const ublas::vector<double> &w, 
                                                const std::vector<int> &initial_permutation, 
                                                const std::vector< ublas::matrix<double> > &cut_list,
                                                const bool &relax) {
            IloEnv env = model.getEnv();
            std::vector<int> pi(n);
            std::iota (std::begin(pi), std::end(pi), 1); // Fill with 0, 1, ..., n.
            // IMPORTANT: FIX BIGM TYPE TO 1, SINCE TYPE 2 IS NOT GIVING CORRECT RESULTS FOR SOME INSTANCES
            double bigM = calculate_bigM(m, n, P_bar, P_hat, pi, 1);
            // 1. First stage variables and constraints
            NumVarMatrix d(env, n);  // d is bin => here and now variable, first stage
            for (IloInt i = 0; i < n; i++) {
                if(relax) {
                    d[i] = IloNumVarArray(env, n, 0, 1);  // Float
                } else {
                    d[i] = IloNumVarArray(env, n, 0, 1, IloNumVar::Bool);  // d is bin
                }
                std::stringstream sd;
                sd << "d(" << i << ")";
                d[i].setNames(sd.str().c_str());
            }
            // 2. Second stage variables and constraints
            for(int cut_count = 0; cut_count < cut_list.size(); cut_count++) {
                ublas::matrix<double> sep = cut_list[cut_count];  // the matrix containing current cut
                LOG(INFO) << "Processing cut: " << sep;
                //cout << "Processing cut: " << sep << "\n";
                // S[r, i] >= 0 : start time of job i on machine r => [M, N]
                NumVarMatrix S(env, m);  
                // q[r, i, k] >= 0 : Surplus variables => [M, N, N]
                NumVar3Matrix q(env, m);
                for (IloInt r = 0; r < m; r++) {
                    S[r] = IloNumVarArray(env, n, 0, IloInfinity);  // S[r, i] >= 0
                    std::stringstream ss;
                    ss << "S_" << cut_count << "(" << r << ")";
                    S[r].setNames(ss.str().c_str());
                    // https://www.ibm.com/support/pages/sample-create-and-use-multi-dimensional-ilonumvararray
                    q[r] = NumVarMatrix(env, n);
                    for (IloInt i = 0; i < n; i++) {
                        q[r][i] = IloNumVarArray(env, n, 0, IloInfinity);  // q[r, i, k] >= 0
                        std::stringstream sq;
                        sq << "q_" << cut_count << "(" << r << ", " << i << ")";
                        q[r][i].setNames(sq.str().c_str());
                    }
                }
                // add the constraints
                // (1) new_S[r, i] + P_bar[r, i] + (P_hat[r, i] * sep[r, i]) <= new_S[r+1, i]   for r in 1..m-1, i in 1..n
                for (IloInt r = 0; r < m - 1; r++) {
                    for (IloInt i = 0; i < n; i++) {
                        IloExpr expr1(env);  // (1)
                        expr1 += S[r][i] + P_bar(r+1, i+1) + (P_hat(r+1, i+1) * sep(r+1, i+1));
                        model.add(expr1 <= S[r+1][i]);
                        expr1.end();
                    }
                }
                // (2) new_S[r, i] - new_S[r, k] + P * D[i, k] - (P_bar[r, k] + (P_hat[r, k] * sep[r, k])) == new_q[r, i, k],  for r in 1..m, i in 1..n-1, k=i+1..n
                // (3) P - (P_bar[r, i] + P_hat[r, i] * sep[r, i]) - (P_bar[r, k] + (P_hat[r, k] * sep[r, k])) >= new_q[r, i, k],  for r in 1..m, i in 1..n-1, k=i+1..n
                for (IloInt r = 0; r < m; r++) {
                    for (IloInt i = 0; i < n - 1; i++) {
                        for (IloInt k = i + 1; k < n; k++) {
                            IloExpr expr1(env);  // (2)
                            expr1 += S[r][i] - S[r][k] + bigM * d[i][k] - (P_bar(r+1, k+1) + (P_hat(r+1, k+1) * sep(r+1, k+1)));
                            model.add(expr1 == q[r][i][k]);
                            expr1.end();
                            IloExpr expr2(env);  // (3)
                            expr2 += bigM - (P_bar(r+1, i+1) + P_hat(r+1, i+1) * sep(r+1, i+1)) - (P_bar(r+1, k+1) + (P_hat(r+1, k+1) * sep(r+1, k+1)));
                            model.add(expr2 >= q[r][i][k]);
                            expr2.end();
                        }
                    }
                }
                // add the objective function
                // (4) obj >= sum(w[i] * (new_S[m, i] + P_bar[m, i] + (P_hat[m, i] * sep[m, i])) for i in Jobs)
                IloExpr expr0(env);
                for (IloInt i = 0; i < n; i++) {
                    // w[i] * (new_S[m, i] + P_bar[m, i] + (P_hat[m, i]*z_sep[m, i]))
                    expr0 += w(i + 1) * (S[m - 1][i] + P_bar(m, i+1) + P_hat(m, i+1) * sep(m, i+1));
                }
                model.add(obj >= expr0);
                expr0.end();
            }
            model.add(IloMinimize(env, obj));
            // set mip start using the given permutation
            set_warm_start_solution(cplex, model, n, initial_permutation, d);
            return d;
        }

        void LiaoYouPFSP_WCT_Model::fix_partial_permutation_in_model (IloModel &model, NumVarMatrix &d, IloRangeArray &c,
                                                                const std::vector<int> &seq_1, 
                                                                const int &n) {
            IloEnv env = model.getEnv();
            for(int pos_k = 0; pos_k < seq_1.size(); pos_k++) {  // for all jobs k in position pos_k in seq_1
                int k = seq_1[pos_k] - 1;  // job number must be zero-based in the MIP model !
                // d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
                // Assign d[i][k] = 1 to all jobs i that come before job k in permutation seq_1
                for (int pos_i = 0; pos_i < pos_k; pos_i++) {  // for each position pos_i < pos_k
                    // job i is in position pos_i, and comes before job k
                    int i = seq_1[pos_i] - 1;  // job number must be zero-based in the MIP model !
                    c.add(d[i][k] == 1);
                    c.add(d[k][i] == 0);
                }
            }
            model.add(c);
        }

        /**
         * Calculates the completion time matrix C, given a problem instance (m x n), processing time matrices P_bar and
         * P_hat, a permutation pi and budget scenario m_dev.
         * @param m number of machines
         * @param n number of jobs
         * @param P_bar matrix of nominal processing times
         * @param P_hat matrix of processing time variations
         * @param pi the job permutation (job number starts at 1)
         * @param m_dev processing time scenario (i.e., for each machine r, which jobs will have their proc. time deviated )
         * @return completion time matrix C (zero-based)
         */
        ublas::matrix<double> LiaoYouPFSP_WCT_Model::calculate_cmax_given_scenario(int m, int n, const ublas::matrix<double> &P_bar,
                                                            const ublas::matrix<double> &P_hat, const std::vector<int> &pi,
                                                            const std::vector<std::vector<int>> &m_dev) {
            ublas::matrix<double> P_sum(P_bar);
            for (int i = 0; i < m; i++) {
                for (int j : m_dev[i]) {
                    if (j > 0) {
                        P_sum(i + 1, j) += P_hat(i + 1, j);
                    }
                }
            }
            std::vector<int> perm;
            for(int x = 0; x < pi.size(); x++) {
                perm.push_back(pi[x]);
            }
            ublas::matrix<double> C = PFSProblem::makespan(perm, m, P_sum);
            return C;
        }

        /**
         * Calculates the BigM value to be used in one of the Rob-PFSP-WCT MILP models, according to the type:
         * - type == 0 : standard upper bound for BigM (naive);
         * - type == 2 : Calculate the Makespan (Cmax) for the worst-case scenario (on all machines r, all jobs j have their
         *                  processing times deviated to their maximum value : P_bar(r,j) + P_hat(r,j)).
         * @param m number of machines
         * @param n number of jobs
         * @param P_bar matrix of nominal processing times
         * @param P_hat matrix of processing time variations
         * @return bigM value.
         */
        double LiaoYouPFSP_WCT_Model::calculate_bigM(int m, int n, const ublas::matrix<double> &P_bar,
                                                const ublas::matrix<double> &P_hat,
                                                const std::vector<int> &pi, const unsigned& bigM_type = 1) {
            std::vector<std::vector<int>> worst_Gamma_scenario(m, std::vector<int>());
            std::vector<std::vector<int>> best_Gamma_scenario(m, std::vector<int>());
            for (int r = 0; r < m; r++) {
                for (int i = 0; i < n; i++) {
                    worst_Gamma_scenario[r].push_back(i + 1);
                }
            }
            ublas::matrix<double> max_C = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario);
            ublas::matrix<double> min_C = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, best_Gamma_scenario);
            double max_value = calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario)(m,n) + 1;
            // std::cout << "max_value = " << max_value << "\n";
            // std::cout << "bigM_type = " << bigM_type << "\n";
            ublas::matrix<double> bigM = ublas::matrix<double>(m, n, max_value);
            if(bigM_type <= 1) {
                std::priority_queue<double> q;
                max_value = 0.0;
                for (int r = 0; r < m; r++) {
                    for (int i = 0; i < n; i++) {
                        // max_value += P_bar(r + 1, i + 1) + P_hat(r + 1, i + 1);
                        q.push(P_bar(r + 1, i + 1) + P_hat(r + 1, i + 1));
                    }
                }
                int k = m + n - 1; // number of indices we need
                for (int i = 0; i < k; ++i) {
                    max_value += q.top();
                    q.pop();
                }                
                return max_value + 1;
            } else {  // if(bigM_type == 2) {
                return calculate_cmax_given_scenario(m, n, P_bar, P_hat, pi, worst_Gamma_scenario)(m,n) + 1;
            }
        }

        /**
         * Set MIP start solution using the Upper Bound Solution given by permutation perm.
         */
        void LiaoYouPFSP_WCT_Model::set_warm_start_solution(IloCplex &cplex, IloModel &model, const int &n, 
                const std::vector<int>& perm, NumVarMatrix &d) {
            IloEnv env = model.getEnv();
            IloNumVarArray startVar(env);
            IloNumArray startVal(env);
            for(int pos_k = 0; pos_k < perm.size(); pos_k++) {  // for all jobs k in position pos_k in perm
                int k = perm[pos_k] - 1;  // job number must be zero-based in the MIP model !
                // d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
                // Assign d[i][k] = 1 to all jobs i that come before job k in permutation seq_1
                for (int pos_i = 0; pos_i < pos_k; pos_i++) {  // for each position pos_i < pos_k
                    // job i is in position pos_i, and comes before job k
                    int i = perm[pos_i] - 1;  // job number must be zero-based in the MIP model !
                    // i = 1..n-1 , k = i+1..n
                    /// c.add(d[i][k] == 1);
                    if((i < n - 1) && (k > i)) {
                        startVar.add(d[i][k]);
                        startVal.add(1);
                    }
                    /// c.add(d[k][i] == 0);
                    if((k < n - 1) && (i > k)) {
                        startVar.add(d[k][i]);
                        startVal.add(0);
                    }
                }
            }
            cplex.addMIPStart(startVar, startVal);
            startVal.end();
            startVar.end();
            LOG(INFO) << "Liao-You PFSP-WCT MP warm start solution set.";
        }

        std::vector<int> LiaoYouPFSP_WCT_Model::getSolutionAsPermutation(IloCplex &cplex, IloModel &model, 
                const int &n, NumVarMatrix &d) {
            ublas::matrix<int> d_matrix = ublas::zero_matrix<double>(n, n);
            // i = 1..n-1 , k = i+1..n
            for (IloInt i = 0; i < n - 1; i++) {
                for (IloInt k = i + 1; k < n; k++) {
                    // d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
                    if(cplex.getValue(d[i][k]) >= 0.95) {
                        d_matrix(i, k) = 1;
                    }
                }
            }
            return topological_sort_by_dfs(n, d_matrix);
        }

        /**
         * Return a [toplogical sort](https://en.wikipedia.org/wiki/Topological_sorting) of a directed 
         * graph `g` as a vector of vertices in topological order.
         * https://www.boost.org/doc/libs/1_55_0/libs/graph/doc/topological_sort.html
         */
        std::vector<int> LiaoYouPFSP_WCT_Model::topological_sort_by_dfs(const int &n, const ublas::matrix<int> &d_matrix) {
            typedef adjacency_list< boost::vecS, boost::vecS, directedS, property<vertex_color_t, default_color_type> > Graph;
            typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef std::pair<std::size_t,std::size_t> Pair;
            // Assemble a graph s.t. if edge (u,v) appears in the graph, then v comes before u in the ordering
            std::vector<Pair> edge_list;
            cout << "Graph edges: ";
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < n; k++) {
                    // d[i][k] = 1, if job i is scheduled any time before job k; 0, otherwise
                    if(d_matrix(i, k) == 1) {
                        edge_list.push_back(Pair(i, k));
                        cout << "i = " << i << ", k = " << k << "; ";
                    }
                }
            }
            cout << "\n";
            Graph G(n);
            for (std::size_t j = 0; j < edge_list.size(); ++j) {
                add_edge(edge_list[j].first, edge_list[j].second, G);
            }
            //Graph G(edge_list.begin(), edge_list.end(), n);
            typedef std::vector< Vertex > VertexOrder;
            VertexOrder vertex_order;
            topological_sort(G, std::back_inserter(vertex_order));

            std::vector<int> permutation;
            cout << "A topological ordering: ";
            // get the property map for vertex indices
            typedef property_map<Graph, vertex_index_t>::type IndexMap;
            IndexMap index = get(vertex_index, G);
            for (VertexOrder::reverse_iterator ii=vertex_order.rbegin(); ii!=vertex_order.rend(); ++ii) {
                cout << index[*ii] << " ";
                permutation.push_back(index[*ii] + 1);
            }
            cout << endl;
            return permutation;
        }

        std::vector<int> LiaoYouPFSP_WCT_Model::topological_sort_by_dfs_v2(const int &nv, const ublas::matrix<int> &g) {
            // graph `g` as a vector of vertices in topological order.
            // nv = number of vertices
            ublas::vector<int> vcolor = ublas::zero_vector<int>(nv);
            std::vector<int> verts;
            int w = 0;
            int u = 0;
            for(int v = 0; v < nv; v++) {  // v in 1:nv
                if (vcolor(v) != 0)  continue;
                std::stack<int> S;  // S = Vector{Int64}([v])
                S.push(v);
                vcolor(v) = 1;
                while(!S.empty()) {
                    u = S.top();
                    w = 0;
                    for(int n = 0; n < nv; n++) { //  n in 1:nv  # out_neighbors(g, u)
                        if(g(u, n) == 1) {
                            if(vcolor(n) == 1) {
                                throw std::runtime_error("The input graph contains at least one loop.");
                            } else if(vcolor(n) == 0) {
                                w = n;
                                break;
                            }
                        }
                    }
                    if(w != 0) {
                        vcolor(w) = 1;
                        S.push(w);
                    } else {
                        vcolor(u) = 2;
                        verts.push_back(u);
                        S.pop();
                    }
                }
            }
            std::reverse(verts.begin(), verts.end());
            return verts;
        }

    } /* namespace hybrid */
} /* namespace robust */
