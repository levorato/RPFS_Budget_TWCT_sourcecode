//
// Created by mlevorato on 5/1/20.
//

#ifndef PFSP_ROBPFSINSTANCE_WCT_H
#define PFSP_ROBPFSINSTANCE_WCT_H

#include <stdexcept>
#include <memory>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/multi_array.hpp>

#include "../deterministic/PFSProblem.h"
#include "RobPFSInstance.h"
#include "../GlobalTypes.h"
#include "../PFSP_Parameters.h"
#include "./budget/PFSP_WCT_Budget_WorstCase.h"
#include "./budget/PFSPBudgetScenario.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    class RobPFSInstance_WCT : public RobPFSInstance {
    public:
        /**
         * Build a Robust PFS Instance for the budget constraints problem variation (Ying, 2015).
         * @param instance_name
         * @param M
         * @param N
         * @param pp_bar
         * @param pp_hat
         */
        RobPFSInstance_WCT(const string& instance_name, const unsigned &M, const unsigned &N,
                            const boost::numeric::ublas::matrix<double> &pp_bar,
                            const boost::numeric::ublas::matrix<double> &pp_hat,
                            const boost::numeric::ublas::vector<double> &p_w,
                            const Robust_Parameters& p_rob_params) :
                RobPFSInstance(instance_name, M, N, pp_bar, pp_hat, p_rob_params), w(p_w), num_runs(0),
                E(ublas::zero_matrix<double>(N + 2, M + 2)), R(ublas::zero_matrix<double>(N + 2, M + 2)),
                l_kl(nullptr), p() {
            double gamma = ((m * n * rob_params.budget_gamma[0]) / 100.0);
            p = calculate_p_matrix(m, n, p_bar, p_hat, gamma);
        }

        /*******************************************************************************
         * PUBLIC METHOD calcTotalCosts()
        ******************************************************************************/
        double calcTotalCosts(std::vector<Job> &pi, int nUsedJobs, unsigned recalculate_obj) {
            if(robust_type == RobustType::budget) {
                return calculate_worst_case_wct(pi, nUsedJobs, recalculate_obj);  // use approximate DP worst-case value?
            } else {  // adrs // TODO INVOKE KOUVELIS ADRS CONTINUOUS/DISCRETE WORST-CASE HERE
                throw std::invalid_argument("Invalid robust_type argument.");
            }
        }

        double calculate_objective(const std::vector<int> &pi) {
            std::vector<Job> perm;
            for (int j : pi)
                perm.push_back(Job(j - 1, 0));
            return calcTotalCosts(perm, perm.size(), 0);
        }

        /**
         * Calculates the worst-case Weighted Completion Time, given a permutation pi and the number of jobs nUsedJobs.
         * By default, the calculation is based on the DP worst-case procedure, which does not yield optimum results.
         * Therefore, every t iterations (parameter rob_params.obj_update_freq), the objective function of the incumbent
         * is updated according to the exact approach, based on the (slower) MILP worst-case procedure.
         * @param pi
         * @param nUsedJobs
         * @param calculation_type Frequency for recalculating the metaheuristic objective function (approximating the obj):
         *        0 == always recalculate, i.e., use worst-case MIP (exact solution);
                  1 == no recalculation (use DP v1.0 all the time). Only use worst-case MIP (exact solution) to give final result;
                  2 == no recalculation (use DP v2.0 all the time). Only use worst-case MIP (exact solution) to give final result;
                  3 == obj recalculated at the end of each MH iteration (use DP v1.0 otherwise);
                  4 == obj recalculated at the end of each MH iteration (use DP v2.0 otherwise);
                  5 == obj recalculated at the end of each MH iteration and also at each VND improvement (use DP v1.0 otherwise);
                  6 == obj recalculated at the end of each MH iteration and also at each VND improvement (use DP v2.0 otherwise).
         * @return objective function value (approximated or exact).
         */
        double calculate_worst_case_wct(std::vector<Job> &pi, int nUsedJobs, unsigned calculation_type,
                                        bool verbose = false,
                                        const std::pair<double,double> cutoff_value = std::pair<double,double>()) {
            num_runs++;
            // Use deterministic objective validation if gamma == 0 or gamma == 100 (saves time)
            if((rob_params.budget_gamma[0] == 0) || (rob_params.budget_gamma[0] == 100)) {
                std::vector<int> perm;
                for (int i = 0; i < n; i++)
                    perm.push_back(pi[i].id);
                return PFSProblem::calculate_wct(perm, m, nUsedJobs, p, w);
            } else {
                if (calculation_type == 0) {
                    PFSPBudgetScenario scenario = PFSP_WCT_Budget_WorstCase::worst_case_wct_mip(nUsedJobs, p_bar, p_hat,
                                                                                                rob_params, pi, w,
                                                                                                (num_runs == 1), false,
                                                                                                0, verbose, cutoff_value);
                    return scenario.value;
                } else if ((calculation_type == 1) || (calculation_type == 3) || (calculation_type == 5)) {
                    return PFSP_WCT_Budget_WorstCase::worst_case_wct_dp(nUsedJobs, p_bar, p_hat, rob_params, pi, w);
                } else {  // 2, 4, or 6
                    return PFSP_WCT_Budget_WorstCase::worst_case_wct_dp_v2(nUsedJobs, p_bar, p_hat, rob_params, pi, w);
                }
            }
        }

        /*******************************************************************************
         * PUBLIC METHOD improveByShiftingJobToLeft()
         * k is the position of the job on the right extreme.
         * Note: pi and wct are in_out parameters.
         * Returns true if the resulting makespan has changed.
         ******************************************************************************/
        bool improveByShiftingJobToLeft(std::vector<Job> &pi, double &wct, int k, long &num_visits,
                                        bool first_improvement, unsigned force_recalculate_obj) {
            int bestPosition = k;
            double minWCT = std::numeric_limits<double>::max();
            double newWCT = std::numeric_limits<double>::max();

            // Calculate bestPosition (0...k) and minWCT (mVector)
            for (int i = k; i >= 0; i--) {
                num_visits++;
                std::vector<Job> s(pi);
                // Shift left each job in position k to each possible position i in {0, ..., k - 1}
                Job auxJob = s[k];
                for (int j = k; j > i; j--)
                    s[j] = s[j - 1];
                s[i] = auxJob;
                double newWCT = calculate_worst_case_wct(s, n, force_recalculate_obj);  // use approximate DP worst-case value
                // TIE ISSUE #2 - In case of tie, do swap
                if (newWCT <= minWCT) {
                    minWCT = newWCT;
                    bestPosition = i;
                    if(first_improvement && newWCT < minWCT)  break;
                }
            }
            // Update solution with bestPosition and minWCT
            if (bestPosition < k) // if i == k do nothing
            {
                Job auxJob = pi[k];
                for (int i = k; i > bestPosition; i--)
                    pi[i] = pi[i - 1];
                pi[bestPosition] = auxJob;
            }
            wct = minWCT;  // return value (in-out parameter)
            return true;
        }

        double totalTime() const {
            double T = 0;
            if(robust_type == RobustType::budget) {
                for (unsigned i = 1; i <= m; i++)
                    for (unsigned j = 1; j <= n; j++)
                        T += p_bar(i, j);
            } else {  // adrs
                if(isDiscrete()) {
                    boost::numeric::ublas::matrix<double> p_all(m + 1, n + 1, 0.0);
                    for (boost::numeric::ublas::matrix<double> p : p_lambda) {
                        for (unsigned i = 1; i <= m; i++)
                            for (unsigned j = 1; j <= n; j++)
                                p_all(i, j) += p(i, j);
                    }
                    for (unsigned i = 1; i <= m; i++)
                        for (unsigned j = 1; j <= n; j++)
                            T += p_all(i, j) / lambda;
                } else {  // continuous
                    for (unsigned i = 1; i <= m; i++)
                        for (unsigned j = 1; j <= n; j++)
                            T += (p_over(i, j) - p_under(i, j)) / 2;
                }
            }
            return T;
        }

        /**
         * Given a scenario, calculate the weighted sum of completion times (WCT) objective.
         * The scenario parameter contains the list of jobs j that will deviate to
         * the worst-case value (P_bar[i, j] + P_hat[i, j]) in machine M_i, according to the
         * budgeted uncertainty set.
         * @param pi permutation
         * @param scenario list of jobs which will deviated to the worst-case processing time on each machine.
         * @return weighted sum of completion times (WCT)
         */
        double calculate_wct_given_scenario(const Permutation &pi, PFSPBudgetScenario &scenario) {
            boost::numeric::ublas::matrix<Time> p_time = scenario.getProcessingTimeMatrix(this);
            return PFSProblem::calculate_wct(pi, m, n, p_time, w);
        }

        string objective_name() {
            return string("WCT");
        }

        ublas::matrix<double> calculate_p_matrix(const int &m, const int &n, const ublas::matrix<double> &p_bar,
                                                        const ublas::matrix<double> &p_hat, const double &gamma) {
            ublas::matrix<double> p_time(p_bar);
            // obtain the gamma largest processing times p_hat and oscillate them
            std::vector< std::tuple<int, int, int> > v;
            for(int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                    v.emplace_back(p_hat(i, j), i, j);
                }
            }
            // Using sort() function to sort by 1st element of tuple, in non-ascending order
            std::sort(v.begin(), v.end(), [](const std::tuple<int, int, int> &i, const std::tuple<int, int, int> &j) {
                return i > j;
            });
            int max_gamma = (int)floor(gamma);
            for(int x = 0; x < max_gamma; x++) {
                int i = std::get<1>(v[x]);
                int j = std::get<2>(v[x]);
                p_time(i, j) += p_hat(i, j);
            }
            return p_time;
        }

        /**
         * l_kl(j) or P_j(i, h) : shorthand notation for partial sums of processing times over several
         * adjacent machines (between machines k and l, inclusive : k <= u <= l).
         * Note: calculation of l_kl is based solely on nominal processing times matrix p_bar, if Gamma == 0.
         * @param instance
         * @return
         */
        void initialize_l_kl_matrix(const Robust_Parameters &rob_params) {
            // IMPORTANT: Initialization of l_kl auxiliary vectors used by lower bounds
            // For each machine pair (k, l)
            // Sum of processing times between machines k + 1 and l - 1 : Create a 3D array
            ublas::matrix<double> p_time = this->p_bar;
            if(rob_params.budget_gamma[0] == 100) {  // oscillate the processing time of all operations
                p_time += this->p_hat;
            }
            l_kl = std::make_shared<boost::multi_array<double, 3>>(boost::extents[n + 2][m + 2][m + 2]);
            std::fill_n(l_kl->data(), l_kl->num_elements(), 0);
            for (int j = 1; j <= n; j++) {
                for (int k = 1; k <= m; k++) {
                    for (int l = k; l <= m; l++) {
                        for (int u = k; u <= l; u++) {
                            (*l_kl)[j][k][l] += p_time(u, j);
                        }
                    }
                }
            }
        }

        boost::multi_array<double, 3>& get_l_kl_matrix() {
            if(l_kl.get() == nullptr) {
                initialize_l_kl_matrix(rob_params);
            }
            return *l_kl;
        }

        const ublas::matrix<double>& get_oscillated_p_matrix() const {
            return p;
        }

        /** Vector of job weights */
        boost::numeric::ublas::vector<double> w;
        /** Counter which controls the number of times the objective function was called */
        long num_runs;

        // The following data structures are used for Lower Bound calculation in the combinatorial branch-and-bound
        std::shared_ptr<boost::multi_array<double, 3>> l_kl;
        // p is the matrix of processing times obtained by oscillating the gamma operations with smallest ptimes
        ublas::matrix<double> p;
        // Etk is an underestimate of the earliest start time of the t-th job on machine k
        ublas::matrix<double> E;
        // R(t, k) is an underestimate of the earliest completion time of the t-th job on machine k
        ublas::matrix<double> R;
    };
}

#endif //PFSP_ROBPFSINSTANCE_WCT_H
