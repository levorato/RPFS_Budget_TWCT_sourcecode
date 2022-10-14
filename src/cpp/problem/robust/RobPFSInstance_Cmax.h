//
// Created by mlevorato on 7/26/19.
//

#ifndef FLOWSHOP_SOLVER_ROBPFSINSTANCE_H
#define FLOWSHOP_SOLVER_ROBPFSINSTANCE_H

#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/multi_array.hpp>

#include "../deterministic/PFSProblem.h"
#include "RobPFSInstance.h"
#include "../GlobalTypes.h"
#include "../PFSP_Parameters.h"
#include "budget/PFSP_Cmax_Budget_WorstCase.h"
#include "./budget/PFSPBudgetScenario.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    class RobPFSInstance_Cmax : public RobPFSInstance {
    public:
        /**
         * Build a Robust PFS Instance for the PFSP ADRS Problem, discrete scenarios variation (Kouvelis et al., 2000).
         * @param instance_name
         * @param lbd
         * @param M
         * @param N
         */
        RobPFSInstance_Cmax(const string& instance_name, const unsigned &p_lambda, const unsigned &M, const unsigned &N,
                            const Robust_Parameters& p_rob_params) :
                RobPFSInstance(instance_name, p_lambda, M, N, p_rob_params) {

        }

        /**
         * Build a Robust PFS Instance for the PFSP ADRS Problem, continuous interval variation (Kouvelis et al., 2000).
         * @param instance_name
         * @param lbd
         * @param M
         * @param N
         */
        RobPFSInstance_Cmax(const string& instance_name, const unsigned &M, const unsigned &N,
                            const Robust_Parameters& p_rob_params) :
                RobPFSInstance(instance_name, M, N, p_rob_params) {
        }

        /**
         * Build a Robust PFS Instance for the budget constraints problem variation (Ying, 2015).
         * @param instance_name
         * @param M
         * @param N
         * @param pp_bar
         * @param pp_hat
         */
        RobPFSInstance_Cmax(const string& instance_name, const unsigned &M, const unsigned &N,
                            const boost::numeric::ublas::matrix<double> &pp_bar,
                            const boost::numeric::ublas::matrix<double> &pp_hat,
                            const Robust_Parameters& p_rob_params) :
                RobPFSInstance(instance_name, M, N, pp_bar, pp_hat, p_rob_params)  {
        }
        ~RobPFSInstance_Cmax() {   }

        /*******************************************************************************
         * PUBLIC METHOD calcTotalCosts()
        ******************************************************************************/
        double calcTotalCosts(std::vector<Job> &pi, int nUsedJobs, unsigned recalculate_obj = 0) const {
            if(robust_type == RobustType::budget) {
                return PFSP_Cmax_Budget_WorstCase::worst_case_cmax_dp(nUsedJobs, p_bar, p_hat, rob_params, pi);
            } else {  // adrs // TODO INVOKE KOUVELIS ADRS CONTINUOUS/DISCRETE WORST-CASE HERE
                throw std::invalid_argument("Invalid robust_type argument.");
            }
        }

        double calculate_objective(const std::vector<int> &pi) const {
            std::vector<Job> perm;
            for (int j : pi)
                perm.push_back(Job(j - 1, 0));
            return calcTotalCosts(perm, perm.size(), 0);
        }

        /*******************************************************************************
         * PUBLIC METHOD improveByShiftingJobToLeft()
         * k is the position of the job on the right extreme.
         *
         * Note: pi and makespan are in_out parameters.
         * Returns true if the resulting makespan has changed.
         ******************************************************************************/
        bool improveByShiftingJobToLeft(std::vector<Job> &pi, double &makespan, int k, long &num_visits,
                bool first_improvement, unsigned force_recalculate_obj) {
            int bestPosition = k;
            double minMakespan = std::numeric_limits<double>::max();
            double newMakespan = std::numeric_limits<double>::max();

            // Calculate bestPosition (0...k) and minMakespan (mVector)
            for (int i = k; i >= 0; i--) {
                num_visits++;
                std::vector<Job> s(pi);
                // Shift left each job in position k to each possible position i in {0, ..., k - 1}
                Job auxJob = s[k];
                for (int j = k; j > i; j--)
                    s[j] = s[j - 1];
                s[i] = auxJob;
                newMakespan = PFSP_Cmax_Budget_WorstCase::worst_case_cmax_dp(n, p_bar, p_hat, rob_params, s);
                // TIE ISSUE #2 - In case of tie, do swap
                if (newMakespan <= minMakespan) {
                    minMakespan = newMakespan;
                    bestPosition = i;
                    if(first_improvement && newMakespan < minMakespan)  break;
                }
            }
            // Update solution with bestPosition and minMakespan
            if (bestPosition < k) // if i == k do nothing
            {
                Job auxJob = pi[k];
                for (int i = k; i > bestPosition; i--)
                    pi[i] = pi[i - 1];
                pi[bestPosition] = auxJob;
            }
            makespan = minMakespan;
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
         * Given a scenario, calculate the makespan (Cmax) objective.
         * The scenario parameter contains the list of jobs j that will deviate to
         * the worst-case value (P_bar[i, j] + P_hat[i, j]) in machine M_i, according to the
         * budgeted uncertainty set.
         * @param pi permutation
         * @param scenario list of jobs which will deviated to the worst-case processing time on each machine.
         * @return makespan
         */
        double calculate_cmax_given_scenario(const Permutation &pi, PFSPBudgetScenario &scenario) {
            boost::numeric::ublas::matrix<Time> p_time = scenario.getProcessingTimeMatrix(this);
            return PFSProblem::calculateCmax(pi, m, n, p_time);
        }

        string objective_name() {
            return string("Cmax");
        }
    };

}

#endif //FLOWSHOP_SOLVER_ROBPFSINSTANCE_H
