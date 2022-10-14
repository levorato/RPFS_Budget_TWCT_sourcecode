//
// Created by Mario Costa Levorato Junior on 2019-02-01.
//

#ifndef FLOWSHOP_SOLVER_GRASP_RANDNEHT_H
#define FLOWSHOP_SOLVER_GRASP_RANDNEHT_H

#include "Job.h"
#include "SchedulingSolution.h"
#include "Inputs.h"
#include "Test.h"
#include "Randomness.h"
#include "RandomUtil.h"
#include "../../../util/include/TimeDateUtil.h"
#include <list>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <functional>
#include <algorithm>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o


// Number of possible alpha values in the adaptive approach
#define NUM_ALPHA_VALUES 4

namespace problem {
namespace pfsp {

    using namespace util;
/**
 * NEH heuristic with Taillard accelerations.
 *
 */
    template<class ProblemInstance>
    class RandNEHT {
    private:
        Test aTest;     // Test to run (including parameters and instance's name)
        Inputs<ProblemInstance> inputs;  // Instance inputs
        int nJobs;      // #Jobs
        int nMachines;  // #Machines
        std::vector<int> positions; // Array of randomly selected positions
        Job nextJob;
        // ************* Adaptive alpha attributes *************************************
        // Possible alpha values to draw from
        boost::numeric::ublas::vector<double> v_alpha;
        // Vector containing the probability of selecting (at random) a specific value of alpha
        std::vector<double> p_alpha;
        // the probability distribution for the alpha_i values
        std::discrete_distribution<> d_alpha;
        // sum of all solutions found using α = α_i
        boost::numeric::ublas::vector<double> sum_sol;
        // number of all solutions found using α = α_i
        boost::numeric::ublas::vector<long> num_sol;
        // auxiliary vector used in the adaptive approach probability calculation
        std::vector<double> q;
        // Statistic for the construction phase
        std::vector<double> alpha_history;
        std::vector< std::vector<double> > p_history;
        std::vector< std::vector<double> > q_history;
        // time limit measurement
        boost::timer::cpu_timer timer;
        double startTime;
        double elapsedTime;

    public:
        RandNEHT(const Test &test, const Inputs<ProblemInstance> &inputData) : aTest(test), inputs(inputData),
                                                nJobs(inputData.getNumberOfJobs()),
                                                nMachines(inputData.getNumberOfMachines()),
                                                positions(inputData.getNumberOfJobs(), 0),
                                                nextJob(0, 0), v_alpha(NUM_ALPHA_VALUES, -1.0),
                                                p_alpha(NUM_ALPHA_VALUES, 0.0), q(NUM_ALPHA_VALUES, 0.0),
                                                sum_sol(NUM_ALPHA_VALUES, 0.0), num_sol(NUM_ALPHA_VALUES, 0),
                                                d_alpha(p_alpha.begin(), p_alpha.end()), alpha_history(),
                                                timer(), startTime(0.0), elapsedTime(0.0) {
            timer.start();  // Time measurement
            startTime = TimeDateUtil::calculate_time_spent_ns(timer);
        }

        /**
         * Pseudo-random alpha value selector, based on a spinning wheel.
         * Uses a memory (function parameter) that, for each value of alpha, stores the sum of lower bound values and the
         * total number of times alpha_i was chosen to assume its worst-case proc. time value.
         * This method helps in the construction phase of the metaheuristic, executing the *adaptive* approach
         * when selecting the value of the alpha parameter.
         *  • Dynamic approach: The dynamic strategy often initializes the α value to a random
            value at each iteration of the GRASP metaheuristic. For instance, a uniform
            distribution may be used in the range [0.5, 0.9], or a decreasing non-uniform
            distribution.
            • Adaptive approach: In this strategy, the value of α is self-tuned. The value is updated
            automatically during the search function of the memory of the search.
               * Example 2.42 * Self-tuning in reactive GRASP. 
            In reactive GRASP, the value of the parameter α is periodically updated according to the quality
            of the obtained solutions. At each iteration, the parameter value is selected from a discrete
            set of possible values \Phi = {α_1, . . . , α_m}. The probability associated with each α_i is
            initialized to a same value p[i] = 1/m, i ∈ [1, . . . ,m]. The adaptive initialization in the
            reactive GRASP consists in updating these probabilities according to the quality of solutions
            obtained for each α_i. Let z* be the incumbent solution and avg_sol[i] the average value of all solutions
            found using α = α_i. Then, the probability p[i] for each value of α is updated as follows:
                p[i] = q[i] / sum( j = 1..m, q[j] ) , i ∈ [1, . . . ,m]
            where q[j] = z* / avg_sol[i]. 
            Hence, larger values of p[i] correspond to more suitable values for the parameter α[i].
            @return the index i of alpha_i, drawn from a uniform distribution.
         */
        inline unsigned calculate_adaptive_alpha(const double &z_star) {
            if(v_alpha(0) == -1.0) {  // initialize the alpha probability vector
                v_alpha(0) = 0.2;  v_alpha(1) = 0.4;  v_alpha(2) = 0.6;  v_alpha(3) = 0.8;
                //for(int i = 0; i < NUM_ALPHA_VALUES; i++) {
                //    p_alpha[i] = 1.0 / NUM_ALPHA_VALUES;  // p[i] = 1/m, i ∈ [1, . . . ,m]
                //}
                p_alpha[0] = 0.1;  p_alpha[1] = 0.1;  p_alpha[2] = 0.4;  p_alpha[3] = 0.4;
                // update the probability distribution
                std::discrete_distribution<>::param_type dd(p_alpha.begin(), p_alpha.end());
                d_alpha.param(dd);
            } else {
                update_adaptive_alpha_probability(z_star);
            }
            vector<double> p = d_alpha.probabilities();
            // Using the probabilities p_alpha(i), draw an index of alpha_i in the range {0, ..., NUM_ALPHA_VALUES - 1}
            int i = d_alpha(rng);
            return i;
        }

        /**
         * Update the probability of alpha_i being chosen in the adaptive approach.
         * Let z* be the incumbent solution and avg_sol[i] the average value of all solutions
         *  found using α = α_i. Then, the probability p[i] for each value of α is updated as follows:
         *      p[i] = q[i] / sum( j = 1..m, q[j] ) , i ∈ [1, . . . ,m]
         *  where q[i] = z* / avg_sol[i].
         * @param z_star value of the incumbent solution
         */
        inline void update_adaptive_alpha_probability(const double &z_star) {
            bool all_values_visited = true;
            for (int i = 0; i < NUM_ALPHA_VALUES; i++) {
                if(num_sol(i) == 0) {
                    all_values_visited = false;
                }
            }
            if(all_values_visited) {
                double sum_q = 0.0;
                for (int i = 0; i < NUM_ALPHA_VALUES; i++) {
                    q[i] = z_star / (sum_sol(i) / num_sol(i));
                    sum_q += q[i];
                }
                for (int i = 0; i < NUM_ALPHA_VALUES; i++) {
                    p_alpha[i] = q[i] / sum_q;  // p[i] = 1/m, i ∈ [1, . . . ,m]
                }
                // update the probability distribution
                std::discrete_distribution<>::param_type dd(p_alpha.begin(), p_alpha.end());
                d_alpha.param(dd);
            }
        }

        SchedulingSolution<ProblemInstance> build_solution_RCL(std::vector<Job> &effList, const double &incumbent_value,
                const bool &adaptive) {
            // 0. Define a new solution
            SchedulingSolution<ProblemInstance> currentSol(nJobs, nMachines, inputs);
            Randomness<ProblemInstance> random(aTest, inputs);
            double alpha = 0.0;
            unsigned alpha_idx = 0;
            if(adaptive) {
                alpha_idx = calculate_adaptive_alpha(incumbent_value);
                alpha = v_alpha(alpha_idx);
                alpha_history.push_back(alpha);
                q_history.push_back(q);
                p_history.push_back(p_alpha);
            } else {  // alpha lies in the interval [beta1, beta2]
                alpha = aTest.getBeta1() + aTest.getRandomStream().nextDouble() * (aTest.getBeta2() - aTest.getBeta1());
            }
            std::list<Job> s;
            std::copy( effList.begin(), effList.end(), std::back_inserter( s ) );  // linked list of Jobs
            std::vector< std::list<Job>::iterator > job_ptr(nJobs + 1, s.end());
            for (int i = 0; i < nJobs; i++) {
                std::vector<Job> RCL;
                elapsedTime = ElapsedTime::calcElapsed(startTime, TimeDateUtil::calculate_time_spent_ns(timer));
                if(elapsedTime >= aTest.getMaxTime()) {  // Check for time limit
                    alpha = 1.0;  // switch to fully random (faster) construction if time limit is exceeded
                }
                if(alpha < 1.0) {
                    for (std::list<Job>::iterator it = s.begin(); it != s.end(); ++it) {  // For each Job in the list s
                        // Add nextJob j to the end of currentSol (partial solution)
                        currentSol.jobs[i] = *it;
                        it->cost = currentSol.calcTotalCosts(i + 1, aTest.get_force_recalculate_obj());
                        job_ptr[it->getId()] = it;
                    }
                    s.sort( [](const Job &a, const Job &b) {return a.cost < b.cost; });
                    unsigned element_count = 0;
                    unsigned RCL_size = std::max((int)ceil(alpha * s.size()), 1);
                    for (Job j : s) {
                        RCL.push_back(j);
                        element_count++;
                        if(element_count >= RCL_size)  break;
                    }
                } else {
                    for (std::list<Job>::iterator it = s.begin(); it != s.end(); ++it) {  // For each Job in the list s
                        Job j = *it;
                        job_ptr[j.getId()] = it;
                        RCL.push_back(j);
                    }
                }
                int list_size = RCL.size();
                int index = random.getRandomPosition(list_size, "u");
                if(list_size == 0 || index >= list_size){
                    cerr << "[RandNEHT] Error in GRASP construction!\n";
                    break;
                }
                nextJob = RCL[index];
                currentSol.jobs[i] = nextJob;
                s.erase(job_ptr[nextJob.getId()]);
                //std::fill_n(job_cost.begin(), job_cost.size(), 0);  // reset the job_cost array
                //std::fill_n(job_ptr.begin(), job_ptr.size(), s.end());
            }
            currentSol.setCosts(currentSol.calcTotalCosts(nJobs, aTest.get_force_recalculate_obj()));
            // Update adaptive alpha statistics
            if(adaptive) {
                sum_sol(alpha_idx) += currentSol.getCosts();
                num_sol(alpha_idx)++;
            }
            return currentSol;
        }

        SchedulingSolution<ProblemInstance> solve(std::vector<Job> &effList, const bool &useRandomSelection,
                const double &incumbent_value = 0.0, const bool &use_adaptive_construction = false) {
            // 0. Define a new solution
            SchedulingSolution<ProblemInstance> currentSol(nJobs, nMachines, inputs);
            Randomness<ProblemInstance> random(aTest, inputs);
            if(useRandomSelection && (use_adaptive_construction || (aTest.getBeta1() < 1.0))) {  // RCL approach
                return build_solution_RCL(effList, incumbent_value, use_adaptive_construction);
            } else {
                // 1. Calculate the array of randomly selected positions of jobs in effList
                if (!useRandomSelection) // classical NEH solution
                    for (int i = 0; i < nJobs; i++)
                        positions[i] = i;
                else {  // fully randomized NEHT approach (== GRASP alpha = 1.0)
                    positions = random.calcPositionsArray(aTest.getDistribution()); // Randomized NEH solution
                }
                // 2. Insert the first job in the solution (not an empty solution anymore)
                nextJob = effList[positions[0]];
                currentSol.jobs[0] = nextJob;
                // 3. Complete the solution with the remaining jobs
                for (int i = 1; i < nJobs; i++) {
                    // Add nextJob to the end of currentSol (partial solution)
                    nextJob = effList[positions[i]];
                    currentSol.jobs[i] = nextJob;
                    elapsedTime = ElapsedTime::calcElapsed(startTime, TimeDateUtil::calculate_time_spent_ns(timer));
                    if(elapsedTime < aTest.getMaxTime()) {  // Check for time limit
                        // Try to improve currentSol by shifting nextJob to the left
                        currentSol.improveByShiftingJobToLeft(i, aTest.is_first_improvement(), aTest.get_force_recalculate_obj());
                    }
                }
                return currentSol;
            }
        }

        const std::vector<double> &get_alpha_history() const {
            return alpha_history;
        }

        const std::vector< std::vector<double> > &get_alpha_quality() const {
            return q_history;
        }

        const std::vector< std::vector<double> > &get_alpha_probability() const {
            return p_history;
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_RANDNEHT_H
