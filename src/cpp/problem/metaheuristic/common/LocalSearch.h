//
// Created by Mario Costa Levorato Junior on 2019-02-01.
//

#ifndef FLOWSHOP_SOLVER_GRASP_LOCALSEARCH_H
#define FLOWSHOP_SOLVER_GRASP_LOCALSEARCH_H

#include "SchedulingSolution.h"
#include "Inputs.h"
#include "Randomness.h"
#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o
#include <glog/logging.h>
#include <algorithm>
#include "../../../util/include/TimeDateUtil.h"

namespace problem {
namespace pfsp {

    using namespace util;

    /**
     * This class contains all the local search procedures.
     */
    template<class ProblemInstance>
    class LocalSearch {
    private:
        Inputs<ProblemInstance> inputs; // Instance inputs
        Randomness<ProblemInstance> random;

        std::vector<int> positions, positions_nJobs_minus_2; // Arrays of randomly selected positions
        int nJobs; // #Jobs

    public:
        LocalSearch(const Test &test, const Inputs<ProblemInstance> &inputData) : inputs(inputData),
            random(test, inputData), positions(inputData.getNumberOfJobs(), 0),
            positions_nJobs_minus_2(inputData.getNumberOfJobs() - 2, 0),
            nJobs(inputData.getNumberOfJobs()) {
            positions = random.calcPositionsArray("u"); //uniform
            positions_nJobs_minus_2 = RandomUtil::generate_random_vector(0, inputData.getNumberOfJobs() - 2);
        }

        /**
         * Local search routine used by the ILS metaheuristic.
         * Returns the number of visited solutions during local search.
         * @param aSol
         * @return
         */
        inline long globalImprovement(SchedulingSolution<ProblemInstance> &aSol,
                const bool &first_improvement, long &num_improvements, double elapsed_time, double time_limit,
                const unsigned &force_recalculate_obj) {
            bool hasImproved = false;
            double beforeCosts = 0.0;
            long num_visited_solutions = 0;
            boost::timer::cpu_timer timer;
            timer.start();  // Time measurement
            do {
                hasImproved = false;
                beforeCosts = aSol.getCosts();

                num_visited_solutions += randomJobShifting(aSol, first_improvement, elapsed_time, time_limit,
                                                           force_recalculate_obj);

                if (beforeCosts > aSol.getCosts()) {
                    hasImproved = true;
                    num_improvements++;
                    if(first_improvement)  break;
                }
                double elapsed = elapsed_time + TimeDateUtil::calculate_time_spent(timer);
                if(elapsed >= time_limit)  break;
            } while (hasImproved);
            timer.stop();
            return num_visited_solutions;
        }

        inline SchedulingSolution<ProblemInstance> search_neighborhood_by_number(const unsigned &r,
                const SchedulingSolution<ProblemInstance> &S_star, const bool &first_improvement,
                long &num_visited_solutions, long &num_improvements, bool &hasImproved,
                const unsigned &force_recalculate_obj, double elapsed_time, double time_limit) {
            hasImproved = false;
            SchedulingSolution<ProblemInstance> aSol(S_star);
            if(force_recalculate_obj >= 5) {
                aSol.costs = aSol.calcTotalCosts(aSol.getNJobs(), force_recalculate_obj);
            }
            double beforeCosts = aSol.getCosts();
            int d = 2;  // TODO Select the number of swaps via random distribution ?
            switch (r) {
                default:
                case 1:
                    num_visited_solutions += enhancedInsertion(aSol, first_improvement, force_recalculate_obj);
                    break;
                case 2:
                    num_visited_solutions += enhancedSwap(aSol, first_improvement, force_recalculate_obj);
                    break;
                case 3:
                    num_visited_solutions += randomDestructionConstruction(aSol, d, first_improvement,
                            force_recalculate_obj);
                    break;
                case 4:
                    num_visited_solutions += randomJobShifting(aSol, first_improvement, elapsed_time, time_limit,
                                                               force_recalculate_obj);
                    break;
                case 5:
                    num_visited_solutions += partialImprovement(aSol, first_improvement, force_recalculate_obj);
                    break;
                case 6:
                    num_visited_solutions += randomSwapping(aSol, force_recalculate_obj);
                    break;
                case 7:
                    num_visited_solutions += randomInsertion(aSol, force_recalculate_obj);
                    break;
                case 8:
                    num_visited_solutions += adjacentSwap(aSol, d, force_recalculate_obj);
                    break;
            }
            if (beforeCosts > aSol.getCosts()) {
                hasImproved = true;
                num_improvements++;
            }
            return aSol;
        }

        /**
         * Executes the Variable Neighborhood Descent (VND) algorithm.
         * Repeatedly derives the local optimum solution Sol(l) in the l-neighborhood of
         * the current solution Sol.
         * @param mh_sol initial problem solution to be explored by VND local search
         * @param vnd_permutation the permutation of neighborhood numbers to be used by VND (used in parametrization)
         * @param l_size the size of the VND neighborhood
         * @param first_improvement true if first-improvement is enabled (stop search at 1st improvement found)
         * @param time_limit time limit in seconds for the VND procedure
         * @param time_spent_so_far time spent so far by the whole metaheuristic (to be subtracted from time limit)
         * @param mh_iteration current metaheuristic iteration (to be printed on time results csv file)
         * @param time_results stringstream containing the time results CSV content (with each solution improvement found)
         * @param force_recalculate_obj flag that indicates
         * @param random_vnd if true, apply the RVND search (shuffle the vnd_permutation)
         */
        inline long VariableNeighborhoodDescent(SchedulingSolution<ProblemInstance> &mh_sol,
                const std::vector<unsigned> &vnd_perm, const int &l_size,
                const bool &first_improvement, const long &time_limit, const double &time_spent_so_far,
                const long &mh_iteration, stringstream &time_results, long &num_improvements,
                const unsigned &force_recalculate_obj,
                const bool &random_vnd = false) {
            // k is the current neighborhood distance in the local search
            int r = 1, iteration = 0;
            long num_visited_solutions = 0;
            SchedulingSolution<ProblemInstance> S_star(mh_sol); // S* := mh_sol
            if (force_recalculate_obj >= 5) {  // force recalculation of the objective function in VND improvement
                double S_star_cost = mh_sol.forceUpdateCosts(force_recalculate_obj);
                mh_sol.costs = S_star_cost;
                S_star.costs = S_star_cost;
            }
            // LOG(INFO)<< "VND local search (neighborhood size = " << l_size << ") ...\n";
            boost::timer::cpu_timer timer;
            timer.start();  // Time measurement
            double elapsed = 0.0;
            double startTime = TimeDateUtil::calculate_time_spent_ns(timer);
            std::vector<unsigned> vnd_permutation(vnd_perm);
            if (random_vnd) {  // shuffle the vnd neighborhood order
                std::random_shuffle(vnd_permutation.begin(), vnd_permutation.end());
            }
            while ((r <= l_size) && (time_spent_so_far + elapsed < time_limit)) {
                // N := Nl(S*) : Apply a local search in S_star using the r-neighborhood
                bool hasImproved = false;
                SchedulingSolution<ProblemInstance> S_l = search_neighborhood_by_number(vnd_permutation[r], S_star,
                                                                                        first_improvement,
                                                                                        num_visited_solutions,
                                                                                        num_improvements, hasImproved,
                                                                                        force_recalculate_obj,
                                                                                        time_spent_so_far + elapsed,
                                                                                        time_limit);
                if (hasImproved) {
                    if (force_recalculate_obj >= 5) {  // force recalculation of the objective function in VND improvement
                        double S_l_cost = S_l.forceUpdateCosts(force_recalculate_obj);
                        if (S_l_cost < mh_sol.costs) {
                            S_l.costs = S_l_cost;
                            mh_sol = S_l;
                        }
                    }
                    S_star = S_l;
                    r = 1;  // return to the first neighborhood
                    // Registers the best result at time intervals
                    //notifyNewValue(S_star, time_spent_so_far + elapsed, mh_iteration, time_results);
                } else {  // no better result found in neighborhood, go to the next one
                    r++;
                }
                iteration++;
                elapsed = TimeDateUtil::calculate_time_spent(timer);
            }
            VLOG(1) << "VND local search done. Obj = " << S_star.getCosts() <<
                                     ". Time spent: " << (time_spent_so_far + elapsed) << " s";
            if (force_recalculate_obj <= 4) {
                if (S_star.costs < mh_sol.costs) {
                    mh_sol = S_star;
                }
            }
            return num_visited_solutions;
        }

        inline void notifyNewValue(SchedulingSolution<ProblemInstance> &sol, const double& timeSpent,
                                   const int& iteration, stringstream& timeResults) {
            timeResults << fixed << setprecision(4) << timeSpent << "," << sol.costs
                        << "," << (iteration+1) << "\n";
        }

        /*******************************************************************************
        * randomShiftingJob()
        *
        * Tries to improve the solution by taking each job, randomly without repetition,
        * moving it to the last position and improve by shifting it to the left.
        * Returns the number of visited solutions in the local search.
        ******************************************************************************/
        inline long randomJobShifting(SchedulingSolution<ProblemInstance> &aSol, const bool &first_improvement,
                                      double elapsed_time, double time_limit,
                                      const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            long num_visited_solutions = 0;
            boost::timer::cpu_timer timer;
            timer.start();  // Time measurement

            for (int i = 0; i < nJobs - 1; i++) {
                SchedulingSolution<ProblemInstance> aTemp(aSol);
                int j = positions[i];  // obtain random number between 0 and nJobs - 1
                if (j < nJobs - 1) {
                    Job aJob = aTemp.jobs[j];
                    //System.arraycopy(aTemp.jobs, j + 1, aTemp.jobs, j, nJobs - 1 - j);
                    for(; j < nJobs - 1; j++ ) {
                        aTemp.jobs[j] = aTemp.jobs[j + 1];
                    }
                    aTemp.jobs[nJobs - 1] = aJob;
                }
                // TODO FIXME Check if the recalculation of the objective function could be made only at the end
                num_visited_solutions += aTemp.improveByShiftingJobToLeft(nJobs - 1, first_improvement,
                                                                         force_recalculate_obj);
                if(aTemp.getCosts() < aSol.getCosts()) {
                    aSol = aTemp;
                    if (first_improvement) {  // solution improved
                        break;
                    }
                }
                double elapsed = elapsed_time + TimeDateUtil::calculate_time_spent(timer);
                if(elapsed >= time_limit)  break;
            }
            timer.stop();
            return num_visited_solutions;
        }

        /*******************************************************************************
        * enhancedSwap()
        * Swaps two jobs at random and tries to improve the solution by shifting each
        * modified to job to the left.
        * Perturbation method
        *******************************************************************************/
        inline long enhancedSwap(SchedulingSolution<ProblemInstance> &aSol, const bool &first_improvement,
                                 const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            int posA = positions[0];   // random.getRandomPosition(nJobs, "uniform");
            int posB = positions[nJobs - 1];  // random.getRandomPosition(nJobs, "uniform");

            long num_visited_solutions = 0;
            swapJobs(aSol, posA, posB);
            if (posA < posB) {
                num_visited_solutions += aSol.improveByShiftingJobToLeft(posA, first_improvement, force_recalculate_obj);
                num_visited_solutions += aSol.improveByShiftingJobToLeft(posB, first_improvement, force_recalculate_obj);
            } else {
                num_visited_solutions += aSol.improveByShiftingJobToLeft(posB, first_improvement, force_recalculate_obj);
                num_visited_solutions += aSol.improveByShiftingJobToLeft(posA, first_improvement, force_recalculate_obj);
            }
            aSol.setCosts(aSol.calcTotalCosts(nJobs, force_recalculate_obj));
            return num_visited_solutions;
        }

        inline long swapJobs(SchedulingSolution<ProblemInstance> &aSol, const int &posA, const int &posB) {
            Job aux = aSol.jobs[posA];
            aSol.setJob(posA, aSol.jobs[posB]);
            aSol.setJob(posB, aux);
            return 1;
        }

        /**
         * Swaps 2 jobs at random.
         * @param aSol
         * @return
         */
        inline long randomSwapping(SchedulingSolution<ProblemInstance> &aSol,
                                   const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            int posA = positions[0];   // random.getRandomPosition(nJobs, "uniform");
            int posB = positions[nJobs - 1];  // random.getRandomPosition(nJobs, "uniform");
            swapJobs(aSol, posA, posB);
            aSol.setCosts(aSol.calcTotalCosts(nJobs, force_recalculate_obj));
            return 1;
        }

        /**
         * Swaps n pairs of adjancent jobs from the permutation.
         * @param aSol
         * @param nPairs
         * @return
         */
        inline long adjacentSwap(SchedulingSolution<ProblemInstance> &aSol, const int &nPairs,
                                 const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            for (int i = 0; i < nPairs; i++) {
                int pos = positions[i];  // obtain random number between 0 and nJobs - 1, without replacement
                if(pos == nJobs - 1)  continue;
                swapJobs(aSol, pos, pos + 1);
            }
            aSol.setCosts(aSol.calcTotalCosts(nJobs, force_recalculate_obj));
            return 1;
        }

        /**
         * Randomly selects d jobs from the permutation and tries to reinsert each one
         * of them in the last position, then improves by shifting each job to the left.
         * @param aSol
         * @param d
         * @param first_improvement
         * @return
         */
        inline long randomDestructionConstruction(SchedulingSolution<ProblemInstance> &aSol, const int &d,
            const bool &first_improvement, const unsigned &force_recalculate_obj) {
            std::vector<Job> jobList;
            long num_visited_solutions = 0;
            RandomUtil::shuffle_vector(positions);
            for (int i = 0; i < d; i++) {
                int pos = positions[i];  // obtain random number between 0 and nJobs - 1, without replacement
                jobList.push_back(aSol.jobs[pos]);
                //System.arraycopy(aSol.jobs, pos + 1, aSol.jobs, pos, nJobs - 1 - pos);
                for(int j = pos; j < nJobs - 1; j++ ) {
                    aSol.jobs[j] = aSol.jobs[j + 1];
                }
            }

            for (int i = 0; i < d; i++) {
                aSol.jobs[nJobs - d + i] = jobList[i];
                num_visited_solutions += aSol.improveByShiftingJobToLeft(nJobs - d + i, first_improvement, force_recalculate_obj);
            }
            return num_visited_solutions;
        }

        /**
         * Randomly selects one job to be removed from its current position
         * and inserted in a random position.
         */
        inline long randomInsertion(SchedulingSolution<ProblemInstance> &aSol, const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            int inPos = positions[0];   // random.getRandomPosition(nJobs, "uniform");
            int endPos = positions[nJobs - 1];  // random.getRandomPosition(nJobs, "uniform");

            this->insertion(inPos, endPos, aSol);
            return 1;
        }

        /**
         * Randomly selects one job to be removed from its current position and inserted in a random position,
         * then tries to improve the solution by shifting the modified to job to the left.
         */
        inline long enhancedInsertion(SchedulingSolution<ProblemInstance> &aSol, const bool &first_improvement,
                                      const unsigned &force_recalculate_obj) {
            RandomUtil::shuffle_vector(positions);
            int inPos = positions[0];   // random.getRandomPosition(nJobs, "uniform");
            int endPos = positions[1];  // random.getRandomPosition(nJobs, "uniform");
            for(int x = nJobs - 1; x >= 0 && (inPos == endPos || endPos == nJobs - 1); x--) {
                endPos = positions[x];  // random.getRandomPosition(nJobs, "uniform");
            }

            if (inPos < endPos) {
                int aux = endPos;
                endPos = inPos;
                inPos = aux;
            }
            this->insertion(inPos, endPos, aSol);

            long num_visits = aSol.improveByShiftingJobToLeft(endPos, first_improvement, force_recalculate_obj);

            aSol.setCosts(aSol.calcTotalCosts(nJobs, force_recalculate_obj));
            return num_visits;
        }

        inline void insertion(const int &inPos, const int &endPos, SchedulingSolution<ProblemInstance> &aSol) {
            int dif = endPos - inPos;

            Job inJob = aSol.jobs[inPos];
            if (dif > 0) {
                // System.arraycopy(aSol.jobs, inPos + 1, aSol.jobs, inPos, dif);
                for(int j = inPos; j < inPos + dif; j++ ) {
                    aSol.jobs[j] = aSol.jobs[j + 1];
                }
            } else {
                dif = -dif;
                // System.arraycopy(aSol.jobs, endPos, aSol.jobs, endPos + 1, dif);
                for(int j = endPos + dif - 1; j >= endPos; j-- ) {
                    aSol.jobs[j + 1] = aSol.jobs[j];
                }
            }
            aSol.jobs[endPos] = inJob;
        }

        /**
         * Randomly selects a job between positions 2 and (n-1) and tries to improve the solution by shifting
         * the selected job to the left. Repeats until the randomly-selected job matches its sequence number
         * in the permutation.
         */
        inline long partialImprovement(SchedulingSolution<ProblemInstance> &aSol, const bool &first_improvement,
                                       const unsigned &force_recalculate_obj) {
            int endPos = 0;
            Job auxJob = aSol.jobs[endPos];
            long num_visits = 0;
            RandomUtil::shuffle_vector(positions_nJobs_minus_2);

            int count = 0;
            while (auxJob == aSol.jobs[endPos] && count < positions_nJobs_minus_2.size()) {
                endPos = 1 + positions_nJobs_minus_2[count];  // random.getRandomPosition(nJobs - 2, "uniform");
                auxJob = aSol.jobs[endPos];
                num_visits += aSol.improveByShiftingJobToLeft(endPos, first_improvement, force_recalculate_obj);
                count++;
            }
            aSol.setCosts(aSol.calcTotalCosts(nJobs, force_recalculate_obj));
            return num_visits;
        }

    };
}
}

#endif //FLOWSHOP_SOLVER_LOCALSEARCH_H
