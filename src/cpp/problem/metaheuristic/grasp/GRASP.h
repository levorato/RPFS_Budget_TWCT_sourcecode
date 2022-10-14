//
// Created by mlevorato on 11/27/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_H
#define FLOWSHOP_SOLVER_GRASP_H

#include <vector>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o
#include <glog/logging.h>
#include "../common/Job.h"
#include "../common/SchedulingSolution.h"
#include "../common/Construction.h"
#include "../common/LocalSearch.h"
#include "../common/RandNEHT.h"
#include "../common/Inputs.h"
#include "../common/Outputs.h"
#include "../../../util/include/TimeDateUtil.h"

namespace problem {
namespace pfsp {

    using namespace util;

    template<class ProblemInstance>
    class GRASP {
    public:
        /**General Values*/
        static double TIME_FACTOR;  // 0.03

        Test aTest; // Test to run (including parameters and instance's name)
        Inputs<ProblemInstance> inputs; // Instance inputs

        int nJobs; // #Jobs
        int nMachines; // #Machines
        std::vector<Job> effList; // Jobs sorted by processing time'
        SchedulingSolution<ProblemInstance> nehSol;

        RandNEHT<ProblemInstance> nehtAlg; // Randomized NEH with Taillard's accelerations
        Construction<ProblemInstance> constructor; // Randomized NEH with Taillard's accelerations
        LocalSearch<ProblemInstance> locSearch; // Local Search procedures

        boost::timer::cpu_timer timer;
        long startTime;
        double elapsedTime;
        long maxTime;

        GRASP(const Test &test, const Inputs<ProblemInstance> &inputData, const bool &p_incremental_update_obj = true) :
            aTest(test), inputs(inputData),
            nehSol(inputData.getNumberOfJobs(), inputData.getNumberOfMachines(), inputData), effList(),
            startTime(0), elapsedTime(0), maxTime(0),
            nJobs(inputData.getNumberOfJobs()),
            nMachines(inputData.getNumberOfMachines()), nehtAlg(test, inputData), // Rand NEH with Taillard's
            constructor(test, inputData), // Rand NEH with Taillard's
            locSearch(test, inputData), // Local Search procedures
            timer() {
        }

        // non-trivial construction work
        void init() {
            effList = createEffList();
            timer.start();  // Time measurement
            startTime = TimeDateUtil::calculate_time_spent_ns(timer);
            nehSol = nehtAlg.solve(effList, false); // Computation of the NEH
            // solution
            elapsedTime = TimeDateUtil::calculate_time_spent(timer);
            nehSol.setTime(elapsedTime);
        }

        Outputs<ProblemInstance> run() {
            // 0. Create a base solution and set initial time variables
            double maxTime = aTest.getMaxTime();
            double elapsed = 0.0;
            long iteration = 0;  // current GRASP iteration
            long num_visited_solutions = 0;  // number of solutions examined by GRASP (given all neighborhoods)
            long num_improvements = 0;  // number of global improvements to the initial solution value
            long num_improvements_ls = 0;  // number of improvements in the local search
            stringstream timeResults;  // Stringstream containing the best result found at each moment of time.
            LOG(INFO) << "[GRASP] Starting. Time limit = " << maxTime << " s";

            // 3. ITERATIVE DISRUPTION AND LOCAL SEARCH
            SchedulingSolution<ProblemInstance> bestSol(nehSol);
            double bestSolcost = bestSol.calcTotalCosts(nJobs, aTest.get_force_recalculate_obj());
            if(aTest.get_force_recalculate_obj() >= 3) {  // force recalculation of the obj. function inside GRASP
                bestSolcost = bestSol.forceUpdateCosts(aTest.get_force_recalculate_obj());
                bestSol.costs = bestSolcost;
            }
            elapsed = TimeDateUtil::calculate_time_spent(timer);
            notifyNewValue(bestSol, elapsed, 0, timeResults);  // register first solution

            while( elapsed < maxTime )
            {
                iteration++;
                // 3.0 Greedy Randomized Construction of a new solution
                SchedulingSolution<ProblemInstance> currentSol = nehtAlg.solve(effList, true, bestSol.getCosts(),
                        aTest.is_adaptive_construction());
                if(iteration == 1) {
                    elapsed = TimeDateUtil::calculate_time_spent(timer);
                    notifyNewValue(currentSol, elapsed, -iteration, timeResults);  // register first solution
                }
                // Use the Construction class for deterministic GRASP
                //SchedulingSolution<ProblemInstance> currentSol = constructor.solve(effList);

                // 3.1 Local Search
                //num_visited_solutions += locSearch.globalImprovement(currentSol, incremental_update_obj, num_improvements_ls);
                num_visited_solutions += locSearch.VariableNeighborhoodDescent(currentSol, aTest.get_vnd_permutation(),
                        aTest.get_vnd_size(), aTest.is_first_improvement(), maxTime, elapsed, iteration, timeResults,
                        num_improvements_ls, aTest.get_force_recalculate_obj(), aTest.is_RVND());

                // force recalculation of the obj. function ONLY at each GRASP iteration
                double currentSolcost = currentSol.getCosts();
                if((aTest.get_force_recalculate_obj() == 3) || (aTest.get_force_recalculate_obj() == 4)) {
                    currentSolcost = currentSol.forceUpdateCosts(aTest.get_force_recalculate_obj());
                    currentSol.costs = currentSolcost;
                }

                // 3.2 Update Elapsed time
                elapsed = TimeDateUtil::calculate_time_spent(timer);

                // 3.3 Acceptance Criterion
                if (currentSolcost < bestSolcost) {
                    bestSol = currentSol;
                    bestSolcost = currentSolcost;
                    bestSol.costs = currentSolcost;
                    bestSol.setTime(elapsed);
                    notifyNewValue(bestSol, elapsed, iteration, timeResults);  // register new solution
                    num_improvements++;
                }
            }
            if(elapsed >= maxTime) {
                LOG(INFO) << "[GRASP] Time limit! t=" << elapsed << " s; iter = " << iteration;
                LOG(INFO) << "[GRASP] Number of iterations performed: " << iteration;
            }
            if(aTest.get_force_recalculate_obj() <= 2) {  // force recalculation of the obj. function ONLY in the end
                double cost = bestSol.forceUpdateCosts(aTest.get_force_recalculate_obj());
                bestSol.costs = cost;
                elapsed = TimeDateUtil::calculate_time_spent(timer);
                bestSol.setTime(elapsed);
            }
            std::cout << aTest.getInstanceName() << " "
                               << bestSol.getCosts() << " " << bestSol.getTime() << "\n" /* << " " << nRuns*/;

            // 4. Set output
            return Outputs<ProblemInstance>(nehSol, bestSol, iteration, num_visited_solutions, timeResults.str(),
                    elapsed, num_improvements, nehtAlg.get_alpha_history(), nehtAlg.get_alpha_quality(),
                    nehtAlg.get_alpha_probability());
        }

    private:
        inline void notifyNewValue(SchedulingSolution<ProblemInstance> &sol, const double& timeSpent,
                                   const int& iteration, stringstream& timeResults) {
            timeResults << fixed << setprecision(4) << timeSpent << "," << sol.costs
                        << "," << iteration << "\n";
        }

        std::vector<Job> createEffList() {
            std::vector<Job> array = inputs.getJobs();
            // Sort using the compareTo() method of the Job class (TIE ISSUE #1)
            //std::sort (array.begin(), array.end());
            return array;
        }

        void printSolOnScreen(const SchedulingSolution<ProblemInstance> &aGraspSolution, const bool &isNEHSol) {
            if (isNEHSol)
                std::cout << "\n** NEH SchedulingSolution: ";
            else
                std::cout << "\n** NEH = " << nehSol.getCosts() << "; Seed = "
                                 << aTest.getSeed() << "; OBS = " << aGraspSolution.getCosts()
                                 << "; Time = " << aGraspSolution.getTime();
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_GRASP_H
