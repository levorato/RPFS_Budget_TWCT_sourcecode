//
// Created by Mario Costa Levorato Junior on 2019-02-01.
//

#ifndef FLOWSHOP_SOLVER_GRASP_CONSTRUCTION_H
#define FLOWSHOP_SOLVER_GRASP_CONSTRUCTION_H

#include <vector>
#include <list>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "SchedulingSolution.h"
#include "Job.h"
#include "Inputs.h"
#include "Randomness.h"

namespace problem {
namespace pfsp {

    template<class ProblemInstance>
    class Construction {
    private:
        Test aTest;     // Test to run (including parameters and instance's name)
        Inputs<ProblemInstance> inputs;  // Instance inputs
        int nJobs;      // #Jobs
        int nMachines;  // #Machines
        std::vector<int> positions; // Array of randomly selected positions
        Job nextJob;
        // The value of the alpha parameter is found by an on-line adaptive mechanism
        double ALPHA;
        typedef boost::random::uniform_int_distribution<int> UNIFORM_INT_DISTRIBUTION;

    public:
        Construction(const Test &test, const Inputs<ProblemInstance> &inputData, const double &p_alpha = -1.0) :
                                                    aTest(test), inputs(inputData),
                                                    nJobs(inputData.getNumberOfJobs()),
                                                    nMachines(inputData.getNumberOfMachines()),
                                                    positions(inputData.getNumberOfJobs(), 0),
                                                    nextJob(0, 0), ALPHA(p_alpha) {
        }

        SchedulingSolution<ProblemInstance> solve(std::vector<Job> &effList) {
            // 0. Define a new solution
            SchedulingSolution<ProblemInstance> currentSol(nJobs, nMachines, inputs);
            Randomness<ProblemInstance> rand(aTest, inputs);
            std::list<Job> s;
            std::copy( effList.begin(), effList.end(), std::back_inserter( s ) );  // linked list of Jobs
            std::vector<double> job_cost(nJobs + 1, 0);
            std::vector< std::list<Job>::iterator > job_ptr(nJobs + 1, s.end());

            for (int i = 0; i < nJobs; i++) {
                std::vector<Job> RCL;
                if(ALPHA < 1.0) {
                    double min = std::numeric_limits<double>::max(), max = std::numeric_limits<double>::min();
                    for (std::list<Job>::iterator it = s.begin(); it != s.end(); ++it) {  // For each Job in the list s
                        // Add nextJob j to the end of currentSol (partial solution)
                        Job j = *it;
                        currentSol.jobs[i] = j;
                        double cost = currentSol.calcTotalCosts(i + 1);
                        job_cost[j.getId()] = cost;
                        job_ptr[j.getId()] = it;
                        if (cost < min)
                            min = cost;
                        if (cost > max)
                            max = cost;
                    }
                    double th = min + ALPHA * (double) (max - min);
                    std::vector<Job> RCL;
                    for (Job j : s) {
                        if (job_cost[j.getId()] <= th) {
                            RCL.push_back(j);
                        }
                    }
                } else {
                    for (std::list<Job>::iterator it = s.begin(); it != s.end(); ++it) {  // For each Job in the list s
                        Job j = *it;
                        job_ptr[j.getId()] = it;
                        RCL.push_back(j);
                    }
                }
                int list_size = RCL.size();
                assert(list_size > 0);
                int index = rand.getRandomPosition(list_size, "u");
                nextJob = RCL[index];
                currentSol.jobs[i] = nextJob;
                s.erase(job_ptr[nextJob.getId()]);
                //std::fill_n(job_cost.begin(), job_cost.size(), 0);  // reset the job_cost array
                //std::fill_n(job_ptr.begin(), job_ptr.size(), s.end());
            }

            return currentSol;
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_CONSTRUCTION_H
