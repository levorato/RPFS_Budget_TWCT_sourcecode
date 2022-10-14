//
// Created by mlevorato on 11/27/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_INPUTS_H
#define FLOWSHOP_SOLVER_GRASP_INPUTS_H

#include <vector>
#include "Job.h"

namespace problem {
namespace pfsp {

    template<class ProblemInstance>
    class Inputs {
    private:
        int nJobs;      // #Jobs
        int nMachines;  // #Machines
        std::vector<Job> jobs;     // Array of jobs

    public:
        ProblemInstance I;  // class with problem instance data

        Inputs(int nJobsInProblem, int nMachinesInProblem, const ProblemInstance &_I) : nJobs(nJobsInProblem),
                nMachines(nMachinesInProblem), jobs(), I(_I) {
            for (int i = 0; i < nJobs; i++)
                jobs.push_back(Job(i, nMachines));
        }

        Inputs(const Inputs &rhs) : Inputs(rhs.nJobs, rhs.nMachines, rhs.I) { /* copy construction from rhs*/
            /*
            for (int i = 0; i < rhs.getNumberOfJobs(); i++) {
                jobs[i].setTotalProcessingTime(rhs.jobs[i].getTotalProcessingTime());
            } */
        }

        Inputs &operator=(const Inputs &rhs) {
            // check for self-assignment
            if(&rhs == this)
                return *this;
            nJobs = rhs.nJobs;
            nMachines = rhs.nMachines;
            I = rhs.I;
            /*
            for (int i = 0; i < rhs.getNumberOfJobs(); i++) {
                jobs[i].setTotalProcessingTime(rhs.jobs[i].getTotalProcessingTime());
            } */
            return *this;
        }

        void initialize() {

        }

        int getNumberOfJobs() const {
            return nJobs;
        }

        int getNumberOfMachines() const {
            return nMachines;
        }

        const std::vector<Job> getJobs() const {
            return jobs;
        }

        Job getJob(const int &job_number) {
            return jobs[job_number];
        }

        void setJob(const int &job_number, const Job &job) {
            jobs[job_number] = job;
        }

        double totalTime() const {
            return I.totalTime();
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_INPUTS_H
