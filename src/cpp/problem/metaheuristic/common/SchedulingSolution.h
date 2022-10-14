//
// Created by Mario Costa Levorato Junior on 2019-02-01.
//

#ifndef FLOWSHOP_SOLVER_GRASP_SOLUTION_H
#define FLOWSHOP_SOLVER_GRASP_SOLUTION_H

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>

#include "Job.h"
#include "ElapsedTime.h"
#include "Inputs.h"
#include "../../GlobalTypes.h"

namespace problem {
namespace pfsp {
    /**
    * This class represents a solution (sequence of jobs) of the flow shop. It also
    * includes some important methods to be applied over a solution, like the
    * Taillard accelerations to improve a given solution.
    */
    template<class ProblemInstance>
    class SchedulingSolution {
    public:
        double costs; // solution costs = end of processing time for all jobs
        int nJobs; // number of jobs in the problem
        std::vector<problem::pfsp::Job> jobs; // array of jobs in this solution
        int nMachines; // number of machines in the problem
        double time; // elapsed computational time (in seconds)
        Inputs<ProblemInstance> inputs;

        SchedulingSolution(int nJobsInProblem, int nMachinesInProblem, const Inputs<ProblemInstance> &_inputs) :
                costs(0), nJobs(nJobsInProblem), nMachines(nMachinesInProblem), time(0), jobs(),
                inputs(_inputs) {
            for (int i = 0; i < nJobs; i++)
                jobs.push_back(problem::pfsp::Job(i, nMachines));
        }

        SchedulingSolution(int nJobsInProblem, int nMachinesInProblem, const std::vector<problem::pfsp::Job>& p_jobs, const double &p_costs,
                           const double &p_time, const Inputs<ProblemInstance> &_inputs) :
            nJobs(nJobsInProblem), nMachines(nMachinesInProblem), costs(p_costs), time(p_time),
            jobs(p_jobs), inputs(_inputs)
        {

        }

        SchedulingSolution &operator=(const SchedulingSolution &rhs) {
            // check for self-assignment
            if(&rhs == this)
                return *this;
            nJobs = rhs.nJobs;
            nMachines = rhs.nMachines;
            // copy the jobs vector
            jobs = rhs.jobs;

            setCosts(rhs.getCosts());
            setTime(rhs.getTime());
            return *this;
        }

        inline void setCosts(const double &c) {
            costs = c;
        }


        inline void setTime(const double &t) {
            time = t;
        }

        inline void setJob(const int &pos, const problem::pfsp::Job &aJob) {
            jobs[pos] = aJob;
        }

        inline double getCosts() const {
            return costs;
        }

        const inline std::vector<problem::pfsp::Job> getJobs() const {
            return jobs;
        }

        inline double getTime() const {
            return time;
        }

        inline int getNJobs() const {
            return nJobs;
        }

        inline int getNMachines() const {
            return nMachines;
        }

        double calcTotalCosts(int nUsedJobs, const unsigned &force_recalculate_obj) {
            return inputs.I.calcTotalCosts(jobs, nUsedJobs, force_recalculate_obj);
        }

        /**
         * Force the full recalculation of the objective function.
         */
        double forceUpdateCosts(unsigned recalculate_obj) {
            return inputs.I.calcTotalCosts(jobs, nJobs, 0);
        }

        /*******************************************************************************
         * improveByShiftingJobToLeft() This method implements
         * Taillard accelerations where k is the position of the job on the right
         * extreme.
         *
         * This method also updates the solution cost (makespan) if k == nJobs -1
         * If first_improvement is true, interrupts the local search on the first
         * improvement found.
         * Returns the number of visited solutions.
         ******************************************************************************/
        long improveByShiftingJobToLeft(int k, bool first_improvement, const unsigned &force_recalculate_obj) {
            double makespan = this->getCosts();
            long num_visits = 0;
            if(inputs.I.improveByShiftingJobToLeft(jobs, makespan, k, num_visits, first_improvement, force_recalculate_obj)) {
                this->setCosts(makespan);
            }
            return num_visits;
        }

        std::string toString(const bool &printDetails) {
            std::stringstream ss;
            ss << "\r\n";
            ss << "Sol costs: " << this->getCosts() << "\r\n";
            double time = this->getTime();
            int timeInt = (int) round(time);
            ss << "Sol time: " << pfsp::ElapsedTime::calcHMS(timeInt) << " (" << time << " sec.)";
            ss << "\r\n";
            if (printDetails) {
                ss << "List of jobs: \r\n";
                for (int i = 0; i < jobs.size(); i++)
                    ss << jobs[i].getId() << "\r\n";
            }
            return ss.str();
        }

        std::ostream &operator<<(std::ostream &os) {
            os << "time spent: " << time << "\nobjective value: " << costs << std::endl;
            os << "[ ";
            for (int i = 0; i < jobs.size(); i++)
                os << jobs[i].getId();
            os << " ]\r\n";
            return os;
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_GRASP_SOLUTION_H
