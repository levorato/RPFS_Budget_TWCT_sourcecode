//
// Created by mlevorato on 12/5/19.
//

#ifndef FLOWSHOP_SOLVER_PFSPBUDGETSCENARIO_H
#define FLOWSHOP_SOLVER_PFSPBUDGETSCENARIO_H

#include "../../GlobalTypes.h"
#include "../PFSPScenario.h"
#include "../RobPFSInstance.h"
#include <vector>
#include <list>
#include <boost/numeric/ublas/matrix.hpp>
#include <chrono>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> //include all types plus i/o

namespace robust {

    using namespace std;
    using namespace boost;
    //using namespace parameters;

    typedef std::list<int> JobList;

    /**
     * Jobs that have uncertain processing times can be represented in a scenario with two strings of numbers,
     * each consisting of a permutation of n jobs for M1 and M2, respectively.
     * The first floor(T1) number in the first string and the first floor(T2) number in the second string denote
     * the jobs that have uncertain processing times to reach their worst-case values in M1 and M2, respectively.
     * Meanwhile, the number in the floor(T1)+1 position of the first string and the number in the floor(T2)+1
     * position of the second string represent the jobs that have uncertain processing times to be changed by
     * (T1 - floor(T1))*p_hat(1, t_1) in M1 and (T1 - floor(T2))*p_hat(2, t_2) in M2, respectively.
     * The processing times of the remaining jobs in the two strings are certain (ie, p_ij = p_bar_ij).
     */
    class PFSPBudgetScenario : public PFSPScenario {
    public:
        /**
         * Constructor for Ying's metaheuristics.
         * @param p_T
         * @param p_dev
         */
        PFSPBudgetScenario(const std::vector<Time> &p_T, const std::vector<Permutation> &p_dev) :
                PFSPScenario(), T(p_T), dev(p_dev), dev_multiplier(), p_time_calculated(false),
                use_ceil(false) {
        }

        /**
         * Constructor for GRASP and ILS metaheuristics, whenever the worst-case scenario is calculated by
         * one of these metaheuristics.
         * @param p_T
         * @param p_dev
         * @param p_dev_multiplier
         */
        PFSPBudgetScenario(const std::vector<Time> &p_T, const std::vector<Permutation> &p_dev,
                           const std::vector< std::vector<double> > &p_dev_multiplier) :
                PFSPScenario(), T(p_T), dev(p_dev), dev_multiplier(p_dev_multiplier), p_time_calculated(false),
                use_ceil(false) {
        }

        PFSPBudgetScenario(const PFSPBudgetScenario &rhs) : PFSPScenario(rhs.p_time, rhs.value), T(rhs.T), dev(rhs.dev),
                                                            dev_multiplier(rhs.dev_multiplier),
                                                            p_time_calculated(rhs.p_time_calculated),
                                                            use_ceil(rhs.use_ceil) {

        }

        // The indices of both vectors below start at 0 !
        std::vector<double> T;        // vector of Budget parameters T (in % of the number of jobs n)
        std::vector<Permutation> dev; // vector of list of jobs that deviated to worst-case values, for each machine M_i
        // For each Machine M_i, stores the multiplier of how much each job deviates its proc. time
        std::vector<std::vector<double>> dev_multiplier;
        // The 1-index-based matrix below is only filled after the invocation of getProcessingTimeMatrix()
        // boost::numeric::ublas::matrix<double> p_time;  -> already defined in superclass PFSPScenario
        bool p_time_calculated;
        // objective value of the optimal solution for this scenario
        // double value;   -> already defined in superclass PFSPScenario
        bool use_ceil;  // use ceil for budget parameter calculation

        // TODO TROCAR O VECTOR DE JOBS PELA LISTA ENCADEADA
        const boost::numeric::ublas::matrix<Time>& getProcessingTimeMatrix(const RobPFSInstance *instance) {
            if(! p_time_calculated) {
                const boost::numeric::ublas::matrix<double>& p_bar = instance->get_p_bar();
                const boost::numeric::ublas::matrix<double>& p_hat = instance->get_p_hat();
                p_time = boost::numeric::ublas::zero_matrix<Time>(instance->getNumberOfMachines() + 1,
                        instance->getNumberOfJobs() + 1);
                if(dev_multiplier.size() == 0) {  // Ying's way of calculating the processing time deviations
                    for (unsigned i = 1; i <= instance->getNumberOfMachines(); i++) {  // For each machine i
                        // adjust the value of the budget paramter T_i as a percentage of the number of jobs n
                        double T_i = (instance->getNumberOfJobs() * T[i - 1]) / 100;
                        if (use_ceil) T_i = ceil(T_i);
                        unsigned ki = (int) floor(T_i);
                        int job = 0;
                        // The first ki = floor(T[i]) number in the permutation dev[i] denote the jobs that reach their
                        // worst-case values in machine M_i : jobs in [ dev[i][1] ... dev[i][ki] ]
                        for (unsigned k = 1; k <= ki; k++) {  // For each position k in permutation dev[i]
                            job = dev[i - 1][k - 1];  // zero-based index
                            p_time(i, job) = p_bar(i, job) + p_hat(i, job);
                        }
                        // The number in the floor(T[i])+1 position of the permutation dev[i] represent the job that
                        // will have uncertain processing time to be changed by (T[i] - floor(T[i]))*p_hat(i, t_i) in M_i
                        if (ki + 1 <= instance->getNumberOfJobs()) {
                            job = dev[i - 1][ki]; // zero-based index
                            p_time(i, job) = p_bar(i, job) + (T_i - floor(T_i)) * p_hat(i, job);
                        }
                        // The processing times of the remaining jobs in the permutation dev[i] are certain (ie, p_ij = p_bar_ij)
                        for (unsigned k = ki + 2;
                             k <= instance->getNumberOfJobs(); k++) {  // For each position k in permutation dev[i]
                            job = dev[i - 1][k - 1];  // zero-based index
                            p_time(i, job) = p_bar(i, job);
                        }
                    }
                } else {
                    for (unsigned i = 1; i <= instance->getNumberOfMachines(); i++) {  // For each machine i
                        // adjust the value of the budget paramter T_i as a percentage of the number of jobs n
                        for (unsigned k = 0; k <= instance->getNumberOfJobs(); k++) {
                            p_time(i, k) = p_bar(i, k);
                        }
                        for (unsigned k = 1; k <= dev_multiplier[i - 1].size(); k++) {  // For each position k in permutation dev[i]
                            int job = dev[i - 1][k - 1];  // zero-based index
                            double deviation = dev_multiplier[i - 1][k - 1];
                            p_time(i, job) = p_bar(i, job) + deviation * p_hat(i, job);
                        }
                    }
                }
            }
            return p_time;
        }

        /**
         * Updates the processing time matrix after a swap operation between the jobs in positions posA and posB in
         * the budget sequence (scenario.dev) on machine seq_num. Used by Ying's metaheuristics (SA and IG) only.
         * @param instance
         * @return
         */
        void updateProcessingTimeMatrix(const RobPFSInstance *instance,
                const unsigned &seq_num, const unsigned &posA, const unsigned &posB) {
            // In which machine did the swap operation occur? => Machine seq_num in {0, ..., m - 1}
            unsigned i = seq_num + 1;
            double T_i = (instance->n * T[i - 1]) / 100;
            if(use_ceil)  T_i = ceil(T_i);
            // The job in posA has uncertain processing times ==> posA < (floor(T[i])+1) - 1
            int job = dev[i - 1][posA];  // zero-based index
            unsigned ki = (int)floor(T_i);
            if(posA <= ki - 1) {  // The job in posA is *uncertain* and will reach worst-case values in machine M_i
                p_time(i, job) = instance->p_bar(i, job) + instance->p_hat(i, job);
            } else {  // job with uncertain processing time to be changed by (T[i] - floor(T[i]))*p_hat(i, t_i) in M_i
                p_time(i, job) = instance->p_bar(i, job) + (T_i - floor(T_i)) * instance->p_hat(i, job);
            }
            // The processing time of the job in posB in the permutation dev[i] is *certain* (ie, p_ij = p_bar_ij)
            job = dev[i - 1][posB];  // zero-based index
            p_time(i, job) = instance->p_bar(i, job);
        }

        friend std::ostream& operator << (std::ostream& os, const PFSPBudgetScenario& scenario)
        {
            os << "[";
            for(unsigned i = 1; i <= scenario.dev.size(); i++) {  // For each machine i
                os << "dev[" << i << "] = ";
                for(unsigned int j = 0; j < scenario.dev[i - 1].size(); j++) {
                    if(scenario.dev_multiplier.size() > 0) {
                        if(scenario.dev_multiplier[i - 1].size() > 0) {
                            os << "(" << scenario.dev[i - 1][j];
                            os << "; " << scenario.dev_multiplier[i - 1][j] << ") ";
                        } else {
                            os << scenario.dev[i - 1][j] << " ";
                        }
                    } else {
                        os << scenario.dev[i - 1][j] << " ";
                    }
                }
            }
            os << " ]";
            return os;
        }

    };
}

#endif //FLOWSHOP_SOLVER_PFSPBUDGETSCENARIO_H
