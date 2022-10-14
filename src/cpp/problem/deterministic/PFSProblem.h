/*
 * PFSProblem.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_INCLUDE_PFSPROBLEM_H_
#define SRC_CPP_PROBLEM_INCLUDE_PFSPROBLEM_H_

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include "../GlobalTypes.h"
#include "PFSInstance.h"
#include "../PFSSolution.h"
#include "../metaheuristic/common/Job.h"

using namespace std;
using namespace boost;
using namespace problem::pfsp;

class PFSProblem {
public:
	PFSProblem();
	virtual ~PFSProblem();

	static boost::numeric::ublas::matrix<double> makespan(Permutation my_seq, PFSInstance instance);
	static boost::numeric::ublas::matrix<double> makespan(Permutation my_seq, int nb_machines,
	        const boost::numeric::ublas::matrix<double>& p_ij);
    static void makespan(boost::numeric::ublas::matrix<double>& c_ij, Permutation my_seq, int nb_machines,
                         const boost::numeric::ublas::matrix<double>& p_ij);
    static void incremental_makespan(boost::numeric::ublas::matrix<double>& c_ij, Permutation my_seq, int nb_machines,
                                                 const boost::numeric::ublas::matrix<double>& p_ij, const int& j);

    static double calculateCmax(Permutation my_seq, PFSInstance instance);
    static double computeMakespan_alternative(Permutation my_seq, PFSInstance instance);
    static double calculateCmax(Permutation my_seq, int m, int n, const boost::numeric::ublas::matrix<double>& p_ij);
    static double calculateCmax(const std::vector<Job> &pi, int m, int n, const boost::numeric::ublas::matrix<double>& p_ij);

    /**
     * Calculate the weighted sum of job completion times (WCT).
     * @param pi the permutation of jobs
     * @param m number of machines
     * @param n number of jobs
     * @param p_ij processing times of job j, on each machine i
     * @param w job weights, for each job j
     * @return
     */
    static double calculate_wct(const Permutation &pi, const int &m, const int &n,
            const boost::numeric::ublas::matrix<double>& p_ij, const boost::numeric::ublas::vector<double>& w);
};

#endif /* SRC_CPP_PROBLEM_INCLUDE_PFSPROBLEM_H_ */
