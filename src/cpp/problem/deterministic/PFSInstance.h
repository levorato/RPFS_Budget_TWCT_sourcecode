/*
 * PFSInstance.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_INCLUDE_PFSINSTANCE_H_
#define SRC_CPP_PROBLEM_INCLUDE_PFSINSTANCE_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include "../metaheuristic/common/Job.h"

using namespace std;
using namespace problem::pfsp;

class PFSInstance {
public:
    // name which identified the instance
    string name;
    // number of jobs
    unsigned int n;
    // number of machines
    unsigned int m;
    // processing time matrix => indices start at 1 !!!
    boost::numeric::ublas::matrix<double> p_time;
    // best known lower bound
    double lower_bound;
    // best known upper bound
    double upper_bound;

    PFSInstance(const string& instance_name, boost::numeric::ublas::matrix<double> p) : name(instance_name), p_time(p) {
		// p_time matrix is m x n
		m = p_time.size1() - 1;
		n = p_time.size2() - 1;
		upper_bound = std::numeric_limits<double>::max();
		lower_bound = 0.0;
	}
	virtual ~PFSInstance();

	const boost::numeric::ublas::matrix<double>& getTimeMatrix() const {
		return p_time;
	}

	const unsigned int& getNumberOfMachines() {
		return m;
	}

	const unsigned int& getNumberOfJobs() {
		return n;
	}

    PFSInstance generate_reverse_instance() {
		boost::numeric::ublas::matrix<double> p_reversed(p_time.size1(), p_time.size2());
        for(unsigned i = 1, k = m; i <= m; i++, k--) {  // iterate over all columns of p
            for(unsigned j = 1; j <= n; j++) {
                p_reversed(k, j) = p_time(i, j);
            }
        }
        return PFSInstance(name + "_reversed", p_reversed);
	}

    /*******************************************************************************
     * PUBLIC METHOD calcTotalCosts()
    ******************************************************************************/
    double calcTotalCosts(std::vector<Job> &pi, int nUsedJobs, unsigned recalculate_obj = 0) const {
        // indices from p_time start at one, not zero!
        // nUsedJobs = # of jobs in the partially filled solution
        std::vector< std::vector<double> > tcosts(nUsedJobs, std::vector<double>(m, 0));
        for (int column = 0; column < m; column++)
            for (int row = 0; row < nUsedJobs; row++) {
                if (column == 0 && row == 0)
                    tcosts[0][0] = p_time(0 + 1, pi[0].id);
                else if (column == 0)
                    tcosts[row][0] = tcosts[row - 1][0]
                                     + p_time(0 + 1, pi[row].id);
                else if (row == 0)
                    tcosts[0][column] = tcosts[0][column - 1]
                                        + p_time(column + 1, pi[0].id);
                else {
                    double max = std::max(tcosts[row - 1][column],
                                          tcosts[row][column - 1]);
                    tcosts[row][column] = max
                                          + p_time(column + 1, pi[row].id);
                }
            }
        return tcosts[nUsedJobs - 1][m - 1];
    }

    double calculate_objective(const std::vector<int> &pi) const {
        std::vector<Job> perm;
        for (int j : pi)
            perm.push_back(Job(j - 1, 0));
        return calcTotalCosts(perm, perm.size(), 0);
    }

    /*******************************************************************************
     * PUBLIC METHOD improveByShiftingJobToLeft() This method implements
     * Taillard's accelerations where k is the position of the job on the right
     * extreme.
     *
     * This method also updates the solution cost (makespan) if k == nJobs -1
     * Note: pi and makespan are in_out parameters.
     * Returns true if the resulting makespan has changed.
     ******************************************************************************/
    bool improveByShiftingJobToLeft(std::vector<Job> &pi, double &makespan, int k, long &num_visits,
                                    bool first_improvement, unsigned force_recalculate_obj) {
        int bestPosition = k;
        double minMakespan = std::numeric_limits<double>::max();
        double newMakespan = std::numeric_limits<double>::max();

        double maxSum = 0;
        double newSum = 0;

        // Calculate eMatrix
        std::vector< std::vector<double> > eMatrix = calcEMatrix(pi, k);

        // Calculate qMatrix
        std::vector< std::vector<double> > qMatrix = calcQMatrix(pi, k);

        // Calculate fMatrix
        std::vector< std::vector<double> > fMatrix = calcFMatrix(pi, k, eMatrix);

        // Calculate bestPosition (0...k) and minMakespan (mVector)
        for (int i = k; i >= 0; i--) {
            num_visits++;
            maxSum = 0;
            for (int j = 0; j < m; j++) {
                newSum = fMatrix[i][j] + qMatrix[i][j];
                if (newSum > maxSum)
                    maxSum = newSum;
            }
            newMakespan = maxSum;
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
        if (k == n - 1) {
            // this->setCosts(minMakespan);
            makespan = minMakespan;
            return true;
        }
        return false;
    }

    double totalTime() const {
        double T = 0;
        for(unsigned i=1; i<=m; i++)
            for(unsigned j=1; j<=n; j++)
                T += p_time(i, j);
        return T;
    }

private:
    /*******************************************************************************
     * PRIVATE METHOD calcQMatrix()
     ******************************************************************************/
    std::vector< std::vector<double> > calcQMatrix(const std::vector<Job> &pi, const int &k) {
        std::vector< std::vector<double> > q(k + 1, std::vector<double>(m, 0));
        // indices from p_time start at one, not zero!
        for (int i = k; i >= 0; i--) {
            for (int j = m - 1; j >= 0; j--) {
                if (i == k)
                    q[k][j] = 0; // dummy file to make possible fMatrix +
                    // qMatrix
                else if (i == k - 1 && j == m - 1)
                    q[k - 1][m - 1] = p_time(m - 1 + 1, pi[k - 1].id);
                else if (j == m - 1)
                    q[i][m - 1] = q[i + 1][m - 1]
                                          + p_time(m - 1 + 1, pi[i].id);
                else if (i == k - 1)
                    q[k - 1][j] = q[k - 1][j + 1]
                                  + p_time(j + 1, pi[k - 1].id);
                else {
                    int max = std::max(q[i + 1][j], q[i][j + 1]);
                    q[i][j] = max + p_time(j + 1, pi[i].id);
                }
            }
        }
        return q;
    }

    /*******************************************************************************
     * PRIVATE METHOD calcFMatrix()
     ******************************************************************************/
    std::vector< std::vector<double> > calcFMatrix(const std::vector<Job> &pi, const int &k,
            const std::vector< std::vector<double> > &e) {
        std::vector< std::vector<double> > f(k + 1, std::vector<double>(m, 0));
        // indices from p_time start at one, not zero!
        for (int i = 0; i <= k; i++) {
            for (int j = 0; j < m; j++) {
                if (i == 0 && j == 0)
                    f[0][0] = p_time(0 + 1, pi[k].id);
                else if (j == 0)
                    f[i][0] = e[i - 1][0] + p_time(0 + 1, pi[k].id);
                else if (i == 0)
                    f[0][j] = f[0][j - 1] + p_time(j + 1, pi[k].id);
                else {
                    int max = std::max(e[i - 1][j], f[i][j - 1]);
                    f[i][j] = max + p_time(j + 1, pi[k].id);
                }
            }
        }
        return f;
    }

    /*******************************************************************************
     * PRIVATE METHOD calcEMatrix()
     ******************************************************************************/
    std::vector< std::vector<double> > calcEMatrix(const std::vector<Job> &pi, const int &k) {
        std::vector< std::vector<double> > e(k, std::vector<double>(m, 0));
        // indices from p_time start at one, not zero!
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < m; j++) {
                if (i == 0 && j == 0)
                    e[0][0] = p_time(0 + 1, pi[0].id);
                else if (j == 0)
                    e[i][0] = e[i - 1][0] + p_time(0 + 1, pi[i].id);
                else if (i == 0)
                    e[0][j] = e[0][j - 1] + p_time(j + 1, pi[0].id);
                else {
                    int max = std::max(e[i - 1][j], e[i][j - 1]);
                    e[i][j] = max + p_time(j + 1, pi[i].id);
                }
            }
        }
        return e;
    }
};

#endif /* SRC_CPP_PROBLEM_INCLUDE_PFSINSTANCE_H_ */
