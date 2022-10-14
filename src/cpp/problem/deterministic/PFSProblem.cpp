/*
 * PFSProblem.cpp
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#include "PFSProblem.h"
#include "../../util/NumericUtil.h"

#include <algorithm>
#include <glog/logging.h>

using namespace boost;

PFSProblem::PFSProblem() {
	// TODO Auto-generated constructor stub

}

PFSProblem::~PFSProblem() {
	// TODO Auto-generated destructor stub
}

boost::numeric::ublas::matrix<double> PFSProblem::makespan(Permutation my_seq, PFSInstance instance) {
    const boost::numeric::ublas::matrix<double>& p_ij = instance.getTimeMatrix();
	unsigned int nb_machines = instance.getNumberOfMachines();
	return PFSProblem::makespan(my_seq, nb_machines, p_ij);
}


boost::numeric::ublas::matrix<double> PFSProblem::makespan(Permutation my_seq, int nb_machines, const boost::numeric::ublas::matrix<double>& p_ij) {
    boost::numeric::ublas::matrix<double> c_ij = boost::numeric::ublas::zero_matrix<double>(nb_machines + 1, my_seq.size() + 2);
    auto less_eps = [&](const double& a, const double& b) {
        return NumericUtil::TestsetIsLT(a, b);
    };
    if(my_seq.size() > 0) {
        c_ij(1, 1) = p_ij(1, my_seq[0]);
        for (int i = 2; i <= nb_machines; i++) {
            c_ij(i, 1) = c_ij(i - 1, 1) + p_ij(i, my_seq[0]);
        }
        for (int j = 2; j <= my_seq.size(); j++) {
            c_ij(1, j) = c_ij(1, j - 1) + p_ij(1, my_seq[j - 1]);
            for (int i = 2; i <= nb_machines; i++) {
                c_ij(i, j) = std::max(c_ij(i, j - 1), c_ij(i - 1, j), less_eps) + p_ij(i, my_seq[j - 1]);
            }
        }
    }
    //print(c_ij)
    return c_ij;
}

void PFSProblem::makespan(boost::numeric::ublas::matrix<double>& c_ij, Permutation my_seq, int nb_machines,
                          const boost::numeric::ublas::matrix<double>& p_ij) {
    auto less_eps = [&](const double& a, const double& b) {
        return NumericUtil::TestsetIsLT(a, b);
    };
    if(my_seq.size() > 0) {
        c_ij(1, 1) = p_ij(1, my_seq[0]);
        for (int i = 2; i <= nb_machines; i++) {
            c_ij(i, 1) = c_ij(i - 1, 1) + p_ij(i, my_seq[0]);
        }
        for (int j = 2; j <= my_seq.size(); j++) {
            c_ij(1, j) = c_ij(1, j - 1) + p_ij(1, my_seq[j - 1]);
            for (int i = 2; i <= nb_machines; i++) {
                c_ij(i, j) = std::max(c_ij(i, j - 1), c_ij(i - 1, j), less_eps) + p_ij(i, my_seq[j - 1]);
            }
        }
    }
    //print(c_ij)
    //return c_ij;
}

void PFSProblem::incremental_makespan(boost::numeric::ublas::matrix<double>& c_ij, Permutation my_seq, int nb_machines,
                                      const boost::numeric::ublas::matrix<double>& p_ij, const int& j) {
    auto less_eps = [&](const double& a, const double& b) {
        return NumericUtil::TestsetIsLT(a, b);
    };
    if(my_seq.size() > 0) {
        if(j == 1) {
            c_ij(1, 1) = p_ij(1, my_seq[0]);
            for (int i = 2; i <= nb_machines; i++) {
                c_ij(i, 1) = c_ij(i - 1, 1) + p_ij(i, my_seq[0]);
            }
        } else {
            c_ij(1, j) = c_ij(1, j - 1) + p_ij(1, my_seq[j - 1]);
            for (int i = 2; i <= nb_machines; i++) {
                c_ij(i, j) = std::max(c_ij(i, j - 1), c_ij(i - 1, j), less_eps) + p_ij(i, my_seq[j - 1]);
            }
        }
    }
    //print(c_ij)
    //return c_ij;
}

double PFSProblem::calculateCmax(Permutation my_seq, PFSInstance instance) {
    // TODO remover validacao do calculo da FO
    //LOG(INFO) << "[PFSProblem] Cmax_1";
    double cmax_1 = computeMakespan_alternative(my_seq, instance);
    LOG(INFO) << "[PFSProblem] Cmax_1 = " << cmax_1;
    double cmax_2 = makespan(my_seq, instance)(instance.m, my_seq.size());
    LOG(INFO) << "[PFSProblem] Cmax_2 = " << cmax_2;
    if (cmax_1 != cmax_2) {
        LOG(ERROR) << "[PFSProblem] Makespan calculation DOES NOT MATCH!";
    }
    return cmax_2;
}

double PFSProblem::computeMakespan_alternative(Permutation my_seq, PFSInstance instance)
{
    int numberOfJobs = my_seq.size();
    int numberOfMachines = instance.m;
    const boost::numeric::ublas::matrix<double>& p_ij = instance.getTimeMatrix();
    std::vector<std::vector<double>> completeTaskTime(numberOfJobs + 1, std::vector<double>(numberOfMachines + 1, 0));
    auto less_eps = [&](const double& a, const double& b) {
        return NumericUtil::TestsetIsLT(a, b);
    };
    for (int i = 1; i <= numberOfJobs; i++)
    {
        for (int j = 1; j <= numberOfMachines; j++)
        {

            double first = i > 0 ? completeTaskTime[i - 1][j] : 0;
            double second = j > 0 ? completeTaskTime[i][j - 1] : 0;

            double greaterTaskTime = max(first, second, less_eps);
            completeTaskTime[i][j] = greaterTaskTime + p_ij(j, my_seq[i - 1]); // + resultSchedule[i].getTimeJob(j);

        }
    }
    return completeTaskTime.back().back();
}

double PFSProblem::calculateCmax(Permutation my_seq, int m, int n, const boost::numeric::ublas::matrix<double>& p_ij)
{
    return makespan(my_seq, m, p_ij)(m, n);
}

double PFSProblem::calculateCmax(const std::vector<Job> &pi, int m, int n, const boost::numeric::ublas::matrix<double>& p_ij)
{
    auto less_eps = [&](const double& a, const double& b) {
        return NumericUtil::TestsetIsLT(a, b);
    };
    boost::numeric::ublas::matrix<double> c_ij = boost::numeric::ublas::zero_matrix<double>(m + 1, pi.size() + 2);
    for(int j = 1; j < pi.size() + 1; j++) {
        c_ij(0, j) = c_ij(0, j - 1) + p_ij(1, pi[j - 1].id);
    }
    for(int i = 1; i <= m; i++) {
        for(int j = 1; j <= pi.size(); j++) {
            c_ij(i, j) = std::max(c_ij(i - 1, j), c_ij(i, j - 1), less_eps) + p_ij(i, pi[j - 1].id);
        }
    }
    return c_ij(m, n);
}

double PFSProblem::calculate_wct(const Permutation &pi, const int &m, const int &n,
        const boost::numeric::ublas::matrix<double>& p_ij, const boost::numeric::ublas::vector<double>& w)
{
    boost::numeric::ublas::matrix<double> C = makespan(pi, m, p_ij);
    double wct = 0.0;
    for(int j = 1; j <= n; j++) {  // C[m, j] is the final completion time of job in position j
        int k = pi[j - 1];
        wct += C(m, j) * w(k);  // C[m, j] * w[k]
    }
    return wct;
}
