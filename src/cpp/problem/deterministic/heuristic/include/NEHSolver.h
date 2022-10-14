/*
 * NEHSolver.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_NEHSOLVER_H_
#define SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_NEHSOLVER_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <glog/logging.h>
#include <boost/foreach.hpp>

#include "../../PFSInstance.h"
#include "../../../PFSSolution.h"
#include "../../PFSProblem.h"
#include "../../../../util/include/CollectionUtil.h"
#include "../../../../util/include/SchedulingUtil.h"
#include "../../../GlobalTypes.h"

using namespace std;
using namespace boost;

namespace problem {
    namespace pfsp {
        namespace heuristic {

            using namespace util;
            using namespace problem::common;

            class NEHSolver {
            public:
                NEHSolver();

                virtual ~NEHSolver();

                double sum_processing_time_for_job(unsigned int index_job, const boost::numeric::ublas::matrix<double> &data,
                                                   unsigned int nb_machines) {
                    double sum_p = 0.0;
                    for (int i = 1; i <= nb_machines; i++) {
                        sum_p += data(i, index_job);
                    }
                    return sum_p;
                }

                std::vector<int>
                order_neh(const boost::numeric::ublas::matrix<double> &data, unsigned int nb_machines, unsigned int nb_jobs) {
                    std::vector<int> my_seq;
                    for (int j = 0; j < nb_jobs; j++) {
                        my_seq.push_back(j);
                    }
                    // create a vector of sum_processing_time, to be used in the sort operation of my_seq
                    std::vector<double> processing_time_vector_for_job = build_sum_processing_time_vector(data,
                                                                                                          nb_machines,
                                                                                                          nb_jobs);
                    std::sort(my_seq.begin(), my_seq.end(),
                              util::vector_comparator(&processing_time_vector_for_job, false));  // descending order
                    return my_seq;
                }

                std::vector<double> build_sum_processing_time_vector(const boost::numeric::ublas::matrix<double> &data,
                                                                     unsigned int nb_machines, unsigned int nb_jobs) {
                    std::vector<double> processing_time_vector(nb_jobs);
                    for (int i = 1; i <= nb_jobs; i++) {  // i = index_job
                        processing_time_vector[i - 1] = sum_processing_time_for_job(i, data, nb_machines);
                    }
                    return processing_time_vector;
                }

                std::vector<int>
                insertion(const std::vector<int> &sequence, unsigned int index_position, unsigned int value) {
                    // make a copy of the sequence vector --> new_seq = sequence[:]
                    std::vector<int> new_seq(sequence);
                    new_seq.insert(new_seq.begin() + index_position, value);
                    return new_seq;
                }

                PFSSolution solve(PFSInstance instance) {
                    int m = instance.m;
                    int n = instance.n;
                    ublas::matrix<double>& p_time = instance.p_time;
                    int pi[NJOBS+1];  // the permutation pi
                    PFSProblem problem;

                    NEH(n, m, p_time, pi);
                    double c_max_best = SchedulingUtil::Cmax(n, m, instance.p_time, pi);
                    std::vector<int> seq_current;
                    seq_current.assign(pi + 1, pi + n + 1);

                    // call the naive algorithm for result validation
                    /*
                    PFSSolution naive_sol = NEH_naive(instance);
                    if(naive_sol.value != c_max_best) {
                        std::cout << "Opt perm: " << seq_current << " vs Naive perm: " << naive_sol.permutation << "\n";
                        std::cout << "Opt value: " << c_max_best << " vs Naive value: " << naive_sol.value << "\n";
                        double c_max_best2 = problem.makespan(seq_current, instance)(m, n);
                        std::cout << "Opt value 2: " << c_max_best2 << "\n";
                    } */
                    double c_max_best2 = problem.makespan(seq_current, instance)(m, n);
                    assert(c_max_best2 == c_max_best);

                    if (seq_current.size() != n) {
                        LOG(ERROR) << "Invalid permutation size!";
                    }
                    return PFSSolution(seq_current, c_max_best, n);
                }

                void print_solution(std::vector<int> permutation, double cmax) {
                    // TODO implement printing in log file!
                }

                // Nawaz, Enscore, Ham's algorithm for F//Cmax problem,
                //        efficient implementation from Taillard's paper, O(n^2*m)
                // http://new.zsd.iiar.pwr.wroc.pl/files/sources/sax.cpp
                void NEH(int n, int m, ublas::matrix<double>& p, int pi[])
                { int i,j,k,l;
                    long c,cp,s,sum[NJOBS+1];
                    long r[MACHINES+1][NJOBS+1],q[MACHINES+2][NJOBS+2],d[MACHINES+1];

                    for (j=1;j<=n;j++)
                    { pi[j]=j; s=0;                                        // default pi[]
                        for (i=1;i<=m;i++) s+=p(i,j);
                        sum[j]=s;
                    }
                    SchedulingUtil::_sort(1,n,sum,pi);

                    for (i=0;i<=m;i++) r[i][0]=0;                         // r[][0] edge values
                    for (j=0;j<=n;j++) r[0][j]=0;                         // r[0][] edge values
                    d[0]=0;                                               // d[0] edge value

                    for (k=2;k<=n;k++)
                    {
                        for (i=1;i<=m;i++) for (j=1;j<=k;j++)
                                r[i][j]=max(r[i][j-1],r[i-1][j])+p(i,pi[j]);     // set new r[][]

                        for (i=0;i<=m;i++) q[i][k]=0;                       // q[][k] edge values
                        for (j=0;j<=k;j++) q[m+1][j]=0;                     // q[m+1][] edge values
                        for (i=m;i>=1;i--) for (j=k-1;j>=1;j--)
                                q[i][j]=max(q[i][j+1],q[i+1][j])+p(i,pi[j]);     // set new q[][]

                        cp=r[m][k]; i=k;                                    // insert on k
                        for (j=k-1;j>=1;j--)
                        { for (l=1;l<=m;l++) d[l]=max(d[l-1],r[l][j-1])+p(l,pi[k]); // set d[]
                            c=d[1]+q[1][j];
                            for (l=2;l<=m;l++) c=max(c,d[l]+q[l][j]);         // set new cmax
                            if (c<=cp) { cp=c; i=j; }                         // store best location
                        }
                        l=pi[k]; for (j=k;j>i;j--) pi[j]=pi[j-1]; pi[i]=l;  // adjust pi[]
                    }
                }

                problem::common::PFSSolution NEH_naive(PFSInstance instance) {
                    unsigned int nb_machines = instance.m;
                    unsigned int nb_jobs = instance.n;
                    std::vector<int> best_seq;
                    std::vector<int> order_seq = order_neh(instance.p_time, nb_machines, nb_jobs);
                    for(int i = 0; i < order_seq.size(); i++) {
                        order_seq[i]++;
                    }
                    std::vector<int> seq_current;
                    seq_current.push_back(order_seq[0]);  // seq_current = [order_seq[0]];
                    PFSProblem problem;

                    for (int i = 1; i < nb_jobs; i++) {
                        double min_cmax = std::numeric_limits<double>::max();
                        for (int j = 0; j < i + 1; j++) {
                            std::vector<int> tmp_seq = insertion(seq_current, j, order_seq[i]);
                            double cmax_tmp = problem.makespan(tmp_seq, instance)(nb_machines, tmp_seq.size());
                            //print_solution(tmp_seq, cmax_tmp);
                            if (min_cmax > cmax_tmp) {
                                best_seq = tmp_seq;
                                min_cmax = cmax_tmp;
                            }
                        }
                        seq_current = best_seq;
                    }
                    double c_max_best = problem.makespan(seq_current, instance)(nb_machines, nb_jobs);
                    if (seq_current.size() != nb_jobs) {
                        LOG(ERROR) << "Invalid permutation size!";
                    }
                    return PFSSolution(seq_current, c_max_best, nb_jobs);
                }

            };


/*
 * https://github.com/martinWANG2014/TP_Johnson_CDS_NEH/blob/master/neh.py
 import commonFunction
import interface

#####################################################
# algorithm for neh
#####################################################


def sum_processing_time(index_job, data, nb_machines):
    sum_p = 0
    for i in range(nb_machines):
        sum_p += data[i][index_job]
    return sum_p


def order_neh(data, nb_machines, nb_jobs):
    my_seq = []
    for j in range(nb_jobs):
        my_seq.append(j)
    return sorted(my_seq, key=lambda x: sum_processing_time(x, data, nb_machines), reverse=True)


def insertion(sequence, index_position, value):
    new_seq = sequence[:]
    new_seq.insert(index_position, value)
    return new_seq


def neh(data, nb_machines, nb_jobs):
    order_seq = order_neh(data, nb_machines, nb_jobs)
    seq_current = [order_seq[0]]
    for i in range(1, nb_jobs):
        min_cmax = float("inf")
        for j in range(0, i + 1):
            tmp_seq = insertion(seq_current, j, order_seq[i])
            cmax_tmp = commonFunction.makespan(tmp_seq, data, nb_machines)[nb_machines - 1][len(tmp_seq)]
            print(tmp_seq, cmax_tmp)
            if min_cmax > cmax_tmp:
                best_seq = tmp_seq
                min_cmax = cmax_tmp
        seq_current = best_seq
    return seq_current, commonFunction.makespan(seq_current, data, nb_machines)[nb_machines - 1][nb_jobs]


# run NEH
nbm, nbj, p_ij = commonFunction.read_from_file("example3.txt")
seq, cmax = neh(p_ij, nbm, nbj)
print("nbMachines:", nbm)
print("nbJobs:", nbj)
print("data: p_ij, the processing time of jth job on ith machine\n", p_ij)
print("neh:", seq, cmax)
interface.graphic("NEH", seq, nbj, nbm, commonFunction.makespan(seq, p_ij, nbm), p_ij)

 */

// https://stackoverflow.com/questions/29290778/slow-performance-using-stl-in-neh-algorithm
/*
#include <iostream>
#include <map>
#include <sstream>
#include <limits>
#include <vector>
#include <fstream>
#include <list>
#include <algorithm>
#include <numeric>
#include <queue>
#include <memory>


class Task
{
public:
    Task(unsigned int ordNum = 0) :m_ordNum(ordNum) { }
    unsigned int totalTasksTime() const;
    unsigned int ordNum() const { return m_ordNum; }
    unsigned int machinesNum() const { return m_params.size(); }
    unsigned int machineTaskTime(unsigned int machineNum) const { return m_params[machineNum - 1]; }
protected:
    std::vector<unsigned int> m_params;
    unsigned int m_ordNum;
    unsigned int m_totalTasksTime;
};


unsigned int Task::totalTasksTime() const
{
    return m_totalTasksTime;
}


class Instance
{
public:
    Instance() { }
    Instance(const std::string& name) :m_name(name) { }
    const std::string& name() const { return m_name; }
    void name(std::string& newName) { m_name = newName; }
    void neh(std::list<unsigned int>& permutation, unsigned int &totalTime) const;
    const Task* getTask(unsigned int taskNum) const { return &m_tasks[taskNum]; }
private:
    unsigned int calculateTotalTime(const std::list<unsigned int>& permutationList, unsigned int bestTimeFound) const;
    std::vector<Task> m_tasks;
    std::string m_name;
};


typedef std::map<unsigned int, unsigned int> MapIterator;
typedef std::vector<Task>::const_iterator TaskVecIterator;

bool compareTasksPtrBySumTime(const Task* t1, const Task* t2)
{
    unsigned int t1TotalTime = t1->totalTasksTime(), t2TotalTime = t2->totalTasksTime();

    bool w1 = t1TotalTime < t2TotalTime, w2 = t1TotalTime == t2TotalTime && t1->ordNum() > t2->ordNum();

    return w1 || w2;
}

void Instance::neh(std::list<unsigned int>& permutation, unsigned int &totalTime) const
{
    // Sorting tasks after total execution time
    std::list<const Task*> sortedTaskList;
    for (unsigned int i = 0; i < m_tasks.size(); i++)
        sortedTaskList.push_back(&m_tasks[i]);

    sortedTaskList.sort(compareTasksPtrBySumTime);

    while (!sortedTaskList.empty()) //
    {
        const Task* taskPtr = sortedTaskList.back(); sortedTaskList.pop_back();
        unsigned int taskNum = taskPtr->ordNum();

        std::list<unsigned int>::iterator bestPosition = permutation.begin();
        unsigned int bestTotalTime = std::numeric_limits<unsigned int>::max();
        permutation.push_front(taskNum);
        for (std::list<unsigned int>::iterator it = permutation.begin(); // procurando a melhor posicao
             it != permutation.end(); it++)
          {
            unsigned int currentTotalTime = calculateTotalTime(permutation, bestTotalTime);

            if (bestTotalTime > currentTotalTime)
              {
            bestTotalTime = currentTotalTime;
            bestPosition = it;
              }
            // TODO Resolver para o C++ 11
            auto nextIt = it; nextIt++;
            if (nextIt != permutation.end())
              std::swap(*it, *nextIt);

          }

        totalTime = bestTotalTime;
        permutation.insert(bestPosition, taskNum);
        permutation.pop_back();

    }
    std::cout << "Done:" << name() << std::endl;
}

unsigned int Instance::calculateTotalTime(const std::list<unsigned int>& permutationList, unsigned int bestTimeFound) const
{
    unsigned int rows = m_tasks[*permutationList.begin() - 1].machinesNum() + 1, columns = permutationList.size() + 1;
    ublas::matrix<int> matrix(rows, columns);
    unsigned int totalTime = 0;

    for (unsigned int c = 0; c < columns; c++)
        matrix(0, c) = 0;

    for (unsigned int r = 0; r < rows; r++)
        matrix(r, 0) = 0;

    std::list<unsigned int>::const_iterator it = permutationList.begin();
    for (unsigned int c = 1; c < columns; c++)
    {
        unsigned int taskNum = *it;
        for (unsigned int r = 1; r < rows; r++)
            (matrix(r, c) = std::max(matrix(r, c - 1), matrix(r - 1, c)) + m_tasks[taskNum - 1].machineTaskTime(r));// >bestTimeFound;
        //  return std::numeric_limits<unsigned int>::max();

        it++;
    }

    return matrix(rows - 1, columns - 1);
}

class InstanceVector
{
public:
    void neh(std::list< std::list<unsigned int> >& result) const;
    void neh(std::ostream& os) const;
private:
    std::vector<Instance> m_instances;
};



void InstanceVector::neh(std::list< std::list<unsigned int> >& results) const
{
    std::vector<Instance>::const_iterator it;
    for (it = m_instances.begin(); it != m_instances.end(); it++)
    {
        std::list<unsigned int> resultInstance;
        unsigned int totalTimeInstance;
        it->neh(resultInstance, totalTimeInstance);
        results.push_back(resultInstance);
    }
}

void InstanceVector::neh(std::ostream& os) const
{
    std::list< std::list<unsigned int> > results;
    for (std::vector<Instance>::const_iterator it = m_instances.begin();
        it != m_instances.end(); it++)
    {
        std::list<unsigned int> resultInstance;
        unsigned int totalTimeInstance;
        it->neh(resultInstance, totalTimeInstance);
        // TODO Resolver para o C++ 11
        results.push_back(std::move(resultInstance));
    }

    for (std::list< std::list<unsigned int> >::const_iterator it = results.begin();
        it != results.end(); it++)
    {
        for (std::list<unsigned int>::const_iterator itPermutation = it->begin(); itPermutation != it->end(); itPermutation++)
            os << *itPermutation << " ";
        os << std::endl;
    }
}
*/

        } /* namespace heuristic */
    } //* namespace pfsp */
} /* namespace problem */


#endif /* SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_NEHSOLVER_H_ */
