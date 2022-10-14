/*
 * JohnsonSolver.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_JOHNSONSOLVER_H_
#define SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_JOHNSONSOLVER_H_
#pragma once
#include <glog/logging.h>

#include "../../PFSInstance.h"
#include "../../../PFSSolution.h"
#include "../../PFSProblem.h"
#include "HeuristicCommons.h"
#include "../../../../util/include/SchedulingUtil.h"

namespace problem {
    namespace pfsp {
        namespace heuristic {

            using namespace util;
            using namespace problem::common;

            class JohnsonSolver {
            public:
                JohnsonSolver();

                virtual ~JohnsonSolver();

                static PFSSolution solve(PFSInstance instance) {
                    int m = instance.m;
                    int n = instance.n;
                    ublas::matrix<double>& p_time = instance.p_time;
                    return solve(m, n, p_time);
                }

                // johnson's algorithm for F2//Cmax problem
                // http://new.zsd.iiar.pwr.wroc.pl/files/sources/sax.cpp
                static void johnson(int n, long aa[], long bb[], int pi[]) {
                    int a=0,b=n+1,j;

                    for (j=1;j<=n;j++) { if (aa[j]<=bb[j]) pi[++a]=j; else pi[--b]=j; }
                    std::vector<int> test(n + 1), pi1, pi2;
                    // test.assign(pi + 1, pi + n + 1);
                    // std::cout << "Opt vetor pi: " << test << "\n";
                    // std::cout << "a: " << a << "\n";
                    // std::cout << "b: " << b << "\n";
                    if (a>0){   SchedulingUtil::sort(1,a,aa,pi); pi1.assign(pi + 1, pi + a + 1);  }
                    if (b<=n){  SchedulingUtil::_sort(b,n,bb,pi); pi2.assign(pi + b, pi + n + 1); }
                    // std::cout << "Opt vetor pi final: " << pi1 << " " << pi2 << "\n";
                }

                static PFSSolution johnson_naive(PFSInstance instance) {
                    unsigned int nb_machines = instance.m;
                    unsigned int nb_jobs = instance.n;
                    PFSProblem problem;
                    // seq_current = u(data, nb_jobs) + v(data, nb_jobs)
                    std::vector<int> seq_current = HeuristicCommons::u(instance.p_time, nb_jobs);
                    std::vector<int> vet_v = HeuristicCommons::v(instance.p_time, nb_jobs);
                    seq_current.insert(seq_current.end(), vet_v.begin(), vet_v.end());
                    // LOG(INFO) << "Final vector of jobs: " << HeuristicCommons::print(seq_current);

                    double c_max_best = problem.makespan(seq_current, instance)(nb_machines, nb_jobs);
                    if (seq_current.size() != nb_jobs) {
                        LOG(ERROR) << "Invalid permutation size!";
                    }
                    return PFSSolution(seq_current, c_max_best, nb_jobs);
                }

                static PFSSolution solve(int m, int n, const ublas::matrix<double>& p_time,
                        const int& algorithm = 1) {
                    if(algorithm <= 1) {
                        Permutation s(n, 0);
                        std::iota(std::begin(s), std::end(s), 1);
                        auto mid = partition(s.begin(),s.end(),[&p_time](const int& i) { return p_time(1,i)<p_time(2,i);});
                        std::sort(s.begin(),mid,[&p_time](const int& i, const int& j) { return p_time(1,i)<p_time(1,j);});
                        std::sort(mid,s.end(),[&p_time](const int& i, const int& j) { return p_time(2,i)>p_time(2,j);});
                        Time cmax_s = PFSProblem::calculateCmax(s, 2, n, p_time);
                        return PFSSolution(s, cmax_s, n);
                    } else {  // if(algorithm >= 2) {
                        int pi[NJOBS + 1];  // the permutation pi
                        long aa[NJOBS + 1];  // aa[] is the processing time vector for machine 1
                        long bb[NJOBS + 1];  // bb[] is the processing time vector for machine 2
                        for (int job_count = 1; job_count <= n; job_count++) {
                            aa[job_count] = p_time(1, job_count);
                            bb[job_count] = p_time(2, job_count);
                        }
                        johnson(n, aa, bb, pi);
                        double c_max_best = SchedulingUtil::Cmax(n, m, p_time, pi);
                        std::vector<int> seq_current(n);
                        seq_current.assign(pi + 1, pi + n + 1);
                        if (seq_current.size() != n) {
                            LOG(ERROR) << "Invalid permutation size!";
                        }
                        return PFSSolution(seq_current, c_max_best, n);
                    }
                }

                /**
                 * Arrange the remaining (n - l) unassigned jobs in Johnson order.
                 * @param n  number of jobs
                 * @param partial_seq  partial sequence with currently assigned jobs
                 * @param p  processing time matrix
                 * @return complete permutation with all jobs positioned.
                 */
                static Permutation remaining_johnson_order(const int &n, const Permutation &partial_seq,
                                                    const boost::numeric::ublas::matrix<Time>& p) {
                    // A) Find the jobs that are not present in partial_seq and store them in seq_1
                    std::vector<int> S;
                    std::vector<bool> present(n + 1, false);
                    for (int x = 0; x < partial_seq.size(); x++) {
                        int j = partial_seq[x];
                        present[j] = true;
                    }
                    for (int j = 1; j <= n; j++) {
                        if(! present[j]) {
                            S.push_back(j);  // store the unassigned job j
                        }
                    }
                    // B) sort the jobs in S in Johnson order, according to the matrix {p}_{ij}
                    std::vector<int> seq_2, vet_v2;
                    for (int x = 0; x < S.size(); x++) {
                        int j = S[x];
                        if (p(1, j) < p(2, j)) {  seq_2.push_back(j);
                        } else {  vet_v2.push_back(j);  }
                    }
                    std::vector<double> row_vector_0 = util::CollectionUtil::matrix_row_to_vector(p, 1);
                    std::sort(seq_2.begin(), seq_2.end(), util::vector_comparator(&row_vector_0, true));  // ascending order
                    std::vector<double> row_vector_1 = util::CollectionUtil::matrix_row_to_vector(p, 2);
                    std::sort(vet_v2.begin(), vet_v2.end(), util::vector_comparator(&row_vector_1, false));  // descending order
                    seq_2.insert(seq_2.end(), vet_v2.begin(), vet_v2.end());
                    // merge the partial sequence with the remaining jobs (sorted in Johnson order)
                    std::vector<int> seq_1(partial_seq);
                    seq_1.insert(seq_1.end(), seq_2.begin(), seq_2.end());
                    return seq_1;
                }
            };

        } /* namespace heuristic */
    } //* namespace pfsp */
} /* namespace problem */

#endif /* SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_JOHNSONSOLVER_H_ */
