/*
 * PFSSolution.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_INCLUDE_PFSSOLUTION_H_
#define SRC_CPP_PROBLEM_INCLUDE_PFSSOLUTION_H_

#include <iostream>
#include <vector>
#include <algorithm>

#include "GlobalTypes.h"
#include "./metaheuristic/common/Job.h"

using namespace std;

namespace problem {
    namespace common {

        class PFSSolution {
        public:
            PFSSolution() : permutation(), value(0.0), permutation_size(0), lower_bound(0.0), time_spent(0.0), 
                is_optimal(false), job_list(), num_iterations(0), validated(false), gap(100.0)  {

            }

            PFSSolution(Permutation perm, double sol_val, unsigned long perm_size) : permutation(perm), value(sol_val),
                                                                                     permutation_size(perm_size),
                                                                                     lower_bound(sol_val), time_spent(0.0),
                                                                                     is_optimal(false), job_list(),
                                                                                     num_iterations(0), num_improvements(0),
                                                                                     validated(false), gap(100.0) 
            {

            }

            PFSSolution(Permutation perm, double sol_val, unsigned long perm_size, double LB) : permutation(perm), value(sol_val),
                                                                                                permutation_size(perm_size), lower_bound(LB),
                                                                                                time_spent(0.0), is_optimal(false),
                                                                                                job_list(),
                                                                                                num_iterations(0), num_improvements(0),
                                                                                                validated(false), gap(100.0) 
            {

            }

            /**
             * Build a PFSP solution using 2 stacks, one with jobs of the left (head)
             * and the other one with jobs of the right (tail).
             */
            PFSSolution(std::vector<int> head, std::vector<int> tail, unsigned int perm_size) :
                    permutation(perm_size, 0), value(0), permutation_size(perm_size), time_spent(0.0), is_optimal(false) {

                permutation = generate_permutation_from_head_tail(head, tail, perm_size);
            }

            /**
             * Build a PFSP solution using only one stack, with jobs of the left (head).
             */
            PFSSolution(std::vector<int> head, unsigned int perm_size) :
                    permutation(perm_size, 0), value(0), permutation_size(perm_size), time_spent(0.0), is_optimal(false) {
                permutation = head;
            }

            virtual ~PFSSolution();

            const Permutation& getPermutation() const {
                return permutation;
            }

            std::string getPermutationAsString() const {
                stringstream ss;
                ss << "[";
                for(int i : permutation) {
                    ss << i << " ";
                }
                ss << "]";
                return ss.str();
            }

            friend std::ostream& operator<<(std::ostream& os, const PFSSolution& sol);

            static Permutation generate_permutation_from_head_tail(std::vector<int> head,
                                                                   std::vector<int> tail, int permutation_size) {
                Permutation perm(permutation_size, 0);
                int i = 0;
                int size_difference = permutation_size - (head.size() + tail.size());
                for(i = 0; i < head.size(); i++) {
                    perm[i] = head[i];
                }
                for(int j = 0; j < size_difference; j++, i++) {
                    perm[i] = -1;
                }
                for(int j = tail.size() - 1; j >= 0; j--, i++) {
                    perm[i] = tail[j];
                }
                return perm;
            }

            std::vector<pfsp::Job> getJobs() {
                return job_list;
            }

            Permutation permutation;
            double value;  // solution value
            double lower_bound;  // if the solution is obtained from a branch and bound algorithm which has not finished
            long permutation_size;
            double time_spent;
            bool is_optimal;
            std::vector<pfsp::Job> job_list;
            long num_iterations;
            long num_improvements;
            bool validated;
            double gap;

            /**
             * Indicates that the current solution is not complete (i.e. does not have a full permutation).
             */
            bool is_partial_solution() {
                return permutation.size() != permutation_size;
            }

            bool is_approximate() {
                return value != lower_bound;
            }

            void sendToFile(const string &filePath);
        };

        inline std::ostream& operator << (std::ostream& os, const std::vector<int>& v)
        {
            os << "[";
            int neg_count = 0;
            for (typename std::vector<int>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
            {
                if((*ii) < 0) {
                    neg_count++;
                } else {
                    if(neg_count > 0) {
                        os << " [" << neg_count << " blank(s) ]";
                        neg_count = 0;
                    }
                    os << " " << *ii;
                }
            }
            os << " ]";
            return os;
        }
    }
}


#endif /* SRC_CPP_PROBLEM_INCLUDE_PFSSOLUTION_H_ */
