//
// Created by Mario Costa Levorato Junior on 2019-06-18.
//

#ifndef FLOWSHOP_SOLVER_PARTIALPERMUTATION_H
#define FLOWSHOP_SOLVER_PARTIALPERMUTATION_H

#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "../../GlobalTypes.h"

#define MAX_TIME std::numeric_limits<Time>::max()

class PartialPermutation {
public:
    std::vector<t_job> pi;
    unsigned n1,n2;

    PartialPermutation(unsigned n, unsigned k) : pi(k+1), n1(0), n2(k+1) {
        iota(pi.begin(),pi.end(),0);
    }

    PartialPermutation(unsigned n) : PartialPermutation(n,n) {}

    PartialPermutation(PartialPermutation&& other) {
        this->swap(other);
    }

    PartialPermutation(const PartialPermutation& other) {
        pi = other.pi;
        n1 = other.n1;
        n2 = other.n2;
    }

    PartialPermutation& operator=(PartialPermutation other) {
        this->swap(other);
        return *this;
    }

    void swap(PartialPermutation& other) {
        using std::swap;
        pi.swap(other.pi);
        swap(n1,other.n1);
        swap(n2,other.n2);
    }

    void addi(unsigned ji) {
        n1++;
        std::swap(pi[n1],pi[ji]);
    }

    void ddai(unsigned ji) {
        n2--;
        std::swap(pi[n2],pi[ji]);
    }

    bool empty() const {
        if (n1>0)
            return false;
        if (n2<pi.size())
            return false;
        return true;
    }

    unsigned totalJobs() const {
        return pi.size()-1;
    }

    bool full() const {
        return n1+1 >= n2;
    }

    bool determined() const {
        return n1+1 >= n2;
    }

    unsigned numFixed() const {
        return n1+pi.size()-n2;
    }

    unsigned numFree() const {
        return n2-n1-1;
    }

    bool free(t_job j) {
        return std::find(pi.begin()+n1+1,pi.begin()+n2,j) != pi.begin()+n2;
    }
};


#endif //FLOWSHOP_SOLVER_PARTIALPERMUTATION_H
