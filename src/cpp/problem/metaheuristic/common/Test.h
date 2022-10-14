//
// Created by mlevorato on 11/28/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_TEST_H
#define FLOWSHOP_SOLVER_GRASP_TEST_H

#include <string>
#include "../../../util/RandomStream.h"

namespace problem {
namespace pfsp {

    using namespace std;
    using namespace util;

    class Test {
    private:
        string instanceName;
        int maxTime; // maximum computational time (in seconds)
        int nIter; // number of first-level iterations
        string distribution;
        double beta1;
        double beta2;
        int seed;
        RandomStream random_ng; // Random Number Generator

        // Parameters used by VND local search
        std::vector<unsigned> vnd_permutation;
        unsigned vnd_size;
        bool first_improvement;
        bool random_vnd;
        bool adaptive_construction;

        // Parameters used internally by the metaheuristics ***************************************
        // IMPORTANT: boolean that controls if the objective function value should be recalculated
        // (using a specific function) in a specific point of the metaheuristic:
        // 0 == no recalculation;
        // 1 == recalculated at the end of each MH iteration;
        // 2 == recalculated at the end of each MH iteration and also at each VND size change.
        // 3 and 4 == <reserved for future use>
        unsigned force_recalculate_obj;

    public:
        Test(string name, int time, int n, string d, float b1, float b2, int s,
                const std::vector<unsigned> &p_vnd_permutation, const unsigned &p_vnd_size,
                const bool &p_first_improvement, const bool &p_random_vnd, const bool &p_adaptive_construction,
                const unsigned &p_force_recalculate_obj) :
            instanceName(name), maxTime(time), nIter(n), distribution(d), beta1(b1), beta2(b2), seed(s), random_ng(s),
            vnd_permutation(p_vnd_permutation), vnd_size(p_vnd_size), first_improvement(p_first_improvement),
            random_vnd(p_random_vnd), adaptive_construction(p_adaptive_construction),
            force_recalculate_obj(p_force_recalculate_obj) {
        }

        string getInstanceName() const
        {
            return instanceName;
        }

        int getMaxTime() const
        {
            return maxTime;
        }

        int getNIter() const
        {
            return nIter;
        }

        string getDistribution() const
        {
            return distribution;
        }

        double getBeta1() const
        {
            return beta1;
        }

        double getBeta2() const
        {
            return beta2;
        }

        int getSeed() const
        {
            return seed;
        }

        RandomStream& getRandomStream()
        {
            return random_ng;
        }

        std::vector<unsigned> get_vnd_permutation() {
            return vnd_permutation;
        }

        unsigned get_vnd_size() {
            return vnd_size;
        }

        bool is_first_improvement() {
            return first_improvement;
        }

        bool is_RVND() {
            return random_vnd;
        }

        bool is_adaptive_construction() {
            return adaptive_construction;
        }

        unsigned get_force_recalculate_obj() {
            return force_recalculate_obj;
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_TEST_H
