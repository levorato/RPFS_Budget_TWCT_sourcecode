//
// Created by mlevorato on 11/27/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_RANDOMNESS_H
#define FLOWSHOP_SOLVER_GRASP_RANDOMNESS_H

#include "Inputs.h"
#include "Test.h"
#include "RandomUtil.h"

namespace problem {
namespace pfsp {

    using namespace util;

    template<class ProblemInstance>
    class Randomness {
    private:
        Test aTest;
        Inputs<ProblemInstance> inputs;

    public:
        Randomness(Test test, Inputs<ProblemInstance> inputData) : aTest(test), inputs(inputData)
        {
        }

        std::vector<int> calcPositionsArray(const std::string &distribution)
        {
            int nJobs = inputs.getNumberOfJobs();
            if(distribution.at(0) == 'u') {  // uniform
                return RandomUtil::generate_random_vector(0, nJobs);
            }
            std::vector<int> posArray(nJobs, 0);
            std::vector<int> auxArray(nJobs, 0); // array of "pointers" to jobs in effList

            // Reset auxArray
            for(int i = 0; i < nJobs; i++ )
                auxArray[i] = i;

            // Assign new random positions
            for(int i = 0; i < nJobs; i++ )
            {
                int pos = getRandomPosition(nJobs - i, distribution, aTest.getBeta1());
                posArray[i] = auxArray[pos];
                // System.arraycopy(auxArray, pos+1, auxArray, pos, nJobs-i-1-pos);
        		for(int j = pos; j < nJobs - i - 1; j++ ) {
                    auxArray[j] = auxArray[j + 1];
                }
            }
            return posArray;
        }

        int getRandomPosition(int n, const std::string &dist) {
            return getRandomPosition(n, dist, 0.0);
        }

        int getRandomPosition(int n, std::string dist, double beta) // random between 0 and n-1
        {
            int pos = 0;
            // Assuming Uniform Distribution
            pos = RandomUtil::nextInt(n - 1);  // replaced nextInt ?
            return pos;
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_RANDOMNESS_H
