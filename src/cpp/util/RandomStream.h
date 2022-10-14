//
// Created by mlevorato on 11/28/19.
//

#ifndef FLOWSHOP_SOLVER_RANDOMSTREAM_H
#define FLOWSHOP_SOLVER_RANDOMSTREAM_H

#include "random.h"  // includes std::mt19937 rng (randomNumberGenerator) as global variable (extern)
#include <stdexcept>

namespace util {

    /**
     *  Random number generator class.
     */
    class RandomStream {
    public:
        RandomStream(const long &p_seed) : seed(p_seed) {

        }

        std::double_t nextDouble()
        {
            return this->nextDouble(0.0, 1.0);
        }

        std::double_t nextDouble(std::double_t minValue, std::double_t maxValue)
        {
            if (minValue < 0.0 || maxValue < minValue)
            {
                throw std::invalid_argument("minValue and maxValue must be non-negative. maxValue must be > minvalue");
            }
            std::uniform_real_distribution<std::double_t> distribution(minValue, maxValue);
            return distribution(rng);
        }

        int nextInt(int minValue, int maxValue) {
            if (minValue < 0.0 || maxValue < minValue)
            {
                throw std::invalid_argument("minValue and maxValue must be non-negative. maxValue must be > minvalue");
            }
            std::uniform_int_distribution<> distribution(minValue, maxValue);
            return distribution(rng);
        }

    private:
        long seed;
    };
}

#endif //FLOWSHOP_SOLVER_RANDOMSTREAM_H
