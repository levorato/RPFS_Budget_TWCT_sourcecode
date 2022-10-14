/*
 * RandomUtil.h
 *
 *  Created on: 06/05/2014
 *      Author: czt0
 */

#ifndef RANDOMUTIL_H_
#define RANDOMUTIL_H_

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <random>
#include <numeric>
#include "../random.h"

namespace util {

// other options of generator: minstd_rand, mt19937
typedef boost::minstd_rand RandomGeneratorType;

class RandomUtil {
public:
	virtual ~RandomUtil();
	RandomUtil(std::uint_least32_t seed);

	static void setSeed(int seed) {
		// other option for seed: generator.seed(boost::random::random_device()());
        rg.seed(seed);
		isSeeded = true;
	}

	// TODO Change to std rng
	/**
	 * Returns a random integer value uniformly distributed in the set of integer numbers {lowerLimit, ..., upperLimit}.
	 * @param lowerLimit
	 * @param upperLimit
	 * @return
	 * DEPRECATED
	 */
	static int next(int lowerLimit, int upperLimit) {
		// we are supposing the random generator is already seeded here
		rg.seed(boost::random::random_device()());
		boost::uniform_int<> distribution(lowerLimit, upperLimit);
		boost::variate_generator<RandomGeneratorType&, boost::uniform_int<> > LimitedInt(
				rg, distribution);
		return LimitedInt();
	}

	// DEPRECATED
	static double nextDouble_0_1() {
        // we are supposing the random generator is already seeded here
        rg.seed(boost::random::random_device()());
        boost::uniform_real<> distribution(0.0, 1.0);
        boost::variate_generator<RandomGeneratorType&, boost::uniform_real<> > LimitedReal(
                rg, distribution);
        return LimitedReal();
	}

    static std::double_t nextDouble()
	{
		return nextDouble(0.0, 1.0);
	}

    static std::double_t nextDouble(std::double_t minValue, std::double_t maxValue)
	{
		if (minValue < 0.0 || maxValue < minValue)
		{
			throw std::invalid_argument("minValue and maxValue must be non-negative. maxValue must be > minvalue");
		}
		std::uniform_real_distribution<std::double_t> distribution(minValue, maxValue);
		return distribution(rng);
	}

    static std::int32_t nextInt()
	{
		return next(0, std::numeric_limits<std::int32_t>::max());
	}

	/**
	 * Returns a random integer value uniformly distributed in the set of integer numbers {0, ..., maxValue}.
	 * @param maxValue
	 * @return
	 */
    static std::int32_t nextInt(int32_t maxValue)
	{
        std::uniform_int_distribution<> dis(0);
        return (dis(rng) % (maxValue + 1));
		//return next(0, maxValue);
	}

	/**
	 * Returns a random integer value uniformly distributed in the set of integer numbers {minValue, ..., maxValue}.
	 * @param minValue
	 * @param maxValue
	 * @return
	 */
    static std::int32_t nextInt(int32_t minValue, int32_t maxValue)
	{
		if (minValue < 0 || maxValue < minValue)
		{
			throw std::invalid_argument("minValue and maxValue must be non-negative. maxValue must be > minvalue");
		}
        unsigned left = maxValue - minValue + 1;
        std::uniform_int_distribution<> dis(0);
        return minValue + (dis(rng) % left);
        // std::uniform_int_distribution<std::int32_t> distribution(minValue, maxValue);
		// return distribution(rng);
	}

    void nextBytes(std::vector<uint8_t>& buffer)
	{
		for (auto &i : buffer)
		{
			i = static_cast<std::uint8_t>(byteDistribution(rng));
		}
	}

	/**
	 * Generate a shuffled vector of ints with integers in the range [initial_value, initial_value + n).
	 * @param n number of elements in the vector
	 * @return shuffled vector of ints
	 */
	static std::vector<int> generate_random_vector(const int &initial_value, const unsigned &n) {
        std::vector<int> v(n);
        std::iota (std::begin(v), std::end(v), initial_value); // Fill with 1, 2, ..., n
        std::shuffle(v.begin(), v.end(), rng);
        return v;
    }

    /**
     * Shuffle a vector v
     * @param v shuffled vector
     */
    static inline void shuffle_vector(std::vector<int> &v) {
        std::shuffle(v.begin(), v.end(), rng);
	}

	static RandomGeneratorType rg;
	static bool isSeeded;
	std::uniform_int_distribution<std::int32_t> byteDistribution;

};

} /* namespace util */
#endif /* RANDOMUTIL_H_ */
