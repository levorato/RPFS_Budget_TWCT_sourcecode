/*
 * RandomUtil.cpp
 *
 *  Created on: 06/05/2014
 *      Author: czt0
 */

#include "include/RandomUtil.h"

namespace util {

bool RandomUtil::isSeeded = false;
RandomGeneratorType RandomUtil::rg;

// https://codereview.stackexchange.com/questions/83603/random-class-in-c
RandomUtil::RandomUtil(std::uint_least32_t seed)
        : byteDistribution(0, 256)
    {
        setupRandom(seed);
    }

RandomUtil::~RandomUtil() {
	// TODO Auto-generated destructor stub
}

} /* namespace util */
