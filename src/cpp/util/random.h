#ifndef FLOWSHOP_SOLVER_RANDOM_H
#define FLOWSHOP_SOLVER_RANDOM_H

#pragma once

#include <random>

extern std::mt19937 rng;

unsigned long setupRandom(unsigned long seed=0);


#endif //FLOWSHOP_SOLVER_RANDOM_H
