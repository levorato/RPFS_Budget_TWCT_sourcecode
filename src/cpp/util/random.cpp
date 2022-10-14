#include "random.h"

#include <time.h>

#include <fstream>
using namespace std;

mt19937 rng;

#ifdef _MSC_VER 
#define srand48(x) srand((int)(x))
#define drand48() ((double)rand()/RAND_MAX)
#endif

unsigned long setupRandom(unsigned long seed) {
    if (seed==0) {
        seed = time(0);
        ifstream f("/dev/urandom");
        if (f.good()) {
            f.read((char *)(&seed), sizeof(unsigned int));
        }
    }
    rng.seed(seed);
    srand48(seed);
    return seed;
}
