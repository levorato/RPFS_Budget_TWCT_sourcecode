//
// Created by Mario Costa Levorato Junior on 2019-01-15.
//

#ifndef FLOWSHOP_SOLVER_GLOBALTYPES_H
#define FLOWSHOP_SOLVER_GLOBALTYPES_H

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <cstdint>

// IMPORTANT: BRANCH AND BOUND CONVERGENCE AND PRECISION EPS
#include "../util/precision.h"

// -----------------------------------------------------------
typedef std::vector<int> Permutation;
namespace ublas = boost::numeric::ublas;
typedef double Time;
typedef unsigned t_job;
// the following types are used by the ILS procedure ---------
using t_time = double; // uint_fast64_t;
//using t_job = int_fast16_t;
using t_machine = int_fast16_t;
using t_flow_time = int_fast64_t;
using t_size_type = uint_fast32_t;
// ------------------------------------------------------------

enum BranchingType { alternating, smallest_dynamic };

enum NodeSelectionStrategy { dfs, cplex };

enum BranchingDirection { no_dir, prefix, suffix };

enum MIPUsage { no_mip, linear_relax, full_mip };

enum Robust_UB_Type { grasp };

// custom specialization of std::hash can be injected in namespace std
namespace std
{
    template<class T>
    struct hash<std::vector<T>> {
        typedef std::vector<T> argument_type;
        typedef std::size_t result_type;
        result_type operator() (argument_type const& key) const noexcept {
            std::hash<T> hasher;
            size_t result = 0;   // I would still seed this.
            for(size_t i = 0; i < key.size(); ++i) {
                // h ^= std::hash<int>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
                result = (result << 1) ^ hasher(key[i]); // ??
            }
            return result;
        }

        // TODO colocar hash como membro da classe
        // std::hash<T> *hasher;
    };
}

#endif //FLOWSHOP_SOLVER_GLOBALTYPES_H
