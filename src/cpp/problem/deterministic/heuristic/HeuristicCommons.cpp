/*
 * HeuristicCommons.cpp
 *
 *  Created on: 22 de mar de 2018
 *      Author: mlevorato
 */

#include "include/HeuristicCommons.h"
#include <glog/logging.h>
#include <vector>
#include <boost/numeric/ublas/storage.hpp>

using std::vector;

#include <iostream>
#include "../../GlobalTypes.h"

using std::ostream;
using namespace boost;
using namespace std;

namespace problem {
    namespace pfsp {
        namespace heuristic {


    HeuristicCommons::HeuristicCommons() {
        // TODO Auto-generated constructor stub

    }

    HeuristicCommons::~HeuristicCommons() {
        // TODO Auto-generated destructor stub
    }


/**
 * Jobs are sequenced in nondecreasing order (ascending) of processing times.
 * @param data
 * @param nbj
 * @return
 */
    std::vector<int> HeuristicCommons::u(ublas::matrix<double> &data, unsigned int nbj) {
        std::vector<int> my_u;
        for (int j = 1; j <= nbj; j++) {  // Get all jobs j with p_1j <= p_2j
            if (data(1, j) <= data(2, j)) {
                my_u.push_back(j);
            }
        }
        std::vector<double> row_vector = util::CollectionUtil::matrix_row_to_vector(data, 1);
        //std::cout << "Naive vetor u: " << my_u[my_u.size() - 1] << "\n";

        // sequence the jobs in nondecreasing order of processing time
        // LOG(INFO) << "p_time: " << print(row_vector);
        std::sort(my_u.begin(), my_u.end(), util::vector_comparator(&row_vector, true));  // ascending order
        // LOG(INFO) << "Vector after sorting by p_time: " << print(my_u);
        return my_u;
    }

/**
 * Jobs are sequenced in nonincreasing order (descending) of processing times.
 * @param data
 * @param nbj
 * @return
 */
    std::vector<int> HeuristicCommons::v(ublas::matrix<double> &data, unsigned int nbj) {
        std::vector<int> my_v;
        for (int j = 1; j <= nbj; j++) {  // Get all jobs j with p_1j > p_2j
            if (data(1, j) > data(2, j)) {
                my_v.push_back(j);
            }
        }
        std::vector<double> row_vector = util::CollectionUtil::matrix_row_to_vector(data, 2);
        //std::cout << "Naive vetor v: " << my_v[0] << "\n";

        // sequence the jobs in nonincreasing order of processing time
        // LOG(INFO) << "p_time: " << print(row_vector);
        std::sort(my_v.begin(), my_v.end(), util::vector_comparator(&row_vector, false));  // descending order
        //std::cout << "Vector after sorting by p_time: " << print(my_v);
        std::vector<int> ptime_bb;
        //std::cout << "Ptime naive bb: ";
        //for(int x = 0; x < my_v.size(); x++)  std::cout << (row_vector[my_v[x]]) << " ";
        //std::cout << "\n";
        return my_v;
    }

    std::string HeuristicCommons::print(const std::vector<double> &v) {
        std::stringstream ss;
        if (!v.empty()) {
            ss << "[ ";
            std::copy(v.begin(), v.end(), std::ostream_iterator<double>(ss, ", "));
            ss << " ]";
        }
        return ss.str();
    }

    std::string HeuristicCommons::print(const std::vector<int> &v) {
        std::stringstream ss;
        if (!v.empty()) {
            ss << "[ ";
            std::copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, ", "));
            ss << " ]";
        }
        return ss.str();
    }

        } /* namespace heuristic */
    } //* namespace pfsp */
} /* namespace problem */

