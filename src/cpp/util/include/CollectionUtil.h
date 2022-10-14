/*
 * CollectionUtil.h
 *
 *  Created on: 22 de mar de 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_UTIL_INCLUDE_COLLECTIONUTIL_H_
#define SRC_CPP_UTIL_INCLUDE_COLLECTIONUTIL_H_

#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/irange.hpp>
#include "../../problem/GlobalTypes.h"

namespace util {

    using namespace std;
    using namespace boost;


    class CollectionUtil {
    public:
        CollectionUtil();

        virtual ~CollectionUtil();

        static std::vector<double>
        matrix_row_to_vector(const ublas::matrix<double> &matrix, const unsigned int &row_number) {
            std::vector<double> row_as_vector(matrix.size2());
            const ublas::matrix_row<const ublas::matrix<double> > mr(matrix, row_number);
            std::copy(mr.begin(), mr.end(), row_as_vector.begin());
            return row_as_vector;
        }

        static std::vector<double>
        ublas_vector_to_std_vector(const ublas::vector<double> &v) {
            std::vector<double> row_as_vector(v.size());
            std::copy(v.begin(), v.end(), row_as_vector.begin());
            return row_as_vector;
        }

        /**
         * Returns a permutation of jobs in the format {1, 2, ... , n}
         * @return
         */
        static std::vector<int> create_initial_job_list(int n) {
            std::vector<int> v;
            boost::range::push_back(v, boost::irange(1, n + 1));
            return v;
        }
    };


/**
 * Comparison class for a given ublas::vector<double> object
 */
    class vector_comparator {
    public:
        vector_comparator(std::vector<double> *v, bool asc) : vet(v), ascending(asc) {}

        ~vector_comparator() {}

        // TODO implement proper double comparison (precision)!
        bool operator()(unsigned int i, unsigned int j) {
            if (ascending) {
                return ((*vet)[i] < (*vet)[j]);
            } else {
                return ((*vet)[i] > (*vet)[j]);
            }
        }

        std::vector<double> *vet;
        bool ascending;
    };

} /* namespace util */

#endif /* SRC_CPP_UTIL_INCLUDE_COLLECTIONUTIL_H_ */
