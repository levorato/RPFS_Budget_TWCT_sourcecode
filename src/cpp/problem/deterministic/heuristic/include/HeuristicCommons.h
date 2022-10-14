/*
 * HeuristicCommons.h
 *
 *  Created on: 22 de mar de 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_HEURISTICCOMMONS_H_
#define SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_HEURISTICCOMMONS_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include "../../../../util/include/CollectionUtil.h"

using namespace boost;

namespace problem {
	namespace pfsp {
		namespace heuristic {

class HeuristicCommons {
public:
	HeuristicCommons();
	virtual ~HeuristicCommons();

	static std::vector<int> u(boost::numeric::ublas::matrix<double>& data, unsigned int nbj);
	static std::vector<int> v(boost::numeric::ublas::matrix<double>& data, unsigned int nbj);

	static std::string print(const std::vector<double>& v);
    static std::string print(const std::vector<int>& v);
};

		} /* namespace heuristic */
	} //* namespace pfsp */
} /* namespace problem */

#endif /* SRC_CPP_PROBLEM_HEURISTIC_INCLUDE_HEURISTICCOMMONS_H_ */
