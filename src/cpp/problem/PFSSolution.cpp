/*
 * PFSSolution.cpp
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#include "PFSSolution.h"
#include "../util/include/FileUtil.h"

using namespace problem::common;

PFSSolution::~PFSSolution() {
	// TODO Auto-generated destructor stub
}

namespace problem {
    namespace common {
        std::ostream &operator<<(std::ostream &os, const PFSSolution &sol) {
            os << "Objective value: " << sol.value << "\nPermutation: " << sol.permutation;
            return os;
        }

        void PFSSolution::sendToFile(const string &filePath)
        {
            stringstream ss;
            ss << "Objective value: " << value << "\nPermutation: " << permutation << "\n";
            // LOG(INFO) << "[PFSSolution] Saving results file to " << filePath;
            util::FileUtil::saveTextContentsToFile(filePath, ss.str(), false);
        }
    }
}
