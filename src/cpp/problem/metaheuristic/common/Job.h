//
// Created by Mario Costa Levorato Junior on 2019-02-01.
//

#ifndef FLOWSHOP_SOLVER_GRASP_JOB_H
#define FLOWSHOP_SOLVER_GRASP_JOB_H

#include <string>
#include <sstream>
#include <vector>
#include "../../GlobalTypes.h"
#include "RandomUtil.h"

namespace problem {
namespace pfsp {

    using namespace util;

    /**
     * This class represents a job in the flowshop problem.
     */
    class Job {

    public:
        int id; // job ID
        double cost;

        Job(int order, int nMachines) : id(order + 1), cost(0.0)
        {
        }

        Job(int _id) : id(_id), cost(0.0) {

        }

        void setId(int jobId) {
            id = jobId;
        }

        inline int getId() const {
            return id;
        }

        bool operator==(const Job& rhs) const {
            return id == rhs.id;
        }

        bool operator<(const Job& rhs) const {
            return id < rhs.id;
        }

        bool operator<(const int& rhs) const {
            return id < rhs;
        }

        friend std::ostream& operator<<( std::ostream& o, const pfsp::Job& j) {
            return o << j.id;
        }

        std::string toString() {
            std::stringstream ss;
            ss << "\n Job Id: " << this->getId() << " ";
            return ss.str();
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_JOB_H
