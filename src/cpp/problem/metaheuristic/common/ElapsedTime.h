//
// Created by mlevorato on 11/27/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_ELAPSEDTIME_H
#define FLOWSHOP_SOLVER_GRASP_ELAPSEDTIME_H

#include <cmath>
#include <sstream>

namespace problem {
namespace pfsp {

    class ElapsedTime {
    public:
        ElapsedTime() {
        }

        static double calcElapsed(const double &start, const double &end) {
            double elapsed = (end - start) / 1.0e+9;
            return elapsed;
        }

        static std::string calcElapsedHMS(const long &start, const long &end) {
            std::stringstream ss;
            double elapsed = (end - start) / 1.0e+9;
            ss << calcHMS((int) round(elapsed));
            return ss.str();
        }

        static std::string calcHMS(int timeInSeconds) {
            std::stringstream ss;

            int hours, minutes, seconds;
            hours = timeInSeconds / 3600;
            timeInSeconds = timeInSeconds - (hours * 3600);
            minutes = timeInSeconds / 60;
            timeInSeconds = timeInSeconds - (minutes * 60);
            seconds = timeInSeconds;

            ss << hours << "h " << minutes << "m " << seconds << "s";

            return ss.str();
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_ELAPSEDTIME_H
