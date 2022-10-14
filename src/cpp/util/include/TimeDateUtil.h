/*
 * TimeDateUtil.h
 *
 *  Created on: 23/05/2013
 *      Author: Mario Levorato
 */

#ifndef TIMEDATEUTIL_H_
#define TIMEDATEUTIL_H_

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/timer/timer.hpp>
#include <boost/chrono.hpp>

namespace util {

    using boost::timer::cpu_timer;
    using boost::timer::cpu_times;
    using boost::timer::nanosecond_type;
    typedef boost::chrono::duration<double> sec;

class TimeDateUtil {
public:
	TimeDateUtil();
	virtual ~TimeDateUtil();

	static std::string FormatTime(boost::posix_time::ptime now);
	static std::string generateRandomId();
	static std::string getDateAndTime();
	static double calculate_time_spent(boost::timer::cpu_timer &timer);
	static double get_time_spent_in_seconds(boost::timer::cpu_times elapsed_times);
    static double calculate_time_spent_ns(boost::timer::cpu_timer &timer);
};

} /* namespace util */
#endif /* TIMEDATEUTIL_H_ */
