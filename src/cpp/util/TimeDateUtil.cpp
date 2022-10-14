/*
 * TimeDateUtil.cpp
 *
 *  Created on: 23/05/2013
 *      Author: Mario Levorato
 */

#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>
#include "include/TimeDateUtil.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace util {

TimeDateUtil::TimeDateUtil() {
	// TODO Auto-generated constructor stub

}

TimeDateUtil::~TimeDateUtil() {
	// TODO Auto-generated destructor stub
}

std::string TimeDateUtil::FormatTime(boost::posix_time::ptime now)
{
  using namespace boost::posix_time;
  static std::locale loc(std::wcout.getloc(),
                         new wtime_facet(L"%Y%m%d_%H%M%S"));

  std::basic_stringstream<wchar_t> wss;
  wss.imbue(loc);
  wss << now;

  std::wstring ws = wss.str();

  return std::string ( ws.begin(), ws.end() );
}

std::string TimeDateUtil::generateRandomId() {
	using namespace boost::posix_time;
	ptime now = second_clock::local_time();

	std::string time = TimeDateUtil::FormatTime(now);
	boost::uuids::uuid tag = boost::uuids::random_generator()();
	const std::string tmp = boost::lexical_cast<std::string>(tag);
	return time + "_" + tmp;
}

std::string TimeDateUtil::getDateAndTime() {
	using namespace boost::posix_time;
	ptime now = second_clock::local_time();

	return TimeDateUtil::FormatTime(now);
}

double TimeDateUtil::calculate_time_spent(boost::timer::cpu_timer &timer) {
    sec seconds = boost::chrono::nanoseconds(timer.elapsed().wall);
    double time_spent = ((double)seconds.count());
    return time_spent;
}

double TimeDateUtil::calculate_time_spent_ns(boost::timer::cpu_timer &timer) {
    auto nanoseconds = boost::chrono::nanoseconds(timer.elapsed().wall);
    double time_spent = ((double)nanoseconds.count());
    return time_spent;
}

double TimeDateUtil::get_time_spent_in_seconds(boost::timer::cpu_times elapsed_times) {
    sec seconds = boost::chrono::nanoseconds(elapsed_times.wall);
    double time_spent = ((double)seconds.count());
    return time_spent;
}

} /* namespace util */
