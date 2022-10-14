/*
 * ExecutionInfo.h
 *
 *  Created on: 24/04/2014
 *      Author: czt0
 */

#ifndef EXECUTIONINFO_H_
#define EXECUTIONINFO_H_

#include <string>

namespace util {

class ExecutionInfo {
public:
	ExecutionInfo(std::string id, std::string fileid, std::string outfolder);
	virtual ~ExecutionInfo();

	std::string executionId;
	std::string fileId;
	std::string outputFolder;
};

} /* namespace util */
#endif /* EXECUTIONINFO_H_ */
