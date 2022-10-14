//============================================================================
// Name        : CommandLineInterfaceController.h
// Author      : Mario Levorato
// Version     :
// Copyright   : Copyright (c) 2013
// Description : Command Line Interface (CLI) class. Processes program
//               arguments and invokes problem solving.
//============================================================================

#ifndef COMMANDLINEINTERFACECONTROLLER_H_
#define COMMANDLINEINTERFACECONTROLLER_H_

#include <string>
#include <boost/filesystem.hpp>
#include <boost/exception/all.hpp>
#include <boost/program_options.hpp>
#include <boost/any.hpp>

#include "../../problem/PFSP_Parameters.h"

namespace controller {

    using namespace std;
    using namespace parameters;
    namespace fs = boost::filesystem;

    typedef boost::error_info<struct tag_stack_str, std::string> stack_info;

    class CommandLineInterfaceController {
    public:
        CommandLineInterfaceController();

        virtual ~CommandLineInterfaceController();

        string getTimeAndDateAsString();

        int processArgumentsAndExecute(int argc, char *argv[]);

    private:
        void processInputFile(const PFSP_Parameters &pfsp_params,
                                               const BranchBound_Parameters &bb_params,
                                               Robust_Parameters &rob_params,
                                               const UpperBound_Parameters &ub_params,
                                               const bool &robust);

        void PrintVariableMap(const boost::program_options::variables_map vm);

        /*
         * Tests if a specified file exists.
         */
        bool testInputFile(fs::path filePath);

        void readPropertiesFile();

        void initLogging(const string &outputFolder, const string &executionId, char *argv[]);

        static void handler();

        string logSeverity;
    };
} /* namespace controller */

#endif /* COMMANDLINEINTERFACECONTROLLER_H_ */
