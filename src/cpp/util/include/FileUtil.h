/*
 * FileUtil.h
 *
 *  Created on: Apr 3, 2016
 *      Author: mlevorato
 */

#ifndef SRC_UTIL_INCLUDE_FILEUTIL_H_
#define SRC_UTIL_INCLUDE_FILEUTIL_H_

#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdarg>

#include <glog/logging.h>
#include <boost/filesystem.hpp>

namespace util {

    using namespace std;

    class FileUtil {
    public:
        FileUtil();

        virtual ~FileUtil();

        static bool saveTextContentsToFile(const string &file_path, const string &file_contents, bool append = true) {
            boost::filesystem::path filepath(file_path);
            LOG(INFO) << "Saving output text file to '" << filepath.string() << "'.";

            ofstream out(filepath.string().c_str(), append ? ios::out | ios::app : ios::out);
            if (!out) {
                LOG(FATAL) << "Cannot open output result file to: " << filepath.string();
                filepath += ".txt";
                out = ofstream(filepath.string().c_str(), append ? ios::out | ios::app : ios::out);
                if(!out)  return false;
            }
            out << file_contents;
            // Closes the file
            out.close();
            return true;
        }

        static string create_output_folder(const string &outputFolder, const string &outputsolfname = "") {
            boost::filesystem::path rootPath(outputFolder);
            boost::system::error_code returnedError;
            boost::filesystem::create_directories(rootPath, returnedError);
            boost::filesystem::path base_dir = outputFolder;
            if (outputsolfname != "")
                base_dir /= outputsolfname;
            boost::filesystem::create_directories(base_dir, returnedError);
            if (returnedError) {
                LOG(FATAL) << "Cannot create output directory: " << rootPath.string();
                throw std::runtime_error("Cannot create output directory.");
            }
            return base_dir.string();
        }


        /** prints a message */
        static void TestinfoMessage(
                FILE *file,               /**< file stream to print into, or NULL for stdout */
                const char *formatstr,          /**< format string like in printf() function */
                ...                  /**< variable argument list -> format arguments line in printf() function */
        ) {
            // TestmessageVFPrintInfo(scip->messagehdlr, file, formatstr, ap);
            va_list args;
            va_start (args, formatstr);
            vfprintf(file, formatstr, args);
            va_end (args);
        }

        /** prints a message */
        static void TestmessageFPrintInfo(
                FILE *file,               /**< file stream to print into, or NULL for stdout */
                const char *formatstr,          /**< format string like in printf() function */
                ...                 /**< variable argument list -> format arguments line in printf() function */
        ) {
            // TestmessageVFPrintInfo(scip->messagehdlr, file, formatstr, ap);
            va_list args;
            va_start (args, formatstr);
            vfprintf(file, formatstr, args);
            va_end (args);
        }
    };

} /* namespace util */

#endif /* SRC_UTIL_INCLUDE_FILEUTIL_H_ */
