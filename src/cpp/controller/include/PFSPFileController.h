/*
 * PFSPFileController.h
 *
 *  Created on: Mar 12, 2018
 *      Author: mlevorato
 */

#ifndef SRC_CPP_CONTROLLER_INCLUDE_PFSPFILECONTROLLER_H_
#define SRC_CPP_CONTROLLER_INCLUDE_PFSPFILECONTROLLER_H_

#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <string>
#include <cstdio>
#include <sstream>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <glog/logging.h>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "../../problem/deterministic/PFSInstance.h"
#include "../../problem/PFSSolution.h"
#include "../../util/include/FileUtil.h"
#include "../../problem/GlobalTypes.h"
#include "../../util/include/SchedulingUtil.h"
#include "../../problem/robust/RobPFSInstance_Cmax.h"
#include "../../problem/robust/RobPFSInstance_WCT.h"
#include "../../problem/PFSP_Parameters.h"

using namespace boost;
using namespace boost::algorithm;
using namespace std;
namespace fs = boost::filesystem;
using namespace problem::common;

class PFSPFileController {
public:
    PFSPFileController();

    virtual ~PFSPFileController();

    static void writeSolutionToFile(boost::filesystem::path rootPath, problem::common::PFSSolution &solution,
            const double& timeSpent) {
        ostringstream oss;
        oss << "time spent: " << timeSpent << endl;
        oss << "solution status: " << (solution.is_optimal ? "Optimal" : "Not Optimal") << endl;
        oss << "size: " << solution.permutation_size << endl;
        oss << solution << endl;

        LOG(INFO) << oss.str();

        // Saves elapsed time and best solution to output file
        boost::filesystem::path result_file_path(rootPath);
        //result_file_path += "/result.txt";
        util::FileUtil::saveTextContentsToFile(result_file_path.string(), oss.str());
    }

    static PFSInstance readFromFile(const parameters::PFSP_Parameters &pfsp_params) {
        boost::filesystem::path filePath = pfsp_params.filePath;
        int m, n;
        std::ifstream infile(filePath.string().c_str(), std::ios::in | std::ios::binary);
        LOG(INFO) << "Reading input file: " << filePath.string();
        string fileId = filePath.filename().string();

        if (infile) {
            LOG(INFO) << "Reading input file line by line, avoiding full file reading.";
            // captures the first line of the file containing info about number of jobs and machines
            string line;
            std::getline(infile, line);

            char_separator<char> sep2(" \t");
            tokenizer<char_separator<char> > tokens(line, sep2);
            vector<string> vec;
            vec.assign(tokens.begin(), tokens.end());
            // n (jobs) m (machines)
            n = boost::lexical_cast<long>(vec.at(0));
            if (vec.size() == 1) { // only n is given => number of machines (m) is in the next line
                std::getline(infile, line);
                tokenizer<char_separator<char> > tokens2(line, sep2);
                vec.assign(tokens2.begin(), tokens2.end());
                m = boost::lexical_cast<long>(vec.at(0));
            } else {
                m = boost::lexical_cast<long>(vec.at(1));
            }
            double ub = std::numeric_limits<double>::max();
            double lb = 0.0;
            if (vec.size() >= 5) {
                ub = boost::lexical_cast<double>(vec.at(3));
                lb = boost::lexical_cast<double>(vec.at(4));
            }

            // read the matrix of processing times: M x J (machines x jobs)
            boost::numeric::ublas::matrix<double> matrix = boost::numeric::ublas::zero_matrix<double>(m + 1, n + 1);
            unsigned long a = 1;
            while (std::getline(infile, line)) {
                trim(line);
                tokenizer<char_separator<char> > tokens2(line, sep2);
                vec.assign(tokens2.begin(), tokens2.end());
                // cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

                for (unsigned long b = 1; b <= vec.size(); b++) {
                    try {
                        int value = 0;
                        sscanf(vec.at(b - 1).c_str(), "%d", &value);
                        // std::cout << "Adding job processing time (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
                        if (a > m || b > n) {
                            LOG(ERROR)
                                << "Error: invalid processing time. Either number of machines (m) or jobs (n) was exceeded.";
                            cerr
                                    << "Error: invalid processing time. Either number of machines (m) or jobs (n) was exceeded.";
                        } else {
                            matrix(a, b) = value;
                        }
                    } catch (boost::bad_lexical_cast const &) {
                        LOG(FATAL) << "Error: input string was not valid" << std::endl;
                        LOG(FATAL) << "Error: input string was not valid" << std::endl;
                        LOG(ERROR) << "vec.at(b) = " << vec.at(b - 1);
                    }
                }
                a++;
            }
            // LOG(INFO) << "Matrix read: \n" << matrix;
            infile.close();
            LOG(INFO) << "Successfully read input file.";
            LOG(INFO) << "Number of jobs (n) = " << n << ", number of machines (m) = " << m;
            PFSInstance instance(fileId, matrix);
            instance.lower_bound = lb;
            instance.upper_bound = ub;
            LOG(INFO) << "Upper bound = " << instance.upper_bound << "; Lower bound = "
                                    << instance.lower_bound;

            return instance;
        } else {
            LOG(ERROR) << "Failed to read input file!";
            throw "Failed to read input file!";
        }
    }

    static robust::RobPFSInstance_WCT readRobustInstance_newformat(const parameters::PFSP_Parameters &pfsp_params,
                                                                        parameters::Robust_Parameters &rob_params,
                                                                        const string &objective) {
        boost::filesystem::path filePath = pfsp_params.filePath;
        unsigned int m = 0, n = 0;
        std::ifstream infile(filePath.string().c_str(), std::ios::in | std::ios::binary);
        LOG(INFO) << "Reading robust input file (new format)): " << filePath.string();
        string fileId = filePath.filename().string();
        std::string filePath_str = filePath.string();
        boost::algorithm::to_lower(filePath_str);

        if (infile) {
            vector<string> vec;
            char_separator<char> sep2(" \t");
            string line;
            int step = 0;
            // captures the first line of the file containing info about number of jobs and machines
            while (std::getline(infile, line)) {  // For each job j
                trim(line);
                tokenizer<char_separator<char> > tokens(line, sep2);
                vec.assign(tokens.begin(), tokens.end());
                if (vec.size() == 0) continue;
                if ((line.find("!") != std::string::npos) ||
                    (line.find("#") != std::string::npos)) {  // Ignore comment lines
                    continue;
                }
                // n (jobs) m (machines)
                n = boost::lexical_cast<long>(vec.at(0));
                if (vec.size() == 1) { // only n is given => number of machines (m) is in the next line
                    std::getline(infile, line);
                    tokenizer<char_separator<char> > tokens2(line, sep2);
                    vec.assign(tokens2.begin(), tokens2.end());
                    m = boost::lexical_cast<long>(vec.at(0));
                } else {
                    m = boost::lexical_cast<long>(vec.at(1));
                }
                LOG(INFO) << "The number of machines is " << m << ".";
                LOG(INFO) << "The number of jobs is " << n << ".";
                break;
            }
            // vector of processing time weights
            boost::numeric::ublas::vector<double> w = boost::numeric::ublas::zero_vector<double>(n + 1);
            // read the matrix of processing times: M x J (machines x jobs)
            boost::numeric::ublas::matrix<double> P_bar = boost::numeric::ublas::zero_matrix<double>(m + 1,
                                                                                                     n + 1);
            boost::numeric::ublas::matrix<double> P_hat = boost::numeric::ublas::zero_matrix<double>(m + 1,
                                                                                                     n + 1);
            unsigned long i = 0, j = 1;
            step = 2;
            if(objective == "cmax") {
                step = 3;
            }
            while (std::getline(infile, line)) {  // For each job j
                trim(line);
                tokenizer<char_separator<char> > tokens2(line, sep2);
                vec.assign(tokens2.begin(), tokens2.end());
                if(vec.size() == 0)  continue;
                if ((line.find("!") != std::string::npos) || (line.find("#") != std::string::npos)) {  // Ignore comment lines
                    continue;
                }
                if(step == 2) {
                    // Read the vector of job weights (w)
                    double value = 0.0;
                    for(string x : vec) {
                        sscanf(x.c_str(), "%lf", &value);
                        w(j) = value;
                        j++;
                    }
                    if(j > n){  step = 3;  j = 1;  }
                } else if(step == 3) {
                    // Now read the processing time matrices
                    for (i = 0; i < m; i++) {  // For each machine i => get p_bar values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(i).c_str(), "%lf", &value);
                            P_bar(i + 1, j) = value;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(i) = " << vec.at(i);
                        }
                    }
                    j++;
                    if(j > n){  step = 4;  j = 1;  }
                } else if(step == 4) {
                    for (i = 0; i < m; i++) {  // For each machine i => get p_hat values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(i).c_str(), "%lf", &value);
                            P_hat(i + 1, j) = value;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(i) = " << vec.at(i);
                        }
                    }
                    j++;
                    if(j > n)  break;
                }
            }
            LOG(INFO) << "P_bar read: \n" << P_bar;
            LOG(INFO) << "P_hat read: \n" << P_hat;
            LOG(INFO) << "Job weight vector : " << w;
            infile.close();
            LOG(INFO) << "Successfully read (new format) robust input file.";
            robust::RobPFSInstance_WCT instance(fileId, m, n, P_bar, P_hat, w, rob_params);
            // determine the value of the alpha parameter (from file path)
            string alpha_str = filePath.string();
            alpha_str = alpha_str.substr(alpha_str.rfind(boost::filesystem::path::preferred_separator) + 1);
            int idx = alpha_str.find("_wct_inputs");
            if(idx == string::npos) {
                idx = alpha_str.find("_cmax_inputs");
            }
            alpha_str = alpha_str.substr(idx - 2, 2);
            double a = 0.0;
            sscanf(alpha_str.c_str(), "%lf", &a);
            instance.alpha = a;
            instance.robust_type = robust::RobustType::budget;
            LOG(INFO) << "Alpha value: " << a;
            return instance;
        } else {
            LOG(ERROR) << "Failed to read robust input file!";
            throw "Failed to read robust input file!";
        }
    }

    static robust::RobPFSInstance_Cmax readRobustInstance_Cmax_FromFile(const parameters::PFSP_Parameters &pfsp_params,
            parameters::Robust_Parameters &rob_params, const string &file_type) {
        boost::filesystem::path filePath = pfsp_params.filePath;
        unsigned m = 0, n = 0;
        std::ifstream infile(filePath.string().c_str(), std::ios::in | std::ios::binary);
        LOG(INFO) << "Reading robust Cmax input file: " << filePath.string();
        string fileId = filePath.filename().string();
        std::string filePath_str = filePath.string();
        boost::algorithm::to_lower(filePath_str);

        if (infile) {
            LOG(INFO) << "Reading input file line by line, avoiding full file reading.";
            LOG(INFO) << "Input file type is " << file_type;
            vector<string> vec;
            char_separator<char> sep2(" \t");
            string line;
            if (file_type == "ying") {  // Ying Robust PFSP Budget instance file
                // read the matrix of processing times: M x J (machines x jobs)
                boost::numeric::ublas::matrix<double> P_bar = boost::numeric::ublas::zero_matrix<double>(MACHINES + 1, NJOBS + 1);
                boost::numeric::ublas::matrix<double> P_hat = boost::numeric::ublas::zero_matrix<double>(MACHINES + 1, NJOBS + 1);
                unsigned long i = 0, j = 0, read_m_n = 0;
                unsigned int m = 0, n = 0;
                while (std::getline(infile, line)) {  // For each job j
                    trim(line);
                    tokenizer<char_separator<char> > tokens2(line, sep2);
                    vec.assign(tokens2.begin(), tokens2.end());
                    if(vec.size() == 0)  continue;
                    if (line.find("!") != std::string::npos) {  // Ignore comment lines
                        read_m_n = 1;
                        continue;
                    }
                    j++;
                    if(j == 1) {  // use the first line to detect the number of machines
                        if(read_m_n == 1) { // ! n m alpha
                            sscanf(vec.at(0).c_str(), "%u", &n);
                            sscanf(vec.at(1).c_str(), "%u", &m);
                            LOG(INFO) << "The number of jobs is " << n << ".";
                            LOG(INFO) << "The number of machines is " << m << ".";
                            continue;
                        } else {
                            m = (int) vec.size() / 2;
                            LOG(INFO) << "The number of machines is " << m << ".";
                        }
                    }
                    if(read_m_n == 1) {
                        j = 1;
                        read_m_n = 2;
                    }
                    // cout << "Line is: " << line << " vec.size = " << vec.size() << endl;
                    // format : 13	16	1.3	1.6  ==> p_bar_1j p_bar_2j p_hat_1j p_hat_2j
                    unsigned long count = 0;
                    for (i = 0; i < m; i++) {  // For each machine i => get p_bar values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(count).c_str(), "%lf", &value);
                            P_bar(i + 1, j) = value;
                            count++;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(count) = " << vec.at(count);
                        }
                    }
                    for (i = 0; i < m; i++) {  // For each machine i => get p_hat values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(count).c_str(), "%lf", &value);
                            P_hat(i + 1, j) = value;
                            count++;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(count) = " << vec.at(count);
                        }
                    }
                }
                n = j;
                // Reduce matrices size (according to the value of n)
                P_bar.resize(m + 1, n + 1, true);
                P_hat.resize(m + 1, n + 1, true);
                // LOG(INFO) << "P_bar read: \n" << P_bar;
                // LOG(INFO) << "P_hat read: \n" << P_hat;
                infile.close();
                LOG(INFO) << "Number of jobs (n) = " << n << ", number of machines (m) = " << m;
                LOG(INFO) << "Successfully read Ying robust input file.";
                robust::RobPFSInstance_Cmax instance(fileId, m, n, P_bar, P_hat, rob_params);
                // determine the value of the alpha parameter (from file path)
                string alpha_str = filePath.string();
                alpha_str = alpha_str.substr(0, alpha_str.rfind(boost::filesystem::path::preferred_separator));
                alpha_str = alpha_str.substr(alpha_str.rfind(boost::filesystem::path::preferred_separator)+1);
                alpha_str = alpha_str.substr(0, alpha_str.size() - 1);  // remove the \% at the end of the string
                double a = 0.0;
                sscanf(alpha_str.c_str(), "%lf", &a);
                instance.alpha = a;
                instance.robust_type = robust::RobustType::budget;
                return instance;
            } else if (file_type == "kouvelis") {  // Kouvelis Robust PFSP instance
                double alpha, alpha1, alpha2;
                string beta;
                do {
                    // captures the first line of the file containing info about number of jobs and machines
                    std::getline(infile, line);
                    tokenizer<char_separator<char> > tokens(line, sep2);
                    vec.assign(tokens.begin(), tokens.end());
                    // if(vec.size() > 0)  if(vec[0].find("#") != string::npos)  continue;  // ignore line comments
                } while(vec.size() == 0 || vec[0].find("#") != string::npos);  // ignore empty lines
                // n (jobs) m (machines)
                m = boost::lexical_cast<long>(vec.at(0));
                n = boost::lexical_cast<long>(vec.at(1));
                long lambda = boost::lexical_cast<long>(vec.at(3));
                LOG(INFO) << "Number of jobs (n) = " << n << ", number of machines (m) = " << m;
                if(lambda == 0) {  // processing time intervals
                    LOG(INFO) << "Processing time intervals instance.";
                    alpha1 = boost::lexical_cast<double>(vec.at(4));
                    alpha2 = boost::lexical_cast<double>(vec.at(5));
                    beta = vec.at(6);
                    boost::numeric::ublas::matrix<double> P_low = boost::numeric::ublas::zero_matrix<double>(m + 1, n + 1);
                    boost::numeric::ublas::matrix<double> P_high = boost::numeric::ublas::zero_matrix<double>(m + 1, n + 1);
                    // read the matrix of processing times: M x J (machines x jobs)
                    unsigned long i = 0;
                    while (std::getline(infile, line)) {  // For each machine i
                        trim(line);
                        tokenizer<char_separator<char> > tokens2(line, sep2);
                        vec.assign(tokens2.begin(), tokens2.end());
                        if(vec.size() == 0)  continue;  // ignore empty lines
                        for (int b = 0; b < vec.size(); b++) {  // For each job b
                            try {
                                double value = 0.0;
                                sscanf(vec.at(b).c_str(), "%lf", &value);
                                if(i < m) {  // lower bound values
                                    P_low(i + 1, b + 1) = value;
                                } else {  // upper bound values
                                    P_high(i - m + 1, b + 1) = value;
                                }
                            } catch (boost::bad_lexical_cast const &) {
                                LOG(FATAL) << "Error: input string was not valid" << std::endl;
                                LOG(ERROR) << "vec.at(b) = " << vec.at(b);
                            }
                        }
                        i++;
                    }
                    infile.close();
                    LOG(INFO) << "Successfully read Kouvelis time intervals robust input file.";
                    robust::RobPFSInstance_Cmax instance(fileId, m, n, rob_params);
                    //LOG(INFO) << "P_low read: \n" << P_low;
                    //LOG(INFO) << "P_high read: \n" << P_high;
                    instance.p_under = P_low;
                    instance.p_over = P_high;
                    instance.lambda = lambda;
                    instance.alpha1 = alpha1;
                    instance.alpha2 = alpha2;
                    instance.beta = beta;
                    instance.robust_type = robust::RobustType::adrs;
                    return instance;
                } else {  // discrete scenarios
                    LOG(INFO) << "Discrete time scenarios instance, number of scenarios = " << lambda;
                    alpha = boost::lexical_cast<double>(vec.at(4));
                    beta = vec.at(5);
                    robust::RobPFSInstance_Cmax instance(fileId, m, n, rob_params);
                    instance.p_lambda = std::vector< boost::numeric::ublas::matrix<double> >(lambda,
                            boost::numeric::ublas::zero_matrix<double>(m + 1, n + 1));
                    // read the matrix of processing times: M x J (machines x jobs)
                    unsigned long i = 0, lambda_count = 0;
                    while (std::getline(infile, line)) {  // For each machine i
                        trim(line);
                        tokenizer<char_separator<char> > tokens2(line, sep2);
                        vec.assign(tokens2.begin(), tokens2.end());
                        if(vec.size() == 0)  continue;  // ignore empty lines
                        if(vec[0].find("Scenario") != string::npos)  continue;
                        if(i == m) {
                            i = 0;
                            lambda_count++;
                        }
                        for (int b = 0; b < vec.size(); b++) {  // For each job b
                            try {
                                double value = 0.0;
                                sscanf(vec.at(b).c_str(), "%lf", &value);
                                instance.p_lambda[lambda_count](i + 1, b + 1) = value;
                            } catch (boost::bad_lexical_cast const &) {
                                LOG(FATAL) << "Error: input string was not valid" << std::endl;
                                LOG(ERROR) << "vec.at(b) = " << vec.at(b);
                            }
                        }
                        i++;
                    }
                    if(lambda_count + 1 != lambda) {
                        LOG(ERROR) << "Error reading instance: lambda count does not match!";
                    }
                    //for(boost::numeric::ublas::matrix<double> p : instance.p_lambda) {
                    //    LOG(INFO) << "P read: \n" << p;
                    //}
                    infile.close();
                    LOG(INFO) << "Successfully read Kouvelis discrete scenarios robust input file.";
                    instance.lambda = lambda;
                    instance.alpha = alpha;
                    instance.beta = beta;
                    instance.robust_type = robust::RobustType::adrs;
                    return instance;
                }
            } else {  // "newformat"
                robust::RobPFSInstance_WCT instance = readRobustInstance_newformat(pfsp_params, rob_params, "cmax");
                return robust::RobPFSInstance_Cmax(fileId, instance.m, instance.n, instance.p_bar,
                                                   instance.p_hat, rob_params);
            }
            // captures the first line of the file containing info about number of jobs and machines
            std::getline(infile, line);
            tokenizer<char_separator<char> > tokens(line, sep2);
            vec.assign(tokens.begin(), tokens.end());
            // n (jobs) m (machines)
            n = boost::lexical_cast<long>(vec.at(0));
            if (vec.size() == 1) { // only n is given => number of machines (m) is in the next line
                std::getline(infile, line);
                tokenizer<char_separator<char> > tokens2(line, sep2);
                vec.assign(tokens2.begin(), tokens2.end());
                m = boost::lexical_cast<long>(vec.at(0));
            } else {
                m = boost::lexical_cast<long>(vec.at(1));
            }
            double ub = std::numeric_limits<double>::max();
            double lb = 0.0;
            if (vec.size() >= 5) {
                ub = boost::lexical_cast<double>(vec.at(3));
                lb = boost::lexical_cast<double>(vec.at(4));
            }

            // read the matrix of processing times: M x J (machines x jobs)
            boost::numeric::ublas::matrix<double> matrix = boost::numeric::ublas::zero_matrix<double>(m + 1, n + 1);
            unsigned long a = 1;
            while (std::getline(infile, line)) {
                trim(line);
                tokenizer<char_separator<char> > tokens2(line, sep2);
                vec.assign(tokens2.begin(), tokens2.end());
                // cout << "Line is: " << line << " vec.size = " << vec.size() << endl;

                for (unsigned long b = 1; b <= vec.size(); b++) {
                    try {
                        int value = 0;
                        sscanf(vec.at(b - 1).c_str(), "%d", &value);
                        // std::cout << "Adding job processing time (" << a-1 << ", " << b-1 << ") = " << value << std::endl;
                        if (a > m || b > n) {
                            LOG(ERROR)
                                << "Error: invalid processing time. Either number of machines (m) or jobs (n) was exceeded.";
                            cerr
                                    << "Error: invalid processing time. Either number of machines (m) or jobs (n) was exceeded.";
                        } else {
                            matrix(a, b) = value;
                        }
                    } catch (boost::bad_lexical_cast const &) {
                        LOG(FATAL) << "Error: input string was not valid" << std::endl;
                        LOG(FATAL) << "Error: input string was not valid" << std::endl;
                        LOG(ERROR) << "vec.at(b) = " << vec.at(b - 1);
                    }
                }
                a++;
            }
            // LOG(INFO) << "Matrix read: \n" << matrix;
            infile.close();
            LOG(INFO) << "Successfully read input file.";
            LOG(INFO) << "Number of jobs (n) = " << n << ", number of machines (m) = " << m;
            robust::RobPFSInstance_Cmax instance(fileId, m, n, rob_params);
            instance.robust_type = robust::RobustType::budget;
            return instance;
        } else {
            LOG(ERROR) << "Failed to read robust input file!";
            throw "Failed to read robust input file!";
        }
    }

    static robust::RobPFSInstance_WCT readRobustInstance_WCT_FromFile(const parameters::PFSP_Parameters &pfsp_params,
                    parameters::Robust_Parameters &rob_params,
                    const string &file_type) {
        boost::filesystem::path filePath = pfsp_params.filePath;
        unsigned int m = 0, n = 0;
        std::ifstream infile(filePath.string().c_str(), std::ios::in | std::ios::binary);
        LOG(INFO) << "Reading robust WCT input file: " << filePath.string();
        string fileId = filePath.filename().string();
        std::string filePath_str = filePath.string();
        boost::algorithm::to_lower(filePath_str);

        if (infile) {
            LOG(INFO) << "Reading input file line by line, avoiding full file reading.";
            vector<string> vec;
            char_separator<char> sep2(" \t");
            string line;
            if (file_type == "ying") {  // Ying Robust PFSP Budget instance file
                // vector of processing time weights
                boost::numeric::ublas::vector<double> w = boost::numeric::ublas::zero_vector<double>(NJOBS + 1);
                // read the matrix of processing times: M x J (machines x jobs)
                boost::numeric::ublas::matrix<double> P_bar = boost::numeric::ublas::zero_matrix<double>(MACHINES + 1,
                        NJOBS + 1);
                boost::numeric::ublas::matrix<double> P_hat = boost::numeric::ublas::zero_matrix<double>(MACHINES + 1,
                        NJOBS + 1);
                unsigned long i = 0, j = 0;
                m = 0;  n = 0;
                while (std::getline(infile, line)) {  // For each job j
                    trim(line);
                    tokenizer<char_separator<char> > tokens2(line, sep2);
                    vec.assign(tokens2.begin(), tokens2.end());
                    if(vec.size() == 0)  continue;
                    if (line.find("!") != std::string::npos) {  // Ignore comment lines
                        continue;
                    }
                    if(j == 0) {  // use the first line to detect the number of machines
                        // n (jobs) m (machines)
                        n = boost::lexical_cast<long>(vec.at(0));
                        if (vec.size() == 1) { // only n is given => number of machines (m) is in the next line
                            std::getline(infile, line);
                            tokenizer<char_separator<char> > tokens2(line, sep2);
                            vec.assign(tokens2.begin(), tokens2.end());
                            m = boost::lexical_cast<long>(vec.at(0));
                        } else {
                            m = boost::lexical_cast<long>(vec.at(1));
                        }
                        LOG(INFO) << "The number of machines is " << m << ".";
                        LOG(INFO) << "The number of jobs is " << n << ".";
                        // Read the vector of job weights (w)
                        std::getline(infile, line);
                        trim(line);
                        tokenizer<char_separator<char> > tokens3(line, sep2);
                        vec.assign(tokens3.begin(), tokens3.end());
                        double value = 0.0;
                        for(int i = 0; i < n; i++) {
                            sscanf(vec.at(i).c_str(), "%lf", &value);
                            w(i + 1) = value;
                        }
                        j++;
                        continue;
                    }
                    // Now read the processing time matrices
                    // cout << "Line is: " << line << " vec.size = " << vec.size() << endl;
                    // format : 13	16	1.3	1.6  ==> p_bar_1j p_bar_2j p_hat_1j p_hat_2j
                    unsigned long count = 0;
                    for (i = 0; i < m; i++) {  // For each machine i => get p_bar values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(count).c_str(), "%lf", &value);
                            P_bar(i + 1, j) = value;
                            count++;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(count) = " << vec.at(count);
                        }
                    }
                    for (i = 0; i < m; i++) {  // For each machine i => get p_hat values
                        try {
                            double value = 0.0;
                            sscanf(vec.at(count).c_str(), "%lf", &value);
                            P_hat(i + 1, j) = value;
                            count++;
                        } catch (boost::bad_lexical_cast const &) {
                            LOG(FATAL) << "Error: input string was not valid" << std::endl;
                            LOG(ERROR) << "vec.at(count) = " << vec.at(count);
                        }
                    }
                    j++;
                }
                // Reduce matrices size (according to the value of n)
                P_bar.resize(m + 1, n + 1, true);
                P_hat.resize(m + 1, n + 1, true);
                w.resize(n + 1, true);
                 LOG(INFO) << "P_bar read: \n" << P_bar;
                 LOG(INFO) << "P_hat read: \n" << P_hat;
                 LOG(INFO) << "Job weight vector : " << w;
                infile.close();
                LOG(INFO) << "Number of jobs (n) = " << n << ", number of machines (m) = " << m;
                LOG(INFO) << "Successfully read Ying WCT robust input file.";
                robust::RobPFSInstance_WCT instance(fileId, m, n, P_bar, P_hat, w, rob_params);
                // determine the value of the alpha parameter (from file path)
                string alpha_str = filePath.string();
                alpha_str = alpha_str.substr(0, alpha_str.rfind(boost::filesystem::path::preferred_separator));
                alpha_str = alpha_str.substr(alpha_str.rfind(boost::filesystem::path::preferred_separator)+1);
                alpha_str = alpha_str.substr(0, alpha_str.size() - 1);  // remove the \% at the end of the string
                double a = 0.0;
                sscanf(alpha_str.c_str(), "%lf", &a);
                instance.alpha = a;
                instance.robust_type = robust::RobustType::budget;
                return instance;
            } else {  // "newformat"
                return readRobustInstance_newformat(pfsp_params, rob_params, "wct");
            }
        } else {
            LOG(ERROR) << "Failed to read robust WCT input file!";
            throw "Failed to read robust WCT input file!";
        }
    }
};


#endif /* SRC_CPP_CONTROLLER_INCLUDE_PFSPFILECONTROLLER_H_ */
