//
// Created by mlevorato on 11/28/19.
//

#ifndef FLOWSHOP_SOLVER_GRASP_OUTPUTS_H
#define FLOWSHOP_SOLVER_GRASP_OUTPUTS_H

#include <string>
#include <sstream>
#include <glog/logging.h>
#include "FileUtil.h"
#include "SchedulingSolution.h"

namespace problem {
namespace pfsp {

    using namespace std;
    using namespace util;
    using namespace boost;

    template<class ProblemInstance>
    class Outputs {
    private:
        SchedulingSolution<ProblemInstance> nehSol; // NEH solution
        SchedulingSolution<ProblemInstance> ourBestSol; // Best solution provided by GRASP
        long iterations;  // number of performed GRASP iterations
        long num_visited_solutions;  // number of solutions examined by GRASP (given all neighborhoods)
        string time_results;
        double time_spent;  // total time spent in the procedure
        long num_improvements;  // number of improvements in the solution value
        // Constructive phase statistics
        std::vector<double> alpha_history;
        std::vector< std::vector<double> > alpha_quality, alpha_probability;

    public:
        Outputs(SchedulingSolution<ProblemInstance> neh, SchedulingSolution<ProblemInstance> obs,
                const long &p_iterations, const long &p_visited_sols, const string &p_time_results,
                const double &p_time_spent, const long &p_num_improvements, const std::vector<double> &alpha_hist,
                const std::vector< std::vector<double> > &alpha_qual,
                const std::vector< std::vector<double> > &alpha_prob) :
                nehSol(neh), time_spent(p_time_spent),
                ourBestSol(obs), num_improvements(p_num_improvements),
                iterations(p_iterations), num_visited_solutions(p_visited_sols), time_results(p_time_results),
                alpha_history(alpha_hist), alpha_quality(alpha_qual), alpha_probability(alpha_prob)
        {
        }

        SchedulingSolution<ProblemInstance> getOurBestSol() const {
            return ourBestSol;
        }

        SchedulingSolution<ProblemInstance> getNehSol() const {
            return nehSol;
        }

        long get_number_iterations() const {
            return iterations;
        }

        long get_num_visited_solutions() const {
            return num_visited_solutions;
        }

        string get_time_results() const {
            return time_results;
        }

        double get_time_spent() const {
            return time_spent;
        }

        long get_number_improvements() const {
            return num_improvements;
        }

        const vector<double> &getAlphaHistory() const {
            return alpha_history;
        }

        const std::vector< std::vector<double> > &getAlphaQuality() const {
            return alpha_quality;
        }

        const std::vector< std::vector<double> > &getAlphaProbability() const {
            return alpha_probability;
        }

        void sendToFile(const string &filePath)
        {
            stringstream ss;
            ss << "***************************************************";
            ss << "*       RESULTS FROM PFSP GRASP                   *";
            ss << "***************************************************";

            ss << "\r\n";
            ss << "--------------------------------------------";
            ss << "               NEH SchedulingSolution                 ";
            ss << "--------------------------------------------";
            ss << nehSol.toString(true);

            ss << "\r\n";
            ss << "--------------------------------------------";
            ss << "Our best solution (provided by the GRASP)   ";
            ss << "--------------------------------------------";
            ss << ourBestSol.toString(true);

            LOG(INFO) << "[GRASP PFSP] Saving results file to " << filePath;
            util::FileUtil::saveTextContentsToFile(filePath, ss.str(), false);
        }
    };
}
}

#endif //FLOWSHOP_SOLVER_OUTPUTS_H

