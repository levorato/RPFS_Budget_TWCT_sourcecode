//
// Created by mlevorato on 5/1/20.
//

#ifndef PFSP_ROBPFSINSTANCE_H
#define PFSP_ROBPFSINSTANCE_H


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/multi_array.hpp>

#include "../GlobalTypes.h"
#include "../PFSP_Parameters.h"

namespace robust {

    using namespace std;
    using namespace boost;
    using namespace parameters;

    typedef enum RobustType {
        adrs, budget
    } RobustType;

    class RobPFSInstance {
    public:
        /**
         * Build a Robust PFS Instance for the PFSP ADRS Problem, discrete scenarios variation (Kouvelis et al., 2000).
         * @param instance_name
         * @param lbd
         * @param M
         * @param N
         */
        RobPFSInstance(const string& instance_name, const unsigned &p_lambda, const unsigned &M, const unsigned &N,
                            const Robust_Parameters& p_rob_params) :
                name(instance_name), lambda(p_lambda), m(M), n(N), p_under(), p_over(), p_bar(), p_hat(),
                p_lambda(), alpha(0.0), alpha1(0.0), alpha2(0.0), beta(), rob_params(p_rob_params),
                robust_type() {

        }

        /**
         * Build a Robust PFS Instance for the PFSP ADRS Problem, continuous interval variation (Kouvelis et al., 2000).
         * @param instance_name
         * @param lbd
         * @param M
         * @param N
         */
        RobPFSInstance(const string& instance_name, const unsigned &M, const unsigned &N,
                            const Robust_Parameters& p_rob_params) :
                name(instance_name), lambda(0), m(M), n(N), p_under(M + 1, N + 1, 0.0), p_over(M + 1, N + 1, 0.0),
                p_lambda(), alpha(0.0), alpha1(0.0), alpha2(0.0), beta(), p_bar(), p_hat(), rob_params(p_rob_params),
                robust_type() {
        }

        /**
         * Build a Robust PFS Instance for the budget constraints problem variation (Ying, 2015).
         * @param instance_name
         * @param M
         * @param N
         * @param pp_bar
         * @param pp_hat
         */
        RobPFSInstance(const string& instance_name, const unsigned &M, const unsigned &N,
                const boost::numeric::ublas::matrix<double> &pp_bar, const boost::numeric::ublas::matrix<double> &pp_hat,
                const Robust_Parameters& p_rob_params) :
                name(instance_name), lambda(0), m(M), n(N), p_under(pp_bar - pp_hat), p_over(pp_bar + pp_hat),
                p_lambda(), p_bar(pp_bar),
                p_hat(pp_hat), alpha(0.0), alpha1(0.0), alpha2(0.0), beta(), rob_params(p_rob_params), robust_type()  {
        }
        ~RobPFSInstance() {   }

        const unsigned int& getNumberOfMachines() const {
            return m;
        }

        const unsigned int& getNumberOfJobs() const {
            return n;
        }

        const unsigned int& getNumberOfScenarios() const {
            return lambda;
        }

        bool isDiscrete() const {
            return lambda > 0;
        }

        const boost::numeric::ublas::matrix<double>& get_p_bar() const {
            return p_bar;
        }

        const boost::numeric::ublas::matrix<double>& get_p_hat() const {
            return p_hat;
        }

        // name which identified the instance
        string name;
        // number of jobs
        unsigned int n;
        // number of machines
        unsigned int m;
        // CASE A : DISCRETE SCENARIOS
        unsigned int lambda;
        // processing time matrix for discrete scenarios => indices start at 1 !!!
        std::vector< boost::numeric::ublas::matrix<double> > p_lambda;
        string beta;  // and alpha
        // CASE B : PROCESSING TIME INTERVALS MIN AND MAX VALUES
        boost::numeric::ublas::matrix<Time> p_over, p_under;
        double alpha1, alpha2;  // and beta
        // CASE C : PROCESSING TIME INTERVALS AVERAGE (p_bar) AND DEVIATION VALUES (p_hat)
        boost::numeric::ublas::matrix<double> p_bar, p_hat;
        double alpha;
        // enum representing the type of robust problem
        RobustType robust_type;
        // The set of robust parameters associated with the instance
        Robust_Parameters rob_params;
    };

}


#endif //PFSP_ROBPFSINSTANCE_H
