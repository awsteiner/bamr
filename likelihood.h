#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/root_brent_gsl.h>
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_math.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace std::placeholders;

typedef boost::numeric::ublas::vector<double> ubvector;

struct mass_data {
  // std::string to store the full names of stars
  vector<string> name_ns;
  vector<string> name_wd;
  vector<string> name_ms;

  // std::string to store the shortened names of stars
  vector<string> id_ns;
  vector<string> id_wd;
  vector<string> id_ms;

  //std::vector to store the measured NS mass
  vector<double> mass_ns;
  vector<double> mass_wd;
  vector<double> mass_ms;

  //std::vector to store +68% central limits of NS masses
  vector<double> uplim_ns;
  vector<double> uplim_wd;
  vector<double> lim_ms;

  //std::vector to store -68% central limits of NS masses
  vector<double> lowlim_ns;
  vector<double> lowlim_wd;
};

class like {

  private:
    double norm_pdf(double);
    double norm_cdf(double);
    double skew_norm(double, double, double, double);
    double asym_norm(double, double, double);
    double f2solve(double, double &, double &);
    double df2solve(double, double &, double &);
    double calc_par_d(double, double);
    double calc_likelihood_ns(const ubvector &, vec_index &);
    double calc_likelihood_wd(const ubvector &, vec_index &);
    double calc_likelihood_ms(const ubvector &, vec_index &);

  public:
    // Total number of stars in each population
    static const int N_ns=22, N_wd=32, N_ms=16;

    // Set tolerance for small likelihoods
    const double tol = 1.0e-6;

    // std::string vectors to store the full names/identifiers of the stars
    string nsns_id[N_ns],  nswd_id[N_wd],  nsms_id[N_ms];

    // std::string vectors to store the shortened names/identifiers
    string nsns_sid[N_ns], nswd_sid[N_wd], nsms_sid[N_ms];

    /* Arrays to store the raw data in the format: [mi][li, ui] for NS-NS
    and NS-WD, and [mi][si] for NS-MS (where si=li=ui) */
    double nsns[N_ns][3], nswd[N_wd][3], nsms[N_ms][2];

    void load_data();
    void set_params(vec_index &);
    double calc_likelihood(const ubvector &, vec_index &);
};

#endif
