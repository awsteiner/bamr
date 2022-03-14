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

  // std::vector to store the measured NS mass
  vector<double> mass_ns;
  vector<double> mass_wd;
  vector<double> mass_ms;

  // std::vector to store +68% central limits of NS masses
  vector<double> uplim_ns;
  vector<double> uplim_wd;
  vector<double> lim_ms; // Symmetric 68% limits for NS-MS

  // std::vector to store -68% central limits of NS masses
  vector<double> lowlim_ns;
  vector<double> lowlim_wd;

  // Total number of stars in all populations
  size_t n_stars;

  // Function to load population mass data
  void load_data();

};

class likelihood {

  private:

    // PDF of standard normal distribution N(0,1)
    double norm_pdf(double);

    // CDF of standard N(0,1) in terms of erf(x)
    double norm_cdf(double);

    // Skewed Normal PDF [eq. 13, Kiziltan et al. (2013)]
    double skew_norm(double, double, double, double);

    // Asymmetric Normal PDF [eq. 14, Kiziltan et al. (2013)]
    double asym_norm(double, double, double);

    // The function to solve f(x)=0 [see refs/calc.pdf]
    double f2solve(double, double &, double &);
    
    /* Derivative of the function to solve f'(x)
    (for use with root_stef only) */
    double df2solve(double, double &, double &);
    
    /* Solver to calculate scale parameters d for a given 
    asymmetry parameter c of function asym_norm (AN) */
    double calc_par_d(double, double);


  public:
  
    /* Tolerance for small weights, below which weights are
    ignored. This ensures that log-weights do not explode. */
    const double tol = 1.0e-6;
    
    // Function to fill pvi object with all parameter names
    void set_params(vec_index &);

    // Likelihood function for NS-NS
    double get_weight_ns(const ubvector &, vec_index &);
    
    // Likelihood function for NS-WD
    double get_weight_wd(const ubvector &, vec_index &);
    
    // Likelihood function for NS-MS
    double get_weight_ms(const ubvector &, vec_index &);
    
    /* Combined likelihood function for all stars, except 
    GW170817, QLMXBs, PREs, and NICER */
    double get_weight(const ubvector &, vec_index &);

};

#endif
