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
  vector<string> name_lms;
  vector<string> name_hms;

  // std::string to store the shortened names of stars
  vector<string> id_ns;
  vector<string> id_wd;
  vector<string> id_lms;
  vector<string> id_hms;

  // std::vector to store the measured NS mass
  vector<double> mass_ns;
  vector<double> mass_wd;
  vector<double> mass_lms;
  vector<double> mass_hms;

  // std::vector to store +68% central limits of NS masses
  vector<double> uplim_ns;
  vector<double> uplim_wd;
  vector<double> lim_lms; // Symmetric 68% limits for NS-MS
  vector<double> lim_hms;

  // std::vector to store -68% central limits of NS masses
  vector<double> lowlim_ns;
  vector<double> lowlim_wd;

  // Total number of stars in all populations
  size_t n_stars;

  // Function to load population mass data
  void load_data();
};

class likelihood {
  
  public:

  bool debug;
  
  // Constructor to load source data
  likelihood () {
    md.load_data();
    debug=false;
  }

  virtual ~likelihood() {
  } 

  // Object to load source data from class mass_data
  mass_data md;
  
  /* Tolerance for small weights, below which weights are
  ignored. This ensures that log-weights do not explode. */
  const double tol=1.0e-6;
    
  // Counts the total number of population parameters 
  size_t n_params;

  // Counts the total number of distribution parameters
  size_t n_dist_pars;

  // Vector to store the names of all population parameters
  vector<string> par_names;

  // Vector to store the units of all population parameters
  vector<string> par_units;

  /* Skewed Normal PDF: SN(M_star, mean, width, skewness) 
  [eq. 13, Kiziltan et al. (2013)] */
  double skew_norm(double, double, double, double);

  /* Asymmetric Normal PDF: AN(mass-M_star, asym, scale)
  [eq. 14, Kiziltan et al. (2013)] */
  double asym_norm(double, double, double);

  // The function to solve f(x)=0 [see refs/calc.pdf]
    double f_to_solve(double, double &, double &);
    
  /* Derivative of the function to solve f'(x)
  (for use with root_stef only) */
  double df_to_solve(double, double &, double &);
    
  /* Solver to calculate scale parameters d for a given 
  asymmetry parameter c of function asym_norm (AN):
  get_scale(lowlim, highlim) */
  double get_scale(double, double);

  /* Function to fill vectors with names and units of all 
  population parameter */
  void get_params();
    
  /* Function to fill the pvi object with names of all 
  population parameter */
  void set_params(vec_index &);

  // Likelihood function for NS-NS
  double get_weight_ns(const ubvector &, vec_index &, int &);
    
  // Likelihood function for NS-WD
  double get_weight_wd(const ubvector &, vec_index &, int &);
    
  // Likelihood function for NS-MS/HMXB
  double get_weight_hms(const ubvector &, vec_index &, int &);

  // Likelihood function for NS-MS/LMXB
  double get_weight_lms(const ubvector &, vec_index &, int &);
    
  /* Combined likelihood function for all stars, except 
  GW170817, QLMXBs, PREs, and NICER */
  double get_weight(const ubvector &, vec_index &, int &);

  
  private:

  // PDF of standard normal distribution N(0,1)
  double norm_pdf(double);

  // CDF of standard N(0,1) in terms of erf(x)
  double norm_cdf(double);

};

#endif
