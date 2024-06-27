#ifndef NS_POP_H
#define NS_POP_H

#include <fstream>
#include <cmath>
#include <string>
#include <o2scl/funct.h>
#include <o2scl/constants.h>
#include <o2scl/root_brent_gsl.h>
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_math.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace std::placeholders;

typedef boost::numeric::ublas::vector<double> ubvector;

struct pop_data {
  
  /* std::string to store the full names of stars */
  vector<string> name_ns;
  vector<string> name_wd;
  vector<string> name_lx;
  // vector<string> name_hms;

  /* std::string to store the shortened names of stars */
  vector<string> id_ns;
  vector<string> id_wd;
  vector<string> id_lx;
  // vector<string> id_hms;

  /* std::vector to store the measured NS mass */
  vector<double> mass_ns;
  vector<double> mass_wd;
  vector<double> mass_lx;
  // vector<double> mass_hms;

  /* std::vector to store 68% central limits of NS masses */
  vector<double> low_ns;
  vector<double> low_wd;
  vector<double> high_ns;
  vector<double> high_wd;
  vector<double> low_lx;
  vector<double> high_lx;
  // vector<double> low_hx;
  // vector<double> high_hx;

  /* Total number of stars in all populations */
  size_t n_stars;

  /* Function to load population mass data */
  void load_data();

};

class ns_pop {
  
  public:

    bool debug;
  
    /* Constructor to load source data */
    ns_pop() {
      pd.load_data();
      debug=false;
    }

    virtual ~ns_pop() {
    } 

    /* Object to load source data from class pop_data */
    pop_data pd;

    /* Counts the total number of population parameters */
    size_t n_pop_params;

    /* std::vector<string> to store the names and units of 
    population parameters */
    vector<string> par_names;
    vector<string> par_units;

    /* std::vector<double> to store lower and upper limits of 
    priors for population parameters */
    vector<double> par_low;
    vector<double> par_high;
    
    /* std::vector<double> to store the initial points for 
    population parameters */
    vector<double> par_init;

    /* Skewed Normal PDF: SN(M_star, mean, width, skewness) 
    [eq. 13, Kiziltan et al. (2013)] */
    double skewed_norm(double, double, double, double);

    /* Asymmetric Normal PDF: AN(mass-M_star, asym, scale)
    [eq. 14, Kiziltan et al. (2013)] */
    double asym_norm(double, double, double);

    /* Derivative of the Skewed Normal PDF */
    double deriv_sn(int , double , double , double , double );

    /* Derivative of the Asymmetric Normal PDF */
    double deriv_an(double , double , double );

    /* Function to fill vectors with the names and units of 
    population parameter */
    void get_param_info();

    /* Function to fill the pvi object with names of all 
    population parameter */
    void set_params(vec_index &);

    /* Likelihood function for NS-NS */
    double get_weight_ns(const ubvector &, vec_index &, int &);

    /* Likelihood function for NS-WD */
    double get_weight_wd(const ubvector &, vec_index &, int &);

    /* Likelihood function for NS-MS/HMXB */
    double get_weight_hx(const ubvector &, vec_index &, int &);

    /* Likelihood function for NS-MS/LMXB */
    double get_weight_lx(const ubvector &, vec_index &, int &);

    /* Combined likelihood function for all stars, except 
    GW170817, QLMXBs, PREs, and NICER */
    double get_weight(const ubvector &, vec_index &, int &);

    /* Store the values of AN(x|c,d) for all stars in a binary */
    vector<double> an_ns;
    vector<double> an_wd;
    vector<double> an_lx;

    /* Store the values of SN(x|m,s,a) for all stars in a binary*/
    vector<double> sn_ns;
    vector<double> sn_wd;
    vector<double> sn_lx;


  protected:

    /* PDF of standard normal distribution N(0,1) */
    double norm_pdf(double);

    /* CDF of standard N(0,1) in terms of erf(x) */
    double norm_cdf(double);

};

class eqn_solver {

  public:

    /* The function to solve f(x)=0 [see refs/calc.pdf] */
    double f_to_solve(double, double &, double &);
    double f2_to_solve(double, double &, double &);

    /* Derivative of the function to solve f'(x)
    (for use with root_stef only) */
    double df_to_solve(double, double &, double &);

    /* Solver to calculate scale parameter d for a given 
    asymmetry parameter c of function asym_norm (AN):
    get_scale(lowlim, highlim) */
    double get_scale(double, double);

    /* Solver to calculate mass m2 given chirp mass M_chirp 
    and mass m1 */
    double get_m2(double, double);

};

#endif
