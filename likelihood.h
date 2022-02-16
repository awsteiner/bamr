#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_math.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace std::placeholders;

typedef boost::numeric::ublas::vector<double> ubvector;

// std::string vectors to store the full names/identifiers of the stars
string nsns_id[26],  nswd_id[39],  nsms_id[16];

// std::string vectors to store the shortened names/identifiers
string nsns_sid[26], nswd_sid[39], nsms_sid[16];

/* Arrays to store the raw data in the format: [mi][li, ui] for NS-NS
and NS-WD, and [mi][si] for NS-MS (where si=li=ui) */
double nsns[26][3], nswd[39][3], nsms[16][2];

class like {
private:
  double norm_pdf(double);
  double norm_cdf(double);
  double skew_norm(double, double, double, double);
  double asym_norm(double, double, double);
  double f2solve(double, double &, double &);
  double df2solve(double, double &, double &);
  double calc_d(double, double);
  double wgt_ns(double, double, double, ubvector);
  double wgt_wd(double, double, double, ubvector);
  double wgt_ms(double, double, double, ubvector);
public:
  double wgt_all(const ubvector &, vec_index &);
  void set_params(vec_index &);
  void load_data();
}
