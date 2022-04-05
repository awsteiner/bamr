#include "likelihood.h"
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>

using namespace std;
using namespace o2scl;

int n_times=1;
double m=0.0, w=0.0, s=0.0, M=0.0;

const int n_pts=10;
const double dm = (2.5-0.5)/((double)n_pts);
const double dw = (1.0-0.0)/((double)n_pts);
const double ds = (1.0+1.0)/((double)n_pts);
const double dM = (2.3-1.0)/((double)n_pts);

// The likelihood function for NS-NS (see refs/method.pdf)
double get_weight_ns(size_t nv, const ubvector &x) {
    
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];

  likelihood lk;
  lk.get_params();

  mass_data md;
  md.load_data();

  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;
  
  for (size_t i=0; i<1; i++) {
    mass = md.mass_ns[i]; 
    uplim = md.uplim_ns[i];
    lowlim = md.lowlim_ns[i]; 
    asym = sqrt(uplim/lowlim);
    scale = lk.get_scale(lowlim, uplim);
    M_star = x[i+3];
    wgt_star = lk.asym_norm(mass-M_star, asym, scale) 
      * lk.skew_norm(M_star, mean, width, skewness);
    if (wgt_star==0.0) {
      cout << "Zero weight found in NS" << endl;
      wgt_star = 1.0; // Ignore small likelihoods
    }
    wgt *= wgt_star; 
  }
  return wgt;
}

// The likelihood function for NS-WD (see refs/method.pdf)
double get_weight_wd(size_t nv, const ubvector &x) {
  
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];
  double M_star = M; 

  likelihood lk;
  lk.get_params();

  mass_data md;
  md.load_data();

  double mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<md.id_wd.size(); i++) {
    mass = md.mass_wd[i]; 
    uplim = md.uplim_wd[i];
    lowlim = md.lowlim_wd[i];
    asym = sqrt(uplim/lowlim);
    scale = lk.get_scale(lowlim, uplim);
    wgt_star = lk.asym_norm(mass-M_star, asym, scale) 
      * lk.skew_norm(M_star, mean, width, skewness);
    /*if (wgt_star==0.0) {
      wgt_star = 1.0; // Ignore small likelihoods
    }*/
    wgt *= wgt_star; 
  }
  return wgt;
}

// The likelihood function for NS-MS (see refs/method.pdf)
double get_weight_ms(size_t nv, const ubvector &x) {
  
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];
  double M_star = M; 
  
  likelihood lk;
  lk.get_params();

  mass_data md;
  md.load_data();
  
  double mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<md.id_ms.size(); i++) {
    mass = md.mass_ms[i]; 
    uplim = md.lim_ms[i];
    lowlim = uplim; // Symmetric 68% limits
    asym = sqrt(uplim/lowlim); 
    scale = lk.get_scale(lowlim, uplim);
    wgt_star = lk.asym_norm(mass-M_star, asym, scale) 
      * lk.skew_norm(M_star, mean, width, skewness);
    /*if (wgt_star==0.0) {
      wgt_star = 1.0; // Ignore small likelihoods
    }*/
    wgt *= wgt_star; 
  }
  return wgt;
}

// The combined likelihood function to be calculated
double get_weight(size_t nv, const ubvector &x) {

  double wgt_ns, wgt_wd, wgt_ms, wgt; 

  // Calculate log-likelihood for each population
  wgt_ns = get_weight_ns(nv, x);
  wgt_wd = get_weight_wd(nv, x);
  wgt_ms = get_weight_ms(nv, x);

  // Multiply all likelihoods. Note: This is log-likelihood.
  wgt = wgt_ns * wgt_wd * wgt_ms;
  
  // Return the log-likelihood
  return wgt;
}

int main(void) {

  likelihood lk;
  lk.get_params();

  mass_data md;
  md.load_data();

  ofstream file; // file_ns, file_wd, file_ms;
  // file_ns.open("L_ns.dat");
  // file_wd.open("L_wd.dat");
  // file_ms.open("L_ms.dat");
  file.open("test.dat");
 
  double res=0.0,    err=0.0;
  double res_ns=0.0, err_ns=0.0;
  double res_wd=0.0, err_wd=0.0;
  double res_ms=0.0, err_ms=0.0;
  
  size_t n_dpars = 3;
  size_t n_mpars = md.id_ns.size();
  size_t n_pars = n_dpars + n_mpars;

  mcarlo_vegas<> gm;
  ubvector a(4), b(4);

  a[0]=0.5;  b[0]=2.5;
  a[1]=0.0;  b[1]=1.0;
  a[2]=-1.0; b[2]=1.0;
  a[3]=1.55;  b[3]=1.56;

  /*for (size_t i=0; i<n_mpars; i++) {
    a[i+n_dpars] = 1.4;
    b[i+n_dpars] = 1.6;
  }*/

  multi_funct f = get_weight;
  multi_funct f_ns = get_weight_ns;
  multi_funct f_wd = get_weight_wd;
  multi_funct f_ms = get_weight_ms;
 
  gm.n_points=100000;

  // cout << "n_params = " << lk.n_params << endl;
  // cout << "n_dist_pars = " << lk.n_dist_pars << endl;
  // cout << "id_ns.size() = " << md.id_ns.size() << endl;

  /* for (int i=0; i<n_pts; i++) {
    gm.minteg_err(f, 3, a, b, res, err);
    M += dM;
    file << M << "\t" << res << endl;
    cout << "i=" << i << "\t M = " << M << "\t L(M) = " << res << endl;
  } */
  
  gm.minteg_err(f_ns, 4, a, b, res_ns, err_ns);
  
  // gm.minteg_err(f_wd, 3, a, b, res_wd, err_wd);
  // gm.minteg_err(f_ms, 3, a, b, res_ms, err_ms);
  
  // cout << "res_all = " << res << "\t err_all = " << err << endl;
  
  cout << "res_ns = " << res_ns << "\t err_ns = " << err_ns << endl;
  
  // cout << "res_wd = " << res_wd << "\t err_wd = " << err_wd << endl;
  // cout << "res_ms = " << res_ms << "\t err_ms = " << err_ms << endl;
 
  // file_ns.close();
  // file_wd.close();
  // file_ms.close();
  
  file.close();
  return 0;
}