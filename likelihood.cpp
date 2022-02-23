#include "likelihood.h"

// PDF of standard normal distribution N(0,1)
double like::norm_pdf(double x) {
  return exp(-0.5*x*x) / sqrt(2.0*M_PI);
}

// CDF of standard N(0,1) in terms of erf(x)
double like::norm_cdf(double x) {
  return 0.5 * (1.0 + erf(x/sqrt(2.0)));
}

// Skewed Normal PDF [eq. 13, refs/kiziltan13]
double like::skew_norm(double x, double mean, double width, double asym) {
  return 2.0 * norm_pdf((x-mean)/width)
    * norm_cdf((x-mean)*asym/width) / width;
}

// Asymmetric Normal PDF [eq. 14, refs/kiziltan13]
double like::asym_norm(double x, double c, double d) {
  double a = 2.0 / (d*(c+1.0/c));
  if (x>=0.0)
    return a * norm_pdf(x/(c*d));
  else
    return a * norm_pdf(c*x/d);
}

// This is the function to solve [see refs/calc.pdf]
double like::f2solve(double x, double &l, double &u) {
  double c = sqrt(u/l);
  return c*c*erf(u/(sqrt(2.0)*c*x)) - erf(-c*l/(sqrt(2.0)*x))
    - 0.68*(c*c+1.0);
}

// Derivative of the function to solve (for use with root_stef)
double like::df2solve(double x, double &l, double &u) {
  double c = sqrt(u/l);
  double a = sqrt(2.0/M_PI) * c / x / x;
  return a*l*exp(-pow(u/(sqrt(2.0)*c*x), 2.0))
    - a*u*exp(-pow(c*l/sqrt(2.0)/x, 2.0));
}

// The solver that calculates parameters dj, given cj = sqrt(uj/lj)
double like::calc_d(double l, double u) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  // Only print something out if one of the tests fails
  t.set_output_level(1);
  
  // The solver, specifying the function type: funct<double>
  root_brent_gsl<> solver;
  
  like c;
  
  /* This is the code that allows specification of class member
     functions as functions to solve. This approach avoids the use of
     static variables and functions and multiple inheritance at the
     expense of a little overhead. We need to provide the address of
     an instantiated object and the address of the member function. */
  funct f2 = bind(mem_fn<double(double, double &, double &)>
		  (&like::f2solve), &c, _1, ref(l), ref(u));
  
  /* funct df2 = bind(mem_fn<double(double, double &, double &)>
     (&like::df2solve), &c, _1, ref(l), ref(u)); */
  
  // The root is bracketted in [x1, x2]
  double x1=0.0, x2=1.0;
  
  /* The value verbose=1 prints out iteration information
     and verbose=2 requires a keypress between iterations. */
  solver.verbose=0;
  
  solver.solve_bkt(x1, x2, f2);
  
  // Obtain and summarize test results
  // t.report();
  
  return x1;
}

// The likelihood function for NS-NS (see refs/method.pdf)
double like::calc_likelihood_ns(const ubvector &pars, vec_index &pvi) {

  double mean = pars[pvi["mean_ns"]];
  double width = pars[pvi["width_ns"]];
  double asym = pars[pvi["asym_ns"]];

  double mj, lj, uj, cj, dj, Lj, L=1.0;

  for (size_t j=0; j<N_ns; j++) {
    mj = nsns[j][0]; 
    uj = nsns[j][1];
    lj = nsns[j][2]; 
    cj = sqrt(uj/lj);
    dj = calc_d(lj, uj);
    double Mj = pars[pvi[((string)"M_")+nsns_sid[j]]];
    Lj = asym_norm(mj-Mj, cj, dj) * skew_norm(Mj, mean, width, asym);
    if (Lj<tol) Lj = 1.0; // Ignore small likelihoods
    L *= Lj; 
  }
  return L;
}

// The likelihood function for NS-WD (see refs/method.pdf)
double like::calc_likelihood_wd(const ubvector &pars, vec_index &pvi) {
  
  double mean = pars[pvi["mean_wd"]];
  double width = pars[pvi["width_wd"]];
  double asym = pars[pvi["asym_wd"]];
  
  double mj, lj, uj, cj, dj, Lj, L=1.0;
  
  for (size_t j=0; j<N_wd; j++) {
    mj = nswd[j][0];
    uj = nswd[j][1];
    lj = nswd[j][2];
    cj = sqrt(uj/lj);
    dj = calc_d(lj, uj);
    double Mj = pars[pvi[((string)"M_")+nswd_sid[j]]];
    Lj = asym_norm(mj-Mj, cj, dj) * skew_norm(Mj, mean, width, asym);
    if (Lj<tol) Lj = 1.0; // Ignore small likelihoods
    L *= Lj; 
  }
  return L;
}

// The likelihood function for NS-MS (see refs/method.pdf)
double like::calc_likelihood_ms(const ubvector &pars, vec_index &pvi) {
  
  double mean = pars[pvi["mean_ms"]];
  double width = pars[pvi["width_ms"]];
  double asym = pars[pvi["asym_ms"]];

  double mj, lj, uj, cj, dj, Lj, L=1.0;

  for (size_t j=0; j<N_ms; j++) {
    mj = nsms[j][0];
    uj = nsms[j][1];
    lj = uj; // Symmetric 68% limits
    cj = sqrt(uj/lj); 
    dj = calc_d(lj, uj);
    double Mj = pars[pvi[((string)"M_")+nsms_sid[j]]];
    Lj = asym_norm(mj-Mj, cj, dj) * skew_norm(Mj, mean, width, asym);
    if (Lj<tol) Lj = 1.0; // Ignore small likelihoods
    L *= Lj; 
  }
  return L;
}

/** \brief Set the vec_index object with the parameters from
    the mass data.
    
    This function will be called by bamr to fill the \c pvi
    object with the all parameters from the data set.
*/
void like::set_params(vec_index &pvi) {
  
// Start with NS-NS parameters
  pvi.append("mean_ns");
  pvi.append("width_ns");
  pvi.append("asym_ns");
  for(size_t i=0; i<N_ns; i++) {
    string mass_par=((string)"M_")+nsns_sid[i];
    pvi.append(mass_par);
  }
// Next, fill in NS-WD parameters
  pvi.append("mean_wd");
  pvi.append("width_wd");
  pvi.append("asym_wd");
  for(size_t i=0; i<N_wd; i++) {
    string mass_par=((string)"M_")+nswd_sid[i];
    pvi.append(mass_par);
  }
// Finally, fill in NS-MS parameters
  pvi.append("mean_ms");
  pvi.append("width_ms");
  pvi.append("asym_ms");
  for(size_t i=0; i<N_ms; i++) {
    string mass_par=((string)"M_")+nsms_sid[i];
    pvi.append(mass_par);
  }
  return;
}

/// The combined likelihood function to be calculated
double like::calc_likelihood(const ubvector &pars, vec_index &pvi) {
  /*
  default_random_engine seed;
  uniform_real_distribution<double> fmean(0.5, 2.5);
  uniform_real_distribution<double> fsigma(0.0, 1.0);
  uniform_real_distribution<double> falpha(-1.0, 1.0);
  uniform_real_distribution<double> fmass(1.0, 2.3);
  */

  double L_ns, L_wd, L_ms, L;
  
  this->load_data(); // Load source data 

  // Calculate likelihood for each population
  L_ns = calc_likelihood_ns(pars, pvi);
  L_wd = calc_likelihood_wd(pars, pvi);
  L_ms = calc_likelihood_ms(pars, pvi);

  // Multiply all likelihoods. Note: This is not log-likelihood.
  L = L_ns * L_wd * L_ms;
  
  return log(L);
}
