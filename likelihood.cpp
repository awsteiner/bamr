#include "likelihood.h"


// PDF of standard normal distribution N(0,1)
double likelihood::norm_pdf(double x) {
  return exp(-0.5*x*x) / sqrt(2.0*M_PI);
}


// CDF of standard normal N(0,1) in terms of erf(x)
double likelihood::norm_cdf(double x) {
  return 0.5 * (1.0 + erf(x/sqrt(2.0)));
}


// Skewed Normal PDF 
double likelihood::skew_norm(double x, double mean, double width, 
    double skewness) {
  return 2.0 * norm_pdf((x-mean)/width)
    * norm_cdf((x-mean)*skewness/width) / width;
}


// Asymmetric Normal PDF 
double likelihood::asym_norm(double x, double c, double d) {
  double a = 2.0 / (d*(c+1.0/c));
  if (x>=0.0) return a * norm_pdf(x/(c*d));
  else return a * norm_pdf(c*x/d);
}


// This is the function to solve 
double eqn_solver::f_to_solve(double x, double &l, double &u) {
  double c = sqrt(u/l);
  return c*c*erf(u/(sqrt(2.0)*c*x)) - erf(-c*l/(sqrt(2.0)*x))
    - 0.68*(c*c+1.0);
}


// Derivative of the function to solve 
double eqn_solver::df_to_solve(double x, double &l, double &u) {
  double c = sqrt(u/l);
  double a = sqrt(2.0/M_PI) * c / x / x;
  return a*l*exp(-pow(u/(sqrt(2.0)*c*x), 2.0))
    - a*u*exp(-pow(c*l/sqrt(2.0)/x, 2.0));
}


// Solver to calculates parameters d, given c
double eqn_solver::get_scale(double l, double u) {
  
  cout.setf(ios::scientific);
  
  // The solver, specifying the function type: funct<double>
  root_brent_gsl<> solver;
  
  /* This is the code that allows specification of class member
     functions as functions to solve. We need to provide the address of
     an instantiated object and the address of the member function. */
  eqn_solver es;
  funct f = bind(mem_fn<double(double, double &, double &)>
		  (&eqn_solver::f_to_solve), &es, _1, ref(l), ref(u));
  
  /* funct df2 = bind(mem_fn<double(double, double &, double &)>
     (&likelihood::df_to_solve), &c, _1, ref(l), ref(u)); */
  
  // The root is bracketted in [x1, x2]
  double x1=0.0, x2=1.0;
  
  /* The value verbose=1 prints out iteration information
     and verbose=2 requires a keypress between iterations. */
  solver.verbose=0;
  
  solver.solve_bkt(x1, x2, f); 
  // cout << "f(x) = " << f(x1) << endl;
  
  return x1;
}


// The likelihood function for NS-NS (see refs/method.pdf)
double likelihood::get_weight_ns(const ubvector &pars, vec_index &pvi,
                                 int &ret) {

  double mean = pars[pvi["mean_NS"]];
  double log10_width = pars[pvi["log10_width_NS"]];
  double width = pow(10.0, log10_width);
  double skewness = pars[pvi["skewness_NS"]];
  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, log_wgt=0.0;
  eqn_solver es;
  
  if (debug) {
    cout << "index name mass(data) asym scale M_star(param) "
         << "mean width skewness wgt AN SN" << endl;
  }
  for (size_t i=0; i<md.id_ns.size(); i++) {
    mass = md.mass_ns[i]; 
    uplim = md.uplim_ns[i];
    lowlim = md.lowlim_ns[i]; 
    asym = sqrt(uplim/lowlim);
    scale = es.get_scale(lowlim, uplim);
    M_star = pars[pvi[string("M_")+md.id_ns[i]]];
    wgt_star = asym_norm(mass-M_star, asym, scale) 
      * skew_norm(M_star, mean, width, skewness);
    
    if (debug) {
      cout << "NS: " << i << " " << md.id_ns[i] << " "
           << mass << " " << asym << " " << scale << " " << M_star << " "
           << mean << " " << width << " " << skewness << " "
           << wgt_star; 
      cout << " " << asym_norm(mass-M_star, asym, scale)  << " "
           << skew_norm(M_star, mean, width, skewness) << endl;
    }
    if (wgt_star==0.0) {
      ret = 1;
      break;
    }
    log_wgt += log(wgt_star); 
  }
  if (debug) cout << "NS: log_wgt = " << log_wgt << endl;
  if (ret==0) return log_wgt;
  else return 0.0;
}


// The likelihood function for NS-WD (see refs/method.pdf)
double likelihood::get_weight_wd(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  
  double mean = pars[pvi["mean_WD"]];
  double log10_width = pars[pvi["log10_width_WD"]];
  double width = pow(10.0, log10_width);
  double skewness = pars[pvi["skewness_WD"]];
  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, log_wgt=0.0;
  eqn_solver es;

  if (debug) {
    cout << "index name mass(data) asym scale M_star(param) "
         << "mean width skewness wgt AN SN" << endl;
  }
  for (size_t i=0; i<md.id_wd.size(); i++) {
    mass = md.mass_wd[i]; 
    uplim = md.uplim_wd[i];
    lowlim = md.lowlim_wd[i];
    asym = sqrt(uplim/lowlim);
    scale = es.get_scale(lowlim, uplim);
    M_star = pars[pvi[string("M_")+md.id_wd[i]]];
    wgt_star = asym_norm(mass-M_star, asym, scale) 
      * skew_norm(M_star, mean, width, skewness);
    
    if (debug) {
      cout << "WD: " << i << " " << md.id_wd[i] << " "
           << mass << " " << asym << " " << scale << " " << M_star << " "
           << mean << " " << width << " " << skewness << " "
           << wgt_star; 
      cout << " " << asym_norm(mass-M_star, asym, scale)  << " "
           << skew_norm(M_star, mean, width, skewness) << endl;
    }
    if (wgt_star==0.0) {
      ret = 1;
      break;
    }
    log_wgt += log(wgt_star); 
  }
  if (debug) {
    cout << "WD: " << log_wgt << endl;
    exit(-1);
  }
  if (ret==0) return log_wgt;
  else return 0.0;
}


// The likelihood function for NS-MS/HMXB
double likelihood::get_weight_hms(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  
  double mean = pars[pvi["mean_HMS"]];
  double log10_width = pars[pvi["log10_width_HMS"]];
  double width = pow(10.0, log10_width);
  double skewness = pars[pvi["skewness_HMS"]];
  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, log_wgt=0.0;
  eqn_solver es;

  for (size_t i=0; i<md.id_hms.size(); i++) {
    mass = md.mass_hms[i]; 
    uplim = md.lim_hms[i];
    lowlim = uplim; // Symmetric 68% limits
    asym = sqrt(uplim/lowlim); 
    scale = es.get_scale(lowlim, uplim);
    M_star = pars[pvi[string("M_")+md.id_hms[i]]];
    wgt_star = asym_norm(mass-M_star, asym, scale) 
      * skew_norm(M_star, mean, width, skewness);
    if (wgt_star==0.0) {
      ret = 1;
      break;
    }
    log_wgt += log(wgt_star); 
  }
  if (debug) {
    cout << "MS/HMXB: " << log_wgt << endl;
    exit(-1);
  }
  if (ret==0) return log_wgt;
  else return 0.0;
}


// The likelihood function for NS-MS/LMXB 
double likelihood::get_weight_lms(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  
  double mean = pars[pvi["mean_LMS"]];
  double log10_width = pars[pvi["log10_width_LMS"]];
  double width = pow(10.0, log10_width);
  double skewness = pars[pvi["skewness_LMS"]];
  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, log_wgt=0.0;
  eqn_solver es;

  for (size_t i=0; i<md.id_lms.size(); i++) {
    mass = md.mass_lms[i]; 
    uplim = md.lim_lms[i];
    lowlim = uplim; // Symmetric 68% limits
    asym = sqrt(uplim/lowlim); 
    scale = es.get_scale(lowlim, uplim);
    M_star = pars[pvi[string("M_")+md.id_lms[i]]];
    wgt_star = asym_norm(mass-M_star, asym, scale) 
      * skew_norm(M_star, mean, width, skewness);
    if (wgt_star==0.0) {
      ret = 1;
      break;
    }
    log_wgt += log(wgt_star); 
  }
  if (debug) {
    cout << "MS/LMXB: " << log_wgt << endl;
    exit(-1);
  }
  if (ret==0) return log_wgt;
  else return 0.0;
}


// The combined likelihood function to be calculated
double likelihood::get_weight(const ubvector &pars, vec_index &pvi,
                              int &ret) {

  double wgt_ns, wgt_wd, wgt_hms, wgt_lms, wgt; 

  // Calculate log-likelihood for each population
  wgt_ns = get_weight_ns(pars, pvi, ret);
  wgt_wd = get_weight_wd(pars, pvi, ret);
  wgt_hms = get_weight_hms(pars, pvi, ret);
  wgt_lms = get_weight_lms(pars, pvi, ret);

  // Multiply all likelihoods. Note: This is log-likelihood.
  wgt = wgt_ns + wgt_wd + wgt_hms + wgt_lms;
  
  // Return the log-likelihood
  return wgt;
}


void likelihood::get_params() {

  /* Note: md.load_data() must be called before calling this function */

  // double 🐖=1.0;
  
  // Fill names and units of distribution parameters
  par_names.push_back("mean_NS");
  par_units.push_back("Msun");
  par_names.push_back("log10_width_NS");
  par_units.push_back("");
  par_names.push_back("skewness_NS");
  par_units.push_back("");
  par_names.push_back("mean_WD");
  par_units.push_back("Msun");
  par_names.push_back("log10_width_WD");
  par_units.push_back("");
  par_names.push_back("skewness_WD");
  par_units.push_back("");
  par_names.push_back("mean_HMS");
  par_units.push_back("Msun");
  par_names.push_back("log10_width_HMS");
  par_units.push_back("");
  par_names.push_back("skewness_HMS");
  par_units.push_back("");
  par_names.push_back("mean_LMS");
  par_units.push_back("Msun");
  par_names.push_back("log10_width_LMS");
  par_units.push_back("");
  par_names.push_back("skewness_LMS");
  par_units.push_back("");

  n_dist_pars = par_names.size();

  // Fill names and units of mass parameters
  for(size_t i=0; i<md.id_ns.size(); i++) {
    par_names.push_back(string("M_")+md.id_ns[i]);
    par_units.push_back("Msun");
  }
  for(size_t i=0; i<md.id_wd.size(); i++) {
    par_names.push_back(string("M_")+md.id_wd[i]);
    par_units.push_back("Msun");
  }
  for(size_t i=0; i<md.id_hms.size(); i++) {
    par_names.push_back(string("M_")+md.id_hms[i]);
    par_units.push_back("Msun");
  }
  for(size_t i=0; i<md.id_lms.size(); i++) {
    par_names.push_back(string("M_")+md.id_lms[i]);
    par_units.push_back("Msun");
  }

  n_params = par_names.size();

  return;
}


/** \brief Set the vec_index object with the parameters from
    the mass data.
    
    This function will be called by bamr to fill the \c pvi
    object with the all parameters from the data set.
*/
void likelihood::set_params(vec_index &pvi) {

  // Fill in NS-NS parameters
  pvi.append("mean_NS");
  pvi.append("width_NS");
  pvi.append("asym_NS");
  for(size_t i=0; i<md.id_ns.size(); i++) {
    string mass_par=string("M_")+md.id_ns[i];
    pvi.append(mass_par);
  }
  // Fill in NS-WD parameters
  pvi.append("mean_WD");
  pvi.append("width_WD");
  pvi.append("asym_WD");
  for(size_t i=0; i<md.id_wd.size(); i++) {
    string mass_par=string("M_")+md.id_wd[i];
    pvi.append(mass_par);
  }
  // Fill in NS-MS (HMXBs & LMXBs) parameters
  pvi.append("mean_HMS");
  pvi.append("width_HMS");
  pvi.append("asym_HMS");
  for(size_t i=0; i<md.id_hms.size(); i++) {
    string mass_par=string("M_")+md.id_hms[i];
    pvi.append(mass_par);
  }
  pvi.append("mean_LMS");
  pvi.append("width_LMS");
  pvi.append("asym_LMS");
  for(size_t i=0; i<md.id_lms.size(); i++) {
    string mass_par=string("M_")+md.id_lms[i];
    pvi.append(mass_par);
  }
  return;
}
