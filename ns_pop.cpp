#include "ns_pop.h"


// PDF of standard normal distribution N(0,1)
double ns_pop::norm_pdf(double x) {
  return exp(-0.5*x*x)/sqrt(2.0*M_PI);
}

// CDF of standard normal N(0,1) in terms of erf(x)
double ns_pop::norm_cdf(double x) {
  return 0.5*(1.0+erf(x/sqrt(2.0)));
}

// Skewed Normal PDF 
double ns_pop::skewed_norm(double x, double m, double s, double a) {
  return 2.0*norm_pdf((x-m)/s)*norm_cdf((x-m)*a/s)/s;
}

// Asymmetric Normal PDF 
double ns_pop::asym_norm(double x, double c, double d) {
  double k = 2.0/(d*(c+1.0/c));
  if (x<0.0) return k*norm_pdf(c*x/d);
  else return k*norm_pdf(x/c/d);
}

double ns_pop::deriv_sn(const int i_wrt, double x, double m, 
                        double s, double a) {
  double k=1.0/s/sqrt(2.0*M_PI);
  double ef=exp(-0.5*pow((x-m)/s,2.0));
  double p=(x-m)*a/s/sqrt(2.0);
  double t1=a*sqrt(2.0/M_PI)*exp(-p*p)/s;
  double t2=(x-m)*(1.0+erf(p))/s/s;

  if (i_wrt==0) return k*ef*(t1-t2);
  if (i_wrt==1) return k*ef*(t2-t1);
  if (i_wrt==2) return k*ef*(t2-t1)*(x-m)/s;
  if (i_wrt==3) return k*ef*exp(-p*p);
}

double ns_pop::deriv_an(double x, double c, double d) {
  double k=2.0/(d*(c+1.0/c));
  if (x<0.0) return -k*x*norm_pdf(c*x/d)*c*c/d/d;
  else return -k*x*norm_pdf(x/(c*d))/c/c/d/d;
}


// The likelihood function for NS-NS (see refs/method.pdf)
double ns_pop::get_weight_ns(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  ret=0;
  double mean=pars[pvi["mean_NS"]];
  double log10_width=pars[pvi["log10_width_NS"]];
  double width=pow(10.0, log10_width);
  double skew=pars[pvi["skewness_NS"]];
  double log_wgt=0.0;
  eqn_solver es;
  
  if (this->debug) {
    cout << "index name mass(data) asym scale M_star(param) "
         << "mean width skewness wgt AN SN" << endl;
  }
  for (size_t i=0; i<pd.id_ns.size(); i++) {
    double mass=pd.mass_ns[i];
    double high=pd.high_ns[i];
    double low=pd.low_ns[i]; 
    double asym=sqrt(high/low);
    double scale=es.get_scale(low, high);
    double M_star=pars[pvi[string("M_")+pd.id_ns[i]]];
    an_ns[i]=asym_norm(mass-M_star, asym, scale);
    sn_ns[i]=skewed_norm(M_star, mean, width, skew);
    double wgt_star=an_ns[i]*sn_ns[i];
    
    if (this->debug) {
      cout << "NS: " << i << " " << pd.id_ns[i] << " "
           << mass << " " << asym << " " << scale << " " 
           << M_star << " " << mean << " " << width << " " 
           << skew << " " << wgt_star;
      cout << " " << an_ns[i]  << " " << sn_ns[i] << endl;
    }
    if (wgt_star<=0.0) {
      /* Record index i (via ret) for bookkeeping (scr_out): 
      30 is added to avoid ret=0 when wgt_star=0 */
      ret=30+i; 
      return 0.0;
    }
    log_wgt+=log(wgt_star); 
  }
  if (this->debug) cout << "NS: log_wgt = " << log_wgt << endl;
  return log_wgt;
}


// The likelihood function for NS-WD (see refs/method.pdf)
double ns_pop::get_weight_wd(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  ret=0;
  double mean=pars[pvi["mean_WD"]];
  double log10_width=pars[pvi["log10_width_WD"]];
  double width=pow(10.0, log10_width);
  double skew=pars[pvi["skewness_WD"]];
  double log_wgt=0.0;
  eqn_solver es;

  if (this->debug) {
    cout << "index name mass(data) asym scale M_star(param) "
         << "mean width skewness wgt AN SN" << endl;
  }
  for (size_t i=0; i<pd.id_wd.size(); i++) {
    double mass=pd.mass_wd[i];
    double high=pd.high_wd[i];
    double low=pd.low_wd[i];
    double asym=sqrt(high/low);
    double scale=es.get_scale(low, high);
    double M_star=pars[pvi[string("M_")+pd.id_wd[i]]];
    an_wd[i]=asym_norm(mass-M_star, asym, scale);
    sn_wd[i]=skewed_norm(M_star, mean, width, skew);
    double wgt_star=an_wd[i]*sn_wd[i];
    
    if (this->debug) {
      cout << "WD: " << i << " " << pd.id_wd[i] << " "
           << mass << " " << asym << " " << scale << " " 
           << M_star << " " << mean << " " << width << " " 
           << skew << " " << wgt_star; 
      cout << " " << an_wd[i]  << " " << sn_wd[i] << endl;
    }
    if (wgt_star<=0.0) {
      ret=60+i;
      return 0.0;
    }
    log_wgt+=log(wgt_star); 
  }
  if (this->debug) cout << "WD: " << log_wgt << endl;
  return log_wgt;
}


// The likelihood function for NS-MS/LMXB 
double ns_pop::get_weight_lx(const ubvector &pars, vec_index &pvi,
                                 int &ret) {
  ret=0;
  double mean=pars[pvi["mean_LMS"]];
  double log10_width=pars[pvi["log10_width_LMS"]];
  double width=pow(10.0, log10_width);
  double skew=pars[pvi["skewness_LMS"]];
  double log_wgt=0.0;
  eqn_solver es;

  if (this->debug) {
    cout << "index name mass(data) asym scale M_star(param) "
         << "mean width skewness wgt AN SN" << endl;
  }
  for (size_t i=0; i<pd.id_lx.size(); i++) {
    double mass=pd.mass_lx[i]; 
    double low=pd.low_lx[i];
    double high=pd.high_lx[i];
    double asym=sqrt(high/low);
    double scale=es.get_scale(low, high);
    double M_star=pars[pvi[string("M_")+pd.id_lx[i]]];
    an_lx[i]=asym_norm(mass-M_star, asym, scale);
    sn_lx[i]=skewed_norm(M_star, mean, width, skew);
    double wgt_star=an_lx[i]*sn_lx[i];
    
    if (this->debug) {
      cout << "LMXB: " << i << " " << pd.id_lx[i] << " "
           << mass << " " << asym << " " << scale << " " 
           << M_star << " " << mean << " " << width << " " 
           << skew << " " << wgt_star; 
      cout << " " << an_lx[i]  << " " << sn_lx[i] << endl;
    }
    if (wgt_star<=0.0) {
      ret=100+i;
      return 0.0;
    }
    log_wgt+=log(wgt_star); 
  }
  if (this->debug) {
    cout << "LMXB: " << log_wgt << endl;
    exit(-1);
  }
  return log_wgt;
}


// The combined likelihood function to be calculated
double ns_pop::get_weight(const ubvector &pars, vec_index &pvi,
                              int &ret) {

  double wgt_ns, wgt_wd, wgt_hms, wgt_lx, wgt; 

  // Calculate log-likelihood for each population
  wgt_ns=get_weight_ns(pars, pvi, ret);
  wgt_wd=get_weight_wd(pars, pvi, ret);
  // wgt_hms=get_weight_hx(pars, pvi, ret);
  wgt_lx=get_weight_lx(pars, pvi, ret);

  // Multiply all likelihoods. Note: This is log-likelihood.
  wgt=wgt_ns+wgt_wd+wgt_lx; // + wgt_hms
  
  // Return the log-likelihood
  return wgt;
}


void ns_pop::get_param_info() {

  /* Note: pd.load_data() must be called before calling this function */

  // double 🐖=1.0;
  
  // Fill names, units, and initial points of distribution parameters
  par_names.push_back("mean_NS");
  par_units.push_back("Msun");
  par_init.push_back(1.3);
  
  par_names.push_back("log10_width_NS");
  par_units.push_back("");
  par_init.push_back(-0.7);
  
  par_names.push_back("skewness_NS");
  par_units.push_back("");
  par_init.push_back(0.0);
  
  par_names.push_back("mean_WD");
  par_units.push_back("Msun");
  par_init.push_back(1.7);
  
  par_names.push_back("log10_width_WD");
  par_units.push_back("");
  par_init.push_back(-0.5);
  
  par_names.push_back("skewness_WD");
  par_units.push_back("");
  par_init.push_back(0.0);
  
  /* par_names.push_back("mean_HMS");
  par_units.push_back("Msun");
  par_init.push_back(1.5);
  
  par_names.push_back("log10_width_HMS");
  par_units.push_back("");
  par_init.push_back(-0.5);
  
  par_names.push_back("skewness_HMS");
  par_units.push_back("");
  par_init.push_back(0.0); */
  
  par_names.push_back("mean_LMS");
  par_units.push_back("Msun");
  par_init.push_back(1.5);
  
  par_names.push_back("log10_width_LMS");
  par_units.push_back("");
  par_init.push_back(-0.5);
  
  par_names.push_back("skewness_LMS");
  par_units.push_back("");
  par_init.push_back(0.0);

  // Fill names, units, and initial points for mass parameters
  for(size_t i=0; i<pd.id_ns.size(); i++) {
    par_names.push_back(string("M_")+pd.id_ns[i]);
    par_units.push_back("Msun");
    par_init.push_back(pd.mass_ns[i]);
  }
  for(size_t i=0; i<pd.id_wd.size(); i++) {
    par_names.push_back(string("M_")+pd.id_wd[i]);
    par_units.push_back("Msun");
    par_init.push_back(pd.mass_wd[i]);
  }
  /* for(size_t i=0; i<pd.id_hms.size(); i++) {
    par_names.push_back(string("M_")+pd.id_hms[i]);
    par_units.push_back("Msun");
    par_init.push_back(pd.mass_hms[i]);
  } */
  for(size_t i=0; i<pd.id_lx.size(); i++) {
    par_names.push_back(string("M_")+pd.id_lx[i]);
    par_units.push_back("Msun");
    par_init.push_back(pd.mass_lx[i]);
  }

  // Set priors for distribution parameters
  for(size_t i=0;i<3;i++) {
    par_low.push_back(0.5);
    par_high.push_back(2.5);
    par_low.push_back(-6.0);
    par_high.push_back(0.0);
    par_low.push_back(-1.0);
    par_high.push_back(1.0);
  }

  // Set priors for mass parameters
  for (size_t i=0; i<pd.n_stars; i++) {
    par_low.push_back(1.0);
    par_high.push_back(2.5);
  }
  
  n_pop_params=par_names.size();

  return;
}


/** \brief Set the vec_index object with the parameters from
    the mass data.
    
    This function will be called by bamr to fill the \c pvi
    object with the all parameters from the data set.
*/
void ns_pop::set_params(vec_index &pvi) {

  // Fill in DNS parameters
  pvi.append("mean_NS");
  pvi.append("log10_width_NS");
  pvi.append("skewness_NS");
  for(size_t i=0; i<pd.id_ns.size(); i++) {
    string mass_par=string("M_")+pd.id_ns[i];
    pvi.append(mass_par);
  }
  // Fill in NS-WD parameters
  pvi.append("mean_WD");
  pvi.append("log10_width_WD");
  pvi.append("skewness_WD");
  for(size_t i=0; i<pd.id_wd.size(); i++) {
    string mass_par=string("M_")+pd.id_wd[i];
    pvi.append(mass_par);
  }
  // Fill in LMXB parameters
  pvi.append("mean_LMS");
  pvi.append("log10_width_LMS");
  pvi.append("skewness_LMS");
  for(size_t i=0; i<pd.id_lx.size(); i++) {
    string mass_par=string("M_")+pd.id_lx[i];
    pvi.append(mass_par);
  }
  return;
}


// This is the function to solve for d, given c
double eqn_solver::f_to_solve(double x, double &l, double &u) {
  double c=sqrt(u/l);
  return c*c*erf(u/(sqrt(2.0)*c*x)) - erf(-c*l/(sqrt(2.0)*x))
    - 0.68*(c*c+1.0);
}

// This is the function to solve for m2, given q
double eqn_solver::f2_to_solve(double x, double &M_chirp, double &m1) {
  double k=pow(m1/M_chirp, 5.0);
  return k*pow(x,3.0)-x-1.0;
}


// Derivative of the function to solve 
double eqn_solver::df_to_solve(double x, double &l, double &u) {
  double c=sqrt(u/l);
  double a=sqrt(2.0/M_PI) * c / x / x;
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
  funct f=bind(mem_fn<double(double, double &, double &)>
		  (&eqn_solver::f_to_solve), &es, _1, ref(l), ref(u));
  
  /* funct df2 = bind(mem_fn<double(double, double &, double &)>
     (&eqn_solver::df_to_solve), &c, _1, ref(l), ref(u)); */
  
  // The root is bracketted in [x1, x2]
  double x1=0.0, x2=1.0;
  
  /* The value verbose=1 prints out iteration information
     and verbose=2 requires a keypress between iterations. */
  solver.verbose=0;
  
  solver.solve_bkt(x1, x2, f); 
  // cout << "f(x) = " << f(x1) << endl;
  
  return x1;
}

// Solver to calculates mass2, given M_chirp and m1
double eqn_solver::get_m2(double M_chirp, double m1) {
  
  cout.setf(ios::scientific);
  
  // The solver, specifying the function type: funct<double>
  root_brent_gsl<> solver;

  eqn_solver es;
  funct f=bind(mem_fn<double(double, double &, double &)>
		  (&eqn_solver::f2_to_solve), &es, _1, ref(M_chirp), ref(m1));
  
  // The root is bracketted in [x1, x2]
  double x1=0.0, x2=3.0;
  
  solver.solve_bkt(x1, x2, f); 
  
  // Note: x1=q=m2/m1
  return m1*x1;
}