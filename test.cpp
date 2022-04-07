#include "test.h"

// The likelihood function for NS-NS (see refs/method.pdf)
double t_class::get_weight_ns(size_t n_vars, const ubvector &x) {
    
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];

  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;
  
  for (size_t i=0; i<mdat.id_ns.size(); i++) {
    mass = mdat.mass_ns[i]; 
    uplim = mdat.uplim_ns[i];
    lowlim = mdat.lowlim_ns[i]; 
    asym = sqrt(uplim/lowlim);
    scale = like.get_scale(lowlim, uplim);
    M_star = x[i+3];
    wgt_star = like.asym_norm(mass-M_star, asym, scale) 
      * like.skew_norm(M_star, mean, width, skewness);
    /*if (wgt_star==0.0) {
      cout << "Zero weight found in NS" << endl;
      wgt_star = 1.0; // Ignore small likelihoods
    }*/
    wgt *= wgt_star; 
    // cout << "mass=" << mass << "\t M_star=" << M_star << endl;
    // cout << "asym=" << asym << "\t scale=" << scale << endl;
    // cout << "AN=" << like.asym_norm(mass-M_star, asym, scale) << endl;
    // cout << "SN=" << like.skew_norm(M_star, mean, width, skewness) << endl;
  }
  return wgt;
}

// The likelihood function for NS-WD (see refs/method.pdf)
double t_class::get_weight_wd(size_t n_vars, const ubvector &x) {
  
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];

  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<mdat.id_wd.size(); i++) {
    mass = mdat.mass_wd[i]; 
    uplim = mdat.uplim_wd[i];
    lowlim = mdat.lowlim_wd[i];
    asym = sqrt(uplim/lowlim);
    M_star = x[i+3];
    scale = like.get_scale(lowlim, uplim);
    wgt_star = like.asym_norm(mass-M_star, asym, scale) 
      * like.skew_norm(M_star, mean, width, skewness);
    /*if (wgt_star==0.0) {
      wgt_star = 1.0; // Ignore small likelihoods
    }*/
    wgt *= wgt_star; 
  }
  return wgt;
}

// The likelihood function for NS-MS (see refs/method.pdf)
double t_class::get_weight_ms(size_t n_vars, const ubvector &x) {
  
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];
  
  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<mdat.id_ms.size(); i++) {
    mass = mdat.mass_ms[i]; 
    uplim = mdat.lim_ms[i];
    lowlim = uplim; // Symmetric 68% limits
    asym = sqrt(uplim/lowlim); 
    scale = like.get_scale(lowlim, uplim);
    M_star = x[i+3];
    wgt_star = like.asym_norm(mass-M_star, asym, scale) 
      * like.skew_norm(M_star, mean, width, skewness);
    /*if (wgt_star==0.0) {
      wgt_star = 1.0; // Ignore small likelihoods
    }*/
    wgt *= wgt_star; 
  }
  return wgt;
}

// The combined likelihood function to be calculated
double t_class::get_weight(size_t n_vars, const ubvector &x) {

  double wgt_ns, wgt_wd, wgt_ms, wgt; 

  // Calculate log-likelihood for each population
  wgt_ns = get_weight_ns(n_vars, x);
  wgt_wd = get_weight_wd(n_vars, x);
  wgt_ms = get_weight_ms(n_vars, x);

  // Multiply all likelihoods. Note: This is log-likelihood.
  wgt = wgt_ns * wgt_wd * wgt_ms;
  
  // Return the log-likelihood
  return wgt;
}


void t_class::set_dist_limits() {

  double y, mean, width, skewness, M_star;
  double m_min=1.2, m_max=1.8;
  double w_min=0.4, w_max=1.0;
  double s_min=-1.0, s_max=1.0;
  double M_min=1.0, M_max=2.4;
  double dx = 0.01;

  for (M_star=M_min; M_star<=M_max; M_star+=dx) {
    for (mean=m_min; mean<=m_max; mean+=dx) {
      for (width=w_min; width<=w_max; width+=dx) {
        for (skewness=s_min; skewness<=s_max; skewness+=dx) {
          y = like.skew_norm(M_star, mean, width, skewness);
          cout << fixed << setprecision(4);
          /*cout << "mu=" << mu << "\t sigma=" << sigma << "\t alpha=" << alpha
            << "\t M=" << M << endl; */
          if (y<tol) {
            cout << "y<tol at mu=" << mean << ",\t sigma=" << width 
              << ",\t alpha=" << skewness << "\t (M=" << M_star << ")" << endl;
          }
        }
      }
    }
  }
}

/*
void t_class::test_dist_limits() {

} */


void t_class::set_mass_limits() {
  
  bool found_xmin = false;
  double x, y, x_min, x_max;
  double mass, lowlim, highlim, asym, scale; 
  double M_min=1.0, M_max=2.4, dx=1.0e-4;

  for (size_t i=0; i<mdat.id_ns.size(); i++) {
    mass = mdat.mass_ns[i];
    lowlim = mdat.lowlim_ns[i];
    highlim = mdat.uplim_ns[i];
    asym = sqrt(highlim/lowlim);
    scale = like.get_scale(lowlim, highlim);
    for (x=M_min; x<=M_max; x+=dx) {
      y = like.asym_norm(mass-x, asym, scale);
      if (y>tol) {
        if (!found_xmin) {
          x_min=x;
          found_xmin=true;
        }
        x_max=x;
        if (x+dx>M_max) {
          low_m_ns.push_back(x_min);
          high_m_ns.push_back(x_max);
        }
      }
      else if (found_xmin) {
        low_m_ns.push_back(x_min);
        high_m_ns.push_back(x_max);
        break;
      }
    }
    found_xmin=false;
  }
  for (size_t i=0; i<mdat.id_wd.size(); i++) {
    mass = mdat.mass_wd[i];
    lowlim = mdat.lowlim_wd[i];
    highlim = mdat.uplim_wd[i];
    asym = sqrt(highlim/lowlim);
    scale = like.get_scale(lowlim, highlim);
    for (x=M_min; x<=M_max; x+=dx) {
      y = like.asym_norm(mass-x, asym, scale);
      if (y>tol) {
        if (!found_xmin) {
          x_min=x;
          found_xmin=true;
        }
        x_max=x;
        if (x+dx>M_max) {
          low_m_wd.push_back(x_min);
          high_m_wd.push_back(x_max);
        }
      }
      else if (found_xmin) {
        low_m_wd.push_back(x_min);
        high_m_wd.push_back(x_max);
        break;
      }
    }
    found_xmin=false;
  }
  for (size_t i=0; i<mdat.id_ms.size(); i++) {
    mass = mdat.mass_ms[i];
    lowlim = mdat.lim_ms[i];
    highlim = lowlim;
    asym = sqrt(highlim/lowlim);
    scale = like.get_scale(lowlim, highlim);
    for (x=M_min; x<=M_max; x+=dx) {
      y = like.asym_norm(mass-x, asym, scale);
      if (y>tol) {
        if (!found_xmin) {
          x_min=x;
          found_xmin=true;
        }
        x_max=x;
        if (x+dx>M_max) {
          low_m_ms.push_back(x_min);
          high_m_ms.push_back(x_max);
        }
      }
      else if (found_xmin) {
        low_m_ms.push_back(x_min);
        high_m_ms.push_back(x_max);
        break;
      }
    }
    found_xmin=false;
  }
}


void t_class::test_mass_limits() {

  for (size_t i=0; i<low_m_ns.size(); i++) {
    cout << fixed << setprecision(4);
    cout << "NS: i=" << i << "  low = " << low_m_ns[i] 
      << "\t high = " << high_m_ns[i];
    if (mdat.mass_ns[i]>low_m_ns[i] && mdat.mass_ns[i]<high_m_ns[i]) {
      cout << "\t low<mass<high" << endl; 
    }
    else cout << "\t Unbounded mass for " << mdat.id_ns[i] << endl;
  }
  for (size_t i=0; i<low_m_wd.size(); i++) {
    cout << fixed << setprecision(4);
    cout << "WD: i=" << i << "  low = " << low_m_wd[i] 
      << "\t high = " << high_m_wd[i];
    if (mdat.mass_wd[i]>low_m_wd[i] && mdat.mass_wd[i]<high_m_wd[i]) {
      cout << "\t low<mass<high" << endl; 
    }
    else cout << "\t Unbounded mass for " << mdat.id_wd[i] << endl;
  }
  for (size_t i=0; i<low_m_ms.size(); i++) {
    cout << fixed << setprecision(4);
    cout << "MS: i=" << i << "  low = " << low_m_ms[i] 
      << "\t high = " << high_m_ms[i];
    if (mdat.mass_ms[i]>low_m_ms[i] && mdat.mass_ms[i]<high_m_ms[i]) {
      cout << "\t low<mass<high" << endl; 
    }
    else cout << "\t Unbounded mass for " << mdat.id_ms[i] << endl;
  }
}


int main(void) {

  ofstream file; // file_ns, file_wd, file_ms;
  // file_ns.open("L_ns.dat");
  // file_wd.open("L_wd.dat");
  // file_ms.open("L_ms.dat");
  // file.open("test.dat");
 
  double res=0.0,    err=0.0;
  double res_ns=0.0, err_ns=0.0;
  double res_wd=0.0, err_wd=0.0;
  double res_ms=0.0, err_ms=0.0;
  
  t_class tc;
  likelihood like;
  mass_data mdat;
  mdat.load_data();

  tc.set_mass_limits();
  // tc.test_mass_limits();

  tc.set_dist_limits();
  
  mcarlo_vegas<> gm;

  size_t n_pars_ns = 3+mdat.id_ns.size();
  size_t n_pars_wd = 3+mdat.id_wd.size();
  size_t n_pars_ms = 3+mdat.id_ms.size();
  
  ubvector a_ns(n_pars_ns), b_ns(n_pars_ns);
  ubvector a_wd(n_pars_wd), b_wd(n_pars_wd);
  ubvector a_ms(n_pars_ms), b_ms(n_pars_ms);

  a_ns[0]=1.2;  b_ns[0]=1.8;
  a_ns[1]=0.5;  b_ns[1]=1.0;
  a_ns[2]=-1.0; b_ns[2]=1.0;
  for (size_t i=0; i<mdat.id_ns.size(); i++) {
    a_ns[i+3] = tc.low_m_ns[i];
    b_ns[i+3] = tc.high_m_ns[i];
    if (a_ns[i]>=b_ns[i]) {cout << "NS: Invalid limits!" << endl;}
  }
  a_wd[0]=1.2;  b_wd[0]=1.8;
  a_wd[1]=0.5;  b_wd[1]=1.0;
  a_wd[2]=-1.0; b_wd[2]=1.0;
  for (size_t i=0; i<mdat.id_wd.size(); i++) {
    a_wd[i+3] = tc.low_m_wd[i];
    b_wd[i+3] = tc.high_m_wd[i];
    if (a_wd[i]>=b_wd[i]) {cout << "WD: Invalid limits!" << endl;}
  }
  a_ms[0]=1.2;  b_ms[0]=1.8;
  a_ms[1]=0.5;  b_ms[1]=1.0;
  a_ms[2]=-1.0; b_ms[2]=1.0;
  for (size_t i=0; i<mdat.id_ms.size(); i++) {
    a_ms[i+3] = tc.low_m_ms[i];
    b_ms[i+3] = tc.high_m_ms[i];
    if (a_ms[i]>=b_ms[i]) {cout << "MS: Invalid limits!" << endl;}
  }

  multi_funct f = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight), &tc, _1, _2);
  multi_funct f_ns = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_ns), &tc, _1, _2);
  multi_funct f_wd = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_wd), &tc, _1, _2);
  multi_funct f_ms = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_ms), &tc, _1, _2);
 
  gm.n_points=100000;

  // cout << "n_params = " << like.n_params << endl;
  // cout << "n_dist_pars = " << like.n_dist_pars << endl;
  // cout << "id_ns.size() = " << mdat.id_ns.size() << endl;

  double x, y, m, c, d, u, l, mu, sigma, alpha, M;
  int k;

  /*cout << "Enter k=" ;
  cin >> k; */
  
  /*for (mu=0.5; mu<=2.5; mu+=0.1) {
    for (sigma=0.0; sigma<=1.0; sigma+=0.1) {
      for (alpha=-1.0; alpha<=1.0; alpha+=0.1) {
        for (M=0.0; M<=2.5; M+=0.1) {
          y = like.skew_norm(M, mu, sigma, alpha);
          cout << "mu=" << mu << "\t sigma=" << sigma << "\t alpha=" << alpha
            << "\t M=" << M << endl;
          if (y<tc.tol) {cout << "y = " << y << " < tol at mu=" << mu 
          << ",\t sigma=" << sigma << ",\t alpha=" << alpha << ",\t M=" << M << endl;}
        }
      }
    }
  } */

  gm.minteg_err(f_ns, n_pars_ns, a_ns, b_ns, res_ns, err_ns);
  gm.minteg_err(f_wd, n_pars_wd, a_wd, b_wd, res_wd, err_wd);
  gm.minteg_err(f_ms, n_pars_ms, a_ms, b_ms, res_ms, err_ms);
  
  cout << scientific;
  cout << "res_ns = " << res_ns << "\t err_ns = " << err_ns << endl;
  cout << "res_wd = " << res_wd << "\t err_wd = " << err_wd << endl;
  cout << "res_ms = " << res_ms << "\t err_ms = " << err_ms << endl;
 
  // file_ns.close();
  // file_wd.close();
  // file_ms.close();
  // file.close();

  return 0;
}
