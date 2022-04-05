#include "test.h"

// The likelihood function for NS-NS (see refs/method.pdf)
double t_class::get_weight_ns(size_t n_vars, const ubvector &x) {
    
  double mean = x[0];
  double width = x[1];
  double skewness = x[2];

  double M_star, mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;
  
  for (size_t i=0; i<1; i++) {
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
  double M_star = M; 

  double mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<mdat.id_wd.size(); i++) {
    mass = mdat.mass_wd[i]; 
    uplim = mdat.uplim_wd[i];
    lowlim = mdat.lowlim_wd[i];
    asym = sqrt(uplim/lowlim);
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
  double M_star = M; 
  
  double mass, lowlim, uplim, asym, scale, wgt_star, wgt=1.0;

  for (size_t i=0; i<mdat.id_ms.size(); i++) {
    mass = mdat.mass_ms[i]; 
    uplim = mdat.lim_ms[i];
    lowlim = uplim; // Symmetric 68% limits
    asym = sqrt(uplim/lowlim); 
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

void t_class::set_limits() {

}

int main(void) {

  ofstream file; // file_ns, file_wd, file_ms;
  // file_ns.open("L_ns.dat");
  // file_wd.open("L_wd.dat");
  // file_ms.open("L_ms.dat");
  file.open("test.dat");
 
  double res=0.0,    err=0.0;
  double res_ns=0.0, err_ns=0.0;
  double res_wd=0.0, err_wd=0.0;
  double res_ms=0.0, err_ms=0.0;
  
  mass_data mdat;
  mdat.load_data();

  size_t n_dpars = 3;
  size_t n_mpars = mdat.id_ns.size();
  size_t n_pars = n_dpars + n_mpars;

  mcarlo_vegas<> gm;
  ubvector a(4), b(4);

  a[0]=0.5;  b[0]=2.5;
  a[1]=0.0;  b[1]=1.0;
  a[2]=-1.0; b[2]=1.0;
  a[3]=mdat.mass_ns[0]-0.155;  b[3]=mdat.mass_ns[0];

  /*for (size_t i=0; i<n_mpars; i++) {
    a[i+n_dpars] = 1.4;
    b[i+n_dpars] = 1.6;
  }*/

  t_class tc;

  multi_funct f = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight), &tc, _1, _2);
  multi_funct f_ns = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_ns), &tc, _1, _2);
  multi_funct f_wd = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_wd), &tc, _1, _2);
  multi_funct f_ms = bind(mem_fn<double(size_t, const ubvector &)>
    (&t_class::get_weight_ms), &tc, _1, _2);
 
  gm.n_points=1000000;

  // cout << "n_params = " << like.n_params << endl;
  // cout << "n_dist_pars = " << like.n_dist_pars << endl;
  // cout << "id_ns.size() = " << mdat.id_ns.size() << endl;

  double x, y, c, d, u, l;
  int k;

  // cout << "Enter k=" ;
  // cin >> k;

  /*for (x=0.0; x<=2.5; x+=0.001) {
    u = mdat.uplim_ns[k];
    l = mdat.lowlim_ns[k]; 
    c = sqrt(u/l);
    d = like.get_scale(l, u);
    y = like.asym_norm(x, c, d);
    file << x << "\t" << y << endl;
    cout << "x = " << x << "\t AN(x) = " << y << endl;
    M += dM;
  } 
  cout << "scale = " << d << endl; */

  // gm.minteg_err(f, 3, a, b, res, err);
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