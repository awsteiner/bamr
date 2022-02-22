#include "likelihood.h"


  /* We'll use this to count the number of function
     evaluations required by the integration routines */
  // int nf;

void like::load_data() {
  data_t dt;
  dt.mass.push_back(1.4);
  /* Data from NS-NS/DNS (some with asymmetric 68% limits)
     [source: Alsing et al. 2018] */
  nsns_id[0] = "J0453+1559 pulsar"; nsns_sid[0] = "J0453p";
  nsns[0][0] = 1.559; nsns[0][1] = 0.004; nsns[0][2] = 0.004;

  nsns_id[1] = "J0453+1559 comp."; nsns_sid[1] = "J0453c";
  nsns[1][0] = 1.174; nsns[1][1] = 0.004; nsns[1][2] = 0.004;

  nsns_id[2] = "J1906+0746 pulsar"; nsns_sid[2] = "J1906p";
  nsns[2][0] = 1.291; nsns[2][1] = 0.011; nsns[2][2] = 0.011;

  nsns_id[3] = "J1906+0746 comp."; nsns_sid[3] = "J1906c";
  nsns[3][0] = 1.322; nsns[3][1] = 0.011; nsns[3][2] = 0.011;

  nsns_id[4] = "B1534+12 pulsar"; nsns_sid[4] = "B1534p"; 
  nsns[4][0] = 1.3332; nsns[4][1] = 0.001; nsns[4][2] = 0.001;

  nsns_id[5] = "B1534+12 comp."; nsns_sid[5] = "B1534c";
  nsns[5][0] = 1.3452; nsns[5][1] = 0.001; nsns[5][2] = 0.001;

  nsns_id[6] = "B1913+16 pulsar"; nsns_sid[6] = "B1913p";
  nsns[6][0] = 1.4398; nsns[6][1] = 0.0002; nsns[6][2] = 0.0002;

  nsns_id[7] = "B1913+16 comp."; nsns_sid[7] = "B1913c";
  nsns[7][0] = 1.3886; nsns[7][1] = 0.0002; nsns[7][2] = 0.0002;

  nsns_id[8] = "B2127+11C pulsar"; nsns_sid[8] = "B2127p";
  nsns[8][0] = 1.358; nsns[8][1] = 0.01; nsns[8][2] = 0.01;

  nsns_id[9] = "B2127+11C comp."; nsns_sid[9] = "B2127c";
  nsns[9][0] = 1.354; nsns[9][1] = 0.01; nsns[9][2] = 0.01;

  nsns_id[10] = "J0737-3039A"; nsns_sid[10] = "J0737A";
  nsns[10][0] = 1.3381; nsns[10][1] = 0.0007; nsns[10][2] = 0.0007;

  nsns_id[11] = "J0737-3039B"; nsns_sid[11] = "J0737B";
  nsns[11][0] = 1.2489; nsns[11][1] = 0.0007; nsns[11][2] = 0.0007;

  nsns_id[12] = "J1756-2251 pulsar"; nsns_sid[12] = "J1756p";
  nsns[12][0] = 1.312; nsns[12][1] = 0.017; nsns[12][2] = 0.017;

  nsns_id[13] = "J1756-2251 comp."; nsns_sid[13] = "J1756c";
  nsns[13][0] = 1.258; nsns[13][1] = 0.017; nsns[13][2] = 0.017;

  nsns_id[14] = "J1807-2500B pulsar"; nsns_sid[14] = "J1807p";
  nsns[14][0] = 1.3655; nsns[14][1] = 0.0021; nsns[14][2] = 0.0021;

  nsns_id[15] = "J1807-2500B comp."; nsns_sid[15] = "J1807c";
  nsns[15][0] = 1.2064; nsns[15][1] = 0.002; nsns[15][2] = 0.002;

  // [source: Kiziltan et al. 2013]
  nsns_id[16] = "J1518+4904 pulsar"; nsns_sid[16] = "J1518p";
  nsns[16][0] = 1.56; nsns[16][1] = 0.13; nsns[16][2] = 0.44;

  nsns_id[17] = "J1518+4904 comp."; nsns_sid[17] = "J1518c";
  nsns[17][0] = 1.05; nsns[17][1] = 0.45; nsns[17][2] = 0.11;

  nsns_id[18] = "J1811-1736 pulsar"; nsns_sid[18] = "J1811p";
  nsns[18][0] = 1.56; nsns[18][1] = 0.24; nsns[18][2] = 0.45;

  nsns_id[19] = "J1811-1736 comp."; nsns_sid[19] = "J1811c";
  nsns[19][0] = 1.12; nsns[19][1] = 0.47; nsns[19][2] = 0.13;

  nsns_id[20] = "J1829+2456 pulsar"; nsns_sid[20] = "J1829p";
  nsns[20][0] = 1.20; nsns[20][1] = 0.12; nsns[20][2] = 0.46;

  nsns_id[21] = "J1829+2456 comp."; nsns_sid[21] = "J1829c";
  nsns[21][0] = 1.40; nsns[21][1] = 0.46; nsns[21][2] = 0.12;

  /* Data from NS-WD (some with asymmetric 68% limits)
     [source: Alsing et al. 2018, Kiziltan et al. 2013] */
  nswd_id[0] = "J2045+3633"; nswd_sid[0] = "J2045";
  nswd[0][0] = 1.33; nswd[0][1] = 0.3; nswd[0][2] = 0.3;

  nswd_id[1] = "J2053+4650"; nswd_sid[1] = "J2053";
  nswd[1][0] = 1.40; nswd[1][1] = 0.21; nswd[1][2] = 0.21;

  nswd_id[2] = "J1713+0747"; nswd_sid[2] = "J1713"; 
  nswd[2][0] = 1.35; nswd[2][1] = 0.07; nswd[2][2] = 0.07;

  nswd_id[3] = "B1855+09"; nswd_sid[3] = "B1855";
  nswd[3][0] = 1.37; nswd[3][1] = 0.13; nswd[3][2] = 0.13;

  nswd_id[4] = "J0751+1807"; nswd_sid[4] = "J0751";
  nswd[4][0] = 1.72; nswd[4][1] = 0.07; nswd[4][2] = 0.07;

  nswd_id[5] = "J1141-6545"; nswd_sid[5] = "J1141";
  nswd[5][0] = 1.27; nswd[5][1] = 0.01; nswd[5][2] = 0.01;

  nswd_id[6] = "J1738+0333"; nswd_sid[6] = "J1738";
  nswd[6][0] = 1.47; nswd[6][1] = 0.07; nswd[6][2] = 0.07;

  nswd_id[7] = "J1614-2230"; nswd_sid[7] = "J1614";
  nswd[7][0] = 1.908; nswd[7][1] = 0.016; nswd[7][2] = 0.016;

  nswd_id[8] = "J0348+0432"; nswd_sid[8] = "J0348";
  nswd[8][0] = 2.01; nswd[8][1] = 0.04; nswd[8][2] = 0.04;

  nswd_id[9] = "J2222-0137"; nswd_sid[9] = "J2222";
  nswd[9][0] = 1.76; nswd[9][1] = 0.06; nswd[9][2] = 0.06;

  nswd_id[10] = "J2234+0611"; nswd_sid[10] = "J2234";
  nswd[10][0] = 1.393; nswd[10][1] = 0.013; nswd[10][2] = 0.013;

  nswd_id[11] = "J1949+3106"; nswd_sid[11] = "J1949";
  nswd[11][0] = 1.47; nswd[11][1] = 0.43; nswd[11][2] = 0.43;

  nswd_id[12] = "J1012+5307"; nswd_sid[12] = "J1012";
  nswd[12][0] = 1.83; nswd[12][1] = 0.11; nswd[12][2] = 0.11;

  nswd_id[13] = "J0437-4715"; nswd_sid[13] = "J0437";
  nswd[13][0] = 1.44; nswd[13][1] = 0.07; nswd[13][2] = 0.07;

  nswd_id[14] = "J1909-3744"; nswd_sid[14] = "J1909";
  nswd[14][0] = 1.48; nswd[14][1] = 0.03; nswd[14][2] = 0.03;

  nswd_id[15] = "J1802-2124"; nswd_sid[15] = "J1802";
  nswd[15][0] = 1.24; nswd[15][1] = 0.11; nswd[15][2] = 0.11;

  nswd_id[16] = "J1911-5958A"; nswd_sid[16] = "J1911";
  nswd[16][0] = 1.34; nswd[16][1] = 0.08; nswd[16][2] = 0.08;

  nswd_id[17] = "J2043+1711"; nswd_sid[17] = "J2043";
  nswd[17][0] = 1.38; nswd[17][1] = 0.13; nswd[17][2] = 0.13;

  nswd_id[18] = "J0337+1715"; nswd_sid[18] = "J0337";
  nswd[18][0] = 1.4378; nswd[18][1] = 0.0013; nswd[18][2] = 0.0013;

  nswd_id[19] = "J1946+3417"; nswd_sid[19] = "J1946";
  nswd[19][0] = 1.828; nswd[19][1] = 0.022; nswd[19][2] = 0.022;

  nswd_id[20] = "J1918-0642"; nswd_sid[20] = "J1918";
  nswd[20][0] = 1.29; nswd[20][1] = 0.1; nswd[20][2] = 0.1;

  nswd_id[21] = "J1600-3053"; nswd_sid[21] = "J1600";
  nswd[21][0] = 2.3; nswd[21][1] = 0.7; nswd[21][2] = 0.7;

  // [source: Kiziltan et al. 2013]
  nswd_id[22] = "J0621+1002"; nswd_sid[22] = "J0621";
  nswd[22][0] = 1.70; nswd[22][1] = 0.10; nswd[22][2] = 0.17;

  nswd_id[23] = "B2303+46"; nswd_sid[23] = "B2303";
  nswd[23][0] = 1.38; nswd[23][1] = 0.06; nswd[23][2] = 0.10;

  nswd_id[24] = "J0024-7204H"; nswd_sid[24] = "J0024";
  nswd[24][0] = 1.48; nswd[24][1] = 0.03; nswd[24][2] = 0.06;

  nswd_id[25] = "J0514-4002A"; nswd_sid[25] = "J0514";
  nswd[25][0] = 1.49; nswd[25][1] = 0.04; nswd[25][2] = 0.27;

  nswd_id[26] = "B1516+02B"; nswd_sid[26] = "B1516";
  nswd[26][0] = 2.10; nswd[26][1] = 0.19; nswd[26][2] = 0.19;

  nswd_id[27] = "J1748-2446I"; nswd_sid[27] = "J1748I";
  nswd[27][0] = 1.91; nswd[27][1] = 0.02; nswd[27][2] = 0.10;

  nswd_id[28] = "J1748-2446J"; nswd_sid[28] = "J1748J";
  nswd[28][0] = 1.79; nswd[28][1] = 0.02; nswd[28][2] = 0.10;

  nswd_id[29] = "B1802-07"; nswd_sid[29] = "B1802";
  nswd[29][0] = 1.26; nswd[29][1] = 0.08; nswd[29][2] = 0.17;

  nswd_id[30] = "B1911-5958A"; nswd_sid[30] = "B1911";
  nswd[30][0] = 1.40; nswd[30][1] = 0.16; nswd[30][2] = 0.10;

  // [source: Cromartie et al. 2020]
  nswd_id[31] = "J0740+6620"; nswd_sid[31] = "J0740";
  nswd[31][0] = 2.14; nswd[31][1] = 0.20; nswd[31][2] = 0.18;

  /* Data from NS-MS (some from X-ray/Optical, all with symmetric 68% limits)
     [source: Alsing et al. 2018] */
  nsms_id[0]  = "4U1700-377"; nsms_sid[0] = "4U1700";
  nsms[0][0]  = 1.96;  nsms[0][1] = 0.19; 

  nsms_id[1]  = "Cyg X-2"; nsms_sid[1] = "CygX2";
  nsms[1][0]  = 1.71;  nsms[1][1] = 0.21; 

  nsms_id[2]  = "SMC X-1"; nsms_sid[2] = "SMCX1";
  nsms[2][0]  = 1.21;  nsms[2][1] = 0.12; 
  
  nsms_id[3]  = "Cen X-3"; nsms_sid[3] = "CenX3";
  nsms[3][0]  = 1.57;  nsms[3][1] = 0.16; 
  
  nsms_id[4]  = "XTE J2123-058"; nsms_sid[4] = "XTEJ2123";
  nsms[4][0]  = 1.53;  nsms[4][1] = 0.42; 
  
  nsms_id[5]  = "4U 1822-371"; nsms_sid[5] = "4U1822";
  nsms[5][0]  = 1.96;  nsms[5][1] = 0.36; 
  
  nsms_id[6]  = "OAO 1657-415"; nsms_sid[6] = "OAO1657";
  nsms[6][0]  = 1.74;  nsms[6][1] = 0.3;  
  
  nsms_id[7]  = "J013236.7+303228"; nsms_sid[7] = "J013236";
  nsms[7][0]  = 2.0;   nsms[7][1] = 0.4;  
  
  nsms_id[8]  = "Vela X-1"; nsms_sid[8] = "VelaX1";
  nsms[8][0]  = 2.12;  nsms[8][1] = 0.16; 
  
  nsms_id[9]  = "4U1538-522"; nsms_sid[9] = "4U1538";
  nsms[9][0]  = 1.02;  nsms[9][1] = 0.17; 
  
  nsms_id[10] = "LMC X-4"; nsms_sid[10] = "LMCX4";
  nsms[10][0] = 1.57;  nsms[10][1] = 0.11; 
  
  nsms_id[11] = "Her X-1"; nsms_sid[11] = "HerX1";
  nsms[11][0] = 1.073; nsms[11][1] = 0.36; 
  
  nsms_id[12] = "2S 0921-630"; nsms_sid[12] = "2S0921";
  nsms[12][0] = 1.44;  nsms[12][1] = 0.1;  
  
  nsms_id[13] = "EXO 1722-363"; nsms_sid[13] = "EXO1722";
  nsms[13][0] = 1.91;  nsms[13][1] = 0.45; 
  
  nsms_id[14] = "SAX J1802.7-2017"; nsms_sid[14] = "SAXJ1802";
  nsms[14][0] = 1.57;  nsms[14][1] = 0.25; 
  
  nsms_id[15] = "XTE J1855-026"; nsms_sid[15] = "XTEJ1855";
  nsms[15][0] = 1.41;  nsms[15][1] = 0.24; 
  
  /* Data from NS-MS (all with symmetric 68% limits)
     [source: Alsing et al. 2018]
     nsms_rd[0][0] = 1.58;   nsms_rd[0][1] = 0.34; 
     nsms_id[0] = "J0045-7319";
     nsms_rd[1][0] = 1.71;   nsms_rd[1][1] = 0.16;  
     nsms_id[1] = "J1023+0038";
     nsms_rd[2][0] = 1.666;  nsms_rd[2][1] = 0.01;  
     nsms_id[2] = "J1903+0327"; */
}

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


int main() {
  like c;
  c.calc_likelihood();
  return 0;
}
