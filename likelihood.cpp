#include "likelihood.h"
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace std::placeholders;

typedef boost::numeric::ublas::vector<double> ubvector;

class cls {

public:

  /* We'll use this to count the number of function
     evaluations required by the integration routines */
  int nf;

  // PDF of standard normal distribution N(0,1)
  double phi(double x) {
    return exp(-0.5*x*x) / sqrt(2.0*M_PI);
  }

  // CDF of standard N(0,1) in terms of erf(x)
  double Phi(double x) {
    return 0.5 * (1.0 + erf(x/sqrt(2.0)));
  }

  // Skewed Normal PDF [eq. 13, refs/kiziltan13]
  double SN(double x, double mu, double sigma, double alpha) {
    return 2.0 * phi((x-mu)/sigma) * Phi((x-mu)*alpha/sigma) / sigma;
  }

  // Asymmetric Normal PDF [eq. 14, refs/kiziltan13]
  double AN(double x, double c, double d) {
    double a = 2.0 / (d*(c+1.0/c));
    if (x>=0.0)
      return a * phi(x/(c*d));
    else
      return a * phi(c*x/d);
  }

  // This is the function to solve [see refs/calc.pdf]
  double f2solve(double x, double &l, double &u) {
    double c = sqrt(u/l);
    return c*c*erf(u/(sqrt(2.0)*c*x)) - erf(-c*l/(sqrt(2.0)*x))
      - 0.68*(c*c+1.0);
  }

  // Derivative of the function to solve (for use with root_stef)
  double df2solve(double x, double &l, double &u) {
    double c = sqrt(u/l);
    double a = sqrt(2.0/M_PI) * c / x / x;
    return a*l*exp(-pow(u/(sqrt(2.0)*c*x), 2.0))
      - a*u*exp(-pow(c*l/sqrt(2.0)/x, 2.0));
  }

  // The solver that calculates parameters d_j, given c_j = sqrt(u_j/l_j)
  double calculate_d(double l, double u) {
    cout.setf(ios::scientific);

    test_mgr t;
    // Only print something out if one of the tests fails
    t.set_output_level(1);

    // The solver, specifying the function type: funct<double>
    root_brent_gsl<> solver;

    cls c;

    /* This is the code that allows specification of class member
       functions as functions to solve. This approach avoids the use of
       static variables and functions and multiple inheritance at the
       expense of a little overhead. We need to provide the address of
       an instantiated object and the address of the member function. */
    funct f2 = bind(mem_fn<double(double, double &, double &)>
			 (&cls::f2solve), &c, _1, ref(l), ref(u));

    /* funct df2 = bind(mem_fn<double(double, double &, double &)>
       (&cls::df2solve), &c, _1, ref(l), ref(u)); */

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
  double likelihood_ns(double mu, double sigma, double alpha, ubvector M_j) {
    double m_j, l_j, u_j, c_j, d_j, w_j, L_j, L=1.0;
    for (int j=0; j<26; j++) {
      m_j = nsns[j][0]; u_j = nsns[j][1];
      l_j = nsns[j][2]; c_j = sqrt(u_j/l_j);
      d_j = calculate_d(l_j, u_j);
      w_j = m_j - M_j.at(j);
      L_j = AN(w_j, c_j, d_j) * SN(M_j, mu, sigma, alpha);
      if (L_j<=1.0e-6) L_j = 1.0; // Ignore small likelihoods
      L *= L_j; // Note: This is not log-likelihood
    }
    return L;
  }

  // The likelihood function for NS-WD (see refs/method.pdf)
  double likelihood_wd(double mu, double sigma, double alpha, ubvector M_j) {
    double m_j, l_j, u_j, c_j, d_j, w_j, L_j, L=1.0;
    for (int j=0; j<38; j++) {
      m_j = nswd[j][0]; u_j = nswd[j][1];
      l_j = nswd[j][2]; c_j = sqrt(u_j/l_j);
      d_j = calculate_d(l_j, u_j);
      w_j = m_j - M_j.at(j);
      L_j = AN(w_j, c_j, d_j) * SN(M_j, mu, sigma, alpha);
      if (L_j<=1.0e-6) L_j = 1.0;
      L *= L_j;
    }
    return L;
  }

  // The likelihood function for NS-MS (see refs/method.pdf)
  double likelihood_ms(double mu, double sigma, double alpha, ubvector M_j) {
    double m_j, l_j, u_j, c_j, d_j, w_j, L_j, L=1.0;
    for (int j=0; j<16; j++) {
      m_j = nsms[j][0]; u_j = nsms[j][1];
      l_j = u_j; c_j = sqrt(u_j/l_j); // Symmetric 68% limits
      d_j = calculate_d(l_j, u_j);
      w_j = m_j - M_j.at(j);
      L_j = AN(w_j, c_j, d_j) * SN(M_j, mu, sigma, alpha);
      if (L_j<=1.0e-6) L_j = 1.0;
      L *= L_j;
    }
    return L;
  }

  /*
   */

  /** \brief Set the vec_index object with the parameters from
      the mass data

      This function will be called by bamr to fill the \c pvi
      object with the mass parameters from the data set.
  */
  void set_params(vec_index &pvi) {

    for(size_t i=0;i<26;i++) {
      string mass_par=((string)"M_")+nsns_id[i];
      pvi.append(mass_par);
    }
    
    return;
  }
  
  /// The combined likelihood function to be calculated
  double likelihood(const ubvector &pars, vec_index &pvi) {

    double mu_NSNS=pars[pvi["mu_NSNS"]];
    //double M_1913a=pars[pvi["M_1913a"]];

    for(size_t i=0;i<26;i++) {
      double mass_par=pars[pvi[((string)"M_")+nsns_id[i]]];
      double temp=SN(nsns[i][0],mass_par);
    }

    default_random_engine seed;
    uniform_real_distribution<double> fmean(0.5, 2.5);
    uniform_real_distribution<double> fsigma(0.0, 1.0);
    uniform_real_distribution<double> falpha(-1.0, 1.0);
    uniform_real_distribution<double> fmass(1.0, 2.3);

    int j, k, n=94;

    //ubvector M_j;
    //M_j.reserve(n);

    double mu, sigma, alpha, L_ns, L_wd, L_ms, L=1.0;

    get_data(); // Load data from likelihood.h

    for (k=1; k<=1e5; k++) {
      mu = fmean(seed);
      sigma = fsigma(seed);
      alpha = falpha(seed);

      for(j=0; j<94; j++) {
	M_j.push_back(fmass(seed));
      }
      // Calculate likelihood for each population
      L_ns = likelihood_ns(mu, sigma, alpha, M_j);
      L_wd = likelihood_wd(mu, sigma, alpha, M_j);
      L_ms = likelihood_ms(mu, sigma, alpha, M_j);
      // Multiply the likelihoods
      L *= L_ns * L_wd * L_ms;
    }
  }
  return log(L);
};


int main() {
  cls c;
  c.likelihood();
  return 0;
}
