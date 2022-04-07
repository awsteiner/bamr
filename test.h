#include "likelihood.h"
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>

using namespace std;
using namespace o2scl;

class t_class {
    
    public:

    t_class() {
        mdat.load_data();
        like.get_params();
    }

    virtual ~t_class() {};

    likelihood like;
    mass_data mdat;
    
    vector<double> low_m_ns;
    vector<double> high_m_ns;
    vector<double> low_m_wd;
    vector<double> high_m_wd;
    vector<double> low_m_ms;
    vector<double> high_m_ms;

    double low_mean_ns;
    double high_mean_ns;
    double low_mean_wd;
    double high_mean_wd;
    double low_mean_ms;
    double high_mean_ms;

    double low_width_ns;
    double high_width_ns;
    double low_width_wd;
    double high_width_wd;
    double low_width_ms;
    double high_width_ms;
    
    double low_skewness_ns;
    double high_skewness_ns;
    double low_skewness_wd;
    double high_skewness_wd;
    double low_skewness_ms;
    double high_skewness_ms;


    // double mu=0.0, sigma=0.0, alpha=0.0;

    static const int n_pts=100;
    const double dm = (2.5-0.5)/n_pts;
    const double dw = (1.0-0.0)/n_pts;
    const double ds = (1.0+1.0)/n_pts;
    const double dM = (2.3-1.0)/n_pts;
    const double tol = 1.0e-10;

    double get_weight(size_t, const ubvector &);
    double get_weight_ns(size_t, const ubvector &);
    double get_weight_wd(size_t, const ubvector &);
    double get_weight_ms(size_t, const ubvector &);

    void set_mass_limits();
    void set_dist_limits();
    void test_mass_limits();
    void test_dist_limits();
};