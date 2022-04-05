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
    
    vector<double> low;
    vector<double> high;

    double m=0.0, w=0.0, s=0.0, M=1.0;

    static const int n_pts=100;
    const double dm = (2.5-0.5)/n_pts;
    const double dw = (1.0-0.0)/n_pts;
    const double ds = (1.0+1.0)/n_pts;
    const double dM = (2.3-1.0)/n_pts;

    double get_weight(size_t, const ubvector &);
    double get_weight_ns(size_t, const ubvector &);
    double get_weight_wd(size_t, const ubvector &);
    double get_weight_ms(size_t, const ubvector &);

    void set_limits();
    void mc_integrate();
};