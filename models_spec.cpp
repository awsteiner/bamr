#include "models_spec.h"

#include "bint.h"

#ifdef BAMR_MPI
#include <mpi.h>
#endif

#include <o2scl/hdf_file.h>
#include <o2scl/lu.h>
#include <o2scl/columnify.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/anneal_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace o2scl_linalg;
using namespace bamr;

void qmc_spec_ligo::compute_eos(const ubvector &e, int &success,
				std::ofstream &scr_out, model_data &dat) {
  return;
}

