/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2020, Andrew W. Steiner
  
  This file is part of Bamr.
  
  Bamr is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Bamr is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Bamr. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
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

