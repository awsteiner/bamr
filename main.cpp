/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2017, Andrew W. Steiner
  
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
#ifndef BAMR_NO_MPI
#include <mpi.h>
#endif

#include "models.h"
#include "bamr_class.h"
#include "mcmc_bamr.h"

using namespace std;
using namespace bamr;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  // Set error handler for this thread
  o2scl::err_hnd_cpp error_handler;
  o2scl::err_hnd=&error_handler;
  
#ifndef BAMR_NO_MPI
  // Init MPI
  MPI_Init(&argc,&argv);
#endif

  // Main bamr object
  bamr::mcmc_bamr b(1);

  // Set starting time
#ifdef O2SCL_MPI
  b.mpi_start_time=MPI_Wtime();
#else
  b.mpi_start_time=time(0);
#endif

  b.setup_cli();
  
  // Set command-line args
  b.cl_args.clear();
  for(int i=0;i<argc;i++) {
    b.cl_args.push_back(argv[i]);
  }

  b.cl.prompt="bamr> ";
  b.cl.run_auto(argc,argv);

#ifndef BAMR_NO_MPI
  // Finalize MPI
  MPI_Finalize();
#endif

  return 0;
}
