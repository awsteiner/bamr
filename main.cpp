/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2016, Andrew W. Steiner
  
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
#include "bamr.h"

using namespace bamr;

int main(int argc, char *argv[]) {

#ifndef BAMR_NO_MPI
  // Init MPI
  MPI_Init(&argc,&argv);
#endif

  // Settings object
  settings set;

  // Main bamr object
  bamr::bamr_class b(set);

  // Run!
  b.run(argc,argv);
  
#ifndef BAMR_NO_MPI
  // Finalize MPI
  MPI_Finalize();
#endif

  return 0;
}
