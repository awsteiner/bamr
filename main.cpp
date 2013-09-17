/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013, Andrew W. Steiner
  
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
#include <mpi.h>

#include "bamr.h"

int main(int argc, char *argv[]) {

  // ---------------------------------------
  // Init MPI
  
  MPI_Init(&argc,&argv);

  // ---------------------------------------
  // Main bamr object 
  
  o2scl::bamr b;
  b.run(argc,argv);
  
  // ---------------------------------------
  // Finalize MPI

  MPI_Finalize();

  return 0;
}
