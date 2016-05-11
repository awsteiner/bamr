/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016, Andrew W. Steiner
  
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
#include "mcmc.h"

class model {

  virtual void init() {
    return;
  };
    
  /** \brief Desc
   */
  size_t nparams;

  /** \brief Desc
   */
  virtual double compute_point(ubvector &pars, std::ofstream &scr_out,
			       int &success, ubvector &dat)=0;
    
  /** \brief Desc
   */
  virtual void first_point(ubvector &pars) {
    ubvector low(nparams), high(nparams);
    low_limits(low);
    high_limits(high);
    for(size_t i=0;i<nparams;i++) {
      pars[i]=(low[i]+high[i])/2.0;
    }
    return;
  }
    
  /// \name Functions for MCMC parameters
  //@{
  /** \brief Set the lower boundaries for all the parameters
   */
  virtual void low_limits(ubvector &pars)=0;

  /** \brief Set the upper boundaries for all the parameters
   */
  virtual void high_limits(ubvector &pars)=0;

  /// Return the name of parameter with index \c i
  virtual std::string param_name(size_t i) {
  }

  /// Return the unit of parameter with index \c i
  virtual std::string param_unit(size_t i) {
  }
  //@}


};


int main(void) {

  // ---------------------------------------
  // Init MPI
  
#ifndef BAMR_NO_MPI
  MPI_Init(&argc,&argv);
#endif

  // ---------------------------------------
  // Main bamr object 
  
  mcmc_class<ubvector,model> mcmc;
  mcmc.run(argc,argv);
  
  // ---------------------------------------
  // Finalize MPI

#ifndef BAMR_NO_MPI
  MPI_Finalize();
#endif

  return;
}

#endif
