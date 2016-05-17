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

using namespace std;
using namespace mcmc_namespace;

typedef boost::numeric::ublas::vector<double> ubvector;

/** \brief A very simple function to simulate
 */
class model_ex : public mcmc_namespace::model_demo {
  
public:

  model_ex() {
    nparams=1;
  }
  
  /** \brief Compute the function
   */
  virtual double compute_point(ubvector &pars, std::ofstream &scr_out,
			       int &success, ubvector &dat) {
    success=0;
    return exp(-pars[0]*pars[0]/2.0);
  }
    
  /// \name Functions for MCMC parameters
  //@{
  /** \brief Set the lower boundaries for all the parameters
   */
  virtual void low_limits(ubvector &pars) {
    pars[0]=-5.0;
    return;
  }

  /** \brief Set the upper boundaries for all the parameters
   */
  virtual void high_limits(ubvector &pars) {
    pars[0]=4.0;
    return;
  }
  
  /// Set up the parameter names
  virtual void param_names(std::vector<std::string> &names) {
    names.resize(1);
    names[0]="x";
    return;
  }
  
  /// Set up the parameter units
  virtual void param_units(std::vector<std::string> &units) {
    units.resize(1);
    return;
  }
  //@}

};

int main(int argc, char *argv[]) {

#ifndef MCMC_NO_MPI
  // Init MPI
  MPI_Init(&argc,&argv);
#endif

  // ---------------------------------------
  
  // Model object
  std::shared_ptr<model_ex> m(new model_ex);
  
  // MCMC object 
  mcmc_class<ubvector,model_ex> mcmc(m);

  // Run!
  mcmc.run(argc,argv);
  
  // ---------------------------------------

#ifndef MCMC_NO_MPI
  // Finalize MPI
  MPI_Finalize();
#endif

  return 0;
}
