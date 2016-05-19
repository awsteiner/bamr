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
  
  /** \brief Compute the value of the function to simulate
      at the parameter point \c pars
  */
  virtual double compute_point(ubvector &pars, std::ofstream &scr_out,
			       int &success, ubvector &dat) {
    success=0;
    return exp(-pars[0]*pars[0]/2.0);
  }
  
  /** \brief Specify the initial point
   */
  virtual void initial_point(ubvector &pars) {
    pars[0]=-0.01;
    return;
  }
  
  /** \brief Set parameter information
   */
  virtual void get_param_info(std::vector<std::string> &names,
			      std::vector<std::string> &units,
			      ubvector &low, ubvector &high) {
    names[0]="x";
    units[0]="";
    low[0]=-5.0;
    high[0]=5.0;
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
  mcmc_base<ubvector,model_ex> mcmc2(m);

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
