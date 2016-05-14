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

/** \brief Desc
 */
class model_ex : public mcmc_namespace::default_model {
  
public:
  
  /** \brief Desc
   */
  virtual double compute_point(ubvector &pars, std::ofstream &scr_out,
			       int &success, ubvector &dat) {
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
    pars[0]=5.0;
    return;
  }

  /// Return the name of parameter with index \c i
  virtual std::string param_name(size_t i) {
    return "x";
  }

  /// Return the unit of parameter with index \c i
  virtual std::string param_unit(size_t i) {
    return "";
  }
  //@}

};

int main(int argc, char *argv[]) {

  // ---------------------------------------
  // Init MPI
  
#ifndef BAMR_NO_MPI
  MPI_Init(&argc,&argv);
#endif

  // ---------------------------------------
  // Main bamr object 

  std::shared_ptr<model_ex> m(new model_ex);
  mcmc_class<ubvector,model_ex> mcmc(m);
  mcmc.run(argc,argv);
  
  // ---------------------------------------
  // Finalize MPI

#ifndef BAMR_NO_MPI
  MPI_Finalize();
#endif

  return 0;
}
