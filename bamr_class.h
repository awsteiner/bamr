/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
/** \file bamr_class.h
    \brief Definition of main bamr class
*/
#ifndef BAMR_CLASS_H
#define BAMR_CLASS_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#ifndef BAMR_NO_MPI
#include <mpi.h>
#endif

#include <o2scl/rng_gsl.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/cholesky.h>
#include <o2scl/mcmc_para.h>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "nstar_cold2.h"
#include "models.h"

/** \brief Main namespace
    
    The bamr namespace which holds all classes and functions.
    
    This file is documented in bamr.h .
*/
namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  typedef std::function<int(size_t,const ubvector &, double &,
			    model_data &)> point_funct;
  
  typedef std::function<int(const ubvector &,double,
			    std::vector<double> &,model_data &)> fill_funct;

  /** \brief Desc
   */
  class bamr_class {
    
  public:

    /// The Schwarzchild radius in km
    double schwarz_km;
    
    /// Neutron star data
    ns_data nsd;

    /// Pointer to settings object
    settings *setp;
    
    /// Model object
    std::shared_ptr<model> mod;

    /// Number of parameters
    size_t n_params;

    bamr_class() {
      schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;
    }

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual int compute_point(const ubvector &pars, std::ofstream &scr_out, 
				 double &weight, model_data &dat);

    /** \brief Fill vector in <tt>line</tt> with data from the
	current Monte Carlo point
    */
    virtual int fill(const ubvector &pars, double weight, 
		      std::vector<double> &line, model_data &dat);
    
  };

}

#endif
