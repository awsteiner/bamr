/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2018, Andrew W. Steiner
  
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

#ifdef BAMR_MPI
#include <mpi.h>
#endif

#include "nstar_cold2.h"
#include "models.h"

/** \brief Main namespace
    
    The bamr namespace which holds all classes and functions.
    
    This file is documented in bamr_class.h .
*/
namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  typedef std::function<int(size_t,const ubvector &, double &,
			    model_data &)> point_funct;
  
  typedef std::function<int(const ubvector &,double,
			    std::vector<double> &,model_data &)> fill_funct;

  /** \brief Compute neutron star structure for each MCMC point

      There will be a number of instances of this class equal to the
      number of OpenMP threads, all of which share the same neutron
      star data (in \ref nsd) and settings (in \ref set). They will
      each have their own model object, stored in \ref mod, which is
      set in \ref mcmc_bamr::set_model() .
  */
  class bamr_class {
    
  public:

    /// The Schwarzchild radius in km
    double schwarz_km;
    
    /// Pointer to neutron star data
    std::shared_ptr<const ns_data> nsd;
    
    /// Pointer to settings object
    std::shared_ptr<const settings> set;
    
    /// Model object
    std::shared_ptr<model> mod;

    /// Model type string
    std::string model_type;
    
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
