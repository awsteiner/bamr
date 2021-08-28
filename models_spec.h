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
#ifndef MODELS_SPEC_H
#define MODELS_SPEC_H

#include <iostream>

#include <gsl/gsl_qrng.h>

#include <o2scl/eos_had_schematic.h>
#include <o2scl/poly.h>
#include <o2scl/eos_had_schematic.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/root_cern.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/interpm_idw.h>

#include "nstar_cold2.h"
#include "models.h"
#include "invert_tov.h"
#include "Faddeeva.h"

namespace bamr {

  class qmc_spec_ligo : public qmc_threep {

  public:

    qmc_threep_ligo(std::shared_ptr<const settings> s,
		    std::shared_ptr<const ns_data> n) {
    }
    
    /// Parameter for transition density
    o2scl::cli::parameter_double p_rho_trans;
    
    /** \brief Set parameter information
     */
    virtual void get_param_info(std::vector<std::string> &names,
				std::vector<std::string> &units,
				ubvector &low, ubvector &high) {
      
      names={"a","alpha","param_S","param_L","index1","trans1",
	     "index2","trans2","index3","M_chirp","q"};

      units={"MeV","","MeV","MeV","","1/fm^4","","1/fm^4","","Msun",""};
  
      low.resize(n_eos_params+nsd->n_sources);
      // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
      // this should allow S values as small as 28.5
      low[0]=12.5;
      // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
      low[1]=0.47;
      low[2]=29.5;
      low[3]=30.0;
  
      low[4]=0.2;
      low[5]=0.75;
      low[6]=0.2;
      low[7]=0.75;
      low[8]=0.2;

      low[9]=1.0;
      low[10]=0.0;
    
      high.resize(n_eos_params+nsd->n_sources);
      high[0]=13.5;
      high[1]=0.53;
      high[2]=36.1;
      high[3]=70.0;

      high[4]=8.0;
      high[5]=8.0;
      high[6]=8.0;
      high[7]=8.0;
      high[8]=8.0;
      high[9]=1.5;
      high[10]=2.5;
    
      // Go to the parent which takes care of the data-related
      // parameters
      model::get_param_info(names,units,low,high);

    }
    
    /** \brief Specify the initial point
     */
    virtual void initial_point(ubvector &params) {
      return;
    }

    /** \brief Setup model parameters */
    virtual void setup_params(o2scl::cli &cl) {
      p_rho_trans.d=&rho_trans;
      p_rho_trans.help="Transition from neutron matter to polytropes.";
      cl.par_list.insert(std::make_pair("rho_trans",&p_rho_trans));
      
      return;
    }

    /** \brief Remove model parameters */
    virtual void remove_params(o2scl::cli &cl) {
      size_t i=cl.par_list.erase("rho_trans");
      if (i!=1) {
	O2SCL_ERR("Failed to erase parameter 'rho_trans'.",o2scl::exc_esanity);
      }
      return;
    }
    
    /** \brief Copy model parameters */
    virtual void copy_params(model &m) {
      // Dynamic casts to references throw exceptions when they fail
      // while dynamic casts to pointers return null pointers when
      // they fail.
      qmc_threep_ligo &tp=dynamic_cast<qmc_threep_ligo &>(m);
      rho_trans=tp.rho_trans;
      return;
    }

    virtual void compute_eos(const ubvector &e, int &success,
			     std::ofstream &scr_out, model_data &dat);
    
  };
  
}

#endif
