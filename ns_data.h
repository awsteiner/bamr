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
/** \file ns_data.h
    \brief Definition of \ref ns_data
*/
#ifndef NS_DATA_H
#define NS_DATA_H

#include <iostream>

#ifndef BAMR_NO_MPI
#include <mpi.h>
#endif

#include <o2scl/nstar_cold.h>
#include <o2scl/eos_had_schematic.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/cli.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_io.h>
#include <o2scl/mcmc_para.h>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "nstar_cold2.h"
#include "settings.h"

namespace bamr {
  
  /** \brief Desc
   */
  class ns_data {

  public:

    ns_data() {
      n_sources=0;
    }      
    
    /// \name Input neutron star data
    //@{
    /// Input probability distributions
    std::vector<o2scl::table3d> source_tables;

    /// The names for each source
    std::vector<std::string> source_names;

    /// The names of the table in the data file
    std::vector<std::string> table_names;

    /// File names for each source
    std::vector<std::string> source_fnames;

    /// Slice names for each source
    std::vector<std::string> slice_names;

    /// The initial set of neutron star masses
    std::vector<double> init_mass_fracs;

    /** \brief The number of sources
     */
    size_t n_sources;

    /** \brief Add a data distribution to the list
     */
    virtual int add_data(std::vector<std::string> &sv, bool itive_com);

    /** \brief Load input probability distributions
     */
    virtual void load_mc(std::ofstream &scr_out, int mpi_nprocs,
			 int mpi_rank, settings &set);
    //@}

  };
  
}

#endif
