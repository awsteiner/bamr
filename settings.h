/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2020, Andrew W. Steiner
  
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
/** \file settings.h
    \brief Definition of \ref settings
*/
#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "nstar_cold2.h"

namespace bamr {
  
  /** \brief Settings object
   */
  class settings {

  public:

    settings() {
      grid_size=100;
      debug_load=false;
      min_max_mass=2.0;
      min_mass=1.0;
      debug_star=false;
      debug_eos=false;
      baryon_density=true;
      exit_mass=10.0;
      input_dist_thresh=0.0;
      use_crust=true;
      best_detail=false;
      inc_baryon_mass=false;
      norm_max=false;
      mvsr_pr_inc=1.1;
      nb_low=0.04;
      nb_high=1.24;
      e_low=0.3;
      e_high=10.0;
      m_low=0.2;
      m_high=3.0;
      verbose=0;
      in_m_min=1.0;
      in_m_max=3.0;
      in_r_min=5.0;
      in_r_max=18.0;
      addl_quants=false;
      mass_switch=0;
      compute_cthick=false;
      crust_from_L=false;
      mpi_load_debug=false;
      data_dir="data";
      apply_emu=false;
      //emu_train="";
      mmax_deriv = false;
    }
    
    /// \name Parameter objects for the 'set' command
    //@{
    o2scl::cli::parameter_double p_min_max_mass;
    o2scl::cli::parameter_double p_exit_mass;
    o2scl::cli::parameter_double p_input_dist_thresh;
    o2scl::cli::parameter_double p_min_mass;
    o2scl::cli::parameter_int p_grid_size;
    o2scl::cli::parameter_int p_mass_switch;
    o2scl::cli::parameter_int p_verbose;
    o2scl::cli::parameter_bool p_debug_star;
    o2scl::cli::parameter_bool p_debug_load;
    o2scl::cli::parameter_bool p_debug_eos;
    o2scl::cli::parameter_bool p_baryon_density;
    o2scl::cli::parameter_bool p_use_crust;
    o2scl::cli::parameter_bool p_inc_baryon_mass;
    o2scl::cli::parameter_bool p_norm_max;
    o2scl::cli::parameter_bool p_compute_cthick;
    o2scl::cli::parameter_bool p_addl_quants;
    o2scl::cli::parameter_bool p_crust_from_L;
    o2scl::cli::parameter_bool p_mpi_load_debug;
    o2scl::cli::parameter_bool p_apply_intsc;
    o2scl::cli::parameter_bool p_cached_intsc;
    o2scl::cli::parameter_double p_nb_low;
    o2scl::cli::parameter_double p_nb_high;
    o2scl::cli::parameter_double p_e_low;
    o2scl::cli::parameter_double p_e_high;
    o2scl::cli::parameter_double p_m_low;
    o2scl::cli::parameter_double p_m_high;
    o2scl::cli::parameter_double p_mvsr_pr_inc;
    o2scl::cli::parameter_bool p_prior_q;
    o2scl::cli::parameter_bool p_prior_eta;
    o2scl::cli::parameter_bool p_prior_delm;
    o2scl::cli::parameter_string p_data_dir;
    o2scl::cli::parameter_bool p_apply_emu;
    o2scl::cli::parameter_bool p_couple_threads;
    o2scl::cli::parameter_string p_emu_train;
    o2scl::cli::parameter_bool p_mmax_deriv;
    //@}

    /// Verbosity parameter
    int verbose;
    
    /** \name Limits on mass and radius from source data files

	These are automatically computed in load_mc() as the
	smallest rectangle in the \f$ (M,R) \f$ plane which
	encloses all of the user-specified source data
    */
    //@{
    double in_m_min;
    double in_m_max;
    double in_r_min;
    double in_r_max;
    //@}

    /// \name Other parameters accessed by 'set' and 'get'
    //@{
    /// Number of bins for all histograms (default 100)
    int grid_size;

    /// Pressure increment for the M vs. R curve (default 1.1)
    double mvsr_pr_inc;

    /** \brief If true, normalize the data distributions so that the
	max is one, otherwise, normalize so that the integral is one
	(default false)
    */
    bool norm_max;

    /// If true, use the default crust (default true)
    bool use_crust;

    /** \brief If true, output debug information about the input data 
	files (default false)
    */
    bool debug_load;

    /// If true, output stellar properties for debugging (default false)
    bool debug_star;
    
    /// If true, output equation of state for debugging (default false)
    bool debug_eos;

    /// If true, compute the baryon density (default true)
    bool baryon_density;

    /// The lower threshold for the input distributions (default 0.0)
    double input_dist_thresh;

    /// The upper mass threshold (default 10.0)
    double exit_mass;

    /// If true, debug MPI file loading 
    bool mpi_load_debug;
    
    /** \brief Minimum mass allowed for any of the individual neutron
	stars (default 1.0)
    */
    double min_mass;
  
    /// Minimum allowed maximum mass (default 2.0)
    double min_max_mass;

    /** \brief If true, output more detailed information about the 
	best point (default false)
    */
    bool best_detail;
    
    /** \brief If true, output information about the baryon mass
	as well as the gravitational mass (default false)
    */
    bool inc_baryon_mass;
    /// If true, compute crust thicknesses (default false)
    bool compute_cthick;
    /// If true (default false)
    bool addl_quants;
    /** \brief If true, compute a crust consistent with current 
	value of L

	Only works if \ref use_crust, \ref baryon_density, and
	\ref compute_cthick are true and the model 
	provides S and L.
    */
    bool crust_from_L;
    int mass_switch;

    /** \brief If true, include emulator from sklearn
     */
    bool apply_emu;

    /** \brief Desc
     */
    std::string emu_train;
    
    bool couple_threads;
    //@}

    /** \brief If true, use the eta prior
     */
    bool prior_eta;

    /** \brief If true, use the q prior
     */
    bool prior_q;

    /** \brief If true, use the delta_m prior
     */
    bool prior_delm;

    /** \name Histogram limits
     */
    //@{
    double nb_low;
    double nb_high;
    double e_low;
    double e_high;
    double m_low;
    double m_high;
    //@}

    /** \brief If true, apply an intrinic scattering correction
     */
    bool apply_intsc;

    /** \brief Desc
     */
    bool cached_intsc;

    /** \brief Desc
     */
    std::string data_dir;
    
    /** \brief Add parameters to the \ref o2scl::cli object
     */
    void setup_cli(o2scl::cli &cl) {
      
      // ---------------------------------------
      // Set parameters
      
      p_grid_size.i=&grid_size;
      p_grid_size.help="Grid size (default 100).";
      cl.par_list.insert(std::make_pair("grid_size",&p_grid_size));
      
      p_min_max_mass.d=&min_max_mass;
      p_min_max_mass.help=((std::string)"Minimum maximum mass ")
	+"(in solar masses, default 2.0).";
      cl.par_list.insert(std::make_pair("min_max_mass",&p_min_max_mass));
      
      p_min_mass.d=&min_mass;
      p_min_mass.help=
	((std::string)"Minimum possible mass for any of individual ")+
	"neutron stars in solar masses. The default is 1.0 solar masses.";
      cl.par_list.insert(std::make_pair("min_mass",&p_min_mass));
      
      p_exit_mass.d=&exit_mass;
      p_exit_mass.help=((std::string)"Upper limit on maximum mass ")+
	"(default 10.0). When the maximum mass is larger than this value, "+
	"the current point in the parameter space is output to 'cout' and "+
	"execution is aborted. This is sometimes useful in debugging the "+
	"initial guess.";
      cl.par_list.insert(std::make_pair("exit_mass",&p_exit_mass));
      
      p_input_dist_thresh.d=&input_dist_thresh;
      p_input_dist_thresh.help=((std::string)"Input distribution ")+
	"threshold. This is the artificial lower limit for the "+
	"(renormalized) probability of a (R,M) pair as reported by the "+
	"data file. If the weight is smaller than or equal to this value, "+
	"an exception is thrown. Changing this value is sometimes "+
	"useful to gracefully avoid zero probabilities in the input "+
	"data files. The default is 0.";
      cl.par_list.insert(std::make_pair("input_dist_thresh",
					&p_input_dist_thresh));
      
      p_debug_star.b=&debug_star;
      p_debug_star.help=((std::string)"If true, output stellar properties ")+
	"to file with suffix '_scr' at each point (default false).";
      cl.par_list.insert(std::make_pair("debug_star",&p_debug_star));
      
      p_norm_max.b=&norm_max;
      p_norm_max.help=((std::string)"If true, normalize by max probability ")+
	"or if false, normalize by total integral (default false). If "+
	"the Bayes factor is to be computed, this should be false.";
      cl.par_list.insert(std::make_pair("norm_max",&p_norm_max));
      
      p_debug_load.b=&debug_load;
      p_debug_load.help=((std::string)"If true, output info on loaded data ")+
	"(default false).";
      cl.par_list.insert(std::make_pair("debug_load",&p_debug_load));
      
      p_debug_eos.b=&debug_eos;
      p_debug_eos.help=((std::string)"If true, output initial equation ")+
	"of state to file 'debug_eos.o2' and abort (default false).";
      cl.par_list.insert(std::make_pair("debug_eos",&p_debug_eos));
      
      p_baryon_density.b=&baryon_density;
      p_baryon_density.help=((std::string)"If true, compute baryon density ")+
	"and associated profiles (default true).";
      cl.par_list.insert(std::make_pair("baryon_density",
					&p_baryon_density));
      
      p_use_crust.b=&use_crust;
      p_use_crust.help=((std::string)"If true, use the default crust ")+
	"(default true).";
      cl.par_list.insert(std::make_pair("use_crust",&p_use_crust));
      
      p_inc_baryon_mass.b=&inc_baryon_mass;
      p_inc_baryon_mass.help=((std::string)"If true, compute the baryon ")+
	"mass (default false)";
      cl.par_list.insert(std::make_pair("inc_baryon_mass",
					&p_inc_baryon_mass));
      
      p_mvsr_pr_inc.d=&mvsr_pr_inc;
      p_mvsr_pr_inc.help=((std::string)"The multiplicative pressure ")+
	"increment for the TOV solver (default 1.1). This can be "+
	"set to a value closer to 1.0 to make the mass-radius curve "+
	"more accurate.";
      cl.par_list.insert(std::make_pair("mvsr_pr_inc",&p_mvsr_pr_inc));
      
      // --------------------------------------------------------
      
      p_nb_low.d=&nb_low;
      p_nb_low.help=
	"Smallest baryon density grid point in 1/fm^3 (default 0.04).";
      cl.par_list.insert(std::make_pair("nb_low",&p_nb_low));
      
      p_nb_high.d=&nb_high;
      p_nb_high.help=
	"Largest baryon density grid point in 1/fm^3 (default 1.24).";
      cl.par_list.insert(std::make_pair("nb_high",&p_nb_high));
      
      p_e_low.d=&e_low;
      p_e_low.help=
	"Smallest energy density grid point in 1/fm^4 (default 0.3).";
      cl.par_list.insert(std::make_pair("e_low",&p_e_low));
      
      p_e_high.d=&e_high;
      p_e_high.help=
	"Largest energy density grid point in 1/fm^4 (default 10.0).";
      cl.par_list.insert(std::make_pair("e_high",&p_e_high));
      
      p_m_low.d=&m_low;
      p_m_low.help="Smallest mass grid point in Msun (default 0.2).";
      cl.par_list.insert(std::make_pair("m_low",&p_m_low));
      
      p_m_high.d=&m_high;
      p_m_high.help="Largest mass grid point in Msun (default 3.0).";
      cl.par_list.insert(std::make_pair("m_high",&p_m_high));
      
      p_crust_from_L.b=&crust_from_L;
      p_crust_from_L.help=((std::string)"If true, compute the core-crust ")+
	"transition density from L and S (requires a model which "+
	"provides these quantities). This also requires baryon_density, "+
	"compute_cthick, and use_crust are true (default false).";
      cl.par_list.insert(std::make_pair("crust_from_L",&p_crust_from_L));
      
      p_compute_cthick.b=&compute_cthick;
      p_compute_cthick.help=((std::string)"If true, compute the ")+
	"thickness of the crust (default false). This also requires "+
	"that baryon_density and use_crust are true.";
      cl.par_list.insert(std::make_pair("compute_cthick",&p_compute_cthick));

      p_addl_quants.b=&addl_quants;
      p_addl_quants.help=((std::string)"If true, compute moment of ")+
	"inertia, binding energy, and tidal deformability (default "+
	"false). This requires that baryon_mass is true.";
      cl.par_list.insert(std::make_pair("addl_quants",&p_addl_quants));

      p_mpi_load_debug.b=&mpi_load_debug;
      p_mpi_load_debug.help="";
      cl.par_list.insert(std::make_pair("mpi_load_debug",&p_mpi_load_debug));

      p_mass_switch.i=&mass_switch;
      p_mass_switch.help="";
      cl.par_list.insert(std::make_pair("mass_switch",&p_mass_switch));

      p_verbose.i=&verbose;
      p_verbose.help=((std::string)"This controls verbose output ")+
	"not already managed by 'mcmc_verbose'. Currently, no output "+
	"is performed unless verbose>=2, in which case some progress "+
	"on evaluating each point is written to cout. This value "+
	"defaults to zero.";
      cl.par_list.insert(std::make_pair("verbose",&p_verbose));

      p_prior_q.b=&prior_q;
      p_prior_q.help="";
      cl.par_list.insert(std::make_pair("prior_q",&p_prior_q));

      p_prior_eta.b=&prior_eta;
      p_prior_eta.help="";
      cl.par_list.insert(std::make_pair("prior_eta",&p_prior_eta));

      p_prior_delm.b=&prior_delm;
      p_prior_delm.help="";
      cl.par_list.insert(std::make_pair("prior_delm",&p_prior_delm));

      p_apply_intsc.b=&apply_intsc;
      p_apply_intsc.help="help";
      cl.par_list.insert(std::make_pair("apply_intsc",&p_apply_intsc));
      
      p_cached_intsc.b=&cached_intsc;
      p_cached_intsc.help="help";
      cl.par_list.insert(std::make_pair("cached_intsc",&p_cached_intsc));
      
      p_data_dir.str=&data_dir;
      p_data_dir.help="help";
      cl.par_list.insert(std::make_pair("data_dir",&p_data_dir));
      
      p_couple_threads.b=&couple_threads;
      p_couple_threads.help="help";
      cl.par_list.insert(std::make_pair("couple_threads",&p_couple_threads));
      
      p_apply_emu.b=&apply_emu;
      p_apply_emu.help="Activate emulator";
      cl.par_list.insert(std::make_pair("apply_emu",&p_apply_emu));

      // AWS: section for mmax_deriv
      p_mmax_deriv.b=&mmax_deriv;
      p_mmax_deriv.help="help";
      cl.par_list.insert(std::make_pair("mmax_deriv",&p_mmax_deriv));

      p_emu_train.str=&emu_train;
      p_emu_train.help="help";
      cl.par_list.insert(std::make_pair("emu_train",&p_emu_train));
      
      return;
    }
    
  };

}

#endif
