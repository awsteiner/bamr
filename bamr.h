/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
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
/** \file bamr.h
    \brief Definition of main bamr class
*/
#ifndef BAMR_H
#define BAMR_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#ifndef NO_MPI
#include <mpi.h>
#endif

#include <o2scl/rng_gsl.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/cholesky.h>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "mcmc.h"
#include "nstar_cold2.h"
#include "models.h"

/** \brief Main namespace
    
    The bamr namespace which holds all classes and functions.
    
    This file is documented in bamr.h .
*/
namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;

  /** \brief Statistical analysis of EOS from M and R constraints

      \comment
      \note Right now the EOS is rejected if the pressure decreases
      with increasing density at any density, even if it happens at a
      density which is larger than the central density of the maximum
      mass star.
      1/6/16 - I originally thought this was a problem, but
      actually there is no problem here, as EOSs can always be
      fixed to ensure that pressures always increase. 
      \endcomment

      \todo It's not clear if successive calls of the mcmc command
      really work. For now, one may have ensure the program exits
      after each mcmc() run. 

      \comment
      \todo Fix issue of block_counter giving confusing output if
      there are too few MCMC points between each block. 
      2/2/16 - I think this is fixed now
      \endcomment

      \todo More testing

      \todo Better documentation

      \todo Help with plots

      \future Allow non-tabulated data specified as a function?

      \future The code internally stores two copies of the model
      objects, which makes things a bit easier to handle MC
      rejections, but also requires models to have this funny
      copy_params() function to copy model parameters between model
      objects. There's probably a better way to do this.
  */
  class bamr_class : public mcmc_namespace::mcmc_class<model_data,model> {

#ifdef O2SCL_NEVER_DEFINED
    
  protected:

    /// \name Parameter objects for the 'set' command
    //@{
    o2scl::cli::parameter_double p_min_max_mass;
    o2scl::cli::parameter_double p_exit_mass;
    o2scl::cli::parameter_double p_input_dist_thresh;
    o2scl::cli::parameter_double p_min_mass;
    o2scl::cli::parameter_int p_grid_size;
    o2scl::cli::parameter_bool p_debug_star;
    o2scl::cli::parameter_bool p_debug_line;
    o2scl::cli::parameter_bool p_debug_load;
    o2scl::cli::parameter_bool p_debug_eos;
    o2scl::cli::parameter_bool p_baryon_density;
    o2scl::cli::parameter_bool p_use_crust;
    o2scl::cli::parameter_bool p_inc_baryon_mass;
    o2scl::cli::parameter_bool p_norm_max;
    o2scl::cli::parameter_double p_nb_low;
    o2scl::cli::parameter_double p_nb_high;
    o2scl::cli::parameter_double p_e_low;
    o2scl::cli::parameter_double p_e_high;
    o2scl::cli::parameter_double p_m_low;
    o2scl::cli::parameter_double p_m_high;
    o2scl::cli::parameter_double p_mvsr_pr_inc;
    o2scl::cli::parameter_string p_prefix;
    //@}

    //@}

    /// A string indicating which model is used, set in \ref set_model().
    std::string model_type;
    
    /// Number of parameters, set in \ref set_model()
    size_t nparams;
    
    //@}

    /// \name Main functions called from the command-line interface
    //@{
    /** \brief Set the model for the EOS to use
     */
    virtual int set_model(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the first point in the parameter space
     */
    virtual int set_first_point(std::vector<std::string> &sv, bool itive_com);

    /** \brief Prepare for MCMC using Metropolis-Hastings
     */
    virtual int hastings(std::vector<std::string> &sv, bool itive_com);
			 
    /** \brief Add a data distribution to the list
     */
    virtual int add_data(std::vector<std::string> &sv, bool itive_com);

    /** \brief Run a MCMC simulation

	In order to save the code from copying results over if a step
	is rejected, the Metropolis steps alternate between two kinds.
	In the first kind of step, \ref modp is the previous point and
	\ref modp2 is the new point to be considered, and in the
	second kind of step, \ref modp2 is the old point and \ref modp
	is the new point to be considered. This two-step procedure is
	also why there are two TOV solvers, i.e. \ref ts and \ref ts2.
    */
    virtual int mcmc(std::vector<std::string> &sv, bool itive_com);
    //@}

    /// \name Internal functions 
    //@{
    /// Setup column names and units for data table
    virtual void table_names_units(std::string &s, std::string &u);

    /** \brief Make any necessary preparations for the mcmc() function

	This is called by \ref mcmc(). If the return value is non-zero
	then it is assumed that the calculation fails and mcmc()
	returns.

	This function does nothing and returns zero by default.
    */
    virtual int mcmc_init();

    /** \brief Initialize the expectation value objects

	This function is called by mcmc().
    */
    virtual void init_grids_table(ubvector &low, ubvector &high);

    /** \brief Further preparations of the EOS before calling the TOV 
	solver
	
	This function is empty by default.
    */
    virtual void prepare_eos(ubvector &pars, model_data &dat, int &success);

    /** \brief Add a measurement
     */
    void add_measurement(ubvector &pars, double weight, model_data &dat,
			 bool new_meas, size_t n_meas);

    /** \brief Fill vector in <tt>line</tt> with data from the
	current Monte Carlo point
    */
    virtual void fill_line(ubvector &pars, double weight, model_data &dat,
			   bool new_meas, size_t n_meas,
			   std::vector<double> &line);
    
    /** \brief Write initial data to HDF file
     */
    virtual void first_update(o2scl_hdf::hdf_file &hf, model &modp);
    
    /** \brief Write histogram data to files with prefix \c fname
     */
    virtual void update_files(ubvector &current);
    
    /** \brief Set up the 'cli' object

	This function just adds the four commands and the 'set' parameters
    */
    virtual void setup_cli();
    
    /** \brief Load input probability distributions (called by mcmc())
     */
    virtual void load_mc();
  
  protected:
    
    /// Output the "best" EOS obtained so far (called by mcmc())
    void output_best(ubvector &best, double w_best, model_data &dat);

  public:

    bamr_class();

    virtual ~bamr_class();

    /// Desc
    void run(int argc, char *argv[]) {
      mcmc_class::run(argc,argv);
      return;
    }

    
#endif
    
  };

}

#endif
