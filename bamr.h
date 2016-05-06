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
#include "entry.h"
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
  class bamr_class : public mcmc_class<model_data,model> {
    
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

    /// \name Grids
    //@{
    o2scl::uniform_grid<double> nb_grid;
    o2scl::uniform_grid<double> e_grid;
    o2scl::uniform_grid<double> m_grid;
    //@}

    /// \name Other parameters accessed by 'set' and 'get'
    //@{
    /// Pressure increment for the M vs. R curve (default 1.1)
    double mvsr_pr_inc;

    /** \brief If true, normalize the data distributions so that the
	max is one, otherwise, normalize so that the integral is one
	(default true)
    */
    bool norm_max;

    /// If true, use the default crust (default true)
    bool use_crust;

    /** \brief If true, output debug information about the input data 
	files (default false)
    */
    bool debug_load;

    /** \brief If true, output each line of the table as it's stored
	(default false)
    */
    bool debug_line;
    
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

    /** \brief Minimum mass allowed for any of the individual neutron
	stars (default 0.8)
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
    //@}

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

    /** \brief The number of sources
     */
    size_t nsources;

    /// The initial set of neutron star masses
    std::vector<double> first_mass;
    //@}

    /// \name Parameter limits
    //@{
    entry low, high;
    //@}

    /// Number of bins for all histograms (default 100)
    int grid_size;

    /// A string indicating which model is used, set in \ref set_model().
    std::string model_type;
    
    /// True if the model provides S and L
    bool has_esym;

    /// True if the model has an EOS
    bool has_eos;

    /// Number of parameters, set in \ref set_model()
    size_t nparams;
    
    /// Schwarzchild radius (set in constructor)
    double schwarz_km;
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
    /// The function which selects the next neutron star masses
    virtual void select_mass(entry &e_current, entry &e_next, double mmax);

    /// Setup column names and units for data table
    virtual void table_names_units(std::string &s, std::string &u);

    /** \brief Make any necessary preparations for the mcmc() function

	This is called by \ref mcmc(). If the return value is non-zero
	then it is assumed that the calculation fails and mcmc()
	returns.

	This function does nothing and returns zero by default.
    */
    virtual int mcmc_init();

    /** \brief Decide to accept or reject the step 
     */
    virtual bool make_step(double w_current, double w_next, bool debug,
			   bool warm_up, int iteration,
			   double q_current, double q_next, double z);

    /** \brief Initialize the expectation value objects

	This function is called by mcmc().
    */
    virtual void init_grids_table(entry &low, entry &high);

    /** \brief Further preparations of the EOS before calling the TOV 
	solver
	
	This function is empty by default.
    */
    virtual void prepare_eos(entry &e, eos_tov &dat, int &success);

    /** \brief Add a measurement
     */
    virtual void add_measurement
      (entry &e, std::shared_ptr<o2scl::table_units<> > tab_eos,
       std::shared_ptr<o2scl::table_units<> > tab_mvsr,
       double weight, bool new_meas, size_t n_meas, ubvector &weights);

    /** \brief Fill vector in <tt>line</tt> with data from the
	current Monte Carlo point
    */
    virtual void fill_line
      (entry &e, std::shared_ptr<o2scl::table_units<> > tab_eos,
       std::shared_ptr<o2scl::table_units<> > tab_mvsr,
       double weight, bool new_meas, size_t n_meas, ubvector &weights,
       std::vector<double> &line);
    
    /** \brief Write initial data to HDF file
     */
    virtual void first_update(o2scl_hdf::hdf_file &hf, model &modp);
    
    /** \brief Write histogram data to files with prefix \c fname
     */
    virtual void update_files(model &modp, entry &e_current);
    
    /** \brief Set up the 'cli' object

	This function just adds the four commands and the 'set' parameters
    */
    virtual void setup_cli();
    
    /** \brief Load input probability distributions (called by mcmc())
     */
    virtual void load_mc();
  
    /** \brief Compute \f$ P(M|D) \f$ from parameters in entry \c e
	
	Called by mcmc().
    */
    virtual double compute_weight(entry &e, eos_tov &dat, int &success,
				  ubvector &wgts, bool warm_up);

  public:
    
    /** \brief Tabulate EOS and then use in cold_nstar
	
	Called by compute_weight().

	\todo Temporarily made public for drdp project hack
    */
    virtual void compute_star(entry &e, eos_tov &dat, int &success);
    
  protected:
    
    /// EOS interpolation object for TOV solver
    o2scl::eos_tov_interp teos;
    
    /// Output the "best" EOS obtained so far (called by mcmc())
    virtual void output_best
      (entry &e_best, double w_best,
       std::shared_ptr<o2scl::table_units<> > tab_eos,
       std::shared_ptr<o2scl::table_units<> > tab_mvsr,
       ubvector &wgts);
    //@}

  public:

    bamr_class();

    virtual ~bamr_class();

    /// \name Return codes for each point
    //@{
    static const int ix_success=0;
    static const int ix_mr_outside=1;
    static const int ix_r_outside=2;
    static const int ix_zero_wgt=3;
    static const int ix_press_dec=4;
    static const int ix_nb_problem=5;
    static const int ix_nb_problem2=6;
    static const int ix_crust_unstable=7;
    static const int ix_mvsr_failed=8;
    static const int ix_tov_failure=9;
    static const int ix_small_max=10;
    static const int ix_tov_conv=11;
    static const int ix_mvsr_table=12;
    static const int ix_acausal=13;
    static const int ix_acausal_mr=14;
    static const int ix_param_mismatch=15;
    static const int ix_neg_pressure=16;
    static const int ix_no_eos_table=17;
    static const int ix_eos_solve_failed=18;
    //@}

  };

}

#endif
