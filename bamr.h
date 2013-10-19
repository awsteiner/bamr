/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2013, Andrew W. Steiner
  
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

#include <o2scl/cold_nstar.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/hist.h>
#include <o2scl/shared_ptr.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "misc.h"
#include "entry.h"
#include "models.h"

namespace bamr {

  /** \brief Statistical analysis of EOS from M and R constraints

      \todo It's not clear if successive calls of the mcmc command
      really work. For now, one may have ensure the program exits
      after each mcmc() run. 
      \todo More testing
      \todo Compute baryonic masses
      \todo Better documentation
      \todo Help with plots
      \todo Currently warm_up is changed to false only during 
      the first accepted step after <tt>iteration > n_warm_up</tt>,
      not immediately after <tt>iteration > n_warm_up</tt>.
      \future Allow the user to control how often the output 
      file is updated (currently every 10 successful MH steps or 
      when a block is finished)
      \future Allow non-tabulated data specified as a function
      \future Get rid of the 'best' file
  */
  class bamr_class {
    
  protected:

    /// \name Parameter objects for the 'set' command
    //@{
    o2scl::cli::parameter_double p_max_time;
    o2scl::cli::parameter_double p_min_max_mass;
    o2scl::cli::parameter_double p_step_fac;
    o2scl::cli::parameter_double p_exit_mass;
    o2scl::cli::parameter_double p_input_dist_thresh;
    o2scl::cli::parameter_double p_min_mass;
    o2scl::cli::parameter_int p_warm_up;
    o2scl::cli::parameter_int p_grid_size;
    o2scl::cli::parameter_int p_user_seed;
    o2scl::cli::parameter_int p_max_iters;
    o2scl::cli::parameter_bool p_norm_max;
    o2scl::cli::parameter_bool p_debug_star;
    o2scl::cli::parameter_bool p_debug_load;
    o2scl::cli::parameter_bool p_debug_eos;
    o2scl::cli::parameter_bool p_output_next;
    o2scl::cli::parameter_bool p_baryon_density;
    o2scl::cli::parameter_bool p_use_crust;
    o2scl::cli::parameter_string p_first_point_file;
    o2scl::cli::parameter_double p_nb_low;
    o2scl::cli::parameter_double p_nb_high;
    o2scl::cli::parameter_double p_e_low;
    o2scl::cli::parameter_double p_e_high;
    o2scl::cli::parameter_double p_m_low;
    o2scl::cli::parameter_double p_m_high;
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
    /** \brief If true, output debug information about the input data 
	files (default false)
    */
    bool debug_load;

    /** \brief If true, normalize the data distributions so that the
	max is one, otherwise, normalize so that the integral is one
    */
    bool norm_max;

    /// Maximum number of iterations (default 0)
    int max_iters;

    /// If true, use the default crust (default true)
    bool use_crust;

    /// If true, output stellar properties for debugging
    bool debug_star;
    
    /// If true, output equation of state for debugging 
    bool debug_eos;

    /// If true, output next point
    bool output_next;

    /// If true, compute the baryon density
    bool baryon_density;

    /// MCMC stepsize factor (default 15)
    double step_fac;

    /// The lower threshold for the input distributions
    double input_dist_thresh;

    /// An upper mass threshold
    double exit_mass;

    /// Minimum mass allowed for any of the individual neutron stars
    double min_mass;
  
    /// Minimum allowed maximum mass
    double min_max_mass;

    /// Number of warm up steps (successful steps not iterations)
    int n_warm_up;

    /** \brief Time in seconds (3600 seconds is one hour, default is
	86400 seconds or 1 day)
    */
    double max_time;

    /// If non-zero, use as the seed for the random number generator
    int user_seed;

    /** \brief If true, output more detailed information about the 
	best point
    */
    bool best_detail;
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

    /// \name MPI properties
    //@{
    /// The MPI processor rank
    int mpi_rank;

    /// The MPI number of processors
    int mpi_nprocs;

    /// The MPI starting time
    double mpi_start_time;
    //@}
    
    /// \name Input data
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
    //@}

    /// \name Parameter limits
    //@{
    entry low, high;
    //@}

    /// The first point in the parameter space
    ubvector first_point;

    /// The initial set of neutron star masses
    std::vector<double> first_mass;

    /// \name Other variables
    //@{
    /// The file containing the initial point
    std::string first_point_file;

    /// Main data table for Markov chain
    o2scl::table_units<> tc;
    
    /// The number of Metropolis steps which succeeded
    size_t mh_success;

    /// The number of Metropolis steps which failed
    size_t mh_failure;

    /// Total number of mcmc iterations
    size_t mcmc_iterations;

    /// Store the full Markov chain
    std::vector<double> markov_chain;

#ifdef BAMR_READLINE
    /// Command-line interface
    o2scl::cli_readline cl;
#else
    /// Command-line interface
    o2scl::cli cl;
#endif

    /// If true, then parameter names have been written to the HDF file
    bool first_file_update;

    /// Number of bins for all histograms (default 100)
    int grid_size;

    /// Number of Markov chain segments
    size_t n_chains;

    /// Number of chains
    size_t chain_size;
  
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

    /// The model for the EOS
    model *modp;

    /// The second copy of the model for the EOS
    model *modp2;

    /// Random number generator
    o2scl::rng_gsl gr;
  
    /// The screen output file
    std::ofstream scr_out;
    //@}

    /// \name Main functions called from the command-line interface
    //@{
    /** \brief Set the model for the EOS to use
     */
    virtual int set_model(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the first point in the parameter space
     */
    virtual int set_first_point(std::vector<std::string> &sv, bool itive_com);

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
    virtual void select_mass(entry &e_current, entry &e_next, double mmax,
			     bool &bad_step);

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
			   bool warm_up, int iteration);

    /** \brief Initialize the expectation value objects

	This function is called by mcmc().
    */
    virtual void init_grids_table(entry &low, entry &high);

    /** \brief Further preparations of the EOS before calling the TOV 
	solver
	
	This function is empty by default.
    */
    virtual void prepare_eos(entry &e, model &modref, o2scl::tov_solve *tsr, 
			     bool &success);

    /** \brief Add a measurement
     */
    virtual void add_measurement
      (entry &e, o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_eos,
       o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_mvsr,
       double weight, bool new_meas, size_t n_meas, ubvector &weights);

    /** \brief Fill vector in <tt>line</tt> with data from the
	current Monte Carlo point
     */
    virtual void fill_line
      (entry &e, o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_eos,
       o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_mvsr,
       double weight, bool new_meas, size_t n_meas, ubvector &weights,
       std::vector<double> &line);
    
    /** \brief Write initial data to HDF file
     */
    virtual void first_update(o2scl_hdf::hdf_file &hf, model &modp);
    
    /** \brief Write histogram data to files with prefix \c fname
     */
    virtual void update_files(std::string fname_prefix, 
			      model &modp, entry &e_current);
    
    /** \brief Set up the 'cli' object

	This function just adds the four commands and the 'set' parameters
    */
    virtual void setup_cli();
    
    /** \brief Load input probability distributions (called by mcmc())
     */
    virtual void load_mc();
  
    /** \brief Compute P(M|D) from parameters in entry \c e
	
	Called by mcmc().
    */
    virtual double compute_weight(entry &e, model &modref, 
				  o2scl::tov_solve *ts, bool &success,
				  ubvector &wgts, bool warm_up);
    
    /** \brief Tabulate EOS and then use in cold_nstar

	Called by compute_weight().
    */
    virtual void compute_star(entry &e, model &modref, o2scl::tov_solve *ts, 
			     bool &success);

    /// Output the "best" EOS obtained so far (called by mcmc())
    virtual void output_best
      (std::string fname_prefix, entry &e_best, double w_best,
       o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_eos,
       o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_mvsr,
       ubvector &wgts);
    //@}

    /// \name TOV solver objects
    //@{
    /// EOS interpolation object for TOV solver
    o2scl::tov_interp_eos teos;

    /// Pointer to a TOV solver
    o2scl::tov_solve *ts;
  
    /// Second pointer to a TOV solver
    o2scl::tov_solve *ts2;

    /// Default TOV solver
    o2scl::tov_solve def_ts;
  
    /// Second default TOV solver
    o2scl::tov_solve def_ts2;
    //@}

    /** \brief The arguments sent to the command-line
     */
    std::vector<std::string> run_args;

  public:

    bamr_class();

    virtual ~bamr_class();

    /// Main wrapper for parsing command-line arguments
    virtual void run(int argc, char *argv[]);

  };

}

#endif
