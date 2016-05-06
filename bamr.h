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

#include "nstar_cold2.h"
#include "entry.h"
#include "models.h"

/** \brief Main namespace
    
    The bamr namespace which holds all classes and functions.
    
    This file is documented in bamr.h .
*/
namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;

  class mcmc_functions {
    virtual void compute_point(ubvector &pars);
    virtual void first_point(ubvector &pars)=0;
    virtual void low_limits(ubvector &pars)=0;
    virtual void high_limits(ubvector &pars)=0;
  };
  
  /** \brief Desc
   */
  template<class data_t> class mcmc_class {
    
  public:

    /// Desc
    size_t nparams;

    /// Desc
    std::vector<std::string> param_names;

    /// Desc
    bool use_smove;
  
    /// Desc
    size_t nwalk;

    /// \name Member data for the Metropolis-Hastings step
    //@{
    /// A Gaussian probability distribution
    o2scl::prob_dens_gaussian pdg;
    
    /// If true, then use Metropolis-Hastings with a multivariate Gaussian
    int hg_mode;
    
    /// The Cholesky decomposition of the covariance matrix
    ubmatrix hg_chol;
    
    /// The inverse of the covariance matrix
    ubmatrix hg_covar_inv;
    
    /// The normalization factor
    double hg_norm;
    
    /// The location of the peak
    ubvector hg_best;
    //@}

    /** \brief Error handler for each thread
     */
    o2scl::err_hnd_cpp error_handler;

    /** \brief Prefix for output filenames
     */
    std::string prefix;
    
    /// \name Parameter objects for the 'set' command
    //@{
    o2scl::cli::parameter_double p_max_time;
    o2scl::cli::parameter_double p_step_fac;
    o2scl::cli::parameter_int p_n_warm_up;
    o2scl::cli::parameter_int p_grid_size;
    o2scl::cli::parameter_int p_user_seed;
    o2scl::cli::parameter_int p_max_iters;
    o2scl::cli::parameter_int p_file_update_iters;
    o2scl::cli::parameter_int p_max_chain_size;
    o2scl::cli::parameter_bool p_output_next;
    o2scl::cli::parameter_bool p_use_smove;
    o2scl::cli::parameter_string p_prefix;
    //@}

    /** \brief The number of MCMC successes between file updates
	(default 40)
    */
    int file_update_iters;

    /** \brief Maximum size of Markov chain (default 10000)
     */
    int max_chain_size;
    
    /// Maximum number of iterations (default 0)
    int max_iters;

    /// MCMC stepsize factor (default 15.0)
    double step_fac;

    /** \brief Number of warm up steps (successful steps not
	iterations) (default 0)
	
	\note Not to be confused with <tt>warm_up</tt>, which is 
	a boolean local variable in some functions not an int.
    */
    int n_warm_up;

    /** \brief Time in seconds (default is 86400 seconds or 1 day)
     */
    double max_time;

    /** \brief If non-zero, use as the seed for the random number 
	generator (default 0)
    */
    int user_seed;

    /// \name MPI properties
    //@{
    /// The MPI processor rank
    int mpi_rank;

    /// The MPI number of processors
    int mpi_nprocs;

    /// The MPI starting time
    double mpi_start_time;
    //@}
    
    /// The number of Metropolis steps which succeeded
    size_t mh_success;

    /// The number of Metropolis steps which failed
    size_t mh_failure;

    /// Total number of mcmc iterations
    size_t mcmc_iterations;

    /// Number of Markov chain segments
    size_t n_chains;

    /// Number of chains
    size_t chain_size;
  
#ifdef BAMR_READLINE
    /// Command-line interface
    o2scl::cli_readline cl;
#else
    /// Command-line interface
    o2scl::cli cl;
#endif

    /// Main data table for Markov chain
    o2scl::table_units<> tc;
    
    /// Random number generator
    o2scl::rng_gsl gr;
  
    /// The screen output file
    std::ofstream scr_out;

    /// If true, scr_out has been opened
    bool file_opened;

    /// If true, output next point (default true)
    bool output_next;

    /** \brief The arguments sent to the command-line
     */
    std::vector<std::string> cl_args;

    /// Vector of data objects
    std::vector<data_t> data_arr;
  
    /** \brief Set up the 'cli' object
      
	This function just adds the four commands and the 'set' parameters
    */
    virtual void setup_cli() {

      // ---------------------------------------
      // Set parameters
    
      p_max_time.d=&max_time;
      p_max_time.help=
	"Maximum run time in seconds (default 86400 sec or 1 day).";
      cl.par_list.insert(std::make_pair("max_time",&p_max_time));
    
      p_step_fac.d=&step_fac;
      p_step_fac.help=
	((std::string)"MCMC step factor. The step size for each ")+
	"variable is taken as the difference between the high and low "+
	"limits divided by this factor (default 15.0). This factor can "+
	"be increased if the acceptance rate is too small, but care must "+
	"be taken, e.g. if the conditional probability is multimodal. If "+
	"this step size is smaller than 1.0, it is reset to 1.0 .";
      cl.par_list.insert(std::make_pair("step_fac",&p_step_fac));

      p_n_warm_up.i=&n_warm_up;
      p_n_warm_up.help=((std::string)"Minimum number of warm up iterations ")+
	"(default 0).";
      cl.par_list.insert(std::make_pair("n_warm_up",&p_n_warm_up));

      p_file_update_iters.i=&file_update_iters;
      p_file_update_iters.help=
	((std::string)"Number of MCMC successes between ")+
	"file upates (default 10, minimum value 1).";
      cl.par_list.insert(std::make_pair("file_update_iters",
					&p_file_update_iters));

      p_user_seed.i=&user_seed;
      p_user_seed.help=((std::string)"Seed for multiplier for random number ")+
	"generator. If zero is given (the default), then mcmc() uses "+
	"time(0) to generate a random seed.";
      cl.par_list.insert(std::make_pair("user_seed",&p_user_seed));

      p_max_iters.i=&max_iters;
      p_max_iters.help=((std::string)"If non-zero, limit the number of ")+
	"iterations to be less than the specified number (default zero).";
      cl.par_list.insert(std::make_pair("max_iters",&p_max_iters));

      p_output_next.b=&output_next;
      p_output_next.help=((std::string)"If true, output next point ")+
	"to the '_scr' file before calling TOV solver (default true).";
      cl.par_list.insert(std::make_pair("output_next",&p_output_next));

      p_use_smove.b=&use_smove;
      p_use_smove.help="";
      cl.par_list.insert(std::make_pair("use_smove",&p_use_smove));

      p_max_chain_size.i=&max_chain_size;
      p_max_chain_size.help=
	((std::string)"Maximum Markov chain size (default ")+
	"10000).";
      cl.par_list.insert(std::make_pair("max_chain_size",&p_max_chain_size));

      p_prefix.str=&prefix;
      p_prefix.help="Output file prefix (default 'bamr').";
      cl.par_list.insert(std::make_pair("prefix",&p_prefix));

      return;
    }    
    
    /// Desc
    virtual int mcmc(std::vector<std::string> &sv, bool itive_com) {
      
      // Make sure that first_update() is called when necessary
      first_file_update=false;
      
#ifndef NO_MPI
      mpi_start_time=MPI_Wtime();
#else
      mpi_start_time=time(0);
#endif
      
      if (file_opened==false) {
	// Open main output file
	scr_out.open((prefix+"_"+std::to_string(mpi_rank)+"_scr").c_str());
	scr_out.setf(ios::scientific);
	file_opened=true;
	scr_out << "Opened main file in command 'mcmc'." << endl;
      }
      
      // Fix file_update_iters if necessary
      if (file_update_iters<1) {
	scr_out << "Parameter 'file_update_iters' less than 1. Set equal to 1."
		<< endl;
	file_update_iters=1;
      }
      
      if (max_chain_size<1) {
	O2SCL_ERR("Parameter 'max_chain_size' must be larger than 1.",
		  exc_einval);
      }
      
      // Fix step_fac if it's too small
      if (step_fac<1.0) {
	step_fac=1.0;
	scr_out << "Fixed 'step_fac' to 1.0." << endl;
      }
      
      // Run init() function (have to make sure to do this after opening
      // scr_out). 
      int iret=mcmc_init();
      if (iret!=0) return iret;
      
      // Set RNG seed
      unsigned long int seed=time(0);
      if (user_seed!=0) {
	seed=user_seed;
      }
      seed*=(mpi_rank+1);
      gr.set_seed(seed);
      pdg.set_seed(seed);
      scr_out << "Using seed " << seed 
	      << " for processor " << mpi_rank+1 << "/" 
	      << mpi_nprocs << "." << endl;
      scr_out.precision(12);
      scr_out << " Start time: " << mpi_start_time << endl;
      scr_out.precision(6);
      
      // First MC point
      std::vector<ubvector> current(1);
      std::vector<double> w_current(1);
      std::vector<bool> step_flags(1);
      data_arr.resize(2);
      if (use_smove) {
	current.resize(nwalk);
	w_current.resize(nwalk);
	step_flags.resize(nwalk);
	data_arr.resize(2*nwalk);
	for(size_t i=0;i<nwalk;i++) {
	  step_flags[i]=false;
	}
      }

#ifndef NO_MPI
  int buffer=0, tag=0;
  if (mpi_nprocs>1 && mpi_rank>0) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
  }
#endif

  if (first_point_file.length()>0) {
  
    if (first_point_type==fp_last) {

      // Read file 
      scr_out << "Reading last point from file '" << first_point_file
	      << "'." << endl;
      hdf_file hf;
      hf.open(first_point_file);
      
      // Read table
      size_t file_n_chains;
      hf.get_szt("n_chains",file_n_chains);
      std::string chain_name=std::string("markov_chain")+
	o2scl::szttos(file_n_chains-1);
      table_units<> file_tab;
      hdf_input(hf,file_tab,chain_name);
      size_t last_line=file_tab.get_nlines()-1;
      
      // Get parameters
      for(size_t i=0;i<nparams;i++) {
	string pname=((string)"param_")+param_name(i);
	current[0][i]=file_tab.get(pname,last_line);
	scr_out << "Parameter named " << param_name(i) << " " 
		<< current[0][i] << endl;
      }
      
      // Finish up
      scr_out << endl;
      hf.close();

    } else if (first_point_type==fp_best) {

      ubvector best_point;
      hdf_file hf;
      hf.open(first_point_file);
      hf.getd_vec("best_point",best_point);
      hf.close();
      scr_out << "Reading best point from file '" << first_point_file
	      << "'." << endl;
      for(size_t i=0;i<nparams;i++) {
	current[0][i]=best_point[i];
	scr_out << "Parameter " << i << " : " << current[0][i] << endl;
      }
      scr_out << "Best weight: " << best_point[nparams] << endl;
      scr_out << endl;

    } else {

      // Read file 
      scr_out << "Reading " << first_point_type << "th point from file '" 
	      << first_point_file
	      << "'." << endl;
      hdf_file hf;
      hf.open(first_point_file);
      
      // Read table
      size_t file_n_chains, row=first_point_type;
      hf.get_szt("n_chains",file_n_chains);
      
      table_units<> file_tab;
      for(size_t k=0;k<file_n_chains;k++) {
	std::string chain_name=std::string("markov_chain")+o2scl::szttos(k);
	hdf_input(hf,file_tab,chain_name);
	if (file_tab.get_nlines()>row) {
	  k=file_n_chains;
	} else {
	  row-=file_tab.get_nlines();
	}
      }
      if (row>=file_tab.get_nlines()) {
	scr_out << "Couldn't find point " << first_point_type 
		<< " in file. Using last point." << endl;
	row=file_tab.get_nlines()-1;
      }
      
      // Get parameters
      for(size_t i=0;i<nparams;i++) {
	string pname=((string)"param_")+modp->param_name(i);
	current[0][i]=file_tab.get(pname,row);
	scr_out << "Parameter named " << param_name(i) << " " 
		<< current[0][i] << endl;
      }
      
      // Finish up
      scr_out << endl;
      hf.close();
    }

  } else if (first_point.size()>0) {
    
    scr_out << "First point from command-line." << endl;
    for(size_t i=0;i<nparams;i++) {
      current[0][i]=first_point[i];
      scr_out << current[0][i] << endl;
    }
    scr_out << endl;

  } else {
    
    scr_out << "First point from default." << endl;
    data_arr[0].modp->first_point(e_current);

  }

#ifndef NO_MPI
  if (mpi_nprocs>1 && mpi_rank<mpi_nprocs-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,tag,MPI_COMM_WORLD);
  }
#endif

  scr_out << "First point: " << e_current << endl;

  // Determine initial masses

  for(size_t i=0;i<nsources;i++) {
    e_current.mass[i]=first_mass[i];
    e_current.rad[i]=0.0;
  }

  // Set lower and upper bounds for parameters
  low.allocate(nparams,nsources);
  high.allocate(nparams,nsources);
  data_arr[0].modp->low_limits(low);
  data_arr[0].modp->high_limits(high);

  n_chains=0;

  // Vector containing weights
  ubvector wgts(nsources);

  // Entry objects (Must be after read_input() since nsources is set
  // in that function.)
  entry e_next(nparams,nsources);
  entry e_best(nparams,nsources);

  // Load data
  load_mc();

  // Set lower and upper bounds for masses and radii from load_mc()
  for(size_t i=0;i<nsources;i++) {
    low.mass[i]=in_m_min;
    high.mass[i]=in_m_max;
    low.rad[i]=in_r_min;
    high.rad[i]=in_r_max;
  }

  // Prepare statistics
  init_grids_table(low,high);

  // Weights for each entry
  double w_current=0.0, w_next=0.0, w_best=0.0;

  // Warm-up flag, not to be confused with 'n_warm_up', i.e. the
  // number of warm_up iterations.
  bool warm_up=true;
  if (n_warm_up==0) warm_up=false;

  // Keep track of successful and failed MH moves
  mh_success=0;
  mh_failure=0;

  // ---------------------------------------------------
  // Compute initial point and initial weights

  int suc;
  double q_current=0.0, q_next=0.0;

  if (use_smove) {
    // Stretch-move steps

    size_t ij_best=0;

    // Initialize each walker in turn
    for(size_t ij=0;ij<nwalk;ij++) {

      size_t init_iters=0;
      bool done=false;

      while (!done) {

	// Begin with the intial point
	entry e_first(nparams,nsources);
	data_arr[ij].modp->first_point(e_first);

	// Make a perturbation from the initial point
	for(size_t ik=0;ik<nparams;ik++) {
	  do {
	    e_curr_arr[ij].params[ik]=e_first.params[ik]+
	      (gr.random()*2.0-1.0)*(high.params[ik]-low.params[ik])/100.0;
	  } while (e_curr_arr[ij].params[ik]>=high.params[ik] ||
		   e_curr_arr[ij].params[ik]<=low.params[ik]);
	}
	for(size_t ik=0;ik<nsources;ik++) {
	  do {
	    e_curr_arr[ij].mass[ik]=first_mass[ik]+
	      (gr.random()*2.0-1.0)*(high.mass[ik]-low.mass[ik])/100.0;
	  } while (e_curr_arr[ij].mass[ik]>=high.mass[ik] ||
		   e_curr_arr[ij].mass[ik]<=low.mass[ik]);
	  e_curr_arr[ij].rad[ik]=0.0;
	}
	
	// Compute the weight
	w_curr_arr[ij]=compute_weight(e_curr_arr[ij],data_arr[ij],
				      suc,wgts,warm_up);
	scr_out << "SM Init: " << ij << " "
		<< e_curr_arr[ij] << " " << w_curr_arr[ij] << endl;

	// Keep track of the best point and the best index
	if (ij==0) {
	  w_best=w_curr_arr[0];
	} else if (w_curr_arr[ij]>w_best) {
	  ij_best=ij;
	  w_best=w_curr_arr[ij];
	}

	// Increment iteration count
	init_iters++;

	// If we have a good point, stop the loop
	if (w_curr_arr[ij]>0.0) {
	  done=true;
	  ret_codes[suc]++;
	} else if (init_iters>1000) {
	  scr_out << "Failed to construct initial walkers." << endl;
	  return 1;
	}
      }

      // For the initial point for this walker, add it
      // to the result table
      {
	shared_ptr<table_units<> > tab_eos;
	shared_ptr<table_units<> > tab_mvsr;
	tab_eos=data_arr[ij].modp->cns.get_eos_results();
	tab_mvsr=data_arr[ij].ts.get_results();
	if (warm_up==false) {
	  // Add the initial point if there's no warm up
	  add_measurement(e_curr_arr[ij],tab_eos,tab_mvsr,w_curr_arr[ij],
			  true,mh_success,wgts);
	}
      }
    }

    // Output the best initial walker if necessary
    {
      e_best=e_curr_arr[ij_best];
      shared_ptr<table_units<> > tab_eos;
      shared_ptr<table_units<> > tab_mvsr;
      tab_eos=data_arr[ij_best].modp->cns.get_eos_results();
      tab_mvsr=data_arr[ij_best].ts.get_results();
      output_best(e_best,w_best,tab_eos,tab_mvsr,wgts);
    }

  } else {
    // Normal or Metropolis-Hastings steps

    // Compute weight for initial point
    w_current=compute_weight(e_current,data_arr[0],suc,wgts,warm_up);
    ret_codes[suc]++;
    scr_out << "Initial weight: " << w_current << endl;
    if (w_current<=0.0) {
      for(size_t i=0;i<nsources;i++) {
	scr_out << i << " " << wgts[i] << endl;
      }
      scr_out << "Initial weight zero. Aborting." << endl;
      exit(-1);
    }

    // Compute the initial Hastings proposal weight
    if (hg_mode>0) {
      q_current=approx_like(e_current);
    }

    // Add measurement to output table and output best 
    {
      shared_ptr<table_units<> > tab_eos;
      shared_ptr<table_units<> > tab_mvsr;
      tab_eos=data_arr[0].modp->cns.get_eos_results();
      tab_mvsr=data_arr[0].ts.get_results();
      if (warm_up==false) {
	// Add the initial point if there's no warm up
	add_measurement(e_current,tab_eos,tab_mvsr,w_current,
			true,mh_success,wgts);
      }
      e_best=e_current;
      w_best=w_current;
      output_best(e_current,w_current,tab_eos,tab_mvsr,wgts);
    }
    
  }

  // ---------------------------------------------------

  // Initialize radii of e_next to zero
  for(size_t k=0;k<nsources;k++) {
    e_next.rad[k]=0.0;
  }

  // Keep track of total number of different points in the parameter
  // space that are considered (some of them may not result in TOV
  // calls because, for example, the EOS was acausal).
  mcmc_iterations=0;

  // Switch between two different Metropolis steps
  bool first_half=true;

  // Main loop
  bool main_done=false;

  // The MCMC is arbitrarily broken up into 20 'blocks', making
  // it easier to keep track of progress and ensure file updates
  size_t block_counter=0;

  while (!main_done) {

    // Walker to move for smove
    size_t ik=0;
    double smove_z=0.0;
    
    // ---------------------------------------------------
    // Select next point
    
    if (use_smove) {

      // Choose walker to move
      ik=mcmc_iterations % nwalk;
      
      bool in_bounds;
      size_t step_iters=0;
      
      do {

	in_bounds=true;
	
	// Choose jth walker
	size_t ij;
	do {
	  ij=((size_t)(gr.random()*((double)nwalk)));
	} while (ij==ik || ij>=nwalk);
	
	// Select z 
	double p=gr.random();
	double a=step_fac;
	smove_z=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;
	
	// Create new trial point
	for(size_t i=0;i<nparams;i++) {
	  e_next.params[i]=e_curr_arr[ij].params[i]+smove_z*
	    (e_curr_arr[ik].params[i]-e_curr_arr[ij].params[i]);
	  if (e_next.params[i]>=high.params[i] ||
	      e_next.params[i]<=low.params[i]) {
	    in_bounds=false;
	  }
	}
	for(size_t i=0;i<nsources;i++) {
	  e_next.mass[i]=e_curr_arr[ij].mass[i]+smove_z*
	    (e_curr_arr[ik].mass[i]-e_curr_arr[ij].mass[i]);
	  if (e_next.mass[i]>=high.mass[i] ||
	      e_next.mass[i]<=low.mass[i]) {
	    in_bounds=false;
	  }
	}
	
	step_iters++;
	if (step_iters==1000) {
	  scr_out << "Failed to find suitable step at point 1." << endl;
	  cerr << "Failed to find suitable step at point 1." << endl;
	  return 2;
	}

      } while (in_bounds==false);

    } else if (hg_mode>0) {
      
      // Make a Metropolis-Hastings step based on previous data
      
      size_t nv=e_next.np+nsources;
      ubvector hg_temp(nv), hg_z(nv);
      
      bool out_of_range;
      int hg_it=0;
      
      do {
	
	for(size_t k=0;k<nv;k++) {
	  hg_z[k]=pdg.sample();
	}
	hg_temp=prod(hg_chol,hg_z);
	for(size_t k=0;k<nv;k++) {
	  if (k<e_next.np) {
	    e_next.params[k]=hg_best[k]+hg_temp[k];
	  } else {
	    e_next.mass[k-e_next.np]=hg_best[k]+hg_temp[k];
	  }
	}
	
	out_of_range=false;
	for(size_t k=0;k<e_next.np;k++) {
	  if (e_next.params[k]<low.params[k] ||
	      e_next.params[k]>high.params[k]) {
	    out_of_range=true;
	  }
	}
	for(size_t k=0;k<nsources;k++) {
	  if (e_next.mass[k]<low.mass[k] ||
	      e_next.mass[k]>high.mass[k]) {
	    out_of_range=true;
	  }
	}

	hg_it++;
	if (hg_it>1000) {
	  O2SCL_ERR("Sanity check in hg step.",exc_esanity);
	}

      } while (out_of_range==true);

      q_next=approx_like(e_next);

    } else {

      // Make a step, ensure that we're in bounds and that
      // the masses are not too large
      for(size_t k=0;k<e_next.np;k++) {
	
	e_next.params[k]=e_current.params[k]+(gr.random()*2.0-1.0)*
	  (high.params[k]-low.params[k])/step_fac;
	
	// If it's out of range, redo step near boundary
	if (e_next.params[k]<low.params[k]) {
	  e_next.params[k]=low.params[k]+gr.random()*
	    (high.params[k]-low.params[k])/step_fac;
	} else if (e_next.params[k]>high.params[k]) {
	  e_next.params[k]=high.params[k]-gr.random()*
	    (high.params[k]-low.params[k])/step_fac;
	}
	
	if (e_next.params[k]<low.params[k] || 
	    e_next.params[k]>high.params[k]) {
	  O2SCL_ERR("Sanity check in parameter step.",exc_esanity);
	}
      }
      
      if (nsources>0) {
	// Just use a large value (1.0e6) here since we don't yet
	// know the maximum mass
	select_mass(e_current,e_next,1.0e6);
      }
      
    }

    // End of select next point
    // ---------------------------------------------------

    // Output the next point
    if (output_next) {
      scr_out << "Iteration, next: " << mcmc_iterations << " " 
	      << e_next << endl;
    }
      
    // ---------------------------------------------------
    // Compute next weight

    if (use_smove) {
      if (step_flags[ik]==false) {
	w_next=compute_weight(e_next,data_arr[ik+nwalk],suc,wgts,warm_up);
      } else {
	w_next=compute_weight(e_next,data_arr[ik],suc,wgts,warm_up);
      }
    } else {
      if (first_half) {
	w_next=compute_weight(e_next,data_arr[1],suc,wgts,warm_up);
      } else {
	w_next=compute_weight(e_next,data_arr[0],suc,wgts,warm_up);
      }
    }
    ret_codes[suc]++;

    // ---------------------------------------------------
    
    // Test to ensure new point is good
    if (suc==ix_success) {

      // Test radii
      for(size_t i=0;i<nsources;i++) {
	if (e_next.rad[i]>high.rad[i] || e_next.rad[i]<low.rad[i]) {
	  scr_out << "Rejected: Radius out of range." << endl;
	  suc=ix_r_outside;
	  i=nsources;
	}
      }
	
      // Ensure non-zero weight
      if (w_next==0.0) {
	scr_out << "Rejected: Zero weight." << endl;
	suc=ix_zero_wgt;
      }
	
    }

    bool force_file_update=false;

    // If the new point is still good, compare with
    // the Metropolis algorithm
    if (suc==ix_success) {

      if (debug) {
	cout << first_half << " Next: " 
	     << e_next.params[0] << " " << w_next << endl;
      }
      
      bool accept;
      if (use_smove) {
	accept=make_step(w_curr_arr[ik],w_next,debug,warm_up,
			 mcmc_iterations,q_current,q_next,smove_z);
      } else {
	accept=make_step(w_current,w_next,debug,warm_up,
			 mcmc_iterations,q_current,q_next,0.0);
      }
      
      if (accept) {

	mh_success++;

	// Add measurement from new point
	shared_ptr<table_units<> > tab_eos;
	shared_ptr<table_units<> > tab_mvsr;

	if (use_smove) {
	  if (step_flags[ik]==false) {
	    tab_eos=data_arr[ik+nwalk].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[ik+nwalk].ts.get_results();
	  } else {
	    tab_eos=data_arr[ik].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[ik].ts.get_results();
	  }
	} else {
	  if (first_half) {
	    tab_eos=data_arr[1].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[1].ts.get_results();
	  } else {
	    tab_eos=data_arr[0].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[0].ts.get_results();
	  }
	}
	tab_eos->set_interp_type(itp_linear);
	tab_mvsr->set_interp_type(itp_linear);

	// Store results from new point
	if (!warm_up) {
	  add_measurement(e_next,tab_eos,tab_mvsr,w_next,true,
			  mh_success,wgts);
	  if (debug) {
	    cout << first_half << " Adding new: " 
		 << e_next.params[0] << " " << w_next << " "
		 << tab_mvsr->max("gm") << endl;
	  }
	}

	// Output the new point
	scr_out << "MC Acc: " << mh_success << " " << e_next << " " 
		<< w_next << endl;
	  
	// Keep track of best point
	if (w_next>w_best) {
	  e_best=e_next;
	  w_best=w_next;
	  output_best(e_best,w_best,tab_eos,tab_mvsr,wgts);
	  force_file_update=true;
	}

	// Prepare for next point
	if (use_smove) {
	  e_curr_arr[ik]=e_next;
	  w_curr_arr[ik]=w_next;
	} else {
	  e_current=e_next;
	  w_current=w_next;
	}

	// Flip "first_half" parameter
	if (use_smove) {
	  step_flags[ik]=!(step_flags[ik]);
	} else {
	  first_half=!(first_half);
	  if (debug) cout << "Flip: " << first_half << endl;
	}
	  
      } else {
	    
	// Point was rejected

	mh_failure++;

	shared_ptr<table_units<> > tab_eos;
	shared_ptr<table_units<> > tab_mvsr;
	
	if (use_smove) {
	  if (step_flags[ik]==false) {
	    tab_eos=data_arr[ik].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[ik].ts.get_results();
	  } else {
	    tab_eos=data_arr[ik+nwalk].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[ik+nwalk].ts.get_results();
	  }
	} else {
	  if (first_half) {
	    tab_eos=data_arr[0].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[0].ts.get_results();
	  } else {
	    tab_eos=data_arr[1].modp->cns.get_eos_results();
	    tab_mvsr=data_arr[1].ts.get_results();
	  }
	}
	tab_eos->set_interp_type(itp_linear);
	tab_mvsr->set_interp_type(itp_linear);

	// Repeat measurement of old point
	if (!warm_up) {
	  if (use_smove) {
	    add_measurement(e_curr_arr[ik],tab_eos,tab_mvsr,
			    w_curr_arr[ik],false,mh_success,wgts);
	    if (debug) {
	      cout << step_flags[ik] << " Adding old: "
		   << e_curr_arr[ik].params[0] << " " << w_curr_arr[ik] << " "
		   << tab_mvsr->max("gm") << endl;
	    }
	  } else {
	    add_measurement(e_current,tab_eos,tab_mvsr,
			    w_current,false,mh_success,wgts);
	    if (debug) {
	      cout << first_half << " Adding old: "
		   << e_current.params[0] << " " << w_current << " "
		   << tab_mvsr->max("gm") << endl;
	    }
	  }
	}

	// Output the old point
	if (use_smove) {
	  scr_out << "MC Rej: " << mh_success << " " << e_curr_arr[ik]
		  << " " << w_curr_arr[ik] << " " << w_next << endl;
	} else {
	  scr_out << "MC Rej: " << mh_success << " " << e_current 
		  << " " << w_current << " " << w_next << endl;
	}

	// Keep track of best point
	if (w_next>w_best) {
	  e_best=e_next;
	  w_best=w_next;
	  output_best(e_best,w_best,tab_eos,tab_mvsr,wgts);
	  force_file_update=true;
	  scr_out << "Best point with rejected step: " << w_next << " " 
		  << w_best << endl;
	}

      }
	  
      // End of "if (suc==ix_success)"
    }

    // ---------------------------------------------------------------

    // After the warm-up is over, the calculation is abritrarily
    // broken up into 20 blocks. The purpose of these blocks is simply
    // to allow easier tracking of progress and to force periodic file
    // updates.

    // This section determines if we have finished a block or if we
    // have finished the full calculation. Note that the value of
    // mcmc_iterations isn't incremented until later.
    
    if (warm_up==false) {

      // If 'max_iters' is zero, then presume we're running over
      // a fixed time interval
      if (max_iters==0) {

	// Determine time elapsed
#ifndef NO_MPI
	double elapsed=MPI_Wtime()-mpi_start_time;
#else
	double elapsed=time(0)-mpi_start_time;
#endif
	// Force a file update when we've finished a block or if
	// we've run out of time
	if (elapsed>max_time/((double)20)*((double)(block_counter+1)) ||
	    (elapsed>max_time && block_counter==19)) {
	  force_file_update=true;
	  block_counter++;
	  scr_out << "Finished block " << block_counter << " of 20." << endl;
	}

	// Output elapsed time every 10 iterations. The value of
	// mcmc_iterations isn't increased until later.
	if ((mcmc_iterations+1)%10==0) {
	  scr_out << "Elapsed time: " << elapsed << " of " << max_time
		  << " seconds" << endl;
	}
	
	if (elapsed>max_time) {
	  main_done=true;
	}

      } else {

	// Otherwise, 'max_iters' is non-zero, so we're completing
	// a fixed number of iterations.

	// Force a file update when we've finished a block
	if (((int)mcmc_iterations)+1>max_iters*(((int)block_counter)+1)/20) {
	  force_file_update=true;
	  block_counter++;
	  scr_out << "Iteration " << mcmc_iterations+1 << " of " 
		  << max_iters << ", and finished block " 
		  << block_counter << " of 20." << endl;
	}

	if (((int)mcmc_iterations)+1>max_iters) {
	  scr_out << "Iteration count, " << mcmc_iterations 
		  << ", exceed maximum number, " << max_iters << "." << endl;
	  main_done=true;
	}
	
      }
     
    }

    // --------------------------------------------------------------
    // Store a copy of measurements in file if 'force_file_update' is
    // true and for a fixed interval of MCMC successes and if the
    // table is at the maximum size. By default file_default_iters is
    // 10 and so the files are updated for every 10 MCMC successes.
    
    if (!warm_up && (force_file_update ||
		     ((int)tc.get_nlines())==max_chain_size || 
		     (mcmc_iterations+1) % file_update_iters==0)) {
      scr_out << "Updating files." << endl;
      update_files(*(data_arr[0].modp),e_current);
      scr_out << "Done updating files." << endl;
    }

    // --------------------------------------------------------------

    // Increment iteration counter
    mcmc_iterations++;

    // Leave warm_up mode if necessary
    if (((int)mcmc_iterations)>n_warm_up && warm_up==true) {
      warm_up=false;
      scr_out << "Setting warm_up to false. Reset start time." << endl;
#ifndef NO_MPI
      max_time-=MPI_Wtime()-mpi_start_time;
      scr_out << "Resetting max_time to : " << max_time << endl;
      mpi_start_time=MPI_Wtime();
#else
      max_time-=time(0)-mpi_start_time;
      scr_out << "Resetting max_time to : " << max_time << endl;
      mpi_start_time=time(0);
#endif
      scr_out.precision(12);
      scr_out << " Start time: " << mpi_start_time << endl;
      scr_out.precision(6);
    }
    
    // End of main loop
  }
  
  return 0;
    }

    /// Main wrapper for parsing command-line arguments
    virtual void run(int argc, char *argv[]) {
  
      // ---------------------------------------
      // Set error handler for this thread
  
      o2scl::err_hnd=&error_handler;
  
      // ---------------------------------------
      // Process command-line arguments and run
  
      setup_cli();

#ifndef NO_MPI
      // Get MPI rank, etc.
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpi_nprocs);
#endif

      // Process arguments
      for(int i=0;i<argc;i++) {
	cl_args.push_back(argv[i]);
      }

      cl.prompt="mcmc> ";
      cl.run_auto(argc,argv);

      if (file_opened) {
	// Close main output file
	scr_out.close();
      }
 
      return;
    }    

    mcmc_class() {

      // Parameters
      prefix="mcmc";
      file_update_iters=40;
      max_chain_size=10000;
      max_iters=0;
      user_seed=0;
      n_warm_up=0;
      // Default to 24 hours
      max_time=3.6e3*24;
      output_next=true;

      // MC step parameters
      use_smove=false;
      hg_mode=0;
      step_fac=15.0;
      nwalk=10;

      // Initial values
      mpi_nprocs=1;
      mpi_rank=0;
      chain_size=0;
      n_chains=0;
    }
    
  };
  
  /** \brief Desc
   */
  class eos_tov {
    
  public:

    /// The model for the EOS
    model *modp;

    /// Desc
    ubvector rad;
    
    /// TOV solver
    o2scl::tov_solve ts;
    
    eos_tov() {
      modp=0;
      ts.verbose=0;
      ts.set_units("1/fm^4","1/fm^4","1/fm^3");
      ts.err_nonconv=false;
    }

    /** \brief Desc
     */
    eos_tov(const eos_tov &e) {
      modp=0;
    }
    
    /** \brief Desc
     */
    eos_tov &operator=(const eos_tov &e) {
      if (this!=&e) {
	modp=0;
      }
      return *this;
    }
    
  };
  
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
  class bamr_class : public mcmc_class<eos_tov> {
    
  protected:

    /// \name Member data for the Metropolis-Hastings step
    //@{
    /// Return the approximate likelihood
    double approx_like(entry &e);
    //@}

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

    /// \name Other variables
    //@{
    /// The first point in the parameter space
    ubvector first_point;
    
    /// The file containing the initial point
    std::string first_point_file;

    /// \name Integer designating how to set the initial point
    //@{
    int first_point_type;
    static const int fp_unspecified=-1;
    static const int fp_last=-2;
    static const int fp_best=-3;
    //@}

    /// If true, then \ref first_update() has been called
    bool first_file_update;

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
    std::vector<int> ret_codes;
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
