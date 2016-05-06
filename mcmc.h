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
/** \file mcmc.h
    \brief Definition of main mcmc class
*/
#ifndef MCMC_H
#define MCMC_H

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

#ifndef O2SCL_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

/** \brief Main namespace
    
    The mcmc namespace which holds all classes and functions.
    
    This file is documented in mcmc.h .
*/
namespace mcmc_namespace {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  /** \brief Desc
   */
  class default_model {
    
    /** \brief Desc
     */
    size_t nparams;

    /** \brief Desc
     */
    virtual double compute_point(ubvector &pars, std::ofstream &scr_out,
				 int &success, ubvector &dat)=0;
    
    /** \brief Desc
     */
    virtual void first_point(ubvector &pars) {
      ubvector low(nparams), high(nparams);
      low_limits(low);
      high_limits(high);
      for(size_t i=0;i<nparams;i++) {
	pars[i]=(low[i]+high[i])/2.0;
      }
      return;
    }
    
    /// \name Functions for MCMC parameters
    //@{
    /** \brief Set the lower boundaries for all the parameters
     */
    virtual void low_limits(ubvector &pars)=0;

    /** \brief Set the upper boundaries for all the parameters
     */
    virtual void high_limits(ubvector &pars)=0;

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i)=0;

    /// Return the unit of parameter with index \c i
    virtual std::string param_unit(size_t i)=0;
    //@}

    /// \name Functions for model parameters fixed during the MCMC run
    //@{
    /** \brief Setup model parameters */
    virtual void setup_params(o2scl::cli &cl) {
      return;
    }
    
    /** \brief Remove model parameters */
    virtual void remove_params(o2scl::cli &cl) {
      return;
    }
    //@}

  };

  /** \brief Desc
   */
  template<class data_t=ubvector,
    class model_t=default_model> class mcmc_class {
    
  public:
  
    /// \name Member data for the Metropolis-Hastings step
    //@{
    /// Return the approximate likelihood
  double approx_like(ubvector &pars);
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

  std::vector<int> ret_codes;
  
  /// If true, then \ref first_update() has been called
  bool first_file_update;
  
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
  
#ifdef O2SCL_READLINE
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
  
  /// Map of model objects
  std::map<std::string,model_t,std::greater<std::string> > model_arr;

  /// Current model
  std::string curr_model;
  
  /** \brief Set up the 'cli' object
      
      This function just adds the four commands and the 'set' parameters
  */
  void setup_cli() {

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

  /** \brief Desc
   */
  int mcmc_init() {
    return 0;
  };
  
  /** \brief Desc
   */
  virtual void output_best
  (ubvector &best, double w_best, data_t &d);

  /** \brief Desc
   */
  int mcmc(std::vector<std::string> &sv, bool itive_com) {

    // Get model object
    typename std::map<std::string,model_t,
      std::greater<std::string> >::const_iterator model_it=
      model_arr.find(curr_model);
    if (model_it==model_arr.end()) {
      std::cerr << "Model type " << curr_model << " unknown." << std::endl;
      return 2;
    }
    model_t &m=model_it->second;

    // Set number of parameters
    nparams=m.nparams;
    
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
      scr_out.setf(std::ios::scientific);
      file_opened=true;
      scr_out << "Opened main file in command 'mcmc'." << std::endl;
    }
      
    // Fix file_update_iters if necessary
    if (file_update_iters<1) {
      scr_out << "Parameter 'file_update_iters' less than 1. Set equal to 1."
	      << std::endl;
      file_update_iters=1;
    }
      
    if (max_chain_size<1) {
      O2SCL_ERR("Parameter 'max_chain_size' must be larger than 1.",
		o2scl::exc_einval);
    }
      
    // Fix step_fac if it's too small
    if (step_fac<1.0) {
      step_fac=1.0;
      scr_out << "Fixed 'step_fac' to 1.0." << std::endl;
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
	    << mpi_nprocs << "." << std::endl;
    scr_out.precision(12);
    scr_out << " Start time: " << mpi_start_time << std::endl;
    scr_out.precision(6);
      
    // First MC point
    std::vector<ubvector> current(1);
    std::vector<double> w_current(1);
    std::vector<bool> step_flags(1);
    data_arr.resize(2);
    w_current[0]=0.0;

    // For stretch-moves, allocate for each walker
    if (use_smove) {
      current.resize(nwalk);
      w_current.resize(nwalk);
      step_flags.resize(nwalk);
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<nwalk;i++) {
	step_flags[i]=false;
	w_current[i]=0.0;
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
		<< "'." << std::endl;
	o2scl_hdf::hdf_file hf;
	hf.open(first_point_file);
      
	// Read table
	size_t file_n_chains;
	hf.get_szt("n_chains",file_n_chains);
	std::string chain_name=std::string("markov_chain")+
	  o2scl::szttos(file_n_chains-1);
	o2scl::table_units<> file_tab;
	hdf_input(hf,file_tab,chain_name);
	size_t last_line=file_tab.get_nlines()-1;
      
	// Get parameters
	for(size_t i=0;i<nparams;i++) {
	  std::string pname=((std::string)"param_")+m.param_name(i);
	  current[0][i]=file_tab.get(pname,last_line);
	  scr_out << "Parameter named " << m.param_name(i) << " " 
		  << current[0][i] << std::endl;
	}
      
	// Finish up
	scr_out << std::endl;
	hf.close();

      } else if (first_point_type==fp_best) {
	
	std::vector<double> best_point;
	o2scl_hdf::hdf_file hf;
	hf.open(first_point_file);
	hf.getd_vec("best_point",best_point);
	hf.close();
	scr_out << "Reading best point from file '" << first_point_file
		<< "'." << std::endl;
	for(size_t i=0;i<nparams;i++) {
	  current[0][i]=best_point[i];
	  scr_out << "Parameter " << i << " : " << current[0][i] << std::endl;
	}
	scr_out << "Best weight: " << best_point[nparams] << std::endl;
	scr_out << std::endl;

      } else {

	// Read file 
	scr_out << "Reading " << first_point_type << "th point from file '" 
		<< first_point_file
		<< "'." << std::endl;
	o2scl_hdf::hdf_file hf;
	hf.open(first_point_file);
      
	// Read table
	size_t file_n_chains, row=first_point_type;
	hf.get_szt("n_chains",file_n_chains);
      
	o2scl::table_units<> file_tab;
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
		  << " in file. Using last point." << std::endl;
	  row=file_tab.get_nlines()-1;
	}
      
	// Get parameters
	for(size_t i=0;i<nparams;i++) {
	  std::string pname=((std::string)"param_")+m.param_name(i);
	  current[0][i]=file_tab.get(pname,row);
	  scr_out << "Parameter named " << m.param_name(i) << " " 
		  << current[0][i] << std::endl;
	}
      
	// Finish up
	scr_out << std::endl;
	hf.close();
      }

    } else if (first_point.size()>0) {
    
      scr_out << "First point from command-line." << std::endl;
      for(size_t i=0;i<nparams;i++) {
	current[0][i]=first_point[i];
	scr_out << current[0][i] << std::endl;
      }
      scr_out << std::endl;

    } else {
    
      scr_out << "First point from default." << std::endl;
      m.first_point(current[0]);

    }

#ifndef NO_MPI
    if (mpi_nprocs>1 && mpi_rank<mpi_nprocs-1) {
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,tag,MPI_COMM_WORLD);
    }
#endif

    scr_out << "First point: ";
    vector_out(scr_out,current[0],true);

    // Determine initial masses
    std::cout << "fixme." << std::endl;
    exit(-1);

    //for(size_t i=0;i<nsources;i++) {
    //e_current.mass[i]=first_mass[i];
    //e_current.rad[i]=0.0;
    //}

    // Set lower and upper bounds for parameters
    ubvector low(nparams), high(nparams);
    m.low_limits(low);
    m.high_limits(high);

    n_chains=0;

    // Entry objects (Must be after read_input() since nsources is set
    // in that function.)
    ubvector next(nparams), best(nparams);

    // Weights for each entry
    double w_next=0.0, w_best=0.0;

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
	  ubvector first(nparams);
	  m.first_point(first);

	  // Make a perturbation from the initial point
	  for(size_t ik=0;ik<nparams;ik++) {
	    do {
	      current[ij][ik]=first[ik]+
		(gr.random()*2.0-1.0)*(high[ik]-low[ik])/100.0;
	    } while (current[ij][ik]>=high[ik] ||
		     current[ij][ik]<=low[ik]);
	  }
	
	  // Compute the weight
	  w_current[ij]=m.compute_point(current[ij],scr_out,suc,data_arr[ij]);
	  scr_out << "SM Init: " << ij << " "
		  << current[ij] << " " << w_current[ij] << std::endl;

	  // Keep track of the best point and the best index
	  if (ij==0) {
	    w_best=w_current[0];
	  } else if (w_current[ij]>w_best) {
	    ij_best=ij;
	    w_best=w_current[ij];
	  }

	  // Increment iteration count
	  init_iters++;

	  // If we have a good point, stop the loop
	  if (w_current[ij]>0.0) {
	    done=true;
	    ret_codes[suc]++;
	  } else if (init_iters>1000) {
	    scr_out << "Failed to construct initial walkers." << std::endl;
	    return 1;
	  }
	}

	// For the initial point for this walker, add it
	// to the result table
	if (warm_up==false) {
	  // Add the initial point if there's no warm up
	  add_measurement(current[ij],w_current[ij],data_arr[ij],
			  true,mh_success);
	}
      }

      // Output the best initial walker if necessary
      {
	best=current[ij_best];
	output_best(best,w_best,data_arr[ij_best]);
      }

    } else {
      // Normal or Metropolis-Hastings steps

      // Compute weight for initial point
      w_current[0]=m.compute_point(current[0],scr_out,suc,data_arr[0]);
      ret_codes[suc]++;
      scr_out << "Initial weight: " << w_current << std::endl;
      if (w_current[0]<=0.0) {
	scr_out << "Initial weight zero. Aborting." << std::endl;
	exit(-1);
      }

      // Compute the initial Hastings proposal weight
      if (hg_mode>0) {
	q_current=approx_like(current[0]);
      }

      // Add measurement to output table and output best 
      {
	if (warm_up==false) {
	  // Add the initial point if there's no warm up
	  add_measurement(current[0],w_current,data_arr[0],
			  true,mh_success);
	}
	best=current[0];
	w_best=w_current[0];
	output_best(current[0],w_current[0],data_arr[0]);
      }
    
    }

    // ---------------------------------------------------

    // Initialize radii of next to zero
    for(size_t k=0;k<nsources;k++) {
      next.rad[k]=0.0;
    }

    // Keep track of total number of different points in the parameter
    // space that are considered (some of them may not result in TOV
    // calls because, for example, the EOS was acausal).
    mcmc_iterations=0;

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
	    next[i]=current[ij][i]+smove_z*
	      (current[ik][i]-current[ij][i]);
	    if (next[i]>=high[i] ||
		next[i]<=low[i]) {
	      in_bounds=false;
	    }
	  }
	  for(size_t i=0;i<nsources;i++) {
	    next.mass[i]=current[ij].mass[i]+smove_z*
	      (current[ik].mass[i]-current[ij].mass[i]);
	    if (next.mass[i]>=high.mass[i] ||
		next.mass[i]<=low.mass[i]) {
	      in_bounds=false;
	    }
	  }
	
	  step_iters++;
	  if (step_iters==1000) {
	    scr_out << "Failed to find suitable step at point 1." << std::endl;
	    cerr << "Failed to find suitable step at point 1." << std::endl;
	    return 2;
	  }

	} while (in_bounds==false);

      } else if (hg_mode>0) {
      
	// Make a Metropolis-Hastings step based on previous data
      
	size_t nv=next.np+nsources;
	ubvector hg_temp(nv), hg_z(nv);
      
	bool out_of_range;
	int hg_it=0;
      
	do {
	
	  for(size_t k=0;k<nv;k++) {
	    hg_z[k]=pdg.sample();
	  }
	  hg_temp=prod(hg_chol,hg_z);
	  for(size_t k=0;k<nv;k++) {
	    if (k<next.np) {
	      next[k]=hg_best[k]+hg_temp[k];
	    } else {
	      next.mass[k-next.np]=hg_best[k]+hg_temp[k];
	    }
	  }
	
	  out_of_range=false;
	  for(size_t k=0;k<next.np;k++) {
	    if (next[k]<low[k] ||
		next[k]>high[k]) {
	      out_of_range=true;
	    }
	  }
	  for(size_t k=0;k<nsources;k++) {
	    if (next.mass[k]<low.mass[k] ||
		next.mass[k]>high.mass[k]) {
	      out_of_range=true;
	    }
	  }

	  hg_it++;
	  if (hg_it>1000) {
	    O2SCL_ERR("Sanity check in hg step.",exc_esanity);
	  }

	} while (out_of_range==true);

	q_next=approx_like(next);

      } else {

	// Make a step, ensure that we're in bounds and that
	// the masses are not too large
	for(size_t k=0;k<next.np;k++) {
	
	  next[k]=current[0][k]+(gr.random()*2.0-1.0)*
	    (high[k]-low[k])/step_fac;
	
	  // If it's out of range, redo step near boundary
	  if (next[k]<low[k]) {
	    next[k]=low[k]+gr.random()*
	      (high[k]-low[k])/step_fac;
	  } else if (next[k]>high[k]) {
	    next[k]=high[k]-gr.random()*
	      (high[k]-low[k])/step_fac;
	  }
	
	  if (next[k]<low[k] || 
	      next[k]>high[k]) {
	    O2SCL_ERR("Sanity check in parameter step.",exc_esanity);
	  }
	}
      
	if (nsources>0) {
	  // Just use a large value (1.0e6) here since we don't yet
	  // know the maximum mass
	  select_mass(current[0],next,1.0e6);
	}
      
      }

      // End of select next point
      // ---------------------------------------------------

      // Output the next point
      if (output_next) {
	scr_out << "Iteration, next: " << mcmc_iterations << " " 
		<< next << std::endl;
      }
      
      // ---------------------------------------------------
      // Compute next weight

      if (use_smove) {
	if (step_flags[ik]==false) {
	  w_next=compute_weight(next,data_arr[ik+nwalk],suc,wgts,warm_up);
	} else {
	  w_next=compute_weight(next,data_arr[ik],suc,wgts,warm_up);
	}
      } else {
	if (step_flags[0]) {
	  w_next=compute_weight(next,data_arr[1],suc,wgts,warm_up);
	} else {
	  w_next=compute_weight(next,data_arr[0],suc,wgts,warm_up);
	}
      }
      ret_codes[suc]++;

      // ---------------------------------------------------
    
      // Test to ensure new point is good
      if (suc==ix_success) {

	// Test radii
	for(size_t i=0;i<nsources;i++) {
	  if (next.rad[i]>high.rad[i] || next.rad[i]<low.rad[i]) {
	    scr_out << "Rejected: Radius out of range." << std::endl;
	    suc=ix_r_outside;
	    i=nsources;
	  }
	}
	
	// Ensure non-zero weight
	if (w_next==0.0) {
	  scr_out << "Rejected: Zero weight." << std::endl;
	  suc=ix_zero_wgt;
	}
	
      }

      bool force_file_update=false;

      // If the new point is still good, compare with
      // the Metropolis algorithm
      if (suc==ix_success) {

	if (debug) {
	  cout << step_flags[0] << " Next: " 
	       << next[0] << " " << w_next << std::endl;
	}
      
	bool accept;
	if (use_smove) {
	  accept=make_step(w_current[ik],w_next,debug,warm_up,
			   mcmc_iterations,q_current,q_next,smove_z);
	} else {
	  accept=make_step(w_current,w_next,debug,warm_up,
			   mcmc_iterations,q_current,q_next,0.0);
	}
      
	if (accept) {

	  mh_success++;

	  // Store results from new point
	  if (!warm_up) {
	    if (use_smove) {
	      if (step_flags[ik]==false) {
		add_measurement(next,data_arr[ik+nwalk],w_next,true,
				mh_success,wgts);
	      } else {
		add_measurement(next,data_arr[ik],w_next,true,
				mh_success,wgts);
	      }
	    } else {
	      if (step_flags[0]==false) {
		add_measurement(next,data_arr[1],w_next,true,
				mh_success,wgts);
	      } else {
		add_measurement(next,data_arr[0],w_next,true,
				mh_success,wgts);
	      }
	    }
	    if (debug) {
	      cout << step_flags[0] << " Adding new: " 
		   << next[0] << " " << w_next << " "
		   << tab_mvsr->max("gm") << std::endl;
	    }
	  }

	  // Output the new point
	  scr_out << "MC Acc: " << mh_success << " " << next << " " 
		  << w_next << std::endl;
	  
	  // Keep track of best point
	  if (w_next>w_best) {
	    best=next;
	    w_best=w_next;
	    output_best(best,w_best,tab_eos,tab_mvsr,wgts);
	    force_file_update=true;
	  }

	  // Prepare for next point
	  if (use_smove) {
	    current[ik]=next;
	    w_current[ik]=w_next;
	  } else {
	    current[0]=next;
	    w_current=w_next;
	  }

	  // Flip "first_half" parameter
	  if (use_smove) {
	    step_flags[ik]=!(step_flags[ik]);
	  } else {
	    step_flags[0]=!(step_flags[0]);
	    if (debug) cout << "Flip: " << step_flags[0] << std::endl;
	  }
	  
	} else {
	    
	  // Point was rejected

	  mh_failure++;


	  // Repeat measurement of old point
	  if (!warm_up) {
	    if (use_smove) {
	      if (step_flags[ik]==false) {
		add_measurement(current[ik],data_arr[ik],
				w_current[ik],false,mh_success,wgts);
	      } else {
		add_measurement(current[ik],data_arr[ik+nwalk],
				w_current[ik],false,mh_success,wgts);
	      }
	      if (debug) {
		cout << step_flags[ik] << " Adding old: "
		     << current[ik][0] << " " << w_current[ik] << " "
		     << tab_mvsr->max("gm") << std::endl;
	      }
	    } else {
	      if (step_flags[0]==false) {
		add_measurement(current[0],data_arr[0],
				w_current,false,mh_success,wgts);
	      } else {
		add_measurement(current[0],data_arr[1],
				w_current,false,mh_success,wgts);
	      }
	      if (debug) {
		cout << step_flags[0] << " Adding old: "
		     << current[0][0] << " " << w_current << " "
		     << tab_mvsr->max("gm") << std::endl;
	      }
	    }
	  }

	  // Output the old point
	  if (use_smove) {
	    scr_out << "MC Rej: " << mh_success << " " << current[ik]
		    << " " << w_current[ik] << " " << w_next << std::endl;
	  } else {
	    scr_out << "MC Rej: " << mh_success << " " << current[0] 
		    << " " << w_current << " " << w_next << std::endl;
	  }

	  // Keep track of best point
	  if (w_next>w_best) {
	    best=next;
	    w_best=w_next;
	    output_best(best,w_best,tab_eos,tab_mvsr,wgts);
	    force_file_update=true;
	    scr_out << "Best point with rejected step: " << w_next << " " 
		    << w_best << std::endl;
	  }

	}
	  
	// End of "if (suc==ix_success)"
      }

      // ---------------------------------------------------------------

      // After the warm-up is over, the calculation is abritrarily
      // broken up into 20 blocks. The purpose of these blocks is
      // simply to allow easier tracking of progress and to force
      // periodic file updates.

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
	    scr_out << "Finished block " << block_counter << " of 20." << std::endl;
	  }

	  // Output elapsed time every 10 iterations. The value of
	  // mcmc_iterations isn't increased until later.
	  if ((mcmc_iterations+1)%10==0) {
	    scr_out << "Elapsed time: " << elapsed << " of " << max_time
		    << " seconds" << std::endl;
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
		    << block_counter << " of 20." << std::endl;
	  }

	  if (((int)mcmc_iterations)+1>max_iters) {
	    scr_out << "Iteration count, " << mcmc_iterations 
		    << ", exceed maximum number, " << max_iters << "." << std::endl;
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
	scr_out << "Updating files." << std::endl;
	update_files(*(data_arr[0].modp),current[0]);
	scr_out << "Done updating files." << std::endl;
      }

      // --------------------------------------------------------------

      // Increment iteration counter
      mcmc_iterations++;

      // Leave warm_up mode if necessary
      if (((int)mcmc_iterations)>n_warm_up && warm_up==true) {
	warm_up=false;
	scr_out << "Setting warm_up to false. Reset start time." << std::endl;
#ifndef NO_MPI
	max_time-=MPI_Wtime()-mpi_start_time;
	scr_out << "Resetting max_time to : " << max_time << std::endl;
	mpi_start_time=MPI_Wtime();
#else
	max_time-=time(0)-mpi_start_time;
	scr_out << "Resetting max_time to : " << max_time << std::endl;
	mpi_start_time=time(0);
#endif
	scr_out.precision(12);
	scr_out << " Start time: " << mpi_start_time << std::endl;
	scr_out.precision(6);
      }
    
      // End of main loop
    }
  
    return 0;
  }

  /// Main wrapper for parsing command-line arguments
  void run(int argc, char *argv[]) {
  
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

  // End of mcmc_namespace namespace
}

#endif
