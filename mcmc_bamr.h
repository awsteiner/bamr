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
/** \file mcmc_bamr.h
    \brief Definition of main bamr class
*/
#ifndef MCMC_BAMR_H
#define MCMC_BAMR_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#ifndef BAMR_NO_MPI
#include <mpi.h>
#endif

#include <o2scl/hdf_file.h>
#include <o2scl/mcmc_para.h>

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#include "bamr_class.h"

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
  class mcmc_bamr :
    public o2scl::mcmc_para_cli<point_funct,fill_funct,model_data,ubvector> {

  protected:

    /// A string indicating which model is used, set in \ref set_model().
    std::string model_type;
    
    /// Vector of \ref bamr_class objects (one for each OpenMP thread)
    std::vector<bamr_class> bc_arr;

    /// Settings object
    std::shared_ptr<settings> set;

    /// Data object
    std::shared_ptr<ns_data> nsd;

    /// Number of OpenMP threads
    size_t n_threads;
    
    /// \name Main functions called from the command-line interface
    //@{
    /** \brief Set the model for the EOS to use
     */
    virtual int set_model(std::vector<std::string> &sv, bool itive_com);
    //@}

    /** \brief Write initial data to HDF file
     */
    virtual void file_header(o2scl_hdf::hdf_file &hf);
    
    /** \brief Make any necessary preparations for the mcmc() function

	This is called by \ref mcmc(). If the return value is non-zero
	then it is assumed that the calculation fails and mcmc()
	returns.

	This function does nothing and returns zero by default.
    */
    virtual int mcmc_init();

    /** \brief Add a data distribution to the list
     */
    virtual int add_data(std::vector<std::string> &sv, bool itive_com);

    /** \brief Perform the MCMC simulation
     */
    virtual int mcmc_func(std::vector<std::string> &sv, bool itive_com);
    
  public:

    /** \brief The command-line interface object
     */
#ifdef BAMR_READLINE
    o2scl::cli_readline cl;
#else
    o2scl::cli cl;
#endif

    /** \brief Create a \ref mcmc_bamr object with the specified
	number of OpenMP threads (default 1)
     */
    mcmc_bamr(size_t n_omp_threads=1) {
      model_type="";
      n_threads=n_omp_threads;
      bc_arr.resize(n_threads);
      set=std::shared_ptr<settings>(new settings);
      nsd=std::shared_ptr<ns_data>(new ns_data);
      for(size_t i=0;i<n_threads;i++) {
	bc_arr[i].set=set;
	bc_arr[i].nsd=nsd;
      }
    }
    
    virtual ~mcmc_bamr() {
    }
    
    /** \brief Set up the 'cli' object

	This function adds three commands (mcmc, model, add-data) and
	the 'set' parameters.
    */
    virtual void setup_cli();
    
  };

}

#endif
