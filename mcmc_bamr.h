/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
  and Andrew W. Steiner
  
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

#ifdef BAMR_MPI
#include <mpi.h>
#endif

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

#include <o2scl/hdf_file.h>
#include <o2scl/mcmc_para.h>
#include <o2scl/kde_python.h>

#ifdef O2SCL_NEVER_DEFINED
#include "emulator_bamr.h"
#endif

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

      \todo There's a problem with testing debug_eos and debug_star
      because of the complications of calling exit(0) inside
      an OpenMP/MPI parallel region. For this reason test1 and
      test2 are temporarily disabled. I think MPI_Barrier() may
      help this problem, but this fix hasn't been tested and
      may not completely fix it.

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
//#ifdef BAMR_KDE
/*  class mcmc_bamr :
    public o2scl::mcmc_para_cli
  <point_funct,fill_funct,model_data,ubvector,
   mcmc_stepper_mh<point_funct,model_data,
                   ubvector,ubmatrix,
                   prob_cond_mdim_indep<>>>
*/
//#else

#ifdef ANDREW
  class mcmc_bamr :
    public o2scl::mcmc_para_emu
  <point_funct,fill_funct,model_data,ubvector>
#else
  class mcmc_bamr :
    public o2scl::mcmc_para_cli
  <point_funct,fill_funct,model_data,ubvector,
    mcmc_stepper_hmc<point_funct,model_data,ubvector>
#endif
   
//#endif
   
  {
  protected:

    /** \brief The input tensor for the KDE proposal distribution
     */
    o2scl::tensor<> ten_in;

    /// The KDE proposal distribution
    std::shared_ptr<o2scl::kde_python<ubvector>> kp;
    
    /// The Gaussian proposal distribution
    std::shared_ptr<o2scl::prob_dens_mdim_gaussian<>> gpp;
    
    /** \brief If true, use index2 to take derivative of M_max
     */
    bool dv_index2;

#ifdef O2SCL_NEVER_DEFINED
    
    /** \brief Train file name for python emulator
     */
    std::string emu_train;

    PyObject *train_modFile;
    PyObject *train_tParam_Names;
    PyObject *train_trainClass;
    PyObject *train_instance;
    PyObject *train_trainMthd;
    PyObject *train_pArgs;

    /// The number of sources
    PyObject *addtl_sources;
    
    PyObject *train_res;
    PyObject *train_pTemp;
    PyObject *train_temp;

    bool py_train;

    /** \brief Train the emulator based on data from \c file_name
        using parameter names listed in \c names.
        
        This is called in mcmc_bamr::mcmc_func().
    */
    int train(std::string file_name, std::vector<std::string> &names);

    /** \brief Calculate posteriors from the emulated points.
     */
    virtual int emu_points(std::vector<std::string> &sv,
                           bool itive_com);
    
    virtual int emu_train2(std::vector<std::string> &sv,
                           bool itive_com);
#endif

    /** \brief A string specifying the sampling algorithm for MCMC
        Current options are: 'kde', 'kde_sklearn', 'gauss' and 'hmc'
    */
    std::string mcmc_method;

    /// A string indicating which model is used, set in \ref set_model().
    std::string model_type;

    /** \brief The index of parameters

        This is created in \ref mcmc_func().
     */
    vec_index pvi;

    /** \brief Vector of \ref bamr_class objects (one for each OpenMP
        thread)

        This is currently a pointer because it makes it a lot easier
        to replace these pointers with children. A shared_ptr might be
        better, but I've had problems implementing vector<shared_ptr>
        correctly.
    */
    std::vector<bamr_class *> bc_arr;

#ifdef O2SCL_NEVER_DEFINED
    std::vector<emulator_bamr> eb_arr;
#endif

    /** \brief The \ref bamr::settings object 
        (shared by instances of \ref bamr_class)
    */
    std::shared_ptr<settings> set;

    /** \brief The \ref bamr::ns_data object
        (shared by instances of \ref bamr_class)
    */
    std::shared_ptr<ns_data> nsd;

    /// \name Main functions called from the command-line interface
    //@{
    /** \brief Set the MCMC sampling method
     */
    virtual int set_method(std::vector<std::string> &sv,
                                bool itive_com);
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
    */
    virtual int mcmc_init();

    /** \brief Add a data distribution to the list
     */
    virtual int add_data(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Add a data distribution to the list
     */
    virtual int add_data_alt(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Perform the MCMC simulation
     */
    virtual int mcmc_func(std::vector<std::string> &sv, bool itive_com);

#ifdef ANDREW
    /// Wrapper to the point function which uses the emulator
    virtual int point_wrapper(size_t it, size_t np, const ubvector &p,
                              double &log_wgt, model_data &dat);
#endif
    
    /** \brief Set the number of OpenMP threads
     */
    virtual int threads(std::vector<std::string> &sv, bool itive_com);

    /** \brief Use the last point in a specifed file for the 
        initial point
    */
    virtual int initial_point_last(std::vector<std::string> &sv,
                                   bool itive_com);
    
    /** \brief Use the highest likelihood point in the specified file
        for the initial point
    */
    virtual int initial_point_best(std::vector<std::string> &sv,
                                   bool itive_com);

    /// Desc
    virtual int combine_files(std::vector<std::string> &sv,
                                   bool itive_com);
    
    /** \brief Read previous results from a file
     */
    virtual int read_prev_results_mb(std::vector<std::string> &sv,
                                     bool itive_com);
    
  public:

    /** \brief The command-line interface object
     */
#ifdef BAMR_READLINE
    o2scl::cli_readline cl;
#else
    o2scl::cli cl;
#endif
    
    /** \brief Create a \ref mcmc_bamr object 
     */
    mcmc_bamr();
    
    virtual ~mcmc_bamr() {
    }
    
    /** \brief Set up the 'cli' object

        This function adds three commands (mcmc, model, add-data) and
        the 'set' parameters.
    */
    virtual void setup_cli_mb();
    
  };

}

#endif
