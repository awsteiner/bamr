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
/** \file bamr_class.h
    \brief Definition of main bamr class
*/
#ifndef BAMR_CLASS_H
#define BAMR_CLASS_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <Python.h>

#ifdef BAMR_MPI
#include <mpi.h>
#endif

#include "nstar_cold2.h"
#include "models.h"
#include "filters.h"

/** \brief Main namespace
    
    The bamr namespace which holds all classes and functions.
    
    This file is documented in bamr_class.h .
*/
namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  typedef std::function<int(size_t,const ubvector &, double &,
                            model_data &)> point_funct;
  
  typedef std::function<int(const ubvector &,double,
                            std::vector<double> &,model_data &)> fill_funct;

  typedef std::function<int(size_t,ubvector &,point_funct &,
                            ubvector &, model_data &, bool &)> deriv_funct;

  /** \brief Compute neutron star structure for each MCMC point

      There will be a number of instances of this class equal to the
      number of OpenMP threads, all of which share the same neutron
      star data (in \ref nsd) and settings (in \ref set). They will
      each have their own model object, stored in \ref mod, which is
      set in \ref mcmc_bamr::set_model() .
  */
  class bamr_class {

  public:
    /*
      PyObject *train_modFile;
      PyObject *train_tParam_Names;
      PyObject *train_trainClass;
      PyObject *train_instance;
      PyObject *train_trainMthd;
      PyObject *train_pArgs;
      PyObject *addtl_sources;
      PyObject *train_res;
      PyObject *train_pTemp;
      PyObject *train_temp;
      PyObject *test_show;
      PyObject *test_vals;
      PyObject *target_cols;
      PyObject *target_pred;
    */
    
    /** \brief If true, include emulator from sklearn
     */
#ifdef O2SCL_NEVER_DEFINED
    bool apply_emu;
#endif
    
    /** \brief If true, use index2 to take derivative of M_max
     */
    bool dv_index2;

    /** \brief Train file name for python emulator
     */
    //std::string emu_train;    

    //bool py_train;

    void setup_filters();
    
    // -------------------------------------------------------
    // New bint variables

    /** \brief A counter of intrinsic scatter calculations
        for debugging
     */
    int intsc_counter;

    /** \brief The number of OpenMP threads (set where?)
     */
    int n_threads;

#ifdef BAMR_FFTW3    
    /** \brief The filter objects (one for each thread)
     */
    std::vector<filters::Filter *> flt;
#endif
    
    /** \brief Desc
     */
    o2scl::tensor_grid<> fft_data[22];
    
    /// Copy of the original NS measurement data
    std::vector<o2scl::table3d> source_tables_is;

    /// Copy of the original NS measurement data
    std::vector<o2scl::table3d> source_tables_alt_is;

    // -------------------------------------------------------

    /// The Schwarzchild radius in km
    double schwarz_km;
    
    /// Pointer to neutron star data
    std::shared_ptr<ns_data> nsd;
    
    /** \brief Pointer to settings object
        
        AWS: 4/5/2020 changed from "const settings" to "settings"
        to enable the python interface
    */
    std::shared_ptr<settings> set;
    
    /// Model object
    std::shared_ptr<model> mod;

    /// Model type string
    std::string model_type;

    /// MCMC method type
    std::string mcmc_method;

    /// Initial point file name
    std::string init_file;

    /// The index of parameters
    o2scl::vec_index pvi;
    
    /// Vector to store log-weights to be passed to table
    std::vector<double> wgt_pop;

    /// Vector to store the output quantities for GW190425
    double wgt_gw19;
    std::vector<double> fsn_gw19;
    std::vector<double> mass_gw19;

    /// Vector to store the output quantities for GW170817
    double wgt_gw17;
    std::vector<double> fsn_gw17;
    std::vector<double> mass_gw17;

    /// Vector to store the output quantities for EM sources
    std::vector<double> wgt_em;
    std::vector<double> fsn_em;

    bamr_class() {
      schwarz_km=o2scl_const::schwarzchild_radius_f<double>()/1.0e3;
      n_threads=1;
    }

    // Empty destructor to make sure its virtual
    virtual ~bamr_class() {
    }

    /** \brief Compute the EOS corresponding to parameters in 
        \c e and put output in \c tab_eos
    */
    virtual int compute_point(const ubvector &pars, std::ofstream &scr_out, 
                              double &log_wgt, model_data &dat);

    virtual int compute_point_ext(const ubvector &pars, std::ofstream &scr_out, 
                                  double &log_wgt, model_data &dat);

    virtual int compute_deriv(ubvector &, point_funct &,
                              ubvector &, model_data &, bool &);

    virtual int numeric_deriv(size_t, ubvector &, point_funct &, 
                              double &, double &, model_data &);
    
    virtual int compute_gw17(const ubvector &, double &, model_data &);
    virtual int compute_gw19(const ubvector &, double &, model_data &);
    virtual int compute_ems(size_t, const ubvector &, double &, model_data &);
    virtual int compute_dist(size_t, const ubvector &, double &, model_data &);
    
    /** \brief Fill vector in <tt>line</tt> with data from the
        current Monte Carlo point
    */
    virtual int fill(const ubvector &pars, double weight, 
                     std::vector<double> &line, model_data &dat);

  protected:

    /// If true, initialize the atm parameters
    bool init_eval=true;

    /// If true, atm parameters are held fixed
    bool atms_fixed=false;

    /// Vector to store atm parameters
    std::vector<bool> atms;

    /// Vector to store gradients at current point
    ubvector grad2;

    /// Function for Random Walk over atm parameters
    virtual void compute_atms(const ubvector &, model_data &);
    
  };

}

#endif
