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

  /** \brief Compute neutron star structure for each MCMC point

      There will be a number of instances of this class equal to the
      number of OpenMP threads, all of which share the same neutron
      star data (in \ref nsd) and settings (in \ref set). They will
      each have their own model object, stored in \ref mod, which is
      set in \ref mcmc_bamr::set_model() .
  */
  class bamr_class {

  public:

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

    /** \brief If true, include emulator from sklearn
     */
    //bool apply_emu;

    /** \brief If true, use index2 to take derivative of M_max
     */
    bool dv_index2;

    /** \brief Train file name for python emulator
     */
    //std::string emu_train;    

    bool py_train;

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
    
    /** \brief The LIGO data
     */
    o2scl::tensor_grid<> ligo_data_table;

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

    /** \brief Desc
     */
    class py_param_info {
      
    public:

      int np;
      std::vector<double> low;
      std::vector<double> high;
      std::vector<std::string> names;
      std::vector<std::string> units;
      std::vector<int> name_counts;
      std::vector<int> unit_counts;
      std::vector<char> name_c;
      std::vector<char> unit_c;
    };

    /// Desc
    py_param_info ppi;
    
    /// Desc
    vec_index pvi;
    
    /// Vector to store log_weights to be passed to table
    vector<double> pop_weights(4);

    bamr_class() {
      schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;
    }

    // Empty destructor to make sure its virtual
    virtual ~bamr_class() {
    }

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual int compute_point(const ubvector &pars, std::ofstream &scr_out, 
			      double &log_wgt, model_data &dat);
    
    /** \brief Fill vector in <tt>line</tt> with data from the
	current Monte Carlo point
    */
    virtual int fill(const ubvector &pars, double weight, 
		     std::vector<double> &line, model_data &dat);
    
  };

}

extern "C" {

  /** \brief Create a \ref bamr_class and \ref model_data object and
      return the associated pointers
  */
  void create_pointers(char *model_name, void *&bcp2, void *&mdp2,
		       void *&nsd2, void *&setp2, char *data_dir,
		       int verbose);

  /** \brief Set a parameter
   */
  void set_parameter(void *bcp2, void *setp2, char *param_name, double val);

  /** \brief Set a string parameter
   */
  void set_parameter_string(void *bcp2, void *setp2, char *param_name,
			    char *val);

  /** \brief Desc
   */
  int init(void *bcp2, void *mdp2, void *nsd2, void *setp2,
	   int *np, int *&name_counts, char *&names,
	   int *&unit_counts, char *&units,
	   double *&low, double *&high);
  
  /** \brief Compute a point using the parameters given in \c vals
   */
  int compute_point(void *bcp2, void *mdp2, int nv, double *vals,
		    double *log_wgt);

  /** \brief Summarize tables
   */
  void summarize_tables(void *mdp2);
  
  /** \brief Get a column from the M-R table
   */
  int get_mvsr_column(void *mdp2, char *col_name, int &n, double *&ptr);

  /** \brief Get a column from the source table
   */
  int get_source_column(void *mdp2, char *col_name, int &n, double *&ptr);
  
  /** \brief Get a column from the reinterpolated table
   */
  int get_grid_column(void *mdp2, char *col_name, int &n, double *&ptr);

  /** \brief Get a column from the EOS table
   */
  int get_eos_column(void *mdp2, char *col_name, int &n, double *&ptr);

  /** \brief Get a constant from the M-R table
   */
  double get_mvsr_constant(void *mdp2, char *con_name);

  /** \brief Get a constant from the EOS table
   */
  double get_eos_constant(void *mdp2, char *con_name);

  /** \brief Free the memory associated with the \ref bamr_class
      and \ref model_data objects
   */
  void destroy_pointers(void *bcp2, void *mdp2);

  /** \brief Add a M-R data set
   */
  void add_data(void *nsd2, char *name, char *fname, char *slice,
		double mass_frac, char *table);

  /** \brief Add a M-R data set
   */
  void add_data_alt(void *nsd2, char *name, char *fname,
		    char *fname_alt, char *slice,
		    double mass_frac, char *table);
  
}

#endif
