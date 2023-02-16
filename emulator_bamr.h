/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
  and Andrew W. Steiner
  
  This file is part of Bamr.
  
  Bamr is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published 
  by the Free Software Foundation; either version 3 of the License, 
  or (at your option) any later version.
  
  Bamr is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Bamr. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef EMULATOR_BAMR_H
#define EMULATOR_BAMR_H

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <Python.h>
#include <o2scl/hdf_file.h>
#include <o2scl/mcmc_para.h>
#include "emulator.h"
#include "bamr_class.h"

#ifdef BAMR_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

#ifdef BAMR_MPI
#include <mpi.h>
#endif

namespace bamr {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef std::function<int(size_t, const ubvector &, double &,
                            model_data &)> point_funct;
  typedef std::function<int(const ubvector &, double,
                            std::vector<double> &,
                            model_data &)> fill_funct;
  
  /** \brief A specialized emulator for this code
      
      Note that this emulator has a data type of array<double,ndat>,
      because that is what the MCMC is using, but internally em1 uses a
      data type of vector<double>, since not all output fields are
      always emulated.
  */
  class emulator_bamr : public o2scl::emulator_unc<model_data,
                                      model_data, ubvector> {  
  public:
    
    int mpi_rank;
    int mpi_size;
    
    /** \brief The internal generic emulator 
    */
    o2scl::emulator_interpm_idw_table<std::vector<double>,ubvector> em1;
    
    /** \brief The second internal generic emulator 
    */
    o2scl::emulator_interpm_krige_table<std::vector<double>,ubvector> em2;
    
    o2scl::emulator_python<std::vector<double>,ubvector> em3;
    
    /** \brief Pointer to the data_eval object for function evaluations 
    */
    bamr_class *bcp;

    std::ofstream *sop;
    
    /** \brief List of parameters and output quantities (including log_wgt) 
    */
    std::vector<std::string> list;
    
    /** \brief Tables from which to train 
    */
    o2scl::table_units<> table;
    
    /** \brief List of filenames containing training tables 
    */
    std::vector<std::string> files;
    
    /** \brief Number of parameters 
    */
    size_t np;
    
    /** \brief Number of output quantities for em1 
    */
    size_t nout;
    
    /** \brief Mapping of parameters to indices 
    */
    o2scl::vec_index pvi;
    
    /** \brief Mapping of output data to array indices 
    */
    o2scl::vec_index dvi;
    
    /** \brief Random number generator 
    */
    o2scl::rng<> r;
    
    emulator_bamr();
    
    /** \brief Train the emulator with the data in \c tab_train 
    */
    void train(o2scl::table_units<> &tab_train, o2scl::vec_index &pvii,
               o2scl::vec_index &dvii, bamr_class *bcpi,
               std::ofstream *sopi);
    
    /** \brief Evaluate the emulator, and or the full function if necessary 
    */
    virtual int eval(size_t n, const ubvector &p, double &log_wgt,
                     model_data &dat);
    
    /** \brief Evaluate the emulator, and or the full function if necessary 
    */
    virtual int eval_unc(size_t n, const ubvector &p, double &log_wgt,
                         double &lw_unc, model_data &dat,
                         model_data &dat_unc);
  };
}
#endif