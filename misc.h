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
#ifndef MISC_H
#define MISC_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/cold_nstar.h>
#include <o2scl/hist.h>
#include <o2scl/hist_2d.h>
#include <o2scl/expect_val.h>

#ifndef DOXYGEN
namespace o2scl {
#endif

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

  /** \brief A simplified version of cold_nstar 
      
      This simplified version only computes the energy density and 
      pressure rather than the \o2 version which computes several
      extra quantities. This class is part of \ref model.
  */
  class cold_nstar2 : public cold_nstar {

  public:

    /// Compute the core EOS 
    int calc_eos(double &n1, double &e1, double np_0=0.0);

  };

#ifndef DOXYGEN
}
#endif

#endif
