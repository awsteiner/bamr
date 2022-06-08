/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Andrew W. Steiner
  
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
/** \file nstar_cold2.h
    \brief Definition of nstar_cold2
*/
#ifndef MISC_H
#define MISC_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/nstar_cold.h>

namespace bamr {

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

  /** \brief A simplified version of nstar_cold 
      
      This simplified version only computes the energy density and 
      pressure rather than the \o2 version which computes several
      extra quantities. This class is part of \ref bamr::model.
  */
  class nstar_cold2 : public o2scl::nstar_cold {

  public:

    /// Compute the core EOS 
    int calc_eos(double np_0=0.0);

  };

}

#endif
