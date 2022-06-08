/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2022, Joonas Nattila and Andrew W. Steiner
  
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
#ifndef FILTERS_H
#define FILTERS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>

#include <boost/numeric/ublas/matrix.hpp>

#include "fftw3.h"

#include <o2scl/table3d.h>

namespace filters {
  
  using ubmatrix = boost::numeric::ublas::matrix<double>;

  /** \brief General filtering using fftw for image convolution
      
      Because of the fftw3 usage this class is internally heavily
      relying on C code instead of C++.
  */
  class Filter {

  private:

    /// x image size
    int Nx;
    /// y image size
    int Ny;
    /// z image size
    int Nz;
    
    /// Circular indexing
    inline int circular(int x, int M) {
      if (x<0)    return x+M;
      if (x >= M) return x-M;
      return x;
    }

    /** \brief Internal circular indexing,
     
	NOTE: this is different from the rest of the code but done like this
	in order to appear consistent with the fttw examples.
    */
    inline int index(int i, int j) {
      i = circular(i, Nx);
      j = circular(j, Ny);

      assert(i >= 0 && i < Nx);
      assert(j >= 0 && j < Ny);

      return i + j*Nx;
    }

    /// Neighbor-padded arrays for convolution
    fftw_complex *arr;

    /// Actual (zero-padded) convolution kernel
    fftw_complex *kernel;

    /// Fftw3 transform plans
    //@{
    fftw_plan p_kernel, p_forw_arr, p_back_arr;
    //@}

  public:

    /// Build filter assuming information from 3x3x1 tiles
  Filter(int NxMesh, int NyMesh) : Nx( 3*NxMesh ), Ny( 3*NyMesh ), Nz( 1 ) {
      
      // allocate
      //arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx*Ny*Nz);
      arr = (fftw_complex*) fftw_alloc_complex(Nx*Ny*Nz);

      for(int j = 0 ; j < Ny ; ++j) {
	for(int i = 0 ; i < Nx ; ++i) {
	  arr[ index(i,j) ][0] = 0.0;
	  arr[ index(i,j) ][1] = 0.0;
	}
      }

      //kernel = (fftw_complex*) fftw_malloc(sizeof
      // (fftw_complex) * Nx*Ny*Nz);
	  
      kernel = (fftw_complex*) fftw_alloc_complex(Nx*Ny*Nz);

      // in-place complex to complex transform
      // TODO: change to real2complex and complex2real plans
      // TODO: use fftw guru interface to apply same plan
      // for all arrays (in same direction)
	  
      p_kernel    = fftw_plan_dft_2d
	( Nx,Ny, kernel, kernel, FFTW_FORWARD,  FFTW_MEASURE);

      p_forw_arr   = fftw_plan_dft_2d
	( Nx,Ny, arr, arr, FFTW_FORWARD,  FFTW_MEASURE );

      p_back_arr   = fftw_plan_dft_2d
	( Nx,Ny, arr, arr, FFTW_BACKWARD, FFTW_MEASURE );

      /*
	fftw_print_plan(p_forw_arr);

	fftw_print_plan(p_back_arr);
      */
    }

    /// Explicit destructor for fftw arrays and plans
    virtual ~Filter() {
      fftw_free(arr);
      fftw_free(kernel);

      fftw_destroy_plan(p_kernel);
      fftw_destroy_plan(p_forw_arr);
      fftw_destroy_plan(p_back_arr);
    }

    /// zero-padded centralized index
    //
    //NOTE:
    // The FFT is relative to the middle and in its scale there are
    // negative points. In the memory the points are 0...n-1, but the
    // FFT treats them as -ceil(n/2)...floor(n/2), where 0 is
    // -ceil(n/2) and n-1 is floor(n/2)
    //
    //NOTE: this is not how fftw3 does its indexing anymore!
    //
    //inline virtual std::tuple<int,int> zero_padded_index(int i, int j) 
    //{
    //  int iff = ceil(height/2);
    //  int jff = ceil(width /2);
    //  
    //  return std::make_tuple(iff, jff);
    //}

    /// Zero-wrapped index convention used by fftw3
    inline virtual std::tuple<int,int> zero_wrapped_index(int i, int j) {

      int iff = i >= 0 ? i : Nx + i;
      int jff = j >= 0 ? j : Ny + j;

      return std::make_tuple(iff, jff);
    }

    /** \brief Initialize kernel 

	In practice this is just a complex way to set kernel[0,0] = 1,
	but we do it like this for clarity as the syntax is used later on.
    */
    virtual void init_kernel() {

      double val;

      // kernel size (even number because of initialization)
      //int knx = height/3;
      //int kny = width /3;

      // halo region size
      int h1 = Nx/2; // division is floor(x/y) automatically
      int w1 = Ny/2;
      //int h2 = (height + 2 - 1)/2; // ceil(x/y)
      //int w2 = (width  + 2 - 1)/2;

      int wi,wj;

      for(int j=-w1; j<=w1 ; ++j) {
	for(int i =-h1; i<=h1 ; ++i) {
	  auto zindx = zero_wrapped_index(i,j);
	  wi = std::get<0>(zindx);
	  wj = std::get<1>(zindx);
          
	  val = 0.0;
	  if ((i ==  0) && (j ==  0)) val = 1.0;

	  //if ((i ==  1) && (j ==  0)) val = 1.0;
	  //if ((i ==  0) && (j ==  1)) val = 1.0;
	  //if ((i ==  1) && (j ==  1)) val = 1.0;
	  //if ((i == -1) && (j ==  0)) val = 1.0;
	  //if ((i ==  0) && (j == -1)) val = 1.0;
	  //if ((i == -1) && (j == -1)) val = 1.0;
	  //if ((i == 1 ) && (j == -1)) val = 1.0;
	  //if ((i ==-1 ) && (j ==  1)) val = 1.0;
        
	  kernel[ index(wi,wj) ][0] = val; // real part
	  kernel[ index(wi,wj) ][1] = 0.0; // complex part
	}
      }
    }

    /// Apply 3-point digital filter directly
    virtual void direct_convolve_3point() {

      // 3-point digital filter
      std::vector<double> coeffs = {{ 1., 2., 1.,
				      2., 4., 2.,
				      1., 2., 1. }};
      // normalize
      double norm = 0.0;
      for(auto& coeff : coeffs) norm += coeff;
      for(auto& coeff : coeffs) coeff /= norm;

      std::vector<double> image1, image2, image3;
      image1.resize(Nx*Ny*Nz); 

      for (int j=0; j < Ny;  ++j) {
	for (int i=0; i < Nx; ++i) {
	  image1[ index(i,j) ] = arr[ index(i,j) ][0];
	}
      }

      direct_convolve(image1.data(), coeffs.data(), sqrt( coeffs.size() ) );

      // normalize and copy back
      norm = 1.0;
      for (int j=0; j < Ny;  ++j) {
	for (int i=0; i < Nx; ++i) {
	  arr[ index(i,j) ][0] = image1[ index(i,j) ]/norm;
	}
      }

      return;
    }

    /** \brief initialize 3-point digital filter into kernel array for FFT
	transformations
    */
    virtual void init_3point_kernel(int times) {

      // kernel
      //int K = 3; // three point kernel

      // 3-point digital filter
      std::vector<double> coeffs = {{ 1., 2., 1.,
				      2., 4., 2.,
				      1., 2., 1. }};

      double norm = 0.0;
      for(auto& coeff : coeffs) norm += coeff;
      for(auto& coeff : coeffs) coeff /= norm;


      // create temporary Real number image array
      // NOTE: can not easily copy kernel into pure real part due to
      // interleaved nature
      std::vector<double> image;
      image.resize(Nx*Ny*Nz); 

      for (int j=0; j < Ny; ++j) 
	for (int i=0; i < Nx; ++i)
	  image[ index(i,j) ] = kernel[ index(i,j) ][0];

      // perform convolution N times
      for(int N=0; N < times; N++) {
	direct_convolve(image.data(), coeffs.data(), sqrt(coeffs.size()) );
      }
    
      // normalize
      // norm = 0.0;
      //for (int i=0; i < height; ++i)
      //  for (int j=0; j < width;  ++j) 
      //    norm += image[ index(i,j) ];
      norm = 1.0;

      // normalize and copy back
      for (int j=0; j < Ny;  ++j) {
	for (int i=0; i < Nx; ++i) {
	  kernel[ index(i,j) ][0] = image[ index(i,j) ]/norm;
	}
      }

      return;
    }

    /** \brief Direct circular convolve two arrays
	
	NOTE: assumes height x width for image size
	
	in, out are m x n images (integer data) K is the kernel size
	(KxK) - currently needs to be an odd number, e.g. 3 coeffs[K][K]
	is a 2D array of integer coefficients scale is a scaling factor
	to normalise the filter gain
    */
    void direct_convolve(double* image, double* kernel_loc, int K) {

      // out array
      double data;
      std::vector<double> out;
      out.resize(Nx*Ny*Nz);

      //for (int j = K/2; j < width -K/2; ++j) // iterate through image
      // iterate through circular image
      for (int j=0; j < Ny; ++j) {

	//for (int i = K/2; i < height - K/2; ++i) // iterate through image
	// iterate through circular image
	for (int i=0; i < Nx; ++i) {
	  
	  // sum will be the sum of input data * coeff terms
	  double sum = 0.0; 
	  
	  // convolution of single point
	  for (int jj = -K/2; jj <= K/2; ++jj) {
	    // iterate over kernel
	    for (int ii = -K/2; ii <= K/2; ++ii) {
	      data = image[ index(i+ii, j+jj) ]; 
	      double coeff = kernel_loc[ (ii + K/2)*K + (jj + K/2) ];
	      
	      sum += data * coeff;
	    }
	  }
	  // sum of convolution products and store in output
	  out[ index(i,j) ] = sum; 
	}
	// end of conv
      } 
      
      // copy (real part) back
      for (int j=0; j < Ny;  ++j) {
	for (int i=0; i < Nx; ++i) {
	  image[ index(i,j) ] = out[ index(i,j) ];
	}
      }

      return;
    } 

    /// Initialize Gaussian kernel given sigma x and sigma y
    virtual void init_gaussian_kernel(double sigx, double sigy) {

      // kernel size (even number because of initialization)
      //int knx = height/3;
      //int kny = width /3;

      // halo region size
      int h1 = Nx/2; // division is floor(x/y) automatically
      int w1 = Ny/2;
      int h2 = (Nx + 2 - 1)/2; // ceil(x/y)
      int w2 = (Ny + 2 - 1)/2;

      int wi,wj;
      double val;

      double sum = 0.0;
      for(int j=-w1; j<w2; ++j) {
	for(int i =-h1; i<h2; ++i) {
	  auto zindx = zero_wrapped_index(i,j);
	  wi = std::get<0>(zindx);
	  wj = std::get<1>(zindx);

	  assert(wi >= 0 && wi < Nx);
	  assert(wj >= 0 && wj < Ny);

	  val = 1.0;
	  val *= exp( -0.5*((double)(i*i))/sigx/sigx);
	  val *= exp( -0.5*((double)(j*j))/sigy/sigy);

	  kernel[ index(wi,wj) ][0] = val; // real part
	  kernel[ index(wi,wj) ][1] = 0.0; // complex part
	  sum += val;
	}
      }

      // Normalize
      
      if (false) {
	for(int j=-w1; j<w2; ++j) {
	  for(int i =-h1; i<h2; ++i) {
	    auto zindx = zero_wrapped_index(i,j);
	    wi = std::get<0>(zindx);
	    wj = std::get<1>(zindx);
	    kernel[ index(wi,wj) ][0] /= sum;
	  }
	}
      }

      // AWS: I added this on 9/17/18, and now as of 10/3/18 I
      // don't think it's necessary, but it doesn't harm 
      // anything so I'm leaving it in for now.
      
      if (true) {
	for(int j = 0 ; j < Ny ; ++j) {
	  for(int i = 0 ; i < Nx ; ++i) {
	    arr[ index(i,j) ][0] = 0.0;
	    arr[ index(i,j) ][1] = 0.0;
	  }
	}
      }
      
      return;
    }

    /// initialize 2d box sinc filter (in frequency space)
    /*
      virtual void init_sinc_kernel(double X, double Y) 
      {
      double val, u,v;
      int h2 = floor(height/2);
      int w2 = floor(width /2);
      int wi, wj;

      double sum = 0.0;
      for(int i =-h2; i<=h2 ; ++i) {
      for(int j=-w2; j<=w2 ; ++j) {
      auto zindx = zero_wrapped_index(i,j);
      wi = std::get<0>(zindx);
      wj = std::get<1>(zindx);

      u = (double)i;
      v = (double)j;

      val = 0.0;
      if ((abs(i) < height/3) || (abs(j) < width/3)){
      val = X*sin(pi*X*u)/(pi*X*u) * 
      Y*sin(pi*Y*v)/(pi*Y*v);
      }
      //val = (X*sin(pi*X*u)/(pi*X*u)) * (Y*sin(pi*Y*v)/(pi*Y*v));
      if ((i == 0) && (j == 0)) val = 4.0*pi*pi*X*Y; // center to avoid 0/0
      else if (i == 0) val = 2.0*pi*X*Y*sin(pi*Y*v)/(pi*Y*v);
      else if (j == 0) val = 2.0*pi*X*Y*sin(pi*X*u)/(pi*X*u);

      // windowing
      val *= 0.42 - 0.5*cos(2*pi*3.0/height) + 0.08*cos(4.0*pi*3.0/height);
      val *= 0.42 - 0.5*cos(2*pi*3.0/width ) + 0.08*cos(4.0*pi*3.0/width);

      kernel[ index(wi,wj) ][0] = val; // real part
      kernel[ index(wi,wj) ][1] = 0.0; // complex part
      sum += val;
      }
      }

      // normalize 
      for(int i =0; i<height ; ++i)
      for(int j=0; j<width ; ++j)
      kernel[ index(i,j) ][0] /= sum;

      }
    */

    /** \brief Low-pass filter in frequency space
	
	uses cutoff to describe how many array elements are filtered
    */
    virtual void init_lowpass_fft_kernel(int cutoff) {

      // kernel size (even number because of initialization)
      //int knx = height/3;
      //int kny = width /3;

      // halo region size
      int h1 = Nx/2; // division is floor(x/y) automatically
      int w1 = Ny/2;
      int h2 = (Nx + 2 - 1)/2; // ceil(x/y)
      int w2 = (Ny + 2 - 1)/2;

      int wi,wj;
      double val;

      double sum = 0.0;
      for(int i =-h1; i<h2; ++i) {
	for(int j=-w1; j<w2; ++j) {
	  auto zindx = zero_wrapped_index(i,j);
	  wi = std::get<0>(zindx);
	  wj = std::get<1>(zindx);

	  assert(wi >= 0 && wi < Nx);
	  assert(wj >= 0 && wj < Ny);

	  val = 0.0;
	  // circle (Bessel)
	  if (sqrt(i*i + j*j) < cutoff) val = 1.0; 

	  // box (2D Sinc)
	  //if ((sqrt(i*i) < cutoff) && (sqrt(j*j) < cutoff)) val = 1.0; 

	  kernel[ index(wi,wj) ][0] = val; // real part
	  kernel[ index(wi,wj) ][1] = 0.0; // complex part
	  sum += val;
	}
      }
      return;
    }

    /// Normalize fft transformation
    void normalize() {

      for(int k = 0 ; k < Nz ; ++k) {
	for(int j = 0 ; j < Ny ; ++j) {
	  for(int i  = 0 ; i < Nx ; ++i) {
	    arr[ index(i,j) ][0] /= Nx*Ny*Nz;
	  }
	}
      }
      return;
    }

    /// FFT kernel (once is enough)
    virtual void fft_kernel() {
      fftw_execute(p_kernel);
    }

    /// FFT image forward
    virtual void fft_image_forward() {
      fftw_execute(p_forw_arr);
      return;
    }

    /// FFT image backwards and normalize 
    virtual void fft_image_backward() { 
      fftw_execute(p_back_arr);
      normalize();
      return;
    }

    /// Multiply kernel and image
    void apply_kernel() {

      double x1, y1, x2, y2, x3, y3;
      double u, v;
      for(int j = 0 ; j < Ny ; ++j) {
	for(int i  = 0 ; i < Nx ; ++i) {

	  x1 = arr[ index(i,j) ][0];
	  y1 = arr[ index(i,j) ][1];

	  u = kernel[ index(i,j) ][0];
	  v = kernel[ index(i,j) ][1];

	  arr[ index(i,j) ][0] = x1*u - y1*v;
	  arr[ index(i,j) ][1] = x1*v + y1*u;

	}
      }

      return;
    }

    /// Set image from \ref o2scl::table3d object
    void set_image(const o2scl::table3d &in,
		   const std::string &slice_name) {

      const auto& data = in.get_slice(slice_name);
      int Nx_loc = in.get_nx();
      int Ny_loc = in.get_ny();

      int indx;
      int s = 0; // z-dimension; collapsed

      //for(int s=0; r<Nz; s++) {
      for(int r=0; r<Ny_loc; r++) {
	for(int q=0; q<Nx_loc; q++) {
          
	  indx = index(q, r);
	  arr[ indx ][0] = data(q,r);
	}
      }
      return;
    }


    /// Get image back from filter class
    void get_image(o2scl::table3d& out,
		   const std::string& slice_name) {

      auto& data = out.get_slice(slice_name);
      int Nx_loc = out.get_nx();
      int Ny_loc = out.get_ny();

      int indx;
      int s = 0; // z-dimension; collapsed

      //for(int s=0; r<Nz; s++) {
      for(int r=0; r<Ny_loc; r++) {
	for(int q=0; q<Nx_loc; q++) {
          
	  indx = index(q, r);
	  data(q,r) = arr[ indx ][0];
	}
      }
      return;
    }

    /// \name Utility functions for debugging
    //@{
    void set_kernel(std::vector<double>& in) {
      if(in.size() != (size_t)Nx*Ny*Nz) {
	std::cout << "error in size!\n";
      }

      int q = 0;
      for(int j = 0 ; j < Ny ; ++j) {
	for(int i = 0 ; i < Nx ; ++i, q++) {
	  kernel[ index(i,j) ][0] = in[q];
	}
      }
      return;
    }

    /** \brief Desc
     */
    void set_image(std::vector<double>& in) {

      if(in.size() != (size_t)Nx*Ny*Nz) {
	std::cout << "error in size!\n";
      }

      int q = 0;
      for(int j = 0 ; j < Ny ; ++j) {
	for(int i = 0 ; i < Nx ; ++i, q++) {
	  arr[ index(i,j) ][0] = in[q];
	}
      }
      
      return;
    }

    /** \brief Desc
     */
    std::vector<double> get_kernel() {
      std::vector<double> ret;

      for(int j = 0 ; j < Ny ; ++j) {
	for(int i = 0 ; i < Nx ; ++i) {  
	  ret.push_back( kernel[ index(i,j) ][0] );
	}
      }

      return ret;
    }

    /** \brief Desc
     */
    std::vector<double> get_image() {
      std::vector<double> ret;

      for(int j = 0 ; j < Ny ; ++j) {
	for(int i = 0 ; i < Nx ; ++i) {  
	  ret.push_back( arr[ index(i,j) ][0] );
	}
      }

      return ret;
    }
    //@}

  };

  // end of namespace filters
} 

#endif
