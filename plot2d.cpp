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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <TEllipse.h>

#include <o2scl/table.h>
#include <o2scl/table3d.h>
#include <o2scl/format_float.h>
#include <o2scl/hdf_file.h>
#include <o2scl/vec_stats.h>
#include <o2scl/contour.h>
#include <o2scl/expect_val.h>
#include <o2scl/hist_2d.h>
#include <o2scl/graph.h>
#include <o2scl/interp.h>
#include <o2scl/cli_readline.h>

// For hist_set_ev
#include "misc.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_graph;

typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

/** \brief Define a new histogram object which accesses the
    weights through <tt>operator()</tt>
*/
class hist_2d_matrix : public hist_2d {

public:

  // Return a reference to weight at <tt>(i,j)</tt>
  double &operator()(size_t i, size_t j) {
    return get_wgt_i(i,j);
  }

};

/** \brief Set of one-dimensional histograms over a grid in x
 */
class hist_set_ev_x {

public:

  /// Histogram
  hist_2d_matrix h;
  /// Count of measurements for each histogram in the set
  ubvector_size_t count;
  /// Expectation value object
  matrix_ev ev;

  /// Record hist in ev object
  void update() {
    ev.add(h);
    h.clear_wgts();
    //count.set_all(0);
    return;
  }

};

/** \brief Desc
 */
class plot2d_class {

public:

  /// Desc
  double xscale;
  /// Desc
  double one_sigma;
  /// Desc
  double two_sigma;
  /// Desc
  double three_sigma;
  /// Desc
  double yscale;
  /// Desc
  bool logx;
  /// Desc
  bool logy;
  /// Desc
  string out_file;
  /// Desc
  string xlabel;
  /// Desc
  string ylabel;
  /// Desc
  string title;
  /// Desc
  std::vector<double> ellvec;
  /// Desc
  double min_total;

  /// \name Desc
  //@{
  static const int no_contours=0;
  static const int x_contours=1;
  static const int y_contours=2;
  static const int xy_contours=3;
  int contours;
  //@}

  /// \name Desc
  //@{
  static const int norm_none=0;
  static const int norm_col_unity=1;
  static const int norm_col_weight=2;
  static const int norm_row_weight=3;
  static const int norm_max_unity=4;
  int normalize;
  //@}

  plot2d_class() {
    xscale=1.0;
    yscale=1.0;
    logx=false;
    logy=false;
    contours=x_contours;
    min_total=0.0;
    normalize=norm_col_weight;
    xlabel="x-axis Label";
    ylabel="y-axis Label";
    title="title";

    one_sigma=gsl_sf_erf(1.0/sqrt(2.0));
    two_sigma=gsl_sf_erf(2.0/sqrt(2.0));
    three_sigma=gsl_sf_erf(3.0/sqrt(2.0));
  }

  /** \brief Desc
   */
  int contour_fun(std::vector<std::string> &sv, bool itive_com) {
    if (sv[1]=="none") {
      contours=no_contours;
    } else if (sv[1]=="x") {
      contours=x_contours;
    } else if (sv[1]=="y") {
      contours=y_contours;
    } else if (sv[1]=="xy") {
      contours=xy_contours;
    } else {
      cout << "Incorrect argument to contours." << endl;
      return gsl_efailed;
    }
    return 0;
  }

  /** \brief Desc
   */
  int normal_fun(std::vector<std::string> &sv, bool itive_com) {
    if (sv[1]=="none") {
      normalize=norm_none;
    } else if (sv[1]=="max") {
      normalize=norm_max_unity;
    } else if (sv[1]=="colu") {
      normalize=norm_col_unity;
    } else if (sv[1]=="col") {
      normalize=norm_col_weight;
    } else if (sv[1]=="row") {
      normalize=norm_row_weight;
    } else {
      cout << "Incorrect argument to normalize." << endl;
      return gsl_efailed;
    }
    return 0;
  }

  /** \brief Desc
   */
  int ell(std::vector<std::string> &sv, bool itive_com) {
    ellvec.push_back(stod(sv[1]));
    ellvec.push_back(stod(sv[2]));
    ellvec.push_back(stod(sv[3]));
    ellvec.push_back(stod(sv[4]));
    return 0;
  }

  /** \brief Desc
   */
  int show(std::vector<std::string> &sv, bool itive_com) {
    
    // object name
    string name=sv[1];
    
    // file list
    string file=sv[2];

    int ix=stoi(sv[3]);
    int jx=stoi(sv[4]);

    // -----------------------------------------------------------------
    // Read data file

    hdf_file hf;
    hist_set_ev_x data;
    size_t nh=0;

    cout << "Opening file named: " << file << endl;
    hf.open(file);
    
    cout << "Reading histogram named: " << name+".h" << endl;
    hdf_input(hf,data.h,name+".h");
    nh=data.h.size_x();
    cout << "Reading vector_ev object named: " << name+".ev" << endl;
    hdf_input(hf,data.ev,name+".ev");
    
    cout << "Closing HDF file." << endl;
    hf.close();
    
    matrix_ev &mev=data.ev;
    hist_2d_matrix &h=data.h;
    const tensor3 &mevt=mev.get_data();

    double total=0.0;
    for(size_t i=0;i<nh;i++) {
      for(size_t j=0;j<nh;j++) {
	total+=h.get_wgt_i(i,j);
      }
    }
    cout << "Histogram total: " << total << endl;

    size_t n_blocks, n_per_block;
    size_t i_block, i_current_block;
    mev.get_blocks(n_blocks,n_per_block);
    mev.get_block_indices(i_block,i_current_block);
    cout << "n_blocks, n_per_block: " << n_blocks << " " 
	 << n_per_block << endl;
    cout << "i_block, i_current_block: " << i_block << " " 
	 << i_current_block << endl;
    
    cout << "x: " << ix << " " << h.get_x_rep_i(ix) << "\t";
    cout << "y: " << jx << " " << h.get_y_rep_i(jx) << endl;
    for(size_t i=0;i<20;i++) {
      cout.width(2);
      cout << i << " ";
      if (mevt.get(((size_t)ix),((size_t)jx),i)>=0.0) {
	cout << " " << mevt.get(((size_t)ix),((size_t)jx),i) << endl;
      } else {
	cout << mevt.get(((size_t)ix),((size_t)jx),i) << endl;
      }
    }

    return 0;
  }

  /** \brief Desc
   */
  int row(std::vector<std::string> &sv, bool itive_com) {
    
    // object name
    string name=sv[1];
    
    // file list
    string file=sv[2];

    int ix=stoi(sv[3]);

    // -----------------------------------------------------------------
    // Read data file

    hdf_file hf;
    hist_set_ev_x data;
    size_t nh=0;

    cout << "Opening file: " << file << endl;
    hf.open(file);
    
    cout << "Reading hist: " << name+".h" << endl;
    hdf_input(hf,data.h,name+".h");
    nh=data.h.size_x();
    cout << "Reading vector_ev: " << name+".ev" << endl;
    hdf_input(hf,data.ev,name+".ev");
    
    cout << "Closing HDF file: " << endl;
    hf.close();

    matrix_ev &mev=data.ev;
    hist_2d_matrix &h=data.h;
    const tensor3 &mevt=mev.get_data();
    
    size_t n_blocks, n_per_block;
    size_t i_block, i_current_block;
    mev.get_blocks(n_blocks,n_per_block);
    mev.get_block_indices(i_block,i_current_block);
    cout << "n_blocks, n_per_block: " << n_blocks << " " 
	 << n_per_block << endl;
    cout << "i_block, i_current_block: " << i_block << " " 
	 << i_current_block << endl;
    
    cout << "Row: " << ix << " " << h.get_x_rep_i(ix) << endl;
    
    ubvector total(nh);
    for(size_t jx=0;jx<nh;jx++) {
      total[jx]=0.0;
      for(size_t i=0;i<i_block;i++) {
	total[jx]+=mevt.get(((size_t)ix),((size_t)jx),i);
      }
    }
    cout.setf(ios::showpos);
    for(size_t jx=0;jx<nh;jx++) {
      if (total[jx]>0.0) {
	cout << h.get_y_rep_i(jx) << " ";
      }
    }
    cout << endl;
    for(size_t i=0;i<i_block;i++) {
      for(size_t jx=0;jx<nh;jx++) {
	if (total[jx]>0.0) {
	  cout << mevt.get(((size_t)ix),((size_t)jx),i) << " ";
	}
      }
      cout << endl;
    }
    cout.unsetf(ios::showpos);

    return 0;
  }

  /** \brief Desc
   */
  int tabulate(std::vector<std::string> &sv, bool itive_com) {
    
    // object name
    string name=sv[1];
    
    // file list
    vector<string> files;

    // Form list of data files
    for(size_t i=2;i<sv.size();i++) files.push_back(sv[i]);
    size_t nf=files.size();

    // data
    vector<hist_set_ev_x> data;
    data.resize(nf);
    
    // Averages
    vector<ubmatrix> avg, std_dev, avg_err;
    avg.resize(nf);
    std_dev.resize(nf);
    avg_err.resize(nf);

    // -----------------------------------------------------------------
    // Read data files

    hdf_file hf;
    size_t nh=0;

    for(size_t i=0;i<nf;i++) {
      cout << "Opening file named: " << files[i] << endl;
      hf.open(files[i]);
      
      cout << "Reading histogram." << endl;
      hdf_input(hf,data[i].h,name+".h");
      nh=data[0].h.size_x();
      cout << "Histogram size: " << nh << endl;

      double total=0.0;
      for(size_t ii=0;ii<nh;ii++) {
	for(size_t j=0;j<nh;j++) {
	  total+=data[0].h.get_wgt_i(ii,j);
	}
      }
      cout << "Histogram total: " << total << endl;

      cout << "Reading vector_ev object." << endl;
      hdf_input(hf,data[i].ev,name+".ev");

      cout << "Closing HDF file." << endl;
      hf.close();

      if (true) {
	matrix_ev &mev=data[i].ev;
	const tensor3 &mev_data=mev.get_data();
	size_t n_blocks, n_per_block;
	size_t i_block, i_current_block;
	mev.get_blocks(n_blocks,n_per_block);
	mev.get_block_indices(i_block,i_current_block);
	cout << "n_blocks, n_per_block: " << n_blocks << " " 
	     << n_per_block << endl;
	cout << "i_block, i_current_block: " << i_block << " " 
	     << i_current_block << endl;
	for(size_t ib=0;ib<n_blocks;ib++) {
	  double total=0.0;
	  for(size_t j=0;j<nh;j++) {
	    for(size_t k=0;k<nh;k++) {
	      total+=mev_data.get(j,k,ib);
	    }
	  }
	  cout << "Total for block " << ib << ": " << total << endl;
	}
      }
	

      cout << "Computing file stats." << endl;
      avg[i].resize(nh,nh);
      std_dev[i].resize(nh,nh);
      avg_err[i].resize(nh,nh);
      data[i].ev.current_avg(avg[i],std_dev[i],avg_err[i]);
      //cout << "Min and max of average for this file: " 
      //<< *std::min_element(avg[i].begin(),avg[i].end()) << " " 
      //<< *std::max_element(avg[i].begin(),avg[i].end()) << endl;
    }

    // -----------------------------------------------------------------
    // Compute weighted average over all files

    cout << "Averaging over all files." << endl;

    ubmatrix tavg(nh,nh);

    for(size_t i=0;i<nh;i++) {
      for(size_t j=0;j<nh;j++) {
	// Construct a vector of the averages for each file
	std::vector<double> xa, xe;
	for(size_t k=0;k<nf;k++) {
	  xa.push_back(avg[k](i,j));
	  xe.push_back(avg_err[k](i,j));
	}
	// If there isn't information on errors, just use the
	// unweighted mean
	if (*std::max_element(xe.begin(),xe.end())==0.0) {
	  tavg(i,j)=vector_mean(nf,xa);
	} else {
	  // Weighted average
	  for(size_t k=0;k<nf;k++) {
	    // There's no error information for this file, so just
	    // use the maximum error over all files 
	    if (xe[k]==0.0) xe[k]=*max_element(xe.begin(),xe.end());
	    xe[k]=fabs(1.0/xe[k]/xe[k]);
	  }
	  tavg(i,j)=wvector_mean(nf,xa,xe);
	}
      }
    }

    // -----------------------------------------------------------------
    // Normalize the full histogram to ensure the maximum value is one

    if (normalize==norm_max_unity) {
      double max=0.0;
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  if (tavg(i,j)>max) max=tavg(i,j);
	}
      }
      if (max>0.0) {
	for(size_t i=0;i<nh;i++) {
	  for(size_t j=0;j<nh;j++) {
	    tavg(i,j)/=max;
	  }
	}
      }
    } 

    // -----------------------------------------------------------------
    // Normalize each column to ensure the maximum value is 1

    if (normalize==norm_col_unity) {
      for(size_t j=0;j<nh;j++) {
	double max=0.0;
	for(size_t i=0;i<nh;i++) {
	  if (tavg(i,j)>max) max=tavg(i,j);
	}
	if (max>0.0) {
	  for(size_t i=0;i<nh;i++) {
	    tavg(i,j)/=max;
	  }
	}
      }
    } 

    // -----------------------------------------------------------------
    // Normalize each column to ensure it's maximum value is equal to
    // the total number of measurements in that column

    if (normalize==norm_col_weight) {
      // Total for each column
      ubvector total(nh);
      // Max in each column
      ubvector max(nh);
      fill(total.begin(),total.end(),0.0);
      fill(max.begin(),max.end(),0.0);
      // Compute total and max in each column
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  total[j]+=tavg(i,j);
	  if (tavg(i,j)>max[j]) max[j]=tavg(i,j);
	}
      }
      // Do normalization
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  if (*max_element(total.begin(),total.end())>0.0 && max[j]>0.0) {
	    tavg(i,j)*=total[j]/max[j]/
	      (*max_element(total.begin(),total.end()));
	  }
	}
      }
    }

    // -----------------------------------------------------------------
    // Normalize each row to ensure it's maximum value is equal to
    // the total number of measurements in that row

    if (normalize==norm_row_weight) {
      // Total for each row
      ubvector total(nh);
      // Max in each row
      ubvector max(nh);
      fill(total.begin(),total.end(),0.0);
      fill(max.begin(),max.end(),0.0);
      // Compute total and max in each row
      for(size_t i=0;i<nh;i++) {
	for(size_t j=0;j<nh;j++) {
	  total[i]+=tavg(i,j);
	  if (tavg(i,j)>max[i]) max[i]=tavg(i,j);
	}
      }
      // Do normalization
      for(size_t i=0;i<nh;i++) {
	for(size_t j=0;j<nh;j++) {
	  if ((*max_element(total.begin(),total.end()))>0.0 && max[i]>0.0) {
	    tavg(i,j)*=total[i]/max[i]/
	      (*max_element(total.begin(),total.end()));
	  }
	}
      }
    }

    std::vector<double> xg, yg;
    for(size_t i=0;i<nh;i++) {
      xg.push_back(data[0].h.get_x_rep_i(i));
      yg.push_back(data[0].h.get_y_rep_i(i));
    }
    table3d tnew;
    tnew.set_xy("x",nh,xg,"y",nh,yg);
    tnew.new_slice("weight");
    for(size_t i=0;i<nh;i++) {
      for(size_t j=0;j<nh;j++) {
	tnew.set(i,j,"weight",tavg(i,j));
      }
    }
    
    if (out_file.length()>0) {
      hf.open_or_create(out_file);
    } else {
      hf.open_or_create("tab.o2");
    }
    hdf_output(hf,tnew,"tab");
    hf.close();

    return 0;
  }

  /** \brief Desc
   */
  int plot(std::vector<std::string> &sv, bool itive_com) {
    
    // object name
    string name=sv[1];
    
    // file list
    vector<string> files;

    // Form list of data files
    for(size_t i=2;i<sv.size();i++) files.push_back(sv[i]);
    size_t nf=files.size();

    // data
    vector<hist_set_ev_x> data;
    data.resize(nf);
    
    // Averages
    vector<ubmatrix> avg, std_dev, avg_err;
    avg.resize(nf);
    std_dev.resize(nf);
    avg_err.resize(nf);

    // -----------------------------------------------------------------
    // Read data files

    hdf_file hf;
    size_t nh=0;

    for(size_t i=0;i<nf;i++) {
      cout << "Opening file named: " << files[i] << endl;
      hf.open(files[i]);
      
      cout << "Reading histogram." << endl;
      hdf_input(hf,data[i].h,name+".h");
      nh=data[0].h.size_x();
      cout << "Histogram size: " << nh << endl;

      double total=0.0;
      for(size_t ii=0;ii<nh;ii++) {
	for(size_t j=0;j<nh;j++) {
	  total+=data[0].h.get_wgt_i(ii,j);
	}
      }
      cout << "Histogram total: " << total << endl;

      cout << "Reading vector_ev object." << endl;
      hdf_input(hf,data[i].ev,name+".ev");

      cout << "Closing HDF file." << endl;
      hf.close();

      if (true) {
	matrix_ev &mev=data[i].ev;
	const tensor3 &mev_data=mev.get_data();
	size_t n_blocks, n_per_block;
	size_t i_block, i_current_block;
	mev.get_blocks(n_blocks,n_per_block);
	mev.get_block_indices(i_block,i_current_block);
	cout << "n_blocks, n_per_block: " << n_blocks << " " 
	     << n_per_block << endl;
	cout << "i_block, i_current_block: " << i_block << " " 
	     << i_current_block << endl;
	for(size_t ib=0;ib<n_blocks;ib++) {
	  double total=0.0;
	  for(size_t j=0;j<nh;j++) {
	    for(size_t k=0;k<nh;k++) {
	      total+=mev_data.get(j,k,ib);
	    }
	  }
	  cout << "Total for block " << ib << ": " << total << endl;
	}
      }
	

      cout << "Computing file stats." << endl;
      avg[i].resize(nh,nh);
      std_dev[i].resize(nh,nh);
      avg_err[i].resize(nh,nh);
      data[i].ev.current_avg(avg[i],std_dev[i],avg_err[i]);
      //cout << "Min and max of average for this file: " 
      //<< *std::min_element(avg[i].begin(),avg[i].end()) << " " 
      //<< *std::max_element(avg[i].begin(),avg[i].end()) << endl;
    }

    // -----------------------------------------------------------------
    // Compute weighted average over all files

    cout << "Averaging over all files." << endl;

    ubmatrix tavg(nh,nh);

    for(size_t i=0;i<nh;i++) {
      for(size_t j=0;j<nh;j++) {
	// Construct a vector of the averages for each file
	std::vector<double> xa, xe;
	for(size_t k=0;k<nf;k++) {
	  xa.push_back(avg[k](i,j));
	  xe.push_back(avg_err[k](i,j));
	}
	// If there isn't information on errors, just use the
	// unweighted mean
	if (*std::max_element(xe.begin(),xe.end())==0.0) {
	  tavg(i,j)=vector_mean(nf,xa);
	} else {
	  // Weighted average
	  for(size_t k=0;k<nf;k++) {
	    // There's no error information for this file, so just
	    // use the maximum error over all files 
	    if (xe[k]==0.0) xe[k]=*max_element(xe.begin(),xe.end());
	    xe[k]=fabs(1.0/xe[k]/xe[k]);
	  }
	  tavg(i,j)=wvector_mean(nf,xa,xe);
	}
      }
    }

    // -----------------------------------------------------------------
    // Normalize the full histogram to ensure the maximum value is one

    if (normalize==norm_max_unity) {
      double max=0.0;
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  if (tavg(i,j)>max) max=tavg(i,j);
	}
      }
      if (max>0.0) {
	for(size_t i=0;i<nh;i++) {
	  for(size_t j=0;j<nh;j++) {
	    tavg(i,j)/=max;
	  }
	}
      }
    } 

    // -----------------------------------------------------------------
    // Normalize each column to ensure the maximum value is 1

    if (normalize==norm_col_unity) {
      for(size_t j=0;j<nh;j++) {
	double max=0.0;
	for(size_t i=0;i<nh;i++) {
	  if (tavg(i,j)>max) max=tavg(i,j);
	}
	if (max>0.0) {
	  for(size_t i=0;i<nh;i++) {
	    tavg(i,j)/=max;
	  }
	}
      }
    } 

    // -----------------------------------------------------------------
    // Normalize each column to ensure it's maximum value is equal to
    // the total number of measurements in that column

    if (normalize==norm_col_weight) {
      // Total for each column
      ubvector total(nh);
      // Max in each column
      ubvector max(nh);
      fill(total.begin(),total.end(),0.0);
      fill(max.begin(),max.end(),0.0);
      // Compute total and max in each column
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  total[j]+=tavg(i,j);
	  if (tavg(i,j)>max[j]) max[j]=tavg(i,j);
	}
      }
      // Do normalization
      for(size_t j=0;j<nh;j++) {
	for(size_t i=0;i<nh;i++) {
	  if (*max_element(total.begin(),total.end())>0.0 && max[j]>0.0) {
	    tavg(i,j)*=total[j]/max[j]/
	      (*max_element(total.begin(),total.end()));
	  }
	}
      }
    }

    // -----------------------------------------------------------------
    // Normalize each row to ensure it's maximum value is equal to
    // the total number of measurements in that row

    if (normalize==norm_row_weight) {
      // Total for each row
      ubvector total(nh);
      // Max in each row
      ubvector max(nh);
      fill(total.begin(),total.end(),0.0);
      fill(max.begin(),max.end(),0.0);
      // Compute total and max in each row
      for(size_t i=0;i<nh;i++) {
	for(size_t j=0;j<nh;j++) {
	  total[i]+=tavg(i,j);
	  if (tavg(i,j)>max[i]) max[i]=tavg(i,j);
	}
      }
      // Do normalization
      for(size_t i=0;i<nh;i++) {
	for(size_t j=0;j<nh;j++) {
	  if ((*max_element(total.begin(),total.end()))>0.0 && max[i]>0.0) {
	    tavg(i,j)*=total[i]/max[i]/
	      (*max_element(total.begin(),total.end()));
	  }
	}
      }
    }

    // -----------------------------------------------------------------
    // contours
    
    std::vector<double> c_ab, c_hi_1, c_hi_2, c_lo_1, c_lo_2, c_mid;

    if (contours==x_contours) {

      cout << "Computing contours giving limits on x-axis." << endl;
      cout << "yval    total   -2sig    -1sig    peak   +1sig   +2sig" << endl;
      cout.precision(4);
      int verbose=1;

      // Compute contours
      for(size_t iy=0;iy<nh;iy++) {
	//cout.width(2);
	//cout << iy << " ";
    
	// Create two vectors to interpolate from
	std::vector<double> yx, yy;
    
	// Fill remaining entries
	for(size_t ix=0;ix<nh;ix++) {
	  yx.push_back(data[0].h.get_x_rep_i(ix));
	  if (ix==0 || ix==nh-1) {
	    yy.push_back(0.0);
	  } else {
	    yy.push_back(tavg(ix,iy));
	  }
	}

	// Compute total;
	double total=vector_integ_linear(yx.size(),yx,yy);
	cout << data[0].h.get_y_rep_i(iy)*yscale << " ";
	cout << total << " ";

	// Only add points if total is large enough
	if (total>min_total) {

	  std::vector<double> locs;
	  double lev;

	  vector_invert_enclosed_sum(one_sigma*total,yx.size(),yx,yy,lev);
	  vector_find_level(lev,yx.size(),yx,yy,locs);
	  double xlo1=0.0;
	  double xhi1=0.0;
	  bool suc=true;
	  if (locs.size()>0) {
	    xlo1=*min_element(locs.begin(),locs.end());
	    xhi1=*max_element(locs.begin(),locs.end());
	  } else {
	    suc=false;
	  }
          
	  vector_invert_enclosed_sum(two_sigma*total,yx.size(),yx,yy,lev);
	  vector_find_level(lev,yx.size(),yx,yy,locs);
	  double xlo2=0.0;
	  double xhi2=0.0;
	  if (locs.size()>0) {
	    xlo2=*min_element(locs.begin(),locs.end());
	    xhi2=*max_element(locs.begin(),locs.end());
	  } else {
	    suc=false;
	  }

	  if (suc) {
	    c_ab.push_back(data[0].h.get_y_rep_i(iy));
	    double yy_max=*max_element(yy.begin(),yy.end());
	    size_t iy=0;
	    for(size_t i=0;i<yy.size();i++) {
	      if (yy[i]==yy_max) iy=i;
	    }
	    c_mid.push_back(yx[iy]);
	    c_lo_1.push_back(xlo1);
	    c_lo_2.push_back(xlo2);
	    c_hi_1.push_back(xhi1);
	    c_hi_2.push_back(xhi2);
	    cout << xlo2*xscale << " " << xlo1*xscale << " " 
		 << (xlo1+xhi1)/2.0*xscale << " " 
		 << xhi1*xscale << " " << xhi2*xscale << endl;
	  } else {
	    cout << endl;
	  }
	} else {
	  cout << endl;
	}
      }

      if (c_ab.size()==0) {
	cout << "No contour information." << endl;
	return gsl_efailed;
      }

      if (false) {
	// Reformat for a LaTeX table
	format_float ff;
	ff.set_pad_zeros(true);
	ff.latex_mode();
    
	interp_o2scl<std::vector<double> > si;
    
	for(double mtab=1.0;mtab<2.21;mtab+=0.1) {
	  ff.set_sig_figs(2);
	  cout << ff.convert(mtab) << " & ";
	  ff.set_sig_figs(4);
	  cout << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_lo_2)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_lo_1)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_mid)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_hi_1)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_hi_2))
	       << " \\\\" << endl;
	}
      }

      cout.precision(6);

    } else if (contours==y_contours) {

      cout.precision(4);
      int verbose=1;

      // Compute "contours"
      for(size_t ix=0;ix<nh-1;ix++) {
	//cout.width(2);
	//cout << ix << " ";
    
	// Create two vectors to interpolate from
	std::vector<double> yx, yy;
    
	// Fill remaining entries
	for(size_t iy=0;iy<nh;iy++) {
	  yx.push_back(data[0].h.get_y_rep_i(iy));
	  if (iy==0 || iy==nh-1) {
	    yy.push_back(0.0);
	  } else {
	    yy.push_back(tavg(ix,iy));
	  }
	}

	// Compute total;
	double total=vector_integ_linear(yx.size(),yx,yy);
	cout << data[0].h.get_x_rep_i(ix)*xscale << " ";
	cout << total << " ";

	// Only add points if total is large enough
	if (total>min_total) {

	  std::vector<double> locs;
	  double lev;

	  vector_invert_enclosed_sum(one_sigma*total,yx.size(),yx,yy,lev);
	  vector_find_level(lev,yx.size(),yx,yy,locs);
	  double xlo1=*min_element(locs.begin(),locs.end());
	  double xhi1=*max_element(locs.begin(),locs.end());
          
	  vector_invert_enclosed_sum(two_sigma*total,yx.size(),yx,yy,lev);
	  vector_find_level(lev,yx.size(),yx,yy,locs);
	  double xlo2=*min_element(locs.begin(),locs.end());
	  double xhi2=*max_element(locs.begin(),locs.end());

	  c_ab.push_back(data[0].h.get_x_rep_i(ix));
	  double yy_max=*max_element(yy.begin(),yy.end());
	  size_t iy=0;
	  for(size_t i=0;i<yy.size();i++) {
	    if (yy[i]==yy_max) iy=i;
	  }
	  c_mid.push_back(yx[iy]);
	  c_lo_1.push_back(xlo1);
	  c_lo_2.push_back(xlo2);
	  c_hi_1.push_back(xhi1);
	  c_hi_2.push_back(xhi2);
	  cout << xlo2*yscale << " " << xlo1*yscale << " " 
	       << (xlo1+xhi1)/2.0*yscale << " " 
	       << xhi1*yscale << " " << xhi2*yscale << endl;
	} else {
	  cout << endl;
	}
      }

      if (c_ab.size()==0) {
	cout << "No contour information." << endl;
	return gsl_efailed;
      }

      if (false) {
	// Reformat for a LaTeX table
	format_float ff;
	ff.set_pad_zeros(true);
	ff.latex_mode();
    
	interp_o2scl<std::vector<double> > si;
	
	for(double mtab=1.0;mtab<2.21;mtab+=0.1) {
	  ff.set_sig_figs(2);
	  cout << ff.convert(mtab) << " & ";
	  ff.set_sig_figs(4);
	  cout << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_lo_2)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_lo_1)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_mid)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_hi_1)) << " & "
	       << ff.convert(si.interp(mtab,c_ab.size(),c_ab,c_hi_2))
	       << " \\\\" << endl;
	}
      }

    }

    // -----------------------------------------------------------------
    
    TApplication *theApp=new TApplication("App",0,NULL);
    TCanvas *c1;

    {
      string wtitle="t1", wtitle2="t2", xtitle="", ytitle="";
      double prmar=0.05;
      double ptmar=0.05;
      double plmar=0.11;
      double pbmar=0.1;

      double xlo=0.0;
      double xhi=-1.0;
      double ylo=0.0;
      double yhi=-1.0;

      bool logz=false;
    
      int verbose=2;

      // Create the canvas
      c1=new TCanvas(wtitle.c_str(),wtitle2.c_str(),600,0,680,640);
      c1->SetFillColor(10);
      TPad *p1=new TPad("t3","",0.0,0.0,1.0,1.0);
      p1->SetFillColor(10);
      p1->SetTopMargin(ptmar);
      p1->SetRightMargin(0.19);
      p1->SetLeftMargin(plmar);
      p1->SetBottomMargin(pbmar);

      if (logx) p1->SetLogx();
      if (logy) p1->SetLogy();

      p1->Draw();
      p1->cd();
      
      size_t nx=data[0].h.size_x(), ny=data[0].h.size_y();

      // -----------------------------------------------------------------
      // Set colors

      //double lcolr=1, lcolg=0.5, lcolb=0.5;
      double lcolr=1, lcolg=1, lcolb=1;
      double hcolr=1, hcolg=0, hcolb=0;

      for(size_t i=0;i<40;i++) {
	TColor *colort=(TColor *)(gROOT->GetListOfColors()->At(11+i));
	double cr=lcolr+((double)i)*(hcolr-lcolr)/39.0;
	double cg=lcolg+((double)i)*(hcolg-lcolg)/39.0;
	double cb=lcolb+((double)i)*(hcolb-lcolb)/39.0;
	if (false && i==0) {
	  cr=1.0;
	  cg=1.0;
	  cb=1.0;
	}
	colort->SetRGB(cr,cg,cb);
      }

      // -----------------------------------------------------------------
      // Determine the plot limits
    
      cout << "Computing plot limits." << endl;

      // Initialize to avoid warnings about uninit'ed vars
      double left=0.0, right=0.0, top=0.0, bottom=0.0;

      left=data[0].h.get_x_low_i(0)*xscale;
      right=data[0].h.get_x_high_i(nx-1)*xscale;
      bottom=data[0].h.get_y_low_i(0)*yscale;
      top=data[0].h.get_y_high_i(ny-1)*yscale;

      // -----------------------------------------------------------------
      // Compute z range
    
      double zmin, zmax;
      zmin=tavg(0,0);
      zmax=tavg(0,0);
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  if (tavg(i,j)<zmin) zmin=tavg(i,j);
	  if (tavg(i,j)>zmax) zmax=tavg(i,j);
	}
      }

      // -----------------------------------------------------------------
      // 

      if (verbose>1) {
	cout << "X range: " << left << " " << right << endl;
	cout << "Y range: " << bottom << " " << top << endl;
	cout << "Z range: " << zmin << " " << zmax << endl;
      }

      // -----------------------------------------------------------------
      // Draw axes

      TH1 *th1=p1->DrawFrame(left,bottom,right,top);
      th1->GetXaxis()->SetLabelFont(132);
      th1->GetYaxis()->SetLabelFont(132);
      th1->GetXaxis()->CenterTitle(kTRUE);
      th1->GetYaxis()->CenterTitle(kTRUE);
      th1->GetXaxis()->SetTitleFont(132);
      th1->GetYaxis()->SetTitleFont(132);
      if (xtitle.length()>0) th1->GetXaxis()->SetTitle(xtitle.c_str());
      if (ytitle.length()>0) th1->GetYaxis()->SetTitle(ytitle.c_str());
    
    
      // -----------------------------------------------------------------
      // Construct box boundaries
    
      ubvector xleft(nx), xright(nx);
      for(size_t i=0;i<nx;i++) {
	xleft[i]=data[0].h.get_x_low_i(i)*xscale;
	xright[i]=data[0].h.get_x_high_i(i)*xscale;
      }
    
      ubvector ybottom(ny), ytop(ny);
      for(size_t i=0;i<ny;i++) {
	ybottom[i]=data[0].h.get_y_low_i(i)*yscale;
	ytop[i]=data[0].h.get_y_high_i(i)*yscale;
      }
      
      // -----------------------------------------------------------------
      // Plot data

      for(size_t i=0;i<nx-1;i++) {
	for(size_t j=0;j<ny-1;j++) {
	
	  // From 11 to 50 inclusive
	  double cd;
	  if (logz) {
	    cd=(log10(tavg(i,j))-log10(zmin))/
	      (log10(zmax)-log10(zmin))*39.0;
	  } else {
	    cd=(tavg(i,j)-zmin)/(zmax-zmin)*39.0;
	  }
	  int color=((int)(cd+11));
	
	  TBox *tb=new TBox(xleft[i],ybottom[j],xright[i],ytop[j]);
	  tb->SetFillColor(color);
	  tb->Draw();

	}
      }

      // -----------------------------------------------------------------
      // Replot axes

      p1->cd();
    
      TGaxis *aleft, *aright, *atop, *abottom;
      if (logy) {
	aright=new TGaxis(right,bottom,right,top,bottom,top,510,"+G");
	aleft=new TGaxis(left,bottom,left,top,bottom,top,510,"-G");
      } else {
	aright=new TGaxis(right,bottom,right,top,bottom,top,510,"+");
	aleft=new TGaxis(left,bottom,left,top,bottom,top,510,"-");
      }
      if (logx) {
	atop=new TGaxis(left,top,right,top,left,right,510,"-G");
	abottom=new TGaxis(left,bottom,right,bottom,left,right,510,"+G");
      } else {
	atop=new TGaxis(left,top,right,top,left,right,510,"-");
	abottom=new TGaxis(left,bottom,right,bottom,left,right,510,"+");
      }
      aleft->SetLabelSize(0.0);
      aright->SetLabelSize(0.0);
      atop->SetLabelSize(0.0);
      abottom->SetLabelSize(0.0);
      aleft->Draw();
      aright->Draw();
      atop->Draw();
      abottom->Draw();

      // -----------------------------------------------------------------
      // Plot z-axis legend

      // Horizontal location of the legend
      // FIXME: These need to be recomputed if logx is true!
      double left3=right+(right-left)*0.02;
      double right3=right+(right-left)*0.09;

      size_t nzl=40;
      for(size_t i=0;i<nzl;i++) {
      
	// Bottom and top margins of this box
	double bottom3, top3;
	if (logy) {
	  bottom3=bottom*pow(top/bottom,((double)i)/nzl);
	  top3=bottom*pow(top/bottom,((double)(i+1))/nzl);
	} else {
	  bottom3=((double)i)/nzl*(top-bottom)+bottom;
	  top3=((double)(i+1))/nzl*(top-bottom)+bottom;
	}

	// Compute color
	int color;
	if (logz) {
	  double value=
	    pow(10.0,((double)i+1)/nzl*(log10(zmax)-log10(zmin))+log10(zmin));
	  double cd=(log10(value)-log10(zmin))/(log10(zmax)-log10(zmin))*39.0;
	  color=((int)(cd+11));
	} else {
	  color=((int)((i+1)*39.0/((double)nzl)+11.0));
	}
      
	// Draw the box
	TBox *tb=new TBox(left3,bottom3,right3,top3);
	tb->SetFillColor(color);
	tb->Draw();
      }

      // Plot the RHS axis for the legend
      TGaxis *a1;
      if (logz) {
	a1=new TGaxis(right3,bottom,right3,top,zmin,zmax,510,"+LG");
      } else {
	a1=new TGaxis(right3,bottom,right3,top,zmin,zmax,510,"+L");
      }
      a1->SetLabelFont(132);
      a1->SetLabelSize(0.04);
      a1->SetLabelOffset(0.01);
      a1->SetTickSize(0.13);
      a1->CenterTitle(kTRUE);
      a1->Draw();

      TLine *l1=new TLine(left3,bottom,left3,top);
      l1->Draw();
      TLine *l2=new TLine(left3,bottom,right3,bottom);
      l2->Draw();
      TLine *l3=new TLine(left3,top,right3,top);
      l3->Draw();
      
      // -----------------------------------------------------------------
      // Plot axis labels
      
      TLatex tt;
      tt.SetTextAlign(22);
      tt.SetTextFont(132);
      
      double xlx=(left+right)/2.0;
      double xly=bottom-(top-bottom)/15.0;
      if (logx) {
	xlx=sqrt(left*right);
      }
      if (logy) {
	xly=bottom/pow(top/bottom,1.0/15.0);
      }
      tt.DrawLatex(xlx,xly,xlabel.c_str());

      double ylx=left-(right-left)/9.0;
      double yly=(top+bottom)/2.0;
      if (logx) {
	ylx=left/pow(right/left,1.0/9.0);
      }
      if (logy) {
	yly=sqrt(top*bottom);
      }

      tt.SetTextAngle(90);
      tt.DrawLatex(ylx,yly,ylabel.c_str());
      tt.SetTextAngle(0);
      tt.DrawLatex(left+(right-left)/5.0,bottom+(top-bottom)/5.0,
		   title.c_str());
    }
  
    // -----------------------------------------------------------------
    // Plot contour lines

    if (contours==x_contours) {
      //TGraph *gmid=new TGraph(c_ab.size());
      TGraph *glo1=new TGraph(c_ab.size());
      TGraph *glo2=new TGraph(c_ab.size());
      TGraph *ghi1=new TGraph(c_ab.size());
      TGraph *ghi2=new TGraph(c_ab.size());
      for(size_t i=0;i<c_ab.size();i++) {
	//gmid->SetPoint(i,c_mid[i],c_ab[i]);
	glo1->SetPoint(i,c_lo_1[i]*xscale,c_ab[i]*yscale);
	ghi1->SetPoint(i,c_hi_1[i]*xscale,c_ab[i]*yscale);
	glo2->SetPoint(i,c_lo_2[i]*xscale,c_ab[i]*yscale);
	ghi2->SetPoint(i,c_hi_2[i]*xscale,c_ab[i]*yscale);
      }
      glo1->SetLineStyle(2);
      ghi1->SetLineStyle(2);
      //gmid->Draw();
      glo1->Draw();
      glo2->Draw();
      ghi1->Draw();
      ghi2->Draw();
    } else if (contours==y_contours) {
      //TGraph *gmid=new TGraph(c_ab.size());
      TGraph *glo1=new TGraph(c_ab.size());
      TGraph *glo2=new TGraph(c_ab.size());
      TGraph *ghi1=new TGraph(c_ab.size());
      TGraph *ghi2=new TGraph(c_ab.size());
      for(size_t i=0;i<c_ab.size();i++) {
	//gmid->SetPoint(i,c_mid[i],c_ab[i]);
	glo1->SetPoint(i,c_ab[i]*xscale,c_lo_1[i]*yscale);
	ghi1->SetPoint(i,c_ab[i]*xscale,c_hi_1[i]*yscale);
	glo2->SetPoint(i,c_ab[i]*xscale,c_lo_2[i]*yscale);
	ghi2->SetPoint(i,c_ab[i]*xscale,c_hi_2[i]*yscale);
      }
      glo1->SetLineStyle(2);
      ghi1->SetLineStyle(2);
      //gmid->Draw();
      glo1->Draw();
      glo2->Draw();
      ghi1->Draw();
      ghi2->Draw();
    } else if (contours==xy_contours) {
      
      ubvector integx(101), integy(101);
      for(size_t i=0;i<101;i++) {
 	integx[i]=((double)i)/100.0;
      }
      for(size_t k=0;k<101;k++) {
 	integy[k]=0.0;
 	for(size_t i=0;i<nh;i++) {
 	  for(size_t j=0;j<nh;j++) {
 	    if (tavg(i,j)>integx[k]) integy[k]+=tavg(i,j);
 	  }
 	}
      }
      double total=integy[0], one_sig=0.0, two_sig=0.0;
      for(size_t i=0;i<100;i++) {
 	if (integy[i]>two_sigma*total && integy[i+1]<two_sigma*total) {
 	  two_sig=(integx[i]+integx[i+1])/2.0;
 	  i=100;
 	}
      }
      for(size_t i=0;i<100;i++) {
 	if (integy[i]>one_sigma*total && integy[i+1]<one_sigma*total) {
 	  one_sig=(integx[i]+integx[i+1])/2.0;
 	  i=100;
 	}
      }
      
      if (one_sig>0.0 && two_sig>0.0) {
 	contour co;
 	ubvector xg(nh), yg(nh);
 	for(size_t i=0;i<nh;i++) {
 	  xg[i]=data[0].h.get_x_rep_i(i);
 	  yg[i]=data[0].h.get_y_rep_i(i);
 	}
 	co.set_data(nh,nh,xg,yg,tavg);
 	double levels[2]={one_sig,two_sig};
 	co.set_levels(2,levels);
 	
 	vector<contour_line> conts;
 	co.calc_contours(conts);
 	size_t nc=conts.size();
 	
 	cout << "Number of contours: " << nc << endl;
	double minx=100.0, miny=100.0, maxx=0.0, maxy=0.0;
 	
 	for(size_t i=0;i<nc;i++) {
 	  TGraph *gr1=new TGraph(conts[i].x.size());
 	  for(size_t j=0;j<conts[i].x.size();j++) {
 	    gr1->SetPoint(j,conts[i].x[j]*xscale,conts[i].y[j]*yscale);
	    if (conts[i].x[j]<minx) minx=conts[i].x[j];
	    if (conts[i].x[j]>maxx) maxx=conts[i].x[j];
	    if (conts[i].y[j]<miny) miny=conts[i].y[j];
	    if (conts[i].y[j]>maxy) maxy=conts[i].y[j];
 	  }
 	  if (conts[i].level>sqrt(one_sig*two_sig)) gr1->SetLineStyle(2);
 	  gr1->Draw();
 	}

	cout << "x min, max: " << minx << " " << maxx << endl;
	cout << "y min, max: " << miny << " " << maxy << endl;

      } else {
 	cout << "Error in contours." << endl;
      }
    }

    // -----------------------------------------------------------------

    for(size_t i=0;i<ellvec.size();i+=4) {
      cout << ellvec[i] << " " << ellvec[i+1] << " " 
	   << ellvec[i+2] << " " << ellvec[i+3] << endl;
	
      TEllipse *tel1=new TEllipse(ellvec[i],ellvec[i+1],ellvec[i+2],
				  ellvec[i+3]);
      tel1->SetFillStyle(0);
      tel1->Draw();
    }
  
    // -----------------------------------------------------------------

    c1->cd();
    c1->Update();
  
    // Output to file or screen
    if (out_file.length()>0) {
      c1->Print(out_file.c_str());
    } else {
      theApp->Run(kTRUE);
    }

    return 0;
  }

  /** \brief Desc
   */
  void run(int argc, char *argv[]) {
  
    // ---------------------------------------
    // Specify command-line option object
    
    cli_readline cl;
    cl.prompt="plot2d_class> ";
    cl.gnu_intro=false;
    
    // ---------------------------------------
    // Set options
    
    static const int nopt=7;
    comm_option_s options[nopt]={
      {'p',"plot","Short desc.",2,-1,"<object> <file1> [file2 file3 ...]",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::plot),
       cli::comm_option_both},
      {'t',"tabulate","Short desc.",2,-1,"<object> <file1> [file2 file3 ...]",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::tabulate),
       cli::comm_option_both},
      {'s',"show","Short desc.",4,4,"<object> <file1>",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::show),
       cli::comm_option_both},
      {'r',"row","Short desc.",3,3,"<object> <file1>",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::row),
       cli::comm_option_both},
      {'c',"contours","Short desc.",1,1,"<type>",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::contour_fun),
       cli::comm_option_both},
      {'n',"normalize","Short desc.",1,1,"<type>",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::normal_fun),
       cli::comm_option_both},
      {'e',"ell","Short desc.",4,4,"<x> <y> <xw> <yw>",
       ((string)"Long ")+"desc.",
       new comm_option_mfptr<plot2d_class>(this,&plot2d_class::ell),
       cli::comm_option_both}
    };
    cl.set_comm_option_vec(nopt,options);

    // ---------------------------------------
    // Set parameters
    
    cli::parameter_double p_xscale;
    p_xscale.d=&xscale;
    p_xscale.help="";
    cl.par_list.insert(make_pair("xscale",&p_xscale));

    cli::parameter_double p_min_total;
    p_min_total.d=&min_total;
    p_min_total.help="";
    cl.par_list.insert(make_pair("min_total",&p_min_total));

    cli::parameter_double p_yscale;
    p_yscale.d=&yscale;
    p_yscale.help="";
    cl.par_list.insert(make_pair("yscale",&p_yscale));

    cli::parameter_string p_xlabel;
    p_xlabel.str=&xlabel;
    p_xlabel.help="";
    cl.par_list.insert(make_pair("xlabel",&p_xlabel));

    cli::parameter_string p_ylabel;
    p_ylabel.str=&ylabel;
    p_ylabel.help="";
    cl.par_list.insert(make_pair("ylabel",&p_ylabel));

    cli::parameter_bool p_logx;
    p_logx.b=&logx;
    p_logx.help="";
    cl.par_list.insert(make_pair("logx",&p_logx));

    cli::parameter_bool p_logy;
    p_logy.b=&logy;
    p_logy.help="";
    cl.par_list.insert(make_pair("logy",&p_logy));

    cli::parameter_string p_out_file;
    p_out_file.str=&out_file;
    p_out_file.help="";
    cl.par_list.insert(make_pair("out_file",&p_out_file));

    cli::parameter_string p_title;
    p_title.str=&title;
    p_title.help="";
    cl.par_list.insert(make_pair("title",&p_title));

    // ---------------------------------------
    // Process command-line arguments and run
    
    cl.run_auto(argc,argv);
    
    return;
  }

};

int main(int argc, char *argv[]) {
  cout.setf(ios::scientific);
  plot2d_class pc;
  pc.run(argc,argv);
  return 0;
}
