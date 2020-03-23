/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2019, Andrew W. Steiner
  
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

#include "process.h"

#include <o2scl/mcarlo_miser.h>
#include <o2scl/mcarlo_vegas.h>
#include <o2scl/multi_funct.h>

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

process::process() : one_sigma(gsl_sf_erf(1.0/sqrt(2.0))),
		     two_sigma(gsl_sf_erf(2.0/sqrt(2.0))),
		     three_sigma(gsl_sf_erf(3.0/sqrt(2.0))) {
  verbose=1;
  xscale=1.0;
  yscale=1.0;
  xset=false;
  user_xlow=0.0;
  user_xhigh=0.0;
  yset=false;
  user_ylow=0.0;
  user_yhigh=0.0;
  errors=true;
  hist_size_int=100;
  logx=false;
  logy=false;
  logz=false;
  line_start=0;
  ff.latex_mode();
  ff.set_sig_figs(4);
  ff.set_pad_zeros(true);
  ff.set_exp_limits(-5,5);
  n_blocks=1;
  grid_init=0;
  weights_col="";
}

void process::swap(double &w1, size_t &i1, ubvector &x1, 
		    double &w2, size_t &i2, ubvector &x2) {
  double temp=w1; w1=w2; w2=temp;
  size_t itemp=i1; i1=i2; i2=itemp;
  o2scl::vector_swap<ubvector,ubvector,double>(x1,x2);
  return;
}

double process::dist_sq(const ubvector &x, size_t row) {
  double dist=0.0;
  for(size_t i=0;i<x.size();i++) {
    dist+=pow((x[i]-bfactor_data[i][row])/scale[i],2.0);
  }
  return dist;
}

double process::func(size_t n, const ubvector &x) {

  // Initialize from first three points
  size_t i0=0, i1=1, i2=2;
  double min0=dist_sq(x,0);
  double min1=dist_sq(x,1);
  double min2=dist_sq(x,2);
  ubvector x0(n), x1(n), x2(n);
  for(size_t i=0;i<n;i++) {
    x0[i]=bfactor_data[i][0];
    x1[i]=bfactor_data[i][1];
    x2[i]=bfactor_data[i][2];
  }
  if (min2<min1) swap(min1,i1,x1,min2,i2,x2);
  if (min1<min0) swap(min0,i0,x0,min1,i1,x1);
  if (min2<min1) swap(min1,i1,x1,min2,i2,x2);
  
  // Sort through the rest of the chain
  for(size_t i=3;i<bfactor_data[0].size();i++) {
    double thisd=dist_sq(x,i);
    if (thisd<min2) {
      min2=thisd;
      i2=i;
      for(size_t j=0;j<n;j++) {
	x2[j]=bfactor_data[j][i];
      }
      if (min2<min1) swap(min1,i1,x1,min2,i2,x2);
      if (min1<min0) swap(min0,i0,x0,min1,i1,x1);
    }
  }

  // Return the final result
  double norm=1.0/min0+1.0/min1+1.0/min2;
  double ret=(bfactor_data[bfactor_data.size()-1][i0]/min0+
	      bfactor_data[bfactor_data.size()-1][i1]/min1+
	      bfactor_data[bfactor_data.size()-1][i2]/min2)/norm;

  if (false) {
    cout << i0 << " " << i1 << " " << i2 << " "
	 << bfactor_data[bfactor_data.size()-1][i0] << " "
	 << bfactor_data[bfactor_data.size()-1][i1] << " "
	 << bfactor_data[bfactor_data.size()-1][i2] << endl;
    cout << min0 << " " << min1 << " " << min2 << " " << norm << " " 
	 << ret << endl;
    char ch;
    cin >> ch;
  }
    
  return ret;
}

int process::set_params_limits(std::vector<std::string> &sv, bool itive_com) {
  size_t np=(sv.size()-1)/3;
  std::cout << np << std::endl;
  if (np==0) {
    cerr << "No parameters specified." << endl;
  }
  x_low.resize(np);
  x_high.resize(np);
  x_params.resize(np);
  for(size_t j=0;j<np;j++) {
    x_params[j]=sv[j*3+1];
    x_low[j]=function_to_double(sv[j*3+2]);
    x_high[j]=function_to_double(sv[j*3+3]);
    cout << "Added parameter: " << x_params[j] << " with lower limit "
	 << x_low[j] << " and upper limit " << x_high[j] << endl;
  }
  return 0;
}

int process::bfactor(std::vector<std::string> &sv, bool itive_com) {

  if (x_params.size()==0) {
    cerr << "Parameters must be set first." << endl;
    return 1;
  }
  
  // file list
  vector<string> files;

  // Form list of data files
  for(size_t i=1;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();

  // Number of parameters and scales
  size_t n_params=x_params.size();
  cout << n_params << " parameters." << endl;
  scale.resize(n_params);
  for(size_t i=0;i<n_params;i++) {
    scale[i]=fabs(x_high[i]-x_low[i]);
    cout << "Parameter " << i << " " << x_params[i] << " "
	 << x_low[i] << " " << x_high[i] << " " << scale[i] << endl;
  }
  
  // Setup bfactor_data with empty columns, adding one more column for
  // the weights
  vector<double> empty;
  for(size_t i=0;i<n_params;i++) {
    bfactor_data.push_back(empty);
  }
  bfactor_data.push_back(empty);
  
  // Read data
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;

    // Open file
    cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);

    // Read table
    std::string tab_name;
    table_units<> tab;
    hdf_input(hf,tab,tab_name);

    // Parse parameters into bfactor_data
    for(size_t ell=0;ell<n_params;ell++) {
      for(size_t k=0;k<tab.get_nlines();k++) {
	bfactor_data[ell].push_back(tab.get(x_params[ell],k));
      }
    }
    
    // Parse weights into bfactor_data
    for(size_t k=0;k<tab.get_nlines();k++) {
      bfactor_data[bfactor_data.size()-1].push_back
	(exp(tab.get("log_wgt",k)));
    }

    // Go to next file
  }

  {
    mcarlo_miser<> gm;
    multi_funct mff=std::bind(std::mem_fn<double(size_t,const ubvector &)>
			      (&process::func),this,std::placeholders::_1,
			      std::placeholders::_2);
    
    gm.n_points=100000;
    double res, err;
    cout << "Computing integral:" << endl;
    gm.minteg_err(mff,n_params,x_low,x_high,res,err);
    
    double hc_vol=1.0;
    for(size_t k=0;k<n_params;k++) {
      hc_vol*=fabs(x_high[k]-x_low[k]);
    }
    
    cout << "Hypercube volume: " << hc_vol << endl;
    cout << "Integral, error: " << res << " " << err << endl;
  }
  {
    mcarlo_vegas<> gm;
    multi_funct mff=std::bind(std::mem_fn<double(size_t,const ubvector &)>
			      (&process::func),this,std::placeholders::_1,
			      std::placeholders::_2);
    
    gm.n_points=100000;
    double res, err;
    cout << "Computing integral:" << endl;
    gm.minteg_err(mff,n_params,x_low,x_high,res,err);
    
    double hc_vol=1.0;
    for(size_t k=0;k<n_params;k++) {
      hc_vol*=fabs(x_high[k]-x_low[k]);
    }
    
    cout << "Hypercube volume: " << hc_vol << endl;
    cout << "Integral, error: " << res << " " << err << endl;
  }

  return 0;
}
	
int process::auto_corr(std::vector<std::string> &sv, bool itive_com) {

  // Setup histogram size
  size_t hist_size=(size_t)hist_size_int;
  if (hist_size_int<=0) {
    cout << "Histogram size <=0, using 100." << endl;
    hist_size=100;
  }
    
  // column name
  string name=sv[1];

  // (output file is sv[2])

  // Parse file list
  vector<string> files;
  for(size_t i=3;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();

  // Read data

  vector<double> values;

  for(size_t i=0;i<nf;i++) {
    hdf_file hf;

    // Open file
    cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);

    // Get number of chains
    size_t n_chains;
    hf.get_szt("n_chains",n_chains);
    cout << n_chains << " separate chains." << endl;

    size_t line_counter=0;

    // Read each chain
    for(size_t j=0;j<n_chains;j++) {

      // Read table
      std::string tab_name="markov_chain_"+szttos(j);
      table_units<> tab;
      hdf_input(hf,tab,tab_name);
      cout << "Read table " << tab_name << endl;

      // Parse table into values and weights
      for(size_t k=0;k<tab.get_nlines();k++) {
	if (((int)line_counter)>=line_start) {
	  values.push_back(tab.get(name,k));
	}
	line_counter++;
      }
      // Go to next chain
    }
    // Go to next file
  }

  // Store the autocorrelation data
  size_t kmax=values.size()/3;
  vector<double> kvec(kmax-1), acvec(kmax-1), acerr(kmax-1), acfit(kmax-1);
  for(size_t k=1;k<kmax;k++) {
    kvec[k-1]=((double)k);
    acvec[k-1]=vector_lagk_autocorr(values.size(),values,k);
    acerr[k-1]=fabs(acvec[k-1])/10.0;
    acfit[k-1]=0.0;
  }

  // Put the autocorrelation data into a table
  table<> tab;
  tab.line_of_names("k ac err fit");
  tab.inc_maxlines(kvec.size());
  tab.set_nlines(kvec.size());
  tab.swap_column_data("k",kvec);
  tab.swap_column_data("ac",acvec);
  tab.swap_column_data("err",acerr);
  tab.swap_column_data("fit",acfit);

  // Create the output file
  hdf_file hf;
  hf.open_or_create(sv[2]);
  hdf_output(hf,tab,"ac");
  hf.close();
      
  return 0;
}

int process::xlimits(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<3) {
    // If there are less than 2 arguments, reset the y limits
    if (verbose>0) {
      cout << "Setting 'xset' to false." << endl;
    }
    xset=false;
    return 0;
  }
    
  if (sv.size()==3) {
    // If there are 2 arguments, use them to set the y limits
    user_xlow=o2scl::function_to_double(sv[1]);
    user_xhigh=o2scl::function_to_double(sv[2]);
    xset=true;
    if (verbose>0) {
      cout << "X limits are " << user_xlow << " and " 
	   << user_xhigh << " ." << endl;
    }
    return 0;
  }

  // If there are three arguments, presume they're stored in
  // a specified file as named double-precision numbers
  string file=sv[1];
  string low_name=sv[2];
  string high_name=sv[3];

  hdf_file hf;
  if (verbose>1) {
    cout << "Reading file named '" << file << "'." << endl;
  }
  hf.open(file);
  xset=true;
  hf.getd(low_name,user_xlow);
  hf.getd(high_name,user_xhigh);
  hf.close();

  if (verbose>0) {
    cout << "X limits are " << user_xlow << " and " 
	 << user_xhigh << " ." << endl;
  }
      
  return 0;
}

int process::ylimits(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<3) {
    // If there are less than 2 arguments, reset the y limits
    if (verbose>0) {
      cout << "Setting 'yset' to false." << endl;
    }
    yset=false;
    return 0;
  }

  if (sv.size()==3) {
    // If there are 2 arguments, use them to set the y limits
    user_ylow=o2scl::function_to_double(sv[1]);
    user_yhigh=o2scl::function_to_double(sv[2]);
    yset=true;
    if (verbose>0) {
      cout << "Y limits are " << user_ylow << " and " 
	   << user_yhigh << " ." << endl;
    }
    return 0;
  }

  // If there are three arguments, presume they're stored in
  // a specified file as named double-precision numbers
  string file=sv[1];
  string low_name=sv[2];
  string high_name=sv[3];

  hdf_file hf;
  if (verbose>1) {
    cout << "Reading file named '" << file << "'." << endl;
  }
  hf.open(file);
  yset=true;
  hf.getd(low_name,user_ylow);
  hf.getd(high_name,user_yhigh);
  hf.close();

  if (verbose>0) {
    cout << "Y limits are " << user_ylow << " and " 
	 << user_yhigh << " ." << endl;
  }

  return 0;
}
  
int process::mass_sel(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<4) {
    cerr << "Not enough arguments to 'mass-sel'." << endl;
    return 1;
  }
  
  // mass value
  double m_val=o2scl::stod(sv[1]);
  string suffix=sv[2];
  if (verbose>0) {
    cout << "Mass value " << m_val << " and suffix " << suffix << endl;
  }
  
  // Form list of data files
  vector<string> files;
  for(size_t i=3;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
  
  // ------------------------------------------------------------
  // Read all of the data files in order
      
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;
    
    // Open file with write access
    if (verbose>0) cout << "Opening file: " << files[i] << endl;
    hf.open(files[i],true);

    int addl_quants;
    hf.geti("addl_quants",addl_quants);
    if (addl_quants>0 && verbose>0) {
      cout << "Value 'addl_quants' is true." << endl;
    }

    // Reading mass grid

    uniform_grid<> m_grid;
    hdf_input(hf,m_grid,"m_grid");
    vector<double> m_vec;
    m_grid.vector(m_vec);
    
    // Get number of chains
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    if (verbose>0) {
      if (n_chains>1) {
	cout << n_chains << " separate chains." << endl;
      } else {
	cout << "1 chain." << endl;
      }
    }

    // Read each chain
    for(size_t j=0;j<n_chains;j++) {

      // Read table
      std::string tab_name="markov_chain_"+szttos(j);
      table_units<> tab;
      hdf_input(hf,tab,tab_name);

      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows_func(constraint);
	size_t nlines_new=tab.get_nlines();
	if (verbose>0) {
	  cout << "Applied constraint \"" << constraint
	       << "\" and went from " << nlines_old << " to "
	       << nlines_new << " lines." << endl;
	}
      }

      tab.new_column(((std::string)"R")+suffix);
      tab.set_unit(((std::string)"R")+suffix,"km");
      tab.new_column(((std::string)"PM")+suffix);
      tab.set_unit(((std::string)"PM")+suffix,"1/fm^4");
      if (addl_quants>0) {
	tab.new_column(((std::string)"I")+suffix);
	tab.set_unit(((std::string)"I")+suffix,"Msun*km^2");
	tab.new_column(((std::string)"MB")+suffix);
	tab.set_unit(((std::string)"MB")+suffix,"Msun");
	tab.new_column(((std::string)"BE")+suffix);
	tab.set_unit(((std::string)"BE")+suffix,"Msun");
      }
      
      // Parse table into values and weights
      for(size_t k=0;k<tab.get_nlines();k++) {

	vector<double> R_y, PM_y, I_y, MB_y, BE_y;
	for(size_t hh=0;hh<m_grid.get_npoints();hh++) {
	  R_y.push_back(tab.get(((std::string)"R_")+o2scl::szttos(hh),k));
	  PM_y.push_back(tab.get(((std::string)"PM_")+o2scl::szttos(hh),k));
	  if (addl_quants>0) {
	    I_y.push_back(tab.get(((std::string)"I_")+o2scl::szttos(hh),k));
	    MB_y.push_back(tab.get(((std::string)"MB_")+o2scl::szttos(hh),k));
	    BE_y.push_back(tab.get(((std::string)"BE_")+o2scl::szttos(hh),k));
	  }
	}

	interp<vector<double> > oi;
	oi.set_type(itp_linear);
	tab.set(((std::string)"R")+suffix,k,
		oi.eval(m_val,m_grid.get_npoints(),m_vec,R_y));
	tab.set(((std::string)"PM")+suffix,k,
		oi.eval(m_val,m_grid.get_npoints(),m_vec,PM_y));
	if (addl_quants>0) {
	  tab.set(((std::string)"I")+suffix,k,
		  oi.eval(m_val,m_grid.get_npoints(),m_vec,I_y));
	  tab.set(((std::string)"MB")+suffix,k,
		  oi.eval(m_val,m_grid.get_npoints(),m_vec,MB_y));
	  tab.set(((std::string)"BE")+suffix,k,
		  oi.eval(m_val,m_grid.get_npoints(),m_vec,BE_y));
	}
      }

      if (verbose>0) {
	cout << "Table " << tab_name << " lines: " 
	     << tab.get_nlines();
	cout << endl;
      }

      // Output table to file
      hdf_output(hf,tab,tab_name);
      
      // Go to next chain
    }

    // Close file
    hf.close();
    
    // Go to next file
  }
  if (verbose>0) {
    cout << "Done." << endl;
  }

  return 0;
}

int process::hist(std::vector<std::string> &sv, bool itive_com) {
    
  // Setup histogram size
  size_t hist_size=(size_t)hist_size_int;
  if (hist_size_int<=0) {
    cout << "Histogram size <=0, using 100." << endl;
    hist_size=100;
  }
    
  // column name
  string name=sv[1];

  // (output file is sv[2])

  // Form list of data files
  vector<string> files;
  for(size_t i=3;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
  
  
  // Storage for all of the values and weights
  vector<double> values, weights;

  // ------------------------------------------------------------
  // Read all of the data files in order
      
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;

    // Open file
    if (verbose>0) cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);

    // Get number of chains
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    if (verbose>0) {
      if (n_chains>1) {
	cout << n_chains << " separate chains." << endl;
      } else {
	cout << "1 chain." << endl;
      }
    }

    size_t line_counter=0;

    // Read each chain
    for(size_t j=0;j<n_chains;j++) {

      // Read table
      std::string tab_name="markov_chain_"+szttos(j);
      table_units<> tab;
      hdf_input(hf,tab,tab_name);

      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows_func(constraint);
	size_t nlines_new=tab.get_nlines();
	if (verbose>0) {
	  cout << "Applied constraint \"" << constraint
	       << "\" and went from " << nlines_old << " to "
	       << nlines_new << " lines." << endl;
	}
      }

      // Parse table into values and weights
      for(size_t k=0;k<tab.get_nlines();k++) {
	if (((int)line_counter)>=line_start) {
	  values.push_back(tab.get(name,k));

	  if (weights_col.size()>0) {
	    weights.push_back(tab.get(weights_col,k));
	  } else {
	    weights.push_back(1.0);
	  }

	}
	line_counter++;
      }

      if (verbose>0) {
	cout << "Table " << tab_name << " lines: " 
	     << tab.get_nlines();
	if (line_start>0) {
	  cout << " skipping: " << line_start;
	}
	cout << endl;
      }

      // Go to next chain
    }
    // Go to next file
  }
  if (verbose>0) {
    cout << "Done with reading files.\n" << endl;
  }

  // ------------------------------------------------------------

  // ------------------------------------------------------------
  // Get limits for grid
      
  double min, max;
  if (!xset) {
    min=*min_element(values.begin(),values.end());
    max=*max_element(values.begin(),values.end());
    min*=xscale;
    max*=xscale;
  } else {
    min=user_xlow;
    max=user_xhigh;
  }
  if (verbose>0) {
    cout << "xsize, min, max: " << values.size() << " " 
	 << min << " " << max << endl;
  }

  // ------------------------------------------------------------
  // Double check limits in the case of a logarithmic grid
      
  if (logx) {
    if (max<=0.0) {
      O2SCL_ERR2("Both min and max negative with logarithmic ",
		 "grid in hist().",exc_efailed);
    }
    // Reset 'min' to be smallest positive value
    if (min<=0.0) {
      double min_new=max;
      for(size_t i=0;i<values.size();i++) {
	if (values[i]>0.0 && values[i]<min_new) min_new=values[i];
      }
      cout << "Resetting 'min' from " << min << " to " 
	   << min_new << "." << endl;
      min=min_new;
    }
  }

  // ------------------------------------------------------------
  // Create expval_scalar objects

  if (n_blocks<=0) {
    cout << "Variable 'n_blocks' <= 0. Setting to 20." << endl;
    n_blocks=20;
  }
  vector<expval_scalar> sev(hist_size);
  for(size_t i=0;i<hist_size;i++) {
    sev[i].set_blocks(n_blocks,1);
  }

  // ------------------------------------------------------------
  // Setup grid and histogram. We make the grid slightly bigger
  // than the min-max range to ensure finite precision errors
  // don't give problems with elements outside the histogram.
      
  o2scl::hist h;
  uniform_grid<double> grid;
  if (logx) {
    grid=uniform_grid_log_end<double>(min*(1.0-1.0e-6),
				      max*(1.0+1.0e-6),hist_size);
    h.set_bin_edges(grid);
  } else {
    grid=uniform_grid_end<double>(min-fabs(min)*1.0e-6,
				  max+fabs(max)*1.0e-6,hist_size);
    h.set_bin_edges(grid);
  }
      
  // ------------------------------------------------------------
  // Fill expval_scalar objects using histogram
      
  size_t block_size=values.size()/((size_t)n_blocks);
  for(int i=0;i<n_blocks;i++) {
    if (verbose>1) {
      cout << "Block " << i << endl;
    }
    for(size_t j=0;j<block_size;j++) {
      double val=values[i*block_size+j]*xscale;
      if (val>min && val<max) {
	h.update(val,weights[i*block_size+j]);
      }
    }
    // Update expval_scalar object with histogram for this block
    for(size_t j=0;j<hist_size;j++) {
      sev[j].add(h[j]);
    }
    h.clear_wgts();
  }

  // ------------------------------------------------------------
  // Setup table
      
  table<> t;
  if (errors) {
    t.line_of_names("reps avgs errs plus minus");
  } else {
    t.line_of_names("reps avgs errs");
  }
  t.set_nlines(hist_size);

  // ------------------------------------------------------------
  // Compute values for table
      
  double avg, std_dev, avg_err;
  for(size_t i=0;i<hist_size;i++) {
    t.set("reps",i,h.get_rep_i(i));
    sev[i].current_avg(avg,std_dev,avg_err);
    t.set("avgs",i,avg);
    t.set("errs",i,avg_err);
  }

  // Set endpoints to zero for vector_invert_enclosed_sum()
  t.set("avgs",0,0.0);
  t.set("avgs",hist_size-1,0.0);
    
  // Compute total integral
  double total=vector_integ_xy_interp(hist_size,t["reps"],t["avgs"],
				      itp_linear);
  if (verbose>0) cout << "Total integral: " << total << endl;
  
  std::vector<double> locs;
  double lev;

  // Get one-sigma error ranges
  vector_invert_enclosed_sum(one_sigma*total,hist_size,
			     t["reps"],t["avgs"],lev);

  vector_find_level(lev,hist_size,t["reps"],t["avgs"],locs);
  double xlo1=locs[0], xhi1=locs[0];
  for(size_t i=0;i<locs.size();i++) {
    if (locs[i]<xlo1) xlo1=locs[i];
    if (locs[i]>xhi1) xhi1=locs[i];
  }

  // Get two-sigma error ranges
  vector_invert_enclosed_sum(two_sigma*total,hist_size,
			     t["reps"],t["avgs"],lev);
  vector_find_level(lev,hist_size,t["reps"],t["avgs"],locs);
  double xlo2=locs[0], xhi2=locs[0];
  for(size_t i=0;i<locs.size();i++) {
    if (locs[i]<xlo2) xlo2=locs[i];
    if (locs[i]>xhi2) xhi2=locs[i];
  }

  // Get peak position
  double xpeak=vector_max_quad_loc<std::vector<double>,double>
    (hist_size,t["reps"],t["avgs"]);
      
  if (verbose>0) {
    cout << name << " -2,-1,0,+1,+2 sigma: " << endl;
    if (false) {
      cout.unsetf(ios::scientific);
      cout.precision(4);
      cout << xlo2 << " & " << xlo1 << " & " 
	   << xpeak << " & " << xhi1 << " & " << xhi2 << " \\\\" << endl;
      cout.precision(6);
      cout.setf(ios::scientific);
    }
    cout << ff.convert(xlo2) << " & " << ff.convert(xlo1) << " & " 
	 << ff.convert(xpeak) << " & " << ff.convert(xhi1) << " & " 
	 << ff.convert(xhi2) << " \\\\" << endl;
  }

  // Setup min and max y values (distinct from 'min' and
  // 'max' for the grid above)
  double xmin=t.min("reps");
  double xmax=t.max("reps");
  double ymin=t.min("avgs");
  double ymax=t.max("avgs");

  // Create columns for errors
  if (errors) {
    for(size_t i=0;i<hist_size;i++) {
      if (t.get("avgs",i)>=0.0) {
	if (t.get("avgs",i)-t.get("errs",i)>=0.0) {
	  t.set("minus",i,t.get("avgs",i)-t.get("errs",i));
	} else {
	  t.set("minus",i,0.0);
	}
	t.set("plus",i,t.get("avgs",i)+t.get("errs",i));
      } else {
	t.set("minus",i,t.get("avgs",i)-t.get("errs",i));
	if (t.get("avgs",i)+t.get("errs",i)>0.0) {
	  t.set("plus",i,0.0);
	} else {
	  t.set("plus",i,t.get("avgs",i)+t.get("errs",i));
	}
      }
    }
    double plus_max=t.max("plus");
    double minus_min=t.min("minus");
    if (plus_max>ymax) ymax=plus_max;
    if (minus_min<ymin) ymin=minus_min;
  }

  if (verbose>0) {
    cout << "Writing to file " << sv[2] << endl;
  }
  hdf_file hf;
  hf.open_or_create(sv[2]);
  hdf_output(hf,t,"hist_table");
  hdf_output(hf,grid,"xgrid");
  hf.setd("xmin",xmin);
  hf.setd("xmax",xmax);
  if (logx) {
    hf.seti("logx",1);
  } else {
    hf.seti("logx",0);
  }
  hf.setd("ymin",ymin);
  hf.setd("ymax",ymax);
  hf.setd("xlo2",xlo2);
  hf.setd("xhi2",xhi2);
  hf.setd("xlo1",xlo1);
  hf.setd("xhi1",xhi1);
  hf.setd("xpeak",xpeak);
  hf.close();

  return 0;
}

int process::curve_set(std::vector<std::string> &sv, bool itive_com) {

  string infile=sv[1];

  string grid_name=sv[2];
  
  double low=o2scl::stod(sv[3]);
  double high=o2scl::stod(sv[4]);
  
  // prefix name
  string prefix=sv[5];

  string outfile=sv[6];

  // ------------------------------------------------------------
  // Read all of the data files in order
  
  hdf_file hf;
  
  // Open file
  if (verbose>0) cout << "Opening file: " << infile << endl;
  hf.open(infile);
  
  // Read table
  std::string tab_name="markov_chain_0";
  table_units<> tab;
  hdf_input(hf,tab,tab_name);
  hf.close();
  
  if (constraint.size()>0) {
    size_t nlines_old=tab.get_nlines();
    tab.delete_rows_func(constraint);
    size_t nlines_new=tab.get_nlines();
    if (verbose>0) {
      cout << "Applied constraint \"" << constraint
	   << "\" and went from " << nlines_old << " to "
	   << nlines_new << " lines." << endl;
    }
  }
  
  if (verbose>0) {
    cout << "Table " << tab_name << " lines: " 
	 << tab.get_nlines();
  }
  
  if (verbose>0) {
    cout << "Done with reading files.\n" << endl;
  }

  // ------------------------------------------------------------

  size_t n_rows_old=tab.get_nlines();
  table_units<> t_new;
  t_new.new_column(grid_name);
  for(size_t i=0;i<n_rows_old;i++) {
    t_new.new_column(prefix+"_"+o2scl::szttos(i));
  }
  uniform_grid_end<double> ug(low,high,99);

  for(size_t j=0;j<100;j++) {
    vector<double> line;
    line.push_back(ug[j]);
    for(size_t k=0;k<n_rows_old;k++) {
      line.push_back(tab.get(prefix+"_"+o2scl::szttos(j),k));
    }
    t_new.line_of_data(line.size(),line);
  }

  hf.open(infile);
  tab_name="markov_chain_0";
  hdf_output(hf,t_new,tab_name);
  hf.close();
  
  return 0;
}

int process::hist2(std::vector<std::string> &sv, bool itive_com) {
    
  // Setup histogram size
  size_t hist_size=(size_t)hist_size_int;
  if (hist_size_int<=0) {
    cout << "Histogram size <=0, using 100." << endl;
    hist_size=100;
  }
    
  // Column names
  string xname=sv[1];
  string yname=sv[2];

  // (output file is sv[3])

  // List of data files
  vector<string> files;
  for(size_t i=4;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
    
  // Ouptut table
  table_units<> tab;

  // Temporary column storage
  vector<double> valuesx, valuesy, weights;

  // data
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;
    cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    cout << n_chains << " chains." << endl;

    size_t line_counter=0;

    for(size_t j=0;j<n_chains;j++) {

      std::string tab_name="markov_chain_"+szttos(j);
      hdf_input(hf,tab,tab_name);
      cout << "Read table " << tab_name << endl;

      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows_func(constraint);
	size_t nlines_new=tab.get_nlines();
	if (verbose>0) {
	  cout << "Applied constraint \"" << constraint
	       << "\" and went from " << nlines_old << " to "
	       << nlines_new << " lines." << endl;
	}
      }

      for(size_t k=0;k<tab.get_nlines();k++) {
	if (((int)line_counter)>=line_start) {
	  valuesx.push_back(tab.get(xname,k));
	  valuesy.push_back(tab.get(yname,k));

	  if (weights_col.size()>0) {
	    weights.push_back(tab.get(weights_col,k));
	  } else {
	    weights.push_back(1.0);
	  }
	  
	}
	line_counter++;
      }
    }
  }

  // Form x and y limits
  double xmin, xmax;
  if (!xset) {
    xmin=*min_element(valuesx.begin(),valuesx.end());
    xmax=*max_element(valuesx.begin(),valuesx.end());
    xmin*=xscale;
    xmax*=xscale;
    cout << xname << " xsize, xmin, xmax: " << valuesx.size() << " " 
	 << xmin << " " << xmax << endl;
  } else {
    xmin=user_xlow;
    xmax=user_xhigh;
  }

  double ymin, ymax;
  if (!yset) {
    ymin=*min_element(valuesy.begin(),valuesy.end());
    ymax=*max_element(valuesy.begin(),valuesy.end());
    ymin*=yscale;
    ymax*=yscale;
    cout << yname << " ysize, ymin, ymax: " << valuesy.size() << " " 
	 << ymin << " " << ymax << endl;
  } else {
    ymin=user_ylow;
    ymax=user_yhigh;
  }

  // Create expval_scalar objects
  if (n_blocks<=0) n_blocks=20;
  vector<expval_scalar> sev(hist_size*hist_size);
  for(size_t i=0;i<hist_size*hist_size;i++) {
    sev[i].set_blocks(n_blocks,1);
  }

  // Grid and histogram
  uniform_grid_end<double> xgrid(xmin*(1.0-1.0e-6),
				 xmax*(1.0+1.0e-6),hist_size);
  uniform_grid_end<double> ygrid(ymin*(1.0-1.0e-6),
				 ymax*(1.0+1.0e-6),hist_size);
  hist_2d h;
  h.set_bin_edges(xgrid,ygrid);
    
  // Fill expval_scalar objects using histogram
  size_t block_size=valuesx.size()/20;
  for(size_t i=0;((int)i)<n_blocks;i++) {
    cout << "Block " << i << endl;
    for(size_t j=0;j<block_size;j++) {
      double vx=valuesx[i*block_size+j]*xscale;
      double vy=valuesy[i*block_size+j]*yscale;
      if (vx<xmax && vx>xmin && vy<ymax && vy>ymin) {
	h.update(vx,vy,weights[i*block_size+j]);
      }
    }
    for(size_t j=0;j<hist_size;j++) {
      for(size_t k=0;k<hist_size;k++) {
	sev[j*hist_size+k].add(h.get_wgt_i(j,k));
      }
    }
    h.clear_wgts();
  }
    
  // Create x and y grids
  ubvector xreps(hist_size), yreps(hist_size);
  for(size_t i=0;i<hist_size;i++) {
    xreps[i]=h.get_x_rep_i(i);
    yreps[i]=h.get_y_rep_i(i);
  }

  // Create table
  table3d t3d;
  t3d.set_xy(xname,xreps.size(),xreps,yname,yreps.size(),yreps);
  t3d.new_slice("avgs");
    
  // Collect averages into table
  double std_dev, avg_err, avg, smax=0.0;
  for(size_t i=0;i<hist_size;i++) {
    for(size_t j=0;j<hist_size;j++) {
      sev[i*hist_size+j].current_avg(avg,std_dev,avg_err);
      if (avg>smax) smax=avg;
      t3d.set(i,j,"avgs",avg);
    }
  }

  // Renormalize so that maximum value is 1.0 (important for
  // contours below)
  for(size_t i=0;i<hist_size;i++) {
    for(size_t j=0;j<hist_size;j++) {
      t3d.set(i,j,"avgs",t3d.get(i,j,"avgs")/smax);
    }
  }

  cout << "Writing to file " << sv[3] << endl;
  hdf_file hf;
  hf.open_or_create(sv[3]);
  hdf_output(hf,(const table3d &)t3d,"hist2_table");

  // ------------------------------------------
  // Construct x and y values for interpolating the double integral
  
  ubvector integx(101), integy(101);
  for(size_t i=0;i<101;i++) {
    integx[i]=((double)i)/100.0;
  }
  for(size_t k=0;k<101;k++) {
    integy[k]=0.0;
    for(size_t i=0;i<hist_size;i++) {
      for(size_t j=0;j<hist_size;j++) {
	if (t3d.get(i,j,"avgs")>integx[k]) integy[k]+=t3d.get(i,j,"avgs");
      }
    }
  }

  // Get the total 
  double total=integy[0];

  // Create contour object
  contour co;
  
  // Set the contour object data
  ubvector xg(hist_size), yg(hist_size);
  for(size_t i=0;i<hist_size;i++) {
    xg[i]=t3d.get_grid_x(i);
    yg[i]=t3d.get_grid_y(i);
  }
  co.set_data(hist_size,hist_size,xg,yg,t3d.get_slice("avgs"));
    
  // ------------------------------------------------------------------
  // Compute contour lines for all of the requested levels. We could
  // compute all the levels at once with the contour class, but here
  // we separate the calculation so that it is easy to separate the
  // contour lines for each requested level.
  
  for(size_t ic=0;ic<cont_levels.size();ic++) {
    
    string sig_level;

    if (cont_levels[ic]==one_sigma) {
      sig_level="1s";
    } else if(cont_levels[ic]==two_sigma) {
      sig_level="2s";
    } else if(cont_levels[ic]==three_sigma) {
      sig_level="3s";
    } else {
      sig_level=o2scl::szttos(ic);
    }

    vector<contour_line> conts;
    size_t nc=0;
      
    ubvector levels(1);
    
    // Compute the function value associated with the contour levels
    levels[0]=0.0;
    for(size_t i=0;i<100;i++) {
      if (integy[i]>cont_levels[ic]*total &&
	  integy[i+1]<cont_levels[ic]*total) {
	levels[0]=(integx[i]+integx[i+1])/2.0;
	i=100;
      }
    }
    if (levels[0]==0.0) {
      cout << "Failed to find any contour levels for level " 
  	   << cont_levels[ic] << " and function value " << levels[0] << endl;
    }

    // Set the levels
    co.set_levels(levels.size(),levels);
    
    // Compute the contours
    co.calc_contours(conts);
    nc=conts.size();
    cout << "Number of contours: " << nc << endl;
    
    // If those levels were found, write the associated contours
    // to the HDF5 file
    hf.setd(((string)"level_")+sig_level,levels[0]);
    hf.set_szt(((string)"n_contours_")+sig_level,nc);
    if (nc>0) {
      hdf_output(hf,conts,((string)"contours_")+sig_level);
    }
  }
  hf.close();

  return 0;
}

int process::hist_set(std::vector<std::string> &sv, bool itive_com) {
    
  // Setup histogram size
  size_t hist_size=(size_t)hist_size_int;
  if (hist_size_int<=0) {
    cout << "Histogram size <=0, using 100." << endl;
    hist_size=100;
  }

  string type=sv[1];
  double index_low, index_high;
  o2scl::stod_nothrow(sv[2],index_low);
  o2scl::stod_nothrow(sv[3],index_high);
  string set_prefix=sv[4];
  int renorm=o2scl::stoi(sv[5]);
  // (output file is sv[6])
  // (input files are sv[7] ... )

  // file list
  vector<string> files;

  // Form list of data files
  for(size_t i=7;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
      
  // The value of 'grid_size' from the bamr output file
  size_t grid_size=100;

  // Storage for all of the data
  vector<double> weights;
  vector<vector<double> > values_ser;

  // The grid determined by the boundaries specified in low_name
  // and high_name by the user (unscaled)
  uniform_grid<double> index_grid;

  // Count the data over each grid
  ubvector count;

  // Series min and series max, just determined by the min and
  // max values over all points in all chains
  double ser_min=0.0, ser_max=0.0;

  if ((logx && type==((string)"x")) ||
      (logy && type==((string)"y"))) {
    index_grid=uniform_grid_log_end<double>
      (index_low,index_high,grid_size-1);
  } else {
    index_grid=uniform_grid_end<double>
      (index_low,index_high,grid_size-1);
  }
  
  // ------------------------------------------------------------
  // Read the data
      
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;
    if (verbose>0) {
      cout << "Opening file: " << files[i] << endl;
    }
    hf.open(files[i]);

    // If we're reading the first file, obtain the grid information
    if (i==0) {
      values_ser.resize(grid_size);
      count.resize(index_grid.get_npoints());
    }

    // Count the total number of lines over all
    // the individual tables in the file
    size_t line_counter=0;

    // Read chain from file
    table_units<> tab;
    std::string tab_name="";
    hdf_input(hf,tab,tab_name);
    if (verbose>0) {
      cout << "Read table with "
	   << tab.get_nlines() << " lines." << endl;
    }
    
    // Apply constraint to this table
    if (constraint.size()>0) {
      size_t nlines_old=tab.get_nlines();
      tab.delete_rows_func(constraint);
      size_t nlines_new=tab.get_nlines();
      if (verbose>0) {
	cout << "Applied constraint \"" << constraint
	     << "\" and went from " << nlines_old << " to "
	     << nlines_new << " lines." << endl;
      }
    }
    
    // Read each line in the current chain
    for(size_t k=0;k<tab.get_nlines();k++) {
      
      if (((int)line_counter)>=line_start) {
	
	// For each column in the series
	for(size_t ell=grid_init;ell<grid_size;ell++) {
	  
	  // Column name and value
	  string col=set_prefix+"_"+szttos(ell);
	  double val=tab.get(col,k);
	  
	  // Add the value even if it's zero (blank)
	  values_ser[ell].push_back(val);
	  
	  // But if it's less than or equal to zero, don't count it
	  // for min, max, or count. This is important because we
	  // don't want to count past the M-R curve or the end of the
	  // EOS.
	  if (val>0.0) {
	    count[ell]++;
	    // The first time, initialize to the first non-zero point
	    if (ser_min==0.0) {
	      ser_min=val;
	      ser_max=val;
	    }
	    if (val<ser_min) ser_min=val;
	    if (val>ser_max) ser_max=val;
	  }
	  // Next column
	}
	if (weights_col.size()>0) {
	  weights.push_back(tab.get(weights_col,k));
	} else {
	  weights.push_back(1.0);
	}
      }
      line_counter++;
      
      // Next table line (k)
    }

    if (verbose>0) {
      cout << "Read file " << files[i] << "\n  and found min, max: "
	   << ser_min << " " << ser_max << " (not yet rescaled) " << endl;
    }
    
    // Next file (i)
  }
  if (verbose>0) {
    cout << "Done reading files." << endl;
  }

  // Rescale series limits, and set from user-specified values
  // if present
  if (type==((string)"x")) {
    ser_min*=yscale;
    ser_max*=yscale;
    if (verbose>0 && fabs(yscale-1.0)<1.0e-6) {
      cout << "Rescaling y by " << yscale << endl;
    }
  } else {
    ser_min*=xscale;
    ser_max*=xscale;
    if (verbose>0 && fabs(xscale-1.0)<1.0e-6) {
      cout << "Rescaling x by " << xscale << endl;
    }
  }
  if (type==((string)"x")) {
    if (yset) {
      ser_min=user_ylow;
      ser_max=user_yhigh;
      cout << "Set y limits to " << ser_min << " and " << ser_max << endl;
    }
  } else {
    if (xset) {
      ser_min=user_xlow;
      ser_max=user_xhigh;
      cout << "Set x limits to " << ser_min << " and " << ser_max << endl;
    }
  }

  if (verbose>0) {
    cout << "Total lines: " << weights.size() << " " 
	 << values_ser[0].size() << " ser_min: " << ser_min << " ser_max: " 
	 << ser_max << endl;
  }

  // Create expval_scalar objects
  if (n_blocks<=0) {
    n_blocks=20;
  }
  cout << "Using " << n_blocks << " blocks." << endl;
  vector<expval_scalar> sev(hist_size*grid_size);
  for(size_t i=0;i<hist_size*grid_size;i++) {
    sev[i].set_blocks(n_blocks,1);
  }

  uniform_grid<double> ser_grid;
  if ((logy && type==((string)"x")) ||
      (logx && type==((string)"y"))) {
    ser_grid=uniform_grid_log_end<double>(ser_min,ser_max,hist_size);
  } else {
    ser_grid=uniform_grid_end<double>(ser_min,ser_max,hist_size);
  }

  // Histograms to compute contour lines. The histogram 'h' stores
  // the values for the current block, and 'hsum' stores the
  // values over all blocks.
  o2scl::hist h, hsum;
  h.set_bin_edges(ser_grid);
  hsum.set_bin_edges(ser_grid);

  // Vectors to store contour results. Each additional contour level
  // requires two vectors: a low and high vector.
  table<> cont_res;
  cont_res.new_column("grid");
  for(size_t j=0;j<cont_levels.size()*2;j++) {
    cont_res.new_column(((std::string)"c")+o2scl::szttos(j));
  }

  // Fill expval_scalar objects using histogram
  size_t block_size=weights.size()/20;
  for(size_t k=grid_init;k<grid_size;k++) {

    string col=set_prefix+"_"+szttos(k);
    cout.precision(3);
    cout << "Col ";
    cout.width(4);
    cout << col << " ix: ";
    if (type==((string)"x")) {
      cout << index_grid[k]*xscale << " cnt: ";
    } else {
      cout << index_grid[k]*yscale << " cnt: ";
    }
    cout << count[k] << " ";

    // Add up the entries for the histograms
    for(size_t i=0;((int)i)<n_blocks;i++) {
      for(size_t j=0;j<block_size;j++) {
	double val=values_ser[k][i*block_size+j];
	if (type==((string)"x")) val*=yscale;
	else val*=xscale;
	if (val>ser_min && val<ser_max) {
	  h.update(val,weights[i*block_size+j]);
	  hsum.update(val,weights[i*block_size+j]);
	}
      }
      // If requested, renormalize the histogram so that the maximum
      // value is 1
      if (renorm>=1) {
	double max=h.get_max_wgt();
	if (max==0.0) {
	  std::cout << std::endl;
	  std::cout << "Histogram has no entries." << std::endl;
	}
	for(size_t j=0;j<hist_size;j++) {
	  h.set_wgt_i(j,h.get_wgt_i(j)/max);
	}
      }
      for(size_t j=0;j<hist_size;j++) {
	sev[j*grid_size+k].add(h.get_wgt_i(j));
      }
      h.clear_wgts();
    }

    // Copy results from hsum histogram
    vector<double> cont_wgt, cont_grid;
    for(size_t i=0;i<hsum.size();i++) {
      cont_wgt.push_back(hsum.get_wgt_i(i));
      cont_grid.push_back(hsum.get_rep_i(i));
    }
	
    // Compute contour levels
    for(size_t j=0;j<cont_levels.size();j++) {
      if (count[k]>0.0) {
	std::vector<double> locs;
	cont_wgt[0]=0.0;
	cont_wgt[cont_wgt.size()-1]=0.0;
	if (o2scl::vector_sum_double(cont_wgt)>0.0) {
	  vector_region_fracint(cont_wgt.size(),cont_grid,cont_wgt,
				cont_levels[j],locs);
	  if (locs.size()==0) {
	    cout << "<failed> ";
	    if (false) {
	      cout << cont_levels[j] << endl;
	      for(size_t i=0;i<cont_wgt.size();i++) {
		cout << i << " " << cont_grid[i] << " "
		     << cont_wgt[i] << endl;
	      }
	    }
	  } else {
	    double lmin=vector_min_value<vector<double>,double>(locs);
	    double lmax=vector_max_value<vector<double>,double>(locs);
	    cout << lmin << " " << lmax << " ";
	    if (cont_res.get_nlines()<=k) {
	      cont_res.set_nlines(k+1);
	      if (type==((string)"x")) {
		cont_res.set("grid",k,index_grid[k]*xscale);
	      } else {
		cont_res.set("grid",k,index_grid[k]*yscale);
	      }
	    }
	    cont_res.set(2*j+1,k,lmin);
	    cont_res.set(2*j+2,k,lmax);
	  }
	}
      }
    }
    
    hsum.clear_wgts();
    cout << endl;
    cout.precision(6);

    // End of loop over grid index
  }
      
  // Create both grids for the table3d object
  ubvector reps(hist_size);
  for(size_t i=0;i<hist_size;i++) {
    reps[i]=h.get_rep_i(i);
  }
  ubvector index_grid_vec;
  index_grid_vec.resize(index_grid.get_npoints());
  vector_grid(index_grid,index_grid_vec);
  if (type==((string)"x")) {
    for(size_t i=0;i<index_grid_vec.size();i++) {
      index_grid_vec[i]*=xscale;
    }
  } else {
    for(size_t i=0;i<index_grid_vec.size();i++) {
      index_grid_vec[i]*=yscale;
    }
  }

  // Create the table3d object
  table3d t3d;
  if (type==((string)"x")) {
    t3d.set_xy("x",index_grid_vec.size(),index_grid_vec,
	       "y",reps.size(),reps);
  } else {
    t3d.set_xy("x",reps.size(),reps,
	       "y",index_grid_vec.size(),index_grid_vec);
  }
  t3d.new_slice("avgs");
    
  // Collect averages into the table3d slice
  double std_dev, avg_err, avg, smax=0.0;
  for(size_t j=0;j<hist_size;j++) {
    for(size_t k=grid_init;k<grid_size;k++) {
      sev[j*grid_size+k].current_avg(avg,std_dev,avg_err);
      if (avg>smax) smax=avg;
      if (type==((string)"x")) {
	t3d.set(k,j,"avgs",avg);
      } else {
	t3d.set(j,k,"avgs",avg);
      }
    }
  }

  // Renormalize by the maximum value
  for(size_t j=0;j<t3d.get_nx();j++) {
    for(size_t k=0;k<t3d.get_ny();k++) {
      t3d.set(j,k,"avgs",t3d.get(j,k,"avgs")/smax);
    }
  }

  // Perform file output
  cout << "Writing table to file " << sv[6] << endl;
  hdf_file hf;
  hf.open_or_create(sv[6]);
  hdf_output(hf,(const table3d &)t3d,"hist_set");
  hdf_output(hf,index_grid,"index_grid");
  hdf_output(hf,cont_res,"cont_res");
  hf.close();

  return 0;
}

int process::contours(std::vector<std::string> &sv, bool itive_com) {
  cont_levels.clear();
  for(size_t i=1;i<sv.size();i++) {
    if (sv[i]==((string)"1s")) {
      cont_levels.push_back(one_sigma);
    } else if (sv[i]==((string)"2s")) {
      cont_levels.push_back(two_sigma);
    } else if (sv[i]==((string)"3s")) {
      cont_levels.push_back(three_sigma);
    } else {
      cont_levels.push_back(o2scl::function_to_double(sv[i]));
    }
  }
  return 0;
}

int process::combine(std::vector<std::string> &sv, bool itive_com) {
    
  // Thinning factor
  size_t thin_factor=stoszt(sv[1]);
  if (thin_factor==0) thin_factor=1;

  // Output file
  string out_file=sv[2];

  // file list
  vector<string> files;

  // Form list of data files
  for(size_t i=3;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
  
  table_units<> full;
  
  /*
    double e_low, e_high, m_low, m_high, nb_low, nb_high;
    ubvector low, high;
    size_t nparams, nsources, lgrid_size, file_row_start;
    vector<string> param_names, source_names;
  */
  //size_t file_row_start;
  
  if (line_start>0) {
    cout << "Skipping " << line_start
	 << " lines at the beginning of every file." << endl;
  }
  
  // data
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;
    
    /*
    // Get number of chains
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    if (n_chains>1) {
    cout << n_chains << " separate chains." << endl;
    } else {
    cout << "1 chain." << endl;
    }
    */
    
    size_t line_counter=0;
    
    //file_row_start=full.get_nlines();
    
    // Read each chain
    //for(size_t j=0;j<n_chains;j++) {
    
    // Read table
    cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);
    std::string tab_name="markov_chain_0";
    table_units<> tab;
    hdf_input(hf,tab,tab_name);
    hf.close();
    
    // Set up columns in fill table over all files
    if (i==0) {
      for(size_t ell=0;ell<tab.get_ncolumns();ell++) {
	full.new_column(tab.get_column_name(ell));
      }
    }
    
    // Get misc parameters
    /*
      if (i==0 && j==0) {
      hf.getd("e_low",e_low);
      hf.getd("e_high",e_high);
      hf.getd("nb_low",nb_low);
      hf.getd("nb_high",nb_high);
      hf.getd("m_low",m_low);
      hf.getd("m_high",m_high);
      hf.getd_vec_copy("low",low);
      hf.getd_vec_copy("high",high);
      hf.get_szt("grid_size",lgrid_size);
      hf.get_szt("n_params",nparams);
      hf.get_szt("n_sources",nsources);
      //hf.gets_vec("param_names",param_names);
      hf.gets_vec("source_names",source_names);
      }
    */
    
    // Add data to table for this file
    
    /*
      if (j==0) {
      
      // For the first chain, we want to skip according to line_start
      for(size_t k=((size_t)line_start);k<tab.get_nlines();k+=thin_factor) {
      vector<double> line;
      for(size_t ell=0;ell<tab.get_ncolumns();ell++) {
      line.push_back(tab.get(tab.get_column_name(ell),k));
      }
      full.line_of_data(line.size(),line);
      }
      
      cout << "Added table " << tab_name << " lines: " 
      << tab.get_nlines();
      if (line_start>0) {
      cout << " skipped: " << line_start;
      }
      cout << endl;
      
      } else {
    */
    
    // For the remaining chains we don't skip any, we just
    // start at the first row
    for(size_t k=line_start;k<tab.get_nlines();k+=thin_factor) {
      vector<double> line;
      for(size_t ell=0;ell<tab.get_ncolumns();ell++) {
	line.push_back(tab.get(tab.get_column_name(ell),k));
      }
      full.line_of_data(line.size(),line);
    }
    
    cout << "Added table " << tab_name << " lines: " 
	 << tab.get_nlines() << endl;
    //}
    
    // Go to next chain
    //}
    
    /*	
    // Load weight column for this file from full table
    vector<double> subweights;
    for(size_t k=file_row_start;k<full.get_nlines();k++) {
    subweights.push_back(exp(full.get("log_wgt",k)));
    }
    
    // Compute autocorrelation
    size_t kmax=subweights.size()/3;
    ubvector kvec(kmax-1), acvec(kmax-1), acerr(kmax-1);
    for(size_t k=1;k<kmax;k++) {
    kvec[k-1]=((double)k);
    acvec[k-1]=vector_lagk_autocorr(subweights.size(),subweights,k);
    acerr[k-1]=fabs(acvec[k-1])/10.0;
    }
    */
    
    cout << "Current total number of lines: "
	 << full.get_nlines() << endl;
    
    // Go to next file
  }
  
  cout << "Writing table to file " << out_file << endl;
  hdf_file hf;
  hf.open_or_create(out_file);
  string type;
  if (hf.find_object_by_name("markov_chain_0",type)==0) {
    O2SCL_ERR("Refusing to overwrite file with data in process 'combine'.",
	      exc_efailed);
  }
  /*
  hf.setd("e_low",e_low);
  hf.setd("e_high",e_high);
  hf.setd("nb_low",nb_low);
  hf.setd("nb_high",nb_high);
  hf.setd("m_low",m_low);
  hf.setd("m_high",m_high);
  hf.set_szt("grid_size",lgrid_size);
  hf.set_szt("n_chains",1);
  hf.setd_vec_copy("low",low);
  hf.setd_vec_copy("high",high);
  hf.set_szt("n_params",nparams);
  hf.set_szt("n_sources",nsources);
  //hf.sets_vec("param_names",param_names);
  hf.sets_vec("source_names",source_names);
  */
  hdf_output(hf,full,"markov_chain_0");
  hf.close();

  return 0;
}
    
void process::setup_cli() {

    // ---------------------------------------
    // Setup CLI readline history

#ifdef O2SCL_READLINE
    char *hd=getenv("HOME");
    std::string histfile;
    if (hd) {
      histfile=((std::string)hd)+"/.process_hist";
      cl.set_histfile(histfile);
    }
#endif
    
  // ---------------------------------------
  // Set options
  
  static const int nopt=12;
  comm_option_s options[nopt]={
    {'x',"xlimits","Set histogram limits for first variable",0,3,
     "<low-value high-value> or <file> <low-name> <high-name> or <>",
     ((string)"In the case of no arguments, set the limits ")+
     "to the default. In the case of two arguments, take the two "+
     "specified numbers as the limits, in the case of three "+
     "arguments, assume that the first argument is the bamr output "+
     "file, and the second and third are the names of values "+
     "in the output file which should be set as limits.",
     new comm_option_mfptr<process>(this,&process::xlimits),
     cli::comm_option_both},
    {'y',"ylimits","Set histogram limits for the second variable",0,3,
     "<low-value high-value> or <file> <low-name> <high-name> or <>",
     ((std::string)"In the case of no arguments, set the limits ")+
     "to the default. In the case of two arguments, take the two "+
     "specified numbers as the limits, in the case of three "+
     "arguments, assume that the first argument is the bamr output "+
     "file, and the second and third are the names of values "+
     "in the output file which should be set as limits.",
     new comm_option_mfptr<process>(this,&process::ylimits),
     cli::comm_option_both},
    {'C',"combine","Combine two bamr output files.",2,-1,
     "<thin factor> <output filename> <file1> <file2> [file 3...]",
     ((string)"Long ")+"desc.",
     new comm_option_mfptr<process>(this,&process::combine),
     cli::comm_option_both},
    {0,"contours","Specify contour levels",-1,-1,
     "[level 1] [level 2] ...",
     ((string)"Specify a list of contour levels for the 'hist2' and ")+
     "'hist-set' commands. The levels can either be numerical values "+
     "or they can be 1s, 2s, or 3s, corresponding to 1-, 2-, and 3-"+
     "sigma confidence limits, respectively.",
     new comm_option_mfptr<process>(this,&process::contours),
     cli::comm_option_both},
    {0,"mass-sel","Create new columns from a specified mass.",3,-1,
     "<mass> <suffix> <file 1> [file 2] ...",
     ((string)"Long ")+"desc.",
     new comm_option_mfptr<process>(this,&process::mass_sel),
     cli::comm_option_both},
    {0,"hist2","Create a histogram from two columns of MCMC data.",4,-1,
     "<x> <y> <out_file> <file1> [file2 file 3...]",
     ((string)"Using the 'markov_chain' objects in a bamr output ")+
     "file, construct a two-dimensional histogram from the "+
     "specified columns and store the histogram in 'out_file'.",
     new comm_option_mfptr<process>(this,&process::hist2),
     cli::comm_option_both},
    {0,"hist","Create a histogram from one column of MCMC data.",3,-1,
     "<column> <out_file> <file1> [file2 file 3...]",
     ((string)"Using the 'markov_chain' objects in a bamr output ")+
     "file, construct a histogram from the specified column, and "+
     "store the histogram in 'out_file'.",
     new comm_option_mfptr<process>(this,&process::hist),
     cli::comm_option_both},
    {0,"auto-corr","Create a table of autocorrelation data",3,-1,
     "<column> <out_file> <file1> [file2 file 3...]",
     ((string)"Using the 'markov_chain' objects in a bamr output ")+
     "file, create a table containing autocorrelation "+
     "information about the specified column, and "+
     "store the table in 'out_file'.",
     new comm_option_mfptr<process>(this,&process::auto_corr),
     cli::comm_option_both},
    {0,"hist-set","Create an ensemble of of 1-d histograms.",6,-1,
     "<direction> <low> <high> <set> <renorm> <out_file> <file1> [file2...]",
     ((string)"Using the 'markov_chain' objects in a bamr output ")+
     "file, create an ensemble of 1-d histograms from a set "+
     "of columns in each chain. Typical uses are: "+
     "\'process -hist-set y m_low m_high R 0 out.o2 x_0_out\', "+
     "\'process -hist-set x e_low e_high P 0 out.o2 x_0_out\', and "+
     "\'process -hist-set x nb_low nb_high P 0 out.o2 x_0_out\'. ",
     new comm_option_mfptr<process>(this,&process::hist_set),
     cli::comm_option_both},
    {0,"curve-set","",6,-1,
     "<in file> <grid name> <low> <high> <prefix> <out file>","",
     new comm_option_mfptr<process>(this,&process::curve_set),
     cli::comm_option_both},
    {'b',"bfactor","Compute the Bayes factor.",1,-1,
     "<x name> <y name> <file1> [file2 file 3...]",
     ((string)"Long ")+"desc.",
     new comm_option_mfptr<process>(this,&process::bfactor),
     cli::comm_option_both},
    {0,"set-params-limits","Set parameters and limits",-1,-1,
     "<1: name low high> <2: name low high> ... <n: name low high>",
     ((string)"Long ")+"desc.",
     new comm_option_mfptr<process>(this,&process::set_params_limits),
     cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);

  // ---------------------------------------
  // Set parameters
    
  p_xscale.d=&xscale;
  p_xscale.help="Scale parameter for first variable.";
  cl.par_list.insert(make_pair("xscale",&p_xscale));
    
  p_yscale.d=&yscale;
  p_yscale.help="Scale parameter for second variable.";
  cl.par_list.insert(make_pair("yscale",&p_yscale));

  p_errors.b=&errors;
  p_errors.help="If true, output more error information.";
  cl.par_list.insert(make_pair("errors",&p_errors));

  p_logx.b=&logx;
  p_logx.help="If true, use a logarithmic x-axis (default false).";
  cl.par_list.insert(make_pair("logx",&p_logx));

  p_logy.b=&logy;
  p_logy.help="If true, use a logarithmic y-axis (default false).";
  cl.par_list.insert(make_pair("logy",&p_logy));

  p_logz.b=&logz;
  p_logz.help="If true, use a logarithmic z-axis (default false).";
  cl.par_list.insert(make_pair("logz",&p_logz));

  p_hist_size.i=&hist_size_int;
  p_hist_size.help="Histogram size (default 100).";
  cl.par_list.insert(make_pair("hist_size",&p_hist_size));

  p_n_blocks.i=&n_blocks;
  p_n_blocks.help="Number of blocks (default 1).";
  cl.par_list.insert(make_pair("n_blocks",&p_n_blocks));

  p_grid_init.i=&grid_init;
  p_grid_init.help="Initial column number for grid (default 0).";
  cl.par_list.insert(make_pair("grid_init",&p_grid_init));

  p_line_start.i=&line_start;
  p_line_start.help=((string)"Number of initial rows to skip over ")+
    "when reading MCMC tables.";
  cl.par_list.insert(make_pair("line_start",&p_line_start));

  p_verbose.i=&verbose;
  p_verbose.help="Verbosity parameter (default 1).";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  p_constraint.str=&constraint;
  p_constraint.help="Constraint to apply (default is an empty string).";
  cl.par_list.insert(make_pair("constraint",&p_constraint));

  p_weights.str=&weights_col;
  p_weights.help="Column with weights.";
  cl.par_list.insert(make_pair("weights",&p_weights));

  cl.prompt="process> ";

  return;
}

void process::run(int argc, char *argv[]) {
  
  // ---------------------------------------
  // Process command-line arguments and run

  setup_cli();
  
  cl.run_auto(argc,argv);

  return;

}

