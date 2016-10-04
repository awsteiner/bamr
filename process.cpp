/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
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

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

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
      std::string tab_name="markov_chain"+szttos(j);
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

/** \brief Set limits for the x-axis
 */
int process::xlimits(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<3) {
    if (verbose>0) {
      cout << "Setting 'xset' to false." << endl;
    }
    xset=false;
    return 0;
  }
    
  if (sv.size()==3) {
    user_xlow=o2scl::function_to_double(sv[1]);
    user_xhigh=o2scl::function_to_double(sv[2]);
    xset=true;
    if (verbose>0) {
      cout << "X limits are " << user_xlow << " and " 
	   << user_xhigh << " ." << endl;
    }
    return 0;
  }

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

/** \brief Set limits for the y-axis
 */
int process::ylimits(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<3) {
    if (verbose>0) {
      cout << "Setting 'yset' to false." << endl;
    }
    yset=false;
    return 0;
  }

  if (sv.size()==3) {
    user_ylow=o2scl::function_to_double(sv[1]);
    user_yhigh=o2scl::function_to_double(sv[2]);
    yset=true;
    if (verbose>0) {
      cout << "Y limits are " << user_ylow << " and " 
	   << user_yhigh << " ." << endl;
    }
    return 0;
  }

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
  
/** \brief Create a histogram from a specified column

    \future This function isn't that efficient. It first reads all
    of the data into one long vector and then reparses the long
    vector into a set of \ref o2scl::expval_scalar objects. The
    averages and errors for the \ref o2scl::expval_scalar objects
    are stored into a table, and the table is written to the
    output file. This could be improved.
*/
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
      std::string tab_name="markov_chain"+szttos(j);
      table_units<> tab;
      hdf_input(hf,tab,tab_name);

      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows(constraint);
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
	  weights.push_back(tab.get("mult",k));
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
  double total=vector_integ_linear(hist_size,t["reps"],t["avgs"]);
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

/** \brief Create a two-dimensional histogram from two
    user-specified columns
*/
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

      std::string tab_name="markov_chain"+szttos(j);
      hdf_input(hf,tab,tab_name);
      cout << "Read table " << tab_name << endl;

      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows(constraint);
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
	  weights.push_back(tab.get("mult",k));
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

  // ------------------------------------------
  // Compute contour lines
    
  vector<contour_line> conts;
  size_t nc=0;
    
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
  ubvector levels(cont_levels.size());

  // Compute the function values associated with the contour levels
  for(size_t k=0;k<cont_levels.size();k++) {
    levels[k]=0.0;
    for(size_t i=0;i<100;i++) {
      if (integy[i]>cont_levels[k]*total && 
	  integy[i+1]<cont_levels[k]*total) {
	levels[k]=(integx[i]+integx[i+1])/2.0;
	i=100;
      }
    }
    if (levels[k]==0.0) {
      cout << "Failed to find contour level for level " 
	   << levels[k] << endl;
    }
  }

  // If those levels were found, plot the associated contours
  contour co;
      
  // Set the contour object data
  ubvector xg(hist_size), yg(hist_size);
  for(size_t i=0;i<hist_size;i++) {
    xg[i]=t3d.get_grid_x(i);
    yg[i]=t3d.get_grid_y(i);
  }
  co.set_data(hist_size,hist_size,xg,yg,t3d.get_slice("avgs"));
      
  // Set the levels
  co.set_levels(levels.size(),levels);
      
  // Compute the contours
  co.calc_contours(conts);
  nc=conts.size();
  cout << "Number of contours: " << nc << endl;

  cout << "Writing to file " << sv[3] << endl;
  hdf_file hf;
  hf.open_or_create(sv[3]);
  hdf_output(hf,t3d,"hist2_table");
  hf.setd_vec_copy("levels",levels);
  hf.set_szt("n_contours",nc);
  if (nc>0) {
    hdf_output(hf,conts,"contours");
  }
  hf.close();

  return 0;
}

/** \brief Create a set of histograms from a set of columns in
    the bamr MCMC output
*/
int process::hist_set(std::vector<std::string> &sv, bool itive_com) {
    
  // Setup histogram size
  size_t hist_size=(size_t)hist_size_int;
  if (hist_size_int<=0) {
    cout << "Histogram size <=0, using 100." << endl;
    hist_size=100;
  }

  string type=sv[1];
  string low_name=sv[2];
  string high_name=sv[3];
  string set_prefix=sv[4];

  // (output file is sv[5])

  // file list
  vector<string> files;

  // Form list of data files
  for(size_t i=6;i<sv.size();i++) files.push_back(sv[i]);
  size_t nf=files.size();
      
  // The value of 'grid_size' from the bamr output file
  size_t grid_size=0;

  // Storage for all of the data
  vector<double> weights;
  vector<vector<double> > values_ser;

  // The grid determined by the boundaries specified in low_name
  // and high_name by the user (unscaled)
  double index_low, index_high;
  uniform_grid<double> index_grid;

  // Count the data over each grid
  ubvector count;

  // Series min and series max, just determined by the min and
  // max values over all points in all chains
  double ser_min=0.0, ser_max=0.0;

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
      hf.get_szt("grid_size",grid_size);
      hf.getd(low_name,index_low);
      hf.getd(high_name,index_high);
      if (verbose>0) {
	cout << "grid_size, index_low, index_high: " 
	     << grid_size << " " << index_low << " " << index_high << endl;
      }
      if ((logx && type==((string)"x")) ||
	  (logy && type==((string)"y"))) {
	index_grid=uniform_grid_log_end<double>
	  (index_low,index_high,grid_size-1);
      } else {
	index_grid=uniform_grid_end<double>
	  (index_low,index_high,grid_size-1);
      }
      values_ser.resize(grid_size);
      count.resize(index_grid.get_npoints());
    }

    // Obtain the number of chains in this file
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    if (verbose>0) {
      if (n_chains==1) {
	cout << n_chains << " separate chain." << endl;
      } else {
	cout << n_chains << " separate chains." << endl;
      }
    }

    // Count the total number of lines over all
    // the individual tables in the file
    size_t line_counter=0;

    // Process each chain in turn
    for(size_t j=0;j<n_chains;j++) {
	  
      // Read chain from file
      table_units<> tab;
      std::string tab_name="markov_chain"+szttos(j);
      hdf_input(hf,tab,tab_name);
      if (verbose>0) {
	cout << "Read table " << tab_name << " lines: " 
	     << tab.get_nlines() << endl;
      }

      // Apply constraint to this table
      if (constraint.size()>0) {
	size_t nlines_old=tab.get_nlines();
	tab.delete_rows(constraint);
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
	  for(size_t ell=0;ell<grid_size;ell++) {

	    // Column name and value
	    string col=set_prefix+"_"+szttos(ell);
	    double val=tab.get(col,k);

	    // Add the value even if it's zero (blank)
	    values_ser[ell].push_back(val);

	    // But if it's less than or equal to zero, don't count
	    // it for min, max, or count. This is important
	    // because we don't want to count past the M-R curve
	    // or the end of the EOS.
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
	  weights.push_back(tab.get("mult",k));
	}
	line_counter++;

	// Next table line (k)
      }

      // Next chain (j)
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
  } else {
    ser_min*=xscale;
    ser_max*=xscale;
  }
  if (type==((string)"x")) {
    if (yset) {
      ser_min=user_ylow;
      ser_max=user_yhigh;
    }
  } else {
    if (xset) {
      ser_min=user_xlow;
      ser_max=user_xhigh;
    }
  }

  if (verbose>0) {
    cout << "Total lines: " << weights.size() << " " 
	 << values_ser[0].size() << " ser_min: " << ser_min << " ser_max: " 
	 << ser_max << endl;
  }

  // Create expval_scalar objects
  if (n_blocks<=0) n_blocks=20;
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

  // Vectors to store contour results. Each additional
  // contour level requires two vectors: a low and high
  // vector.
  ubmatrix cont_res(cont_levels.size()*2,grid_size);

  // Fill expval_scalar objects using histogram
  size_t block_size=weights.size()/20;
  for(size_t k=0;k<grid_size;k++) {

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
	  vector_region_parint(cont_wgt.size(),cont_grid,cont_wgt,
			       cont_levels[j],locs);
	  double lmin=vector_min_value<vector<double>,double>(locs);
	  double lmax=vector_max_value<vector<double>,double>(locs);
	  cout << lmin << " " << lmax << " ";
	  cont_res(2*j,k)=lmin;
	  cont_res(2*j+1,k)=lmax;
	} else {
	  // If there's not enough data, just report zero
	  cont_res(2*j,k)=0.0;
	  cont_res(2*j+1,k)=0.0;
	}
      } else {
	// If there's not enough data, just report zero
	cont_res(2*j,k)=0.0;
	cont_res(2*j+1,k)=0.0;
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
    for(size_t k=0;k<grid_size;k++) {
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
  cout << "Writing table to file " << sv[5] << endl;
  hdf_file hf;
  hf.open_or_create(sv[5]);
  hdf_output(hf,t3d,"hist_set");
  hdf_output(hf,index_grid,"index_grid");
  hf.setd_mat_copy("cont_res",cont_res);
      
  hf.close();

  return 0;
}

/** \brief Specify which contour levels to use
 */
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

/** \brief Combine several <tt>bamr</tt> output files
 */
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

  double e_low, e_high, m_low, m_high, nb_low, nb_high;
  ubvector low, high;
  size_t nparams, nsources, lgrid_size, file_row_start;
  vector<string> param_names, source_names;

  // data
  for(size_t i=0;i<nf;i++) {
    hdf_file hf;

    // Open file
    cout << "Opening file: " << files[i] << endl;
    hf.open(files[i]);

    // Get number of chains
    size_t n_chains;
    hf.get_szt_def("n_chains",1,n_chains);
    if (n_chains>1) {
      cout << n_chains << " separate chains." << endl;
    } else {
      cout << "1 chain." << endl;
    }

    size_t line_counter=0;

    file_row_start=full.get_nlines();

    // Read each chain
    for(size_t j=0;j<n_chains;j++) {

      // Read table
      std::string tab_name="markov_chain"+szttos(j);
      table_units<> tab;
      cout << "H1." << endl;
      hdf_input(hf,tab,tab_name);
      cout << "H2." << endl;
	
      // Set up columns in fill table over all files
      if (i==0 && j==0) {
	for(size_t ell=0;ell<tab.get_ncolumns();ell++) {
	  full.new_column(tab.get_column_name(ell));
	}
      }

      // Get misc parameters
      if (i==0 && j==0) {
	hf.getd("e_low",e_low);
	hf.getd("e_high",e_high);
	hf.getd("nb_low",nb_low);
	hf.getd("nb_high",nb_high);
	hf.getd("m_low",m_low);
	hf.getd("m_high",m_high);
	//hf.getd_vec_copy("low",low);
	//hf.getd_vec_copy("high",high);
	hf.get_szt("grid_size",lgrid_size);
	hf.get_szt("n_params",nparams);
	hf.get_szt("n_sources",nsources);
	//hf.gets_vec("param_names",param_names);
	hf.gets_vec("source_names",source_names);
      }
	  
      // Add data to table for this file
      
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

	// For the remaining chains we don't skip any, we just
	// start at the first row
	for(size_t k=0;k<tab.get_nlines();k+=thin_factor) {
	  vector<double> line;
	  for(size_t ell=0;ell<tab.get_ncolumns();ell++) {
	    line.push_back(tab.get(tab.get_column_name(ell),k));
	  }
	  full.line_of_data(line.size(),line);
	}
	
	cout << "Added table " << tab_name << " lines: " 
	     << tab.get_nlines() << endl;
      }

      // Go to next chain
    }
	
	
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
	
    cout << "full.nlines=" << full.get_nlines() << endl;

    // Go to next file
  }

  cout << "Writing table to file " << out_file << endl;
  hdf_file hf;
  hf.open_or_create(out_file);
  string type;
  if (hf.find_group_by_name("markov_chain0",type)==0) {
    O2SCL_ERR("Refusing to overwrite file with data in process 'combine'.",
	      exc_efailed);
  }
  hf.setd("e_low",e_low);
  hf.setd("e_high",e_high);
  hf.setd("nb_low",nb_low);
  hf.setd("nb_high",nb_high);
  hf.setd("m_low",m_low);
  hf.setd("m_high",m_high);
  hf.set_szt("grid_size",lgrid_size);
  hf.set_szt("n_chains",1);
  //hf.setd_vec_copy("low",low);
  //hf.setd_vec_copy("high",high);
  hf.set_szt("nparams",nparams);
  hf.set_szt("nsources",nsources);
  hf.sets_vec("param_names",param_names);
  hf.sets_vec("source_names",source_names);
  hdf_output(hf,full,"markov_chain0");
  hf.close();

  return 0;
}
//@}
    
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
  n_blocks=0;
}

void process::setup_cli() {

  // ---------------------------------------
  // Set options
  
  static const int nopt=8;
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
    {0,"hist-set","Create an ensemble of of 1-d histograms.",3,-1,
     "<direction> <low> <high> <set> <out_file> <file1> [file2...]",
     ((string)"Using the 'markov_chain' objects in a bamr output ")+
     "file, create an ensemble of 1-d histograms from a set "+
     "of columns in each chain. Typical uses are: "+
     "\'process -hist-set y m_low m_high R out.o2 x_0_out\', "+
     "\'process -hist-set x e_low e_high P out.o2 x_0_out\', and "+
     "\'process -hist-set x nb_low nb_high P out.o2 x_0_out\'. ",
     new comm_option_mfptr<process>(this,&process::hist_set),
     cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);

  // ---------------------------------------
  // Set parameters
    
  p_xscale.d=&xscale;
  p_xscale.help="Scale parameter for first variable";
  cl.par_list.insert(make_pair("xscale",&p_xscale));
    
  p_yscale.d=&yscale;
  p_yscale.help="Scale parameter for second variable";
  cl.par_list.insert(make_pair("yscale",&p_yscale));

  p_errors.b=&errors;
  p_errors.help="If true, output more error information";
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
  p_n_blocks.help="Number of blocks (default 20).";
  cl.par_list.insert(make_pair("n_blocks",&p_n_blocks));

  p_line_start.i=&line_start;
  p_line_start.help=((string)"Number of initial rows to skip over ")+
    "when reading MCMC tables.";
  cl.par_list.insert(make_pair("line_start",&p_line_start));

  p_verbose.i=&verbose;
  p_verbose.help="Verbosity (default 1).";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  p_constraint.str=&constraint;
  p_constraint.help="Constraint (default is an empty string).";
  cl.par_list.insert(make_pair("constraint",&p_constraint));

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

