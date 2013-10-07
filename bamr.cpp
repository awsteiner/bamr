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
#include <mpi.h>

#include <o2scl/vector.h>

#include "bamr.h"

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
  
bamr::bamr() {

  nparams=0;
  nsources=0;
  first_file_update=false;
  hist_size=100;
  model_type="";

  // Default parameter values

  first_point_file="";
  debug_load=false;
  step_fac=15.0;
  schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;
  // Default to 24 hours
  max_time=3.6e3*24;
  n_warm_up=500;
  in_file="default.in";
  // Minimum allowed maximum mass
  min_max_mass=2.0;
  // Minimum neutron star mass
  min_mass=0.8;
  debug_star=false;
  debug_eos=false;
  output_next=true;
  baryon_density=true;
  exit_mass=10.0;
  input_dist_thresh=0.0;
  use_crust=true;
  user_seed=0;
  best_detail=false;
  max_iters=0;
  chain_size=0;
  n_chains=0;
  norm_max=true;
  has_eos=true;

  in_m_min=0.8;
  in_m_max=3.0;
  in_r_min=5.0;
  in_r_max=18.0;
  
  // -----------------------------------------------------------
  // Grid limits
  
  nb_low=0.04;
  nb_high=1.24;
  
  e_low=0.3;
  e_high=10.0;

  // We set the mass minimum to 0.2 so that we get a full mass vs.
  // radius curve, even though the data doesn't go that far down
  //r_low=5.0;
  //r_high=18.0;
  m_low=0.2;
  m_high=3.0;
  
  modp=0;
  modp2=0;

  // -----------------------------------------------------------
  // Set up TOV solvers

  teos.verbose=0;

  def_ts.verbose=0;
  def_ts.set_units("1/fm^4","1/fm^4","1/fm^3");
  def_ts.set_eos(teos);
  def_ts.err_nonconv=false;
  ts=&def_ts;

  def_ts2.verbose=0;
  def_ts2.set_units("1/fm^4","1/fm^4","1/fm^3");
  def_ts2.set_eos(teos);
  def_ts2.err_nonconv=false;
  ts2=&def_ts2;
}

bamr::~bamr() {
  if (modp!=0) delete modp;
  if (modp2!=0) delete modp2;
}

void bamr::table_names_units(std::string &s, std::string &u) {
  s="N mult ";
  u+=". . ";
  s+="weight ";
  if (norm_max) {
    u+=". ";
  } else {
    u+=((string)"1/km^")+szttos(nsources)+"/Msun^"+
      szttos(nsources)+" ";
  }
  for(size_t i=0;i<nsources;i++) {
    s+=((string)"wgt_")+source_names[i]+" ";
    if (norm_max) {
      u+=". ";
    } else {
      u+="1/km/Msun ";
    }
  }
  for(size_t i=0;i<nparams;i++) {
    s+=((string)"param_")+modp->param_name(i)+" ";
    u+=modp->param_unit(i)+" ";
  }
  for(size_t i=0;i<nsources;i++) {
    s+=((string)"R_")+source_names[i]+" ";
    u+="km ";
  }
  for(size_t i=0;i<nsources;i++) {
    s+=((string)"M_")+source_names[i]+" ";
    u+="Msun ";
  }
  if (has_eos) {
    for(size_t i=0;i<hist_size;i++) {
      s+=((string)"P_")+szttos(i)+" ";
      u+="1/fm^4 ";
    }
  }
  for(size_t i=0;i<hist_size;i++) {
    s+=((string)"R_")+szttos(i)+" ";
    u+="km ";
    if (has_eos) {
      s+=((string)"PM_")+szttos(i)+" ";
      u+="1/fm^4 ";
    }
  }
  if (has_eos) {
    if (baryon_density) {
      for(size_t i=0;i<hist_size;i++) {
	s+=((string)"Pnb_")+szttos(i)+" ";
	u+="1/fm^4 ";
	s+=((string)"EoA_")+szttos(i)+" ";
	u+="MeV ";
      }
    }
    if (has_esym) {
      s+="S L ";
      u+="MeV MeV ";
    }
    s+="M_max P_max e_max ";
    u+="Msun 1/fm^4 1/fm^4 ";
    if (baryon_density) {
      s+="nb_max ";
      u+="1/fm^3 ";
    }
    for(size_t i=0;i<nsources;i++) {
      s+=((string)"ce_")+source_names[i]+" ";
      u+="1/fm^4 ";
    }
    if (baryon_density) {
      for(size_t i=0;i<nsources;i++) {
	s+=((string)"cnb_")+source_names[i]+" ";
	u+="1/fm^3 ";
      }
    }
  }

  return;
}

void bamr::init_grids_table(entry &low, entry &high) {
  
  scr_out << "In init()." << std::endl;

  if (low.np==0 || high.np==0) {
    O2SCL_ERR("No parameters in bamr::init().",gsl_einval);
  }

  // -----------------------------------------------------------
  // Make grids

  nb_grid=uniform_grid_end<double>(nb_low,nb_high,hist_size);
  e_grid=uniform_grid_end<double>(e_low,e_high,hist_size);
  m_grid=uniform_grid_end<double>(m_low,m_high,hist_size);

  e_hist.set_bin_edges(e_grid);
  m_hist.set_bin_edges(m_grid);
  nb_hist.set_bin_edges(nb_grid);

  // -----------------------------------------------------------
  // Init table

  std::string s, u;
  table_names_units(s,u);
  tc.line_of_names(s);

  {
    size_t ctr=0;
    std::string unit;
    std::istringstream is(u);
    while(is >> unit) {
      if (unit!=((string)".")) {
	tc.set_unit(tc.get_column_name(ctr),unit);
      }
      ctr++;
    } 
    if (ctr!=tc.get_ncolumns()) {
      O2SCL_ERR("Column/unit alignment in bamr::init_grids_table().",
		gsl_esanity);
    }
  }

  return;
}

void bamr::fill_line
(entry &e, o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_eos,
 o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_mvsr,
 double weight, bool new_meas, size_t n_meas, ubvector &wgts,
 std::vector<double> &line) {

  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0;

  if (has_eos) {

    // The highest baryon density in the EOS table
    nbmax2=tab_eos->max("nb");
    // The central energy density in the maximum mass configuration
    emax=tab_mvsr->max("ed");
    // The central pressure in the maximum mass configuration
    pmax=tab_mvsr->max("pr");
    // The maximum mass
    mmax=tab_mvsr->max("gm");

    if (baryon_density) {
      // The central baryon density in the maximum mass configuration
      nbmax=tab_mvsr->get("nb",tab_mvsr->lookup("gm",mmax));
    }

  } else {
    // Need to set mmax for no EOS models to figure out how 
    // high up we should go for the radius grid 
    mmax=3.0;
  }

  line.push_back(n_meas);
  line.push_back(1.0);
  line.push_back(weight);
  for(size_t i=0;i<nsources;i++) {
    line.push_back(wgts[i]);
  }
  for(size_t i=0;i<nparams;i++) {
    line.push_back(e.params[i]);
  }
  for(size_t i=0;i<nsources;i++) {
    line.push_back(e.rad[i]);
  }
  for(size_t i=0;i<nsources;i++) {
    line.push_back(e.mass[i]);
  }
  if (has_eos) {
    for(size_t i=0;i<hist_size;i++) {
      double eval=e_hist.get_rep_i(i);
      // Make sure the energy density from the grid
      // isn't beyond the table limit
      double emax2=tab_eos->max("ed");
      if (eval<emax2) {
	double pres_temp=tab_eos->interp("ed",eval,"pr");
	if (pres_temp<pmax) {
	  line.push_back(pres_temp);
	} else {
	  line.push_back(0.0);
	}
      } else {
	line.push_back(0.0);
      }
    }
  }
  for(size_t i=0;i<hist_size;i++) {
    double mval=m_hist.get_rep_i(i);
    if (mval<mmax) {
      line.push_back(tab_mvsr->interp("gm",mval,"r"));
      if (has_eos) {
	line.push_back(tab_mvsr->interp("gm",mval,"pr"));
      }
    } else {
      line.push_back(0.0);
      if (has_eos) {
	line.push_back(0.0);
      }
    }
  }
  if (has_eos) {
    if (baryon_density) {
      for(size_t i=0;i<hist_size;i++) {
	double nbval=nb_hist.get_rep_i(i);
	if (nbval<nbmax2) {
	  double pres_temp=tab_eos->interp("nb",nbval,"pr");
	  if (pres_temp<pmax) {
	    line.push_back(pres_temp);
	  } else {
	    line.push_back(0.0);
	  }
	  double eval2=tab_eos->interp("nb",nbval,"ed");
	  double eoa_val2=eval2/nbval-939.0/o2scl_const::hc_mev_fm;
	  line.push_back(eoa_val2);
	} else {
	  line.push_back(0.0);
	  line.push_back(0.0);
	}
      }
    }
    if (has_esym) {
      line.push_back(tab_eos->get_constant("S"));
      line.push_back(tab_eos->get_constant("L"));
    }
    line.push_back(mmax);
    line.push_back(pmax);
    line.push_back(emax);
    if (baryon_density) line.push_back(nbmax);
    for(size_t i=0;i<nsources;i++) {
      line.push_back(tab_mvsr->interp("gm",e.mass[i],"ed"));
    }
    if (baryon_density) {
      for(size_t i=0;i<nsources;i++) {
	line.push_back(tab_mvsr->interp("gm",e.mass[i],"nb"));
      }
    }
  }

  return;
}

void bamr::add_measurement
(entry &e, o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_eos,
 o2scl::o2_shared_ptr<o2scl::table_units<> >::type tab_mvsr,
 double weight, bool new_meas, size_t n_meas, ubvector &wgts) {
  
  // Test to see if we need to add a new line of data or
  // increment the weight on the previous line
  if (tc.get_nlines()==0 || new_meas==true) {
    // weight!=tc.get("weight",tc.get_nlines()-1)) {
    
    std::vector<double> line;
    fill_line(e,tab_eos,tab_mvsr,weight,new_meas,n_meas,wgts,line);
    
    // Done adding values, check size and add to table
    if (line.size()!=tc.get_ncolumns()) {
      scr_out << "line.size(): " << line.size() << endl;
      scr_out << "tc.get_ncolumns(): " << tc.get_ncolumns() << endl;
      for(size_t i=0;i<line.size() && i<tc.get_ncolumns();i++) {
	scr_out << line[i] << " " << tc.get_column_name(i) << endl;
      }
      O2SCL_ERR("Alignment problem.",gsl_efailed);
    }
    tc.line_of_data(line.size(),line);
    
  } else if (tc.get_nlines()>0) {
    tc.set("mult",tc.get_nlines()-1,
	   tc.get("mult",tc.get_nlines()-1)+1.0);
  }
  
  return;
}

void bamr::first_update(hdf_file &hf, model &modp) {

  vector<string> param_names;
  for(size_t i=0;i<nparams;i++) {
    param_names.push_back(modp.param_name(i));
  }
  hf.sets_vec("param_names",param_names);

  hf.sets_vec("source_names",source_names);
  hf.sets_vec("source_fnames",source_fnames);
  hf.sets_vec("slice_names",slice_names);

  hf.set_szt("hist_size",hist_size);
  hf.set_szt("nparams",nparams);
  hf.set_szt("nsources",nsources);
  hf.sets("model",model_type);
  hf.sets("in_file",in_file);
  hf.setd("max_time",max_time);
  hf.seti("user_seed",user_seed);
  hf.seti("n_warm_up",n_warm_up);
  hf.setd("min_mass",min_mass);
  hf.setd("exit_mass",exit_mass);
  hf.setd("min_max_mass",min_max_mass);
  hf.setd("step_fac",step_fac);
  hf.setd("input_dist_thresh",input_dist_thresh);
  hf.seti("use_crust",use_crust);
  hf.seti("baryon_density",baryon_density);
  hf.seti("max_iters",max_iters);
  hf.seti("debug_load",debug_load);
  hf.seti("debug_eos",debug_eos);
  hf.seti("debug_star",debug_star);
  hf.seti("best_detail",best_detail);
  hf.seti("output_next",output_next);
  hf.setd("nb_low",nb_low);
  hf.setd("nb_high",nb_high);
  hf.setd("e_low",e_low);
  hf.setd("e_high",e_high);
  hf.setd("m_low",m_low);
  hf.setd("m_high",m_high);
    
  std::vector<double> low_vec, high_vec;
  for(size_t i=0;i<nparams;i++) {
    low_vec.push_back(low.params[i]);
    high_vec.push_back(high.params[i]);
  }
  for(size_t i=0;i<nsources;i++) {
    low_vec.push_back(low.mass[i]);
    high_vec.push_back(high.mass[i]);
  }
  for(size_t i=0;i<nsources;i++) {
    low_vec.push_back(low.rad[i]);
    high_vec.push_back(high.rad[i]);
  }
  hf.setd_vec("low",low_vec);
  hf.setd_vec("high",high_vec);

  return;
}

void bamr::update_files(string fname_prefix, model &modp, entry &e_current) {

  hdf_file hf;

  // File for parameter histograms
  hf.open_or_create(fname_prefix+"_out");
    
  // First time, output some initial quantities
  if (first_file_update==false) {
    first_update(hf,modp);
    first_file_update=true;
  }

  hf.set_szt("mh_success",mh_success_cnt);
  hf.set_szt("mh_failure",mh_failure_cnt);

  // Store full Markov chain
  if (n_chains==0) n_chains++;
  hf.set_szt("n_chains",n_chains);
  string ch_name="markov_chain"+szttos(n_chains-1);
  hdf_output(hf,tc,ch_name);
  if (tc.get_nlines()>10000) {
    tc.clear_data();
    n_chains++;
    // Store the new empty chain in the HDF5 file just in case we
    // stop before the next call to update_files().
    string ch_name="markov_chain"+szttos(n_chains-1);
    hdf_output(hf,tc,ch_name);
  }

  hf.close();

  return;
}

void bamr::load_mc() {

  double tot, max;

  string name;
    
  if (nsources>0) {

    source_tables.resize(nsources);
    
    scr_out << "\nInput data files: " << endl;
    
    if (norm_max) {
      scr_out << "Normalizing maximum probability to 1." << endl;
    } else {
      scr_out << "Normalizing integral of distribution to 1." << endl;
    }

    scr_out << "File name total max P(10,1.4)" << endl;

#ifdef BAMR_MPI_LOAD

    bool mpi_load_debug=true;
    int buffer=0, tag=0;
    
    // Choose which file to read first for this rank
    int filestart=0;
    if (mpi_rank>mpi_nprocs-((int)nsources) && mpi_rank>0) {
      filestart=mpi_nprocs-mpi_rank;
    }
    if (mpi_load_debug) {
      scr_out << "Variable 'filestart' is " << filestart << " for rank "
	      << mpi_rank << "." << endl;
    }
    
    // Loop through all files
    for(int k=0;k<((int)nsources);k++) {
      
      // For k=0, we choose some ranks to begin reading, the others
      // have to wait. For k>=1, all ranks have to wait their turn.
      if (k>0 || (mpi_rank>0 && mpi_rank<=mpi_nprocs-((int)nsources))) {
	int prev=mpi_rank-1;
	if (prev<0) prev+=mpi_nprocs;
	if (mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " waiting for " 
		  << prev << "." << endl;
	}
	MPI_Recv(&buffer,1,MPI_INT,prev,tag,MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
      }
      
      // Determine which file to read next
      int file=filestart+k;
      if (file>=((int)nsources)) file-=nsources;

      if (mpi_load_debug) {
	scr_out << "Rank " << mpi_rank << " reading file " 
		<< file << "." << endl;
      }

      // 
      {
	hdf_file hf;
	hf.open(source_fnames[file]);
	hdf_input(hf,source_tables[file]);
	hf.close();
      }
      
      // Update input limits
      if (file==0) {
	in_r_min=source_tables[file].get_grid_x(0);
	in_r_max=source_tables[file].get_grid_x(source_tables[file].get_nx()-1);
	in_m_min=source_tables[file].get_grid_y(0);
	in_m_max=source_tables[file].get_grid_y(source_tables[file].get_ny()-1);
      } else {
	if (in_r_min>source_tables[file].get_grid_x(0)) {
	  in_r_min=source_tables[file].get_grid_x(0);
	}
	if (in_r_max<source_tables[file].get_grid_x
	    (source_tables[file].get_nx()-1)) {
	  in_r_max=source_tables[file].get_grid_x
	    (source_tables[file].get_nx()-1);
	}
	if (in_m_min>source_tables[file].get_grid_y(0)) {
	  in_m_min=source_tables[file].get_grid_y(0);
	}
	if (in_m_max<source_tables[file].get_grid_y
	    (source_tables[file].get_ny()-1)) {
	  in_m_max=source_tables[file].get_grid_y
	    (source_tables[file].get_ny()-1);
	}
      }

      // Renormalize
      tot=0.0;
      max=0.0;
      for(size_t i=0;i<source_tables[file].get_nx();i++) {
	for(size_t j=0;j<source_tables[file].get_ny();j++) {
	  tot+=source_tables[file].get(i,j,slice_names[file]);
	  if (source_tables[file].get(i,j,slice_names[file])>max) {
	    max=source_tables[file].get(i,j,slice_names[file]);
	  }
	}
      }
      for(size_t i=0;i<source_tables[file].get_nx();i++) {
	for(size_t j=0;j<source_tables[file].get_ny();j++) {
	  if (norm_max) {
	    source_tables[file].set
	      (i,j,slice_names[file],
	       source_tables[file].get(i,j,slice_names[file])/max);
	  } else {
	    source_tables[file].set
	      (i,j,slice_names[file],
	       source_tables[file].get(i,j,slice_names[file])/tot);
	  }
	}
      }

      if (debug_load) {
	cout << source_fnames[file] << endl;
	for(size_t i=0;i<source_tables[file].get_nx();i++) {
	  cout << i << " " << source_tables[file].get_grid_x(i) << endl;
	}
	for(size_t j=0;j<source_tables[file].get_ny();j++) {
	  cout << j << " " << source_tables[file].get_grid_y(j) << endl;
	}
	for(size_t i=0;i<source_tables[file].get_nx();i++) {
	  for(size_t j=0;j<source_tables[file].get_ny();j++) {
	    cout << source_tables[file].get(i,j,slice_names[file]) << " ";
	  }
	  cout << endl;
	}
      }

      scr_out.setf(ios::left);
      scr_out.width(25);
      scr_out << source_fnames[file] << " ";
      scr_out.width(6);
      scr_out << source_names[file] << " " << tot << " " << max << " ";
      scr_out.unsetf(ios::left);
      scr_out << source_tables[file].interp(10.0,1.4,slice_names[file]) << endl;
      
      // Send a message, unless the rank is the last one to read a
      // file.
      if (k<((int)nsources)-1 || mpi_rank<mpi_nprocs-((int)nsources)) {
	int next=mpi_rank+1;
	if (next>=mpi_nprocs) next-=mpi_nprocs;
	if (mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " sending to " 
		  << next << "." << endl;
	}
	MPI_Send(&buffer,1,MPI_INT,next,tag,MPI_COMM_WORLD);
      }
      
    }
    
#else

    for(size_t k=0;k<nsources;k++) {

      {
	hdf_file hf;
	hf.open(source_fnames[k]);
	hdf_input(hf,source_tables[k]);
	hf.close();
      }
      
      // Update input limits
      if (k==0) {
	in_r_min=source_tables[k].get_grid_x(0);
	in_r_max=source_tables[k].get_grid_x(source_tables[k].get_nx()-1);
	in_m_min=source_tables[k].get_grid_y(0);
	in_m_max=source_tables[k].get_grid_y(source_tables[k].get_ny()-1);
      } else {
	if (in_r_min>source_tables[k].get_grid_x(0)) {
	  in_r_min=source_tables[k].get_grid_x(0);
	}
	if (in_r_max<source_tables[k].get_grid_x
	    (source_tables[k].get_nx()-1)) {
	  in_r_max=source_tables[k].get_grid_x(source_tables[k].get_nx()-1);
	}
	if (in_m_min>source_tables[k].get_grid_y(0)) {
	  in_m_min=source_tables[k].get_grid_y(0);
	}
	if (in_m_max<source_tables[k].get_grid_y
	    (source_tables[k].get_ny()-1)) {
	  in_m_max=source_tables[k].get_grid_y(source_tables[k].get_ny()-1);
	}
      }

      // Renormalize
      tot=0.0;
      max=0.0;
      for(size_t i=0;i<source_tables[k].get_nx();i++) {
	for(size_t j=0;j<source_tables[k].get_ny();j++) {
	  tot+=source_tables[k].get(i,j,slice_names[k]);
	  if (source_tables[k].get(i,j,slice_names[k])>max) {
	    max=source_tables[k].get(i,j,slice_names[k]);
	  }
	}
      }
      for(size_t i=0;i<source_tables[k].get_nx();i++) {
	for(size_t j=0;j<source_tables[k].get_ny();j++) {
	  if (norm_max) {
	    source_tables[k].set(i,j,slice_names[k],
				 source_tables[k].get(i,j,slice_names[k])/max);
	  } else {
	    source_tables[k].set(i,j,slice_names[k],
				 source_tables[k].get(i,j,slice_names[k])/tot);
	  }
	}
      }

      if (debug_load) {
	cout << source_fnames[k] << endl;
	for(size_t i=0;i<source_tables[k].get_nx();i++) {
	  cout << i << " " << source_tables[k].get_grid_x(i) << endl;
	}
	for(size_t j=0;j<source_tables[k].get_ny();j++) {
	  cout << j << " " << source_tables[k].get_grid_y(j) << endl;
	}
	for(size_t i=0;i<source_tables[k].get_nx();i++) {
	  for(size_t j=0;j<source_tables[k].get_ny();j++) {
	    cout << source_tables[k].get(i,j,slice_names[k]) << " ";
	  }
	  cout << endl;
	}
      }

      scr_out.setf(ios::left);
      scr_out.width(25);
      scr_out << source_fnames[k] << " ";
      scr_out.width(6);
      scr_out << source_names[k] << " " << tot << " " << max << " ";
      scr_out.unsetf(ios::left);
      scr_out << source_tables[k].interp(10.0,1.4,slice_names[k]) << endl;

    }

#endif

    scr_out << endl;
  }

  if (in_m_min<min_mass) in_m_min=min_mass;

  // Merciless hack
  //in_m_max=2.5;
  
  scr_out << "M limits: (" 
	  << in_m_min << "," << in_m_max << ")" << endl;
  scr_out << "R limits: ("
	  << in_r_min << "," << in_r_max << ")" << endl;
  
  return;
}
  
void bamr::prepare_eos(entry &e, model &modref, tov_solve *tsr, 
		       bool &success) {
  return;
}

void bamr::compute_star(entry &e, model &modref, tov_solve *tsr, 
			bool &success) {
  
  success=true;

  // Compute the EOS first
  bool test_eos_fail=false;
  if (has_eos) {
    modref.compute_eos(e,test_eos_fail,scr_out);
  }

  // Ensure we're using linear interpolation
  o2_shared_ptr<table_units<> >::type tab_eos=modref.cns.get_eos_results();

  if (has_eos) {
    tab_eos->set_interp_type(itp_linear);

    // Check that pressure is increasing. If we're using a crust EOS,
    // choose 0.6 as an arbitrary low-density cutoff which corresponds
    // to about n_B=0.12 fm^{-3}, and check the low-density part below
    // instead.
  
    for(size_t i=0;(tab_eos->get_nlines()>0 && i<tab_eos->get_nlines()-1);i++) {
      if ((!use_crust || tab_eos->get("ed",i)>0.6) && 
	  tab_eos->get("pr",i+1)<tab_eos->get("pr",i)) {
	scr_out << "Rejected: Pressure decreasing." << endl;
	scr_out << "ed=" << tab_eos->get("ed",i) 
		<< " pr=" << tab_eos->get("pr",i) << endl;
	scr_out << "ed=" << tab_eos->get("ed",i+1) 
		<< " pr=" << tab_eos->get("pr",i+1) << endl;
	success=false;
	return;
      }
    }
  }

  // If it failed, stop
  if (test_eos_fail==true) {
    success=false;
    return;
  }

  // If requested, compute the baryon density automatically
  if (has_eos && baryon_density) {

    // Obtain the baryon density calibration point from the model

    double n1, e1;
    modref.baryon_density_point(n1,e1);
    
    if (n1<=0.0 && e1<=0.0) {
      O2SCL_ERR2("Computing the baryon density requires one ",
		 "calibration point in bamr::compute_star().",gsl_einval);
    }

    // Compute inverse of gibbs energy density, 'igb'
    tab_eos->new_column("igb");
    for(size_t i=0;i<tab_eos->get_nlines();i++) {

      if (tab_eos->get("ed",i)+tab_eos->get("pr",i)<=0.0 ||
	  !gsl_finite(tab_eos->get("ed",i)) ||
	  !gsl_finite(tab_eos->get("pr",i))) {
	scr_out << "Inverse Gibbs not finite." << endl;
	scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	scr_out << "ed pr" << endl;
	for(size_t i=0;i<tab_eos->get_nlines();i++) {
	  scr_out << i << " "
		  << tab_eos->get("ed",i) << " "
		  << tab_eos->get("pr",i) << endl;
	}
	O2SCL_ERR("Inverse Gibbs not finite.",gsl_efailed);
      }
      tab_eos->set("igb",i,1.0/(tab_eos->get("ed",i)+tab_eos->get("pr",i)));
    }

    // Compute integral of 'igb' relative to ed='e1', called 'iigb'
    tab_eos->new_column("iigb");
    for(size_t i=0;i<tab_eos->get_nlines();i++) {
      if (e1<=tab_eos->get("ed",i)) {
	double val=tab_eos->integ("ed",e1,tab_eos->get("ed",i),"igb");
	if (!gsl_finite(val)) {
	  scr_out << "Baryon integral not finite." << endl;
	  scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	  scr_out << "ed pr" << endl;
	  for(size_t i=0;i<tab_eos->get_nlines();i++) {
	    scr_out << i << " "
		    << tab_eos->get("ed",i) << " "
		    << tab_eos->get("pr",i) << endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",gsl_efailed);
	}
	tab_eos->set("iigb",i,val);
      } else {
	double val=-tab_eos->integ("ed",tab_eos->get("ed",i),e1,"igb");
	if (!gsl_finite(val)) {
	  scr_out << "Baryon integral not finite (2)." << endl;
	  scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	  scr_out << "ed pr" << endl;
	  for(size_t i=0;i<tab_eos->get_nlines();i++) {
	    scr_out << i << " "
		    << tab_eos->get("ed",i) << " "
		    << tab_eos->get("pr",i) << endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",gsl_efailed);
	}
	tab_eos->set("iigb",i,val);
      }
    }

    // Compute normalization constant
    double Anb=n1/exp(tab_eos->interp("ed",e1,"iigb"));
    if (!gsl_finite(Anb) || Anb<0.0) {
      scr_out << "Baryon density normalization problem." << endl;
      success=false;
      return;
    }

    // Now compute baryon density

    tab_eos->new_column("nb");

    for(size_t i=0;i<tab_eos->get_nlines();i++) {      

      // If the density is too low, then just use zero baryon density.
      // This pressure (10^{-5} fm^{-4}) corresponds to a baryon 
      // density of about 3e-3 fm^{-3}.

      if (use_crust==false && tab_eos->get("pr",i)<1.0e-5) {

	tab_eos->set("nb",i,0.0);

      } else {
	
	double nbt=Anb*exp(tab_eos->get("iigb",i));
	if (!gsl_finite(nbt)) {
	  scr_out << "Baryon density normalization problem." << endl;
	  success=false;
	  return;
	  /*
	    cout << "Baryon density not finite." << endl;
	    cout << "i,ed,pr,igb,iigb: " << i << " " << Anb << " "
	    << tab_eos->get("ed",i) << " "
	    << tab_eos->get("pr",i) << " " << tab_eos->get("igb",i) << " "
	    << tab_eos->get("iigb",i) << endl;
	    exit(-1);
	  */
	} 
	tab_eos->set("nb",i,nbt);
      }
    }

    // End of loop 'if (has_eos && baryon_density) {' 
  }

  o2_shared_ptr<table_units<> >::type tab_mvsr;

  if (has_eos) {

    // Perform any additional necessary EOS preparations
    prepare_eos(e,modref,tsr,success);
    if (success==false) {
      return;
    }

    // Read the EOS into the tov_eos object. We don't specify a baryon
    // number density column here, even though it might be provided,
    // just to make the TOV solver a bit faster.
    tab_eos->set_unit("ed","1/fm^4");
    tab_eos->set_unit("pr","1/fm^4");
    teos.read_table(*tab_eos,"ed","pr");
    
    if (use_crust) {
    
      if (false) {
      
	scr_out << "Checking crust-core: " << endl;
      
	double plo, pt, phi;
	teos.get_transition(plo,pt,phi);
      
	double ed_last=0.0;
	for(double pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	  scr_out << pr << " ";
	
	  if (pr<plo) {
	    double edlo,nblo;
	    teos.get_eden_low(pr,edlo,nblo);
	    scr_out << edlo << " ";
	  } else {
	    scr_out << 0.0 << " ";
	  }
	
	  double ed, nb;
	  int ip=teos.get_eden_full(pr,ed,nb);
	  scr_out << ed << " ";
	  if (ip==tov_interp_eos::itrans) {
	    scr_out << "0 ";
	  } else if (ip==tov_interp_eos::icrust) {
	    scr_out << "- ";
	  } else {
	    scr_out << "+ ";
	  }
	
	  if (pr>phi) {
	    double edhi,nbhi;
	    teos.get_eden_high(pr,edhi,nbhi);
	    scr_out << edhi << " ";
	  } else {
	    scr_out << 0.0 << " ";
	  }
	
	  if (baryon_density) {
	    scr_out << tab_eos->interp("pr",pr,"nb");
	  }
	  scr_out << endl;
	}
	scr_out << endl;
      }

      double ed_last=0.0;
      // This range corresponds to between about n_B=0.01 and 0.17
      // fm^{-3}
      for(double pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	double ed, nb;
	teos.get_eden_full(pr,ed,nb);
	if (ed_last>1.0e-20 && ed<ed_last) {
	  scr_out << "Stability problem near crust-core transition." << endl;
	  if (has_esym) {
	    scr_out << "S=" << tab_eos->get_constant("S")*hc_mev_fm 
		    << " L=" << tab_eos->get_constant("L")*hc_mev_fm << endl;
	  }
	  scr_out << "ed_last,pr,ed: " << ed_last << " " << pr 
		  << " " << ed << endl;
	  scr_out << endl;
	  scr_out << "pr ed phase edlo edhi" << endl;
	  double edlo, edhi;
	  for(pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	    int ip=teos.get_eden_full(pr,ed,nb);
	    teos.get_eden_low(pr,edlo,nb);
	    teos.get_eden_high(pr,edhi,nb);
	    scr_out << pr << " " << ed;
	    scr_out << " " << ip << " " << edlo << " " << edhi << endl;
	  }
	  scr_out << endl;
	  success=false;
	  return;
	}
	ed_last=ed;
      }

      // End of 'if (use_crust)'
    }

    // If necessary, output debug information (We want to make sure this
    // is after tov_eos::read_table() so that we can debug the
    // core-crust transition)
    if (debug_eos) {
      hdf_file hfde;

      hfde.open_or_create("debug_eos.o2");

      // Output the core EOS
      hdf_output(hfde,*tab_eos,"eos");

      // Output core and crust EOS as reported by the tov_interp_eos object
      table_units<> full_eos;
      full_eos.line_of_names("ed pr");
      full_eos.set_unit("ed","1/fm^4");
      full_eos.set_unit("pr","1/fm^4");
      for(double pr=1.0e-20;pr<10.0;pr*=1.05) {
	double ed, nb;
	teos.get_eden_full(pr,ed,nb);
	double line[2]={ed,pr};
	full_eos.line_of_data(2,line);
	if (pr<1.0e-4) pr*=1.2;
      }
      hdf_output(hfde,full_eos,"full_eos");

      hfde.close();
      if (!debug_star) {
	scr_out << "Automatically exiting since 'debug_eos' is true." << endl;
	exit(-1);
      }
    }

    // Solve for M vs. R curve
    int ret=tsr->mvsr();
    if (ret!=0) {
      scr_out << "M vs. R failed." << endl;
      success=false;
      return;
    }
    tab_mvsr=tsr->get_results();
    tab_mvsr->set_interp_type(itp_linear);
  
    // If the EOS is sufficiently stiff, the TOV solver will output
    // gibberish, i.e. masses and radii equal to zero, especially at the
    // higher pressures. This rules these EOSs out, as they would likely
    // be acausal anyway. If this happens frequently, it might signal
    // a problem. 
    size_t ir=tab_mvsr->get_nlines()-1;
    if ((*tab_mvsr)["gm"][ir]<1.0e-10 ||
	(*tab_mvsr)["gm"][ir-1]<1.0e-10) {
      scr_out << "TOV failure fix." << endl;
      success=false;
      return;
    }

    // Check that maximum mass is large enough,
    double mmax=tab_mvsr->max("gm");
    if (mmax<min_max_mass) {
      scr_out << "Maximum mass too small: " << mmax << " < "
	      << min_max_mass << "." << endl;
      success=false;
      return;
    }

    // Check the radius of the maximum mass star
    size_t ix_max=tab_mvsr->lookup("gm",mmax);
    if (tab_mvsr->get("r",ix_max)>1.0e4) {
      scr_out << "TOV convergence problem: " << endl;
      success=false;
      return;
    }

    // Check that all stars have masses below the maximum mass
    // of the current M vs R curve
    bool mass_fail=false;
    for(size_t i=0;i<nsources;i++) {
      if (e.mass[i]>mmax) mass_fail=true;
    }
    if (true) {
      if (mass_fail==true) {
	bool bad_step;
	scr_out << "Adjusting masses for M_{max} = " << mmax << endl;
	scr_out << "Old entry: " << e << endl;
	select_mass(e,e,mmax,bad_step);
	mass_fail=false;
	for(size_t i=0;i<nsources;i++) {
	  if (e.mass[i]>mmax) mass_fail=true;
	}
	if (mass_fail || bad_step) {
	  scr_out << "Failure with new mass update." << endl;
	  O2SCL_ERR("Failure in new mass update.",gsl_efailed);
	}
	scr_out << "New entry: " << e << endl;
      }
    } else {
      if (mass_fail==true) {
	scr_out << "Rejected: Mass of object larger than maximum. " << endl;
	success=false;
	return;
      }
    }

    // Remove table entries with pressures above the maximum pressure
    double prmax=tab_mvsr->get("pr",
			       tab_mvsr->lookup("gm",tab_mvsr->max("gm")));
    tab_mvsr->delete_rows(((string)"pr>")+dtos(prmax));
  
    // Make sure that the M vs. R curve generated enough data. This
    // is not typically an issue.
    if (tab_mvsr->get_nlines()<10) {
      scr_out << "M vs. R failed to generate lines." << endl;
      success=false;
      return;
    }

    // Compute speed of sound squared
    tab_mvsr->deriv("ed","pr","dpde");

  } else {

    // If there's no EOS, then the model object gives the M-R curve
    tab_mvsr=tsr->get_results();
    modref.compute_mr(e,scr_out,tab_mvsr,success);
    if (success==false) {
      return;
    }

  }
  
  // Output M vs. R curve
  if (debug_star) {
    hdf_file hfds;
    hfds.open_or_create("debug_star.o2");
    hdf_output(hfds,*tab_mvsr,"mvsr");
    hfds.close();
    scr_out << "Automatically exiting since 'debug_star' is true." << endl;
    exit(-1);
  }

  // Compute the radius for each source
  for(size_t i=0;i<nsources;i++) {
    e.rad[i]=tab_mvsr->interp("gm",e.mass[i],"r");
  }

  // Check causality
  if (has_eos) {
    
    for(size_t i=0;i<tab_mvsr->get_nlines();i++) {
      if ((*tab_mvsr)["dpde"][i]>1.0) {
	scr_out.precision(4);
	scr_out << "Rejected: Acausal."<< endl;
	scr_out << "ed_max="
		<< tab_mvsr->max("ed") << " ed_bad="
		<< (*tab_mvsr)["ed"][i] << " pr_max=" 
		<< tab_mvsr->max("pr") << " pr_bad=" 
		<< (*tab_mvsr)["pr"][i] << endl;
	scr_out.precision(6);
	success=false;
	return;
      }
    }

  } else {

    for(size_t i=0;i<nsources;i++) {
      if (e.rad[i]<2.94*schwarz_km/2.0*e.mass[i]) {
	scr_out << "Source " << source_names[i] << " acausal." << endl;
	success=false;
	return;
      }
    }

  }

  return;
}

double bamr::compute_weight(entry &e, model &modref, tov_solve *tsr, 
			    bool &success, ubvector &wgts, bool warm_up) {
			     
  // Compute the M vs R curve and return if it failed
  bool compute_star_success;
  compute_star(e,modref,tsr,compute_star_success);
  if (compute_star_success==false) {
    success=false;
    return 0.0;
  }
    
  bool mr_fail=false;
  for(size_t i=0;i<nsources;i++) {
    if (e.mass[i]<in_m_min || e.mass[i]>in_m_max ||
	e.rad[i]<in_r_min || e.rad[i]>in_r_max) {
      mr_fail=true;
    }
  }

  if (mr_fail==true) {
    scr_out << "Rejected: Mass or radius outside range." << endl;
    if (nsources>0) {
      scr_out.precision(2);
      scr_out.setf(ios::showpos);
      for(size_t i=0;i<nsources;i++) {
	scr_out << e.mass[i] << " ";
      }
      scr_out << endl;
      for(size_t i=0;i<nsources;i++) {
	scr_out << e.rad[i] << " ";
      }
      scr_out << endl;
      scr_out.precision(6);
      scr_out.unsetf(ios::showpos);
    }
    success=false;
    return 0.0;
  }

  success=true;
  double ret=1.0;

  o2_shared_ptr<table_units<> >::type tab_mvsr=tsr->get_results();
  tab_mvsr->set_interp_type(itp_linear);
  double m_max_current=tab_mvsr->max("gm");

  // -----------------------------------------------
  // Compute the weights for each source

  if (debug_star) scr_out << "Name M R Weight" << endl;

  for(size_t i=0;i<nsources;i++) {
	
    // Double check that current M and R is in the range of
    // the provided input data
    if (e.rad[i]<source_tables[i].get_x_data()[0] ||
	e.rad[i]>source_tables[i].get_x_data()[source_tables[i].get_nx()-1] ||
	e.mass[i]<source_tables[i].get_y_data()[0] ||
	e.mass[i]>source_tables[i].get_y_data()[source_tables[i].get_ny()-1]) {
      wgts[i]=0.0;
    } else {
      // If it is, compute the weight
      wgts[i]=source_tables[i].interp(e.rad[i],e.mass[i],slice_names[i]);
    }

    // If the weight is lower than the threshold, set it equal
    // to the threshold
    if (wgts[i]<input_dist_thresh) wgts[i]=input_dist_thresh;

    // Include the weight for this source 
    ret*=wgts[i];

    if (debug_star) {
      scr_out << source_names[i] << " " << e.mass[i] << " " 
	      << e.rad[i] << " " << wgts[i] << endl;
    }
	
    // Go to the next source
  }

  if (debug_star) scr_out << endl;
      
  // -----------------------------------------------
  // Exit if the current maximum mass is too large

  if (m_max_current>exit_mass) {
    scr_out.setf(ios::scientific);
    scr_out << "Exiting because maximum mass larger than 'exit_mass'." 
	    << endl;
    scr_out << "e,ret: " << e << " " << ret << endl;
    exit(-1);
  }

  return ret;
}

void bamr::read_input(entry &e_current) {
    
  bool debug=true;

  string stmp;
  
  // Open file
  if (debug) {
    scr_out << "Opening input file named: " 
	    << in_file.c_str() << endl;
  }
  ifstream fin(in_file.c_str());
    
  // Read number of sources, and allocate e_current 

  fin >> nsources;
  e_current.allocate(nparams,nsources);

  // Read parameter names and initial guesses

  if (debug) {
    scr_out << "Looking for " << nparams 
	    << " parameter names and initial guesses." << endl;
  }
  for(size_t i=0;i<e_current.np;i++) {
    fin >> stmp;
    fin >> e_current.params[i];
    if (debug) {
      scr_out << i << " " << e_current.params[i] << endl;
    }
  }

  // Read information for each source

  if (nsources>0) {

    fin >> stmp;
    fin >> stmp;
    fin >> stmp;
    fin >> stmp;
    fin >> stmp;
      
    source_names.resize(nsources);
    source_fnames.resize(nsources);
    slice_names.resize(nsources);
      
    for(size_t i=0;i<nsources;i++) {
      fin >> source_names[i];
      fin >> source_fnames[i];
      fin >> slice_names[i];
      fin >> e_current.mass[i];
      if (debug) {
	scr_out << i << " ";
	scr_out.width(6);
	scr_out << source_names[i] << " ";
	scr_out.width(30);
	scr_out.setf(ios::left);
	scr_out << source_fnames[i] << " ";
	scr_out.unsetf(ios::left);
	scr_out.width(6);
	scr_out << slice_names[i] << " " << e_current.mass[i] << endl;
      }
    }

  }

  fin.close();

  return;
}

int bamr::set_model(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()<2) {
    cout << "Model name not given." << endl;
    return gsl_efailed;
  }
  model_type=sv[1];
  if (modp!=0) {
    delete modp;
    modp=0;
  }
  if (modp2!=0) {
    delete modp2;
    modp2=0;
  }
  if (sv[1]==((string)"twop")) {
    modp=new two_polytropes;
    modp2=new two_polytropes;
    nparams=8;
    has_esym=true;
    has_eos=true;
    // Cannot output to scr_out since scr_out isn't set until mcmc()
    // is called. 
    // scr_out << "Selected two polytropes." << endl;
  } else if (sv[1]==((string)"altp")) {
    modp=new alt_polytropes;
    modp2=new alt_polytropes;
    nparams=8;
    has_esym=true;
    has_eos=true;
    //scr_out << "Selected improved polytropes." << endl;
  } else if (sv[1]==((string)"fixp")) {
    modp=new fixed_pressure;
    modp2=new fixed_pressure;
    nparams=8;
    has_esym=true;
    has_eos=true;
    //scr_out << "Selected fixed pressure." << endl;
  } else if (sv[1]==((string)"qstar")) {
    modp=new quark_star;
    modp2=new quark_star;
    nparams=4;
    has_esym=false;
    has_eos=true;
    //scr_out << "Selected strange quark star model." << endl;
  } else if (sv[1]==((string)"genq")) {
    modp=new generic_quarks;
    modp2=new generic_quarks;
    nparams=9;
    has_esym=true;
    has_eos=true;
    //scr_out << "Selected generic quarks." << endl;
  } else {
    cout << "Model unknown." << endl;
    nparams=0;
    model_type="";
    has_esym=true;
    has_eos=true;
    return gsl_efailed;
  }

  return 0;
}

void bamr::output_best(string fname_prefix, entry &e_best, double w_best,
		       o2_shared_ptr<table_units<> >::type tab_eos,
		       o2_shared_ptr<table_units<> >::type tab_mvsr,
		       ubvector &wgts) {
  
  scr_out << "Best: " << e_best << " " << w_best << endl;

  // "Best" point in parameter space and its associated weight
  string fname_best=fname_prefix+"_best";
  ofstream fout(fname_best.c_str(),ios::app);
  fout.setf(ios::scientific);
  fout.precision(12);
  fout << "Best: " << e_best << " Weights: ";
  for(size_t i=0;i<nsources;i++) {
    fout << wgts[i] << " ";
  }
  fout << w_best << endl;
  fout.close();
  
  string fname_best_out=fname_prefix+"_out";
  hdf_file hf;
  hf.open_or_create(fname_best_out);
  std::vector<double> best_point;
  for(size_t i=0;i<nparams;i++) {
    best_point.push_back(e_best.params[i]);
  }
  for(size_t i=0;i<nsources;i++) {
    best_point.push_back(e_best.mass[i]);
  }
  best_point.push_back(w_best);
  hf.setd_vec("best_point",best_point);
  hf.close();

  if (best_detail) {

    // "Best" EOS
    string fname_best_eos=fname_prefix+"_best_eos";
    hdf_file hf;
    hf.open_or_create(fname_best_eos);
    hdf_output(hf,*tab_eos,"best_eos");
    hf.close();
    
    // "Best" M vs. R curve
    string fname_best_mvsr=fname_prefix+"_best_mvsr";
    hf.open_or_create(fname_best_mvsr);
    hdf_output(hf,*tab_mvsr,"best_mvsr");
    hf.close();
  }

  return;
}

int bamr::mcmc_init() {
  return 0;
}

bool bamr::make_step(double w_current, double w_next, bool debug,
		     bool warm_up, int iteration) {
  
  double r=gr.random();

  bool accept=false;

  // Metropolis algorithm
  if (r<w_next/w_current) accept=true;

  if (debug) {
    cout << "Metropolis: " << r << " " << w_next/w_current << " " 
	 << accept << endl;
  }
  
  return accept;
}

void bamr::select_mass(entry &e_current, entry &e_next, double mmax,
		       bool &bad_step) {

  // Step for masses
  for(size_t k=0;k<nsources;k++) {
    size_t step_count=0;
    do {
      bad_step=false;
      if (e_current.mass[k]<mmax) {
	e_next.mass[k]=e_current.mass[k]+(gr.random()*2.0-1.0)*
	  (high.mass[k]-low.mass[k])/step_fac;
      } else {
	e_next.mass[k]=low.mass[k]+gr.random()*(mmax-low.mass[k]);
      }
      if (e_next.mass[k]<low.mass[k] || e_next.mass[k]>high.mass[k] || 
	  e_next.mass[k]>mmax) {
	bad_step=true;
      }
      if (bad_step) {
	scr_out << "Bad step: " << k << " " << low.mass[k] << " "
		<< e_current.mass[k] << " " << e_next.mass[k] << " "
		<< mmax << " " << high.mass[k] << endl;
      }
    } while (bad_step && step_count<100);
    if (step_count==100) {
      scr_out << "Too many steps in parameter space failed." << endl;
      bad_step=true;
      return;
    }
  }

  return;
}

int bamr::mcmc(std::vector<std::string> &sv, bool itive_com) {

  int iret=mcmc_init();
  if (iret!=0) return iret;

  // User-specified filename prefix
  if (sv.size()<2) {
    cout << "No filename prefix given in bamr::mcmc()." << endl;
    return gsl_efailed;
  }
  string fname_prefix=sv[1];
    
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_nprocs);
  mpi_start_time=MPI_Wtime();

  // Update filename with processor rank
  fname_prefix+=((string)"_")+itos(mpi_rank);
    
  // Open main output file
  scr_out.open((fname_prefix+"_scr").c_str());
  scr_out.setf(ios::scientific);
  
  // Check model
  if (model_type.length()==0 || modp==0 || modp2==0) {
    scr_out << "Model not set." << endl;
    return gsl_efailed;
  }

  // Low-density crust
  if (use_crust) {
    teos.default_low_dens_eos();

    // Get the transition density from the crust
    double plo, pt, phi;
    teos.get_transition(plo,pt,phi);
    // We set the transition density a bit lower (because by default
    // it's the largest density in the crust EOS) and then add a 
    // small width
    teos.transition_mode=tov_interp_eos::smooth_trans;
    teos.set_transition(pt/1.2,1.2);

  } else {
    teos.no_low_dens_eos();
  }

  bool debug=false;
  if (debug) cout.setf(ios::scientific);
  if (debug) cout.setf(ios::showpos);
    
  // Set RNG seed
  unsigned long int seed=time(0);
  if (user_seed!=0) seed=user_seed;
  seed*=(mpi_rank+1);
  gr.set_seed(seed);
  scr_out << "Using seed " << seed 
	  << " for processor " << mpi_rank+1 << "/" 
	  << mpi_nprocs << "." << endl;
  scr_out.precision(12);
  scr_out << " Start time: " << mpi_start_time << endl;
  scr_out.precision(6);

  // Read input file
  entry e_current;
  read_input(e_current);

  // Set lower and upper bounds for parameters
  low.allocate(nparams,nsources);
  high.allocate(nparams,nsources);
  modp->low_limits(low);
  modp->high_limits(high);

  n_chains=0;

  // Vector containing weights
  ubvector wgts(nsources);

  // Entry objects (Must be after read_input() since nsources is set
  // in that function.)
  entry e_next(nparams,nsources);
  entry e_best(nparams,nsources);

  // Load data
  load_mc();

  // If necessary, read initial guess
  if (first_point_file.length()>0) {
    std::vector<double> best_point;
    hdf_file hf;
    hf.open(first_point_file);
    hf.getd_vec("best_point",best_point);
    hf.close();
    scr_out << "Reading best point from file '" << first_point_file
	    << "'." << endl;
    for(size_t i=0;i<nparams;i++) {
      e_current.params[i]=best_point[i];
      scr_out << e_current.params[i] << endl;
    }
    scr_out << endl;
    for(size_t i=nparams;i<nsources+nparams;i++) {
      e_current.mass[i-nparams]=best_point[i];
      scr_out << e_current.mass[i-nparams] << endl;
    }
    scr_out << endl;
  }

  // Set lower and upper bounds for masses and radii from load_mc()
  for(size_t i=0;i<nsources;i++) {
    low.mass[i]=in_m_min;
    high.mass[i]=in_m_max;
    low.rad[i]=in_r_min;
    high.rad[i]=in_r_max;
  }

  // Prepare statistics
  init_grids_table(low,high);

  // Weights for each entry
  double w_current, w_next, w_best=0.0;

  // Warm-up flag
  bool warm_up=true;
  if (n_warm_up==0) warm_up=false;

  // Keep track of successful and failed MH moves
  mh_success_cnt=0;
  mh_failure_cnt=0;
  
  // Compute initial weight
  bool suc;
  w_current=compute_weight(e_current,*modp,ts,suc,wgts,warm_up);
  scr_out << "Initial weight: " << w_current << endl;
  if (w_current<=0.0) {
    for(size_t i=0;i<nsources;i++) {
      scr_out << i << " " << wgts[i] << endl;
    }
    scr_out << "Initial weight zero. Aborting." << endl;
    exit(-1);
  }
  
  {
    o2_shared_ptr<table_units<> >::type tab_eos;
    o2_shared_ptr<table_units<> >::type tab_mvsr;
    tab_eos=modp->cns.get_eos_results();
    tab_mvsr=ts->get_results();
    if (warm_up==false) {
      // Add the initial point if there's no warm up
      add_measurement(e_current,tab_eos,tab_mvsr,w_current,true,
		      mh_success_cnt,wgts);
    }
    e_best=e_current;
    w_best=w_current;
    output_best(fname_prefix,e_current,w_current,tab_eos,tab_mvsr,wgts);
  }

  // Initialize radii of e_next to zero
  for(size_t k=0;k<nsources;k++) {
    e_next.rad[k]=0.0;
  }

  // Keep track of total number of different points in the parameter
  // space that are considered (some of them may not result in TOV
  // calls because, for example, the EOS was acausal).
  int iteration=0;

  // Switch between two different Metropolis steps
  bool first_half=true;

  // Main loop
  bool main_done=false;
  size_t block_counter=0;
  while (!main_done) {

    // Make a step, ensure that we're in bounds and that
    // the masses are not too large
    bool bad_step;
    for(size_t k=0;k<e_next.np;k++) {
      size_t step_count=0;
      do {
	bad_step=false;
	e_next.params[k]=e_current.params[k]+(gr.random()*2.0-1.0)*
	  (high.params[k]-low.params[k])/step_fac;
	if (e_next.params[k]<low.params[k] || 
	    e_next.params[k]>high.params[k]) {
	  bad_step=true;
	}
	if (bad_step) {
	  scr_out << "Bad step: " << k << " " << low.params[k] << " "
		  << e_current.params[k] << " " << e_next.params[k] << " "
		  << high.params[k] << endl;
	}
      } while (bad_step && step_count<100);
      if (step_count==1e2) {
	scr_out << "Too many steps in parameter space failed." << endl;
	return gsl_efailed;
      }
    }
    
    select_mass(e_current,e_next,high.mass[0],bad_step);
    if (bad_step) return gsl_efailed;

    // Output the next point
    if (output_next) {
      scr_out << "Iteration, next: " << iteration << " " << e_next << endl;
    }
      
    // Compute next weight
    if (first_half) {
      w_next=compute_weight(e_next,*modp2,ts2,suc,wgts,warm_up);
    } else {
      w_next=compute_weight(e_next,*modp,ts,suc,wgts,warm_up);
    }
      
    // Test to ensure new point is good
    if (suc==true) {

      // Test radii
      for(size_t i=0;i<nsources;i++) {
	if (e_next.rad[i]>high.rad[i] || e_next.rad[i]<low.rad[i]) {
	  scr_out << "Rejected: Radius out of range." << endl;
	  suc=false;
	  i=nsources;
	}
      }
	
      // Ensure non-zero weight
      if (w_next==0.0) {
	scr_out << "Rejected: Zero weight." << endl;
	suc=false;
      }
	
    }

    // If the new point is still good, compare with
    // the Metropolis algorithm
    if (suc==true) {

      if (debug) {
	cout << first_half << " Next: " 
	     << e_next.params[0] << " " << w_next << endl;
      }
      
      bool accept=make_step(w_current,w_next,debug,warm_up,iteration);
      
      if (accept) {

	mh_success_cnt++;

	if (iteration>n_warm_up && warm_up==true) {
	  warm_up=false;
	  scr_out << "Setting warm_up to false. Reset start time." << endl;
	  if (true) {
	    max_time-=MPI_Wtime()-mpi_start_time;
	    scr_out << "Resetting max_time to : " << max_time << endl;
	  }
	  mpi_start_time=MPI_Wtime();
	  scr_out.precision(12);
	  scr_out << " Start time: " << mpi_start_time << endl;
	  scr_out.precision(6);
	}
	    
	// Add measurement from new point
	o2_shared_ptr<table_units<> >::type tab_eos;
	o2_shared_ptr<table_units<> >::type tab_mvsr;

	if (first_half) {
	  tab_eos=modp2->cns.get_eos_results();
	  tab_mvsr=ts2->get_results();
	} else {
	  tab_eos=modp->cns.get_eos_results();
	  tab_mvsr=ts->get_results();
	}
	tab_eos->set_interp_type(itp_linear);
	tab_mvsr->set_interp_type(itp_linear);

	// Store results from new point
	if (!warm_up) {
	  add_measurement(e_next,tab_eos,tab_mvsr,w_next,true,
			  mh_success_cnt,wgts);
	  if (debug) {
	    cout << first_half << " Adding new: " 
		 << e_next.params[0] << " " << w_next << " "
		 << tab_mvsr->max("gm") << endl;
	  }
	}

	// Output the new point
	scr_out << "MH Acc: " << mh_success_cnt << " " << e_next << " " 
		<< w_next << endl;
	  
	// Keep track of best point
	if (w_next>w_best) {
	  e_best=e_next;
	  w_best=w_next;
	  output_best(fname_prefix,e_best,w_best,tab_eos,tab_mvsr,wgts);
	}

	// Prepare for next point
	e_current=e_next;
	w_current=w_next;

	// Flip "first_half" parameter
	first_half=!(first_half);
	if (debug) cout << "Flip: " << first_half << endl;
	  
      } else {
	    
	// Point was rejected

	mh_failure_cnt++;

	o2_shared_ptr<table_units<> >::type tab_eos;
	o2_shared_ptr<table_units<> >::type tab_mvsr;
	
	if (first_half) {
	  tab_eos=modp->cns.get_eos_results();
	  tab_mvsr=ts->get_results();
	} else {
	  tab_eos=modp2->cns.get_eos_results();
	  tab_mvsr=ts2->get_results();
	}
	tab_eos->set_interp_type(itp_linear);
	tab_mvsr->set_interp_type(itp_linear);

	// Repeat measurement of old point
	if (!warm_up) {
	  add_measurement(e_current,tab_eos,tab_mvsr,w_current,false,
			  mh_success_cnt,wgts);
	  if (debug) {
	    cout << first_half << " Adding old: " 
		 << e_current.params[0] << " " << w_current << " "
		 << tab_mvsr->max("gm") << endl;
	  }
	}

	// Output the old point 
	scr_out << "MH Rej: " << mh_success_cnt << " " << e_current 
		<< " " << w_current << " " << w_next << endl;

	// Keep track of best point
	if (w_next>w_best) {
	  e_best=e_next;
	  w_best=w_next;
	  output_best(fname_prefix,e_best,w_best,tab_eos,tab_mvsr,wgts);
	  //cerr << "Best point with rejected step?" << endl;
	  //exit(-1);
	}

      }
	  
      // End of "if (suc==true)"
    }

    bool force_file_update=false;

    if (warm_up==false && iteration%10==9) {
    
      if (max_iters==0) {

	// If necessary, collect a block of histograms
	double elapsed=MPI_Wtime()-mpi_start_time;
	if (elapsed>max_time/((double)20)*((double)(block_counter+1)) ||
	    (elapsed>max_time && block_counter==19)) {
	  force_file_update=true;
	  block_counter++;
	  scr_out << "Finished block " << block_counter << " of 20." << endl;
	}
	
	// Check time
	scr_out << "Elapsed time: " << elapsed << " " << max_time << endl;
	if (elapsed>max_time) {
	  main_done=true;
	}

      } else {

	if (iteration+1>max_iters*(((int)block_counter)+1)/20) {
	  force_file_update=true;
	  block_counter++;
	  scr_out << "Iteration " << iteration+1 << " of " 
		  << max_iters << ", and finished block " 
		  << block_counter << " of 20." << endl;
	}

	if (iteration+1>max_iters) {
	  scr_out << "Iteration count, " << iteration 
		  << ", exceed maximum number, " << max_iters << "." << endl;
	  main_done=true;
	}
	
      }
     
    }

    // Store a copy of measurements in file
    if (!warm_up && (force_file_update || mh_success_cnt%10==9)) {
      update_files(fname_prefix,*modp,e_current);
    }
    
    // Increment iteration counter
    iteration++;

    // End of main loop
  }
   
  scr_out.close();
 
  return 0;
}

void bamr::setup_cli() {

  // ---------------------------------------
  // Set options
    
  static const int nopt=3;
  comm_option_s options[nopt]={
    {'m',"mcmc","Perform the Markov Chain Monte Carlo simulation.",
     1,1,"<filename prefix>",((string)"This is the main part of ")+
     "the code which performs the simulation. Make sure to set the "+
     "model first using the 'model' command and the input file "+
     "by setting the 'in_file' variable using the 'set' command. "+
     "The one required argument to mcmc is the prefix for the "+
     "output files.",
     new comm_option_mfptr<bamr>(this,&bamr::mcmc),
     cli::comm_option_both},
    {'o',"model","Choose model.",
     1,1,"<model name>",((string)"Choose the EOS parameterization model. ")+
     "Typical values are 'twop', 'impp', 'fixp', and 'genq'.",
     new comm_option_mfptr<bamr>(this,&bamr::set_model),
     cli::comm_option_both},
  };
  cl.set_comm_option_vec(nopt,options);

  // ---------------------------------------
  // Set parameters
    
  p_max_time.d=&max_time;
  p_max_time.help="Maximum run time in seconds (default 86400 sec or 1 day).";
  cl.par_list.insert(make_pair("max_time",&p_max_time));

  p_min_max_mass.d=&min_max_mass;
  p_min_max_mass.help=((string)"Minimum maximum mass ")
    +"(in solar masses, default 1.66).";
  cl.par_list.insert(make_pair("min_max_mass",&p_min_max_mass));

  p_min_mass.d=&min_mass;
  p_min_mass.help=((string)"Minimum possible mass for any of individual ")+
    "neutron stars in solar masses. The default is 0.8 solar masses.";
  cl.par_list.insert(make_pair("min_mass",&p_min_mass));

  p_step_fac.d=&step_fac;
  p_step_fac.help=((string)"MCMC step factor. The step size for each ")+
    "variable is taken as the difference between the high and low "+
    "limits divided by this factor (default 15.0). This factor can "+
    "be increased if the acceptance rate is too small, but care must "+
    "be taken, e.g. if the conditional probability is multimodal.";
  cl.par_list.insert(make_pair("step_fac",&p_step_fac));

  p_exit_mass.d=&exit_mass;
  p_exit_mass.help=((string)"Upper limit on maximum mass ")+
    "(default 10.0). When the maximum mass is larger than this value, "+
    "the current point in the parameter space is output to 'cout' and "+
    "execution is aborted. This is sometimes useful in debugging the "+
    "initial guess.";
  cl.par_list.insert(make_pair("exit_mass",&p_exit_mass));

  p_input_dist_thresh.d=&input_dist_thresh;
  p_input_dist_thresh.help=((string)"Input distribution threshold. ")+
    "This is the artificial lower limit for the (renormalized) "+
    "probability of a (R,M) pair as reported by the data file. If the "+
    "weight is smaller than or equal to this value, an exception is "+
    "thrown. Changing this value is sometimes "+
    "useful to gracefully avoid zero probabilities in the input "+
    "data files. The default is 0.";
  cl.par_list.insert(make_pair("input_dist_thresh",&p_input_dist_thresh));

  p_warm_up.i=&n_warm_up;
  p_warm_up.help=((string)"Minimum number of warm up iterations ")+
    "(default 500, can be zero).";
  cl.par_list.insert(make_pair("warm_up",&p_warm_up));

  p_user_seed.i=&user_seed;
  p_user_seed.help=((string)"Seed for multiplier for random number ")+
    "generator. If zero is given (the default), then mcmc() uses "+
    "time(0) to generate a random seed.";
  cl.par_list.insert(make_pair("user_seed",&p_user_seed));

  p_max_iters.i=&max_iters;
  p_max_iters.help=((string)"If non-zero, limit the number of ")+
    "iterations to be less than the specified number (default zero).";
  cl.par_list.insert(make_pair("max_iters",&p_max_iters));

  p_debug_star.b=&debug_star;
  p_debug_star.help=((string)"If true, output stellar properties ")+
    "to file with suffix '_scr' at each point (default false).";
  cl.par_list.insert(make_pair("debug_star",&p_debug_star));

  p_norm_max.b=&norm_max;
  p_norm_max.help=((string)"If true, normalize by max probability ")+
    "or if false, normalize by total integral (default true).";
  cl.par_list.insert(make_pair("norm_max",&p_norm_max));

  p_debug_load.b=&debug_load;
  p_debug_load.help=((string)"If true, output info on loaded data ")+
    "(default false).";
  cl.par_list.insert(make_pair("debug_load",&p_debug_load));
  
  p_debug_eos.b=&debug_eos;
  p_debug_eos.help=((string)"If true, output initial equation of state ")+
    "to file 'debug_eos.o2' and abort (default false).";
  cl.par_list.insert(make_pair("debug_eos",&p_debug_eos));

  p_output_next.b=&output_next;
  p_output_next.help=((string)"If true, output next point ")+
    "to the '_scr' file before calling TOV solver (default true).";
  cl.par_list.insert(make_pair("output_next",&p_output_next));

  p_baryon_density.b=&baryon_density;
  p_baryon_density.help=((string)"If true, compute baryon density ")+
    "and store histograms with baryon densities "+
    "(default true).";
  cl.par_list.insert(make_pair("baryon_density",&p_baryon_density));

  p_use_crust.b=&use_crust;
  p_use_crust.help=((string)"If true, use the default crust (default ")+
    "true).";
  cl.par_list.insert(make_pair("use_crust",&p_use_crust));

  p_in_file.str=&in_file;
  p_in_file.help="Input file name (default \"default.in\").";
  cl.par_list.insert(make_pair("in_file",&p_in_file));

  p_first_point_file.str=&first_point_file;
  p_first_point_file.help="First point file name (default \"\").";
  cl.par_list.insert(make_pair("first_point_file",&p_first_point_file));

  // --------------------------------------------------------

  p_nb_low.d=&nb_low;
  p_nb_low.help="Smallest baryon density grid point in 1/fm^3 (default 0.04).";
  cl.par_list.insert(make_pair("nb_low",&p_nb_low));

  p_nb_high.d=&nb_high;
  p_nb_high.help="Largest baryon density grid point in 1/fm^3 (default 1.04).";
  cl.par_list.insert(make_pair("nb_high",&p_nb_high));

  p_e_low.d=&e_low;
  p_e_low.help="Smallest energy density grid point in 1/fm^4 (default 0.3).";
  cl.par_list.insert(make_pair("e_low",&p_e_low));

  p_e_high.d=&e_high;
  p_e_high.help="Largest energy density grid point in 1/fm^4 (default 10.0).";
  cl.par_list.insert(make_pair("e_high",&p_e_high));
  
  p_m_low.d=&m_low;
  p_m_low.help="Smallest mass grid point in Msun (default 0.2).";
  cl.par_list.insert(make_pair("m_low",&p_m_low));

  p_m_high.d=&m_high;
  p_m_high.help="Largest mass grid point in Msun (default 3.0).";
  cl.par_list.insert(make_pair("m_high",&p_m_high));

  // --------------------------------------------------------
  
  return;
}

void bamr::run(int argc, char *argv[]) {
  
  // ---------------------------------------
  // Process command-line arguments and run
  
  setup_cli();
  
  cl.prompt="bamr> ";
  cl.run_auto(argc,argv);
    
  return;
}

