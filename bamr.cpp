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
#include "bamr.h"

#include <o2scl/vector.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

bamr_class::bamr_class() {

  prefix="bamr";
  model_type="";
  has_eos=true;
  schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;

  in_m_min=0.8;
  in_m_max=3.0;
  in_r_min=5.0;
  in_r_max=18.0;

  // Default parameter values

  grid_size=100;
  debug_load=false;
  // Minimum allowed maximum mass
  min_max_mass=2.0;
  // Minimum neutron star mass
  min_mass=0.8;
  debug_star=false;
  debug_line=false;
  debug_eos=false;
  baryon_density=true;
  exit_mass=10.0;
  input_dist_thresh=0.0;
  use_crust=true;
  best_detail=false;
  inc_baryon_mass=false;
  norm_max=true;
  mvsr_pr_inc=1.1;

  // -----------------------------------------------------------
  // Grid limits
  
  nb_low=0.04;
  nb_high=1.24;
  
  e_low=0.3;
  e_high=10.0;

  m_low=0.2;
  m_high=3.0;
  
  teos.verbose=0;
}

bamr_class::~bamr_class() {
}

void bamr_class::table_names_units(std::string &s, std::string &u) {
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
    s+=((string)"param_")+data_arr[0].modp->param_name(i)+" ";
    u+=data_arr[0].modp->param_unit(i)+" ";
  }

  // It is important here that all of these columns which store values
  // over a grid are either always positive or always negative,
  // because the code reports zero in the fill_line() function for
  // values beyond the end of the EOS or the M-R curve. 
  for(size_t i=0;i<nsources;i++) {
    s+=((string)"Rns_")+source_names[i]+" ";
    u+="km ";
  }
  for(size_t i=0;i<nsources;i++) {
    s+=((string)"Mns_")+source_names[i]+" ";
    u+="Msun ";
  }
  if (has_eos) {
    for(int i=0;i<grid_size;i++) {
      s+=((string)"P_")+szttos(i)+" ";
      u+="1/fm^4 ";
    }
  }
  for(int i=0;i<grid_size;i++) {
    s+=((string)"R_")+szttos(i)+" ";
    u+="km ";
    if (has_eos) {
      s+=((string)"PM_")+szttos(i)+" ";
      u+="1/fm^4 ";
    }
  }
  if (has_eos) {
    if (baryon_density) {
      for(int i=0;i<grid_size;i++) {
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
    s+="R_max M_max P_max e_max ";
    u+="km Msun 1/fm^4 1/fm^4 ";
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

void bamr_class::init_grids_table(ubvector &low, ubvector &high) {
  
  if (low.np==0 || high.np==0) {
    O2SCL_ERR("No parameters in bamr_class::init().",exc_einval);
  }

  // -----------------------------------------------------------
  // Make grids

  nb_grid=uniform_grid_end<double>(nb_low,nb_high,grid_size-1);
  e_grid=uniform_grid_end<double>(e_low,e_high,grid_size-1);
  m_grid=uniform_grid_end<double>(m_low,m_high,grid_size-1);

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
      O2SCL_ERR("Column/unit alignment in bamr_class::init_grids_table().",
		exc_esanity);
    }
  }

  return;
}

void bamr_class::fill_line
(ubvector &e, std::shared_ptr<o2scl::table_units<> > tab_eos,
 std::shared_ptr<o2scl::table_units<> > tab_mvsr,
 double weight, bool new_meas, size_t n_meas, ubvector &wgts,
 std::vector<double> &line) {

  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (has_eos) {

    // The highest baryon density in the EOS table
    nbmax2=tab_eos->max("nb");
    // The central energy density in the maximum mass configuration
    emax=tab_mvsr->max("ed");
    // The central pressure in the maximum mass configuration
    pmax=tab_mvsr->max("pr");
    // The maximum mass
    mmax=tab_mvsr->get_constant("new_max");
    // The radius of the maximum mass star
    rmax=tab_mvsr->get_constant("new_r_max");
    
    if (baryon_density) {
      // The central baryon density in the maximum mass configuration
      nbmax=tab_mvsr->get_constant("new_nb_max");
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
    for(int i=0;i<grid_size;i++) {
      double eval=e_grid[i];
      // Make sure the energy density from the grid
      // isn't beyond the table limit
      double emax2=tab_eos->max("ed");
      if (eval<emax2) {
	double pres_temp=tab_eos->interp("ed",eval,"pr");
	//if (pres_temp<pmax) {
	line.push_back(pres_temp);
	//} else {
	//line.push_back(0.0);
	//}
      } else {
	line.push_back(0.0);
      }
    }
  }

  // It is important here that all of these columns which store values
  // over a grid are either always positive or always negative,
  // because the code reports zero in the fill_line() function for
  // values beyond the end of the EOS or the M-R curve. 
  for(int i=0;i<grid_size;i++) {
    double mval=m_grid[i];
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
      for(int i=0;i<grid_size;i++) {
	double nbval=nb_grid[i];
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
    line.push_back(rmax);
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

void bamr_class::add_measurement
(ubvector &e, std::shared_ptr<o2scl::table_units<> > tab_eos,
 std::shared_ptr<o2scl::table_units<> > tab_mvsr,
 double weight, bool new_meas, size_t n_meas, ubvector &wgts) {

  // Test to see if we need to add a new line of data or
  // increment the weight on the previous line
  if (tc.get_nlines()==0 || new_meas==true) {
    
    std::vector<double> line;
    fill_line(e,tab_eos,tab_mvsr,weight,new_meas,n_meas,wgts,line);
    
    // Done adding values, check size and add to table
    if (line.size()!=tc.get_ncolumns()) {
      scr_out << "line.size(): " << line.size() << endl;
      scr_out << "tc.get_ncolumns(): " << tc.get_ncolumns() << endl;
      for(size_t i=0;i<line.size() && i<tc.get_ncolumns();i++) {
	scr_out << line[i] << " " << tc.get_column_name(i) << endl;
      }
      O2SCL_ERR("Alignment problem.",exc_efailed);
    }
    tc.line_of_data(line.size(),line);

    if (debug_line) {
      vector<string> sc_in, sc_out;
      for(size_t k=0;k<line.size();k++) {
	sc_in.push_back(tc.get_column_name(k)+": "+o2scl::dtos(line[k]));
      }
      o2scl::screenify(line.size(),sc_in,sc_out);
      for(size_t k=0;k<sc_out.size();k++) {
	cout << sc_out[k] << endl;
      }
      cout << "Press a key and enter to continue." << endl;
      char ch;
      cin >> ch;
    }
    
  } else if (tc.get_nlines()>0) {
    tc.set("mult",tc.get_nlines()-1,
	   tc.get("mult",tc.get_nlines()-1)+1.0);
  }
  
  return;
}

void bamr_class::first_update(hdf_file &hf, model &modp) {

  vector<string> param_names;
  for(size_t i=0;i<nparams;i++) {
    param_names.push_back(modp.param_name(i));
  }
  hf.sets_vec("param_names",param_names);

  hf.sets_vec("source_names",source_names);
  hf.sets_vec("source_fnames",source_fnames);
  hf.sets_vec("slice_names",slice_names);

  hf.set_szt("grid_size",grid_size);
  hf.set_szt("nparams",nparams);
  hf.set_szt("nsources",nsources);
  hf.sets("model",model_type);
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
  hf.seti("debug_line",debug_line);
  hf.seti("debug_eos",debug_eos);
  hf.seti("debug_star",debug_star);
  hf.seti("best_detail",best_detail);
  hf.seti("file_update_iters",file_update_iters);
  hf.seti("inc_baryon_mass",inc_baryon_mass);
  hf.seti("output_next",output_next);
  hf.setd("nb_low",nb_low);
  hf.setd("nb_high",nb_high);
  hf.setd("e_low",e_low);
  hf.setd("e_high",e_high);
  hf.setd("m_low",m_low);
  hf.setd("m_high",m_high);
  hf.seti("first_point_type",first_point_type);
  hf.sets("first_point_file",first_point_file);
  hf.setd_vec_copy("first_point",first_point);

  hdf_output(hf,nb_grid,"nb_grid");
  hdf_output(hf,e_grid,"e_grid");
  hdf_output(hf,m_grid,"m_grid");
    
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

  hf.sets_vec("cl_args",cl_args);

  return;
}

void bamr_class::update_files(model &modp, ubvector &e_current) {

  hdf_file hf;

  // Open main update file
  hf.open_or_create(prefix+"_"+std::to_string(mpi_rank)+"_out");
    
  // First time, output some initial quantities
  if (first_file_update==false) {
    first_update(hf,modp);
    first_file_update=true;
  }

  hf.set_szt("mh_success",mh_success);
  hf.set_szt("mh_failure",mh_failure);
  hf.set_szt("mcmc_iterations",mcmc_iterations);
  hf.seti_vec("ret_codes",ret_codes);
  
  // Store Markov chain
  if (n_chains==0) n_chains++;
  hf.set_szt("n_chains",n_chains);
  string ch_name="markov_chain"+szttos(n_chains-1);
  hdf_output(hf,tc,ch_name);
  if (((int)tc.get_nlines())==max_chain_size) {
    tc.clear_data();
    n_chains++;

    // 03/04/16 - I don't think this is necessary
    // Store the new empty chain in the HDF5 file just in case we
    // stop before the next call to update_files().
    //string ch_name="markov_chain"+szttos(n_chains-1);
    //hdf_output(hf,tc,ch_name);
  }

  hf.close();

  return;
}

void bamr_class::load_mc() {

  double tot, max;

  string name;
    
  if (nsources>0) {

    source_tables.resize(nsources);

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

      hdf_file hf;
      hf.open(source_fnames[file]);
      if (table_names[file].length()>0) {
	hdf_input(hf,source_tables[file],table_names[file]);
      } else {
	hdf_input(hf,source_tables[file]);
      }
      hf.close();
      
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
      
      hdf_file hf;
      hf.open(source_fnames[k]);
      if (table_names[k].length()>0) {
	hdf_input(hf,source_tables[k],table_names[k]);
      } else {
	hdf_input(hf,source_tables[k]);
      }
      hf.close();
    }
    
#endif
    
    scr_out << "\nInput data files: " << endl;
    
    if (norm_max) {
      scr_out << "Normalizing maximum probability to 1." << endl;
    } else {
      scr_out << "Normalizing integral of distribution to 1." << endl;
    }

    scr_out << "File                          name   total        "
	    << "max          P(10,1.4)" << endl;

    for(size_t k=0;k<nsources;k++) {
      
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
      scr_out.width(29);
      string stempx=source_fnames[k].substr(0,29);
      scr_out << stempx << " ";
      scr_out.width(6);
      scr_out << source_names[k] << " " << tot << " " << max << " ";
      scr_out.unsetf(ios::left);
      scr_out << source_tables[k].interp(10.0,1.4,slice_names[k]) << endl;
      
    }
    
    scr_out << endl;
  }

  if (in_m_min<min_mass) in_m_min=min_mass;
  
#ifdef AWS_HACK
  in_m_min=1.3;
  in_m_max=1.5;
#endif
  
  scr_out << "M limits: (" 
	  << in_m_min << "," << in_m_max << ")" << endl;
  scr_out << "R limits: ("
	  << in_r_min << "," << in_r_max << ")" << endl;
  scr_out << endl;
  
  return;
}
  
void bamr_class::prepare_eos(ubvector &e, eos_tov &dat, int &success) {
  return;
}

void bamr_class::compute_star(ubvector &e, eos_tov &dat, int &success) {
			      
  
  success=ix_success;

  // Compute the EOS first
  if (has_eos) {
    dat.modp->compute_eos(e,success,scr_out);
    if (success!=ix_success) return;
  }
  
  // Ensure we're using linear interpolation
  shared_ptr<table_units<> > tab_eos=dat.modp->cns.get_eos_results();

  if (has_eos) {
    tab_eos->set_interp_type(itp_linear);

    // Check that pressure is increasing. If we're using a crust EOS,
    // choose 0.6 as an arbitrary low-density cutoff which corresponds
    // to about n_B=0.12 fm^{-3}, and check the low-density part below
    // instead.
  
    for(size_t i=0;(tab_eos->get_nlines()>0 && i<tab_eos->get_nlines()-1);
	i++) {
      if ((!use_crust || tab_eos->get("ed",i)>0.6) && 
	  tab_eos->get("pr",i+1)<tab_eos->get("pr",i)) {
	scr_out << "Rejected: Pressure decreasing." << endl;
	scr_out << "ed=" << tab_eos->get("ed",i) 
		<< " pr=" << tab_eos->get("pr",i) << endl;
	scr_out << "ed=" << tab_eos->get("ed",i+1) 
		<< " pr=" << tab_eos->get("pr",i+1) << endl;
	success=ix_press_dec;
	return;
      }
    }
  }

  // If requested, compute the baryon density automatically
  if (has_eos && baryon_density && !tab_eos->is_column("nb")) {
    
    // Obtain the baryon density calibration point from the model
    double n1, e1;
    dat.modp->baryon_density_point(n1,e1);
    
    if (n1<=0.0 && e1<=0.0) {
      O2SCL_ERR2("Computing the baryon density requires one ",
		 "calibration point in bamr_class::compute_star().",exc_einval);
    }

    // Compute inverse of gibbs energy density, 'igb'

    tab_eos->new_column("igb");
    tab_eos->set_unit("igb","fm^4");

    for(size_t i=0;i<tab_eos->get_nlines();i++) {

      if (tab_eos->get("ed",i)+tab_eos->get("pr",i)<=0.0 ||
	  !std::isfinite(tab_eos->get("ed",i)) ||
	  !std::isfinite(tab_eos->get("pr",i))) {
	scr_out << "Inverse Gibbs not finite." << endl;
	scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	scr_out << "ed pr" << endl;
	for(size_t i=0;i<tab_eos->get_nlines();i++) {
	  scr_out << i << " "
		  << tab_eos->get("ed",i) << " "
		  << tab_eos->get("pr",i) << endl;
	}
	O2SCL_ERR("Inverse Gibbs not finite.",exc_efailed);
      }
      tab_eos->set("igb",i,1.0/(tab_eos->get("ed",i)+tab_eos->get("pr",i)));
    }

    // Compute integral of 'igb' relative to ed='e1', called 'iigb'

    tab_eos->new_column("iigb");

    for(size_t i=0;i<tab_eos->get_nlines();i++) {
      if (e1<=tab_eos->get("ed",i)) {
	double val=tab_eos->integ("ed",e1,tab_eos->get("ed",i),"igb");
	if (!std::isfinite(val)) {
	  scr_out << "Baryon integral not finite." << endl;
	  scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	  scr_out << "ed pr" << endl;
	  for(size_t i=0;i<tab_eos->get_nlines();i++) {
	    scr_out << i << " "
		    << tab_eos->get("ed",i) << " "
		    << tab_eos->get("pr",i) << endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",exc_efailed);
	}
	tab_eos->set("iigb",i,val);
      } else {
	double val=-tab_eos->integ("ed",tab_eos->get("ed",i),e1,"igb");
	if (!std::isfinite(val)) {
	  scr_out << "Baryon integral not finite (2)." << endl;
	  scr_out << "n1=" << n1 << " e1=" << e1 << endl;
	  scr_out << "ed pr" << endl;
	  for(size_t i=0;i<tab_eos->get_nlines();i++) {
	    scr_out << i << " "
		    << tab_eos->get("ed",i) << " "
		    << tab_eos->get("pr",i) << endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",exc_efailed);
	}
	tab_eos->set("iigb",i,val);
      }
    }

    // Compute normalization constant
    double Anb=n1/exp(tab_eos->interp("ed",e1,"iigb"));
    if (!std::isfinite(Anb) || Anb<0.0) {
      scr_out << "Baryon density normalization problem." << endl;
      success=ix_nb_problem;
      return;
    }

    // Now compute baryon density

    tab_eos->new_column("nb");
    tab_eos->set_unit("nb","1/fm^3");

    for(size_t i=0;i<tab_eos->get_nlines();i++) {      

      // If the density is too low, then just use zero baryon density.
      // This pressure (10^{-5} fm^{-4}) corresponds to a baryon 
      // density of about 3e-3 fm^{-3}.

      if (use_crust==false && tab_eos->get("pr",i)<1.0e-5) {

	tab_eos->set("nb",i,0.0);

      } else {
	
	double nbt=Anb*exp(tab_eos->get("iigb",i));
	if (!std::isfinite(nbt)) {
	  scr_out << "Baryon density normalization problem (2)." << endl;
	  success=ix_nb_problem2;
	  return;
	} 
	tab_eos->set("nb",i,nbt);
      }
    }
    
    // End of loop 'if (has_eos && baryon_density && 
    // !tab_eos->is_column("nb")) {' 
  }

  shared_ptr<table_units<> > tab_mvsr;

  if (has_eos) {

    // Perform any additional necessary EOS preparations
    prepare_eos(e,dat,success);
    if (success!=ix_success) {
      return;
    }

    // Read the EOS into the tov_eos object.
    if (baryon_density && inc_baryon_mass) {
      tab_eos->set_unit("ed","1/fm^4");
      tab_eos->set_unit("pr","1/fm^4");
      tab_eos->set_unit("nb","1/fm^3");
      teos.read_table(*tab_eos,"ed","pr","nb");
    } else {
      tab_eos->set_unit("ed","1/fm^4");
      tab_eos->set_unit("pr","1/fm^4");
      teos.read_table(*tab_eos,"ed","pr");
    }
    
    if (use_crust) {
    
      double ed_last=0.0;
      // This range corresponds to between about n_B=0.01 and 0.17
      // fm^{-3}
      for(double pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	double ed, nb;
	teos.ed_nb_from_pr(pr,ed,nb);
	if (ed_last>1.0e-20 && ed<ed_last) {
	  scr_out << "Stability problem near crust-core transition." << endl;
	  if (has_esym) {
	    scr_out << "S=" << tab_eos->get_constant("S")*hc_mev_fm 
		    << " L=" << tab_eos->get_constant("L")*hc_mev_fm << endl;
	  }
	  scr_out << "Energy decreased with increasing pressure "
		  << "at pr=" << pr << endl;
	  scr_out << endl;
	  scr_out << "Full EOS near transition: " << endl;
	  scr_out << "pr ed" << endl;
	  for(pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	    teos.ed_nb_from_pr(pr,ed,nb);
	    scr_out << pr << " " << ed << endl;
	  }
	  scr_out << endl;
	  success=ix_crust_unstable;
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
	teos.get_eden_user(pr,ed,nb);
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
    dat.ts.princ=mvsr_pr_inc;
    int info=dat.ts.mvsr();
    if (info!=0) {
      scr_out << "M vs. R failed: info=" << info << endl;
      success=ix_mvsr_failed;
      return;
    }
    tab_mvsr=dat.ts.get_results();
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
      success=ix_tov_failure;
      return;
    }

    // Check that maximum mass is large enough,
    double mmax=tab_mvsr->max("gm");
    if (mmax<min_max_mass) {
      scr_out << "Maximum mass too small: " << mmax << " < "
	      << min_max_mass << "." << endl;
      success=ix_small_max;
      return;
    }

    // Check the radius of the maximum mass star
    size_t ix_max=tab_mvsr->lookup("gm",mmax);
    if (tab_mvsr->get("r",ix_max)>1.0e4) {
      scr_out << "TOV convergence problem: " << endl;
      success=ix_tov_conv;
      return;
    }

    // Check that all stars have masses below the maximum mass
    // of the current M vs R curve
    bool mass_fail=false;
    for(size_t i=0;i<nsources;i++) {
      if (e.mass[i]>mmax) mass_fail=true;
    }

    // If a star is too massive, readjust it's mass accordingly
    if (mass_fail==true) {
      scr_out << "Adjusting masses for M_{max} = " << mmax << endl;
      scr_out << "Old entry: " << e << endl;
      select_mass(e,e,mmax);
      scr_out << "New entry: " << e << endl;
    }
    
    tab_mvsr->add_constant
      ("new_max",vector_max_quad<vector<double>,double>
       (tab_mvsr->get_nlines(),(*tab_mvsr)["r"],(*tab_mvsr)["gm"]));
    /*
      if (true) {
      size_t ix=vector_max_index<vector<double>,double>
      (tab_mvsr->get_nlines(),(*tab_mvsr)["gm"]);
      if (ix!=0 && ix<tab_mvsr->get_nlines()-1) {
      scr_out << tab_mvsr->get("gm",ix-1) << " ";
      scr_out << tab_mvsr->get("r",ix-1) << " ";
      scr_out << tab_mvsr->get("nb",ix-1) << endl;
      scr_out << tab_mvsr->get("gm",ix) << " ";
      scr_out << tab_mvsr->get("r",ix) << " ";
      scr_out << tab_mvsr->get("nb",ix) << endl;
      scr_out << tab_mvsr->get("gm",ix+1) << " ";
      scr_out << tab_mvsr->get("r",ix+1) << " ";
      scr_out << tab_mvsr->get("nb",ix+1) << endl;
      }
      }
    */
    
    tab_mvsr->add_constant
      ("new_r_max",vector_max_quad_loc<vector<double>,double>
       (tab_mvsr->get_nlines(),(*tab_mvsr)["r"],(*tab_mvsr)["gm"]));

    if (baryon_density) {
      
      double nb1=tab_mvsr->get("nb",tab_mvsr->lookup("gm",mmax));
      if (nb1>0.01) {
	double nb2=vector_max_quad_loc<vector<double>,double>
	  (tab_mvsr->get_nlines(),(*tab_mvsr)["nb"],(*tab_mvsr)["gm"]);
	tab_mvsr->add_constant("new_nb_max",nb2);
	if (nb2>5.0) {
	  scr_out << "nb_check: " << nb2 << " " << nb1 << endl;
	  for(size_t i=0;i<tab_mvsr->get_nlines();i++) {
	    scr_out << i << " " << tab_mvsr->get("gm",i) << " "
		    << tab_mvsr->get("r",i) << " " 
		    << tab_mvsr->get("nb",i) << endl;
	  }
	}
      } else {
	tab_mvsr->add_constant("new_nb_max",nb1);
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
      success=ix_mvsr_table;
      return;
    }

    // Compute speed of sound squared
    tab_mvsr->deriv("ed","pr","dpde");

  } else {

    // If there's no EOS, then the model object gives the M-R curve
    tab_mvsr=dat.ts.get_results();
    dat.modp->compute_mr(e,scr_out,tab_mvsr,success);
    if (success!=ix_success) {
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
	success=ix_acausal;
	return;
      }
    }

  } else {

    for(size_t i=0;i<nsources;i++) {
      if (e.rad[i]<2.94*schwarz_km/2.0*e.mass[i]) {
	scr_out << "Source " << source_names[i] << " acausal." << endl;
	success=ix_acausal_mr;
	return;
      }
    }

  }

  return;
}

double bamr_class::compute_weight(ubvector &e, eos_tov &dat, 
				  int &success, ubvector &wgts, bool warm_up) {
			     
  // Compute the M vs R curve and return if it failed
  compute_star(e,dat,success);
  if (success!=ix_success) {
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
    success=ix_mr_outside;
    return 0.0;
  }

  success=ix_success;
  double ret=1.0;

  shared_ptr<table_units<> > tab_mvsr=dat.ts.get_results();
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
    scr_out << "Exiting because maximum mass (" << m_max_current 
	    << ") larger than exit_mass (" << exit_mass << ")." 
	    << endl;
    scr_out.precision(12);
    scr_out << "e,ret: " << e << " " << ret << endl;
    scr_out.precision(6);
    exit(-1);
  }

  return ret;
}

int bamr_class::add_data(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<5) {
    cout << "Not enough arguments given to 'add-data'." << endl;
    return exc_efailed;
  }

  source_names.push_back(sv[1]);
  source_fnames.push_back(sv[2]);
  slice_names.push_back(sv[3]);
  first_mass.push_back(o2scl::stod(sv[4]));
  if (sv.size()==6) {
    table_names.push_back(sv[5]);
  } else {
    table_names.push_back("");
  }

  nsources++;

  return 0;
}

int bamr_class::hastings(std::vector<std::string> &sv, 
			 bool itive_com) {

  bool debug=true;

  if (file_opened==false) {
    // Open main output file
    scr_out.open((prefix+"_"+std::to_string(mpi_rank)+"_scr").c_str());
    scr_out.setf(ios::scientific);
    file_opened=true;
    scr_out << "Opened main file in command 'hastings'." << endl;
  }

  if (sv.size()<2) {
    cout << "No arguments given to 'hastings'." << endl;
    return exc_efailed;
  }

  if (model_type.length()==0) {
    cout << "No model selected in 'hastings'." << endl;
    return exc_efailed;
  }

#ifndef BAMR_NO_MPI
  int buffer=0, tag=0;
  if (mpi_nprocs>1 && mpi_rank>0) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
  }
#endif
  
  // Read the data file
  std::string fname=sv[1];
  scr_out << "Opening file " << fname << " for hastings." << endl;
  hdf_file hf;
  hf.open(fname);
  table_units<> file_tab;
  hdf_input(hf,file_tab,"markov_chain0");
  hf.close();
  scr_out << "Done opening file " << fname << " for hastings." << endl;

#ifndef BAMR_NO_MPI
  if (mpi_nprocs>1 && mpi_rank<mpi_nprocs-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,tag,MPI_COMM_WORLD);
  }
#endif

  // Create a new column equal to mult times weight
  file_tab.function_column("mult*weight","mwgt");
  
  // Remove
  double max_mwgt=file_tab.max("mwgt");
  if (debug) scr_out << "lines: " << file_tab.get_nlines() << endl;
  file_tab.add_constant("max_mwgt",max_mwgt);
  file_tab.delete_rows("mwgt<0.1*max_mwgt");
  if (debug) scr_out << "lines: " << file_tab.get_nlines() << endl;
  
  // The total number of variables
  size_t nv=nparams+nsources;
  if (debug) {
    scr_out << nparams << " parameters and " << nsources << " sources."
	 << endl;
  }
  hg_best.resize(nv);
  
  // Find the average values
  for(size_t i=0;i<nparams;i++) {
    string str_i=((string)"param_")+data_arr[0].modp->param_name(i);
    hg_best[i]=wvector_mean(file_tab.get_nlines(),file_tab[str_i],
			    file_tab["mwgt"]);
  }
  for(size_t i=0;i<nsources;i++) {
    string str_i=((string)"Mns_")+source_names[i];
    hg_best[i+nparams]=wvector_mean(file_tab.get_nlines(),file_tab[str_i],
				    file_tab["mwgt"]);
  }
  
  // Construct the covariance matrix
  ubmatrix covar(nv,nv);
  for(size_t i=0;i<nparams;i++) {
    string str_i=((string)"param_")+data_arr[0].modp->param_name(i);
    for(size_t j=i;j<nparams;j++) {
      string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
      covar(i,j)=wvector_covariance(file_tab.get_nlines(),
				    file_tab[str_i],file_tab[str_j],
				    file_tab["mult"]);
      if (debug) {
	scr_out << "Covar: " << i << " " << j << " "
	     << covar(i,j) << endl;
      }
      covar(j,i)=covar(i,j);
    }
    for(size_t j=0;j<nsources;j++) {
      string str_j=((string)"Mns_")+source_names[j];
      covar(i,j+nparams)=wvector_covariance(file_tab.get_nlines(),
					    file_tab[str_i],file_tab[str_j],
					    file_tab["mult"]);
      if (debug) {
	scr_out << "Covar: " << i << " " << j+nparams << " "
	     << covar(i,j+nparams) << endl;
      }
      covar(j+nparams,i)=covar(i,j+nparams);
    }
  }
  for(size_t i=0;i<nsources;i++) {
    string str_i=((string)"Mns_")+source_names[i];
    for(size_t j=i;j<nsources;j++) {
      string str_j=((string)"Mns_")+source_names[j];
      covar(i+nparams,j+nparams)=wvector_covariance(file_tab.get_nlines(),
				    file_tab[str_i],file_tab[str_j],
				    file_tab["mult"]);
      if (debug) {
	scr_out << "Covar: " << i+nparams << " " << j+nparams << " "
	     << covar(i+nparams,j+nparams) << endl;
      }
      covar(j+nparams,i+nparams)=covar(i+nparams,j+nparams);
    }
  }

  // Perform the Cholesky decomposition
  hg_chol=covar;
  o2scl_linalg::cholesky_decomp(nv,hg_chol);

  // Find the inverse
  hg_covar_inv=hg_chol;
  o2scl_linalg::cholesky_invert<ubmatrix>(nv,hg_covar_inv);
  
  // Force hg_chol to be lower triangular
  for(size_t i=0;i<nv;i++) {
    for(size_t j=0;j<nv;j++) {
      if (i<j) hg_chol(i,j)=0.0;
    }
  }

  // Compute the normalization, weighted by the likelihood function
  hg_norm=1.0;
  size_t step=file_tab.get_nlines()/20;
  if (step<1) step=1;
  double renorm=0.0;
  double wgt_sum=0.0;
  for(size_t i=0;i<file_tab.get_nlines();i+=step) {
    ubvector e(nparams,nsources);
    for(size_t j=0;j<nparams;j++) {
      string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
      e.params[j]=file_tab.get(str_j,i);
    }
    for(size_t j=0;j<nsources;j++) {
      string str_j=((string)"Mns_")+source_names[j];
      e.mass[j]=file_tab.get(str_j,i);
    }
    double wgt=file_tab.get("mult",i)*file_tab.get("weight",i);
    double rat=wgt/approx_like(e);
    renorm+=wgt*wgt/approx_like(e);
    if (debug) {
      scr_out << wgt << " " << approx_like(e) << " " << rat << endl;
    }
    wgt_sum+=wgt;
  }
  renorm/=((double)wgt_sum);
  hg_norm*=renorm;
  if (debug) {
    scr_out << "New normalization: " << hg_norm << endl;
  }

  step=file_tab.get_nlines()/20;
  if (step<1) step=1;
  for(size_t i=0;i<file_tab.get_nlines();i+=step) {
    ubvector e(nparams,nsources);
    for(size_t j=0;j<nparams;j++) {
      string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
      e.params[j]=file_tab.get(str_j,i);
    }
    for(size_t j=0;j<nsources;j++) {
      string str_j=((string)"Mns_")+source_names[j];
      e.mass[j]=file_tab.get(str_j,i);
    }
    double wgt=file_tab.get("mult",i)*file_tab.get("weight",i);
    double rat=wgt/approx_like(e);
    if (debug) {
      scr_out << wgt << " " << approx_like(e) << " " << rat << endl;
    }
  }
  hg_mode=1;

  return 0;
}

int bamr_class::set_first_point(std::vector<std::string> &sv, 
				bool itive_com) {

  if (sv.size()<2) {
    cout << "No arguments given to 'first-point'." << endl;
    return exc_efailed;
  }

  if (sv[1]==((string)"values")) {

    first_point.resize(sv.size()-1);
    for(size_t i=2;i<sv.size();i++) {
      first_point[i-2]=o2scl::function_to_double(sv[i]);
    }
    first_point_type=fp_unspecified;

  } else if (sv[1]==((string)"prefix")) {
  
    first_point_type=fp_last;
    first_point_file=sv[2]+((std::string)"_")+std::to_string(mpi_rank)+"_out";

  } else if (sv[1]==((string)"last")) {
    first_point_type=fp_last;
    first_point_file=sv[2];
  } else if (sv[1]==((string)"best")) {
    first_point_type=fp_best;
    first_point_file=sv[2];
  } else if (sv[1]==((string)"N")) {
    first_point_type=o2scl::stoi(sv[2]);
    first_point_file=sv[3];
  }

  return 0;
}

int bamr_class::set_model(std::vector<std::string> &sv, bool itive_com) {
  // We cannot use scr_out here because it isn't set until the call
  // to mcmc().
  if (sv.size()<2) {
    cerr << "Model name not given." << endl;
    return exc_efailed;
  }
  model_type=sv[1];
  if (data_arr.size()>0) {
    data_arr[0].modp->remove_params(cl);
    for(size_t i=0;i<data_arr.size();i++) {
      delete data_arr[i].modp;
    }
    data_arr.clear();
  }
  if (sv[1]==((string)"twop")) {
    nparams=8;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new two_polytropes;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new two_polytropes;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"altp")) {
    nparams=8;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new alt_polytropes;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new alt_polytropes;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"fixp")) {
    nparams=8;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new fixed_pressure;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new fixed_pressure;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"qstar")) {
    nparams=4;
    has_esym=false;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new quark_star;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new quark_star;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"genq")) {
    nparams=9;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new generic_quarks;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new generic_quarks;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"qmc")) {
    nparams=7;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new qmc_neut;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new qmc_neut;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"qmc_threep")) {
    nparams=9;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new qmc_threep;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new qmc_threep;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"qmc_fixp")) {
    nparams=8;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new qmc_fixp;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new qmc_fixp;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else if (sv[1]==((string)"qmc_twolines")) {
    nparams=8;
    has_esym=true;
    has_eos=true;

    if (use_smove) {
      data_arr.resize(2*nwalk);
      for(size_t i=0;i<2*nwalk;i++) {
	data_arr[i].modp=new qmc_twolines;
	data_arr[i].ts.set_eos(teos);
      }
    } else {
      data_arr.resize(2);
      for(size_t i=0;i<2;i++) {
	data_arr[i].modp=new qmc_twolines;
	data_arr[i].ts.set_eos(teos);
      }
    }

  } else {
    cerr << "Model unknown." << endl;
    nparams=0;
    model_type="";
    has_esym=true;
    has_eos=true;
    return exc_efailed;
  }
  data_arr[0].modp->setup_params(cl);

  return 0;
}

void bamr_class::output_best(ubvector &e_best, double w_best,
			     shared_ptr<table_units<> > tab_eos,
			     shared_ptr<table_units<> > tab_mvsr,
			     ubvector &wgts) {
  
  scr_out << "Best: " << e_best << " " << w_best << endl;

  string fname_best_out=prefix+"_"+std::to_string(mpi_rank)+"_out";
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

  if (best_detail) {

    // "Best" EOS
    hdf_output(hf,*tab_eos,"best_eos");
    
    // "Best" M vs. R curve
    hdf_output(hf,*tab_mvsr,"best_mvsr");

  }

  hf.close();

  return;
}

int bamr_class::mcmc_init() {
  return 0;
}

void bamr_class::setup_cli() {

  // ---------------------------------------
  // Set options
    
  static const int nopt=5;
  comm_option_s options[nopt]={
    {'m',"mcmc","Perform the Markov Chain Monte Carlo simulation.",
     0,0,"",((string)"This is the main part of ")+
     "the code which performs the simulation. Make sure to set the "+
     "model first using the 'model' command first.",
     new comm_option_mfptr<bamr_class>(this,&bamr_class::mcmc),
     cli::comm_option_both},
    {'o',"model","Choose model.",
     1,1,"<model name>",((string)"Choose the EOS parameterization model. ")+
     "Possible values are 'twop', 'altp', 'fixp', 'genq', 'qstar', "+
     "'qmc', 'qmc_threep' ,'qmc_fixp', and 'qmc_twolines'. A "+
     "model must be chosen before a MCMC run.",
     new comm_option_mfptr<bamr_class>(this,&bamr_class::set_model),
     cli::comm_option_both},
    {'a',"add-data","Add data source to the list.",
     4,5,"<name> <file> <slice> <initial mass> [obj name]",
     ((string)"Specify data as a table3d object in a HDF5 file. ")+
     "The string <name> is the name used, <file> is the filename, "+
     "<slice> is the name of the slice in the table3d object, "+
     "<initial mass> is the initial mass for the first point, and "+
     "[obj name] is the optional name of table3d object in <file>. "+
     "If [obj name] is not specified, then the first table3d object "+
     "is used.",new comm_option_mfptr<bamr_class>(this,&bamr_class::add_data),
     cli::comm_option_both},
    {'f',"first-point","Set the starting point in the parameter space",
     1,-1,"<mode> [...]",
     ((string)"Mode can be one of 'best', 'last', 'N', or 'values'. ")+
     "If mode is 'best', then it uses the point with the largest "+
     "weight and the second argument specifies the file. If mode is "+
     "'last' then it uses the last point and the second argument "+
     "specifies the file. If mode is 'N' then it uses the Nth point, "+
     "the second argument specifies the value of N and the third "+
     "argument specifies the file. If mode is 'values', then the remaining "+
     "arguments specify all the parameter values. On the command-line, "+
     "enclose negative values in quotes and parentheses, i.e. \"(-1.00)\" "+
     "to ensure they do not get confused with other options.",
     new comm_option_mfptr<bamr_class>(this,&bamr_class::set_first_point),
     cli::comm_option_both},
    {'s',"hastings","Specify distribution for M-H step",
     1,1,"<filename>",
     ((string)"Desc. ")+"Desc2.",
     new comm_option_mfptr<bamr_class>(this,&bamr_class::hastings),
     cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);

  // ---------------------------------------
  // Set parameters
    
  p_grid_size.i=&grid_size;
  p_grid_size.help="Grid size (default 100).";
  cl.par_list.insert(make_pair("grid_size",&p_grid_size));

  p_min_max_mass.d=&min_max_mass;
  p_min_max_mass.help=((string)"Minimum maximum mass ")
    +"(in solar masses, default 2.0).";
  cl.par_list.insert(make_pair("min_max_mass",&p_min_max_mass));

  p_min_mass.d=&min_mass;
  p_min_mass.help=((string)"Minimum possible mass for any of individual ")+
    "neutron stars in solar masses. The default is 0.8 solar masses.";
  cl.par_list.insert(make_pair("min_mass",&p_min_mass));

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
  
  p_debug_line.b=&debug_line;
  p_debug_line.help=((string)"If true, output each line as its stored ")+
    "(default false).";
  cl.par_list.insert(make_pair("debug_line",&p_debug_line));
  
  p_debug_eos.b=&debug_eos;
  p_debug_eos.help=((string)"If true, output initial equation of state ")+
    "to file 'debug_eos.o2' and abort (default false).";
  cl.par_list.insert(make_pair("debug_eos",&p_debug_eos));

  p_baryon_density.b=&baryon_density;
  p_baryon_density.help=((string)"If true, compute baryon density ")+
    "and associated profiles (default true).";
  cl.par_list.insert(make_pair("baryon_density",&p_baryon_density));

  p_use_crust.b=&use_crust;
  p_use_crust.help=((string)"If true, use the default crust (default ")+
    "true).";
  cl.par_list.insert(make_pair("use_crust",&p_use_crust));

  p_inc_baryon_mass.b=&inc_baryon_mass;
  p_inc_baryon_mass.help=((string)"If true, compute the baryon mass ")+
    "(default false)";
  cl.par_list.insert(make_pair("inc_baryon_mass",&p_inc_baryon_mass));

  p_mvsr_pr_inc.d=&mvsr_pr_inc;
  p_mvsr_pr_inc.help=((string)"The multiplicative pressure increment for ")+
    "the TOV solver (default 1.1).";
  cl.par_list.insert(make_pair("mvsr_pr_inc",&p_mvsr_pr_inc));

  // --------------------------------------------------------

  p_nb_low.d=&nb_low;
  p_nb_low.help="Smallest baryon density grid point in 1/fm^3 (default 0.04).";
  cl.par_list.insert(make_pair("nb_low",&p_nb_low));

  p_nb_high.d=&nb_high;
  p_nb_high.help="Largest baryon density grid point in 1/fm^3 (default 1.24).";
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
  
  mcmc_class::setup_cli();
  
  return;
}

