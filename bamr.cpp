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

#ifdef NEVER_DEFINED

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

#endif

int bamr_class::set_model(std::vector<std::string> &sv, bool itive_com) {
  // We cannot use scr_out here because it isn't set until the call
  // to mcmc().
  if (sv.size()<2) {
    cerr << "Model name not given." << endl;
    return exc_efailed;
  }
  if (model_type==sv[1]) {
    cerr << "Model already set to " << sv[1] << endl;
    return 0;
  }
  mod->remove_params(cl);
  model_type=sv[1];
  if (sv[1]==((string)"twop")) {
    std::shared_ptr<model> mnew(new two_polytropes(set));
    mod=mnew;
  } else if (sv[1]==((string)"altp")) {
    std::shared_ptr<model> mnew(new alt_polytropes(set));
    mod=mnew;
  } else if (sv[1]==((string)"fixp")) {
    std::shared_ptr<model> mnew(new fixed_pressure(set));
    mod=mnew;
  } else if (sv[1]==((string)"qstar")) {
    std::shared_ptr<model> mnew(new quark_star(set));
    mod=mnew;
  } else if (sv[1]==((string)"genq")) {
    std::shared_ptr<model> mnew(new generic_quarks(set));
    mod=mnew;
  } else if (sv[1]==((string)"qmc")) {
    std::shared_ptr<model> mnew(new qmc_neut(set));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_threep")) {
    std::shared_ptr<model> mnew(new qmc_threep(set));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_fixp")) {
    std::shared_ptr<model> mnew(new qmc_fixp(set));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_twolines")) {
    std::shared_ptr<model> mnew(new qmc_twolines(set));
    mod=mnew;
  } else {
    cerr << "Model unknown." << endl;
    return exc_efailed;
  }
  mod->setup_params(cl);

  return 0;
}

void bamr_class::setup_cli() {

  // ---------------------------------------
  // Set options
    
  static const int nopt=4;
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
    /*    {'a',"add-data","Add data source to the list.",
     4,5,"<name> <file> <slice> <initial mass> [obj name]",
     ((string)"Specify data as a table3d object in a HDF5 file. ")+
     "The string <name> is the name used, <file> is the filename, "+
     "<slice> is the name of the slice in the table3d object, "+
     "<initial mass> is the initial mass for the first point, and "+
     "[obj name] is the optional name of table3d object in <file>. "+
     "If [obj name] is not specified, then the first table3d object "+
     "is used.",new comm_option_mfptr<bamr_class>(this,&bamr_class::add_data),
     cli::comm_option_both},
    */
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

  // --------------------------------------------------------
  
  mcmc_class::setup_cli();
  
  return;
}

