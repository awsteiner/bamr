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

void bamr_class::table_names_units(std::string &s, std::string &u) {

  model &m=*this->mod;
  
  mcmc_class::table_names_units(s,u);
  
  if (set.norm_max) {
    u+=". ";
  } else {
    u+=((std::string)"1/km^")+std::to_string(m.nsources)+"/Msun^"+
      std::to_string(m.nsources)+" ";
  }
  for(size_t i=0;i<m.nsources;i++) {
    s+=((std::string)"wgt_")+m.source_names[i]+" ";
    if (set.norm_max) {
      u+=". ";
    } else {
      u+="1/km/Msun ";
    }
  }
  
  // It is important here that all of these columns which store values
  // over a grid are either always positive or always negative,
  // because the code reports zero in the fill_line() function for
  // values beyond the end of the EOS or the M-R curve. 
  for(size_t i=0;i<m.nsources;i++) {
    s+=((std::string)"Rns_")+m.source_names[i]+" ";
    u+="km ";
  }
  for(size_t i=0;i<m.nsources;i++) {
    s+=((std::string)"Mns_")+m.source_names[i]+" ";
    u+="Msun ";
  }
  
  if (m.has_eos) {
    for(int i=0;i<set.grid_size;i++) {
      s+=((string)"P_")+std::to_string(i)+" ";
      u+="1/fm^4 ";
    }
  }
  
  for(int i=0;i<set.grid_size;i++) {
    s+=((string)"R_")+std::to_string(i)+" ";
    u+="km ";
    if (m.has_eos) {
      s+=((string)"PM_")+std::to_string(i)+" ";
      u+="1/fm^4 ";
    }
  }
  if (m.has_eos) {
    if (set.baryon_density) {
      for(int i=0;i<set.grid_size;i++) {
	s+=((string)"Pnb_")+std::to_string(i)+" ";
	u+="1/fm^4 ";
	s+=((string)"EoA_")+std::to_string(i)+" ";
	u+="MeV ";
      }
    }
    if (m.has_esym) {
      s+="S L ";
      u+="MeV MeV ";
    }
    s+="R_max M_max P_max e_max ";
    u+="km Msun 1/fm^4 1/fm^4 ";
    if (set.baryon_density) {
      s+="nb_max ";
      u+="1/fm^3 ";
    }
    for(size_t i=0;i<m.nsources;i++) {
      s+=((string)"ce_")+m.source_names[i]+" ";
      u+="1/fm^4 ";
    }
    if (set.baryon_density) {
      for(size_t i=0;i<m.nsources;i++) {
	s+=((string)"cnb_")+m.source_names[i]+" ";
	u+="1/fm^3 ";
      }
    }
  }

  return;
}

void bamr_class::fill_line(ubvector &pars, double weight, model_data &dat,
			   std::vector<double> &line) {

  mcmc_class::fill_line(pars,weight,dat,line);

  model &m=*this->mod;
  
  shared_ptr<table_units<> > tab_eos=dat.eos;
  shared_ptr<table_units<> > tab_mvsr=dat.mvsr;
  
  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (m.has_eos) {

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
    
    if (set.baryon_density) {
      // The central baryon density in the maximum mass configuration
      nbmax=tab_mvsr->get_constant("new_nb_max");
    }

  } else {
    // Need to set mmax for no EOS models to figure out how 
    // high up we should go for the radius grid 
    mmax=3.0;
  }

  for(size_t i=0;i<m.nsources;i++) {
    line.push_back(dat.wgts[i]);
  }
  for(size_t i=0;i<nparams;i++) {
    line.push_back(pars[i]);
  }
  for(size_t i=0;i<m.nsources;i++) {
    line.push_back(dat.rad[i]);
  }
  if (m.has_eos) {
    for(int i=0;i<set.grid_size;i++) {
      double eval=m.e_grid[i];
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
  for(int i=0;i<set.grid_size;i++) {
    double mval=m.m_grid[i];
    if (mval<mmax) {
      line.push_back(tab_mvsr->interp("gm",mval,"r"));
      if (m.has_eos) {
	line.push_back(tab_mvsr->interp("gm",mval,"pr"));
      }
    } else {
      line.push_back(0.0);
      if (m.has_eos) {
	line.push_back(0.0);
      }
    }
  }
  if (m.has_eos) {
    if (set.baryon_density) {
      for(int i=0;i<set.grid_size;i++) {
	double nbval=m.nb_grid[i];
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
    if (m.has_esym) {
      line.push_back(tab_eos->get_constant("S"));
      line.push_back(tab_eos->get_constant("L"));
    }
    line.push_back(rmax);
    line.push_back(mmax);
    line.push_back(pmax);
    line.push_back(emax);
    if (set.baryon_density) line.push_back(nbmax);
    for(size_t i=0;i<m.nsources;i++) {
      line.push_back(tab_mvsr->interp("gm",pars[nparams-m.nsources+i],"ed"));
    }
    if (set.baryon_density) {
      for(size_t i=0;i<m.nsources;i++) {
	line.push_back(tab_mvsr->interp("gm",pars[nparams-m.nsources+i],"nb"));
      }
    }
  }

  return;
}

void bamr_class::first_update(o2scl_hdf::hdf_file &hf) {

  mcmc_class::first_update(hf);

  model &m=*this->mod;
  
  hf.sets_vec("source_names",m.source_names);
  hf.sets_vec("source_fnames",m.source_fnames);
  hf.sets_vec("slice_names",m.slice_names);

  hf.set_szt("grid_size",set.grid_size);
  hf.set_szt("nsources",m.nsources);
  hf.sets("model",model_type);
  hf.setd("min_mass",set.min_mass);
  hf.setd("exit_mass",set.exit_mass);
  hf.setd("min_max_mass",set.min_max_mass);
  //hf.setd("input_dist_thresh",input_dist_thresh);
  hf.seti("use_crust",set.use_crust);
  hf.seti("baryon_density",set.baryon_density);
  hf.seti("debug_load",set.debug_load);
  hf.seti("debug_eos",set.debug_eos);
  hf.seti("debug_star",set.debug_star);
  hf.seti("inc_baryon_mass",set.inc_baryon_mass);
  hf.setd("nb_low",set.nb_low);
  hf.setd("nb_high",set.nb_high);
  hf.setd("e_low",set.e_low);
  hf.setd("e_high",set.e_high);
  hf.setd("m_low",set.m_low);
  hf.setd("m_high",set.m_high);

  hdf_output(hf,m.nb_grid,"nb_grid");
  hdf_output(hf,m.e_grid,"e_grid");
  hdf_output(hf,m.m_grid,"m_grid");
    
  return;
}

int bamr_class::mcmc_init() {

  // -----------------------------------------------------------
  // Make grids

  model &m=*this->mod;
  
  m.nb_grid=uniform_grid_end<double>(set.nb_low,set.nb_high,set.grid_size-1);
  m.e_grid=uniform_grid_end<double>(set.e_low,set.e_high,set.grid_size-1);
  m.m_grid=uniform_grid_end<double>(set.m_low,set.m_high,set.grid_size-1);

  return 0;
}

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
  };
  cl.set_comm_option_vec(nopt,options);

  // --------------------------------------------------------
  
  mcmc_class::setup_cli();
  
  return;
}

