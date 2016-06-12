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
  
  mcmc_bamr::table_names_units(s,u);

  /*
    if (set.norm_max) {
    u+=". ";
    } else {
    u+=((std::string)"1/km^")+std::to_string(nsd.nsources)+"/Msun^"+
    std::to_string(nsd.nsources)+" ";
    }
  */
  
  for(size_t i=0;i<nsd.nsources;i++) {
    s+=((std::string)"wgt_")+nsd.source_names[i]+" ";
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
  for(size_t i=0;i<nsd.nsources;i++) {
    s+=((std::string)"Rns_")+nsd.source_names[i]+" ";
    u+="km ";
  }

  for(size_t i=0;i<nsd.nsources;i++) {
    s+=((std::string)"Mns_")+nsd.source_names[i]+" ";
    u+="km ";
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
    for(size_t i=0;i<nsd.nsources;i++) {
      s+=((string)"ce_")+nsd.source_names[i]+" ";
      u+="1/fm^4 ";
    }
    if (set.baryon_density) {
      for(size_t i=0;i<nsd.nsources;i++) {
	s+=((string)"cnb_")+nsd.source_names[i]+" ";
	u+="1/fm^3 ";
      }
      s+="gm_nb1 r_nb1 ";
      u+="Msun km ";
      s+="gm_nb2 r_nb2 ";
      u+="Msun km ";
      s+="gm_nb3 r_nb3 ";
      u+="Msun km ";
      s+="gm_nb4 r_nb4 ";
      u+="Msun km ";
      s+="gm_nb5 r_nb5 ";
      u+="Msun km ";
    }
    if (set.compute_cthick) {
      if (set.nt_corr) {
        s+="nt prt ";
        u+="1/fm^3 1/fm^4 ";
      }
      for(int i=0;i<set.grid_size;i++) {
        s+=((string)"ct06_")+std::to_string(i)+" ";
        u+="km ";
        s+=((string)"ct08_")+std::to_string(i)+" ";
        u+="km ";
        s+=((string)"ct10_")+std::to_string(i)+" ";
        u+="km ";
      }
    }
  }
  if (set.addl_quants) {
    for(int i=0;i<set.grid_size;i++) {
      s+=((string)"Mb_")+std::to_string(i)+" ";
      u+="Msun ";
      s+=((string)"be_")+std::to_string(i)+" ";
      u+="Msun ";
      s+=((string)"I_")+std::to_string(i)+" ";
      u+="Msun*km^2 ";
      s+=((string)"lambda_")+std::to_string(i)+" ";
      u+=". ";
    }
  }

  return;
}

void bamr_class::fill_line(ubvector &pars, double weight, 
			   std::vector<double> &line, model_data &dat) {
  
  std::cout << "bfl1" << std::endl;
  
  mcmc_bamr::fill_line(pars,weight,line,dat);
  
  std::cout << "bfl2" << std::endl;
  
  model &m=*this->mod;
  
  size_t nparams=this->param_names.size();
  
  std::cout << "bfl3" << std::endl;

  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (m.has_eos) {

    // The highest baryon density in the EOS table
    nbmax2=dat.eos->max("nb");
    // The central energy density in the maximum mass configuration
    emax=dat.mvsr->max("ed");
    // The central pressure in the maximum mass configuration
    pmax=dat.mvsr->max("pr");
    // The maximum mass
    mmax=dat.mvsr->get_constant("new_max");
    // The radius of the maximum mass star
    rmax=dat.mvsr->get_constant("new_r_max");
    
    if (set.baryon_density) {
      // The central baryon density in the maximum mass configuration
      nbmax=dat.mvsr->get_constant("new_nb_max");
    }

  } else {
    // Need to set mmax for no EOS models to figure out how 
    // high up we should go for the radius grid 
    mmax=3.0;
  }

  for(size_t i=0;i<nsd.nsources;i++) {
    line.push_back(dat.wgts[i]);
  }
  for(size_t i=0;i<nsd.nsources;i++) {
    line.push_back(dat.rad[i]);
  }
  for(size_t i=0;i<nsd.nsources;i++) {
    line.push_back(dat.mass[i]);
  }
  
  if (m.has_eos) {
    for(int i=0;i<set.grid_size;i++) {
      double eval=m.e_grid[i];
      // Make sure the energy density from the grid
      // isn't beyond the table limit
      double emax2=dat.eos->max("ed");
      if (eval<emax2) {
	double pres_temp=dat.eos->interp("ed",eval,"pr");
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
      line.push_back(dat.mvsr->interp("gm",mval,"r"));
      if (m.has_eos) {
	line.push_back(dat.mvsr->interp("gm",mval,"pr"));
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
	  double pres_temp=dat.eos->interp("nb",nbval,"pr");
	  if (pres_temp<pmax) {
	    line.push_back(pres_temp);
	  } else {
	    line.push_back(0.0);
	  }
	  double eval2=dat.eos->interp("nb",nbval,"ed");
	  double eoa_val2=eval2/nbval-939.0/o2scl_const::hc_mev_fm;
	  line.push_back(eoa_val2);
	} else {
	  line.push_back(0.0);
	  line.push_back(0.0);
	}
      }
    }
    if (m.has_esym) {
      line.push_back(dat.eos->get_constant("S"));
      line.push_back(dat.eos->get_constant("L"));
    }
    line.push_back(rmax);
    line.push_back(mmax);
    line.push_back(pmax);
    line.push_back(emax);
    if (set.baryon_density) line.push_back(nbmax);
    for(size_t i=0;i<nsd.nsources;i++) {
      line.push_back(dat.mvsr->interp("gm",pars[nparams-nsd.nsources+i],"ed"));
    }
    if (set.baryon_density) {
      for(size_t i=0;i<nsd.nsources;i++) {
	line.push_back(dat.mvsr->interp
		       ("gm",pars[nparams-nsd.nsources+i],"nb"));
      }
    }
  }
  
  if (set.baryon_density) {
    line.push_back(dat.mvsr->get_constant("gm_nb1"));
    line.push_back(dat.mvsr->get_constant("r_nb1"));
    line.push_back(dat.mvsr->get_constant("gm_nb2"));
    line.push_back(dat.mvsr->get_constant("r_nb2"));
    line.push_back(dat.mvsr->get_constant("gm_nb3"));
    line.push_back(dat.mvsr->get_constant("r_nb3"));
    line.push_back(dat.mvsr->get_constant("gm_nb4"));
    line.push_back(dat.mvsr->get_constant("r_nb4"));
    line.push_back(dat.mvsr->get_constant("gm_nb5"));
    line.push_back(dat.mvsr->get_constant("r_nb5"));
  }
  if (set.compute_cthick) {
    if (set.nt_corr) {
      line.push_back(dat.eos->get_constant("nt"));
      line.push_back(dat.eos->get_constant("prt"));
    }
    for(int i=0;i<set.grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax && dat.mvsr->is_column("r2")) {
	double rval=dat.mvsr->interp("gm",mval,"r");
	line.push_back(rval-dat.mvsr->interp("gm",mval,"r0"));
	line.push_back(rval-dat.mvsr->interp("gm",mval,"r1"));
	line.push_back(rval-dat.mvsr->interp("gm",mval,"r2"));
      } else {
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
      }
    }
  }
  /*
  if (epja_mode) {
    double mmax=dat.mvsr->get_constant("new_max");
    for(int i=0;i<grid_size;i++) {
      double mval=m_grid[i];
      if (mval<mmax) {
	double bm=dat.mvsr->interp("gm",mval,"bm");
	double rad=dat.mvsr->interp("gm",mval,"r");
	// rjw is km^4, so dividing by km/Msun gives Msun*km^2
	double I=dat.mvsr->interp("gm",mval,"rjw")/3.0/schwarz_km;
	line.push_back(bm);
	line.push_back((bm-mval)/mval);
	// Make unitless by dividing by G^2
	line.push_back(I/mval/mval/mval/schwarz_km/schwarz_km*4.0);
      } else {
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
      }
    }
    }
  */
  return;
}

void bamr_class::first_update(o2scl_hdf::hdf_file &hf) {

  mcmc_bamr::first_update(hf);

  model &m=*this->mod;
  
  hf.sets_vec("source_names",nsd.source_names);
  hf.sets_vec("source_fnames",nsd.source_fnames);
  hf.sets_vec("slice_names",nsd.slice_names);

  hf.set_szt("grid_size",set.grid_size);
  hf.set_szt("nsources",nsd.nsources);
  hf.sets("model",model_type);
  hf.setd("min_mass",set.min_mass);
  hf.setd("exit_mass",set.exit_mass);
  hf.setd("min_max_mass",set.min_max_mass);
  hf.setd("input_dist_thresh",set.input_dist_thresh);
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

  std::cout << "Hx." << std::endl;
  mcmc_bamr::mcmc_init();
  std::cout << "Hy." << std::endl;
  
  model &m=*this->mod;
  
  if (set.compute_cthick && !set.baryon_density) {
    scr_out << "Cannot use 'compute_cthick=true' with "
	    << "'baryon_density=false'." << endl;
    return exc_efailed;
  }
  if (set.nt_corr && (!set.compute_cthick || !set.baryon_density)) {
    scr_out << "Cannot use 'nt_corr=true' with "
	    << "'compute_cthick=false' or "
	    << "'baryon_density=false'." << endl;
    return exc_efailed;
  }
  if (set.nt_corr && !m.has_esym) {
    scr_out << "Cannot use 'nt_corr=true' with a model which does not "
	    << "provide S and L." << endl;
    return exc_efailed;
  }
  if (set.crust_from_L && (!m.has_esym || !set.use_crust ||
			   !set.baryon_density || !set.compute_cthick)) {
    scr_out << "Cannot use 'crust_from_L=true' with a model which does not "
	    << "provide S and L or with 'use_crust=false' or with "
	    << "'baryon_density=false or with 'compute_cthick=false'." 
	    << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------
  // Make grids

  m.nb_grid=uniform_grid_end<double>(set.nb_low,set.nb_high,set.grid_size-1);
  m.e_grid=uniform_grid_end<double>(set.e_low,set.e_high,set.grid_size-1);
  m.m_grid=uniform_grid_end<double>(set.m_low,set.m_high,set.grid_size-1);

  // -----------------------------------------------------------
  // Load data

  m.mpi_rank=mpi_rank;
  m.mpi_nprocs=mpi_nprocs;
  m.load_mc(this->scr_out);

  // -----------------------------------------------------------
  // Prepare data objects

  std::cout << "H0 " << data_arr.size() << " " << nsd.nsources << std::endl;
  for(size_t i=0;i<data_arr.size();i++) {
    data_arr[i].rad.resize(nsd.nsources);
    data_arr[i].mass.resize(nsd.nsources);
    data_arr[i].wgts.resize(nsd.nsources);
  }
  std::cout << "H1." << std::endl;

  // -----------------------------------------------------------
  // Prepare crust

  if (set.use_crust) {
    m.teos.default_low_dens_eos();
    
    // Get the transition density from the crust
    double pt, pw;
    m.teos.get_transition(pt,pw);
    // We set the transition density a bit lower (because by default
    // it's the largest pressure in the crust EOS) and then add a 
    // small width
    m.teos.transition_mode=eos_tov_interp::smooth_trans;
    m.teos.set_transition(pt/1.2,1.2);
    
  } else {
    m.teos.no_low_dens_eos();
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
  if (model_type==sv[1]) {
    cerr << "Model already set to " << sv[1] << endl;
    return 0;
  }
  if (model_type.length()>0) {
    mod->remove_params(cl);
  }
  if (sv[1]==((string)"twop")) {
    std::shared_ptr<model> mnew(new two_polytropes(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"altp")) {
    std::shared_ptr<model> mnew(new alt_polytropes(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"fixp")) {
    std::shared_ptr<model> mnew(new fixed_pressure(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"qstar")) {
    std::shared_ptr<model> mnew(new quark_star(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"genq")) {
    std::shared_ptr<model> mnew(new generic_quarks(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"qmc")) {
    std::shared_ptr<model> mnew(new qmc_neut(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_threep")) {
    std::shared_ptr<model> mnew(new qmc_threep(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_fixp")) {
    std::shared_ptr<model> mnew(new qmc_fixp(set,nsd));
    mod=mnew;
  } else if (sv[1]==((string)"qmc_twolines")) {
    std::shared_ptr<model> mnew(new qmc_twolines(set,nsd));
    mod=mnew;
  } else {
    cerr << "Model unknown." << endl;
    return exc_efailed;
  }
  model_type=sv[1];
  mod->setup_params(cl);

  return 0;
}

int bamr_class::mcmc_func(std::vector<std::string> &sv, bool itive_com) {

  std::cout << "H7." << std::endl;

  std::vector<std::string> names;
  std::vector<std::string> units;

  ubvector low;
  ubvector high;
  mod->get_param_info(names,units,low,high);

  size_t nparams=names.size();
  std::cout << "H7b " << nparams << std::endl;
  ubvector init(nparams);
  mod->initial_point(init);
  
  int success;
  bamr::point_funct mf=std::bind
    (std::mem_fn<double(const ubvector &,ofstream &,int &,model_data &)>
     (&model::compute_point),mod,
     std::placeholders::_2,std::ref(scr_out),std::ref(success),
     std::placeholders::_3);
  bamr::measure_funct mt=std::bind
    (std::mem_fn<int(const ubvector &,double,size_t,bool,model_data &)>
     (&mcmc_bamr::add_line),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
     std::placeholders::_5);
  
  std::cout << "H8." << std::endl;

  this->mcmc(names.size(),init,low,high,mf,mt);
  
  return 0;
}

int bamr_class::add_data(std::vector<std::string> &sv, bool itive_com) {
  nsd.add_data(sv,itive_com);
  return 0;
}

void bamr_class::setup_cli() {
  
  mcmc_bamr::setup_cli();

  set.setup_cli(cl);
  
  // ---------------------------------------
  // Set options
    
  static const int nopt=3;
  comm_option_s options[nopt]={
    {'m',"mcmc","Perform the Markov Chain Monte Carlo simulation.",
     0,0,"",((std::string)"This is the main part of ")+
     "the code which performs the simulation. Make sure to set the "+
     "model first using the 'model' command first.",
     new o2scl::comm_option_mfptr<bamr_class>(this,&bamr_class::mcmc_func),
     o2scl::cli::comm_option_both},
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
  };
  cl.set_comm_option_vec(nopt,options);

  // --------------------------------------------------------
  
  return;
}

