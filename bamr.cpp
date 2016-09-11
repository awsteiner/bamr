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

int bamr_class::fill(const ubvector &pars, double weight, 
		      std::vector<double> &line, model_data &dat) {
  
  model &m=*this->mod;
  
  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (m.has_eos) {

    // The central energy density in the maximum mass configuration
    emax=dat.mvsr->max("ed");
    // The central pressure in the maximum mass configuration
    pmax=dat.mvsr->max("pr");
    // The maximum mass
    mmax=dat.mvsr->get_constant("m_max");
    // The radius of the maximum mass star
    rmax=dat.mvsr->get_constant("r_max");
    
    if (set.baryon_density) {

      // The highest baryon density in the EOS table
      nbmax2=dat.eos->max("nb");
      // The central baryon density in the maximum mass configuration
      nbmax=dat.mvsr->get_constant("nb_max");
      
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
      double val=dat.mvsr->interp
	("gm",pars[this->nparams-nsd.nsources+i],"ed");
      line.push_back(val);
    }
    if (set.baryon_density) {
      for(size_t i=0;i<nsd.nsources;i++) {
	line.push_back(dat.mvsr->interp
		       ("gm",pars[this->nparams-nsd.nsources+i],"nb"));
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
    line.push_back(dat.eos->get_constant("nt"));
    line.push_back(dat.eos->get_constant("prt"));
    for(int i=0;i<set.grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax) {
	double rval=dat.mvsr->interp("gm",mval,"r");
	line.push_back(rval-dat.mvsr->interp("gm",mval,"r0"));
      } else {
	line.push_back(0.0);
      }
    }
  }
  if (set.addl_quants) {
    double mmax=dat.mvsr->get_constant("m_max");
    for(int i=0;i<set.grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax) {
	double bm=dat.mvsr->interp("gm",mval,"bm");
	double rad=dat.mvsr->interp("gm",mval,"r");
	// rjw is km^4, so dividing by km/Msun gives Msun*km^2
	double I=dat.mvsr->interp("gm",mval,"rjw")/3.0/schwarz_km;
	line.push_back(bm);
	line.push_back((bm-mval)/mval);
	// Make unitless by dividing by G^2
	double I_bar=I/mval/mval/mval/schwarz_km/schwarz_km*4.0;
	line.push_back(I_bar);

#ifdef NEVER_DEFINED
	// Relation for lambda given I from Kent Yagi based
	// on Yagi and Yunes, Science, (2014)
	double a0=-210.327;
	double a1=481.472;
	double a2=-464.271;
	double a3=241.564;
	double a4=-70.6449;
	double a5=10.9775;
	double a6=-0.707066;

	double li=log(I_bar);
	double li2=li*li;
	double li3=li*li2;
	double li4=li*li3;
	double li5=li*li4;
	double li6=li*li5;
	  
	double lambda_bar=exp(a0+a1*li+a2*li2+a3*li3+a4*li4+a5*li5+a6*li6);
	double lambda=lambda_bar*pow(mval,5.0);
	cout << I_bar << " " << lambda_bar << endl;
	exit(-1);
#endif
	
      } else {
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
      }
    }
  }
  
  return o2scl::success;
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
  hf.setd_vec_copy("low",this->low_copy);
  hf.setd_vec_copy("high",this->high_copy);

  hdf_output(hf,m.nb_grid,"nb_grid");
  hdf_output(hf,m.e_grid,"e_grid");
  hdf_output(hf,m.m_grid,"m_grid");
    
  return;
}

int bamr_class::mcmc_init() {

  if (this->verbose>=2) {
    std::cout << "Start bamr_class::mcmc_init()." << std::endl;
  }
  
  model &m=*this->mod;
  
  mcmc_bamr::mcmc_init();

  if (this->file_opened==false) {
    // Open main output file
    this->scr_out.open((this->prefix+"_"+
			o2scl::itos(this->mpi_rank)+"_scr").c_str());
    this->scr_out.setf(std::ios::scientific);
    this->file_opened=true;
    this->scr_out << "Opened main file in function 'bamr_class::mcmc_init()'."
		  << std::endl;
  }
  
  // -----------------------------------------------------------
  // Make sure the settings are consistent

  // Does inc_baryon_mass also need baryon_density?
  if (set.inc_baryon_mass && !set.baryon_density) {
    scr_out << "Cannot use inc_baryon_mass=true with "
	    << "baryon_density=false." << endl;
    return exc_efailed;
  }
  if (set.compute_cthick && (!set.baryon_density || !set.use_crust)) {
    scr_out << "Cannot use compute_cthick=true with "
	    << "baryon_density=false or use_crust=false." << endl;
    return exc_efailed;
  }
  if (set.crust_from_L && (!m.has_esym || !set.use_crust ||
			   !set.baryon_density)) {
    scr_out << "crust_from_L: " << set.crust_from_L << std::endl;
    scr_out << "has_esym: " << m.has_esym << std::endl;
    scr_out << "use_crust: " << set.use_crust << std::endl;
    scr_out << "baryon_density: " << set.baryon_density << std::endl;
    scr_out << "Cannot use crust_from_L=true with a model which does not "
	    << "provide S and L\nor with use_crust=false or with "
	    << "baryon_density=false." << endl;
    return exc_efailed;
  }
  if (set.addl_quants && !set.inc_baryon_mass) {
    scr_out << "Cannot do additional quantities without including "
	    << "baryon mass." << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------
  // Add columns to table

  for(size_t i=0;i<nsd.nsources;i++) {
    this->tab->new_column(((std::string)"wgt_")+nsd.source_names[i]);
    if (!set.norm_max) {
      this->tab->set_unit(((std::string)"wgt_")+nsd.source_names[i],
			  "1/km/Msun");
    }
  }
  
  // It is important here that all of these columns which store values
  // over a grid are either always positive or always negative,
  // because the code reports zero in the fill_line() function for
  // values beyond the end of the EOS or the M-R curve. 
  for(size_t i=0;i<nsd.nsources;i++) {
    this->tab->new_column(((std::string)"Rns_")+nsd.source_names[i]);
    this->tab->set_unit(((std::string)"Rns_")+nsd.source_names[i],
			"km");
  }
  
  for(size_t i=0;i<nsd.nsources;i++) {
    this->tab->new_column(((std::string)"Mns_")+nsd.source_names[i]);
    this->tab->set_unit(((std::string)"Mns_")+nsd.source_names[i],
			"Msun");
  }
  
  if (m.has_eos) {
    for(int i=0;i<set.grid_size;i++) {
      this->tab->new_column(((string)"P_")+o2scl::itos(i));
      this->tab->set_unit(((string)"P_")+o2scl::itos(i),
			  "1/fm^4");
    }
  }
  
  for(int i=0;i<set.grid_size;i++) {
    this->tab->new_column(((string)"R_")+o2scl::itos(i));
    this->tab->set_unit(((string)"R_")+o2scl::itos(i),
			"km");
    if (m.has_eos) {
      this->tab->new_column(((string)"PM_")+o2scl::itos(i));
      this->tab->set_unit(((string)"PM_")+o2scl::itos(i),
			  "1/fm^4");
    }
  }
  if (m.has_eos) {
    if (set.baryon_density) {
      for(int i=0;i<set.grid_size;i++) {
	this->tab->new_column(((string)"Pnb_")+o2scl::itos(i));
	this->tab->set_unit(((string)"Pnb_")+o2scl::itos(i),
			    "1/fm^4");
	this->tab->new_column(((string)"EoA_")+o2scl::itos(i));
	this->tab->set_unit(((string)"EoA_")+o2scl::itos(i),
			    "MeV");
      }
    }
    if (m.has_esym) {
      this->tab->new_column("S");
      this->tab->set_unit("S","1/fm");
      this->tab->new_column("L");
      this->tab->set_unit("L","1/fm");
    }
    this->tab->new_column("R_max");
    this->tab->set_unit("R_max","km");
    this->tab->new_column("M_max");
    this->tab->set_unit("M_max","Msun");
    this->tab->new_column("P_max");
    this->tab->set_unit("P_max","1/fm^4");
    this->tab->new_column("e_max");
    this->tab->set_unit("e_max","1/fm^4");
    if (set.baryon_density) {
      this->tab->new_column("nb_max");
      this->tab->set_unit("nb_max","1/fm^3");
    }
    for(size_t i=0;i<nsd.nsources;i++) {
      this->tab->new_column(((string)"ce_")+nsd.source_names[i]);
      this->tab->set_unit(((string)"ce_")+nsd.source_names[i],
			  "1/fm^4");
    }
    if (set.baryon_density) {
      for(size_t i=0;i<nsd.nsources;i++) {
	this->tab->new_column(((string)"cnb_")+nsd.source_names[i]);
	this->tab->set_unit(((string)"cnb_")+nsd.source_names[i],
			      "1/fm^3");
      }
      this->tab->new_column("gm_nb1");
      this->tab->set_unit("gm_nb1","Msun");
      this->tab->new_column("r_nb1");
      this->tab->set_unit("r_nb1","km");
      this->tab->new_column("gm_nb2");
      this->tab->set_unit("gm_nb2","Msun");
      this->tab->new_column("r_nb2");
      this->tab->set_unit("r_nb2","km");
      this->tab->new_column("gm_nb3");
      this->tab->set_unit("gm_nb3","Msun");
      this->tab->new_column("r_nb3");
      this->tab->set_unit("r_nb3","km");
      this->tab->new_column("gm_nb4");
      this->tab->set_unit("gm_nb4","Msun");
      this->tab->new_column("r_nb4");
      this->tab->set_unit("r_nb4","km");
      this->tab->new_column("gm_nb5");
      this->tab->set_unit("gm_nb5","Msun");
      this->tab->new_column("r_nb5");
      this->tab->set_unit("r_nb5","km");
    }
    if (set.compute_cthick) {
      this->tab->new_column("nt");
      this->tab->set_unit("nt","1/fm^3");
      this->tab->new_column("prt");
      this->tab->set_unit("prt","1/fm^4");
      for(int i=0;i<set.grid_size;i++) {
        this->tab->new_column(((string)"ct_")+o2scl::itos(i));
        this->tab->set_unit(((string)"ct_")+o2scl::itos(i),"km");
      }
    }
  }
  if (set.addl_quants) {
    for(int i=0;i<set.grid_size;i++) {
      this->tab->new_column(((string)"Mb_")+o2scl::itos(i));
      this->tab->set_unit(((string)"Mb_")+o2scl::itos(i),"Msun");
      this->tab->new_column(((string)"be_")+o2scl::itos(i));
      this->tab->set_unit(((string)"be_")+o2scl::itos(i),"Msun");
      this->tab->new_column(((string)"I_")+o2scl::itos(i));
      this->tab->set_unit(((string)"I_")+o2scl::itos(i),
			    "Msun*km^2");
      //this->tab->new_column(((string)"lambda_")+o2scl::itos(i));
    }
  }

  // -----------------------------------------------------------
  // Make grids

  m.nb_grid=uniform_grid_end<double>(set.nb_low,set.nb_high,set.grid_size-1);
  m.e_grid=uniform_grid_end<double>(set.e_low,set.e_high,set.grid_size-1);
  m.m_grid=uniform_grid_end<double>(set.m_low,set.m_high,set.grid_size-1);

  // -----------------------------------------------------------
  // Load data

  nsd.load_mc(this->scr_out,mpi_nprocs,mpi_rank,set);

  // -----------------------------------------------------------
  // Prepare data objects

  for(size_t i=0;i<data_arr.size();i++) {
    data_arr[i].rad.resize(nsd.nsources);
    data_arr[i].mass.resize(nsd.nsources);
    data_arr[i].wgts.resize(nsd.nsources);
  }

  if (this->verbose>=2) {
    std::cout << "End bamr_class::mcmc_init()." << std::endl;
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

int bamr_class::compute_point(const ubvector &pars, std::ofstream &scr_out, 
			      double &weight, model_data &dat) {
  return mod->compute_point(pars,scr_out,weight,dat);
}

int bamr_class::mcmc_func(std::vector<std::string> &sv, bool itive_com) {

  if (model_type.length()==0) {
    cerr << "Model not set in 'mcmc' command." << endl;
    return 1;
  }
  
  std::vector<std::string> names;
  std::vector<std::string> units;

  ubvector low;
  ubvector high;
  // Get names and units for parameters from model (which also
  // automatically includes nuisance variables for the data points)
  // The other columns and units are specified in mcmc_init()
  // function manually using a call to table::new_column().
  mod->get_param_info(names,units,low,high);
  set_names_units(names,units);

  ubvector init(names.size());
  mod->initial_point(init);
  
  bamr::point_funct mf=std::bind
    (std::mem_fn<int(const ubvector &,ofstream &,double &,model_data &)>
     (&bamr_class::compute_point),this,
     std::placeholders::_2,std::ref(scr_out),std::placeholders::_3,
     std::placeholders::_4);
  bamr::fill_funct mt=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     model_data &)>
     (&bamr_class::fill),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);

  set.verbose=this->verbose;
  this->mcmc(init.size(),init,low,high,mf,mt);
  
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

