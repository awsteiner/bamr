/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2018, Andrew W. Steiner
  
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
#include "mcmc_bamr.h"

#include <o2scl/vector.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

mcmc_bamr::mcmc_bamr(size_t n_omp_threads) {
  model_type="";
  this->n_threads=n_omp_threads;
  bc_arr.resize(this->n_threads);
  set=std::make_shared<settings>();
  nsd=std::make_shared<ns_data>();
  for(size_t i=0;i<this->n_threads;i++) {
    bc_arr[i]=new bamr_class;
    bc_arr[i]->set=set;
    bc_arr[i]->nsd=nsd;
  }
}

void mcmc_bamr::file_header(o2scl_hdf::hdf_file &hf) {

  mcmc_para_cli::file_header(hf);
  
  model &m=*(bc_arr[0]->mod);
  
  hf.sets_vec("source_names",nsd->source_names);
  hf.sets_vec("source_fnames",nsd->source_fnames);
  hf.sets_vec("slice_names",nsd->slice_names);

  hf.set_szt("grid_size",set->grid_size);
  hf.set_szt("n_sources",nsd->n_sources);
  hf.sets("model",model_type);
  hf.setd("min_mass",set->min_mass);
  hf.setd("exit_mass",set->exit_mass);
  hf.setd("min_max_mass",set->min_max_mass);
  hf.setd("input_dist_thresh",set->input_dist_thresh);
  hf.seti("use_crust",set->use_crust);
  hf.seti("baryon_density",set->baryon_density);
  hf.seti("debug_load",set->debug_load);
  hf.seti("debug_eos",set->debug_eos);
  hf.seti("debug_star",set->debug_star);
  hf.seti("inc_baryon_mass",set->inc_baryon_mass);
  hf.seti("addl_quants",set->addl_quants);
  hf.setd("nb_low",set->nb_low);
  hf.setd("nb_high",set->nb_high);
  hf.setd("e_low",set->e_low);
  hf.setd("e_high",set->e_high);
  hf.setd("m_low",set->m_low);
  hf.setd("m_high",set->m_high);

  hdf_output(hf,m.nb_grid,"nb_grid");
  hdf_output(hf,m.e_grid,"e_grid");
  hdf_output(hf,m.m_grid,"m_grid");

  return;
}

int mcmc_bamr::mcmc_init() {

  if (this->verbose>=2) {
    std::cout << "(rank " << this->mpi_rank
	      << ") Start mcmc_bamr::mcmc_init()." << std::endl;
  }
  
  model &m=*(bc_arr[0]->mod);
  
  // This ensures enough space for all the
  // default return values in models.h
  this->ret_value_counts.resize(21);

  // Copy parameter values to all of the model objects
  for(size_t i=1;i<bc_arr.size();i++) {
    model &m2=*(bc_arr[i]->mod);
    m.copy_params(m2);
  }
  
  mcmc_para_cli::mcmc_init();

  // -----------------------------------------------------------
  // Make sure the settings are consistent

  // Does inc_baryon_mass also need baryon_density?
  if (set->inc_baryon_mass && !set->baryon_density) {
    scr_out << "Cannot use inc_baryon_mass=true with "
	    << "baryon_density=false." << endl;
    return exc_efailed;
  }
  if (set->compute_cthick && (!set->baryon_density || !set->use_crust)) {
    scr_out << "Cannot use compute_cthick=true with "
	    << "baryon_density=false or use_crust=false." << endl;
    return exc_efailed;
  }
  if (set->crust_from_L && (!m.has_esym || !set->use_crust ||
			   !set->baryon_density)) {
    scr_out << "crust_from_L: " << set->crust_from_L << std::endl;
    scr_out << "has_esym: " << m.has_esym << std::endl;
    scr_out << "use_crust: " << set->use_crust << std::endl;
    scr_out << "baryon_density: " << set->baryon_density << std::endl;
    scr_out << "Cannot use crust_from_L=true with a model which does not "
	    << "provide S and L\nor with use_crust=false or with "
	    << "baryon_density=false." << endl;
    return exc_efailed;
  }
  if (set->addl_quants && !set->inc_baryon_mass) {
    scr_out << "Cannot do additional quantities without including "
	    << "baryon mass." << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------
  // Add columns to table

  for(size_t i=0;i<nsd->n_sources;i++) {
    this->table->new_column(((std::string)"wgt_")+nsd->source_names[i]);
    if (!set->norm_max) {
      this->table->set_unit(((std::string)"wgt_")+nsd->source_names[i],
			  "1/km/Msun");
    }
  }
  
  // It is important here that all of these columns which store values
  // over a grid are either always positive or always negative,
  // because the code reports zero in the fill_line() function for
  // values beyond the end of the EOS or the M-R curve. 
  for(size_t i=0;i<nsd->n_sources;i++) {
    this->table->new_column(((std::string)"Rns_")+nsd->source_names[i]);
    this->table->set_unit(((std::string)"Rns_")+nsd->source_names[i],
			"km");
  }
  
  for(size_t i=0;i<nsd->n_sources;i++) {
    this->table->new_column(((std::string)"Mns_")+nsd->source_names[i]);
    this->table->set_unit(((std::string)"Mns_")+nsd->source_names[i],
			"Msun");
  }
  
  if (m.has_eos) {
    for(int i=0;i<set->grid_size;i++) {
      this->table->new_column(((string)"P_")+o2scl::itos(i));
      this->table->set_unit(((string)"P_")+o2scl::itos(i),
			  "1/fm^4");
    }
  }
  
  for(int i=0;i<set->grid_size;i++) {
    this->table->new_column(((string)"R_")+o2scl::itos(i));
    this->table->set_unit(((string)"R_")+o2scl::itos(i),
			"km");
    if (m.has_eos) {
      this->table->new_column(((string)"PM_")+o2scl::itos(i));
      this->table->set_unit(((string)"PM_")+o2scl::itos(i),
			  "1/fm^4");
    }
  }
  if (m.has_eos) {
    if (set->baryon_density) {
      for(int i=0;i<set->grid_size;i++) {
	this->table->new_column(((string)"Pnb_")+o2scl::itos(i));
	this->table->set_unit(((string)"Pnb_")+o2scl::itos(i),
			    "1/fm^4");
	this->table->new_column(((string)"EoA_")+o2scl::itos(i));
	this->table->set_unit(((string)"EoA_")+o2scl::itos(i),
			    "MeV");
      }
    }
    if (m.has_esym) {
      this->table->new_column("S");
      this->table->set_unit("S","1/fm");
      this->table->new_column("L");
      this->table->set_unit("L","1/fm");
    }
    this->table->new_column("R_max");
    this->table->set_unit("R_max","km");
    this->table->new_column("M_max");
    this->table->set_unit("M_max","Msun");
    this->table->new_column("P_max");
    this->table->set_unit("P_max","1/fm^4");
    this->table->new_column("e_max");
    this->table->set_unit("e_max","1/fm^4");
    if (set->baryon_density) {
      this->table->new_column("nb_max");
      this->table->set_unit("nb_max","1/fm^3");
    }
    for(size_t i=0;i<nsd->n_sources;i++) {
      this->table->new_column(((string)"ce_")+nsd->source_names[i]);
      this->table->set_unit(((string)"ce_")+nsd->source_names[i],
			  "1/fm^4");
    }
    if (set->baryon_density) {
      for(size_t i=0;i<nsd->n_sources;i++) {
	this->table->new_column(((string)"cnb_")+nsd->source_names[i]);
	this->table->set_unit(((string)"cnb_")+nsd->source_names[i],
			      "1/fm^3");
      }
      this->table->new_column("gm_nb1");
      this->table->set_unit("gm_nb1","Msun");
      this->table->new_column("r_nb1");
      this->table->set_unit("r_nb1","km");
      this->table->new_column("gm_nb2");
      this->table->set_unit("gm_nb2","Msun");
      this->table->new_column("r_nb2");
      this->table->set_unit("r_nb2","km");
      this->table->new_column("gm_nb3");
      this->table->set_unit("gm_nb3","Msun");
      this->table->new_column("r_nb3");
      this->table->set_unit("r_nb3","km");
      this->table->new_column("gm_nb4");
      this->table->set_unit("gm_nb4","Msun");
      this->table->new_column("r_nb4");
      this->table->set_unit("r_nb4","km");
      this->table->new_column("gm_nb5");
      this->table->set_unit("gm_nb5","Msun");
      this->table->new_column("r_nb5");
      this->table->set_unit("r_nb5","km");
    }
    if (set->compute_cthick) {
      this->table->new_column("nt");
      this->table->set_unit("nt","1/fm^3");
      this->table->new_column("Pt");
      this->table->set_unit("Pt","1/fm^4");
      for(int i=0;i<set->grid_size;i++) {
        this->table->new_column(((string)"CT_")+o2scl::itos(i));
        this->table->set_unit(((string)"CT_")+o2scl::itos(i),"km");
      }
    }
  }
  if (set->addl_quants) {
    for(int i=0;i<set->grid_size;i++) {
      this->table->new_column(((string)"MB_")+o2scl::itos(i));
      this->table->set_unit(((string)"MB_")+o2scl::itos(i),"Msun");
      this->table->new_column(((string)"BE_")+o2scl::itos(i));
      this->table->set_unit(((string)"BE_")+o2scl::itos(i),"Msun");
      this->table->new_column(((string)"I_")+o2scl::itos(i));
      this->table->set_unit(((string)"I_")+o2scl::itos(i),
			    "Msun*km^2");
      //this->table->new_column(((string)"lambda_")+o2scl::itos(i));
    }
  }
  if (nsd->source_fnames_alt.size()>0) {
    for(size_t i=0;i<nsd->n_sources;i++) {
      this->table->new_column(((std::string)"alt_")+o2scl::szttos(i));
    }
  }

  // -----------------------------------------------------------
  // Make grids

  for(size_t i=0;i<n_threads;i++) {
    bc_arr[i]->mod->nb_grid=uniform_grid_end<double>
      (set->nb_low,set->nb_high,set->grid_size-1);
    bc_arr[i]->mod->e_grid=uniform_grid_end<double>
      (set->e_low,set->e_high,set->grid_size-1);
    bc_arr[i]->mod->m_grid=uniform_grid_end<double>
      (set->m_low,set->m_high,set->grid_size-1);
  }

  // -----------------------------------------------------------
  // Load data

  nsd->load_mc(this->scr_out,mpi_size,mpi_rank,set);

  // -----------------------------------------------------------
  // Prepare data objects

  for(size_t i=0;i<data_arr.size();i++) {
    data_arr[i].rad.resize(nsd->n_sources);
    data_arr[i].mass.resize(nsd->n_sources);
    data_arr[i].wgts.resize(nsd->n_sources);
  }

  if (this->verbose>=2) {
    std::cout << "(rank " << this->mpi_rank
	      << ") End mcmc_bamr::mcmc_init()." << std::endl;
  }

  return 0;
}

int mcmc_bamr::set_model(std::vector<std::string> &sv, bool itive_com) {
  
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
    bc_arr[0]->mod->remove_params(cl);
  }
  if (sv[1]==((string)"twop")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new two_polytropes(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"altp")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new alt_polytropes(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"fixp")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new fixed_pressure(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"qstar")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new quark_star(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"genq")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new generic_quarks(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"qmc")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new qmc_neut(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"qmc_threep")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new qmc_threep(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"qmc_fixp")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new qmc_fixp(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"qmc_twolines")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new qmc_twolines(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else {
    cerr << "Model unknown." << endl;
    return exc_efailed;
  }
  model_type=sv[1];
  bc_arr[0]->mod->setup_params(cl);
  return 0;
}

int mcmc_bamr::mcmc_func(std::vector<std::string> &sv, bool itive_com) {

  if (model_type.length()==0) {
    cerr << "Model not set in 'mcmc' command." << endl;
    return 1;
  }
  
  std::vector<std::string> names;
  std::vector<std::string> units;

  ubvector low;
  ubvector high;
  // Get upper and lower parameter limits and also the column names
  // and units for the data table (which also automatically includes
  // nuisance variables for the data points). The other columns and
  // units are specified in mcmc_init() function manually using a call
  // to table::new_column().
  bc_arr[0]->mod->get_param_info(names,units,low,high);
  set_names_units(names,units);

  // Get the parameter initial values for this model 
  ubvector init(names.size());
  bc_arr[0]->mod->initial_point(init);
  this->initial_points.clear();
  for(size_t i=0;i<n_threads;i++) {
    this->initial_points.push_back(init);
  }

  vector<bamr::point_funct> pfa(n_threads);
  vector<bamr::fill_funct> ffa(n_threads);
  for(size_t i=0;i<n_threads;i++) {
    pfa[i]=std::bind
      (std::mem_fn<int(const ubvector &,ofstream &,double &,model_data &)>
       (&bamr_class::compute_point),bc_arr[i],std::placeholders::_2,
       std::ref(scr_out),std::placeholders::_3,std::placeholders::_4);
    ffa[i]=std::bind
      (std::mem_fn<int(const ubvector &,double,vector<double> &,
		       model_data &)>
       (&bamr_class::fill),bc_arr[i],std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
  }
  
  this->mcmc(init.size(),low,high,pfa,ffa);
  
  return 0;
}

int mcmc_bamr::add_data(std::vector<std::string> &sv, bool itive_com) {
  nsd->add_data(sv,itive_com);
  return 0;
}

int mcmc_bamr::add_data_alt(std::vector<std::string> &sv, bool itive_com) {
  nsd->add_data_alt(sv,itive_com);
  return 0;
}

void mcmc_bamr::setup_cli() {
  
  mcmc_para_cli::setup_cli(cl);

  set->setup_cli(cl);
  
  // ---------------------------------------
  // Set options
    
  static const int nopt=4;
  comm_option_s options[nopt]={
    {'m',"mcmc","Perform the Markov Chain Monte Carlo simulation.",
     0,0,"",((std::string)"This is the main part of ")+
     "the code which performs the simulation. Make sure to set the "+
     "model first using the 'model' command first.",
     new o2scl::comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::mcmc_func),
     o2scl::cli::comm_option_both},
    {'o',"model","Choose model.",
     1,1,"<model name>",((string)"Choose the EOS parameterization model. ")+
     "Possible values are 'twop', 'altp', 'fixp', 'genq', 'qstar', "+
     "'qmc', 'qmc_threep' ,'qmc_fixp', and 'qmc_twolines'. A "+
     "model must be chosen before a MCMC run.",
     new comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::set_model),
     cli::comm_option_both},
    {'a',"add-data","Add data source to the list.",
     4,5,"<name> <file> <slice> <initial mass> [obj name]",
     ((string)"Specify data as a table3d object in a HDF5 file. ")+
     "The string <name> is the name used, <file> is the filename, "+
     "<slice> is the name of the slice in the table3d object, "+
     "<initial mass> is the initial mass for the first point, and "+
     "[obj name] is the optional name of table3d object in <file>. "+
     "If [obj name] is not specified, then the first table3d object "+
     "is used.",new comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::add_data),
     cli::comm_option_both},
    {0,"add-data-alt","Add data source to the list.",
     5,6,"<name> <file> <alt file> <slice> <initial mass> [obj name]",
     ((string)"Specify data as a table3d object in two HDF5 files. ")+
     "The string <name> is the name used, <file> and <alt file> are "+
     "the filenames, "+
     "<slice> is the name of the slice in the table3d object, "+
     "<initial mass> is the initial mass for the first point, and "+
     "[obj name] is the optional name of table3d object in <file>. "+
     "If [obj name] is not specified, then the first table3d object "+
     "is used.",new comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::add_data_alt),
     cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);

  // --------------------------------------------------------
  
  return;
}

