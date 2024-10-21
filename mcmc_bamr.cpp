/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2024, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
  and Andrew W. Steiner
  
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
#include <o2scl/interpm_idw.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/interpm_python.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf; // For I/O with HDF files
using namespace o2scl_const; // For pi, pi^2, etc.
using namespace bamr;

mcmc_bamr::mcmc_bamr() {
  model_type="";
  set=std::make_shared<settings>();
  nsd=std::make_shared<ns_data>();

  bc_arr.resize(1);
  bc_arr[0]=new bamr_class;
  bc_arr[0]->set=set;
  bc_arr[0]->nsd=nsd;

}

#ifdef O2SCL_NEVER_DEFINED

int mcmc_bamr::train(std::string file_name, std::vector<std::string> &names) {
  
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.append('./')");

  // Todo: check to see if threading really works
  //PyEval_InitThreads();
  Py_DECREF(PyImport_ImportModule("threading"));

  // train and test file names
  string train_file = file_name;
  string test_file = "test_data";

  // Import python module
  train_modFile = PyImport_ImportModule("emu");
  if (train_modFile == 0) {
    PyErr_Print();
    std::exit(1);
  }

  // Copy parameter names to python module. This does not
  // currently include the alt_ parameters for the atmosphere
  train_tParam_Names=PyList_New(names.size());
  for(size_t i=0; i<names.size(); i++){
    PyList_SetItem(train_tParam_Names, i, 
                   PyUnicode_FromString(names[i].c_str()));
  }

  // Python class object
  train_trainClass = PyObject_GetAttrString(train_modFile, "modGpr");
  assert(train_trainClass != 0);

  // Create an instance of the modGpr class
  if(PyCallable_Check(train_trainClass)) {
    train_instance = PyObject_CallObject(train_trainClass, 0);
  }
  assert(train_instance != 0);

  if(nsd->n_sources == 0 && !set->apply_intsc){
    addtl_sources = PyLong_FromSize_t(0);
  }
  if(nsd->n_sources>0 && !set->apply_intsc){
    addtl_sources = PyLong_FromSize_t(nsd->n_sources);
  }
  if(nsd->n_sources>0 && set->apply_intsc){
    addtl_sources = PyLong_FromSize_t(nsd->n_sources);
  } 

  // Python arguments for the modGpr::modTrain() function
  train_pArgs = PyTuple_Pack(4, 
                             PyUnicode_FromString(train_file.c_str()),
                             train_tParam_Names, train_tParam_Names,
                             addtl_sources);

  train_trainMthd = PyObject_GetAttrString(train_instance, "modTrain");
  
  // Call Python training function
  if (PyCallable_Check(train_trainMthd)) {
    PyObject_CallObject(train_trainMthd, train_pArgs);
  }

  return 0;
}

int mcmc_bamr::emu_train2(std::vector<std::string> &sv, bool itive_com) {

  vector<string> files;
  for(size_t k=1;k<sv.size();k++) {
    files.push_back(sv[k]);
  }

  eb_arr.resize(n_threads);

  for(size_t k=0;k<n_threads;k++) {
    table_units<> tab;
    //eb_arr[k].train(tab,bc_arr[i],&scr_out);
  }
    
  return 0;
}

int mcmc_bamr::emu_points(std::vector<std::string> &sv, bool itive_com) {


  if(sv.size()<2){
    cout << "Need an emulated output filename." << endl;
  }
  if(sv.size()<3){
    cout << "Need an posterior output filename." << endl;
  }
  string emu_file = sv[1];
  string post_out = sv[2];
  
  // Initial row number
  size_t init_row = 0;
  if(sv.size()<4){
    cout << "Computing postesrior from the first row." << endl;
  }else{
    init_row = o2scl::stoszt(sv[3]);
    cout << "Computing postesrior from row number " << init_row << endl;
  }

  
  // initialize the grids and columns
  mcmc_init();

  // Read emulated file to table
  o2scl::table_units<> emu_init_table;
  o2scl::table_units<> out_table;      
  hdf_file hf_emu;
  hf_emu.open(emu_file);
  hf_emu.get_szt("n_params",this->n_params);
  hdf_input(hf_emu,emu_init_table,"markov_chain_0");
  hf_emu.close();
  cout << "Emulated file copied to table" << endl;

  model_data test_point;
  double log_wgt;
  ubvector emu_pars(n_params);

  model &m=*(bc_arr[0]->mod);
  bamr_class &bc=dynamic_cast<bamr_class &>(*(bc_arr[0]));
  if (set->apply_intsc) {
    bc.setup_filters();
  }

  // Open or create output file
  hdf_file hf_out;

  size_t nrows = emu_init_table.get_nlines();

  int pthread=0;
  bool set_col =false;
  std::clock_t start_time = std::clock();

  
  // Check start time, which can be used to update file after some interval
  // cout << "Start time : " << start_time << endl;
  
  for(size_t i=init_row; i<nrows; i++){

    //cout << "working on row : " << i  << endl;

    // copy parameter values to use in bamr_class::compute_point()
    for(size_t j=5;j<5+n_params;j++) {
      emu_pars(j-5) = emu_init_table.get(emu_init_table.get_column_name(j), i);
    }

    // Compute point success status
    size_t iret = bc.compute_point(emu_pars, scr_out, log_wgt, test_point);

    if(iret==0){

      // copy row from tables in model_data 
      ubvector temp_mvsr_row;
      test_point.mvsr.get_row(i, temp_mvsr_row);
      ubvector temp_eos_row;
      test_point.eos.get_row(i, temp_eos_row);
      ubvector temp_gridt_row;
      test_point.gridt.get_row(i, temp_gridt_row);
      /*
        ubvector temp_sourcet_row;
        test_point.sourcet.get_row(i, temp_sourcet_row);
      */


      cout << "mvsr table size : " << temp_mvsr_row.size() << endl;
      cout << "eos table size : " << temp_eos_row.size() << endl;
      cout << "gridt table size : " << temp_gridt_row.size() << endl;
      //cout << "eos table size : " << temp_sourcet_row.size() << endl;


      exit(0);
      // copy data done.

      vector<string> cols;
      vector<double> col_vals;

      cout << "compute_point return status : " << iret << endl;
      cout << "predicted log_wgt : " <<
        emu_init_table.get(emu_init_table.get_column_name(4), i) << endl;
      cout << "compute_point log_wgt : " << log_wgt << endl;

    std:string temp_const;
      double temp_val;

      for(size_t j=0; j<test_point.mvsr.get_nconsts(); j++){
        test_point.mvsr.get_constant(j, temp_const, temp_val);
        cols.push_back(temp_const);
        col_vals.push_back(temp_val);
      }

      for(size_t j=0; j<test_point.eos.get_nconsts(); j++){
        test_point.eos.get_constant(j, temp_const, temp_val);
        cols.push_back(temp_const);
        col_vals.push_back(temp_val);
      }
      

      // Create/check column names in the new table
      if(set_col==false){
        for(size_t j=0; j<cols.size(); j++){
          out_table.new_column(cols[j]);
        }
        // add M-R grid columns
        for(int ij=0;ij<set->grid_size;ij++) {
          out_table.new_column(((string)"R_")+o2scl::itos(ij));
        }
        /*
        // add EOS grid
        for(int i=0;i<set->grid_size;i++) {
        out_table.new_column(((string)"P_")+o2scl::itos(i));
        }
        */
        set_col=true;
      }

      // Interpolate M-R grid
      test_point.mvsr.set_interp_type(itp_linear);
      for(int ik=0;ik<set->grid_size;ik++) {
        col_vals.push_back(test_point.mvsr.interp("gm",(ik+1)*
                                                  (3.0-0.02)/100.0,"r"));
      }
      
      /*
        for(int i=0;i<set->grid_size;i++) {
        double eval = m.e_grid[i];
        double pres_temp=test_point.eos.interp("ed",eval,"pr");
        col_vals.push_back(pres_temp);
        }
      */    
      // cout << out_table.get_ncolumns() << " " << col_vals.size() << endl;
      out_table.line_of_data(col_vals);
      double duration = (std::clock()-start_time)/(double) CLOCKS_PER_SEC;
      cout << "duration : "<< duration << endl;
      if(duration-60.0 > 0.0 && duration-60.0 < 10.0){
        hf_out.open_or_create(post_out);
        hdf_output(hf_out, out_table, "emulated");
        hf_out.close();
      }

    }else{
      continue;
    }
  }
  
  hf_out.open_or_create(post_out);
  hdf_output(hf_out, out_table, "emulated");
  hf_out.close();
  
  return 0;
}

#endif

int mcmc_bamr::threads(std::vector<std::string> &sv, bool itive_com) {
  
  if (sv.size()==1) {
    cerr << "Number of threads not specified in 'threads'." << endl;
    return 1;                                                          
  }

  if (model_type.length()>0) {
    cerr << "Threads must be set before model." << endl;
    return 2;
  }
  
  size_t n_threads_old=n_threads;
  for(size_t i=0;i<n_threads_old;i++) {
    delete bc_arr[i];
  }
  n_threads=o2scl::stoszt(sv[1]);
  
  bc_arr.resize(n_threads);
  for(size_t i=0;i<n_threads;i++) {
    bc_arr[i]=new bamr_class;
    bc_arr[i]->set=set;
    bc_arr[i]->nsd=nsd;
    bc_arr[i]->n_threads=n_threads;
  }
  
  return 0;
}
  
void mcmc_bamr::file_header(o2scl_hdf::hdf_file &hf) {

  mcmc_para_cli::file_header(hf);
  
  model &m=*(bc_arr[0]->mod);
  
  hf.sets_vec_copy("source_names",nsd->source_names);
  hf.sets_vec_copy("source_fnames",nsd->source_fnames);
  hf.sets_vec_copy("slice_names",nsd->slice_names);

  hf.set_szt("grid_size",set->grid_size);
  hf.set_szt("n_sources",nsd->n_sources);
  hf.sets("method",mcmc_method);
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

  if (bc_arr.size()<1) {
    O2SCL_ERR("Object bc_arr invalid.",o2scl::exc_esanity);
  }
  model &m=*(bc_arr[0]->mod);
  
  // This ensures enough space for all the
  // default return values in models.h
  this->ret_value_counts.resize(this->n_threads);
  for(size_t it=0;it<this->n_threads;it++) {
    // The size must be at least (total # of error codes + 1)
    this->ret_value_counts[it].resize(31);
  }

  // Copy parameter values to all of the model objects
  for(size_t i=1;i<bc_arr.size();i++) {
    model &m2=*(bc_arr[i]->mod);
    m.copy_params(m2);
  }
  
  mcmc_para_cli::mcmc_init();

  // Enable/diable storing rejected MCMC points
  /*if (set->use_emulator) this->store_rejects=true;
  else this->store_rejects=false;*/
  
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

  if (set->emu_tov==false) {

    // -----------------------------------------------------------
    // Add columns to table

    /* These columns are redundant because the output table also 
    contains log_wgt_sources */
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
        this->table->new_column(((string)"cs2_")+o2scl::itos(i));
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
                                "1/fm");
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
      if (set->mmax_deriv) {
        this->table->new_column("dpdM");
        this->table->new_column("M_max2");
        this->table->set_unit("M_max2","Msun");
        if (model_type==((string)"tews_threep_ligo")) {
          this->table->set_unit("dpdM","Msun");
        } else if (model_type==((string)"tews_fixp_ligo")) {
          this->table->set_unit("dpdM","Msun*fm^4");
        } else if (model_type==((string)"new_poly")) {
          this->table->set_unit("dpdM","Msun");
        } else if (model_type==((string)"new_lines")) {
          this->table->set_unit("dpdM","Msun");
        }
      }
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
        this->table->new_column(((string)"I_bar_")+o2scl::itos(i));

        this->table->new_column(((string)"Lambda_bar_")+o2scl::itos(i));
      }
    }

    if (nsd->source_fnames_alt.size()>0) {
      for(size_t i=0;i<nsd->n_sources;i++) {
        this->table->new_column(((string)"atm_")+nsd->source_names[i]);
      }
    }

    if (set->inc_ligo) {
      this->table->new_column("M_chirp_gw17");
      this->table->set_unit("M_chirp_gw17","Msun");
      this->table->new_column("m1_gw17");
      this->table->set_unit("m1_gw17","Msun");
      this->table->new_column("m2_gw17");
      this->table->set_unit("m2_gw17","Msun");
      this->table->new_column("R1");
      this->table->set_unit("R1","km");
      this->table->new_column("R2");
      this->table->set_unit("R2","km");
      this->table->new_column("I1");
      this->table->set_unit("I1","Msun*km^2");
      this->table->new_column("I2");
      this->table->set_unit("I2","Msun*km^2");
      this->table->new_column("I_bar1");
      this->table->new_column("I_bar2");
      this->table->new_column("Lambda1");
      this->table->new_column("Lambda2");
      this->table->new_column("Lambdat");
      this->table->new_column("del_Lambdat");    
      this->table->new_column("prob_gw17");
      this->table->new_column("eta");
      this->table->new_column("m2_gw19");
      this->table->set_unit("m2_gw19","Msun");
      this->table->new_column("prob_gw19");
    }
  }
  
  if (nsd->n_sources>0){
    for(size_t i=0;i<nsd->n_sources;i++) {
      this->table->new_column(((string)"log_wgt_")+
                              nsd->source_names[i]);
    }
  }
  if (set->inc_pop) {
    this->table->new_column("log_wgt_NS");
    this->table->new_column("log_wgt_WD");
    this->table->new_column("log_wgt_LMS");
    this->table->new_column("log_wgt_pop");
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
  // Setup filters

  for(size_t i=0;i<n_threads;i++) {
    bamr_class &bc=dynamic_cast<bamr_class &>(*(bc_arr[i]));
    if (set->apply_intsc) {
      bc.setup_filters();
    }
  }

  // Read FFT cache
  
  if (set->cached_intsc) {
    for(size_t i=0;i<n_threads;i++) {
      bamr_class &bc=dynamic_cast<bamr_class &>(*(bc_arr[i]));
      hdf_file hfx;
      for(size_t ii=0;ii<nsd->n_sources;ii++) {
        string fname=((string)"data/cache/tg_")+szttos(ii)+"_0";
        hfx.open(fname);
        hdf_input(hfx,bc.fft_data[ii*2],"tg");
        hfx.close();
        fname=((string)"data/cache/tg_")+szttos(ii)+"_1";
        hfx.open(fname);
        hdf_input(hfx,bc.fft_data[ii*2+1],"tg");
        hfx.close();
      }
    }
  }

  if (this->verbose>=2) {
    std::cout << "(rank " << this->mpi_rank
              << ") End mcmc_bamr::mcmc_init()." << std::endl;
  }

  return 0;
}

int mcmc_bamr::set_method(std::vector<std::string> &sv, bool itive_com) {
  
  if (sv.size()<2) {
    cerr << "MCMC method not given." << endl;
    return exc_efailed;
  }
  if (mcmc_method==sv[1]) {
    cerr << "Method already set to " << sv[1] << endl;
    return 0;
  }
  mcmc_method=sv[1];
  bc_arr[0]->mcmc_method=sv[1];
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
  } else if (sv[1]==((string)"tews_threep_ligo")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new tews_threep_ligo(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"tews_fixp_ligo")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new tews_fixp_ligo(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"new_poly")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new new_poly(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"new_lines")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new new_lines(set,nsd));
      bc_arr[i]->mod=mnew;
      bc_arr[i]->model_type=sv[1];
    }
  } else if (sv[1]==((string)"prx")) {
    for(size_t i=0;i<n_threads;i++) {
      std::shared_ptr<model> mnew(new prx(set,nsd));
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

int mcmc_bamr::combine_files(std::vector<std::string> &sv,
			     bool itive_com) {

  if (sv.size()<3) {
    O2SCL_ERR("Not enough args combine_files.",
	      o2scl::exc_einval);
  }
  size_t n_files=sv.size()-2;
  string file_final=sv[sv.size()-1];
  table_units<> t_final;

  for(size_t i_file=0;i_file<n_files;i_file++) {
  
    // Read file
    hdf_file hf;
    hf.open(sv[i_file+1]);
    table_units<> t0;
    hdf_input(hf,t0);
    
    if (t0.is_column("walker")==false || t0.is_column("thread")==false) {
      O2SCL_ERR("No walker or thread columns.",
		o2scl::exc_einval);
    }
    
    // Remove empty rows from table.
    t0.delete_rows_func("mult<0.5");
  
    // Compute number of walkers and threads
    size_t n_walker=((size_t)(t0.max("walker")+1.00001));
    size_t n_thread=((size_t)(t0.max("thread")+1.00001));
    cout << "Determined file has " << n_walker << " walkers, " << n_thread
	 << " threads, and " << t0.get_nlines() << " lines." << endl;
    if (t0.get_nlines()<n_walker*n_thread) {
      O2SCL_ERR("Not enough lines in table.",
		o2scl::exc_einval);
    }
    
    // Remove first point
    vector<size_t> list;
    list.push_back(0);
    for(size_t j=1;j<n_walker*n_thread;j++) {
      list.push_back(j);
    }
    t0.delete_rows_list(list);
    cout << "Table now has " << t0.get_nlines() << " lines." << endl;
    
    // Remove last point
    list.clear();
    for(size_t j=0;j<n_walker;j++) {
      for(size_t k=0;k<n_thread;k++) {
	for(int i=t0.get_nlines()-1;i>=0;i--) {
	  if (fabs(t0.get("walker",i)-((double)j))<1.0e-4 &&
	      fabs(t0.get("thread",i)-((double)k))<1.0e-4) {
	    list.push_back(i);
	    cout << "Deleting row " << i << " for thread "
		 << k << " and walker " << j << endl;
	    i=-1;
	  }
	}
      }
    }
    t0.delete_rows_list(list);
    cout << "Table now has " << t0.get_nlines() << " lines." << endl;

    cout << "Copying " << t0.get_nlines() << " lines to final table." << endl;
    
    if (i_file==0) {
      
      t_final=t0;
      
    } else {

      if (t0.get_ncolumns()!=t_final.get_ncolumns()) {
	O2SCL_ERR("Different number of columns.",
		  o2scl::exc_einval);
      }

      size_t t_final_orig=t_final.get_nlines();
      t_final.set_nlines(t_final.get_nlines()+t0.get_nlines());
      for(size_t i=0;i<t0.get_nlines();i++) {
	for(size_t j=0;j<t0.get_ncolumns();j++) {
	  t_final.set(j,t_final_orig+i,t0.get(j,i));
	}
      }
    }
    
  }

  cout << "Selecting independent samples from " << t_final.get_nlines()
       << " lines." << endl;

  // Compute the autocorrelation length from log_wgt
  if (t_final.get_nlines()!=t_final.get_maxlines()) {
    t_final.set_maxlines(t_final.get_nlines());
  }
  
  /*
    AWS: commenting these out temporarily because they depend on
    a more recent o2scl
    
  std::vector<double> ac, ftom;
  o2scl::vector_autocorr_vector_fftw_mult(t_final["log_wgt"],
                                          t_final["mult"],ac);

  size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom,3);
  cout << "Autocorrelation length is: " << ac_len << endl;

  // Create a separate table of statistically independent samples
  table_units<> indep;
  copy_table_thin_mcmc(ac_len,t_final,indep,"mult",3);
  
  cout << "Found " << indep.get_nlines() << " independent samples." << endl;
  
  hdf_file hfx;
  hfx.open_or_create(file_final);
  hdf_output(hfx,indep,"markov_chain_0");
  hfx.close();
  */
  
  return 0;
}

int mcmc_bamr::initial_point_last(std::vector<std::string> &sv,
                                  bool itive_com) {

  if (sv.size()<2) {
    cerr << "Need a filename for initial_point_last()." << endl;
    return 1;
  }

  if (model_type.length()<2) {
    cerr << "Model not specified in initial_point_last()." << endl;
    return 2;
  }
      
  size_t np;
  size_t np_ligo=nsd->n_ligo_params;
  size_t np_eos=bc_arr[0]->mod->n_eos_params;
  size_t np_src=nsd->n_sources;
  size_t np_pop=nsd->pop.n_pop_params;
  string fname=sv[1];
  size_t pos=fname.find("<rank>");

  if (pos!=std::string::npos) {
    fname.replace(pos,6,o2scl::itos(mpi_rank));
  }

  bc_arr[0]->init_file=fname;

  // Determine the number of parameters according to the
  // settings
  if (set->inc_ligo && set->inc_pop) {
    np=np_eos+np_ligo+np_src+np_pop;
  } else if (set->inc_ligo && !set->inc_pop) {
    np=np_eos+np_ligo+np_src;
  } else if (!set->inc_ligo && set->inc_pop) {
    np=np_eos+np_src+np_pop;
  } else {
    np=np_eos+np_src;
  }

  cout << "mcmc_bamr::initial_point_last set: " << np << " parameters."
       << endl;
  this->initial_points_file_last(fname,np);
  
  return 0;
}

int mcmc_bamr::initial_point_best(std::vector<std::string> &sv,
                                  bool itive_com) {
  
  if (sv.size()<2) {
    cerr << "Need a filename for initial_point_best()." << endl;
    return 1;
  }
      
  if (model_type.length()<2) {
    cerr << "Model not specified in initial_point_best()." << endl;
    return 2;
  }
  
  size_t np;
  size_t np_ligo=nsd->n_ligo_params;
  size_t np_eos=bc_arr[0]->mod->n_eos_params;
  size_t np_src=nsd->n_sources;
  size_t np_pop=nsd->pop.n_pop_params;
  string fname=sv[1];
  size_t pos=fname.find("<rank>");
  
  if (pos!=std::string::npos) {
    fname.replace(pos,6,o2scl::itos(mpi_rank));
  }

  if (set->inc_ligo && set->inc_pop) {
    np=np_eos+np_ligo+np_src+np_pop;
  }
  else if (set->inc_ligo && !set->inc_pop) {
    np=np_eos+np_ligo+np_src;
  }
  else if (!set->inc_ligo && set->inc_pop) {
    np=np_eos+np_src+np_pop;
  }
  else {
    np=np_eos+np_src;
  }

  this->initial_points_file_best(fname, np);
  
  return 0;
}

int mcmc_bamr::read_prev_results_mb(std::vector<std::string> &sv,
                                    bool itive_com) {

  O2SCL_ERR("Not implemented yet.",o2scl::exc_eunimpl);
  
  if (sv.size()<2) {
    cerr << "Need a filename for read_prev_results_mb()." << endl;
    return 1;
  }

  if (model_type.length()<2) {
    cerr << "Model not specified in read_prev_results_mb()." << endl;
    return 2;
  }
  
  model &m=*(bc_arr[0]->mod);
  size_t np=m.n_eos_params+nsd->n_sources;
  
  // Ensure that multiple threads aren't reading from the 
  // filesystem at the same time
#ifdef BAMR_MPI
  int tag=0, buffer=0;
  if (mpi_size>1 && mpi_rank>=1) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
             tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
#endif
  
  string fname=sv[1];
  size_t pos=fname.find("<rank>");
  if (pos!=std::string::npos) {
    fname.replace(pos,6,o2scl::itos(mpi_rank));
  }
  cout << "Rank " << mpi_rank
       << " is reading previous results from " << fname << " ." << endl;
  hdf_file hf;
  hf.open(fname);
  mcmc_para_table::read_prev_results(hf,np);
  hf.close();
  
#ifdef BAMR_MPI
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
             tag,MPI_COMM_WORLD);
  }
#endif

  return 0;
}

#ifdef ANDREW
int mcmc_bamr::point_wrapper(size_t it, size_t np, const ubvector &p,
                             double &log_wgt, model_data &dat) {

  // We have to run the exact code a couple of times to make sure the
  // model data object is filled with data before we run the
  // emulator
  
  if (dat.eos.get_nlines()>0 && dat.mvsr.get_nlines()>0 &&
      dat.gridt.get_nlines()>0 && dat.mvsr.is_constant("R_max") &&
      dat.mvsr.is_constant("P_max") && n_retrain>0) {
    
    //std::cout << "Point emu. " << &dat << std::endl;
    ubvector out(1);
    emu[it]->eval(p,out);
    log_wgt=out[0];
    //std::cout << "Point emu done." << std::endl;
    
  } else {
    
    std::cout << "Point exact. " << &dat << std::endl;
    int ret=((*func_ptr)[it])(np,p,log_wgt,dat);
    std::cout << "Point exact done ret=" << ret << std::endl;
    return ret;
    
  }

  return 0;
}

#endif

int mcmc_bamr::mcmc_func(std::vector<std::string> &sv, bool itive_com) {

  if (model_type.length()==0) {
    cerr << "Model not set in 'mcmc' command." << endl;
    return 1;
  }

  std::vector<std::string> names;
  std::vector<std::string> units;

  vector<double> low, high;
  // Get upper and lower parameter limits and also the column names
  // and units for the data table (which also automatically includes
  // nuisance variables for the data points). The other columns and
  // units are specified in mcmc_init() function manually using a call
  // to table::new_column().
  bc_arr[0]->mod->get_param_info(names,units,low,high); 

  nsd->data_params(names,units,low,high,set);
  
  if (set->apply_intsc) {

    for(size_t i=0;i<nsd->n_sources;i++) {
      names.push_back(((string)"log10_is_")+nsd->source_names[i]);
      units.push_back("");
      low.push_back(-2.0);
      high.push_back(2.0);
    }
  }
  
  // Population parameters
  if (set->inc_pop) {
    
    ns_pop &pop=nsd->pop;

    // Set names, units, low, high for population parameters
    for (size_t i=0; i<pop.n_pop_params; i++) {
      names.push_back(pop.par_names[i]);
      units.push_back(pop.par_units[i]);
      low.push_back(pop.par_low[i]);
      high.push_back(pop.par_high[i]);
    }

  }

#ifdef O2SCL_NEVER_DEFINED
 if (set->apply_emu) {    
    for(size_t i=0;i<nsd->n_sources;i++) {
      names.push_back(((string)"atm_")+o2scl::szttos(i));
      units.push_back("");
      low.push_back(0.0);
      high.push_back(1.0);
    }
  }
#endif

  // Send names and units to o2scl
  set_names_units(names,units);
  
  // Set initial points if they have not already been set by the user
  if (this->initial_points.size()==0) {
    
    // Get the parameter initial values for this model 
    vector<double> init;
    
    bc_arr[0]->mod->initial_point(init);

    nsd->initial_point(set,init);
    
    if (set->apply_intsc) {
      for (size_t i=0; i<nsd->n_sources; i++) {
        init[i+bc_arr[0]->mod->n_eos_params+nsd->n_sources]=-0.5;
      }
    }
    /* Note:
       bc_arr[0]->mod->n_eos_params+nsd->n_sources=23
       bc_arr[0]->mod->n_eos_params=12 (9+3)
       nsd->n_sources=11 
    */

    // Set initial points for the population parameters
    if (set->inc_pop) {
      
      ns_pop &pop=nsd->pop;

      for (size_t i=0; i<pop.n_pop_params; i++) {
        init.push_back(pop.par_init[i]);
      }
    }
    
    // AWS: 3/20/18: I changed this part, because we want the MCMC class
    // to be able to expect that different entries in the initial_points
    // array
    this->initial_points.clear();
    
    ubvector init2(init.size());
    vector_copy(init,init2);
    this->initial_points.push_back(init2);
    
    if (this->verbose>1) {
      cout << "Summary of default initial point: " << endl;
      cout << "Sizes of names, units, low, high, and init: "
           << names.size() << " " << units.size() << " "
           << low.size() << " " << high.size() << " "
           << init.size() << endl;
      cout << "Parameter index, name, unit, low, init, high: " << endl;
      for(size_t j=0;j<names.size();j++) {
        cout.width(3);
        cout << j << " ";
        cout.width(18);
        cout << names[j] << " ";
        cout.width(6);
        cout.setf(ios::left);
        cout << units[j] << " ";
        cout.unsetf(ios::left);
        cout.setf(ios::showpos);
        cout << low[j] << " " << init[j] << " " << high[j];
        cout.unsetf(ios::showpos);
        if (init[j]<low[j]) {
          cout << " L";
        }
        if (init[j]>high[j]) {
          cout << " H";
        }
        cout << endl;
      }
    }
    
  } else {

    if (this->verbose>1) {
      cout << "Parameters index, name, unit, low, high: " << endl;
      for(size_t j=0;j<names.size();j++) {
        cout.width(3);
        cout << j << " ";
        cout.setf(ios::left);
        cout.width(18);
        cout << names[j];
        cout << " ";
        cout.width(6);
        cout << units[j] << " ";
        cout.unsetf(ios::left);
        cout.setf(ios::showpos);
        cout << low[j] << " " << initial_points[0][j] << " "
             << high[j] << endl;
        cout.unsetf(ios::showpos);
      }
    }
  }

#ifdef O2SCL_NEVER_DEFINED
  if (set->emu_aws) {
    if (eb_arr.size()!=n_threads) {
      cerr << "Not enough emulators." << endl;
    }
  }
#endif  

  vector<bamr::point_funct> pfa(n_threads);
  vector<bamr::fill_funct> ffa(n_threads);
  for(size_t i=0;i<n_threads;i++) {
    if (this->store_rejects==true) { 
      pfa[i]=std::bind
        (std::mem_fn<int(const ubvector &,ofstream &,double &,model_data &)>
         (&bamr_class::compute_point_ext),bc_arr[i],std::placeholders::_2,
         std::ref(scr_out),std::placeholders::_3,std::placeholders::_4);

#ifdef O2SCL_NEVER_DEFINED
    if (set->emu_aws) {
      pfa[i]=std::bind
        (std::mem_fn<int(size_t n,const ubvector &,double &,model_data &)>
         (&emulator_bamr::eval),eb_arr[i],std::placeholders::_1,
         std::placeholders::_2,std::placeholders::_3,
         std::placeholders::_4);
    }
#endif

    } else {
      pfa[i]=std::bind
        (std::mem_fn<int(const ubvector &,ofstream &,double &,model_data &)>
         (&bamr_class::compute_point),bc_arr[i],std::placeholders::_2,
         std::ref(scr_out),std::placeholders::_3,std::placeholders::_4);
    }
    ffa[i]=std::bind
      (std::mem_fn<int(const ubvector &,double,vector<double> &,
                       model_data &)>
       (&bamr_class::fill),bc_arr[i],std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
  }

#ifdef O2SCL_NEVER_DEFINED  
  if (set->apply_emu) {
    cout << "Applying train function." << endl;

    // train the module
    int pinfo=train(set->emu_train, names);
    if(pinfo != 0){
      cout << "Training Failed. " << endl;
      exit(-1);
    }

    // Copy trained method to bint classes
    for(size_t i=0;i<n_threads;i++){
      bamr_class &bc=dynamic_cast<bamr_class &>(*(bc_arr[i]));

      // copy pyobject to bint class
      //bc.emu_train=emu_train;
      bc.train_modFile=train_modFile;
      bc.train_trainClass=train_trainClass;
      bc.train_instance=train_instance;
      bc.train_trainMthd=train_trainMthd;
      bc.train_tParam_Names=train_tParam_Names;
      bc.addtl_sources=addtl_sources;
    }

    // Delete unnecessary PyObjects
    Py_DECREF(train_modFile);
    Py_DECREF(train_instance);
    Py_DECREF(train_trainClass);
  }
#endif

  //-------------------------------------------------------------------

  // Note that kde_python doesn't work with n_threads>1

#ifdef ANDREW

  if (true) {
    
    this->n_retrain=1000;
    this->emu_file="interp";
    this->show_emu=1;
    this->max_train_size=10000;
    this->test_emu_file="test_emu.o2";

    this->emu.resize(1);

    int intp=3;
    if (intp==1) {
      
      // Set up the shared pointer to the interpolation object
      std::shared_ptr<interpm_idw<boost::numeric::ublas::vector<double>,
				  o2scl::const_matrix_view_table<>,
				  o2scl::matrix_view_table<>>> ii
	(new interpm_idw<boost::numeric::ublas::vector<double>,
	 o2scl::const_matrix_view_table<>,
	 o2scl::matrix_view_table<>>);
      this->emu[0]=ii;
      ii->n_extra=5;
      
    } else if (intp==2) {

      std::shared_ptr<interpm_krige_optim<>> iko(new interpm_krige_optim<>);
      typedef const const_matrix_row_gen
	<o2scl::const_matrix_view_table<>> mat_x_row_t;
         
      // Setup the multidimensional covariance object
      vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
      vmfrn.resize(1);
      std::shared_ptr<mcovar_funct_rbf_noise<
	ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
				    mat_x_row_t>);
      vmfrn[0]=mfrn;
      mfrn->len.resize(names.size());
      
      iko->verbose=2;
      iko->def_mmin.verbose=1;
      iko->full_min=true;
      vector<double> len_list={0.1,0.3};
      vector<double> l10_list={-15,-13};
      vector<vector<double>> ptemp;
      for(size_t j=0;j<names.size();j++) {
	ptemp.push_back(len_list);
      }
      ptemp.push_back(l10_list);
      vector<vector<vector<double>>> param_lists;
      param_lists.push_back(ptemp);
      
      iko->set_covar(vmfrn,param_lists);

      this->emu[0]=iko;
      
    } else if (intp==3) {

      std::shared_ptr<interpm_python<boost::numeric::ublas::vector<double>,
				     o2scl::const_matrix_view_table<>,
				     o2scl::matrix_view_table<>>> ip
	(new interpm_python<boost::numeric::ublas::vector<double>,
	 o2scl::const_matrix_view_table<>,
	 o2scl::matrix_view_table<>>("o2sclpy","set_data_str","eval","eval_unc",
				     "interpm_sklearn_gp",
				     ((std::string)"verbose=1,")+
				     "transform_in=quant,"+
                                     "normalize_y=True",1));
      this->emu[0]=ip;
      
    } else {

      std::shared_ptr<interpm_python<boost::numeric::ublas::vector<double>,
				     o2scl::const_matrix_view_table<>,
				     o2scl::matrix_view_table<>>> ip
	(new interpm_python<boost::numeric::ublas::vector<double>,
	 o2scl::const_matrix_view_table<>,
	 o2scl::matrix_view_table<>>("o2sclpy","set_data_str","eval","eval",
				     "interpm_tf_dnn",
				     ((std::string)"verbose=1,")+
				     "transform_in=quant,"+
                                     "transform_out=quant,"+
                                     "hlayers=[200,400,200]",1));
      this->emu[0]=ip;
      
    }
  }
  
  // #ifdef BAMR_KDE
  if (mcmc_method==string("kde") ||
      mcmc_method==string("kde_sklearn") ||
      mcmc_method==string("gauss")) {
    
    // Copy the table data to a tensor for use in kde_python.
    // We need a copy for each thread because kde_python takes
    // over the tensor data.
    
    // AWS: I'm leaving this for Anik to change
    
#ifdef BAMR_MPI
    // Get MPI rank, etc.
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
      
    // Ensure that multiple MPI ranks aren't reading from the
    // filesystem at the same time
    int tag=0, buffer=0;
    if (mpi_size>1 && mpi_rank>=1) {
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
         tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif

    // Input table object
    o2scl::table_units<> tab_in;

    hdf_file hf;

    // Set input file name based on the model type
    string fname;

    if (model_type==string("new_lines")) {
      if (set->model_dpdm==1) fname="out/ml_train";
      else fname="out/nl_train";
    }
    else if (model_type==string("new_poly")) {
      if (set->model_dpdm==1) fname="out/mp_train";
      else fname="out/np_train";
    }
    
    // Fill input data
    size_t n_pars=names.size(); 
    cout << "KDE is reading the input table of " << n_pars 
         << " parameters from file " << fname << endl;
    hf.open(fname);
    hdf_input(hf,tab_in);
    hf.close();

    // Print column names
    /*for (size_t i=0; i<n_pars; i++) {
      cout << tab_in.get_column_name(i+5) << " ";
    }
    cout << endl;*/

    // Copy input data to the tensor object
    vector<size_t> in_size={tab_in.get_nlines(),n_pars};
    ten_in.resize(2,in_size);
    for (size_t j=0; j<tab_in.get_nlines(); j++) {
      vector<size_t> ix;
      for (size_t i=0; i<n_pars; i++) {
	      ix={j,i};
	      ten_in.get(ix)=tab_in.get(i+5,j);
      }
    }
    
    std::shared_ptr<mcmc_stepper_mh<point_funct,
                                    model_data,ubvector,
                                    ubmatrix,prob_cond_mdim_indep<>>>
      mh_stepper(new mcmc_stepper_mh<point_funct,model_data,
                 ubvector,ubmatrix,prob_cond_mdim_indep<>>);
    stepper=mh_stepper;
    
    // Train the KDE
    if (mcmc_method==string("kde")) {

      // Weights can be set here, but they are presumed to be
      // the same if this vector is empty
      vector<double> weights;

      kp=std::shared_ptr<kde_python<ubvector>>(new kde_python<ubvector>);
      kp->set_function("o2sclpy", ten_in, weights,
		       "verbose=0", "kde_scipy");
      
      // Setting the KDE as the base distribution for the independent
      // conditional probability. The kde_python class does not work
      // for more than one OpenMP thread.
      mh_stepper->proposal.resize(1);
      mh_stepper->proposal[0].set_base(kp);
    
    } else if (mcmc_method==string("kde_sklearn")) {
      
      kp=std::shared_ptr<kde_python<ubvector>>(new kde_python<ubvector>);
      uniform_grid_log_end<double> ug(1.0e-3,1.0e3,99);
      vector<double> bw_array;
      ug.vector(bw_array);
      kp->set_function("o2sclpy", ten_in, bw_array,
		       "verbose=0", "kde_sklearn");
      
      // Setting the KDE as the base distribution for the independent
      // conditional probability. The kde_python class does not work
      // for more than one OpenMP thread.
      mh_stepper->proposal.resize(1);
      mh_stepper->proposal[0].set_base(kp);
      
    } else {

      std::cout << "Setting up Gaussian:" << std::endl;
      ubvector std(n_pars), avg(n_pars);
      cout << "j param,avg,std: " << endl;
      for(size_t j=0;j<n_pars;j++) {
	      avg[j]=vector_mean(tab_in.get_nlines(),tab_in[j+5]);
	      std[j]=vector_stddev(tab_in.get_nlines(),tab_in[j+5]);
	      cout << "param,avg,stddev: " << j << " " << avg[j] 
	           << " " << std[j] << endl;
      }
      
      ubmatrix covar(n_pars,n_pars);
      for(size_t i=0;i<n_pars;i++) {
	      for(size_t j=0;j<n_pars;j++) {
	        if (i==j) {
	          covar(i,j)=std[j]*std[j];
	          //covar(i,j)/=var_dec_factor;
	        } else {
	          covar(i,j)=vector_covariance(tab_in.get_nlines(),tab_in[i+5],
					    tab_in[j+5]);
	          covar(i,j)/=2.0;
	        }
	      }
      }
      
      gpp=std::shared_ptr<prob_dens_mdim_gaussian<>>
	    (new prob_dens_mdim_gaussian<>);
      gpp->set_covar(n_pars,avg,covar);
      gpp->pdg.set_seed(mpi_rank*clock());
      
      // Setting the KDE as the base distribution for the independent
      // conditional probability. This code may need to be changed
      // for more than one OpenMP thread.
      mh_stepper->proposal.resize(1);
      mh_stepper->proposal[0].set_base(gpp);
    
    }
    
#ifdef BAMR_MPI
    // Send a message to the next MPI rank
    if (mpi_size>1 && mpi_rank<mpi_size-1) {
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
         tag,MPI_COMM_WORLD);
    }
#endif

  }
  
#endif

  if (mcmc_method==string("hmc")) {

#ifdef BAMR_MPI
    // Ensure that multiple MPI ranks aren't reading from the
    // filesystem at the same time
    int tag=0, buffer=0;
    if (mpi_size>1 && mpi_rank>=1) {
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
         tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif
    
    std::shared_ptr<mcmc_stepper_hmc<point_funct,
                                     model_data,ubvector>> hmc_stepper
      (new mcmc_stepper_hmc<point_funct,model_data,ubvector>);
    stepper=hmc_stepper;
    
    size_t np=names.size();
    
    hmc_stepper->auto_grad.resize(np);
    for (size_t i=0; i<np; i++) {
      hmc_stepper->auto_grad[i]=false;
    }

    hmc_stepper->hmc_step.resize(np);

    // Scale the step sizes
    for (size_t i=0; i<34; i++) {
      hmc_stepper->hmc_step[i]=1.0e-2*(high[i]-low[i]);
    }
    
    // Mass parameters
    pop_data &pd=nsd->pd;
    for (size_t i=0; i<pd.n_stars; i++) {
      double width=2.0*min(pd.lo_nsp[i],pd.hi_nsp[i]);
      hmc_stepper->hmc_step[34+i]=1.0e-2*width;
    }

    hmc_stepper->traj_length=1;

    vector<bamr::deriv_funct> gfa(n_threads);
    using namespace std::placeholders;
    for (size_t i=0; i<n_threads; i++) {
      gfa[i]=std::bind(std::mem_fn<int(ubvector &,
        point_funct &,ubvector &,model_data &,bool &)>
        (&bamr_class::compute_deriv),bc_arr[i],_2,_3,_4,_5,_6);
    }

    hmc_stepper->set_gradients(gfa);

#ifdef BAMR_MPI
    // Send a message to the next MPI rank
    if (mpi_size>1 && mpi_rank<mpi_size-1) {
      int tag=0, buffer=0;
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
         tag,MPI_COMM_WORLD);
    }
#endif

  }

  // ---------------------------------------

  for(size_t j=0;j<names.size();j++) {
    pvi.append(names[j]);
  }
  // Copy the pvi object from mcmc_bamr to bamr_class so it can
  // be used in compute_point()
  for(size_t j=0;j<bc_arr.size();j++) {
    bc_arr[j]->pvi=pvi;
  }
  
  // Perform the MCMC simulation
  ubvector low2(low.size()), high2(high.size());
  vector_copy(low,low2);
  vector_copy(high,high2);

  std::vector<model_data> dat_arr(2*this->n_walk*this->n_threads);

  if (verbose>2) {
    cout << "In mcmc_bamr::mcmc_func(): Going to mcmc_fill()." << endl;
  }
#ifdef ANDREW
  this->mcmc_emu(names.size(),low2,high2,pfa,ffa,dat_arr);
#else
  this->mcmc_fill(names.size(),low2,high2,pfa,ffa,dat_arr);
#endif

#ifdef O2SCL_NEVER_DEFINED  
  if (set->apply_emu) {
    Py_Finalize();
  }
#endif

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

void mcmc_bamr::setup_cli_mb() {

  mcmc_para_cli::setup_cli(cl);

  set->setup_cli(cl);
  
  // ---------------------------------------
  // Set options
    
  static const int nopt=10; // nopt=12 with commented out 2 options
  comm_option_s options[nopt]=
    {
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
      {'t',"method","Choose MCMC method.",
        1,1,"<method type>",((string)"Choose the MCMC sampling method. ")+
        "Possible values are 'kde' and 'hmc'.",
        new comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::set_method),
        cli::comm_option_both},
      {0,"threads","Specify number of OpenMP threads",
       1,1,"<number>",((string)"The threads command must be ")+
       "before any model selection, any changes to the data, and "+
       "any changes to settings.",
       new comm_option_mfptr<mcmc_bamr>(this,&mcmc_bamr::threads),
       cli::comm_option_both},
      {'a',"add-data","Add data source to the list.",
       4,5,"<name> <file> <slice> <initial mass> [obj name]",
       ((string)"Specify data as a table3d object in a HDF5 file. ")+
       "The string <name> is the name used, <file> is the filename, "+
       "<slice> is the name of the slice in the table3d object, "+
       "<initial mass> is the initial mass for the first point, and "+
       "[obj name] is the optional name of table3d object in <file>. "+
       "If [obj name] is not specified, then the first table3d object "+
       "is used.",new comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::add_data),
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
       "is used.",new comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::add_data_alt),
       cli::comm_option_both},
      {0,"initial-point-last","Set initial point from file.",1,1,
       "<filename>",((string)"Note that this command must be ")+
       "after the model selection and, after the setting of inc_pop "+
       " and inc_ligo, "+
       "and after the specification of the data, so that it can "+
       "determine the correct number of parameters to obtain from "+
       "the file.",
       new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::initial_point_last),
       o2scl::cli::comm_option_both},
      {0,"initial-point-best","Set initial point from file.",1,1,
       "<filename>",((string)"Note that this command must be ")+
       "after the model selection and, after the setting of inc_pop "+
       " and inc_ligo, "+
       "and after the specification of the data, so that it can "+
       "determine the correct number of parameters to obtain from "+
       "the file.",
       new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::initial_point_best),
       o2scl::cli::comm_option_both},
      {0,"combine-files","Desc.",-1,-1,"<file 1> <file 2> ... <target file>","",
       new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::combine_files),
       o2scl::cli::comm_option_both},
      {0,"read-prev-results","Read previous results from file (unfinished).",
       1,1,"<filename>","Long. desc.",
       new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::read_prev_results_mb),
       o2scl::cli::comm_option_both} // Insert a comma to add more options
      /*{0,"emu-points","emu-points help.",
       2,3,"<input filename> <output filename> <strting row number>",
       "Long description.",new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::emu_points),
       o2scl::cli::comm_option_both},
      {0,"emu-train","emu-train help.",
       -1,-1,"<training file 1> [training file 2] ...",
       "Long description.",new o2scl::comm_option_mfptr<mcmc_bamr>
       (this,&mcmc_bamr::emu_train2),
       o2scl::cli::comm_option_both}*/
    };
  cl.set_comm_option_vec(nopt,options);

  // --------------------------------------------------------
  
  return;
}

