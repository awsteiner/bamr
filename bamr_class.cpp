/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
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
#include "bamr_class.h"

#include <o2scl/vector.h>
#include <o2scl/hdf_io.h>

#include "likelihood.h"

using namespace std;
using namespace o2scl;
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

void bamr_class::setup_filters() {

#ifdef BAMR_FFTW3

  flt.resize(n_threads*nsd->n_sources);
  
  // Input and output table references for convenience
  std::vector<o2scl::table3d> &in=nsd->source_tables;
  std::vector<o2scl::table3d> &in_alt=nsd->source_tables_alt;
  std::vector<o2scl::table3d> &out=source_tables_is;
  std::vector<o2scl::table3d> &out_alt=source_tables_alt_is;
  
  // Copy the original tables over if we're running this
  // code for the first time
  out.resize(nsd->n_sources*n_threads);
  out_alt.resize(nsd->n_sources*n_threads);
  for (size_t i=0;i<nsd->n_sources*n_threads;i++) {
    out[i]=in[i % nsd->n_sources];
    out_alt[i]=in_alt[i % nsd->n_sources];
  }
  
  int mpi_rank=0, mpi_size=1;
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

  for(int k=0;k<n_threads;k++) {
    for(size_t j=0;j<nsd->n_sources;j++) {
      size_t Nx = in[j].get_nx();
      size_t Ny = in[j].get_ny();
      flt[k*nsd->n_sources+j]=new filters::Filter(Nx,Ny);
    }
  }

#ifdef BAMR_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
             tag,MPI_COMM_WORLD);
  }
#endif

#endif

  return;
}

int bamr_class::fill(const ubvector &pars, double weight, 
		     std::vector<double> &line, model_data &dat) {

  if (set->apply_emu) {
    return 0;
  } else {

    for(size_t i=0;i<nsd->n_sources;i++) {
      line.push_back(dat.sourcet.get("wgt",i));
    }
    for(size_t i=0;i<nsd->n_sources;i++) {
      line.push_back(dat.sourcet.get("R",i));
    }
    for(size_t i=0;i<nsd->n_sources;i++) {
      line.push_back(dat.sourcet.get("M",i));
    }
    
    // Reference to model object for convenience
    model &m=*this->mod;

    if (m.has_eos) {
      for(int i=0;i<set->grid_size;i++) {
        line.push_back(dat.gridt.get("P",i));
        line.push_back(dat.gridt.get("cs2",i));
      }
    }

    for(int i=0;i<set->grid_size;i++) {
      line.push_back(dat.gridt.get("R",i));
      line.push_back(dat.gridt.get("PM",i));
    }
    
    if (m.has_eos) {
      if (set->baryon_density) {
        for(int i=0;i<set->grid_size;i++) {
          line.push_back(dat.gridt.get("Pnb",i));
          line.push_back(dat.gridt.get("EoA",i));
        }
      }
      
      if (m.has_esym) {
        line.push_back(dat.eos.get_constant("S"));
        line.push_back(dat.eos.get_constant("L"));
      }
      
      line.push_back(dat.mvsr.get_constant("R_max"));
      line.push_back(dat.mvsr.get_constant("M_max"));
      if (set->mmax_deriv) {
        line.push_back(dat.eos.get_constant("dpdM"));
      }
      line.push_back(dat.mvsr.get_constant("P_max"));
      line.push_back(dat.mvsr.get_constant("e_max"));
      if (set->baryon_density) {
        line.push_back(dat.mvsr.get_constant("nb_max"));
      }
      
      for(size_t i=0;i<nsd->n_sources;i++) {
        line.push_back(dat.sourcet.get("ce",i));
      }
      if (set->baryon_density) {
        for(size_t i=0;i<nsd->n_sources;i++) {
          line.push_back(dat.sourcet.get("cnb",i));
        }
      }
    }
    
    if (set->baryon_density) {
      line.push_back(dat.mvsr.get_constant("gm_nb1"));
      line.push_back(dat.mvsr.get_constant("r_nb1"));
      line.push_back(dat.mvsr.get_constant("gm_nb2"));
      line.push_back(dat.mvsr.get_constant("r_nb2"));
      line.push_back(dat.mvsr.get_constant("gm_nb3"));
      line.push_back(dat.mvsr.get_constant("r_nb3"));
      line.push_back(dat.mvsr.get_constant("gm_nb4"));
      line.push_back(dat.mvsr.get_constant("r_nb4"));
      line.push_back(dat.mvsr.get_constant("gm_nb5"));
      line.push_back(dat.mvsr.get_constant("r_nb5"));
    }
    
    if (set->compute_cthick) {
      line.push_back(dat.eos.get_constant("nt"));
      line.push_back(dat.eos.get_constant("prt"));
      for(int i=0;i<set->grid_size;i++) {
        line.push_back(dat.gridt.get("CT",i));
      }
    }

    if (set->addl_quants) {
      for(int i=0;i<set->grid_size;i++) {
        line.push_back(dat.gridt.get("MB",i));
        line.push_back(dat.gridt.get("BE",i));
        line.push_back(dat.gridt.get("I",i));
        line.push_back(dat.gridt.get("I_bar",i));
        line.push_back(dat.gridt.get("Lambda_bar",i));
      }
    }
    
    if (nsd->source_fnames_alt.size()>0) {
      for(size_t i=0;i<nsd->n_sources;i++) {
        line.push_back(dat.sourcet.get("atm",i));
      }
    }

    if (set->inc_ligo) {
      line.push_back(dat.eos.get_constant("M_chirp"));
      line.push_back(dat.eos.get_constant("m1"));
      line.push_back(dat.eos.get_constant("m2"));
      line.push_back(dat.eos.get_constant("R1"));
      line.push_back(dat.eos.get_constant("R2"));
      line.push_back(dat.eos.get_constant("I1"));
      line.push_back(dat.eos.get_constant("I2"));
      line.push_back(dat.eos.get_constant("I_bar1"));
      line.push_back(dat.eos.get_constant("I_bar2"));
      line.push_back(dat.eos.get_constant("Lambda1"));
      line.push_back(dat.eos.get_constant("Lambda2"));
      line.push_back(dat.eos.get_constant("Lambdat"));
      line.push_back(dat.eos.get_constant("del_Lambdat"));    
      line.push_back(dat.eos.get_constant("ligo_prob"));
      line.push_back(dat.eos.get_constant("eta"));
    }
    
    if (nsd->n_sources>0) {
      for(size_t i=0;i<nsd->n_sources;i++) {
        if (dat.eos.is_constant(((std::string)"log_wgt_")+
                                nsd->source_names[i])){
          line.push_back(dat.eos.get_constant(((std::string)"log_wgt_")+
                                              nsd->source_names[i]));
          
        } else {
          line.push_back(-800);
        }
      }
    }
    
    if (set->use_population) {
      std::cout << "XZ: " << pop_weights.size() << std::endl;
      for (size_t i=0; i<pop_weights.size(); i++) {
        line.push_back(pop_weights[i]);
      }
    }
    return o2scl::success;
  }  
  
}

int bamr_class::compute_point(const ubvector &pars, std::ofstream &scr_out, 
			      double &log_wgt, model_data &dat) {

  log_wgt=0.0;

  int iret;

  if (set->apply_emu) {

    // create vector for emulator prediction
    ubvector test_pars;
    
    // copy mcmc param values
    test_pars = pars;
    
    // update emulator parameter vector with H or He atm values
    if (nsd->n_sources>0) {
      
      test_pars.resize(pars.size()+nsd->n_sources);
      
      for(size_t i=0; i<pars.size(); i++){
        test_pars[i] = pars[i];
      }
      /* 
         MCMC paprmeter vector contains moddel params and $mf_$'s
         from the sources. We calculate "atm" values from the $mf_$'s
         and pass the additional atm values to the "emy.py". "emu.py"
         was trained with "atm" columns with the mcmc_params. To 
         emulate a point we need to update the "atm" values.
      */      
      for(size_t i=(pars.size()-nsd->n_sources); i<pars.size(); i++){
        double atm=pars[i]*1.0e8-((double)((int)(pars[i]*1.0e8)));
        if(atm<2/3){
          test_pars[pars.size()] = 0;
        } else {
          test_pars[pars.size()] = 1;
        }
      }
    }

    // Create new pylist from param_vals
    test_vals = PyList_New(test_pars.size());
    for(size_t i=0; i<test_pars.size(); i++){
      PyList_SetItem(test_vals, i, PyFloat_FromDouble(test_pars[i]));
    }

    /* 
       As a test, call emu.py:modGpr:show().
    */

    if (PyCallable_Check(train_trainMthd)) {
      target_pred=PyObject_CallObject
        (train_trainMthd, 
         PyTuple_Pack(4,PyUnicode_FromString(set->emu_train.c_str()),
                      train_tParam_Names,test_vals,addtl_sources));
    }

    // Finally, set the value of log_wgt equal to the value returned
    // by the python emulator
    
    // Check current MPI rank
    int mpi_rank = 0;
#ifdef BAMR_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#endif
    
    // Check current OpenMP thread
    int pthread=0;
#ifdef O2SCL_OPENMP      
    pthread=omp_get_thread_num();
#endif

    // Prediction vector copied from the python list stored in
    // target_pred
    ubvector preds;
    preds.resize(PyList_Size(target_pred));
    
    for (long int i=0; i < PyList_Size(target_pred); i++) {
      PyObject *pTarget = PyList_GetItem(target_pred, i);
      preds[i] = PyFloat_AsDouble(pTarget);
    }
    
    log_wgt = preds[0];
    if (false) {
      cout << "Emulated log_wgt by rank "<< mpi_rank
           <<" and thread "<< pthread <<
        " : " << log_wgt << endl;
    }

    double pred_Mmax = preds[2];
    if (pred_Mmax < 2.0) {
      iret=1;
    }

    double pred_e_max=preds[3];

    // Check speed of sound causal limit
    for (size_t i=0;i<100;i++) {
      double e_i=mod->e_grid[i];
      if (e_i<pred_e_max) {
        if (preds[preds.size()-(i+1)] > 1.0) {
          iret=1;
        }
      }
    }

    iret = 0;
    
  } else {

    // Reference to model object for convenience
    model &m=*this->mod;
    
    // Compute the M vs R curve and return a non-zero value if it failed
    m.compute_star(pars,scr_out,iret,dat,model_type);
    if (iret!=0) {
      log_wgt=0.0;
      return iret;
    }

    // Calculate likelihood if using mass data from populations
    if (set->use_population) {
      
      likelihood &like = nsd->pop_like;

      std::cout << "XY: " << pop_weights.size() << endl;
      if (pop_weights.size()==0) pop_weights.resize(5); 
      pop_weights[0] = like.get_weight_ns(pars, pvi, iret);
      if (iret!=0) {
        log_wgt=0.0;
        return iret;
      }
      pop_weights[1] = like.get_weight_wd(pars, pvi, iret);
      if (iret!=0) {
        log_wgt=0.0;
        return iret;
      }
      pop_weights[2] = like.get_weight_hms(pars, pvi, iret);
      if (iret!=0) {
        log_wgt=0.0;
        return iret;
      }
      pop_weights[3] = like.get_weight_lms(pars, pvi, iret);
      if (iret!=0) {
        log_wgt=0.0;
        return iret;
      }

      for (int i=0; i<4; i++) pop_weights[4] += pop_weights[i];
      log_wgt += pop_weights[4];
      
      if (iret==0) {
        cout << "Final pop result: ";
        vector_out(cout, pop_weights, true);
      }

    }

    // ----------------------------------------------------------------
    // Exit early if the mass and radius for any of the masses or radii
    // are out of range
	  
    for(size_t i=0;i<nsd->n_sources;i++) {
      double mass=dat.sourcet.get("M",i);
      double rad=dat.sourcet.get("R",i);
      if (mass<set->in_m_min || mass>set->in_m_max || 
          rad<set->in_r_min || rad>set->in_r_max) {
        scr_out << "Rejected: Mass or radius outside range." << std::endl;
        scr_out << "M limits: " << set->in_m_min << " "
                << set->in_m_max << std::endl;
        scr_out << "R limits: " << set->in_r_min << " "
                << set->in_r_max << std::endl;
        if (nsd->n_sources>0) {
          scr_out.precision(2);
          scr_out.setf(ios::showpos);
          scr_out << "M ";
          for(size_t j=0;j<nsd->n_sources;j++) {
            scr_out << dat.sourcet.get("M",j) << " ";
          }
          scr_out << std::endl;
          scr_out << "R ";
          for(size_t j=0;j<nsd->n_sources;j++) {
            scr_out << dat.sourcet.get("R",j) << " ";
          }
          scr_out << std::endl;
          scr_out.precision(6);
          scr_out.unsetf(ios::showpos);
        }
        log_wgt=0.0;
        return m.ix_mr_outside;
      }
    }

    // -----------------------------------------------
    // Determine the atm parameter
	  
    
    for (size_t i=0;i<nsd->n_sources;i++) {
	    
      // Determine H or He from mass parameter
      double mf;
      if (set->inc_ligo) {
        mf=pars[i+mod->n_eos_params+3];
      } else {
        mf=pars[i+mod->n_eos_params];
      }
      double d_atm=mf*1.0e8-((double)((int)(mf*1.0e8)));
      if (d_atm<2.0/3.0) {
        dat.sourcet.set("atm",i,0.0);
      } else {
        dat.sourcet.set("atm",i,1.0);
      }
    }

    if (set->apply_intsc==false) {

      // -----------------------------------------------
      // Compute the weights for each source
	    
      dat.mvsr.set_interp_type(o2scl::itp_linear);
	    
      double m_max_current=dat.mvsr.max("gm");
	    
      if (set->verbose>=2) scr_out << "Name M R Weight" << std::endl;
	    
      for(size_t i=0;i<nsd->n_sources;i++) {
	      
        double mass=dat.sourcet.get("M",i);
        double rad=dat.sourcet.get("R",i);
        bool atm=false;
        if (dat.sourcet.get("atm",i)>0.5) atm=true;
	      
        // Double check that current M and R is in the range of
        // the provided input data
        if (rad<nsd->source_tables[i].get_x_data()[0] ||
            rad>nsd->source_tables[i].get_x_data()
            [nsd->source_tables[i].get_nx()-1] ||
            mass<nsd->source_tables[i].get_y_data()[0] ||
            mass>nsd->source_tables[i].get_y_data()
            [nsd->source_tables[i].get_ny()-1]) {
		
          dat.sourcet.set("wgt",i,0.0);
		
        } else {
		
          // If M and R are in range, compute the weight
		
          if (nsd->source_fnames_alt.size()>0) {
		  
            // Compute alternate probability from an insignificant bit
            // in the mass 
		  
            if (atm==false) {
              dat.sourcet.set("wgt",i,
                              nsd->source_tables[i].interp
                              (rad,mass,nsd->slice_names[i]));
            } else {
              dat.sourcet.set("wgt",i,
                              nsd->source_tables_alt[i].interp
                              (rad,mass,nsd->slice_names[i]));
            }
		  
          } else {
            dat.sourcet.set("wgt",i,
                            nsd->source_tables[i].interp
                            (rad,mass,nsd->slice_names[i]));
          }
		
          // If the weight is lower than the threshold, set it equal
          // to the threshold
          if (dat.sourcet.get("wgt",i)<set->input_dist_thresh) {
            dat.sourcet.set("wgt",i,set->input_dist_thresh);
          }
		
        }
	      
        // If the data gives a zero weight, just return a factor
        // of 1e8 smaller than the peak value
        if (dat.sourcet.get("wgt",i)<=0.0) {
          dat.sourcet.set("wgt",i,
                          o2scl::matrix_max_value<ubmatrix,double>
                          (nsd->source_tables[i].get_slice
                           (nsd->slice_names[i]))/1.0e8);
        }
	      
        // If the weight is zero, then return failure
        if (dat.sourcet.get("wgt",i)<=0.0) {
          scr_out << "Weight zero for source " << i << " "
                  << nsd->source_names[i] << " with mass " << mass
                  << " and radius " << rad
                  << " with atm=" << atm << endl;
          return m.ix_mr_outside;
        }
	      
        // Include the weight for this source
        log_wgt+=log(dat.sourcet.get("wgt",i));

        // Update each weight into output table
        dat.eos.add_constant(((std::string)"log_wgt_")+nsd->source_names[i]
          ,log(dat.sourcet.get("wgt",i)));
	      
        if (set->verbose>=2) {
          scr_out.width(10);
          scr_out << nsd->source_names[i] << " "
                  << mass << " " 
                  << rad << " " << dat.sourcet.get("wgt",i) << std::endl;
        }
	      
        // Go to the next source
      }

      if (set->debug_star) scr_out << std::endl;
	    
      // -----------------------------------------------
      // Exit if the current maximum mass is too large
	    
      if (m_max_current>set->exit_mass) {
        scr_out.setf(ios::scientific);
        scr_out << "Exiting because maximum mass (" << m_max_current 
                << ") larger than exit_mass (" << set->exit_mass << ")." 
                << std::endl;
        scr_out.precision(12);
        vector_out(scr_out,pars);
        scr_out << " " << log_wgt << std::endl;
        scr_out.precision(6);
        exit(0);
      }
	    
      if (set->verbose>=2) {
        cout << "End model::compute_point()." << endl;
      }
	    
      if (iret!=m.ix_success) {
        // We shouldn't be returning a non-zero value if success is
        // non-zero, so we double check this here
        O2SCL_ERR("Sanity check for success flag in model::compute_point.",
                  o2scl::exc_esanity);
      }

	    
    } else {

      // Apply intrinsic scatter

      // Input and output table references for convenience
      std::vector<o2scl::table3d> &in=nsd->source_tables;
      std::vector<o2scl::table3d> &in_alt=nsd->source_tables_alt;
      std::vector<o2scl::table3d> &out=source_tables_is;
      std::vector<o2scl::table3d> &out_alt=source_tables_alt_is;

      for (size_t i=0;i<nsd->n_sources;i++) {

        bool atm=false;
        if (dat.sourcet.get("atm",i)>0.5) atm=true;
	      
        int ithread=0;
#ifdef O2SCL_OPENMP      
        ithread=omp_get_thread_num();
        if (omp_get_num_threads()!=n_threads) {
          O2SCL_ERR("Thread count doesn't match.",o2scl::exc_esanity);
        }
#endif
        int fix=ithread*nsd->n_sources+i;

        bool make_tensor_files=false;
        if (make_tensor_files) {
          // Create the tensor files for 'cached_intsc'
	        
          for (size_t iik=0;iik<2;iik++) {
	          
            tensor_grid3<> tg3;
            size_t sz[3]={in[i].get_nx(),in[i].get_ny(),21};
            tg3.resize(3,sz);
            vector<double> gp;
            for(size_t i3=0;i3<in[i].get_nx();i3++) {
              gp.push_back(in[i].get_grid_x(i3));
            }
            for(size_t i3=0;i3<in[i].get_ny();i3++) {
              gp.push_back(in[i].get_grid_y(i3));
            }
            for (double iix=-2.0;iix<2.01;iix+=0.2) {
              gp.push_back(iix);
            }
            tg3.set_grid_packed(gp);
	          
            size_t iiz=0;
	          
#ifdef BAMR_FFTW3
            for (double iix=-2.0;iix<2.01;iix+=0.2) {
	            
              double intrsc=pow(10.0,iix);
              double dx = in[i].get_grid_x(1)-in[i].get_grid_x(0);
              double dy = in[i].get_grid_y(1)-in[i].get_grid_y(0);
	            
              // radius direction
              double sigx = intrsc/dx*5.0;
              // mass direction
              double sigy = intrsc/dy*1.4;  
	            
              cout << "Filtering " << nsd->source_names[i]
                   << " log10(intsig): "
                   << pars[i+mod->n_eos_params+nsd->n_sources]
                   << " with intsig: " << intrsc
                   << "\n  sigx: " << sigx << " sigy: " << sigy << std::endl;

              // Set kernel
              flt[fix]->init_gaussian_kernel(sigx,sigy);
              flt[fix]->fft_kernel(); 
	            
              // Set input image
              if (iik==0) {
                flt[fix]->set_image(in[i],nsd->slice_names[i]);
              } else {
                flt[fix]->set_image(in_alt[i],nsd->slice_names[i]);
              }
	            
              // Convolve
              flt[fix]->fft_image_forward();
              flt[fix]->apply_kernel();
              flt[fix]->fft_image_backward();
	            
              // Read image back
              if (iik==0) {
                flt[fix]->get_image(out[fix],nsd->slice_names[i]);
              } else {
                flt[fix]->get_image(out_alt[fix],nsd->slice_names[i]);
              }
	            
              if (iik==0) {
                ubmatrix &outs=out[fix].get_slice(nsd->slice_names[i]);
                cout << "outs: " << intrsc << " " << iik << " "
                     << outs(0,0) << endl;
                for(size_t i3=0;i3<in[i].get_nx();i3++) {
                  for(size_t j3=0;j3<in[i].get_nx();j3++) {
                    tg3.set(i3,j3,iiz,outs(i3,j3));
                  }
                }
              } else {
                ubmatrix &outs=out_alt[fix].get_slice(nsd->slice_names[i]);
                cout << "outs: " << intrsc << " " << iik << " "
                     << outs(0,0) << endl;
                for(size_t i3=0;i3<in[i].get_nx();i3++) {
                  for(size_t j3=0;j3<in[i].get_nx();j3++) {
                    tg3.set(i3,j3,iiz,outs(i3,j3));
                  }
                }
              }

              iiz++;
            }
#endif
		  
            hdf_file hfx;
            string fnamex=set->data_dir+"/cache/tg_"+o2scl::szttos(i)+
              "_"+o2scl::szttos(iik);
            hfx.open_or_create(fnamex);
            hdf_output(hfx,tg3,"tg");
            hfx.close();
          }
	        
          if (i==nsd->n_sources-1) {
            exit(-1);
          }
          // End of 'if (make_tensor_files)'
        }
	     
        // Calculate 2D intrinsic scatter
        // value is normalized with (12km / 1.5Msun)
        // NOTE: implicitly assumes uniform grid
        double intrsc = pow(10.0,pars[i+mod->n_eos_params+nsd->n_sources]);
        double dx = in[i].get_grid_x(1)-in[i].get_grid_x(0);
        double dy = in[i].get_grid_y(1)-in[i].get_grid_y(0);

        // radius direction
        double sigx = intrsc/dx*5.0;
        // mass direction
        double sigy = intrsc/dy*1.4;  

        if (set->verbose>1) {
          scr_out << "Filtering " << nsd->source_names[i]
                  << " log10(intsig): "
                  << pars[i+mod->n_eos_params+nsd->n_sources]
                  << " with intsig: " << intrsc
                  << "\n  sigx: " << sigx << " sigy: " << sigy << std::endl;
        }

#ifdef BAMR_FFTW3
        if (set->cached_intsc==false) {
          // Set kernel
          flt[fix]->init_gaussian_kernel(sigx,sigy);
          flt[fix]->fft_kernel(); 
	        
          // Set input image
          if (atm==false) {
            flt[fix]->set_image(in[i],nsd->slice_names[i]);
          } else {
            flt[fix]->set_image(in_alt[i],nsd->slice_names[i]);
          }
	        
          // Convolve
          flt[fix]->fft_image_forward();
          flt[fix]->apply_kernel();
          flt[fix]->fft_image_backward();
	        
          // Read image back
          if (atm==false) {
            flt[fix]->get_image(out[fix],nsd->slice_names[i]);
          } else {
            flt[fix]->get_image(out_alt[fix],nsd->slice_names[i]);
          }

          if (true) {
            if (atm==false) {
              ubmatrix &outs=out[fix].get_slice(nsd->slice_names[i]);
              double sum=matrix_sum<ubmatrix,double>(outs);
              // Renormalize
              for(size_t j=0;j<outs.size1();j++) {
                for(size_t k=0;k<outs.size2();k++) {
                  outs(j,i)/=sum;
                }
              }
            } else {
              ubmatrix &outs=out_alt[fix].get_slice(nsd->slice_names[i]);
              double sum=matrix_sum<ubmatrix,double>(outs);
              // Renormalize
              for(size_t j=0;j<outs.size1();j++) {
                for(size_t k=0;k<outs.size2();k++) {
                  outs(j,k)/=sum;
                }
              }
            }
          }
	      
          // End of 'if (cached_intsc==false)'
        }
#endif

        // End of loop over sources
      }

      for(size_t i=0;i<nsd->n_sources;i++) {
	      
        double mass=dat.sourcet.get("M",i);
        double rad=dat.sourcet.get("R",i);
	      
        if (mass<set->in_m_min || mass>set->in_m_max || 
            rad<set->in_r_min || rad>set->in_r_max) {
		
          scr_out << "Rejected: Mass or radius outside range." << std::endl;
          scr_out << "M limits: " << set->in_m_min << " "
                  << set->in_m_max << std::endl;
          scr_out << "R limits: " << set->in_r_min << " "
                  << set->in_r_max << std::endl;
          if (nsd->n_sources>0) {
            scr_out.precision(2);
            scr_out.setf(ios::showpos);
            scr_out << "M ";
            for(size_t j=0;j<nsd->n_sources;j++) {
              scr_out << dat.sourcet.get("M",j) << " ";
            }
            scr_out << std::endl;
            scr_out << "R ";
            for(size_t j=0;j<nsd->n_sources;j++) {
              scr_out << dat.sourcet.get("R",j) << " ";
            }
            scr_out << std::endl;
            scr_out.precision(6);
            scr_out.unsetf(ios::showpos);
          }
          log_wgt=0.0;
          iret=m.ix_mr_outside;
        }
      }

      if (iret==0) {
	      
        log_wgt=0.0;
	      
        dat.mvsr.set_interp_type(o2scl::itp_linear);
        double m_max_current=dat.mvsr.max("gm");
	      
        // -----------------------------------------------
        // Compute the weights for each source
	      
        if (set->verbose>=2) scr_out << "Name M R Weight" << std::endl;
	      
        for(size_t i=0;i<nsd->n_sources;i++) {

          double mass=dat.sourcet.get("M",i);
          double rad=dat.sourcet.get("R",i);
          bool atm=false;
          if (dat.sourcet.get("atm",i)>0.5) atm=true;
		
          int ithread=0;
#ifdef O2SCL_OPENMP      
          ithread=omp_get_thread_num();
          if (omp_get_num_threads()!=n_threads) {
            O2SCL_ERR("Thread count doesn't match.",o2scl::exc_esanity);
          }
#endif
          int fix=ithread*nsd->n_sources+i;
	        
          // Double check that current M and R is in the range of
          // the provided input data
          if (rad<out[fix].get_x_data()[0] ||
              rad>out[fix].get_x_data()[out[fix].get_nx()-1] ||
              mass<out[fix].get_y_data()[0] ||
              mass>out[fix].get_y_data()[out[fix].get_ny()-1]) {

            dat.sourcet.set("wgt",i,0.0);
	          
          } else {

            if (atm==false) {
              if (set->cached_intsc) {

                ubvector iu(3);
                iu[0]=rad;
                iu[1]=mass;
                iu[2]=pars[i+mod->n_eos_params+nsd->n_sources];
                dat.sourcet.set("wgt",i,fft_data[i*2].interp_linear(iu));
              } else {
                dat.sourcet.set("wgt",i,
                                out[fix].interp
                                (rad,mass,nsd->slice_names[i]));
              }
            } else {
              if (set->cached_intsc) {
                ubvector iu(3);
		      
                iu[0]=rad;
                iu[1]=mass;
                iu[2]=pars[i+mod->n_eos_params+nsd->n_sources];
                dat.sourcet.set("wgt",i,
                                fft_data[(i*2)+1].interp_linear(iu));
              } else {
                dat.sourcet.set("wgt",i,
                                out_alt[fix].interp
                                (rad,mass,nsd->slice_names[i]));
              }
            }
	          
            // If the weight is lower than the threshold, set it equal
            // to the threshold
            if (dat.sourcet.get("wgt",i)<set->input_dist_thresh) {
              dat.sourcet.set("wgt",i,set->input_dist_thresh);
            }
	          
          }
	        
          // If the data gives a zero weight, just return a factor
          // of 1e8 smaller than the peak value
          if (dat.sourcet.get("wgt",i)<=0.0) {
            typedef boost::numeric::ublas::matrix<double> ubmatrix;
            dat.sourcet.set("wgt",i,
                            o2scl::matrix_max_value<ubmatrix,double>
                            (out[fix].get_slice(nsd->slice_names[i]))/1.0e8);
          }
		
          // If the weight is zero, then return failure
          if (dat.sourcet.get("wgt",i)<=0.0) {
            scr_out << "Weight zero for source " << i << " "
                    << nsd->source_names[i]
                    << " with mass " << mass << " and radius "
                    << rad << " with atm=" << atm << endl;
            iret=m.ix_mr_outside;
          }
	        
          // Include the weight for this source
          log_wgt+=log(dat.sourcet.get("wgt",i));

        // Update each weight into output table
        dat.eos.add_constant(((std::string)"log_wgt_")+nsd->source_names[i]
          ,log(dat.sourcet.get("wgt",i)));
	        
          if (set->verbose>=2) {
            scr_out.width(10);
            scr_out << nsd->source_names[i] << " "
                    << mass << " " 
                    << rad << " " << dat.sourcet.get("wgt",i) << std::endl;
          }
	        
          // Go to the next source
        }

        if (!std::isfinite(log_wgt)) {
          scr_out << "IS weight not finite." << endl;
          iret=m.ix_mr_outside;
          log_wgt=0.0;
        }
	      
      }
	    
      if (set->debug_star) scr_out << std::endl;

      // End of 'if (apply_intsc)'
    }
    
    // Add Tews et al. probability to the log likelihood
    if (iret==0 && (model_type==((string)"tews_threep_ligo") ||
                    model_type==((string)"tews_fixp_ligo"))) {
      log_wgt+=dat.eos.get_constant("tews");
    }
    
    // Section for additional LIGO constraints 
    if (iret==0 && set->inc_ligo) {
            
      double M_chirp_det=0.0, q=0.0, z_cdf; 
      double M_chirp, z, m1=0.0, m2=0.0;
      M_chirp_det=pars[m.n_eos_params];
      q=pars[m.n_eos_params+1];
      z_cdf=pars[m.n_eos_params+2];
      prob_dens_gaussian pdg(0.0099,0.0009);
      z=pdg.invert_cdf(z_cdf);
      M_chirp=M_chirp_det/(1.0+z);
      dat.eos.add_constant("M_chirp",M_chirp);
      m1=M_chirp*pow(1.0+q,0.2)/pow(q,0.6);
      m2=M_chirp*pow(q,0.4)*pow(1.0+q,0.2);

      double Mmax=dat.mvsr.max("gm");
      
      if (m1>Mmax || m2>Mmax || m1<m2) {
        
        log_wgt=0.0;
        iret=1;
        
      } else {
        
        // radii
        double R1=dat.mvsr.interp("gm",m1,"r");
        double R2=dat.mvsr.interp("gm",m2,"r");
        // I's in Msun*km^2
        double I1=dat.mvsr.interp("gm",m1,"rjw")/3.0/schwarz_km;
        double I2=dat.mvsr.interp("gm",m2,"rjw")/3.0/schwarz_km;
        // To compute I_bar, divide by G^2*M^3
        double G=schwarz_km/2.0;
        double I_bar1=I1/G/G/m1/m1/m1;
        double I_bar2=I2/G/G/m2/m2/m2;
        dat.eos.add_constant("m1",m1);
        dat.eos.add_constant("m2",m2);
        dat.eos.add_constant("R1",R1);
        dat.eos.add_constant("R2",R2);
        dat.eos.add_constant("I1",I1);
        dat.eos.add_constant("I2",I2);
        dat.eos.add_constant("I_bar1",I_bar1);
        dat.eos.add_constant("I_bar2",I_bar2);

        // Jim's fit from Steiner, Lattimer, and Brown (2016)
        double b0=-30.5395;
        double b1=38.3931;
        double b2=-16.3071;
        double b3=3.36972;
        double b4=-0.26105;
            
        double li=log(I_bar1);
        double li2=li*li;
        double li3=li*li2;
        double li4=li*li3;
        double li5=li*li4;
        double li6=li*li5;

        double Lambda1=exp(b0+b1*li+b2*li2+b3*li3+b4*li4);

        li=log(I_bar2);
        li2=li*li;
        li3=li*li2;
        li4=li*li3;
        li5=li*li4;
        li6=li*li5;
        double Lambda2=exp(b0+b1*li+b2*li2+b3*li3+b4*li4);

        dat.eos.add_constant("Lambda1",Lambda1);
        dat.eos.add_constant("Lambda2",Lambda2);
        double Lambdat=16.0/13.0*((m1+12.0*m2)*pow(m1,4.0)*Lambda1+
                                  (m2+12.0*m1)*pow(m2,4.0)*Lambda2)/
          pow(m1+m2,5.0);
        dat.eos.add_constant("Lambdat",Lambdat);

        double eta=(m1*m2)/((m1+m2)*(m1+m2));
        dat.eos.add_constant("eta",eta);
        
        double eta2=eta*eta, eta3=eta2*eta;
        double del_Lambdat=0.5*
          ((sqrt(1.0-4.0*eta)*(1.0-(13272.0/1319.0)*eta+
                               (8944.0/1319.0)*eta2)*(Lambda1+Lambda2))+
           ((1.0-(15910.0/1319.0)*eta+(32850.0/1319.0)*eta2+
             (3380.0/1319.0)*eta3)*(Lambda1-Lambda2)));
        dat.eos.add_constant("del_Lambdat",del_Lambdat);
        
        // Add the LIGO log-likelihood
        double prob_data=-800.0;

        // Check the ligo_data_table for new or old data
        if (true) {
                
          ubvector lin_v(3);
          lin_v[0]=M_chirp_det;
          lin_v[1]=q;
          lin_v[2]=Lambdat;

          double prob=nsd->ligo_data_table.interp_linear(lin_v);
          // If the point is outside of the range specified
          // in the data file, give it a very small probability
          for(size_t jj=0;jj<3;jj++) {
            if (lin_v[jj]<nsd->ligo_data_table.get_grid(jj,0) ||
                lin_v[jj]>nsd->ligo_data_table.get_grid
                (jj,nsd->ligo_data_table.get_size(jj)-1)) {
              scr_out << "LIGO quantity " << jj
                      << " out of range: " << endl;
              size_t n_ligo=nsd->ligo_data_table.get_size(jj);
              scr_out << lin_v[jj] << " "
                      << nsd->ligo_data_table.get_grid(jj,0) << " "
                      << nsd->ligo_data_table.get_grid(jj,n_ligo-1) 
                      << endl;
              prob=-800.0;
            }
          }
          prob_data=prob;
          
        }

        dat.eos.add_constant("ligo_prob",prob_data);
        log_wgt+=(prob_data);
      }
      
      // End of section for additional LIGO constraints
    }


    // If the gridt table has not yet been initialized perform that
    // initialization

    if (dat.gridt.get_ncolumns()==0) {
      
      dat.gridt.set_nlines(set->grid_size);
      dat.gridt.new_column("m_grid");
      for(int i=0;i<set->grid_size;i++) {
        dat.gridt.set("m_grid",i,m.m_grid[i]);
      }
      
      dat.gridt.new_column("R");
      if (m.has_eos) {
        dat.gridt.new_column("e_grid");
        for(int i=0;i<set->grid_size;i++) {
          dat.gridt.set("e_grid",i,m.e_grid[i]);
        }
        dat.gridt.new_column("P");
        dat.gridt.new_column("cs2");
        dat.gridt.new_column("PM");
      }

      if (set->baryon_density) {
        dat.gridt.line_of_names("nb_grid Pnb EoA");
        for(int i=0;i<set->grid_size;i++) {
          dat.gridt.set("nb_grid",i,m.nb_grid[i]);
        }
      }

      if (set->compute_cthick) {
        dat.gridt.new_column("CT");
      }
      if (set->addl_quants) {
        dat.gridt.line_of_names("MB BE I I_bar Lambda_bar");
      }
    }

    size_t n_params=pars.size();
	  
    double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

    if (m.has_eos) {

      // The central energy density in the maximum mass configuration
      emax=dat.mvsr.max("ed");
      // The central pressure in the maximum mass configuration
      pmax=dat.mvsr.max("pr");

      dat.mvsr.add_constant("P_max",pmax);
      dat.mvsr.add_constant("e_max",emax);

      // The maximum mass
      mmax=dat.mvsr.get_constant("M_max");
      // The radius of the maximum mass star
      rmax=dat.mvsr.get_constant("R_max");

      if (set->baryon_density) {

        // The highest baryon density in the EOS table
        nbmax2=dat.eos.max("nb");
	      
        // The central baryon density in the maximum mass configuration
        nbmax=dat.mvsr.get_constant("nb_max");
	      
        dat.mvsr.add_constant("nb_max",nbmax);
	      
      }

    } else {
      // Need to set mmax for no EOS models to figure out how 
      // high up we should go for the radius grid 
      mmax=3.0;
    }
	  
    if (m.has_eos) {
      for(int i=0;i<set->grid_size;i++) {
        double eval=m.e_grid[i];
        // Make sure the energy density from the grid isn't beyond the
        // last energy density computed by the EOS model (but still
        // include energy densities larger than the maximum energy
        // density of the maximum mass star)
        double emax2=dat.eos.max("ed");
        if (eval<emax2) {
          double pres_temp=dat.eos.interp("ed",eval,"pr");
          double cs2_temp=dat.eos.interp("ed",eval,"cs2");
          dat.gridt.set("P",i,pres_temp);
          dat.gridt.set("cs2",i,cs2_temp);
        } else {
          dat.gridt.set("P",i,0.0);
          dat.gridt.set("cs2",i,0.0);
        }
      }
    }

    // It is important here that all of these columns which store values
    // over a grid are either always positive or always negative,
    // because the code reports zero in the fill_line() function for
    // values beyond the end of the EOS or the M-R curve. 
    for(int i=0;i<set->grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax) {
        dat.gridt.set("R",i,dat.mvsr.interp("gm",mval,"r"));
        if (m.has_eos) {
          dat.gridt.set("PM",i,dat.mvsr.interp("gm",mval,"pr"));
        }
      } else {
        dat.gridt.set("R",i,0.0);
        if (m.has_eos) {
          dat.gridt.set("PM",i,0.0);
        }
      }
    }
    if (m.has_eos) {
      if (set->baryon_density) {
        for(int i=0;i<set->grid_size;i++) {
          double nbval=m.nb_grid[i];
          if (nbval<nbmax2) {
            double pres_temp=dat.eos.interp("nb",nbval,"pr");
            if (pres_temp<pmax) {
              dat.gridt.set("Pnb",i,pres_temp);
            } else {
              dat.gridt.set("Pnb",i,0.0);
            }
            double eval2=dat.eos.interp("nb",nbval,"ed");
            double eoa_val2=eval2/nbval-939.0/o2scl_const::hc_mev_fm;
            dat.gridt.set("EoA",i,eoa_val2);
          } else {
            dat.gridt.set("Pnb",i,0.0);
            dat.gridt.set("EoA",i,0.0);
          }
        }
      }

      for(size_t i=0;i<nsd->n_sources;i++) {
        double val=dat.mvsr.interp
          ("gm",pars[n_params-nsd->n_sources+i],"ed");
        dat.sourcet.set("ce",i,val);
      }
      if (set->baryon_density) {
        for(size_t i=0;i<nsd->n_sources;i++) {
          double val2=dat.mvsr.interp
            ("gm",pars[n_params-nsd->n_sources+i],"nb");
          dat.sourcet.set("cnb",i,val2);
        }
      }
    }
	  
    if (set->compute_cthick) {
      for(int i=0;i<set->grid_size;i++) {
        double mval=m.m_grid[i];
        if (mval<mmax) {
          double rval=dat.mvsr.interp("gm",mval,"r");
          // Compute the crust thickness by subtracting the radius
          // of the crust core transition from the full radius
          dat.gridt.set("CT",i,rval-dat.mvsr.interp("gm",mval,"r0"));
        } else {
          dat.gridt.set("CT",i,0.0);
        }
      }
    }
	  
    if (set->addl_quants) {
      for(int i=0;i<set->grid_size;i++) {
        double mval=m.m_grid[i];

        if (mval<mmax) {
		
          // Baryonic mass
          double bm=dat.mvsr.interp("gm",mval,"bm");
          dat.gridt.set("MB",i,bm);

          // Binding energy
          dat.gridt.set("BE",i,bm-mval);

          // Moment of inertia
          double rad=dat.mvsr.interp("gm",mval,"r");
          // rjw is km^4, so dividing by km/Msun gives Msun*km^2
          double I=dat.mvsr.interp("gm",mval,"rjw")/3.0/schwarz_km;
          dat.gridt.set("I",i,I);

          // To compute I_bar, divide by G^2*M^3
          double I_bar=I*4.0/schwarz_km/schwarz_km/mval/mval/mval;
          dat.gridt.set("I_bar",i,I_bar);

          // Jim's fit from Steiner, Lattimer, and Brown (2016)
          double b0=-30.5395;
          double b1=38.3931;
          double b2=-16.3071;
          double b3=3.36972;
          double b4=-0.26105;
	    
          double li=log(I_bar);
          double li2=li*li;
          double li3=li*li2;
          double li4=li*li3;
		
          double Lambda_bar=exp(b0+b1*li+b2*li2+b3*li3+b4*li4);

          dat.gridt.set("Lambda_bar",i,Lambda_bar);
		
        } else {

          dat.gridt.set("MB",i,0.0);
          dat.gridt.set("BE",i,0.0);
          dat.gridt.set("I",i,0.0);
          dat.gridt.set("I_bar",i,0.0);
          dat.gridt.set("Lambda_bar",i,0.0);
		
        }
      }
	  
    }
  }
  return iret;
}

