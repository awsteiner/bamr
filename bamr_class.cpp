/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2026, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
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

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace bamr;
using namespace std::placeholders;

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

  model &m=*this->mod;
  
  if (!set->emu_tov) {

    for(size_t i=0;i<nsd->n_sources;i++) {
      line.push_back(dat.sourcet.get("R",i));
    }
    for(size_t i=0;i<nsd->n_sources;i++) {
      line.push_back(dat.sourcet.get("M",i));
    }

    if (m.has_eos) {
      for(int i=0;i<set->grid_size;i++) {
        line.push_back(dat.gridt.get("P",i));
        line.push_back(dat.gridt.get("cs2",i));
      }
    }

    for(int i=0;i<set->grid_size;i++) {
      line.push_back(dat.gridt.get("R",i));
      if (m.has_eos) {
        line.push_back(dat.gridt.get("PM",i));
      }
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
        line.push_back(dat.eos.get_constant("M_max2"));
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
      line.push_back(dat.eos.get_constant("M_chirp_gw17"));
      line.push_back(dat.eos.get_constant("m1_gw17"));
      line.push_back(dat.eos.get_constant("m2_gw17"));
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
      line.push_back(dat.eos.get_constant("eta"));
      line.push_back(dat.eos.get_constant("m2_gw19"));
      line.push_back(dat.eos.get_constant("log_wgt_gw17"));
      line.push_back(dat.eos.get_constant("log_wgt_gw19"));
      if (set->inc_pop) {
        line.push_back(dat.eos.get_constant("log_SN_gw17"));
        line.push_back(dat.eos.get_constant("log_SN_gw19"));
      }
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
      if (set->inc_pop) {
        for (size_t i=0; i<nsd->n_sources; i++) {
          if (nsd->source_names[i]!=string("0030")) {
            line.push_back(dat.eos.get_constant(std::string("log_SN_")+
                                                nsd->source_names[i]));
          }
        }
      }
    }
    
    if (set->inc_pop) {
      for (size_t i=0; i<wgt_pop.size(); i++) {
        line.push_back(wgt_pop[i]);
      }
    }

    if (m.has_eos) {
      if (set->mmax_deriv) {
        line.push_back(dat.eos.get_constant("log_dpdM"));
      }
    }
    
  } 
  
  return o2scl::success;
  
}

int bamr_class::compute_point(const ubvector &pars, std::ofstream &scr_out, 
                              double &log_wgt, model_data &dat) {

  log_wgt=0.0;
  int iret;

  // Reference to model object for convenience
  model &m=*this->mod;

  if (!set->emu_tov) {
    
    // Compute the M vs R curve and return a non-zero value if it failed
    m.compute_star(pars,scr_out,iret,dat,model_type);
    
    if (iret!=m.ix_success) {
      if (set->verbose>=2) {
        cout << "models::compute_star() failure:"
             << " ix_return=" << iret << endl;
      }
      log_wgt=0.0;
      return iret;
    }

    // If likelihood is also a function of M_max, multiply by dpdM
    if (set->mmax_deriv && set->model_dpdm) {
      log_wgt+=dat.eos.get_constant("log_dpdM");
    }
    
    // -----------------------------------------------
    // Determine the atm parameter

    if (mcmc_method!=string("hmc")) {
      for (size_t i=0; i<nsd->n_sources; i++) {
        // Determine H or He from mass parameter
        double mf;
        if (set->inc_ligo) {
          mf=pars[i+mod->n_eos_params+nsd->n_ligo_params];
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
    }

    if (wgt_em.size()!=nsd->n_sources) wgt_em.resize(nsd->n_sources);
    if (fsn_em.size()!=nsd->n_sources) fsn_em.resize(nsd->n_sources);

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
              
        // Include the weight for this source
        wgt_em[i]=dat.sourcet.get("wgt",i);

        /* If population is included, calculate the skewed normal (SN) 
        PDF for the sources: QLMXBs, PREs, and NICER */
        if (set->inc_pop) {
          if (nsd->source_names[i]!=string("0030")) {
            ns_pop &pop=nsd->pop;
            double mean=pars[pvi["mean_LMS"]];
            double width=pow(10.0, pars[pvi["log10_width_LMS"]]);
            double skew=pars[pvi["skewness_LMS"]];

            double mf;
            if (set->inc_ligo) {
              mf=pars[i+mod->n_eos_params+nsd->n_ligo_params];
            } else {
              mf=pars[i+mod->n_eos_params];
            }
          
            double m_em=mf*m_max_current;
            fsn_em[i]=pop.skewed_norm(m_em,mean,width,skew);
            dat.eos.add_constant(string("log_SN_")+nsd->source_names[i],
                                 log(fsn_em[i]));
          } else {
            fsn_em[i]=1.0;
          }
        }

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
        cout << "bamr_class::compute_point() exited:"
             << " M_max > exit_mass" << endl;
        exit(0);
      }
            
      if (iret!=m.ix_success) {
        // We shouldn't be returning a non-zero value if success is
        // non-zero, so we double check this here
        O2SCL_ERR("Sanity check for success flag in model::compute_point.",
                  o2scl::exc_esanity);
      }
            
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

  if (iret==0 && set->verbose>=1) {
    cout << "bamr_class::compute_point() success:"
         << " log_wgt=" << log_wgt << endl;
  }

  return iret;
}

int bamr_class::compute_point_ext(const ubvector &pars, std::ofstream &scr_out, 
                                  double &log_wgt, model_data &dat) {
  int ret=compute_point(pars, scr_out, log_wgt, dat);
  if (ret!=0) log_wgt=-800.0-((double)ret);
  return 0;  
}


void bamr_class::compute_atms(const ubvector &pars, model_data &dat) {
  return;
}


int bamr_class::compute_gw17(const ubvector &pars, double &log_wgt,
                             model_data &dat) {
  return o2scl::success;
}


int bamr_class::compute_gw19(const ubvector &pars, double &log_wgt, 
                             model_data &dat) {
  return o2scl::success;
}


int bamr_class::compute_ems(size_t ix, const ubvector &pars, 
                           double &log_wgt, model_data &dat) {
  return o2scl::success;
}


int bamr_class::compute_dist(size_t ix, const ubvector &pars, 
                             double &log_wgt, model_data &dat) {
  return o2scl::success;
}


int bamr_class::numeric_deriv(size_t ix, const ubvector &x2, point_funct &pf,
                         double &pfx, double &g, model_data &dat) {
  return o2scl::success;
}


int bamr_class::compute_deriv(const ubvector &pars, point_funct &pf,
                              ubvector &grad, model_data &dat) {
  return 0;
}