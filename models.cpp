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
#include "models.h"

#include "bamr.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace bamr;

model::model(settings &s, ns_data &n) : set(s), nsd(n) {
  
  cns.nb_start=0.01;
  ts.verbose=0;
  ts.set_units("1/fm^4","1/fm^4","1/fm^3");
  ts.err_nonconv=false;
  ts.set_eos(teos);
      
  has_eos=true;
  schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;
      
  teos.verbose=0;

  // Transition density parameters
  nt_a.set_center(0.1327);
  nt_a.set_sigma(0.0024);
  nt_b.set_center(-0.0898);
  nt_b.set_sigma(0.0055);
  nt_c.set_center(0.0228);
  nt_c.set_sigma(0.0028);

  // If crust_from_L is true and the transition density is outside
  // this range, then we reject and go to the next point
  nt_low=0.02;
  nt_high=0.14;

  // -----------------------------------------------------------
  // Prepare crust

  if (set.use_crust) {
    teos.default_low_dens_eos();
    
    // Get the transition density from the crust
    double pt, pw;
    teos.get_transition(pt,pw);
    // We set the transition density a bit lower (because by default
    // it's the largest pressure in the crust EOS) and then add a
    // small width. This ensures the pressure will be increasing and
    // prefers the core EOS over the crust EOS near the transition.
    teos.transition_mode=eos_tov_interp::smooth_trans;
    teos.set_transition(pt/1.2,1.2);
    
  } else {
    teos.no_low_dens_eos();
  }
  
}

void ns_data::load_mc(std::ofstream &scr_out, int mpi_nprocs, int mpi_rank,
		      settings &set) {
      
  double tot, max;
      
  std::string name;

  if (source_names.size()!=source_fnames.size() ||
      source_names.size()!=init_mass_fracs.size()) {
    O2SCL_ERR("Incorrect input data sizes.",o2scl::exc_esanity);
  }
  
  if (n_sources>0) {
    
    if (set.verbose>=2) {
      cout << "bamr: Loading " << n_sources << " data files." << endl;
    }
    
    source_tables.resize(n_sources);
    
#ifdef BAMR_MPI_LOAD

    bool mpi_load_debug=false;
    int buffer=0, tag=0;
    
    // Choose which file to read first for this rank
    int filestart=0;
    if (mpi_rank>mpi_nprocs-((int) n_sources) && mpi_rank>0) {
      filestart=mpi_nprocs-mpi_rank;
    }
    if (mpi_load_debug) {
      scr_out << "Variable 'filestart' is " << filestart << " for rank "
	      << mpi_rank << "." << std::endl;
    }
    
    // Loop through all files
    for(int k=0;k<((int) n_sources);k++) {
      
      // For k=0, we choose some ranks to begin reading, the others
      // have to wait. For k>=1, all ranks have to wait their turn.
      if (k>0 || (mpi_rank>0 && mpi_rank<=mpi_nprocs-((int) n_sources))) {
	int prev=mpi_rank-1;
	if (prev<0) prev+=mpi_nprocs;
	if (mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " waiting for " 
		  << prev << "." << std::endl;
	}
	MPI_Recv(&buffer,1,MPI_INT,prev,tag,MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
      }
      
      // Determine which file to read next
      int file=filestart+k;
      if (file>=((int) n_sources)) file-= n_sources;

      if (mpi_load_debug) {
	scr_out << "Rank " << mpi_rank << " reading file " 
		<< file << "." << std::endl;
      }

      o2scl_hdf::hdf_file hf;
      hf.open(source_fnames[file]);
      if (table_names[file].length()>0) {
	hdf_input(hf,source_tables[file],table_names[file]);
      } else {
	hdf_input(hf,source_tables[file]);
      }
      hf.close();
      
      // Send a message, unless the rank is the last one to read a
      // file.
      if (k<((int)n_sources)-1 || mpi_rank<mpi_nprocs-((int)n_sources)) {
	int next=mpi_rank+1;
	if (next>=mpi_nprocs) next-=mpi_nprocs;
	if (mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " sending to " 
		  << next << "." << std::endl;
	}
	MPI_Send(&buffer,1,MPI_INT,next,tag,MPI_COMM_WORLD);
      }
      
    }
    
#else
    
    for(size_t k=0;k<n_sources;k++) {
      
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

    scr_out << "\nInput data files: " << std::endl;
    
    if (set.norm_max) {
      scr_out << "Normalizing maximum probability to 1." << std::endl;
    } else {
      scr_out << "Normalizing integral of distribution to 1." << std::endl;
    }

    scr_out << "File                          name   total        "
	    << "max          P(10,1.4)" << std::endl;

    for(size_t k=0;k<n_sources;k++) {
      scr_out << "Here: " << k << endl;
      
      // Update input limits
      if (k==0) {
	set.in_r_min=source_tables[k].get_grid_x(0);
	set.in_r_max=source_tables[k].get_grid_x
	  (source_tables[k].get_nx()-1);
	set.in_m_min=source_tables[k].get_grid_y(0);
	set.in_m_max=source_tables[k].get_grid_y
	  (source_tables[k].get_ny()-1);
      } else {
	if (set.in_r_min>source_tables[k].get_grid_x(0)) {
	  set.in_r_min=source_tables[k].get_grid_x(0);
	}
	if (set.in_r_max<source_tables[k].get_grid_x
	    (source_tables[k].get_nx()-1)) {
	  set.in_r_max=source_tables[k].get_grid_x
	    (source_tables[k].get_nx()-1);
	}
	if (set.in_m_min>source_tables[k].get_grid_y(0)) {
	  set.in_m_min= source_tables[k].get_grid_y(0);
	}
	if (set.in_m_max<source_tables[k].get_grid_y
	    (source_tables[k].get_ny()-1)) {
	  set.in_m_max= source_tables[k].get_grid_y
	    (source_tables[k].get_ny()-1);
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
	  if (set.norm_max) {
	     source_tables[k].set
	      (i,j,slice_names[k],source_tables[k].get
	       (i,j,slice_names[k])/max);		   
	  } else {
	     source_tables[k].set
	      (i,j,slice_names[k],source_tables[k].get
	       (i,j,slice_names[k])/tot);		   
	  }
	}
      }

      if (set.debug_load) {
	std::cout << source_fnames[k] << std::endl;
	for(size_t i=0;i<source_tables[k].get_nx();i++) {
	  std::cout << i << " " << source_tables[k].get_grid_x(i)
		    << std::endl;
	}
	for(size_t j=0;j<source_tables[k].get_ny();j++) {
	  std::cout << j << " " << source_tables[k].get_grid_y(j)
		    << std::endl;
	}
	for(size_t i=0;i<source_tables[k].get_nx();i++) {
	  for(size_t j=0;j<source_tables[k].get_ny();j++) {
	    std::cout << source_tables[k].get(i,j,slice_names[k])
		      << " ";
	  }
	  std::cout << std::endl;
	}
      }

      scr_out.setf(std::ios::left);
      scr_out.width(29);
      std::string stempx=source_fnames[k].substr(0,29);
      scr_out << stempx << " ";
      scr_out.width(6);
      scr_out << source_names[k] << " " << tot << " " << max << " ";
      scr_out.unsetf(std::ios::left);
      scr_out << source_tables[k].interp(10.0,1.4,slice_names[k])
	      << std::endl;
      
    }
    
    scr_out << std::endl;
  } else {
    if (set.verbose>=2) {
      cout << "bamr: no sources." << endl;
    }
  }    

  if (set.in_m_min<set.min_mass) set.in_m_min=set.min_mass;
  
  scr_out << "M limits: (" 
	  << set.in_m_min << "," << set.in_m_max << ")" << std::endl;
  scr_out << "R limits: ("
	  << set.in_r_min << "," << set.in_r_max << ")" << std::endl;
  
  return;
}
  
int ns_data::add_data(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<5) {
    std::cerr << "Not enough arguments given to 'add-data'." << std::endl;
    return o2scl::exc_efailed;
  }
      
  source_names.push_back(sv[1]);
  source_fnames.push_back(sv[2]);
  slice_names.push_back(sv[3]);
  init_mass_fracs.push_back(o2scl::stod(sv[4]));
  if (sv.size()==6) {
    table_names.push_back(sv[5]);
  } else {
    table_names.push_back("");
  }
      
  n_sources++;

  return 0;
}
    
void model::compute_star(const ubvector &pars, std::ofstream &scr_out, 
			 int &ret, model_data &dat) {

  if (set.verbose>=2) {
    cout << "Start model::compute_star()." << endl;
  }
  
  ret=ix_success;

  if (has_eos) {
    
    // ---------------------------------------------------------------
    // Compute the EOS
  
    if (set.verbose>=2) {
      cout << "Going to model::compute_eos()." << endl;
    }
    compute_eos(pars,ret,scr_out,dat);
    if (ret!=ix_success) return;
    if (set.verbose>=2) {
      cout << "Back from model::compute_eos()." << endl;
    }

    // Ensure we're using linear interpolation
    dat.eos->set_interp_type(o2scl::itp_linear);

    // ---------------------------------------------------------------
    // Check that pressure is increasing. If we're using a crust EOS,
    // choose 0.6 as an arbitrary low-density cutoff which corresponds
    // to about n_B=0.12 fm^{-3}, and check the low-density part below
    // instead.
  
    for(size_t i=0;(dat.eos->get_nlines()>0 &&
		    i<dat.eos->get_nlines()-1);i++) {
      if ((!set.use_crust || dat.eos->get("ed",i)>0.6) && 
	  dat.eos->get("pr",i+1)<dat.eos->get("pr",i)) {
	scr_out << "Rejected: Pressure decreasing." << std::endl;
	scr_out << "ed=" << dat.eos->get("ed",i) 
		<< " pr=" << dat.eos->get("pr",i) << std::endl;
	scr_out << "ed=" << dat.eos->get("ed",i+1) 
		<< " pr=" << dat.eos->get("pr",i+1) << std::endl;
	ret=ix_press_dec;
	return;
      }
    }
  }

  // ---------------------------------------------------------------
  // If requested, compute the baryon density automatically
  
  if (has_eos && set.baryon_density && !dat.eos->is_column("nb")) {
    
    // Obtain the baryon density calibration point from the model
    if (nb_n1<=0.0 && nb_e1<=0.0) {
      O2SCL_ERR2("Computing the baryon density requires one ",
		 "calibration point in compute_star().",
		 o2scl::exc_einval);
    }

    // Compute inverse of gibbs energy density, 'igb'

    dat.eos->new_column("igb");
    dat.eos->set_unit("igb","fm^4");

    for(size_t i=0;i<dat.eos->get_nlines();i++) {

      if (dat.eos->get("ed",i)+dat.eos->get("pr",i)<=0.0 ||
	  !std::isfinite(dat.eos->get("ed",i)) ||
	  !std::isfinite(dat.eos->get("pr",i))) {
	scr_out << "Inverse Gibbs not finite." << std::endl;
	scr_out << "n1=" << nb_n1 << " e1=" << nb_e1 << std::endl;
	scr_out << "ed pr" << std::endl;
	for(size_t i=0;i<dat.eos->get_nlines();i++) {
	  scr_out << i << " "
		  << dat.eos->get("ed",i) << " "
		  << dat.eos->get("pr",i) << std::endl;
	}
	O2SCL_ERR("Inverse Gibbs not finite.",o2scl::exc_efailed);
      }
      dat.eos->set("igb",i,1.0/(dat.eos->get("ed",i)+dat.eos->get("pr",i)));
    }

    // Compute integral of 'igb' relative to ed='e1', called 'iigb'

    dat.eos->new_column("iigb");
    
    for(size_t i=0;i<dat.eos->get_nlines();i++) {
      if (nb_e1<=dat.eos->get("ed",i)) {
	double val=dat.eos->integ("ed",nb_e1,dat.eos->get("ed",i),"igb");
	if (!std::isfinite(val)) {
	  scr_out << "Baryon integral not finite." << std::endl;
	  scr_out << "n1=" << nb_n1 << " e1=" << nb_e1 << std::endl;
	  scr_out << "ed pr" << std::endl;
	  for(size_t i=0;i<dat.eos->get_nlines();i++) {
	    scr_out << i << " "
		    << dat.eos->get("ed",i) << " "
		    << dat.eos->get("pr",i) << std::endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",o2scl::exc_efailed);
	}
	dat.eos->set("iigb",i,val);
      } else {
	double val=-dat.eos->integ("ed",dat.eos->get("ed",i),nb_e1,"igb");
	if (!std::isfinite(val)) {
	  scr_out << "Baryon integral not finite (2)." << std::endl;
	  scr_out << "n1=" << nb_n1 << " e1=" << nb_e1 << std::endl;
	  scr_out << "ed pr" << std::endl;
	  for(size_t i=0;i<dat.eos->get_nlines();i++) {
	    scr_out << i << " "
		    << dat.eos->get("ed",i) << " "
		    << dat.eos->get("pr",i) << std::endl;
	  }
	  O2SCL_ERR("Baryon integral not finite.",o2scl::exc_efailed);
	}
	dat.eos->set("iigb",i,val);
      }
    }

    // Compute normalization constant
    double Anb=nb_n1/exp(dat.eos->interp("ed",nb_e1,"iigb"));
    if (!std::isfinite(Anb) || Anb<0.0) {
      scr_out << "Baryon density normalization problem." << std::endl;
      ret=ix_nb_problem;
      return;
    }

    // Now compute baryon density

    dat.eos->new_column("nb");
    dat.eos->set_unit("nb","1/fm^3");

    for(size_t i=0;i<dat.eos->get_nlines();i++) {      

      // If the density is too low, then just use zero baryon density.
      // This pressure (10^{-5} fm^{-4}) corresponds to a baryon 
      // density of about 3e-3 fm^{-3}.

      if (set.use_crust==false && dat.eos->get("pr",i)<1.0e-5) {

	dat.eos->set("nb",i,0.0);

      } else {
	
	double nbt=Anb*exp(dat.eos->get("iigb",i));
	if (!std::isfinite(nbt)) {
	  scr_out << "Baryon density normalization problem (2)."
		  << std::endl;
	  ret=ix_nb_problem2;
	  return;
	} 
	dat.eos->set("nb",i,nbt);
      }
    }
    
    // End of loop 'if (has_eos && baryon_density && 
    // !dat.eos->is_column("nb")) {' 
  }

  // ---------------------------------------------------------------

  if (has_eos) {
    
    // ---------------------------------------------------------------
    // Compute crust transition density and set crust (if necessary)

    if (set.compute_cthick && set.baryon_density) {

      // If crust_from_L is true, compute the pressure at the number density
      // specified by the correlation. Otherwise, just use 0.08

      if (set.crust_from_L) {

	if (set.verbose>=2) {
	  std::cout << "Starting crust from L." << std::endl;
	}
	
	double asamp, bsamp, csamp;
	if (false) {
	  double f1=pars[0]*1.0e5-floor(pars[0]*1.0e5);
	  if (f1<=0.0) f1=1.0e-5;
	  else if (f1>=1.0) f1=1.0-1.0e-5;
	  asamp=nt_a.invert_cdf(f1);

	  double f2=pars[1]*1.0e5-floor(pars[1]*1.0e5);
	  f2=0.5;
	  if (f2<=0.0) f2=1.0e-5;
	  else if (f2>=1.0) f2=1.0-1.0e-5;
	  bsamp=nt_b.invert_cdf(f2);

	  double f3=pars[2]*1.0e5-floor(pars[2]*1.0e5);
	  if (f3<=0.0) f3=1.0e-5;
	  else if (f3>=1.0) f3=1.0-1.0e-5;
	  csamp=nt_c.invert_cdf(f3);
	} else {
	  asamp=nt_a.mean();
	  bsamp=nt_b.mean();
	  csamp=nt_c.mean();
	}
      
	double St=dat.eos->get_constant("S")*o2scl_const::hc_mev_fm;
	double Lt=dat.eos->get_constant("L")*o2scl_const::hc_mev_fm;
	double nt=(asamp+bsamp*(Lt/70.0)+csamp*(Lt*Lt/4900.0))*(St/30.0);

	if (nt<nt_low || nt>nt_high) {
	  scr_out << "Transition density, " << nt << ", out of range." << endl;
	  ret=ix_trans_invalid;
	  return;
	}
	dat.eos->add_constant("nt",nt);
	
	double prt=dat.eos->interp("nb",nt,"pr");
	dat.eos->add_constant("prt",prt);
	
	// Add the transition pressure to the tov_solve object
	
	if (ts.pr_list.size()>0) {
	  ts.pr_list.clear();
	}
	ts.pr_list.push_back(prt);
	
	// Set the crust and it's transition pressure
      
	if (dat.eos->get_constant("S")*hc_mev_fm<28.0 || 
	    dat.eos->get_constant("S")*hc_mev_fm>38.0 || 
	    dat.eos->get_constant("L")*hc_mev_fm<25.0 ||
	    dat.eos->get_constant("L")*hc_mev_fm>115.0 ||
	    dat.eos->get_constant("L")*hc_mev_fm>
	    dat.eos->get_constant("S")*hc_mev_fm*5.0-65.0) {
	  scr_out << "S or L out of range" << endl;
	  ret=ix_SL_invalid;
	  return;
	}
      
	// This function expects its argument in MeV
	teos.ngl13_low_dens_eos2
	  (dat.eos->get_constant("S")*hc_mev_fm,
	   dat.eos->get_constant("L")*hc_mev_fm,nt,"");
      
	// Set the transition pressure and width. Note that
	// the ngl13 EOS is in units of Msun/km^3, so we 
	// convert prt to those units
	teos.transition_mode=eos_tov_interp::match_line;
	teos.set_transition(o2scl_settings.get_convert_units().convert
			    ("1/fm^4","Msun/km^3",prt),1.4);

	if (set.verbose>=2) {
	  std::cout << "Done with crust from L." << std::endl;
	}
	
      } else {

	// Otherwise, if we're not determining the crust from L, then
	// compute the crust thickness based on a density of 0.08
	// fm^{-3}

	double nt=0.08;
	dat.eos->add_constant("nt",nt);
	double prt=dat.eos->interp("nb",0.08,"pr");
	dat.eos->add_constant("prt",prt);
	if (ts.pr_list.size()>0) {
	  ts.pr_list.clear();
	}
	ts.pr_list.push_back(prt);
	
      }
    
    } else {

      if (ts.pr_list.size()>0) {
	ts.pr_list.clear();
      }

    }

    // ---------------------------------------------------------------
    // Read the EOS into the tov_eos object
    
    table_units<> &teos_temp=*(dat.eos);
    if (set.baryon_density && set.inc_baryon_mass) {
      dat.eos->set_unit("ed","1/fm^4");
      dat.eos->set_unit("pr","1/fm^4");
      dat.eos->set_unit("nb","1/fm^3");
      if (set.verbose>=2) {
	std::cout << "Going to read_table() (with nb)." << std::endl;
      }
      teos.read_table(teos_temp,"ed","pr","nb");
      if (set.verbose>=2) {
	std::cout << "Done in read_table() (with nb)." << std::endl;
      }
    } else {
      dat.eos->set_unit("ed","1/fm^4");
      dat.eos->set_unit("pr","1/fm^4");
      if (set.verbose>=2) {
	std::cout << "Going to read_table() (without nb)." << std::endl;
      }
      teos.read_table(teos_temp,"ed","pr");
      if (set.verbose>=2) {
	std::cout << "Done in read_table() (without nb)." << std::endl;
      }
    }

    // ---------------------------------------------------------------
    // We checked that the pressure of the core EOS was increasing
    // earlier. Here we also double check that the EOS is increasing
    // near the crust-core transition.
    
    if (set.use_crust) {
    
      double ed_last=0.0;
      // This range corresponds to between about n_B=0.01 and 0.17
      // fm^{-3}
      for(double pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	double ed, nb;
	teos.ed_nb_from_pr(pr,ed,nb);
	if (ed_last>1.0e-20 && ed<ed_last) {
	  scr_out << "Stability problem near crust-core transition."
		  << std::endl;
	  if (has_esym) {
	    scr_out << "S=" << dat.eos->get_constant("S")*hc_mev_fm 
		    << " L=" << dat.eos->get_constant("L")*hc_mev_fm
		    << std::endl;
	  }
	  scr_out << "Energy decreased with increasing pressure "
		  << "at pr=" << pr << std::endl;
	  scr_out << std::endl;
	  scr_out << "Full EOS near transition: " << std::endl;
	  scr_out << "pr ed" << std::endl;
	  for(pr=1.0e-4;pr<2.0e-2;pr*=1.1) {
	    teos.ed_nb_from_pr(pr,ed,nb);
	    scr_out << pr << " " << ed << std::endl;
	  }
	  scr_out << std::endl;
	  ret=ix_crust_unstable;
	  return;
	}
	ed_last=ed;
      }

      // End of 'if (set.use_crust)'
    }

    // ---------------------------------------------------------------
    // If necessary, output debug information (We want to make sure
    // this is after tov_eos::read_table() so that we can debug the
    // core-crust transition)
    
    if (set.debug_eos) {
      o2scl_hdf::hdf_file hfde;

      hfde.open_or_create("debug_eos.o2");

      // Output the core EOS
      hdf_output(hfde,*dat.eos,"eos");

      // Output core and crust EOS as reported by the tov_interp_eos object
      o2scl::table_units<> full_eos;
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
      if (!set.debug_star) {
	scr_out << "Automatically exiting since 'debug_eos' is true."
		<< std::endl;
	exit(0);
      }
    }

    // ---------------------------------------------------------------
    // Solve for M vs. R curve
    
    double m_max=0.0;
    ts.princ=set.mvsr_pr_inc;
    ts.set_table(dat.mvsr);
    if (set.addl_quants) {
      ts.ang_vel=true;
      ts.calc_gpot=true;
    } else {
      ts.ang_vel=false;
      ts.calc_gpot=false;
    }
    if (set.verbose>=2) {
      std::cout << "Going to TOV." << std::endl;
    }
    int info=ts.mvsr();
    if (set.verbose>=2) {
      cout << "Done with TOV." << endl;
    }
    if (info!=0) {
      scr_out << "M vs. R failed: info=" << info << std::endl;
      ret=ix_mvsr_failed;
      return;
    }
    dat.mvsr->set_interp_type(o2scl::itp_linear);

    // ---------------------------------------------------------------
    // Add baryon density to M vs. R table if it's not already there

    if (set.baryon_density && !set.inc_baryon_mass) {
      dat.mvsr->add_col_from_table(*dat.eos,"pr","nb","pr");
      dat.mvsr->set_unit("nb","1/fm^3");
    }
    
    // ---------------------------------------------------------------
    // Check that maximum mass is large enough. Note that the
    // mass-radius curve is not differentiable near M_{max}, so the
    // best way to increase the accuracy here is to make
    // set.mvsr_pr_inc smaller.
    
    m_max=dat.mvsr->max("gm");
    dat.mvsr->add_constant("m_max",m_max);
    if (m_max<set.min_max_mass) {
      scr_out << "Maximum mass too small: " << m_max << " < "
	      << set.min_max_mass << "." << std::endl;
      ret=ix_small_max;
      return;
    }

    // ---------------------------------------------------------------
    // If the EOS is sufficiently stiff, the TOV solver will output
    // gibberish, i.e. masses and radii equal to zero, especially at
    // the higher pressures. This rules these EOSs out, as they would
    // likely be acausal anyway. If this happens frequently, it might
    // signal a problem.
    
    size_t ir=dat.mvsr->get_nlines()-1;
    if ((*dat.mvsr)["gm"][ir]<1.0e-10 ||
	(*dat.mvsr)["gm"][ir-1]<1.0e-10) {
      scr_out << "TOV failure fix." << std::endl;
      ret=ix_tov_failure;
      return;
    }

    // ---------------------------------------------------------------
    // Check the radius of the maximum mass star
    
    size_t ix_max=dat.mvsr->lookup("gm",m_max);
    double r_max=dat.mvsr->get("r",ix_max);
    dat.mvsr->add_constant("r_max",r_max);
    if (r_max>1.0e4) {
      scr_out << "TOV convergence problem: " << std::endl;
      ret=ix_tov_conv;
      return;
    }
    
    // ---------------------------------------------------------------

    // Compute the central baryon density in the maximum mass star
    if (set.baryon_density) {
      double nb_max=dat.mvsr->get("nb",ix_max);
      dat.mvsr->add_constant("nb_max",nb_max);
    }
    
    // Remove table entries with pressures above the maximum pressure
    dat.mvsr->set_nlines(ix_max+1);
  
    // Make sure that the M vs. R curve generated enough data. This
    // is not typically an issue.
    if (dat.mvsr->get_nlines()<10) {
      scr_out << "M vs. R failed to generate lines." << std::endl;
      ret=ix_mvsr_table;
      return;
    }

    // Compute the masses and radii for each source
    for(size_t i=0;i<nsd.n_sources;i++) {
      if (set.mass_switch==0) {
	dat.mass[i]=m_max*pars[this->n_eos_params+i];
      } else if (set.mass_switch==1) {
	dat.mass[i]=0.4*pars[this->n_eos_params+i]+1.3;
      } else {
	dat.mass[i]=0.2*pars[this->n_eos_params+i]+1.3;
      }
      dat.rad[i]=dat.mvsr->interp("gm",dat.mass[i],"r");
    }

    // End of loop 'if (has_eos)'
  } else {

    // Compute mass-radius curve directly
    compute_mr(pars,ret,scr_out,dat);
    if (ret!=ix_success) {
      return;
    }

    // Compute the masses and radii for each source
    for(size_t i=0;i<nsd.n_sources;i++) {
      dat.mass[i]=pars[this->n_eos_params+i];
      dat.rad[i]=dat.mvsr->interp("gm",dat.mass[i],"r");
    }

  }
  
  // Output M vs. R curve
  if (set.debug_star) {
    o2scl_hdf::hdf_file hfds;
    hfds.open_or_create("debug_star.o2");
    hdf_output(hfds,*dat.mvsr,"mvsr");
    hfds.close();
    scr_out << "Automatically exiting since 'debug_star' is true."
	    << std::endl;
    exit(0);
  }

  // -----------------------------------------------------------------
  // Check causality. Note that we have to do this after the rows for
  // the unstable branch have been removed from the mass-radius table.
  
  if (has_eos) {
    
    // Compute speed of sound squared
    dat.mvsr->deriv("ed","pr","dpde");

    for(size_t i=0;i<dat.mvsr->get_nlines();i++) {
      if ((*dat.mvsr)["dpde"][i]>1.0) {
	scr_out.precision(4);
	scr_out << "Rejected: Acausal."<< std::endl;
	scr_out << "ed_max="
		<< dat.mvsr->max("ed") << " ed_bad="
		<< (*dat.mvsr)["ed"][i] << " pr_max=" 
		<< dat.mvsr->max("pr") << " pr_bad=" 
		<< (*dat.mvsr)["pr"][i] << std::endl;
	scr_out.precision(6);
	ret=ix_acausal;
	return;
      }
    }
    
  } else {

    for(size_t i=0;i<nsd.n_sources;i++) {
      if (dat.rad[i]<2.94*schwarz_km/2.0*dat.mass[i]) {
	scr_out << "Source " << nsd.source_names[i] << " acausal."
		<< std::endl;
	ret=ix_acausal_mr;
	return;
      }
    }

  }
  
  // ---------------------------------------------------------------
  // Compute M, R for fixed central baryon densities

  if (has_eos && set.baryon_density) {
    double ed1=dat.eos->interp("nb",0.16,"ed");
    dat.mvsr->add_constant("gm_nb1",dat.mvsr->interp("ed",ed1,"gm"));
    dat.mvsr->add_constant("r_nb1",dat.mvsr->interp("ed",ed1,"r"));
    double ed2=dat.eos->interp("nb",0.32,"ed");
    dat.mvsr->add_constant("gm_nb2",dat.mvsr->interp("ed",ed2,"gm"));
    dat.mvsr->add_constant("r_nb2",dat.mvsr->interp("ed",ed2,"r"));
    double ed3=dat.eos->interp("nb",0.48,"ed");
    dat.mvsr->add_constant("gm_nb3",dat.mvsr->interp("ed",ed3,"gm"));
    dat.mvsr->add_constant("r_nb3",dat.mvsr->interp("ed",ed3,"r"));
    double ed4=dat.eos->interp("nb",0.64,"ed");
    dat.mvsr->add_constant("gm_nb4",dat.mvsr->interp("ed",ed4,"gm"));
    dat.mvsr->add_constant("r_nb4",dat.mvsr->interp("ed",ed4,"r"));
    double ed5=dat.eos->interp("nb",0.80,"ed");
    dat.mvsr->add_constant("gm_nb5",dat.mvsr->interp("ed",ed5,"gm"));
    dat.mvsr->add_constant("r_nb5",dat.mvsr->interp("ed",ed5,"r"));
  }

  if (set.verbose>=2) {
    cout << "End model::compute_star()." << endl;
  }

  return;
}

int model::compute_point(const ubvector &pars, std::ofstream &scr_out, 
			 double &log_weight, model_data &dat) {

  if (set.verbose>=2) {
    cout << "Start model::compute_point()." << endl;
  }

  // Compute the M vs R curve and return if it failed
  int iret;
  compute_star(pars,scr_out,iret,dat);
  if (iret!=ix_success) {
    log_weight=0.0;
    return iret;
  }
  
  for(size_t i=0;i<nsd.n_sources;i++) {
    if (dat.mass[i]<set.in_m_min || dat.mass[i]>set.in_m_max || 
	dat.rad[i]<set.in_r_min || dat.rad[i]>set.in_r_max) {
      scr_out << "Rejected: Mass or radius outside range." << std::endl;
      scr_out << "M limits: " << set.in_m_min << " "
	      << set.in_m_max << std::endl;
      scr_out << "R limits: " << set.in_r_min << " "
	      << set.in_r_max << std::endl;
      if (nsd.n_sources>0) {
	scr_out.precision(2);
	scr_out.setf(ios::showpos);
	scr_out << "M ";
	for(size_t i=0;i<nsd.n_sources;i++) {
	  scr_out << dat.mass[i] << " ";
	}
	scr_out << std::endl;
	scr_out << "R ";
	for(size_t i=0;i<nsd.n_sources;i++) {
	  scr_out << dat.rad[i] << " ";
	}
	scr_out << std::endl;
	scr_out.precision(6);
	scr_out.unsetf(ios::showpos);
      }
      log_weight=0.0;
      return ix_mr_outside;
    }
  }

  log_weight=0.0;
      
  dat.mvsr->set_interp_type(o2scl::itp_linear);
  double m_max_current=dat.mvsr->max("gm");

  // -----------------------------------------------
  // Compute the weights for each source
      
  if (set.debug_star) scr_out << "Name M R Weight" << std::endl;
  
  for(size_t i=0;i<nsd.n_sources;i++) {
	
    // Double check that current M and R is in the range of
    // the provided input data
    if (dat.rad[i]<nsd.source_tables[i].get_x_data()[0] ||
	dat.rad[i]>nsd.source_tables[i].get_x_data()
	[nsd.source_tables[i].get_nx()-1] ||
	dat.mass[i]<nsd.source_tables[i].get_y_data()[0] ||
	dat.mass[i]>nsd.source_tables[i].get_y_data()
	[nsd.source_tables[i].get_ny()-1]) {
      dat.wgts[i]=0.0;
    } else {
      // If it is, compute the weight
      dat.wgts[i]=nsd.source_tables[i].interp
	(dat.rad[i],dat.mass[i],nsd.slice_names[i]);
				      
    }
	
    // If the weight is lower than the threshold, set it equal
    // to the threshold
    if (dat.wgts[i]<set.input_dist_thresh) {
      dat.wgts[i]=set.input_dist_thresh;
    }
    // Include the weight for this source
    if (dat.wgts[i]<=0.0) {
      scr_out << "Weight zero for source " << nsd.source_names[i] << endl;
      return ix_mr_outside;
    }
    log_weight+=log(dat.wgts[i]);
	
    if (set.debug_star) {
      scr_out << nsd.source_names[i] << " "
	      << dat.mass[i] << " " 
	      << dat.rad[i] << " " << dat.wgts[i] << std::endl;
    }
    if (set.verbose>=2) {
      cout << nsd.source_names[i] << " "
	   << dat.mass[i] << " " 
	   << dat.rad[i] << " " << dat.wgts[i] << std::endl;
    }
    
    // Go to the next source
  }
      
  if (set.debug_star) scr_out << std::endl;
      
  // -----------------------------------------------
  // Exit if the current maximum mass is too large
      
  if (m_max_current>set.exit_mass) {
    scr_out.setf(ios::scientific);
    scr_out << "Exiting because maximum mass (" << m_max_current 
	    << ") larger than exit_mass (" << set.exit_mass << ")." 
	    << std::endl;
    scr_out.precision(12);
    vector_out(scr_out,pars);
    scr_out << " " << log_weight << std::endl;
    scr_out.precision(6);
    exit(0);
  }

  if (set.verbose>=2) {
    cout << "End model::compute_point()." << endl;
  }

  if (iret!=ix_success) {
    // We shouldn't be returning a non-zero value if success is
    // non-zero, so we double check this here
    O2SCL_ERR("Sanity check for success flag in model::compute_point.",
	      o2scl::exc_esanity);
  }
  return o2scl::success;
}

void two_polytropes::setup_params(o2scl::cli &cl) {
  p_kin_sym.d=&se.a;
  p_kin_sym.help="Kinetic part of symmetry energy.";
  cl.par_list.insert(make_pair("kin_sym",&p_kin_sym));

  return;
}

void two_polytropes::copy_params(model &m) {
  // Dynamic casts throw exceptions when they fail
  two_polytropes &tp=dynamic_cast<two_polytropes &>(m);
  se.a=tp.se.a;
  
  return;
}

void two_polytropes::remove_params(o2scl::cli &cl) {
  size_t i=cl.par_list.erase("kin_sym");
  if (i!=1) {
    O2SCL_ERR("Failed to erase parameter 'kin_sym'.",o2scl::exc_esanity);
  }
  return;
}

two_polytropes::two_polytropes(settings &s, ns_data &n) : model(s,n) {

  se.kpp=0.0;
  se.n0=0.16;
  se.eoa=-16.0/hc_mev_fm;
  se.a=17.0/hc_mev_fm;
    
  // We include muons by default, but they rarely appear at low-density
  cns.include_muons=true;

  this->n_eos_params=8;
  this->has_esym=true;
}

void two_polytropes::get_param_info(std::vector<std::string> &names,
				    std::vector<std::string> &units,
				    ubvector &low, ubvector &high) {

  names={"comp","kprime","esym","gamma","trans1","index1",
	 "trans2","index2"};
  
  units={"1/fm","1/fm","1/fm","","1/fm^4","","1/fm^4",""};

  low.resize(n_eos_params+nsd.n_sources);
  low[0]=180.0/hc_mev_fm;
  low[1]=-1000.0/hc_mev_fm;
  low[2]=28.0/hc_mev_fm;
  low[3]=0.0;
  // The value 0.75 fm^{-4} is about the energy density of nuclear
  // matter
  low[4]=0.75;
  low[5]=0.2;
  low[6]=0.75;
  low[7]=0.2;
  
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=300.0/hc_mev_fm;
  // FSU gold is -280 MeV or so
  high[1]=-200.0/hc_mev_fm;
  high[2]=38.0/hc_mev_fm;
  high[3]=1.2;
  // The value of high.trans1 has to be small enough because we
  // don't want to evaluate the schematic EOS to too high of a
  // density.
  high[4]=3.0;
  high[5]=1.5;
  high[6]=8.0;
  high[7]=2.0;

  // Go to the parent which takes care of the data-related
  // parameters
  model::get_param_info(names,units,low,high);
  
  return;
}

void two_polytropes::initial_point(ubvector &params) {
  params[0]=1.0;
  params[1]=-3.0;
  params[2]=0.165;
  params[3]=0.644;
  params[4]=1.51;
  params[5]=0.576;
  params[6]=4.60;
  params[7]=1.21;
  model::initial_point(params);
  return;
}

void two_polytropes::compute_eos(const ubvector &params, int &ret,
				 ofstream &scr_out, model_data &dat) {

  ret=ix_success;
  if (params[4]>params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    ret=ix_param_mismatch;
    return;
  }
  
  // Set hadronic EOS from ubvector information
  se.comp=params[0];
  se.kprime=params[1];
  se.b=params[2]-se.a;
  se.gamma=params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.set_eos_table(dat.eos);
  cns.calc_eos();
  dat.eos->set_interp_type(itp_linear);

  dat.eos->add_constant("S",params[2]);
  dat.eos->add_constant("L",se.fesym_slope(0.16));
  
  // Transition densities
  double ed1=params[4];
  double ed2=params[6];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=dat.eos->interp("ed",ed1,"nb");
  double pr1=dat.eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,1.0+1.0/params[5]);

  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    ret=ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=dat.eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<dat.eos->get_nlines();i++) {
    if ((*dat.eos)["ed"][i]>ed1) {
      dat.eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    ret=ix_no_eos_table;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=ed1;ed<ed2-0.001;ed+=0.05) {
    double pr=coeff1*pow(ed,1.0+1.0/params[5]);
    double nb=nb1*pow(ed/ed1,1.0+params[5])/
      pow((ed+pr)/(ed1+pr1),params[5]);
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    ret=ix_no_eos_table;
    return;
  }
  
  // Boundary baryon density and pressure
  double pr2=coeff1*pow(ed2,1.0+1.0/params[5]);
  double nb2=nb1*pow(ed2/ed1,1.0+params[5])/
    pow((ed2+pr2)/(ed1+pr1),params[5]);

  // Determine 2nd polytrope coefficient
  double coeff2=pr2/pow(ed2,1.0+1.0/params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    ret=ix_neg_pressure;
    return;
  }

  // Add second polytrope to table
  for(double ed=ed2;ed<=10.0;ed+=0.05) {
    double pr=coeff2*pow(ed,1.0+1.0/params[7]);
    double nb=nb2*pow(ed/ed2,1.0+params[7])/
      pow((ed+pr)/(ed2+pr2),params[7]);
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }
  
  return;
}

void alt_polytropes::get_param_info(std::vector<std::string> &names,
				    std::vector<std::string> &units,
				    ubvector &low, ubvector &high) {

  two_polytropes::get_param_info(names,units,low,high);

  low[5]=1.5;
  low[7]=0.0;

  high[5]=6.0;
  high[7]=3.0;
  
  names={"comp","kprime","esym","gamma","trans1","exp1",
	 "trans2","exp2"};
  
  units={"1/fm","1/fm","1/fm","","1/fm^4","","1/fm^4",""};
  
  return;
}

void alt_polytropes::initial_point(ubvector &params) {
  params[0]=1.0;
  params[1]=-2.66;
  params[2]=0.165;
  params[3]=0.66;
  params[4]=1.48;
  params[5]=2.913;
  params[6]=4.066;
  params[7]=1.80;
  model::initial_point(params);
  return;
}

void alt_polytropes::compute_eos(const ubvector &params, int &ret,
				 ofstream &scr_out, model_data &dat) {
  
  ret=ix_success;
  if (params[4]>params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    ret=ix_param_mismatch;
    return;
  }

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from ubvector information
  se.comp=params[0];
  se.kprime=params[1];
  se.b=params[2]-se.a;
  se.gamma=params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.set_eos_table(dat.eos);
  cns.calc_eos();
  dat.eos->set_interp_type(itp_linear);

  dat.eos->add_constant("S",params[2]);
  dat.eos->add_constant("L",se.fesym_slope(0.16));

  // Transition densities
  double ed1=params[4];
  double ed2=params[6];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=dat.eos->interp("ed",ed1,"nb");
  double pr1=dat.eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,params[5]);
    
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    ret=ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=dat.eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<dat.eos->get_nlines();i++) {
    if ((*dat.eos)["ed"][i]>ed1) {
      dat.eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    ret=ix_no_eos_table;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=ed1;ed<ed2-0.001;ed+=0.05) {
    double pr=coeff1*pow(ed,params[5]);
    double nb=nb1*pow(ed/ed1,params[5]/(params[5]-1.0))*
      pow((ed+pr)/(ed1+pr1),1.0/(1.0-params[5]));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    ret=ix_no_eos_table;
    return;
  }

  // Boundary baryon density and pressure
  double pr2=coeff1*pow(ed2,params[5]);
  double nb2=nb1*pow(ed2/ed1,params[5]/(params[5]-1.0))*
    pow((ed2+pr2)/(ed1+pr1),1.0/(1.0-params[5]));
  
  // Determine 2nd polytrope coefficient
  double coeff2=pr2/pow(ed2,params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    ret=ix_neg_pressure;
    return;
  }

  // Add second polytrope to table
  for(double ed=ed2;ed<=10.0;ed+=0.05) {
    double pr=coeff2*pow(ed,params[7]);
    double nb=nb2*pow(ed/ed2,params[7]/(params[7]-1.0))*
      pow((ed+pr)/(ed2+pr2),1.0/(1.0-params[7]));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  return;
}
  
void fixed_pressure::get_param_info(std::vector<std::string> &names,
				    std::vector<std::string> &units,
				    ubvector &low, ubvector &high) {

  two_polytropes::get_param_info(names,units,low,high);

  low[4]=0.0;
  low[5]=0.0;
  low[6]=0.0;
  low[7]=0.0;
  
  high[4]=0.3;
  // This parameter is the pressure at 3 fm^{-4} minus the pressure at
  // 2 fm^{-4}, thus from causality cannot be larger than 1 fm^{-4}.
  high[5]=1.0;
  // This parameter is the pressure at 5 fm^{-4} minus the pressure at
  // 3 fm^{-4}, thus from causality cannot be larger than 2 fm^{-4}.
  high[6]=2.0;
  // This parameter is the pressure at 7 fm^{-4} minus the pressure at
  // 5 fm^{-4}. It is not always limited by causality, because the
  // central energy density is sometimes smaller than 5 fm^{-4}. We
  // allow it to be larger than 2.0 fm^{-4}, but in practice larger
  // values are relatively rare.
  high[7]=2.5;

  names[4]="pres1";
  names[5]="pres2";
  names[6]="pres3";
  names[7]="pres4";
  
  units[4]="1/fm^4";
  units[5]="1/fm^4";
  units[6]="1/fm^4";
  units[7]="1/fm^4";

  return;
}

void fixed_pressure::initial_point(ubvector &params) {
  params[0]=1.0;
  params[1]=-2.5;
  params[2]=0.165;
  params[3]=0.8;
  params[4]=0.024;
  params[5]=0.74;
  params[6]=0.60;
  params[7]=1.84;
  model::initial_point(params);
  return;
}

void fixed_pressure::compute_eos(const ubvector &params, int &ret,
				 ofstream &scr_out, model_data &dat) {

  ret=ix_success;

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from ubvector information
  se.comp=params[0];
  se.kprime=params[1];
  se.b=params[2]-se.a;
  se.gamma=params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.set_eos_table(dat.eos);
  cns.calc_eos();
  dat.eos->set_interp_type(itp_linear);

  dat.eos->add_constant("S",params[2]);
  dat.eos->add_constant("L",se.fesym_slope(0.16));
  
  // Compute boundary energy density, baryon density and pressure
  double ed_last=1.0;
  double nb_last=dat.eos->interp("ed",1.0,"nb");
  double pr_last=dat.eos->interp("ed",1.0,"pr");
    
  // Remove extra rows from EOS near saturation
  for(size_t i=0;i<dat.eos->get_nlines();i++) {
    if ((*dat.eos)["ed"][i]>1.0) {
      dat.eos->delete_row(i);
      i=0;
    }
  }

  // Compute pressures on grid, ed=2, 3, 5, 7 fm^{-4}
  double p2=pr_last+params[4];
  double p3=p2+params[5];
  double p5=p3+params[6];

  // Computes slopes (squared sound speeds)
  double cs2_1=params[4]/(2.0-1.0);
  double cs2_2=params[5]/(3.0-2.0);
  double cs2_3=params[6]/(5.0-3.0);
  double cs2_4=params[7]/(7.0-5.0);

  // Add 1st high-density EOS
  for(double ed=1.0;ed<2.0-1.0e-4;ed+=0.1) {
    double pr=pr_last+params[4]*(ed-1.0)/(2.0-1.0);
    double nb=nb_last*pow((ed+pr_last+cs2_1*(ed-ed_last))/
			  (ed_last+pr_last),1.0/(1.0+cs2_1));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }
  double nb2=nb_last*pow((2.0+pr_last+cs2_1*(2.0-ed_last))/
			 (ed_last+pr_last),1.0/(1.0+cs2_1));
  
  // Add 2nd high-density EOS
  for(double ed=2.0;ed<3.0-1.0e-4;ed+=0.1) {
    double pr=p2+params[5]*(ed-2.0)/(3.0-2.0);
    double nb=nb2*pow((ed+p2+cs2_2*(ed-2.0))/(2.0+p2),1.0/(1.0+cs2_2));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }
  double nb3=nb2*pow((3.0+p2+cs2_2*(3.0-2.0))/(2.0+p2),1.0/(1.0+cs2_2));

  // Add 3rd high-density EOS
  for(double ed=3.0;ed<5.0-1.0e-4;ed+=0.1) {
    double pr=p3+params[6]*(ed-3.0)/(5.0-3.0);
    double nb=nb3*pow((ed+p3+cs2_3*(ed-3.0))/(3.0+p3),1.0/(1.0+cs2_3));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }
  double nb5=nb3*pow((5.0+p3+cs2_3*(5.0-3.0))/(3.0+p3),1.0/(1.0+cs2_3));

  // Add 4th high-density EOS
  for(double ed=5.0;ed<10.0-1.0e-4;ed+=0.2) {
    double pr=p5+params[7]*(ed-5.0)/(7.0-5.0);
    double nb=nb5*pow((ed+p5+cs2_4*(ed-5.0))/(5.0+p5),1.0/(1.0+cs2_4));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  return;
}
  
void generic_quarks::get_param_info(std::vector<std::string> &names,
				    std::vector<std::string> &units,
				    ubvector &low, ubvector &high) {

  names={"comp","kprime","esym","gamma","trans1","exp1",
	 "trans2","a2","a4"};

  units={"1/fm","1/fm","1/fm","","1/fm^4","","1/fm^4","1/fm^2",""};
  
  low.resize(n_eos_params+nsd.n_sources);
  low[0]=180.0/hc_mev_fm;
  low[1]=-1000.0/hc_mev_fm;
  low[2]=28.0/hc_mev_fm;
  //low[3]=0.2;
  low[3]=0.0;
  // 0.75 is about the energy density of nuclear matter
  low[4]=0.75;
  low[5]=2.0;
  low[6]=0.75;
  // a2
  low[7]=-0.3;
  // a4
  low[8]=0.045;
    
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=300.0/hc_mev_fm;
  // FSU gold is -280 or so
  high[1]=-200.0/hc_mev_fm;
  high[2]=38.0/hc_mev_fm;
  high[3]=1.2;
  // The value of high.trans1 has to be small enough because we
  // don't want to evaluate the schematic EOS to too high of a
  // density.
  high[4]=3.0;
  high[5]=12.0;
  high[6]=4.5;
  // a2
  high[7]=0.3;
  // a4
  high[8]=0.08;

  return;
}

void generic_quarks::initial_point(ubvector &params) {
  params[0]=1.19;
  params[1]=-2.52;
  params[2]=0.188;
  params[3]=0.357;
  params[4]=1.86;
  params[5]=5.70;
  params[6]=2.29;
  params[7]=0.1907;
  params[8]=0.0796;
  model::initial_point(params);
  return;
}

void generic_quarks::compute_eos(const ubvector &params, int &ret,
				 ofstream &scr_out, model_data &dat) {

  ret=ix_success;

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from ubvector information
  se.comp=params[0];
  se.kprime=params[1];
  se.b=params[2]-se.a;
  se.gamma=params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.set_eos_table(dat.eos);
  cns.calc_eos();
  dat.eos->set_interp_type(itp_linear);

  dat.eos->add_constant("S",params[2]);
  dat.eos->add_constant("L",se.fesym_slope(0.16));

  // Double check that the table is non-empty (we have to do this
  // especially for the size_t index in the for loop below)
  if (dat.eos->get_nlines()==0) {
    O2SCL_ERR("Table empty in generic quarks.",exc_efailed);
  }

  // Transition between nuclear part and polytrope
  double ed1=params[4];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=dat.eos->interp("ed",ed1,"nb");
  double pr1=dat.eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,params[5]);
  
  // Store the transition pressure for later use
  dat.eos->add_constant("pr_pt",pr1);
  
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    ret=ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=dat.eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<dat.eos->get_nlines();i++) {
    if ((*dat.eos)["ed"][i]>ed1) {
      dat.eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    ret=ix_no_eos_table;
    return;
  }
  
  // Transition between polytrope and quarks
  double ed_trans=params[6];
  double pr_trans=coeff1*pow(ed_trans,params[5]);
  double nb_trans=nb1*pow(ed_trans/ed1,params[5]/(params[5]-1.0))*
    pow((ed_trans+pr_trans)/(ed1+pr1),1.0/(1.0-params[5]));

  // Store the transition pressure for later use
  dat.eos->add_constant("pr_q",pr_trans);

  // Add first polytrope to table
  for(double ed=ed1;ed<=ed_trans;ed+=0.05) {
    double pr=coeff1*pow(ed,params[5]);
    double nb=nb1*pow(ed/ed1,params[5]/(params[5]-1.0))*
      pow((ed+pr)/(ed1+pr1),1.0/(1.0-params[5]));
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (dat.eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    ret=ix_no_eos_table;
    return;
  }

  // Quark EOS parameters
  double a2=params[7];
  double a4=params[8];
  
  // Coefficients of quadratic
  double quad_b=a2/2.0/a4;
  double quad_c=-(ed_trans+pr_trans)/4.0/a4;

  // Roots of quadratic
  double musq1=(-quad_b+sqrt(quad_b*quad_b-4.0*quad_c))/2.0;
  double musq2=(-quad_b-sqrt(quad_b*quad_b-4.0*quad_c))/2.0;

  // Select correct root
  double mu_start=0.0;
  if (musq1>0.0) {
    if (musq2<0.0) {
      mu_start=sqrt(musq1);
    } else {
      if (musq1>musq2) {
	mu_start=sqrt(musq1);
      } else {
	mu_start=sqrt(musq2);
      }
    }
  } else {
    if (musq2>0.0) {
      mu_start=sqrt(musq2);
    } else {
      scr_out << "Rejected: Neither mu^2 is positive." << endl;
      ret=ix_eos_solve_failed;
      return;
    }
  }
    
  // Compute maximum chemical potential and bag constant
  double mu_end=10.0;
  if (mu_end<2.0*mu_start) mu_end=2.0*mu_start;
  double musq=mu_start*mu_start;
  double bag=ed_trans-a2*musq-3.0*a4*musq*musq;

  // Correction factor for baryon density. Selecting a curve which is
  // continuous in pressure and energy density isn't enough because
  // that only specifies the baryon density up to an arbitrary scale.
  // We need to select the EOS which has the right baryon density to
  // match with the rest of the EOS above. Instead of doing this, we
  // simply renormalize the baryon density after the fact, using the
  // baryon density from the polytrope at the transition between the
  // polytrope and the quark matter EOS.
  double nb_corr=nb_trans/((2.0*a2*mu_start+4.0*a4*mu_start*musq)/3.0);

  // Loop through quark EOS. We ensure to start slightly larger mu
  // than delta_mu to ensure we don't get two grid points with the
  // same energy density (this causes problems with interpolation
  // later).
  double delta_mu=(mu_end-mu_start)/500;
  for(double mu=mu_start+delta_mu;mu<=mu_end;mu+=delta_mu) {
    musq=mu*mu;
      
    // Ensure pressure is increasing
    double dPde=(a2+2.0*a4*musq)/(a2+6.0*a4*musq);
    if (dPde<0.0) {
      scr_out << "Rejected: dPdeps<0.0 in quarks." << endl;
      ret=ix_acausal;
      return;
    }

    double ed=bag+a2*musq+3.0*a4*musq*musq;
    double pr=-bag+a2*musq+a4*musq*musq;
    double nb=(2.0*a2*mu+4.0*a4*mu*musq)/3.0*nb_corr;
    double line[3]={ed,pr,nb};
    dat.eos->line_of_data(3,line);
  }

  return;
}

int quark_star::pressure(size_t nv, const ubvector &x, ubvector &y) {
  double pr, mu=x[0];
  
  double mu2=mu*mu, mu4=mu2*mu2;
  
  double a2=-3.0*(ms*ms-4.0*Delta*Delta)/4.0/o2scl_const::pi2;
  double a4=3.0*(1.0-c)/4.0/o2scl_const::pi2;
    
  pr=a4*mu4+a2*mu2-B;
  pr-=3.0*pow(ms,4.0)/8.0/o2scl_const::pi2*log(ms/2.0/mu);
    
  y[0]=pr;
  return 0;
}

double quark_star::pressure2(double mu) {
  double pr;
    
  double mu2=mu*mu, mu4=mu2*mu2;
    
  double a2=-3.0*(ms*ms-4.0*Delta*Delta)/4.0/o2scl_const::pi2;
  double a4=3.0*(1.0-c)/4.0/o2scl_const::pi2;
    
  pr=a4*mu4+a2*mu2-B;
  pr-=3.0*pow(ms,4.0)/8.0/o2scl_const::pi2*log(ms/2.0/mu);

  return pr;
}

void quark_star::get_param_info(std::vector<std::string> &names,
				std::vector<std::string> &units,
				ubvector &low, ubvector &high) {

  names={"B","c","Delta","ms"};

  units={"1/fm","","1/fm","1/fm"};
  
  low.resize(n_eos_params+nsd.n_sources);
  // B
  low[0]=-10.0;
  // c
  low[1]=0.0;
  // Delta
  low[2]=0.0;
  // ms
  low[3]=0.75;
    
  high.resize(n_eos_params+nsd.n_sources);
  // B
  high[0]=10.0;
  // c
  high[1]=0.4;
  // Delta
  high[2]=1.0;
  // ms
  high[3]=2.5;

  return;
}

void quark_star::initial_point(ubvector &params) {
  params[0]=0.2446;
  params[1]=0.0740;
  params[2]=0.00289;
  params[3]=0.0754;
  model::initial_point(params);
  return;
}

void quark_star::compute_eos(const ubvector &params, int &ret,
			     std::ofstream &scr_out, model_data &dat) {
  
  ret=ix_success;

  B=params[0];
  c=params[1];
  Delta=params[2];
  ms=params[3];

  // Compute chemical potential at zero pressure
  ubvector x(1), y(1);

  // First, go over a wide range of chemical potentials to get
  // a good initial guess
  double dmu=0.1, mu_0=dmu, mu_1=0.0;
  for(x[0]=7.0;x[0]>=dmu-1.0e-6;x[0]-=dmu) {
    pressure(1,x,y);
    double yhigh=y[0];
    x[0]-=dmu;
    pressure(1,x,y);
    double ylow=y[0];
    x[0]+=dmu;
    if (yhigh>0.0 && ylow<0.0) {
      mu_0=x[0];
      x[0]=0.0;
      mu_1=x[0]-dmu;
    }
  }

  if (mu_0<dmu+1.0e-6) {
    ret=ix_eos_solve_failed;
    scr_out << "No zero pressure solution." << std::endl;
    return;
  }

  // Then call the root finder
  x[0]=mu_0;
  mm_funct11 fmf=std::bind(std::mem_fn<int(size_t,const ubvector &,ubvector &)>
			   (&quark_star::pressure),
			   this,std::placeholders::_1,std::placeholders::_2,
			   std::placeholders::_3);
  gmh.err_nonconv=false;
  int solve_ret=gmh.msolve(1,x,fmf);
  if (solve_ret!=0) {
    scr_out << "Solver failed in qstar." << std::endl;
    ret=ix_eos_solve_failed;
    return;
  }
  mu_0=x[0];
  //grb.solve_bkt(mu_1,mu_0,fmf);
    
  dat.eos->clear_table();
  dat.eos->new_column("ed");
  dat.eos->new_column("pr");
  dat.eos->new_column("nb");
  dat.eos->set_unit("ed","1/fm^4");
  dat.eos->set_unit("pr","1/fm^4");
  dat.eos->set_unit("nb","1/fm^3");

  double ed_last=0.0;
  double mu_max=mu_0*2.0;
  double pi2=o2scl_const::pi2;

  for(double mu=mu_0;mu<mu_max;mu+=(mu_max-mu_0)/100.0) {
      
    double mu2=mu*mu, mu3=mu2*mu, mu4=mu2*mu2;
    double ms2=ms*ms, ms3=ms2*ms, ms4=ms2*ms2;
      
    double a2=-3.0*(ms2-4.0*Delta*Delta)/4.0/pi2;
    double a4=3.0*(1.0-c)/4.0/pi2;

    double pr=a4*mu4+a2*mu2-B-3.0*ms4/8.0/pi2*log(ms/2.0/mu)+
      0.375*ms4/pi2/mu;
    double dpdmu=3.0*mu3*(1-c)/pi2-1.5/pi2*mu*(ms2-4.0*Delta*Delta);
    double ed=-pr+dpdmu*mu;
    double nb=dpdmu/3.0;

    // Check that energy density is increasing
    if (ed<ed_last) {
      scr_out << "Energy density decreasing in quark_star." << std::endl;
      ret=ix_acausal;
      return;
    }

    // First point
    if (fabs(mu-mu_0)<1.0e-6) {

      // For baryon density computation
      nb_n1=nb;
      nb_e1=ed;

      // Double check that the energy per baryon at zero density
      // is less than iron
      if (ed/nb>931.0/o2scl_const::hc_mev_fm) {
	scr_out << "Not absolutely stable." << std::endl;
	ret=ix_param_mismatch;
	return;
      }

      // Fix the first point at exactly zero pressure
      double line[3]={ed,0.0,nb};
      dat.eos->line_of_data(3,line);
	
    } else {
	
      double line[3]={ed,pr,nb};
      dat.eos->line_of_data(3,line);

    }

    ed_last=ed;
  }

  return;
}

// --------------------------------------------------------------

qmc_neut::qmc_neut(settings &s, ns_data &n) :
  model(s,n) {
  
  rho0=0.16;

  // Set sigma for Gaussian distribution
  pdg.set_sigma(1.0);

  double ratio_data[11][3]={
    {5.067731e-01,8.448435e-01,1.508267e-02},
    {1.013546e+00,8.892959e-01,4.387342e-02},
    {1.520319e+00,9.280464e-01,5.807839e-02},
    {2.027092e+00,9.529515e-01,5.861344e-02},
    {2.533866e+00,9.682471e-01,5.491887e-02},
    {3.040639e+00,9.776334e-01,5.004978e-02},
    {3.547412e+00,9.832998e-01,4.514579e-02},
    {4.054185e+00,9.865683e-01,4.078075e-02},
    {4.560958e+00,9.883023e-01,3.734590e-02},
    {5.067731e+00,9.890893e-01,3.506759e-02},
    {5.574504e+00,9.934692e-01,3.632158e-02}
  };
  
  // Recast data
  ed_corr.resize(11);
  pres_corr.resize(11);
  pres_err.resize(11);
  for(size_t i=0;i<11;i++) {
    ed_corr[i]=ratio_data[i][0];
    pres_corr[i]=ratio_data[i][1];
    pres_err[i]=ratio_data[i][2];
  }

  // Set interpolation objects with columns from ratio matrix
  si.set(11,ed_corr,pres_corr,itp_linear);
  si_err.set(11,ed_corr,pres_err,itp_linear);
  
  rho_trans=0.48;

  this->n_eos_params=7;
  this->has_esym=true;
}

qmc_neut::~qmc_neut() {
}

void qmc_neut::get_param_info(std::vector<std::string> &names,
			      std::vector<std::string> &units,
			      ubvector &low, ubvector &high) {

  names={"a","alpha","b","beta","index1","trans1","index2"};

  units={"MeV","","MeV","","","1/fm^4",""};
  
  low.resize(n_eos_params+nsd.n_sources);
  low[0]=12.7;
  low[1]=0.48;
  low[2]=1.0;
  low[3]=2.1;
  low[4]=0.2;
  low[5]=2.0;
  low[6]=0.2;
    
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=13.3;
  high[1]=0.52;
  high[2]=5.0;
  high[3]=2.5;
  high[4]=4.0;
  high[5]=8.0;
  high[6]=4.0;
    
  // Go to the parent which takes care of the data-related
  // parameters
  model::get_param_info(names,units,low,high);

  return;
}

void qmc_neut::initial_point(ubvector &params) {

  params[0]=1.276936e+01;
  params[1]=5.043647e-01;
  params[2]=4.584098e+00;
  params[3]=2.323736e+00;
  params[4]=0.576;
  params[5]=4.60;
  params[6]=1.21;

  model::initial_point(params);
  return;
}

void qmc_neut::compute_eos(const ubvector &params, int &ret,
			   ofstream &scr_out, model_data &dat) {

  ret=ix_success;
  
  // Hack to start with a fresh table
  dat.eos->clear_table();
  dat.eos->line_of_names("ed pr");
  dat.eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=params[0];
  double alpha=params[1];
  double b=params[2];
  double beta=params[3];
    
  double L=3.0*(a*alpha+b*beta);
  dat.eos->add_constant("S",(a+b+16.0)/hc_mev_fm);
  dat.eos->add_constant("L",L/hc_mev_fm);

  double ed=0.0, pr=0.0;

  double gauss=pdg();
  if (fabs(gauss)>3.0) gauss=0.0;

  for(double rho=0.02;rho<rho_trans+0.001;rho+=0.01) {
    double rho1=rho/rho0;
    double rho1a=pow(rho1,alpha);
    double rho1b=pow(rho1,beta);
    double ene=a*rho1a+b*rho1b;
    ed=rho*(ene/hc_mev_fm+o2scl_settings.get_convert_units().convert
            ("kg","1/fm",o2scl_mks::mass_neutron));
    pr=rho*(a*alpha*rho1a+b*beta*rho1b)/hc_mev_fm;
      
    // Correct the pressure by a factor to correct for
    // neutron -> neutron star matter
    if (true) {
      if (ed<0.5) {
        pr*=0.8;
      } else if (ed<5.6) {
        double fact=gauss*si_err.eval(ed)+si.eval(ed);
        if (fact>1.0) fact=1.0;
        pr*=fact;
      }
    }
      
    double line[2]={ed,pr};
    dat.eos->line_of_data(2,line);
  }

  // Set values for the computation of the baryon density
  // from the last point of the QMC parameterization
  nb_n1=rho_trans;
  nb_e1=ed;

  // Check that the last Monte Carlo energy density
  // isn't greater than the transition between the
  // polytropes
  if (params[5]<ed) {
    scr_out << "First polytrope doesn't appear." << endl;
    ret=ix_param_mismatch;
    return;
  }

  // Polytropic index at higher densities
  double index1=params[4];

  // Compute coefficient given index
  double exp=1.0+1.0/index1;
  double coeff=pr/pow(ed,exp);

  // Compute stepsize in energy density
  double delta_ed=(params[5]-ed)/50.01;

  // Ensure first point is not the same as the last point
  // (this causes problems for the baryon density)
  ed+=delta_ed;

  // Add first polytrope to table
  for(;ed<params[5];ed+=delta_ed) {
    double line[2]={ed,coeff*pow(ed,exp)};
    dat.eos->line_of_data(2,line);
  }

  // Compute second coefficient given index
  double index2=params[6];
  double exp2=1.0+1.0/index2;

  // Compute coefficient
  double pr_last=coeff*pow(params[5],exp);
  double coeff2=pr_last/pow(params[5],exp2);

  // Add second polytrope to table
  delta_ed=(10.0-params[5])/50.01;
  for(ed=params[5];ed<10.0;ed+=delta_ed) {
    double line[2]={ed,coeff2*pow(ed,exp2)};
    dat.eos->line_of_data(2,line);
  }

  return;
}

// --------------------------------------------------------------

qmc_threep::qmc_threep(settings &s, ns_data &n) :
  model(s,n) {
  
  rho0=0.16;
  rho_trans=0.16;

  this->n_eos_params=9;
  this->has_esym=true;
}

qmc_threep::~qmc_threep() {
}

void qmc_threep::get_param_info(std::vector<std::string> &names,
				std::vector<std::string> &units,
				ubvector &low, ubvector &high) {

  names={"a","alpha","param_S","param_L","index1","trans1","index2","trans2",
	 "index3"};

  units={"MeV","","MeV","MeV","","1/fm^4","","1/fm^4",""};
  
  low.resize(n_eos_params+nsd.n_sources);
  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  low[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  low[1]=0.47;
  low[2]=29.5;
  low[3]=30.0;
  
  low[4]=0.2;
  low[5]=0.75;
  low[6]=0.2;
  low[7]=0.75;
  low[8]=0.2;
    
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=13.5;
  high[1]=0.53;
  high[2]=36.1;
  high[3]=70.0;

  high[4]=8.0;
  high[5]=8.0;
  high[6]=8.0;
  high[7]=8.0;
  high[8]=8.0;
    
  // Go to the parent which takes care of the data-related
  // parameters
  model::get_param_info(names,units,low,high);

  return;
}

void qmc_threep::initial_point(ubvector &params) {

  params[0]=13.0;
  params[1]=0.5;
  params[2]=32.0;
  params[3]=50.0;
  params[4]=0.5;
  params[5]=2.0;
  params[6]=0.5;
  params[7]=2.5;
  params[8]=1.0;

  model::initial_point(params);
  return;
}

void qmc_threep::compute_eos(const ubvector &params, int &ret,
			     ofstream &scr_out, model_data &dat) {

  bool debug=false;

  ret=ix_success;
  
  // Hack to start with a fresh table
  dat.eos->clear_table();
  dat.eos->line_of_names("ed pr");
  dat.eos->set_interp_type(itp_linear);
  //dat.eos->set_unit("ed","1/fm^4");
  //dat.eos->set_unit("pr","1/fm^4");
  //dat.eos->set_unit("nb","1/fm^3");
  //dat.eos->line_of_names("ed pr nb");

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=params[0];
  double alpha=params[1];
  double Stmp=params[2];
  double Ltmp=params[3];

  /*
    This is based on limits from two lines, as in Jim and I's EPJA
    review. In (S,L) space, the lower line is (29,0) to (35,55),
    and the upper line is (26.5,0) to (33.5,100)
  */
  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range: " << Stmp << " " << Ltmp << endl;
    scr_out << 9.17*Stmp-266.0 << " " << 14.3*Stmp-379.0 << endl;
    ret=ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    ret=ix_param_mismatch;
    return;
  }
  if (debug) {
    scr_out << "b=" << b << " beta=" << beta << endl;
  }

  dat.eos->add_constant("S",Stmp/hc_mev_fm);
  dat.eos->add_constant("L",Ltmp/hc_mev_fm);

  double index1=params[4];
  double exp1=1.0+1.0/index1;
  double trans1=params[5];
  double index2=params[6];
  double exp2=1.0+1.0/index2;
  double trans2=params[7];
  double index3=params[8];
  double exp3=1.0+1.0/index3;

  double ed=0.0, pr=0.0, ed_last=0.0, pr_last=0.0, nb_last=0.0;

  double rho;
  for(rho=0.02;rho<rho_trans;rho+=0.01) {
    double rho1=rho/rho0;
    double rho1a=pow(rho1,alpha);
    double rho1b=pow(rho1,beta);
    double ene=a*rho1a+b*rho1b;
    ed=rho*(ene/hc_mev_fm+o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron));
    pr=rho*(a*alpha*rho1a+b*beta*rho1b)/hc_mev_fm;
    
    double line[2]={ed,pr};
    dat.eos->line_of_data(2,line);
    // double line[3]={ed,pr,rho};
    // dat.eos->line_of_data(3,line);
    if (debug) scr_out << line[0] << " " << line[1] << endl;
    ed_last=ed;
    pr_last=pr;
    nb_last=rho;
  }

  // Set values for the computation of the baryon density
  // from the last point of the QMC parameterization
  nb_n1=nb_last;
  nb_e1=ed_last;

  // Check that the transition densities are ordered
  if (ed_last>trans1 || trans1>trans2) {
    scr_out << "Transition densities misordered." << endl;
    scr_out << ed_last << " " << trans1 << " " << trans2 << endl;
    ret=ix_param_mismatch;
    return;
  }

  // Compute coefficient given index
  double coeff1=pr_last/pow(ed_last,exp1);

  // Compute stepsize in energy density
  double delta_ed=(trans1-ed_last)/30.01;

  // Add first polytrope to table
  for(ed=ed_last+delta_ed;ed<trans1;ed+=delta_ed) {
    pr=coeff1*pow(ed,exp1);
    double line[2]={ed,pr};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_threep (2): " << line[0] << " "
	   << line[1] << endl;
      cerr << index1 << " " << exp1 << endl;
      cerr << coeff1 << " " << exp1 << " " << ed_last << " " << trans1 << endl;
      cerr << ed << " " << pr << " " << nb_n1 << " " << nb_e1 << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) scr_out << line[0] << " " << line[1] << endl;
    ed_last=ed;
    pr_last=pr;
  }

  // Compute second coefficient given index
  double coeff2=pr_last/pow(ed_last,exp2);

  // Add second polytrope to table
  delta_ed=(trans2-trans1)/20.01;
  for(ed=trans1;ed<trans2;ed+=delta_ed) {
    pr=coeff2*pow(ed,exp2);
    double line[2]={ed,pr};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_threep (3): " << line[0] << " "
	   << line[1] << endl;
      cerr << index2 << " " << exp2 << endl;
      cerr << trans1 << " " << trans2 << " " << coeff2 << " "
	   << pr_last << " " << ed_last << " " << exp2 << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) scr_out << line[0] << " " << line[1] << endl;
    ed_last=ed;
    pr_last=pr;
  }

  // Compute third coefficient given index
  double coeff3=pr_last/pow(ed_last,exp3);

  // Add third polytrope to table
  delta_ed=(10.0-trans2)/20.01;
  for(ed=trans2;ed<10.0;ed+=delta_ed) {
    pr=coeff3*pow(ed,exp3);
    double line[2]={ed,pr};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_threep (4): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) scr_out << line[0] << " " << line[1] << endl;
  }

  return;
}

// --------------------------------------------------------------

qmc_fixp::qmc_fixp(settings &s, ns_data &n) :
  model(s,n) {
  
  nb0=0.16;
  nb_trans=0.16;

  ed1=2.0;
  ed2=3.0;
  ed3=5.0;
  ed4=7.0;

  this->n_eos_params=8;
  this->has_esym=true;
}

qmc_fixp::~qmc_fixp() {
}

void qmc_fixp::get_param_info(std::vector<std::string> &names,
			      std::vector<std::string> &units,
			      ubvector &low, ubvector &high) {

  names={"a","alpha","param_S","param_L","pres1","pres2","pres3","pres4"};

  units={"MeV","","MeV","MeV","1/fm^4","1/fm^4","1/fm^4","1/fm^4"};
  
  low.resize(n_eos_params+nsd.n_sources);
  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  low[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  low[1]=0.47;
  low[2]=29.5;
  low[3]=30.0;
  
  low[4]=0.0;
  low[5]=0.0;
  low[6]=0.0;
  low[7]=0.0;
    
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=13.5;
  high[1]=0.53;
  high[2]=36.1;
  high[3]=70.0;

  // These parameters are limited by causality, but if the user
  // changes the values of ed1, ed2, ed3, and ed4, then the upper
  // limits change accordingly. To make things easier, we just choose
  // relatively large values for these upper limits for now.
  high[4]=0.3;
  high[5]=1.5;
  high[6]=2.5;
  high[7]=2.5;
    
  // Go to the parent which takes care of the data-related
  // parameters
  model::get_param_info(names,units,low,high);

  return;
}

void qmc_fixp::initial_point(ubvector &params) {

  params[0]=1.276936e+01;
  params[1]=5.043647e-01;
  params[2]=30.0;
  params[3]=40.0;
  params[4]=0.014;
  params[5]=0.74;
  params[6]=0.60;
  params[7]=1.84;

  model::initial_point(params);
  return;
}

void qmc_fixp::compute_eos(const ubvector &params, int &ret,
			   ofstream &scr_out, model_data &dat) {

  ret=ix_success;
  bool debug=false;
  
  // Hack to start with a fresh table
  dat.eos->clear_table();
  dat.eos->line_of_names("ed pr");
  dat.eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=params[0];
  double alpha=params[1];
  double Stmp=params[2];
  double Ltmp=params[3];

  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range." << endl;
    ret=ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta || b<0.5) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    ret=ix_param_mismatch;
    return;
  }
  
  dat.eos->add_constant("S",Stmp/hc_mev_fm);
  dat.eos->add_constant("L",Ltmp/hc_mev_fm);

  double ed=0.0, pr=0.0, ed_trans=0.0, pr_trans=0.0;

  if (debug) {
    cout.setf(ios::scientific);
    cout << "a,alpha,b,beta: "
	 << a << " " << alpha << " " << b << " " << beta << endl;
    cout << endl;
    cout << "EOS below saturation:" << endl;
    cout << "ed           pr" << endl;
  }
  for(double nb=0.02;nb<nb_trans+1.0e-6;nb+=0.001) {
    double nb1=nb/nb0;
    double nb1a=pow(nb1,alpha);
    double nb1b=pow(nb1,beta);
    double ene=a*nb1a+b*nb1b;
    ed=nb*(ene/hc_mev_fm+o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_neutron));
    pr=nb*(a*alpha*nb1a+b*beta*nb1b)/hc_mev_fm;
      
    double line[2]={ed,pr};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_fixp (4): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << ed << " " << pr << endl;
    ed_trans=ed;
    pr_trans=pr;
  }
  if (debug) cout << endl;

  // Set values for the computation of the baryon density
  // from the last point of the QMC parameterization
  nb_n1=nb_trans;
  nb_e1=ed_trans;

  if (ed_trans>ed1) {
    scr_out << "Transition densities misordered." << endl;
    scr_out << ed_trans << " " << ed1 << endl;
    ret=ix_param_mismatch;
    return;
  }

  // Compute pressures on grid, p1 is the pressure at ed1
  double p1=pr_trans+params[4];
  // Variable p2 is the pressure at ed2
  double p2=p1+params[5];
  // Variable p3 is the pressure at ed3
  double p3=p2+params[6];
  
  // Add 1st high-density EOS
  if (debug) {
    cout << "First line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  double delta_ed=(ed1-ed_trans)/20.0;
  for(double ed=ed_trans+delta_ed;ed<ed1-1.0e-4;ed+=delta_ed) {
    double line[2]={ed,pr_trans+params[4]*(ed-ed_trans)/(ed1-ed_trans)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_fixp (5): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 2nd high-density EOS
  if (debug) {
    cout << "Second line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed1;ed<ed2-1.0e-4;ed+=(ed2-ed1)/10.0) {
    double line[2]={ed,p1+params[5]*(ed-ed1)/(ed2-ed1)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_fixp (6): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 3rd high-density EOS
  if (debug) {
    cout << "Third line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed2;ed<ed3-1.0e-4;ed+=(ed3-ed2)/10.0) {
    double line[2]={ed,p2+params[6]*(ed-ed2)/(ed3-ed2)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_fixp (7): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 4th high-density EOS
  if (debug) {
    cout << "Fourth line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed3;ed<10.0-1.0e-4;ed+=(ed4-ed3)/10.0) {
    double line[2]={ed,p3+params[7]*(ed-ed3)/(ed4-ed3)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_fixp (8): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) {
    cout << "Exiting since debug in qmc_fixp::compute_eos() is true."
	 << endl;
    exit(-1);
  }

  return;
}

// --------------------------------------------------------------

qmc_twolines::qmc_twolines(settings &s, ns_data &n) :model(s,n) {
  nb0=0.16;
  nb_trans=0.16;
}

qmc_twolines::~qmc_twolines() {
}

void qmc_twolines::get_param_info(std::vector<std::string> &names,
				  std::vector<std::string> &units,
				  ubvector &low, ubvector &high) {

  names={"a","alpha","param_S","param_L","pres1","ed1","pres2","ed2"};
  
  units={"MeV","","MeV","MeV","1/fm^4","1/fm^4","1/fm^4","1/fm^4"};

  low.resize(n_eos_params+nsd.n_sources);
  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  low[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  low[1]=0.47;
  low[2]=29.5;
  low[3]=30.0;
  
  low[4]=0.0;
  low[5]=0.0;
  low[6]=0.0;
  low[7]=0.0;
    
  high.resize(n_eos_params+nsd.n_sources);
  high[0]=13.5;
  high[1]=0.53;
  high[2]=36.1;
  high[3]=70.0;

  high[4]=0.3;
  high[5]=1.5;
  high[6]=2.5;
  high[7]=2.5;
      
  return;
}

void qmc_twolines::initial_point(ubvector &params) {

  params[0]=1.276936e+01;
  params[1]=5.043647e-01;
  params[2]=30.0;
  params[3]=40.0;
  params[4]=0.3;
  params[5]=1.0;
  params[6]=0.8;
  params[7]=1.84;

  model::initial_point(params);
  return;
}

void qmc_twolines::compute_eos(const ubvector &params, int &ret,
			       ofstream &scr_out, model_data &dat) {

  ret=ix_success;
  bool debug=false;
  
  // Hack to start with a fresh table
  dat.eos->clear_table();
  dat.eos->line_of_names("ed pr");
  dat.eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=params[0];
  double alpha=params[1];
  double Stmp=params[2];
  double Ltmp=params[3];

  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range." << endl;
    ret=ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta || b<0.5) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    ret=ix_param_mismatch;
    return;
  }
  
  dat.eos->add_constant("S",Stmp/hc_mev_fm);
  dat.eos->add_constant("L",Ltmp/hc_mev_fm);

  double ed=0.0, pr=0.0, ed_last=0.0, nb_last=0.0, pr_last=0.0;

  double ed1=params[5];
  double ed2=params[7];
  if (ed1>ed2) {
    scr_out << "Transition densities misordered." << endl;
    scr_out << ed1 << " " << ed2 << endl;
    ret=ix_param_mismatch;
    return;
  }

  if (debug) {
    cout.setf(ios::scientific);
    cout << a << " " << alpha << " " << b << " " << beta << endl;
    cout << endl;
  }
  bool done=false;
  for(double nb=0.02;done==false;nb+=0.001) {
    double nb1=nb/nb0;
    double nb1a=pow(nb1,alpha);
    double nb1b=pow(nb1,beta);
    double ene=a*nb1a+b*nb1b;
    ed=nb*(ene/hc_mev_fm+o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_neutron));
    pr=nb*(a*alpha*nb1a+b*beta*nb1b)/hc_mev_fm;

    if (ed>ed1) {
      done=true;
    } else {
      double line[2]={ed,pr};
      if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
	cerr << "Problem in qmc_twolines (4): " << line[0] << " "
	     << line[1] << endl;
	exit(-1);
      }
      dat.eos->line_of_data(2,line);
      if (debug) cout << ed << " " << pr << endl;
      ed_last=ed;
      pr_last=pr;
      nb_last=nb;
      
      // Set values for the computation of the baryon density
      // from the last point of the QMC parameterization
      nb_n1=nb;
      nb_e1=ed;
    }
  }
  if (debug) cout << endl;

  // Compute pressures at end points
  double p2=pr_last+params[4];
  double p3=p2+params[6];
  
  // Add 1st high-density EOS
  if (debug) {
    cout << "First line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  double delta_ed=(ed2-ed_last)/100.0;
  for(double ed=ed_last+delta_ed;ed<ed2-1.0e-4;ed+=delta_ed) {
    double line[2]={ed,pr_last+params[4]*(ed-ed_last)/(ed2-ed_last)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_twolines (5): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;
  
  // Add 2nd high-density EOS
  if (debug) {
    cout << "Second segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed2;ed<10.0-1.0e-4;ed+=(10.0-ed2)/20.0) {
    double line[2]={ed,p2+params[6]*(ed-ed2)/(10.0-ed2)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in qmc_twolines (6): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    dat.eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) {
    cout << endl;
    exit(-1);
  }

  return;
}

