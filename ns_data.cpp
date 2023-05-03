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
#include "ns_data.h"

#include <o2scl/hdf_io.h>

#ifdef BAMR_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace bamr;

void ns_data::load_mc(std::ostream &scr_out, int mpi_size, int mpi_rank,
		      std::shared_ptr<settings> set) {

  double tot, max;

  std::string name;

  if (source_names.size()!=source_fnames.size() ||
      source_names.size()!=init_mass_fracs.size()) {
    O2SCL_ERR("Incorrect input data sizes.",o2scl::exc_esanity);
  }
  
  // If requested, add the LIGO data
  if (set->inc_ligo) {
    hdf_file hfx;
    hfx.open("data/ligo/gw170817_kde.o2");
    hdf_input(hfx,ligo_data_table,name);
    hfx.close();
    hfx.open("data/ligo/190425_mass1.o2");
    hdf_input(hfx,ligo_gw19,name);
    hfx.close();
  }
  
  if (n_sources>0) {

    if (set->verbose>=2) {
      cout << "bamr: Loading " << n_sources << " data files " 
	   << "with rank " << mpi_rank << " and size " 
	   << mpi_size << endl;
    }
    
    source_tables.resize(n_sources);
    if (source_fnames_alt.size()>0) {
      source_tables_alt.resize(n_sources);
    }
    
#ifdef BAMR_MPI

    int buffer=0, tag=0;
    
    // Choose which file to read first for this rank
    int filestart=0;
    if (mpi_rank>mpi_size-((int)n_sources) && mpi_rank>0) {
      filestart=mpi_size-mpi_rank;
    }
    if (set->mpi_load_debug) {
      scr_out << "Variable 'filestart' is " << filestart << " for rank "
	      << mpi_rank << "." << std::endl;
    }

    for(int k=0;k<((int)n_sources);k++) {
      
      // For k=0, we choose some ranks to begin reading, the others
      // have to wait. For k>=1, all ranks have to wait their turn.
      if (mpi_size>1 && (k>0 || (mpi_rank>0 &&
				 mpi_rank<=mpi_size-((int)n_sources)))) {
	int prev=mpi_rank-1;
	if (prev<0) prev+=mpi_size;
	if (set->mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " waiting for " 
		  << prev << "." << std::endl;
	}
	MPI_Recv(&buffer,1,MPI_INT,prev,tag,MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
      }
      
      // Determine which file to read next
      int file=filestart+k;
      if (file>=((int)n_sources)) file-=n_sources;

      if (set->mpi_load_debug) {
	scr_out << "Rank " << mpi_rank << " reading file " 
		<< file << "." << std::endl;
      } else {

	o2scl_hdf::hdf_file hf;
	hf.open(source_fnames[file]);
	if (table_names[file].length()>0) {
	  hdf_input(hf,source_tables[file],table_names[file]);
	} else {
          std::string name_temp;
	  hdf_input_n(hf,source_tables[file],name_temp);
	}
	source_tables[file].set_interp_type(o2scl::itp_linear);
	hf.close();
      
	if (source_fnames_alt.size()>0) {
	  o2scl_hdf::hdf_file hf2;
	  hf2.open(source_fnames_alt[file]);
	  if (table_names[file].length()>0) {
            hdf_input(hf2,source_tables_alt[file],table_names[file]);
	  } else {
            std::string name_temp;
	    hdf_input_n(hf2,source_tables_alt[file],name_temp);
	  }
	  source_tables_alt[file].set_interp_type(o2scl::itp_linear);
	  hf2.close();
	}
	
      }
      
      // Send a message, unless the rank is the last one to read a
      // file.
      if (mpi_size>1 && (k<((int)n_sources)-1 ||
                         mpi_rank<mpi_size-((int)n_sources))) {
	int next=mpi_rank+1;
	if (next>=mpi_size) next-=mpi_size;
	if (set->mpi_load_debug) {
	  scr_out << "Rank " << mpi_rank << " sending to " 
		  << next << "." << std::endl;
	}
	MPI_Send(&buffer,1,MPI_INT,next,tag,MPI_COMM_WORLD);
      }
      
    }

    if (set->mpi_load_debug==true) {
      scr_out << "Exiting since mpi_load_debug is true." << endl;
      cout << "Exiting since mpi_load_debug is true." << endl;
      // Ensure all the debug information is output before
      // we call exit();
      MPI_Barrier(MPI_COMM_WORLD);
      exit(-1);
    }
    
#else
    
    for(size_t k=0;k<n_sources;k++) {
      
      hdf_file hf;
      hf.open(source_fnames[k]);
      if (table_names[k].length()>0) {
	hdf_input(hf,source_tables[k],table_names[k]);
      } else {
        std::string name_temp;
	hdf_input_n(hf,source_tables[k],name_temp);
      }
      source_tables[k].set_interp_type(o2scl::itp_linear);
      hf.close();

      if (source_fnames_alt.size()>0) {
	o2scl_hdf::hdf_file hf2;
	hf2.open(source_fnames_alt[k]);
	if (table_names[k].length()>0) {
	  hdf_input(hf2,source_tables_alt[k],table_names[k]);
	} else {
          std::string name_temp;
          hdf_input_n(hf2,source_tables_alt[k],name_temp);
	}
	source_tables_alt[k].set_interp_type(o2scl::itp_linear);
	hf2.close();
      }
      
    }
    
#endif

    scr_out << "\nInput data files: " << std::endl;
    
    if (set->norm_max) {
      scr_out << "Normalizing maximum probability to 1." << std::endl;
    } else {
      scr_out << "Normalizing integral of distribution to 1." << std::endl;
    }

    scr_out << "File                          name   total        "
	    << "max          P(10,1.4)" << std::endl;

    for(size_t k=0;k<n_sources;k++) {
      
      // Update input limits
      if (k==0) {
	set->in_r_min=source_tables[k].get_grid_x(0);
	set->in_r_max=source_tables[k].get_grid_x
	  (source_tables[k].get_nx()-1);
	set->in_m_min=source_tables[k].get_grid_y(0);
	set->in_m_max=source_tables[k].get_grid_y
	  (source_tables[k].get_ny()-1);
      } else {
	if (set->in_r_min>source_tables[k].get_grid_x(0)) {
	  set->in_r_min=source_tables[k].get_grid_x(0);
	}
	if (set->in_r_max<source_tables[k].get_grid_x
	    (source_tables[k].get_nx()-1)) {
	  set->in_r_max=source_tables[k].get_grid_x
	    (source_tables[k].get_nx()-1);
	}
	if (set->in_m_min>source_tables[k].get_grid_y(0)) {
	  set->in_m_min= source_tables[k].get_grid_y(0);
	}
	if (set->in_m_max<source_tables[k].get_grid_y
	    (source_tables[k].get_ny()-1)) {
	  set->in_m_max= source_tables[k].get_grid_y
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
	  if (set->norm_max) {
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

      if (source_fnames_alt.size()>0) {
	tot=0.0;
	max=0.0;
	for(size_t i=0;i<source_tables_alt[k].get_nx();i++) {
	  for(size_t j=0;j<source_tables_alt[k].get_ny();j++) {
	    tot+=source_tables_alt[k].get(i,j,slice_names[k]);
	    if (source_tables_alt[k].get(i,j,slice_names[k])>max) {
	      max=source_tables_alt[k].get(i,j,slice_names[k]);
	    }
	  }
	}
	for(size_t i=0;i<source_tables_alt[k].get_nx();i++) {
	  for(size_t j=0;j<source_tables_alt[k].get_ny();j++) {
	    if (set->norm_max) {
	      source_tables_alt[k].set
		(i,j,slice_names[k],source_tables_alt[k].get
		 (i,j,slice_names[k])/max);
	    } else {
	      source_tables_alt[k].set
		(i,j,slice_names[k],source_tables_alt[k].get
		 (i,j,slice_names[k])/tot);
	    }
	  }
	}
      }	    

      if (set->debug_load) {
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
      scr_out << source_fnames[k].substr(0,29) << " ";
      scr_out.width(6);
      scr_out << source_names[k] << " " << tot << " " << max << " ";
      scr_out.unsetf(std::ios::left);
      scr_out << source_tables[k].interp(10.0,1.4,slice_names[k])
	      << std::endl;
      
      if (source_fnames_alt.size()>0) {
	scr_out.setf(std::ios::left);
	scr_out.width(29);
	scr_out << source_fnames_alt[k].substr(0,29) << " ";
	scr_out.width(6);
	scr_out << source_names[k] << " " << tot << " " << max << " ";
	scr_out.unsetf(std::ios::left);
	scr_out << source_tables_alt[k].interp(10.0,1.4,slice_names[k])
		<< std::endl;
      }
      
    }
    
    scr_out << std::endl;
  } else {
    if (set->verbose>=2) {
      cout << "bamr: no sources." << endl;
    }
  }    

  if (set->in_m_min<set->min_mass) set->in_m_min=set->min_mass;
  
  scr_out << "M limits: (" 
	  << set->in_m_min << "," << set->in_m_max << ")" << std::endl;
  scr_out << "R limits: ("
	  << set->in_r_min << "," << set->in_r_max << ")" << std::endl;

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

  cout << "Xere." << endl;

  return 0;
}
    
int ns_data::add_data_alt(std::vector<std::string> &sv, bool itive_com) {
      
  if (sv.size()<6) {
    std::cerr << "Not enough arguments given to 'add-data-alt'." << std::endl;
    return o2scl::exc_efailed;
  }

  if (source_fnames.size()!=source_fnames_alt.size()) {
    std::cerr << "When using 'add-data-alt' all objects must have an "
	      << "alternate data set." << std::endl;
    return o2scl::exc_efailed;
  }

  source_names.push_back(sv[1]);
  source_fnames.push_back(sv[2]);
  source_fnames_alt.push_back(sv[3]);
  slice_names.push_back(sv[4]);
  init_mass_fracs.push_back(o2scl::stod(sv[5]));
  if (sv.size()>=7) {
    table_names.push_back(sv[6]);
  } else {
    table_names.push_back("");
  }

  n_sources++;

  return 0;
}

void ns_data::data_params(std::vector<std::string> &names,
                          std::vector<std::string> &units,
                          std::vector<double> &low,
                          std::vector<double> &high,
                          std::shared_ptr<settings> set) {

  if (set->inc_ligo) {
    names.push_back("M_chirp_det");
    names.push_back("q");
    names.push_back("z_cdf");
    names.push_back("m1_gw19");
    units.push_back("Msun");
    units.push_back("");
    units.push_back("");
    units.push_back("Msun");
    low.push_back(1.1971);
    low.push_back(0.0);
    low.push_back(0.0);
    low.push_back(0.5);
    high.push_back(1.1979);
    high.push_back(1.0);
    high.push_back(1.0);
    high.push_back(2.5);
  }
  
  for(size_t i=0;i<n_sources;i++) {
    names.push_back("mf_"+source_names[i]);
    units.push_back("");
    low.push_back(0.0);
    high.push_back(1.0);
  }

  return;
}

void ns_data::initial_point(std::shared_ptr<settings> set,
                            std::vector<double> &init) {

  if (set->inc_ligo) {
    init.push_back(1.1975);
    init.push_back(0.6);
    init.push_back(0.5);
    init.push_back(1.82);
  }
  for(size_t i=0;i<n_sources;i++) {
    init.push_back(0.7);
  }
  
  return;
}

