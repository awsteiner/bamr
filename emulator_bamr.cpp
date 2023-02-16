/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2022, Mohammad Al-Mamun, Mahmudul Hasan Anik, 
  and Andrew W. Steiner
  
  This file is part of Bamr.
  
  Bamr is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published 
  by the Free Software Foundation; either version 3 of the License, 
  or (at your option) any later version.
  
  Bamr is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Bamr. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include "emulator_bamr.h"
#include <o2scl/vector.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf; 
using namespace o2scl_const; 
using namespace bamr;

/** \brief The interpolation estimate objects
 */

emulator_bamr::emulator_bamr() {
  
  mpi_rank=0;
  mpi_size=1;
  
#ifdef BAMR_MPI    
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
#endif
  
  unsigned long int seed=time(0);
  seed*=(mpi_rank+1);
  r.set_seed(seed);
}

void emulator_bamr::train(o2scl::table_units<> &tab_train,
                          o2scl::vec_index &pvii,
                          o2scl::vec_index &dvii,
                          bamr_class *bcpi,
                          std::ofstream *sopi) {

  bcp=bcpi;
  sop=sopi;
  dvi=dvii;
  pvi=pvii;
    
  list.clear();
  for(size_t j=0;j<pvi.size();j++) {
    list.push_back(pvi[j]);
  }
  /* 
     M_i, R_i, wgt_i for n_sources [2: ix_mr_outside]
     M_max=dat.mvsr.max("gm") [10: ix_small_max]
     m_max2 (add col to table) [10: ix_small_max]
     ed_max 
     cs2_i for i in (0,99) [13: ix_acausal]
     dpdM [exclude since ix_infinite never happens]
     ed_last (add col to table) [which one? there are 3 in NP]
  */
  list.push_back("a");
  list.push_back("alpha");
  list.push_back("param_S");
  list.push_back("param_L");
  list.push_back("exp1");
  list.push_back("trans1");
  list.push_back("exp2");
  list.push_back("trans2");
  list.push_back("exp3");
  list.push_back("M_chirp_det");
  list.push_back("q");
  list.push_back("z_cdf");
  list.push_back("mf_6304");
  list.push_back("mf_6397");
  list.push_back("mf_M13");
  list.push_back("mf_M28");
  list.push_back("mf_M30");
  list.push_back("mf_wCen");
  list.push_back("mf_X7");
  list.push_back("mf_1810b");
  list.push_back("mf_1724b");
  list.push_back("mf_1702");
  list.push_back("mf_0030");
  list.push_back("mean_NS");
  list.push_back("log10_width_NS");
  list.push_back("skewness_NS");
  list.push_back("mean_WD");
  list.push_back("log10_width_WD");
  list.push_back("skewness_WD");
  list.push_back("mean_LMS");
  list.push_back("log10_width_LMS");
  list.push_back("skewness_LMS");
  list.push_back("M_J0453p");
  list.push_back("M_J0453c");
  list.push_back("M_J1906p");
  list.push_back("M_J1906c");
  list.push_back("M_B1534p");
  list.push_back("M_B1534c");
  list.push_back("M_B1913p");
  list.push_back("M_B1913c");
  list.push_back("M_B2127p");
  list.push_back("M_B2127c");
  list.push_back("M_J0737A");
  list.push_back("M_J0737B");
  list.push_back("M_J1756p");
  list.push_back("M_J1756c");
  list.push_back("M_J1807p");
  list.push_back("M_J1807c");
  list.push_back("M_J1518p");
  list.push_back("M_J1518c");
  list.push_back("M_J1811p");
  list.push_back("M_J1811c");
  list.push_back("M_J1829p");
  list.push_back("M_J1829c");
  list.push_back("M_J2045");
  list.push_back("M_J2053");
  list.push_back("M_J1713");
  list.push_back("M_B1855");
  list.push_back("M_J0751");
  list.push_back("M_J1141");
  list.push_back("M_J1738");
  list.push_back("M_J1614");
  list.push_back("M_J0348");
  list.push_back("M_J2222");
  list.push_back("M_J2234");
  list.push_back("M_J1949");
  list.push_back("M_J1012");
  list.push_back("M_J0437");
  list.push_back("M_J1909");
  list.push_back("M_J1802");
  list.push_back("M_J1911");
  list.push_back("M_J2043");
  list.push_back("M_J0337");
  list.push_back("M_J1946");
  list.push_back("M_J1918");
  list.push_back("M_J1600");
  list.push_back("M_J0621");
  list.push_back("M_B2303");
  list.push_back("M_J0024");
  list.push_back("M_J0514");
  list.push_back("M_B1516");
  list.push_back("M_J1748I");
  list.push_back("M_J1748J");
  list.push_back("M_B1802");
  list.push_back("M_B1911");
  list.push_back("M_J0740");
  list.push_back("M_CygX2");
  list.push_back("M_XTEJ2123");
  list.push_back("M_4U1822");
  list.push_back("M_HerX1");
  list.push_back("M_2S0921");
  
  list.push_back("log_wgt");
  
  // list.push_back("log_wgt_NS");
  // list.push_back("log_wgt_WD");
  // list.push_back("log_wgt_LMS");
  // list.push_back("M_max");
  // list.push_back("R_43");

  cout << "Training column list (size " << list.size() << "): "; 
  o2scl::vector_out(std::cout,list,true);

  np=pvi.size();
  nout=list.size()-pvi.size();

  table.clear();
  for(size_t j=0;j<list.size();j++) {
    table.new_column(list[j]);
  }
  
  for(size_t k=0;k<files.size();k++) {
    
#ifdef BAMR_MPI    
    // Ensure that multiple MPI ranks are not writing to the 
    // filesystem at the same time
    int tag=0, buffer=0;
    if (mpi_size>1 && mpi_rank>=1) {
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
               tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#endif
    
    std::cout << "Rank " << mpi_rank << " reading training table "
              << "with index " << k << " from file "
              << files[k] << std::endl;
    table_units<> tab_k;
    hdf_file hf;
    hf.open(files[k]);
    hdf_input(hf,tab_k);
    hf.close();
    cout << "Rank " << mpi_rank << ": table has "
         << tab_k.get_nlines() << " lines." << endl;
    
#ifdef BAMR_MPI
    // Send a message to the next MPI rank
    if (mpi_size>1 && mpi_rank<mpi_size-1) {
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
               tag,MPI_COMM_WORLD);
    }
#endif

    for(size_t j=0;j<list.size();j++) {
      if (tab_k.is_column(list[j])==false) {
        cout << "Table in file " << files[k] << " does not have "
             << "column " << list[j] << "." << endl;
        exit(-1);
      }
    }

    // Delete any rows with a small log_wgt, which are emulated, or
    // have a value of "mult" which is zero

    if (tab_k.is_column("emulated")) {
      if (tab_k.is_column("mult")) {
        tab_k.delete_rows_func(((string)"emulated>0.5 || ")+
                               "log_wgt<(-700) || abs(mult)<0.5");
      } else {
        tab_k.delete_rows_func("emulated>0.5 || log_wgt<(-700)");
      }
    } else {
      if (tab_k.is_column("mult")) {
        tab_k.delete_rows_func("log_wgt<(-700) || abs(mult)<0.5");
      } else {
        tab_k.delete_rows_func("log_wgt<(-700)");
      }
    }
    cout << "Rank " << mpi_rank << ": table now has "
         << tab_k.get_nlines() << " lines." << endl;

    // Add this table to the combined table
    for(size_t i=0;i<tab_k.get_nlines();i++) {
      vector<double> line;
      for(size_t j=0;j<list.size();j++) {
        line.push_back(tab_k.get(list[j],i));
      }
      table.line_of_data(line.size(),line);
    }
    
  }

  cout << "Rank " << mpi_rank << ": combined table has "
       << table.get_nlines() << " lines." << endl;
  
  // Go through the combined table and delete nearly equal rows
  double tol_abs=1.0e-12;
  double tol_rel=1.0e-12;

  vector<size_t> row_list;
  
  for(size_t i=0;i<table.get_nlines();i++) {
    
    if (i%500==499) {
      std::cout << "Rank " << mpi_rank
                << " progress: i+1= " << i+1 << " of "
                << table.get_nlines() << endl;
    }
    // Check for duplicates
    for(size_t j=i+1;j<table.get_nlines();j++) {
      bool match=true;
      for(size_t k=0;k<np && match==true;k++) {
        if (fabs(table.get(list[k],i))>tol_abs ||
            fabs(table.get(list[k],j))>tol_abs) {
          if (fabs(table.get(list[k],i)-table.get(list[k],j))/
              fabs(table.get(list[k],i)+table.get(list[k],j))>tol_rel) {
            match=false;
          }
        }
      }
      if (match==true) {
        row_list.push_back(j);
        if (false) {
          std::cout << "Match between rows " << i << " and " << j
                    << " " << table.get(list[0],i)
                    << " " << table.get(list[0],j)
                    << " " << table.get(list[1],i)
                    << " " << table.get(list[1],j)
                    << std::endl;
        }
      }
    }
  }
  
  table.delete_rows_list(row_list);

  cout << "Rank " << mpi_rank << ": combined table now has "
       << table.get_nlines() << " lines." << endl;
  table.summary(&std::cout);

#ifdef BAMR_MPI    
  // Ensure that multiple MPI ranks are not writing to the 
  // filesystem at the same time
  int tag=0, buffer=0;
  if (mpi_size>1 && mpi_rank>=1) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
             tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
#endif
  
  std::cout << "Rank " << mpi_rank << " writing final training table."
            << std::endl;
  hdf_file hf2;
  hf2.open_or_create(((string)"train_")+o2scl::itos(mpi_rank)+"_out");
  hdf_output(hf2,table,"train");
  hf2.close();
  
#ifdef BAMR_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
             tag,MPI_COMM_WORLD);
  }
#endif
  
  em1.set(np,nout,0,table,list);
  
  // module (*.py), class_name (emu), train_func, point_func, n_pars, file
  
  // em3.set("anik","emu","train","eval","log_wgt",np,table,list);

  cout << "Rank " << mpi_rank << " done setting emulator." << endl;

  return;
}
  
int emulator_bamr::eval(size_t n, const ubvector &p, double &log_wgt,
                        model_data &dat) {

  if (false) {
    std::cout << "p: ";
    o2scl::vector_out(std::cout,p,true);
  }
    
  // Show that we can use the o2scl emulator or the full function
  std::vector<double> x(nout);
  double log_wgt_unc;
  std::vector<double> x_unc(nout);
  int em1_ret=em1.eval_unc(n,p,log_wgt,log_wgt_unc,x,x_unc);

  if (false) {
    std::cout << "x: ";
    o2scl::vector_out(std::cout,x,true);
  }
  std::cout << "log_wgt, unc, Mns_max: " << em1_ret << " " << log_wgt << " "
            << log_wgt_unc << " " << x[8] << std::endl;

  /*
  // Clear the dat array
  for(size_t i=0;i<ndat;i++) {
  dat[i]=0.0;
  }
    
  // Translate the emulated data into the 'dat' array used
  // by the data_eval point function. Skip k=0 because
  // it's already stored in log_wgt
  for(size_t k=1;k<nout;k++) {
  if (false) {
  std::cout << "Mapping: " << k << " " << list[np+k] << std::endl;
  }
  dat[dvi[list[np+k]]]=x[k];
  }

  dat[dvi["emulated"]]=1.0;
  dat[dvi["log_wgt_unc"]]=log_wgt_unc;
  */

  if (false) {
    std::cout << "p: ";
    o2scl::vector_out(std::cout,p,true);
  }

  double xrand=r.random();
  if (log_wgt_unc>1.0e-2 && 
      (log_wgt+2*log_wgt_unc>-40.0 || x[8]>1.947 || xrand<0.02)) {
    if (log_wgt+2*log_wgt_unc>-40.0) {
      std::cout << "High log_wgt." << std::endl;
    }
    if (xrand<0.02) {
      std::cout << "Random." << std::endl;
    }
    if (x[8]>1.947) {
      std::cout << "High M_max." << std::endl;
    }
      
    double log_wgt_old=log_wgt;

    int iret=bcp->compute_point(p,*sop,log_wgt,dat);

    //dat[dvi["log_wgt_unc"]]=fabs(log_wgt_old-log_wgt);
      
    std::cout << "log_wgt_old, log_wgt: " << log_wgt_old << " "
              << log_wgt << std::endl;
      
    if (false) {
      std::cout << "Character:" << std::endl;
      char ch;
      std::cin >> ch;
    }

    return iret;
  }
    
  return 0;
}

int emulator_bamr::eval_unc(size_t n, const ubvector &p, double &log_wgt,
                            double &lw_unc, model_data &dat,
                            model_data &dat_unc) {
  return eval(n,p,log_wgt,dat);
}

