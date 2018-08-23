/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  hdf_file hf;
  std::string name;
  
  {
    cout << "test_debug_eos: " << endl;
    hf.open("data_temp/debug_eos.o2");
    table_units<> eos;
    name="eos";
    hdf_input(hf,eos,name);
    hf.close();

    // Check causality
    for(size_t j=0;j<eos.get_nlines()-1;j++) {
      if (j==0) {
	t.set_output_level(2);
      } else {
	t.set_output_level(1);
      }
      t.test_gen(eos.get("ed",j+1)>eos.get("ed",j),"ed");
      t.test_gen(eos.get("pr",j+1)>eos.get("pr",j),"pr");
    }
      
    t.set_output_level(2);
    cout << endl;
  }

  {
    cout << "test_debug_star: " << endl;
    hf.open("data_temp/debug_star.o2");
    table_units<> mvsr;
    name="";
    hdf_input(hf,mvsr,name);
    hf.close();

    // Check maximum mass
    t.test_gen(mvsr.max("gm")>2.0,"M_max");
    cout << endl;
  }

  cout << "test_data: " << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_data_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check average radius for one of the stars
    t.test_rel(vector_mean(mcmc["Rns_1608"]),11.5,1.0,"R_1608");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
  }
  cout << endl;

  cout << "test_nodata: " << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_nodata_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check average radius
    t.test_rel(vector_mean(mcmc["R_43"]),11.5,1.0,"R_43");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_cthick:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_cthick_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check crust thickness
    t.test_rel(vector_mean(mcmc["CT_50"]),1.0,1.0,"CT_50");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
    
  }
  cout << endl;

  cout << "test_fixp:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/fixp_nodata_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check fixp EOS parameter
    t.test_rel(vector_mean(mcmc["pres1"]),1.0,1.0,"pres1");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
    
  }
  cout << endl;

  cout << "test_qt:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/qt_nodata_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check fixp EOS parameter
    t.test_rel(vector_mean(mcmc["index3"]),5.0,1.0,"index3");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
    
  }
  cout << endl;

  cout << "test_qf:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/qf_nodata_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check fixp EOS parameter
    t.test_rel(vector_mean(mcmc["pres1"]),1.0,1.0,"pres1");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
    
  }
  cout << endl;

  cout << "test_warmup:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_warmup_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    mcmc.delete_rows_func("mult<0.5");
    cout << mcmc.get_nlines() << " " << n_accept[0] << " "
	 << n_accept[1] << endl;
    //t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+2,"iters table");
    
  }
  cout << endl;

  cout << "test_ai:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_ai_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_addl:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_addl_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check moment of inertia
    t.test_rel(vector_mean(mcmc["I_43"]),70.0,10.0,"I_43");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_crustL:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_crustL_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check crust thickness parameter
    t.test_rel(vector_mean(mcmc["CT_50"]),1.0,1.0,"CT_50");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_storej:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_storej_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_tableseq:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_tableseq_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  cout << "test_rejtab:" << endl;
  for(size_t irank=0;irank<2;irank++) {

    hf.open(((string)"data_temp/twop_rejtab_")+o2scl::szttos(irank)+
	    "_out");
    table_units<> mcmc;
    vector<size_t> n_accept, n_reject;
    name="";
    hdf_input(hf,mcmc,name);
    hf.get_szt_vec("n_accept",n_accept);
    hf.get_szt_vec("n_reject",n_reject);
    hf.close();

    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==100,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==100,"iters thread 1");
    t.test_gen(mcmc.get_nlines()==n_accept[0]+n_accept[1]+1,"iters table");
    
  }
  cout << endl;

  return 0;
}
