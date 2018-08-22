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
  t.set_output_level(1);

  hdf_file hf;
  std::string name;
  
  {
    hf.open("data_temp/debug_eos.o2");
    table_units<> eos;
    name="eos";
    hdf_input(hf,eos,name);
    hf.close();

    // Check causality
    for(size_t j=0;j<eos.get_nlines()-1;j++) {
      t.test_gen(eos.get("ed",j+1)>eos.get("ed",j),"ed");
      t.test_gen(eos.get("pr",j+1)>eos.get("pr",j),"pr");
    }
      
  }

  {
    hf.open("data_temp/debug_star.o2");
    table_units<> mvsr;
    name="";
    hdf_input(hf,mvsr,name);
    hf.close();

    // Check maximum mass
    t.test_gen(mvsr.max("gm")<2.0,"M_max");
  }

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
    t.test_rel(vector_mean(mcmc["R_1608"]),11.5,1.0,"R_1608");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    
  }

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
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    
  }

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
    t.test_rel(vector_mean(mcmc["ct_50"]),1.0,1.0,"ct_50");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    
  }

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
    t.test_rel(vector_mean(mcmc["X"]),1.0,1.0,"X");
    // Check total iterations
    t.test_gen(n_accept[0]+n_reject[0]==300,"iters thread 0");
    t.test_gen(n_accept[1]+n_reject[1]==300,"iters thread 1");
    
  }

  return 0;
}
