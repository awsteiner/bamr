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
#include "nstar_cold2.h"

using namespace std;
using namespace o2scl;
using namespace bamr;

int nstar_cold2::calc_eos(double np_0) {

  int ret;
  
  if (eos_set==false) {
    O2SCL_ERR("EOS not set in nstar_cold2::calc_eos().",exc_efailed);
    return exc_efailed;
  }

  /* This code was in because saturation was used to calibrate the nb
     calculation. The nb calculation is now in the table, so there's
     probably no need for this exception here, but we leave it for now
     just in case.
  */
  if (nb_end<0.16) {
    O2SCL_ERR("nstar_cold2 doesn't support small nb_ends.",exc_efailed);
    return exc_efailed;
  }

  eost->clear();

  eost->line_of_names("ed pr nb");
  eost->set_unit("ed","1/fm^4");
  eost->set_unit("pr","1/fm^4");
  eost->set_unit("nb","1/fm^3");
  
  double x;
  // Initial guess for proton number density
  if (fabs(np_0)<1.0e-12) x=nb_start/3.0;
  else x=np_0;
  
  funct sf=std::bind(std::mem_fn<double(double)>
		       (&nstar_cold2::solve_fun),
		       this,std::placeholders::_1);

  bool success=true, n1_set=false;
  double y;

  for(barn=nb_start;barn<=nb_end+dnb/10.0;barn+=dnb) {

    ret=rp->solve(x,sf);
    y=solve_fun(x);
    
    if (ret!=0 || fabs(y)>solver_tol) {
      success=false;
    }

    if (include_muons) {
      h=hb+e+mu;
    } else {
      h=hb+e;
    }

    double line[3]={h.ed,h.pr,barn};
    eost->line_of_data(3,line);

  }

  if (success==false) {
    O2SCL_ERR("Solving for EOS failed in calc_eos().",exc_efailed);
  }

  return 0;
}

