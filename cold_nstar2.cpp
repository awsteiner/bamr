/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2014, Andrew W. Steiner
  
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
#include "cold_nstar2.h"

using namespace std;
using namespace o2scl;
using namespace bamr;

int cold_nstar2::calc_eos(double &n1, double &e1, double np_0) {

  int ret;
  
  if (eos_set==false) {
    O2SCL_ERR("EOS not set in cold_nstar2::calc_eos().",exc_efailed);
    return exc_efailed;
  }

  if (nb_end<0.16) {
    O2SCL_ERR("cold_nstar2 doesn't support small nb_ends.",exc_efailed);
    return exc_efailed;
  }
  
  eost->clear_table();
  eost->line_of_names("ed pr");
  eost->set_unit("ed","1/fm^4");
  eost->set_unit("pr","1/fm^4");
  
  double x;
  // Initial guess for proton number density
  if (fabs(np_0)<1.0e-12) x=nb_start/3.0;
  else x=np_0;
  
  funct_mfptr<cold_nstar2> sf(this,&cold_nstar2::solve_fun);
  
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
    
    double line[2]={h.ed,h.pr};
    eost->line_of_data(2,line);
    
    // Set calibration for baryon density near the
    // nuclear saturation density
    if (np->n+pp->n>=0.01599 && n1_set==false) {
      n1=np->n+pp->n;
      e1=h.ed;
      n1_set=true;
    }

  }

  if (success==false) {
    O2SCL_ERR_RET("Solving for EOS failed in calc_eos().",exc_efailed);
  }

  return 0;
}

