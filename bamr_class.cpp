/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
// For I/O with HDF files
using namespace o2scl_hdf;
// For pi, pi^2, etc.
using namespace o2scl_const;
using namespace bamr;

int bamr_class::fill(const ubvector &pars, double weight, 
		     std::vector<double> &line, model_data &dat) {
  
  model &m=*this->mod;

  size_t n_params=pars.size();
  
  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (m.has_eos) {

    // The central energy density in the maximum mass configuration
    emax=dat.mvsr->max("ed");
    // The central pressure in the maximum mass configuration
    pmax=dat.mvsr->max("pr");
    // The maximum mass
    mmax=dat.mvsr->get_constant("m_max");
    // The radius of the maximum mass star
    rmax=dat.mvsr->get_constant("r_max");

    if (set->baryon_density) {

      // The highest baryon density in the EOS table
      nbmax2=dat.eos->max("nb");
      // The central baryon density in the maximum mass configuration
      nbmax=dat.mvsr->get_constant("nb_max");
      
    }

  } else {
    // Need to set mmax for no EOS models to figure out how 
    // high up we should go for the radius grid 
    mmax=3.0;
  }

  for(size_t i=0;i<nsd->n_sources;i++) {
    line.push_back(dat.wgts[i]);
  }
  for(size_t i=0;i<nsd->n_sources;i++) {
    line.push_back(dat.rad[i]);
  }
  for(size_t i=0;i<nsd->n_sources;i++) {
    line.push_back(dat.mass[i]);
  }
  
  if (m.has_eos) {
    for(int i=0;i<set->grid_size;i++) {
      double eval=m.e_grid[i];
      // Make sure the energy density from the grid
      // isn't beyond the table limit
      double emax2=dat.eos->max("ed");
      if (eval<emax2) {
	double pres_temp=dat.eos->interp("ed",eval,"pr");
	//if (pres_temp<pmax) {
	line.push_back(pres_temp);
	//} else {
	//line.push_back(0.0);
	//}
      } else {
	line.push_back(0.0);
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
      line.push_back(dat.mvsr->interp("gm",mval,"r"));
      if (m.has_eos) {
	line.push_back(dat.mvsr->interp("gm",mval,"pr"));
      }
    } else {
      line.push_back(0.0);
      if (m.has_eos) {
	line.push_back(0.0);
      }
    }
  }
  if (m.has_eos) {
    if (set->baryon_density) {
      for(int i=0;i<set->grid_size;i++) {
	double nbval=m.nb_grid[i];
	if (nbval<nbmax2) {
	  double pres_temp=dat.eos->interp("nb",nbval,"pr");
	  if (pres_temp<pmax) {
	    line.push_back(pres_temp);
	  } else {
	    line.push_back(0.0);
	  }
	  double eval2=dat.eos->interp("nb",nbval,"ed");
	  double eoa_val2=eval2/nbval-939.0/o2scl_const::hc_mev_fm;
	  line.push_back(eoa_val2);
	} else {
	  line.push_back(0.0);
	  line.push_back(0.0);
	}
      }
    }
    if (m.has_esym) {
      line.push_back(dat.eos->get_constant("S"));
      line.push_back(dat.eos->get_constant("L"));
    }
    line.push_back(rmax);
    line.push_back(mmax);
    line.push_back(pmax);
    line.push_back(emax);
    if (set->baryon_density) line.push_back(nbmax);
    for(size_t i=0;i<nsd->n_sources;i++) {
      double val=dat.mvsr->interp
	("gm",pars[n_params-nsd->n_sources+i],"ed");
      line.push_back(val);
    }
    if (set->baryon_density) {
      for(size_t i=0;i<nsd->n_sources;i++) {
	double val2=dat.mvsr->interp
	  ("gm",pars[n_params-nsd->n_sources+i],"nb");
	line.push_back(val2);
      }
    }
  }
  
  if (set->baryon_density) {
    line.push_back(dat.mvsr->get_constant("gm_nb1"));
    line.push_back(dat.mvsr->get_constant("r_nb1"));
    line.push_back(dat.mvsr->get_constant("gm_nb2"));
    line.push_back(dat.mvsr->get_constant("r_nb2"));
    line.push_back(dat.mvsr->get_constant("gm_nb3"));
    line.push_back(dat.mvsr->get_constant("r_nb3"));
    line.push_back(dat.mvsr->get_constant("gm_nb4"));
    line.push_back(dat.mvsr->get_constant("r_nb4"));
    line.push_back(dat.mvsr->get_constant("gm_nb5"));
    line.push_back(dat.mvsr->get_constant("r_nb5"));
  }
  if (set->compute_cthick) {
    line.push_back(dat.eos->get_constant("nt"));
    line.push_back(dat.eos->get_constant("prt"));
    for(int i=0;i<set->grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax) {
	double rval=dat.mvsr->interp("gm",mval,"r");
	line.push_back(rval-dat.mvsr->interp("gm",mval,"r0"));
      } else {
	line.push_back(0.0);
      }
    }
  }
  if (set->addl_quants) {
    double mmax=dat.mvsr->get_constant("m_max");
    for(int i=0;i<set->grid_size;i++) {
      double mval=m.m_grid[i];
      if (mval<mmax) {
	double bm=dat.mvsr->interp("gm",mval,"bm");
	double rad=dat.mvsr->interp("gm",mval,"r");
	// rjw is km^4, so dividing by km/Msun gives Msun*km^2
	double I=dat.mvsr->interp("gm",mval,"rjw")/3.0/schwarz_km;
	line.push_back(bm);
	line.push_back((bm-mval)/mval);
	// Make unitless by dividing by G^2
	double I_bar=I/mval/mval/mval/schwarz_km/schwarz_km*4.0;
	line.push_back(I_bar);

#ifdef NEVER_DEFINED
	// Relation for lambda given I from Kent Yagi based
	// on Yagi and Yunes, Science, (2014)
	double a0=-210.327;
	double a1=481.472;
	double a2=-464.271;
	double a3=241.564;
	double a4=-70.6449;
	double a5=10.9775;
	double a6=-0.707066;

	double li=log(I_bar);
	double li2=li*li;
	double li3=li*li2;
	double li4=li*li3;
	double li5=li*li4;
	double li6=li*li5;
	  
	double lambda_bar=exp(a0+a1*li+a2*li2+a3*li3+a4*li4+a5*li5+a6*li6);
	double lambda=lambda_bar*pow(mval,5.0);
	cout << I_bar << " " << lambda_bar << endl;
	exit(-1);
#endif
	
      } else {
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
      }
    }
  }

  return o2scl::success;
}

int bamr_class::compute_point(const ubvector &pars, std::ofstream &scr_out, 
			      double &weight, model_data &dat) {
  return mod->compute_point(pars,scr_out,weight,dat);
}

