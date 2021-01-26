/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2019, Andrew W. Steiner
  
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

int bamr_class::fill_emu(const ubvector &pars, double weight, 
         std::vector<double> &line, model_data &dat) {

  return o2scl::success;
}

int bamr_class::fill(const ubvector &pars, double weight, 
		     std::vector<double> &line, model_data &dat) {

  model &m=*this->mod;

  size_t n_params=pars.size();
  
  double nbmax2=0.0, emax=0.0, pmax=0.0, nbmax=0.0, mmax=0.0, rmax=0.0;

  if (m.has_eos) {

    // The central energy density in the maximum mass configuration
    cout << "bamr_class: check emax." << endl;
    emax=dat.mvsr.max("ed");
    // The central pressure in the maximum mass configuration
    pmax=dat.mvsr.max("pr");
    // The maximum mass
    mmax=dat.mvsr.get_constant("m_max");
    // The radius of the maximum mass star
    rmax=dat.mvsr.get_constant("r_max");

    if (set->baryon_density) {

      // The highest baryon density in the EOS table
      nbmax2=dat.eos.max("nb");
      // The central baryon density in the maximum mass configuration
      nbmax=dat.mvsr.get_constant("nb_max");
      
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
      double emax2=dat.eos.max("ed");
      if (eval<emax2) {
	double pres_temp=dat.eos.interp("ed",eval,"pr");
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
      line.push_back(dat.mvsr.interp("gm",mval,"r"));
      if (m.has_eos) {
	line.push_back(dat.mvsr.interp("gm",mval,"pr"));
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
	  double pres_temp=dat.eos.interp("nb",nbval,"pr");
	  if (pres_temp<pmax) {
	    line.push_back(pres_temp);
	  } else {
	    line.push_back(0.0);
	  }
	  double eval2=dat.eos.interp("nb",nbval,"ed");
	  double eoa_val2=eval2/nbval-939.0/o2scl_const::hc_mev_fm;
	  line.push_back(eoa_val2);
	} else {
	  line.push_back(0.0);
	  line.push_back(0.0);
	}
      }
    }
    if (m.has_esym) {
      line.push_back(dat.eos.get_constant("S"));
      line.push_back(dat.eos.get_constant("L"));
    }
    line.push_back(rmax);
    line.push_back(mmax);
    line.push_back(pmax);
    line.push_back(emax);
    if (set->baryon_density) line.push_back(nbmax);
    for(size_t i=0;i<nsd->n_sources;i++) {
      double val=dat.mvsr.interp
	("gm",pars[n_params-nsd->n_sources+i],"ed");
      line.push_back(val);
    }
    if (set->baryon_density) {
      for(size_t i=0;i<nsd->n_sources;i++) {
	double val2=dat.mvsr.interp
	  ("gm",pars[n_params-nsd->n_sources+i],"nb");
	line.push_back(val2);
      }
    }
  }
  
  if (set->baryon_density) {
    cout << "bamr_class:: check gms." << endl;
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
      double mval=m.m_grid[i];
      if (mval<mmax) {
	double rval=dat.mvsr.interp("gm",mval,"r");
	line.push_back(rval-dat.mvsr.interp("gm",mval,"r0"));
      } else {
	line.push_back(0.0);
      }
    }
  }
  if (set->addl_quants) {
    for(int i=0;i<set->grid_size;i++) {
      double mval=m.m_grid[i];

      if (mval<mmax) {
	
	// Baryonic mass
	double bm=dat.mvsr.interp("gm",mval,"bm");
	line.push_back(bm);

	// Binding energy
	line.push_back(bm-mval);

	// Moment of inertia
	double rad=dat.mvsr.interp("gm",mval,"r");
	// rjw is km^4, so dividing by km/Msun gives Msun*km^2
	double I=dat.mvsr.interp("gm",mval,"rjw")/3.0/schwarz_km;
	line.push_back(I);

	// To compute I_bar, divide by G^2*M^3
	double I_bar=I*4.0/schwarz_km/schwarz_km/mval/mval/mval;
	line.push_back(I_bar);

	/*
	// Relation for lambda given I from Kent Yagi based
	// on Yagi and Yunes, Science, (2014)
	double a0=-210.327;
	double a1=481.472;
	double a2=-464.271;
	double a3=241.564;
	double a4=-70.6449;
	double a5=10.9775;
	double a6=-0.707066;
	*/
	
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

	line.push_back(Lambda_bar);
	
      } else {
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
	line.push_back(0.0);
      }
    }
  }
  if (nsd->source_fnames_alt.size()>0) {
    for(size_t i=0;i<nsd->n_sources;i++) {
      // Compute alternate probability from an insignificant bit
      // in the mass 
      double alt=dat.mass[i]*1.0e8-((double)((int)(dat.mass[i]*1.0e8)));
      if (alt<2.0/3.0) {
	line.push_back(0.0);
      } else {
	line.push_back(1.0);
      }
    }
  }

  return o2scl::success;
}

int bamr_class::compute_point(const ubvector &pars, std::ofstream &scr_out, 
			      double &weight, model_data &dat) {

  // Compute the M vs R curve and return if it failed
  int iret;
  mod->compute_star(pars,scr_out,iret,dat);
  if (iret!=0) {
    weight=0.0;
    return iret;
  }
  
  return mod->compute_point(pars,scr_out,weight,dat);
}

int bamr_class::compute_point_emu(const ubvector &pars, std::ofstream &scr_out, 
            double &weight, model_data &dat) {

  return mod->compute_point_emu(pars,scr_out,weight,dat);
}

void create_pointers(void *&bcp2,
		     void *&mdp2) {

  bamr::bamr_class *bcp=new bamr::bamr_class;
  bcp2=(void *)bcp;
  cout << "Creating settings object." << endl;
  bcp->set=std::make_shared<settings>();
  cout << "Creating ns_data object." << endl;
  bcp->nsd=std::make_shared<ns_data>();
  cout << "Creating model." << endl;
  std::shared_ptr<model> mnew(new two_polytropes(bcp->set,bcp->nsd));
  cout << "Setting model." << endl;
  bcp->mod=mnew;
  bcp->model_type="twop";
  cout << "Creating model_data object." << endl;
  bamr::model_data *mdp=new bamr::model_data;
  mdp2=(void *)mdp;

  cout << "Calling get_param_info()." << endl;
  std::vector<std::string> names;
  std::vector<std::string> units;

  ubvector low;
  ubvector high;
  bcp->mod->get_param_info(names,units,low,high);

  cout << "Getting initial point." << endl;
  ubvector init(names.size());
  bcp->mod->initial_point(init);

  cout.setf(ios::scientific);
  o2scl::vector_out(cout,init,true);

  ofstream fout;//("bamr.py.out");
  
  cout << "Calling compute_point()." << endl;
  double weight;
  int ret=bcp->compute_point(init,fout,weight,*mdp);
  cout << "ret,weight: " << ret << " " << weight << endl;

  if (ret==0) {
    table_units<> &mvsr=mdp->mvsr;
    table_units<> &eos=mdp->eos;
    cout << "M_max: " << mvsr.max("gm") << endl;
  }

  //fout.close();
  
  return;
}

void py_compute_point(void *bcp2, void *mdp2,
		      int nv, double *vals) {
  
  bamr::bamr_class *bcp=(bamr::bamr_class *)bcp2;
  bamr::model_data *mdp=(bamr::model_data *)mdp2;
  
  ofstream fout;
  ubvector point(nv);
  for(int j=0;j<nv;j++) {
    point[j]=vals[j];
    cout << j << " " << point[j] << endl;
  }
  double weight;
  int ret=bcp->compute_point(point,fout,weight,*mdp);
  cout << "ret,weight: " << ret << " " << weight << endl;
  
  return;
}

int get_mvsr_column(void *mdp2, char *col_name, int &n, double *&ptr) {
  
  bamr::model_data *mdp=(bamr::model_data *)mdp2;
  n=mdp->mvsr.get_nlines();
  std::string stmp=col_name;
  if (mdp->mvsr.is_column(stmp)==false) {
    return 2;
  }
  const std::vector<double> &col=mdp->mvsr.get_column(stmp);
  ptr=(double *)&col[0];
  return 0;
}

void destroy_pointers(void *bcp2, void *mdp2) {
  bamr::bamr_class *bcp=(bamr::bamr_class *)bcp2;
  bamr::model_data *mdp=(bamr::model_data *)mdp2;
  delete bcp;
  delete mdp;
  return;
}
