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
#include "models.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;
using namespace bamr;

two_polytropes::two_polytropes() {
  se.kpp=0.0;
  se.n0=0.16;
  se.eoa=-16.0/hc_mev_fm;
  se.a=17.0/hc_mev_fm;
    
  neut.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  prot.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_proton),2.0);

  cns.set_n_and_p(neut,prot);
  // We include muons by default, but they rarely appear at low-density
  cns.include_muons=true;
}

void two_polytropes::low_limits(entry &e) {
    
  e.params[0]=180.0/hc_mev_fm;
  e.params[1]=-1000.0/hc_mev_fm;
  e.params[2]=28.0/hc_mev_fm;
  e.params[3]=0.0;
  // The value 0.75 fm^{-4} is about the energy density of nuclear
  // matter
  e.params[4]=0.75;
  e.params[5]=0.2;
  e.params[6]=0.75;
  e.params[7]=0.2;

  return;
}

void two_polytropes::high_limits(entry &e) {
    
  e.params[0]=300.0/hc_mev_fm;
  // FSU gold is -280 MeV or so
  e.params[1]=-200.0/hc_mev_fm;
  e.params[2]=38.0/hc_mev_fm;
  e.params[3]=1.2;
  // The value of high.trans1 has to be small enough because we
  // don't want to evaluate the schematic EOS to too high of a
  // density.
  e.params[4]=3.0;
  e.params[5]=1.5;
  e.params[6]=8.0;
  e.params[7]=2.0;

  return;
}

string two_polytropes::param_name(size_t i) {
  if (i==0) return "comp";
  else if (i==1) return "kprime";
  else if (i==2) return "esym";
  else if (i==3) return "gamma";
  else if (i==4) return "trans1";
  else if (i==5) return "index1";
  else if (i==6) return "trans2";
  return "index2";
}

string two_polytropes::param_unit(size_t i) {
  if (i==0) return "1/fm";
  else if (i==1) return "1/fm";
  else if (i==2) return "1/fm";
  else if (i==3) return ".";
  else if (i==4) return "1/fm^4";
  else if (i==5) return ".";
  else if (i==6) return "1/fm^4";
  return ".";
}

void two_polytropes::first_point(entry &e) {
  e.params[0]=1.0;
  e.params[1]=-3.0;
  e.params[2]=0.165;
  e.params[3]=0.644;
  e.params[4]=1.51;
  e.params[5]=0.576;
  e.params[6]=4.60;
  e.params[7]=1.21;
  return;
}

void two_polytropes::compute_eos(entry &e, bool &fail, ofstream &scr_out) {

  fail=false;
  if (e.params[4]>e.params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    fail=true;
    return;
  }

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos(nb_n1,nb_e1);
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // What does this do? 1/31/11 - It appears not to 
  // get called frequently
  size_t nl=tab_eos->get_nlines();
  for(size_t i=0;nl>0 && i<nl-1;i++) {
    if ((*tab_eos)["ed"][i]>(*tab_eos)["ed"][i+1]) {
      tab_eos->set("ed",i+1,(*tab_eos)["ed"][i]*2.0);
      scr_out << "Pressure grid fix." << endl;
    }
  }

  // Determine 1st polytrope coefficient
  double pr1=tab_eos->interp("ed",e.params[4],"pr");
  double coeff1=pr1/pow(e.params[4],1.0+1.0/e.params[5]);

  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    fail=true;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<e.params[4]) {
    scr_out << "Gap between low-density EOS and polytrope " << endl;
    exit(-1);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>e.params[4]) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    fail=true;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=e.params[4];ed<e.params[6]-0.001;ed+=0.05) {
    double line[2]={ed,coeff1*pow(ed,1.0+1.0/e.params[5])};
    tab_eos->line_of_data(2,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    fail=true;
    return;
  }

  // Determine 2nd polytrope coefficient
  double pr2=tab_eos->interp("ed",e.params[6],"pr");
  double coeff2=pr2/pow(e.params[6],1.0+1.0/e.params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    fail=true;
    return;
  }

  // Add second polytrope to table
  for(double ed=e.params[6];ed<=10.0;ed+=0.05) {
    double line[2]={ed,coeff2*pow(ed,1.0+1.0/e.params[7])};
    tab_eos->line_of_data(2,line);
  }

  return;
}

void alt_polytropes::low_limits(entry &e) {
  two_polytropes::low_limits(e);
  e.params[5]=0.0;
  e.params[7]=0.0;
  return;
}

void alt_polytropes::high_limits(entry &e) {
  two_polytropes::high_limits(e);
  e.params[5]=6.0;
  e.params[7]=3.0;
  return;
}

string alt_polytropes::param_name(size_t i) {
  if (i==5) return "exp1";
  else if (i==7) return "exp2";
  return two_polytropes::param_name(i);
}

string alt_polytropes::param_unit(size_t i) {
  if (i==5) return ".";
  else if (i==7) return ".";
  return two_polytropes::param_unit(i);
}

void alt_polytropes::first_point(entry &e) {
  e.params[0]=1.0;
  e.params[1]=-2.66;
  e.params[2]=0.165;
  e.params[3]=0.66;
  e.params[4]=1.48;
  e.params[5]=2.913;
  e.params[6]=4.066;
  e.params[7]=1.80;
  return;
}

void alt_polytropes::compute_eos(entry &e, bool &fail, ofstream &scr_out) {
  
  fail=false;
  if (e.params[4]>e.params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    fail=true;
    return;
  }

  schematic_eos &se=this->se;
  cold_nstar2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos(nb_n1,nb_e1);
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // What does this do? 1/31/11 - It appears not to 
  // get called frequently
  size_t nl=tab_eos->get_nlines();
  for(size_t i=0;nl>0 && i<nl-1;i++) {
    if ((*tab_eos)["ed"][i]>(*tab_eos)["ed"][i+1]) {
      tab_eos->set("ed",i+1,(*tab_eos)["ed"][i]*2.0);
      scr_out << "Pressure grid fix." << endl;
    }
  }

  // Determine 1st polytrope coefficient
  double pr1=tab_eos->interp("ed",e.params[4],"pr");
  double coeff1=pr1/pow(e.params[4],e.params[5]);
    
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    fail=true;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<e.params[4]) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    exit(-1);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>e.params[4]) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    fail=true;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=e.params[4];ed<e.params[6]-0.001;ed+=0.05) {
    double line[2]={ed,coeff1*pow(ed,e.params[5])};
    tab_eos->line_of_data(2,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    fail=true;
    return;
  }

  // Determine 2nd polytrope coefficient
  double pr2=tab_eos->interp("ed",e.params[6],"pr");
  double coeff2=pr2/pow(e.params[6],e.params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    fail=true;
    return;
  }

  // Add second polytrope to table
  for(double ed=e.params[6];ed<=10.0;ed+=0.05) {
    double line[2]={ed,coeff2*pow(ed,e.params[7])};
    tab_eos->line_of_data(2,line);
  }

  return;
}
  
void fixed_pressure::low_limits(entry &e) {
  two_polytropes::low_limits(e);
  e.params[4]=0.0;
  e.params[5]=0.0;
  e.params[6]=0.0;
  e.params[7]=0.0;
  return;
}

void fixed_pressure::high_limits(entry &e) {
  two_polytropes::high_limits(e);
  e.params[4]=0.3;
  e.params[5]=1.5;
  e.params[6]=2.5;
  e.params[7]=2.5;
  return;
}

string fixed_pressure::param_name(size_t i) {
  if (i==4) return "pres1";
  else if (i==5) return "pres2";
  else if (i==6) return "pres3";
  else if (i==7) return "pres4";
  return two_polytropes::param_name(i);
}

string fixed_pressure::param_unit(size_t i) {
  if (i==4) return "1/fm^4";
  else if (i==5) return "1/fm^4";
  else if (i==6) return "1/fm^4";
  else if (i==7) return "1/fm^4";
  return two_polytropes::param_unit(i);
}

void fixed_pressure::first_point(entry &e) {
  e.params[0]=1.0;
  e.params[1]=-2.5;
  e.params[2]=0.165;
  e.params[3]=0.8;
  e.params[4]=0.024;
  e.params[5]=0.74;
  e.params[6]=0.60;
  e.params[7]=1.84;
  return;
}

void fixed_pressure::compute_eos(entry &e, bool &fail, ofstream &scr_out) {

  fail=false;

  schematic_eos &se=this->se;
  cold_nstar2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos(nb_n1,nb_e1);
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // Compute boundary energy density and pressure
  double ed_last=1.0;
  double pr_last=tab_eos->interp("ed",1.0,"pr");
    
  // Remove extra rows from EOS near saturation
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>1.0) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Compute pressures on grid, ed=2, 3, 5, 7 fm^{-4}
  double p2=pr_last+e.params[4];
  double p3=p2+e.params[5];
  double p5=p3+e.params[6];

  // Add 1st high-density EOS
  for(double ed=1.0;ed<2.0-1.0e-4;ed+=0.1) {
    double line[2]={ed,pr_last+e.params[4]*(ed-1.0)/(2.0-1.0)};
    tab_eos->line_of_data(2,line);
  }

  // Add 2nd high-density EOS
  for(double ed=2.0;ed<3.0-1.0e-4;ed+=0.1) {
    double line[2]={ed,p2+e.params[5]*(ed-2.0)/(3.0-2.0)};
    tab_eos->line_of_data(2,line);
  }

  // Add 3rd high-density EOS
  for(double ed=3.0;ed<5.0-1.0e-4;ed+=0.1) {
    double line[2]={ed,p3+e.params[6]*(ed-3.0)/(5.0-3.0)};
    tab_eos->line_of_data(2,line);
  }

  // Add 4th high-density EOS
  for(double ed=5.0;ed<10.0-1.0e-4;ed+=0.2) {
    double line[2]={ed,p5+e.params[7]*(ed-5.0)/(7.0-5.0)};
    tab_eos->line_of_data(2,line);
  }

  return;
}
  
void generic_quarks::low_limits(entry &e) {
    
  e.params[0]=180.0/hc_mev_fm;
  e.params[1]=-1000.0/hc_mev_fm;
  e.params[2]=28.0/hc_mev_fm;
  //e.params[3]=0.2;
  e.params[3]=0.0;
  // 0.75 is about the energy density of nuclear matter
  e.params[4]=0.75;
  e.params[5]=2.0;
  e.params[6]=0.75;
  // a2
  e.params[7]=-0.3;
  // a4
  e.params[8]=0.045;

  return;
}

void generic_quarks::high_limits(entry &e) {
    
  e.params[0]=300.0/hc_mev_fm;
  // FSU gold is -280 or so
  e.params[1]=-200.0/hc_mev_fm;
  e.params[2]=38.0/hc_mev_fm;
  e.params[3]=1.2;
  // The value of high.trans1 has to be small enough because we
  // don't want to evaluate the schematic EOS to too high of a
  // density.
  e.params[4]=3.0;
  e.params[5]=12.0;
  e.params[6]=4.5;
  // a2
  e.params[7]=0.3;
  // a4
  e.params[8]=0.08;

  return;
}

string generic_quarks::param_name(size_t i) {
  if (i==0) return "comp";
  else if (i==1) return "kprime";
  else if (i==2) return "esym";
  else if (i==3) return "gamma";
  else if (i==4) return "trans1";
  else if (i==5) return "exp1";
  else if (i==6) return "trans2";
  else if (i==7) return "a2";
  return "a4";
}

string generic_quarks::param_unit(size_t i) {
  if (i==0) return "1/fm";
  else if (i==1) return "1/fm";
  else if (i==2) return "1/fm";
  else if (i==3) return ".";
  else if (i==4) return "1/fm^4";
  else if (i==5) return ".";
  else if (i==6) return "1/fm^4";
  else if (i==7) return "1/fm^2";
  return ".";
}

void generic_quarks::first_point(entry &e) {
  e.params[0]=1.19;
  e.params[1]=-2.52;
  e.params[2]=0.188;
  e.params[3]=0.357;
  e.params[4]=1.86;
  e.params[5]=5.70;
  e.params[6]=2.29;
  e.params[7]=0.1907;
  e.params[8]=0.0796;
  return;
}

void generic_quarks::compute_eos(entry &e, bool &fail, ofstream &scr_out) {

  fail=false;

  schematic_eos &se=this->se;
  cold_nstar2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos(nb_n1,nb_e1);
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // Double check that the table is non-empty (we have to do this
  // especially for the size_t index in the for loop below)
  if (tab_eos->get_nlines()==0) {
    O2SCL_ERR("Table empty in generic quarks.",exc_efailed);
  }

  // Determine 1st polytrope coefficient
  for(size_t i=0;i<tab_eos->get_nlines()-1;i++) {
    if ((*tab_eos)["ed"][i]>(*tab_eos)["ed"][i+1]) {
      tab_eos->set("ed",i+1,(*tab_eos)["ed"][i]*2.0);
      scr_out << "Pressure grid fix." << flush;
    }
  }
  double pr1=tab_eos->interp("ed",e.params[4],"pr");
  double coeff1=pr1/pow(e.params[4],e.params[5]);
  
  tab_eos->add_constant("pr_pt",pr1);
  
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    fail=true;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<e.params[4]) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    exit(-1);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>e.params[4]) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    fail=true;
    return;
  }

  // Add first polytrope to table
  for(double ed=e.params[4];ed<=e.params[6];ed+=0.05) {
    double line[2]={ed,coeff1*pow(ed,e.params[5])};
    tab_eos->line_of_data(2,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    fail=true;
    return;
  }

  double a2=e.params[7];
  double a4=e.params[8];
  double ed_trans=e.params[6];
  double pr_trans=coeff1*pow(ed_trans,e.params[5]);

  tab_eos->add_constant("pr_q",pr_trans);

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
      fail=true;
      return;
    }
  }
    
  // Compute maximum chemical potential and bag constant
  double mu_end=10.0;
  if (mu_end<2.0*mu_start) mu_end=2.0*mu_start;
  double musq=mu_start*mu_start;
  double bag=ed_trans-a2*musq-3.0*a4*musq*musq;

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
      fail=true;
      return;
    }
      
    double line[2]={bag+a2*musq+3.0*a4*musq*musq,
		    -bag+a2*musq+a4*musq*musq};
    tab_eos->line_of_data(2,line);
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

void quark_star::low_limits(entry &e) {

  // B
  e.params[0]=-10.0;
  // c
  e.params[1]=0.0;
  // Delta
  e.params[2]=0.0;
  // ms
  e.params[3]=0.75;
    
  return;
}
  
void quark_star::high_limits(entry &e) {

  // B
  e.params[0]=10.0;
  // c
  e.params[1]=0.4;
  // Delta
  e.params[2]=1.0;
  // ms
  e.params[3]=2.5;

  return;
}

std::string quark_star::param_name(size_t i) {
  if (i==0) return "B";
  else if (i==1) return "c";
  else if (i==2) return "Delta";
  return "ms";
}

std::string quark_star::param_unit(size_t i) {
  if (i==0) return "1/fm";
  else if (i==1) return ".";
  else if (i==2) return "1/fm";
  return "1/fm";
}
  
void quark_star::first_point(entry &e) {
  e.params[0]=0.2446;
  e.params[1]=0.0740;
  e.params[2]=0.00289;
  e.params[3]=0.0754;
  return;
}

void quark_star::compute_eos(entry &e, bool &fail, std::ofstream &scr_out) {
  
  fail=false;

  B=e.params[0];
  c=e.params[1];
  Delta=e.params[2];
  ms=e.params[3];

  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
    
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
    fail=true;
    scr_out << "No zero pressure solution." << std::endl;
    return;
  }

  // Then call the root finder
  x[0]=mu_0;
  mm_funct_mfptr<quark_star> fmf(this,&quark_star::pressure);
  gmh.err_nonconv=false;
  int ret=gmh.msolve(1,x,fmf);
  if (ret!=0) {
    scr_out << "Solver failed in qstar." << std::endl;
    fail=true;
    return;
  }
  mu_0=x[0];
  //grb.solve_bkt(mu_1,mu_0,fmf);
    
  tab_eos->clear_table();
  tab_eos->new_column("ed");
  tab_eos->new_column("pr");
  tab_eos->new_column("nb");
  tab_eos->set_unit("ed","1/fm^4");
  tab_eos->set_unit("pr","1/fm^4");
  tab_eos->set_unit("nb","1/fm^3");

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
      fail=true;
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
	fail=true;
	return;
      }

      // Fix the first point at exactly zero pressure
      double line[3]={ed,0.0,nb};
      tab_eos->line_of_data(3,line);
	
    } else {
	
      double line[3]={ed,pr,nb};
      tab_eos->line_of_data(3,line);

    }

    ed_last=ed;
  }

  return;
}
