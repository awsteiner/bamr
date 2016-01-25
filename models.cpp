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
    O2SCL_ERR("Failed to erase parameter.",o2scl::exc_esanity);
  }
  return;
}

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

void two_polytropes::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;
  if (e.params[4]>e.params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    success=bamr_class::ix_param_mismatch;
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
  cns.calc_eos();
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));
  
  // Transition densities
  double ed1=e.params[4];
  double ed2=e.params[6];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=tab_eos->interp("ed",ed1,"nb");
  double pr1=tab_eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,1.0+1.0/e.params[5]);

  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    success=bamr_class::ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>ed1) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=ed1;ed<ed2-0.001;ed+=0.05) {
    double pr=coeff1*pow(ed,1.0+1.0/e.params[5]);
    double nb=nb1*pow(ed/ed1,1.0+e.params[5])/
      pow((ed+pr)/(ed1+pr1),e.params[5]);
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }
  
  // Boundary baryon density and pressure
  double pr2=coeff1*pow(ed2,1.0+1.0/e.params[5]);
  double nb2=nb1*pow(ed2/ed1,1.0+e.params[5])/
    pow((ed2+pr2)/(ed1+pr1),e.params[5]);

  // Determine 2nd polytrope coefficient
  double coeff2=pr2/pow(ed2,1.0+1.0/e.params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    success=bamr_class::ix_neg_pressure;
    return;
  }

  // Add second polytrope to table
  for(double ed=ed2;ed<=10.0;ed+=0.05) {
    double pr=coeff2*pow(ed,1.0+1.0/e.params[7]);
    double nb=nb2*pow(ed/ed2,1.0+e.params[7])/
      pow((ed+pr)/(ed2+pr2),e.params[7]);
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }
  
  return;
}

void alt_polytropes::low_limits(entry &e) {
  two_polytropes::low_limits(e);
  e.params[5]=1.5;
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

void alt_polytropes::compute_eos(entry &e, int &success, ofstream &scr_out) {
  
  success=bamr_class::ix_success;
  if (e.params[4]>e.params[6]) {
    scr_out << "Rejected: Transition densities misordered." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos();
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // Transition densities
  double ed1=e.params[4];
  double ed2=e.params[6];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=tab_eos->interp("ed",ed1,"nb");
  double pr1=tab_eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,e.params[5]);
    
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    success=bamr_class::ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>ed1) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }

  // Add first polytrope to table. The shift of 0.001 is
  // important to ensure that we don't have two points
  // at the same energy density
  for(double ed=ed1;ed<ed2-0.001;ed+=0.05) {
    double pr=coeff1*pow(ed,e.params[5]);
    double nb=nb1*pow(ed/ed1,e.params[5]/(e.params[5]-1.0))*
      pow((ed+pr)/(ed1+pr1),1.0/(1.0-e.params[5]));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }

  // Boundary baryon density and pressure
  double pr2=coeff1*pow(ed2,e.params[5]);
  double nb2=nb1*pow(ed2/ed1,e.params[5]/(e.params[5]-1.0))*
    pow((ed2+pr2)/(ed1+pr1),1.0/(1.0-e.params[5]));
  
  // Determine 2nd polytrope coefficient
  double coeff2=pr2/pow(ed2,e.params[7]);

  if (coeff2<0.0 || pr2<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (2)." << endl;
    success=bamr_class::ix_neg_pressure;
    return;
  }

  // Add second polytrope to table
  for(double ed=ed2;ed<=10.0;ed+=0.05) {
    double pr=coeff2*pow(ed,e.params[7]);
    double nb=nb2*pow(ed/ed2,e.params[7]/(e.params[7]-1.0))*
      pow((ed+pr)/(ed2+pr2),1.0/(1.0-e.params[7]));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
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

void fixed_pressure::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos();
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));
  
  // Compute boundary energy density, baryon density and pressure
  double ed_last=1.0;
  double nb_last=tab_eos->interp("ed",1.0,"nb");
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

  // Computes slopes (squared sound speeds)
  double cs2_1=e.params[4]/(2.0-1.0);
  double cs2_2=e.params[5]/(3.0-2.0);
  double cs2_3=e.params[6]/(5.0-3.0);
  double cs2_4=e.params[7]/(7.0-5.0);

  // Add 1st high-density EOS
  for(double ed=1.0;ed<2.0-1.0e-4;ed+=0.1) {
    double pr=pr_last+e.params[4]*(ed-1.0)/(2.0-1.0);
    double nb=nb_last*pow((ed+pr_last+cs2_1*(ed-ed_last))/
			  (ed_last+pr_last),1.0/(1.0+cs2_1));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }
  double nb2=nb_last*pow((2.0+pr_last+cs2_1*(2.0-ed_last))/
			 (ed_last+pr_last),1.0/(1.0+cs2_1));
  
  // Add 2nd high-density EOS
  for(double ed=2.0;ed<3.0-1.0e-4;ed+=0.1) {
    double pr=p2+e.params[5]*(ed-2.0)/(3.0-2.0);
    double nb=nb2*pow((ed+p2+cs2_2*(ed-2.0))/(2.0+p2),1.0/(1.0+cs2_2));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }
  double nb3=nb2*pow((3.0+p2+cs2_2*(3.0-2.0))/(2.0+p2),1.0/(1.0+cs2_2));

  // Add 3rd high-density EOS
  for(double ed=3.0;ed<5.0-1.0e-4;ed+=0.1) {
    double pr=p3+e.params[6]*(ed-3.0)/(5.0-3.0);
    double nb=nb3*pow((ed+p3+cs2_3*(ed-3.0))/(3.0+p3),1.0/(1.0+cs2_3));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }
  double nb5=nb3*pow((5.0+p3+cs2_3*(5.0-3.0))/(3.0+p3),1.0/(1.0+cs2_3));

  // Add 4th high-density EOS
  for(double ed=5.0;ed<10.0-1.0e-4;ed+=0.2) {
    double pr=p5+e.params[7]*(ed-5.0)/(7.0-5.0);
    double nb=nb5*pow((ed+p5+cs2_4*(ed-5.0))/(5.0+p5),1.0/(1.0+cs2_4));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
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

void generic_quarks::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;

  eos_had_schematic &se=this->se;
  nstar_cold2 &cns=this->cns;

  // Set hadronic EOS from entry information
  se.comp=e.params[0];
  se.kprime=e.params[1];
  se.b=e.params[2]-se.a;
  se.gamma=e.params[3];
  se.kpp=0.0;

  // Compute low-density eos
  cns.nb_end=0.6;
  cns.set_eos(se);
  cns.calc_eos();
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->set_interp_type(itp_linear);

  tab_eos->add_constant("S",e.params[2]);
  tab_eos->add_constant("L",se.fesym_slope(0.16));

  // Double check that the table is non-empty (we have to do this
  // especially for the size_t index in the for loop below)
  if (tab_eos->get_nlines()==0) {
    O2SCL_ERR("Table empty in generic quarks.",exc_efailed);
  }

  // Transition between nuclear part and polytrope
  double ed1=e.params[4];

  // Boundary baryon density and pressure by interpolating
  // the table
  double nb1=tab_eos->interp("ed",ed1,"nb");
  double pr1=tab_eos->interp("ed",ed1,"pr");

  // Determine 1st polytrope coefficient
  double coeff1=pr1/pow(ed1,e.params[5]);
  
  // Store the transition pressure for later use
  tab_eos->add_constant("pr_pt",pr1);
  
  if (coeff1<0.0 || pr1<0.0) {
    scr_out << "Rejected: Negative polytrope coefficient or "
	    << "matching pressure (1)." << endl;
    success=bamr_class::ix_neg_pressure;
    return;
  }

  // Double check that there is no gap in density between
  // the low-density EOS and the first polytrope
  double ed_last=tab_eos->max("ed");
  if (ed_last<ed1) {
    scr_out << "Gap between low-density EOS and polytrope. " << endl;
    O2SCL_ERR("Gap between low-density EOS and polytrope.",exc_efailed);
  }

  // Remove rows beyond 1st transition
  for(size_t i=0;i<tab_eos->get_nlines();i++) {
    if ((*tab_eos)["ed"][i]>ed1) {
      tab_eos->delete_row(i);
      i=0;
    }
  }

  // Check that low-density EOS has statistics
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (1)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }
  
  // Transition between polytrope and quarks
  double ed_trans=e.params[6];
  double pr_trans=coeff1*pow(ed_trans,e.params[5]);
  double nb_trans=nb1*pow(ed_trans/ed1,e.params[5]/(e.params[5]-1.0))*
    pow((ed_trans+pr_trans)/(ed1+pr1),1.0/(1.0-e.params[5]));

  // Store the transition pressure for later use
  tab_eos->add_constant("pr_q",pr_trans);

  // Add first polytrope to table
  for(double ed=ed1;ed<=ed_trans;ed+=0.05) {
    double pr=coeff1*pow(ed,e.params[5]);
    double nb=nb1*pow(ed/ed1,e.params[5]/(e.params[5]-1.0))*
      pow((ed+pr)/(ed1+pr1),1.0/(1.0-e.params[5]));
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
  }

  // Check that matching didn't fail
  if (tab_eos->get_nlines()<3) {
    scr_out << "Rejected: Polytrope fit failed (2)." << endl;
    success=bamr_class::ix_no_eos_table;
    return;
  }

  // Quark EOS parameters
  double a2=e.params[7];
  double a4=e.params[8];
  
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
      success=bamr_class::ix_eos_solve_failed;
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
      success=bamr_class::ix_acausal;
      return;
    }

    double ed=bag+a2*musq+3.0*a4*musq*musq;
    double pr=-bag+a2*musq+a4*musq*musq;
    double nb=(2.0*a2*mu+4.0*a4*mu*musq)/3.0*nb_corr;
    double line[3]={ed,pr,nb};
    tab_eos->line_of_data(3,line);
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

void quark_star::compute_eos(entry &e, int &success, std::ofstream &scr_out) {
  
  success=bamr_class::ix_success;

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
    success=bamr_class::ix_eos_solve_failed;
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
  int ret=gmh.msolve(1,x,fmf);
  if (ret!=0) {
    scr_out << "Solver failed in qstar." << std::endl;
    success=bamr_class::ix_eos_solve_failed;
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
      success=bamr_class::ix_acausal;
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
	success=bamr_class::ix_param_mismatch;
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

// --------------------------------------------------------------

qmc_neut::qmc_neut() {
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
}

qmc_neut::~qmc_neut() {
}

void qmc_neut::low_limits(entry &e) {

  e.params[0]=12.7;
  e.params[1]=0.48;
  e.params[2]=1.0;
  e.params[3]=2.1;
  e.params[4]=0.2;
  e.params[5]=2.0;
  e.params[6]=0.2;
    
  return;
}

void qmc_neut::high_limits(entry &e) {
    
  e.params[0]=13.3;
  e.params[1]=0.52;
  e.params[2]=5.0;
  e.params[3]=2.5;
  e.params[4]=4.0;
  e.params[5]=8.0;
  e.params[6]=4.0;
    
  return;
}

string qmc_neut::param_name(size_t i) {
  if (i==0) return "a";
  else if (i==1) return "alpha";
  else if (i==2) return "b";
  else if (i==3) return "beta";
  else if (i==4) return "index1";
  else if (i==5) return "trans1";
  return "index2";
}

string qmc_neut::param_unit(size_t i) {
  if (i==0) return "MeV";
  else if (i==1) return ".";
  else if (i==2) return "MeV";
  else if (i==3) return ".";
  else if (i==4) return ".";
  else if (i==5) return "1/fm^4";
  return ".";
}

void qmc_neut::first_point(entry &e) {

  e.params[0]=1.276936e+01;
  e.params[1]=5.043647e-01;
  e.params[2]=4.584098e+00;
  e.params[3]=2.323736e+00;
  e.params[4]=0.576;
  e.params[5]=4.60;
  e.params[6]=1.21;

  return;
}

void qmc_neut::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;
  
  // Hack to start with a fresh table
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->clear_table();
  tab_eos->line_of_names("ed pr");
  tab_eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=e.params[0];
  double alpha=e.params[1];
  double b=e.params[2];
  double beta=e.params[3];
    
  double L=3.0*(a*alpha+b*beta);
  tab_eos->add_constant("S",(a+b+16.0)/hc_mev_fm);
  tab_eos->add_constant("L",L/hc_mev_fm);

  double ed=0.0, pr=0.0;

  double gauss=pdg.sample();
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
    tab_eos->line_of_data(2,line);
  }

  // Set values for the computation of the baryon density
  // from the last point of the QMC parameterization
  nb_n1=rho_trans;
  nb_e1=ed;

  // Check that the last Monte Carlo energy density
  // isn't greater than the transition between the
  // polytropes
  if (e.params[5]<ed) {
    scr_out << "First polytrope doesn't appear." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  // Polytropic index at higher densities
  double index1=e.params[4];

  // Compute coefficient given index
  double exp=1.0+1.0/index1;
  double coeff=pr/pow(ed,exp);

  // Compute stepsize in energy density
  double delta_ed=(e.params[5]-ed)/50.01;

  // Ensure first point is not the same as the last point
  // (this causes problems for the baryon density)
  ed+=delta_ed;

  // Add first polytrope to table
  for(;ed<e.params[5];ed+=delta_ed) {
    double line[2]={ed,coeff*pow(ed,exp)};
    tab_eos->line_of_data(2,line);
  }

  // Compute second coefficient given index
  double index2=e.params[6];
  double exp2=1.0+1.0/index2;

  // Compute coefficient
  double pr_last=coeff*pow(e.params[5],exp);
  double coeff2=pr_last/pow(e.params[5],exp2);

  // Add second polytrope to table
  delta_ed=(10.0-e.params[5])/50.01;
  for(ed=e.params[5];ed<10.0;ed+=delta_ed) {
    double line[2]={ed,coeff2*pow(ed,exp2)};
    tab_eos->line_of_data(2,line);
  }

  return;
}

// --------------------------------------------------------------

qmc_threep::qmc_threep() {
  rho0=0.16;
  rho_trans=0.16;
}

qmc_threep::~qmc_threep() {
}

void qmc_threep::low_limits(entry &e) {

  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  e.params[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  e.params[1]=0.47;
  e.params[2]=29.5;
  e.params[3]=30.0;
  
  e.params[4]=0.2;
  e.params[5]=0.75;
  e.params[6]=0.2;
  e.params[7]=0.75;
  e.params[8]=0.2;
    
  return;
}

void qmc_threep::high_limits(entry &e) {
    
  e.params[0]=13.5;
  e.params[1]=0.53;
  e.params[2]=36.1;
  e.params[3]=70.0;

  e.params[4]=8.0;
  e.params[5]=8.0;
  e.params[6]=8.0;
  e.params[7]=8.0;
  e.params[8]=8.0;
    
  return;
}

string qmc_threep::param_name(size_t i) {
  if (i==0) return "a";
  else if (i==1) return "alpha";
  else if (i==2) return "S";
  else if (i==3) return "L";
  else if (i==4) return "index1";
  else if (i==5) return "trans1";
  else if (i==6) return "index2";
  else if (i==7) return "trans2";
  return "index3";
}

string qmc_threep::param_unit(size_t i) {
  if (i==0) return "MeV";
  else if (i==1) return ".";
  else if (i==2) return "MeV";
  else if (i==3) return "MeV";
  else if (i==4) return ".";
  else if (i==5) return "1/fm^4";
  else if (i==6) return ".";
  else if (i==7) return "1/fm^4";
  return ".";
}

void qmc_threep::first_point(entry &e) {

  e.params[0]=13.0;
  e.params[1]=0.5;
  e.params[2]=32.0;
  e.params[3]=50.0;
  e.params[4]=0.5;
  e.params[5]=2.0;
  e.params[6]=0.5;
  e.params[7]=2.5;
  e.params[8]=1.0;

  return;
}

void qmc_threep::compute_eos(entry &e, int &success, ofstream &scr_out) {

  bool debug=false;

  success=bamr_class::ix_success;
  
  // Hack to start with a fresh table
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->clear_table();
  tab_eos->line_of_names("ed pr");
  tab_eos->set_interp_type(itp_linear);
  //tab_eos->set_unit("ed","1/fm^4");
  //tab_eos->set_unit("pr","1/fm^4");
  //tab_eos->set_unit("nb","1/fm^3");
  //tab_eos->line_of_names("ed pr nb");

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=e.params[0];
  double alpha=e.params[1];
  double Stmp=e.params[2];
  double Ltmp=e.params[3];

  /*
    This is based on limits from two lines, as in Jim and I's EPJA
    review. In (S,L) space, the lower line is (29,0) to (35,55),
    and the upper line is (26.5,0) to (33.5,100)
   */
  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range: " << Stmp << " " << Ltmp << endl;
    scr_out << 9.17*Stmp-266.0 << " " << 14.3*Stmp-379.0 << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }
  if (debug) {
    scr_out << "b=" << b << " beta=" << beta << endl;
  }

  tab_eos->add_constant("S",Stmp/hc_mev_fm);
  tab_eos->add_constant("L",Ltmp/hc_mev_fm);

  double index1=e.params[4];
  double exp1=1.0+1.0/index1;
  double trans1=e.params[5];
  double index2=e.params[6];
  double exp2=1.0+1.0/index2;
  double trans2=e.params[7];
  double index3=e.params[8];
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
    tab_eos->line_of_data(2,line);
    // double line[3]={ed,pr,rho};
    // tab_eos->line_of_data(3,line);
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
    success=bamr_class::ix_param_mismatch;
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
      cerr << coeff1 << " " << exp1 << " " << ed_last << " " << trans1 << endl;
      cerr << ed << " " << pr << " " << nb_n1 << " " << nb_e1 << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
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
      cerr << "Problem in model (3): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
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
      cerr << "Problem in model (4): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) scr_out << line[0] << " " << line[1] << endl;
  }

  return;
}

// --------------------------------------------------------------

qmc_fixp::qmc_fixp() {
  nb0=0.16;
  nb_trans=0.16;

  ed1=2.0;
  ed2=3.0;
  ed3=5.0;
  ed4=7.0;
}

qmc_fixp::~qmc_fixp() {
}

void qmc_fixp::low_limits(entry &e) {

  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  e.params[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  e.params[1]=0.47;
  e.params[2]=29.5;
  e.params[3]=30.0;
  
  e.params[4]=0.0;
  e.params[5]=0.0;
  e.params[6]=0.0;
  e.params[7]=0.0;
    
  return;
}

void qmc_fixp::high_limits(entry &e) {
    
  e.params[0]=13.5;
  e.params[1]=0.53;
  e.params[2]=36.1;
  e.params[3]=70.0;

  e.params[4]=0.3;
  e.params[5]=1.5;
  e.params[6]=2.5;
  e.params[7]=2.5;
    
  return;
}

string qmc_fixp::param_name(size_t i) {
  if (i==0) return "a";
  else if (i==1) return "alpha";
  else if (i==2) return "S";
  else if (i==3) return "L";
  else if (i==4) return "pres1";
  else if (i==5) return "pres2";
  else if (i==6) return "pres3";
  return "pres4";
}

string qmc_fixp::param_unit(size_t i) {
  if (i==0) return "MeV";
  else if (i==1) return ".";
  else if (i==2) return "MeV";
  else if (i==3) return "MeV";
  else if (i==4) return "1/fm^4";
  else if (i==5) return "1/fm^4";
  else if (i==6) return "1/fm^4";
  return "1/fm^4";
}

void qmc_fixp::first_point(entry &e) {

  e.params[0]=1.276936e+01;
  e.params[1]=5.043647e-01;
  e.params[2]=30.0;
  e.params[3]=40.0;
  e.params[4]=0.014;
  e.params[5]=0.74;
  e.params[6]=0.60;
  e.params[7]=1.84;

  return;
}

void qmc_fixp::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;
  bool debug=true;
  
  // Hack to start with a fresh table
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->clear_table();
  tab_eos->line_of_names("ed pr");
  tab_eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=e.params[0];
  double alpha=e.params[1];
  double Stmp=e.params[2];
  double Ltmp=e.params[3];

  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta || b<0.5) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }
  
  tab_eos->add_constant("S",Stmp/hc_mev_fm);
  tab_eos->add_constant("L",Ltmp/hc_mev_fm);

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
      cerr << "Problem in model (4): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
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
    success=bamr_class::ix_param_mismatch;
    return;
  }

  // Compute pressures on grid, p1 is the pressure at ed1
  double p1=pr_trans+e.params[4];
  // Variable p2 is the pressure at ed2
  double p2=p1+e.params[5];
  // Variable p3 is the pressure at ed3
  double p3=p2+e.params[6];
  
  // Add 1st high-density EOS
  if (debug) {
    cout << "First line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  double delta_ed=(ed1-ed_trans)/20.0;
  for(double ed=ed_trans+delta_ed;ed<ed1-1.0e-4;ed+=delta_ed) {
    double line[2]={ed,pr_trans+e.params[4]*(ed-ed_trans)/(ed1-ed_trans)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (5): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 2nd high-density EOS
  if (debug) {
    cout << "Second line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed1;ed<ed2-1.0e-4;ed+=(ed2-ed1)/10.0) {
    double line[2]={ed,p1+e.params[5]*(ed-ed1)/(ed2-ed1)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (6): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 3rd high-density EOS
  if (debug) {
    cout << "Third line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed2;ed<ed3-1.0e-4;ed+=(ed3-ed2)/10.0) {
    double line[2]={ed,p2+e.params[6]*(ed-ed2)/(ed3-ed2)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (7): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 4th high-density EOS
  if (debug) {
    cout << "Fourth line segment: " << endl;
    cout << "ed           pr" << endl;
  }
  for(double ed=ed3;ed<10.0-1.0e-4;ed+=(ed4-ed3)/10.0) {
    double line[2]={ed,p3+e.params[7]*(ed-ed3)/(ed4-ed3)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (8): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
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

qmc_twolines::qmc_twolines() {
  rho0=0.16;
  rho_trans=0.16;
}

qmc_twolines::~qmc_twolines() {
}

void qmc_twolines::low_limits(entry &e) {

  // The paper gives 12.7-13.4, we enlarge this to 12.5 to 13.5, and
  // this should allow S values as small as 28.5
  e.params[0]=12.5;
  // The paper gives 0.475 to 0.514, we enlarge this to 0.47 to 0.53
  e.params[1]=0.47;
  e.params[2]=29.5;
  e.params[3]=30.0;
  
  e.params[4]=0.0;
  e.params[5]=0.0;
  e.params[6]=0.0;
  e.params[7]=0.0;
    
  return;
}

void qmc_twolines::high_limits(entry &e) {
    
  e.params[0]=13.5;
  e.params[1]=0.53;
  e.params[2]=36.1;
  e.params[3]=70.0;

  e.params[4]=0.3;
  e.params[5]=1.5;
  e.params[6]=2.5;
  e.params[7]=2.5;
    
  return;
}

string qmc_twolines::param_name(size_t i) {
  if (i==0) return "a";
  else if (i==1) return "alpha";
  else if (i==2) return "S";
  else if (i==3) return "L";
  else if (i==4) return "pres1";
  else if (i==5) return "pres2";
  else if (i==6) return "pres3";
  return "pres4";
}

string qmc_twolines::param_unit(size_t i) {
  if (i==0) return "MeV";
  else if (i==1) return ".";
  else if (i==2) return "MeV";
  else if (i==3) return "MeV";
  else if (i==4) return "1/fm^4";
  else if (i==5) return "1/fm^4";
  else if (i==6) return "1/fm^4";
  return "1/fm^4";
}

void qmc_twolines::first_point(entry &e) {

  e.params[0]=1.276936e+01;
  e.params[1]=5.043647e-01;
  e.params[2]=30.0;
  e.params[3]=40.0;
  e.params[4]=0.014;
  e.params[5]=0.74;
  e.params[6]=0.60;
  e.params[7]=1.84;

  return;
}

void qmc_twolines::compute_eos(entry &e, int &success, ofstream &scr_out) {

  success=bamr_class::ix_success;
  bool debug=false;
  
  // Hack to start with a fresh table
  o2_shared_ptr<table_units<> >::type tab_eos=cns.get_eos_results();
  tab_eos->clear_table();
  tab_eos->line_of_names("ed pr");
  tab_eos->set_interp_type(itp_linear);

  // Add the QMC calculations over the suggested range, but go a
  // little bit lower in density to make sure we extend all the way
  // down to the crust

  double a=e.params[0];
  double alpha=e.params[1];
  double Stmp=e.params[2];
  double Ltmp=e.params[3];

  if (Ltmp<9.17*Stmp-266.0 || Ltmp>14.3*Stmp-379.0) {
    scr_out << "L out of range." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  double b=Stmp-16.0-a;
  double beta=(Ltmp/3.0-a*alpha)/b;
  if (b<=0.0 || beta<=0.0 || alpha>beta || b<0.5) {
    scr_out << "Parameter b=" << b << " or beta=" 
	    << beta << " out of range." << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }
  
  tab_eos->add_constant("S",Stmp/hc_mev_fm);
  tab_eos->add_constant("L",Ltmp/hc_mev_fm);

  double ed=0.0, pr=0.0, ed_last=0.0, nb_last=0.0, pr_last=0.0;

  double ed_trans=e.params[5];
  double ed1=e.params[7];
  if (ed_trans>ed1) {
    scr_out << "Transition densities misordered." << endl;
    scr_out << ed_trans << " " << ed1 << endl;
    success=bamr_class::ix_param_mismatch;
    return;
  }

  if (debug) {
    cout.setf(ios::scientific);
    cout << a << " " << alpha << " " << b << " " << beta << endl;
    cout << endl;
  }
  bool done=false;
  for(double rho=0.02;done==false;rho+=0.001) {
    double rho1=rho/rho0;
    double rho1a=pow(rho1,alpha);
    double rho1b=pow(rho1,beta);
    double ene=a*rho1a+b*rho1b;
    ed=rho*(ene/hc_mev_fm+o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_neutron));
    pr=rho*(a*alpha*rho1a+b*beta*rho1b)/hc_mev_fm;

    if (ed>ed_trans) {
      done=true;
    } else {
      double line[2]={ed,pr};
      if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
	cerr << "Problem in model (4): " << line[0] << " "
	     << line[1] << endl;
	exit(-1);
      }
      tab_eos->line_of_data(2,line);
      if (debug) cout << ed << " " << pr << endl;
      ed_last=ed;
      pr_last=pr;
      nb_last=rho;
      
      // Set values for the computation of the baryon density
      // from the last point of the QMC parameterization
      nb_n1=rho;
      nb_e1=ed;
    }
  }
  if (debug) cout << endl;

  // Compute pressures on grid, ed=2.0, 3.0, 5, 7 fm^{-4}
  double p2=pr_last+e.params[4];
  double p3=p2+e.params[6];
  
  // Add 1st high-density EOS
  double delta_ed=(ed1-ed_last)/20.0;
  for(double ed=ed_last+delta_ed;ed<ed1-1.0e-4;ed+=delta_ed) {
    double line[2]={ed,p2+e.params[4]*(ed-ed_last)/(ed1-ed_last)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (5): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  // Add 2nd high-density EOS
  for(double ed=ed1;ed<10.0-1.0e-4;ed+=(10.0-ed1)/10.0) {
    double line[2]={ed,p3+e.params[5]*(ed-ed1)/(10.0-ed1)};
    if (!gsl_finite(line[0]) || !gsl_finite(line[1])) {
      cerr << "Problem in model (6): " << line[0] << " "
	   << line[1] << endl;
      exit(-1);
    }
    tab_eos->line_of_data(2,line);
    if (debug) cout << line[0] << " " << line[1] << endl;
  }
  if (debug) cout << endl;

  return;
}

