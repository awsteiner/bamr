/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2013, Andrew W. Steiner
  
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
#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/table.h>
#include <o2scl/table3d.h>
#include <o2scl/format_float.h>
#include <o2scl/vec_stats.h>
#include <o2scl/contour.h>
#include <o2scl/hist.h>
#include <o2scl/hdf_io.h>
#include <o2scl/expect_val.h>
#include <o2scl/graph.h>
#include <o2scl/interp.h>
#include <o2scl/cli_readline.h>
#include <o2scl/vector.h>

// For hist_1d_ev
#include "misc.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_graph;

typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

/** \brief Desc
 */
class plot_class {
  
public:

  /// Desc
  double scale;
  /// Desc
  string out_file;
  /// Desc
  string xlabel;
  /// Desc
  string title;
  
  plot_class() {
    scale=1.0;
    xlabel="x-axis Label";
    title="title";
  }
  
  /** \brief Desc
   */
  int plot(std::vector<std::string> &sv, bool itive_com) {

    // object name
    string name=sv[1];

    // file list
    vector<string> files;

    // Form list of data files
    for(size_t i=2;i<sv.size();i++) files.push_back(sv[i]);
    size_t nf=files.size();

    // data
    vector<hist_ev_1d> data;
    data.resize(nf);
    
    // Averages
    ubmatrix avg, std_dev, avg_err;
    size_t m_block=0, m_per_block=0;

    // Read data files
    hdf_file hf;
    size_t nh=0;
    for(size_t i=0;i<nf;i++) {
      cout << "Opening file: " << files[i] << endl;
      hf.open(files[i]);
      
      cout << "Reading hist: " << name+".h" << endl;
      hdf_input(hf,data[i].h,name+".h");
      if (i==0) {
	nh=data[0].h.size();
	avg.resize(nf,nh);
	std_dev.resize(nf,nh);
	avg_err.resize(nf,nh);
      }
      cout << "Reading vector_ev: " << endl;
      hdf_input(hf,data[i].ev,name+".ev");
      cout << "Reading count: " << endl;
      hf.get_szt(name+".count",data[i].count);
      
      cout << "Closing HDF file: " << endl;
      hf.close();

      cout << "Compute file stats: " << endl;
      ubmatrix_row avg_row(avg,i);
      ubmatrix_row std_dev_row(std_dev,i);
      ubmatrix_row avg_err_row(avg_err,i);
      data[i].ev.current_avg_stats(avg_row,std_dev_row,
				   avg_err_row,m_block,m_per_block);
    }

    // Compute weighted average over all files
    std::vector<double> full_avg(nh);
    for(size_t i=0;i<nh;i++) {
      ubmatrix_column avg_col(avg,i);
      ubmatrix_column avge_col(avg_err,i);
      full_avg[i]=wvector_mean(nf,avg_col,avge_col);
    }

    // Get reps from hist
    std::vector<double> reps;
    for(size_t i=0;i<data[0].h.size();i++) {
      reps.push_back(data[0].h.get_rep_i(i));
    }

    // Compute extremal values for x-axis
    double xmin=reps[0]*scale, xmax=reps[0]*scale;
    for(size_t i=1;i<reps.size();i++) {
      if (reps[i]*scale<xmin) xmin=reps[i]*scale;
      if (reps[i]*scale>xmax) xmax=reps[i]*scale;
    }

    // Renormalize y-axis to peak = 1
    
    double norm=full_avg[0];
    for(size_t i=1;i<full_avg.size();i++) {
      if (full_avg[i]>norm) norm=full_avg[i];
    }
    for(size_t i=0;i<full_avg.size();i++) full_avg[i]/=norm;

    // Set endpoints to zero for vector_invert_enclosed_sum()
    full_avg[0]=0.0;
    full_avg[full_avg.size()-1]=0.0;
    
    // Compute total integral
    double total=vector_integ_linear(data[0].h.size(),reps,full_avg);
    cout << "Total: " << total << endl;
  
    std::vector<double> locs;
    double lev;

    // Get one-sigma error ranges
    vector_invert_enclosed_sum(0.68*total,reps.size(),reps,full_avg,lev);
    vector_find_level(lev,reps.size(),reps,full_avg,locs);
    double xlo1=locs[0], xhi1=locs[0];
    for(size_t i=0;i<locs.size();i++) {
      if (locs[i]<xlo1) xlo1=locs[i];
      if (locs[i]>xhi1) xhi1=locs[i];
    }

    // Get two-sigma error ranges
    vector_invert_enclosed_sum(0.95*total,reps.size(),reps,full_avg,lev);
    vector_find_level(lev,reps.size(),reps,full_avg,locs);
    double xlo2=locs[0], xhi2=locs[0];
    for(size_t i=0;i<locs.size();i++) {
      if (locs[i]<xlo2) xlo2=locs[i];
      if (locs[i]>xhi2) xhi2=locs[i];
    }
    
    cout << xlo2*scale << " " << xlo1*scale << " " 
	 << (xlo1+xhi1)/2.0*scale << " "
	 << xhi1*scale << " " << xhi2*scale << endl;

    TApplication *theApp=new TApplication("App",0,NULL);
    TCanvas *c1;
    TPad *p1;
    TH1 *th1;
    TGaxis *a1, *b1;
    double ymax=1.0;
    new_graph_ticks(c1,p1,th1,"c1","d1","e1",xmin-(xmax-xmin)/40.0,
		    -ymax/40.0,xmax+(xmax-xmin)/40.0,ymax*1.1,a1,b1);
    p1->SetRightMargin(0.02);
    p1->SetTopMargin(0.02);
    p1->SetLeftMargin(0.11);
    p1->SetBottomMargin(0.11);
    p1->Draw();

    TGraph *gr1=new TGraph(reps.size());
    for(size_t i=0;i<reps.size();i++) {
      gr1->SetPoint(i,reps[i]*scale,full_avg[i]);
    }
    gr1->SetName(name.c_str());
    gr1->Draw();

    TLatex tt;
    tt.SetTextAlign(22);
    tt.SetTextFont(132);
    tt.DrawLatex((xmin+xmax)/2.0,ymax*(-0.12),xlabel.c_str());
    tt.DrawLatex(xmin+(xmax-xmin)/4.0,0.8,title.c_str());
    tt.SetTextAngle(90);
    tt.DrawLatex(xmin-(xmax-xmin)/(8.5),ymax*0.5,"MC Weight");
    tt.SetTextAngle(0);

    c1->Draw();
  
    // Output to file or screen
    if (out_file.length()>0) {
      c1->Print(out_file.c_str());
    } else {
      theApp->Run(kTRUE);
    }

    return 0;
  }

  /** \brief Desc
   */
  void run(int argc, char *argv[]) {
  
    // ---------------------------------------
    // Specify command-line option object
    
    cli_readline cl;
    cl.prompt="bayes_mvsr> ";
    cl.gnu_intro=false;
    
    // ---------------------------------------
    // Set options
    
    static const int nopt=1;
    comm_option_s options[nopt]={
      {'p',"plot","Short desc.",2,-1,"[file]",((string)"Long ")+"desc.",
       new comm_option_mfptr<plot_class>(this,&plot_class::plot),
       cli::comm_option_both}
    };
    cl.set_comm_option_vec(nopt,options);

    // ---------------------------------------
    // Set parameters
    
    cli::parameter_double p_scale;
    p_scale.d=&scale;
    p_scale.help="";
    cl.par_list.insert(make_pair("scale",&p_scale));

    cli::parameter_string p_xlabel;
    p_xlabel.str=&xlabel;
    p_xlabel.help="";
    cl.par_list.insert(make_pair("xlabel",&p_xlabel));

    cli::parameter_string p_out_file;
    p_out_file.str=&out_file;
    p_out_file.help="";
    cl.par_list.insert(make_pair("out_file",&p_out_file));

    cli::parameter_string p_title;
    p_title.str=&title;
    p_title.help="";
    cl.par_list.insert(make_pair("title",&p_title));

    // ---------------------------------------
    // Process command-line arguments and run
    
    cl.run_auto(argc,argv);
    
    return;
  }

};

int main(int argc, char *argv[]) {
  cout.setf(ios::scientific);
  plot_class pc;
  pc.run(argc,argv);
  return 0;
}
