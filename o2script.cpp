#include "ns_pop.h"

using namespace std;
using namespace o2scl;

int main() {
    ofstream file;

/*  ns_pop nsp;
    pop_data pd;
    vector<string> eos_pars;

    eos_pars.push_back("a");
    eos_pars.push_back("alpha");
    eos_pars.push_back("param_S");
    eos_pars.push_back("param_L");
    eos_pars.push_back("exp1");
    eos_pars.push_back("trans1");
    eos_pars.push_back("exp2");
    eos_pars.push_back("trans2");
    eos_pars.push_back("exp3");
    eos_pars.push_back("M_chirp_det");
    eos_pars.push_back("q");
    eos_pars.push_back("z_cdf");

    pd.load_data();
    nsp.get_param_info();
    
    file.open("o2plot.txt");

    // Trace plots
    file << "# Trace plots:" << endl;
    for (size_t i=0; i<eos_pars.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -ytitle " 
        << eos_pars[i] << " -xtitle steps -plot1 "
        << eos_pars[i] << " -show &" << endl;
    }
    for (size_t i=0; i<nsp.par_names.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -ytitle " 
        << nsp.par_names[i] << " -xtitle steps -plot1 "
        << nsp.par_names[i] << " -show &" << endl;
    }
    
    // Histogran plots
    file << "# Histogram plots:" << endl;
    for (size_t i=0; i<eos_pars.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -xtitle "
            << eos_pars[i] << " -ytitle frequency -to-hist "
            << eos_pars[i] << " $binsz -plot -show &" << endl;
    }
    for (size_t i=0; i<nsp.par_names.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -xtitle "
            << nsp.par_names[i] << " -ytitle frequency -to-hist "
            << nsp.par_names[i] << " $binsz -plot -show &" << endl;
    }

    // Likelihood plots
    file << "# Likelihood plots:" << endl;
    for (size_t i=0; i<eos_pars.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -function 1 n -xtitle "
            << eos_pars[i] << " -ytitle log_wgt -scatter " << eos_pars[i] 
            << " log_wgt n log_wgt -show &" << endl;
    }
    for (size_t i=0; i<nsp.par_names.size(); i++) {
        file << "o2graph -read $file markov_chain_0 -function 1 n -xtitle "
            << nsp.par_names[i] << " -ytitle log_wgt -scatter " 
            << nsp.par_names[i] << " log_wgt n log_wgt -show &" << endl;
    }*/

    file.open("o2plot.txt"); 
    file << "o2 -read pop.o2 -function ";
    double x=0.03;
    for (int i=0; i<100; i++){
        file << "\"(1/sqrt(2*3.1416)/width_NS)*exp(-0.5*(("
            << x << "-mean_NS)/width_NS)^2)*(1+erf(("
            << x << "-mean_NS)*(skewness_NS)/width_NS/sqrt(2)))\" "
            << "SN_NS_" << i << " -function ";
        x+=0.03;
    }
    x=0.03;
    for (int i=0; i<100; i++) {
        file << "\"(1/sqrt(2*3.1416)/width_WD)*exp(-0.5*(("
            << x << "-mean_WD)/width_WD)^2)*(1+erf(("
            << x << "-mean_WD)*(skewness_WD)/width_WD/sqrt(2)))\" "
            << "SN_WD_" << i << " -function ";
        x+=0.03;
    }
    x=0.03;
    for (int i=0; i<100; i++) {
        if (i<99) {
            file << "\"(1/sqrt(2*3.1416)/width_LM)*exp(-0.5*(("
            << x << "-mean_LM)/width_LM)^2)*(1+erf(("
            << x << "-mean_LM)*(skewness_LM)/width_LM/sqrt(2)))\" "
            << "SN_LM_" << i << " -function ";
            x+=0.03;    
        }
        else {
            file << "\"(1/sqrt(2*3.1416)/width_LM)*exp(-0.5*(("
            << x << "-mean_LM)/width_LM)^2)*(1+erf(("
            << x << "-mean_LM)*(skewness_LM)/width_LM/sqrt(2)))\" "
            << "SN_LM_" << i << " ";
        }
    }
    file << "-internal mdist.o2" << endl;
    file.close();
    return 0;
}
