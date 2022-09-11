#include "ns_pop.h"

using namespace std;
using namespace o2scl;

int main() {
    ofstream file;
    ns_pop nsp;
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
    }

    file.close();

    return 0;
}
