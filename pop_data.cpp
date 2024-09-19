/* This file containes NS mass data from all populations with 
their associated 68% central limits */

#include "ns_pop.h"

using namespace std;

void pop_data::load_data() {
  
  /* Data from NS-NS/DNS (some with asymmetric 68% limits)
     [source: Alsing et al. 2018] */
  
  vector<double> low_ns, high_ns;

  name_ns.push_back("J0453+1559 pulsar"); // 34
  id_ns.push_back("J0453p");
  mass_ns.push_back(1.559);
  high_ns.push_back(0.004);
  low_ns.push_back(0.004);

  name_ns.push_back("J0453+1559 comp."); // 35
  id_ns.push_back("J0453c");
  mass_ns.push_back(1.174);
  high_ns.push_back(0.004);
  low_ns.push_back(0.004);

  name_ns.push_back("J1906+0746 pulsar"); // 36
  id_ns.push_back("J1906p");
  mass_ns.push_back(1.291);
  high_ns.push_back(0.011);
  low_ns.push_back(0.011);

  name_ns.push_back("J1906+0746 comp."); // 37
  id_ns.push_back("J1906c");
  mass_ns.push_back(1.322);
  high_ns.push_back(0.011);
  low_ns.push_back(0.011);

  name_ns.push_back("B1534+12 pulsar"); // 38
  id_ns.push_back("B1534p");
  mass_ns.push_back(1.3332);
  high_ns.push_back(0.001);
  low_ns.push_back(0.001);

  name_ns.push_back("B1534+12 comp."); // 39
  id_ns.push_back("B1534c");
  mass_ns.push_back(1.3452);
  high_ns.push_back(0.001);
  low_ns.push_back(0.001);

  name_ns.push_back("B1913+16 pulsar"); // 40
  id_ns.push_back("B1913p");
  mass_ns.push_back(1.4398);
  high_ns.push_back(0.0002);
  low_ns.push_back(0.0002);

  name_ns.push_back("B1913+16 comp."); // 41
  id_ns.push_back("B1913c");
  mass_ns.push_back(1.3886);
  high_ns.push_back(0.0002);
  low_ns.push_back(0.0002);

  name_ns.push_back("B2127+11C pulsar"); // 42
  id_ns.push_back("B2127p");
  mass_ns.push_back(1.358);
  high_ns.push_back(0.01);
  low_ns.push_back(0.01);

  name_ns.push_back("B2127+11C comp."); // 43
  id_ns.push_back("B2127c");
  mass_ns.push_back(1.354);
  high_ns.push_back(0.01);
  low_ns.push_back(0.01);

  name_ns.push_back("J0737-3039A"); // 44
  id_ns.push_back("J0737A");
  mass_ns.push_back(1.3381);
  high_ns.push_back(0.0007);
  low_ns.push_back(0.0007);

  name_ns.push_back("J0737-3039B"); // 45
  id_ns.push_back("J0737B");
  mass_ns.push_back(1.2489);
  high_ns.push_back(0.0007);
  low_ns.push_back(0.0007);

  name_ns.push_back("J1756-2251 pulsar"); // 46
  id_ns.push_back("J1756p");
  mass_ns.push_back(1.312);
  high_ns.push_back(0.017);
  low_ns.push_back(0.017);

  name_ns.push_back("J1756-2251 comp."); // 47
  id_ns.push_back("J1756c");
  mass_ns.push_back(1.258);
  high_ns.push_back(0.017);
  low_ns.push_back(0.017);

  name_ns.push_back("J1807-2500B pulsar"); // 48
  id_ns.push_back("J1807p");
  mass_ns.push_back(1.3655);
  high_ns.push_back(0.0021);
  low_ns.push_back(0.0021);

  name_ns.push_back("J1807-2500B comp."); // 49
  id_ns.push_back("J1807c");
  mass_ns.push_back(1.2064);
  high_ns.push_back(0.002);
  low_ns.push_back(0.002);

  // [source: Kiziltan et al. 2013]
  name_ns.push_back("J1518+4904 pulsar"); // 50
  id_ns.push_back("J1518p");
  mass_ns.push_back(1.56);
  high_ns.push_back(0.13);
  low_ns.push_back(0.44);

  name_ns.push_back("J1518+4904 comp."); // 51
  id_ns.push_back("J1518c");
  mass_ns.push_back(1.05);
  high_ns.push_back(0.45);
  low_ns.push_back(0.11);

  name_ns.push_back("J1811-1736 pulsar"); // 52
  id_ns.push_back("J1811p");
  mass_ns.push_back(1.56);
  high_ns.push_back(0.24);
  low_ns.push_back(0.45);

  name_ns.push_back("J1811-1736 comp."); // 53
  id_ns.push_back("J1811c");
  mass_ns.push_back(1.12);
  high_ns.push_back(0.47);
  low_ns.push_back(0.13);

  name_ns.push_back("J1829+2456 pulsar"); // 54
  id_ns.push_back("J1829p");
  mass_ns.push_back(1.20);
  high_ns.push_back(0.12);
  low_ns.push_back(0.46);

  name_ns.push_back("J1829+2456 comp."); // 55
  id_ns.push_back("J1829c");
  mass_ns.push_back(1.40);
  high_ns.push_back(0.46);
  low_ns.push_back(0.12);

  /* Data from NS-WD (some with asymmetric 68% limits)
     [source: Alsing et al. 2018, Kiziltan et al. 2013] */

  vector<double> low_wd, high_wd;

  name_wd.push_back("J2045+3633"); // 56
  id_wd.push_back("J2045");
  mass_wd.push_back(1.33);
  high_wd.push_back(0.3);
  low_wd.push_back(0.3);

  name_wd.push_back("J2053+4650"); // 57
  id_wd.push_back("J2053");
  mass_wd.push_back(1.40);
  high_wd.push_back(0.21);
  low_wd.push_back(0.21);

  name_wd.push_back("J1713+0747"); // 58
  id_wd.push_back("J1713");
  mass_wd.push_back(1.35);
  high_wd.push_back(0.07);
  low_wd.push_back(0.07);

  name_wd.push_back("B1855+09"); // 59
  id_wd.push_back("B1855");
  mass_wd.push_back(1.37);
  high_wd.push_back(0.13); 
  low_wd.push_back(0.13);

  name_wd.push_back("J0751+1807"); // 60
  id_wd.push_back("J0751");
  mass_wd.push_back(1.72);
  high_wd.push_back(0.07);
  low_wd.push_back(0.07);

  name_wd.push_back("J1141-6545"); // 61
  id_wd.push_back("J1141");
  mass_wd.push_back(1.27);
  high_wd.push_back(0.01);
  low_wd.push_back(0.01);

  name_wd.push_back("J1738+0333"); // 62
  id_wd.push_back("J1738");
  mass_wd.push_back(1.47);
  high_wd.push_back(0.07);
  low_wd.push_back(0.07);

  name_wd.push_back("J1614-2230"); // 63
  id_wd.push_back("J1614");
  mass_wd.push_back(1.908);
  high_wd.push_back(0.016);
  low_wd.push_back(0.016);

  name_wd.push_back("J0348+0432"); // 64
  id_wd.push_back("J0348");
  mass_wd.push_back(2.01);
  high_wd.push_back(0.04);
  low_wd.push_back(0.04);

  name_wd.push_back("J2222-0137"); // 65
  id_wd.push_back("J2222");
  mass_wd.push_back(1.76);
  high_wd.push_back(0.06);
  low_wd.push_back(0.06);

  name_wd.push_back("J2234+0611"); // 66
  id_wd.push_back("J2234");
  mass_wd.push_back(1.393);
  high_wd.push_back(0.013);
  low_wd.push_back(0.013);

  name_wd.push_back("J1949+3106"); // 67
  id_wd.push_back("J1949");
  mass_wd.push_back(1.47);
  high_wd.push_back(0.43);
  low_wd.push_back(0.43);

  name_wd.push_back("J1012+5307"); // 68
  id_wd.push_back("J1012");
  mass_wd.push_back(1.83);
  high_wd.push_back(0.11);
  low_wd.push_back(0.11);

  name_wd.push_back("J0437-4715"); // 69
  id_wd.push_back("J0437");
  mass_wd.push_back(1.44);
  high_wd.push_back(0.07);
  low_wd.push_back(0.07);

  name_wd.push_back("J1909-3744"); // 70
  id_wd.push_back("J1909");
  mass_wd.push_back(1.48);
  high_wd.push_back(0.03);
  low_wd.push_back(0.03);

  name_wd.push_back("J1802-2124"); // 71
  id_wd.push_back("J1802");
  mass_wd.push_back(1.24);
  high_wd.push_back(0.11);
  low_wd.push_back(0.11);

  name_wd.push_back("J1911-5958A"); // 72
  id_wd.push_back("J1911");
  mass_wd.push_back(1.34);
  high_wd.push_back(0.08);
  low_wd.push_back(0.08);

  name_wd.push_back("J2043+1711"); // 73
  id_wd.push_back("J2043");
  mass_wd.push_back(1.38);
  high_wd.push_back(0.13);
  low_wd.push_back(0.13);

  name_wd.push_back("J0337+1715"); // 74
  id_wd.push_back("J0337");
  mass_wd.push_back(1.4378);
  high_wd.push_back(0.0013);
  low_wd.push_back(0.0013);

  name_wd.push_back("J1946+3417"); // 75
  id_wd.push_back("J1946");
  mass_wd.push_back(1.828);
  high_wd.push_back(0.022);
  low_wd.push_back(0.022);

  name_wd.push_back("J1918-0642"); // 76
  id_wd.push_back("J1918");
  mass_wd.push_back(1.29);
  high_wd.push_back(0.1);
  low_wd.push_back(0.1);

  name_wd.push_back("J1600-3053"); // 77
  id_wd.push_back("J1600");
  mass_wd.push_back(2.3);
  high_wd.push_back(0.7);
  low_wd.push_back(0.7);

  // [source: Kiziltan et al. 2013]
  name_wd.push_back("J0621+1002"); // 78
  id_wd.push_back("J0621");
  mass_wd.push_back(1.70);
  high_wd.push_back(0.10);
  low_wd.push_back(0.17);

  name_wd.push_back("B2303+46"); // 79
  id_wd.push_back("B2303");
  mass_wd.push_back(1.38);
  high_wd.push_back(0.06);
  low_wd.push_back(0.10);

  name_wd.push_back("J0024-7204H"); // 80
  id_wd.push_back("J0024");
  mass_wd.push_back(1.48);
  high_wd.push_back(0.03);
  low_wd.push_back(0.06);

  name_wd.push_back("J0514-4002A"); // 81
  id_wd.push_back("J0514");
  mass_wd.push_back(1.49);
  high_wd.push_back(0.04);
  low_wd.push_back(0.27);

  name_wd.push_back("B1516+02B"); // 82
  id_wd.push_back("B1516");
  mass_wd.push_back(2.10);
  high_wd.push_back(0.19);
  low_wd.push_back(0.19);

  name_wd.push_back("J1748-2446I"); // 83
  id_wd.push_back("J1748I");
  mass_wd.push_back(1.91);
  high_wd.push_back(0.02);
  low_wd.push_back(0.10);

  name_wd.push_back("J1748-2446J"); // 84
  id_wd.push_back("J1748J");
  mass_wd.push_back(1.79);
  high_wd.push_back(0.02);
  low_wd.push_back(0.10);

  name_wd.push_back("B1802-07"); // 85
  id_wd.push_back("B1802");
  mass_wd.push_back(1.26);
  high_wd.push_back(0.08);
  low_wd.push_back(0.17);

  name_wd.push_back("B1911-5958A"); // 86
  id_wd.push_back("B1911");
  mass_wd.push_back(1.40);
  high_wd.push_back(0.16);
  low_wd.push_back(0.10);

  // [source: Cromartie et al. 2020]
  name_wd.push_back("J0740+6620"); // 87
  id_wd.push_back("J0740");
  mass_wd.push_back(2.14);
  high_wd.push_back(0.20);
  low_wd.push_back(0.18);

  /* Data from NS-MS (2 types of X-ray/optical binaries: HMXB
  and LMXB (includes QLMXBs), all with symmetric 68% limits) 
  [source: Alsing et al. 2018] */
  
  // High Mass X-ray Binaries (HMXBs)
  /* name_hx.push_back("4U1700-377"); id_hx.push_back("4U1700");
  mass_hx.push_back(1.96); low_hx.push_back(0.19); 

  name_hx.push_back("SMC X-1"); id_hx.push_back("SMCX1");
  mass_hx.push_back(1.21); low_hx.push_back(0.12); 
  
  name_hx.push_back("Cen X-3"); id_hx.push_back("CenX3");
  mass_hx.push_back(1.57); low_hx.push_back(0.16); 

  name_hx.push_back("OAO 1657-415"); id_hx.push_back("OAO1657");
  mass_hx.push_back(1.74); low_hx.push_back(0.3);  
  
  name_hx.push_back("J013236.7+303228"); id_hx.push_back("J013236");
  mass_hx.push_back(2.0); low_hx.push_back(0.4);  
  
  name_hx.push_back("Vela X-1"); id_hx.push_back("VelaX1");
  mass_hx.push_back(2.12); low_hx.push_back(0.16); 
  
  name_hx.push_back("4U1538-522"); id_hx.push_back("4U1538");
  mass_hx.push_back(1.02); low_hx.push_back(0.17); 
  
  name_hx.push_back("LMC X-4"); id_hx.push_back("LMCX4");
  mass_hx.push_back(1.57); low_hx.push_back(0.11);
  
  name_hx.push_back("EXO 1722-363"); id_hx.push_back("EXO1722");
  mass_hx.push_back(1.91); low_hx.push_back(0.45); 
  
  name_hx.push_back("SAX J1802.7-2017"); id_hx.push_back("SAXJ1802");
  mass_hx.push_back(1.57); low_hx.push_back(0.25); 
  
  name_hx.push_back("XTE J1855-026"); id_hx.push_back("XTEJ1855");
  mass_hx.push_back(1.41); low_hx.push_back(0.24); */

  /* Data from NS-MS (2 types of X-ray/optical binaries: HMXB
  and LMXB (includes QLMXBs), all with symmetric 68% limits) 
  [source: Alsing et al. 2018] */

  vector<double> low_lx, high_lx;

  // Low Mass X-ray Binaries (LMXBs)
  name_lx.push_back("Cyg X-2"); // 88
  id_lx.push_back("CygX2");
  mass_lx.push_back(1.71);
  low_lx.push_back(0.21);
  
  name_lx.push_back("XTE J2123-058"); // 89
  id_lx.push_back("XTEJ2123");
  mass_lx.push_back(1.53);
  low_lx.push_back(0.42);
  
  name_lx.push_back("4U 1822-371"); // 90
  id_lx.push_back("4U1822");
  mass_lx.push_back(1.96);
  low_lx.push_back(0.36);
  
  name_lx.push_back("Her X-1"); // 91
  id_lx.push_back("HerX1");
  mass_lx.push_back(1.073);
  low_lx.push_back(0.36);
  
  name_lx.push_back("2S 0921-630"); // 92
  id_lx.push_back("2S0921");
  mass_lx.push_back(1.44);
  low_lx.push_back(0.1);

  // Symmetric 68% limits for LMXBs
  for (size_t i=0; i<low_lx.size(); i++) {
    high_lx.push_back(low_lx[i]);
  }

  // Fill input data and compute hyperparameters 
  // for all NS populations
  
  eqn_solver es;
  
  for (size_t i=0; i<name_ns.size(); i++) {
    double c=sqrt(high_ns[i]/low_ns[i]);
    double d=es.get_scale(low_ns[i], high_ns[i]);
    id_nsp.push_back(id_ns[i]);
    mass_nsp.push_back(mass_ns[i]);
    asym_ns.push_back(c);
    asym_nsp.push_back(c);
    scale_ns.push_back(d);
    scale_nsp.push_back(d);
    lo_nsp.push_back(low_ns[i]);
    hi_nsp.push_back(high_ns[i]);
  }
  for (size_t i=0; i<name_wd.size(); i++) {
    double c=sqrt(high_wd[i]/low_wd[i]);
    double d=es.get_scale(low_wd[i], high_wd[i]);
    id_nsp.push_back(id_wd[i]);
    mass_nsp.push_back(mass_wd[i]);
    asym_wd.push_back(c);
    asym_nsp.push_back(c);
    scale_wd.push_back(d);
    scale_nsp.push_back(d);
    lo_nsp.push_back(low_wd[i]);
    hi_nsp.push_back(high_wd[i]);
  }
  for (size_t i=0; i<name_lx.size(); i++) {
    double c=sqrt(high_lx[i]/low_lx[i]);
    double d=es.get_scale(low_lx[i], high_lx[i]);
    id_nsp.push_back(id_lx[i]);
    mass_nsp.push_back(mass_lx[i]);
    asym_lx.push_back(c);
    asym_nsp.push_back(c);
    scale_lx.push_back(d);
    scale_nsp.push_back(d);
    lo_nsp.push_back(low_lx[i]);
    hi_nsp.push_back(high_lx[i]);
  }
  
  // Count the total number of stars in all populations
  n_dns=id_ns.size();
  n_nswd=id_wd.size();
  n_lmxb=id_lx.size();

  // Currently, 59 = 22 + 32 + 5
  this->n_stars=n_dns+n_nswd+n_lmxb;
}
