/* This file containes NS mass data from all populations with 
their associated 68% central limits */

#include "likelihood.h"

using namespace std;

void mass_data::load_data() {
  
  /* Data from NS-NS/DNS (some with asymmetric 68% limits)
     [source: Alsing et al. 2018] */
  name_ns.push_back("J0453+1559 pulsar"); id_ns.push_back("J0453p");
  mass_ns.push_back(1.559); uplim_ns.push_back(0.004); 
  lowlim_ns.push_back(0.004);

  name_ns.push_back("J0453+1559 comp."); id_ns.push_back("J0453c");
  mass_ns.push_back(1.174); uplim_ns.push_back(0.004); 
  lowlim_ns.push_back(0.004);

  name_ns.push_back("J1906+0746 pulsar"); id_ns.push_back("J1906p");
  mass_ns.push_back(1.291); uplim_ns.push_back(0.011); 
  lowlim_ns.push_back(0.011);

  name_ns.push_back("J1906+0746 comp."); id_ns.push_back("J1906c");
  mass_ns.push_back(1.322); uplim_ns.push_back(0.011); 
  lowlim_ns.push_back(0.011);

  name_ns.push_back("B1534+12 pulsar"); id_ns.push_back("B1534p"); 
  mass_ns.push_back(1.3332); uplim_ns.push_back(0.001); 
  lowlim_ns.push_back(0.001);

  name_ns.push_back("B1534+12 comp."); id_ns.push_back("B1534c");
  mass_ns.push_back(1.3452); uplim_ns.push_back(0.001); 
  lowlim_ns.push_back(0.001);

  name_ns.push_back("B1913+16 pulsar"); id_ns.push_back("B1913p");
  mass_ns.push_back(1.4398); uplim_ns.push_back(0.0002); 
  lowlim_ns.push_back(0.0002);

  name_ns.push_back("B1913+16 comp."); id_ns.push_back("B1913c");
  mass_ns.push_back(1.3886); uplim_ns.push_back(0.0002); 
  lowlim_ns.push_back(0.0002);

  name_ns.push_back("B2127+11C pulsar"); id_ns.push_back("B2127p");
  mass_ns.push_back(1.358); uplim_ns.push_back(0.01); 
  lowlim_ns.push_back(0.01);

  name_ns.push_back("B2127+11C comp."); id_ns.push_back("B2127c");
  mass_ns.push_back(1.354); uplim_ns.push_back(0.01); 
  lowlim_ns.push_back(0.01);

  name_ns.push_back("J0737-3039A"); id_ns.push_back("J0737A");
  mass_ns.push_back(1.3381); uplim_ns.push_back(0.0007); 
  lowlim_ns.push_back(0.0007);

  name_ns.push_back("J0737-3039B"); id_ns.push_back("J0737B");
  mass_ns.push_back(1.2489); uplim_ns.push_back(0.0007); 
  lowlim_ns.push_back(0.0007);

  name_ns.push_back("J1756-2251 pulsar"); id_ns.push_back("J1756p");
  mass_ns.push_back(1.312); uplim_ns.push_back(0.017); 
  lowlim_ns.push_back(0.017);

  name_ns.push_back("J1756-2251 comp."); id_ns.push_back("J1756c");
  mass_ns.push_back(1.258); uplim_ns.push_back(0.017); 
  lowlim_ns.push_back(0.017);

  name_ns.push_back("J1807-2500B pulsar"); id_ns.push_back("J1807p");
  mass_ns.push_back(1.3655); uplim_ns.push_back(0.0021); 
  lowlim_ns.push_back(0.0021);

  name_ns.push_back("J1807-2500B comp."); id_ns.push_back("J1807c");
  mass_ns.push_back(1.2064); uplim_ns.push_back(0.002); 
  lowlim_ns.push_back(0.002);

  // [source: Kiziltan et al. 2013]
  name_ns.push_back("J1518+4904 pulsar"); id_ns.push_back("J1518p");
  mass_ns.push_back(1.56); uplim_ns.push_back(0.13); 
  lowlim_ns.push_back(0.44);

  name_ns.push_back("J1518+4904 comp."); id_ns.push_back("J1518c");
  mass_ns.push_back(1.05); uplim_ns.push_back(0.45); 
  lowlim_ns.push_back(0.11);

  name_ns.push_back("J1811-1736 pulsar"); id_ns.push_back("J1811p");
  mass_ns.push_back(1.56); uplim_ns.push_back(0.24); 
  lowlim_ns.push_back(0.45);

  name_ns.push_back("J1811-1736 comp."); id_ns.push_back("J1811c");
  mass_ns.push_back(1.12); uplim_ns.push_back(0.47); 
  lowlim_ns.push_back(0.13);

  name_ns.push_back("J1829+2456 pulsar"); id_ns.push_back("J1829p");
  mass_ns.push_back(1.20); uplim_ns.push_back(0.12); 
  lowlim_ns.push_back(0.46);

  name_ns.push_back("J1829+2456 comp."); id_ns.push_back("J1829c");
  mass_ns.push_back(1.40); uplim_ns.push_back(0.46); 
  lowlim_ns.push_back(0.12);

  /* Data from NS-WD (some with asymmetric 68% limits)
     [source: Alsing et al. 2018, Kiziltan et al. 2013] */
  name_wd.push_back("J2045+3633"); id_wd.push_back("J2045");
  mass_wd.push_back(1.33); uplim_wd.push_back(0.3); 
  lowlim_wd.push_back(0.3);

  name_wd.push_back("J2053+4650"); id_wd.push_back("J2053");
  mass_wd.push_back(1.40); uplim_wd.push_back(0.21); 
  lowlim_wd.push_back(0.21);

  name_wd.push_back("J1713+0747"); id_wd.push_back("J1713"); 
  mass_wd.push_back(1.35); uplim_wd.push_back(0.07); 
  lowlim_wd.push_back(0.07);

  name_wd.push_back("B1855+09"); id_wd.push_back("B1855");
  mass_wd.push_back(1.37); uplim_wd.push_back(0.13); 
  lowlim_wd.push_back(0.13);

  name_wd.push_back("J0751+1807"); id_wd.push_back("J0751");
  mass_wd.push_back(1.72); uplim_wd.push_back(0.07); 
  lowlim_wd.push_back(0.07);

  name_wd.push_back("J1141-6545"); id_wd.push_back("J1141");
  mass_wd.push_back(1.27); uplim_wd.push_back(0.01); 
  lowlim_wd.push_back(0.01);

  name_wd.push_back("J1738+0333"); id_wd.push_back("J1738");
  mass_wd.push_back(1.47); uplim_wd.push_back(0.07); 
  lowlim_wd.push_back(0.07);

  name_wd.push_back("J1614-2230"); id_wd.push_back("J1614");
  mass_wd.push_back(1.908); uplim_wd.push_back(0.016); 
  lowlim_wd.push_back(0.016);

  name_wd.push_back("J0348+0432"); id_wd.push_back("J0348");
  mass_wd.push_back(2.01); uplim_wd.push_back(0.04); 
  lowlim_wd.push_back(0.04);

  name_wd.push_back("J2222-0137"); id_wd.push_back("J2222");
  mass_wd.push_back(1.76); uplim_wd.push_back(0.06); 
  lowlim_wd.push_back(0.06);

  name_wd.push_back("J2234+0611"); id_wd.push_back("J2234");
  mass_wd.push_back(1.393); uplim_wd.push_back(0.013); 
  lowlim_wd.push_back(0.013);

  name_wd.push_back("J1949+3106"); id_wd.push_back("J1949");
  mass_wd.push_back(1.47); uplim_wd.push_back(0.43); 
  lowlim_wd.push_back(0.43);

  name_wd.push_back("J1012+5307"); id_wd.push_back("J1012");
  mass_wd.push_back(1.83); uplim_wd.push_back(0.11); 
  lowlim_wd.push_back(0.11);

  name_wd.push_back("J0437-4715"); id_wd.push_back("J0437");
  mass_wd.push_back(1.44); uplim_wd.push_back(0.07); 
  lowlim_wd.push_back(0.07);

  name_wd.push_back("J1909-3744"); id_wd.push_back("J1909");
  mass_wd.push_back(1.48); uplim_wd.push_back(0.03); 
  lowlim_wd.push_back(0.03);

  name_wd.push_back("J1802-2124"); id_wd.push_back("J1802");
  mass_wd.push_back(1.24); uplim_wd.push_back(0.11); 
  lowlim_wd.push_back(0.11);

  name_wd.push_back("J1911-5958A"); id_wd.push_back("J1911");
  mass_wd.push_back(1.34); uplim_wd.push_back(0.08); 
  lowlim_wd.push_back(0.08);

  name_wd.push_back("J2043+1711"); id_wd.push_back("J2043");
  mass_wd.push_back(1.38); uplim_wd.push_back(0.13); 
  lowlim_wd.push_back(0.13);

  name_wd.push_back("J0337+1715"); id_wd.push_back("J0337");
  mass_wd.push_back(1.4378); uplim_wd.push_back(0.0013); 
  lowlim_wd.push_back(0.0013);

  name_wd.push_back("J1946+3417"); id_wd.push_back("J1946");
  mass_wd.push_back(1.828); uplim_wd.push_back(0.022); 
  lowlim_wd.push_back(0.022);

  name_wd.push_back("J1918-0642"); id_wd.push_back("J1918");
  mass_wd.push_back(1.29); uplim_wd.push_back(0.1); 
  lowlim_wd.push_back(0.1);

  name_wd.push_back("J1600-3053"); id_wd.push_back("J1600");
  mass_wd.push_back(2.3); uplim_wd.push_back(0.7); 
  lowlim_wd.push_back(0.7);

  // [source: Kiziltan et al. 2013]
  name_wd.push_back("J0621+1002"); id_wd.push_back("J0621");
  mass_wd.push_back(1.70); uplim_wd.push_back(0.10); 
  lowlim_wd.push_back(0.17);

  name_wd.push_back("B2303+46"); id_wd.push_back("B2303");
  mass_wd.push_back(1.38); uplim_wd.push_back(0.06); 
  lowlim_wd.push_back(0.10);

  name_wd.push_back("J0024-7204H"); id_wd.push_back("J0024");
  mass_wd.push_back(1.48); uplim_wd.push_back(0.03); 
  lowlim_wd.push_back(0.06);

  name_wd.push_back("J0514-4002A"); id_wd.push_back("J0514");
  mass_wd.push_back(1.49); uplim_wd.push_back(0.04); 
  lowlim_wd.push_back(0.27);

  name_wd.push_back("B1516+02B"); id_wd.push_back("B1516");
  mass_wd.push_back(2.10); uplim_wd.push_back(0.19); 
  lowlim_wd.push_back(0.19);

  name_wd.push_back("J1748-2446I"); id_wd.push_back("J1748I");
  mass_wd.push_back(1.91); uplim_wd.push_back(0.02); 
  lowlim_wd.push_back(0.10);

  name_wd.push_back("J1748-2446J"); id_wd.push_back("J1748J");
  mass_wd.push_back(1.79); uplim_wd.push_back(0.02); 
  lowlim_wd.push_back(0.10);

  name_wd.push_back("B1802-07"); id_wd.push_back("B1802");
  mass_wd.push_back(1.26); uplim_wd.push_back(0.08); 
  lowlim_wd.push_back(0.17);

  name_wd.push_back("B1911-5958A"); id_wd.push_back("B1911");
  mass_wd.push_back(1.40); uplim_wd.push_back(0.16); 
  lowlim_wd.push_back(0.10);

  // [source: Cromartie et al. 2020]
  name_wd.push_back("J0740+6620"); id_wd.push_back("J0740");
  mass_wd.push_back(2.14); uplim_wd.push_back(0.20); 
  lowlim_wd.push_back(0.18);

  /* Data from NS-MS (some from X-ray/Optical, all with symmetric 
  68% limits) [source: Alsing et al. 2018] */
  name_ms.push_back("4U1700-377"); id_ms.push_back("4U1700");
  mass_ms.push_back(1.96); lim_ms.push_back(0.19); 

  name_ms.push_back("Cyg X-2"); id_ms.push_back("CygX2");
  mass_ms.push_back(1.71); lim_ms.push_back(0.21); 

  name_ms.push_back("SMC X-1"); id_ms.push_back("SMCX1");
  mass_ms.push_back(1.21); lim_ms.push_back(0.12); 
  
  name_ms.push_back("Cen X-3"); id_ms.push_back("CenX3");
  mass_ms.push_back(1.57); lim_ms.push_back(0.16); 
  
  name_ms.push_back("XTE J2123-058"); id_ms.push_back("XTEJ2123");
  mass_ms.push_back(1.53); lim_ms.push_back(0.42); 
  
  name_ms.push_back("4U 1822-371"); id_ms.push_back("4U1822");
  mass_ms.push_back(1.96); lim_ms.push_back(0.36); 
  
  name_ms.push_back("OAO 1657-415"); id_ms.push_back("OAO1657");
  mass_ms.push_back(1.74); lim_ms.push_back(0.3);  
  
  name_ms.push_back("J013236.7+303228"); id_ms.push_back("J013236");
  mass_ms.push_back(2.0); lim_ms.push_back(0.4);  
  
  name_ms.push_back("Vela X-1"); id_ms.push_back("VelaX1");
  mass_ms.push_back(2.12); lim_ms.push_back(0.16); 
  
  name_ms.push_back("4U1538-522"); id_ms.push_back("4U1538");
  mass_ms.push_back(1.02); lim_ms.push_back(0.17); 
  
  name_ms.push_back("LMC X-4"); id_ms.push_back("LMCX4");
  mass_ms.push_back(1.57); lim_ms.push_back(0.11); 
  
  name_ms.push_back("Her X-1"); id_ms.push_back("HerX1");
  mass_ms.push_back(1.073); lim_ms.push_back(0.36); 
  
  name_ms.push_back("2S 0921-630"); id_ms.push_back("2S0921");
  mass_ms.push_back(1.44); lim_ms.push_back(0.1);  
  
  name_ms.push_back("EXO 1722-363"); id_ms.push_back("EXO1722");
  mass_ms.push_back(1.91); lim_ms.push_back(0.45); 
  
  name_ms.push_back("SAX J1802.7-2017"); id_ms.push_back("SAXJ1802");
  mass_ms.push_back(1.57); lim_ms.push_back(0.25); 
  
  name_ms.push_back("XTE J1855-026"); id_ms.push_back("XTEJ1855");
  mass_ms.push_back(1.41); lim_ms.push_back(0.24); 
  
  // Count the total number of stars in all populations
  this->n_stars = id_ns.size() + id_wd.size() + id_ms.size();

}
