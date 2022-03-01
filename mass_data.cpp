/* This file containes NS mass data from all populations with 
their associated 68% central limits */

#include "likelihood.h"

using namespace std;

void like::load_data() {
  
  mass_data md;
  
  /* Data from NS-NS/DNS (some with asymmetric 68% limits)
     [source: Alsing et al. 2018] */
  md.name_ns.push_back("J0453+1559 pulsar"); md.id_ns.push_back("J0453p");
  md.mass_ns.push_back(1.559); md.uplim_ns.push_back(0.004); md.lowlim_ns.push_back(0.004);

  md.name_ns.push_back("J0453+1559 comp."); md.id_ns.push_back("J0453c");
  md.mass_ns.push_back(1.174); md.uplim_ns.push_back(0.004); md.lowlim_ns.push_back(0.004);

  md.name_ns.push_back("J1906+0746 pulsar"); md.id_ns.push_back("J1906p");
  md.mass_ns.push_back(1.291); md.uplim_ns.push_back(0.011); md.lowlim_ns.push_back(0.011);

  md.name_ns.push_back("J1906+0746 comp."); md.id_ns.push_back("J1906c");
  md.mass_ns.push_back(1.322); md.uplim_ns.push_back(0.011); md.lowlim_ns.push_back(0.011);

  md.name_ns.push_back("B1534+12 pulsar"); md.id_ns.push_back("B1534p"); 
  md.mass_ns.push_back(1.3332); md.uplim_ns.push_back(0.001); md.lowlim_ns.push_back(0.001);

  md.name_ns.push_back("B1534+12 comp."); md.id_ns.push_back("B1534c");
  md.mass_ns.push_back(1.3452); md.uplim_ns.push_back(0.001); md.lowlim_ns.push_back(0.001);

  md.name_ns.push_back("B1913+16 pulsar"); md.id_ns.push_back("B1913p");
  md.mass_ns.push_back(1.4398); md.uplim_ns.push_back(0.0002); md.lowlim_ns.push_back(0.0002);

  md.name_ns.push_back("B1913+16 comp."); md.id_ns.push_back("B1913c");
  md.mass_ns.push_back(1.3886); md.uplim_ns.push_back(0.0002); md.lowlim_ns.push_back(0.0002);

  md.name_ns.push_back("B2127+11C pulsar"); md.id_ns.push_back("B2127p");
  md.mass_ns.push_back(1.358); md.uplim_ns.push_back(0.01); md.lowlim_ns.push_back(0.01);

  md.name_ns.push_back("B2127+11C comp."); md.id_ns.push_back("B2127c");
  md.mass_ns.push_back(1.354); md.uplim_ns.push_back(0.01); md.lowlim_ns.push_back(0.01);

  md.name_ns.push_back("J0737-3039A"); md.id_ns.push_back("J0737A");
  md.mass_ns.push_back(1.3381); md.uplim_ns.push_back(0.0007); md.lowlim_ns.push_back(0.0007);

  md.name_ns.push_back("J0737-3039B"); md.id_ns.push_back("J0737B");
  md.mass_ns.push_back(1.2489); md.uplim_ns.push_back(0.0007); md.lowlim_ns.push_back(0.0007);

  md.name_ns.push_back("J1756-2251 pulsar"); md.id_ns.push_back("J1756p");
  md.mass_ns.push_back(1.312); md.uplim_ns.push_back(0.017); md.lowlim_ns.push_back(0.017);

  md.name_ns.push_back("J1756-2251 comp."); md.id_ns.push_back("J1756c");
  md.mass_ns.push_back(1.258); md.uplim_ns.push_back(0.017); md.lowlim_ns.push_back(0.017);

  md.name_ns.push_back("J1807-2500B pulsar"); md.id_ns.push_back("J1807p");
  md.mass_ns.push_back(1.3655); md.uplim_ns.push_back(0.0021); md.lowlim_ns.push_back(0.0021);

  md.name_ns.push_back("J1807-2500B comp."); md.id_ns.push_back("J1807c");
  md.mass_ns.push_back(1.2064); md.uplim_ns.push_back(0.002); md.lowlim_ns.push_back(0.002);

  // [source: Kiziltan et al. 2013]
  md.name_ns.push_back("J1518+4904 pulsar"); md.id_ns.push_back("J1518p");
  md.mass_ns.push_back(1.56); md.uplim_ns.push_back(0.13); md.lowlim_ns.push_back(0.44);

  md.name_ns.push_back("J1518+4904 comp."); md.id_ns.push_back("J1518c");
  md.mass_ns.push_back(1.05); md.uplim_ns.push_back(0.45); md.lowlim_ns.push_back(0.11);

  md.name_ns.push_back("J1811-1736 pulsar"); md.id_ns.push_back("J1811p");
  md.mass_ns.push_back(1.56); md.uplim_ns.push_back(0.24); md.lowlim_ns.push_back(0.45);

  md.name_ns.push_back("J1811-1736 comp."); md.id_ns.push_back("J1811c");
  md.mass_ns.push_back(1.12); md.uplim_ns.push_back(0.47); md.lowlim_ns.push_back(0.13);

  md.name_ns.push_back("J1829+2456 pulsar"); md.id_ns.push_back("J1829p");
  md.mass_ns.push_back(1.20); md.uplim_ns.push_back(0.12); md.lowlim_ns.push_back(0.46);

  md.name_ns.push_back("J1829+2456 comp."); md.id_ns.push_back("J1829c");
  md.mass_ns.push_back(1.40); md.uplim_ns.push_back(0.46); md.lowlim_ns.push_back(0.12);

  /* Data from NS-WD (some with asymmetric 68% limits)
     [source: Alsing et al. 2018, Kiziltan et al. 2013] */
  md.name_wd.push_back("J2045+3633"); md.id_wd.push_back("J2045");
  md.mass_wd.push_back(1.33); md.uplim_wd.push_back(0.3); md.lowlim_wd.push_back(0.3);

  md.name_wd.push_back("J2053+4650"); md.id_wd.push_back("J2053");
  md.mass_wd.push_back(1.40); md.uplim_wd.push_back(0.21); md.lowlim_wd.push_back(0.21);

  md.name_wd.push_back("J1713+0747"); md.id_wd.push_back("J1713"); 
  md.mass_wd.push_back(1.35); md.uplim_wd.push_back(0.07); md.lowlim_wd.push_back(0.07);

  md.name_wd.push_back("B1855+09"); md.id_wd.push_back("B1855");
  md.mass_wd.push_back(1.37); md.uplim_wd.push_back(0.13); md.lowlim_wd.push_back(0.13);

  md.name_wd.push_back("J0751+1807"); md.id_wd.push_back("J0751");
  md.mass_wd.push_back(1.72); md.uplim_wd.push_back(0.07); md.lowlim_wd.push_back(0.07);

  md.name_wd.push_back("J1141-6545"); md.id_wd.push_back("J1141");
  md.mass_wd.push_back(1.27); md.uplim_wd.push_back(0.01); md.lowlim_wd.push_back(0.01);

  md.name_wd.push_back("J1738+0333"); md.id_wd.push_back("J1738");
  md.mass_wd.push_back(1.47); md.uplim_wd.push_back(0.07); md.lowlim_wd.push_back(0.07);

  md.name_wd.push_back("J1614-2230"); md.id_wd.push_back("J1614");
  md.mass_wd.push_back(1.908); md.uplim_wd.push_back(0.016); md.lowlim_wd.push_back(0.016);

  md.name_wd.push_back("J0348+0432"); md.id_wd.push_back("J0348");
  md.mass_wd.push_back(2.01); md.uplim_wd.push_back(0.04); md.lowlim_wd.push_back(0.04);

  md.name_wd.push_back("J2222-0137"); md.id_wd.push_back("J2222");
  md.mass_wd.push_back(1.76); md.uplim_wd.push_back(0.06); md.lowlim_wd.push_back(0.06);

  md.name_wd.push_back("J2234+0611"); md.id_wd.push_back("J2234");
  md.mass_wd.push_back(1.393); md.uplim_wd.push_back(0.013); md.lowlim_wd.push_back(0.013);

  md.name_wd.push_back("J1949+3106"); md.id_wd.push_back("J1949");
  md.mass_wd.push_back(1.47); md.uplim_wd.push_back(0.43); md.lowlim_wd.push_back(0.43);

  md.name_wd.push_back("J1012+5307"); md.id_wd.push_back("J1012");
  md.mass_wd.push_back(1.83); md.uplim_wd.push_back(0.11); md.lowlim_wd.push_back(0.11);

  md.name_wd.push_back("J0437-4715"); md.id_wd.push_back("J0437");
  md.mass_wd.push_back(1.44); md.uplim_wd.push_back(0.07); md.lowlim_wd.push_back(0.07);

  md.name_wd.push_back("J1909-3744"); md.id_wd.push_back("J1909");
  md.mass_wd.push_back(1.48); md.uplim_wd.push_back(0.03); md.lowlim_wd.push_back(0.03);

  md.name_wd.push_back("J1802-2124"); md.id_wd.push_back("J1802");
  md.mass_wd.push_back(1.24); md.uplim_wd.push_back(0.11); md.lowlim_wd.push_back(0.11);

  md.name_wd.push_back("J1911-5958A"); md.id_wd.push_back("J1911");
  md.mass_wd.push_back(1.34); md.uplim_wd.push_back(0.08); md.lowlim_wd.push_back(0.08);

  md.name_wd.push_back("J2043+1711"); md.id_wd.push_back("J2043");
  md.mass_wd.push_back(1.38); md.uplim_wd.push_back(0.13); md.lowlim_wd.push_back(0.13);

  md.name_wd.push_back("J0337+1715"); md.id_wd.push_back("J0337");
  md.mass_wd.push_back(1.4378); md.uplim_wd.push_back(0.0013); md.lowlim_wd.push_back(0.0013);

  md.name_wd.push_back("J1946+3417"); md.id_wd.push_back("J1946");
  md.mass_wd.push_back(1.828); md.uplim_wd.push_back(0.022); md.lowlim_wd.push_back(0.022);

  md.name_wd.push_back("J1918-0642"); md.id_wd.push_back("J1918");
  md.mass_wd.push_back(1.29); md.uplim_wd.push_back(0.1); md.lowlim_wd.push_back(0.1);

  md.name_wd.push_back("J1600-3053"); md.id_wd.push_back("J1600");
  md.mass_wd.push_back(2.3); md.uplim_wd.push_back(0.7); md.lowlim_wd.push_back(0.7);

  // [source: Kiziltan et al. 2013]
  md.name_wd.push_back("J0621+1002"); md.id_wd.push_back("J0621");
  md.mass_wd.push_back(1.70); md.uplim_wd.push_back(0.10); md.lowlim_wd.push_back(0.17);

  md.name_wd.push_back("B2303+46"); md.id_wd.push_back("B2303");
  md.mass_wd.push_back(1.38); md.uplim_wd.push_back(0.06); md.lowlim_wd.push_back(0.10);

  md.name_wd.push_back("J0024-7204H"); md.id_wd.push_back("J0024");
  md.mass_wd.push_back(1.48); md.uplim_wd.push_back(0.03); md.lowlim_wd.push_back(0.06);

  md.name_wd.push_back("J0514-4002A"); md.id_wd.push_back("J0514");
  md.mass_wd.push_back(1.49); md.uplim_wd.push_back(0.04); md.lowlim_wd.push_back(0.27);

  md.name_wd.push_back("B1516+02B"); md.id_wd.push_back("B1516");
  md.mass_wd.push_back(2.10); md.uplim_wd.push_back(0.19); md.lowlim_wd.push_back(0.19);

  md.name_wd.push_back("J1748-2446I"); md.id_wd.push_back("J1748I");
  md.mass_wd.push_back(1.91); md.uplim_wd.push_back(0.02); md.lowlim_wd.push_back(0.10);

  md.name_wd.push_back("J1748-2446J"); md.id_wd.push_back("J1748J");
  md.mass_wd.push_back(1.79); md.uplim_wd.push_back(0.02); md.lowlim_wd.push_back(0.10);

  md.name_wd.push_back("B1802-07"); md.id_wd.push_back("B1802");
  md.mass_wd.push_back(1.26); md.uplim_wd.push_back(0.08); md.lowlim_wd.push_back(0.17);

  md.name_wd.push_back("B1911-5958A"); md.id_wd.push_back("B1911");
  md.mass_wd.push_back(1.40); md.uplim_wd.push_back(0.16); md.lowlim_wd.push_back(0.10);

  // [source: Cromartie et al. 2020]
  md.name_wd.push_back("J0740+6620"); md.id_wd.push_back("J0740");
  md.mass_wd.push_back(2.14); md.uplim_wd.push_back(0.20); md.lowlim_wd.push_back(0.18);

  /* Data from NS-MS (some from X-ray/Optical, all with symmetric 68% limits)
     [source: Alsing et al. 2018] */
  md.name_ms.push_back("4U1700-377"); md.id_ms.push_back("4U1700");
  md.mass_ms.push_back(1.96); md.lim_ms.push_back(0.19); 

  md.name_ms.push_back("Cyg X-2"); md.id_ms.push_back("CygX2");
  md.mass_ms.push_back(1.71); md.lim_ms.push_back(0.21); 

  md.name_ms.push_back("SMC X-1"); md.id_ms.push_back("SMCX1");
  md.mass_ms.push_back(1.21); md.lim_ms.push_back(0.12); 
  
  md.name_ms.push_back("Cen X-3"); md.id_ms.push_back("CenX3");
  md.mass_ms.push_back(1.57); md.lim_ms.push_back(0.16); 
  
  md.name_ms.push_back("XTE J2123-058"); md.id_ms.push_back("XTEJ2123");
  md.mass_ms.push_back(1.53); md.lim_ms.push_back(0.42); 
  
  md.name_ms.push_back("4U 1822-371"); md.id_ms.push_back("4U1822");
  md.mass_ms.push_back(1.96); md.lim_ms.push_back(0.36); 
  
  md.name_ms.push_back("OAO 1657-415"); md.id_ms.push_back("OAO1657");
  md.mass_ms.push_back(1.74); md.lim_ms.push_back(0.3);  
  
  md.name_ms.push_back("J013236.7+303228"); md.id_ms.push_back("J013236");
  md.mass_ms.push_back(2.0); md.lim_ms.push_back(0.4);  
  
  md.name_ms.push_back("Vela X-1"); md.id_ms.push_back("VelaX1");
  md.mass_ms.push_back(2.12); md.lim_ms.push_back(0.16); 
  
  md.name_ms.push_back("4U1538-522"); md.id_ms.push_back("4U1538");
  md.mass_ms.push_back(1.02); md.lim_ms.push_back(0.17); 
  
  md.name_ms.push_back("LMC X-4"); md.id_ms.push_back("LMCX4");
  md.mass_ms.push_back(1.57); md.lim_ms.push_back(0.11); 
  
  md.name_ms.push_back("Her X-1"); md.id_ms.push_back("HerX1");
  md.mass_ms.push_back(1.073); md.lim_ms.push_back(0.36); 
  
  md.name_ms.push_back("2S 0921-630"); md.id_ms.push_back("2S0921");
  md.mass_ms.push_back(1.44); md.lim_ms.push_back(0.1);  
  
  md.name_ms.push_back("EXO 1722-363"); md.id_ms.push_back("EXO1722");
  md.mass_ms.push_back(1.91); md.lim_ms.push_back(0.45); 
  
  md.name_ms.push_back("SAX J1802.7-2017"); md.id_ms.push_back("SAXJ1802");
  md.mass_ms.push_back(1.57); md.lim_ms.push_back(0.25); 
  
  md.name_ms.push_back("XTE J1855-026"); md.id_ms.push_back("XTEJ1855");
  md.mass_ms.push_back(1.41); md.lim_ms.push_back(0.24); 
  
  /* Data from NS-MS (all with symmetric 68% limits)
     [source: Alsing et al. 2018]
     nsms_rd[0][0] = 1.58;   nsms_rd[0][1] = 0.34; 
     nsms_id[0] = "J0045-7319");
     nsms_rd[1][0] = 1.71;   nsms_rd[1][1] = 0.16;  
     nsms_id[1] = "J1023+0038");
     nsms_rd[2][0] = 1.666;  nsms_rd[2][1] = 0.01;  
     nsms_id[2] = "J1903+0327"); */
}