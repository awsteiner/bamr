#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <string>
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/constants.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>
#include <boost/numeric/ublas/vector.hpp>
#include <gsl/gsl_math.h>

using namespace std;

// Arrays to store the names/identifiers of the stars
string nsns_id[26], nswd_id[39], nsms_id[16];

/* Arrays to store the raw data in the format: [mi][li, ui] for NS-NS
and NS-WD, and [mi][si] for NS-MS (where si=li=ui) */
double nsns[26][3], nswd[39][3], nsms[16][2];

void get_data() {
  /* Data from NS-NS/DNS (some with asymmetric 68% limits)
     [source: Alsing et al. 2018] */
  nsns_id[0] = "J0453+1559 pulsar";
  nsns[0][0] = 1.559; nsns[0][1] = 0.004; nsns[0][2] = 0.004;

  nsns_id[1] = "J0453+1559 comp.";
  nsns[1][0] = 1.174; nsns[1][1] = 0.004; nsns[1][2] = 0.004;

  nsns_id[2] = "J1906+0746 pulsar";
  nsns[2][0] = 1.291; nsns[2][1] = 0.011; nsns[2][2] = 0.011;

  nsns_id[3] = "J1906+0746 comp.";
  nsns[3][0] = 1.322; nsns[3][1] = 0.011; nsns[3][2] = 0.011;

  nsns_id[4] = "B1534+12 pulsar";
  nsns[4][0] = 1.3332; nsns[4][1] = 0.001; nsns[4][2] = 0.001;

  nsns_id[5] = "B1534+12 comp.";
  nsns[5][0] = 1.3452; nsns[5][1] = 0.001; nsns[5][2] = 0.001;

  nsns_id[6] = "B1913+16 pulsar";
  nsns[6][0] = 1.4398; nsns[6][1] = 0.0002; nsns[6][2] = 0.0002;

  nsns_id[7] = "B1913+16 comp.";
  nsns[7][0] = 1.3886; nsns[7][1] = 0.0002; nsns[7][2] = 0.0002;

  nsns_id[8] = "B2127+11C pulsar";
  nsns[8][0] = 1.358; nsns[8][1] = 0.01; nsns[8][2] = 0.01;

  nsns_id[9] = "B2127+11C comp.";
  nsns[9][0] = 1.354; nsns[9][1] = 0.01; nsns[9][2] = 0.01;

  nsns_id[10] = "J0737-3039A";
  nsns[10][0] = 1.3381; nsns[10][1] = 0.0007; nsns[10][2] = 0.0007;

  nsns_id[11] = "J0737-3039B";
  nsns[11][0] = 1.2489; nsns[11][1] = 0.0007; nsns[11][2] = 0.0007;

  nsns_id[12] = "J1756-2251 pulsar";
  nsns[12][0] = 1.312; nsns[12][1] = 0.017; nsns[12][2] = 0.017;

  nsns_id[13] = "J1756-2251 comp.";
  nsns[13][0] = 1.258; nsns[13][1] = 0.017; nsns[13][2] = 0.017;

  nsns_id[14] = "J1807-2500B pulsar";
  nsns[14][0] = 1.3655; nsns[14][1] = 0.0021; nsns[14][2] = 0.0021;

  nsns_id[15] = "J1807-2500B comp.";
  nsns[15][0] = 1.2064; nsns[15][1] = 0.002; nsns[15][2] = 0.002;

  // [source: Kiziltan et al. 2013]
  nsns_id[16] = "J1518+4904 pulsar";
  nsns[16][0] = 1.56; nsns[16][1] = 0.13; nsns[16][2] = 0.44;

  nsns_id[17] = "J1518+4904 comp.";
  nsns[17][0] = 1.05; nsns[17][1] = 0.45; nsns[17][2] = 0.11;

  nsns_id[18] = "J1756-2251 pulsar";
  nsns[18][0] = 1.40; nsns[18][1] = 0.02; nsns[18][2] = 0.03;

  nsns_id[19] = "J1756-2251 comp.";
  nsns[19][0] = 1.18; nsns[19][1] = 0.03; nsns[19][2] = 0.02;

  nsns_id[20] = "J1811-1736 pulsar";
  nsns[20][0] = 1.56; nsns[20][1] = 0.24; nsns[20][2] = 0.45;

  nsns_id[21] = "J1811-1736 comp.";
  nsns[21][0] = 1.12; nsns[21][1] = 0.47; nsns[21][2] = 0.13;

  nsns_id[22] = "J1829+2456 pulsar";
  nsns[22][0] = 1.20; nsns[22][1] = 0.12; nsns[22][2] = 0.46;

  nsns_id[23] = "J1829+2456 comp.";
  nsns[23][0] = 1.40; nsns[23][1] = 0.46; nsns[23][2] = 0.12;

  nsns_id[24] = "J1906+0746 pulsar";
  nsns[24][0] = 1.248; nsns[24][1] = 0.018; nsns[24][2] = 0.018;

  nsns_id[25] = "J1906+0746 comp.";
  nsns[25][0] = 1.365; nsns[25][1] = 0.018; nsns[25][2] = 0.018;

  /* Data from NS-WD (some with asymmetric 68% limits)
     [source: Alsing et al. 2018, Kiziltan et al. 2013] */
  nswd_id[0] = "J2045+3633";
  nswd[0][0] = 1.33; nswd[0][1] = 0.3; nswd[0][2] = 0.3;

  nswd_id[1] = "J2053+4650";
  nswd[1][0] = 1.40; nswd[1][1] = 0.21; nswd[1][2] = 0.21;

  nswd_id[2] = "J1713+0747";
  nswd[2][0] = 1.35; nswd[2][1] = 0.07; nswd[2][2] = 0.07;

  nswd_id[3] = "B1855+09";
  nswd[3][0] = 1.37; nswd[3][1] = 0.13; nswd[3][2] = 0.13;

  nswd_id[4] = "J0751+1807";
  nswd[4][0] = 1.72; nswd[4][1] = 0.07; nswd[4][2] = 0.07;

  nswd_id[5] = "J1141-6545";
  nswd[5][0] = 1.27; nswd[5][1] = 0.01; nswd[5][2] = 0.01;

  nswd_id[6] = "J1738+0333";
  nswd[6][0] = 1.47; nswd[6][1] = 0.07; nswd[6][2] = 0.07;

  nswd_id[7] = "J1614-2230";
  nswd[7][0] = 1.908; nswd[7][1] = 0.016; nswd[7][2] = 0.016;

  nswd_id[8] = "J0348+0432";
  nswd[8][0] = 2.01; nswd[8][1] = 0.04; nswd[8][2] = 0.04;

  nswd_id[9] = "J2222-0137";
  nswd[9][0] = 1.76; nswd[9][1] = 0.06; nswd[9][2] = 0.06;

  nswd_id[10] = "J2234+0611";
  nswd[10][0] = 1.393; nswd[10][1] = 0.013; nswd[10][2] = 0.013;

  nswd_id[11] = "J1949+3106";
  nswd[11][0] = 1.47; nswd[11][1] = 0.43; nswd[11][2] = 0.43;

  nswd_id[12] = "J1012+5307";
  nswd[12][0] = 1.83; nswd[12][1] = 0.11; nswd[12][2] = 0.11;

  nswd_id[13] = "J0437-4715";
  nswd[13][0] = 1.44; nswd[13][1] = 0.07; nswd[13][2] = 0.07;

  nswd_id[14] = "J1909-3744";
  nswd[14][0] = 1.48; nswd[14][1] = 0.03; nswd[14][2] = 0.03;

  nswd_id[15] = "J1802-2124";
  nswd[15][0] = 1.24; nswd[15][1] = 0.11; nswd[15][2] = 0.11;

  nswd_id[16] = "J1911-5958A";
  nswd[16][0] = 1.34; nswd[16][1] = 0.08; nswd[16][2] = 0.08;

  nswd_id[17] = "J2043+1711";
  nswd[17][0] = 1.38; nswd[17][1] = 0.13; nswd[17][2] = 0.13;

  nswd_id[18] = "J0337+1715";
  nswd[18][0] = 1.4378; nswd[18][1] = 0.0013; nswd[18][2] = 0.0013;

  nswd_id[19] = "J1946+3417";
  nswd[19][0] = 1.828; nswd[19][1] = 0.022; nswd[19][2] = 0.022;

  nswd_id[20] = "J1918-0642";
  nswd[20][0] = 1.29; nswd[20][1] = 0.1; nswd[20][2] = 0.1;

  nswd_id[21] = "J1600-3053";
  nswd[21][0] = 2.3; nswd[21][1] = 0.7; nswd[21][2] = 0.7;

  // [source: Kiziltan et al. 2013]
  nswd_id[22] = "J0437-4715";
  nswd[22][0] = 1.76; nswd[22][1] = 0.20; nswd[22][2] = 0.20;

  nswd_id[23] = "J0621+1002";
  nswd[23][0] = 1.70; nswd[23][1] = 0.10; nswd[23][2] = 0.17;

  nswd_id[24] = "J0751+1807";
  nswd[24][0] = 1.26; nswd[24][1] = 0.14; nswd[24][2] = 0.14;

  nswd_id[25] = "J1012+5307";
  nswd[25][0] = 1.64; nswd[25][1] = 0.22; nswd[25][2] = 0.22;

  nswd_id[26] = "J1614-2230";
  nswd[26][0] = 1.97; nswd[26][1] = 0.04; nswd[26][2] = 0.04;

  nswd_id[27] = "J1713+0747";
  nswd[27][0] = 1.53; nswd[27][1] = 0.08; nswd[27][2] = 0.06;

  nswd_id[28] = "B1855+09";
  nswd[28][0] = 1.57; nswd[28][1] = 0.12; nswd[28][2] = 0.11;

  nswd_id[29] = "J1909-3744";
  nswd[29][0] = 1.438; nswd[29][1] = 0.024; nswd[29][2] = 0.024;

  nswd_id[30] = "B2303+46";
  nswd[30][0] = 1.38; nswd[30][1] = 0.06; nswd[30][2] = 0.10;

  nswd_id[31] = "J0024-7204H";
  nswd[31][0] = 1.48; nswd[31][1] = 0.03; nswd[31][2] = 0.06;

  nswd_id[32] = "J0514-4002A";
  nswd[32][0] = 1.49; nswd[32][1] = 0.04; nswd[32][2] = 0.27;

  nswd_id[33] = "B1516+02B";
  nswd[33][0] = 2.10; nswd[33][1] = 0.19; nswd[33][2] = 0.19;

  nswd_id[34] = "J1748-2446I";
  nswd[34][0] = 1.91; nswd[34][1] = 0.02; nswd[34][2] = 0.10;

  nswd_id[35] = "J1748-2446J";
  nswd[35][0] = 1.79; nswd[35][1] = 0.02; nswd[35][2] = 0.10;

  nswd_id[36] = "B1802-07";
  nswd[36][0] = 1.26; nswd[36][1] = 0.08; nswd[36][2] = 0.17;

  nswd_id[37] = "B1911-5958A";
  nswd[37][0] = 1.40; nswd[37][1] = 0.16; nswd[37][2] = 0.10;

  // [source: Cromartie et al. 2020]
  nswd_id[38] = "J0740+6620";
  nswd[38][0] = 2.14; nswd[38][1] = 0.20; nswd[38][2] = 0.18;

  /* Data from x-ray/optical (all with symmetric 68% limits)
     [source: Alsing et al. 2018] */
  nsms[0][0]  = 1.96;  nsms[0][1]  = 0.19; nsms_id[0]  = "4U1700-377";
  nsms[1][0]  = 1.71;  nsms[1][1]  = 0.21; nsms_id[1]  = "Cyg X-2";
  nsms[2][0]  = 1.21;  nsms[2][1]  = 0.12; nsms_id[2]  = "SMC X-1";
  nsms[3][0]  = 1.57;  nsms[3][1]  = 0.16; nsms_id[3]  = "Cen X-3";
  nsms[4][0]  = 1.53;  nsms[4][1]  = 0.42; nsms_id[4]  = "XTE J2123-058";
  nsms[5][0]  = 1.96;  nsms[5][1]  = 0.36; nsms_id[5]  = "4U 1822-371";
  nsms[6][0]  = 1.74;  nsms[6][1]  = 0.3;  nsms_id[6]  = "OAO 1657-415";
  nsms[7][0]  = 2.0;   nsms[7][1]  = 0.4;  nsms_id[7]  = "J013236.7+303228";
  nsms[8][0]  = 2.12;  nsms[8][1]  = 0.16; nsms_id[8]  = "Vela X-1";
  nsms[9][0]  = 1.02;  nsms[9][1]  = 0.17; nsms_id[9]  = "4U1538-522";
  nsms[10][0] = 1.57;  nsms[10][1] = 0.11; nsms_id[10] = "LMC X-4";
  nsms[11][0] = 1.073; nsms[11][1] = 0.36; nsms_id[11] = "Her X-1";
  nsms[12][0] = 1.44;  nsms[12][1] = 0.1;  nsms_id[12] = "2S 0921-630";
  nsms[13][0] = 1.91;  nsms[13][1] = 0.45; nsms_id[13] = "EXO 1722-363";
  nsms[14][0] = 1.57;  nsms[14][1] = 0.25; nsms_id[14] = "SAX J1802.7-2017";
  nsms[15][0] = 1.41;  nsms[15][1] = 0.24; nsms_id[15] = "XTE J1855-026";

  /* Data from NS-MS (all with symmetric 68% limits)
     [source: Alsing et al. 2018]
  nsms_rd[0][0] = 1.58;   nsms_rd[0][1] = 0.34;  nsms_id[0] = "J0045-7319";
  nsms_rd[1][0] = 1.71;   nsms_rd[1][1] = 0.16;  nsms_id[1] = "J1023+0038";
  nsms_rd[2][0] = 1.666;  nsms_rd[2][1] = 0.01;  nsms_id[2] = "J1903+0327"; */
}
