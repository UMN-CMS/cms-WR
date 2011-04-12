#ifndef _HNUANALINPUT_H
#define _HNUANALINPUT_H

const char* finalVarHist = "hNu/cut8_mWRmass/mWR";

const double mcnormlumipbinv = 36.1;

/*******************************************************
 * BACKGROUND YIELDS: Number of events (and fractional
 *   errors) estimated after all selections performed
 *   per background type per invpb of luminosity
 *******************************************************/

// XSEC from CMS 194 ± 72 (stat) ± 24 (syst) ± 21 (lumi), plus W->mu nu
//const double ttb_yieldpb_mc = 194.3*(0.11*0.11)*(3818+6)/193317;

//const double tt_yieldpb_mc = 1.56/mcnormlumipbinv, tt_yield_ferror = 0.40; // fractional error
//const double zj_yieldpb_mc = 0.86/mcnormlumipbinv, zj_yield_ferror = 0.10;
//const double qcd_yieldpb_data = 0.0026;
  
//const double wj_yield_mc = lumi*53711.0*1/1.11e7;   // No W+jet after fix of Jet ID

// Apr. 11:
const double tt_yieldpb_mc = 1.03/mcnormlumipbinv,  tt_yield_ferror = 0.21; // fractional error
const double zj_yieldpb_mc = 0.85/mcnormlumipbinv,  zj_yield_ferror = 0.20;
const double vv_yieldpb_mc = 0.037/mcnormlumipbinv, vv_yield_ferror = 0.24;
const double tw_yieldpb_mc = 0.030/mcnormlumipbinv, tw_yield_ferror = 0.24;
const double wj_yieldpb_mc = 0.022/mcnormlumipbinv, wj_yield_ferror = 0.24;

const double qcd_yieldpb_data = 0.0026;

/*****************************************
 * ROOT "expo" fit slope parameters:
 ****************************************/

//const double tt_expfit_p1 = -5.08e-3; // Jan 24
//const double zj_expfit_p1 = -4.35e-3; // Jan 24
//const double qcd_expfit_p1 = -0.003923; // Jan 24

const double tt_expfit_p1 = -5.53e-3; // Apr 11
const double zj_expfit_p1 = -4.63e-3;  // Apr 11
const double wj_expfit_p1 = -6.62e-3;  // Apr 4
const double vv_expfit_p1 = -4.69e-3;  // Apr 4
const double tw_expfit_p1 = -8.17e-3;  // Apr 4
//const double qcd_expfit_p1= -8.05e-3;  // Apr 4


#endif // _HNUANALINPUT_H
