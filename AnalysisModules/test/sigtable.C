#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

#include "hnuUtils.C"

using namespace std;

const double optmwrcuts[] = {
  //760, 800, 880,1000,1080,1160,1240          // significance-optimized 1000-1600
  //520,600,640,720,800,880,960,1040,1120,1200 // added low masspoints 700-900
  560,640,720,760,800,840,840,880,880,880      // switched to CLfactor optimization
};

#define sqr(_x) ((_x)*(_x))

const double othuncsqr = 
  sqr(.005)+ // trigeff
  sqr(.02)+  // muid
  sqr(.03)+  // isr/fsr
  sqr(.04);  // pdf

//======================================================================

double calcErrFromJESorMES(TH1 *hsignm, TH1 *hsighi, TH1 *hsiglo)
{
  const double optmwrcut = 520;

  int    lobin    = hsignm->FindBin(optmwrcut);
  int    hibin    = hsignm->FindBin(1+hsignm->GetNbinsX());
      
  double sumsignm = hsignm->Integral(lobin,hibin);
  double sumsighi = hsighi->Integral(lobin,hibin);
  double sumsiglo = hsiglo->Integral(lobin,hibin);

  double sigfrcdevhi = (sumsighi-sumsignm)/sumsignm;
  double sigfrcdevlo = (sumsiglo-sumsignm)/sumsignm;

  return max(fabs(sigfrcdevhi),fabs(sigfrcdevlo));
}                                                 // calcErrFromJESorMES

//======================================================================

void calcSig(TH1 *sighist,double cutgev,
	     double& sumsig,double& mcstaterr)
{
  int lobin   = sighist->FindBin(cutgev);
  int hibin   = sighist->FindBin(1+sighist->GetNbinsX());
  sumsig      = sighist->Integral(lobin,hibin);
  mcstaterr   = 1./sqrt(sighist->GetEntries());
}                                                             // calcSig

//======================================================================

void
printTable(map<string,wTH1*>&  m_wth1)
{
  FILE *systerrfp = fopen("sigsysterr.h","w");

  printf("MWR   MNu  Cut     Ssig   Serr  JES(%%) MES(%%) MC(%%) Tot(%%)\n");

  fprintf(systerrfp,"float fracerrPerMasspt[][16] = {\n");
  fprintf(systerrfp,"  /* MNu = ");
  for (int mnu=100;mnu<=1600;mnu+=100)
    fprintf (systerrfp,"\t%d",mnu);
  fprintf(systerrfp,"*/\n");
  for (int i=0; i<10; i++) {
    double mwr = 700. + i*100.;
    cout<<"=============================================================="<<endl;

    fprintf(systerrfp,"/*MWR=%4.0f*/{ ", mwr);

    double optmwrcut = optmwrcuts[i];

    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",      mwr,mnu ); string signom(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_jeshi",mwr,mnu ); string sigjhi(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_jeslo",mwr,mnu ); string sigjlo(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_meshi",mwr,mnu ); string sigmhi(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_meslo",mwr,mnu ); string sigmlo(s);

      assert (m_wth1[signom]);
      assert (m_wth1[sigjhi]);
      assert (m_wth1[sigjlo]);
      assert (m_wth1[sigmhi]);
      assert (m_wth1[sigmlo]);

      TH1   *hsignom  = m_wth1[signom]->histo();
      TH1   *hsigjhi  = m_wth1[sigjhi]->histo();
      TH1   *hsigjlo  = m_wth1[sigjlo]->histo();
      TH1   *hsigmhi  = m_wth1[sigmhi]->histo();
      TH1   *hsigmlo  = m_wth1[sigmlo]->histo();

      double jesunc = calcErrFromJESorMES(hsignom,hsigjhi,hsigjlo);
      double mesunc = calcErrFromJESorMES(hsignom,hsigmhi,hsigmlo);
      double totsig, mcstatunc;
      calcSig(hsignom,optmwrcut,totsig,mcstatunc);

      double totunc = sqrt(sqr(jesunc)+sqr(mesunc)+sqr(mcstatunc)+othuncsqr);
      double toterr = totsig*totunc;

      printf("%4.0f %4.0f %4.0f | ",mwr,mnu,optmwrcut);
      printf("%6.2f %6.2f %5.1f%% %5.1f%% %5.1f%% %5.1f%%\n",
	     totsig,toterr,jesunc*100,mesunc*100,mcstatunc*100,totunc*100);

      fprintf(systerrfp,"%.3g, ", totunc);

    } // mnu loop

    for (double mnu=mwr; mnu<=1600; mnu+=100)
      fprintf(systerrfp,"\t-1,");
    fprintf(systerrfp,"},\n");

  } // mwr loop

  fprintf(systerrfp,"};\n");
  fclose(systerrfp);
}                                                          // printTable

//======================================================================

void sigtable()
{
  map<string,wTH1 *> m_wth1;

  // preload map so I don't have to do finds and inserts
  //
  for (double mwr=700; mwr<=1600; mwr+=1000) {
    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);
      m_wth1[name] = NULL;
    }
  }

  getHistosFromRE("jesmesfinals.root",".*",m_wth1);

  printTable(m_wth1);
}
