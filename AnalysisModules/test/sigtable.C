#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"

#include "optiUtils.C"

using namespace std;

const double optmwrcuts[] = {
  //760, 800, 880,1000,1080,1160,1240          // significance-optimized 1000-1600
  //520,600,640,720,800,880,960,1040,1120,1200 // added low masspoints 700-900
  560,640,720,760,800,840,840,880,880,880      // switched to CLfactor optimization
};

#define sqr(_x) ((_x)*(_x))

const double othuncsqr = sqr(.01)+sqr(.005)+sqr(.02)+sqr(.04); // MES,trigeff,muid,pdf

//======================================================================

double calcErrFromJES(TH1 *hsignm, TH1 *hsighi, TH1 *hsiglo)
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
}                                                      // calcErrFromJES

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
  printf("MWR   MNu  Cut     Ssig   Serr  JES(%%)  MC(%%) Tot(%%)\n");

  for (int i=0; i<10; i++) {
    double mwr = 700. + i*100.;
    cout<<"=============================================================="<<endl;

    double optmwrcut = optmwrcuts[i];

    for (double mnu=100; mnu<mwr; mnu+=100) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );    string signm(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_hi",mwr,mnu ); string sighi(s);
      sprintf( s, "WR%.0f_nuRmu%.0f_lo",mwr,mnu ); string siglo(s);

      assert (m_wth1[signm]);
      assert (m_wth1[sighi]);
      assert (m_wth1[siglo]);

      TH1   *hsignm   = m_wth1[signm]->histo();
      TH1   *hsighi   = m_wth1[sighi]->histo();
      TH1   *hsiglo   = m_wth1[siglo]->histo();

      double jesunc = calcErrFromJES(hsignm,hsighi,hsiglo);
      double totsig, mcstatunc;
      calcSig(hsignm,optmwrcut, totsig,mcstatunc);

      double totunc = sqrt(sqr(jesunc)+sqr(mcstatunc)+othuncsqr);
      double toterr = totsig*totunc;

      printf("%4.0f %4.0f %4.0f | ",mwr,mnu,optmwrcut);
      printf("%6.2f %6.2f %5.1f%% %5.1f%% %5.1f%%\n",
	     totsig,toterr,jesunc*100,mcstatunc*100,totunc*100);
    } // mnu loop
  } // mwr loop
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
