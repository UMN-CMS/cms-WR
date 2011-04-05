#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "tdrstyle.C"
#include "TVectorD.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"

#include "optiUtils.C"
#include "CLfast.C"

using namespace std;

//======================================================================
// Internal types

struct optVars_t {
  optVars_t():
    cutgev(0),ttint(0),zjint(0),wjint(0),vvint(0),twint(0),qcdint(0),
    sumback(0),sumsig(0),sigeff(0),signif(0),clfactor(9e99) {}
  double cutgev,ttint,zjint,wjint,vvint,twint,qcdint,sumback,sumsig,sigeff,signif,clfactor;
};

struct mWRscanPerMasspt_t {
  double mwr,mnu;       // masspoint
  int minclfactbin;     // index into v_vars for CLs-optimized variables
  optVars_t v_vars[27]; // (1600-520)/40 = 27
  optVars_t minclfact;  // variables yielding minimized CLs factor
  optVars_t maxsignif;  // variables yielding maximized significance (s/sqrt(s+b))
};

struct mLLscanPerMasspt_t {
  double mwr,mnu;
  int minclfactbin;
  optVars_t v_vars[30]; // (200-120)/4 + (600-200)/40 = 30
  optVars_t minclfact;
  optVars_t maxsignif;
};

struct optFits_t {
  optFits_t():tt(0),zj(0),wj(0),vv(0),tw(0),qcd(0) {}
  TF1 *tt,*zj,*wj,*vv,*tw,*qcd;
};

//======================================================================
// global/static variables and constants

//const double zjet_scale = 1.57;
//const double zjet_scale = 1.46;
const double zjet_scale = 1.15;
const double qcd_mllcutsurvivalweightfactor = 0.1;

// Bryan uses p0*exp(p1*x)
const double bryansqcdp0 = 0.1149;     // constant
const double bryansqcdp1 = -0.008046;  // slope

// I/Root uses expo(p0 + p1*x)
//
const double myqcdp0 = log(qcd_mllcutsurvivalweightfactor*bryansqcdp0);

const double wjp0 = -1.729 + log(0.022/0.14);
const double wjp1 = -0.006617; // from chi2 fit of mWR histo after cut4_LLJJpt

const double binwidthGeV=40.0;
const double zoombinwidthGeV=4.0;

//const int rebinMWRval = 5;
const int rebinMWRval = 2;
const int rebinMLLval = 1;

const int nmassptsMWR = 105;
const int nmassptsMLL = 84;

//======================================================================

void fitndraw(wTH1 *wth1,const char *fitoption)
{
  gStyle->SetOptStat(110010);
  gStyle->SetOptFit(1);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.03);

  gPad->SetLogy(1);
  gPad->SetRightMargin(0.05);
  wth1->SetStats(true,0,.6,.7,.99,.99);
  wth1->histo()->GetXaxis()->SetNdivisions(505);
  wth1->Draw("HIST");
  gPad->Update();
  wth1->DrawFits("same",fitoption);
  gPad->Update();
  wth1->DrawStats();
  gPad->Update();
}                                                            // fitndraw

//======================================================================
inline
double pow10floor(double x) {
  return (pow(10.,std::floor(log10(x))));
}
double pow10ceil(double x) {
  return (pow(10.,std::ceil(log10(x))));
}

void drawGraph(double mwr, double mnu, const char *fitfn,
	       double fitrangemin,double fitrangemax,
	       const char *xtitle,
	       const TVectorD& vxc, const TVectorD& vyc,
	       const TVectorD& vxb, const TVectorD& vyb,
	       const TVectorD& vxs, const TVectorD& vys,
	       const TVectorD& vxf, const TVectorD& vyf,
	       double& optgev, double& optgevlo, double& optgevhi)
{
  char s[80];
  sprintf( s, "WR%.0f",mwr );
  string name(s);

  TGraph *grs = new TGraph(vxs,vys);  grs->SetMarkerStyle(24);
  TGraph *grb = new TGraph(vxb,vyb);  grb->SetMarkerStyle(25);

  TVectorD newvyf(vyf);
  newvyf *= 0.2; // significance/5

  double sbmin=9e99, sigclfmin=9e99;
  double sbmax=9e-99,sigclfmax=9e-99;;
  for (int i=0; i<vxb.GetNoElements(); i++) {
    sbmin =min(sbmin,min(vyb[i],vys[i]));
    sbmax =max(sbmax,max(vyb[i],vys[i]));
    sigclfmin=min(sigclfmin,min(vyc[i],newvyf[i]));
    sigclfmax=max(sigclfmax,max(vyc[i],newvyf[i]));
  }
  //grs->GetYaxis()->SetRangeUser( pow10floor(sbmin),
                                 //pow10ceil (sbmax) );

  grs->GetYaxis()->SetRangeUser( 0, 1.3*ceil (sbmax) );

  //grs->GetXaxis()->SetRangeUser(100,250);
  //grs->GetYaxis()->SetRangeUser( 0.1, 10 );

  //TLegend *leg = new TLegend(.2,.45,.6,.6);
  TLegend *leg = new TLegend(.3,.8,.75,.99);

  leg->SetTextSize(.04);
  leg->SetTextFont(42);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetLineWidth(1);
  leg->SetNColumns(1);

  string title=name+"; "+string(xtitle);
  grs->SetTitle(title.c_str());
  grs->GetXaxis()->SetNdivisions(505);
  grs->GetYaxis()->SetTitle("# Events");
  grs->Draw("AP"); grb->Draw("P same");
  gPad->Update();

  // Draw alternative y axis on the right for clfactor
  float rightmin = 0.0; // 0.8*sigclfmin;
  float rightmax = 1.0; // 1.3*sigclfmax;
  float scale    = (gPad->GetUymax()-gPad->GetUymin())/(rightmax-rightmin);

#if 0
  printf ("%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
	  sigclfmin,sigclfmax,rightmin,rightmax,
	  gPad->GetUymin(),gPad->GetUymax(),scale);
#endif
  TVectorD newvyc(vyc);
  newvyc -= rightmin;
  newvyc *=  scale;
  newvyc += gPad->GetUymin();

  TGraph *grc = new TGraph(vxc,newvyc);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  grc->SetMarkerStyle(29);
  grc->SetMarkerColor(kRed);
  grc->SetMarkerSize(1.5);
  //grc->Draw("P");

  newvyf -= rightmin;
  newvyf *=  scale;
  newvyf += gPad->GetUymin();

  TGraph *grf = new TGraph(vxf,newvyf);
  grf->SetMarkerStyle(34);
  grf->SetMarkerColor(kBlue);
  //grf->SetMarkerSize(1.5);
  grf->Draw("P");

  // Fit CL factor to input fit function, fitrange,
  //   find the 10% error bars on the optimal GeV cut
  //
  TF1 *f1 = new TF1((name+"_clfit").c_str(),fitfn,fitrangemin,fitrangemax);
   
  grc->Fit(f1,"RQ");
  optgev         = f1->GetMinimumX(fitrangemin,fitrangemax);
  double optval  = f1->Eval(optgev);
  optgevhi       = max(optgev+40.,f1->GetX(1.10*optval,optgev,fitrangemax));
  optgevlo       = min(optgev-40.,f1->GetX(1.10*optval,fitrangemin,optgev));

  TLine *cuthiline = new TLine(); cuthiline->SetLineColor(kRed);
  TLine *cutloline = new TLine(); cutloline->SetLineColor(kRed);

  //cutloline->DrawLine(optgevlo,0,optgevlo,optval*1.10);
  //cuthiline->DrawLine(optgevhi,0,optgevhi,optval*1.10);
  gPad->Update();

  //leg->AddEntry(grc,"CL factor","P");
  leg->AddEntry(grf,"Significance/5","P");
  leg->AddEntry(grb,"#Sigma Background","P");

  sprintf(s,"(W_{R},N_{R})=(%.0f,%.0f)",mwr,mnu);
  leg->AddEntry(grs,s,"P");
  gPad->Update();
  gPad->SetLogy(0);
  leg->Draw();
  gPad->SetRightMargin(0.15);
  gPad->SetTicks(1,0); // suppress tickmarks on the right from the main graph
  gPad->Update();

  //draw an axis on the right side as the last thing.
  TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
			    gPad->GetUxmax(), gPad->GetUymax(),
			    rightmin,rightmax,510,"+L");
  axis->SetTitleFont(42);
  axis->SetTitleSize(.06);
  axis->SetTitleOffset(1.2);
  axis->SetLabelFont(42);
  axis->SetLabelSize(.05);
  axis->SetLabelOffset(0.007);
  axis->SetLineColor(kBlue); //kRed);
  axis->SetTitle("Signif./5"); // CL factor, Signif./5");
  axis->Draw();
  gPad->Update();
}                                                           // drawGraph

//======================================================================

void drawMWRcutGeVgraph(TVectorD& vx, const TVectorD& vy,
			const TVectorD& veyl, const TVectorD& veyh)
{
  TCanvas *c1 = new TCanvas("c1","c1",1100,300);

  TVectorD zerov(vx.GetNoElements()); // for x errors
  zerov.Zero();

  vx += 1;

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(vx,vy,zerov,zerov,veyl,veyh);

  gr->GetYaxis()->SetTitle("Opt. MWR Cut (GeV)");
  gr->GetYaxis()->SetTitleOffset(0.6);
  gr->SetMarkerStyle(7);
  gr->Draw("AP");
  gPad->SetLeftMargin(0.07);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.2);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  c1->Update();
  TAxis *xax = gr->GetXaxis();
  xax->Set(107,0,106);
  int ibin=2;
  for( double mwr=700; mwr<=1600; mwr+=100 ) {
    for( double mnu=100; mnu<mwr; mnu+=100 ) {
      char s[80];
      sprintf( s, "%.0f,%.0f",mwr,mnu );
      xax->SetBinLabel(ibin++,s);
    }
  }

  c1->Update();
  
  c1->SaveAs("optMWRcutvalPerMasspt.png");

}                                                  // drawMWRcutGeVgraph

//======================================================================

void drawMLLcutGeVgraph(TVectorD& vx, const TVectorD& vy,
			const TVectorD& veyl, const TVectorD& veyh)
{
  TCanvas *c1 = new TCanvas("c1","c1",1200,300);

  TVectorD zerov(vx.GetNoElements()); // for x errors
  zerov.Zero();

  vx += 1;

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(vx,vy,zerov,zerov,veyl,veyh);

  //gr->GetXaxis()->SetRangeUser(0,42);
  gr->GetYaxis()->SetRangeUser(100,300);
  gr->GetYaxis()->SetTitle("Opt. MLL Cut (GeV)");
  gr->GetYaxis()->SetTitleOffset(0.6);
  gr->SetMarkerStyle(7);
  gr->Draw("AP");
  gPad->SetLeftMargin(0.07);
  gPad->SetRightMargin(0.02);
  gPad->SetBottomMargin(0.2);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  c1->Update();
#if 1
  TAxis *xax = gr->GetXaxis();
  cout << vx.GetNrows() << endl;
  cout << vy.GetNrows() << endl;
  cout << xax->GetNbins() << endl;
  xax->Set(83,0,82);
  int ibin=2;
  for( double mwr=1000; mwr<=1600; mwr+=100 ) {
    for( double mnu=100; mnu<mwr; mnu+=100, ibin+=1 /* 2 */ ) {
      char s[80];
      sprintf( s, "%.0f,%.0f",mwr,mnu );
      xax->SetBinLabel(ibin,s);
    }
  }
#endif
  c1->Update();
  
  c1->SaveAs("optMLLcutvalPerMasspt.png");
}                                                  // drawMLLcutGeVgraph

//======================================================================

void loadMWRFits(map<string,wTH1*>&  m_wth1,
		 optFits_t& fits,
		 bool write2file=false)
{
  cout << "===> loadMWRFits <===" << endl;

  wTH1 *ttbar = m_wth1["ttjets_m4"];fits.tt = new TF1("ttfit","expo",650,2000); ttbar->loadFitFunction(fits.tt);
  wTH1 *zjets = m_wth1["zjets_m4"]; fits.zj = new TF1("zjfit","expo",650,2000); zjets->loadFitFunction(fits.zj);
  wTH1 *vv    = m_wth1["vv_m4"];    fits.vv = new TF1("vvfit","expo",650,2000); vv   ->loadFitFunction(fits.vv);
  wTH1 *tw    = m_wth1["tw_m4"];    fits.tw = new TF1("twfit","expo",650,2000); tw   ->loadFitFunction(fits.tw);
#if 0 // see below
  wTH1 *wjets = m_wth1["wjets_m4"]; fits.wj = new TF1("wjfit","expo",560,2000); wjets->loadFitFunction(fits.wj);
#endif

  ttbar->histo()->Rebin(rebinMWRval);
  vv   ->histo()->Rebin(rebinMWRval);
  tw   ->histo()->Rebin(rebinMWRval);

  // Scale zjets according to Z-peak fit:
  zjets->histo()->Scale(zjet_scale);
  zjets->histo()->Rebin(rebinMWRval);

  // QCD estimated from data, see notes for fit parameters above
  //
  fits.qcd = new TF1("qcdfit","expo",650,2000);
  fits.qcd->SetParameter(0,myqcdp0);
  fits.qcd->SetParameter(1,bryansqcdp1);

  // Wjets is so bad on statistics, even the first cut leaves
  //  barely anything to fit. So we take the fit parameters for
  //  that cut level and then scale it by the additional whack
  //  that the final cuts produce, similar to QCD.
  //
  fits.wj = new TF1("wjfit","expo",650,2000);
  fits.wj->SetParameter(0,wjp0);
  fits.wj->SetParameter(1,wjp1);

  TCanvas *c1=new TCanvas("bgfits","bgfits",700,700);
  c1->Divide(2,2);

  c1->cd(1); fitndraw(zjets,"LL"); // "LL" actually not done correctly by ROOT
  c1->cd(2); fitndraw(ttbar,"LL"); // (the fit parameter errors are all wrong)
  c1->cd(3); fitndraw(vv,"LL");    // but the fits themselves look good for the most part
  c1->cd(4); fitndraw(tw,"LL");    // and so do the fit parameter values
  //c1->cd(5); fitndraw(wjets,"");   // only two points in the distribution, do chi2

  c1->SaveAs("optMWRscanBGfits.png");

  if (write2file) {
    FILE *fp = fopen("optMWRfitpars.h","w");

    const char *fmt = "const double %3sfitpars[] = {%10.4f, %10.4e, %8.4g, %3d, %6.1f, %6.1f};\n";
    fprintf(fp,"// arrays contain:                fitp0,      fitp1,   fitchi2, fitndf, fitxmin, fitxmax\n");
    fprintf(fp,fmt,"tt",
	    fits.tt->GetParameter(0),fits.tt->GetParameter(1),
	    fits.tt->GetChisquare(),fits.tt->GetNDF(),
	    fits.tt->GetXmin(),fits.tt->GetXmax());
    fprintf(fp,fmt,"zj",
	    fits.zj->GetParameter(0),fits.zj->GetParameter(1),
	    fits.zj->GetChisquare(),fits.zj->GetNDF(),
	    fits.zj->GetXmin(),fits.zj->GetXmax());
    fprintf(fp,fmt,"vv",
	    fits.vv->GetParameter(0),fits.vv->GetParameter(1),
	    fits.vv->GetChisquare(),fits.vv->GetNDF(),
	    fits.vv->GetXmin(),fits.vv->GetXmax());
    fprintf(fp,fmt,"tw",
	    fits.tw->GetParameter(0),fits.tw->GetParameter(1),
	    fits.tw->GetChisquare(),fits.tw->GetNDF(),
	    fits.tw->GetXmin(),fits.tw->GetXmax());
    fprintf(fp,fmt,"wj",
	    fits.wj->GetParameter(0),fits.wj->GetParameter(1),
	    fits.wj->GetChisquare(),fits.wj->GetNDF(),
	    fits.wj->GetXmin(),fits.wj->GetXmax());
    fprintf(fp,fmt,"qcd",
	    fits.qcd->GetParameter(0),fits.qcd->GetParameter(1),
	    fits.qcd->GetChisquare(),fits.qcd->GetNDF(),
	    fits.qcd->GetXmin(),fits.qcd->GetXmax());
    
    fclose(fp);
  }
}                                                         // loadMWRFits

//======================================================================

void loadMLLFits(map<string,wTH1*>&  m_wth1,
		 optFits_t& fits,
		 bool write2file=false)
{
  cout << "===> loadMLLFits <===" << endl;

  wTH1 *ttbar = m_wth1["ttjets_m2"];
  wTH1 *zjets = m_wth1["zjets_m2"];
  wTH1 *wjets = m_wth1["wjets_m2"];
  wTH1 *vv    = m_wth1["vv_m2"];
  wTH1 *tw    = m_wth1["tw_m2"];

  fits.tt = new TF1("ttfit","expo",160,1000);
  fits.zj = new TF1("zjfit","expo",160,2000);
  fits.wj = new TF1("wjfit","expo",120,300); // note!
  fits.vv = new TF1("vvfit","expo",160,1200);
  fits.tw = new TF1("twfit","expo",160,1000);

  ttbar->loadFitFunction(fits.tt);
  zjets->loadFitFunction(fits.zj);
  wjets->loadFitFunction(fits.wj);
  vv   ->loadFitFunction(fits.vv);
  tw   ->loadFitFunction(fits.tw);

  ttbar->histo()->Rebin(rebinMLLval);
  vv   ->histo()->Rebin(rebinMLLval);
  tw   ->histo()->Rebin(rebinMLLval);

  // Scale zjets according to Z-peak fit:
  zjets->histo()->Scale(zjet_scale);
  zjets->histo()->Rebin(rebinMLLval);

  TCanvas *c1=new TCanvas("bgfits","bgfits",700,700);
  c1->Divide(3,2);

  fitndraw(tw,"LL");
  fitndraw(zjets,"LL");
  fitndraw(ttbar,"LL");
  fitndraw(vv,"LL");
  fitndraw(wjets,"");

  printf("%-9s%10s%10s%10s%10s%10s\n", "Sample","p0","p1","chi2/ndof","xmin","xmax");
  printf("%-9s%10.4f%10.4f%10.4g/%d%10.0f%10.0f\n","TTJets",
	 fits.tt->GetParameter(0),fits.tt->GetParameter(1),
	 fits.tt->GetChisquare(),fits.tt->GetNDF(),
	 fits.tt->GetXmin(),fits.tt->GetXmax());
  printf("%-9s%10.4f%10.4f%10.4g/%d%10.0f%10.0f\n","ZJets",
	 fits.zj->GetParameter(0),fits.zj->GetParameter(1),
	 fits.zj->GetChisquare(),fits.zj->GetNDF(),
	 fits.zj->GetXmin(),fits.zj->GetXmax());
  printf("%-9s%10.4f%10.4f%10.4g/%d%10.0f%10.0f\n","VV",
	 fits.vv->GetParameter(0),fits.vv->GetParameter(1),
	 fits.vv->GetChisquare(),fits.vv->GetNDF(),
	 fits.vv->GetXmin(),fits.vv->GetXmax());
  printf("%-9s%10.4f%10.4f%10.4g/%d%10.0f%10.0f\n","TW",
	 fits.tw->GetParameter(0),fits.tw->GetParameter(1),
	 fits.tw->GetChisquare(),fits.tw->GetNDF(),
	 fits.tw->GetXmin(),fits.tw->GetXmax());
  printf("%-9s%10.4f%10.4f%10.4g/%d%10.0f%10.0f\n","WJets",
	 fits.wj->GetParameter(0),fits.wj->GetParameter(1),
	 fits.wj->GetChisquare(),fits.wj->GetNDF(),
	 fits.wj->GetXmin(),fits.wj->GetXmax());

  // reset the ranges for background estimation!
  fits.tt->SetRange(120,1000);
  fits.zj->SetRange(120,1000);
  fits.wj->SetRange(120,1000);
  fits.vv->SetRange(120,1000);
  fits.tw->SetRange(120,1000);

  if (write2file) {
    FILE *fp = fopen("optMLLfitpars.h","w");

    const char *fmt = "const double %3sfitpars[] = {%10.4f, %10.4e, %8.4g, %3d, %6.1f, %6.1f};\n";
    fprintf(fp,"// arrays contain:                fitp0,      fitp1,   fitchi2, fitndf, fitxmin, fitxmax\n");

    fprintf(fp,fmt,"tt",
	    fits.tt->GetParameter(0),fits.tt->GetParameter(1),
	    fits.tt->GetChisquare(),fits.tt->GetNDF(),
	    fits.tt->GetXmin(),fits.tt->GetXmax());
    fprintf(fp,fmt,"zj",
	    fits.zj->GetParameter(0),fits.zj->GetParameter(1),
	    fits.zj->GetChisquare(),fits.zj->GetNDF(),
	    fits.zj->GetXmin(),fits.zj->GetXmax());
    fprintf(fp,fmt,"vv",
	    fits.vv->GetParameter(0),fits.vv->GetParameter(1),
	    fits.vv->GetChisquare(),fits.vv->GetNDF(),
	    fits.vv->GetXmin(),fits.vv->GetXmax());
    fprintf(fp,fmt,"tw",
	    fits.tw->GetParameter(0),fits.tw->GetParameter(1),
	    fits.tw->GetChisquare(),fits.tw->GetNDF(),
	    fits.tw->GetXmin(),fits.tw->GetXmax());
    fprintf(fp,fmt,"wj",
	    fits.wj->GetParameter(0),fits.wj->GetParameter(1),
	    fits.wj->GetChisquare(),fits.wj->GetNDF(),
	    fits.wj->GetXmin(),fits.wj->GetXmax());

    fclose(fp);
  }
}                                                         // loadMLLFits

//======================================================================

void calcBack(const optFits_t& fits,optVars_t& vars,double rebinval)
{
  vars.ttint   = fits.tt ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.zjint   = fits.zj ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.wjint   = fits.wj ->Integral(vars.cutgev,2000)/binwidthGeV;
  vars.vvint   = fits.vv ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  vars.twint   = fits.tw ->Integral(vars.cutgev,2000)/(binwidthGeV*rebinval);
  if (fits.qcd)
    vars.qcdint= fits.qcd->Integral(vars.cutgev,2000)/binwidthGeV;
	
  vars.sumback = vars.ttint+vars.zjint+vars.wjint+vars.vvint+vars.twint+vars.qcdint;

}                                                            // calcBack

//======================================================================

void calcSig(TH1 *sighist,optVars_t& vars)
{
  int lobin     = sighist->FindBin(vars.cutgev);
  int hibin     = sighist->FindBin(1+sighist->GetNbinsX());
  vars.sumsig   = sighist->Integral(lobin,hibin);
  vars.sigeff   = sighist->GetEntries()/10000.;
  vars.signif   = vars.sumsig/sqrt(vars.sumsig+vars.sumback);
  vars.clfactor = CLfast_goalSeek(vars.sumsig,vars.sumback); // this takes a while

}                                                             // calcSig

//======================================================================

inline
void putvars(const char *fmt, double mwr, double mnu, const optVars_t& v)
{
  printf(fmt,mwr,mnu,v.cutgev,v.ttint,v.zjint,v.wjint,v.vvint,v.twint,
	 v.qcdint,v.sumback,v.sumsig,v.signif,v.clfactor);
}                                                            // putvars

//======================================================================

void scanMWRCutValsSave2tree(map<string,wTH1*>&  m_wth1,
			     const optFits_t& fits)
{
  cout << "===> scanMWRCutValsSave2tree <===" << endl;

  mWRscanPerMasspt_t scanpt;

  TFile f("mwrscantree.root","RECREATE");
  TTree tree("mytree","MWRscanTree");
  tree.Branch("scanpt",  &scanpt,  32000, 1);

  // calc all backgrounds up front; the background values
  // within scanpt will be set once and written multiple
  // times with different associated signal values to the tree.
  //
  for( int i=0; i<27; i++ ) {
    scanpt.v_vars[i].cutgev=(double)(520+binwidthGeV*i);
    calcBack(fits,scanpt.v_vars[i],rebinMWRval);
  }

  for( scanpt.mwr=700; scanpt.mwr<=1600; scanpt.mwr+=100 ) {
    for( scanpt.mnu=100; scanpt.mnu<scanpt.mwr; scanpt.mnu+=100 ) {

      // get signal histo for this masspoint from the input map
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",scanpt.mwr,scanpt.mnu );
      string name(s);

      if (!m_wth1[name]) {
	cerr << "Didn't find " << name << endl;
	continue;
      }
      cout<<"Scanning cut values for "<<name<<"..."<<flush;

      TH1 *sigh = m_wth1[name]->histo();

      optVars_t minclfact;
      optVars_t maxsignif;

      for(int i=0; i<27; i++ ) {
	optVars_t& vars = scanpt.v_vars[i];
	if( vars.cutgev >= scanpt.mwr ) break;
	calcSig(sigh,vars);
	if (vars.signif   > maxsignif.signif)    maxsignif = vars;
	if (vars.clfactor < minclfact.clfactor)  { minclfact = vars; scanpt.minclfactbin = i; }
      }

      scanpt.maxsignif = maxsignif;
      scanpt.minclfact = minclfact;
      
      tree.Fill();

      cout<<"done."<<endl;
    } // mnu loop
  } // mwr loop

  tree.Write();
  f.Close();
  cout << "Tree written." << endl;
}                                             // scanMWRCutValsSave2tree

//======================================================================

void scanMWRCutValsReadTreePutResults()
{
  cout << "===> scanMWRCutValsReadTreePutResults <===" << endl;

  TFile f("mwrscantree.root");

  TTree *tree = (TTree *)f.Get("mytree");
  if (!tree) {
    exit(-1);
  }
  mWRscanPerMasspt_t *scanpt = 0;

  tree->SetBranchAddress("scanpt", &scanpt);

  cout<<"=================================================================";
  cout<<"================================================================="<<endl;

  printf("%-9s%9s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
"mWR(GeV)","mNU(GeV)","optcutgev","ttbar","zjets","wjets","VV","tW","QCD","Sum BG","Sum Sig.","signif","clfactor");

  cout<<"=================================================================";
  cout<<"================================================================="<<endl;

  // vectors for asymmetric error graph of optimized cut values per masspt
  TVectorD vxclf(nmassptsMWR),vyclf(nmassptsMWR),veylclf(nmassptsMWR),veyhclf(nmassptsMWR);

  double oldmwr = 0;

  TCanvas *c1;
  for (int i=0; i<nmassptsMWR; i++) {

    tree->GetEntry(i);

    //if (scanpt->mwr >= 1000) break;

    char s[80];

    sprintf( s, "WR%.0f_nuRmu%.0f",scanpt->mwr,scanpt->mnu );
    //cout << "Read "<< string(s) << endl;

    TVectorD vxb(27),vyb(27);
    TVectorD vxc(27),vyc(27);
    TVectorD vxs(27),vys(27);
    TVectorD vxf(27),vyf(27);
    int j;
    for( j=0; j<27; j++ ) {
      optVars_t& vars = scanpt->v_vars[j];
      if( vars.cutgev >= scanpt->mwr ) break;
      vxb[j] = vars.cutgev; vyb[j] = vars.sumback;
      vxs[j] = vars.cutgev; vys[j] = vars.sumsig;
      vxc[j] = vars.cutgev; vyc[j] = vars.clfactor;
      vxf[j] = vars.cutgev; vyf[j] = vars.signif;
    }

    // build optimized cutgev-by-masspoint graph
    //
    vxclf[i]=i;
    vyclf[i]=scanpt->minclfact.cutgev;

#if 0
    // cutgev error bars based on 5% variation on clfactor
    //
    int hibin=scanpt->minclfactbin;
    for (int k=scanpt->minclfactbin+1; k<j; k++) {
      double clfhi=scanpt->v_vars[k].clfactor;
      hibin=k;
      if (clfhi > 1.05*scanpt->minclfact.clfactor) break;
    }
    veyhclf[i]=vxc[hibin]-scanpt->minclfact.cutgev;

    int lobin=scanpt->minclfactbin;
    for (int k=scanpt->minclfactbin-1; k>=0; k--) {
      double clflo=scanpt->v_vars[k].clfactor;
      lobin=k;
      if (clflo > 1.05*scanpt->minclfact.clfactor) break;
    }
    veylclf[i]=scanpt->minclfact.cutgev-vxc[lobin];
#endif
    double optgev, optgevlo,optgevhi;
    int polorder = min(6,2*((int)(scanpt->mwr) - 600)/100);
    sprintf (s,"pol%d",polorder);
    string fitfn(s);

    if (scanpt->mwr != oldmwr) {
      if (oldmwr>1.0) {
	c1->SaveAs((string(c1->GetName())+".png").c_str());
	c1->SaveAs((string(c1->GetName())+".eps").c_str());
	c1->SaveAs((string(c1->GetName())+".C").c_str());
      }
      char s[80];
      sprintf( s, "WR%.0f",scanpt->mwr );
      string name(s);

      int nrows = ((int)scanpt->mwr/100)-1; nrows = (nrows%4) ? (1+nrows/4) : nrows/4;

      c1 = new TCanvas(name.c_str(),name.c_str(),4*300,min(750,(nrows*250)));// 500,450);
      c1->Divide(4,nrows);
      oldmwr = scanpt->mwr;
    }
    c1->cd((int)scanpt->mnu/100);

    drawGraph(scanpt->mwr,scanpt->mnu, 
	      fitfn.c_str(),520,scanpt->mwr,
	      "M(W_{R}) Cut Value (GeV)",
	      vxc.GetSub(0,j-1),vyc.GetSub(0,j-1),
	      vxb.GetSub(0,j-1),vyb.GetSub(0,j-1),
	      vxs.GetSub(0,j-1),vys.GetSub(0,j-1),
	      vxf.GetSub(0,j-1),vyf.GetSub(0,j-1),
	      optgev,optgevlo,optgevhi);

    if ((scanpt->mwr == 1300) &&
	(scanpt->mnu == 700) ) { // save separately as an example
      sprintf( s, "WR%.0f_nuRmu%.0f",scanpt->mwr,scanpt->mnu );
      string name(s);
      TCanvas *c2 = new TCanvas(name.c_str(),name.c_str(),500,450);
      drawGraph(scanpt->mwr,scanpt->mnu,
		fitfn.c_str(),520,scanpt->mwr,
		"M(W_{R}) Cut Value (GeV)",
		vxc.GetSub(0,j-1),vyc.GetSub(0,j-1),
		vxb.GetSub(0,j-1),vyb.GetSub(0,j-1),
		vxs.GetSub(0,j-1),vys.GetSub(0,j-1),
		vxf.GetSub(0,j-1),vyf.GetSub(0,j-1),
	      optgev,optgevlo,optgevhi);
      c2->SaveAs((name+".png").c_str());
    }

    // now doing a fit, use fit results
    vyclf  [i]=optgev;
    veyhclf[i]=optgevhi-optgev;
    veylclf[i]=optgev-optgevlo;

    putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)%7.3f\n",
	    scanpt->mwr,scanpt->mnu,scanpt->maxsignif);
    putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)\n\n",
	    scanpt->mwr,scanpt->mnu,scanpt->minclfact);

#if 0
    //printf("#%5s%6s%7s%10s%10s\n", "mwr","mnu","mllcut", "sumbg","sumsig");
    //for( int i=0; i<27; i++ )
    printf("%6.0f%6.0f%7.0f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
	   scanpt->mwr,scanpt->mnu,scanpt->v_vars[0].cutgev,
	   scanpt->v_vars[0].ttint,scanpt->v_vars[0].zjint,scanpt->v_vars[0].wjint,
	   scanpt->v_vars[0].vvint,scanpt->v_vars[0].twint,scanpt->v_vars[0].qcdint);
#endif
  } // tree entry/masspoint loop

  // save last canvas
  if (oldmwr) {
    c1->SaveAs((string(c1->GetName())+".png").c_str());
    c1->SaveAs((string(c1->GetName())+".eps").c_str());
    c1->SaveAs((string(c1->GetName())+".C").c_str());
  }

  drawMWRcutGeVgraph(vxclf,vyclf,veylclf,veyhclf);

  delete tree;
}                                    // scanMWRCutValsReadTreePutResults

//======================================================================
// almost the same as above, but different cut optimization range
//

void scanMLLCutValsSave2tree(map<string,wTH1*>&  m_wth1,
			     const optFits_t& fits)
{
  mLLscanPerMasspt_t scanpt;

  TFile f("mllscantree.root","RECREATE");
  TTree tree("mytree","MLLscanTree");
  tree.Branch("scanpt",  &scanpt,  32000, 1);


  // same number of points for all mass points,
  // so calc backgrounds up front
  //
  for( int i=0; i<30; i++ ) {
    scanpt.v_vars[i].cutgev=(i < 20) ?  // points below 200, use zoom histos (finer binning)
      (double)(120+zoombinwidthGeV*i) :
      (double)(200+binwidthGeV*(i-20));
    calcBack(fits,scanpt.v_vars[i],rebinMLLval);
  }

  for( scanpt.mwr=1000; scanpt.mwr<=1600; scanpt.mwr+=100 ) {
    for( scanpt.mnu=100; scanpt.mnu<scanpt.mwr; scanpt.mnu+=100 ) {

      // get signal histo for this masspoint from the input map
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",scanpt.mwr,scanpt.mnu );
      string name(s);

      if (!m_wth1[name]) {
	cerr << "Didn't find " << name << endl;
	continue;
      }
      cout<<"Scanning cut values for "<<name<<"..."<<flush;

      TH1 *sigh = m_wth1[name]->histo();

      // get zoomed signal histo
      name += "_zoom";
      if( !m_wth1[name] ) {
	cerr << "Didn't find " << name << endl;
	continue;
      }

      TH1 *sighz = m_wth1[name]->histo();

      optVars_t minclfact;
      optVars_t maxsignif;

      for(int i=0; i<30; i++ ) {
	optVars_t& vars = scanpt.v_vars[i];
	(i<20) ? calcSig(sighz,vars) : calcSig(sigh,vars);
	if (vars.signif   > maxsignif.signif)    maxsignif = vars;
	if (vars.clfactor < minclfact.clfactor)  { minclfact = vars; scanpt.minclfactbin = i; }
      }

      scanpt.maxsignif = maxsignif;
      scanpt.minclfact = minclfact;
      
      tree.Fill();

      cout<<"done."<<endl;
    } // mnu loop
  } // mwr loop

  tree.Write();
  f.Close();
  cout << "Tree written." << endl;
}                                             // scanMLLCutValsSave2tree

//======================================================================

void scanMLLCutValsReadTreePutResults()
{
  TFile f("mllscantree.root");

  TTree *tree = (TTree *)f.Get("mytree");
  if (!tree) {
    exit(-1);
  }
  mLLscanPerMasspt_t *scanpt = 0;

  tree->SetBranchAddress("scanpt", &scanpt);

  cout<<"=================================================================";
  cout<<"================================================================="<<endl;

  printf("%-9s%9s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
 "mWR(GeV)","mNU(GeV)","optcutgev","ttbar","zjets","wjets","VV","tW","QCD","Sum BG","Sum Sig.","signif","clfactor");

  cout<<"=================================================================";
  cout<<"================================================================="<<endl;

  // vectors for asymmetric error graph of optimized cut values per masspt
  TVectorD vxclf(nmassptsMLL),vyclf(nmassptsMLL),veylclf(nmassptsMLL),veyhclf(nmassptsMLL);

  for (int i=0; i<nmassptsMLL; i++) {

    tree->GetEntry(i);

    if (scanpt->mwr >= 1400) break;

    char s[80];

    sprintf( s, "WR%.0f_nuRmu%.0f",scanpt->mwr,scanpt->mnu );
    //cout << "Read "<< string(s) << endl;

    TVectorD vxb(30),vyb(30);
    TVectorD vxc(30),vyc(30);
    TVectorD vxs(30),vys(30);
    TVectorD vxf(30),vyf(30);

    for( int j=0; j<30; j++ ) {
      optVars_t& vars = scanpt->v_vars[j];
      vxb[j] = vars.cutgev; vyb[j] = vars.sumback;
      vxc[j] = vars.cutgev; vyc[j] = vars.clfactor;
      vxs[j] = vars.cutgev; vys[j] = vars.sumsig;
      vxf[j] = vars.cutgev; vyf[j] = vars.signif;
    }

    // build optimized cutgev-by-masspoint graph
    //
    vxclf[i]=i;
    vyclf[i]=scanpt->minclfact.cutgev;

#if 0
    // cutgev error bars based on 5% variation on clfactor
    //
    int hibin=scanpt->minclfactbin;
    for (int k=scanpt->minclfactbin+1; k<30; k++) {
      double clfhi=scanpt->v_vars[k].clfactor;
      hibin=k;
      if (clfhi > 1.05*scanpt->minclfact.clfactor) break;
    }
    veyhclf[i]=vxc[hibin]-scanpt->minclfact.cutgev;

    int lobin=scanpt->minclfactbin;
    for (int k=scanpt->minclfactbin-1; k>=0; k--) {
      double clflo=scanpt->v_vars[k].clfactor;
      lobin=k;
      if (clflo > 1.05*scanpt->minclfact.clfactor) break;
    }
    veylclf[i]=scanpt->minclfact.cutgev-vxc[lobin];
#endif

    double optgev, optgevlo, optgevhi;
    drawGraph(scanpt->mwr,scanpt->mnu,
	      "pol6",120,400,
	      "M_{#mu#mu} Cut Value (GeV)",
	      vxc,vyc,vxb,vyb,vxs,vys,vxf,vyf,
	      optgev, optgevlo, optgevhi);

    if ((scanpt->mwr == 1100) &&
	(scanpt->mnu == 500) ) { // save separately as an example
      sprintf( s, "WR%.0f_nuRmu%.0f",scanpt->mwr,scanpt->mnu );
      string name(s);
      TCanvas *c2 = new TCanvas(name.c_str(),name.c_str(),500,450);
      drawGraph(scanpt->mwr,scanpt->mnu,
		"pol6",120,400,
		"M_{#mu#mu} Cut Value (GeV)",
		vxc,vyc,vxb,vyb,vxs,vys,vxf,vyf,
		optgev, optgevlo, optgevhi);

      c2->SaveAs((name+".png").c_str());
      break;
    }

    // now doing a fit, use fit results
    vyclf  [i]=optgev;
    veyhclf[i]=optgevhi-optgev;
    veylclf[i]=optgev-optgevlo;

    putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)%7.3f\n",
	    scanpt->mwr,scanpt->mnu,scanpt->maxsignif);
    putvars("%-9.0f%9.0f%10.0f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f(*)\n\n",
	    scanpt->mwr,scanpt->mnu,scanpt->minclfact);
#if 0
    printf("#%5s%6s%7s%10s%10s\n", "mwr","mnu","mllcut", "sumbg","sumsig");
    for( int i=0; i<30; i++ )
      printf("%6.0f%6.0f%7.0f%10.4f%10.4f\n",mwr,mnu,vxs[i],vyb[i],vys[i]);
#endif
  } // tree entry/masspoint loop

  drawMLLcutGeVgraph(vxclf,vyclf,veylclf,veyhclf);

  delete tree;
}                                    // scanMLLCutValsReadTreePutResults

//======================================================================

void opticut( int mode )
{
  map<string,wTH1 *> m_wth1;
  optFits_t fits;

  setTDRStyle();
  gROOT->ForceStyle();

  // preload map so I don't have to do finds and inserts
  //
  for( double mwr=700; mwr<=1600; mwr+=1000 ) {
    for( double mnu=100; mnu<mwr; mnu+=100 ) {
      char s[80];
      sprintf( s, "WR%.0f_nuRmu%.0f",mwr,mnu );
      string name(s);
      m_wth1[name] = NULL;

      name += "_zoom";
      m_wth1[name] = NULL;
    }
  }

  string input;
  switch( mode ) {
    // ==========> MWR SCANS <==========
  case 0:
    {
      getHistosFromRE("optimMWR.root",".*",m_wth1);
      bool write2file=true;
      loadMWRFits(m_wth1,fits,write2file);
    }
    break;
  case 1:
    cout << "This'll take a while, are you sure?" << endl;
    cin >> input;
    if (tolower(input[0]) == 'y') {
      getHistosFromRE("optimMWR.root",".*",m_wth1);
      loadMWRFits(m_wth1,fits);
      scanMWRCutValsSave2tree(m_wth1,fits);
    }
    break;
  case 2:
    scanMWRCutValsReadTreePutResults();
    break;
#if 0
  case 3:
    getHistosFromRE("optimMWR.root",".*",m_wth1);
    printDanilaTable(m_wth1,fits);    // For when a set of MWR cuts has been decided on
    break;
#endif
    // ==========> MLL SCANS <==========
  case 3:
    {
      getHistosFromRE("optimMLL.root",".*",m_wth1);
      bool write2file=true;
      loadMLLFits(m_wth1,fits,write2file);
    }
    break;
  case 4:
    cout << "This'll take a while, are you sure?" << endl;
    cin >> input;
    if (tolower(input[0]) == 'y') {
      getHistosFromRE("optimMLL.root",".*",m_wth1);
      loadMLLFits(m_wth1,fits);
      scanMLLCutValsSave2tree(m_wth1,fits);
    }
    break;
  case 5:
    scanMLLCutValsReadTreePutResults();
    break;
  default:
    cerr << "Unknown mode " << mode << endl;
    exit(-1);
  }
}
