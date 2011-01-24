#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "THStack.h"

using namespace std;

double computeChi2(TH1 *h1, TH1 *h2,int ifirst,int ilast)
{
  double chi2=0.;
  for (int ibin=ifirst; ibin<=ilast; ibin++) {
    double bin1 = h1->GetBinContent(ibin);
    double bin2 = h2->GetBinContent(ibin);
    double bin1e = h1->GetBinError(ibin);
    double bin2e = h2->GetBinError(ibin);
    chi2 += (bin1-bin2)*(bin1-bin2)/((bin1e*bin1e)+(bin2e*bin2e));
    //cout << bin1 << " " << bin2 << "; ";
  }
  //cout << endl;
  return chi2;
}

const double luminvpb = 33.9;
#if 0
// MC->MC matching, reco level (sucks)
//const int rebinx=5;
const int rebinx=0;
const double ymin    = 5e-5;
const double ymax    = 20.;
const double h1sf    = luminvpb * 3740/1647472.0;
const double nevents = 92291.;
const string format  = "#sigma = %4.2f+/-%4.2fpb\n";
const string str1    = "Summer10 Z+Jets Madgraph";
const string str2    = "Z#rightarrow#mu#mu+Jets Madgraph, UMN";
void myChi2Fit(const char *filename1="summer10_zjets_7tev_madgraph_start36_v10-v2_hnu-anal.root",
	       const char *filename2="madgraph_zmumujets_m180_reco385_hnu-anal.root",
	       const char *bckgrndfn="",
	       const char *path     ="hNu/LLJJptcuts/mMuMu",
	       double minMLL=200.0,
	       double maxMLL=600.0,
	       double ctr_xsecpb=0.5,
	       double xsecincpb=.02)
#elif 0
// MC->MC matching, gen level (much better)
const int rebinx=5;
const double ymin    = 1e-2;
const double ymax    = 1e5;
const double h1sf    = 1.0;
const double nevents = 92291.;
const string format  = "#sigma = %4.2f+/-%4.2fpb\n";
const string str1    = "Summer10 Z+Jets Madgraph";
const string str2    = "Z#rightarrow#mu#mu+Jets Madgraph, UMN";
void myChi2Fit(const char *filename1="grendanal_summer10_zjets_7tev_madgraph.root",
	       const char *filename2="grendanal_madgraph_zmumujets_m10_reco385.root",
	       const char *bckgrndfn="",
	       const char *path     ="gr/ptmin/mMuMu",
	       double minMLL=200.0,
	       double maxMLL=600.0,
	       double ctr_xsecpb=55.0,
	       double xsecincpb=0.5)
#elif 0
// Pythia MC->data matching
const int    rebinx  = 0;
const double ymin    = 5e-2;
const double ymax    = 500;
const double h1sf    = 1.0;
const double nevents = 1647472.0;
const string format  = "#sigma = %4.0f+/-%4.0fpb\n";
const string str1    = "/Mu, Run2010A+B, 33.9/pb";
const string str2    = "Summer10 Z+Jets Madgraph";
void myChi2Fit(const char *filename1="bryansNewMuSkim_Run2010AandBcombinedNov19JSON_hnu-anal.root",
	       const char *filename2="summer10_zjets_7tev_madgraph_start36_v10-v2_hnu-anal.root",
	       const char *bckgrndfn="",
	       const char *path     ="hNu/LLJJptcuts/mMuMuZoom",
	       double minMLL=84.0,
	       double maxMLL=98.0,
	       double ctr_xsecpb=3740.,
	       double xsecincpb=30.)
#else
// Alpgen binned MC->data matching
const int    rebinx  = 0;
const double ymin    = 5e-2;
const double ymax    = 500;
const double h1sf    = 1.0;
const double nevents = luminvpb;
const string format  = "scale factor = %4.2f+/-%4.2f\n";
const string str1    = "/Mu, Run2010A+B, 33.9/pb";
const string str2    = "Fall10 Z+Jets Alpgen";

const string str3    = "Summer10 TTbar";
const double bkxsec  = 94.3;
const double bknev   = 632010;
const double h3sf    = luminvpb*bkxsec/bknev;

void myChi2Fit(const char *filename1="data.root",
	       const char *filename2="fall10zjets.root",
	       const char *bckgrndfn="",
	       const char *path     ="hNu/cut4_LLJJpt/mMuMuZoom",
	       double minMLL=84.0,
	       double maxMLL=98.0,
	       double ctr_xsecpb=1.25,
	       double xsecincpb=0.01)
#endif
{
  //==================================================
  // Open Files, prep the histograms
  //==================================================

  TFile *f1 = new TFile(filename1);
  if (f1->IsZombie()) {
    cerr << "File failed to open, " << string(filename1) << endl;
    return;
  }
  TFile *f2 = new TFile(filename2);
  if (f2->IsZombie()) {
    cerr << "File failed to open, " << string(filename2) << endl;
    return;
  }

  TH1 *h1ref = (TH1 *)f1->Get(path);
  if (!h1ref) {
    cerr << "couldn't find " << path << " in " << string(filename1) << endl;
    return;
  }

  if (rebinx) h1ref->Rebin(rebinx);
  h1ref->Scale(h1sf);

  TH1 *h2mc = (TH1 *)f2->Get(path);
  if (!h2mc) {
    cerr << "couldn't find " << path << " in " << string(filename2) << endl;
    return;
  }

  if (rebinx) h2mc->Rebin(rebinx);

  TH1 *h3bk = NULL;
  if (strlen(bckgrndfn)) {
    TFile *f3 = new TFile(bckgrndfn);
    if (f3->IsZombie()) {
      cerr << "File failed to open, " << string(filename2) << endl;
      return;
    }
    h3bk = (TH1 *)f3->Get(path);
    if (!h3bk) {
      cerr << "couldn't find " << path << " in " << string(bckgrndfn) << endl;
      return;
    }
    if (rebinx) h3bk->Rebin(rebinx);
    h3bk->Scale(h3sf);
  }

  //============================================================
  // Scan over npoints around ctr_xsecpb with step xsecincpb
  // and calculate chi2, minimize.
  //============================================================

  //const double ctr_xsecpb = (2204.0/0.44);

  int istart=-12;
  int iend=12;
  int npoints=iend-istart+1;
  TVectorD vx(npoints), vy(npoints);

  h2mc->GetXaxis()->SetRangeUser(minMLL,maxMLL);
  int ifirst=h2mc->GetXaxis()->GetFirst();
  int ilast=h2mc->GetXaxis()->GetLast();
  //int ndf=ilast-ifirst+1;

  double minchi2=1e99;
  double xsec4min;
  for (int i=istart; i<=iend; i++) {
    double xsecpb      = ctr_xsecpb + (double)i*xsecincpb;
    double scalefactor = luminvpb * xsecpb/nevents;

    TH1 *testh = (TH1 *)h2mc->Clone();
    testh->Scale(scalefactor);
    if (h3bk) testh->Add(h3bk);
    double chi2 = computeChi2(testh,h1ref,ifirst,ilast);
    vx[i-istart] = xsecpb;
    vy[i-istart] = chi2;

    if (chi2 < minchi2) {
      minchi2=chi2;
      xsec4min=xsecpb;
    }
    printf ("%4d\t%7.1f\t%5.3e\t%6.2f\n",
	    i, xsecpb, scalefactor, chi2);

    delete testh;
  }

  cout << "xsec4min = "<<xsec4min<<endl;

  //============================================================
  // Plot the results
  //============================================================

  gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  c1->Divide(2);

  //----------------------------------------
  // Pad 1, plot the chi2 parabola and error estimate

  c1->cd(1);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);

  //gStyle->SetFillColor(10);
  TGraph *gr = new TGraph(vx,vy);
  gr->SetTitle("#chi^{2} Minimization of Z+jets MC to Data");
  //gr->GetXaxis()->SetTitle("Z+jets #sigma (pb)");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  gr->GetYaxis()->SetTitleOffset(1.6);
  gr->SetMarkerStyle(4);

  gr->Draw("AP");
  TF1 *parab = new TF1("parab","[1]+[2]*pow((x-[0]),2.)",
		       ctr_xsecpb-(istart*xsecincpb),
		       ctr_xsecpb+(iend*xsecincpb));
  parab->SetParameters(xsec4min,minchi2,1.);
  TFitResultPtr r  = gr->Fit(parab,"RPS");
  Double_t    par0 = r->Value(0); // retrieve the value for the parameter 0
  Double_t    par1 = r->Value(1); // retrieve the value for the parameter 1
  Double_t    par2 = r->Value(2); // retrieve the value for the parameter 2
  //Double_t    chi2 = r->Chi2();   // to retrieve the fit chi2
  //Double_t    err0 = r->Error(0); // retrieve the error for the parameter 0
  //Double_t    err1 = r->Error(1); // retrieve the error for the parameter 1
  //Double_t    err2 = r->Error(2); // retrieve the error for the parameter 2

  // reset minchi2 and xsec4min from the fit
  xsec4min = par0;
  minchi2  = par1;
  double xsecerror = sqrt(1./par2);

  TLine *lo = new TLine(xsec4min-xsecerror,minchi2,
			xsec4min-xsecerror,minchi2+1);
  lo->SetLineColor(2);
  lo->Draw("same");
  TLine *hi = new TLine(xsec4min+xsecerror,minchi2,
			xsec4min+xsecerror,minchi2+1);
  hi->SetLineColor(2);
  hi->Draw("same");
  gPad->Update();

  TLatex *plabel = new TLatex();
  plabel -> SetNDC();
  plabel -> DrawLatex(.3, .8, Form(format.c_str(),xsec4min,xsecerror));

  //----------------------------------------
  // Pad 2, plot the histos overlaid

  c1->cd(2);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.1);
  gPad->SetLogy(1);
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);

  TLegend *leg = new TLegend(.4,.75,.98,.9);

  double scalefactor = luminvpb * xsec4min/nevents;

  THStack *st = new THStack("stack","stack");
  if (h3bk) st->Add(h3bk);
  TH1 *testh = (TH1 *)h2mc->Clone("fred");
  testh->Scale(scalefactor);
  st->Add(testh);

  TH1 *stkh = (TH1 *)testh->Clone("ginger");
  stkh->Reset();

  //st->GetXaxis()->SetRangeUser(70., 115.);
  //st->GetYaxis()->SetRangeUser(.01, 100.);

  stkh->SetStats(0);
  stkh->GetYaxis()->SetLimits(ymin,ymax);
  stkh->GetYaxis()->SetRangeUser(ymin,ymax);
  stkh->GetXaxis()->SetTitle(h2mc->GetXaxis()->GetTitle());//"M(#mu#mu) (GeV)");
  stkh->GetYaxis()->SetTitle(Form("dN/%d GeV",rebinx ? (rebinx*10) : 10));
  stkh->GetYaxis()->SetTitleOffset(1.3);
  //stkh->SetFillColor(5);
  //stkh->SetLineColor(2);

  st->SetHistogram(stkh);
  st->SetMaximum(stkh->GetYaxis()->GetXmax());
  st->SetMinimum(stkh->GetYaxis()->GetXmin());

  st->Draw("A HIST"); //"F");

  if (h3bk) leg->AddEntry(h3bk,str3.c_str());
  leg->AddEntry(testh,str2.c_str());

  leg->AddEntry(h1ref,str1.c_str(),"LE");
  h1ref->SetStats(0);
  h1ref->Draw("E SAME");

  leg->Draw("same");
  leg->SetFillColor(10);
  leg->SetTextSize(.035);
  leg->Print();
  gPad->Update();
}
