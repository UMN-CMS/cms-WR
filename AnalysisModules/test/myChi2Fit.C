#include <iostream>
#include <string>
#include <vector>
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

#include "tdrstyle.C"
#include "drawStandardTexts.C"

using namespace std;

//======================================================================
// Got this from
// http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html

void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ",
	      bool include_delimiters=false)
{
  tokens.clear();

  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  if (include_delimiters && lastPos>0)
    tokens.push_back(str.substr(0,lastPos));

  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    lastPos = str.find_first_not_of(delimiters, pos);

    if (include_delimiters && pos!=string::npos) {
      tokens.push_back(str.substr(pos, lastPos-pos));
    } //else skip delimiters.

    // Find next delimiter
    pos = str.find_first_of(delimiters, lastPos);
  }
}                                                            // Tokenize

//======================================================================

TH1 *getHisto(const string& path)
{
  TH1 *h(0);
  vector<string> v_tokens;
  Tokenize(path,v_tokens,":");
  if (v_tokens.size() != 2) return h;

  TFile *f = new TFile(v_tokens[0].c_str());
  if (f->IsZombie()) {
    cerr << "File failed to open, " << string(v_tokens[0]) << endl;
  } else
    h = (TH1 *)f->Get(v_tokens[1].c_str());

  return h;
}

//======================================================================

double computeChi2(TH1 *h1, TH1 *h2,int ifirst,int ilast)
{
  double chi2=0.;
  for (int ibin=ifirst; ibin<=ilast; ibin++) {
    double bin1  = h1->GetBinContent(ibin);
    double bin2  = h2->GetBinContent(ibin);
    double bin1e = h1->GetBinError  (ibin);
    double bin2e = h2->GetBinError  (ibin);
    chi2 += (bin1-bin2)*(bin1-bin2)/((bin1e*bin1e)+(bin2e*bin2e));
    //cout << bin1 << " " << bin2 << "; ";
  }
  //cout << endl;
  return chi2;
}

const double luminvpb = 36.145;
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
const string str1    = "/Mu, Run2010A+B, 36.145/pb";
const string str2    = "Summer10 Z+Jets Madgraph";
void myChi2Fit(const char *filename1="bryansNewMuSkim_Run2010AandBcombinedNov19JSON_hnu-anal.root",
	       const char *filename2="summer10_zjets_7tev_madgraph_start36_v10-v2_hnu-anal.root",
	       const char *bckgrndfn="",
	       const char *path     ="hNu/LLJJptcuts/mMuMuZoom",
	       double minMLL=84.0,
	       double maxMLL=98.0,
	       double ctr_xsecpb=3740.,
	       double xsecincpb=30.)
#elif 1
// Alpgen binned MC->data matching
const int    rebinx  = 2;
const double ymin    = 5e-2;
const double ymax    = 1000;
const double h1sf    = 1.0;
const double h2presf = 1.28;
const double nevents = luminvpb;
const string format  = "multiplier = %4.2f+/-%4.2f\n";
const string str1    = "/Mu, Run2010A+B, 36.145/pb";
const string str2    = "Fall10 Z+Jets Alpgen";

const string str3    = "Other background";
//const double bkxsec  = 167;
//const double bknev   = 1165716;
//const double h3sf    = luminvpb*bkxsec/bknev;
const double h3sf    = 1.0;

void myChi2Fit(double ctr_xsecpb=1.11,
	       double xsecincpb=0.01,
	       const string& refpath="data.root:hNu/cut5_Vertex/mMuMuZoom",
	       const string& path2scale="fall10zjetsLO.root:hNu/cut5_Vertex/mMuMuZoom",
	       const string& otherbckgrnd="zpeakbackgrnd.root:sumback_mMuMuZoom",
	       double minMLL=72.0,
	       double maxMLL=112.0)
//double minMLL=84.0,
//double maxMLL=100.0)
#else
// Alpgen binned MC->data matching for Electron channel
const int    rebinx  = 2;
const double ymin    = 5e-2;
const double ymax    = 1000;
const double h1sf    = 1.0;
const double h2presf = 0.94;
const double nevents = luminvpb;
const string format  = "multiplier = %4.2f+/-%4.2f\n";
const string str1    = "Run2010A+B, 36.145/pb";
const string str2    = "Fall10 Z+Jets Alpgen";

const string str3    = "Other background";
const double h3sf    = 1.0;

void myChi2Fit(double ctr_xsecpb=0.95, // 1.02,
	       double xsecincpb=0.01,
	       const string& refpath="fit-input-Mll-10.root:Data",
	       const string& path2scale="fit-input-Mll-10.root:Background;2",
	       const string& otherbckgrnd="fit-input-Mll-10.root:Background;1",
	       //const string& otherbckgrnd="",
	       double minMLL=70.0,
	       double maxMLL=110.0)
//double minMLL=84.0,
//double maxMLL=100.0)
#endif
{
  setTDRStyle();

  //==================================================
  // Open Files, prep the histograms
  //==================================================

  TH1 *h1ref = getHisto(refpath);
  if (!h1ref) {
    cerr << "couldn't find " << refpath << endl;
    return;
  }
  h1ref->UseCurrentStyle();

  if (rebinx) h1ref->Rebin(rebinx);
  h1ref->Scale(h1sf);

  TH1 *h2mc = getHisto(path2scale);
  if (!h2mc) {
    cerr << "couldn't find " << path2scale << endl;
    return;
  }

  h2mc->Scale(h2presf);

  h2mc->UseCurrentStyle();

  if (rebinx) h2mc->Rebin(rebinx);

  TH1 *h3bk = NULL;
  if (otherbckgrnd.size()) {
    h3bk = getHisto(otherbckgrnd);
    if (!h3bk) {
      cerr << "couldn't find " << otherbckgrnd << endl;
      return;
    }
    if (rebinx) h3bk->Rebin(rebinx);
    h3bk->Scale(h3sf);
    h3bk->UseCurrentStyle();
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

  TAxis *xax = h2mc->GetXaxis();
  xax->SetRangeUser(40,150); //(minMLL,maxMLL);
  int ifirst=xax->FindBin(minMLL);
  int ilast=xax->FindBin(maxMLL);
  ilast--;
  int ndf=ilast-ifirst+1;

  cout << "Minimizing over interval " << xax->GetBinLowEdge(ifirst);
  cout << ","                         << xax->GetBinUpEdge(ilast) << endl;
  double minchi2=1e99;
  double xsec4min=0;
  printf ("%4s\t%7s\t%9s\t%6s\n","i","xsecpb","sf","chi2");
  for (int i=istart; i<=iend; i++) {
    double xsecpb      = ctr_xsecpb + (double)i*xsecincpb;
    double scalefactor = luminvpb * xsecpb/nevents;

    TH1 *testh = (TH1 *)h2mc->Clone();
    testh->Scale(scalefactor);
    if (h3bk) testh->Add(h3bk);
    //testh->Sumw2();
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

  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  c1->Divide(2);

  //----------------------------------------
  // Pad 1, plot the chi2 parabola and error estimate

  c1->cd(1);
  //gPad->SetRightMargin(0.05);
  //gPad->SetLeftMargin(0.15);

  //gStyle->SetFillColor(10);
  TGraph *gr = new TGraph(vx,vy);
  gr->UseCurrentStyle();
  //gr->SetTitle("#chi^{2} Minimization of Z+jets MC to Data");
  gr->GetXaxis()->SetTitle("Multiplier for NNLO Z+jets #sigma");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  //gr->GetYaxis()->SetTitleOffset(1.6);
  //gr->SetMarkerStyle(4);

  gr->Draw("AP");
  gPad->Update();
  gStyle->SetOptFit(0);

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
  gPad->SetRightMargin(0.05);
  //gPad->SetTopMargin(0.1);
  gPad->SetLogy(1);
  //gPad->SetGridx(1);
  //gPad->SetGridy(1);

  TLegend *leg = new TLegend(.4,.8,.95,.95);

  double sf4min = luminvpb * xsec4min/nevents;

  TH1 *mcScaled = (TH1 *)h2mc->Clone("fred");
  mcScaled->Scale(sf4min);
  mcScaled->SetFillColor(5);
  //mcScaled->Sumw2();
  //TH1 *sum = (TH1 *)mcScaled->Clone("nancy");
  //if (h3bk) sum->Add(h3bk);

  THStack *sum = new THStack("stack","stack");
  if (h3bk) {
    //h3bk->Sumw2();
    h3bk->SetFillColor(4);
    sum->Add(h3bk);
  }
  sum->Add(mcScaled);

  sum->Draw("HIST F");

  sum->GetXaxis()->SetRangeUser(70., 110.);
  //sum->GetYaxis()->SetRangeUser(.01, 100.);

  //sum->SetStats(0);
  sum->GetYaxis()->SetLimits(ymin,ymax);
  sum->GetYaxis()->SetRangeUser(ymin,ymax);
  sum->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  //sum->GetXaxis()->SetTitle("M(ee) (GeV)");
  sum->GetYaxis()->SetTitle(Form("dN/%d GeV",rebinx ? (rebinx*4) : 4));
  //sum->GetYaxis()->SetTitleOffset(1.3);
  //sum->SetLineColor(2);

  //sum->Draw();
  gPad->Update();

  leg->AddEntry(h1ref,str1.c_str(),"PE");

  leg->AddEntry(mcScaled,str2.c_str(),"F");
  if (h3bk) leg->AddEntry(h3bk,str3.c_str(),"F");

  h1ref->SetStats(0);
  h1ref->Draw("E SAME");

  leg->Draw("same");
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  //leg->SetTextSize(.035);
  leg->Print();
  gPad->Update();

  drawStandardText("CMS Preliminary 2010",.77,.77);
}
