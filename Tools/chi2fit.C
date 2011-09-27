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

#include "hnuUtils.C" // getHisto

using namespace std;

//======================================================================

double computeChi2(TH1 *h1mc, TH1 *h2dt,int ifirst,int ilast)
{
  double chi2=0.;
  for (int ibin=ifirst; ibin<=ilast; ibin++) {
    double bin1  = h1mc->GetBinContent(ibin);
    double bin2  = h2dt->GetBinContent(ibin);
    // In case all the histogram manipulation screws up the errors
    // Take MC uncertainty as sqrt(bin contents)
    // double bin1e = h1mc->GetBinError  (ibin);
    double bin1e = sqrt(h1mc->GetBinContent(ibin));
    double bin2e = h2dt->GetBinError  (ibin);

    // if data statistics are too low, take the estimate of the error
    // from the MC number;
    //
    if (bin2<4.0)
      bin2e = max(bin2e,sqrt(bin1));

    chi2 += (bin1-bin2)*(bin1-bin2)/((bin1e*bin1e)+(bin2e*bin2e));
    //printf("%3d %5.3g %5.3g %5.3g %5.3g %5.3g\n",ibin,bin1,bin2,bin1e,bin2e,chi2);
  }

  return chi2;
}

//======================================================================

const double h1sf = 1.0;
void chi2fit(
	     const string& refpath="data.root:hNu/cut4_Mu1HighPt/mMuMuZoom",
	     const string& path2scale="zjets.root:hNu2011Z30/cut6_Mu1HighPt/mMuMuZoom",
	     const string& otherbckgrnd="zjetsBG.root:sumback_mMuMu_cut6",
	     const bool zfit=true,
	     const bool makeplots=false,
	     const string& plottag="test",
	     double ctr_xsecpb=1.20,
	     double xsecincpb=0.01,
	     double minMLL=70.0,
	     double maxMLL=110.0,
	     double luminvpb=204.0,
	     const int rebinx=4,
	     const double h2sf=1.0,
	     const double h3sf=1.0,
	     const double plotxmin=70.,
	     const double plotxmax=110.0,
	     const double plotymin=0.01,
	     const double plotymax=500.0
	     )
{
  setTDRStyle();

  //==================================================
  // Open Files, prep the histograms
  //==================================================
  std::cout << "Files/Histograms: " << std::endl ; 
  std::cout << refpath << std::endl ; 
  std::cout << path2scale << std::endl ; 
  std::cout << otherbckgrnd << std::endl ; 

  const double nevents = luminvpb;

  char s[100] ;
  sprintf(s,"Data, %i pb^{-1}",(int)luminvpb) ;
  const std::string str1(s) ; 

  const string format   = "multiplier = %4.2f+/-%4.2f\n";
  const string aformat  = "multiplier = %4.2f^{+%4.2f}_{-%4.2f}\n";
  const string str2     = (zfit) ? ("Spring11 Z+Jets (Alpgen)"):("Spring11 t#bar{t}+Jets (Madgraph)") ; 
  const string str3     = "Other background";

  const double zbins[13] = {0,70,74,78,82,86,90,94,98,102,106,110,400} ; 
  const double tbins[9]  = {0,20,60,100,140,180,220,260,300} ; 

  TH1 *h1ref = getHisto(refpath);
  if (!h1ref) {
    cerr << "couldn't find " << refpath << endl;
    return;
  }
  h1ref->UseCurrentStyle();

  // kluge
  if (rebinx) {
      if (rebinx == 999) {
          std::string newName = "tmpHisto" ;
          TH1D* tmpHisto = (TH1D*) h1ref->Clone() ;
          if (zfit) h1ref = (TH1D*) tmpHisto->Rebin(12,newName.c_str(),zbins) ;
          else      h1ref = (TH1D*) tmpHisto->Rebin(8,newName.c_str(),tbins) ;
      } else { 
          std::cout << "Currently kluged...no rescaling" << std::endl ; 
          // h1ref->Rebin(rebinx);
      }
  }
  h1ref->Scale(h1sf);

  TH1 *h2mc = getHisto(path2scale);
  if (!h2mc) {
    cerr << "couldn't find " << path2scale << endl;
    return;
  }

  for (int i=1; i<=h2mc->GetNbinsX(); i++) { 
    std::cout << i << ": " << h2mc->GetBinContent(i) << std::endl ; 
  }
  h2mc->Scale(h2sf);

  for (int i=1; i<=h2mc->GetNbinsX(); i++) { 
    std::cout << i << ": " << h2mc->GetBinContent(i) << std::endl ; 
  }
  h2mc->UseCurrentStyle();

  if (rebinx) {
      if (rebinx == 999) {
          std::string newName = "tmpHisto" ;
          TH1D* tmpHisto = (TH1D*) h2mc->Clone() ;
          if (zfit) h2mc = (TH1D*) tmpHisto->Rebin(12,newName.c_str(),zbins) ;
          else      h2mc = (TH1D*) tmpHisto->Rebin(8,newName.c_str(),tbins) ;
      } else { 
          h2mc->Rebin(rebinx);
      }
  }

  for (int i=1; i<=h2mc->GetNbinsX(); i++) { 
    std::cout << i << ": " << h2mc->GetBinContent(i) << std::endl ; 
  }

  TH1 *h3bk = NULL;
  if (otherbckgrnd.size()) {
    h3bk = getHisto(otherbckgrnd);
    if (!h3bk) {
      cerr << "couldn't find " << otherbckgrnd << endl;
      return;
    }
    if (rebinx) {
      if (rebinx == 999) {
          std::string newName = "tmpHisto" ;
          TH1D* tmpHisto = (TH1D*) h3bk->Clone() ;
          if (zfit) h3bk = (TH1D*) tmpHisto->Rebin(12,newName.c_str(),zbins) ;
          else      h3bk = (TH1D*) tmpHisto->Rebin(8,newName.c_str(),tbins) ;
      } else {
        h3bk->Rebin(rebinx);
      }
    }
    h3bk->Scale(h3sf);
    h3bk->UseCurrentStyle();
  }

  //============================================================
  // Scan over npoints around ctr_xsecpb with step xsecincpb
  // and calculate chi2, minimize.
  //============================================================

  int istart=-20;
  int iend=20;
  int npoints=iend-istart+1;
  TVectorD vx(npoints), vy(npoints);

  TAxis *xax = h2mc->GetXaxis();
  int ifirst=xax->FindBin(minMLL);
  int ilast=xax->FindBin(maxMLL);
  ilast--;
  // int ndf=ilast-ifirst+1;

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

  TGraph *gr = new TGraph(vx,vy);
  gr->UseCurrentStyle();

  if (zfit) gr->GetXaxis()->SetTitle("Multiplier for NNLO Z+jets #sigma");
  else      gr->GetXaxis()->SetTitle("Multiplier for t#bar{t}+jets #sigma");
  gr->GetYaxis()->SetTitle("#chi^{2}");
  gr->GetYaxis()->SetTitleOffset(1.1);

  gr->Draw("AP");
  gPad->Update();
  gStyle->SetOptFit(0);

  // Special calculation for asymmetric errors
  double minx = plotxmax ; double miny = plotymax ; 
  for (int i=0; i<gr->GetN(); i++) { 
      double x, y ; 
      gr->GetPoint(i,x,y) ;
      if ( y < miny ) { minx = x ; miny = y ; } 
  }
  double xloerror = plotxmin ; double xhierror = plotxmax ;
  double dlomin = 100. ; double dhimin = 100. ; 
  for (int i=0; i<gr->GetN(); i++) { 
      double x, y ; 
      gr->GetPoint(i,x,y) ;
      if ( x < minx && fabs(y-miny-1.0)<dlomin ) { dlomin = fabs(y-miny-1.0) ; xloerror = x ; } 
      if ( x > minx && fabs(y-miny-1.0)<dhimin ) { dhimin = fabs(y-miny-1.0) ; xhierror = x ; } 
  }
  
  
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
  // plabel -> DrawLatex(.3, .8, Form(format.c_str(),xsec4min,xsecerror));

  TLatex *palabel = new TLatex();
  palabel -> SetNDC();
  // palabel -> DrawLatex(.35, .8, Form(aformat.c_str(),minx,(minx-xloerror),(xhierror-minx)));

  if ( fabs( 2*minx - xloerror - xhierror ) < 0.02 )
      plabel -> DrawLatex(.3, .8, Form(format.c_str(),xsec4min,xsecerror));
  else
      palabel -> DrawLatex(.35, .8, Form(aformat.c_str(),minx,(xhierror-minx),(minx-xloerror)));
  
  //----------------------------------------
  // Pad 2, plot the histos overlaid

  c1->cd(2);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy(1);

  TLegend *leg = new TLegend(.30,.8,.95,.95);

  double sf4min = luminvpb * xsec4min/nevents;

  TH1 *mcScaled = (TH1 *)h2mc->Clone("fred");
  mcScaled->Scale(sf4min);
  mcScaled->SetFillColor(5);
  //mcScaled->Sumw2();
  //TH1 *sum = (TH1 *)mcScaled->Clone("nancy");
  //if (h3bk) sum->Add(h3bk);

  // mcScaled->Draw() ; 
  // return ; 

  THStack *sum = new THStack("stack","stack");
  if (h3bk) {
    //h3bk->Sumw2();
    h3bk->SetFillColor(4);
    sum->Add(h3bk);
  }
  sum->Add(mcScaled);

  TH1D* dummyHisto = new TH1D("dummyHisto","axis only",100,plotxmin,plotxmax) ; 
  dummyHisto->SetMaximum(plotymax) ; 
  dummyHisto->SetMinimum(plotymin) ;
  if (zfit) dummyHisto->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
  else      dummyHisto->GetXaxis()->SetTitle("M(e#mu) [GeV]");
  int rebinValue ; 
  if ( zfit ) { 
    if ( rebinx && rebinx != 999 ) rebinValue = rebinx * 4 ; 
    else rebinValue = 4 ; 
  } else { 
    rebinValue = 40 ; 
  }
    
  dummyHisto->GetYaxis()->SetTitle(Form("dN/%d GeV",rebinValue));
  dummyHisto->Draw() ; 

  sum->SetMaximum(plotymax) ; 
  sum->SetMinimum(plotymin) ; 
  sum->Draw("HIST SAME F");

  // Stacks are *so* annoying - can't do this until it's drawn
  //sum->GetXaxis()->SetTitle("M(ee) (GeV)");
  // sum->GetXaxis()->SetTitle("M(#mu#mu) (GeV)");
  // sum->GetYaxis()->SetTitle(Form("dN/%d GeV",rebinx ? (rebinx*4) : 4));

  // sum->GetYaxis()->SetRangeUser(plotymin,plotymax);

  // sum->GetXaxis()->SetRangeUser(plotxmin,plotxmax);
  // gPad->Update();

  // return ; 

  leg->AddEntry(h1ref,str1.c_str(),"PE");

  leg->AddEntry(mcScaled,str2.c_str(),"F");
  if (h3bk) leg->AddEntry(h3bk,str3.c_str(),"F");

  h1ref->SetStats(0);
  h1ref->Draw("E SAME");

  leg->Draw("same");
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  leg->Print();
  gPad->Update();

  // drawStandardText("CMS Preliminary",.77,.77);
  drawStandardText("CMS Preliminary",.285,.97);

  // bool printResults = true ; 
  if ( makeplots ) {

      std::string pdfName = "plots42x/znormChi2-" + plottag + ".pdf" ; 
      std::string epsName = "plots42x/znormChi2-" + plottag + ".eps" ; 
      std::string pngName = "plots42x/znormChi2-" + plottag + ".png" ; 
      std::string cName   = "plots42x/znormChi2-" + plottag + ".C" ; 
      
      c1->Print(pdfName.c_str()) ; 
      c1->Print(epsName.c_str()) ; 
      c1->Print(pngName.c_str()) ; 
      c1->Print(cName.c_str()) ; 
  }
}
