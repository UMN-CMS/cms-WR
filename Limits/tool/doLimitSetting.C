#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooHistPdf.h"
#include "RooGlobalFunc.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestResult.h"
#include "TCanvas.h"
#include "TFile.h"
#include "doLimitSetting.hh"

using namespace RooFit;
using namespace RooStats;

void doLimitSetting(TFile* dataf, TFile* signal, int ntoys, LimitPointStruct& info) {
  double lumi=info.lumi; // pb^-1

  const char* finalVarHist="hNu/masscut/mWR";

  const double sig_scale_factor=info.xsec;
  const double ttf_yield_mc=lumi*94.3*611/632010;
  const double zj_yield_mc=lumi*2400.0*35/1647472;

  RooWorkspace w("example");
  w.factory("mwr[400,2000]"); // observable is called mwr
  w.defineSet("obs","mwr");

  TH1* dataHist=(TH1*)(dataf->Get(finalVarHist)->Clone("dataHist"));

  TH1* signalShape=(TH1*)(signal->Get(finalVarHist)->Clone("ss"));
  TH1* signalNorm=(TH1*)(signal->Get("hNu/njet"));

  //  signalShape->Print();
  double signal_eff=signalShape->Integral()*1.0/signalNorm->GetEntries();
  
  signalShape->Scale(1.0/signalShape->Integral());

  RooDataHist* data=new RooDataHist("data","data",RooArgList(*w.var("mwr")),dataHist) ;

  data->Print();

  RooDataHist* signals=new RooDataHist("signals","signals",RooArgList(*w.var("mwr")),signalShape) ;
  RooHistPdf* signalPdf=new RooHistPdf("sig_pdf","sig_pdf",RooArgList(*w.var("mwr")),*signals);

  // Gaussian signal
  w.import(*signalPdf);
  //  w.import(*datardh);
		       //  w.factory("Gaussian::sig_pdf(mwr,sig_mean[600],sig_sigma[10])");

  // ttbar background (from fit)
  w.factory("Exponential::ttb_exp(mwr,-5e-3)");
  w.factory("Uniform::ttb_lin(mwr)");
  w.factory("SUM::bkg_ttb(6.78*ttb_exp,0.336*ttb_lin)");

  w.factory("Exponential::bkg_zj(mwr,-5.2e-3)");
  
  w.factory("SUM::bkg_pdf(ttf_yield[20,0,300]*bkg_ttb,zj_yield[1,0,100]*bkg_zj)");

  // total model with signal and background yields as parameters
  w.factory("SUM::main_pdf(sig_yield[20,0,3000]*sig_pdf,bkg_scale[1,0,100]*bkg_pdf)");

  w.var("ttf_yield")->setVal(611*ttf_yield_mc);
  w.var("zj_yield")->setVal(35*zj_yield_mc);
  double sig_yield_expected=signal_eff*sig_scale_factor*lumi;

  w.var("ttf_yield")->setConstant();
  w.var("zj_yield")->setConstant();
  w.var("bkg_scale")->setConstant();

  // The model for the control sample that constrains the background.
  //w.factory("Gaussian::control_pdf(control_meas[50],bkg_yield,10.)");

  // The total model including the main measurement and the control sample
  //w.factory("PROD::main_with_control(main_pdf,control_pdf)");

  // choose which pdf you want to use
  RooAbsPdf* pdfToUse = w.pdf("main_pdf"); // only use main measurement
  //RooAbsPdf* pdfToUse = w.pdf("main_with_control"); // also include control sample
  
  // define sets for reference later
  w.defineSet("poi","sig_yield");
  //  w.defineSet("nuis","bkg_yield");

  //  RooDataSet* data = w.pdf("main_pdf")->generate(*w.set("obs"),RooFit::Extended());
  data->Print();

  // D E F I N E  N U L L  &  A L T E R N A T I V E   H Y P O T H E S E S  
  ModelConfig b_model("B_model", &w);
  b_model.SetPdf(*pdfToUse);
  b_model.SetObservables(*w.set("obs"));
  b_model.SetParametersOfInterest(*w.set("poi"));
  //  b_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(0.0);
  b_model.SetSnapshot(*w.set("poi"));  

  ModelConfig sb_model("S+B_model", &w);
  sb_model.SetPdf(*pdfToUse);
  sb_model.SetObservables(*w.set("obs"));
  sb_model.SetParametersOfInterest(*w.set("poi"));
  // sb_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(sig_yield_expected);
  sb_model.SetSnapshot(*w.set("poi"));

  
  SimpleLikelihoodRatioTestStat slrts(*b_model.GetPdf(),*sb_model.GetPdf());
  slrts.SetNullParameters(*b_model.GetSnapshot());
  slrts.SetAltParameters(*sb_model.GetSnapshot());

  ToyMCSampler toymcsampler(slrts, ntoys);

  HybridCalculator myH2(*data,sb_model, b_model, &toymcsampler);

  //  myH2.ForcePriorNuisanceNull(*w.pdf("control_pdf")); 
  //  myH2.ForcePriorNuisanceAlt(*w.pdf("control_pdf")); 

  HypoTestResult *res = myH2.GetHypoTest();
  //  res->Print(0);

  // here we get the expected limit
  SamplingDistribution* sd_sb=res->GetAltDistribution();
  SamplingDistribution* sd_b=res->GetNullDistribution();

  double mean_sb=sd_sb->InverseCDFInterpolate(0.5);
  double mean_b=sd_b->InverseCDFInterpolate(0.5);

  double expected_sb=sd_sb->Integral(-10000,mean_b);
  double expected_b=sd_b->Integral(mean_sb,10000);
  
  //  printf("DETAIL : mean_sb=%f  mean_b=%f\n", mean_sb,mean_b);
  printf("SUMMARY  : LUMI=%.1f ipb MW=%.0f GeV  MNu=%.0f GeV  XSEC=%.2e pb \n",info.lumi,info.mwr,info.mnu,sig_scale_factor);
  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",expected_sb,res->CLsplusb(),expected_b,1.0-res->CLb());
  info.cl_sb_obs = res->CLsplusb();
  info.cl_sb_exp = expected_sb;
  info.cl_b_obs = 1.0-res->CLb();
  info.cl_b_exp = expected_b;

  info.cl_sb_exp_p1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.3413));
  info.cl_sb_exp_m1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.3413));
  info.cl_sb_exp_p2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.4773));
  info.cl_sb_exp_m2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.4773));

  printf("SUMMARY3 : EXP(s+b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_sb_exp_m2s,info.cl_sb_exp_m1s,info.cl_sb_exp_p1s,info.cl_sb_exp_p2s);

  info.cl_b_exp_p1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.3413),10000);
  info.cl_b_exp_m1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.3413),10000);
  info.cl_b_exp_p2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.4773),10000);
  info.cl_b_exp_m2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.4773),10000);
  printf("SUMMARY3 : EXP(b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_b_exp_m2s,info.cl_b_exp_m1s,info.cl_b_exp_p1s,info.cl_b_exp_p2s);

  
#ifndef MAKE_MAIN
  TCanvas* c2=new TCanvas("c2","c2",900,900);
  c2->Divide(2,2);
  c2->cd(1);
  //  data->Draw("mwr");

  c2->cd(2);
  dataHist->Draw("HIST");

  c2->cd(3);
  // and make a plot. number of bins is optional (default: 100)
  HypoTestPlot *plot = new HypoTestPlot(*res, 80); 
  plot->Draw();

  c2->cd(4);
  RooPlot* mesframe = w.var("mwr")->frame() ;

  w.pdf("bkg_pdf")->plotOn(mesframe);
  w.pdf("main_pdf")->plotOn(mesframe);
  data->plotOn(mesframe);


  mesframe->Draw();

#endif

  delete dataHist;
  delete signalShape;
  delete signals;
  delete data;

}

void testLimitSetting(TFile* f0, TFile* f1) {
  LimitPointStruct lps;
  lps.xsec=0.7;
  lps.lumi=7.;

  doLimitSetting(f0,f1,400,lps);
}
