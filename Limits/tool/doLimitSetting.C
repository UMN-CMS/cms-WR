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

  const char* finalVarHist="hNu/diLmasscut/mWR";

  const double sig_scale_factor=info.xsec;
  const double ttb_yield_mc=lumi*194.3*(0.11*0.11)*(3818+6)/193317; // XSEC from CMS 194 ± 72 (stat) ± 24 (syst) ± 21 (lumi), plus W->mu nu
  //  const double zj_yield_mc=lumi*5454.0*37/1647472;  // XSEC from P. Dudero
  const double zj_yield_mc=lumi*(4.5/7)*(3353+18)/92291;  //M180 mumu (Match to data)
  
  const double wj_yield_mc=lumi*53711.0*1/1.11e7;   // change from 3 to 1 with higher MLL cut
  const double qcd_yield_data=lumi*0.06514/100/0.003923*exp(520*(-0.003923));

  RooWorkspace w("example");
  w.factory("mwr[520,2000]"); // observable is called mwr
  w.defineSet("obs","mwr");

  TH1* signalShape=(TH1*)(signal->Get(finalVarHist)->Clone("ss"));
  TH1* signalNorm=(TH1*)(signal->Get("hNu/njet"));

  //  signalShape->Print();
  double signal_eff=signalShape->Integral()*1.0/signalNorm->GetEntries();
  
  signalShape->Scale(1.0/signalShape->Integral());

  TH1* dataHist=(TH1*)(dataf->Get(finalVarHist)->Clone("dataHist"));


  info.data=dataHist->Integral(int(dataHist->FindBin(520)),int(dataHist->FindBin(2001)));
  RooDataHist* data=new RooDataHist("data","data",RooArgList(*w.var("mwr")),dataHist) ;

  //  data->Print();

  RooDataHist* signals=new RooDataHist("signals","signals",RooArgList(*w.var("mwr")),signalShape) ;
  RooHistPdf* signalPdf=new RooHistPdf("sig_pdf","sig_pdf",RooArgList(*w.var("mwr")),*signals);

  // Gaussian signal
  w.import(*signalPdf);
  //  w.import(*datardh);
		       //  w.factory("Gaussian::sig_pdf(mwr,sig_mean[600],sig_sigma[10])");

  // ttbar background (from fit)
  w.factory("Exponential::ttb_exp(mwr,-5.22e-3)");
  w.factory("Uniform::ttb_lin(mwr)");
  w.factory("SUM::bkg_ttb(8.34*ttb_exp,1.192*ttb_lin)");

  w.factory("Exponential::bkg_zj(mwr,-3.64e-3)");

  w.factory("Exponential::bkg_wj(mwr,-3e-3)");

  w.factory("Exponential::bkg_qcd(mwr,-0.003923)");
  
  w.factory("SUM::bkg_pdf(ttb_yield[20,0,300]*bkg_ttb,zj_yield[1,0,100]*bkg_zj,wj_yield[1,0,100]*bkg_wj,qcd_yield[1,0,30]*bkg_qcd)");

  // total model with signal and background yields as parameters
  w.factory("SUM::main_pdf(sig_yield[20,0,3000]*sig_pdf,bkg_scale[1,0,1000]*bkg_pdf)");

  w.var("ttb_yield")->setVal(ttb_yield_mc);
  w.var("zj_yield")->setVal(zj_yield_mc);
  w.var("wj_yield")->setVal(wj_yield_mc);
  w.var("qcd_yield")->setVal(qcd_yield_data);
  double sig_yield_expected=signal_eff*sig_scale_factor*lumi;
  info.background=(wj_yield_mc+zj_yield_mc+ttb_yield_mc+qcd_yield_data);
  w.var("bkg_scale")->setVal(info.background);


  w.var("ttb_yield")->setConstant();
  w.var("zj_yield")->setConstant();
  w.var("wj_yield")->setConstant();
  w.var("qcd_yield")->setConstant();
  w.var("bkg_scale")->setConstant();


  info.signal=sig_yield_expected;
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
  //  data->Print();

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
  toymcsampler.SetGenerateBinned(true);

  HybridCalculator myH2(*data,sb_model, b_model, &toymcsampler);

  //  myH2.ForcePriorNuisanceNull(*w.pdf("control_pdf")); 
  //  myH2.ForcePriorNuisanceAlt(*w.pdf("control_pdf")); 

  HypoTestResult *res = myH2.GetHypoTest();
  //  res->Print(0);

  // here we get the expected limit
  SamplingDistribution* sd_sb=res->GetAltDistribution();
  SamplingDistribution* sd_b=res->GetNullDistribution();

  double xobs=res->GetTestStatisticData();
  xobs+=0.0001; // important for handling the "exactly zero data" case, irrelevant for all other purposes
  double alt_obs_sb=sd_sb->Integral(-10000,xobs);
  double alt_obs_b=sd_b->Integral(xobs,10000);

  double mean_sb=sd_sb->InverseCDFInterpolate(0.5);
  double mean_b=sd_b->InverseCDFInterpolate(0.5);

  double expected_sb=sd_sb->Integral(-10000,mean_b);
  double expected_b=sd_b->Integral(mean_sb,10000);
  
  //  printf("DETAIL : mean_sb=%f  mean_b=%f\n", mean_sb,mean_b);
  printf("SUMMARY  : LUMI=%.1f ipb MW=%.0f GeV  MNu=%.0f GeV  XSEC=%.2e pb \n",info.lumi,info.mwr,info.mnu,sig_scale_factor);
  //  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",expected_sb,res->CLsplusb(),expected_b,1.0-res->CLb());
  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",expected_sb,alt_obs_sb,expected_b,alt_obs_b);
  //  info.cl_sb_obs = res->CLsplusb();
  info.cl_sb_obs = alt_obs_sb;
  info.cl_sb_exp = expected_sb;
  //  info.cl_b_obs = 1.0-res->CLb();
  info.cl_b_obs = alt_obs_b;
  info.cl_b_exp = expected_b;

  info.cl_sb_exp_p1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.3413));
  info.cl_sb_exp_m1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.3413)+0.0001);
  info.cl_sb_exp_p2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.4773));
  info.cl_sb_exp_m2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.4773)+0.0001);

  printf("SUMMARY3 : EXP(s+b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_sb_exp_m2s,info.cl_sb_exp_m1s,info.cl_sb_exp_p1s,info.cl_sb_exp_p2s);

  info.cl_b_exp_p1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.3413),10000);
  info.cl_b_exp_m1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.3413),10000);
  info.cl_b_exp_p2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.4773),10000);
  info.cl_b_exp_m2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.4773),10000);
  printf("SUMMARY3 : EXP(b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_b_exp_m2s,info.cl_b_exp_m1s,info.cl_b_exp_p1s,info.cl_b_exp_p2s);
  printf("SUMMARY4: %.2f background / %.2f (%.2f%%) signal / %d data\n",info.background,info.signal,signal_eff*100, int(info.data));


  
#ifndef MAKE_MAIN
  TCanvas* c2=new TCanvas("c2","c2",900,900);
  c2->Divide(1,2);

  c2->cd(1);
  // and make a plot. number of bins is optional (default: 100)
  HypoTestPlot *plot = new HypoTestPlot(*res, 80); 
  plot->Draw();

  c2->cd(2);
  RooPlot* mesframe = w.var("mwr")->frame() ;

  data->plotOn(mesframe);
  //  w.pdf("bkg_pdf")->plotOn(mesframe);
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(sig_yield_expected+info.background));
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(sig_yield_expected+info.background),Components("bkg_pdf"),LineStyle(kDashed));

  mesframe->Draw();

#endif

  delete dataHist;
  delete signalShape;
  delete signals;
  delete data;

}

void testLimitSetting(TFile* f0, TFile* f1, double xsec=1.0) {
  LimitPointStruct lps;
  lps.xsec=xsec;
  lps.lumi=34.0;

  doLimitSetting(f0,f1,5000,lps);
}
