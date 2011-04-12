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
#include "TROOT.h"
#include "TFile.h"
#include "doLimitSetting.hh"
#include "hnuanalinput.h"

#include <algorithm>
#include <functional> // function std::less

using namespace RooFit;
using namespace RooStats;


//================================================================================

void setupWorkspace(RooWorkspace& w,
		    struct BaseInfo& baseinfo,
		    TH1* signalShape,
		    bool doSyst)
{
  char temp[1024]; // for sprintf's

  double lumipbinv=baseinfo.lumi; // pb^-1

  const double sig_scale_factor=baseinfo.xsec;

  // Scale background yield to the requested lumi:
  //
  const double tt_yield_mc   = lumipbinv * tt_yieldpb_mc;
  const double zj_yield_mc   = lumipbinv * zj_yieldpb_mc;
  //const double wj_yield_mc   = lumipbinv * wj_yieldpb_mc;
  //const double vv_yield_mc   = lumipbinv * vv_yieldpb_mc;
  //const double tw_yield_mc   = lumipbinv * tw_yieldpb_mc;
  //const double qcd_yield_data = lumipbinv * qcd_yieldpb_data;

  w.factory("mwr[520,2000]"); // observable is called mwr
  w.defineSet("obs","mwr");

  RooDataHist signals("signals","signals",RooArgList(*w.var("mwr")),signalShape) ;
  RooHistPdf signalPdf("sig_pdf","sig_pdf",RooArgList(*w.var("mwr")),signals);

  // Gaussian signal
  w.import(signalPdf); // signalPdf can go out of scope now

  // ttbar background (from fit)
  sprintf(temp,"Exponential::bkg_tt(mwr,%f)", tt_expfit_p1);  w.factory(temp);
  sprintf(temp,"Exponential::bkg_zj(mwr,%f)", zj_expfit_p1);  w.factory(temp);
  //sprintf(temp,"Exponential::bkg_wj(mwr,%f)", wj_expfit_p1);  w.factory(temp);
  //sprintf(temp,"Exponential::bkg_vv(mwr,%f)", vv_expfit_p1);  w.factory(temp);
  //sprintf(temp,"Exponential::bkg_tw(mwr,%f)", tw_expfit_p1);  w.factory(temp);
  //sprintf(temp,"Exponential::bkg_qcd(mwr,%f)", qcd_expfit_p1);  w.factory(temp);
  
  // total model with signal and background yields as parameters ( *beware when changing luminosity* )

  w.factory(
"SUM::main_pdf(\
sig_yield[20,0,3000]*sig_pdf,\
tt_yield[0.1,1e-5,10]*bkg_tt,\
zj_yield[0.1,1e-5,10]*bkg_zj)");
//vv_yield[0.1,1e-5,10]*bkg_vv)");
//tw_yield[0.1,1e-5,10]*bkg_tw)");
  //wj_yield[0.1,1e-5,10]*bkg_wj,		\
  //qcd_yield[0.1,1e-5,10]*bkg_qcd)");

  if (doSyst) {
    // note the "+1" required for the defintion of K! ( K - 1   = [fractional error] )
    sprintf(temp,"Lognormal::tt_xsec(tt_estimate[%f],tt_yield,%f)",tt_yield_mc,tt_yield_ferror+1); w.factory(temp);
    sprintf(temp,"Lognormal::zj_xsec(zj_estimate[%f],zj_yield,%f)",zj_yield_mc,zj_yield_ferror+1); w.factory(temp);
    //sprintf(temp,"Lognormal::wj_xsec(wj_estimate[%f],wj_yield,%f)",wj_yield_mc,wj_yield_ferror+1); w.factory(temp);
    //sprintf(temp,"Lognormal::vv_xsec(vv_estimate[%f],vv_yield,%f)",vv_yield_mc,vv_yield_ferror+1); w.factory(temp);
    //sprintf(temp,"Lognormal::tw_xsec(tw_estimate[%f],tw_yield,%f)",tw_yield_mc,tw_yield_ferror+1); w.factory(temp);
  }

  w.var("tt_yield")->setVal(tt_yield_mc);
  w.var("zj_yield")->setVal(zj_yield_mc);
  //w.var("wj_yield")->setVal(wj_yield_mc);
  //w.var("vv_yield")->setVal(vv_yield_mc);
  //w.var("tw_yield")->setVal(tw_yield_mc);
  //w.var("qcd_yield")->setVal(qcd_yield_data);

  baseinfo.signal=baseinfo.signal_eff*sig_scale_factor*lumipbinv;

  //  printf("Components: %f %f %f \n",signal_eff,sig_scale_factor,lumipbinv);
  baseinfo.background=(zj_yield_mc+
		       tt_yield_mc);
		       //wj_yield_mc+
		       //vv_yield_mc);
		       //tw_yield_mc);
		       //qcd_yield_data);
  //  w.var("bkg_scale")->setVal(info.background);

  if (!doSyst) {
    w.var("tt_yield")->setConstant();
    w.var("zj_yield")->setConstant();
    //w.var("wj_yield")->setConstant();
    //w.var("tw_yield")->setConstant();
    //w.var("vv_yield")->setConstant();
  }
  // w.var("qcd_yield")->setConstant();
  //  w.var("bkg_scale")->setConstant();

  // The model for the control sample that constrains the background.
  //w.factory("Gaussian::control_pdf(control_meas[50],bkg_yield,10.)");
  
  if (doSyst) {
    //w.factory("PROD::control_pdf(tt_xsec,zj_xsec,wj_xsec,vv_xsec,tw_xsec)");
    w.factory("PROD::control_pdf(tt_xsec,zj_xsec)");

    //sprintf(temp,"Gaussian::control_pdf(tt_estimate[%f],tt_yield,%f)",tt_yield_mc,tt_yield_mc*tt_yield_ferror);
    //w.factory(temp);
    // The total model including the main measurement and the control sample

    w.factory("PROD::main_with_control(main_pdf,control_pdf)");
  }
    
  // define sets for reference later
  w.defineSet("poi","sig_yield");
  if (doSyst) {
    //w.defineSet("nuis","tt_yield,zj_yield,wj_yield,tw_yield,vv_yield");
    w.defineSet("nuis","tt_yield,zj_yield");
  }
}                                                                // setupWorkspace

//================================================================================

void calc_CLsb_CLb(SamplingDistribution *sd_sb,
		   SamplingDistribution *sd_b,
		   struct CLInfo& cl_sb,
		   struct CLInfo& cl_b,
		   double xobs)
{
  double mean_sb = sd_sb->InverseCDFInterpolate(0.5); printf("mean_sb = %f\n", mean_sb);
  double mean_b  =  sd_b->InverseCDFInterpolate(0.5); printf("mean_b  = %f\n", mean_b);

  //  cl_sb_obs = res->CLsplusb();
  //  cl_b_obs  = 1.0-res->CLb();

  cl_sb.obs     = sd_sb->Integral(-10000,xobs);
  cl_sb.exp     = sd_sb->Integral(-10000,mean_b);

  cl_b.obs      =  sd_b->Integral(xobs,10000);
  cl_b.exp      =  sd_b->Integral(mean_sb,10000);

  cl_sb.exp_p1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.3413));
  cl_sb.exp_m1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.3413)+0.0001);
  cl_sb.exp_p2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.4773));
  cl_sb.exp_m2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.4773)+0.0001);

  cl_b.exp_p1s  =  sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.3413),10000);
  cl_b.exp_m1s  =  sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.3413),10000);
  cl_b.exp_p2s  =  sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.4773),10000);
  cl_b.exp_m2s  =  sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.4773),10000);

}                                                                 // calc_CLsb_CLb

//================================================================================

void calc_CLs(const SamplingDistribution *sd_sb,
	      const SamplingDistribution *sd_b,
	      struct CLInfo& cls)
{
  // Gets "ntoys" samples from the PDFs; each element is 
  // a value of the test statistic
  //
  std::vector<double> sampl_sb=sd_sb->GetSamplingDistribution();
  std::vector<double> sampl_b =sd_b ->GetSamplingDistribution();

  printf("i\tb\t\ts+b\n");
  for (size_t i=0; i<sampl_b.size(); i++) {
    printf("%u\t%f\t%f\n",i,sampl_b[i],sampl_sb[i]);
  }

  // Already pre-sorted, but what the 7734
  std::sort(sampl_sb.begin(),sampl_sb.end(),std::less<double>());
  std::sort(sampl_b.begin(),sampl_b.end(),std::less<double>());
  int nsamp=int(sampl_b.size());

  //  double median=sampl_b[sampl_b.size()/2];
  double clmed=-1,clmed_m2s=-1,clmed_m1s=-1,clmed_p1s=-1,clmed_p2s=-1;

  int p1s_pos = int(nsamp*(0.5+0.3413)+0.5)+1; // plus  1 sigma index position in a sample
  int m1s_pos = int(nsamp*(0.5-0.3413)+0.5)+1; // minus 1 sigma index position in a sample
  int p2s_pos = int(nsamp*(0.5+0.4773)+0.5)+1; // plus  2 sigma index position in a sample
  int m2s_pos = int(nsamp*(0.5-0.4773)+0.5)+1; // minus 2 sigma index position in a sample

  printf ("m2s_pos = %d\n",m2s_pos);
  printf ("m1s_pos = %d\n",m1s_pos);
  printf ("p1s_pos = %d\n",p1s_pos);
  printf ("p2s_pos = %d\n",p2s_pos);

  float val=sampl_b[nsamp/2]; // if this ain't the mean, it's pretty close!

  printf ("val = %f\n",val);

  // clmed, etc. get values between 0 and 2
  for (int i=0; i<nsamp; i++) {
    if (sampl_sb[i]>=val              && clmed    <0){clmed    =std::max(0.0f,float(i)/(nsamp/2));printf("CLMED %d %f %f %f\n",
												       i,clmed,sampl_sb[i],val);}
    if (sampl_sb[i]>=sampl_b[m1s_pos] && clmed_m1s<0) clmed_m1s=std::max(0.0f,float(i)/(nsamp/2));
    if (sampl_sb[i]>=sampl_b[m2s_pos] && clmed_m2s<0) clmed_m2s=std::max(0.0f,float(i)/(nsamp/2));
    if (sampl_sb[i]>=sampl_b[p1s_pos] && clmed_p1s<0) clmed_p1s=std::max(0.0f,float(i)/(nsamp/2));
    if (sampl_sb[i]>=sampl_b[p2s_pos] && clmed_p2s<0) clmed_p2s=std::max(0.0f,float(i)/(nsamp/2));
  }

  printf("%f %f %f %f %f \n",sampl_b[m2s_pos],sampl_b[m1s_pos],val,sampl_b[p1s_pos],sampl_b[p2s_pos]);
  if (clmed    ==-1.0 || clmed     > 1) clmed    =1.0;
  if (clmed_m1s==-1.0 || clmed_m1s > 1) clmed_m1s=1.0;
  if (clmed_m2s==-1.0 || clmed_m2s > 1) clmed_m2s=1.0;
  if (clmed_p1s==-1.0 || clmed_p1s > 1) clmed_p1s=1.0;
  if (clmed_p2s==-1.0 || clmed_p2s > 1) clmed_p2s=1.0;
  
  if (clmed>1.0 && clmed<1.001) clmed=1.0;
  if (clmed<0.0) clmed=0.0;

  cls.exp=clmed;
  cls.exp_m2s=clmed_m2s;
  cls.exp_m1s=clmed_m1s;
  cls.exp_p1s=clmed_p1s;
  cls.exp_p2s=clmed_p2s;
  
  // for (int i=0; i<20; i++) 
  //  printf("%3d %f %f \n",i,sampl_sb[i],sampl_b[i]);
}                                                                     // calc_CLs

//================================================================================

void dumpInfo(const LimitPointStruct& info)
{
  printf("SUMMARY  : LUMI=%.1f ipb MW=%.0f GeV  MNu=%.0f GeV  XSEC=%.2e pb \n",
	 info.base.lumi,info.base.mwr,info.base.mnu,info.base.xsec);
  
  //  printf("DETAIL : mean_sb=%f  mean_b=%f\n", mean_sb,mean_b);
  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",
	 info.cl_sb.exp,info.cl_b.obs,info.cl_b.exp,info.cl_b.obs);

  printf("SUMMARY3 : EXP(s+b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",
	 info.cl_sb.exp_m2s,info.cl_sb.exp_m1s,info.cl_sb.exp_p1s,info.cl_sb.exp_p2s);

  printf("SUMMARY3 : EXP(b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",
	 info.cl_b.exp_m2s,
	 info.cl_b.exp_m1s,
	 info.cl_b.exp_p1s,
	 info.cl_b.exp_p2s);

  printf("SUMMARY3 : CLS %f | %f %f %f %f %f\n",
	 info.cls.obs,
	 info.cls.exp,
	 info.cls.exp_p2s,
	 info.cls.exp_p1s,
	 info.cls.exp_m1s,
	 info.cls.exp_m2s);

  printf("SUMMARY4: %.2f background / %.2f (%.2f%%) signal / %d data\n",
	 info.base.background,
	 info.base.signal,
	 info.base.signal_eff*100, 
	 int(info.base.data));
}                                                                      // dumpInfo

//================================================================================

#ifndef MAKE_MAIN
void drawPlot(const RooWorkspace& w,
	      HypoTestResult *res,
	      RooDataHist* data,
	      double norm)
{
  gROOT->SetStyle("Plain");
  TCanvas* c2=new TCanvas("c2","c2",900,450);
  //c2->Divide(1,2);

  //c2->cd(1);
  // and make a plot. number of bins is optional (default: 100)
  HypoTestPlot *plot = new HypoTestPlot(*res, 80); 
  plot->Draw();

#if 0
  c2->cd(2);
  RooPlot* mesframe = w.var("mwr")->frame() ;

  data->plotOn(mesframe);
  //  w.pdf("bkg_pdf")->plotOn(mesframe);
  //w.pdf("main_pdf")->plotOn(mesframe);
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(norm));
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(norm),
			    Components("bkg_tt"),LineStyle(kDashed));

  mesframe->Draw();
#endif
}                                                                      // drawPlot
#endif // !MAKE_MAIN

//================================================================================

void doLimitSetting(TFile* dataf, TFile* signal, int ntoys, LimitPointStruct& info, bool doSyst) {

  // I N I T   D A T A ,   B A C K G R O U N D ,   S I G N A L   H I S T O S

  // extract histograms from input files
  //
  TH1* signalShape=(TH1*)(signal->Get(finalVarHist)->Clone("ss"));
  TH1* signalNorm=(TH1*)(signal->Get("hNu/njet"));

  //  signalShape->Print();
  info.base.signal_eff=signalShape->Integral()*1.0/signalNorm->GetEntries();

  signalShape->Scale(1.0/signalShape->Integral());

  TH1* dataHist=(TH1*)(dataf->Get(finalVarHist)->Clone("dataHist"));
  info.base.data=dataHist->Integral(int(dataHist->FindBin(520)),
				    int(dataHist->FindBin(2001)));

  // S E T U P   T H E   W O R K S P A C E

  RooWorkspace w("mumujj");

  setupWorkspace(w,info.base,signalShape,doSyst);

  RooDataHist* data=new RooDataHist("data","data",RooArgList(*w.var("mwr")),dataHist) ;

  // choose which pdf you want to use
  RooAbsPdf* pdfToUse=w.pdf("main_pdf"); // only use main measurement

  if (doSyst) pdfToUse=w.pdf("main_with_control");

  //  RooDataSet* data = w.pdf("main_pdf")->generate(*w.set("obs"),RooFit::Extended());
  //  data->Print();

  // D E F I N E  N U L L  &  A L T E R N A T I V E   H Y P O T H E S E S  

  ModelConfig b_model("B_model", &w); // background only, null hypothesis
  b_model.SetPdf(*pdfToUse);
  b_model.SetObservables(*w.set("obs"));
  b_model.SetParametersOfInterest(*w.set("poi"));
  if (doSyst)
    b_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(0.0);
  b_model.SetSnapshot(*w.set("poi"));  

  ModelConfig sb_model("S+B_model", &w); // signal+background, alternative hypothesis
  sb_model.SetPdf(*pdfToUse);
  sb_model.SetObservables(*w.set("obs"));
  sb_model.SetParametersOfInterest(*w.set("poi"));
  if (doSyst)
    sb_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(info.base.signal);
  sb_model.SetSnapshot(*w.set("poi"));

  
  SimpleLikelihoodRatioTestStat slrts(*b_model.GetPdf(),*sb_model.GetPdf());
  slrts.SetNullParameters(*b_model.GetSnapshot());
  slrts.SetAltParameters(*sb_model.GetSnapshot());

  ToyMCSampler toymcsampler(slrts, ntoys);
  toymcsampler.SetGenerateBinned(true);

  HybridCalculator myH2(*data,sb_model, b_model, &toymcsampler);

  if (doSyst) {
    myH2.ForcePriorNuisanceNull(*w.pdf("control_pdf")); 
    myH2.ForcePriorNuisanceAlt(*w.pdf("control_pdf")); 
  }

  // R U N   T H E   H Y P O T H E S I S   T E S T ,   C A L C   R E S U L T S

  HypoTestResult *res = myH2.GetHypoTest();
  //  res->Print(0);

  // here we get the expected limit
  SamplingDistribution* sd_sb=res->GetAltDistribution();
  SamplingDistribution* sd_b=res->GetNullDistribution();

  //  printf("Sizes: %d %d %f\n",sd_b->GetSize(),sd_sb->GetSize());

  double xobs=res->GetTestStatisticData();
  xobs+=0.0001; // important for handling the "exactly zero data" case, irrelevant for all other purposes

  calc_CLsb_CLb(sd_sb,sd_b,
		info.cl_sb,info.cl_b,
		xobs);

  calc_CLs(sd_sb,sd_b,info.cls);

  info.cls.obs=res->CLs();

  dumpInfo(info);
  
#ifndef MAKE_MAIN
  drawPlot(w,res,data,info.base.signal+info.base.background);
#endif

  delete sd_sb;
  delete sd_b;
  delete res;

  delete dataHist;
  delete signalShape;
  delete data;
}

void testLimitSetting(TFile* f0, TFile* f1, double xsec=1.0, bool doSyst=true) {
  LimitPointStruct lps;
  lps.base.xsec=xsec;
  lps.base.lumi=36.1;

  doLimitSetting(f0,f1,2000,lps,doSyst);
}
