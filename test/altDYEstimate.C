#include "TStyle.h"
#include "TPaveText.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>

/**
 * steps to use this macro:
 * 1. process data, TTBar and DY MC with analysis.cpp with low dilepton mass sideband requirements and DY MLL scale factors
 * 2. process DY MC with analysis.cpp with signal region requirements and DY MLL scale factors
 * 3. modify the TString dir and treeName to point to correct file directory and tree, and update the root file name such that the correct files are used
 * 4. launch root with "root -l -b"
 * 5. execute ".L altDYEstimate.C"
 * 6. execute "altDYEstimate()" to calculate the scale factors, and compare the data driven DY estimate to DY MC in the signal region
 *
 * the scale factors which interpolate between MLLJJ bins in the sideband and MLLJJ bins in the signal region are saved as plots showing
 * the bin widths and scale factor values.
 *
 */

//std::string chiSquaredNdofString(TF1 * fit);
void fillHisto(TChain * chain, Selector *myEvent, TH1 * h);
void makePlot(TH1 *histo, TString xAxisLabel, TString yAxisLabel, TString plotTitle, std::string outputPlotFileName, TString canvName, Float_t maxYval, Float_t minYval);
void altDYEstimate(){

  ///change dir, treeName and the root file names listed below to read the correct files
  TString dir = "../analysisCppOutputRootFiles/";
  TString treeName = "Tree_Iter0";

  //sideband region files
  TChain * chain_DataMu_CR = new TChain(treeName);
  chain_DataMu_CR->Add(dir+"selected_tree_data_lowdileptonsidebandMuMu.root");	//mu channel data in sideband
  TChain * chain_DYMu_CR = new TChain(treeName);
  chain_DYMu_CR->Add(dir+"selected_tree_DYAMC_lowdileptonsidebandMuMu_withMllWeight.root");	//mu channel DY MC in sideband
  TChain * chain_DataEle_CR = new TChain(treeName);
  chain_DataEle_CR->Add(dir+"selected_tree_data_lowdileptonsidebandEE.root");	//ele channel data in sideband
  TChain * chain_DYEle_CR = new TChain(treeName);
  chain_DYEle_CR->Add(dir+"selected_tree_DYAMC_lowdileptonsidebandEE_withMllWeight.root");	//ele channel DY MC in sideband
  TChain * chain_TTMu_CR = new TChain(treeName);
  chain_TTMu_CR->Add(dir+"selected_tree_TT_lowdileptonsidebandMuMu.root");	//mu channel TTBar MC in sideband
  TChain * chain_TTEle_CR = new TChain(treeName);
  chain_TTEle_CR->Add(dir+"selected_tree_TT_lowdileptonsidebandEE.root");	//ele channel TTBar MC in sideband

  //signal region files
  TChain * chain_DYMu_SR = new TChain(treeName);
  chain_DYMu_SR->Add(dir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");	//mu channel DY MC in signal region
  TChain * chain_DYEle_SR = new TChain(treeName);
  chain_DYEle_SR->Add(dir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");	//ele channel DY MC in signal region

 
  Selector myEvent_DataMu_CR;
  Selector myEvent_DYMu_CR;
  Selector myEvent_DataEle_CR;
  Selector myEvent_DYEle_CR;
  Selector myEvent_DYMu_SR;
  Selector myEvent_DYEle_SR;
  Selector myEvent_TTMu_CR;
  Selector myEvent_TTEle_CR;
  
  myEvent_DataMu_CR.SetBranchAddresses(chain_DataMu_CR);
  myEvent_DYMu_CR.SetBranchAddresses(chain_DYMu_CR);
  myEvent_DataEle_CR.SetBranchAddresses(chain_DataEle_CR);
  myEvent_DYEle_CR.SetBranchAddresses(chain_DYEle_CR);
  myEvent_DYMu_SR.SetBranchAddresses(chain_DYMu_SR);
  myEvent_DYEle_SR.SetBranchAddresses(chain_DYEle_SR);
  myEvent_TTMu_CR.SetBranchAddresses(chain_TTMu_CR);
  myEvent_TTEle_CR.SetBranchAddresses(chain_TTEle_CR);
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///now declare and fill individual histos, and declare ratio histos made by dividing two or more filled histos
 
  /*variable bin widths*/
  Float_t bins[] = {600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 3000};	//from low MLL sideband
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1 *h_WR_mass_DataEleCR = new TH1F("h_WR_mass_DataEleCR","",binnum, bins);
  TH1 *h_WR_mass_DataMuCR = new TH1F("h_WR_mass_DataMuCR","",binnum, bins);
  TH1 *h_WR_mass_DYEleCR = new TH1F("h_WR_mass_DYEleCR","",binnum, bins);
  TH1 *h_WR_mass_DYMuCR = new TH1F("h_WR_mass_DYMuCR","",binnum, bins);
  TH1 *h_WR_mass_DYEleSR = new TH1F("h_WR_mass_DYEleSR","",binnum, bins);
  TH1 *h_WR_mass_DYMuSR = new TH1F("h_WR_mass_DYMuSR","",binnum, bins);
  TH1 *h_WR_mass_TTEleCR = new TH1F("h_WR_mass_TTEleCR","",binnum, bins);
  TH1 *h_WR_mass_TTMuCR = new TH1F("h_WR_mass_TTMuCR","",binnum, bins);
  /**/

  /*
  //fixed bin widths, just for testing
  Int_t fixedWidthBins = 8;
  TH1 *h_WR_mass_DataEleCR = new TH1F("h_WR_mass_DataEleCR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_DataMuCR = new TH1F("h_WR_mass_DataMuCR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_DYEleCR = new TH1F("h_WR_mass_DYEleCR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_DYMuCR = new TH1F("h_WR_mass_DYMuCR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_DYEleSR = new TH1F("h_WR_mass_DYEleSR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_DYMuSR = new TH1F("h_WR_mass_DYMuSR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_TTEleCR = new TH1F("h_WR_mass_TTEleCR","", fixedWidthBins, 600, 3000);
  TH1 *h_WR_mass_TTMuCR = new TH1F("h_WR_mass_TTMuCR","", fixedWidthBins, 600, 3000);
  */

  fillHisto(chain_DataEle_CR, &myEvent_DataEle_CR, h_WR_mass_DataEleCR);
  fillHisto(chain_DataMu_CR, &myEvent_DataMu_CR, h_WR_mass_DataMuCR);
  fillHisto(chain_DYEle_CR, &myEvent_DYEle_CR, h_WR_mass_DYEleCR);
  fillHisto(chain_DYMu_CR, &myEvent_DYMu_CR, h_WR_mass_DYMuCR);
  fillHisto(chain_DYEle_SR, &myEvent_DYEle_SR, h_WR_mass_DYEleSR);
  fillHisto(chain_DYMu_SR, &myEvent_DYMu_SR, h_WR_mass_DYMuSR);
  fillHisto(chain_TTEle_CR, &myEvent_TTEle_CR, h_WR_mass_TTEleCR);
  fillHisto(chain_TTMu_CR, &myEvent_TTMu_CR, h_WR_mass_TTMuCR);
  
  h_WR_mass_DataEleCR->Sumw2();
  h_WR_mass_DataMuCR->Sumw2();
  h_WR_mass_DYEleCR->Sumw2();
  h_WR_mass_DYMuCR->Sumw2();
  h_WR_mass_DYEleSR->Sumw2();
  h_WR_mass_DYMuSR->Sumw2();
  h_WR_mass_TTEleCR->Sumw2();
  h_WR_mass_TTMuCR->Sumw2();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///now fill the histos which show the ratio of two or more quantities
  
  //use this ratio histo to extrapolate DY MC from the signal region to control region
  TH1 *h_ratio_WR_mass_DYEleSRtoCR = (TH1*) h_WR_mass_DYEleSR->Clone();	//DY MC in SR div by DY MC in control region
  h_ratio_WR_mass_DYEleSRtoCR->Divide( ((TH1*) h_WR_mass_DYEleCR->Clone()) );
  
  //if any bin entries in the ratio are negative, set them to zero
  Int_t nbins = h_ratio_WR_mass_DYEleSRtoCR->GetNbinsX();//nbins is the same for all histos
  for(Int_t i=1; i<=nbins ;i++){
	  Double_t binContent = h_ratio_WR_mass_DYEleSRtoCR->GetBinContent(i);
	  if(binContent < 0.0) h_ratio_WR_mass_DYEleSRtoCR->SetBinContent(i, 0.0);
  }
  h_ratio_WR_mass_DYEleSRtoCR->Sumw2();

  //use this ratio histo to extrapolate DY MC from the signal region to control region
  TH1 *h_ratio_WR_mass_DYMuSRtoCR = (TH1*) h_WR_mass_DYMuSR->Clone();	//DY MC in SR div by DY MC in control region
  h_ratio_WR_mass_DYMuSRtoCR->Divide( ((TH1*) h_WR_mass_DYMuCR->Clone()) );
  
  //if any bin entries in the ratio are negative, set them to zero
  for(Int_t i=1; i<=nbins ;i++){
	  Double_t binContent = h_ratio_WR_mass_DYMuSRtoCR->GetBinContent(i);
	  if(binContent < 0.0) h_ratio_WR_mass_DYMuSRtoCR->SetBinContent(i, 0.0);
  }
  h_ratio_WR_mass_DYMuSRtoCR->Sumw2();


  //use this ratio histo to account for the difference in magnitude between DY and TTBar in the sideband (data represents the sum of them)
  //in each bin DY is between 90 and 100 percent of the total background
  TH1 *h_ratio_WR_mass_EleCR_DYtoDYplusTT = (TH1*) h_WR_mass_DYEleCR->Clone();	//DY MC in CR div by DY + TT MC in CR
  
  TH1 *h_ratio_WR_mass_EleCR_DYplusTT = (TH1*) h_WR_mass_DYEleCR->Clone();
  h_ratio_WR_mass_EleCR_DYplusTT->Add( ((TH1*) h_WR_mass_TTEleCR->Clone()) );
  
  h_ratio_WR_mass_EleCR_DYtoDYplusTT->Divide(h_ratio_WR_mass_EleCR_DYplusTT);
  h_ratio_WR_mass_EleCR_DYtoDYplusTT->Sumw2();

  //use this ratio histo to account for the difference in magnitude between DY and TTBar in the sideband (data represents the sum of them)
  //in each bin DY is between 90 and 100 percent of the total background
  TH1 *h_ratio_WR_mass_MuCR_DYtoDYplusTT = (TH1*) h_WR_mass_DYMuCR->Clone();	//DY MC in CR div by DY + TT MC in CR
   
  TH1 *h_ratio_WR_mass_MuCR_DYplusTT = (TH1*) h_WR_mass_DYMuCR->Clone();
  h_ratio_WR_mass_MuCR_DYplusTT->Add( ((TH1*) h_WR_mass_TTMuCR->Clone()) );
  
  h_ratio_WR_mass_MuCR_DYtoDYplusTT->Divide(h_ratio_WR_mass_MuCR_DYplusTT);
  h_ratio_WR_mass_MuCR_DYtoDYplusTT->Sumw2();


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///now multiply data in the control region CR by the ratio hist representing DY MC in SR over DY MC in CR, and
  ///the ratio hist representing the fraction of DY in the total background estimated by data in the CR
  ///to determine the DY background in the SR using collision data in the CR
  
  TH1* h_WR_mass_DataDrivenDYEleSR = (TH1*) h_WR_mass_DataEleCR->Clone();
  h_WR_mass_DataDrivenDYEleSR->Multiply( h_ratio_WR_mass_EleCR_DYtoDYplusTT );	///<this accounts for DY not constituting the entire bkgnd in CR, there is some TTBar
  h_WR_mass_DataDrivenDYEleSR->Multiply( h_ratio_WR_mass_DYEleSRtoCR );	///<this applies the DY SR to CR scale factors to data
  h_WR_mass_DataDrivenDYEleSR->Sumw2();

  TH1* h_WR_mass_DataDrivenDYMuSR = (TH1*) h_WR_mass_DataMuCR->Clone();
  h_WR_mass_DataDrivenDYMuSR->Multiply( h_ratio_WR_mass_MuCR_DYtoDYplusTT );	///<this accounts for DY not constituting the entire bkgnd in CR, there is some TTBar
  h_WR_mass_DataDrivenDYMuSR->Multiply( h_ratio_WR_mass_DYMuSRtoCR );	///<this applies the DY SR to CR scale factors to data
  h_WR_mass_DataDrivenDYMuSR->Sumw2();


  //normalize bin contents to bin widths
  for(Int_t i=1; i<=nbins; i++){
	  Double_t oldBinContent = h_WR_mass_DataDrivenDYMuSR->GetBinContent(i);
	  Double_t oldBinErrors = h_WR_mass_DataDrivenDYMuSR->GetBinError(i);
	  Double_t binWidth = h_WR_mass_DataDrivenDYMuSR->GetBinWidth(i);
	  h_WR_mass_DataDrivenDYMuSR->SetBinContent(i, oldBinContent/binWidth);
	  h_WR_mass_DataDrivenDYMuSR->SetBinError(i, oldBinErrors/binWidth);

	  oldBinContent = h_WR_mass_DYMuSR->GetBinContent(i);
	  oldBinErrors = h_WR_mass_DYMuSR->GetBinError(i);
	  binWidth = h_WR_mass_DYMuSR->GetBinWidth(i);
	  h_WR_mass_DYMuSR->SetBinContent(i, oldBinContent/binWidth);
	  h_WR_mass_DYMuSR->SetBinError(i, oldBinErrors/binWidth);

	  oldBinContent = h_WR_mass_DataDrivenDYEleSR->GetBinContent(i);
	  oldBinErrors = h_WR_mass_DataDrivenDYEleSR->GetBinError(i);
	  binWidth = h_WR_mass_DataDrivenDYEleSR->GetBinWidth(i);
	  h_WR_mass_DataDrivenDYEleSR->SetBinContent(i, oldBinContent/binWidth);
	  h_WR_mass_DataDrivenDYEleSR->SetBinError(i, oldBinErrors/binWidth);

	  oldBinContent = h_WR_mass_DYEleSR->GetBinContent(i);
	  oldBinErrors = h_WR_mass_DYEleSR->GetBinError(i);
	  binWidth = h_WR_mass_DYEleSR->GetBinWidth(i);
	  h_WR_mass_DYEleSR->SetBinContent(i, oldBinContent/binWidth);
	  h_WR_mass_DYEleSR->SetBinError(i, oldBinErrors/binWidth);

  }
  /*
  ///check bin contents
  std::cout<<"electron channel"<<std::endl;
  for(Int_t i=1; i<=nbins; i++){
	  std::cout<<"bin # "<< i<<"\t"<< h_WR_mass_DataDrivenDYEleSR->GetBinContent(i)<<std::endl;

  }
  std::cout<<"\t"<<std::endl;
  std::cout<<"muon channel"<<std::endl;
  for(Int_t i=1; i<=nbins; i++){
	  std::cout<<"bin # "<< i<<"\t"<< h_WR_mass_DataDrivenDYMuSR->GetBinContent(i)<<std::endl;
  }
  */

  gStyle->SetOptStat("");
  ///plot scale factors which interpolate between DY MC in the sideband to DY MC in the signal region
  makePlot(h_ratio_WR_mass_DYMuSRtoCR, "M_{MuMuJJ} [GeV]", "DY M_{MuMuJJ} Signal Region / Control Region", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}", "DYMCMu_SRtoCR_scaleFactors", "DYMuSRtoCR", -1, -2);
  makePlot(h_ratio_WR_mass_DYEleSRtoCR, "M_{EEJJ} [GeV]", "DY M_{EEJJ} Signal Region / Control Region", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}", "DYMCEle_SRtoCR_scaleFactors", "DYEleSRtoCR", -1, -2);

  ///plot scale factors which reflect the fraction of the total MC background in the sideband which comes from DY
  makePlot(h_ratio_WR_mass_MuCR_DYtoDYplusTT, "M_{MuMuJJ} [GeV]", "Control Region DY / DY + TT", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}", "MC_CRMu_DYoverDYplusTT_scaleFactors", "MuCRDYoverDYplusTT", 1.2, 0.6);
  makePlot(h_ratio_WR_mass_EleCR_DYtoDYplusTT, "M_{EEJJ} [GeV]", "Control Region DY / DY + TT", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}", "MC_CREle_DYoverDYplusTT_scaleFactors", "EleCRDYoverDYplusTT", 1.2, 0.6);

  ///plot the net scale factors which are applied to data in the sideband to estimate DY in the signal region
  TH1* MuNetSF = (TH1*) h_ratio_WR_mass_DYMuSRtoCR->Clone();
  MuNetSF->Multiply(h_ratio_WR_mass_MuCR_DYtoDYplusTT);
  makePlot(MuNetSF, "M_{MuMuJJ} [GeV]", "Scale Factor applied to Mu Data", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}", "MuData_scaleFactors", "SFappliedToMuData", 0.2, -0.1);

  TH1* EleNetSF = (TH1*) h_ratio_WR_mass_DYEleSRtoCR->Clone();
  EleNetSF->Multiply(h_ratio_WR_mass_EleCR_DYtoDYplusTT);
  makePlot(EleNetSF, "M_{EEJJ} [GeV]", "Scale Factor applied to Ele Data", "CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}","EleData_scaleFactors","SFappliedToEleData", 0.3, -0.1);


  ///plot DY SR estimate from data overlaid on DY MC in SR
  TCanvas* mycanv_DYMu = new TCanvas( "mycanv_DYMu", "", 0, 0, 600, 600 ) ;
  h_WR_mass_DataDrivenDYMuSR->GetXaxis()->SetTitle("M_{MuMuJJ} [GeV]");
  h_WR_mass_DataDrivenDYMuSR->SetMarkerSize(1.5);
  h_WR_mass_DataDrivenDYMuSR->SetMarkerStyle(21);
  h_WR_mass_DYMuSR->SetTitle("CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}");
  h_WR_mass_DYMuSR->GetYaxis()->SetTitle("Events normalized to bin width");
  h_WR_mass_DYMuSR->GetXaxis()->SetTitle("M_{MuMuJJ} [GeV]");
  h_WR_mass_DYMuSR->SetLineColor(kRed);
  h_WR_mass_DYMuSR->SetLineWidth(3);
  h_WR_mass_DYMuSR->SetMaximum(1.6);
  h_WR_mass_DYMuSR->Draw("histo");
  h_WR_mass_DataDrivenDYMuSR->Draw("epsame");
  
  TLegend *leg_MuMu = new TLegend( 0.6, 0.60, 0.9, 0.90 );
  leg_MuMu->AddEntry( h_WR_mass_DataDrivenDYMuSR, "Data Driven" );
  leg_MuMu->AddEntry( h_WR_mass_DYMuSR, "MC" );
  leg_MuMu->Draw();
  mycanv_DYMu->Print(("DYMuSR_data_driven_and_MC_estimate.pdf"));
  mycanv_DYMu->Print(("DYMuSR_data_driven_and_MC_estimate.png"));

  TCanvas* mycanv_DYEle = new TCanvas( "mycanv_DYEle", "", 0, 0, 600, 600 ) ;
  h_WR_mass_DataDrivenDYEleSR->GetXaxis()->SetTitle("M_{MuMuJJ} [GeV]");
  h_WR_mass_DataDrivenDYEleSR->SetMarkerSize(1.5);
  h_WR_mass_DataDrivenDYEleSR->SetMarkerStyle(21);
  h_WR_mass_DYEleSR->SetTitle("CMS Private   #surds = 13 TeV  #int lumi = 2.6 fb^{-1}");
  h_WR_mass_DYEleSR->GetXaxis()->SetTitle("M_{EEJJ} [GeV]");
  h_WR_mass_DYEleSR->GetYaxis()->SetTitle("Events normalized to bin width");
  h_WR_mass_DYEleSR->SetLineColor(kRed);
  h_WR_mass_DYEleSR->SetLineWidth(3);
  h_WR_mass_DYEleSR->SetMaximum(1.6);
  h_WR_mass_DYEleSR->Draw("histo");
  h_WR_mass_DataDrivenDYEleSR->Draw("epsame");
  
  TLegend *leg_EE = new TLegend( 0.6, 0.60, 0.9, 0.90 );
  leg_EE->AddEntry( h_WR_mass_DataDrivenDYEleSR, "Data Driven" );
  leg_EE->AddEntry( h_WR_mass_DYEleSR, "MC" );
  leg_EE->Draw();
  mycanv_DYEle->Print(("DYEleSR_data_driven_and_MC_estimate.pdf"));
  mycanv_DYEle->Print(("DYEleSR_data_driven_and_MC_estimate.png"));


}//end altDYEstimate()

void fillHisto(TChain * chain, Selector *myEvent, TH1 * h){

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
    if(myEvent->WR_mass > 600.) 
      h->Fill(myEvent->WR_mass,myEvent->weight);
  }
  std::cout<<"histo named\t"<< h->GetName() <<"\thas integral\t"<< h->Integral() <<std::endl;
}//end fillHisto()


void makePlot(TH1 *histo, TString xAxisLabel, TString yAxisLabel, TString plotTitle, std::string outputPlotFileName, TString canvName, Float_t maxYval, Float_t minYval){
  	TCanvas* canv = new TCanvas( canvName, canvName, 0, 0, 600, 600 );
	histo->GetXaxis()->SetTitle(xAxisLabel);
	histo->SetMarkerSize(1.5);
	histo->SetMarkerStyle(21);
	histo->SetTitle(plotTitle);
	histo->GetYaxis()->SetTitle(yAxisLabel);
	histo->SetTitleOffset(1.3, "y");
	if(maxYval > 0.) histo->SetMaximum(maxYval);
	if(minYval > -1.) histo->SetMinimum(minYval);
	histo->Draw("ep");

	canv->Print( (outputPlotFileName + ".pdf").c_str() );
	canv->Print( (outputPlotFileName + ".png").c_str() );

}//end makePlot()

/*
//call this fxn once TF1 is fitted to distribution
std::string chiSquaredNdofString(TF1 * fit){
  std::string tempchiSqd = "#chi^{2}  =  ";
  std::string chiSqdVal = to_string(fit->GetChisquare());
  std::string ndof = to_string(fit->GetNDF());
  tempchiSqd += chiSqdVal.substr(0,4);
  tempchiSqd += " / ";
  tempchiSqd += ndof.substr(0,2);
  return tempchiSqd;
}
*/

