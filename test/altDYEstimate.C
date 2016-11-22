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
 * this macro 
 *
 */

//#define PRINTRATIOS

//std::string chiSquaredNdofString(TF1 * fit);
void fillHisto(TChain * chain, Selector *myEvent, TH1 * h);
void altDYEstimate(){

  TString dir = "../analysisCppOutputRootFiles/";
  TString treeName = "Tree_Iter0";

  //sideband region files
  TChain * chain_DataMu_CR = new TChain(treeName);
  chain_DataMu_CR->Add(dir+"selected_tree_data_lowdileptonsidebandMuMu.root");
  TChain * chain_DYMu_CR = new TChain(treeName);
  chain_DYMu_CR->Add(dir+"selected_tree_DYAMC_lowdileptonsidebandMuMu_withMllWeight.root");
  TChain * chain_DataEle_CR = new TChain(treeName);
  chain_DataEle_CR->Add(dir+"selected_tree_data_lowdileptonsidebandEE.root");
  TChain * chain_DYEle_CR = new TChain(treeName);
  chain_DYEle_CR->Add(dir+"selected_tree_DYAMC_lowdileptonsidebandEE_withMllWeight.root");
  TChain * chain_TTMu_CR = new TChain(treeName);
  chain_TTMu_CR->Add(dir+"selected_tree_TT_lowdileptonsidebandMuMu.root");
  TChain * chain_TTEle_CR = new TChain(treeName);
  chain_TTEle_CR->Add(dir+"selected_tree_TT_lowdileptonsidebandEE.root");

  //signal region files
  TChain * chain_DYMu_SR = new TChain(treeName);
  chain_DYMu_SR->Add(dir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");
  TChain * chain_DYEle_SR = new TChain(treeName);
  chain_DYEle_SR->Add(dir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");

 
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
  //Float_t bins[] = {600, 625, 652, 683, 718, 760, 812, 877, 975, 1160, 2000, 4000 };	//from emu sideband
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
  //fixed bin widths
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
  ///now fill ratio histos
  
  //use this ratio histo to extrapolate DY MC from the signal region to control region
  TH1 *h_ratio_WR_mass_DYEleSRtoCR = (TH1*) h_WR_mass_DYEleSR->Clone();	//DY MC in SR div by DY MC in control region
  h_ratio_WR_mass_DYEleSRtoCR->Divide( ((TH1*) h_WR_mass_DYEleCR->Clone()) );
  //Divide(h_ratio_WR_mass_DYEleSRtoCR, ((TH1*) h_WR_mass_DYEleCR->Clone()));
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
  //Divide(h_ratio_WR_mass_DYMuSRtoCR, ((TH1*) h_WR_mass_DYMuCR->Clone()));
  //if any bin entries in the ratio are negative, set them to zero
  for(Int_t i=1; i<=nbins ;i++){
	  Double_t binContent = h_ratio_WR_mass_DYMuSRtoCR->GetBinContent(i);
	  if(binContent < 0.0) h_ratio_WR_mass_DYMuSRtoCR->SetBinContent(i, 0.0);
  }
  h_ratio_WR_mass_DYMuSRtoCR->Sumw2();


  //use this ratio histo to account for the difference in magnitude between DY and TTBar (data represents the sum of them)
  TH1 *h_ratio_WR_mass_EleCR_DYtoDYplusTT = (TH1*) h_WR_mass_DYEleCR->Clone();	//DY MC in CR div by DY + TT MC in CR
  for(Int_t i=1; i<=nbins ; i++){
	  Double_t numer = h_ratio_WR_mass_EleCR_DYtoDYplusTT->GetBinContent(i);
	  Double_t denom = (h_WR_mass_DYEleCR->GetBinContent(i)) + (h_WR_mass_TTEleCR->GetBinContent(i));
	  h_ratio_WR_mass_EleCR_DYtoDYplusTT->SetBinContent(i, numer/denom);
  }
  h_ratio_WR_mass_EleCR_DYtoDYplusTT->Sumw2();

  //use this ratio histo to account for the difference in magnitude between DY and TTBar (data represents the sum of them)
  TH1 *h_ratio_WR_mass_MuCR_DYtoDYplusTT = (TH1*) h_WR_mass_DYMuCR->Clone();	//DY MC in CR div by DY + TT MC in CR
  for(Int_t i=1; i<=nbins; i++){
	  Double_t numer = h_ratio_WR_mass_MuCR_DYtoDYplusTT->GetBinContent(i);
	  Double_t denom = (h_WR_mass_DYMuCR->GetBinContent(i)) + (h_WR_mass_TTMuCR->GetBinContent(i));
	  h_ratio_WR_mass_MuCR_DYtoDYplusTT->SetBinContent(i, numer/denom);
  }
  h_ratio_WR_mass_MuCR_DYtoDYplusTT->Sumw2();


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///now multiply data in the control region CR by the ratio hist representing DY MC in SR over DY MC in CR, and
  ///the ratio hist representing the fraction of DY in the total background estimated in the CR
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
  
  ///plot DY SR estimate from data overlaid on DY MC in SR
  gStyle->SetOptStat("");
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


  /*
  Float_t mumuEmuSF = 0.657, eeEmuSF = 0.414;	///<used for plotting and chi^2 calculation

  std::string longMuMuEmuSF = to_string(mumuEmuSF);
  std::string shortMuMuEmuSF = longMuMuEmuSF.substr(0,4);
  std::string longEEEmuSF = to_string(eeEmuSF);
  std::string shortEEEmuSF = longEEEmuSF.substr(0,4);
  //gStyle->SetOptFit(100);	///<only show chiSquared divided by nDOF
  gStyle->SetOptFit(0);	///<show nothing
  TChain * chain_EMu = new TChain("Tree_Iter0");
  TChain * chain_EE = new TChain("Tree_Iter0");
  TChain * chain_MuMu = new TChain("Tree_Iter0");
  TChain * chain_EMuData = new TChain("Tree_Iter0");
 
  TString dir = "../analysisCppOutputRootFiles/";
  chain_EMu->Add(dir+"selected_tree_TT_flavoursidebandEMu.root");
  chain_EE->Add(dir+"selected_tree_TT_signal_eeEE.root");
  chain_MuMu->Add(dir+"selected_tree_TT_signal_mumuMuMu.root");
  chain_EMuData->Add(dir+"selected_tree_data_flavoursidebandEMu.root");
 
  Selector myEvent_EMu;
  Selector myEvent_EMuData;
  Selector myEvent_EE;
  Selector myEvent_MuMu;

  myEvent_EMu.SetBranchAddresses(chain_EMu);
  myEvent_EMuData.SetBranchAddresses(chain_EMuData);
  myEvent_EE.SetBranchAddresses(chain_EE);
  myEvent_MuMu.SetBranchAddresses(chain_MuMu);

  Float_t bins[] = { 410, 450, 500, 550, 600, 625, 652, 683, 718, 760, 812, 877, 975, 1160, 2000 };	//for xMax of 2000
  //Float_t bins[] = { 200, 400, 450, 500, 550, 600, 625, 652, 683, 718, 760, 812, 877, 975, 1160, 2000, 6000 };	//for xMax of 6000
  
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;

  //fixed bin width MLLJJ plots with standard domain
  TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",21,530,2000);
  TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",21,530,2000);
  TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",21,530,2000);
  TH1F *h_WR_mass_EMuData = new TH1F("h_WR_mass_EMuData","",21,530,2000);
  
  //variable bin width MLLJJ plots
  //TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",binnum, bins);
  //TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",binnum, bins);
  //TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",binnum, bins);
  //TH1F *h_WR_mass_EMuData = new TH1F("h_WR_mass_EMuData","",binnum, bins);
  
  fillHisto(chain_EMu, &myEvent_EMu, h_WR_mass_EMu);
  fillHisto(chain_EMuData, &myEvent_EMuData, h_WR_mass_EMuData);
  fillHisto(chain_EE, &myEvent_EE, h_WR_mass_EE);
  fillHisto(chain_MuMu, &myEvent_MuMu, h_WR_mass_MuMu);
  
  Double_t MuMuMCIntegral = h_WR_mass_MuMu->Integral();
  Double_t EEMCIntegral = h_WR_mass_EE->Integral();
  Double_t EMuMCIntegral = h_WR_mass_EMu->Integral();

  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"ttbarMC EE/EMu ratio =\t"<< (EEMCIntegral/EMuMCIntegral) << std::endl;
  std::cout<<"ttbarMC MuMu/EMu ratio =\t"<< (MuMuMCIntegral/EMuMCIntegral) << std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::string ttScaleFactorFile = "../data/2015-v1/ttScaleFactors.txt";
  ofstream writeToTTFile(ttScaleFactorFile.c_str(), ofstream::trunc);
  writeToTTFile << "#Channel\tSF" << std::endl;
  writeToTTFile << "EE\t" << (EEMCIntegral/EMuMCIntegral) << std::endl;
  writeToTTFile << "MuMu\t" << (MuMuMCIntegral/EMuMCIntegral) << std::endl;

  ///use this title for all plots
  TString stdTitle = "CMS Private           #surds = 13 TeV #int lumi = 2.6 fb^{-1}";
  h_WR_mass_EMu->SetTitle(stdTitle);
  h_WR_mass_EMuData->SetTitle(stdTitle);
  h_WR_mass_EE->SetTitle(stdTitle);
  h_WR_mass_MuMu->SetTitle(stdTitle);


  TCanvas* mycanvas_EE = new TCanvas( "mycanvas_EE", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_EE->SetLineColor(kRed);
  h_WR_mass_EE->DrawNormalized("same histo");
  TLegend *leg_EE = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_EE->AddEntry( h_WR_mass_EMu, "EMu" );
  leg_EE->AddEntry( h_WR_mass_EE, "EE" );
  leg_EE->Draw();
  mycanvas_EE->Print(("flavor_EE_fixedbinwidth.pdf"));
  mycanvas_EE->Print(("flavor_EE_fixedbinwidth.png"));

  TCanvas* mycanvas_MuMu = new TCanvas( "mycanvas_MuMu", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_MuMu->SetLineColor(kRed);
  h_WR_mass_MuMu->DrawNormalized("same histo");
  TLegend *leg_MuMu = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_MuMu->AddEntry( h_WR_mass_EMu, "EMu" );
  leg_MuMu->AddEntry( h_WR_mass_MuMu, "MuMu" );
  leg_MuMu->Draw();
  mycanvas_MuMu->Print(("flavor_MuMu_fixedbinwidth.pdf"));
  mycanvas_MuMu->Print(("flavor_MuMu_fixedbinwidth.png"));

#ifdef PRINTRATIOS
  std::cout<<"number of overflow events in h_WR_mass_EE histo =\t"<< h_WR_mass_EE->GetBinContent(h_WR_mass_EE->GetNbinsX() + 1) <<std::endl;
  std::cout<<"number of overflow events in h_WR_mass_MuMu histo =\t"<< h_WR_mass_MuMu->GetBinContent(h_WR_mass_MuMu->GetNbinsX() + 1) <<std::endl;
  std::cout<<"number of overflow events in h_WR_mass_EMu histo =\t"<< h_WR_mass_EMu->GetBinContent(h_WR_mass_EMu->GetNbinsX() + 1) <<std::endl;
#endif
  ///dont show anything in histo stats box for the two plots with TF1 curves
  gStyle->SetOptStat("");
  TH1F *h_ratio_EE = (TH1F*)h_WR_mass_EE->Clone();
  TH1F *h_ratio_MuMu = (TH1F*)h_WR_mass_MuMu->Clone();
  h_ratio_EE->Divide(h_WR_mass_EMu);
  h_ratio_EE->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_ratio_EE->GetYaxis()->SetRangeUser(0.31,0.59);
  h_ratio_EE->GetYaxis()->SetTitle("ratio M_{EEJJ} / M_{EMuJJ}");
  h_ratio_EE->SetTitleOffset(1.55,"Y");
  h_ratio_EE->SetTitle(stdTitle);
  h_ratio_EE->SetFillColor(kWhite);
  h_ratio_MuMu->Divide(h_WR_mass_EMu);
  h_ratio_MuMu->SetTitle(stdTitle);
  h_ratio_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_ratio_MuMu->GetYaxis()->SetRangeUser(0.51,0.79);
  h_ratio_MuMu->GetYaxis()->SetTitle("ratio M_{MuMuJJ} / M_{EMuJJ}");
  h_ratio_MuMu->SetTitleOffset(1.55,"Y");
  h_ratio_MuMu->SetFillColor(kWhite);
  
  TCanvas* mycanvas_ratio_EE = new TCanvas( "mycanvas_ratio_EE", "", 0, 0, 600, 600 ) ;
  //TPaveText* chiSqdBoxEE = new TPaveText(1500.,0.54,2000.,0.58);
  TPaveText* chiSqdBoxEE = new TPaveText(475.,0.54,1975.,0.585);	///< for xmax 2000
  chiSqdBoxEE->SetFillColorAlpha(kWhite, 1.0);
  TF1 *f_EE = new TF1("f_EE","[0]",600,1500);
  f_EE->FixParameter(0,eeEmuSF);
  h_ratio_EE->Fit("f_EE");
  chiSqdBoxEE->AddText( TString( chiSquaredNdofString(f_EE) ) );
  chiSqdBoxEE->AddText( TString("ratio = " + shortEEEmuSF) );
  h_ratio_EE->Draw();
  chiSqdBoxEE->Draw("same");
  f_EE->SetLineColor(kBlue);
  f_EE->Draw("same");
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth.png"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_largeXmax.pdf"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_largeXmax.png"));
  mycanvas_ratio_EE->SetLogx(1);
  //chiSqdBoxEE->DrawPave(500.,0.54,1000.,0.58,4,"same");
  //chiSqdBoxEE->DrawPave(300.,0.5,1100.,0.54,4,"same");	//for xmax 2000
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth_logx_largeXmax.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth_logx_largeXmax.png"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx_largeXmax.pdf"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx_largeXmax.png"));


  TCanvas* mycanvas_ratio_MuMu = new TCanvas( "mycanvas_ratio_MuMu", "", 0, 0, 600, 600 ) ;
  //TPaveText* chiSqdBoxMuMu = new TPaveText(1500.,0.73,2000.,0.79);
  TPaveText* chiSqdBoxMuMu = new TPaveText(475.,0.74,1975.,0.785);	///< for xmax 2000
  chiSqdBoxMuMu->SetFillColorAlpha(kWhite, 1.0);
  TF1 *f_MuMu = new TF1("f_MuMu","[0]",600,1500);
  f_MuMu->FixParameter(0,mumuEmuSF);
  h_ratio_MuMu->Fit("f_MuMu");
  chiSqdBoxMuMu->AddText( TString( chiSquaredNdofString(f_MuMu) ) );
  chiSqdBoxMuMu->AddText( TString("ratio = " + shortMuMuEmuSF) );
  h_ratio_MuMu->Draw();
  f_MuMu->SetLineColor(kBlue);
  chiSqdBoxMuMu->Draw("same");
  f_MuMu->Draw("same");
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth.png"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_largeXmax.pdf"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_largeXmax.png"));
  mycanvas_ratio_MuMu->SetLogx(1);
  //chiSqdBoxMuMu->DrawPave(300.,0.74,1100.,0.79,4,"same");
  //chiSqdBoxMuMu->DrawPave(300.,0.70,1100.,0.74,4,"same");	//for xmax 2000
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth_logx.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth_logx.png"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx_largeXmax.pdf"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx_largeXmax.png"));


  gStyle->SetOptStat("nemr");
  TCanvas* canvMuMuEMu = new TCanvas("canvMuMuEMu","",600,600);
  canvMuMuEMu->cd();
  TLegend * legMuMuEMu = new TLegend(0.72,0.6,0.98,0.8);
  legMuMuEMu->AddEntry(h_WR_mass_EMu,"EMu");
  legMuMuEMu->AddEntry(h_WR_mass_MuMu,"MuMu");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_MuMu->Draw("Psame");
  legMuMuEMu->Draw();
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region_fixedbinwidth.pdf","recreate");
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region_fixedbinwidth.png","recreate");

  TCanvas* canvMuMuEMuData = new TCanvas("canvMuMuEMuData","",600,600);
  canvMuMuEMuData->cd();
  TLegend * legMuMuEMuData = new TLegend(0.72,0.6,0.98,0.8);
  legMuMuEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legMuMuEMuData->AddEntry(h_WR_mass_MuMu,"MuMu MC");
  h_WR_mass_MuMu->Draw("histo");
  h_WR_mass_EMuData->Scale(mumuEmuSF);
  h_WR_mass_EMuData->SetMarkerStyle(2);
  h_WR_mass_EMuData->SetMarkerSize(2);
  h_WR_mass_EMuData->Draw("Psame");
  legMuMuEMuData->Draw();
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_fixedbinwidth.pdf","recreate");
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_fixedbinwidth.png","recreate");


  TCanvas* canvEEEMu = new TCanvas("canvEEEMu","",600,600);
  canvEEEMu->cd();
  TLegend * legEEEMu = new TLegend(0.72,0.6,0.98,0.8);
  legEEEMu->AddEntry(h_WR_mass_EMu,"EMu");
  legEEEMu->AddEntry(h_WR_mass_EE,"EE");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_EE->Draw("Psame");
  legEEEMu->Draw();
  canvEEEMu->SaveAs("emujj_and_eejj_signal_region_fixedbinwidth.pdf","recreate");
  canvEEEMu->SaveAs("emujj_and_eejj_signal_region_fixedbinwidth.png","recreate");

  TCanvas* canvEEEMuData = new TCanvas("canvEEEMuData","",600,600);
  canvEEEMuData->cd();
  TLegend * legEEEMuData = new TLegend(0.72,0.6,0.98,0.8);
  legEEEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legEEEMuData->AddEntry(h_WR_mass_EE,"EE MC");
  h_WR_mass_EE->Draw("histo");
  h_WR_mass_EMuData->Scale(1/mumuEmuSF);	///<undo the scaling which was done earlier
  h_WR_mass_EMuData->SetMarkerStyle(2);
  h_WR_mass_EMuData->SetMarkerSize(2);
  h_WR_mass_EMuData->Scale(eeEmuSF);
  h_WR_mass_EMuData->Draw("Psame");
  legEEEMuData->Draw();
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_fixedbinwidth.pdf","recreate");
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_fixedbinwidth.png","recreate");

  //TFile f("flavor_fits.root","RECREATE");
  //h_ratio_EE->Write();
  //h_ratio_MuMu->Write();
  //f_EE->Write();
  //f_MuMu->Write();


  ///rescale the EEJJ and MuMuJJ MC and compare it to the EMu data and EMuJJ MC
  ///this plot with variable bin widths SHOULD NOT be used. If it is ever decided to use this plot with
  ///variable bin widths, the bin contents and bin errors must be divided by the bin widths.
  gStyle->SetOptStat("");
  TCanvas* canvEMuDataTwoRescaledMCs = new TCanvas("canvEMuDataTwoRescaledMCs","",600,600);
  canvEMuDataTwoRescaledMCs->cd();
  //TLegend * legEMuDataTwoRescaledMCs = new TLegend(0.6,0.15,0.94,0.45);	//for variable bin widths
  TLegend * legEMuDataTwoRescaledMCs = new TLegend(0.6,0.55,0.94,0.9);	//for fixed bin widths
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EMuData,"EMu Data");	
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EE,"Rescaled TTBar EE MC");
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_MuMu,"Rescaled TTBar MuMu MC");
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EMu,"TTBar EMu MC");
  legEMuDataTwoRescaledMCs->SetTextSize(0.027);
  h_WR_mass_EE->Scale(1/eeEmuSF);
  h_WR_mass_MuMu->Scale(1/mumuEmuSF);
  h_WR_mass_EMuData->Scale(1/eeEmuSF);	///<undo the scaling which was done earlier
  h_WR_mass_EE->SetLineColor(kRed);
  h_WR_mass_EE->SetLineWidth(3);
  h_WR_mass_EMu->SetLineColor(kBlue);
  h_WR_mass_EMu->SetLineWidth(3);
  h_WR_mass_MuMu->SetLineColor(kBlack);
  h_WR_mass_MuMu->SetLineWidth(3);
  h_WR_mass_EMuData->SetMarkerStyle(20);
  h_WR_mass_EMuData->SetMarkerSize(1);
  h_WR_mass_EMuData->SetMarkerColor(kBlack);
  h_WR_mass_EMuData->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EE->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");


  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1", "p1", 0, 0.25, 1, 1, 0);
  p1->Draw();
  TPad* p2 = new TPad("p2", "p2", 0, 0.1, 1, 0.25 + eps, 0);
  p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);
  p1->cd();
  h_WR_mass_EMuData->SetStats(1);
  h_WR_mass_MuMu->SetStats(1);
  TH1F *ratio_Data_EE = (TH1F*) h_WR_mass_EMuData->Clone();
  TH1F *ratio_Data_EMu = (TH1F*) h_WR_mass_EMuData->Clone();
  TH1F *ratio_Data_MuMu = (TH1F*) h_WR_mass_EMuData->Clone();
  h_WR_mass_EMuData->Draw("ep");
  h_WR_mass_EMu->Draw("histo same");
  h_WR_mass_EE->Draw("histo same");
  h_WR_mass_MuMu->Draw("histo same");
  h_WR_mass_EMuData->Draw("epsame");
  //TString ytitle = "Events/(";
  //ytitle += (h_WR_mass_EMuData->GetXaxis()->GetNbins());
  //ytitle += ")";
  TString ytitle = "Events";
  h_WR_mass_EMuData->GetYaxis()->SetTitle(ytitle.Data());
  h_WR_mass_EMuData->Draw("epsame");
  
  legEMuDataTwoRescaledMCs->Draw();
  canvEMuDataTwoRescaledMCs->cd();
  p2->cd();	///<change to ratio TPad
  Float_t minYratio = 0.41, maxYratio = 1.99;
  Float_t xLabelSize = 0.25, xTitleSize = 0.25;
  ratio_Data_EE->SetTitle(""), ratio_Data_MuMu->SetTitle(""), ratio_Data_EMu->SetTitle("");
  ratio_Data_EE->GetXaxis()->SetTitleSize(xTitleSize), ratio_Data_MuMu->GetXaxis()->SetTitleSize(xTitleSize), ratio_Data_EMu->GetXaxis()->SetTitleSize(xTitleSize);
  ratio_Data_EE->SetLabelSize(xLabelSize,"x"), ratio_Data_MuMu->SetLabelSize(xLabelSize,"x"), ratio_Data_EMu->SetLabelSize(xLabelSize,"x");
  ratio_Data_EE->Sumw2();
  ratio_Data_EE->SetStats(0);
  ratio_Data_EMu->Sumw2();
  ratio_Data_EMu->SetStats(0);
  ratio_Data_MuMu->Sumw2();
  ratio_Data_MuMu->SetStats(0);
  ratio_Data_EE->Divide(h_WR_mass_EE);
  ratio_Data_EE->SetMarkerStyle(20);
  ratio_Data_EE->SetMarkerColor(kRed);
  ratio_Data_EE->SetLabelSize(0.25, "y");
  ratio_Data_EE->GetYaxis()->SetRangeUser(minYratio, maxYratio);
  ratio_Data_EE->GetYaxis()->SetNdivisions(505);
  ratio_Data_EE->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  ratio_Data_EE->SetLabelSize(0.25,"x");

  ratio_Data_EMu->Divide(h_WR_mass_EMu);
  ratio_Data_EMu->SetMarkerStyle(22);
  ratio_Data_EMu->SetMarkerColor(kBlue);
  ratio_Data_EMu->SetLabelSize(0.25, "y");
  ratio_Data_EMu->GetYaxis()->SetRangeUser(minYratio, maxYratio);
  ratio_Data_EMu->GetYaxis()->SetNdivisions(505);
  ratio_Data_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");

  ratio_Data_MuMu->Divide(h_WR_mass_MuMu);
  ratio_Data_MuMu->SetMarkerStyle(21);
  ratio_Data_MuMu->SetMarkerColor(kBlack);
  ratio_Data_MuMu->SetLabelSize(0.25, "y");
  ratio_Data_MuMu->GetYaxis()->SetRangeUser(minYratio, maxYratio);
  ratio_Data_MuMu->GetYaxis()->SetNdivisions(505);
  ratio_Data_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");

  ratio_Data_EE->Draw("p");
  ratio_Data_EMu->Draw("p");
  ratio_Data_MuMu->Draw("p");
  float xmax = ratio_Data_MuMu->GetXaxis()->GetXmax();
  float xmin = ratio_Data_MuMu->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1", xmin, xmax);
  ratio_Data_EE->Draw("p");
  ratio_Data_EMu->Draw("psame");
  ratio_Data_MuMu->Draw("psame");
  f1->Draw("same");
  canvEMuDataTwoRescaledMCs->cd();
  canvEMuDataTwoRescaledMCs->Update();
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_fixedbinwidth.pdf","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_fixedbinwidth.png","recreate");

  p1->SetLogy();
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_fixedbinwidth.pdf","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_fixedbinwidth.png","recreate");
  canvEMuDataTwoRescaledMCs->Close();
  */

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
}

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

