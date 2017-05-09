#include "TStyle.h"
#include "TPaveText.h"
#include "TH1F.h"
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
//#include "tdrStyle.C"
#include <cstdio>
#include <memory>

/**
 * this macro runs on EMu data and TTBar MC minitrees which have been processed by analysis.cpp
 * with -c EMu or -c EE or -c MuMu and signal region requirements.  This macro can be modified, and
 * should be run by itself.  Currently it is not used in wrValidation.sh.
 *
 * This macro creates plots comparing ttbar MC to itself and to EMu data.  No plots use stacked histograms,
 * only overlaid curves.  It also calculates and saves the ttbar MC EE/EMu and MuMu/EMu scale factors to a file in data/2015-v1/.
 */

//#define PRINTRATIOS
#define PRINT

std::string chiSquaredNdofString(TF1 * fit);
void fillHisto(TChain * chain, Selector *myEvent, TH1F * h);
void flavorSideband(){

  Float_t mumuEmuSF = 0.6592, eeEmuSF = 0.4315;	///<used for plotting and chi^2 calculation

  std::string longMuMuEmuSF = to_string(mumuEmuSF);
  std::string shortMuMuEmuSF = longMuMuEmuSF.substr(0,5);
  std::string longEEEmuSF = to_string(eeEmuSF);
  std::string shortEEEmuSF = longEEEmuSF.substr(0,5);
  gStyle->SetOptFit(0);	///<show nothing
  TChain * chain_EMu = new TChain("Tree_Iter0");
  TChain * chain_TopW_EMu = new TChain("Tree_Iter0");
  TChain * chain_EE = new TChain("Tree_Iter0");
  TChain * chain_MuMu = new TChain("Tree_Iter0");
  TChain * chain_EMuData = new TChain("Tree_Iter0");
 
  TString dir = "../analysisCppOutputRootFilesNewHltSf/";
  chain_EMu->Add(dir+"selected_tree_TT_flavoursidebandEMuEE.root");
  chain_EMu->Add(dir+"selected_tree_topW_flavoursidebandEMuEE.root");
  chain_EE->Add(dir+"selected_tree_TT_signal_eeEE.root");
  chain_EE->Add(dir+"selected_tree_topW_signal_eeEE.root");
  chain_MuMu->Add(dir+"selected_tree_TT_signal_mumuMuMu.root");
  chain_MuMu->Add(dir+"selected_tree_topW_signal_mumuMuMu.root");
  chain_EMuData->Add(dir+"selected_tree_data_flavoursidebandEMuEE.root");
  chain_TopW_EMu->Add(dir+"selected_tree_topW_flavoursidebandEMuEE.root");
 
  Selector myEvent_EMu;
  Selector myEvent_EMuData;
  Selector myEvent_TopW_EMu;
  Selector myEvent_EE;
  Selector myEvent_MuMu;

  myEvent_EMu.SetBranchAddresses(chain_EMu);
  myEvent_TopW_EMu.SetBranchAddresses(chain_TopW_EMu);
  myEvent_EMuData.SetBranchAddresses(chain_EMuData);
  myEvent_EE.SetBranchAddresses(chain_EE);
  myEvent_MuMu.SetBranchAddresses(chain_MuMu);

  Float_t bins[] = { 410, 450, 500, 550, 600, 625, 652, 683, 718, 760, 812, 877, 975, 1160, 1999 };	//for xMax of 2000
  //Float_t bins[] = { 200, 400, 450, 500, 550, 600, 625, 652, 683, 718, 760, 812, 877, 975, 1160, 2000, 6000 };	//for xMax of 6000
  //Float_t bins[] = { 1400, 1600, 1800, 2000, 2200, 2400, 2600};
  
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;

  ////fixed bin width MLLJJ plots with standard domain
  //TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",12,400,6000);
  //TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",12,400,6000);
  //TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",12,400,6000);
  //TH1F *h_WR_mass_EMuData = new TH1F("h_WR_mass_EMuData","",12,400,6000);
  
  ////variable bin width MLLJJ plots
  TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",binnum, bins);
  TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",binnum, bins);
  TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",binnum, bins);
  TH1F *h_WR_mass_EMuData = new TH1F("h_WR_mass_EMuData","",binnum, bins);

  fillHisto(chain_EMu, &myEvent_EMu, h_WR_mass_EMu);
  fillHisto(chain_EMuData, &myEvent_EMuData, h_WR_mass_EMuData);
  fillHisto(chain_EE, &myEvent_EE, h_WR_mass_EE);
  fillHisto(chain_MuMu, &myEvent_MuMu, h_WR_mass_MuMu);

 
  ///use this title for all plots
  TString stdTitle = "CMS Preliminary               2.6 fb^{-1} (13 TeV)";
  h_WR_mass_EMu->SetTitle(stdTitle);
  h_WR_mass_EMuData->SetTitle(stdTitle);
  h_WR_mass_EE->SetTitle(stdTitle);
  h_WR_mass_MuMu->SetTitle(stdTitle);


  /**/
  TCanvas* mycanvas_EE = new TCanvas( "mycanvas_EE", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_EE->SetLineColor(kRed);
  h_WR_mass_EE->DrawNormalized("same histo");
  TLegend *leg_EE = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_EE->AddEntry( h_WR_mass_EMu, "TTBar+TopW EMu" );
  leg_EE->AddEntry( h_WR_mass_EE, "TTBar+TopW EE" );
  leg_EE->Draw();
  mycanvas_EE->Print(("flavor_EE_fixedbinwidth.pdf"));
  mycanvas_EE->Print(("flavor_EE_fixedbinwidth.png"));

  TCanvas* mycanvas_MuMu = new TCanvas( "mycanvas_MuMu", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_MuMu->SetLineColor(kRed);
  h_WR_mass_MuMu->DrawNormalized("same histo");
  TLegend *leg_MuMu = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_MuMu->AddEntry( h_WR_mass_EMu, "TTBar+TopW EMu" );
  leg_MuMu->AddEntry( h_WR_mass_MuMu, "TTBar+TopW MuMu" );
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
  h_ratio_EE->GetYaxis()->SetTitle("M_{EEJJ} / M_{EMuJJ}    ");
  h_ratio_EE->SetTitleOffset(1.5,"Y");
  h_ratio_EE->SetTitle(stdTitle);
  h_ratio_EE->SetFillColor(kWhite);
  h_ratio_MuMu->Divide(h_WR_mass_EMu);
  h_ratio_MuMu->SetTitle(stdTitle);
  h_ratio_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  h_ratio_MuMu->GetYaxis()->SetRangeUser(0.51,0.79);
  h_ratio_MuMu->GetYaxis()->SetTitle("M_{MuMuJJ} / M_{EMuJJ}    ");
  h_ratio_MuMu->SetTitleOffset(1.5,"Y");
  h_ratio_MuMu->SetFillColor(kWhite);

  //TDR style formatting
  gStyle->SetOptTitle(1);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(.15);
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.1, "Y");
  gStyle->SetTitleOffset(1.0, "Z");
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(1);
  gStyle->SetHistLineWidth(3);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  TCanvas* mycanvas_ratio_EE = new TCanvas( "mycanvas_ratio_EE", "", 0, 0, 600, 600 ) ;
  //TPaveText* chiSqdBoxEE = new TPaveText(1500.,0.54,2000.,0.58);
  TPaveText* chiSqdBoxEE = new TPaveText(475.,0.54,1975.,0.585);	///< for xmax 2000
  chiSqdBoxEE->SetFillColorAlpha(kWhite, 1.0);
  chiSqdBoxEE->SetTextFont(42);
  //chiSqdBoxEE->SetTextSize(0.05);
  TF1 *f_EE = new TF1("f_EE","[0]",400,1999);
  f_EE->FixParameter(0,eeEmuSF);
  h_ratio_EE->Fit("f_EE");
  chiSqdBoxEE->AddText( "Simulation" );
  chiSqdBoxEE->AddText( TString( chiSquaredNdofString(f_EE) ) );
  chiSqdBoxEE->AddText( TString("SF = " + shortEEEmuSF) );
  chiSqdBoxEE->SetTextColor(1);
  h_ratio_EE->SetTitleFont(42);
  //h_ratio_EE->SetTitleColor(1);
  //h_ratio_EE->SetTitleTextColor(1);
  //h_ratio_EE->SetTitleFillColor(10);
  //h_ratio_EE->SetTitleFontSize(0.05);
  h_ratio_EE->Draw("ep");
  chiSqdBoxEE->Draw("same");
  f_EE->SetLineColor(kBlack);
  f_EE->Draw("same");
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth.png"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_largeXmax.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_largeXmax.png"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth.pdf"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth.png"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth.C"));
  mycanvas_ratio_EE->SetLogx(1);
  //chiSqdBoxEE->DrawPave(500.,0.54,1000.,0.58,4,"same");
  //chiSqdBoxEE->DrawPave(300.,0.5,1100.,0.54,4,"same");	//for xmax 2000
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth_logx_largeXmax.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_fixedbinwidth_logx_largeXmax.png"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx_largeXmax.pdf"));
  //mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx_largeXmax.png"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx.pdf"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE_variablebinwidth_logx.png"));


  TCanvas* mycanvas_ratio_MuMu = new TCanvas( "mycanvas_ratio_MuMu", "", 0, 0, 600, 600 ) ;
  //TPaveText* chiSqdBoxMuMu = new TPaveText(1500.,0.73,2000.,0.79);
  TPaveText* chiSqdBoxMuMu = new TPaveText(475.,0.74,1975.,0.785);	///< for xmax 2000
  chiSqdBoxMuMu->SetFillColorAlpha(kWhite, 1.0);
  chiSqdBoxMuMu->SetTextFont(42);
  TF1 *f_MuMu = new TF1("f_MuMu","[0]",400,1999);
  f_MuMu->FixParameter(0,mumuEmuSF);
  h_ratio_MuMu->Fit("f_MuMu");
  chiSqdBoxMuMu->AddText( "Simulation" );
  chiSqdBoxMuMu->AddText( TString( chiSquaredNdofString(f_MuMu) ) );
  chiSqdBoxMuMu->AddText( TString("SF = " + shortMuMuEmuSF) );
  chiSqdBoxMuMu->SetTextColor(1);
  h_ratio_MuMu->Draw();
  f_MuMu->SetLineColor(kBlack);
  chiSqdBoxMuMu->Draw("same");
  f_MuMu->Draw("same");
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth.png"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_largeXmax.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_largeXmax.png"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth.pdf"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth.png"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth.C"));
  mycanvas_ratio_MuMu->SetLogx(1);
  //chiSqdBoxMuMu->DrawPave(300.,0.74,1100.,0.79,4,"same");
  //chiSqdBoxMuMu->DrawPave(300.,0.70,1100.,0.74,4,"same");	//for xmax 2000
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth_logx.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_fixedbinwidth_logx.png"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx_largeXmax.pdf"));
  //mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx_largeXmax.png"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx.pdf"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu_variablebinwidth_logx.png"));


  gStyle->SetOptStat("");
  TCanvas* canvMuMuEMu = new TCanvas("canvMuMuEMu","",600,600);
  canvMuMuEMu->cd();
  TLegend * legMuMuEMu = new TLegend(0.72,0.6,0.98,0.8);
  legMuMuEMu->AddEntry(h_WR_mass_EMu,"TTBar+TopW EMu MC");
  legMuMuEMu->AddEntry(h_WR_mass_MuMu,"TTBar+TopW MuMu MC");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_MuMu->Draw("Psame");
  legMuMuEMu->Draw();
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region_fixedbinwidth.pdf","recreate");
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region_fixedbinwidth.png","recreate");
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region_fixedbinwidth.C","recreate");
 
  TCanvas* canvMuMuEMuData = new TCanvas("canvMuMuEMuData","",600,600);
  canvMuMuEMuData->cd();
  h_WR_mass_EMuData->SetMarkerStyle(20);
  h_WR_mass_EMuData->SetMarkerSize(1);
  h_WR_mass_EMuData->SetMarkerColor(kBlack);
  h_WR_mass_MuMu->SetLineWidth(3);
  h_WR_mass_MuMu->SetLineColor(kBlack);
  h_WR_mass_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  //h_WR_mass_MuMu->GetYaxis()->SetTitle("Events");	///fixed bin widths
  h_WR_mass_MuMu->GetYaxis()->SetTitle("Events/GeV");	///variable bin widths
  TLegend * legMuMuEMuData = new TLegend(0.72,0.65,0.98,0.9);
  legMuMuEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legMuMuEMuData->AddEntry(h_WR_mass_MuMu,"TTBar+TopW MuMu MC");
  h_WR_mass_EMuData->Scale(mumuEmuSF);
 
  //for ratio plot
  Double_t epsMuMu = 0.001;
  TPad* p1MuMu = new TPad("p1MuMu", "p1MuMu", 0, 0.25, 1, 1, 0);
  p1MuMu->Draw();
  TPad* p2MuMu = new TPad("p2MuMu", "p2MuMu", 0, 0.1, 1, 0.25 + epsMuMu, 0);
  p2MuMu->Draw();
  p1MuMu->SetBottomMargin(0);
  p2MuMu->SetTopMargin(0);
  p1MuMu->cd();
  h_WR_mass_EMuData->SetStats(1);
  h_WR_mass_MuMu->SetStats(1);
  TH1F *ratio_RescaledData_MuMu = (TH1F*) h_WR_mass_EMuData->Clone();
  h_WR_mass_MuMu->Draw("histo");
  h_WR_mass_EMuData->Draw("Psame");
  legMuMuEMuData->Draw();
  
  canvMuMuEMuData->cd();
  p2MuMu->cd();	///<change to ratio TPad
  //Float_t minYratio = 0.41, maxYratio = 1.99;	///fixed bin width
  Float_t minYratio = 0.61, maxYratio = 1.35;	///variable bin width
 
  Float_t xLabelSize = 0.25, xTitleSize = 0.25;
  ratio_RescaledData_MuMu->SetTitle("");
  ratio_RescaledData_MuMu->GetXaxis()->SetTitleSize(xTitleSize);
  ratio_RescaledData_MuMu->SetLabelSize(xLabelSize,"x");
  ratio_RescaledData_MuMu->Sumw2();
  ratio_RescaledData_MuMu->SetStats(0);

  ratio_RescaledData_MuMu->Divide(h_WR_mass_MuMu);
#ifdef PRINT
  Int_t nBinsMuMu = ratio_RescaledData_MuMu->GetNbinsX();
  for(Int_t n=1; n<=nBinsMuMu ; n++){
	  std::cout<< "in ratio of rescaled emujj data to MuMu MC, bin # " << n << " has ratio=\t" << ratio_RescaledData_MuMu->GetBinContent(n) << std::endl;
  }//end loop over bins in ratio plot
#endif 
  ratio_RescaledData_MuMu->SetMarkerStyle(20);
  ratio_RescaledData_MuMu->SetMarkerColor(kBlack);
  ratio_RescaledData_MuMu->SetLabelSize(0.25, "y");
  ratio_RescaledData_MuMu->GetYaxis()->SetRangeUser(minYratio, maxYratio);
  ratio_RescaledData_MuMu->GetYaxis()->SetNdivisions(505);
  ratio_RescaledData_MuMu->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");

  ratio_RescaledData_MuMu->Draw("p");
  float xmax = ratio_RescaledData_MuMu->GetXaxis()->GetXmax();
  float xmin = ratio_RescaledData_MuMu->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1", xmin, xmax);
  f1->Draw("same");
  canvMuMuEMuData->cd();
  canvMuMuEMuData->Update();
 
  //canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_fixedbinwidth.pdf","recreate");
  //canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_fixedbinwidth.png","recreate");
  //canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_fixedbinwidth.C","recreate");
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_variablebinwidth_overflowsInLastBin.pdf","recreate");
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_variablebinwidth_overflowsInLastBin.png","recreate");
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region_variablebinwidth_overflowsInLastBin.C","recreate");


  TCanvas* canvEEEMu = new TCanvas("canvEEEMu","",600,600);
  canvEEEMu->cd();
  TLegend * legEEEMu = new TLegend(0.72,0.6,0.98,0.8);
  legEEEMu->AddEntry(h_WR_mass_EMu,"EMu");
  legEEEMu->AddEntry(h_WR_mass_EE,"EE");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_EE->Draw("Psame");
  legEEEMu->Draw();
  //canvEEEMu->SaveAs("emujj_and_eejj_signal_region_fixedbinwidth.pdf","recreate");
  //canvEEEMu->SaveAs("emujj_and_eejj_signal_region_fixedbinwidth.png","recreate");


  TCanvas* canvEEEMuData = new TCanvas("canvEEEMuData","",600,600);
  canvEEEMuData->cd();
  h_WR_mass_EE->SetLineColor(kBlack);
  h_WR_mass_EE->SetLineWidth(3);
  h_WR_mass_EE->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  //h_WR_mass_EE->GetYaxis()->SetTitle("Events");	///fixed bin widths
  h_WR_mass_EE->GetYaxis()->SetTitle("Events/GeV");	///variable bin widths
  h_WR_mass_EMuData->Scale(1/(mumuEmuSF));	///<undo the scaling which was done earlier
  h_WR_mass_EMuData->Scale(eeEmuSF);
  TLegend * legEEEMuData = new TLegend(0.72,0.65,0.98,0.9);
  legEEEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legEEEMuData->AddEntry(h_WR_mass_EE,"TTBar+TopW EE MC");
  
  //for ratio plot
  Double_t epsEE = 0.001;
  TPad* p1EE = new TPad("p1EE", "p1EE", 0, 0.25, 1, 1, 0);
  p1EE->Draw();
  TPad* p2EE = new TPad("p2EE", "p2EE", 0, 0.1, 1, 0.25 + epsEE, 0);
  p2EE->Draw();
  p1EE->SetBottomMargin(0);
  p2EE->SetTopMargin(0);
  p1EE->cd();
  h_WR_mass_EMuData->SetStats(1);
  h_WR_mass_EE->SetStats(1);
  TH1F *ratio_RescaledData_EE = (TH1F*) h_WR_mass_EMuData->Clone();
  
  h_WR_mass_EE->Draw("histo");
  h_WR_mass_EMuData->Draw("Psame");
  legEEEMuData->Draw();
 
  canvEEEMuData->cd();
  p2EE->cd();	///<change to ratio TPad
  ratio_RescaledData_EE->SetTitle("");
  ratio_RescaledData_EE->GetXaxis()->SetTitleSize(xTitleSize);
  ratio_RescaledData_EE->SetLabelSize(xLabelSize,"x");
  ratio_RescaledData_EE->Sumw2();
  ratio_RescaledData_EE->SetStats(0);

  ratio_RescaledData_EE->Divide(h_WR_mass_EE);
#ifdef PRINT
  Int_t nBinsEE = ratio_RescaledData_EE->GetNbinsX();
  for(Int_t n=1; n<=nBinsEE ; n++){
	  std::cout<< "in ratio of rescaled emujj data to EE MC, bin # " << n << " has ratio=\t" << ratio_RescaledData_EE->GetBinContent(n) << std::endl;
  }//end loop over bins in ratio plot
#endif 
  ratio_RescaledData_EE->SetMarkerStyle(20);
  ratio_RescaledData_EE->SetMarkerColor(kBlack);
  ratio_RescaledData_EE->SetLabelSize(0.25, "y");
  ratio_RescaledData_EE->GetYaxis()->SetRangeUser(minYratio, maxYratio);
  ratio_RescaledData_EE->GetYaxis()->SetNdivisions(505);
  ratio_RescaledData_EE->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");

  ratio_RescaledData_EE->Draw("p");
  f1->Draw("same");
  canvEEEMuData->cd();
  canvEEEMuData->Update();
  //canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_fixedbinwidth.pdf","recreate");
  //canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_fixedbinwidth.png","recreate");
  //canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_fixedbinwidth.C","recreate");
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_variablebinwidth_overflowsInLastBin.pdf","recreate");
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_variablebinwidth_overflowsInLastBin.png","recreate");
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region_variablebinwidth_overflowsInLastBin.C","recreate");


  ///rescale the EEJJ and MuMuJJ MC and compare it to the EMu data and EMuJJ MC
  ///this plot with variable bin widths SHOULD NOT be used. If it is ever decided to use this plot with
  ///variable bin widths, the bin contents and bin errors must be divided by the bin widths.
  gStyle->SetOptStat("");
  TCanvas* canvEMuDataTwoRescaledMCs = new TCanvas("canvEMuDataTwoRescaledMCs","",600,600);
  canvEMuDataTwoRescaledMCs->cd();
  TLegend * legEMuDataTwoRescaledMCs = new TLegend(0.5,0.45,0.9,0.9);
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EMuData,"EMu Data");	
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EE,"Rescaled TTBar+TopW EE MC");
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_MuMu,"Rescaled TTBar+TopW MuMu MC");
  legEMuDataTwoRescaledMCs->AddEntry(h_WR_mass_EMu,"TTBar+TopW EMu MC");
  legEMuDataTwoRescaledMCs->SetTextSize(0.027);
  h_WR_mass_EE->Scale(1/(eeEmuSF));
  h_WR_mass_MuMu->Scale(1/(mumuEmuSF));
  h_WR_mass_EMuData->Scale(1/(eeEmuSF));	///<undo the scaling which was done earlier
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
  //TString ytitle = "Events";	//for fixed bin widths
  TString ytitle = "Events/GeV";	//for variable bin widths
  h_WR_mass_EMuData->GetYaxis()->SetTitle(ytitle.Data());
  h_WR_mass_EMuData->Draw("epsame");
  
  legEMuDataTwoRescaledMCs->Draw();
  canvEMuDataTwoRescaledMCs->cd();
  p2->cd();	///<change to ratio TPad
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
  ratio_Data_EE->Draw("p");
  ratio_Data_EMu->Draw("psame");
  ratio_Data_MuMu->Draw("psame");
  f1->Draw("same");
  canvEMuDataTwoRescaledMCs->cd();
  canvEMuDataTwoRescaledMCs->Update();
  //canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_fixedbinwidth.pdf","recreate");
  //canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_fixedbinwidth.png","recreate");
  //canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_fixedbinwidth.C","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_variablebinwidth.pdf","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_variablebinwidth.png","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_variablebinwidth.C","recreate");

  p1->SetLogy();
  //canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_fixedbinwidth.pdf","recreate");
  //canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_fixedbinwidth.png","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_variablebinwidth.pdf","recreate");
  canvEMuDataTwoRescaledMCs->SaveAs("emujj_data_and_MC_and_rescaled_eejj_and_mumujj_MC_signal_region_log_variablebinwidth.png","recreate");
  canvEMuDataTwoRescaledMCs->Close();

  /**/
}

void fillHisto(TChain * chain, Selector *myEvent, TH1F * h){

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	if(myEvent->sublead_lepton_pt < 53.) continue;
    if(myEvent->WR_mass > 600. && myEvent->dilepton_mass > 200.) 
      h->Fill(myEvent->WR_mass,myEvent->weight);
  }
  //Int_t lowBinOne = h->FindBin(1450.), lowBinTwo = h->FindBin(1700.), lowBinThree = h->FindBin(1850.), lowBinFour = h->FindBin(2150.);
  //Int_t highBinOne = h->FindBin(2000.), highBinTwo = h->FindBin(2200.), highBinThree = h->FindBin(2550.), highBinFour = h->FindBin(3050.);
  std::cout<<"\t"<<std::endl;
  std::cout<<"histo named\t"<< h->GetName() <<"\thas integral\t"<< h->Integral() <<std::endl;
  //std::cout<<"histo named\t"<< h->GetName() <<"\thas in 1.8 TeV window\t"<< h->Integral(lowBinOne, highBinOne) <<std::endl;
  //std::cout<<"histo named\t"<< h->GetName() <<"\thas in 2.0 TeV window\t"<< h->Integral(lowBinTwo, highBinTwo) <<std::endl;
  //std::cout<<"histo named\t"<< h->GetName() <<"\thas in 2.2 TeV window\t"<< h->Integral(lowBinThree, highBinThree) <<std::endl;
  //std::cout<<"histo named\t"<< h->GetName() <<"\thas in 2.4 TeV window\t"<< h->Integral(lowBinFour, highBinFour) <<std::endl;
  std::cout<<"\t"<<std::endl;


  /*for variable bin widths*/
  //rescale bin contents by bin widths
  Int_t nBins = h->GetNbinsX();
  for(Int_t j=1; j<=nBins; j++){
	  //include the overflows in the very last bin shown on the plot
	  if(j==nBins){
		  Double_t origBinContents = h->GetBinContent(j);
		  Double_t overflowContents = h->GetBinContent(j+1);
		  h->SetBinContent(j, origBinContents+overflowContents);
	  }//end work to include overflows in last bin shown on plot
	  //in each bin, divide the bin contents by the bin width
	  Double_t oldBinContents = h->GetBinContent(j);
	  Double_t oldBinErrors = h->GetBinError(j);
	  Double_t binWidth = h->GetBinWidth(j);
	  h->SetBinContent(j, oldBinContents/binWidth);
	  h->SetBinError(j, oldBinErrors/binWidth);
  }//end loop over bins in histo
  /**/
}

//call this fxn once TF1 is fitted to distribution
std::string chiSquaredNdofString(TF1 * fit){
  std::string tempchiSqd = "#chi^{2}  =  ";
  std::string chiSqdVal = to_string(fit->GetChisquare());
  std::string ndof = to_string(fit->GetNDF()-1);
  tempchiSqd += chiSqdVal.substr(0,4);
  tempchiSqd += " / ";
  tempchiSqd += ndof.substr(0,2);
  std::cout<< tempchiSqd << std::endl;
  return tempchiSqd;
}

