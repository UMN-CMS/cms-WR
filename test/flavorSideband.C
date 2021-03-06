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

std::string chiSquaredNdofString(TF1 * fit);
void fillHisto(TChain * chain, Selector *myEvent, TH1F * h);
void flavorSideband(){

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
}

void fillHisto(TChain * chain, Selector *myEvent, TH1F * h){

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
    if(myEvent->WR_mass > 600. && myEvent->dilepton_mass > 200.) 
      h->Fill(myEvent->WR_mass,myEvent->weight);
  }
  std::cout<<"histo named\t"<< h->GetName() <<"\thas integral\t"<< h->Integral() <<std::endl;
}

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

