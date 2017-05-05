#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include "../src/Selector.cc"
#include "../src/SelectorHist.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>

#define useDYMAD
//#define doDataDrivenTT
#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

/**
 * this macro is designed to read several TChains, representing collision data and bkgnd MC, apply cuts, and plot
 * data as black dots and all bkgnd MC as one stacked histogram with several fill colors.
 *
 * This macro should be used by itself on minitrees processed by analysis.cpp with -c EE or -c MuMu
 * and signal region requirements (i.e. only ignoreDyScaleFactors false passed to analysis.cpp on the cmd line).
 *
 */

//switch Selector tag here, and everything else will change accordingly
Selector::tag_t channel = Selector::EE;

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void quickLowFourObjMassMiniPlotter(){

  TChain * chain_DY = new TChain("Tree_Iter0","DY");
  TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
#ifdef doDataDrivenTT
  chain_ttbar = new TChain("Tree_Iter0","TTData");
#endif 
  TChain * chain_WJets = new TChain("Tree_Iter0","WJets");
  TChain * chain_WZ = new TChain("Tree_Iter0","Diboson");
  TChain * chain_ZZ = new TChain("Tree_Iter0","Other");
  TChain * chain_data = new TChain("Tree_Iter0","Data");

  TString localDir = "../analysisCppOutputRootFilesNewHltSf/";
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0;
  switch (channel) {
  case Selector::EE:
#ifndef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");
#endif

#ifdef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYMadInclAndHT_signal_eeEE_withMllWeight.root");
#endif

#ifndef doDataDrivenTT
	tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_eeEE.root");
#endif
#ifdef doDataDrivenTT
	tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
#endif
	wjets = chain_WJets->Add(localDir+"selected_tree_WInclAndHT_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
  	zz = chain_ZZ->Add(localDir+"selected_tree_topW_signal_eeEE.root");
    data = chain_data->Add(localDir+"selected_tree_data_signal_eeEE.root");
    break;
  case Selector::MuMu:
#ifndef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");
#endif

#ifdef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYMadInclAndHT_signal_mumuMuMu_withMllWeight.root");
#endif

#ifndef doDataDrivenTT
	tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_mumuMuMu.root");
#endif
#ifdef doDataDrivenTT
   	tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuMuMu.root");
#endif
	wjets = chain_WJets->Add(localDir+"selected_tree_WInclAndHT_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");
	zz = chain_ZZ->Add(localDir+"selected_tree_topW_signal_mumuMuMu.root");
    data = chain_data->Add(localDir+"selected_tree_data_signal_mumuMuMu.root");
    break;
  default:
    std::cout << "Unknown tag" << std::endl;
  }

  std::cout<<"data = "<< data <<"\tdy = "<< dy << std::endl;
  if(data==0 || dy==0 || tt==0 || wjets==0 || wz==0 || zz==0) exit(-1);

  Selector myEvent_DY;
  Selector myEvent_ttbar;
  Selector myEvent_WJets;
  Selector myEvent_WZ;
  Selector myEvent_ZZ;
  Selector myEvent_data;

  myEvent_DY.SetBranchAddresses(chain_DY);
  myEvent_ttbar.SetBranchAddresses(chain_ttbar);
  myEvent_WJets.SetBranchAddresses(chain_WJets);
  myEvent_WZ.SetBranchAddresses(chain_WZ);
  myEvent_ZZ.SetBranchAddresses(chain_ZZ);
  myEvent_data.SetBranchAddresses(chain_data);

  std::vector<TH1F*> hs_DY;
  MakeHistos(chain_DY, &myEvent_DY, &hs_DY);
  std::vector<TH1F*> hs_ttbar;
  MakeHistos(chain_ttbar, &myEvent_ttbar, &hs_ttbar);
  std::vector<TH1F*> hs_WJets;
  MakeHistos(chain_WJets, &myEvent_WJets, &hs_WJets);
  std::vector<TH1F*> hs_WZ;
  MakeHistos(chain_WZ, &myEvent_WZ, &hs_WZ);
  std::vector<TH1F*> hs_ZZ;
  MakeHistos(chain_ZZ, &myEvent_ZZ, &hs_ZZ);

  std::vector<TH1F*> hs_data;
  MakeHistos(chain_data, &myEvent_data, &hs_data);

  unsigned int nPlots = hs_DY.size();

  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WJets[i],hs_WZ[i],hs_ZZ[i],hs_data[i],xtitles[i],fnames[i]);
  }
  
}//end quickLowFourObjMassMiniPlotter()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",50,0,700);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",50,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",50,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",50,0,700);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",50,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",50,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",50,0,700);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",50,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",50,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",50,0,700);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",50,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",50,-3.15,3.15);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  //TH1F *h_WR_mass = new TH1F("h_WR_mass","",30,150,800);
  //TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",30,150,550);

  /*variable bin widths  dont forget to enable the bin division by its width in the lines below*/
  Float_t bins[] = { 150, 300, 350, 375, 400, 425, 450, 475, 500, 525, 550, 600, 800};
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1F *h_WR_mass = new TH1F("h_WR_mass_varBins","",binnum, bins);

  Float_t mllBins[] = { 155, 200, 210, 225, 240, 255, 275, 300, 330, 375, 550};
  Int_t  mllBinnum = sizeof(mllBins)/sizeof(Float_t) - 1;
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass_varBins","",mllBinnum, mllBins);
  /**/


  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  Float_t ttScaleFactor = 1.0;
  TString chainTitle(chain->GetTitle());
  /*not used previously, but left here for reference
   * if( chainTitle.EqualTo("TTMC") ){
	  ttScaleFactor = (channel == Selector::MuMu) ? 0.958 : 0.954;	//to account for slightly higher number of rescaled ttbar MC events relative to rescaled emu data evts
  }*/
  if( chainTitle.EqualTo("TTData") ){
	  ttScaleFactor = (channel == Selector::MuMu) ? 0.655 : 0.416;	//to rescale emu data evts to estimates of ttbar in electron and muon channels
  }

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	
	if(myEvent->WR_mass > 600.) continue;	//skip evts with four obj mass above 600 GeV
	if(myEvent->dilepton_mass < 200.) continue; //skip evts with dilepton mass below 200 GeV
	if(myEvent->sublead_lepton_pt < 53.) continue;
	
	h_lepton_pt0->Fill(myEvent->lead_lepton_pt,(myEvent->weight)*ttScaleFactor);
    h_lepton_pt1->Fill(myEvent->sublead_lepton_pt,(myEvent->weight)*ttScaleFactor);
    h_lepton_eta0->Fill(myEvent->lead_lepton_eta,(myEvent->weight)*ttScaleFactor);
    h_lepton_eta1->Fill(myEvent->sublead_lepton_eta,(myEvent->weight)*ttScaleFactor);
    h_lepton_phi0->Fill(myEvent->lead_lepton_phi,(myEvent->weight)*ttScaleFactor);
    h_lepton_phi1->Fill(myEvent->sublead_lepton_phi,(myEvent->weight)*ttScaleFactor);

    h_jet_pt0->Fill(myEvent->lead_jet_pt,(myEvent->weight)*ttScaleFactor);
    h_jet_pt1->Fill(myEvent->sublead_jet_pt,(myEvent->weight)*ttScaleFactor);
    h_jet_eta0->Fill(myEvent->lead_jet_eta,(myEvent->weight)*ttScaleFactor);
    h_jet_eta1->Fill(myEvent->sublead_jet_eta,(myEvent->weight)*ttScaleFactor);
    h_jet_phi0->Fill(myEvent->lead_jet_phi,(myEvent->weight)*ttScaleFactor);
    h_jet_phi1->Fill(myEvent->sublead_jet_phi,(myEvent->weight)*ttScaleFactor);
      
    h_WR_mass->Fill(myEvent->WR_mass,(myEvent->weight)*ttScaleFactor);
    h_dilepton_mass->Fill(myEvent->dilepton_mass,(myEvent->weight)*ttScaleFactor);
    h_nPV->Fill(myEvent->nPV,(myEvent->weight)*ttScaleFactor);
  }

  hs->push_back(h_lepton_pt0);
  hs->push_back(h_lepton_pt1);
  hs->push_back(h_jet_pt0);
  hs->push_back(h_jet_pt1);
  hs->push_back(h_lepton_eta0);
  hs->push_back(h_lepton_eta1);
  hs->push_back(h_jet_eta0);
  hs->push_back(h_jet_eta1);
  hs->push_back(h_lepton_phi0);
  hs->push_back(h_lepton_phi1);
  hs->push_back(h_jet_phi0);
  hs->push_back(h_jet_phi1);
  hs->push_back(h_WR_mass);
  hs->push_back(h_dilepton_mass);
  hs->push_back(h_nPV);

  /*shouldn't do this for all histos, just the variable bin width MLLJJ and MLL plots*/
  //normalize histo bins when using variable bin widths
  unsigned int max = hs->size();
  for(unsigned int i=0; i<max; i++){
	  TString histName = (hs->at(i))->GetName();
	  if(histName.Contains("varBins")){
		  //get num bins in histo i
		  Int_t nBins = (hs->at(i))->GetNbinsX();
		  for(Int_t j=1; j<=nBins; j++){
			  if(j==nBins){
				  Double_t origBinContents = (hs->at(i))->GetBinContent(j);
				  Double_t overflowContents = (hs->at(i))->GetBinContent(j+1);
				  std::cout<<"overflow contents =\t"<< overflowContents << std::endl;	//sanity check
				  (hs->at(i))->SetBinContent(j, origBinContents+overflowContents);
			  }//end work to include overflows in last bin shown on plot

			  //in each bin, divide the bin contents by the bin width
			  Double_t oldBinContents = (hs->at(i))->GetBinContent(j);
			  Double_t oldBinErrors = (hs->at(i))->GetBinError(j);
			  Double_t binWidth = (hs->at(i))->GetBinWidth(j);
			  (hs->at(i))->SetBinContent(j, oldBinContents/binWidth);
			  (hs->at(i))->SetBinError(j, oldBinErrors/binWidth);
		  }//end loop over bins in histo

	  }//end filter which selects histos which have variable bin widths

  }//end loop over histos in vector
  /**/

}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  TString ttbarLegEntry = "TTBar MC";
#ifdef doDataDrivenTT
  ttbarLegEntry = "TT Data Driven";
#endif
  TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.80 ) ; 
  leg->AddEntry( hs_ttbar, ttbarLegEntry ) ;
 
#ifndef useDYMAD
  leg->AddEntry( hs_DY, "DY AMCNLO" ) ;
#endif
#ifdef useDYMAD
  leg->AddEntry( hs_DY, "Z/#gamma*+jets" ) ;
#endif 
  leg->AddEntry( hs_ZZ, "TopW MC" ) ; 
  leg->AddEntry( hs_WZ, "Diboson" ) ; 
  leg->AddEntry( hs_WJets, "W+jets" ) ; 
  leg->AddEntry( hs_data, "Data");
  leg->SetFillColor( kWhite ) ; 

  hs_data->Sumw2();
  hs_ttbar->Sumw2();
  hs_WJets->Sumw2();
  hs_WZ->Sumw2();
  hs_ZZ->Sumw2();
  hs_DY->Sumw2();
  
  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 900, 900 ) ;
  //mycanvas->cd();	//only needed when no ratio plot is drawn
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WJets->SetFillColor(kBlue);
  hs_WZ->SetFillColor(kCyan);
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_WZ);
  th->Add(hs_WJets);
  th->Add(hs_ZZ);
  th->Add(hs_DY);
  th->Add(hs_ttbar);
  hs_data->SetMarkerStyle(20);


  Float_t labelSize = 0.25;
  /*for ratio plot*/
  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  /**/
  hs_data->SetStats(0);

  //reset vertical scale in MLLJJ and MLL plots to harmonize plot scales between lepton channels
  hs_data->SetMinimum(0.007);
  th->SetMinimum(0.007);


  TH1F *ratio = (TH1F*)hs_data->Clone();
  //th->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  //hs_data->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  th->SetTitle("CMS Private       2.6 fb^{-1} (13 TeV)");
  hs_data->SetTitle("CMS Private       2.6 fb^{-1} (13 TeV)");
  hs_data->Draw("ep");
  th->Draw("histo same");
  hs_data->Draw("epsame");

  TString ytitle = "Events/(";
  ytitle += (th->GetXaxis()->GetNbins());
  ytitle += ")";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  hs_data->GetYaxis()->SetTitle(ytitle.Data());

  /*if variable bin widths are used*/
  if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetYaxis()->SetTitle("Events/GeV"), hs_data->GetYaxis()->SetTitle("Events/GeV   ");
  else if(fname.EqualTo("Mll")) hs_data->GetXaxis()->SetTitle("M_{LL} [GeV]"), th->GetXaxis()->SetTitle("M_{LL} [GeV]"), th->GetYaxis()->SetTitle("Events/GeV"), hs_data->GetYaxis()->SetTitle("Events/GeV   ");
   /**/
  hs_data->Draw("epsame");
  mycanvas->Update();
 
  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")) ratio->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  else if(fname.EqualTo("Mll")) ratio->GetXaxis()->SetTitle("M_{LL} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(labelSize+0.03);
  ratio->SetLabelSize(labelSize - 0.07,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();	//for ratio plot
  ratio->Sumw2();
  ratio->SetStats(0);

  hs_ttbar->Add(hs_WJets);
  hs_ttbar->Add(hs_WZ);
  hs_ttbar->Add(hs_ZZ);
  hs_ttbar->Add(hs_DY);

  ratio->Divide(hs_ttbar);
  ratio->SetMarkerStyle(21);
  ratio->GetYaxis()->SetTitle("data/MC   ");
  ratio->GetYaxis()->SetTitleSize(0.19);
  ratio->GetYaxis()->SetTitleOffset(0.2);
  //std::cout<<"y axis ratio title size=\t" << ratio->GetYaxis()->GetTitleSize()<<std::endl;
  //std::cout<<"y axis ratio title offset=\t" << ratio->GetYaxis()->GetTitleOffset()<<std::endl;
  ratio->SetLabelSize(labelSize - 0.07,"y");
  Float_t maxYratioRange = 1.98;
  //if(channel == Selector::EE) maxYratioRange = 2.4;
  ratio->GetYaxis()->SetRangeUser(0.5,maxYratioRange);
  ratio->GetYaxis()->SetNdivisions(505);
  /*for ratio plot*/
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  mycanvas->Update();
  /**/

  TString fn = fname;
#ifndef useDYMAD
  if(channel == Selector::EE) fn += "_EEChannelDataAndBkgndMC_DYAMC_VariableBinWidthWithRatio_lowFourObjMass";
  if(channel == Selector::MuMu) fn += "_MuMuChannelDataAndBkgndMC_DYAMC_VariableBinWidthWithRatio_lowFourObjMass";
#endif

#ifdef useDYMAD
  if(channel == Selector::EE) fn += "_EEChannelDataAndBkgndMC_DYMadHTAndIncl_VariableBinWidthWithRatio_lowFourObjMass";
  if(channel == Selector::MuMu) fn += "_MuMuChannelDataAndBkgndMC_DYMadHTAndIncl_VariableBinWidthWithRatio_lowFourObjMass";
#endif

  TString ttbarMode = "_TTBarMC";
#ifdef doDataDrivenTT
  ttbarMode = "_TTBarFromData";
#endif
  fn += ttbarMode;

  if(fname.EqualTo("Mlljj") || fname.EqualTo("Mll")){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  mycanvas->Print((fn+".C").Data());
	  p1->SetLogy();	//for ratio plot
	  //mycanvas->SetLogy();	//only needed when no ratio plot is drawn
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
	  mycanvas->Print((fn+"_log.C").Data());
  }

  mycanvas->Close();
}
