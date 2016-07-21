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
// #include "ExoAnalysis/cmsWR/src/Selector.cc"
// #include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "../src/Selector.cc"
#include "../src/SelectorHist.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>


#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

Selector::tag_t channel = Selector::EMu;

/*
 * this macro is designed to read several TChains, representing data and MC, apply no cuts, and plot
 * data as points and all MC as one stacked histogram with several fill colors
 *
 */

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void quickEMuMiniPlotter(){

  TChain * chain_DY = new TChain("Tree_Iter0","DY");
  TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
  TChain * chain_WJets = new TChain("Tree_Iter0","WJets");
  TChain * chain_WZ = new TChain("Tree_Iter0","WZ");
  TChain * chain_ZZ = new TChain("Tree_Iter0","ZZ");
  TChain * chain_data = new TChain("Tree_Iter0","Data");
 
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0;
  switch (channel) {
  case Selector::EMu:
    //dy = chain_DY->Add("../selected_tree_DYAMC_lowdileptonsidebandEMu_withMllWeight.root");
	dy = chain_DY->Add("../selected_tree_DYAMC_signal_eeEE.root");	//temporary work around, DY will not be plotted
	tt = chain_ttbar->Add("../selected_tree_TT_flavoursidebandEMu.root");
    wjets = chain_WJets->Add("../selected_tree_W_flavoursidebandEMu.root");
    wz = chain_WZ->Add("../selected_tree_WZ_flavoursidebandEMu.root");
    zz = chain_ZZ->Add("../selected_tree_ZZ_flavoursidebandEMu.root");
    data = chain_data->Add("../selected_tree_data_flavoursidebandEMu.root");
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
  
}//end quickEMuMiniPlotter()

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

  Float_t bins[] = { 150, 200, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800};	//wider bins work better at high WR mass
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  //TH1F *h_WR_mass = new TH1F("h_WR_mass","",50,0,2500);
  TH1F *h_WR_mass = new TH1F("h_WR_mass","",binnum, bins);
 
  float dilepton_max = 250.;
  if(channel == Selector::EMu)
    dilepton_max = 1000;
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",50,50,dilepton_max);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  TString chainTitle(chain->GetTitle());
  Float_t ttScaleFactor = 1.0;
  if( chainTitle.EqualTo("TTMC") ){
	  ttScaleFactor = 0.959;	//to account for slightly higher number of ttbar EMu MC events relative to EMu data evts
  }

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);

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

  //normalize histo bins
  unsigned int max = hs->size();
  for(unsigned int i=0; i<max; i++){
	  //get num bins in histo i
	  Int_t nBins = (hs->at(i))->GetNbinsX();
	  for(Int_t j=1; j<=nBins; j++){
		  //in each bin, divide the bin contents by the bin width
		  Double_t oldBinContents = (hs->at(i))->GetBinContent(j);
		  Double_t oldBinErrors = (hs->at(i))->GetBinError(j);
		  Double_t binWidth = (hs->at(i))->GetBinWidth(j);
		  (hs->at(i))->SetBinContent(j, oldBinContents/binWidth);
		  (hs->at(i))->SetBinError(j, oldBinErrors/binWidth);
	  }//end loop over bins in histo

  }//end loop over histos in vector
 
}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ; 
  //leg->AddEntry( hs_DY, "DY" ) ; 
  leg->AddEntry( hs_ttbar, "ttbar" ) ;
  leg->AddEntry( hs_WJets, "WJets" ) ; 
  leg->AddEntry( hs_WZ, "WZ" ) ; 
  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
  //leg->AddEntry( histos[2][0], "10 x WR 2600" ) ; 
  leg->AddEntry( hs_data, "Data");
  leg->SetFillColor( kWhite ) ; 

  hs_data->Sumw2();
  hs_ttbar->Sumw2();
  hs_WJets->Sumw2();
  hs_WZ->Sumw2();
  hs_ZZ->Sumw2();
  
  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WJets->SetFillColor(kBlue);
  hs_WZ->SetFillColor(kCyan);
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_WZ);
  th->Add(hs_WJets);
  th->Add(hs_ZZ);
  th->Add(hs_ttbar);
  //th->Add(hs_DY);
  hs_data->SetMarkerStyle(20);

  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Preliminary #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->SetTitle("CMS Preliminary #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->Draw("ep");
  th->Draw("histo same");
  hs_data->Draw("epsame");
  TString ytitle = "Events/(";
  ytitle += (th->GetXaxis()->GetNbins());
  ytitle += ")";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  hs_data->GetYaxis()->SetTitle(ytitle.Data());
  if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetYaxis()->SetTitle("Events/GeV"), hs_data->GetYaxis()->SetTitle("Events/GeV");
 
  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")) ratio->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(0.2);
  ratio->SetLabelSize(0.1,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();
  ratio->Sumw2();
  ratio->SetStats(0);

  //hs_DY->Add(hs_ttbar);
  //hs_DY->Add(hs_WJets);
  //hs_DY->Add(hs_WZ);
  //hs_DY->Add(hs_ZZ);
  hs_ttbar->Add(hs_WJets);
  hs_ttbar->Add(hs_WZ);
  hs_ttbar->Add(hs_ZZ);
  if(fname.EqualTo("Mll") || fname.EqualTo("Mlljj") ){
	  Float_t dataMCratio = (hs_data->Integral()/hs_ttbar->Integral());
	  Float_t dataEntries = hs_data->GetEntries();
	  Float_t mcEntries = (hs_ttbar->GetEntries()) + (hs_WJets->GetEntries()) + (hs_WZ->GetEntries()) + (hs_ZZ->GetEntries());	///<temp fix
	  //Float_t mcEntries = (hs_DY->GetEntries()) + (hs_ttbar->GetEntries()) + (hs_WJets->GetEntries()) + (hs_WZ->GetEntries()) + (hs_ZZ->GetEntries());	///<use this once DY emu is available
	  Float_t integralUnc = (dataEntries/mcEntries)*sqrt((1/dataEntries) + (1/mcEntries));
	  std::cout<< "in EMu channel "<< fname <<" dataOvrMC ratio=\t"<< dataMCratio <<"\t+/-\t"<< integralUnc << std::endl;
  }

  ratio->Divide(hs_ttbar);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(0.1,"y");
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  mycanvas->Update();

  TString fn = fname + "_sidebandEMuChannelNoDyRescaledTTBarMCNoLLJJCut";

  mycanvas->Print((fn+".pdf").Data());
  mycanvas->Print((fn+".png").Data());
  p1->SetLogy();
  mycanvas->Print((fn+"_log.pdf").Data());
  mycanvas->Print((fn+"_log.png").Data());

  mycanvas->Close();
}
