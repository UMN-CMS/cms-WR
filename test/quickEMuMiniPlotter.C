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

//#define doNarrowMlljj
//#define doNarrowLeadLeptEta
#define useDYMAD

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

Selector::tag_t channel = Selector::EMu;
Float_t mcUnweightedEntriesInExcessBin=0;

/**
 * this macro is designed to read several TChains, representing data and MC, apply no cuts, and plot
 * data as points and all MC as one stacked histogram with several fill colors.
 *
 * This macro should be used by itself on minitrees processed by analysis.cpp with -c EMu.
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

  TString localDir = "../analysisCppOutputRootFiles/";
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0;
  switch (channel) {
  case Selector::EMu:
#ifndef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_flavoursidebandEMu.root");
#endif
#ifdef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYMADInclAndHT_flavoursidebandEMu_withMllWeight.root");
#endif
	tt = chain_ttbar->Add(localDir+"selected_tree_TT_flavoursidebandEMu.root");
    wjets = chain_WJets->Add(localDir+"selected_tree_W_flavoursidebandEMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_flavoursidebandEMu.root");
    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_flavoursidebandEMu.root");
    data = chain_data->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
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

  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","#Delta R lead lepton lead jet","#Delta R lead lepton sublead jet","#Delta R sublead lepton lead jet","#Delta R sublead lepton sublead jet","#Delta R lead lepton sublead lepton","#Delta #phi lead lepton lead jet","#Delta #phi lead lepton sublead jet","#Delta #phi sublead lepton lead jet","#Delta #phi sublead lepton sublead jet","#Delta #phi lead lepton sublead lepton","#Delta #eta lead lepton lead jet","#Delta #eta lead lepton sublead jet","#Delta #eta sublead lepton lead jet","#Delta #eta sublead lepton sublead jet","#Delta #eta lead lepton sublead lepton","muon #eta","muon #phi","muon p_{T}","electron #eta","electron #phi","electron p_{T}","unweighted M_{LLJJ}"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","l1_j1_dr","l1_j2_dr","l2_j1_dr","l2_j2_dr","l1_l2_dr","l1_j1_dphi","l1_j2_dphi","l2_j1_dphi","l2_j2_dphi","l1_l2_dphi","l1_j1_deta","l1_j2_deta","l2_j1_deta","l2_j2_deta","l1_l2_deta","mu_eta","mu_phi","mu_pt","ele_eta","ele_phi","ele_pt","unweightedMLLJJ"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WJets[i],hs_WZ[i],hs_ZZ[i],hs_data[i],xtitles[i],fnames[i]);
  }
  
}//end quickEMuMiniPlotter()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  Float_t nBins = 50;

#ifdef doNarrowMlljj
	nBins = 15;	///use coarse binning when making data MC comparisons in narrow MLLJJ window
#endif

#ifdef doNarrowLeadLeptEta
	nBins = 5;
#endif

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",nBins,0,700);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",nBins,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",nBins,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",nBins,0,700);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",nBins,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",nBins,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",nBins,0,700);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",nBins,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",nBins,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",nBins,0,700);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",nBins,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",nBins,-3.15,3.15);

  //TH1F *h_WR_mass = new TH1F("h_WR_mass","",nBins,0,2500);	//fixed bin width

  /**/
  //Float_t bins[] = { 210, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800, 6000};	//show out to 6 TeV without mass cut without overflow
  Float_t bins[] = { 210, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800};	//standard bins without 600 GeV mass cut, with overflow events
  
  ////Float_t bins[] = { 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800, 2500};	//standard bins with 600 GeV mass cut
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1F *h_WR_mass = new TH1F("h_WR_mass","",binnum, bins);
  TH1F *h_WR_mass_unweighted = new TH1F("h_WR_mass_unweighted","",binnum, bins);


  /**/
 
  float dilepton_max = 250.;
  if(channel == Selector::EMu)
    dilepton_max = 1000;
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",nBins,50,dilepton_max);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  TH1F *h_lead_lept_lead_jet_dr = new TH1F("h_lead_lept_lead_jet_dr","",nBins,-5,5);
  TH1F *h_lead_lept_sublead_jet_dr = new TH1F("h_lead_lept_sublead_jet_dr","",nBins,-5,5);
  TH1F *h_sublead_lept_lead_jet_dr = new TH1F("h_sublead_lept_lead_jet_dr","",nBins,-5,5);
  TH1F *h_sublead_lept_sublead_jet_dr = new TH1F("h_sublead_lept_sublead_jet_dr","",nBins,-5,5);
  TH1F *h_lead_lept_sublead_lept_dr = new TH1F("h_lead_lept_sublead_lept_dr","",nBins,-5,5);

  TH1F *h_lead_lept_lead_jet_dphi = new TH1F("h_lead_lept_lead_jet_dphi","",nBins,-6.5,6.5);
  TH1F *h_lead_lept_sublead_jet_dphi = new TH1F("h_lead_lept_sublead_jet_dphi","",nBins,-6.5,6.5);
  TH1F *h_sublead_lept_lead_jet_dphi = new TH1F("h_sublead_lept_lead_jet_dphi","",nBins,-6.5,6.5);
  TH1F *h_sublead_lept_sublead_jet_dphi = new TH1F("h_sublead_lept_sublead_jet_dphi","",nBins,-6.5,6.5);
  TH1F *h_lead_lept_sublead_lept_dphi = new TH1F("h_lead_lept_sublead_lept_dphi","",nBins,-6.5,6.5);

  TH1F *h_lead_lept_lead_jet_deta = new TH1F("h_lead_lept_lead_jet_deta","",nBins,-5.5,5.5);
  TH1F *h_lead_lept_sublead_jet_deta = new TH1F("h_lead_lept_sublead_jet_deta","",nBins,-5.5,5.5);
  TH1F *h_sublead_lept_lead_jet_deta = new TH1F("h_sublead_lept_lead_jet_deta","",nBins,-5.5,5.5);
  TH1F *h_sublead_lept_sublead_jet_deta = new TH1F("h_sublead_lept_sublead_jet_deta","",nBins,-5.5,5.5);
  TH1F *h_lead_lept_sublead_lept_deta = new TH1F("h_lead_lept_sublead_lept_deta","",nBins,-5.5,5.5);

  TH1F *h_muon_pt = new TH1F("h_muon_pt","",nBins,0,700);
  TH1F *h_muon_eta = new TH1F("h_muon_eta","",nBins,-3,3);
  TH1F *h_muon_phi = new TH1F("h_muon_phi","",nBins,-3.15,3.15);
  TH1F *h_ele_pt = new TH1F("h_ele_pt","",nBins,0,700);
  TH1F *h_ele_eta = new TH1F("h_ele_eta","",nBins,-3,3);
  TH1F *h_ele_phi = new TH1F("h_ele_phi","",nBins,-3.15,3.15);

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  TString chainTitle(chain->GetTitle());
  Float_t ttScaleFactor = 1.0;
  if( chainTitle.EqualTo("TTMC") ){
	  ttScaleFactor = 0.959;	//to account for slightly higher number of ttbar EMu MC events relative to EMu data evts
  }

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	//if(myEvent->WR_mass < 600.) continue;	///MLLJJ cut
	if(myEvent->sublead_lepton_pt < 53.) continue;

#ifdef doNarrowMlljj
	if(myEvent->WR_mass < 1350. || myEvent->WR_mass > 1510.) continue;
#endif

#ifdef doNarrowLeadLeptEta
	if(myEvent->WR_mass < 1350. || myEvent->WR_mass > 1510.) continue;
	if(myEvent->lead_lepton_eta < 0. || myEvent->lead_lepton_eta > 1.) continue;
#endif

	if(myEvent->WR_mass >= 1350 && myEvent->WR_mass < 1510 && !(chainTitle.EqualTo("Data")) ) mcUnweightedEntriesInExcessBin++;

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
    h_WR_mass_unweighted->Fill(myEvent->WR_mass,1*ttScaleFactor);
	h_dilepton_mass->Fill(myEvent->dilepton_mass,(myEvent->weight)*ttScaleFactor);
    h_nPV->Fill(myEvent->nPV,(myEvent->weight)*ttScaleFactor);

	TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, leadJetFourMom, subleadJetFourMom;
	leadLeptonFourMom.SetPtEtaPhiE(myEvent->lead_lepton_pt, myEvent->lead_lepton_eta, myEvent->lead_lepton_phi, myEvent->lead_lepton_pt);
	subleadLeptonFourMom.SetPtEtaPhiE(myEvent->sublead_lepton_pt, myEvent->sublead_lepton_eta, myEvent->sublead_lepton_phi, myEvent->sublead_lepton_pt);
	//E set here for jets should be updated, but doesn't affect dR deta and dphi
	leadJetFourMom.SetPtEtaPhiE(myEvent->lead_jet_pt, myEvent->lead_jet_eta, myEvent->lead_jet_phi, myEvent->lead_jet_pt);
	subleadJetFourMom.SetPtEtaPhiE(myEvent->sublead_jet_pt, myEvent->sublead_jet_eta, myEvent->sublead_jet_phi, myEvent->sublead_jet_pt);
	
	h_lead_lept_sublead_lept_dr->Fill(leadLeptonFourMom.DeltaR(subleadLeptonFourMom),(myEvent->weight)*ttScaleFactor);
  	h_lead_lept_lead_jet_dr->Fill(myEvent->dR_leadlepton_leadjet,(myEvent->weight)*ttScaleFactor);
  	h_lead_lept_sublead_jet_dr->Fill(myEvent->dR_leadlepton_subleadjet,(myEvent->weight)*ttScaleFactor);
  	h_sublead_lept_lead_jet_dr->Fill(myEvent->dR_subleadlepton_leadjet,(myEvent->weight)*ttScaleFactor);
  	h_sublead_lept_sublead_jet_dr->Fill(myEvent->dR_subleadlepton_subleadjet,(myEvent->weight)*ttScaleFactor);

	h_lead_lept_sublead_lept_dphi->Fill(leadLeptonFourMom.DeltaPhi(subleadLeptonFourMom),(myEvent->weight)*ttScaleFactor);
  	h_lead_lept_lead_jet_dphi->Fill(leadLeptonFourMom.DeltaPhi(leadJetFourMom),(myEvent->weight)*ttScaleFactor);
  	h_lead_lept_sublead_jet_dphi->Fill(leadLeptonFourMom.DeltaPhi(subleadJetFourMom),(myEvent->weight)*ttScaleFactor);
  	h_sublead_lept_lead_jet_dphi->Fill(subleadLeptonFourMom.DeltaPhi(leadJetFourMom),(myEvent->weight)*ttScaleFactor);
  	h_sublead_lept_sublead_jet_dphi->Fill(subleadLeptonFourMom.DeltaPhi(subleadJetFourMom),(myEvent->weight)*ttScaleFactor);

	h_lead_lept_sublead_lept_deta->Fill(myEvent->lead_lepton_eta - myEvent->sublead_lepton_eta,(myEvent->weight)*ttScaleFactor);
  	h_lead_lept_lead_jet_deta->Fill(myEvent->lead_lepton_eta - myEvent->lead_jet_eta,(myEvent->weight)*ttScaleFactor);
   	h_lead_lept_sublead_jet_deta->Fill(myEvent->lead_lepton_eta - myEvent->sublead_jet_eta,(myEvent->weight)*ttScaleFactor);
   	h_sublead_lept_lead_jet_deta->Fill(myEvent->sublead_lepton_eta - myEvent->lead_jet_eta,(myEvent->weight)*ttScaleFactor);
   	h_sublead_lept_sublead_jet_deta->Fill(myEvent->sublead_lepton_eta - myEvent->sublead_jet_eta,(myEvent->weight)*ttScaleFactor);

	if(myEvent->lead_lepton_r9 == -1){
		h_muon_pt->Fill(myEvent->lead_lepton_pt,(myEvent->weight)*ttScaleFactor);
		h_muon_eta->Fill(myEvent->lead_lepton_eta,(myEvent->weight)*ttScaleFactor);
		h_muon_phi->Fill(myEvent->lead_lepton_phi,(myEvent->weight)*ttScaleFactor);
		h_ele_pt->Fill(myEvent->sublead_lepton_pt,(myEvent->weight)*ttScaleFactor);
		h_ele_eta->Fill(myEvent->sublead_lepton_eta,(myEvent->weight)*ttScaleFactor);
		h_ele_phi->Fill(myEvent->sublead_lepton_phi,(myEvent->weight)*ttScaleFactor);
	}//lead lepton is muon

	if(myEvent->lead_lepton_r9 > -1){
		h_ele_pt->Fill(myEvent->lead_lepton_pt,(myEvent->weight)*ttScaleFactor);
		h_ele_eta->Fill(myEvent->lead_lepton_eta,(myEvent->weight)*ttScaleFactor);
		h_ele_phi->Fill(myEvent->lead_lepton_phi,(myEvent->weight)*ttScaleFactor);
		h_muon_pt->Fill(myEvent->sublead_lepton_pt,(myEvent->weight)*ttScaleFactor);
		h_muon_eta->Fill(myEvent->sublead_lepton_eta,(myEvent->weight)*ttScaleFactor);
		h_muon_phi->Fill(myEvent->sublead_lepton_phi,(myEvent->weight)*ttScaleFactor);
	}//lead lepton is electron

  }//end loop over events

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
  hs->push_back(h_lead_lept_lead_jet_dr);
  hs->push_back(h_lead_lept_sublead_jet_dr);
  hs->push_back(h_sublead_lept_lead_jet_dr);
  hs->push_back(h_sublead_lept_sublead_jet_dr);
  hs->push_back(h_lead_lept_sublead_lept_dr);
  hs->push_back(h_lead_lept_lead_jet_dphi);
  hs->push_back(h_lead_lept_sublead_jet_dphi);
  hs->push_back(h_sublead_lept_lead_jet_dphi);
  hs->push_back(h_sublead_lept_sublead_jet_dphi);
  hs->push_back(h_lead_lept_sublead_lept_dphi);
  hs->push_back(h_lead_lept_lead_jet_deta);
  hs->push_back(h_lead_lept_sublead_jet_deta);
  hs->push_back(h_sublead_lept_lead_jet_deta);
  hs->push_back(h_sublead_lept_sublead_jet_deta);
  hs->push_back(h_lead_lept_sublead_lept_deta);
  hs->push_back(h_muon_eta);
  hs->push_back(h_muon_phi);
  hs->push_back(h_muon_pt);
  hs->push_back(h_ele_eta);
  hs->push_back(h_ele_phi);
  hs->push_back(h_ele_pt);
  hs->push_back(h_WR_mass_unweighted);
	
  /**/
  //normalize histo bins of WR mass histos
  unsigned int max = hs->size();
  for(unsigned int i=0; i<max; i++){
	  //get num bins in histo i
	  Int_t nBins = (hs->at(i))->GetNbinsX();
	  TString histName = hs->at(i)->GetName();
	  if( histName.Contains("WR_mass") ){
		  for(Int_t j=1; j<=nBins; j++){
			  //in each bin, divide the bin contents by the bin width
			  Double_t oldBinContents = (hs->at(i))->GetBinContent(j);
			  Double_t oldBinErrors = (hs->at(i))->GetBinError(j);
			  Double_t binWidth = (hs->at(i))->GetBinWidth(j);
			  (hs->at(i))->SetBinContent(j, oldBinContents/binWidth);
			  (hs->at(i))->SetBinError(j, oldBinErrors/binWidth);
		  }//end loop over bins in histo

	  }
  }//end loop over histos in vector
  /**/ 
}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.60, 0.60, 0.90, 0.90 ) ; 
#ifndef useDYMAD 
  leg->AddEntry( hs_DY, "DY AMCNLO" ) ; 
#endif
#ifdef useDYMAD
  leg->AddEntry( hs_DY, "DYMadHT+Incl" ) ; 
#endif
  leg->AddEntry( hs_ttbar, "TT MC" ) ;
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
  hs_DY->Sumw2();
  
  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 900, 900 ) ;
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WJets->SetFillColor(kBlue);
  hs_WZ->SetFillColor(kCyan);
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_WZ);
  th->Add(hs_WJets);
  th->Add(hs_DY);
  th->Add(hs_ZZ);
  th->Add(hs_ttbar);
  hs_data->SetMarkerStyle(20);

  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->Draw("ep");
  th->Draw("histo same");
  hs_data->Draw("epsame");
  TString ytitle = "Events/(";
  ytitle += (th->GetXaxis()->GetNbins());
  ytitle += ")";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  hs_data->GetYaxis()->SetTitle(ytitle.Data());
  //for variable size bins normalized to bin width
  if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetYaxis()->SetTitle("Events/GeV   "), hs_data->GetYaxis()->SetTitle("Events/GeV   ");
  //if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetYaxis()->SetTitle("Events"), hs_data->GetYaxis()->SetTitle("Events");


#ifdef doNarrowMlljj
	th->GetYaxis()->SetTitle("Events"), hs_data->GetYaxis()->SetTitle("Events");
#endif

#ifdef doNarrowLeadLeptEta
	th->GetYaxis()->SetTitle("Events"), hs_data->GetYaxis()->SetTitle("Events");
#endif


  Float_t labelSize = 0.25;
  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")) ratio->GetXaxis()->SetTitle("M_{EMuJJ} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(labelSize+0.03);
  ratio->SetLabelSize(labelSize - 0.07,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();
  ratio->Sumw2();
  ratio->SetStats(0);

  hs_ttbar->Add(hs_WJets);
  hs_ttbar->Add(hs_WZ);
  hs_ttbar->Add(hs_ZZ);
  hs_ttbar->Add(hs_DY);
  if(fname.EqualTo("Mlljj") ){
	  Float_t dataMCratio = (hs_data->Integral()/hs_ttbar->Integral());
	  Float_t dataEntries = hs_data->GetEntries();
	  Float_t mcEntries = (hs_ttbar->GetEntries()) + (hs_WJets->GetEntries()) + (hs_WZ->GetEntries()) + (hs_ZZ->GetEntries());	///<temp fix
	  //Float_t mcEntries = (hs_DY->GetEntries()) + (hs_ttbar->GetEntries()) + (hs_WJets->GetEntries()) + (hs_WZ->GetEntries()) + (hs_ZZ->GetEntries());	///<use this once DY emu is available
	  Float_t integralUnc = (dataEntries/mcEntries)*sqrt((1/dataEntries) + (1/mcEntries));
	  std::cout<< "in EMu channel "<< fname <<" dataOvrMC ratio=\t"<< dataMCratio <<"\t+/-\t"<< integralUnc << std::endl;
	  std::cout<<" "<<std::endl;
  	  std::cout<<"bin number\t"<< 6 <<"has data bin contents=\t" << hs_data->GetBinContent(6) <<" and bin error =\t"<< hs_data->GetBinError(6) << std::endl;
  	  std::cout<<"bin number\t"<< 6 <<"has ttbar bin contents=\t" << hs_ttbar->GetBinContent(6) <<" and bin error =\t"<< hs_ttbar->GetBinError(6) << std::endl;
	  std::cout<<"bin number\t"<< 7 <<"has data bin contents=\t" << hs_data->GetBinContent(7) <<" and bin error =\t"<< hs_data->GetBinError(7) << std::endl;
  	  std::cout<<"bin number\t"<< 7 <<"has ttbar bin contents=\t" << hs_ttbar->GetBinContent(7) <<" and bin error =\t"<< hs_ttbar->GetBinError(7) << std::endl;

	  //print actual number of evts in each bin, raw unweighted events, and bin lower edge
	  Int_t centBin = 16;
	  std::cout<<"bin num\t"<< centBin <<"\thas\t"<< hs_ttbar->GetBinLowEdge(centBin) <<"\tGeV lower edge and\t"<< (hs_ttbar->GetBinContent(centBin))*(hs_ttbar->GetBinWidth(centBin)) <<" +/- "<< (hs_ttbar->GetBinError(centBin))*(hs_ttbar->GetBinWidth(centBin)) <<" MC events\t" <<"and unweighted entries=\t" << mcUnweightedEntriesInExcessBin <<std::endl;
	  std::cout<<"bin num\t"<< centBin <<"\thas\t"<< hs_data->GetBinLowEdge(centBin) <<"\tGeV lower edge and\t"<< (hs_data->GetBinContent(centBin))*(hs_data->GetBinWidth(centBin)) <<" +/- "<< (hs_data->GetBinError(centBin))*(hs_data->GetBinWidth(centBin)) <<" data events\t" <<"and unweighted entries=\t" << (hs_data->GetBinContent(centBin))*(hs_data->GetBinWidth(centBin)) <<std::endl;
	  
	  std::cout<<"\t"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<<"bin num "<< i <<" has "<< hs_ttbar->GetBinLowEdge(i) <<" GeV lower edge and "<< (hs_ttbar->GetBinContent(i))*(hs_ttbar->GetBinWidth(i)) <<" +/- "<< (hs_ttbar->GetBinError(i))*(hs_ttbar->GetBinWidth(i)) <<" weighted MC events and " << (hs_data->GetBinContent(i))*(hs_data->GetBinWidth(i)) <<" data evts" <<std::endl;
	  }//end loop over bins in MLLJJ
  }

  if(fname.EqualTo("unweightedMLLJJ") ){
	  std::cout<<"\t"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<<"bin num "<< i <<" has "<< hs_ttbar->GetBinLowEdge(i) <<" GeV lower edge and "<< (hs_ttbar->GetBinContent(i))*(hs_ttbar->GetBinWidth(i)) <<" unweighted MC events" <<std::endl;
	  }//end loop over bins in unweightedMLLJJ
  }


  ratio->Divide(hs_ttbar);
  ratio->SetMarkerStyle(21);
  ratio->GetYaxis()->SetTitle("data/MC    ");
  ratio->GetYaxis()->SetTitleSize(0.19);
  ratio->GetYaxis()->SetTitleOffset(0.2);
  ratio->SetLabelSize(labelSize - 0.07,"y");
  ratio->GetYaxis()->SetRangeUser(0.5,2.4);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  mycanvas->Update();

  //TString fn = fname + "_eMuChannelRescaledTTBarMCNoLLJJCutFixedBinWidth";
  //TString fn = fname + "_eMuChannelRescaledTTBarMCNoLLJJCutVariableBinWidthSixTeVMax";
#ifndef useDYMAD
  TString fn = fname + "_eMuChannelRescaledTTBarMC_DYAMC_NoLLJJCutVariableBinWidth";
#endif

#ifdef useDYMAD
  TString fn = fname + "_eMuChannelRescaledTTBarMC_DYMadHTAndIncl_NoLLJJCutVariableBinWidth";
#endif

#ifdef doNarrowMlljj
	fn = fname + "_eMuChannelRescaledTTBarMCNarrowMlljjWindowFixedBinWidth";
#endif

#ifdef doNarrowLeadLeptEta
	fn = fname + "_eMuChannelRescaledTTBarMCNarrowMlljjAndLeadLeptEtaFixedBinWidth";
#endif


  if(fname.EqualTo("Mlljj")){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  mycanvas->Print((fn+".C").Data());
	  p1->SetLogy();
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
	  mycanvas->Print((fn+"_log.C").Data());
  }

  mycanvas->Close();
}
