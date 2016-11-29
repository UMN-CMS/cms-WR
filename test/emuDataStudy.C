#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <string>
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>

/**
 * this macro is a compliment to quickEMuMiniPlotter.C
 * This macro should be run by itself using minitrees which have been processed by analysis.cpp with -c EMu.
 * It will create plots showing data points and multiple curves, nothing is stacked.
 *
 */

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif
TString dir = "../analysisCppOutputRootFiles/";

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t mlljjLowerCut, Float_t mlljjUpperCut);
void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname);
void emuDataStudy()
{
	
	TString treeName = "Tree_Iter0";
	TChain * chain_DYPowhegEE = new TChain(treeName,"extraBinEE");
	TChain * chain_DYMadInclEE = new TChain(treeName,"lowerBinEE");
	TChain * chain_DYAmcInclEE = new TChain(treeName,"upperBinEE");
	TChain * chain_dataEE = new TChain(treeName,"centralBinEE");

	chain_DYPowhegEE->Add(dir+"selected_tree_data_flavoursidebandEMuEE_cpyThree.root");
	chain_DYMadInclEE->Add(dir+"selected_tree_data_flavoursidebandEMuEE_cpyOne.root");
	chain_DYAmcInclEE->Add(dir+"selected_tree_data_flavoursidebandEMuEE_cpyTwo.root");
	chain_dataEE->Add(dir+"selected_tree_data_flavoursidebandEMuEE.root");
	
	Selector myEvent_DYPowhegEE;
	Selector myEvent_DYMadInclEE;
	Selector myEvent_DYAmcInclEE;
	Selector myEvent_dataEE;
	
	myEvent_DYPowhegEE.SetBranchAddresses(chain_DYPowhegEE);
	myEvent_DYMadInclEE.SetBranchAddresses(chain_DYMadInclEE);
	myEvent_DYAmcInclEE.SetBranchAddresses(chain_DYAmcInclEE);
	myEvent_dataEE.SetBranchAddresses(chain_dataEE);

	std::vector<TH1F*> hs_DYPowhegEE;
	MakeHistos(chain_DYPowhegEE, &myEvent_DYPowhegEE, &hs_DYPowhegEE, 0, 6000);
	std::vector<TH1F*> hs_DYMadInclEE;
	MakeHistos(chain_DYMadInclEE, &myEvent_DYMadInclEE, &hs_DYMadInclEE, 1210, 1350);
	std::vector<TH1F*> hs_DYAmcInclEE;
	MakeHistos(chain_DYAmcInclEE, &myEvent_DYAmcInclEE, &hs_DYAmcInclEE, 1510, 1800);
	std::vector<TH1F*> hs_dataEE;
	MakeHistos(chain_dataEE, &myEvent_dataEE, &hs_dataEE, 1350, 1510);

	//now the vectors of TH1F pointers are filled with pointers to histos with nonzero entries
	unsigned int nPlots = hs_dataEE.size();

  	TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","#Delta R lead lepton lead jet","#Delta R lead lepton sublead jet","#Delta R sublead lepton lead jet","#Delta R sublead lepton sublead jet","#Delta R lead lepton sublead lepton","#Delta #phi lead lepton lead jet","#Delta #phi lead lepton sublead jet","#Delta #phi sublead lepton lead jet","#Delta #phi sublead lepton sublead jet","#Delta #phi lead lepton sublead lepton","#Delta #eta lead lepton lead jet","#Delta #eta lead lepton sublead jet","#Delta #eta sublead lepton lead jet","#Delta #eta sublead lepton sublead jet","#Delta #eta lead lepton sublead lepton","muon #eta","muon #phi","muon p_{T}","electron #eta","electron #phi","electron p_{T}"};
 
   	TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","l1_j1_dr","l1_j2_dr","l2_j1_dr","l2_j2_dr","l1_l2_dr","l1_j1_dphi","l1_j2_dphi","l2_j1_dphi","l2_j2_dphi","l1_l2_dphi","l1_j1_deta","l1_j2_deta","l2_j1_deta","l2_j2_deta","l1_l2_deta","mu_eta","mu_phi","mu_pt","ele_eta","ele_phi","ele_pt"};

	int i = 0;
	for(unsigned int i = 0; i < nPlots; i++) {
		std::string s = std::to_string(i);
		drawPlots(hs_DYPowhegEE[i], hs_DYMadInclEE[i], hs_DYAmcInclEE[i], hs_dataEE[i], xtitles[i], fnames[i]);
	}

}//end emuDataStudy()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t mlljjLowerCut, Float_t mlljjUpperCut){

  Float_t nBins = 15;

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",nBins,0,600);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",nBins,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",nBins,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",nBins,0,600);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",nBins,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",nBins,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",nBins,0,800);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",nBins,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",nBins,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",nBins,0,500);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",nBins,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",nBins,-3.15,3.15);

  TH1F *h_WR_mass = new TH1F("h_WR_mass","",nBins,0,2500);	//fixed bin width

  /*
  //Float_t bins[] = { 210, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800, 6000};	//show out to 6 TeV without mass cut without overflow
  Float_t bins[] = { 210, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800};	//standard bins without 600 GeV mass cut, with overflow events
  
  ////Float_t bins[] = { 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1800, 2500};	//standard bins with 600 GeV mass cut
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1F *h_WR_mass = new TH1F("h_WR_mass","",binnum, bins);
  */
 
  float dilepton_max = 1000.;
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

  TH1F *h_muon_pt = new TH1F("h_muon_pt","",nBins,0,600);
  TH1F *h_muon_eta = new TH1F("h_muon_eta","",nBins,-3,3);
  TH1F *h_muon_phi = new TH1F("h_muon_phi","",nBins,-3.15,3.15);
  TH1F *h_ele_pt = new TH1F("h_ele_pt","",nBins,0,600);
  TH1F *h_ele_eta = new TH1F("h_ele_eta","",nBins,-3,3);
  TH1F *h_ele_phi = new TH1F("h_ele_phi","",nBins,-3.15,3.15);

  Long64_t nEntries = chain->GetEntries();
  cout<< nEntries << endl;
  Float_t ttScaleFactor = 1.0;
  
  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	if(myEvent->WR_mass < mlljjLowerCut || myEvent->WR_mass > mlljjUpperCut) continue;

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
	
  /*
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
  */ 
}//end MakeHistos()


void drawPlots(TH1F* hs_DYPowheg, TH1F* hs_DYMadIncl, TH1F* hs_DYAmcIncl, TH1F* hs_data, TString xtitle, TString fname)
{

	//gStyle->SetOptStat("eou");
	gStyle->SetOptStat("");
	TLegend *leg = new TLegend( 0.60, 0.60, 0.90, 0.90 ) ;
	//leg->AddEntry( hs_DYPowheg, "extra" ) ;
	leg->AddEntry( hs_DYMadIncl, "lower" ) ;
	leg->AddEntry( hs_DYAmcIncl, "upper" ) ;
	leg->AddEntry( hs_data, "central");
	leg->SetFillColor( kWhite ) ;

	TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
	mycanvas->cd();
	hs_DYPowheg->SetLineColor(kBlue);
	hs_DYPowheg->SetLineWidth(3);
	hs_DYMadIncl->SetLineColor(kBlack);
	hs_DYMadIncl->SetLineWidth(3);
	hs_DYAmcIncl->SetLineColor(kRed);
	hs_DYAmcIncl->SetLineWidth(3);
	hs_data->SetMarkerStyle(20);
	hs_data->SetMarkerSize(1);
	hs_data->SetMarkerColor(kBlack);

	/*for ratio plot
	Double_t eps = 0.001;
	TPad* p1 = new TPad("p1", "p1", 0, 0.25, 1, 1, 0);
	p1->Draw();
	TPad* p2 = new TPad("p2", "p2", 0, 0.1, 1, 0.25 + eps, 0);
	p2->Draw();
	p1->SetBottomMargin(0);
	p2->SetTopMargin(0);
	p1->cd();
	*/
	hs_data->SetStats(1);
	hs_DYPowheg->SetStats(1);
	TH1F *ratio_Powheg = (TH1F*)hs_data->Clone();
	TH1F *ratio_Mad = (TH1F*)hs_data->Clone();
	TH1F *ratio_Amc = (TH1F*)hs_data->Clone();
	TString plotTitle = "CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}";
	hs_DYPowheg->SetTitle(plotTitle);
	hs_data->SetTitle(plotTitle);
	TString ytitle = "Events";
	//TString ytitle = "Events/(";
	//ytitle += (hs_data->GetXaxis()->GetNbins());
	//ytitle += (hs_data->GetXaxis()->GetBinWidth(2));
	//ytitle += " GeV)";
	hs_DYAmcIncl->GetYaxis()->SetTitle(ytitle.Data());
	//hs_DYPowheg->GetYaxis()->SetTitle(ytitle.Data());
	hs_DYMadIncl->GetYaxis()->SetTitle(ytitle.Data());
	hs_data->GetYaxis()->SetTitle(ytitle.Data());
	
	hs_data->GetYaxis()->SetTitleOffset(1.5);
	hs_DYAmcIncl->GetYaxis()->SetTitleOffset(1.5);
	hs_DYMadIncl->GetYaxis()->SetTitleOffset(1.5);
	//hs_DYPowheg->GetYaxis()->SetTitleOffset(1.5);
	//hs_data->GetYaxis()->SetTitleSize(0.042);
	//hs_DYAmcIncl->GetYaxis()->SetTitleSize(0.042);
	
	hs_data->GetXaxis()->SetTitleSize(0.04);
	//hs_data->SetLabelSize(0.04, "x");
	//hs_data->SetLabelSize(0.025, "y");
	hs_data->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYAmcIncl->GetXaxis()->SetTitleSize(0.04);
	hs_DYAmcIncl->SetLabelSize(0.04, "x");
	hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYPowheg->GetXaxis()->SetTitleSize(0.04);
	hs_DYPowheg->SetLabelSize(0.04, "x");
	hs_DYPowheg->SetLabelSize(0.025, "y");
	hs_DYPowheg->GetXaxis()->SetTitle(xtitle.Data());
	hs_DYMadIncl->GetXaxis()->SetTitleSize(0.04);
	hs_DYMadIncl->SetLabelSize(0.04, "x");
	hs_DYMadIncl->GetXaxis()->SetTitle(xtitle.Data());
	if(fname.EqualTo("Z_pt") ) hs_data->GetXaxis()->SetTitle("Z P_{T} [GeV]"), hs_DYAmcIncl->GetXaxis()->SetTitle("Z P_{T} [GeV]");
	Float_t max = hs_data->GetBinContent(hs_data->GetMaximumBin());
	Float_t maxLower = hs_DYAmcIncl->GetBinContent(hs_DYAmcIncl->GetMaximumBin());
	Float_t maxUpper = hs_DYMadIncl->GetBinContent(hs_DYMadIncl->GetMaximumBin());
	if(max < maxLower) max = maxLower;
	if(max < maxUpper) max = maxUpper;
	if(fname.EqualTo("mu_eta") ) hs_data->SetMaximum(1.5*max), hs_DYAmcIncl->SetMaximum(1.5*max), hs_DYMadIncl->SetMaximum(1.5*max);
	
	hs_data->Draw("ep");
	//hs_DYPowheg->Draw("histo same");	//comment if only plotting AMC
	hs_DYMadIncl->Draw("histo same");	//comment if only plotting AMC
	hs_DYAmcIncl->Draw("histo same");
	hs_data->Draw("epsame");
	
	hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Amc->GetXaxis()->SetTitle(xtitle.Data());
	if(fname.EqualTo("Mll")) ratio_Amc->GetXaxis()->SetTitle("M_{LL} [GeV]");
	if(fname.EqualTo("Z_pt")) ratio_Amc->GetXaxis()->SetTitle("Z P_{T} [GeV]");
	if(fname.EqualTo("nPV")) ratio_Amc->GetXaxis()->SetTitle("nPV");
	if(fname.EqualTo("nPU")) ratio_Amc->GetXaxis()->SetTitle("nPU");

	Float_t labelSize = 0.25;
	ratio_Amc->GetXaxis()->SetTickSize(0.40);
	ratio_Amc->GetXaxis()->SetTitleSize(labelSize);
	ratio_Amc->SetLabelSize(labelSize - 0.07, "x");
	
	ratio_Powheg->GetXaxis()->SetTitle(xtitle.Data());
	ratio_Powheg->GetXaxis()->SetTickSize(0.40);
	ratio_Powheg->GetXaxis()->SetTitleSize(labelSize);
	ratio_Powheg->SetLabelSize(labelSize - 0.07, "x");
	leg->Draw();
	mycanvas->cd();
	//p2->cd();	//for ratio plot
	ratio_Powheg->Sumw2();
	ratio_Powheg->SetStats(0);
	ratio_Mad->Sumw2();
	ratio_Mad->SetStats(0);
	ratio_Amc->Sumw2();
	ratio_Amc->SetStats(0);

	ratio_Powheg->Divide(hs_DYPowheg);
	ratio_Powheg->SetMarkerStyle(20);
	ratio_Powheg->SetMarkerColor(kRed);
	ratio_Powheg->SetLabelSize(labelSize - 0.07, "y");
	ratio_Powheg->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Powheg->GetYaxis()->SetNdivisions(505);

	ratio_Mad->Divide(hs_DYMadIncl);
	ratio_Mad->SetMarkerStyle(21);
	ratio_Mad->SetMarkerColor(kBlack);
	ratio_Mad->SetLabelSize(labelSize - 0.07, "y");
	ratio_Mad->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Mad->GetYaxis()->SetNdivisions(505);

	ratio_Amc->Divide(hs_DYAmcIncl);
	ratio_Amc->SetMarkerStyle(22);
	ratio_Amc->SetMarkerColor(kBlue);
	ratio_Amc->SetLabelSize(labelSize - 0.07, "y");
	ratio_Amc->GetYaxis()->SetRangeUser(0.95, 1.05);
	ratio_Amc->GetYaxis()->SetNdivisions(505);

	/*for ratio plot
	ratio_Mad->Draw("p");	//comment if only plotting AMC
	ratio_Amc->Draw("p");
	//ratio_Powheg->Draw("p");	//comment if only plotting AMC
	float xmax = ratio_Amc->GetXaxis()->GetXmax();
	float xmin = ratio_Amc->GetXaxis()->GetXmin();
	TF1 *f1 = new TF1("f1", "1", xmin, xmax);
	//ratio_Powheg->Draw("p");
	ratio_Mad->Draw("psame");
	ratio_Amc->Draw("p");
	f1->Draw("same");
	mycanvas->cd();
	*/

	TString cuts = "_emu_data_selfComparison_1210to1800_MLLJJ";
	TString fn = fname + cuts;

	mycanvas->Print((fn + ".pdf").Data());
	mycanvas->Print((fn + ".png").Data());

	mycanvas->SetLogy();	//when not plotting ratio plot
	mycanvas->Print((fn + "_log.pdf").Data());
	mycanvas->Print((fn + "_log.png").Data());

	mycanvas->Close();

}
