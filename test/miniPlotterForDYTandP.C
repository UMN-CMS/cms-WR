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
#include <vector>
#include <iostream>
#include <string>
// #include "ExoAnalysis/cmsWR/src/Selector.cc"
// #include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>


#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

//Selector::tag_t channel = Selector::MuMu;
Selector::tag_t channel = Selector::EE;



void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t leadJetPtCut, Float_t leadLeptonPtCut, Float_t subleadLeptonPtCut, Float_t upperMllCut, Float_t lowerMllCut);
void drawPlots(TH1F* hs_DYPowheg,TH1F* hs_DYMadIncl,TH1F* hs_DYAmcIncl,TH1F* hs_data, TString xtitle, TString fname);
void miniPlotterForDYTandP(){

	TString treeName = "treeDyCheck";
  TChain * chain_DYPowheg = new TChain(treeName);
  TChain * chain_DYMadIncl = new TChain(treeName);
  TChain * chain_DYAmcIncl = new TChain(treeName);
  TChain * chain_data = new TChain(treeName);
  
  switch (channel) {
  case Selector::EE:
    chain_DYPowheg->Add("../selected_tree_DYPOWHEGTANDP_dytagandprobeEE.root");
    chain_DYMadIncl->Add("../selected_tree_DYMADINCLTANDP_dytagandprobeEE.root"); // 0 - Electrons
    chain_DYAmcIncl->Add("../selected_tree_DYAMCINCLTANDP_dytagandprobeEE.root");
    chain_data->Add("../selected_tree_dataEETANDP_dytagandprobeEE.root");
    break;
  case Selector::MuMu:
    chain_DYPowheg->Add("../selected_tree_DYPOWHEGTANDP_dytagandprobeMuMu.root");
    chain_DYMadIncl->Add("../selected_tree_DYMADINCLTANDP_dytagandprobeMuMu.root"); // 1 - Muons
    chain_DYAmcIncl->Add("../selected_tree_DYAMCINCLTANDP_dytagandprobeMuMu.root");
    chain_data->Add("../selected_tree_dataMuMuTANDP_dytagandprobeMuMu.root");
    break;
  default:
    std::cout << "Unknown tag" << std::endl;
  }

  Selector myEvent_DYPowheg;
  Selector myEvent_DYMadIncl;
  Selector myEvent_DYAmcIncl;
  Selector myEvent_data;

  myEvent_DYPowheg.SetBranchAddresses(chain_DYPowheg);
  myEvent_DYMadIncl.SetBranchAddresses(chain_DYMadIncl);
  myEvent_DYAmcIncl.SetBranchAddresses(chain_DYAmcIncl);
  myEvent_data.SetBranchAddresses(chain_data);

  Float_t minJetPt = 20;
  std::vector<TH1F*> hs_DYPowheg;
  MakeHistos(chain_DYPowheg, &myEvent_DYPowheg, &hs_DYPowheg, minJetPt, 32, 20, 120, 60);
  std::vector<TH1F*> hs_DYMadIncl;
  MakeHistos(chain_DYMadIncl, &myEvent_DYMadIncl, &hs_DYMadIncl, minJetPt, 32, 20, 120, 60);
  std::vector<TH1F*> hs_DYAmcIncl;
  MakeHistos(chain_DYAmcIncl, &myEvent_DYAmcIncl, &hs_DYAmcIncl, minJetPt, 32, 20, 120, 60);

  std::vector<TH1F*> hs_data;
  MakeHistos(chain_data, &myEvent_data, &hs_data, minJetPt, 32, 20, 120, 60);

  unsigned int nPlots = hs_DYPowheg.size();

  // hs_data[13]->SetLineColor(kRed);
  // hs_data[13]->Draw();
  // hs_DYMadIncl[13]->Draw("same");
  
  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","Z #phi","Z rapidity","Z p_{T}"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","Z_phi","Z_rapidity","Z_pt"};

  TString jetCutString = "_minJetPt_";
  jetCutString += minJetPt;
  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DYPowheg[i],hs_DYMadIncl[i],hs_DYAmcIncl[i],hs_data[i],xtitles[i],fnames[i]+jetCutString);
  }
  
}//end miniPlotterForDYTandP()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs, Float_t leadJetPtCut, Float_t leadLeptonPtCut, Float_t subleadLeptonPtCut, Float_t upperMllCut, Float_t lowerMllCut){

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",50,0,200);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",50,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",50,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",50,0,200);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",50,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",50,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",50,0,300);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",50,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",50,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",50,0,300);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",50,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",50,-3.15,3.15);

  TH1F *h_WR_mass = new TH1F("h_WR_mass","",50,0,2000);
  float dilepton_max = 150.;
  if(channel == Selector::EMu)
    dilepton_max = 1000;
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",50,50,dilepton_max);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  TH1F *h_Z_phi = new TH1F("h_Z_phi","",50,-3.15,3.15);
  TH1F *h_Z_rapidity = new TH1F("h_Z_rapidity","",50,-5.,5.);
  TH1F *h_Z_pt = new TH1F("h_Z_pt","",50,0.,200.);


  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	if(myEvent->dilepton_mass > upperMllCut || myEvent->dilepton_mass < lowerMllCut) continue;
	if(myEvent->lead_lepton_pt < leadLeptonPtCut || myEvent->sublead_lepton_pt < subleadLeptonPtCut) continue;
	if(myEvent->lead_jet_pt < leadJetPtCut) continue;

	TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, zFourMom;
	leadLeptonFourMom.SetPtEtaPhiE(myEvent->lead_lepton_pt, myEvent->lead_lepton_eta, myEvent->lead_lepton_phi, myEvent->lead_lepton_pt);
	subleadLeptonFourMom.SetPtEtaPhiE(myEvent->sublead_lepton_pt, myEvent->sublead_lepton_eta, myEvent->sublead_lepton_phi, myEvent->sublead_lepton_pt);
	zFourMom = leadLeptonFourMom + subleadLeptonFourMom;

	h_Z_pt->Fill(zFourMom.Pt(), myEvent->weight);
	h_Z_phi->Fill(zFourMom.Phi(), myEvent->weight);
	h_Z_rapidity->Fill(zFourMom.Rapidity(), myEvent->weight);

    h_lepton_pt0->Fill(myEvent->lead_lepton_pt,myEvent->weight);
    h_lepton_pt1->Fill(myEvent->sublead_lepton_pt,myEvent->weight);
    h_lepton_eta0->Fill(myEvent->lead_lepton_eta,myEvent->weight);
    h_lepton_eta1->Fill(myEvent->sublead_lepton_eta,myEvent->weight);
    h_lepton_phi0->Fill(myEvent->lead_lepton_phi,myEvent->weight);
    h_lepton_phi1->Fill(myEvent->sublead_lepton_phi,myEvent->weight);

    h_jet_pt0->Fill(myEvent->lead_jet_pt,myEvent->weight);
    h_jet_pt1->Fill(myEvent->sublead_jet_pt,myEvent->weight);
    h_jet_eta0->Fill(myEvent->lead_jet_eta,myEvent->weight);
    h_jet_eta1->Fill(myEvent->sublead_jet_eta,myEvent->weight);
    h_jet_phi0->Fill(myEvent->lead_jet_phi,myEvent->weight);
    h_jet_phi1->Fill(myEvent->sublead_jet_phi,myEvent->weight);
      
    h_WR_mass->Fill(myEvent->WR_mass,myEvent->weight);
    h_dilepton_mass->Fill(myEvent->dilepton_mass,myEvent->weight);
    h_nPV->Fill(myEvent->nPV,myEvent->weight);
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
  hs->push_back(h_Z_phi);
  hs->push_back(h_Z_rapidity);
  hs->push_back(h_Z_pt);


}

void drawPlots(TH1F* hs_DYPowheg,TH1F* hs_DYMadIncl,TH1F* hs_DYAmcIncl,TH1F* hs_data, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ; 
  leg->AddEntry( hs_DYPowheg, "DY Powheg" ) ; 
  leg->AddEntry( hs_DYMadIncl, "DY MAD Incl" ) ;
  leg->AddEntry( hs_DYAmcIncl, "DY AMC Incl" ) ; 
  //leg->AddEntry( histos[2][0], "10 x WR 2600" ) ; 
  leg->AddEntry( hs_data, "Data");
  leg->SetFillColor( kWhite ) ; 


  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  hs_DYPowheg->SetLineColor(kRed);
  hs_DYPowheg->SetLineWidth(3);
  hs_DYMadIncl->SetLineColor(kBlack);
  hs_DYMadIncl->SetLineWidth(3);
  hs_DYAmcIncl->SetLineColor(kBlue);
  hs_DYAmcIncl->SetLineWidth(3);
  hs_data->SetMarkerStyle(20);
  hs_data->SetMarkerSize(1);
  hs_data->SetMarkerColor(kGreen);


  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  hs_data->SetStats(0);
  TH1F *ratio_Powheg = (TH1F*)hs_data->Clone();
  TH1F *ratio_Mad = (TH1F*)hs_data->Clone();
  TH1F *ratio_Amc = (TH1F*)hs_data->Clone();
  hs_DYPowheg->SetTitle("CMS Preliminary");
  hs_data->SetTitle("CMS Preliminary");
  //th->Draw("histo");
  //hs_data->Draw("epsame");
  hs_data->Draw("ep");
  hs_DYPowheg->Draw("histo same");
  hs_DYMadIncl->Draw("histo same");
  hs_DYAmcIncl->Draw("histo same");
  hs_data->Draw("epsame");
  TString ytitle = "Events/(";
  ytitle += (hs_data->GetXaxis()->GetNbins());
  ytitle += ")";
  hs_DYAmcIncl->GetYaxis()->SetTitle(ytitle.Data());
  hs_DYAmcIncl->GetXaxis()->SetTitle(xtitle.Data());

  ratio_Powheg->GetXaxis()->SetTitle(xtitle.Data());
  ratio_Powheg->GetXaxis()->SetTickSize(0.40);
  ratio_Powheg->GetXaxis()->SetTitleSize(0.2);
  ratio_Powheg->SetLabelSize(0.1,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();	///<change to ratio TPad
  ratio_Powheg->Sumw2();
  ratio_Powheg->SetStats(0);
  ratio_Mad->Sumw2();
  ratio_Mad->SetStats(0);
  ratio_Amc->Sumw2();
  ratio_Amc->SetStats(0);

  ratio_Powheg->Divide(hs_DYPowheg);
  ratio_Powheg->SetMarkerStyle(20);
  ratio_Powheg->SetMarkerColor(kRed);
  ratio_Powheg->SetLabelSize(0.1,"y");
  ratio_Powheg->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio_Powheg->GetYaxis()->SetNdivisions(505);

  ratio_Mad->Divide(hs_DYMadIncl);
  ratio_Mad->SetMarkerStyle(21);
  ratio_Mad->SetMarkerColor(kBlack);
  ratio_Mad->SetLabelSize(0.1,"y");
  ratio_Mad->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio_Mad->GetYaxis()->SetNdivisions(505);

  ratio_Amc->Divide(hs_DYAmcIncl);
  ratio_Amc->SetMarkerStyle(22);
  ratio_Amc->SetMarkerColor(kBlue);
  ratio_Amc->SetLabelSize(0.1,"y");
  ratio_Amc->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio_Amc->GetYaxis()->SetNdivisions(505);

  ratio_Mad->Draw("p");
  ratio_Amc->Draw("p");
  ratio_Powheg->Draw("p");
  float xmax = ratio_Powheg->GetXaxis()->GetXmax();
  float xmin = ratio_Powheg->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio_Powheg->Draw("p");
  ratio_Mad->Draw("psame");
  ratio_Amc->Draw("psame");
  f1->Draw("same");
  mycanvas->cd();

  TString fn = "";

  if(channel == Selector::EE)
    fn = "tempPlots/dyTagAndProbeEE/"+fname;
  if(channel == Selector::MuMu)
    //fn = "/publicweb/j/jchaves/WR/plots/miniTree/Selected/MuMuLowDilepton/"+fname;
    fn = "tempPlots/dyTagAndProbeMuMu/"+fname;

  mycanvas->Print((fn+".pdf").Data());
  mycanvas->Print((fn+".png").Data());
  p1->SetLogy();
  mycanvas->Print((fn+"_log.pdf").Data());
  mycanvas->Print((fn+"_log.png").Data());

  mycanvas->Close();
}
