#include "TH1F.h"
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


void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void miniPlotter(){

  TChain * chain_DY = new TChain("t1");
  chain_DY->Add("../selected_tree_DY_flavoursideband.root");
  TChain * chain_ttbar = new TChain("t1");
  chain_ttbar->Add("../selected_tree_ttbar_flavoursideband.root");
  TChain * chain_WZ = new TChain("t1");
  chain_DY->Add("../selected_tree_WZ_flavoursideband.root");
  TChain * chain_ZZ = new TChain("t1");
  chain_DY->Add("../selected_tree_ZZ_flavoursideband.root");

  TChain * chain_data = new TChain("t1");
  chain_data->Add("../selected_tree_data_EMu_flavoursideband.root");

  Selector myEvent_DY;
  Selector myEvent_ttbar;
  Selector myEvent_WZ;
  Selector myEvent_ZZ;
  Selector myEvent_data;

  myEvent_DY.SetBranchAddresses(chain_DY);
  myEvent_ttbar.SetBranchAddresses(chain_ttbar);
  myEvent_WZ.SetBranchAddresses(chain_WZ);
  myEvent_ZZ.SetBranchAddresses(chain_ZZ);
  myEvent_data.SetBranchAddresses(chain_data);

  std::vector<TH1F*> hs_DY;
  MakeHistos(chain_DY, &myEvent_DY, &hs_DY);
  std::vector<TH1F*> hs_ttbar;
  MakeHistos(chain_ttbar, &myEvent_ttbar, &hs_ttbar);
  std::vector<TH1F*> hs_WZ;
  MakeHistos(chain_WZ, &myEvent_WZ, &hs_WZ);
  std::vector<TH1F*> hs_ZZ;
  MakeHistos(chain_ZZ, &myEvent_ZZ, &hs_ZZ);

  std::vector<TH1F*> hs_data;
  MakeHistos(chain_data, &myEvent_data, &hs_data);

  //std::vector<TH1F*> hs_ttbar = MakeHistos(chain_ttbar, &myEvent_ttbar);

  unsigned int nPlots = hs_DY.size();

  //hs_ttbar[0]->Draw();
  
  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WZ[i],hs_ZZ[i],hs_data[i],xtitles[i],fnames[i]);
  }
  
}

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",50,0,2000);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",50,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",50,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",50,0,2000);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",50,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",50,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",50,0,2000);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",50,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",50,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",50,0,2000);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",50,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",50,-3.15,3.15);

  TH1F *h_WR_mass = new TH1F("h_WR_mass","",50,0,6000);
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",50,0,2000);

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);

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

}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.70 ) ; 
  leg->AddEntry( hs_DY, "DY" ) ; 
  leg->AddEntry( hs_ttbar, "ttbar" ) ;
  leg->AddEntry( hs_WZ, "WZ" ) ; 
  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
  //leg->AddEntry( histos[5][0], "WJets" ) ;  
  //leg->AddEntry( histos[2][0], "10 x WR 2600" ) ; 
  leg->AddEntry( hs_data, "Data");
  leg->SetFillColor( kWhite ) ; 


  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WZ->SetFillColor(kCyan);
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_DY);
  th->Add(hs_ttbar);
  th->Add(hs_WZ);
  th->Add(hs_ZZ);
  hs_data->SetMarkerStyle(20);

  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Preliminary");
  hs_data->SetTitle("CMS Preliminary");
  //th->Draw("histo");
  //hs_data->Draw("epsame");
  hs_data->Draw("ep");
  th->Draw("histo same");
  hs_data->Draw("epsame");
  TString ytitle = "Events/(";
  ytitle += (th->GetXaxis()->GetNbins());
  ytitle += ")";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());

  ratio->GetXaxis()->SetTitle(xtitle.Data());
  //ths[icanvas]->GetXaxis()->SetTickSize(1.0);
  //ths[icanvas]->GetXaxis()->SetTitleSize(0.1);
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(0.2);
  ratio->SetLabelSize(0.1,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();
  ratio->Sumw2();
  ratio->SetStats(0);

  hs_DY->Add(hs_ttbar);
  hs_DY->Add(hs_WZ);
  hs_DY->Add(hs_ZZ);

  ratio->Divide(hs_DY);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(0.1,"y");
  ratio->GetYaxis()->SetRangeUser(0.5,2.5);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();


  //TString fn = "plots/"+fname;
  //TString fn = "~/www/plots/WR/miniTTree/Selected/Flavor/"+fname;
  TString fn = "/publicweb/j/jchaves/WR/plots/miniTree/Selected/Flavor/"+fname;

  mycanvas->Print((fn+".pdf").Data());
  mycanvas->Print((fn+".png").Data());
  p1->SetLogy();
  mycanvas->Print((fn+"_log.pdf").Data());
  mycanvas->Print((fn+"_log.png").Data());

  mycanvas->Close();
}
