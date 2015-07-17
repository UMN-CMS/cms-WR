#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
using namespace std;

void Fill_Histo(TH1F* h1, TTree* tree);

void plot(){
  TFile * hfile0 = new TFile("~/work/WR_skims/skim_ttree_dyjets.root");
  TFile * hfile1 = new TFile("~/work/WR_skims/skim_ttree_ttbar.root");
  TFile * hfile2 = new TFile("~/work/WR_skims/skim_ttree_2600.root");

  TTree *tree0 = (TTree*)hfile0->Get("MakeTTree_Muons/t");  
  TTree *tree1 = (TTree*)hfile1->Get("MakeTTree_Muons/t");  
  TTree *tree2 = (TTree*)hfile2->Get("MakeTTree_Muons/t");  

  TH1F *h_Mlljj_0 = new TH1F("h_Mlljj_0","",100,0,4000);
  TH1F *h_Mlljj_1 = new TH1F("h_Mlljj_1","",100,0,4000);
  TH1F *h_Mlljj_2 = new TH1F("h_Mlljj_2","",100,0,4000);

  TH1F *h_Mll_0 = new TH1F("h_Mll_0","",100,0,3000);
  TH1F *h_Mll_1 = new TH1F("h_Mll_1","",100,0,3000);
  TH1F *h_Mll_2 = new TH1F("h_Mll_2","",100,0,3000);

  THStack *th0 = new THStack("th0","");
  THStack *th1 = new THStack("th1","");

  std::vector<TH1F*> histos0(2); // DY
  std::vector<TH1F*> histos1(2); // TTbar
  std::vector<TH1F*> histos2(2); // WR
  std::vector<THStack*> ths(2); // Stacks
  
  histos0[0] = h_Mlljj_0;
  histos0[1] = h_Mll_0;
  histos1[0] = h_Mlljj_1;
  histos1[1] = h_Mll_1;
  histos2[0] = h_Mlljj_2;
  histos2[1] = h_Mll_2;
  ths[0] = th0;
  ths[1] = th1;

  Fill_Histo(histos0,tree0);
  Fill_Histo(histos1,tree1);
  Fill_Histo(histos2,tree2);

  // Scale = xsection*luminosity/events
  for(std::vector<TH1F*>::size_type i = 0; i != histos0.size(); i++){
    histos0[i]->Scale(6104*18.825/28825132);
    histos0[i]->SetFillColor(kYellow);
    histos1[i]->Scale(670.3*18.825/42730273);
    histos1[i]->SetFillColor(kGreen);
    histos2[i]->Scale(10*0.0142*18.825/50000);
    histos2[i]->SetLineColor(kRed);
    ths[i]->Add(histos0[i]);
    ths[i]->Add(histos1[i]);
  }
  
  TCanvas* mycanvas0 = new TCanvas( "mycanvas0", "", 0, 0, 600, 400 ) ;
  ths[0]->Draw();
  ths[0]->GetXaxis()->SetTitle("Mlljj");
  histos2[0]->Draw("same");
  mycanvas0->Print("~/www/plots/WR/skimmed/Mlljj.pdf");
  mycanvas0->Print("~/www/plots/WR/skimmed/Mlljj.png");
  TCanvas* mycanvas1 = new TCanvas( "mycanvas1", "", 0, 0, 600, 400 ) ;
  ths[1]->Draw();
  ths[1]->GetXaxis()->SetTitle("Mll");
  histos2[1]->Draw("same");
  mycanvas1->Print("~/www/plots/WR/skimmed/Mll.pdf");
  mycanvas1->Print("~/www/plots/WR/skimmed/Mll.png");
  

  TLegend *leg = new TLegend( 0.78, 0.50, 0.98, 0.65 ) ; 
  leg->AddEntry( histos0[0], "DY" ) ; 
  leg->AddEntry( histos1[0], "ttbar" ) ; 
  leg->AddEntry( histos2[0], "10 x WR 2600" ) ; 
  leg->SetFillColor( kWhite ) ; 
  leg->Draw(); 

}


void Fill_Histo(std::vector<TH1F*> h1, TTree* tree){  
  
  int nentries = tree->GetEntries();
  Float_t Mlljj;
  Float_t l1_pt;
  Float_t l2_pt;
  Float_t j1_pt;
  Float_t j2_pt;
  Float_t Mll;
  Float_t dR_l1j1;
  Float_t dR_l1j2;
  Float_t dR_l2j1;
  Float_t dR_l2j2;
  tree->SetBranchAddress("Mlljj",&Mlljj);
  tree->SetBranchAddress("leading_lepton_pt",&l1_pt);
  tree->SetBranchAddress("subleading_lepton_pt",&l2_pt);
  tree->SetBranchAddress("leading_jet_pt",&j1_pt);
  tree->SetBranchAddress("subleading_jet_pt",&j2_pt);
  tree->SetBranchAddress("dilepton_mass",&Mll);
  tree->SetBranchAddress("dilepton_mass",&Mll);
  tree->SetBranchAddress("dR_leadLepton_leadJet",&dR_l1j1);
  tree->SetBranchAddress("dR_leadLepton_subleadJet",&dR_l1j2);
  tree->SetBranchAddress("dR_subleadLepton_leadJet",&dR_l2j1);
  tree->SetBranchAddress("dR_subleadLepton_subleadJet",&dR_l2j2);
  for (Int_t ev = 0; ev < nentries; ev++) {
    tree->GetEntry(ev);
    if(Mlljj>0 && l1_pt>60 && l2_pt>40 && j1_pt>40 && j2_pt>40 && Mll>200 && dR_l1j1 > 0.4 && dR_l1j2 > 0.4 && dR_l2j1 > 0.4 && dR_l2j2 > 0.4)
      {
	h1[0]->Fill(Mlljj);
	h1[1]->Fill(Mll);
      }
  }

}
