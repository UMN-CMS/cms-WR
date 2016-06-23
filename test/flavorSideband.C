#include "TStyle.h"
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
#include <string>
// #include "ExoAnalysis/cmsWR/src/Selector.cc"
// #include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>


void fillHisto(TChain * chain, Selector *myEvent, TH1F * h);
void flavorSideband(){

  gStyle->SetOptFit(1111);
  TChain * chain_EMu = new TChain("Tree_Iter0");
  TChain * chain_EE = new TChain("Tree_Iter0");
  TChain * chain_MuMu = new TChain("Tree_Iter0");
  TChain * chain_EMuData = new TChain("Tree_Iter0");
 
  TString dir = "../";
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

  //TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",50,200,4000);
  //TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",50,200,4000);
  //TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",50,200,4000);
  TH1F *h_WR_mass_EMu = new TH1F("h_WR_mass_EMu","",50,200,2000);
  TH1F *h_WR_mass_EE = new TH1F("h_WR_mass_EE","",50,200,2000);
  TH1F *h_WR_mass_MuMu = new TH1F("h_WR_mass_MuMu","",50,200,2000);
  TH1F *h_WR_mass_EMuData = new TH1F("h_WR_mass_EMuData","",50,200,2000);
  
  fillHisto(chain_EMu, &myEvent_EMu, h_WR_mass_EMu);
  fillHisto(chain_EMuData, &myEvent_EMuData, h_WR_mass_EMuData);
  fillHisto(chain_EE, &myEvent_EE, h_WR_mass_EE);
  fillHisto(chain_MuMu, &myEvent_MuMu, h_WR_mass_MuMu);

  TCanvas* mycanvas_EE = new TCanvas( "mycanvas_EE", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_EE->SetLineColor(kRed);
  h_WR_mass_EE->DrawNormalized("same histo");
  TLegend *leg_EE = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_EE->AddEntry( h_WR_mass_EMu, "EMu" );
  leg_EE->AddEntry( h_WR_mass_EE, "EE" );
  leg_EE->Draw();
  mycanvas_EE->Print(("flavor_EE.pdf"));
  mycanvas_EE->Print(("flavor_EE.png"));

  TCanvas* mycanvas_MuMu = new TCanvas( "mycanvas_MuMu", "", 0, 0, 600, 600 ) ;
  h_WR_mass_EMu->DrawNormalized();
  h_WR_mass_MuMu->SetLineColor(kRed);
  h_WR_mass_MuMu->DrawNormalized("same histo");
  TLegend *leg_MuMu = new TLegend( 0.72, 0.50, 0.98, 0.70 );
  leg_MuMu->AddEntry( h_WR_mass_EMu, "EMu" );
  leg_MuMu->AddEntry( h_WR_mass_MuMu, "MuMu" );
  leg_MuMu->Draw();
  mycanvas_MuMu->Print(("flavor_MuMu.pdf"));
  mycanvas_MuMu->Print(("flavor_MuMu.png"));

  TH1F *h_ratio_EE = (TH1F*)h_WR_mass_EE->Clone();
  TH1F *h_ratio_MuMu = (TH1F*)h_WR_mass_MuMu->Clone();
  h_ratio_EE->Divide(h_WR_mass_EMu);
  h_ratio_MuMu->Divide(h_WR_mass_EMu);
  TCanvas* mycanvas_ratio_EE = new TCanvas( "mycanvas_ratio_EE", "", 0, 0, 600, 600 ) ;
  TF1 *f_EE = new TF1("f_EE","[0]",600,2000);
  f_EE->FixParameter(0,0.427);
  h_ratio_EE->Fit("f_EE");
  h_ratio_EE->Draw();
  f_EE->SetLineColor(kBlue);
  f_EE->Draw("same");
  mycanvas_ratio_EE->Print(("flavor_ratio_EE.pdf"));
  mycanvas_ratio_EE->Print(("flavor_ratio_EE.png"));


  TCanvas* mycanvas_ratio_MuMu = new TCanvas( "mycanvas_ratio_MuMu", "", 0, 0, 600, 600 ) ;
  TF1 *f_MuMu = new TF1("f_MuMu","[0]",600,2000);
  f_MuMu->FixParameter(0,0.64);
  h_ratio_MuMu->Fit("f_MuMu");
  h_ratio_MuMu->Draw();
  f_MuMu->SetLineColor(kBlue);
  f_MuMu->Draw("same");
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu.pdf"));
  mycanvas_ratio_MuMu->Print(("flavor_ratio_MuMu.png"));

  TCanvas* canvMuMuEMu = new TCanvas("canvMuMuEMu","",600,600);
  canvMuMuEMu->cd();
  TLegend * legMuMuEMu = new TLegend(0.72,0.6,0.98,0.8);
  legMuMuEMu->AddEntry(h_WR_mass_EMu,"EMu");
  legMuMuEMu->AddEntry(h_WR_mass_MuMu,"MuMu");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_MuMu->Draw("Psame");
  legMuMuEMu->Draw();
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region.pdf","recreate");
  canvMuMuEMu->SaveAs("emujj_and_mumujj_signal_region.png","recreate");

  TCanvas* canvMuMuEMuData = new TCanvas("canvMuMuEMuData","",600,600);
  canvMuMuEMuData->cd();
  TLegend * legMuMuEMuData = new TLegend(0.72,0.6,0.98,0.8);
  legMuMuEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legMuMuEMuData->AddEntry(h_WR_mass_MuMu,"MuMu MC");
  h_WR_mass_MuMu->Draw("histo");
  h_WR_mass_EMuData->Scale(0.64);
  h_WR_mass_EMuData->SetMarkerStyle(2);
  h_WR_mass_EMuData->SetMarkerSize(2);
  h_WR_mass_EMuData->Draw("Psame");
  legMuMuEMuData->Draw();
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region.pdf","recreate");
  canvMuMuEMuData->SaveAs("rescaled_emujj_data_and_mumujj_MC_signal_region.png","recreate");


  TCanvas* canvEEEMu = new TCanvas("canvEEEMu","",600,600);
  canvEEEMu->cd();
  TLegend * legEEEMu = new TLegend(0.72,0.6,0.98,0.8);
  legEEEMu->AddEntry(h_WR_mass_EMu,"EMu");
  legEEEMu->AddEntry(h_WR_mass_EE,"EE");
  h_WR_mass_EMu->Draw("histo");
  h_WR_mass_EE->Draw("Psame");
  legEEEMu->Draw();
  canvEEEMu->SaveAs("emujj_and_eejj_signal_region.pdf","recreate");
  canvEEEMu->SaveAs("emujj_and_eejj_signal_region.png","recreate");

  TCanvas* canvEEEMuData = new TCanvas("canvEEEMuData","",600,600);
  canvEEEMuData->cd();
  TLegend * legEEEMuData = new TLegend(0.72,0.6,0.98,0.8);
  legEEEMuData->AddEntry(h_WR_mass_EMuData,"Rescaled EMu Data");
  legEEEMuData->AddEntry(h_WR_mass_EE,"EE MC");
  h_WR_mass_EE->Draw("histo");
  h_WR_mass_EMuData->Scale(1/0.64);	///<undo the scaling which was done earlier
  h_WR_mass_EMuData->SetMarkerStyle(2);
  h_WR_mass_EMuData->SetMarkerSize(2);
  h_WR_mass_EMuData->Scale(0.427);
  h_WR_mass_EMuData->Draw("Psame");
  legEEEMuData->Draw();
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region.pdf","recreate");
  canvEEEMuData->SaveAs("rescaled_emujj_data_and_eejj_MC_signal_region.png","recreate");

  TFile f("flavor_fits.root","RECREATE");
  h_ratio_EE->Write();
  h_ratio_MuMu->Write();
  f_EE->Write();
  f_MuMu->Write();

}

void fillHisto(TChain * chain, Selector *myEvent, TH1F * h){

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
    if(myEvent->WR_mass > 600. && myEvent->dilepton_mass > 200.) 
      h->Fill(myEvent->WR_mass,myEvent->weight);
  }
  //std::cout<<"histo named\t"<< h->GetName() <<"\thas integral\t"<< h->Integral() <<std::endl;
}
