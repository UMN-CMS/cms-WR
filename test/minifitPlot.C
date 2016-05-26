#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "RooPlot.h"
#include <vector>
#include <iostream>
#include <string>
// #include "ExoAnalysis/cmsWR/src/Selector.cc"
// #include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>

#define DEBUG
//Selector::tag_t channel = Selector::MuMu;
//Selector::tag_t channel = Selector::EE;
Selector::tag_t channel = Selector::EMu;


void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void minifitPlot(){
  
  TString fname = "";

  switch (channel) {
  case Selector::EE:
    fname = "../selected_tree_DYAMC_signal_eeEE.root";
    //fname = "../selected_tree_DYPOWHEG_signal_eeEE.root";
    break;
  case Selector::MuMu:
    fname = "../selected_tree_DYAMC_signal_mumuMuMu.root";
    break;
  case Selector::EMu:    
    //fname = "../selected_tree_TT_flavoursidebandEMu.root";
    fname = "../selected_tree_data_flavoursidebandEMu.root";
    break;
  default:
    std::cout << "Unknown tag" << std::endl;
  }

  TFile *f = new TFile(fname);
  RooDataSet* data = (RooDataSet*)f->Get("permanentWeightedDataSet");
  TTree* t = (TTree*)f->Get("Tree_Iter0");
  TTree* t_central = (TTree*)f->Get("central_value_tree");
  
  float xmin = 600;
  float xmax = 6500;

  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  RooRealVar x("fourObjectMass","fourObjectMass", xmin, xmax) ;
  RooPlot* xframe = x.frame(RooFit::Title("Permanent RooDataSet")) ;
#ifdef DEBUG
  RooRealVar expPower("expPower", "expPower", -1., -0.00005);
  RooExponential expPdf("rescaledExpPdf", "", x, expPower);
  expPower.setVal(-0.004);
  expPdf.fitTo(*data);
  expPower.Print();
  data->plotOn(xframe);
  expPdf.plotOn(xframe);
#endif
#ifndef DEBUG
  data->plotOn(xframe);
#endif
  xframe->Draw() ;
  TCanvas* mycanvas_mini = new TCanvas( "mycanvas_mini", "", 0, 0, 600, 600 ) ;
  TH1F* h = new TH1F("h","",100, xmin, xmax);
  h->Sumw2();
  Float_t WR_mass;
  Float_t weight;
  Float_t Normalization;
  Float_t fitParams[16];
  t->SetBranchAddress("WR_mass", &WR_mass);
  t->SetBranchAddress("weight", &weight);
  t_central->SetBranchAddress("Normalization", &Normalization);
  t_central->SetBranchAddress("FitParameters", &fitParams);

  int nentries = t->GetEntries();
  for(int i=0; i<nentries; i++){
    t->GetEntry(i);
    if(WR_mass>600)
      h->Fill(WR_mass,weight);
  }
  
  TF1* f_exp = new TF1("f_exp","[0]*exp([1]*x)",600,6500);

  t_central->GetEntry(0);
  h->Scale(1/h->Integral());
  h->Scale(Normalization);
  float rooNorm = fitParams[0]/(TMath::Exp(fitParams[0]*xmax)-TMath::Exp(fitParams[0]*xmin));
  int bin_width = (xmax-xmin)/100;
  f_exp->SetParameter(0,Normalization*rooNorm*bin_width);
  f_exp->SetParameter(1,fitParams[0]);
  h->Draw("ep");
  f_exp->SetLineColor(kRed);
  f_exp->Draw("same");

  if(channel == Selector::EMu){
    //mycanvas->Print("~/www/plots/fit_validation/perm_dataset_emu.pdf");
    //mycanvas_mini->Print("~/www/plots/fit_validation/mini_dataset_emu.pdf");
    mycanvas->Print("plots/fit_validation/perm_dataset_emu.pdf");
    mycanvas_mini->Print("plots/fit_validation/mini_dataset_emu.pdf");
  }
  if(channel == Selector::EE){
    mycanvas->Print("plots/fit_validation/perm_dataset_ee.pdf");
    mycanvas_mini->Print("plots/fit_validation/mini_dataset_ee.pdf");
  }
  if(channel == Selector::MuMu){
    mycanvas->Print("plots/fit_validation/perm_dataset_mumu.pdf");
    mycanvas_mini->Print("plots/fit_validation/mini_dataset_mumu.pdf");
  }
}
