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
#include "../src/Selector.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>

//DEBUG draws a fitted blue curve onto the RooDataSet in a second graph
//by default one plot of the fitted data will be shown in a binned histogram, with a TF1 fit applied to the data
#define DEBUG


void combinedMiniFitPlot(){
  
  TString fname = "../selected_tree_MODE_RGN.root";

  TFile *f = new TFile(fname);
  RooDataSet* data = (RooDataSet*)f->Get("permanentWeightedDataSet");
  TTree* t = (TTree*)f->Get("Tree_Iter0");
  TTree* t_central = (TTree*)f->Get("central_value_tree");
  
  float xmin = 600;
  float xmax = 6500;

  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 );
  RooRealVar x("fourObjectMass","fourObjectMass", xmin, xmax);
  RooPlot* xframe = x.frame(200);	///<number of emulated bins shown in RooPlot
  xframe->SetTitle("Unbinned Permanent RooDataSet");
#ifdef DEBUG
  RooRealVar expPower("expPower", "expPower", -1., -0.00005);
  RooExponential expPdf("rescaledExpPdf", "", x, expPower);
  expPower.setVal(-0.004);
  expPdf.fitTo(*data, RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE), RooFit::PrintLevel(-1));
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
    if(WR_mass>600.)
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

  mycanvas->Print("validationPlots/Mlljj_unbinned_dataset_MODE_RGN.pdf");
  mycanvas->Print("validationPlots/Mlljj_unbinned_dataset_MODE_RGN.png");
  mycanvas_mini->Print("validationPlots/Mlljj_binned_dataset_MODE_RGN.pdf");
  mycanvas_mini->Print("validationPlots/Mlljj_binned_dataset_MODE_RGN.png");

}
