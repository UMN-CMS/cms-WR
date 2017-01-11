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

//#define fixedBinWidths
#define variableBinWidths

//switch btwn DY AMCNLO and DY MADHT samples
#define DOAMC

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

//to change from lowdilepton to lowfourobj, simply search for all instances of lowdilepton, and replace them with lowfourobj
//to change from MuMu to EE, switch the Selector channel in the next line
Selector::tag_t channel = Selector::EE;

/**
 * this macro is designed to read several TChains, representing data and MC, apply no cuts, and plot
 * data as points and all MC as one stacked histogram with several fill colors.
 *
 * This macro should be used by itself on minitrees which have been processed with analysis.cpp, -c MuMu
 * or MuMu, and --isLowDiLepton true.
 *
 */

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void quickMiniPlotterNonTandPSidebands(){

  TChain * chain_DY = new TChain("Tree_Iter0","DY");
  TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
  //TChain * chain_ttbar = new TChain("Tree_Iter0","TTData");
  TChain * chain_WJets = new TChain("Tree_Iter0","WJets");
  TChain * chain_WZ = new TChain("Tree_Iter0","WZ");
  TChain * chain_ZZ = new TChain("Tree_Iter0","ZZ");
  TChain * chain_data = new TChain("Tree_Iter0","Data");

  TString localDir = "../analysisCppOutputRootFiles/";
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0;
  switch (channel) {
	  case Selector::MuMu:
#ifdef DOAMC
		  dy = chain_DY->Add(localDir+"selected_tree_DYAMC_lowdileptonsidebandMuMu_withMllWeight.root");
#endif
#ifndef DOAMC
		  dy = chain_DY->Add(localDir+"selected_tree_DYMADHT_lowdileptonsidebandMuMu_withMllWeight.root");
#endif
		  
		  tt = chain_ttbar->Add(localDir+"selected_tree_TT_lowdileptonsidebandMuMu.root");
		  //tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMu.root");
		  
		  wjets = chain_WJets->Add(localDir+"selected_tree_W_lowdileptonsidebandMuMu.root");
		  wz = chain_WZ->Add(localDir+"selected_tree_WZ_lowdileptonsidebandMuMu.root");
		  zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_lowdileptonsidebandMuMu.root");
		  data = chain_data->Add(localDir+"selected_tree_data_lowdileptonsidebandMuMu.root");
		  break;
	  case Selector::EE:
#ifdef DOAMC
		  dy = chain_DY->Add(localDir+"selected_tree_DYAMC_lowdileptonsidebandEE_withMllWeight.root");
#endif
#ifndef DOAMC
		  dy = chain_DY->Add(localDir+"selected_tree_DYMADHT_lowdileptonsidebandEE_withMllWeight.root");
#endif
		  
		  tt = chain_ttbar->Add(localDir+"selected_tree_TT_lowdileptonsidebandEE.root");
		  //tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMu.root");
		  
		  wjets = chain_WJets->Add(localDir+"selected_tree_W_lowdileptonsidebandEE.root");
		  wz = chain_WZ->Add(localDir+"selected_tree_WZ_lowdileptonsidebandEE.root");
		  zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_lowdileptonsidebandEE.root");
		  data = chain_data->Add(localDir+"selected_tree_data_lowdileptonsidebandEE.root");
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

  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","nPU"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","nPU","unweightedMLLJJ"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WJets[i],hs_WZ[i],hs_ZZ[i],hs_data[i],xtitles[i],fnames[i]);
  }
  
}//end quickMiniPlotterNonTandPSidebands()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  Int_t nbins = 20;
  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",nbins,0,500);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",nbins,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",nbins,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",nbins,0,500);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",nbins,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",nbins,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",nbins,0,500);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",nbins,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",nbins,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",nbins,0,500);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",nbins,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",nbins,-3.15,3.15);

  TH1F *h_WR_mass;
  TH1F *h_WR_mass_unweighted;
  
#ifdef fixedBinWidths
  h_WR_mass = new TH1F("h_WR_mass","",nbins,550,1050);
  h_WR_mass_unweighted = new TH1F("h_WR_mass_unweighted","",nbins,550,1050);
#endif

#ifdef variableBinWidths
  //variable bin widths only for WR mass plot
  //Float_t bins[] = { 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1375, 1450, 1550, 1680, 1900, 2500};	//original
  Float_t bins[] = { 150, 200, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1900, 2500};	//wider bins work better at high WR mass, include overflow evts in last bin shown on plot
  //Float_t bins[] = { 150, 200, 250, 300, 350, 400, 450, 525, 600, 675, 755, 850, 950, 1050, 1150, 1250, 1350, 1510, 1640, 1900, 2500, 6000};	//wider bins, out to 6000 GeV
  
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  h_WR_mass = new TH1F("h_WR_mass","",binnum, bins);
  h_WR_mass_unweighted = new TH1F("h_WR_mass_unweighted","",binnum, bins);
#endif 

  float dilepton_max = 210.;
  if(channel == Selector::EMu)
    dilepton_max = 1000;
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",nbins,50,dilepton_max);
  TH1F *h_nPV = new TH1F("h_nPV","",nbins,0,50);
  TH1F *h_nPU = new TH1F("h_nPU","",nbins,0,50);

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  TString chainTitle(chain->GetTitle());
  Float_t ttScaleFactor = 1.0;
  if( chainTitle.EqualTo("TTMC") ){
	  ttScaleFactor = (channel == Selector::MuMu) ? 0.958 : 0.954;	//to account for slightly higher number of rescaled ttbar MC events relative to rescaled emu data evts
  }
  if( chainTitle.EqualTo("TTData") ){
	  ttScaleFactor = (channel == Selector::MuMu) ? 0.657 : 0.414;	//to rescale emu data evts to estimates of ttbar in electron and muon channels
  }

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	//if(myEvent->dR_leadlepton_leadjet < 1.2 || myEvent->dR_leadlepton_subleadjet < 1.2 || myEvent->dR_subleadlepton_leadjet < 1.2 || myEvent->dR_subleadlepton_subleadjet < 1.2) continue;
	//if(myEvent->lead_jet_pt < 100 || myEvent->sublead_jet_pt < 100) continue;


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
    h_nPU->Fill(myEvent->nPU,(myEvent->weight)*ttScaleFactor);
  
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
  hs->push_back(h_nPU);
  hs->push_back(h_WR_mass_unweighted);

#ifdef variableBinWidths
  //normalize histo bins
  unsigned int max = hs->size();
  for(unsigned int i=0; i<max; i++){
	  //get num bins in histo i
	  Int_t nBins = (hs->at(i))->GetNbinsX();
	  for(Int_t j=1; j<=nBins; j++){
		  //include the overflows in the very last bin shown on the plot
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

  }//end loop over histos in vector
#endif

}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  hs_DY->Sumw2();
  hs_ttbar->Sumw2();
  hs_WJets->Sumw2();
  hs_WZ->Sumw2();
  hs_ZZ->Sumw2();
  hs_data->Sumw2();

  TLegend *leg = new TLegend( 0.6, 0.60, 0.90, 0.90 ) ; 
#ifdef DOAMC
  leg->AddEntry( hs_DY, "DY AMCNLO" ) ; 
#endif
#ifndef DOAMC
  leg->AddEntry( hs_DY, "DY MadgraphHT" ) ; 
#endif
  leg->AddEntry( hs_ttbar, "TTbar MC" ) ;
  leg->AddEntry( hs_WJets, "WJets" ) ; 
  leg->AddEntry( hs_WZ, "WZ" ) ; 
  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
  //leg->AddEntry( histos[2][0], "10 x WR 2600" ) ; 
  leg->AddEntry( hs_data, "Data");
  leg->SetFillColor( kWhite ) ; 


  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
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
  th->Add(hs_ttbar);
  th->Add(hs_DY);
  hs_data->SetMarkerStyle(20);

  /*for ratio plot*/
  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  /**/
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->SetTitle("CMS Private   #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  th->Draw("e1histo");
  hs_data->Draw("epsame");
  //TString ytitle = "Events/(";
  //ytitle += (th->GetXaxis()->GetNbins());
  //ytitle += ")";
  TString ytitle = "Events";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")){
	  th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetYaxis()->SetTitle("Events/GeV"), hs_data->GetYaxis()->SetTitle("Events/GeV");
	  Int_t maxBin = hs_data->GetNbinsX();
  }
  th->GetYaxis()->SetTitleOffset(1.2);
  Float_t labelSize = 0.25;
 
  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")) ratio->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(labelSize);
  ratio->SetLabelSize(labelSize - 0.07,"x");
  leg->Draw(); 
  mycanvas->cd();
  p2->cd();	//for ratio plot
  ratio->Sumw2();
  ratio->SetStats(0);

  hs_DY->Add(hs_ttbar);
  hs_DY->Add(hs_WJets);
  hs_DY->Add(hs_WZ);
  hs_DY->Add(hs_ZZ);
  if(fname.EqualTo("Mll") || fname.EqualTo("Mlljj") ){
	  Float_t dataMCratio = (hs_data->Integral()/hs_DY->Integral());
	  Float_t dataEntries = hs_data->GetEntries();
	  Float_t mcEntries = (hs_DY->GetEntries()) + (hs_ttbar->GetEntries()) + (hs_WJets->GetEntries()) + (hs_WZ->GetEntries()) + (hs_ZZ->GetEntries());
	  Float_t integralUnc = (dataEntries/mcEntries)*sqrt((1/dataEntries) + (1/mcEntries));
	  if(channel == Selector::MuMu) std::cout<< "in MuMu channel "<< fname <<" dataOvrMC ratio=\t"<< dataMCratio <<"\t+/-\t"<< integralUnc << std::endl;
  	  if(channel == Selector::EE) std::cout<< "in EE channel "<< fname <<" dataOvrMC ratio=\t"<< dataMCratio <<"\t+/-\t"<< integralUnc << std::endl;
  }

#ifdef variableBinWidths  
  //////////////////////////variable bin widths
  //use this when bin widths are variable
  if(fname.EqualTo("Mlljj") ){
	  std::cout<<"\t"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  std::cout<<"\t"<<std::endl;
	  std::cout<<"MC evts\tData evts"<<std::endl;
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<< (hs_DY->GetBinContent(i))*(hs_DY->GetBinWidth(i)) <<"\t" << (hs_data->GetBinContent(i))*(hs_data->GetBinWidth(i)) <<std::endl;
	  }//end loop over bins in MLLJJ
	  std::cout<<"\t"<<std::endl;
  }
	
  if(fname.EqualTo("unweightedMLLJJ") ){
	  std::cout<<"\t"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  std::cout<<"\t"<<std::endl;
	  std::cout<<"unweighted MC evts"<<std::endl;
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<< (hs_DY->GetBinContent(i))*(hs_DY->GetBinWidth(i)) << std::endl;
	  }//end loop over bins in unweightedMLLJJ
	  
	  std::cout<<"\t"<<std::endl;
  }
  ////////////////////variable bin widths
#endif

#ifdef fixedBinWidths
  //////////////////////////constant bin widths
  //use this when bin widths are constant
  if(fname.EqualTo("Mlljj") ){
	  std::cout<<"\t"<<std::endl;
	  std::cout<<"MC evts\tData evts"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<< hs_DY->GetBinContent(i) <<"\t"<< hs_data->GetBinContent(i) <<std::endl;
	  }//end loop over bins in MLLJJ
	  std::cout<<"\t"<<std::endl;
  }

/**/  
  if(fname.EqualTo("unweightedMLLJJ") ){
	  std::cout<<"\t"<<std::endl;
	  std::cout<<"unweighted MC evts"<<std::endl;
	  Int_t nbins = hs_data->GetNbinsX();
	  for(Int_t i = 1; i<=nbins; i++){
		  std::cout<< hs_DY->GetBinContent(i) <<std::endl;
	  }//end loop over bins in unweightedMLLJJ
	  
	  std::cout<<"\t"<<std::endl;
  }
  /**/
  ////////////////////constant bin widths
#endif



  ratio->Divide(hs_DY);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(labelSize - 0.07,"y");
  Float_t maxYratioRange = 2.95;
  if(channel == Selector::EE) maxYratioRange = 1.9;
  ratio->GetYaxis()->SetRangeUser(0.51,maxYratioRange);	//upper limit 3. for mu, 2. for ele
  ratio->GetYaxis()->SetNdivisions(505);
  /*for ratio plot*/
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  /**/

  TString fn = "";
  //fn = fname + "_100GeVbinsFromMLLJJ700_variablebinwidths_rescaledEMuData_lowdileptonMuMuChannelDyAmc";	//for ratio plot
 
#ifdef variableBinWidths  
  fn = fname + "_100GeVbinsFromMLLJJ700_variablebinwidthsLastBinIncludesOverflow_rescaledTTBarMC_lowdileptonMuMuChannel";	//for ratio plot
  if(channel == Selector::EE) fn = fname + "_100GeVbinsFromMLLJJ700_variablebinwidthsLastBinIncludesOverflow_rescaledTTBarMC_lowdileptonEEChannel";
#endif

#ifdef fixedBinWidths
  fn = fname + "_fixedbinwidths_rescaledTTBarMC_lowdileptonMuMuChannel";	//for ratio plot
  if(channel == Selector::EE) fn = fname + "_fixedbinwidths_rescaledTTBarMC_lowdileptonEEChannel";
#endif 

#ifdef DOAMC
  fn += "DYAMC";
#endif
#ifndef DOAMC
  fn += "DYMADHT";
#endif

  //fn = fname + "_variablebinwidths_rescaledEMuData_lowdileptonMuMuChannelDyMadHt";	//for ratio plot
  //fn = fname + "_noRatio_variablebinwidths_lowdileptonMuMuChannelDyAmc";	//only needed when no ratio plot is drawn
  //fn = fname + "_100GeVbinsFromMLLJJ700_variablebinwidths_onlyDY_MLLJJbtwn700and1400_lowdileptonMuMuChannelDyMadHt";	//for ratio plot
  //fn = fname + "_100GeVbinsFromMLLJJ700_variablebinwidths_onlyDY_MLLJJbtwn700and1400_lowdileptonMuMuChannelDyAmc";	//for ratio plot

  if(fname.EqualTo("Mlljj") ){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  p1->SetLogy();	//for ratio plot
	  //mycanvas->SetLogy();	//only needed when no ratio plot is drawn
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
  }

  mycanvas->Close();
}
