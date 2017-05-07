#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TPaveText.h"
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
#include "../src/Selector.cc"
#include "../src/SelectorHist.cc"
#include "../src/miniTreeEvent.cc"
#include <cstdio>
#include <memory>

#define unblindedData
#define plotRatio
//#define showRescaledRunOneEEJJExcess
//#define showQCD //dont enable this and showRescaledRunOneEEJJExcess or showWR simultaneously
//#define showWR	//show the distribution for a WR signal mass point as a histogram drawn as a line

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

/**
 * this macro is designed to read several TChains, representing WR signal and bkgnd MC, apply no cuts, and plot
 * WR signal MC as a red curve and all bkgnd MC as one stacked histogram with several fill colors.
 *
 * This macro should be used by itself on minitrees processed by analysis.cpp with -c EE or -c MuMu
 * and signal region requirements (i.e. no additional strings passed to analysis.cpp on the cmd line).
 *
 */

//switch Selector tag here, and everything else will change accordingly
Selector::tag_t channel = Selector::MuMu;

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data,TH1F* hs_Other, TString xtitle, TString fname);
void quickSignalRegionMiniPlotter(){

  TChain * chain_DY = new TChain("Tree_Iter0","DYMC");
  //TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
  TChain * chain_ttbar = new TChain("Tree_Iter0","TTData");
  TChain * chain_WJets = new TChain("Tree_Iter0","WJets");
  TChain * chain_WZ = new TChain("Tree_Iter0","Diboson");
  TChain * chain_ZZ = new TChain("Tree_Iter0","Other");
  TChain * chain_Other = new TChain("Tree_Iter0","HighMassWR");
  TChain * chain_data = new TChain("Tree_Iter0","Data");

  TString localDir = "../analysisCppOutputRootFilesNewHltSf/";
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0, other=0;
  switch (channel) {
  case Selector::EE:
	dy = chain_DY->Add(localDir+"selected_tree_DYMadInclAndHT_signal_eeEE_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_eeEE.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_WInclAndHT_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
	other = chain_Other->Add(localDir+"selected_tree_WRtoEEJJ_3000_1500_signal_eeEE.root");

	//use the ZZ chain to show the WRtoEEJJ signal rescaled to match the normalization of the expected RunI EEJJ excess in 2015 operating conditions
	//show nothing if showRescaledRunOneEEJJExcess is not defined
#ifdef showRescaledRunOneEEJJExcess
	zz = chain_ZZ->Add(localDir+"selected_tree_WRtoEEJJ_2000_1000_signal_eeEE.root");
#endif

#ifdef showWR
	zz = chain_ZZ->Add(localDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
	//zz = chain_ZZ->Add(localDir+"selected_tree_WRtoEEJJ_800_400_signal_eeEE.root");
#endif

#ifdef showQCD
	zz = chain_ZZ->Add(localDir+"selected_tree_qcdData_signal_eeEE.root");
#endif

	//there must be a file linked with the data chain, otherwise this macro will not produce any plots.
#ifndef unblindedData
	data = chain_data->Add(localDir+"selected_tree_WRtoEEJJ_1000_500_signal_eeEE.root");
#endif

#ifdef unblindedData	
	data = chain_data->Add(localDir+"selected_tree_data_signal_eeEE.root");
#endif

	break;
  case Selector::MuMu:

	dy = chain_DY->Add(localDir+"selected_tree_DYMadInclAndHT_signal_mumuMuMu_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_mumuMuMu.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_WInclAndHT_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");
	other = chain_Other->Add(localDir+"selected_tree_WRtoMuMuJJ_3000_1500_signal_mumuMuMu.root");

#ifdef showQCD
    zz = chain_ZZ->Add(localDir+"selected_tree_qcdData_signal_mumuMuMu.root");
#endif

#ifdef showWR
	zz = chain_ZZ->Add(localDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
	//zz = chain_ZZ->Add(localDir+"selected_tree_WRtoMuMuJJ_800_400_signal_mumuMuMu.root");
#endif

#ifndef unblindedData	
    data = chain_data->Add(localDir+"selected_tree_WRtoMuMuJJ_1000_500_signal_mumuMuMu.root");
#endif

#ifdef unblindedData
	data = chain_data->Add(localDir+"selected_tree_data_signal_mumuMuMu.root");
#endif


	break;
  default:
    std::cout << "Unknown tag" << std::endl;
  }

  std::cout<<"data = "<< data <<"\tdy = "<< dy << std::endl;
  //if(data==0 || dy==0 || tt==0 || wjets==0 || wz==0 || zz==0 || other==0) exit(-1);

  Selector myEvent_DY;
  Selector myEvent_ttbar;
  Selector myEvent_WJets;
  Selector myEvent_WZ;
  Selector myEvent_ZZ;
  Selector myEvent_Other;
  Selector myEvent_data;

  myEvent_DY.SetBranchAddresses(chain_DY);
  myEvent_ttbar.SetBranchAddresses(chain_ttbar);
  myEvent_WJets.SetBranchAddresses(chain_WJets);
  myEvent_WZ.SetBranchAddresses(chain_WZ);
  myEvent_ZZ.SetBranchAddresses(chain_ZZ);
  myEvent_Other.SetBranchAddresses(chain_Other);
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
  std::vector<TH1F*> hs_Other;
  MakeHistos(chain_Other, &myEvent_Other, &hs_Other);


  std::vector<TH1F*> hs_data;
  MakeHistos(chain_data, &myEvent_data, &hs_data);

  unsigned int nPlots = hs_DY.size();

  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","Z p_{T}","M_{LLJJ} (GeV)","M_{LLJJ} (GeV)","M_{LLJJ} (GeV)","M_{LJJ} using lead lepton [GeV]","M_{LJJ} using sublead lepton [GeV]","M_{LLJJ}","M_{LJJ} using lead lepton [GeV]","M_{LJJ} using sublead lepton [GeV]"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","Z_pt","Mlljj_lowZpt","Mlljj_varBins","Mlljj_2012Bins","Mljj_leadLept","Mljj_subleadLept","Mlljj_thinBins","Mljj_leadLept_varBins","Mljj_subleadLept_varBins"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WJets[i],hs_WZ[i],hs_ZZ[i],hs_data[i],hs_Other[i],xtitles[i],fnames[i]);
  }
  
}//end quickSignalRegionMiniPlotter()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",14,0,700);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",50,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",50,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",10,0,500);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",50,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",50,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",14,0,700);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",50,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",50,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",14,0,700);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",50,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",50,-3.15,3.15);
  TH1F *h_WR_mass = new TH1F("h_WR_mass_fixedBinWidth","",17,600,4000);	//200 GeV wide bins. DO NOT change this binning
  TH1F *h_WR_mass_thinBins = new TH1F("h_WR_mass_thinBins_fixedBinWidth","",68,600,4000);	//50 GeV wide bins
  TH1F *h_Nu_mass_leadLept = new TH1F("h_Nu_mass_leadLept_fixedBinWidth","",11,200,2400);	//200 GeV wide bins
  TH1F *h_Nu_mass_subleadLept = new TH1F("h_Nu_mass_subleadLept_fixedBinWidth","",10,200,2200);	//200 GeV wide bins
 
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",50,150,2000);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  TH1F *h_Z_pt = new TH1F("h_Z_pt", "", 40, 0., 1200.);
  TH1F *h_WR_mass_lowZpt = new TH1F("h_WR_mass_lowZpt","",40,500,3000);

  Float_t runOneBins[] = { 600, 800, 1000, 1200, 1400, 1600, 1800, 2200, 4000};	//MLLJJ bins used in 2012 search
  Int_t  runOneBinnum = sizeof(runOneBins)/sizeof(Float_t) - 1;
  TH1F *h_WR_mass_2012bins = new TH1F("h_WR_mass_2012bins","", runOneBinnum, runOneBins);

  Float_t runTwoBins[] = { 600, 800, 1000, 1200, 1600, 2400};
  Int_t  runTwoBinnum = sizeof(runTwoBins)/sizeof(Float_t) - 1;
  TH1F *h_Nu_mass_leadLept_varBins = new TH1F("h_Nu_mass_leadLept_varBins","", runTwoBinnum, runTwoBins);
  TH1F *h_Nu_mass_subleadLept_varBins = new TH1F("h_Nu_mass_subleadLept_varBins","", runTwoBinnum, runTwoBins);


  Float_t bins[] = { 600, 750, 900, 1050, 1200, 1500, 1800, 2100, 4000};	//MLLJJ variable bins
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
  TH1F *h_WR_mass_varBins = new TH1F("h_WR_mass_varBins","",binnum, bins);

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  TString chainTitle(chain->GetTitle());
  Float_t scaleFactor = 1.0;
  if( chainTitle.EqualTo("TTMC") ){
	  scaleFactor = (channel == Selector::MuMu) ? 0.958 : 0.954;	//to account for slightly higher number of rescaled ttbar MC events relative to rescaled emu data evts
  }
  if( chainTitle.EqualTo("TTData") ){
	  scaleFactor = (channel == Selector::MuMu) ? 0.66 : 0.43;	//to rescale emu data evts to estimates of ttbar in ele and mu channels
  }


  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	if(myEvent->WR_mass < 600.) continue;
	if(myEvent->dilepton_mass < 200.) continue;
	if(myEvent->sublead_lepton_pt < 53.) continue;

	TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, zFourMom;
	leadLeptonFourMom.SetPtEtaPhiE(myEvent->lead_lepton_pt, myEvent->lead_lepton_eta, myEvent->lead_lepton_phi, myEvent->lead_lepton_pt);
	subleadLeptonFourMom.SetPtEtaPhiE(myEvent->sublead_lepton_pt, myEvent->sublead_lepton_eta, myEvent->sublead_lepton_phi, myEvent->sublead_lepton_pt);
	zFourMom = leadLeptonFourMom + subleadLeptonFourMom;

	h_Z_pt->Fill(zFourMom.Pt(), (myEvent->weight)*scaleFactor);

	h_lepton_pt0->Fill(myEvent->lead_lepton_pt,(myEvent->weight)*scaleFactor);
    h_lepton_pt1->Fill(myEvent->sublead_lepton_pt,(myEvent->weight)*scaleFactor);
    h_lepton_eta0->Fill(myEvent->lead_lepton_eta,(myEvent->weight)*scaleFactor);
    h_lepton_eta1->Fill(myEvent->sublead_lepton_eta,(myEvent->weight)*scaleFactor);
    h_lepton_phi0->Fill(myEvent->lead_lepton_phi,(myEvent->weight)*scaleFactor);
    h_lepton_phi1->Fill(myEvent->sublead_lepton_phi,(myEvent->weight)*scaleFactor);

    h_jet_pt0->Fill(myEvent->lead_jet_pt,(myEvent->weight)*scaleFactor);
    h_jet_pt1->Fill(myEvent->sublead_jet_pt,(myEvent->weight)*scaleFactor);
    h_jet_eta0->Fill(myEvent->lead_jet_eta,(myEvent->weight)*scaleFactor);
    h_jet_eta1->Fill(myEvent->sublead_jet_eta,(myEvent->weight)*scaleFactor);
    h_jet_phi0->Fill(myEvent->lead_jet_phi,(myEvent->weight)*scaleFactor);
    h_jet_phi1->Fill(myEvent->sublead_jet_phi,(myEvent->weight)*scaleFactor);
      
    h_WR_mass->Fill(myEvent->WR_mass,(myEvent->weight)*scaleFactor);
    h_dilepton_mass->Fill(myEvent->dilepton_mass,(myEvent->weight)*scaleFactor);
    h_nPV->Fill(myEvent->nPV,(myEvent->weight)*scaleFactor);
    if(zFourMom.Pt() < 150.) h_WR_mass_lowZpt->Fill(myEvent->WR_mass,(myEvent->weight)*scaleFactor);
	h_WR_mass_2012bins->Fill(myEvent->WR_mass,(myEvent->weight)*scaleFactor);
	h_WR_mass_varBins->Fill(myEvent->WR_mass,(myEvent->weight)*scaleFactor);

	TLorentzVector leadJetFourMom, subleadJetFourMom, leadLeptDijetFourMom, subleadLeptDijetFourMom;
	leadJetFourMom.SetPtEtaPhiE(myEvent->lead_jet_pt, myEvent->lead_jet_eta, myEvent->lead_jet_phi, myEvent->lead_jet_pt);
	subleadJetFourMom.SetPtEtaPhiE(myEvent->sublead_jet_pt, myEvent->sublead_jet_eta, myEvent->sublead_jet_phi, myEvent->sublead_jet_pt);
	leadLeptDijetFourMom = leadLeptonFourMom + leadJetFourMom + subleadJetFourMom;
	subleadLeptDijetFourMom = subleadLeptonFourMom + leadJetFourMom + subleadJetFourMom;

	//plot lepton + dijet mass using leading or subleading lepton
	h_Nu_mass_leadLept->Fill(leadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);
  	h_Nu_mass_subleadLept->Fill(subleadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);
	h_Nu_mass_leadLept_varBins->Fill(leadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);
  	h_Nu_mass_subleadLept_varBins->Fill(subleadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);
   
	h_WR_mass_thinBins->Fill(myEvent->WR_mass, (myEvent->weight)*scaleFactor);
  
  }

#ifdef showRescaledRunOneEEJJExcess
  if( chainTitle.EqualTo("Other") && channel == Selector::EE ){
	  //first reset the bin contents of all Meejj plots (fixed bin width which are wide or narrow, or variable bin width) to zero
	  Int_t numWideFixedBins = h_WR_mass->GetNbinsX();
	  for(Int_t n=1; n<=numWideFixedBins+1; n++){h_WR_mass->SetBinContent(n,0.0);}
	  Int_t numThinFixedBins = h_WR_mass_thinBins->GetNbinsX();
	  for(Int_t n=1; n<=numThinFixedBins+1; n++){h_WR_mass_thinBins->SetBinContent(n,0.0);}
	  Int_t numVarWidthBins = h_WR_mass_2012bins->GetNbinsX();
	  for(Int_t n=1; n<=numVarWidthBins+1; n++){h_WR_mass_2012bins->SetBinContent(n,0.0);}
	  
	  //then create a Gaussian (via TF1) with mean 2000, integral equal to 10.0, and 3sigma = 200
	  //the integral of this Gaussian over different ranges will be used to update the bin contents of each Meejj histogram
	  Float_t minX = 1800.;
	  Float_t maxDomain = 400.;
	  Float_t desiredTotalIntegral = 10.;	//this is equal to the integral of the gaussian TF1 over its entire domain
	  Float_t desiredSigma = (maxDomain/2)/3.;
	  Float_t gausRenorm = sqrt(2*3.14159);
	  TF1 *excess = new TF1("excess","gaus(0)",minX,minX+maxDomain);	//a Gaussian has 3 variable parameters, and the first one (norm) will be labeled parameter 0
	  excess->SetParameter(0, desiredTotalIntegral/(gausRenorm*desiredSigma));	//set normalization.
	  excess->SetParameter(1, minX + (maxDomain/2));	//set mean of Gaussian
	  excess->SetParameter(2, desiredSigma);	//set std deviation of Gaussian

	  ///just to show the shape of the Gaussian
	  TCanvas * cGaus = new TCanvas("cGaus","cGaus", 800, 800);
	  cGaus->cd();
	  //make a text box which shows the Gaussian parameters
	  TPaveText *fitBox = new TPaveText(2075,0.04,2200,0.06);
	  fitBox->SetFillColor(kWhite);
	  fitBox->AddText("mean = " + TString( (to_string(excess->GetParameter(1))).substr(0,6) ) );
	  fitBox->AddText(" ");
	  fitBox->AddText("#sigma = " + TString( (to_string(excess->GetParameter(2))).substr(0,4) ) );
	  fitBox->AddText(" ");
	  fitBox->AddText("#int gaus = " + TString( (to_string(excess->Integral(1800.,2200.))).substr(0,4) ) );
	  excess->GetXaxis()->SetTitle("M_{EEJJ} [GeV]");
	  excess->GetYaxis()->SetTitle("Events");
	  excess->GetYaxis()->SetTitleOffset(1.45);
	  excess->SetTitle("");
	  excess->Draw();
	  fitBox->Draw("same");
	  TString gausFileName = "gaussianEmulationOfRunOneExcess";
	  cGaus->Print((gausFileName+".png").Data());
	  cGaus->Print((gausFileName+".pdf").Data());
	  cGaus->Print((gausFileName+".C").Data());

	  //now integrate the Gaussian over the bin ranges, and update the histogram entries accordingly
	  //for all three histograms, loop over all bins and update the bin contents for bins with lower edge
	  //greater than 1799 GeV, and less than 2199 GeV
	  //wide, fixed width bins
	  for(Int_t n=1; n<=numWideFixedBins; n++){
		  Double_t binLowEdge = h_WR_mass->GetBinLowEdge(n);
		  if(binLowEdge > 1799. && binLowEdge < 2199.){
		  Double_t binWidth = h_WR_mass->GetBinWidth(n);
		  h_WR_mass->SetBinContent(n, excess->Integral(binLowEdge, binLowEdge+binWidth));
		  }//end bin selection
	  }//end loop over all bins

	  //wide, variable width bins
	  for(Int_t n=1; n<=numVarWidthBins; n++){
		  Double_t binLowEdge = h_WR_mass_2012bins->GetBinLowEdge(n);
		  if(binLowEdge > 1799. && binLowEdge < 2199.){
		  Double_t binWidth = h_WR_mass_2012bins->GetBinWidth(n);
		  h_WR_mass_2012bins->SetBinContent(n, excess->Integral(binLowEdge, binLowEdge+binWidth));
		  }//end bin selection
	  }//end loop over all bins
	  
	  //thin, fixed width bins
	  for(Int_t n=1; n<=numThinFixedBins; n++){
		  Double_t binLowEdge = h_WR_mass_thinBins->GetBinLowEdge(n);
		  if(binLowEdge > 1799. && binLowEdge < 2199.){
		  Double_t binWidth = h_WR_mass_thinBins->GetBinWidth(n);
		  h_WR_mass_thinBins->SetBinContent(n, excess->Integral(binLowEdge, binLowEdge+binWidth));
		  }//end bin selection
	  }//end loop over all bins
	  
  }//end selection of Other chain and electron channel
#endif


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
  hs->push_back(h_Z_pt);
  hs->push_back(h_WR_mass_lowZpt);
  hs->push_back(h_WR_mass_varBins);
  hs->push_back(h_WR_mass_2012bins);
  hs->push_back(h_Nu_mass_leadLept);
  hs->push_back(h_Nu_mass_subleadLept);
  hs->push_back(h_WR_mass_thinBins);
  hs->push_back(h_Nu_mass_leadLept_varBins);
  hs->push_back(h_Nu_mass_subleadLept_varBins);


  //h_Nu_mass_leadLept_varBins->Fill(leadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);
  //h_Nu_mass_subleadLept_varBins->Fill(subleadLeptDijetFourMom.M(), (myEvent->weight)*scaleFactor);

  //rescale bin contents by bin width, and add overflow evts to last bin in two WR mass histos with variable bin widths
  //also add overflow events to last bin in fixed-bin-width MLLJJ and M_LJJ plots
  unsigned int max = hs->size();
  for(unsigned int i=0; i<max; i++){
	  TString histName = (hs->at(i))->GetName();
	  //std::cout<<"hist name\t"<< histName <<std::endl;
	  if(histName.Contains("fixedBin")){
		  Int_t nBins = (hs->at(i))->GetNbinsX();
		  for(Int_t j=1; j<=nBins; j++){
			  //include the overflows in the very last bin shown on the plot
			  if(j==nBins){
				  Double_t origBinContents = (hs->at(i))->GetBinContent(j);
				  Double_t overflowContents = (hs->at(i))->GetBinContent(j+1);
				  std::cout<<"overflow contents =\t"<< overflowContents << std::endl;	//sanity check
				  (hs->at(i))->SetBinContent(j, origBinContents+overflowContents);
			  }//end work to include overflows in last bin shown on plot
		  }//end loop over bins in histo
	  }//end adding overflow evts to fixed bin width MLLJJ plot

	  if(histName.Contains("varBins") || histName.Contains("2012bins")){
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

	  }//end filter to select the hists with variable bin widths

  }//end loop over histos in vector
 
}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data,TH1F* hs_Other, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.55, 0.5, 0.90, 0.90 ) ; 
  
  //PAS plot legend entries
  leg->AddEntry( hs_DY, "Z/#gamma*+jets" ) ; 
  leg->AddEntry( hs_ttbar, "Top bkgnds from data" ) ;
  leg->AddEntry( hs_WJets, "W+jets" ) ; 
  leg->AddEntry( hs_WZ, "Diboson" ) ;
  leg->AddEntry(hs_Other, "WR signal M_{WR}=3.0 TeV");

  //AN plot legend entries
  //leg->AddEntry( hs_DY, "DY" ) ; 
  //leg->AddEntry( hs_ttbar, "Top Backgrounds Data Driven" ) ;
  //leg->AddEntry( hs_WJets, "W+jets" ) ; 
  //leg->AddEntry( hs_WZ, "Diboson" ) ; 

#ifdef showRescaledRunOneEEJJExcess
  if(channel == Selector::EE){
	  leg->AddEntry( hs_ZZ, "RunI Excess+Backgrounds");
	  //leg->AddEntry( (TObject*)0, "M_{WR}=2.0 TeV M_{Nu}=1.0 TeV","");	//R is illegible on drawn legend if it is a subscript of W
  }
#endif

#ifdef showQCD
  leg->AddEntry( hs_ZZ, "QCD from data");
#endif

#ifdef showWR
  leg->AddEntry( hs_ZZ, "WR signal");
  leg->AddEntry( (TObject*)0, "M_{WR}=2.2 TeV M_{Nu}=1.1 TeV","");
  //leg->AddEntry( (TObject*)0, "M_{WR}=0.8 TeV M_{Nu}=0.4 TeV","");
#endif

#ifndef unblindedData
  //leg->AddEntry(hs_data, "WR signal x 0.3");
  //leg->AddEntry( (TObject*)0, "M_{WR}=1.0 TeV M_{Nu}=0.5 TeV","");
#endif

#ifdef unblindedData 
  leg->AddEntry( hs_data, "Data");
#endif
  leg->AddEntry( (TObject*)0, "Overflows in last bin","");
  leg->SetFillColor( kWhite );

  //NOTE. if showRescaledRunOneEEJJExcess is defined, then the hs_ZZ histo should not be added to the THStack of all SM
  //backgrounds, nor should it be used in any way in the ratio plot

  hs_data->Sumw2();
  hs_ttbar->Sumw2();
  hs_WJets->Sumw2();
  hs_WZ->Sumw2();
  hs_ZZ->Sumw2();
  hs_Other->Sumw2();
  hs_DY->Sumw2();
  
  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
#ifndef plotRatio
  mycanvas->cd();	//only needed when no ratio plot is drawn
#endif
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WJets->SetFillColor(kBlue);
  hs_WZ->SetFillColor(kCyan);
#ifdef showQCD
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_ZZ); 	//add QCD data driven bkgnd
#endif
  th->Add(hs_WZ);
  th->Add(hs_WJets);
  th->Add(hs_DY);
  th->Add(hs_ttbar);
  hs_data->SetMarkerStyle(20);
  hs_data->SetMarkerSize(1);
  hs_data->SetMarkerColor(kBlack);
#ifdef showRescaledRunOneEEJJExcess
  //add all of the expected backgrounds to the Run1 excess which has been rescaled based on integrated lumi and WR cross sxn
  hs_ZZ->Add(hs_DY);
  hs_ZZ->Add(hs_ttbar);
  hs_ZZ->Add(hs_WJets);
  hs_ZZ->Add(hs_WZ);
  hs_ZZ->SetLineColor(kRed);
  hs_ZZ->SetLineWidth(3);
  hs_ZZ->SetFillColor(kWhite);	//kWhite is 100 percent transparent, so it will not block out filled bins drawn behind it
#endif

  hs_Other->SetFillColor(kWhite);
  hs_Other->SetLineColor(kBlue);
  hs_Other->SetLineStyle(7);
  hs_Other->SetLineWidth(3);

#ifdef showWR
  hs_ZZ->SetLineColor(kRed);
  hs_ZZ->SetLineWidth(3);
  //hs_ZZ->Scale(0.1);
  hs_ZZ->SetFillColor(kWhite);	//kWhite is 100 percent transparent, so it will not block out filled bins drawn behind it
#endif

  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0);
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0);
#ifdef plotRatio 
  p1->Draw();
  p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
#endif
  
  
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Preliminary         2.6 fb^{-1} (13 TeV)");
  hs_data->SetTitle("CMS Preliminary         2.6 fb^{-1} (13 TeV)");
  //th->SetTitle("CMS Private         2.6 fb^{-1} (13 TeV)");
  //hs_data->SetTitle("CMS Private         2.6 fb^{-1} (13 TeV)");
  th->Draw("histo");
#ifdef unblindedData
  hs_data->Draw("EPsame");
#endif

#ifndef unblindedData
  //draw a low WR mass signal as a line
  hs_data->SetLineColor(kRed);
  hs_data->SetLineWidth(3);
  hs_data->SetFillColor(kWhite);
  hs_data->Scale(0.3);
  //hs_data->Draw("HIST same");
#endif

#ifdef showRescaledRunOneEEJJExcess
  if(channel == Selector::EE) hs_ZZ->Draw("HIST same");
#endif

  hs_Other->Draw("HIST same");	//draw the distribution found in high mass WR signal events

#ifdef showWR
  hs_ZZ->Draw("HIST same");
#endif

  TString ytitle = "Events";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  hs_data->GetYaxis()->SetTitle(ytitle.Data());
  if(fname.EqualTo("Mlljj") || fname.EqualTo("Mlljj_varBins") || fname.EqualTo("Mlljj_2012Bins") ) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  if(fname.EqualTo("Mlljj_lowZpt")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV] for Z P_{T}<150"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV] for Z P_{T}<150");
  if(fname.EqualTo("Mlljj_varBins") || fname.EqualTo("Mlljj_2012Bins") ) hs_data->GetYaxis()->SetTitle("Events/GeV"), th->GetYaxis()->SetTitle("Events/GeV");
 
  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj") || fname.EqualTo("Mlljj_varBins") || fname.EqualTo("Mlljj_2012Bins") ) ratio->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(0.2);
  ratio->SetLabelSize(0.18,"x");
  ratio->GetYaxis()->SetTitle("data/MC   ");
  ratio->GetYaxis()->SetTitleSize(0.18);
  ratio->GetYaxis()->SetTitleOffset(0.2);

  leg->Draw(); 
#ifndef plotRatio  
  mycanvas->cd();	//no ratio plot
#endif
#ifdef plotRatio
  p2->cd();	//for ratio plot
#endif
  ratio->Sumw2();
  ratio->SetStats(0);

  /* old check*/
  if(fname.EqualTo("Mlljj_2012Bins")){
	  //print number of evts passing all cuts
	  std::cout<<"in MLLJJ distribution bin 7 there are"<<std::endl;
	  std::cout<< hs_DY->GetBinContent(7)*hs_DY->GetBinWidth(7)<<"\t DY background events"<<std::endl;
	  std::cout<< hs_ttbar->GetBinContent(7)*hs_ttbar->GetBinWidth(7) <<"\t top background events"<<std::endl;
	  std::cout<<"data events\t"<< hs_data->GetBinContent(7)*hs_data->GetBinWidth(7) << std::endl;
	  std::cout<<"in e chnl 1000 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(3,11) <<std::endl;	//670 to 1210
	  std::cout<<"in e chnl 2200 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(14,29) <<std::endl;	//1450 to 2450
	  std::cout<<"in e chnl 3600 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(27,47) <<std::endl;	//2300 to 3870
	  //std::cout<< hs_WJets->Integral() <<"\tWJets weighted evts"<<std::endl;
	  //std::cout<< hs_WZ->Integral() <<"\tWZ weighted evts"<<std::endl;
	  //std::cout<< hs_ZZ->Integral() <<"\tZZ weighted evts"<<std::endl;
  }/**/

  hs_ttbar->Add(hs_WJets);
  hs_ttbar->Add(hs_WZ);
  hs_ttbar->Add(hs_DY);

  ratio->Divide(hs_ttbar);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(0.18,"y");
  ratio->GetYaxis()->SetRangeUser(0.,2.4);
  if(channel == Selector::EE) ratio->GetYaxis()->SetRangeUser(0.,1.9);
  if(fname.EqualTo("Mljj_leadLept") || fname.EqualTo("Mljj_subleadLept") ) ratio->GetYaxis()->SetRangeUser(0.,1.6);
  ratio->GetYaxis()->SetNdivisions(505);
  
#ifdef plotRatio 
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  mycanvas->Update();
#endif

  TString fn = fname + "_MWR3000Signal";
  //TString fn = fname;


#ifndef unblindedData
  //if(channel == Selector::EE) fn += "_SignalRegion_EEChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_MWR1000Signal";
  //if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_MWR1000Signal";
  if(channel == Selector::EE) fn += "_SignalRegion_EEChannelBkgndMC_DYMadHTAndIncl_TTBarFromData";
  if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelBkgndMC_DYMadHTAndIncl_TTBarFromData";
#endif

#ifdef unblindedData
  if(channel == Selector::EE) fn += "_SignalRegion_EEChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_WithUnblindedData";
  if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_WithUnblindedData";
#endif

#ifdef plotRatio
  fn += "_withRatio";
#endif


#ifdef showWR
  fn += "_MWR2200MNu1100Signal";
  //fn += "_MWR800MNu400Signal";
#endif

#ifdef showRescaledRunOneEEJJExcess
  if(channel == Selector::EE) fn += "_withRescaledRunOneEEJJExcess";
#endif
  
  //plots with fixed bin widths
  //if(fname.EqualTo("Mlljj") || fname.EqualTo("Mljj_leadLept") || fname.EqualTo("Mljj_subleadLept") ){
  if(fname.EqualTo("Mlljj")){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  mycanvas->Print((fn+".C").Data());
	  th->SetMinimum(0.1);
	  mycanvas->Update();
#ifdef plotRatio
	  p1->SetLogy();	//for ratio plot
#endif
#ifndef plotRatio
	  mycanvas->SetLogy();	//only needed when no ratio plot is drawn
#endif
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
	  mycanvas->Print((fn+"_log.C").Data());
  }


  //plots with variable bin widths (the bin contents are divided by the bin widths)
  if(fname.EqualTo("Mlljj_2012Bins") ){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  mycanvas->Print((fn+".C").Data());
	  if(fname.EqualTo("Mlljj_2012Bins")) th->SetMinimum(0.0005);
	  mycanvas->Update();
#ifdef plotRatio
	  p1->SetLogy();	//for ratio plot
#endif
#ifndef plotRatio
	  mycanvas->SetLogy();	//only needed when no ratio plot is drawn
#endif
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
	  mycanvas->Print((fn+"_log.C").Data());
  }

  mycanvas->Close();
}
