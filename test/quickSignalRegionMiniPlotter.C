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

//#define doOnlyDYMadInclAndHT
#define useDYMAD

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
Selector::tag_t channel = Selector::EE;

void MakeHistos(TChain* chain, Selector *myEvent, std::vector<TH1F*> *hs);
void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname);
void quickSignalRegionMiniPlotter(){

//#ifndef doOnlyDYMadInclAndHT
//  TChain * chain_DY = new TChain("Tree_Iter0","DYMC");
//  //TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
//  TChain * chain_ttbar = new TChain("Tree_Iter0","TTData");
//#endif
//#ifdef doOnlyDYMadInclAndHT
//  TChain * chain_DY = new TChain("Tree_Iter0","DYMAD");
//  TChain * chain_ttbar = new TChain("Tree_Iter0","NOTT");
//#endif

  TChain * chain_DY = new TChain("Tree_Iter0","DYMC");
  //TChain * chain_ttbar = new TChain("Tree_Iter0","TTMC");
  TChain * chain_ttbar = new TChain("Tree_Iter0","TTData");
  TChain * chain_WJets = new TChain("Tree_Iter0","WJets");
  TChain * chain_WZ = new TChain("Tree_Iter0","WZ");
  TChain * chain_ZZ = new TChain("Tree_Iter0","ZZ");
  TChain * chain_data = new TChain("Tree_Iter0","Data");

  TString localDir = "../analysisCppOutputRootFiles/";
  Int_t data=0, dy=0, tt=0, wjets=0, wz=0, zz=0;
  switch (channel) {
  case Selector::EE:
//#ifndef doOnlyDYMadInclAndHT
//	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");
//	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_eeEE.root");
//    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
//	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_eeEE.root");
//    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
//    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
//#endif
//
//#ifdef doOnlyDYMadInclAndHT
//	wjets = chain_DY->Add(localDir+"selected_tree_DYMAD_signal_eeEEInclusiveLowHT.root");
//	dy = chain_ttbar->Add(localDir+"selected_tree_DYMADHT_signal_eeEEOnlyHT100to200.root");
//	tt = chain_WJets->Add(localDir+"selected_tree_DYMADHT_signal_eeEEOnlyHT200to400.root");
//    wz = chain_WZ->Add(localDir+"selected_tree_DYMADHT_signal_eeEEOnlyHT400to600.root");
//    zz = chain_ZZ->Add(localDir+"selected_tree_DYMADHT_signal_eeEEOnlyHT600toInf.root");
//#endif

#ifndef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_eeEE_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_eeEE.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
#endif

#ifdef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYMADInclAndHT_signal_eeEE_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_eeEE.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_eeEE.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_eeEE.root");
    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_eeEE.root");
#endif
	//data = chain_data->Add(localDir+"selected_tree_WRtoEEJJ_1400_700_signal_eeEE.root");
	data = chain_data->Add(localDir+"selected_tree_WRtoEEJJ_2200_1100_signal_eeEE.root");
    //data = chain_data->Add(localDir+"selected_tree_WRtoEEJJ_3200_1600_signal_eeEE.root");
	//data = chain_data->Add(localDir+"selected_tree_DYAMC_signal_eeEE.root");
	//data = chain_data->Add(localDir+"selected_tree_DYPOWHEG_signal_eeEE.root");
	break;
  case Selector::MuMu:
//#ifndef doOnlyDYMadInclAndHT
//	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");
//	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_mumuMuMu.root");
//    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuMuMu.root");
//	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_mumuMuMu.root");
//    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
//    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");
//#endif

//#ifdef doOnlyDYMadInclAndHT
//	wjets = chain_DY->Add(localDir+"selected_tree_DYMAD_signal_mumuMuMuInclusiveLowHT.root");
//	dy = chain_ttbar->Add(localDir+"selected_tree_DYMADHT_signal_mumuMuMuOnlyHT100to200.root");
//	tt = chain_WJets->Add(localDir+"selected_tree_DYMADHT_signal_mumuMuMuOnlyHT200to400.root");
//    wz = chain_WZ->Add(localDir+"selected_tree_DYMADHT_signal_mumuMuMuOnlyHT400to600.root");
//    zz = chain_ZZ->Add(localDir+"selected_tree_DYMADHT_signal_mumuMuMuOnlyHT600toInf.root");
//#endif

#ifndef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYAMC_signal_mumuMuMu_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_mumuMuMu.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuMuMu.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");
#endif

#ifdef useDYMAD
	dy = chain_DY->Add(localDir+"selected_tree_DYMADInclAndHT_signal_mumuMuMu_withMllWeight.root");
	//tt = chain_ttbar->Add(localDir+"selected_tree_TT_signal_mumuMuMu.root");
    tt = chain_ttbar->Add(localDir+"selected_tree_data_flavoursidebandEMuEE.root");
	wjets = chain_WJets->Add(localDir+"selected_tree_W_signal_mumuMuMu.root");
    wz = chain_WZ->Add(localDir+"selected_tree_WZ_signal_mumuMuMu.root");
    zz = chain_ZZ->Add(localDir+"selected_tree_ZZ_signal_mumuMuMu.root");
#endif	
	//data = chain_data->Add(localDir+"selected_tree_DYPOWHEG_signal_mumuMuMu.root");
   	//data = chain_data->Add(localDir+"selected_tree_DYAMC_signal_mumuMuMu.root");
    data = chain_data->Add(localDir+"selected_tree_WRtoMuMuJJ_2200_1100_signal_mumuMuMu.root");
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

  TString xtitles[] = {"leading lepton p_{T}","subleading lepton p_{T}","leading jet p_{T}","subleading jet p_{T}","leading lepton #eta","subleading lepton #eta","leading jet #eta","subleading jet #eta","leading lepton #phi","subleading lepton #phi","leading jet #phi","subleading jet #phi","Mlljj","dilepton mass","nPV","Z p_{T}","M_{LLJJ} (GeV)"};

  TString fnames[] = {"l1_pt","l2_pt","j1_pt","j2_pt","l1_eta","l2_eta","j1_eta","j2_eta","l1_phi","l2_phi","j1_phi","j2_phi","Mlljj","Mll","nPV","Z_pt","Mlljj_lowZpt"};

  int i = 0;
  for(unsigned int i = 0; i < nPlots; i++){
    std::string s = std::to_string(i);
    drawPlots(hs_DY[i],hs_ttbar[i],hs_WJets[i],hs_WZ[i],hs_ZZ[i],hs_data[i],xtitles[i],fnames[i]);
  }
  
}//end quickSignalRegionMiniPlotter()

void MakeHistos(TChain * chain, Selector *myEvent, std::vector<TH1F*> *hs){

  TH1F *h_lepton_pt0 = new TH1F("h_lepton_pt0","",50,0,700);
  TH1F *h_lepton_eta0 = new TH1F("h_lepton_eta0","",50,-3,3);
  TH1F *h_lepton_phi0 = new TH1F("h_lepton_phi0","",50,-3.15,3.15);
  TH1F *h_lepton_pt1 = new TH1F("h_lepton_pt1","",50,0,700);
  TH1F *h_lepton_eta1 = new TH1F("h_lepton_eta1","",50,-3,3);
  TH1F *h_lepton_phi1 = new TH1F("h_lepton_phi1","",50,-3.15,3.15);

  TH1F *h_jet_pt0 = new TH1F("h_jet_pt0","",50,0,700);
  TH1F *h_jet_eta0 = new TH1F("h_jet_eta0","",50,-3,3);
  TH1F *h_jet_phi0 = new TH1F("h_jet_phi0","",50,-3.15,3.15);
  TH1F *h_jet_pt1 = new TH1F("h_jet_pt1","",50,0,700);
  TH1F *h_jet_eta1 = new TH1F("h_jet_eta1","",50,-3,3);
  TH1F *h_jet_phi1 = new TH1F("h_jet_phi1","",50,-3.15,3.15);
  TH1F *h_WR_mass = new TH1F("h_WR_mass","",15,600,3000);
 
  TH1F *h_dilepton_mass = new TH1F("h_dilepton_mass","",50,150,1600);
  TH1F *h_nPV = new TH1F("h_nPV","",100,0,100);

  TH1F *h_Z_pt = new TH1F("h_Z_pt", "", 40, 0., 1200.);
  TH1F *h_WR_mass_lowZpt = new TH1F("h_WR_mass_lowZpt","",40,500,3000);
 

  Long64_t nEntries = chain->GetEntries();

  cout<< nEntries << endl;

  TString chainTitle(chain->GetTitle());
  Float_t scaleFactor = 1.0;
  if( chainTitle.EqualTo("TTMC") ){
	  scaleFactor = (channel == Selector::MuMu) ? 0.958 : 0.954;	//to account for slightly higher number of rescaled ttbar MC events relative to rescaled emu data evts
  }
  if( chainTitle.EqualTo("TTData") ){
	  scaleFactor = (channel == Selector::MuMu) ? 0.657 : 0.414;	//to rescale emu data evts to estimates of ttbar in ele and mu channels
  }

  for(int ev = 0; ev<nEntries; ++ev){
    chain->GetEntry(ev);
	if(myEvent->WR_mass < 600.) continue;
	if(myEvent->dilepton_mass < 200.) continue;

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
  hs->push_back(h_Z_pt);
  hs->push_back(h_WR_mass_lowZpt);

}

void drawPlots(TH1F* hs_DY,TH1F* hs_ttbar,TH1F* hs_WJets,TH1F* hs_WZ,TH1F* hs_ZZ,TH1F* hs_data, TString xtitle, TString fname){

  TLegend *leg = new TLegend( 0.72, 0.50, 0.98, 0.80 ) ; 
//#ifndef doOnlyDYMadInclAndHT
//  leg->AddEntry( hs_DY, "DY AMCNLO" ) ; 
//  leg->AddEntry( hs_ttbar, "TT Data Driven" ) ;
//  leg->AddEntry( hs_WJets, "WJets" ) ; 
//  leg->AddEntry( hs_WZ, "WZ" ) ; 
//  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
//#endif 
//
//#ifdef doOnlyDYMadInclAndHT
//  leg->AddEntry( hs_DY, "DYMadHT100to200" ) ; 
//  leg->AddEntry( hs_ttbar, "DYMadHT200to400" ) ;
//  leg->AddEntry( hs_WJets, "DYMadIncl HT<100" ) ; 
//  leg->AddEntry( hs_WZ, "DYMadHT400to600" ) ; 
//  leg->AddEntry( hs_ZZ, "DYMadHT600toInf" ) ; 
//  leg->AddEntry( hs_data, "DYAMC No MLL Weight");
//#endif

#ifndef useDYMAD
  leg->AddEntry( hs_DY, "DY AMCNLO" ) ; 
  leg->AddEntry( hs_ttbar, "TT Data Driven" ) ;
  leg->AddEntry( hs_WJets, "WJets" ) ; 
  leg->AddEntry( hs_WZ, "WZ" ) ; 
  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
#endif 

#ifdef useDYMAD
  leg->AddEntry( hs_DY, "DYMadHT+Incl" ) ; 
  leg->AddEntry( hs_ttbar, "TT Data Driven" ) ;
  leg->AddEntry( hs_WJets, "WJets" ) ; 
  leg->AddEntry( hs_WZ, "WZ" ) ; 
  leg->AddEntry( hs_ZZ, "ZZ" ) ; 
#endif 


  //leg->AddEntry( hs_data, "W_{R} Signal","l");
  //leg->AddEntry( (TObject*)0, "M_{WR}=2.2 TeV M_{NuR}=1.1 TeV","");
  leg->SetFillColor( kWhite ) ; 

  hs_data->Sumw2();
  hs_ttbar->Sumw2();
  hs_WJets->Sumw2();
  hs_WZ->Sumw2();
  hs_ZZ->Sumw2();
  hs_DY->Sumw2();
  
  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  mycanvas->cd();	//only needed when no ratio plot is drawn
  THStack* th = new THStack();
  hs_DY->SetFillColor(kYellow);
  hs_ttbar->SetFillColor(kGreen);
  hs_WJets->SetFillColor(kBlue);
  hs_WZ->SetFillColor(kCyan);
  hs_ZZ->SetFillColor(kMagenta);
  th->Add(hs_WZ);
  th->Add(hs_ZZ);
  th->Add(hs_WJets);
  th->Add(hs_DY);
  th->Add(hs_ttbar);
  hs_data->SetLineColor(kRed);
  hs_data->SetLineWidth(2);


  /*for ratio plot
  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  */
  hs_data->SetStats(0);
  TH1F *ratio = (TH1F*)hs_data->Clone();
  th->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  hs_data->SetTitle("CMS Private #surds = 13 TeV #int lumi = 2.6 fb^{-1}");
  th->Draw("histo");
//#ifdef doOnlyDYMadInclAndHT
//  hs_data->Draw("histo");	///<draw this first if hs_data has a greater max than th
//#endif 
  //th->Draw("histo same");
  //hs_data->Draw("histo same");
  //TString ytitle = "Events/(";
  //ytitle += (th->GetXaxis()->GetNbins());
  //ytitle += ")";
  TString ytitle = "Events";
  th->GetYaxis()->SetTitle(ytitle.Data());
  th->GetXaxis()->SetTitle(xtitle.Data());
  hs_data->GetYaxis()->SetTitle(ytitle.Data());
  if(fname.EqualTo("Mlljj")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV]"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  if(fname.EqualTo("Mlljj_lowZpt")) hs_data->GetXaxis()->SetTitle("M_{LLJJ} [GeV] for Z P_{T}<150"), th->GetXaxis()->SetTitle("M_{LLJJ} [GeV] for Z P_{T}<150");

  ratio->GetXaxis()->SetTitle(xtitle.Data());
  if(fname.EqualTo("Mlljj")) ratio->GetXaxis()->SetTitle("M_{LLJJ} [GeV]");
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(0.2);
  ratio->SetLabelSize(0.1,"x");
  leg->Draw(); 
  mycanvas->cd();
  //p2->cd();	//for ratio plot
  ratio->Sumw2();
  ratio->SetStats(0);

  if(fname.EqualTo("Mlljj")){
	  //print number of evts passing all cuts
	  std::cout<<"in MLLJJ distribution there are"<<std::endl;
	  std::cout<< hs_DY->Integral() <<"\tDY weighted evts"<<std::endl;
	  std::cout<< hs_ttbar->Integral() <<"\tttbar weighted evts"<<std::endl;
	  std::cout<<"in e chnl 1000 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(3,11) <<std::endl;	//670 to 1210
	  std::cout<<"in e chnl 2200 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(14,29) <<std::endl;	//1450 to 2450
	  std::cout<<"in e chnl 3600 GeV WR window TTBar has this many evts\t"<< hs_ttbar->Integral(27,47) <<std::endl;	//2300 to 3870
	  std::cout<< hs_WJets->Integral() <<"\tWJets weighted evts"<<std::endl;
	  std::cout<< hs_WZ->Integral() <<"\tWZ weighted evts"<<std::endl;
	  std::cout<< hs_ZZ->Integral() <<"\tZZ weighted evts"<<std::endl;
  }

  hs_ttbar->Add(hs_WJets);
  hs_ttbar->Add(hs_WZ);
  hs_ttbar->Add(hs_ZZ);
  hs_ttbar->Add(hs_DY);

  ratio->Divide(hs_ttbar);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(0.1,"y");
  ratio->GetYaxis()->SetRangeUser(0.5,2.0);
  ratio->GetYaxis()->SetNdivisions(505);
  /*for ratio plot
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  mycanvas->Update();
  */

  TString fn = fname;

/*#ifndef doOnlyDYMadInclAndHT
  if(channel == Selector::EE) fn += "_SignalRegion_EEChannelWR2p2TeVandBkgndMC_TTBarFromData_FixedBinWidthNoRatio";
  if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelWR2p2TeVandBkgndMC_TTBarFromData_FixedBinWidthNoRatio";
#endif

#ifdef doOnlyDYMadInclAndHT
  if(channel == Selector::EE) fn += "_SignalRegion_DYMadInclAndDYMadHTBinnedAndDYAMCNLO_EEChannel_FixedBinWidthNoRatio";
  if(channel == Selector::MuMu) fn += "_SignalRegion_DYMadInclAndDYMadHTBinnedAndDYAMCNLO_MuMuChannel_FixedBinWidthNoRatio";
#endif*/

#ifndef useDYMAD
  if(channel == Selector::EE) fn += "_SignalRegion_EEChannelBkgndMC_DYAMC_TTBarFromData_FixedBinWidthNoRatio";
  if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelBkgndMC_DYAMC_TTBarFromData_FixedBinWidthNoRatio";
#endif

#ifdef useDYMAD
  if(channel == Selector::EE) fn += "_SignalRegion_EEChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_FixedBinWidthNoRatio";
  if(channel == Selector::MuMu) fn += "_SignalRegion_MuMuChannelBkgndMC_DYMadHTAndIncl_TTBarFromData_FixedBinWidthNoRatio";
#endif


  if(fname.EqualTo("Mlljj") || fname.EqualTo("Z_pt") || fname.EqualTo("Mll")){
	  mycanvas->Print((fn+".pdf").Data());
	  mycanvas->Print((fn+".png").Data());
	  mycanvas->Print((fn+".C").Data());
	  th->SetMinimum(0.1);
	  mycanvas->Update();
	  //p1->SetLogy();	//for ratio plot
	  mycanvas->SetLogy();	//only needed when no ratio plot is drawn
	  mycanvas->Print((fn+"_log.pdf").Data());
	  mycanvas->Print((fn+"_log.png").Data());
	  mycanvas->Print((fn+"_log.C").Data());
  }

  mycanvas->Close();
}
