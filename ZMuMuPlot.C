/**
 * this is the main program for the analysis
 * please document it here
 */

#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"

#include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"

#include <vector>
#include <string>
#include <fstream>
#include "ExoAnalysis/cmsWR/src/configReader.cc"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

void ZMuMuPlot(){

  char name[100];// to name the tree per iteration of the Toy
  float integratedLumi = 2.52e3; ///\todo should not be hard-coded!!!

  std::vector<std::string> TTchainNames; 
  TTchainNames.push_back("DYToMuMu_powheg_50to120");
  // TTchainNames.push_back("DYToMuMu_powheg_120to200");
  // TTchainNames.push_back("DYToMuMu_powheg_200to400");
  // TTchainNames.push_back("DYToMuMu_powheg_400to800");
  // TTchainNames.push_back("DYToMuMu_powheg_800to1400");
  // TTchainNames.push_back("DYToMuMu_powheg_1400to2300");
  // TTchainNames.push_back("DYToMuMu_powheg_2300to3500");
  // TTchainNames.push_back("DYToMuMu_powheg_3500to4500");
  // TTchainNames.push_back("DYToMuMu_powheg_4500to6000");
  // TTchainNames.push_back("DYToMuMu_powheg_6000toInf");

  std::vector<std::string> TTchainNamesD; 
  TTchainNamesD.push_back("SingleMu_RunC");
  TTchainNamesD.push_back("SingleMu_RunD_v3");
  TTchainNamesD.push_back("SingleMu_RunD_v4");

  // TTchainNamesD.push_back("DoubleEG_RunC");
  // TTchainNamesD.push_back("DoubleEG_RunD_v3");
  // TTchainNamesD.push_back("DoubleEG_RunD_v4");

  configReader myReader("configs/2015-v1.conf");

  //TChain *chain_DY = myReader.getMiniTreeChain(TTchainNames, "miniTree_dytagandprobe");
  //TChain *chain_DATA = myReader.getMiniTreeChain(TTchainNamesD, "miniTree_dytagandprobe");

  TChain * chain_DY = new TChain("miniTree_dytagandprobe/t");
  chain_DY->Add("~/scratch/dy50_120.root");
  TChain * chain_DY_120 = new TChain("miniTree_dytagandprobe/t");
  chain_DY_120->Add("~/scratch/myfile_dy_120.root");
  TChain * chain_DATA = new TChain("miniTree_dytagandprobe/t");
  chain_DATA->Add("~/scratch/singlemuC.root");
  chain_DATA->Add("~/scratch/singlemuD_v3.root");
  chain_DATA->Add("~/scratch/singlemuD_v4.root");

  std::cout << chain_DY->GetEntries() << std::endl;
  std::cout << chain_DY_120->GetEntries() << std::endl;
  std::cout << chain_DATA->GetEntries() << std::endl;

  TH1F *h_Mll_DY = new TH1F("h_Mll_DY","",40,70,110);
  TH1F *h_Mll_DY_120 = new TH1F("h_Mll_DY_120","",40,70,110);
  TH1F *h_Mll_DATA = new TH1F("h_Mll_DATA","",40,70,110);
  
  miniTreeEvent myEventDY;  
  myEventDY.SetBranchAddresses(chain_DY);
  miniTreeEvent myEventDY_120;  
  myEventDY_120.SetBranchAddresses(chain_DY_120);
  miniTreeEvent myEventDATA;  
  myEventDATA.SetBranchAddresses(chain_DATA);

  int njets_any,njets_eq_0,njets_gt_0,njets_gt_1,njets_gt_2,njets_gt_3,njets_gt_4;
  njets_any = njets_eq_0 = njets_gt_0 = njets_gt_1 = njets_gt_2 = njets_gt_3 = njets_gt_4 = 0;
  int njets_any_data,njets_eq_data_0,njets_gt_data_0,njets_gt_data_1,njets_gt_data_2,njets_gt_data_3,njets_gt_data_4;
  njets_any_data = njets_eq_data_0 = njets_gt_data_0 = njets_gt_data_1 = njets_gt_data_2 = njets_gt_data_3 = njets_gt_data_4 = 0;

  float jet_pt_cut = 0.0;
  float lepton1_pt_cut = 60.0;
  float lepton2_pt_cut = 50.0;
  float jet1_pt_cut = 190.0;
  float jet1_pt_cut_low = 99999190.0;
  float jet2_pt_cut = 40.0;
  float jet2_pt_cut_low = 99999110.0;
  
  int dy_nev = chain_DY->GetEntries();
  for(int i = 0; i < dy_nev; i++){
    chain_DY->GetEntry(i);
    if(myEventDY.muons_p4->size() < 2)// || myEventDY.jets_p4->size() != 0)
      continue;
    float weight = myEventDY.weight * myEventDY.PU_reweight * myEventDY.muon_IDSF_central->at(0) * myEventDY.muon_IsoSF_central->at(0) * myEventDY.muon_IDSF_central->at(1) * myEventDY.muon_IsoSF_central->at(1);
    if(myEventDY.muons_p4->at(0).Pt() > lepton1_pt_cut && myEventDY.muons_p4->at(1).Pt()>lepton2_pt_cut){ 
      if(myEventDY.jets_p4->size() > 2){
	if(myEventDY.jets_p4->at(0).Pt() > jet1_pt_cut && myEventDY.jets_p4->at(1).Pt()>jet2_pt_cut){
	  if(myEventDY.jets_p4->at(0).Pt() < jet1_pt_cut_low){
	    h_Mll_DY->Fill((myEventDY.muons_p4->at(0)+myEventDY.muons_p4->at(1)).M(),weight);
	  }
	}
      }
    }
    vector<TLorentzVector> hard_jets;
    for(auto j:*(myEventDY.jets_p4)){
      if(j.Pt()>jet_pt_cut )
	hard_jets.push_back(j);
    }
    if(myEventDY.muons_p4->at(0).Pt() > lepton1_pt_cut && myEventDY.muons_p4->at(1).Pt()>lepton2_pt_cut){
      if((myEventDY.muons_p4->at(0)+myEventDY.muons_p4->at(1)).M() > 80 && (myEventDY.muons_p4->at(0)+myEventDY.muons_p4->at(1)).M() < 100){
	njets_any++;
	if(hard_jets.size() == 0)
	  njets_eq_0++;
	if(hard_jets.size() > 0)
	  njets_gt_0++;
	if(hard_jets.size() > 1)
	  njets_gt_1++;
	if(hard_jets.size() > 2)
	  njets_gt_2++;
	if(hard_jets.size() > 3)
	  njets_gt_3++;
	if(hard_jets.size() > 4)
	  njets_gt_4++;
      }
    }
  }
  int data_nev = chain_DATA->GetEntries();
  for(int i = 0; i < data_nev; i++){
    chain_DATA->GetEntry(i);
    if(myEventDATA.muons_p4->size() < 2)// || myEventDATA.jets_p4->size() != 0)
      continue;
    if(myEventDATA.muons_p4->at(0).Pt() > lepton1_pt_cut && myEventDATA.muons_p4->at(1).Pt()>lepton2_pt_cut){
      if(myEventDATA.jets_p4->size() > 2){
	if(myEventDATA.jets_p4->at(0).Pt() > jet1_pt_cut && myEventDATA.jets_p4->at(1).Pt()>jet2_pt_cut){
	  if(myEventDATA.jets_p4->at(0).Pt() < jet1_pt_cut_low){
	    h_Mll_DATA->Fill((myEventDATA.muons_p4->at(0)+myEventDATA.muons_p4->at(1)).M());
	  }
	}
      }
    }
    vector<TLorentzVector> hard_jets;
    for(auto j:*(myEventDATA.jets_p4)){
      if(j.Pt()>jet_pt_cut )
	hard_jets.push_back(j);
    }      
    if(myEventDATA.muons_p4->at(0).Pt() > lepton1_pt_cut && myEventDATA.muons_p4->at(1).Pt()>lepton2_pt_cut){
      if((myEventDATA.muons_p4->at(0)+myEventDATA.muons_p4->at(1)).M() > 80 && (myEventDATA.muons_p4->at(0)+myEventDATA.muons_p4->at(1)).M() < 100){
	njets_any_data++;
	if(hard_jets.size() == 0)
	  njets_eq_data_0++;
	if(hard_jets.size() > 0)
	  njets_gt_data_0++;
	if(hard_jets.size() > 1)
	  njets_gt_data_1++;
	if(hard_jets.size() > 2)
	  njets_gt_data_2++;
	if(hard_jets.size() > 3)
	  njets_gt_data_3++;
	if(hard_jets.size() > 4)
	  njets_gt_data_4++;
      }
    }
  }

  int dy_120_nev = chain_DY_120->GetEntries();
  for(int i = 0; i < dy_120_nev; i++){
    chain_DY_120->GetEntry(i);
    if(myEventDY_120.muons_p4->size() < 2)// || myEventDY_120.jets_p4->size() != 0)
      continue;
    if(myEventDY_120.muons_p4->at(0).Pt() > 60 && myEventDY_120.muons_p4->at(1).Pt()>50)
      h_Mll_DY_120->Fill((myEventDY_120.muons_p4->at(0)+myEventDY_120.muons_p4->at(1)).M());
  }

  h_Mll_DY->SetLineColor(kRed);
  h_Mll_DY->Scale(1987*2.52e3/2836876);
  h_Mll_DY_120->SetLineColor(kGreen);
  h_Mll_DY_120->Scale(19.3*2.52e3/100000);

  TCanvas* mycanvas = new TCanvas( "mycanvas", "", 0, 0, 600, 600 ) ;
  Double_t eps = 0.001;
  TPad* p1 = new TPad("p1","p1",0,0.25,1,1,0); p1->Draw();
  TPad* p2 = new TPad("p2","p2",0,0.1,1,0.25+eps,0); p2->Draw();
  p1->SetBottomMargin(0);
  p2->SetTopMargin(0);   
  p1->cd();
  TH1F *ratio = (TH1F*)h_Mll_DATA->Clone();
  h_Mll_DATA->Draw("EP");
  h_Mll_DY->Draw("same histo");
  //h_Mll_DY_120->Draw("same");

  TLegend *leg = new TLegend( 0.78, 0.50, 0.98, 0.65 ) ;
  leg->AddEntry( h_Mll_DATA, "DATA" ) ;
  leg->AddEntry( h_Mll_DY, "DY" ) ;
  //leg->AddEntry( h_Mll_DY_120, "120 < DY < 200 " ) ;
  //Float_t sizeleg = leg->GetTextSize()*1.2;
  //leg->SetTextSize( sizeleg ) ;
  leg->SetFillColor( kWhite ) ;
  ratio->GetXaxis()->SetTickSize(0.40);
  ratio->GetXaxis()->SetTitleSize(0.2);
  ratio->SetLabelSize(0.1,"x");
  leg->Draw();
  mycanvas->cd();
  p2->cd();
  ratio->Sumw2();
  ratio->SetStats(0);
  ratio->Divide(h_Mll_DY);
  ratio->SetMarkerStyle(21);
  ratio->SetLabelSize(0.1,"y");
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->Draw("p");
  float xmax = ratio->GetXaxis()->GetXmax();
  float xmin = ratio->GetXaxis()->GetXmin();
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
  ratio->Draw("p");
  f1->Draw("same");
  mycanvas->cd();
  // PDF name: l1_pt,l2_pt,j1_pt,j2_pt
  mycanvas->Print("mll_dataMC_60_50_190_40.pdf");

  cout<<"Int="<<h_Mll_DATA->Integral()<<endl;

  cout<<"Counts"<<endl;
  cout<<"NJets = Any: "<<njets_any_data/(njets_any*1987*2.52e3/2836876)<<endl;
  cout<<"NJets = 0: "<<njets_eq_data_0/(njets_eq_0*1987*2.52e3/2836876)<<endl;
  cout<<"NJets > 0: "<<njets_gt_data_0/(njets_gt_0*1987*2.52e3/2836876)<<endl;
  cout<<"NJets > 1: "<<njets_gt_data_1/(njets_gt_1*1987*2.52e3/2836876)<<endl;
  cout<<"NJets > 2: "<<njets_gt_data_2/(njets_gt_2*1987*2.52e3/2836876)<<endl;
  cout<<"NJets > 3: "<<njets_gt_data_3/(njets_gt_3*1987*2.52e3/2836876)<<endl;
  cout<<"NJets > 4: "<<njets_gt_data_4/(njets_gt_4*1987*2.52e3/2836876)<<endl;

  TFile f("histos.root","RECREATE");
  h_Mll_DATA->Write();
  h_Mll_DY->Write();
  h_Mll_DY_120->Write();


}
