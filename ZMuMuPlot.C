/**
 * this is the main program for the analysis
 * please document it here
 */

#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
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

  TChain *chain_DY = myReader.getMiniTreeChain(TTchainNames, "miniTree_dytagandprobe");
  TChain *chain_DATA = myReader.getMiniTreeChain(TTchainNamesD, "miniTree_dytagandprobe");

  // TChain * chain_DY = new TChain("miniTree_dytagandprobe/t");
  // chain_DY->Add("~/scratch/myfile_dy.root");
  // TChain * chain_DATA = new TChain("miniTree_dytagandprobe/t");
  // chain_DATA->Add("~/scratch/myfile_data.root");

  std::cout << chain_DY->GetEntries() << std::endl;
  std::cout << chain_DATA->GetEntries() << std::endl;

  TH1F *h_Mll_DY = new TH1F("h_Mll_DY","",150,50,200);
  TH1F *h_Mll_DATA = new TH1F("h_Mll_DATA","",150,50,200);
  
  miniTreeEvent myEventDY;  
  myEventDY.SetBranchAddresses(chain_DY);
  miniTreeEvent myEventDATA;  
  myEventDATA.SetBranchAddresses(chain_DATA);

  int dy_nev = chain_DY->GetEntries();
  for(int i = 0; i < dy_nev; i++){
    chain_DY->GetEntry(i);
    if(myEventDY.muons_p4->size() < 2)
      continue;
    //if(myEventDY.muons_p4->at(0).Pt() > 60 && myEventDY.muons_p4->at(1).Pt()>50)
    h_Mll_DY->Fill((myEventDY.muons_p4->at(0)+myEventDY.muons_p4->at(1)).M());
  }
  int data_nev = chain_DATA->GetEntries();
  for(int i = 0; i < data_nev; i++){
    chain_DATA->GetEntry(i);
    if(myEventDATA.muons_p4->size() < 2)
      continue;
    //if(myEventDATA.muons_p4->at(0).Pt() > 60 && myEventDATA.muons_p4->at(1).Pt()>50)
    h_Mll_DATA->Fill((myEventDATA.muons_p4->at(0)+myEventDATA.muons_p4->at(1)).M());
  }


  h_Mll_DY->SetLineColor(kRed);
  h_Mll_DY->Scale(1975*2.52e3/2836876);
  h_Mll_DY->Draw();

  h_Mll_DATA->Draw("same");

  TLegend *leg = new TLegend( 0.78, 0.50, 0.98, 0.65 ) ;
  leg->AddEntry( h_Mll_DATA, "DATA" ) ;
  leg->AddEntry( h_Mll_DY, "DY" ) ;
  //Float_t sizeleg = leg->GetTextSize()*1.2;
  //leg->SetTextSize( sizeleg ) ;
  leg->SetFillColor( kWhite ) ;
  leg->Draw();

}
