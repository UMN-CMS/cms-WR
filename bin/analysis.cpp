/** 
 * this is the main program for the analysis
 * please document it here
 */

#include "TChain.h"
#include "TFile.h"
// #include "TLorentzVector.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
// #include <iostream>
// // #include "ExoAnalysis/cmsWR/interface/FitRooDataSet.h"

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Selector.h"
// #include "../interface/FitRooDataSet.h"
// #include "../src/Selector.cc"
// #include "../src/miniTreeEvent.cc"

#include <vector>
#include <string>
#include <fstream>
// #include "../interface/miniTreeEvent.h"
#include "FitRooDataSet.h"
#include "rooFitFxns.h"
#include "ToyThrower.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

int main(void){

  using namespace RooFit;

  TString mode = "ttbar";
  
  //TChain * chain = new TChain("MiniTTree/t");
  TChain * chain = new TChain("miniTree_signal/t");

  if(mode.EqualTo("ttbar")){
    //chain->Add("ttree.root");   
    chain->Add("~/eos/cms/store/user/shervin/ntuples/TTJets_DiLept_v1_SHv2/ttree_*.root");
  }

  // Plotting trees
  TFile f("selected_tree_"+mode+".root","recreate");
  TTree * t1 = new TTree("t1","");

  miniTreeEvent myEvent;

  myEvent.SetBranchAddresses(chain);
  Selector selEvent;
  selEvent.SetBranches(t1);
  
  Int_t nToys = 1;
  Int_t nEntries = chain->GetEntries();

  int isData= 1;// Fill in with 1 or 0 based on information from the trees
  TRandom3 Rand;
  std::vector<std::string> List_Systematics;
  std::string word;
  std::ifstream Syst_File;
  Syst_File.open("Systematics_To_Be_Evaluated.txt");
  while(!Syst_File.eof()){
       Syst_File>>word;      
       List_Systematics.push_back(word);
      }

  for(int i = 0; i < nToys; i++){
    Rand.SetSeed(i+1);
    RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,6500);
    RooRealVar evtWeight("evtWeight", "evtWeight", -2,2);
    RooArgSet vars(massWR,evtWeight);
    RooDataSet * tempDataSet = new RooDataSet("temp","temp",vars);

    for(int ev = 0; ev<nEntries; ev++){
      chain->GetEntry(ev);
    
      Selector tmp_selEvent(myEvent);
      selEvent = tmp_selEvent;

#ifdef DEBUG
      cout<<"RUN="<<myEvent.run<<endl;
      cout<<"Mu"<<endl;
      for(auto m:*(myEvent.muons_p4))
        cout<<m.Pt()<<" "<<m.Eta()<<endl;
      cout<<"Jet"<<endl;
      for(auto m:*(myEvent.jets_p4))
        cout<<m.Pt()<<" "<<m.Eta()<<endl;
#endif
      ToyThrower(myEvent, Rand, i+1, List_Systematics, isData); 


      // Select events with one good WR candidate
      // Tags:
      // 0 -- EEJJ Channel
      // 1 -- MuMuJJ Channel
      // 2 -- EMuJJ Channel

      int tag = 1;

      if(selEvent.isPassing(tag) && selEvent.dilepton_mass > 200){
	//cout<<"PASS"<<endl;       
    	
	float weight = selEvent.weight * 1.0; // Product of many weights

	selEvent.lead_jet_pt = selEvent.jets[0].p4.Pt();
	selEvent.sublead_jet_pt = selEvent.jets[1].p4.Pt();
	selEvent.lead_jet_eta = selEvent.jets[0].p4.Eta();
	selEvent.sublead_jet_eta = selEvent.jets[1].p4.Eta();
	selEvent.lead_jet_phi = selEvent.jets[0].p4.Phi();
	selEvent.sublead_jet_phi = selEvent.jets[1].p4.Phi();

	if(tag == 0){
	  selEvent.lead_lepton_pt = selEvent.electrons[0].p4.Pt();
	  selEvent.sublead_lepton_pt = selEvent.electrons[1].p4.Pt();
	  selEvent.lead_lepton_eta = selEvent.electrons[0].p4.Eta();
	  selEvent.sublead_lepton_eta = selEvent.electrons[1].p4.Eta();
	  selEvent.lead_lepton_phi = selEvent.electrons[0].p4.Phi();
	  selEvent.sublead_lepton_phi = selEvent.electrons[1].p4.Phi();
	}
	if(tag == 1){
	  selEvent.lead_lepton_pt = selEvent.muons[0].p4.Pt();
	  selEvent.sublead_lepton_pt = selEvent.muons[1].p4.Pt();
	  selEvent.lead_lepton_eta = selEvent.muons[0].p4.Eta();
	  selEvent.sublead_lepton_eta = selEvent.muons[1].p4.Eta();
	  selEvent.lead_lepton_phi = selEvent.muons[0].p4.Phi();
	  selEvent.sublead_lepton_phi = selEvent.muons[1].p4.Phi();
	}
	if(tag == 2){
	  selEvent.lead_lepton_pt = selEvent.electrons[0].p4.Pt();
	  selEvent.sublead_lepton_pt = selEvent.muons[0].p4.Pt();
	  selEvent.lead_lepton_eta = selEvent.electrons[0].p4.Eta();
	  selEvent.sublead_lepton_eta = selEvent.muons[0].p4.Eta();
	  selEvent.lead_lepton_phi = selEvent.electrons[0].p4.Phi();
	  selEvent.sublead_lepton_phi = selEvent.muons[0].p4.Phi();
	}

	massWR.setVal(selEvent.WR_mass);
	evtWeight.setVal(weight);
	tempDataSet->add(vars);

	t1->Fill();

      }

      // cout<<"Mu"<<endl;
      // for(auto m:*(myEvent.muons_p4))
      //   cout<<m.Pt()<<endl;
      // cout<<"NuMu"<<endl;
      // for(auto m:(selEvent.muons))
      //   cout<<m.p4.Pt()<<endl;
      
    }

    t1->Write();

    tempDataSet->Print();

    RooFitResult * tempFitRslt = NULL;
    fitRooDataSet(tempFitRslt, tempDataSet, expPdfRooAbsPdf);

    std::cout<<"Res="<<std::endl;
    expPdfRooAbsPdf->Print();

  }

  return 0;

}
