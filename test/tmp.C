#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
#include <iostream>
#include "ExoAnalysis/cmsWR/interface/FitRooDataSet.h"
#include "ExoAnalysis/cmsWR/src/Selector.cc"
#include "ExoAnalysis/cmsWR/src/miniTreeEvent.cc"
#include "rooFitFxns.C"
//#include "ToyThrower.C"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

void tmp(){

  using namespace RooFit;
  
  TFile * hfile = new TFile("/afs/cern.ch/user/j/jchavesb/public/WR/ttree.root");
  TTree *tree = (TTree*)hfile->Get("MiniTTree/t");  
  
  // Plotting trees
  TFile f("selected_tree.root","recreate");
  TTree * t1 = new TTree("t1","");

  miniTreeEvent myEvent;
  
  myEvent.SetBranchAddresses(tree, myEvent);
  Selector selEvent;
  selEvent.SetBranches(t1);
  
  int nToys = 1;

  for(int i = 0; i < nToys; i++){

    RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,6500);
    RooRealVar evtWeight("evtWeight", "evtWeight", -2,2);
    RooArgSet vars(massWR,evtWeight);
    RooDataSet * tempDataSet = new RooDataSet("temp","temp",vars);

    for(int ev = 0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
    
      Selector tmp_selEvent(myEvent);

      selEvent = tmp_selEvent;

      //ToyThrower(myEvent, ev+1);

      // Select events with one good WR candidate
      // Tags:
      // 0 -- EEJJ Channel
      // 1 -- MuMuJJ Channel
      // 2 -- EMuJJ Channel

      int tag = 1;
      
      if(selEvent.isPassing(tag) && selEvent.dilepton_mass > 200){
	//cout<<"PASS"<<endl;       
	//cout<<selEvent.WR_mass<<endl;
	
	float weight = selEvent.weight * 1.0; // Product of many weights
	Float_t lead_lepton_pt = (tag == 0 || tag == 2) ? selEvent.electrons[0].p4.Pt() : selEvent.muons[0].p4.Pt();
	Float_t sublead_lepton_pt = (tag == 0 || tag == 1) ? selEvent.electrons[1].p4.Pt() : selEvent.muons[1].p4.Pt();
	Float_t lead_lepton_eta;
	Float_t sublead_lepton_eta;
	Float_t lead_lepton_phi;
	Float_t sublead_lepton_phi;
	Float_t lead_jet_pt;
	Float_t sublead_jet_pt;
	Float_t lead_jet_eta;
	Float_t sublead_jet_eta;
	Float_t lead_jet_phi;
	Float_t sublead_jet_phi;

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

    // RooFitResult * tempFitRslt = NULL;
    // fitRooDataSet(tempFitRslt, tempDataSet, expPdfRooAbsPdf);

    // std::cout<<"Res="<<std::endl;
    // expPdfRooAbsPdf->Print();

  }


}
