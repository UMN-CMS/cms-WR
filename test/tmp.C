#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "../interface/miniTreeEvent.h"
#include "../interface/FitRooDataSet.h"
#include "rooFitFxns.C"
#include "ToyThrower.C"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

void tmp(){

  using namespace RooFit;
 
  TFile * hfile = new TFile("/afs/cern.ch/user/j/jchavesb/public/WR/ttree.root");
  TTree *tree = (TTree*)hfile->Get("MiniTTree/t");  
  
  miniTreeEvent myEvent;

  myEvent.SetBranchAddresses(tree, myEvent);
  
  int nToys = 1;
  int isData= 1;// Fill in with 1 or 0 based on information from the trees
  TRandom3* Rand=new TRandom3();
  vector<string> List_Systematics;
  string word;
  ifstream Syst_File;
  Syst_File.open("Systematics_To_Be_Evaluated.txt");
  while(!Syst_File.eof()){
       Syst_File>>word;      
       List_Systematics.push_back(word);
      }

  for(int i = 0; i < nToys; i++){
    Rand->SetSeed(i+1);
    RooRealVar massWR("fourObjectMass", "fourObjectMass", 600,6500);
    RooRealVar evtWeight("evtWeight", "evtWeight", -2,2);
    RooArgSet vars(massWR,evtWeight);
    RooDataSet * tempDataSet = new RooDataSet("temp","temp",vars);

    for(int ev = 0; ev<tree->GetEntries(); ev++){
      tree->GetEntry(ev);
      
      ToyThrower(myEvent, Rand, i+1, List_Systematics, isData); 

      massWR.setVal(1.0);
      evtWeight.setVal(1.0);
      tempDataSet->add(vars);
    
      // for(auto m:*(myEvent.muons_p4))
      //   cout<<m.Pt()<<endl;
    
    }

    tempDataSet->Print();

    RooFitResult * tempFitRslt = NULL;
    fitRooDataSet(tempFitRslt, tempDataSet, expPdfRooAbsPdf);

    std::cout<<"Res="<<std::endl;
    expPdfRooAbsPdf->Print();
  }

  delete Rand;
}
