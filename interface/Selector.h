#ifndef selector_h
#define selector_h

#include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Objects.h"

typedef std::vector<myElectron> myElectronCollection;
typedef std::vector<myMuon> myMuonCollection;
typedef std::vector<myJet> myJetCollection;

class Selector{
public:

  myElectronCollection electrons;
  myMuonCollection muons;
  myJetCollection jets;

  Float_t WR_mass; // this is of Float_t because want to save it into a tree
  Float_t dilepton_mass;

  Float_t weight;
  
  bool isPassing(Int_t tag);
  
  Selector(const miniTreeEvent& myEvent);
  Selector();
  
  void SetBranches(TTree* tree);
  void SetBranchAddresses(TTree* tree);


private:
	
    void Clear();

};



#endif
