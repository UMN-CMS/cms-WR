#ifndef selector_h
#define selector_h

/* #include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h" */
/* #include "ExoAnalysis/cmsWR/interface/Objects.h" */
#include "../interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Objects.h"


class Selector
{
public:

	myElectronCollection electrons;
	myMuonCollection muons;
	myJetCollection jets;

	Float_t WR_mass; // this is of Float_t because want to save it into a tree
	Float_t dilepton_mass;

	Float_t lead_lepton_pt; // flatten the collections for easy plotting
	Float_t sublead_lepton_pt;
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
