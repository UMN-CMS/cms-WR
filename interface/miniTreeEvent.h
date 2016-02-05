#ifndef miniTreeEvent_h
#define miniTreeEvent_h
#include <TLorentzVector.h>
#include <TTree.h>
#include <TChain.h>

class miniTreeEvent
{
public:
	Int_t run;
	Int_t lumi;
	Long64_t event;
	char datasetName[30];

	// RunTime
	// UInt_t run_time;

	std::vector<TLorentzVector> * electrons_p4;
	std::vector<TLorentzVector> * muons_p4;
	std::vector<TLorentzVector> * jets_p4;

	std::vector<Float_t> * jec_uncertainty;
	std::vector<Float_t> * electron_scale;
	std::vector<Float_t> * electron_smearing;
	std::vector<Int_t> * electron_charge;
	std::vector<Int_t> * muon_charge;
	std::vector<Float_t> * muon_IDSF_central;
	std::vector<Float_t> * muon_IsoSF_central;
	std::vector<Float_t> * muon_IDSF_error;
	std::vector<Float_t> * muon_IsoSF_error;

	Float_t nPU;
	Int_t nPV;
	Float_t weight;
	Float_t PU_reweight;

	miniTreeEvent(); ///< default contructor (empty)
	miniTreeEvent(const miniTreeEvent& otherEvent);

	void clear();
	void SetBranches(TTree* tree);
	void SetBranchAddresses(TChain* tree);


private:
	bool _owningMembers;
};


#endif
