#ifndef miniTreeEvent_h
#define miniTreeEvent_h
#include <TLorentzVector.h>
#include <TTree.h>
#include <TChain.h>

/** \class miniTreeEvent miniTreeEvent.h miniTreeEvent.cc
 * This class defines the content of one event in the miniTree produced from microAOD
 * There are also methods to read and to write from a TTree the event
 */
class miniTreeEvent
{
public:

// public members to be filled by your program
	Int_t run;
	Int_t lumi;
	Long64_t event;
	char datasetName[30]; ///< the main dataset name

	// RunTime
	// UInt_t run_time;

	std::vector<TLorentzVector> * electrons_p4;
	std::vector<Float_t> * electron_scale;
	std::vector<Float_t> * electron_smearing;
	std::vector<Int_t> * electron_charge;

	std::vector<TLorentzVector> * muons_p4;
	std::vector<Int_t> * muon_charge;
	std::vector<Float_t> * muon_IDSF_central;
	std::vector<Float_t> * muon_IsoSF_central;
	std::vector<Float_t> * muon_IDSF_error;
	std::vector<Float_t> * muon_IsoSF_error;

	std::vector<TLorentzVector> * jets_p4;
	std::vector<Float_t> * jec_uncertainty;

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
