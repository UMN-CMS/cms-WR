#ifndef selector_h
#define selector_h

/* #include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h" */
/* #include "ExoAnalysis/cmsWR/interface/Objects.h" */
#include "../interface/miniTreeEvent.h"
#include "ExoAnalysis/cmsWR/interface/Objects.h"
#include <stdexcept>


class Selector
{
public:
	enum tag_t {
		EE,
		MuMu,
		EMu
	};

	static tag_t getTag(std::string ch)
	{
		tag_t channel;
		if(ch == "EE")
			channel = Selector::EE;
		else if(ch == "MuMu")
			channel = Selector::MuMu;
		else if(ch == "EMu")
			channel = Selector::EMu;
		else {
			std::cerr << "[ERROR] Channel " << ch << " not recognized" << std::endl;
			throw std::invalid_argument( "Unrecognized channel" );
		}
		return channel;
	}

	std::string datasetName;

	Float_t WR_mass; // this is of Float_t because want to save it into a tree
	Float_t dilepton_mass;

	Float_t zPt;
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
	Float_t lead_jet_jec_unc;
	Float_t sublead_jet_jec_unc;
	Float_t lead_lepton_weight;
	Float_t sublead_lepton_weight;
	Float_t lead_jet_weight;
	Float_t sublead_jet_weight;
	Float_t dR_leadlepton_leadjet;
	Float_t dR_leadlepton_subleadjet;
	Float_t dR_subleadlepton_leadjet;
	Float_t dR_subleadlepton_subleadjet;
	Float_t lead_lepton_r9;
	Float_t sublead_lepton_r9;
	Float_t lead_lepton_IDSF_error;
	Float_t lead_lepton_IsoSF_error;
	Float_t lead_lepton_RecoSF_error;
	Float_t lead_lepton_HltSF_error;
	Float_t lead_lepton_ESmearing_error;
	Float_t lead_lepton_EScaling_error;
	Float_t sublead_lepton_IDSF_error;
	Float_t sublead_lepton_IsoSF_error;
	Float_t sublead_lepton_RecoSF_error;
	Float_t sublead_lepton_HltSF_error;
	Float_t sublead_lepton_ESmearing_error;
	Float_t sublead_lepton_EScaling_error;

	Double_t nPU;	///<positive only for MC
	Int_t nPV;
	Int_t njets;

	Float_t weight; ///< this variable takes into accont the product of the global_event_weight and the single object weights   the cross sxn normalization is applied in analysis.cpp
	Float_t pu_weight;	///< the PU weight, equivalent to the absolute value of the global_event_weight

	bool isPassing(tag_t tag, bool makeHists = false);
	bool isPassingLooseCuts(tag_t tag);
	bool isPassingPreselect(bool makeHists = false);

	Selector(const miniTreeEvent& myEvent);
	Selector();

	void SetBranches(TTree* tree);
	void SetBranchAddresses(TTree* tree);


private:
	myElectronCollection electrons;
	myMuonCollection muons;
	myJetCollection jets;


	Float_t global_event_weight; ///< this weight is mcGenWeight * PU_reweight


	void Clear();
	bool _isPassing;
	bool _isPassingLooseCuts;
};



#endif
