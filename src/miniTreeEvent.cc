// #include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "../interface/miniTreeEvent.h"

miniTreeEvent::miniTreeEvent():
	electrons_p4(new std::vector<TLorentzVector>),
	electron_scale_error(new std::vector<Float_t>),
	electron_smearing_sigma(new std::vector<Float_t>),
	electron_smearing_sigma_phi_up(new std::vector<Float_t>),
	electron_smearing_sigma_phi_down(new std::vector<Float_t>),
	electron_smearing_sigma_rho_up(new std::vector<Float_t>),
	electron_smearing_sigma_rho_down(new std::vector<Float_t>),
	electron_r9(new std::vector<Float_t>),
	electron_charge(new std::vector<Int_t>),
	electron_IDSF_central(new std::vector<Float_t>),
	electron_IDSF_error(new std::vector<Float_t>),
	electron_RecoSF_central(new std::vector<Float_t>),
	electron_RecoSF_error(new std::vector<Float_t>),
	electron_HltSF_central(new std::vector<Float_t>),
	electron_HltSF_error(new std::vector<Float_t>),
	muons_p4(new std::vector<TLorentzVector>),
	muon_charge(new std::vector<Int_t>),
	muon_IDSF_central(new std::vector<Float_t>),
	muon_IsoSF_central(new std::vector<Float_t>),
	muon_IDSF_error(new std::vector<Float_t>),
	muon_IsoSF_error(new std::vector<Float_t>),
	jets_p4(new std::vector<TLorentzVector>),
	jec_uncertainty(new std::vector<Float_t>),
	jetResolution(new std::vector<Float_t>),
	JER_sf(new std::vector<Float_t>),
	JER_sf_up(new std::vector<Float_t>),
	JER_sf_down(new std::vector<Float_t>),
	genJetPt(new std::vector<Float_t>),
	genJetMatch(new std::vector<bool>),
	_owningMembers(true)
{

	clear();
}


miniTreeEvent::miniTreeEvent(const miniTreeEvent& otherEvent):
	electrons_p4(new std::vector<TLorentzVector>),
	electron_scale_error(new std::vector<Float_t>),
	electron_smearing_sigma(new std::vector<Float_t>),
	electron_smearing_sigma_phi_up(new std::vector<Float_t>),
	electron_smearing_sigma_phi_down(new std::vector<Float_t>),
	electron_smearing_sigma_rho_up(new std::vector<Float_t>),
	electron_smearing_sigma_rho_down(new std::vector<Float_t>),
	electron_r9(new std::vector<Float_t>),
	electron_charge(new std::vector<Int_t>),
	electron_IDSF_central(new std::vector<Float_t>),
	electron_IDSF_error(new std::vector<Float_t>),
	electron_RecoSF_central(new std::vector<Float_t>),
	electron_RecoSF_error(new std::vector<Float_t>),
	electron_HltSF_central(new std::vector<Float_t>),
	electron_HltSF_error(new std::vector<Float_t>),
	muons_p4(new std::vector<TLorentzVector>),
	muon_charge(new std::vector<Int_t>),
	muon_IDSF_central(new std::vector<Float_t>),
	muon_IsoSF_central(new std::vector<Float_t>),
	muon_IDSF_error(new std::vector<Float_t>),
	muon_IsoSF_error(new std::vector<Float_t>),
	jets_p4(new std::vector<TLorentzVector>),
	jec_uncertainty(new std::vector<Float_t>),
	jetResolution(new std::vector<Float_t>),
	JER_sf(new std::vector<Float_t>),
	JER_sf_up(new std::vector<Float_t>),
	JER_sf_down(new std::vector<Float_t>),
	genJetPt(new std::vector<Float_t>),
	genJetMatch(new std::vector<bool>),
	_owningMembers(true)
{
	clear();
	*electrons_p4 = *(otherEvent.electrons_p4);
	*muons_p4 = *(otherEvent.muons_p4);
	*jets_p4 = *(otherEvent.jets_p4);
	*jec_uncertainty = *(otherEvent.jec_uncertainty);
	*jetResolution = *(otherEvent.jetResolution);
	*JER_sf = *(otherEvent.JER_sf);
	*JER_sf_up = *(otherEvent.JER_sf_up);
	*JER_sf_down = *(otherEvent.JER_sf_down);
	*genJetPt = *(otherEvent.genJetPt);
	*genJetMatch = *(otherEvent.genJetMatch);
	*electron_scale_error = *(otherEvent.electron_scale_error);
	*electron_smearing_sigma = *(otherEvent.electron_smearing_sigma);
	*electron_smearing_sigma_phi_up = *(otherEvent.electron_smearing_sigma_phi_up);
	*electron_smearing_sigma_phi_down = *(otherEvent.electron_smearing_sigma_phi_down);
	*electron_smearing_sigma_rho_up = *(otherEvent.electron_smearing_sigma_rho_up);
	*electron_smearing_sigma_rho_down = *(otherEvent.electron_smearing_sigma_rho_down);
	*electron_r9 = *(otherEvent.electron_r9);
	*electron_charge = *(otherEvent.electron_charge);
	*electron_IDSF_central = *(otherEvent.electron_IDSF_central);
	*electron_IDSF_error = *(otherEvent.electron_IDSF_error);
	*electron_RecoSF_central = *(otherEvent.electron_RecoSF_central);
	*electron_RecoSF_error = *(otherEvent.electron_RecoSF_error);
	*electron_HltSF_central = *(otherEvent.electron_HltSF_central);
	*electron_HltSF_error = *(otherEvent.electron_HltSF_error);
	*muon_charge = *(otherEvent.muon_charge);
	*muon_IDSF_central = *(otherEvent.muon_IDSF_central);
	*muon_IsoSF_central = *(otherEvent.muon_IsoSF_central);
	*muon_IDSF_error = *(otherEvent.muon_IDSF_error);
	*muon_IsoSF_error = *(otherEvent.muon_IsoSF_error);
	_owningMembers = false;

	run = otherEvent.run;
	lumi = otherEvent.lumi;
	event = otherEvent.event;

	nPU = otherEvent.nPU;
	nPV = otherEvent.nPV;
	weight = otherEvent.weight;
	PU_reweight = otherEvent.PU_reweight;

	sprintf(datasetName, "%s", otherEvent.datasetName);
}


void miniTreeEvent::clear()
{
	run = lumi = event = 0;

	electrons_p4->clear();
	muons_p4->clear();
	jets_p4->clear();

	jec_uncertainty->clear();
	jetResolution->clear();
	JER_sf->clear();
	JER_sf_up->clear();
	JER_sf_down->clear();
	genJetPt->clear();
	genJetMatch->clear();
	electron_scale_error->clear();
	electron_smearing_sigma->clear();
	electron_smearing_sigma_phi_up->clear();
	electron_smearing_sigma_phi_down->clear();
	electron_smearing_sigma_rho_up->clear();
	electron_smearing_sigma_rho_down->clear();
	electron_r9->clear();
	electron_charge->clear();
	electron_IDSF_central->clear();
	electron_IDSF_error->clear();
	electron_RecoSF_central->clear();
	electron_RecoSF_error->clear();
	electron_HltSF_central->clear();
	electron_HltSF_error->clear();
	muon_charge->clear();
	muon_IDSF_central->clear();
	muon_IsoSF_central->clear();
	muon_IDSF_error->clear();
	muon_IsoSF_error->clear();

	nPU = -999.;
	nPV = 0.;
	weight = 0.0;
	PU_reweight = 0.0;

}
miniTreeEvent::~miniTreeEvent()
{
	clear();
	delete electrons_p4;
	delete electron_scale_error;
	delete electron_smearing_sigma;
	delete electron_smearing_sigma_phi_up;
	delete electron_smearing_sigma_phi_down;
	delete electron_smearing_sigma_rho_up;
	delete electron_smearing_sigma_rho_down;
	delete electron_r9;
	delete electron_charge;
	delete electron_IDSF_central;
	delete electron_IDSF_error;
	delete electron_RecoSF_central;
	delete electron_RecoSF_error;
	delete electron_HltSF_central;
	delete electron_HltSF_error;
	delete muons_p4;
	delete muon_charge;
	delete muon_IDSF_central;
	delete muon_IsoSF_central;
	delete muon_IDSF_error;
	delete muon_IsoSF_error;
	delete jets_p4;
	delete jec_uncertainty;
	delete jetResolution;
	delete JER_sf;
	delete JER_sf_up;
	delete JER_sf_down;
	delete genJetPt;
	delete genJetMatch;
}

void miniTreeEvent::SetBranches(TTree* tree)
{
	tree->Branch("run", &run);
	tree->Branch("lumi", &lumi);
	tree->Branch("event", &event);
	tree->Branch("datasetName", datasetName, "datasetName/C");

	tree->Branch("electrons_p4", electrons_p4, 32000, -1);
	tree->Branch("muons_p4", muons_p4, 32000, -1);
	tree->Branch("jets_p4", jets_p4, 32000, -1);

	tree->Branch("jec_uncertainty", jec_uncertainty);
	tree->Branch("jetResolution", jetResolution);
	tree->Branch("JER_sf", JER_sf);
	tree->Branch("JER_sf_up", JER_sf_up);
	tree->Branch("JER_sf_down", JER_sf_down);
	tree->Branch("genJetPt", genJetPt);
	tree->Branch("genJetMatch", genJetMatch);
	tree->Branch("electron_scale_error", electron_scale_error);
	tree->Branch("electron_smearing_sigma", electron_smearing_sigma);
	tree->Branch("electron_smearing_sigma_phi_up", electron_smearing_sigma_phi_up);
	tree->Branch("electron_smearing_sigma_phi_down", electron_smearing_sigma_phi_down);
	tree->Branch("electron_smearing_sigma_rho_up", electron_smearing_sigma_rho_up);
	tree->Branch("electron_smearing_sigma_rho_down", electron_smearing_sigma_rho_down);
	tree->Branch("electron_r9", electron_r9);
	tree->Branch("electron_charge", electron_charge);
	tree->Branch("muon_charge", muon_charge);
	tree->Branch("muon_IDSF_central", muon_IDSF_central);
	tree->Branch("muon_IsoSF_central", muon_IsoSF_central);
	tree->Branch("muon_IDSF_error", muon_IDSF_error);
	tree->Branch("muon_IsoSF_error", muon_IsoSF_error);

	tree->Branch("nPU", &nPU);
	tree->Branch("nPV", &nPV);
	tree->Branch("weight", &weight);
	tree->Branch("PU_reweight", &PU_reweight);

}

void miniTreeEvent::SetBranchAddresses(TChain* tree)
{

	delete electrons_p4;
	delete muons_p4;
	delete jets_p4;

	delete jec_uncertainty;
	delete jetResolution;
	delete JER_sf;
	delete JER_sf_up;
	delete JER_sf_down;
	delete genJetPt;
	delete genJetMatch;
	delete electron_scale_error;
	delete electron_smearing_sigma;
	delete electron_smearing_sigma_phi_up;
	delete electron_smearing_sigma_phi_down;
	delete electron_smearing_sigma_rho_up;
	delete electron_smearing_sigma_rho_down;
	delete electron_r9;
	delete electron_charge;
	delete muon_charge;
	delete muon_IDSF_central;
	delete muon_IsoSF_central;
	delete muon_IDSF_error;
	delete muon_IsoSF_error;

	_owningMembers = false;

	electrons_p4 = 0;
	muons_p4 = 0;
	jets_p4 = 0;

	jec_uncertainty = 0;
	jetResolution = 0;
	JER_sf = 0;
	JER_sf_up = 0;
	JER_sf_down = 0;
	genJetPt = 0;
	genJetMatch = 0;
	electron_scale_error = 0;
	electron_smearing_sigma = 0;
	electron_smearing_sigma_phi_up = 0;
	electron_smearing_sigma_phi_down = 0;
	electron_smearing_sigma_rho_up = 0;
	electron_smearing_sigma_rho_down = 0;
	electron_r9 = 0;
	electron_charge = 0;
	muon_charge = 0;
	muon_IDSF_central = 0;
	muon_IsoSF_central = 0;
	muon_IDSF_error = 0;
	muon_IsoSF_error = 0;

	tree->SetBranchAddress("run", &run);
	tree->SetBranchAddress("lumi", &lumi);
	tree->SetBranchAddress("event", &event);
	tree->SetBranchAddress("datasetName", datasetName);

	tree->SetBranchAddress("electrons_p4", &electrons_p4);
	tree->SetBranchAddress("muons_p4", &muons_p4);
	tree->SetBranchAddress("jets_p4", &jets_p4);

	tree->SetBranchAddress("jec_uncertainty", &jec_uncertainty);
	tree->SetBranchAddress("jetResolution", &jetResolution);
	tree->SetBranchAddress("JER_sf", &JER_sf);
	tree->SetBranchAddress("JER_sf_up", &JER_sf_up);
	tree->SetBranchAddress("JER_sf_down", &JER_sf_down);
	tree->SetBranchAddress("genJetPt", &genJetPt);
	tree->SetBranchAddress("genJetMatch", &genJetMatch);

	tree->SetBranchAddress("electron_scale_error", &electron_scale_error);
	tree->SetBranchAddress("electron_smearing_sigma", &electron_smearing_sigma);
	tree->SetBranchAddress("electron_smearing_sigma_phi_up", &electron_smearing_sigma_phi_up);
	tree->SetBranchAddress("electron_smearing_sigma_phi_down", &electron_smearing_sigma_phi_down);
	tree->SetBranchAddress("electron_smearing_sigma_rho_up", &electron_smearing_sigma_rho_up);
	tree->SetBranchAddress("electron_smearing_sigma_rho_down", &electron_smearing_sigma_rho_down);
	tree->SetBranchAddress("electron_r9", &electron_r9);
	tree->SetBranchAddress("electron_charge", &electron_charge);

	tree->SetBranchAddress("muon_charge", &muon_charge);
	tree->SetBranchAddress("muon_IDSF_central", &muon_IDSF_central);
	tree->SetBranchAddress("muon_IsoSF_central", &muon_IsoSF_central);
	tree->SetBranchAddress("muon_IDSF_error", &muon_IDSF_error);
	tree->SetBranchAddress("muon_IsoSF_error", &muon_IsoSF_error);

	tree->SetBranchAddress("nPU", &nPU);
	tree->SetBranchAddress("nPV", &nPV);
	tree->SetBranchAddress("weight", &weight);
	tree->SetBranchAddress("PU_reweight", &PU_reweight);

}
