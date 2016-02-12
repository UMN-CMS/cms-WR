#include <iostream>
#include <algorithm>
#include <vector>
// #include "ExoAnalysis/cmsWR/interface/Selector.h"
#include "../interface/Selector.h"
#include "DataFormats/Math/interface/deltaR.h"

float dR_TLV(TLorentzVector t1, TLorentzVector t2)
{
	return deltaR(t1.Eta(), t1.Phi(), t2.Eta(), t2.Phi());
}

void goodJets(myJetCollection *evJets, myJetCollection *selJets)
{
	for(auto j : *evJets) {
		if(j.p4.Pt() > 40 && fabs(j.p4.Eta()) < 2.4)
			selJets->push_back(j);
	}
}

void goodEles(myElectronCollection *evEles, myElectronCollection *selEles)
{
	for(auto e : *evEles) {
		if(e.p4.Pt() > 40 && fabs(e.p4.Eta()) < 2.4 && (fabs(e.p4.Eta()) < 1.4222 || fabs(e.p4.Eta()) > 1.566))
			selEles->push_back(e);
	}
}

void goodMuons(myMuonCollection *evMuons, myMuonCollection *selMuons)
{
	for(auto m : *evMuons) {
		if(m.p4.Pt() > 40 && fabs(m.p4.Eta()) < 2.4)
			selMuons->push_back(m);
	}
}


Selector::Selector(const miniTreeEvent& myEvent) :
	WR_mass(-1),
	dilepton_mass(-1)
{
	datasetName = myEvent.datasetName;

	int nele = myEvent.electrons_p4->size();
	for(int i = 0; i < nele; i++) {
		myElectron ele;
		ele.p4 = myEvent.electrons_p4->at(i);
		ele.scale = myEvent.electron_scale->at(i);
		ele.smearing = myEvent.electron_smearing->at(i);
		ele.charge = myEvent.electron_charge->at(i);
		ele.weight = 1.;
		electrons.push_back(ele);
	}
	int nmu = myEvent.muons_p4->size();
	for(int i = 0; i < nmu; i++) {
		myMuon mu;
		mu.p4 = myEvent.muons_p4->at(i);
		mu.IDSF = myEvent.muon_IDSF_central->at(i);
		mu.IsoSF = myEvent.muon_IsoSF_central->at(i);
		mu.IDSF_error = myEvent.muon_IDSF_error->at(i);
		mu.IsoSF_error = myEvent.muon_IsoSF_error->at(i);
		mu.charge = myEvent.muon_charge->at(i);
		mu.weight = mu.IDSF * mu.IsoSF;
		muons.push_back(mu);
	}
	int njet = myEvent.jets_p4->size();
	for(int i = 0; i < njet; i++) {
		myJet jet;
		jet.p4 = myEvent.jets_p4->at(i);
		jet.jec_uncertainty = myEvent.jec_uncertainty->at(i);
		jet.weight = 1.;
		jets.push_back(jet);
	}

	global_event_weight = myEvent.weight * myEvent.PU_reweight;

//	Clear();
}

Selector::Selector()
{
	Clear();
}

void Selector::Clear()
{
	WR_mass = dilepton_mass = weight = 0.0;
}

bool Selector::isPassing(tag_t tag)
{

	_isPassing = false;
	WR_mass = -1;
	TLorentzVector lead_lepton_p4, sublead_lepton_p4, lead_jet_p4, sublead_jet_p4;

	myJetCollection gJets;
	myElectronCollection gEles;
	myMuonCollection gMuons;

	// Basic Kinematic cuts
	goodJets(&jets, &gJets);
	goodEles(&electrons, &gEles);
	goodMuons(&muons, &gMuons);

	std::sort(gJets.begin(), gJets.end(),
	[](myJet const & a, myJet const & b) {
		return a.p4.Pt() > b.p4.Pt();
	});
	std::sort(gEles.begin(), gEles.end(),
	[](myElectron const & a, myElectron const & b) {
		return a.p4.Pt() > b.p4.Pt();
	});
	std::sort(gMuons.begin(), gMuons.end(),
	[](myMuon const & a, myMuon const & b) {
		return a.p4.Pt() > b.p4.Pt();
	});

	jets = gJets;
	electrons = gEles;
	muons = gMuons;

	// Assert at least 2 good jets
	if(jets.size() < 2) {
		return false;
	}

	lead_jet_pt = jets[0].p4.Pt();
	lead_jet_eta = jets[0].p4.Eta();
	lead_jet_phi = jets[0].p4.Phi();
	lead_jet_weight = 1.0;

	sublead_jet_pt = jets[1].p4.Pt();
	sublead_jet_eta = jets[1].p4.Eta();
	sublead_jet_phi = jets[1].p4.Phi();
	lead_jet_weight = 1.0;

	if(tag == 0) { // EEJJ Channel
		// Assert at least 2 good leptons
		if(electrons.size() < 2) {
			return false;
		}

		lead_lepton_p4 = electrons[0].p4;
		sublead_lepton_p4 = electrons[1].p4;

		lead_lepton_weight = electrons[0].weight;
		sublead_lepton_weight = electrons[1].weight;

	} else if(tag == 1) { // MuMuJJ Channel
		// Assert at least 2 good leptons
		if(muons.size() < 2) {
			return false;
		}

		lead_lepton_p4 = muons[0].p4;
		sublead_lepton_p4 = muons[1].p4;

		lead_lepton_weight = muons[0].weight;
		sublead_lepton_weight = muons[1].weight;

	} else if(tag == 2) { // EMuJJ Channel
		// Assert at least 2 good leptons
		if(electrons.size() < 1 || muons.size() < 1) {
			return false;
		}

		// check which is the leading, which the subleading
		if(electrons[0].p4.Pt() > muons[0].p4.Pt()) { // e > mu
			lead_lepton_p4 = electrons[0].p4;
			sublead_lepton_p4 = muons[0].p4;

			lead_lepton_weight = electrons[0].weight;
			sublead_lepton_weight = muons[0].weight;

		} else {

			sublead_lepton_p4 = electrons[0].p4;
			sublead_lepton_weight = electrons[0].weight;

			lead_lepton_p4 = muons[0].p4;
			lead_lepton_weight = muons[0].weight;
		}
	}

	// check eta and pt cuts
	if(lead_lepton_p4.Pt() < 60) return false;
	if(sublead_lepton_p4.Pt() < 50) return false;

	if(dR_TLV(lead_lepton_p4, gJets[0].p4) < 0.4) return false;
	if(dR_TLV(lead_lepton_p4, gJets[1].p4) < 0.4) return false;
	if(dR_TLV(sublead_lepton_p4, gJets[0].p4) < 0.4) return false;
	if(dR_TLV(sublead_lepton_p4, gJets[1].p4) < 0.4) return false;

	lead_lepton_pt = lead_lepton_p4.Pt();
	lead_lepton_eta = lead_lepton_p4.Eta();
	lead_lepton_phi = lead_lepton_p4.Phi();

	sublead_lepton_pt = sublead_lepton_p4.Pt();
	sublead_lepton_eta = sublead_lepton_p4.Eta();
	sublead_lepton_phi = sublead_lepton_p4.Phi();


	// Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
	WR_mass = (lead_lepton_p4 + sublead_lepton_p4 + gJets[0].p4 + gJets[1].p4).M();
	weight = lead_lepton_weight * sublead_lepton_weight * lead_jet_weight * sublead_jet_weight * global_event_weight;

	dilepton_mass = (lead_lepton_p4 + sublead_lepton_p4).M();

	_isPassing = true;
	return _isPassing;

}

void Selector::SetBranches(TTree* tree)
{

	tree->Branch("lead_lepton_pt", &lead_lepton_pt);
	tree->Branch("sublead_lepton_pt", &sublead_lepton_pt);
	tree->Branch("lead_lepton_eta", &lead_lepton_eta);
	tree->Branch("sublead_lepton_eta", &sublead_lepton_eta);
	tree->Branch("lead_lepton_phi", &lead_lepton_phi);
	tree->Branch("sublead_lepton_phi", &sublead_lepton_phi);
	tree->Branch("lead_jet_pt", &lead_jet_pt);
	tree->Branch("sublead_jet_pt", &sublead_jet_pt);
	tree->Branch("lead_jet_eta", &lead_jet_eta);
	tree->Branch("sublead_jet_eta", &sublead_jet_eta);
	tree->Branch("lead_jet_phi", &lead_jet_phi);
	tree->Branch("sublead_jet_phi", &sublead_jet_phi);

	tree->Branch("weight", &weight);
	tree->Branch("WR_mass", &WR_mass);
	tree->Branch("dilepton_mass", &dilepton_mass);

}

void Selector::SetBranchAddresses(TTree* tree)
{

	tree->SetBranchAddress("lead_lepton_pt", &lead_lepton_pt);
	tree->SetBranchAddress("lead_lepton_eta", &lead_lepton_eta);
	tree->SetBranchAddress("lead_lepton_phi", &lead_lepton_phi);
	tree->SetBranchAddress("sublead_lepton_pt", &sublead_lepton_pt);
	tree->SetBranchAddress("sublead_lepton_eta", &sublead_lepton_eta);
	tree->SetBranchAddress("sublead_lepton_phi", &sublead_lepton_phi);

	tree->SetBranchAddress("lead_jet_pt", &lead_jet_pt);
	tree->SetBranchAddress("lead_jet_eta", &lead_jet_eta);
	tree->SetBranchAddress("lead_jet_phi", &lead_jet_phi);
	tree->SetBranchAddress("sublead_jet_pt", &sublead_jet_pt);
	tree->SetBranchAddress("sublead_jet_eta", &sublead_jet_eta);
	tree->SetBranchAddress("sublead_jet_phi", &sublead_jet_phi);

	tree->SetBranchAddress("weight", &weight);
	tree->SetBranchAddress("WR_mass", &WR_mass);
	tree->SetBranchAddress("dilepton_mass", &dilepton_mass);

}

