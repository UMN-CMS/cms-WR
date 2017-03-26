#include <iostream>
#include <algorithm>
#include <vector>
// #include "ExoAnalysis/cmsWR/interface/Selector.h"
#include "../interface/Selector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "../interface/SelectorHist.h"
//#define DEBUGG

float dR_TLV(TLorentzVector t1, TLorentzVector t2)
{
	return deltaR(t1.Eta(), t1.Phi(), t2.Eta(), t2.Phi());
}

//default jet pt cut is 40, only other jet cuts in this file are in isPassingPreselect
void goodJets(myJetCollection *evJets, myJetCollection *selJets)
{
	for(auto j : *evJets) {
		if(j.p4.Pt() > 40 && fabs(j.p4.Eta()) < 2.4)
			selJets->push_back(j);
	}
}

void goodJetsLooseCuts(myJetCollection *evJets, myJetCollection *selJets)
{
	for(auto j : *evJets) {
		if(fabs(j.p4.Eta()) < 2.4 ) selJets->push_back(j);
	}
}

void goodEles(myElectronCollection *evEles, myElectronCollection *selEles)
{
	for(auto e : *evEles) {
		if(e.p4.Pt() > 40 && fabs(e.p4.Eta()) < 2.4 && (fabs(e.p4.Eta()) < 1.4222 || fabs(e.p4.Eta()) > 1.566))
			selEles->push_back(e);
	}
}

void goodElesLooseCuts(myElectronCollection *evEles, myElectronCollection *selEles)
{
	for(auto e : *evEles) {
		if(fabs(e.p4.Eta()) < 2.4 && (fabs(e.p4.Eta()) < 1.4222 || fabs(e.p4.Eta()) > 1.566)) selEles->push_back(e);
	}
}

void goodMuons(myMuonCollection *evMuons, myMuonCollection *selMuons)
{
	for(auto m : *evMuons) {
		if(m.p4.Pt() > 40 && fabs(m.p4.Eta()) < 2.4)
			selMuons->push_back(m);
	}
}

void goodMuonsLooseCuts(myMuonCollection *evMuons, myMuonCollection *selMuons)
{
	for(auto m : *evMuons) {
		if(fabs(m.p4.Eta()) < 2.4) selMuons->push_back(m);
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
		ele.r9 = myEvent.electron_r9->at(i);
		ele.smearing_error = 0.;	///<temporary
		ele.scale_error = 0.;		///<temporary
		ele.IDSF = myEvent.electron_IDSF_central->at(i);
		ele.IDSF_error = myEvent.electron_IDSF_error->at(i);
		ele.RecoSF = myEvent.electron_RecoSF_central->at(i);
		ele.RecoSF_error = myEvent.electron_RecoSF_error->at(i);
		ele.HltSF = myEvent.electron_HltSF_central->at(i);
		ele.HltSF_error = myEvent.electron_HltSF_error->at(i);
		ele.weight = (ele.IDSF) * (ele.RecoSF) * (ele.HltSF);
		ele.passedID = myEvent.electron_passedHEEP->at(i);
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
		mu.passedID = myEvent.muon_passedIDIso->at(i);
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

	nPV = myEvent.nPV;
	nPU = myEvent.nPU;
	global_event_weight = (myEvent.weight > 0 ? 1 : -1) * myEvent.PU_reweight;
#ifdef DEBUGG
	std::cout << "global_event_weight=\t" << global_event_weight << std::endl;
#endif

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

bool Selector::isPassingLooseCuts(tag_t tag)
{
	_isPassingLooseCuts = false;
	WR_mass = -1, lead_lepton_r9 = -1, sublead_lepton_r9 = -1;
	TLorentzVector lead_lepton_p4, sublead_lepton_p4, lead_jet_p4, sublead_jet_p4;
	lead_lepton_IDSF_error = -9;
	lead_lepton_RecoSF_error = -9;
	lead_lepton_HltSF_error = -9;
	lead_lepton_ESmearing_error = -9;
	lead_lepton_EScaling_error = -9;
	sublead_lepton_IDSF_error = -9;
	sublead_lepton_RecoSF_error = -9;
	sublead_lepton_HltSF_error = -9;
	sublead_lepton_ESmearing_error = -9;
	sublead_lepton_EScaling_error = -9;
	lead_lepton_IsoSF_error = -9;
	sublead_lepton_IsoSF_error = -9;

	myJetCollection gJets;
	myElectronCollection gEles;
	myMuonCollection gMuons;

	// Basic Kinematic cuts (only eta cuts applied in good...LooseCuts() )
	goodJetsLooseCuts(&jets, &gJets);
	goodElesLooseCuts(&electrons, &gEles);
	goodMuonsLooseCuts(&muons, &gMuons);

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

	if(tag == EE) { // EEJJ Channel
		// Assert at least 2 good leptons
		if(electrons.size() < 2) {
			return false;
		}

		lead_lepton_p4 = electrons[0].p4;
		sublead_lepton_p4 = electrons[1].p4;

		lead_lepton_weight = electrons[0].weight;
		sublead_lepton_weight = electrons[1].weight;

		lead_lepton_IDSF_error = electrons[0].IDSF_error;
		lead_lepton_RecoSF_error = electrons[0].RecoSF_error;
		lead_lepton_HltSF_error = electrons[0].HltSF_error;
		lead_lepton_ESmearing_error = electrons[0].smearing_error;
		lead_lepton_EScaling_error = electrons[0].scale_error;

		sublead_lepton_IDSF_error = electrons[1].IDSF_error;
		sublead_lepton_RecoSF_error = electrons[1].RecoSF_error;
		sublead_lepton_HltSF_error = electrons[1].HltSF_error;
		sublead_lepton_ESmearing_error = electrons[1].smearing_error;
		sublead_lepton_EScaling_error = electrons[1].scale_error;

		lead_lepton_r9 = electrons[0].r9;
		sublead_lepton_r9 = electrons[1].r9;
	} else if(tag == MuMu) { // MuMuJJ Channel
		// Assert at least 2 good leptons
		if(muons.size() < 2) {
			return false;
		}

		lead_lepton_p4 = muons[0].p4;
		sublead_lepton_p4 = muons[1].p4;

		lead_lepton_IDSF_error = muons[0].IDSF_error;
		lead_lepton_IsoSF_error = muons[0].IsoSF_error;

		sublead_lepton_IDSF_error = muons[1].IDSF_error;
		sublead_lepton_IsoSF_error = muons[1].IsoSF_error;

		lead_lepton_weight = muons[0].weight;
		sublead_lepton_weight = muons[1].weight;
	} else if(tag == EMu) { // EMuJJ Channel
		// Assert at least 2 good leptons
		if(electrons.size() < 1 || muons.size() < 1) {
			return false;
		}

//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
		// check which is the leading, which the subleading
		if(electrons[0].p4.Pt() > muons[0].p4.Pt()) { // e > mu
			lead_lepton_p4 = electrons[0].p4;
			sublead_lepton_p4 = muons[0].p4;

			lead_lepton_IDSF_error = electrons[0].IDSF_error;
			lead_lepton_RecoSF_error = electrons[0].RecoSF_error;
			lead_lepton_HltSF_error = electrons[0].HltSF_error;
			lead_lepton_ESmearing_error = electrons[0].smearing_error;
			lead_lepton_EScaling_error = electrons[0].scale_error;

			sublead_lepton_IDSF_error = muons[0].IDSF_error;
			sublead_lepton_IsoSF_error = muons[0].IsoSF_error;

			lead_lepton_weight = electrons[0].weight;
			sublead_lepton_weight = muons[0].weight;
			lead_lepton_r9 = electrons[0].r9;
		} else {

			sublead_lepton_p4 = electrons[0].p4;
			sublead_lepton_weight = electrons[0].weight;
			sublead_lepton_IDSF_error = electrons[0].IDSF_error;
			sublead_lepton_RecoSF_error = electrons[0].RecoSF_error;
			sublead_lepton_HltSF_error = electrons[0].HltSF_error;
			sublead_lepton_ESmearing_error = electrons[0].smearing_error;
			sublead_lepton_EScaling_error = electrons[0].scale_error;

			lead_lepton_IDSF_error = muons[0].IDSF_error;
			lead_lepton_IsoSF_error = muons[0].IsoSF_error;
			lead_lepton_p4 = muons[0].p4;
			lead_lepton_weight = muons[0].weight;
			sublead_lepton_r9 = electrons[0].r9;
		}
	}

	//lepton pt cuts necessitated by dytagandprobe triggers
	if(lead_lepton_p4.Pt() < 35) return false;
	if(sublead_lepton_p4.Pt() < 35) return false;


	//defaults if no jets are found in the event
	lead_jet_jec_unc = -10;
	sublead_jet_jec_unc = -10;
	lead_jet_pt = -9;
	lead_jet_eta = -6;
	lead_jet_phi = -6;
	lead_jet_weight = 1.0;
	sublead_jet_weight = 1.0;
	sublead_jet_pt = -9;
	sublead_jet_eta = -6;
	sublead_jet_phi = -6;
	dR_leadlepton_leadjet = 9, dR_leadlepton_subleadjet = 9, dR_subleadlepton_leadjet = 9, dR_subleadlepton_subleadjet = 9;

	njets = jets.size();

	if(jets.size() == 1) {
		lead_jet_pt = jets[0].p4.Pt();
		lead_jet_eta = jets[0].p4.Eta();
		lead_jet_phi = jets[0].p4.Phi();
		lead_jet_weight = 1.0;
		lead_jet_jec_unc = jets[0].jec_uncertainty;
		dR_leadlepton_leadjet = dR_TLV(lead_lepton_p4, jets[0].p4);
		dR_subleadlepton_leadjet = dR_TLV(sublead_lepton_p4, jets[0].p4);
	}//one jet in event

	if(jets.size() > 1) {
		lead_jet_pt = jets[0].p4.Pt();
		lead_jet_eta = jets[0].p4.Eta();
		lead_jet_phi = jets[0].p4.Phi();
		lead_jet_weight = 1.0;
		sublead_jet_weight = 1.0;
		lead_jet_jec_unc = jets[0].jec_uncertainty;
		sublead_jet_jec_unc = jets[1].jec_uncertainty;

		sublead_jet_pt = jets[1].p4.Pt();
		sublead_jet_eta = jets[1].p4.Eta();
		sublead_jet_phi = jets[1].p4.Phi();

		dR_leadlepton_leadjet = dR_TLV(lead_lepton_p4, jets[0].p4);
		dR_subleadlepton_leadjet = dR_TLV(sublead_lepton_p4, jets[0].p4);
		dR_leadlepton_subleadjet = dR_TLV(lead_lepton_p4, jets[1].p4);
		dR_subleadlepton_subleadjet = dR_TLV(sublead_lepton_p4, jets[1].p4);
	}//two or more jets in the event


	lead_lepton_pt = lead_lepton_p4.Pt();
	lead_lepton_eta = lead_lepton_p4.Eta();
	lead_lepton_phi = lead_lepton_p4.Phi();

	sublead_lepton_pt = sublead_lepton_p4.Pt();
	sublead_lepton_eta = sublead_lepton_p4.Eta();
	sublead_lepton_phi = sublead_lepton_p4.Phi();


	// Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
	//WR_mass = (lead_lepton_p4 + sublead_lepton_p4 + gJets[0].p4 + gJets[1].p4).M();
	//weight = lead_lepton_weight * sublead_lepton_weight * lead_jet_weight * sublead_jet_weight * global_event_weight;
	WR_mass = -1;
	weight = lead_lepton_weight * sublead_lepton_weight * global_event_weight;
	pu_weight = fabs(global_event_weight);


#ifdef DEBUGG
	std::cout << "weight (global_event_weight times single object weights)=\t" << weight << std::endl;
#endif

	dilepton_mass = (lead_lepton_p4 + sublead_lepton_p4).M();
	zPt = (lead_lepton_p4 + sublead_lepton_p4).Pt();
	zEta = (lead_lepton_p4 + sublead_lepton_p4).Eta();
	zPhi = (lead_lepton_p4 + sublead_lepton_p4).Phi();
	if(dilepton_mass < 60.0 || dilepton_mass > 120.0) return false;
	//if(dilepton_mass < 60.0) return false;	///<for DYMadIncl and DYMadHT study

	_isPassingLooseCuts = true;
	return _isPassingLooseCuts;

}

bool Selector::isPassing(tag_t tag, bool makeHists)
{

	/**/
	enum det_t {
		DET_ENDCAP,
		DET_BARREL,
		DET_GAP,
	} lead_det, sublead_det;

	enum pair_t {
		P_EE,
		P_BB,
		P_EB,
		P_GAP,
	} pair;
	/**/

	_isPassing = false;
	WR_mass = -1, lead_lepton_r9 = -1, sublead_lepton_r9 = -1;
	TLorentzVector lead_lepton_p4, sublead_lepton_p4, lead_jet_p4, sublead_jet_p4;
	lead_lepton_IDSF_error = -9;
	lead_lepton_RecoSF_error = -9;
	lead_lepton_HltSF_error = -9;
	lead_lepton_ESmearing_error = -9;
	lead_lepton_EScaling_error = -9;
	sublead_lepton_IDSF_error = -9;
	sublead_lepton_RecoSF_error = -9;
	sublead_lepton_HltSF_error = -9;
	sublead_lepton_ESmearing_error = -9;
	sublead_lepton_EScaling_error = -9;
	lead_lepton_IsoSF_error = -9;
	sublead_lepton_IsoSF_error = -9;

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

	if (makeHists) sel::hists("njets", 10, 0, 10)->Fill(jets.size());
	// Assert at least 2 good jets
	if(jets.size() < 2) {
		return false;
	}
	if (makeHists) sel::hists("njets_cut", 10, 0, 10)->Fill(jets.size());

	njets = jets.size();

	lead_jet_pt = jets[0].p4.Pt();
	lead_jet_eta = jets[0].p4.Eta();
	lead_jet_phi = jets[0].p4.Phi();
	lead_jet_weight = 1.0;
	sublead_jet_weight = 1.0;
	lead_jet_jec_unc = jets[0].jec_uncertainty;
	sublead_jet_jec_unc = jets[1].jec_uncertainty;

	sublead_jet_pt = jets[1].p4.Pt();
	sublead_jet_eta = jets[1].p4.Eta();
	sublead_jet_phi = jets[1].p4.Phi();
	lead_jet_weight = 1.0;

	if(tag == EE) { // EEJJ Channel
		// Assert at least 2 good leptons
		if (makeHists) sel::hists("nlep", 10, 0, 10)->Fill(electrons.size());
		if(electrons.size() < 2) {
			return false;
		}
		if (makeHists) sel::hists("nlep_cut", 10, 0, 10)->Fill(electrons.size());

		lead_lepton_p4 = electrons[0].p4;
		sublead_lepton_p4 = electrons[1].p4;

		lead_lepton_passedID = electrons[0].passedID;
		sublead_lepton_passedID = electrons[1].passedID;
		if(lead_lepton_passedID > 0 && sublead_lepton_passedID > 0) return false;	//skip evt if both selected leptons passed ID

		lead_lepton_weight = electrons[0].weight;
		sublead_lepton_weight = electrons[1].weight;

		//update lepton weights to account for probability of jet which passes very loose ID
		//to fake a lepton which passes full ID
		if(lead_lepton_passedID > 0){
			//sublead lepton is a jet reconstructed as a lepton

		}
		else if(sublead_lepton_passedID > 0){
			//lead lepton is a jet reconstructed as a lepton

		}
		else{
			//lead and sublead leptons are jets reco'd as leptons

		}

		lead_lepton_IDSF_error = electrons[0].IDSF_error;
		lead_lepton_RecoSF_error = electrons[0].RecoSF_error;
		lead_lepton_HltSF_error = electrons[0].HltSF_error;
		lead_lepton_ESmearing_error = electrons[0].smearing_error;
		lead_lepton_EScaling_error = electrons[0].scale_error;

		sublead_lepton_IDSF_error = electrons[1].IDSF_error;
		sublead_lepton_RecoSF_error = electrons[1].RecoSF_error;
		sublead_lepton_HltSF_error = electrons[1].HltSF_error;
		sublead_lepton_ESmearing_error = electrons[1].smearing_error;
		sublead_lepton_EScaling_error = electrons[1].scale_error;

		lead_lepton_r9 = electrons[0].r9;
		sublead_lepton_r9 = electrons[1].r9;
	} else if(tag == MuMu) { // MuMuJJ Channel
		// Assert at least 2 good leptons
		if (makeHists) sel::hists("nlep", 10, 0, 10)->Fill(muons.size());
		if(muons.size() < 2) {
			return false;
		}
		if (makeHists) sel::hists("nlep_cut", 10, 0, 10)->Fill(muons.size());

		lead_lepton_p4 = muons[0].p4;
		sublead_lepton_p4 = muons[1].p4;

		lead_lepton_passedID = muons[0].passedID;
		sublead_lepton_passedID = muons[1].passedID;
		if(lead_lepton_passedID > 0 && sublead_lepton_passedID > 0) return false;	//skip evt if both selected leptons passed ID

		lead_lepton_IDSF_error = muons[0].IDSF_error;
		lead_lepton_IsoSF_error = muons[0].IsoSF_error;

		sublead_lepton_IDSF_error = muons[1].IDSF_error;
		sublead_lepton_IsoSF_error = muons[1].IsoSF_error;

		lead_lepton_weight = muons[0].weight;
		sublead_lepton_weight = muons[1].weight;

	} else if(tag == EMu) { // EMuJJ Channel
		// Assert at least 2 good leptons
		if (makeHists) sel::hists("nlep", 10, 0, 10)->Fill(muons.size() + electrons.size());
		if(electrons.size() < 1 || muons.size() < 1) {
			return false;
		}
		if (makeHists) sel::hists("nlep_cut", 10, 0, 10)->Fill(muons.size() + electrons.size());

//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
		// check which is the leading, which the subleading
		if(electrons[0].p4.Pt() > muons[0].p4.Pt()) { // e > mu
			lead_lepton_p4 = electrons[0].p4;
			sublead_lepton_p4 = muons[0].p4;

			lead_lepton_IDSF_error = electrons[0].IDSF_error;
			lead_lepton_RecoSF_error = electrons[0].RecoSF_error;
			lead_lepton_HltSF_error = electrons[0].HltSF_error;
			lead_lepton_ESmearing_error = electrons[0].smearing_error;
			lead_lepton_EScaling_error = electrons[0].scale_error;

			sublead_lepton_IDSF_error = muons[0].IDSF_error;
			sublead_lepton_IsoSF_error = muons[0].IsoSF_error;

			lead_lepton_weight = electrons[0].weight;
			sublead_lepton_weight = muons[0].weight;

			lead_lepton_r9 = electrons[0].r9;
		} else { // mu > e

			sublead_lepton_p4 = electrons[0].p4;
			sublead_lepton_weight = electrons[0].weight;
			sublead_lepton_IDSF_error = electrons[0].IDSF_error;
			sublead_lepton_RecoSF_error = electrons[0].RecoSF_error;
			sublead_lepton_HltSF_error = electrons[0].HltSF_error;
			sublead_lepton_ESmearing_error = electrons[0].smearing_error;
			sublead_lepton_EScaling_error = electrons[0].scale_error;

			lead_lepton_IDSF_error = muons[0].IDSF_error;
			lead_lepton_IsoSF_error = muons[0].IsoSF_error;

			lead_lepton_p4 = muons[0].p4;
			lead_lepton_weight = muons[0].weight;

			sublead_lepton_r9 = electrons[0].r9;
		}
	}//end EMu channel
	/**/
	if(fabs(lead_lepton_p4.Eta()) > 1.566) lead_det = DET_ENDCAP;
	else if(fabs(lead_lepton_p4.Eta()) > 1.4222) lead_det = DET_GAP;
	else lead_det = DET_BARREL;
	if(fabs(sublead_lepton_p4.Eta()) > 1.566) sublead_det = DET_ENDCAP;
	else if(fabs(sublead_lepton_p4.Eta()) > 1.4222) sublead_det = DET_GAP;
	else sublead_det = DET_BARREL;

	if(lead_det == DET_GAP || sublead_det == DET_GAP) pair = P_GAP;
	else if( lead_det == DET_BARREL && sublead_det == DET_BARREL) pair = P_BB;
	else if( lead_det == DET_ENDCAP && sublead_det == DET_ENDCAP) pair = P_EE;
	else pair = P_EB;
	/**/

	// check eta and pt cuts
	if (makeHists) sel::hists("lead_lepton_pt", 100, 0, 200)->Fill(lead_lepton_p4.Pt());
	if(lead_lepton_p4.Pt() < 60) return false;	//default is 60
	if (makeHists) sel::hists("lead_lepton_pt_cut", 100, 0, 200)->Fill(lead_lepton_p4.Pt());
	if (makeHists) sel::hists("sublead_lepton_pt", 100, 0, 200)->Fill(sublead_lepton_p4.Pt());
	if(sublead_lepton_p4.Pt() < 53) return false;	//default changed to 53 from 50 on January 26 2017, after request from conveners to be above trigger threshold
	if (makeHists) sel::hists("sublead_lepton_pt_cut", 100, 0, 200)->Fill(sublead_lepton_p4.Pt());

	if (makeHists) sel::hists("dr", 100, 0, 5)->Fill(dR_TLV(lead_lepton_p4, gJets[0].p4));
	if (makeHists) sel::hists("dr", 100, 0, 5)->Fill(dR_TLV(lead_lepton_p4, gJets[1].p4));
	if (makeHists) sel::hists("dr", 100, 0, 5)->Fill(dR_TLV(sublead_lepton_p4, gJets[0].p4));
	if (makeHists) sel::hists("dr", 100, 0, 5)->Fill(dR_TLV(sublead_lepton_p4, gJets[1].p4));
	if (makeHists) sel::hists("dr_count", 1, 0, 1)->Fill(0);
	if(dR_TLV(lead_lepton_p4, gJets[0].p4) < 0.4) return false;
	if(dR_TLV(lead_lepton_p4, gJets[1].p4) < 0.4) return false;
	if(dR_TLV(sublead_lepton_p4, gJets[0].p4) < 0.4) return false;
	if(dR_TLV(sublead_lepton_p4, gJets[1].p4) < 0.4) return false;
	if (makeHists) sel::hists("dr_cut", 100, 0, 5)->Fill(dR_TLV(lead_lepton_p4, gJets[0].p4));
	if (makeHists) sel::hists("dr_cut", 100, 0, 5)->Fill(dR_TLV(lead_lepton_p4, gJets[1].p4));
	if (makeHists) sel::hists("dr_cut", 100, 0, 5)->Fill(dR_TLV(sublead_lepton_p4, gJets[0].p4));
	if (makeHists) sel::hists("dr_cut", 100, 0, 5)->Fill(dR_TLV(sublead_lepton_p4, gJets[1].p4));
	if (makeHists) sel::hists("dr_count_cut", 1, 0, 1)->Fill(0);

	dR_leadlepton_leadjet = dR_TLV(lead_lepton_p4, jets[0].p4);
	dR_subleadlepton_leadjet = dR_TLV(sublead_lepton_p4, jets[0].p4);
	dR_leadlepton_subleadjet = dR_TLV(lead_lepton_p4, jets[1].p4);
	dR_subleadlepton_subleadjet = dR_TLV(sublead_lepton_p4, jets[1].p4);

	lead_lepton_pt = lead_lepton_p4.Pt();
	lead_lepton_eta = lead_lepton_p4.Eta();
	lead_lepton_phi = lead_lepton_p4.Phi();

	sublead_lepton_pt = sublead_lepton_p4.Pt();
	sublead_lepton_eta = sublead_lepton_p4.Eta();
	sublead_lepton_phi = sublead_lepton_p4.Phi();

	// Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
	WR_mass = (lead_lepton_p4 + sublead_lepton_p4 + gJets[0].p4 + gJets[1].p4).M();
	weight = lead_lepton_weight * sublead_lepton_weight * lead_jet_weight * sublead_jet_weight * global_event_weight;
#ifdef DEBUGG
	std::cout << "weight (global_event_weight times single object weights)=\t" << weight << std::endl;
#endif

	pu_weight = fabs(global_event_weight);

	dilepton_mass = (lead_lepton_p4 + sublead_lepton_p4).M();
	zPt = (lead_lepton_p4 + sublead_lepton_p4).Pt();
	zEta = (lead_lepton_p4 + sublead_lepton_p4).Eta();
	zPhi = (lead_lepton_p4 + sublead_lepton_p4).Phi();
	/**/
	if (makeHists) {
		sel::hists("global", 4, 0, 4)->Fill(int(pair));
		if(pair == P_EE) sel::hists("global_EE", 1, 0, 1)->Fill(0);
		if(pair == P_BB) sel::hists("global_BB", 1, 0, 1)->Fill(0);
		if(pair == P_EB) sel::hists("global_EB", 1, 0, 1)->Fill(0);
		if(pair == P_GAP)  sel::hists("global_GAP", 1, 0, 1)->Fill(0);
	}
	/**/
	_isPassing = true;
	return _isPassing;

}

bool Selector::isPassingPreselect(bool makeHists)
{
	//default lepton pT cuts are 54 and 44, default jet pT cut is 30
	int l54 = 0;
	int l44 = 0;
	int j30 = 0;
	for(auto ele : electrons) {
		if (makeHists) sel::hists("preselect_ele_pt", 100, 0, 200)->Fill(ele.p4.Pt());
		if(ele.p4.Pt() > 54) l54 += 1;
		if(ele.p4.Pt() > 44) l44 += 1;
		if(!makeHists && l54 >= 1 && l44 >= 2) break;
	}
	for(auto mu : muons) {
		if (makeHists) sel::hists("preselect_mu_pt", 100, 0, 200)->Fill(mu.p4.Pt());
		if(mu.p4.Pt() > 54) l54 += 1;
		if(mu.p4.Pt() > 44) l44 += 1;
		if(!makeHists && l54 >= 1 && l44 >= 2) break;
	}
	for(auto jet : jets) {
		if (makeHists) sel::hists("preselect_jet_pt", 100, 0, 200)->Fill(jet.p4.Pt());
		if(jet.p4.Pt() > 30) j30 += 1;
		if(!makeHists && j30 >= 2) break;
	}
	if (makeHists) sel::hists("preselect_count", 1, 0, 1)->Fill(0);
	bool passes = ( l54 >= 1 && l44 >= 2 && j30 >= 2);
	/**/
	if (passes && makeHists) {
		sel::hists("preselect_count_cut", 1, 0, 1)->Fill(0);
		for(auto ele : electrons)
			sel::hists("preselect_ele_pt_cut", 100, 0, 200)->Fill(ele.p4.Pt());
		for(auto mu : muons)
			sel::hists("preselect_mu_pt_cut", 100, 0, 200)->Fill(mu.p4.Pt());
		for(auto jet : jets)
			sel::hists("preselect_jet_pt_cut", 100, 0, 200)->Fill(jet.p4.Pt());
	}
	/**/
	return passes;
}



void Selector::SetBranches(TTree* tree)
{

	tree->Branch("zPt", &zPt);
	tree->Branch("zEta", &zEta);
	tree->Branch("zPhi", &zPhi);
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
	tree->Branch("nPV", &nPV);
	tree->Branch("nPU", &nPU);
	tree->Branch("dR_leadlepton_leadjet", &dR_leadlepton_leadjet);
	tree->Branch("dR_leadlepton_subleadjet", &dR_leadlepton_subleadjet);
	tree->Branch("dR_subleadlepton_leadjet", &dR_subleadlepton_leadjet);
	tree->Branch("dR_subleadlepton_subleadjet", &dR_subleadlepton_subleadjet);
	tree->Branch("lead_lepton_r9", &lead_lepton_r9);
	tree->Branch("sublead_lepton_r9", &sublead_lepton_r9);

	tree->Branch("weight", &weight);
	tree->Branch("WR_mass", &WR_mass);
	tree->Branch("dilepton_mass", &dilepton_mass);
	tree->Branch("pu_weight", &pu_weight);
	tree->Branch("njets", &njets);
	tree->Branch("lead_jet_jec_unc", &lead_jet_jec_unc);
	tree->Branch("sublead_jet_jec_unc", &sublead_jet_jec_unc);

	tree->Branch("lead_lepton_IDSF_error", &lead_lepton_IDSF_error);
	tree->Branch("lead_lepton_IsoSF_error", &lead_lepton_IsoSF_error);
	tree->Branch("lead_lepton_RecoSF_error", &lead_lepton_RecoSF_error);
	tree->Branch("lead_lepton_HltSF_error", &lead_lepton_HltSF_error);
	tree->Branch("lead_lepton_ESmearing_error", &lead_lepton_ESmearing_error);
	tree->Branch("lead_lepton_EScaling_error", &lead_lepton_EScaling_error);
	tree->Branch("sublead_lepton_IDSF_error", &sublead_lepton_IDSF_error);
	tree->Branch("sublead_lepton_IsoSF_error", &sublead_lepton_IsoSF_error);
	tree->Branch("sublead_lepton_RecoSF_error", &sublead_lepton_RecoSF_error);
	tree->Branch("sublead_lepton_HltSF_error", &sublead_lepton_HltSF_error);
	tree->Branch("sublead_lepton_ESmearing_error", &sublead_lepton_ESmearing_error);
	tree->Branch("sublead_lepton_EScaling_error", &sublead_lepton_EScaling_error);

}

void Selector::SetBranchAddresses(TTree* tree)
{

	tree->SetBranchAddress("zPt", &zPt);
	tree->SetBranchAddress("zEta", &zEta);
	tree->SetBranchAddress("zPhi", &zPhi);
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
	tree->SetBranchAddress("nPV", &nPV);
	tree->SetBranchAddress("nPU", &nPU);

	tree->SetBranchAddress("dR_leadlepton_leadjet", &dR_leadlepton_leadjet);
	tree->SetBranchAddress("dR_leadlepton_subleadjet", &dR_leadlepton_subleadjet);
	tree->SetBranchAddress("dR_subleadlepton_leadjet", &dR_subleadlepton_leadjet);
	tree->SetBranchAddress("dR_subleadlepton_subleadjet", &dR_subleadlepton_subleadjet);

	tree->SetBranchAddress("lead_lepton_r9", &lead_lepton_r9);
	tree->SetBranchAddress("sublead_lepton_r9", &sublead_lepton_r9);

	tree->SetBranchAddress("weight", &weight);
	tree->SetBranchAddress("WR_mass", &WR_mass);
	tree->SetBranchAddress("dilepton_mass", &dilepton_mass);
	tree->SetBranchAddress("pu_weight", &pu_weight);
	tree->SetBranchAddress("njets", &njets);
	tree->SetBranchAddress("lead_jet_jec_unc", &lead_jet_jec_unc);
	tree->SetBranchAddress("sublead_jet_jec_unc", &sublead_jet_jec_unc);

	tree->SetBranchAddress("lead_lepton_IDSF_error", &lead_lepton_IDSF_error);
	tree->SetBranchAddress("lead_lepton_IsoSF_error", &lead_lepton_IsoSF_error);
	tree->SetBranchAddress("lead_lepton_RecoSF_error", &lead_lepton_RecoSF_error);
	tree->SetBranchAddress("lead_lepton_HltSF_error", &lead_lepton_HltSF_error);
	tree->SetBranchAddress("lead_lepton_ESmearing_error", &lead_lepton_ESmearing_error);
	tree->SetBranchAddress("lead_lepton_EScaling_error", &lead_lepton_EScaling_error);
	tree->SetBranchAddress("sublead_lepton_IDSF_error", &sublead_lepton_IDSF_error);
	tree->SetBranchAddress("sublead_lepton_IsoSF_error", &sublead_lepton_IsoSF_error);
	tree->SetBranchAddress("sublead_lepton_RecoSF_error", &sublead_lepton_RecoSF_error);
	tree->SetBranchAddress("sublead_lepton_HltSF_error", &sublead_lepton_HltSF_error);
	tree->SetBranchAddress("sublead_lepton_ESmearing_error", &sublead_lepton_ESmearing_error);
	tree->SetBranchAddress("sublead_lepton_EScaling_error", &sublead_lepton_EScaling_error);

}

