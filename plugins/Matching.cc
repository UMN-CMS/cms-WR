#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"


class Matching : public edm::EDAnalyzer
{
public:
	explicit Matching(const edm::ParameterSet&);
	void analyze(const edm::Event&, const edm::EventSetup&);

private:
	const edm::InputTag gen_src;
	const edm::InputTag vertex_src;
	const edm::InputTag gen_jet_src;
	const edm::InputTag electron_src;
	const edm::InputTag muon_src;
	const edm::InputTag jet_src;

	TH1F* mu_pre0;
	TH1F* gen_mu_pre0;
	TH1F* mu_pre1;
	TH1F* gen_mu_pre1;
	TH1F* gen_mu_pre2;
	TH1F* Mlljj;
	TH1F* Mlljj_MC;
	TH1F* gMlljj;
	TH1F* gMlljj2;
	TH1F* gjet_pt;
	TH1F* jetindex0;
	TH1F* jetindex1;
	TH2F* jetindex;
	TH1F* rejected_dR1;
	TH1F* rejected_dR2;

};

Matching::Matching(const edm::ParameterSet& cfg)
	: gen_src(cfg.getParameter<edm::InputTag>("gen_src")),
	  vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
	  gen_jet_src(cfg.getParameter<edm::InputTag>("gen_jet_src")),
	  electron_src(cfg.getParameter<edm::InputTag>("electron_src")),
	  muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
	  jet_src(cfg.getParameter<edm::InputTag>("jet_src"))
{
	edm::Service<TFileService> fs;
	mu_pre0 = fs->make<TH1F>("mu_pre0", "", 5, 0, 5);
	gen_mu_pre0 = fs->make<TH1F>("gen_mu_pre0", "", 5, 0, 5);
	mu_pre1 = fs->make<TH1F>("mu_pre1", "", 4, 0, 4);
	gen_mu_pre1 = fs->make<TH1F>("gen_mu_pre1", "", 7, 0, 7);
	gen_mu_pre2 = fs->make<TH1F>("gen_mu_pre2", "", 10, 0, 10);
	Mlljj = fs->make<TH1F>("Mlljj", "", 100, 0, 4500);
	Mlljj_MC = fs->make<TH1F>("Mlljj_MC", "", 100, 0, 4500);
	gMlljj = fs->make<TH1F>("gMlljj", "", 100, 0, 4500);
	gMlljj2 = fs->make<TH1F>("gMlljj2", "", 100, 0, 4500);
	gjet_pt = fs->make<TH1F>("gjet_pt", "", 100, 0, 1500);
	jetindex0 = fs->make<TH1F>("jetindex0", "", 5, 0, 5);
	jetindex1 = fs->make<TH1F>("jetindex1", "", 5, 0, 5);
	jetindex = fs->make<TH2F>("jetindex", "", 5, 0, 5, 5, 0, 5);
	rejected_dR1 = fs->make<TH1F>("rejected_dR1", "", 100, 0, 10);
	rejected_dR2 = fs->make<TH1F>("rejected_dR2", "", 100, 0, 10);
}


void Matching::analyze(const edm::Event& event, const edm::EventSetup& setup)
{
	using namespace std;
	using namespace edm;
	using namespace reco;

	edm::Handle<pat::ElectronCollection> electrons;
	event.getByLabel(electron_src, electrons);
	edm::Handle<pat::MuonCollection> muons;
	event.getByLabel(muon_src, muons);
	edm::Handle<reco::VertexCollection> primary_vertex;
	event.getByLabel(vertex_src, primary_vertex);
	edm::Handle<pat::JetCollection> jets;
	event.getByLabel(jet_src, jets);
	edm::Handle<reco::GenParticleCollection> gen_particles;
	event.getByLabel(gen_src, gen_particles);
	edm::Handle<reco::GenJetCollection> gen_jets;
	event.getByLabel(gen_jet_src, gen_jets);

	std::vector<pat::Electron> eles;
	std::vector<pat::Muon> mus(2);
	std::vector<reco::GenParticle> WR_daughters;
	std::vector<reco::GenParticle> gen_mus(2);
	std::vector<reco::GenParticle> gen_js(2);
	std::vector<reco::GenParticle> gmus_W;
	std::vector<pat::Jet> js(2);
	std::vector<reco::GenJet> gjets(2);

	vector<int> jet_index = {99, 99};

	//cout<<"event"<<endl;
	mu_pre1->Fill(0);
	for(auto g : *gen_particles) {
		if(abs(g.pdgId()) == 9900012 || abs(g.pdgId()) == 9900016)
			mu_pre1->Fill(1);
		if(abs(g.pdgId()) < 15 && (abs(g.mother()->pdgId()) == 9900024 || abs(g.mother()->pdgId()) == 9900014 ) ) {
			//if(g.pt() > 40 && fabs(g.eta()) < 2.4)
			WR_daughters.push_back(g);
		}
	}

	vector<pat::Jet> selected_jets;
	vector<pat::Jet> rejected_jets;
	for(auto j : *jets) {
		if((j.neutralHadronEnergyFraction() < 0.99 && j.neutralEmEnergyFraction() < 0.99 && (j.chargedMultiplicity() + j.neutralMultiplicity()) > 1 && j.muonEnergyFraction() < 0.8) && ((abs(j.eta()) <= 2.4 && j.chargedHadronEnergyFraction() > 0 && j.chargedMultiplicity() > 0 && j.chargedEmEnergyFraction() < 0.99) || abs(j.eta()) > 2.4))
			selected_jets.push_back(j);
		else
			rejected_jets.push_back(j);
	}

	if(WR_daughters.size() == 4) {
		gMlljj->Fill((WR_daughters[0].p4() + WR_daughters[1].p4() + WR_daughters[2].p4() + WR_daughters[3].p4()).M());
		mu_pre1->Fill(2);
		if(abs(WR_daughters[0].pdgId()) != 13 || abs(WR_daughters[1].pdgId()) != 13 )
			cout << "BREAK" << endl;

		double tmp_dR_0 = 0.4;
		double tmp_dR_1 = 0.4;
		for(auto m : *muons) {
			if(deltaR(WR_daughters[0], m) < tmp_dR_0) {
				mus[0] = m;
				tmp_dR_0 = deltaR(WR_daughters[0], m);
				continue;
			}
			if(deltaR(WR_daughters[1], m) < tmp_dR_1) {
				mus[1] = m;
				tmp_dR_1 = deltaR(WR_daughters[1], m);
			}
		}
		tmp_dR_0 = 0.4;
		tmp_dR_1 = 0.4;
		for(auto gj : *gen_jets) {
			gjet_pt->Fill(gj.pt());
			if(gj.pt() < 40 || fabs(gj.eta()) > 2.4)
				continue;
			if(deltaR(WR_daughters[2], gj) < tmp_dR_0) {
				gjets[0] = gj;
				tmp_dR_0 = deltaR(WR_daughters[2], gj);
				continue;
			}
			if(deltaR(WR_daughters[3], gj) < tmp_dR_1) {
				gjets[1] = gj;
				tmp_dR_1 = deltaR(WR_daughters[3], gj);
			}
		}

		if(WR_daughters[0].pt() != 0 && WR_daughters[1].pt() != 0 && gjets[0].pt() != 0 && gjets[1].pt() != 0)
			gMlljj2->Fill((WR_daughters[0].p4() + WR_daughters[1].p4() + gjets[0].p4() + gjets[1].p4()).M());

		tmp_dR_0 = 0.4;
		tmp_dR_1 = 0.4;
		int i = -1;
		for(auto j : selected_jets) {
			i++;
			if(j.pt() > 40 && fabs(j.eta()) < 2.4) {
				if(deltaR(WR_daughters[2], j) < tmp_dR_0) {
					js[0] = j;
					tmp_dR_0 = deltaR(WR_daughters[2], j);
					jet_index[0] = i;
					continue;
				}
				if(deltaR(WR_daughters[3], j) < tmp_dR_1) {
					js[1] = j;
					tmp_dR_1 = deltaR(WR_daughters[3], j);
					jet_index[1] = i;
				}
			}
		}

		//cout<<"Index="<<jet_index[0] <<" "<<	jet_index[1]<<endl;
		jetindex0->Fill(jet_index[0]);
		jetindex1->Fill(jet_index[1]);
		//jetindex->Fill(jet_index[0],jet_index[1]);

		if(mus[0].pt() == 0)
			mu_pre0->Fill(0);
		else if(mus[1].pt() == 0)
			mu_pre0->Fill(1);
		else if(js[0].pt() == 0)
			mu_pre0->Fill(2);
		else if(js[1].pt() == 0)
			mu_pre0->Fill(3);
		else
			mu_pre0->Fill(4);

		if(mus[0].pt() != 0 && mus[1].pt() != 0 && js[0].pt() != 0 && js[1].pt() != 0)
			Mlljj_MC->Fill((mus[0].p4() + mus[1].p4() + js[0].p4() + js[1].p4()).M());

		/*cout<<"Event"<<endl;
		cout<<WR_daughters[0].pt()<<" "<<WR_daughters[1].pt()<<" "<<WR_daughters[2].pt()<<" "<<WR_daughters[3].pt()<<endl;
		cout<<"all jets"<<endl;
		for(auto i: *jets)
		  cout<<i.pt()<<" "<<deltaR(WR_daughters[0],i)<<" "<<deltaR(WR_daughters[1],i)<<endl;
		cout<<"selected jets"<<endl;
		for(auto i: selected_jets)
		  cout<<i.pt()<<" "<<deltaR(WR_daughters[2],i)<<" "<<deltaR(WR_daughters[3],i)<<endl;
		*/

		for(auto rj : rejected_jets)
			if(rj.pt() > 40 && fabs(rj.eta()) < 2.4) {
				rejected_dR1->Fill(deltaR(WR_daughters[0], rj));
				rejected_dR1->Fill(deltaR(WR_daughters[1], rj));
				rejected_dR2->Fill(deltaR(WR_daughters[3], rj));
				rejected_dR2->Fill(deltaR(WR_daughters[4], rj));
			}
	}

	if(muons->size() == 1 && jets->size() == 0)
		gen_mu_pre1->Fill(0);
	else if(muons->size() > 1  && jets->size() == 0)
		gen_mu_pre1->Fill(1);
	else if(muons->size() > 1  && jets->size() == 1)
		gen_mu_pre1->Fill(2);
	else if(muons->size() > 1  && jets->size() > 1)
		gen_mu_pre1->Fill(3);
	else if(muons->size() == 1 && jets->size() == 1)
		gen_mu_pre1->Fill(4);
	else if(muons->size() == 1 && jets->size() > 1)
		gen_mu_pre1->Fill(5);
	else
		gen_mu_pre1->Fill(6);

	vector<pat::Muon> goodmus(2);
	vector<pat::Jet> goodjets(2);

	if(muons->size() > 0) {
		auto mu_0 = muons->at(0);
		auto iso_mu_0 = mu_0.pfIsolationR04();
		double I_0 = 9.9;
		I_0 = (iso_mu_0.sumChargedHadronPt + max(0., iso_mu_0.sumNeutralHadronEt + iso_mu_0.sumPhotonEt - 0.5 * iso_mu_0.sumPUPt)) / mu_0.pt();
		if(mu_0.pt() > 40 && fabs(mu_0.eta()) < 2.4 && mu_0.isHighPtMuon(primary_vertex->at(0)) && I_0 < 0.2) {
			goodmus[0] = mu_0;
			double tmp_dR = 0.4;
			for(auto g : *gen_particles)
				if(deltaR(g, mu_0) < tmp_dR && g.fromHardProcessBeforeFSR() && abs(g.pdgId()) == 13) {
					gen_mus[0] = g;
					tmp_dR = deltaR(g, mu_0);
				}
		}
	}
	if(muons->size() > 1) {
		auto mu_1 = muons->at(1);
		auto iso_mu_1 = mu_1.pfIsolationR04();
		double I_1 = 9.9;
		I_1 = (iso_mu_1.sumChargedHadronPt + max(0., iso_mu_1.sumNeutralHadronEt + iso_mu_1.sumPhotonEt - 0.5 * iso_mu_1.sumPUPt)) / mu_1.pt();
		if(mu_1.pt() > 40 && fabs(mu_1.eta()) < 2.4 && mu_1.isHighPtMuon(primary_vertex->at(0)) && I_1 < 0.2) {
			goodmus[1] = mu_1;
			double tmp_dR = 0.4;
			for(auto g : *gen_particles)
				if(deltaR(g, mu_1) < tmp_dR && g.fromHardProcessBeforeFSR() && abs(g.pdgId()) == 13) {
					gen_mus[1] = g;
					tmp_dR = deltaR(g, mu_1);
				}
		}
	}

	if(selected_jets.size() > 0) {
		auto j_0 = selected_jets.at(0);
		if(j_0.pt() > 40 && fabs(j_0.eta()) < 2.4) {
			double tmp_dR = 0.4;
			goodjets[0] = j_0;
			for(auto g : *gen_particles)
				if(deltaR(g, j_0) < tmp_dR && g.fromHardProcessBeforeFSR() && abs(g.pdgId()) < 7) {
					gen_js[0] = g;
					tmp_dR = deltaR(g, j_0);
				}
			//if(deltaR(gjets[0],j_0) < tmp_dR){
			//gen_js[0] = gjets[0];
			//tmp_dR = deltaR(gjets[0],j_0);
			//}
			//else if(deltaR(gjets[1],j_0) < tmp_dR){
			//gen_js[0] = gjets[1];
			//tmp_dR = deltaR(gjets[1],j_0);
		}
	}

	if(selected_jets.size() > 1) {
		auto j_1 = selected_jets.at(1);
		if(j_1.pt() > 40 && fabs(j_1.eta()) < 2.4) {
			double tmp_dR = 0.4;
			goodjets[1] = j_1;
			for(auto g : *gen_particles)
				if(deltaR(g, j_1) < tmp_dR && g.fromHardProcessBeforeFSR() && abs(g.pdgId()) < 7) {
					gen_js[1] = g;
					tmp_dR = deltaR(g, j_1);
				}
			//if(deltaR(gjets[1],j_1) < tmp_dR){
			//gen_js[1] = gjets[1];
			//tmp_dR = deltaR(gjets[1],j_1);
			//}
			//else if(deltaR(gjets[0],j_1) < tmp_dR){
			//gen_js[1] = gjets[0];
			//tmp_dR = deltaR(gjets[0],j_1);
			//}
		}
	}

	if(muons->size() > 1)
		gen_mu_pre2->Fill(0);
	if(goodmus[0].pt() != 0 && goodmus[1].pt() != 0)
		gen_mu_pre2->Fill(1);
	if(selected_jets.size() > 1)
		gen_mu_pre2->Fill(2);
	if(goodjets[0].pt() != 0 && goodjets[1].pt() != 0)
		gen_mu_pre2->Fill(3);
	if(goodmus[0].pt() != 0 && goodmus[1].pt() != 0 && goodjets[0].pt() != 0 && goodjets[1].pt() != 0)
		gen_mu_pre2->Fill(4);
	if(muons->size() > 1 && muons->at(1).pt() > 40 && fabs(muons->at(0).eta()) < 2.4 && fabs(muons->at(1).eta()) < 2.4 )
		gen_mu_pre2->Fill(5);
	if(muons->size() > 1 && muons->at(1).pt() > 40 && fabs(muons->at(0).eta()) < 2.4 && fabs(muons->at(1).eta()) < 2.4 && muons->at(0).isHighPtMuon(primary_vertex->at(0)) && muons->at(1).isHighPtMuon(primary_vertex->at(0)))
		gen_mu_pre2->Fill(6);

	double I_0 = 9.9;
	double I_1 = 9.9;

	if(muons->size() > 1) {
		auto mu_0 = muons->at(0);
		auto iso_mu_0 = mu_0.pfIsolationR04();
		I_0 = (iso_mu_0.sumChargedHadronPt + max(0., iso_mu_0.sumNeutralHadronEt + iso_mu_0.sumPhotonEt - 0.5 * iso_mu_0.sumPUPt)) / mu_0.pt();
		auto mu_1 = muons->at(1);
		auto iso_mu_1 = mu_1.pfIsolationR04();
		I_1 = (iso_mu_1.sumChargedHadronPt + max(0., iso_mu_1.sumNeutralHadronEt + iso_mu_1.sumPhotonEt - 0.5 * iso_mu_1.sumPUPt)) / mu_1.pt();
	}

	if(muons->size() > 1 && muons->at(1).pt() > 40 && fabs(muons->at(0).eta()) < 2.4 && fabs(muons->at(1).eta()) < 2.4 && muons->at(0).isHighPtMuon(primary_vertex->at(0)) && muons->at(1).isHighPtMuon(primary_vertex->at(0)) && I_0 < 0.2 && I_1 < 0.2)
		gen_mu_pre2->Fill(7);


	if(muons->size() > 1 && jets->size() > 1) {
		if(gen_mus[0].pt() == 0)
			gen_mu_pre0->Fill(0);
		else if(gen_mus[1].pt() == 0)
			gen_mu_pre0->Fill(1);
		else if(gen_js[0].pt() == 0)
			gen_mu_pre0->Fill(2);
		else if(gen_js[1].pt() == 0)
			gen_mu_pre0->Fill(3);
		else
			gen_mu_pre0->Fill(4);
	}

	if(goodmus[0].pt() != 0 && goodmus[1].pt() != 0 && goodjets[0].pt() != 0 && goodjets[1].pt() != 0) {
		Mlljj->Fill((goodmus[0].p4() + goodmus[1].p4() + goodjets[0].p4() + goodjets[1].p4()).M());
		jetindex->Fill(jet_index[0], jet_index[1]);
	}
	/*
	cout<<"event"<<endl;
	cout<<goodjets[0].pt()-js[0].pt()<<" "<<goodjets[0].pt()-js[1].pt()<<endl;
	cout<<goodjets[1].pt()-js[0].pt()<<" "<<goodjets[1].pt()-js[1].pt()<<endl;

	if(mus1)
	  cout<<"Mus1"<<endl;
	if(mus2)
	  cout<<"Mus2"<<endl;

	for(auto g: *gen_particles){
	  if(g.fromHardProcessBeforeFSR())
	    if(abs(g.pdgId()) == 13)
	cout<<g.pdgId()<<endl;
	  if(abs(g.pdgId()) == 9900012 || abs(g.pdgId()) == 9900016)
	    cout<<"FLAVOR FLAV"<<endl;
	}
	*/
}

DEFINE_FWK_MODULE(Matching);
