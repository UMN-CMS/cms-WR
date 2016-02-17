#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class TunePMuonProducer : public edm::EDProducer
{
public:
	TunePMuonProducer(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);

private:
	edm::EDGetToken srcToken_;
	edm::EDGetToken primaryVertexToken_;

};

TunePMuonProducer::TunePMuonProducer(const edm::ParameterSet& cfg)
{

	produces <pat::MuonCollection>();
	srcToken_ = mayConsume<edm::View<pat::Muon> >(cfg.getParameter<edm::InputTag>("src"));
	primaryVertexToken_ = mayConsume<edm::View<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));

}

void TunePMuonProducer::produce(edm::Event& event, const edm::EventSetup&)
{
	edm::Handle<edm::View<pat::Muon> > muons;
	event.getByToken(srcToken_, muons);
	edm::Handle<edm::View<reco::Vertex> > primary_vertex;
	event.getByToken(primaryVertexToken_, primary_vertex);

	std::auto_ptr<pat::MuonCollection> mus(new pat::MuonCollection);
	reco::Candidate::PolarLorentzVector tmp_muon;

	for(auto mu : *muons) {
		if(mu.tunePMuonBestTrack().isAvailable())
			tmp_muon.SetCoordinates(mu.tunePMuonBestTrack()->pt(), mu.tunePMuonBestTrack()->eta(), mu.tunePMuonBestTrack()->phi(), 0.1057);
		else
			tmp_muon.SetCoordinates(0.0, 0.0, 0.0, 0.0);
		if(mu.isHighPtMuon(primary_vertex->at(0)))
			mu.addUserInt("highPtID", 1);
		else
			mu.addUserInt("highPtID", 0);

		mu.setP4(tmp_muon);
		mus->push_back(mu);
	}
	event.put(mus);
}

DEFINE_FWK_MODULE(TunePMuonProducer);
