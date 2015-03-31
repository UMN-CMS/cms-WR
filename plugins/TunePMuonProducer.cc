#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class TunePMuonProducer : public edm::EDProducer
{
public:
	TunePMuonProducer(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);

private:
	const edm::InputTag src;

};

TunePMuonProducer::TunePMuonProducer(const edm::ParameterSet& cfg)
	: src(cfg.getParameter<edm::InputTag>("src"))
{

	produces <pat::MuonCollection>();

}

void TunePMuonProducer::produce(edm::Event& event, const edm::EventSetup&)
{
	edm::Handle<pat::MuonCollection> muons;
	event.getByLabel(src, muons);

	std::auto_ptr<pat::MuonCollection> mus(new pat::MuonCollection);
	reco::Candidate::PolarLorentzVector tmp_muon;

	for(auto mu : *muons) {
		if(mu.tunePMuonBestTrack().isAvailable())
			tmp_muon.SetCoordinates(mu.tunePMuonBestTrack()->pt(), mu.tunePMuonBestTrack()->eta(), mu.tunePMuonBestTrack()->phi(), 0.1057);
		else
			tmp_muon.SetCoordinates(0.0, 0.0, 0.0, 0.0);
		mu.setP4(tmp_muon);
		mus->push_back(mu);
	}

	event.put(mus);
}

DEFINE_FWK_MODULE(TunePMuonProducer);
