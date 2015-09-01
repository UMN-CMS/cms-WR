#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <string>
#include <vector>

class produceIsHighPtMuons : public edm::EDProducer
{
public:
	produceIsHighPtMuons(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);

private:
	const edm::InputTag src;	///<input particle objects
	edm::Handle<std::vector<reco::Vertex> > vertices;
	edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;	///< needed for isHighPtMuon()
	std::string outputCollName;	///<label name of collection made by this producer

};

produceIsHighPtMuons::produceIsHighPtMuons(const edm::ParameterSet& cfg)
	: src(cfg.getParameter<edm::InputTag>("src")),
	outputCollName(cfg.getParameter<std::string>("outputCollectionName"))
{
 	verticesToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));

	produces <pat::MuonCollection>(outputCollName);

}

void produceIsHighPtMuons::produce(edm::Event& event, const edm::EventSetup&)
{
	event.getByToken(verticesToken, vertices);
	edm::Handle<pat::MuonCollection> muons;
	event.getByLabel(src, muons);

	std::auto_ptr<pat::MuonCollection> mus(new pat::MuonCollection);
	//reco::Candidate::PolarLorentzVector tmp_muon;

	for(auto mu : *muons) {
		/*
		if(mu.tunePMuonBestTrack().isAvailable())
			tmp_muon.SetCoordinates(mu.tunePMuonBestTrack()->pt(), mu.tunePMuonBestTrack()->eta(), mu.tunePMuonBestTrack()->phi(), 0.1057);
		else
			tmp_muon.SetCoordinates(0.0, 0.0, 0.0, 0.0);
		mu.setP4(tmp_muon);
		*/
		if( mu.isHighPtMuon(vertices->at(0)) ) mus->push_back(mu);
		
	}
	event.put(mus, outputCollName);
}

DEFINE_FWK_MODULE(produceIsHighPtMuons);
