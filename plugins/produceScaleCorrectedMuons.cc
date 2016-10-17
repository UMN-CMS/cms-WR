#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/PatAlgos/interface/ObjectModifier.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "ExoAnalysis/cmsWR/interface/rochcor2016.h"
//#include "ExoAnalysis/cmsWR/interface/muresolution_run2.h"
#include <string>
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"
#include "Math/PxPyPzM4D.h"
#include <vector>
using namespace std;

class produceScaleCorrectedMuons : public edm::EDProducer
{
public:
	produceScaleCorrectedMuons(const edm::ParameterSet&);
	virtual void produce(edm::Event&, const edm::EventSetup&);
	TLorentzVector Mu_Original;
	int charge;
	float qter = 1.0;
//        std::unique_ptr<pat::ObjectModifier<pat::Muon> > muonModifier;
private:
	//const edm::InputTag src;	///<input particle objects
	edm::EDGetToken srcToken_;
	std::string outputCollName;     ///<label name of collection made by this producer
};

produceScaleCorrectedMuons::produceScaleCorrectedMuons(const edm::ParameterSet& cfg)
	: srcToken_ ( consumes<edm::View<pat::Muon> >(cfg.getParameter<edm::InputTag>("src"))),
	  outputCollName(cfg.getParameter<std::string>("OutputCollectionName"))
{
	produces <pat::MuonCollection>(outputCollName);
}

void produceScaleCorrectedMuons::produce(edm::Event& event, const edm::EventSetup& setup)
{
	edm::Handle<edm::View<pat::Muon> > muons;
	event.getByToken(srcToken_, muons);
	std::auto_ptr<pat::MuonCollection> mus(new pat::MuonCollection);
	rochcor2016 *rmcor = new rochcor2016();
	for(auto mu : *muons) {
		charge = mu.charge();
		Mu_Original.SetPtEtaPhiE(mu.pt(), mu.eta(), mu.phi(), mu.energy());
		if(!event.isRealData())  rmcor->momcor_mc(Mu_Original, charge, 0, qter);
		else      rmcor->momcor_data(Mu_Original, charge, 0, qter);
		reco::Candidate::PolarLorentzVector p4(Mu_Original.Pt(), Mu_Original.Eta(), Mu_Original.Phi(), 0.1057);
		mu.setP4(p4);
		mus->push_back(mu);
	}
	event.put(mus, outputCollName);
	delete rmcor;
}

DEFINE_FWK_MODULE(produceScaleCorrectedMuons);
