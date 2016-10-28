#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

class Weight : public edm::EDAnalyzer
{
public:
	explicit Weight(const edm::ParameterSet&);

private:
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

	edm::EDGetToken evinfoToken_;

	struct tree_t {
		unsigned run;
		unsigned lumi;
		unsigned event;

		float weight;

		tree_t()
		{
			clear();
		}

		void clear()
		{
			run = lumi = event = 0;
			weight = 0;
		}
	};

	TTree* tree;
	tree_t nt;
};

Weight::Weight(const edm::ParameterSet& cfg) :
  evinfoToken_ ( consumes<GenEventInfoProduct>(edm::InputTag("generator")))
{
	edm::Service<TFileService> fs;
	tree = fs->make<TTree>("t", "");
	tree->Branch("run", &nt.run, "run/i");
	tree->Branch("lumi", &nt.lumi, "lumi/i");
	tree->Branch("event", &nt.event, "event/i");

	tree->Branch("weight", &nt.weight, "weight/F");
}

void Weight::analyze(const edm::Event& event, const edm::EventSetup&)
{
	nt.clear();
	nt.run = event.id().run();
	nt.lumi = event.luminosityBlock();
	nt.event = event.id().event();

	edm::Handle<GenEventInfoProduct> evinfo;
	event.getByToken(evinfoToken_, evinfo);
	nt.weight = evinfo->weight();

	tree->Fill();
}

DEFINE_FWK_MODULE(Weight);
