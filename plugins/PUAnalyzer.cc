#include "TTree.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class PUAnalyzer : public edm::EDAnalyzer
{
public:
	explicit PUAnalyzer(const edm::ParameterSet&);

private:
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

	edm::EDGetToken pileUpInfoToken_;

	TH1F * h_pileup;

};

PUAnalyzer::PUAnalyzer(const edm::ParameterSet& cfg):

	pileUpInfoToken_ ( consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo")))

{
	edm::Service<TFileService> fs;

	h_pileup = fs->make<TH1F>("pileup", "pileup", 1000, 0, 100);
}

void PUAnalyzer::analyze(const edm::Event& event, const edm::EventSetup&)
{

	edm::Handle<edm::View<PileupSummaryInfo> > PU_Info;

	if(!event.isRealData()) {
		event.getByToken(pileUpInfoToken_, PU_Info);
		for(auto p : *PU_Info) {
			int BX = p.getBunchCrossing();
			if(BX == 0)
				h_pileup->Fill(p.getTrueNumInteractions());
		}
	}



}

DEFINE_FWK_MODULE(PUAnalyzer);
