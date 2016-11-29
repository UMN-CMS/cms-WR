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
#include "../interface/miniTreeEvent.h"



typedef double JECUnc_t;
typedef edm::ValueMap<JECUnc_t> JECUnc_Map;

class miniTTree : public edm::EDAnalyzer
{
public:
	explicit miniTTree(const edm::ParameterSet&);

private:
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

	edm::EDGetToken electronsMiniAODToken_;
	edm::EDGetToken muonsMiniAODToken_;
	edm::EDGetToken jetsMiniAODToken_;
	edm::EDGetToken pileUpInfoToken_;
	edm::EDGetToken pileUpReweightToken_;
	edm::EDGetToken primaryVertexToken_;
	edm::EDGetToken evinfoToken_;

	edm::EDGetToken  jec_unc_src;
	edm::EDGetToken  jetResolution_src;
	edm::EDGetToken  JERsf_src;
	edm::EDGetToken  JERsf_up_src;
	edm::EDGetToken  JERsf_down_src;
	edm::EDGetToken  genJetPt_src;
	edm::EDGetToken  genJetMatch_src;

	edm::EDGetToken  ele_scale_error_src;
	edm::EDGetToken  ele_smearing_sigma_src;
	edm::EDGetToken  ele_smearing_sigma_phi_up_src;
	edm::EDGetToken  ele_smearing_sigma_phi_down_src;
	edm::EDGetToken  ele_smearing_sigma_rho_up_src;
	edm::EDGetToken  ele_smearing_sigma_rho_down_src;

	edm::EDGetToken  muon_IDSF_central_src;
	edm::EDGetToken  muon_IsoSF_central_src;
	edm::EDGetToken  muon_IDSF_error_src;
	edm::EDGetToken  muon_IsoSF_error_src;

	edm::EDGetToken datasetNameToken_;
	TTree* tree;
	miniTreeEvent myEvent;

};

miniTTree::miniTTree(const edm::ParameterSet& cfg):

	electronsMiniAODToken_   ( consumes<edm::View<pat::Electron> >(cfg.getParameter<edm::InputTag>("electrons_src"))),
	muonsMiniAODToken_ ( consumes<edm::View<pat::Muon> >(cfg.getParameter<edm::InputTag>("muons_src"))),
	jetsMiniAODToken_ ( consumes<edm::View<pat::Jet> >(cfg.getParameter<edm::InputTag>("jets_src"))),
	pileUpInfoToken_ ( consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"))),
	pileUpReweightToken_ ( consumes<float >(cfg.getParameter<edm::InputTag>("PUWeights_src"))),
	primaryVertexToken_ ( consumes<edm::View<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
	evinfoToken_ ( consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	jec_unc_src ( consumes<JECUnc_Map >(cfg.getParameter<edm::InputTag>("jec_unc_src"))),
	jetResolution_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("jetResolution_src"))),
	JERsf_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("JERsf_src"))),
	JERsf_up_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("JERsf_up_src"))),
	JERsf_down_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("JERsf_down_src"))),
	genJetPt_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("genJetPt_src"))),
	genJetMatch_src ( consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("genJetMatch_src"))),
	ele_scale_error_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_scale_error_src"))),
	ele_smearing_sigma_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_smearing_sigma_src"))),
	ele_smearing_sigma_phi_up_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_smearing_sigma_phi_up_src"))),
	ele_smearing_sigma_phi_down_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_smearing_sigma_phi_down_src"))),
	ele_smearing_sigma_rho_up_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_smearing_sigma_rho_up_src"))),
	ele_smearing_sigma_rho_down_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("ele_smearing_sigma_rho_down_src"))),
	muon_IDSF_central_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("muon_IDSF_central_src"))),
	muon_IsoSF_central_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("muon_IsoSF_central_src"))),
	muon_IDSF_error_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("muon_IDSF_error_src"))),
	muon_IsoSF_error_src ( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("muon_IsoSF_error_src"))),
	datasetNameToken_ ( consumes<std::string>(cfg.getParameter<edm::InputTag>("datasetName")))
{
	edm::Service<TFileService> fs;
	tree = fs->make<TTree>("t", "");

	myEvent.SetBranches(tree);

}

void miniTTree::analyze(const edm::Event& event, const edm::EventSetup&)
{
	myEvent.clear();
	myEvent.run = event.id().run();
	myEvent.lumi = event.luminosityBlock();
	myEvent.event = event.id().event();

	edm::Handle<edm::View<pat::Electron> > electrons;
	event.getByToken(electronsMiniAODToken_, electrons);
	edm::Handle<edm::View<pat::Muon> > muons;
	event.getByToken(muonsMiniAODToken_, muons);
	edm::Handle<edm::View<pat::Jet> > jets;
	event.getByToken(jetsMiniAODToken_, jets);

	edm::Handle<JECUnc_Map > jec_unc;
	event.getByToken(jec_unc_src, jec_unc);
	edm::Handle< edm::ValueMap<float> > jetResolution;
	event.getByToken(jetResolution_src, jetResolution);
	edm::Handle< edm::ValueMap<float> > JERsf;
	event.getByToken(JERsf_src, JERsf);
	edm::Handle< edm::ValueMap<float> > JERsf_up;
	event.getByToken(JERsf_up_src, JERsf_up);
	edm::Handle< edm::ValueMap<float> > JERsf_down;
	event.getByToken(JERsf_down_src, JERsf_down);
	edm::Handle< edm::ValueMap<float> > genjetPt;
	event.getByToken(genJetPt_src, genjetPt);
	edm::Handle< edm::ValueMap<bool> > genjetMatch;
	event.getByToken(genJetMatch_src, genjetMatch);

	edm::Handle< edm::ValueMap<float> > ele_scale_error;
	event.getByToken(ele_scale_error_src, ele_scale_error);
	edm::Handle< edm::ValueMap<float> > ele_smearing_sigma;
	event.getByToken(ele_smearing_sigma_src, ele_smearing_sigma);
	edm::Handle< edm::ValueMap<float> > ele_smearing_sigma_phi_up;
	event.getByToken(ele_smearing_sigma_phi_up_src, ele_smearing_sigma_phi_up);
	edm::Handle< edm::ValueMap<float> > ele_smearing_sigma_phi_down;
	event.getByToken(ele_smearing_sigma_phi_up_src, ele_smearing_sigma_phi_down);
	edm::Handle< edm::ValueMap<float> > ele_smearing_sigma_rho_up;
	event.getByToken(ele_smearing_sigma_rho_up_src, ele_smearing_sigma_rho_up);
	edm::Handle< edm::ValueMap<float> > ele_smearing_sigma_rho_down;
	event.getByToken(ele_smearing_sigma_rho_up_src, ele_smearing_sigma_rho_down);


	edm::Handle< edm::ValueMap<float> > muon_IDSF;
	event.getByToken(muon_IDSF_central_src, muon_IDSF);
	edm::Handle< edm::ValueMap<float> > muon_IsoSF;
	event.getByToken(muon_IsoSF_central_src, muon_IsoSF);
	edm::Handle< edm::ValueMap<float> > muon_IDSF_error;
	event.getByToken(muon_IDSF_error_src, muon_IDSF_error);
	edm::Handle< edm::ValueMap<float> > muon_IsoSF_error;
	event.getByToken(muon_IsoSF_error_src, muon_IsoSF_error);

	edm::Handle<GenEventInfoProduct> evinfo;		
	edm::Handle<edm::View<PileupSummaryInfo> > PU_Info;		
	edm::Handle<float > PU_Weights;

	edm::Handle<edm::View<reco::Vertex> > primary_vertex;
	event.getByToken(primaryVertexToken_, primary_vertex);

	edm::Handle<std::string> datasetName;
	event.getByToken(datasetNameToken_, datasetName);

	sprintf(myEvent.datasetName, "%s", datasetName->c_str());

	if(primary_vertex->size() > 0) {
		for(auto pv : *primary_vertex)
			myEvent.nPV++;
	}

	if(!event.isRealData()) {		
	  event.getByToken(evinfoToken_, evinfo);		
	  myEvent.weight = evinfo->weight();		
	  event.getByToken(pileUpInfoToken_, PU_Info);		
	  for(auto p : *PU_Info) {		
	    int BX = p.getBunchCrossing();		
	    if(BX == 0)		
	      myEvent.nPU = p.getTrueNumInteractions();		
	  }		
	}

	for (size_t i = 0; i < electrons->size(); ++i) {
		const auto ele = electrons->ptrAt(i);
		TLorentzVector p4;
		p4.SetPtEtaPhiM(ele->pt(), ele->eta(), ele->phi(), ele->mass());
		myEvent.electrons_p4->push_back(p4);
		myEvent.electron_scale_error->push_back((*ele_scale_error)[ele]);
		myEvent.electron_smearing_sigma->push_back((*ele_smearing_sigma)[ele]);
		myEvent.electron_smearing_sigma_phi_up->push_back((*ele_smearing_sigma_phi_up)[ele]);
		myEvent.electron_smearing_sigma_phi_down->push_back((*ele_smearing_sigma_phi_down)[ele]);
		myEvent.electron_smearing_sigma_rho_up->push_back((*ele_smearing_sigma_rho_up)[ele]);
		myEvent.electron_smearing_sigma_rho_down->push_back((*ele_smearing_sigma_rho_down)[ele]);
		myEvent.electron_charge->push_back(ele->charge());
		myEvent.electron_r9->push_back(ele->full5x5_r9());
	}

	for (size_t i = 0; i < muons->size(); ++i) {
		const auto mu = muons->ptrAt(i);
		TLorentzVector p4;
		p4.SetPtEtaPhiM(mu->pt(), mu->eta(), mu->phi(), mu->mass());
		myEvent.muons_p4->push_back(p4);
		myEvent.muon_charge->push_back(mu->charge());
		myEvent.muon_IDSF_central->push_back((*muon_IDSF)[mu]);
		myEvent.muon_IsoSF_central->push_back((*muon_IsoSF)[mu]);
		myEvent.muon_IDSF_error->push_back((*muon_IDSF_error)[mu]);
		myEvent.muon_IsoSF_error->push_back((*muon_IsoSF_error)[mu]);
	}

	for (size_t i = 0; i < jets->size(); ++i) {
		const auto jet = jets->ptrAt(i);
		TLorentzVector p4;
		p4.SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), jet->mass());
		myEvent.jets_p4->push_back(p4);
		myEvent.jec_uncertainty->push_back((*jec_unc)[jet]);
		myEvent.jetResolution->push_back((*jetResolution)[jet]);
		myEvent.JER_sf->push_back((*JERsf)[jet]);
		myEvent.JER_sf_up->push_back((*JERsf_up)[jet]);
		myEvent.JER_sf_down->push_back((*JERsf_down)[jet]);
		myEvent.genJetPt->push_back((*genjetPt)[jet]);
		myEvent.genJetMatch->push_back((*genjetMatch)[jet]);
	}

	tree->Fill();
}

DEFINE_FWK_MODULE(miniTTree);
