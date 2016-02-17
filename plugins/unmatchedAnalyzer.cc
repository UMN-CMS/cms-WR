// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/unmatchedAnalyzer
// Class:      unmatchedAnalyzer
//
/**\class unmatchedAnalyzer unmatchedAnalyzer.cc doubleElectronTracklessTrigger/unmatchedAnalyzer/plugins/unmatchedAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sean Kalafut
//         Created:  Wed, 15 April 2015
//
//


// system include files
#include <memory>
#include <map>
#include <utility>
#include <cstring>
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/MET.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TTree.h"
//#include "TLorentzVector.h"
#include <TFile.h>
#include <TBranch.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include "TCollection.h"

//#define DEBUG

//
// class declaration
//

class unmatchedAnalyzer : public edm::EDAnalyzer
{
public:
	explicit unmatchedAnalyzer(const edm::ParameterSet&);
	~unmatchedAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	///calculate the dilepton mass using two const_iterators to reco::Candidate objects
	double getDileptonMass(edm::OwnVector<reco::Candidate>::const_iterator& one, edm::OwnVector<reco::Candidate>::const_iterator& two)
	{
		double mass = TMath::Sqrt( 2 * (one->pt()) * (two->pt()) * ( TMath::CosH( (one->eta()) - (two->eta()) ) - TMath::Cos( (one->phi()) - (two->phi()) ) ) );
		return mass;
	}///end getDileptonMass()

	void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection, bool doDileptonMassCut)
	{

#ifdef DEBUG
		std::cout << "checking pt of particles in findLeadingAndSubleading fxn" << std::endl;
#endif

		if(!doDileptonMassCut) {

			for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++) {
#ifdef DEBUG
				std::cout << "pT = \t" << genIt->pt() << std::endl;
#endif
				if(first == collection->end()) first = genIt;
				else {
					if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ) {
						second = first;
						first = genIt;
					} else if( (second == collection->end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4  ) second = genIt;
				}
			}//end loop over reco::Candidate collection

		}///end if(!doDileptonMassCut)

		if(doDileptonMassCut) {

			for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++) {
#ifdef DEBUG
				std::cout << "pT = \t" << genIt->pt() << std::endl;
#endif
				if(first == collection->end()) first = genIt;
				else {
					if(genIt->pt() > first->pt() && getDileptonMass(first, genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ) {
						second = first;
						first = genIt;
					} else if( (second == collection->end() || genIt->pt() > second->pt()) && getDileptonMass(first, genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ) second = genIt;
				}
			}//end loop over reco::Candidate collection

		}///end if(doDileptonMassCut)

#ifdef DEBUG
		std::cout << "leaving findLeadingAndSubleading fxn" << std::endl;
#endif

	}///end findLeadingAndSubleading()

	/** this fxn sifts through a collection of reco::Candidate objects, finds the highest pT object in the collection, and assigns
	 * a pointer to this object to the iterator named iter
	 * a reference to iter is input to this fxn
	 */
	void findHighestPt(edm::OwnVector<reco::Candidate>::const_iterator& iter, edm::Handle<edm::OwnVector<reco::Candidate> > coll)
	{
		for(edm::OwnVector<reco::Candidate>::const_iterator it = coll->begin(); it != coll->end(); it++) {
			if(iter == coll->end()) iter = it;
			else {
				if(it->pt() > iter->pt()) iter = it;
			}
		}///end loop over reco::Candidate objects in the collection named coll

	}///end findHighestPt()

	/// this fxn finds the highest MET object, and assigns it to a const_iterator of type std::vector<pat::MET>
	void findHighestMET(std::vector<pat::MET>::const_iterator & iter, edm::Handle<std::vector<pat::MET> > coll)
	{
		for(std::vector<pat::MET>::const_iterator it = coll->begin(); it != coll->end(); it++) {
			if(iter == coll->end()) iter = it;
			else {
				if(it->et() > iter->et()) iter = it;
			}
		}///end loop over pat::MET objects in the collection named coll

	}///end findHighestMET()


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;



// ----------member data ---------------------------

	std::string tName;
	bool applyDileptonMassCut;
	double minDileptonMass;

///Handles to RECO object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > leptons;
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > jets;
	edm::Handle<GenEventInfoProduct> genEvtInfo;
	edm::Handle<std::vector<reco::Vertex> > vertices;
	edm::Handle<std::vector<pat::MET> > met;

///tokens to input collections
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > leptonsToken;
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > jetsToken;
	edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;
	edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;
	edm::EDGetTokenT<std::vector<pat::MET> > metToken;


	TTree * tree;

	Int_t runNumber;
	ULong64_t evtNumber;

	Int_t nJets;
	Int_t nLeptons;
	Int_t nVertices;
	Int_t nMETs;

	Float_t missingET;

//first element is leading (highest pT) electron
//second element is subleading electron
	Float_t etaEle[2];
	Float_t ptEle[2];
	Float_t phiEle[2];
	Float_t dileptonMass;

	Float_t etaJet[2];
	Float_t ptJet[2];
	Float_t phiJet[2];
	Float_t dijetMass;

///three and four object masses
	Float_t fourObjectMass;
	Float_t leadLeptonThreeObjMass;
	Float_t subleadingLeptonThreeObjMass;

///deltaR between each lepton and both jets
	Float_t dR_leadingLeptonLeadingJet;
	Float_t dR_leadingLeptonSubleadingJet;
	Float_t dR_subleadingLeptonLeadingJet;
	Float_t dR_subleadingLeptonSubleadingJet;
	Float_t dR_leadingLeptonSubleadingLepton;
	Float_t dR_leadingJetSubleadingJet;


///leadingIsHardest = 1 when leading ele pT > subleading ele pT
///leadingIsHardest = 0 when leading ele pT < subleading ele pT
	Int_t leadingIsHardest;

///pT, eta, and phi of the heavy Nu and WR
	Float_t etaHvyNu;
	Float_t ptHvyNu;
	Float_t phiHvyNu;
	Float_t etaWr;
	Float_t ptWr;
	Float_t phiWr;

	Float_t evWeight;	///< weight of the event.  defaults to 1, only changes for bkgnd MC
	Float_t evWeightSign;	///< if the sign of evWeight is negative, then the evt should not be plotted in a histo

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

unmatchedAnalyzer::unmatchedAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName")),
	applyDileptonMassCut(iConfig.getParameter<bool>("doDileptonMassCut")),
	minDileptonMass(iConfig.getParameter<double>("minDileptonMass"))

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>(tName.c_str(), "event kinematic info");

	tree->Branch("etaEle", etaEle, "etaEle[2]/F");
	tree->Branch("ptEle", ptEle, "ptEle[2]/F");
	tree->Branch("phiEle", phiEle, "phiEle[2]/F");
	tree->Branch("dileptonMass", &dileptonMass, "dileptonMass/F");

	tree->Branch("leadingIsHardest", &leadingIsHardest, "leadingIsHardest/I");

	tree->Branch("runNumber", &runNumber, "runNumber/I");
	tree->Branch("evtNumber", &evtNumber, "evtNumber/l");

	tree->Branch("nJets", &nJets, "nJets/I");
	tree->Branch("nLeptons", &nLeptons, "nLeptons/I");
	tree->Branch("nVertices", &nVertices, "nVertices/I");
	tree->Branch("nMETs", &nMETs, "nMETs/I");

	tree->Branch("missingET", &missingET, "missingET/F");

	tree->Branch("etaJet", etaJet, "etaJet[2]/F");
	tree->Branch("ptJet", ptJet, "ptJet[2]/F");
	tree->Branch("phiJet", phiJet, "phiJet[2]/F");
	tree->Branch("dijetMass", &dijetMass, "dijetMass/F");

	tree->Branch("fourObjectMass", &fourObjectMass, "fourObjectMass/F");
	tree->Branch("leadLeptonThreeObjMass", &leadLeptonThreeObjMass, "leadLeptonThreeObjMass/F");
	tree->Branch("subleadingLeptonThreeObjMass", &subleadingLeptonThreeObjMass, "subleadingLeptonThreeObjMass/F");

	tree->Branch("etaHvyNu", &etaHvyNu, "etaHvyNu/F");
	tree->Branch("ptHvyNu", &ptHvyNu, "ptHvyNu/F");
	tree->Branch("phiHvyNu", &phiHvyNu, "phiHvyNu/F");

	tree->Branch("etaWr", &etaWr, "etaWr/F");
	tree->Branch("ptWr", &ptWr, "ptWr/F");
	tree->Branch("phiWr", &phiWr, "phiWr/F");

	tree->Branch("dR_leadingLeptonLeadingJet", &dR_leadingLeptonLeadingJet, "dR_leadingLeptonLeadingJet/F");
	tree->Branch("dR_leadingLeptonSubleadingJet", &dR_leadingLeptonSubleadingJet, "dR_leadingLeptonSubleadingJet/F");
	tree->Branch("dR_subleadingLeptonLeadingJet", &dR_subleadingLeptonLeadingJet, "dR_subleadingLeptonLeadingJet/F");
	tree->Branch("dR_subleadingLeptonSubleadingJet", &dR_subleadingLeptonSubleadingJet, "dR_subleadingLeptonSubleadingJet/F");
	tree->Branch("dR_leadingLeptonSubleadingLepton", &dR_leadingLeptonSubleadingLepton, "dR_leadingLeptonSubleadingLepton/F");
	tree->Branch("dR_leadingJetSubleadingJet", &dR_leadingJetSubleadingJet, "dR_leadingJetSubleadingJet/F");

	tree->Branch("evWeight", &evWeight, "evWeight/F");
	tree->Branch("evWeightSign", &evWeightSign, "evWeightSign/F");

	leptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("leptonsCollection"));
	jetsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("jetsCollection"));
	genEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
	verticesToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));
	metToken = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));

}


unmatchedAnalyzer::~unmatchedAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
unmatchedAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout << "in analyze method of unmatchedAnalyzer class" << std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	evWeight = 1.0;
	evWeightSign = 1.0;

	iEvent.getByToken(leptonsToken, leptons);
	iEvent.getByToken(jetsToken, jets);
	iEvent.getByToken(verticesToken, vertices);

	iEvent.getByToken(genEventInfoToken, genEvtInfo);	///< get evt weights if analyzing MC

	if(genEvtInfo.isValid() ) {
		///real data does not have gen lvl event weights
		evWeight = genEvtInfo->weight();
		if(evWeight < 0) evWeightSign = -1.0;
	}

	///get MET info if it is available (only from MINIAOD files)
	missingET = 0, nMETs = 0;
	iEvent.getByToken(metToken, met);
	if(met.isValid() ) {
		nMETs = met->size();
		///find the highest MET object in the met handle
		std::vector<pat::MET>::const_iterator leadMET = met->end();
		findHighestMET(leadMET, met);
		///pat::MET objects inherit methods defined for reco::Candidate objects, like eta(), et(), and pt()
		if(leadMET != met->end() ) missingET = leadMET->et();
	}

	///get the number of vertices, jets, and leptons in the event
	nVertices = 1;
	if(vertices.isValid() ) {
		///GEN evt root files do not have a collection of reco::Vertex objects
		nVertices = vertices->size();
	}
	nJets = jets->size();
	nLeptons = leptons->size();

	edm::OwnVector<reco::Candidate>::const_iterator leadingLepton = leptons->end(), subleadingLepton = leptons->end();
	edm::OwnVector<reco::Candidate>::const_iterator leadingJet = jets->end(), subleadingJet = jets->end();

	findLeadingAndSubleading(leadingLepton, subleadingLepton, leptons, applyDileptonMassCut);
	findLeadingAndSubleading(leadingJet, subleadingJet, jets, false);

	if(leadingLepton == leptons->end() || subleadingLepton == leptons->end()) return;	///< skip this evt if two leptons are not found
	if(leadingJet == jets->end() || subleadingJet == jets->end()) return;	///< skip this evt if two jets are not found



	///now that the leading and subleading leptons and jets have been found, fill all of the arrays and single Float_ values
	///which will be saved into the tree
#ifdef DEBUG
	std::cout << "leading lepton pt: \t" << leadingLepton->pt() << std::endl;
	std::cout << "leading lepton pt: \t" << leadingLepton->pt() << std::endl;
	std::cout << "subleading lepton pt: \t" << subleadingLepton->pt() << std::endl;
#endif

	etaEle[0] = leadingLepton->eta();
	ptEle[0] = leadingLepton->pt();
	phiEle[0] = leadingLepton->phi();
	etaEle[1] = subleadingLepton->eta();
	ptEle[1] = subleadingLepton->pt();
	phiEle[1] = subleadingLepton->phi();
	if(ptEle[0] > ptEle[1]) {
		leadingIsHardest = 1;
	} else leadingIsHardest = 0;


#ifdef DEBUG
	std::cout << "leading jet pt: \t" << leadingJet->pt() << std::endl;
	std::cout << "subleading jet pt: \t" << subleadingJet->pt() << std::endl;
#endif

	etaJet[0] = leadingJet->eta();
	ptJet[0] = leadingJet->pt();
	phiJet[0] = leadingJet->phi();
	etaJet[1] = subleadingJet->eta();
	ptJet[1] = subleadingJet->pt();
	phiJet[1] = subleadingJet->phi();

	///make TLorentzVector objects for the four GEN objects.  Then use these LorentzVectors to calculate three and four object
	///invariant mass values.
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > l1 = leadingLepton->p4(), l2 = subleadingLepton->p4();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > j1 = leadingJet->p4(), j2 = subleadingJet->p4();
	dileptonMass = (l1 + l2).M(), dijetMass = (j1 + j2).M();
	leadLeptonThreeObjMass = (l1 + j1 + j2).M(), subleadingLeptonThreeObjMass = (l2 + j1 + j2).M();
	fourObjectMass = (l1 + l2 + j1 + j2).M();
	etaHvyNu = (l2 + j1 + j2).Eta(), ptHvyNu = (l2 + j1 + j2).Pt(), phiHvyNu = (l2 + j1 + j2).Phi();
	etaWr = (l1 + l2 + j1 + j2).Eta(), ptWr = (l1 + l2 + j1 + j2).Pt(), phiWr = (l1 + l2 + j1 + j2).Phi();

#ifdef DEBUG
	std::cout << "dilepton mass = \t" << dileptonMass << std::endl;
	std::cout << "dijet mass = \t" << dijetMass << std::endl;
	std::cout << "\t" << std::endl;
#endif

	///now use the individual GEN object eta and phi values to calculate dR between different (lepton, jet) pairs
	dR_leadingLeptonLeadingJet = deltaR(etaEle[0], phiEle[0], etaJet[0], phiJet[0]);
	dR_leadingLeptonSubleadingJet = deltaR(etaEle[0], phiEle[0], etaJet[1], phiJet[1]);
	dR_subleadingLeptonLeadingJet = deltaR(etaEle[1], phiEle[1], etaJet[0], phiJet[0]);
	dR_subleadingLeptonSubleadingJet = deltaR(etaEle[1], phiEle[1], etaJet[1], phiJet[1]);
	dR_leadingLeptonSubleadingLepton = deltaR(etaEle[0], phiEle[0], etaEle[1], phiEle[1]);
	dR_leadingJetSubleadingJet = deltaR(etaJet[0], phiJet[0], etaJet[1], phiJet[1]);

	tree->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
unmatchedAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
unmatchedAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void unmatchedAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(unmatchedAnalyzer);
