// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/hadronAnalyzer
// Class:      hadronAnalyzer
//
/**\class hadronAnalyzer hadronAnalyzer.cc doubleElectronTracklessTrigger/hadronAnalyzer/plugins/hadronAnalyzer.cc

 Description: [one line class summary]

 ///use this class to analyze several hadrons in an event, and save kinematic info which represents the all hadrons (like sumPT)

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
#include "TLorentzVector.h"
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

class hadronAnalyzer : public edm::EDAnalyzer
{
public:
	explicit hadronAnalyzer(const edm::ParameterSet&);
	~hadronAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;

// ----------member data ---------------------------

	std::string tName;

///Handles to RECO object collections
	edm::Handle<edm::OwnVector<reco::Candidate, edm::ClonePolicy<reco::Candidate> > > inputParticles;

///tokens to input collections
	edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputParticlesToken;


	TTree * tree;

	Int_t runNumber;
	ULong64_t evtNumber;

	Int_t nInputParticles;

//first element is leading (highest pT) electron
//second element is subleading electron
	Float_t sumPt;
	Float_t zPt;

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

hadronAnalyzer::hadronAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName"))

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;

	tree = fs->make<TTree>(tName.c_str(), "event kinematic info");

	tree->Branch("sumPt", &sumPt, "sumPt/F");
	tree->Branch("zPt", &sumPt, "sumPt/F");

	tree->Branch("runNumber", &runNumber, "runNumber/I");
	tree->Branch("evtNumber", &evtNumber, "evtNumber/l");

	tree->Branch("nInputParticles", &nInputParticles, "nInputParticles/I");

	inputParticlesToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputCollection"));

}


hadronAnalyzer::~hadronAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
hadronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout << "in analyze method of hadronAnalyzer class" << std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	sumPt = 0, zPt = 0;

	iEvent.getByToken(inputParticlesToken, inputParticles);
	nInputParticles = inputParticles->size();

	for(edm::OwnVector<reco::Candidate>::const_iterator it = inputParticles->begin(); it!=inputParticles->end(); it++){
		sumPt += it->pt();
	}//end loop over particles in input collection

	if(nInputParticles == 2){
		edm::OwnVector<reco::Candidate>::const_iterator one = inputParticles->begin();
		edm::OwnVector<reco::Candidate>::const_iterator two = (inputParticles->begin())++;
	
		TLorentzVector leadLeptonFourMom, subleadLeptonFourMom, zFourMom;
		leadLeptonFourMom.SetPtEtaPhiE(one->pt(), one->eta(), one->phi(), one->pt());
		subleadLeptonFourMom.SetPtEtaPhiE(two->pt(), two->eta(), two->phi(), two->pt());
		zFourMom = leadLeptonFourMom + subleadLeptonFourMom;
		zPt = zFourMom.Pt();
	}

	tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
hadronAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
hadronAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void hadronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(hadronAnalyzer);
