// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/singleParticleAnalyzer
// Class:      singleParticleAnalyzer
// 
/**\class singleParticleAnalyzer singleParticleAnalyzer.cc doubleElectronTracklessTrigger/singleParticleAnalyzer/plugins/singleParticleAnalyzer.cc

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

class singleParticleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit singleParticleAnalyzer(const edm::ParameterSet&);
      ~singleParticleAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	  /** this fxn sifts through a collection of reco::Candidate objects, finds the highest pT object in the collection, and assigns
	   * a pointer to this object to the iterator named iter
	   * a reference to iter is input to this fxn
	   */
	  void findHighestPt(edm::OwnVector<reco::Candidate>::const_iterator& iter, edm::Handle<edm::OwnVector<reco::Candidate> > coll){
		  for(edm::OwnVector<reco::Candidate>::const_iterator it = coll->begin(); it!=coll->end(); it++){
			  if(iter==coll->end()) iter=it;
			  else{
				  if(it->pt() > iter->pt()) iter=it;
			  }
		  }///end loop over reco::Candidate objects in the collection named coll

	  }///end findHighestPt()
	 

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;

// ----------member data ---------------------------

std::string tName;

///Handles to RECO object collections
edm::Handle<edm::OwnVector<reco::Candidate,edm::ClonePolicy<reco::Candidate> > > rightHandWs;
edm::Handle<GenEventInfoProduct> genEvtInfo;
//edm::Handle<std::vector<reco::Vertex> > vertices;

///tokens to input collections
edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > rightHandWsToken;
edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken;
//edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;


TTree * tree;

Int_t runNumber;
ULong64_t evtNumber;

Int_t nWRs;
Int_t status;
//Int_t nVertices;

//first element is leading (highest pT) electron
//second element is subleading electron
Float_t etaWr;
Float_t ptWr;
Float_t phiWr;
Float_t massWr;

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

singleParticleAnalyzer::singleParticleAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName"))

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"event kinematic info");

   tree->Branch("etaWr",&etaWr,"etaWr/F");
   tree->Branch("ptWr",&ptWr,"ptWr/F");
   tree->Branch("phiWr",&phiWr,"phiWr/F");
   tree->Branch("massWr",&massWr,"massWr/F");
   
   tree->Branch("runNumber",&runNumber,"runNumber/I");
   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
 
   tree->Branch("nWRs",&nWRs,"nWRs/I");
   tree->Branch("status",&status,"status/I");
   
   tree->Branch("evWeight",&evWeight,"evWeight/F");
   tree->Branch("evWeightSign",&evWeightSign,"evWeightSign/F");
 
   rightHandWsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("rightHandWsCollection"));
   genEventInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   //verticesToken = consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"));

}


singleParticleAnalyzer::~singleParticleAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
singleParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout<<"in analyze method of singleParticleAnalyzer class"<<std::endl;
#endif

	evtNumber = iEvent.id().event();
	runNumber = iEvent.id().run();
	evWeight = 1.0;
	evWeightSign = 1.0;

	iEvent.getByToken(rightHandWsToken, rightHandWs);
	//iEvent.getByToken(verticesToken, vertices);
	iEvent.getByToken(genEventInfoToken, genEvtInfo);	///< get evt weights if analyzing MC

	if(genEvtInfo.isValid() ){
		///real data does not have gen lvl event weights
		evWeight = genEvtInfo->weight();
		if(evWeight < 0) evWeightSign = -1.0;
	}

	//nVertices = vertices->size();
	nWRs = rightHandWs->size();

	edm::OwnVector<reco::Candidate>::const_iterator theParticle = rightHandWs->begin();

	status = theParticle->status();
	etaWr = theParticle->eta();
	ptWr = theParticle->pt();
	phiWr = theParticle->phi();
	massWr = theParticle->mass();

	tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
singleParticleAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
singleParticleAnalyzer::endJob() 
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void singleParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(singleParticleAnalyzer);
