// -*- C++ -*-
//
// Package:    testFilter/hasNoHighMassWrObjects
// Class:      hasNoHighMassWrObjects
// 
/**\class hasNoHighMassWrObjects hasNoHighMassWrObjects.cc testFilter/hasNoHighMassWrObjects/plugins/hasNoHighMassWrObjects.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  sean kalafut
//         Created:  Tue, 04 Aug 2015 22:09:29 GMT
//
//


// system include files
#include <memory>
#include <cstring>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/Math/interface/LorentzVector.h"

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



//#define DEBUG
//
// class declaration
//

class hasNoHighMassWrObjects : public edm::EDFilter {
   public:
      explicit hasNoHighMassWrObjects(const edm::ParameterSet&);
      ~hasNoHighMassWrObjects();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      // ----------member data ---------------------------
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputLeadLeptonsToken;
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputSubleadLeptonsToken;
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputJetsToken;
	  double maxWrMass;
	  double maxDileptonMass;

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
hasNoHighMassWrObjects::hasNoHighMassWrObjects(const edm::ParameterSet& iConfig):
	maxWrMass(iConfig.getParameter<double>("maxWrMass")),
	maxDileptonMass(iConfig.getParameter<double>("maxDileptonMass"))
{
   //now do what ever initialization is needed
   inputLeadLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputLeadLeptonsCollTag"));
   inputSubleadLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputSubleadLeptonsCollTag"));
   inputJetsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputJetsCollTag"));

}


hasNoHighMassWrObjects::~hasNoHighMassWrObjects()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
hasNoHighMassWrObjects::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;

#ifdef DEBUG
	std::cout<<"in filter method of hasNoHighMassWrObjects class"<<std::endl;
#endif
	
	Handle<edm::OwnVector<reco::Candidate> > inputLeadLeptonsObjectColl;
	iEvent.getByToken(inputLeadLeptonsToken, inputLeadLeptonsObjectColl);
	
	Handle<edm::OwnVector<reco::Candidate> > inputSubleadLeptonsObjectColl;
	iEvent.getByToken(inputSubleadLeptonsToken, inputSubleadLeptonsObjectColl);

	Handle<edm::OwnVector<reco::Candidate> > inputJetsObjectColl;
	iEvent.getByToken(inputJetsToken, inputJetsObjectColl);

	///loop over all possible combinations of LLJJ objects, and check that none have 
	///mass > maxWrMass
	for(edm::OwnVector<reco::Candidate>::const_iterator leptIt=inputLeadLeptonsObjectColl->begin(); leptIt!=inputLeadLeptonsObjectColl->end();
		   leptIt++){
		
		for(edm::OwnVector<reco::Candidate>::const_iterator leptTwoIt=inputSubleadLeptonsObjectColl->begin(); leptTwoIt!=inputSubleadLeptonsObjectColl->end();
				leptTwoIt++){
			if(leptTwoIt==leptIt) continue;

			for(edm::OwnVector<reco::Candidate>::const_iterator jetIt=inputJetsObjectColl->begin(); jetIt!=inputJetsObjectColl->end();
					jetIt++){

				for(edm::OwnVector<reco::Candidate>::const_iterator jetTwoIt=inputJetsObjectColl->begin(); jetTwoIt!=inputJetsObjectColl->end();
						jetTwoIt++){
					if(jetTwoIt==jetIt) continue;
					ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > l1 = leptIt->p4(), l2 = leptTwoIt->p4();
					ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > j1 = jetIt->p4(), j2 = jetTwoIt->p4();
					if( (l1+l2+j1+j2).M() > maxWrMass ) return false;
					if( (l1+l2).M() > maxDileptonMass ) return false;

				}///end loop over second jet
			}///end loop over first jet
		}///end loop over sublead leptons
	}///end loop over lead leptons

	return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
hasNoHighMassWrObjects::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hasNoHighMassWrObjects::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
hasNoHighMassWrObjects::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
hasNoHighMassWrObjects::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
hasNoHighMassWrObjects::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
hasNoHighMassWrObjects::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hasNoHighMassWrObjects::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(hasNoHighMassWrObjects);
