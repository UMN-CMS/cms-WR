// -*- C++ -*-
//
// 

/**\class applyLeptonJetDrCut applyLeptonJetDrCut.cc doubleElectronTracklessTrigger/applyLeptonJetDrCut/plugins/applyLeptonJetDrCut.cc
 Description: [one line class summary]


 Implementation:
     [Notes on implementation]
*/
//
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
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/Common/interface/OwnVector.h"



#include "DataFormats/EgammaReco/interface/SuperCluster.h"
//#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TTree.h"

//#define NOBJ 500
//#define DEBUG

//
// class declaration
//

class applyLeptonJetDrCut : public edm::EDProducer {
   public:
      explicit applyLeptonJetDrCut(const edm::ParameterSet&);
      ~applyLeptonJetDrCut();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	  ///calculate the dilepton mass using two const_iterators to reco::Candidate objects
	  double getDileptonMass(edm::OwnVector<reco::Candidate>::const_iterator& one, edm::OwnVector<reco::Candidate>::const_iterator& two){
		  double mass = TMath::Sqrt( 2*(one->pt())*(two->pt())*( TMath::CosH( (one->eta())-(two->eta()) ) - TMath::Cos( (one->phi())-(two->phi()) ) ) );
		  return mass;
	  }///end getDileptonMass()

	  void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collection,bool doDileptonMassCut){

#ifdef DEBUG
		  std::cout<<"checking pt of particles in findLeadingAndSubleading fxn"<<std::endl;
#endif

		  if(!doDileptonMassCut){

			  for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++){
#ifdef DEBUG
				  std::cout<<"pT = \t"<< genIt->pt() << std::endl;
#endif
				  if(first==collection->end()) first=genIt;
				  else{
					  if(genIt->pt() > first->pt() && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ){
						  second = first;
						  first = genIt;
					  }
					  else if( ( second==collection->end() || genIt->pt() > second->pt() ) && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ) second = genIt;
				  }
			  }//end loop over reco::Candidate collection

		  }///end if(!doDileptonMassCut)

		  if(doDileptonMassCut){

			  for(edm::OwnVector<reco::Candidate>::const_iterator genIt = collection->begin(); genIt != collection->end(); genIt++){
#ifdef DEBUG
				  std::cout<<"pT = \t"<< genIt->pt() << std::endl;
#endif
				  if(first==collection->end()) first=genIt;
				  else{
					  if(genIt->pt() > first->pt() && getDileptonMass(first,genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ){
						  second = first;
						  first = genIt;
					  }
					  else if( (second==collection->end() || genIt->pt() > second->pt()) && getDileptonMass(first,genIt) > minDileptonMass && deltaR(genIt->eta(), genIt->phi(), first->eta(), first->phi() ) > 0.4 ) second = genIt;
				  }
			  }//end loop over reco::Candidate collection

		  }///end if(doDileptonMassCut)

#ifdef DEBUG
		  std::cout<<"leaving findLeadingAndSubleading fxn"<<std::endl;
#endif

	  }///end findLeadingAndSubleading()


   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 
      // ----------member data ---------------------------
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputJetsToken;
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputLeptonsToken;	///< leptons which have passed dilepton mass and earlier cuts

	  std::string outputCollName;
	  std::string outputTwoCollName;
	  double dRSeparation;	///< minimum dR separation btwn a lepton and jet
	  double minDileptonMass;	///< min dilepton mass

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
applyLeptonJetDrCut::applyLeptonJetDrCut(const edm::ParameterSet& iConfig):
	outputCollName(iConfig.getParameter<std::string>("outputLeptonsCollectionName")),
	outputTwoCollName(iConfig.getParameter<std::string>("outputJetsCollectionName")),
	dRSeparation(iConfig.getParameter<double>("minDrSeparation")),
	minDileptonMass(iConfig.getParameter<double>("minDileptonMassCut"))

{
   
   //edm::Service<TFileService> fs;
   
   ///register the input collections	
   inputJetsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputJetsCollTag"));
   inputLeptonsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputLeptonsCollTag"));

   ///register the collections which are added to the event
   produces<edm::OwnVector<reco::Candidate> >(outputCollName);
   produces<edm::OwnVector<reco::Candidate> >(outputTwoCollName);

}


applyLeptonJetDrCut::~applyLeptonJetDrCut()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
applyLeptonJetDrCut::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef DEBUG
   std::cout<<"entered applyLeptonJetDrCut produce method"<<std::endl;
#endif

   Handle<edm::OwnVector<reco::Candidate> > inputJetsColl;
   iEvent.getByToken(inputJetsToken, inputJetsColl);

   Handle<edm::OwnVector<reco::Candidate> > inputLeptonsColl;
   iEvent.getByToken(inputLeptonsToken, inputLeptonsColl);


#ifdef DEBUG
   std::cout<<"made handles to input collections"<<std::endl;
#endif

   ///make two empty output collections to hold leptons jets, and pointers to these collections
   std::auto_ptr<edm::OwnVector<reco::Candidate> > outputObjColl(new edm::OwnVector<reco::Candidate>());	///leptons
   std::auto_ptr<edm::OwnVector<reco::Candidate> > outputObjTwoColl(new edm::OwnVector<reco::Candidate>());	///jets

   ///loop over all input leptons, and for each lepton loop over the two leading input jets 
   ///the input jets have already passed pT, eta, and ID filters
   ///add leptons to the output collection which are separated from all input jets by dR > 0.4
   ///add two leading jets to a separate output collection
   typedef edm::OwnVector<reco::Candidate>::const_iterator recoCandIt;
   recoCandIt leadJet = inputJetsColl->end(), subleadJet = inputJetsColl->end();
   findLeadingAndSubleading(leadJet, subleadJet, inputJetsColl, false);

   if(leadJet == inputJetsColl->end() || subleadJet == inputJetsColl->end()) return;

   ///add the two leading jets to the outputObjTwoColl
   outputObjTwoColl->push_back(*leadJet);
   outputObjTwoColl->push_back(*subleadJet);

   for(recoCandIt lepton=inputLeptonsColl->begin(); lepton!=inputLeptonsColl->end(); lepton++){
	   double leadJetdRseparation = deltaR(leadJet->eta(),leadJet->phi(),lepton->eta(),lepton->phi());
	   double subleadJetdRseparation = deltaR(subleadJet->eta(),subleadJet->phi(),lepton->eta(),lepton->phi());
	   if(leadJetdRseparation > dRSeparation && subleadJetdRseparation > dRSeparation) outputObjColl->push_back(*lepton);
   }///end loop over input leptons

#ifdef DEBUG
   std::cout<<"about to put collections of reco::Candidate objects into root file"<<std::endl;
#endif
  
   ///now put the collection of matched higher level reco::Candidate objects into the event
   iEvent.put(outputObjColl, outputCollName);
   iEvent.put(outputObjTwoColl, outputTwoCollName);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
applyLeptonJetDrCut::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
applyLeptonJetDrCut::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
applyLeptonJetDrCut::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
applyLeptonJetDrCut::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
applyLeptonJetDrCut::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
applyLeptonJetDrCut::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
applyLeptonJetDrCut::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(applyLeptonJetDrCut);
