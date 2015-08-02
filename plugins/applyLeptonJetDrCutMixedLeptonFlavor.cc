// -*- C++ -*-
//
// 

/**\class applyLeptonJetDrCutMixedLeptonFlavor applyLeptonJetDrCutMixedLeptonFlavor.cc doubleElectronTracklessTrigger/applyLeptonJetDrCutMixedLeptonFlavor/plugins/applyLeptonJetDrCutMixedLeptonFlavor.cc
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

class applyLeptonJetDrCutMixedLeptonFlavor : public edm::EDProducer {
   public:
      explicit applyLeptonJetDrCutMixedLeptonFlavor(const edm::ParameterSet&);
      ~applyLeptonJetDrCutMixedLeptonFlavor();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	  ///calculate the dilepton mass using two const_iterators to reco::Candidate objects
	  double getDileptonMass(edm::OwnVector<reco::Candidate>::const_iterator& one, edm::OwnVector<reco::Candidate>::const_iterator& two){
		  double mass = TMath::Sqrt( 2*(one->pt())*(two->pt())*( TMath::CosH( (one->eta())-(two->eta()) ) - TMath::Cos( (one->phi())-(two->phi()) ) ) );
		  return mass;
	  }///end getDileptonMass()

	  void findLeadingAndSubleading(edm::OwnVector<reco::Candidate>::const_iterator& first, edm::Handle<edm::OwnVector<reco::Candidate> > collectionOne, edm::OwnVector<reco::Candidate>::const_iterator& second, edm::Handle<edm::OwnVector<reco::Candidate> > collectionTwo,bool doDileptonMassCut){

#ifdef DEBUG
		  std::cout<<"checking pt of particles in two handled findLeadingAndSubleading fxn"<<std::endl;
#endif

		  ///find the highest pT object in collectionOne by looping over all contents in collectionOne
		  for(edm::OwnVector<reco::Candidate>::const_iterator genItOne = collectionOne->begin(); genItOne != collectionOne->end(); genItOne++){
#ifdef DEBUG
			  std::cout<<"a particle from collectionOne has pT = \t"<< genItOne->pt() << std::endl;
#endif

			  if(first==collectionOne->end()) first=genItOne;
			  else if(genItOne->pt() > first->pt() ) first = genItOne;
		  }//end loop over reco::Candidate objects in collectionOne


		  if(!doDileptonMassCut){
			  ///now find the highest pT object in collectionTwo which is separated from the highest pT object chosen from collectionOne
			  for(edm::OwnVector<reco::Candidate>::const_iterator genItTwo = collectionTwo->begin(); genItTwo != collectionTwo->end(); genItTwo++){
#ifdef DEBUG
				  std::cout<<"a particle from collectionTwo has pT = \t"<< genItTwo->pt() << std::endl;
#endif

				  if(second==collectionTwo->end() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > 0.4 ) second = genItTwo;

				  else if(genItTwo->pt() > second->pt() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > 0.4) second = genItTwo;
			  }///end loop over objects in collectionTwo

		  }///end if(!doDileptonMassCut)

		  if(doDileptonMassCut){
			  ///now find the highest pT object in collectionTwo which is separated from the highest pT object chosen from collectionOne
			  ///and whose dilepton mass, when combined with the object chosen from collectionOne, is above the dileptonMass threshold value
#ifdef DEBUG
			  std::cout<<"inside if(doDileptonMassCut) of two handled findLeadingAndSubleading() fxn"<<std::endl;
#endif

			  for(edm::OwnVector<reco::Candidate>::const_iterator genItTwo = collectionTwo->begin(); genItTwo != collectionTwo->end(); genItTwo++){
#ifdef DEBUG
				  std::cout<<"a particle from collectionTwo has pT = \t"<< genItTwo->pt() << std::endl;
#endif
			  
				  if(second==collectionTwo->end() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > 0.4 && getDileptonMass(first,genItTwo) > minDileptonMass ) second = genItTwo;

				  else if(genItTwo->pt() > second->pt() && deltaR(first->eta(), first->phi(), genItTwo->eta(), genItTwo->phi()) > 0.4 && getDileptonMass(first,genItTwo) > minDileptonMass ) second = genItTwo;
			  }///end loop over objects in collectionTwo with dilepton mass cut applied

		  }///end if(doDileptonMassCut)

#ifdef DEBUG
		  std::cout<<"leaving two handled findLeadingAndSubleading fxn"<<std::endl;
#endif

	  }///end two handled findLeadingAndSubleading()


	  /**this fxn checks that the object pointed to by objIt does not already exist
	   * in the collection pointed at by ptrToObjColl
	   * returns false if objIt already exists in the collection pointed to by ptrToObjColl
	   */
	  bool isNotDuplicate(edm::OwnVector<reco::Candidate>::const_iterator & objIt,
			  std::auto_ptr<edm::OwnVector<reco::Candidate> >& ptrToObjColl){
		  if(ptrToObjColl->size()==0) return true;
		  for(unsigned int i=0; i<ptrToObjColl->size(); i++){
#ifdef DEBUG
			  std::cout<<"about to check if reco::Candidate object has already been added to another collection"<<std::endl;
			  std::cout<<"size of other collection = \t"<< ptrToObjColl->size() <<std::endl;
#endif
			  if(objIt->pt()==(*ptrToObjColl)[i].pt() && objIt->eta()==(*ptrToObjColl)[i].eta() 
					  && objIt->phi()==(*ptrToObjColl)[i].phi()) return false;

		  }///end loop over objects in collection pointed to by ptrToObjColl
		  return true;
	  }///end isNotDuplicate()


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
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputLeptonsOneToken;	///< leptons selected earlier
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > inputLeptonsTwoToken;	///< leptons selected earlier


	  std::string outputCollName;
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
applyLeptonJetDrCutMixedLeptonFlavor::applyLeptonJetDrCutMixedLeptonFlavor(const edm::ParameterSet& iConfig):
	outputCollName(iConfig.getParameter<std::string>("outputJetsCollectionName")),
	dRSeparation(iConfig.getParameter<double>("minDrSeparation")),
	minDileptonMass(iConfig.getParameter<double>("minDileptonMassCut"))

{
   //edm::Service<TFileService> fs;
   
   ///register the input collections	
   inputJetsToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputJetsCollTag"));
   inputLeptonsOneToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputLeptonsOneCollTag"));
   inputLeptonsTwoToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("inputLeptonsTwoCollTag"));

   ///register the collections which are added to the event
   produces<edm::OwnVector<reco::Candidate> >(outputCollName);

}


applyLeptonJetDrCutMixedLeptonFlavor::~applyLeptonJetDrCutMixedLeptonFlavor()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
applyLeptonJetDrCutMixedLeptonFlavor::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef DEBUG
   std::cout<<"entered applyLeptonJetDrCutMixedLeptonFlavor produce method"<<std::endl;
#endif

   Handle<edm::OwnVector<reco::Candidate> > inputJetsColl;
   iEvent.getByToken(inputJetsToken, inputJetsColl);

   Handle<edm::OwnVector<reco::Candidate> > inputLeptonsOneColl;
   iEvent.getByToken(inputLeptonsOneToken, inputLeptonsOneColl);

   Handle<edm::OwnVector<reco::Candidate> > inputLeptonsTwoColl;
   iEvent.getByToken(inputLeptonsTwoToken, inputLeptonsTwoColl);


#ifdef DEBUG
   std::cout<<"made handles to input collections"<<std::endl;
#endif

   ///make an empty output collection to eventually hold jets, and a pointer to this collection
   std::auto_ptr<edm::OwnVector<reco::Candidate> > outputObjColl(new edm::OwnVector<reco::Candidate>());

   ///find the two highest pT, separated leptons whose dilepton mass is > some threshold, and assign const_iterators to these two leptons
   edm::OwnVector<reco::Candidate>::const_iterator leadLeptonOne=inputLeptonsOneColl->end(), leadLeptonTwo=inputLeptonsTwoColl->end();
   findLeadingAndSubleading(leadLeptonOne, inputLeptonsOneColl, leadLeptonTwo, inputLeptonsTwoColl, true);
   if(leadLeptonOne==inputLeptonsOneColl->end() || leadLeptonTwo==inputLeptonsTwoColl->end() ){
	   iEvent.put(outputObjColl, outputCollName);
	   return;
   }

   ///now that the two hardest, separated leptons have been found (one from each input lepton collection)
   ///two const_iterators should be declared and initialized to const_iterators tied to the two hardest leptons
   edm::OwnVector<reco::Candidate>::const_iterator leadLepton, subleadLepton;	///< use these const_iterators for further analysis
   if(leadLeptonOne->pt() > leadLeptonTwo->pt() ){
	   leadLepton = leadLeptonOne;
	   subleadLepton = leadLeptonTwo;
   }
   else{
	   leadLepton = leadLeptonTwo;
	   subleadLepton = leadLeptonOne;
   }
   

   ///for each possible jet object, check that the jet is at least dRSeparation away from the two hardest leptons which
   ///have dilepton mass > some threshold 
   for(edm::OwnVector<reco::Candidate>::const_iterator jetIt=inputJetsColl->begin();jetIt!=inputJetsColl->end(); jetIt++){
	   double dRleadLepton=deltaR(jetIt->eta(),jetIt->phi(),leadLepton->eta(),leadLepton->phi());
	   double dRsubleadLepton=deltaR(jetIt->eta(),jetIt->phi(),subleadLepton->eta(),subleadLepton->phi());
	   
	   if(dRleadLepton > dRSeparation){
		   if(dRsubleadLepton > dRSeparation){
			   outputObjColl->push_back(*jetIt);
		   }///end filter on dRsubleadLepton > dRSeparation
	   }///end filter on dRleadLepton > dRSeparation
   
   }///end loop over reco::Candidate leptons


#ifdef DEBUG
   std::cout<<"about to put collection of reco::Candidate objects into root file"<<std::endl;
#endif
  
   ///now put the collection of matched higher level reco::Candidate objects into the event
   iEvent.put(outputObjColl, outputCollName);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
applyLeptonJetDrCutMixedLeptonFlavor::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
applyLeptonJetDrCutMixedLeptonFlavor::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
applyLeptonJetDrCutMixedLeptonFlavor::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
applyLeptonJetDrCutMixedLeptonFlavor::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
applyLeptonJetDrCutMixedLeptonFlavor::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
applyLeptonJetDrCutMixedLeptonFlavor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
applyLeptonJetDrCutMixedLeptonFlavor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(applyLeptonJetDrCutMixedLeptonFlavor);
