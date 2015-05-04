// -*- C++ -*-
//
// 

/**\class FindHigherLevelMatchedObject FindHigherLevelMatchedObject.cc doubleElectronTracklessTrigger/FindHigherLevelMatchedObject/plugins/FindHigherLevelMatchedObject.cc
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

#define NOBJ 200
//#define DEBUG

//
// class declaration
//

class FindHigherLevelMatchedObject : public edm::EDProducer {
   public:
      explicit FindHigherLevelMatchedObject(const edm::ParameterSet&);
      ~FindHigherLevelMatchedObject();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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
	  /*
	  edm::EDGetTokenT<std::vector<reco::CompositeCandidate>> momToken;
	  edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> momParentOneToken;
      edm::EDGetTokenT<std::vector<reco::RecoEcalCandidate>> momParentTwoToken;
	
   	  std::string daughterOneCollection;
	  std::string daughterTwoCollection;
	  */
 	  
	  ///the point of this producer is to match an object from the lowLevel collection with
	  ///an object in the higherLevel collection.  An example of this could be gen quarks
	  ///(from "genParticles" with |pdgId| < 7), as lowLevel objects, and genJets as higherLevel
	  ///objects.  The matching is done purely by deltaR.  The closest dR match wins!
	  ///
	  ///This producer adds a collection of reco::Candidate objects to each event.
	  ///there could be more than one object in either of these collections per event
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > lowLevelToken;
	  edm::EDGetTokenT<edm::OwnVector<reco::Candidate> > higherLevelToken;

	  std::string outputCollName;
	  double maxDeltaR;
	  std::string tName;

	  TTree * tree;

	  Int_t runNumber;
	  ULong64_t evtNumber;

	  ///the number of higher lvl objects which could be matched
	  Int_t nHigherLevel;

	  ///deltaR between the lower level object and all higher lvl objects which could be matched 
	  Float_t dR_lowerToHigherLvlObj[NOBJ];
	
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
FindHigherLevelMatchedObject::FindHigherLevelMatchedObject(const edm::ParameterSet& iConfig):
	outputCollName(iConfig.getParameter<std::string>("matchedOutputCollectionName")),
	maxDeltaR(iConfig.getParameter<double>("dRforMatching")),
	tName(iConfig.getParameter<std::string>("treeName"))
{
   
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"matching higher level objects to lower level objects");

   tree->Branch("nHigherLevel",&nHigherLevel,"nHigherLevel/I");
   tree->Branch("dR_lowerToHigherLvlObj",dR_lowerToHigherLvlObj,"dR_lowerToHigherLvlObj[nHigherLevel]/F");
 
   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");


	
   ///register the input collections	
   lowLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("lowLevelCollTag"));
   higherLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("higherLevelCollTag"));

   ///register the collections which are added to the event
   produces<edm::OwnVector<reco::Candidate> >(outputCollName);

   /*	
   momToken = consumes<std::vector<reco::CompositeCandidate>>(iConfig.getParameter<edm::InputTag>("zedLabel"));
   momParentOneToken = consumes<std::vector<reco::RecoEcalCandidate>>(iConfig.getParameter<edm::InputTag>("tracklessHltEle"));
   momParentTwoToken = consumes<std::vector<reco::RecoEcalCandidate>>(iConfig.getParameter<edm::InputTag>("trackedHltEle"));
   //daughterOneCollection = iConfig.getParameter<std::string>("tracklessEleCollectionName");
   //daughterTwoCollection = iConfig.getParameter<std::string>("trackedEleCollectionName");
   
   //register the two collections of products - std::vector<edm::Refs to RecoEcalCandidate objects> (daughters of Z)
   produces<reco::RecoEcalCandidateRefVector>(daughterOneCollection);
   produces<reco::RecoEcalCandidateRefVector>(daughterTwoCollection);
   */

}


FindHigherLevelMatchedObject::~FindHigherLevelMatchedObject()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FindHigherLevelMatchedObject::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef DEBUG
   std::cout<<"entered FindHigherLevelMatch produce method"<<std::endl;
#endif

   ///get the evt and run number, and initialize the other vars saved in the tree
   evtNumber = iEvent.id().event();
   runNumber = iEvent.id().run();
   nHigherLevel = 0;
   for(Int_t i=0; i<NOBJ; i++){
	   dR_lowerToHigherLvlObj[i] = -1;
   }

   Handle<edm::OwnVector<reco::Candidate> > lowLevelObjectColl;
   iEvent.getByToken(lowLevelToken, lowLevelObjectColl);

   Handle<edm::OwnVector<reco::Candidate> > higherLevelObjectColl;
   iEvent.getByToken(higherLevelToken, higherLevelObjectColl);

#ifdef DEBUG
   std::cout<<"made handles to input collections"<<std::endl;
#endif

   ///make empty an output collection, and a pointer to this collection
   std::auto_ptr<edm::OwnVector<reco::Candidate> > outputObjColl(new edm::OwnVector<reco::Candidate>());

   /**now look for objects in higherLevelObjectColl which are within a distance maxDeltaR
	 *away from any object in lowerLevelObjectColl
	 */
   for(edm::OwnVector<reco::Candidate>::const_iterator lowIt=lowLevelObjectColl->begin(); lowIt!=lowLevelObjectColl->end();
		   lowIt++){
	   for(edm::OwnVector<reco::Candidate>::const_iterator higherIt=higherLevelObjectColl->begin(); higherIt!=higherLevelObjectColl->end();
			   higherIt++){
		   double dR = deltaR(higherIt->eta(),higherIt->phi(),lowIt->eta(),lowIt->phi());
		   if(isNotDuplicate(higherIt,outputObjColl)){
			   dR_lowerToHigherLvlObj[nHigherLevel] = dR;
			   nHigherLevel++;
		   }///end filter to check if the higher level object already exists in the output object collection 
		   if(dR <= maxDeltaR && isNotDuplicate(higherIt,outputObjColl)){
			   outputObjColl->push_back(*higherIt);
		   }///end if(dR cut is passed && not duplicate entry)

	   }///end loop over reco::Candidate objects in lowLevelObjectColl

   }///end loop over reco::Candidate objects in higherLevelObjectColl

#ifdef DEBUG
   std::cout<<"about to put collection of matched reco::Candidate objects into root file"<<std::endl;
#endif
  
   ///fill the tree's branches
   tree->Fill();
   
   ///now put the collection of matched higher level reco::Candidate objects into the event
   iEvent.put(outputObjColl, outputCollName);
 
   /*
   //std::cout<<"entered daughter producer code"<<std::endl;

   //read collection of reco::CompositeCandidate objects from iEvent
   Handle<std::vector<reco::CompositeCandidate> > momIn;
   iEvent.getByToken(momToken, momIn);

   Handle<std::vector<reco::RecoEcalCandidate> > momParentOneIn;	//for trackless leg RECs
   iEvent.getByToken(momParentOneToken, momParentOneIn);

   Handle<std::vector<reco::RecoEcalCandidate> > momParentTwoIn;	//for tracked leg RECs
   iEvent.getByToken(momParentTwoToken, momParentTwoIn);

   //std::cout<<"made handles to input collections"<<std::endl;

   //create empty output collections, one for each daughter, and pointers to each collection
   std::auto_ptr<reco::RecoEcalCandidateRefVector> daughterOneRefColl(new reco::RecoEcalCandidateRefVector );	//trackless collection
   std::auto_ptr<reco::RecoEcalCandidateRefVector> daughterTwoRefColl(new reco::RecoEcalCandidateRefVector );	//tracked collection


   
   for(std::vector<reco::CompositeCandidate>::const_iterator momIt = momIn->begin(); momIt != momIn->end(); momIt++){
	   //get a Ref to a daughter via momIt->daughter()->masterClone()
	   //then find the matching (pt, eta, phi) object in momParent(One or Two)In handles, and
	   //save a reference to the object into the appropriate output collection via
	   //getRef(momParent(One or Two)In, index number) 
	   //std::cout<<"looping over CompositeCandidate objects"<<std::endl;
	   if((momIt->daughter("tracklessRecoEle"))->hasMasterClone() ){
		   //std::cout<<"found tracklessRecoEle daughter with a master clone"<<std::endl;
		   reco::CandidateBaseRef dauOneRef = (momIt->daughter("tracklessRecoEle"))->masterClone();
		   //std::cout<<"made a reference obj to a trackless daughter"<<std::endl;
		   for(unsigned int h=0; h<momParentOneIn->size(); h++){
			   if(dauOneRef->pt() == (getRef(momParentOneIn, h))->pt() ){
				   if(dauOneRef->eta() == (getRef(momParentOneIn, h))->eta() ){
					   if(dauOneRef->phi() == (getRef(momParentOneIn, h))->phi() ){
						   daughterOneRefColl->push_back( getRef(momParentOneIn, h) );
					   }//end filter on phi

				   }//end filter on eta

			   }//end filter on pt

		   }//end loop over objects in momParentOneIn

	   }//end requirement that a master clone exists

	   if((momIt->daughter("trackedRecoEle"))->hasMasterClone() ){
		   //std::cout<<"found trackedRecoEle daughter with a master clone"<<std::endl;
		   reco::CandidateBaseRef dauTwoRef = (momIt->daughter("trackedRecoEle"))->masterClone();
		   for(unsigned int m=0; m<momParentTwoIn->size(); m++){
			   if(dauTwoRef->pt() == (getRef(momParentTwoIn, m))->pt() ){
				   if(dauTwoRef->eta() == (getRef(momParentTwoIn, m))->eta() ){
					   if(dauTwoRef->phi() == (getRef(momParentTwoIn, m))->phi() ){
						   daughterTwoRefColl->push_back( getRef(momParentTwoIn, m) );
					   }//end filter on phi

				   }//end filter on eta

			   }//end filter on pt

		   }//end loop over objects in momParentTwoIn

	   }//end requirement that a master clone exists
	
   }//end loop over all CompositeCandidate objects in the event


   //std::cout<<"about to put daughter collections into root file"<<std::endl;
   //now put the two collections of Refs to daughter particles into the event
   iEvent.put(daughterOneRefColl, daughterOneCollection);
   iEvent.put(daughterTwoRefColl, daughterTwoCollection);

   */


}

// ------------ method called once each job just before starting event loop  ------------
void 
FindHigherLevelMatchedObject::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FindHigherLevelMatchedObject::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
FindHigherLevelMatchedObject::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
FindHigherLevelMatchedObject::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
FindHigherLevelMatchedObject::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
FindHigherLevelMatchedObject::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FindHigherLevelMatchedObject::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FindHigherLevelMatchedObject);
