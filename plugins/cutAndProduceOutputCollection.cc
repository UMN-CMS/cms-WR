// -*- C++ -*-
//
// 

/**\class cutAndProduceOutputCollection cutAndProduceOutputCollection.cc doubleElectronTracklessTrigger/cutAndProduceOutputCollection/plugins/cutAndProduceOutputCollection.cc
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
#include "FWCore/Utilities/interface/InputTag.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/Common/interface/OwnVector.h"



#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"


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

/// include files made by me
#include "../interface/CutApplicationInfrastructure.h"

//#define DEBUG

//
// class declaration
//

class cutAndProduceOutputCollection : public edm::EDProducer {
   public:
      explicit cutAndProduceOutputCollection(const edm::ParameterSet&);
      ~cutAndProduceOutputCollection();

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

	  void resetCounters(){
		  ///initialize all variables stored in the output tree
		  nMatchedHigherLevel=0;
		  nWithMatch=0;
	  }///end resetCounters()

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 
      // ----------member data ---------------------------
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
	  std::string tName;

	  ///initCutString must be of the form "float,< or > character,char array"
	  std::string initCutString;	///< use this to initialize the name, threshold val, and upper bound status of SlimCutVar object

	  const std::string inputTagStrings;
	  unsigned int numCollections;	///< the number of input object collections

	  SlimCutVar theCut;
	
	  TTree * tree;
	  Int_t runNumber;
	  ULong64_t evtNumber;

	  /**
	   * nMatchedHigherLevel counts the number of higher level objects which are matched to lower level objects
	   * nWithMatch counts the number of lower lvl objects for which a match is found
	   * if there are N lower level objects, then the max value of nWithMatch is N, and the max value of
	   * nMatchedHigherLevel can be greater than N
	   */
	  Int_t nMatchedHigherLevel;
	  Int_t nWithMatch;
	
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
cutAndProduceOutputCollection::cutAndProduceOutputCollection(const edm::ParameterSet& iConfig):
	outputCollName(iConfig.getParameter<std::string>("outputCollectionName")),
	tName(iConfig.getParameter<std::string>("treeName")),
	initCutString(iConfig.getParameter<std::string>("initializeCut")),
	inputTagStrings(iConfig.getParameter<std::string>("listOfInputCollTagStrings")),
	numCollections(iConfig.getParameter<unsigned int>("nCollections"))
{

   theCut.setAttributesFromString(initCutString);
   std::cout<<"SlimCutVar object member vars are initialized to:\t"<< theCut << std::endl;
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"applying cut and saving reco::Candidate objects which pass cut");

   tree->Branch("nMatchedHigherLevel",&nMatchedHigherLevel,"nMatchedHigherLevel/I");
   tree->Branch("nWithMatch",&nWithMatch,"nWithMatch/I");
 
   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

	
   ///register the input collections	
   lowLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("lowLevelCollTag"));
   //higherLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("higherLevelCollTag"));
   //higherLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("listOfInputCollTags[0]"));
   higherLevelToken = consumes<edm::OwnVector<reco::Candidate> >( *(new edm::InputTag(inputTagStrings)) );




   ///register the collections which are added to the event
   produces<edm::OwnVector<reco::Candidate> >(outputCollName);

}


cutAndProduceOutputCollection::~cutAndProduceOutputCollection()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
cutAndProduceOutputCollection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

#ifdef DEBUG
   std::cout<<"entered cutAndProduceOutputCollection produce method"<<std::endl;
#endif

   ///get the evt and run number, and initialize the other vars saved in the tree
   evtNumber = iEvent.id().event();
   runNumber = iEvent.id().run();
   resetCounters();

   Handle<edm::OwnVector<reco::Candidate> > lowLevelObjectColl;
   iEvent.getByToken(lowLevelToken, lowLevelObjectColl);

   Handle<edm::OwnVector<reco::Candidate> > higherLevelObjectColl;
   iEvent.getByToken(higherLevelToken, higherLevelObjectColl);

#ifdef DEBUG
   std::cout<<"made handles to input collections"<<std::endl;
#endif

   ///make empty an output collection, and a pointer to this collection
   std::auto_ptr<edm::OwnVector<reco::Candidate> > outputObjColl(new edm::OwnVector<reco::Candidate>());

   /**
	* look for objects
	*/
   for(edm::OwnVector<reco::Candidate>::const_iterator lowIt=lowLevelObjectColl->begin(); lowIt!=lowLevelObjectColl->end();
		   lowIt++){
	   bool incrementedWithMatch = false;
	   for(edm::OwnVector<reco::Candidate>::const_iterator higherIt=higherLevelObjectColl->begin(); higherIt!=higherLevelObjectColl->end();
			   higherIt++){
		   if(isNotDuplicate(higherIt,outputObjColl)){
			   outputObjColl->push_back(*higherIt);
			   nMatchedHigherLevel++;
			   if(!incrementedWithMatch){
				   nWithMatch++;
				   incrementedWithMatch = true;
			   }
		   }///end if(cut is passed && not duplicate entry)

	   }///end loop over reco::Candidate objects in higherLevelObjectColl
	  
   }///end loop over reco::Candidate objects in lowLevelObjectColl
 
   ///fill the tree branches
   tree->Fill();

#ifdef DEBUG
   std::cout<<"about to put collection of reco::Candidate objects into root file"<<std::endl;
#endif
  
   ///now put the collection of reco::Candidate objects which pass the cut into the event
   iEvent.put(outputObjColl, outputCollName);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
cutAndProduceOutputCollection::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
cutAndProduceOutputCollection::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
cutAndProduceOutputCollection::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
cutAndProduceOutputCollection::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
cutAndProduceOutputCollection::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
cutAndProduceOutputCollection::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
cutAndProduceOutputCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(cutAndProduceOutputCollection);
