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

#define NOBJ 500
//#define DEBUG

//
// class declaration
//

class applyLeptonJetDrCut : public edm::EDProducer {
   public:
      explicit applyLeptonJetDrCut(const edm::ParameterSet&);
      ~applyLeptonJetDrCut();

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
		  nHigherLevel = 0;
		  nLowerLevel = 0;
		  nMatchedHigherLevel=0;
		  nWithMatch=0;
	
		  for(Int_t i=0; i<NOBJ; i++){
			  dR_lowerToHigherLvlObj[i] = -1;
			  tightest_dR_lowerToMatchedHigherLvlObj[i] = 10;
			  ptLowerLevel[i] = -1;
			  etaLowerLevel[i] = -20;
			  phiLowerLevel[i] = -20;
			  ptHigherLevel[i] = -1;
			  etaHigherLevel[i] = -20;
			  phiHigherLevel[i] = -20;
		  }
	
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
	  double maxDeltaR;
	  std::string tName;

	  TTree * tree;

	  Int_t runNumber;
	  ULong64_t evtNumber;

	  ///the number of higher lvl objects which could be matched
	  Int_t nHigherLevel;

	  ///deltaR between the lower level object and all higher lvl objects which could be matched 
	  Float_t dR_lowerToHigherLvlObj[NOBJ];

	  ///pT, eta, and phi of all higher level objects
	  Float_t ptHigherLevel[NOBJ];
	  Float_t etaHigherLevel[NOBJ];
	  Float_t phiHigherLevel[NOBJ];

	  ///pT, eta, and phi of the lower level objects, and the number of lower lvl objects in the evt
	  ///this info will be useful to study the matching efficiency as a fxn of pT, eta, and phi
	  Int_t nLowerLevel;
	  Float_t ptLowerLevel[NOBJ];
	  Float_t etaLowerLevel[NOBJ];
	  Float_t phiLowerLevel[NOBJ];

	  ///smallest dR btwn each lower lvl object and all higher lvl objects which are successfully matched to lower lvl objects
	  Float_t tightest_dR_lowerToMatchedHigherLvlObj[NOBJ];

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
applyLeptonJetDrCut::applyLeptonJetDrCut(const edm::ParameterSet& iConfig):
	outputCollName(iConfig.getParameter<std::string>("matchedOutputCollectionName")),
	maxDeltaR(iConfig.getParameter<double>("dRforMatching")),
	tName(iConfig.getParameter<std::string>("treeName"))
{
   
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName.c_str(),"matching higher level objects to lower level objects");

   tree->Branch("nHigherLevel",&nHigherLevel,"nHigherLevel/I");
   tree->Branch("dR_lowerToHigherLvlObj",dR_lowerToHigherLvlObj,"dR_lowerToHigherLvlObj[nHigherLevel]/F");
   tree->Branch("ptHigherLevel",ptHigherLevel,"ptHigherLevel[nHigherLevel]/F");
   tree->Branch("etaHigherLevel",etaHigherLevel,"etaHigherLevel[nHigherLevel]/F");
   tree->Branch("phiHigherLevel",phiHigherLevel,"phiHigherLevel[nHigherLevel]/F");

   tree->Branch("nLowerLevel",&nLowerLevel,"nLowerLevel/I");
   tree->Branch("ptLowerLevel",ptLowerLevel,"ptLowerLevel[nLowerLevel]/F");
   tree->Branch("etaLowerLevel",etaLowerLevel,"etaLowerLevel[nLowerLevel]/F");
   tree->Branch("phiLowerLevel",phiLowerLevel,"phiLowerLevel[nLowerLevel]/F");
   tree->Branch("tightest_dR_lowerToMatchedHigherLvlObj",tightest_dR_lowerToMatchedHigherLvlObj,"tightest_dR_lowerToMatchedHigherLvlObj[nLowerLevel]/F");
 
   tree->Branch("nMatchedHigherLevel",&nMatchedHigherLevel,"nMatchedHigherLevel/I");
   tree->Branch("nWithMatch",&nWithMatch,"nWithMatch/I");
 
   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

	
   ///register the input collections	
   lowLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("lowLevelCollTag"));
   higherLevelToken = consumes<edm::OwnVector<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("higherLevelCollTag"));

   ///register the collections which are added to the event
   produces<edm::OwnVector<reco::Candidate> >(outputCollName);

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
   std::cout<<"entered FindHigherLevelMatch produce method"<<std::endl;
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

   /**now look for objects in higherLevelObjectColl which are within a distance maxDeltaR
	 *away from any object in lowerLevelObjectColl
	 */
   for(edm::OwnVector<reco::Candidate>::const_iterator lowIt=lowLevelObjectColl->begin(); lowIt!=lowLevelObjectColl->end();
		   lowIt++){
	   bool incrementedWithMatch = false;
	   for(edm::OwnVector<reco::Candidate>::const_iterator higherIt=higherLevelObjectColl->begin(); higherIt!=higherLevelObjectColl->end();
			   higherIt++){
		   double dR = deltaR(higherIt->eta(),higherIt->phi(),lowIt->eta(),lowIt->phi());
		   if(isNotDuplicate(higherIt,outputObjColl)){
			   dR_lowerToHigherLvlObj[nHigherLevel] = dR;
			   ptHigherLevel[nHigherLevel] = higherIt->pt();
			   etaHigherLevel[nHigherLevel] = higherIt->eta();
			   phiHigherLevel[nHigherLevel] = higherIt->phi();
			   nHigherLevel++;
		   }///end filter to check if the higher level object already exists in the output object collection 
		   if(dR <= maxDeltaR && isNotDuplicate(higherIt,outputObjColl)){
			   outputObjColl->push_back(*higherIt);
			   nMatchedHigherLevel++;
			   ///if dR is less than the current value in tightest_dR...[nLowerLevel], then put dR into the array and overwrite the old element
			   if(dR < tightest_dR_lowerToMatchedHigherLvlObj[nLowerLevel]) tightest_dR_lowerToMatchedHigherLvlObj[nLowerLevel] = dR;
			   if(!incrementedWithMatch){
				   nWithMatch++;
				   incrementedWithMatch = true;
			   }
		   }///end if(dR cut is passed && not duplicate entry)

	   }///end loop over reco::Candidate objects in higherLevelObjectColl
	  
	   ptLowerLevel[nLowerLevel] = lowIt->pt();
	   etaLowerLevel[nLowerLevel] = lowIt->eta();
	   phiLowerLevel[nLowerLevel] = lowIt->phi();
	   nLowerLevel++;

   }///end loop over reco::Candidate objects in lowLevelObjectColl
 
   ///fill the tree branches
   tree->Fill();

#ifdef DEBUG
   std::cout<<"about to put collection of matched reco::Candidate objects into root file"<<std::endl;
#endif
  
   ///now put the collection of matched higher level reco::Candidate objects into the event
   iEvent.put(outputObjColl, outputCollName);
 
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
