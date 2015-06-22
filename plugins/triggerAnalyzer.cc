// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/triggerAnalyzer
// Class:      triggerAnalyzer
// 
/**\class triggerAnalyzer triggerAnalyzer.cc doubleElectronTracklessTrigger/triggerAnalyzer/plugins/triggerAnalyzer.cc

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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

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
#define NPATHS 30 

//
// class declaration
//

class triggerAnalyzer : public edm::EDAnalyzer {
   public:
      explicit triggerAnalyzer(const edm::ParameterSet&);
      ~triggerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;


//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

// ----------member data ---------------------------

///Handles to RECO object collections
edm::Handle<edm::TriggerResults> trigResultsHandl;
edm::Handle<pat::TriggerObjectStandAlone> trigObjsHandl;

///tokens to input collections
edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
edm::EDGetTokenT<pat::TriggerObjectStandAlone> trigObjsToken;

///input the tree name, and the names of all HLT paths which should be investigated
std::string tName;
std::string listOfHltPaths;

TTree * tree;

///save these variables to the TTree
Int_t runNumber;
ULong64_t evtNumber;

Int_t numInterestingPaths;
Float_t numPassingInterestingPaths[NPATHS];
TString namesOfInterestingPaths[NPATHS];


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

triggerAnalyzer::triggerAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName")),
	listOfHltPaths(iConfig.getParameter<std::string>("commaSeparatedHltPaths"))

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs;
  
   tree=fs->make<TTree>(tName.c_str(),"HLT path summary");

   tree->Branch("evtNumber",&evtNumber,"evtNumber/l");
   tree->Branch("runNumber",&runNumber,"runNumber/I");

   tree->Branch("numInterestingPaths",&numInterestingPaths,"numInterestingPaths/I");
   tree->Branch("numPassingInterestingPaths",numPassingInterestingPaths,"numPassingInterestingPaths[numInterestingPaths]/F");
   tree->Branch("namesOfInterestingPaths",namesOfInterestingPaths,"namesOfInterestingPaths[numInterestingPaths]/C");


   ///get the tokens
   trigResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigResultsColl"));
   trigObjsToken = consumes<pat::TriggerObjectStandAlone>(iConfig.getParameter<edm::InputTag>("trigObjectStandAloneColl"));
 
}


triggerAnalyzer::~triggerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
triggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	iEvent.getByToken(trigResultsToken, trigResultsHandl);

#ifdef DEBUG
	std::cout<<"in analyze method of triggerAnalyzer class"<<std::endl;
#endif

	if(!trigResultsHandl.isValid()) return;
	const TriggerNames & triggerNames = iEvent.triggerNames(*trigResultsHandl);

#ifdef DEBUG
	std::cout<<"trigResultsHandl is valid, and successfully returned a list of trigger names"<<std::endl;
#endif

	///now triggerNames has the list of all HLT path names which were run in the event (independent of these paths firing or not)
	unsigned int numNames = triggerNames.size();
	for(unsigned int i=0; i<numNames; i++){
		std::cout<<"trigger number \t"<< i <<"\t is named \t"<< triggerNames.triggerName(i) << std::endl;

	}

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
triggerAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
triggerAnalyzer::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
void 
triggerAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
triggerAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
triggerAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
triggerAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triggerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(triggerAnalyzer);
