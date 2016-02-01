// -*- C++ -*-
//
// Package:    cmsWR/determineCutEfficiency
// Class:      determineCutEfficiency
// 
/**\class determineCutEfficiency determineCutEfficiency.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sean Kalafut
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



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


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

class determineCutEfficiency : public edm::EDAnalyzer {
   public:
      explicit determineCutEfficiency(const edm::ParameterSet&);
      ~determineCutEfficiency();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
virtual void beginJob() override;
virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
virtual void endJob() override;



// ----------member data ---------------------------

std::string tName_;	///< tree name, input by user
double tmpWrMass_, tmpNuMass_;	///<placeholders for the values stored in wrMass and nuMass
Float_t totalNumEvtsBeforeCuts = 0, totalNumEvtsAfterCuts = 0;
	
///Handles to evt count collections
edm::Handle<Float_t> evtCountBeforeCuts;
edm::Handle<Float_t> evtCountAfterCuts;

///tokens to input collections
edm::EDGetTokenT<Float_t> evtCountBeforeCutsToken;
edm::EDGetTokenT<Float_t> evtCountAfterCutsToken;

///TTree and its branches
///the branches will be filled in the analyze() method
///and parsed in the endJob() method 
TTree * tree;
Float_t wrMass, nuMass;	///< these will be input by the user, and saved as TTree branches
Float_t cutEfficiency;	///< fraction of evts which pass all cuts

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

determineCutEfficiency::determineCutEfficiency(const edm::ParameterSet& iConfig):
	tName_(iConfig.getParameter<std::string>("treeName")),
	tmpWrMass_(iConfig.getParameter<double>("massOfWR")),
	tmpNuMass_(iConfig.getParameter<double>("massOfNu"))

{
   edm::Service<TFileService> fs;
   
   tree=fs->make<TTree>(tName_.c_str(),"wr and nu mass, and num events before and after all cuts");

   tree->Branch("wrMass",&wrMass,"wrMass/F");
   tree->Branch("nuMass",&nuMass,"nuMass/F");
   //tree->Branch("evtsBeforeCuts",&evtsBeforeCuts,"evtsBeforeCuts/F");
   //tree->Branch("evtsAfterCuts",&evtsAfterCuts,"evtsAfterCuts/F");
   tree->Branch("cutEfficiency",&cutEfficiency,"cutEfficiency/F");
 
   evtCountBeforeCutsToken = consumes<Float_t>(iConfig.getParameter<edm::InputTag>("evtCountBeforeCutsCollection"));
   evtCountAfterCutsToken = consumes<Float_t>(iConfig.getParameter<edm::InputTag>("evtCountAfterCutsCollection"));

}


determineCutEfficiency::~determineCutEfficiency()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
determineCutEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	iEvent.getByToken(evtCountBeforeCutsToken, evtCountBeforeCuts);
	iEvent.getByToken(evtCountAfterCutsToken, evtCountAfterCuts);
	if(evtCountBeforeCuts.isValid() ) totalNumEvtsBeforeCuts+=1.0;
	if(evtCountAfterCuts.isValid() ) totalNumEvtsAfterCuts+=1.0;

}


// ------------ method called once each job just before starting event loop  ------------
void 
determineCutEfficiency::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
determineCutEfficiency::endJob() 
{
#ifdef DEBUG
	std::cout<<"in endJob method"<<std::endl;
#endif

	///now the efficiency can be computed with totalNumEvts before/after all cuts
	wrMass = tmpWrMass_, nuMass = tmpNuMass_;
	cutEfficiency = 0.0;
	if(totalNumEvtsBeforeCuts > 0) cutEfficiency = totalNumEvtsAfterCuts/totalNumEvtsBeforeCuts;
	tree->Fill();	///< this may not work

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
determineCutEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(determineCutEfficiency);
