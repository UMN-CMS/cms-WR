// -*- C++ -*-
//
// Package:    doubleElectronTracklessTrigger/factAndRenormWeightAnalyzer
// Class:      factAndRenormWeightAnalyzer
//
/**\class factAndRenormWeightAnalyzer factAndRenormWeightAnalyzer.cc doubleElectronTracklessTrigger/factAndRenormWeightAnalyzer/plugins/factAndRenormWeightAnalyzer.cc

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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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
//#define NPATHS 30

//
// class declaration
//

class factAndRenormWeightAnalyzer : public edm::EDAnalyzer
{
public:
	explicit factAndRenormWeightAnalyzer(const edm::ParameterSet&);
	~factAndRenormWeightAnalyzer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;


// ----------member data ---------------------------

///Handles to RECO object collections
	edm::Handle<LHEEventProduct> lheHandle;

///tokens to input collections
	edm::EDGetTokenT<LHEEventProduct> lheToken;

///input the tree name
	std::string tName;

	//TTree * tree;

///save these variables to the TTree
	//Int_t runNumber;
	//ULong64_t evtNumber;

	//Int_t numInterestingPaths;
	//Float_t numPassingInterestingPaths[NPATHS];
	//TString namesOfInterestingPaths[NPATHS];


};

//
// constants, enums and typedefs
//
typedef LHEEventProduct::WGT WGT;	//probably dont need this, included in a header file

//
// static data member definitions
//

//
// constructors and destructor
//

factAndRenormWeightAnalyzer::factAndRenormWeightAnalyzer(const edm::ParameterSet& iConfig):
	tName(iConfig.getParameter<std::string>("treeName"))

{
	//now do what ever initialization is needed
	//edm::Service<TFileService> fs;

	//tree = fs->make<TTree>(tName.c_str(), "factorization and renorm scale weights");

	//tree->Branch("evtNumber", &evtNumber, "evtNumber/l");
	//tree->Branch("runNumber", &runNumber, "runNumber/I");

	//tree->Branch("numInterestingPaths", &numInterestingPaths, "numInterestingPaths/I");
	//tree->Branch("numPassingInterestingPaths", numPassingInterestingPaths, "numPassingInterestingPaths[numInterestingPaths]/F");
	//tree->Branch("namesOfInterestingPaths", namesOfInterestingPaths, "namesOfInterestingPaths[numInterestingPaths]/C");

	///get the tokens
	lheToken = mayConsume<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

}


factAndRenormWeightAnalyzer::~factAndRenormWeightAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
factAndRenormWeightAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	iEvent.getByToken(lheToken, lheHandle);

#ifdef DEBUG
	std::cout<<"\t"<<std::endl;
	std::cout<<"\t"<<std::endl;
	std::cout << "in analyze method of factAndRenormWeightAnalyzer class" << std::endl;
#endif

	if(lheHandle.isValid() && lheHandle->weights().size() >= 9){
		std::vector<WGT> weightVals = lheHandle->weights();
		//std::vector<std::string> weightComments = lheHandle->getComments();

#ifdef DEBUG
		std::cout << "lheHandle is valid, and there are at least 9 different weights" << std::endl;
		std::cout<<"\t"<<std::endl;
#endif

		unsigned int elems = lheHandle->comments_size();
		for(unsigned int i = 0; i < elems; i++) {
			std::cout << "weight number  " << i << "  =  " << weightVals[i].wgt <<"  comment  "<< lheHandle->getComment(i) << std::endl;
		}

		unsigned int elemsWgts = weightVals.size();
		for(unsigned int i = 0; i < elemsWgts; i++) {
			std::cout << "weight number  " << i << "  =  " << weightVals[i].wgt << std::endl;
		}


	}//end if handle is valid and has at least 9 different weights

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
factAndRenormWeightAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void
factAndRenormWeightAnalyzer::endJob()
{

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
factAndRenormWeightAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(factAndRenormWeightAnalyzer);
