// -*- C++ -*-
//
// Package:    testFilter/triggerFilter
// Class:      triggerFilter
// 
/**\class triggerFilter triggerFilter.cc testFilter/triggerFilter/plugins/triggerFilter.cc

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

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Common/interface/TriggerNames.h"


//#define DEBUG
//
// class declaration
//

class triggerFilter : public edm::EDFilter {
   public:
      explicit triggerFilter(const edm::ParameterSet&);
      ~triggerFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
std::string targetHltPathName;	///< check that this HLT path fired in the evt
std::string targetHltPathNameTwo;	///< check that this HLT path fired in the evt


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
triggerFilter::triggerFilter(const edm::ParameterSet& iConfig):
	targetHltPathName(iConfig.getParameter<std::string>("checkThisHltPath")),
	targetHltPathNameTwo(iConfig.getParameter<std::string>("alsoCheckThisHltPath"))

{
   //now do what ever initialization is needed

}


triggerFilter::~triggerFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
triggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

#ifdef DEBUG
	std::cout<<"in filter method of triggerFilter class"<<std::endl;
#endif

	edm::TriggerResultsByName resultsByName = iEvent.triggerResultsByName("HLT");

	///check that the relevant trigger was fired
	if(!resultsByName.accept(targetHltPathName) && !resultsByName.accept(targetHltPathNameTwo) ) return false;

	return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
triggerFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
triggerFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
triggerFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
triggerFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
triggerFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
triggerFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
triggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(triggerFilter);
