// -*- C++ -*-
//
// Package:    HeavyNuGenFilter
// Class:      HeavyNuGenFilter
// 
/**\class HeavyNuGenFilter HeavyNuGenFilter.cc HeavyNu/HeavyNuGenFilter/src/HeavyNuGenFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bryan Dahmes
//         Created:  Wed Sep 22 04:49:56 CDT 2010
// $Id: HeavyNuGenFilter.cc,v 1.2 2010/11/08 10:31:52 bdahmes Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

//
// class declaration
//

class HeavyNuGenFilter : public edm::EDFilter {
public:
  explicit HeavyNuGenFilter(const edm::ParameterSet&);
  ~HeavyNuGenFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  std::vector<int> keepIds_;


  // Counters
  std::vector<int> idCounters_;

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
HeavyNuGenFilter::HeavyNuGenFilter(const edm::ParameterSet& iConfig)
  : idCounters_(1000,0) {
  keepIds_ = iConfig.getParameter< std::vector<int32_t> >("keepIds");
}


HeavyNuGenFilter::~HeavyNuGenFilter()
{

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HeavyNuGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  HeavyNuEvent hnuEvent;

  hnuEvent.mc_class=0;

  edm::Handle<reco::GenParticleCollection> genInfo;
  if (iEvent.getByLabel("genParticles",genInfo)) {
    hnuEvent.decayID(*genInfo);
  } 
  
  if (hnuEvent.mc_class<1000) idCounters_[hnuEvent.mc_class]++;
  for (std::vector<int>::const_iterator i=keepIds_.begin(); i!=keepIds_.end(); i++) {
    if (*i==hnuEvent.mc_class) return true;
  }
  return false;

}

// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNuGenFilter::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNuGenFilter::endJob() {
  std::cout << "Results of filtering: " << std::endl ; 
  int n=0;
  for (std::vector<int>::const_iterator i = idCounters_.begin(); i!=idCounters_.end(); i++) {
    if (*i!=0) std::cout << "Saw " << *i << " event of class " << n << std::endl;
    n++;
  }
}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuGenFilter);
