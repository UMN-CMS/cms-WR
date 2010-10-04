// -*- C++ -*-
//
// Package:    MuFilter
// Class:      MuFilter
// 
/**\class MuFilter MuFilter.cc HeavyNu/MuFilter/src/MuFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bryan Dahmes
//         Created:  Wed Sep 22 04:49:56 CDT 2010
// $Id$
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

//
// class declaration
//

class compare {
public:
  template <class T> bool operator() (const T& a, const T& b) { return a.pt() > b.pt() ; } 
};

class MuFilter : public edm::EDFilter {
public:
  explicit MuFilter(const edm::ParameterSet&);
  ~MuFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  double minPt_ ; // minimum pT cut on skimming objects
  std::vector<std::string> triggerList_ ; // List of HLT paths used for filtering

  // Counters
  int nEvt, nMu1, nMu2, nEB, nEE ; 

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
MuFilter::MuFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  minPt_       = iConfig.getParameter<double>("minPt") ;
  triggerList_ = iConfig.getParameter< std::vector<std::string> >("HLTpaths") ; 
}


MuFilter::~MuFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  edm::Handle<pat::MuonCollection> patMuonCollection ; 
  iEvent.getByLabel("patMuons",patMuonCollection) ; 
  if (!patMuonCollection.isValid()) return false ; 
  pat::MuonCollection pMuons = *(patMuonCollection.product()) ;
  
  edm::Handle<reco::SuperClusterCollection> ebSCCollection ; 
  iEvent.getByLabel("correctedHybridSuperClusters",ebSCCollection) ; 
  if (!ebSCCollection.isValid()) return false ; 
  const reco::SuperClusterCollection ebSCs = *(ebSCCollection.product()) ; 

  edm::Handle<reco::SuperClusterCollection> eeSCCollection ; 
  iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower",eeSCCollection) ; 
  if (!eeSCCollection.isValid()) return false ; 
  const reco::SuperClusterCollection eeSCs = *(eeSCCollection.product()) ; 

  // Require at least one muon
  if ( pMuons.size() == 0 ) return false ; 
  nEvt++ ; // Number of events with at least one muon

  // Sort the collections by pT
  std::sort(pMuons.begin(),pMuons.end(),compare()) ; 
  // std::sort(ebSCs.begin(),ebSCs.end(),compare()) ; 
  // std::sort(eeSCs.begin(),eeSCs.end(),compare()) ; 
  
  if ( pMuons.at(0).pt() < minPt_ ) return false ; 
  nMu1++ ; 
  for (unsigned int i=1; i<pMuons.size(); i++) {
    if ( pMuons.at(i).pt() >= minPt_ ) { 
      nMu2++ ; 
      return true ; 
    }
  }
  for (unsigned int i=0; i<ebSCs.size(); i++)
    if ( ebSCs.at(i).energy()/cosh(ebSCs.at(i).eta()) >= minPt_ ) { nEB++ ; return true ; }
  
  for (unsigned int i=0; i<eeSCs.size(); i++) 
    if ( eeSCs.at(i).energy()/cosh(eeSCs.at(i).eta()) >= minPt_ ) { nEE++ ; return true ; }

  return false ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuFilter::beginJob()
{
  nEvt = 0 ; nMu1 = 0 ; nMu2 = 0 ; nEB = 0 ; nEE = 0 ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuFilter::endJob() {
  std::cout << "Results of filtering: " << std::endl ; 
  std::cout << "Saw " << nEvt << " events with at least one muon" << std::endl ; 
  std::cout << "    " << nMu1 << " events one muon above min pT" << std::endl ; 
  std::cout << "    " << nMu2 << " events two muons above min pT" << std::endl ; 
  std::cout << "    " << nEB << " events one muon, EB SC above min pT" << std::endl ; 
  std::cout << "    " << nEE << " events one muon, EE SC above min pT" << std::endl ; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuFilter);
