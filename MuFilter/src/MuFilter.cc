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
// $Id: MuFilter.cc,v 1.3 2011/05/28 18:25:11 bdahmes Exp $
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
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
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
  bool Overlap(double muEta, double muPhi, double xEta, double xPhi) { return ( deltaR(muEta,muPhi,xEta,xPhi) < overlap_ ) ; }
  
  // ----------member data ---------------------------
  double minMuPt_, minEleEt_ ; // minimum Et/pt cut on skimming objects
  edm::InputTag muonTag_, trackTag_, ebTag_, eeTag_ ; 
  double overlap_ ; 

  // Counters
  int nEvt, nMu1, nMu2, nTrk, nEB, nEE ; 

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
  minMuPt_  = iConfig.getParameter<double>("minMuonPt") ;
  minEleEt_ = iConfig.getParameter<double>("minSCEt") ; 

  muonTag_  = iConfig.getParameter<edm::InputTag>("muonTag") ; 
  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag") ; 
  ebTag_    = iConfig.getParameter<edm::InputTag>("ebTag") ; 
  eeTag_    = iConfig.getParameter<edm::InputTag>("eeTag") ; 

  overlap_  = iConfig.getParameter<double>("overlap") ; 
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
  iEvent.getByLabel(muonTag_,patMuonCollection) ; 
  if (!patMuonCollection.isValid()) return false ; 
  pat::MuonCollection pMuons = *(patMuonCollection.product()) ;

  edm::Handle<pat::GenericParticleCollection> patTrackCollection ; 
  iEvent.getByLabel(trackTag_,patTrackCollection) ; 
  if (!patTrackCollection.isValid()) return false ; 
  pat::GenericParticleCollection pTracks = *(patTrackCollection.product()) ;
  
  edm::Handle<reco::SuperClusterCollection> ebSCCollection ; 
  iEvent.getByLabel(ebTag_,ebSCCollection) ; 
  if (!ebSCCollection.isValid()) return false ; 
  const reco::SuperClusterCollection ebSCs = *(ebSCCollection.product()) ; 

  edm::Handle<reco::SuperClusterCollection> eeSCCollection ; 
  iEvent.getByLabel(eeTag_,eeSCCollection) ; 
  if (!eeSCCollection.isValid()) return false ; 
  const reco::SuperClusterCollection eeSCs = *(eeSCCollection.product()) ; 

  // Require at least one muon
  if ( pMuons.size() == 0 ) return false ; 
  nEvt++ ; // Number of events with at least one muon

  // Sort the collections by pT
  std::sort(pMuons.begin(),pMuons.end(),compare()) ; 
  std::sort(pTracks.begin(),pTracks.end(),compare()) ; 

  double muEta = pMuons.at(0).eta() ; 
  double muPhi = pMuons.at(0).phi() ; 
  if ( pMuons.at(0).pt() < minMuPt_ ) return false ; 
  nMu1++ ; 
  for (unsigned int i=1; i<pMuons.size(); i++) {
    if ( pMuons.at(i).pt() >= minMuPt_ ) { nMu2++ ; return true ; }
  }

  // Check for generic tracks with sufficient pT
  for (unsigned int i=0; i<pTracks.size(); i++) 
    if ( pTracks.at(i).pt() >= minMuPt_ &&
	 !Overlap(muEta,muPhi,pTracks.at(i).eta(),pTracks.at(i).phi()) ) { nTrk++ ; return true ; }

  // Check for ECAL clusters
  for (unsigned int i=0; i<ebSCs.size(); i++)
    if ( ebSCs.at(i).energy()/cosh(ebSCs.at(i).eta()) >= minEleEt_ && 
	 !Overlap(muEta,muPhi,ebSCs.at(i).eta(),ebSCs.at(i).phi()) ) { nEB++ ; return true ; }
  for (unsigned int i=0; i<eeSCs.size(); i++) 
    if ( eeSCs.at(i).energy()/cosh(eeSCs.at(i).eta()) >= minEleEt_ && 
	 !Overlap(muEta,muPhi,eeSCs.at(i).eta(),eeSCs.at(i).phi()) ) { nEE++ ; return true ; }

  return false ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuFilter::beginJob()
{
  nEvt = 0 ; nMu1 = 0 ; nMu2 = 0 ; nEB = 0 ; nEE = 0 ; nTrk = 0 ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuFilter::endJob() {
  std::cout << "Results of filtering: " << std::endl ; 
  std::cout << "Saw " << nEvt << " events with at least one muon" << std::endl ; 
  std::cout << "    " << nMu1 << " events with one muon above min pT" << std::endl ; 
  std::cout << "    " << nTrk << " events with one muon, one generic track above min pT" << std::endl ; 
  std::cout << "    " << nMu2 << " events with two muons above min pT" << std::endl ; 
  std::cout << "    " << nEB  << " events with one muon, EB SC above min pT" << std::endl ; 
  std::cout << "    " << nEE  << " events with one muon, EE SC above min pT" << std::endl ; 
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuFilter);
