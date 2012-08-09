// -*- C++ -*-
//
// Package:    MuFilter
// Class:      QCDFilter
// 
/**\class QCDFilter QCDFilter.cc HeavyNu/MuFilter/src/QCDFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bryan Dahmes
//         Created:  Wed Sep 22 04:49:56 CDT 2010
// $Id: QCDFilter.cc,v 1.4 2011/07/28 07:18:07 bdahmes Exp $
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
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.cc"

//
// class declaration
//

class QCDFilter : public edm::EDFilter {
public:
  explicit QCDFilter(const edm::ParameterSet&);
  ~QCDFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  double minMuonPt_, minElectronPt_ ; 
  bool collectMuons_, collectElectrons_ ; 
  int prescale_ ; 
  double cutoff_ ; 

  edm::InputTag muonTag_, electronTag_ ; 

  // Counters
  int nEvt, nMuE1, nMuE2 ; 
  int n40, n50, n60, n80, n100 ; 
  int nLowPt, nLowPtKeep ; 

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
QCDFilter::QCDFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  minMuonPt_     = iConfig.getParameter<double>("minMuonPt") ;
  minElectronPt_ = iConfig.getParameter<double>("minElectronPt") ;

  prescale_ = iConfig.getParameter<int>("lowPtPrescale") ;
  cutoff_   = iConfig.getParameter<double>("lowPtCutoff") ;

  muonTag_     = iConfig.getParameter<edm::InputTag>("muonTag") ; 
  electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag") ; 

  collectMuons_     = iConfig.getParameter<bool>("collectMuons") ; 
  collectElectrons_ = iConfig.getParameter<bool>("collectElectrons") ; 

}


QCDFilter::~QCDFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
QCDFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  if ( collectMuons_ && collectElectrons_ ) 
    throw cms::Exception("Cannot skim muons and electrons simultaneously.  Filter misconfigured.") ;

  edm::Handle<pat::MuonCollection> patMuonCollection ; 
  iEvent.getByLabel(muonTag_,patMuonCollection) ; 
  if (!patMuonCollection.isValid()) { 
    std::cout << "WARNING: Exiting as valid muon collection absent: " << muonTag_ << std::endl ; 
    return false ; 
  }
  pat::MuonCollection pMuons = *(patMuonCollection.product()) ;

  edm::Handle<pat::ElectronCollection> patElectronCollection ; 
  iEvent.getByLabel(electronTag_,patElectronCollection) ; 
  if (!patElectronCollection.isValid()) {      
    std::cout << "WARNING: Exiting as valid electron collection absent: " << electronTag_ << std::endl ; 
    return false ; 
  }
  pat::ElectronCollection pElectrons = *(patElectronCollection.product()) ;

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);

  // Require at least one lepton 
  if ( collectMuons_     &&     pMuons.size() < 1 ) return false ;
  if ( collectElectrons_ && pElectrons.size() < 1 ) return false ;
  nEvt++ ; // Number of events with at least one lepton

  double highestPtLepton = ( collectMuons_ ? minMuonPt_ : minElectronPt_ ) ; 
  std::vector<unsigned int> locs ; 

  if ( collectMuons_ ) { 
    for (unsigned int i=0; i<pMuons.size(); i++) { 
      pat::Muon muonCand = pMuons.at(i) ; 
      if ( muonCand.pt() > minMuonPt_ ) { 
	locs.push_back(i) ; 
	if ( muonCand.pt() > highestPtLepton ) highestPtLepton = muonCand.pt() ; 
      }
    }
  } else { 
    for (unsigned int i=0; i<pElectrons.size(); i++) { 
      pat::Electron electronCand = pElectrons.at(i) ; 

      double eTcorr   = hnu::getElectronEt( electronCand, true ) ; // Using corrected energy
      double eTuncorr = hnu::getElectronEt( electronCand, false ) ; // Using uncorrected energy
      double eT = std::max( eTcorr,eTuncorr ) ; 
      if ( eT > minElectronPt_ ) { 
	locs.push_back(i) ; 
	if ( eT > highestPtLepton ) highestPtLepton = eT ; 
      }
    }
  }

  // Check if any valid candidate found
  if ( locs.size() < 1 ) return false ; 
  nMuE1++ ; 

  // If more than one candidate found, reject event if both high quality leptons
  // Not implemented yet!!
  if ( locs.size() > 1 ) { 
    int nTight = 0 ; 
    for (unsigned int i=0; i<locs.size(); i++) { 
      if ( collectMuons_ ) { 
	if ( hnu::is2012MuTight(pMuons.at(locs.at(i)),pvHandle) ) nTight++ ; 
      } else { 
	// By assuming rho = 0, fewer multi-electron events will be found
	// More conservative for a skim
	if ( hnu::passesHEEP(pElectrons.at(locs.at(i)),40,0.0) ) nTight++ ; 
      }
    }
    if ( nTight > 1 ) nMuE2++ ; 
  }

  // pT counter
  if ( highestPtLepton > 40.0 )  n40++ ; 
  if ( highestPtLepton > 50.0 )  n50++ ; 
  if ( highestPtLepton > 60.0 )  n60++ ; 
  if ( highestPtLepton > 80.0 )  n80++ ; 
  if ( highestPtLepton > 100.0 ) n100++ ; 

  if ( (prescale_ > 1) && (highestPtLepton < cutoff_) ) { // May need to prescale 
    nLowPt++ ; 
    if ( collectMuons_ ) {
      if ( (nLowPt % prescale_) != 0 ) return false ; 
      else                             nLowPtKeep++ ; 
    }
  }

  return true ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
QCDFilter::beginJob()
{
  nEvt = 0 ; nMuE1 = 0 ; nMuE2 = 0 ; 
  n40 = 0 ; n50 = 0 ; n60 = 0 ; n80 = 0 ; n100 = 0 ; 
  nLowPt = 0 ; nLowPtKeep = 0 ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDFilter::endJob() {
  std::cout << "Results of filtering: " << std::endl ; 
  std::cout << "Saw " << nEvt   << " events with at least one " << (collectMuons_ ? "muon" : "electron") << std::endl ; 
  std::cout << "    " << nMuE1  << " events with one " << (collectMuons_ ? "muon" : "electron") 
	    << " above " << (collectMuons_ ? minMuonPt_ : minElectronPt_) << " GeV" << std::endl ; 
  std::cout << "    " << nMuE2  << " events with two tight " << (collectMuons_ ? "muons" : "electrons") 
	    << " above " << (collectMuons_ ? minMuonPt_ : minElectronPt_) << " GeV" << std::endl ; 
  std::cout << "    " << n40  << " events passing with one " << (collectMuons_ ? "muon" : "electron") << ", pT > 40 GeV"  << std::endl ; 
  std::cout << "    " << n50  << " events passing with one " << (collectMuons_ ? "muon" : "electron") << ", pT > 50 GeV"  << std::endl ; 
  std::cout << "    " << n60  << " events passing with one " << (collectMuons_ ? "muon" : "electron") << ", pT > 60 GeV"  << std::endl ; 
  std::cout << "    " << n80  << " events passing with one " << (collectMuons_ ? "muon" : "electron") << ", pT > 80 GeV"  << std::endl ; 
  std::cout << "    " << n100 << " events passing with one " << (collectMuons_ ? "muon" : "electron") << ", pT > 100 GeV"  << std::endl ; 

  if ( collectMuons_ ) { 
    std::cout << "    " << nLowPt     << " events with one muon below cutoff of " << cutoff_ << " GeV" << std::endl ; 
    std::cout << "    " << nLowPtKeep << " events passing with one muon below cutoff of " << cutoff_ << " GeV" << std::endl ; 
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDFilter);
