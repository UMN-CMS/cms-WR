// -*- C++ -*-
//
// Package:    EleFilter
// Class:      EleFilter
// 
/**\class EleFilter EleFilter.cc HeavyNu/EleFilter/src/EleFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Franzoni
//         Created:  Wed May  2 04:49:56 CDT 2012
// $Id: EleFilter.cc,v 1.4 2012/05/02 07:18:07 franzoni Exp $
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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

//
// class declaration
//

class compareEle {
public:
  template <class T> bool operator() (const T& a, const T& b) { return (  a.pt() > b.pt()  ); } 
};

bool compareSC ( const reco::SuperCluster a,  const reco::SuperCluster b ) { return a.energy()/cosh(a.eta()) > b.energy()/cosh(b.eta()) ;  }


class EleFilter : public edm::EDFilter {
public:
  explicit EleFilter(const edm::ParameterSet&);
  ~EleFilter();
  
private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool Overlap(double eleEta, double elePhi, double xEta, double xPhi) { return ( deltaR(eleEta,elePhi,xEta,xPhi) < overlap_ ) ; }
  
  // ----------member data ---------------------------
  double       minElePt_;
  unsigned int minEleNum_;

  double       minSCEt_;
  double       maxSCAbsEta_;
  unsigned int minSCNum_;

  double       massMin_;
  double       massMax_;

  bool         debug_;

  edm::InputTag electronTag_, trackTag_, ebTag_, eeTag_ ; 
  double overlap_ ; 

  // Counters
  int nEvt, nEvt1, nEvt2, nEvt3;

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
EleFilter::EleFilter(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed

  electronTag_  = iConfig.getParameter<edm::InputTag>("electronTag") ; 

  ebTag_        = iConfig.getParameter<edm::InputTag>("ebTag") ; 
  eeTag_        = iConfig.getParameter<edm::InputTag>("eeTag") ; 

  minElePt_     = iConfig.getParameter<double>("minElePt") ;
  minEleNum_    = iConfig.getParameter<double>("minEleNum") ;
	       
  minSCEt_      = iConfig.getParameter<double>("minSCEt") ;
  maxSCAbsEta_  = iConfig.getParameter<double>("maxSCAbsEta") ;
  minSCNum_     = iConfig.getParameter<double>("minSCNum") ;

  overlap_      = iConfig.getParameter<double>("overlap") ; 

  massMin_      = iConfig.getParameter<double>("massMin") ;
  massMax_      = iConfig.getParameter<double>("massMax") ; 

  debug_        = iConfig.getParameter<bool>("debug") ;

  nEvt=nEvt1=nEvt2=nEvt3=0;
}

EleFilter::~EleFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<pat::ElectronCollection> patElectronCollection ; 
  iEvent.getByLabel(electronTag_,patElectronCollection) ; 
  if (!patElectronCollection.isValid()) {      
    std::cout << "no pat electron collection: " << electronTag_ << " found; bailing out." << std::endl;         assert(0);  } 
  pat::ElectronCollection pElectrons = *(patElectronCollection.product()) ;

  edm::Handle<reco::SuperClusterCollection> ebSCCollection ; 
  iEvent.getByLabel(ebTag_,ebSCCollection) ; 
  if (!ebSCCollection.isValid())     {        std::cout << "no valid collection EB SC: " << ebTag_ << " found; bailing out." << std::endl;         assert(0);  } 
  else {  if(debug_) std::cout << "++ debug: found collection EB SC: " << ebTag_ << ", with " << (*(ebSCCollection.product())).size() << " elements" << std::endl;     }
  const reco::SuperClusterCollection ebSCs = *(ebSCCollection.product()) ; 
  
  edm::Handle<reco::SuperClusterCollection> eeSCCollection ; 
  iEvent.getByLabel(eeTag_,eeSCCollection) ; 
  if (!eeSCCollection.isValid()) {        std::cout << "no valid collection EE SC: " << eeTag_ << " found; bailing out." << std::endl;         assert(0);  } 
  else { if(debug_) std::cout << "++ debug: found collection EE SC: " << eeTag_ << ", with " << (*(eeSCCollection.product())).size() << " elements" << std::endl;     }
  const reco::SuperClusterCollection eeSCs = *(eeSCCollection.product()) ; 
  
  
  // Require at least minEleNum_ electron
  if ( pElectrons.size() < minEleNum_ || pElectrons.size()==0 )     {
    if(debug_) std::cout << "++ debug: no gsf found, or less than: " << minEleNum_ << " ( " << pElectrons.size() << " ) ; returning false " << std::endl;
    return false ; 
  }
  else { if(debug_) std::cout<< "++ debug: found: " << pElectrons.size() << " gsf; continuing (more or equal to minum required: " << minEleNum_ << ")" << std::endl;}
  

  // Sort the collections by pT
  std::sort(pElectrons.begin(),pElectrons.end(),compareEle()) ; 
  
  // require at least minEleNum_ electron above minElePt_
  std::vector<pat::Electron> patSelectedElectrons;
  patSelectedElectrons.reserve( pElectrons.size()  );
  for(unsigned int el=0; el<pElectrons.size(); el++){
    if ( pElectrons.at(el).pt() < minElePt_ ) continue;
    patSelectedElectrons.push_back( pElectrons.at(el) );   
  }

  // if not enough selected electrons, fail
  if( patSelectedElectrons.size() <  minEleNum_){
    if(debug_) std::cout << "++ debug: less than: " << minEleNum_ << " electrons found (" << patSelectedElectrons.size() << ") with pt > " << minElePt_ << "; returning false " << std::endl;
    return false;
  }
  else {
    if(debug_) std::cout << "++ debug " << patSelectedElectrons.size() << " electrons found with pt above: " << minElePt_ << "; minimum num gsf satisfied. " << std::endl;
  }

  nEvt1++; // count how many have required electrons

  // if electrons cover nuber of required electrons AND nuber of required SC's as well, no need to check SC's => pass
  if( patSelectedElectrons.size() >= (minEleNum_ + minSCNum_) ) {
    std::cout << "number of electrons (" << patSelectedElectrons.size() << ") exceeds minEle+minSC " << (minEleNum_ + minSCNum_) << " returning true" << std::endl;
    return true;
  }

  // find SC's that pass et/absEta requirements, but are not yet associated to an electron
  std::vector<reco::SuperCluster> selectedSC;
  selectedSC.reserve( eeSCs.size() );

  // EE sc's
  for(unsigned int eesc=0; eesc<eeSCs.size(); eesc++ ){

    if ( eeSCs.at(eesc).energy()/cosh(eeSCs.at(eesc).eta()) < minSCEt_ ) continue;
    if ( fabs( eeSCs.at(eesc).eta() ) > maxSCAbsEta_                   ) continue;
    
    // make sure this SC is not already among the selected electrons
    for(unsigned int sEl=0; sEl<patSelectedElectrons.size(); sEl++) {
      if (! Overlap(eeSCs.at(eesc).eta(), eeSCs.at(eesc).phi(), patSelectedElectrons.at(sEl).eta(), patSelectedElectrons.at(sEl).phi() )  ){
	selectedSC.push_back( eeSCs.at(eesc) );
      }
    }// loop over selected electrons

  }// loop over EE SC
  if(debug_) std::cout << "++ debug: after loop on EE, found: " << selectedSC.size() << " SC's" << std::endl;

  // EB sc's
  for(unsigned int ebsc=0; ebsc<ebSCs.size(); ebsc++ ){

    if ( ebSCs.at(ebsc).energy()/cosh(ebSCs.at(ebsc).eta()) < minSCEt_ ) continue;
    if ( fabs( ebSCs.at(ebsc).eta() ) > maxSCAbsEta_                   ) continue;
    
    // make sure this SC is not already among the selected electrons
    for(unsigned int sEl=0; sEl<patSelectedElectrons.size(); sEl++) {
      if (! Overlap(ebSCs.at(ebsc).eta(), ebSCs.at(ebsc).phi(), patSelectedElectrons.at(sEl).eta(), patSelectedElectrons.at(sEl).phi() )  ){
	selectedSC.push_back( ebSCs.at(ebsc) );
      }
    }// loop over selected electrons

  }// loop over EB SC
  if(debug_)  std::cout << "++ debug: after loop on EB, found: " << selectedSC.size() << " SC's" << std::endl;

  if( selectedSC.size() < minSCNum_) {
    if(debug_) std::cout << "++debug: num of selected SCs extra2electros ( " << selectedSC.size() << ") less than minimum (" << minSCNum_ << ") returning false;" << std::endl;
    return false;
  }
  
  nEvt2++; // count how many have required SC's

  // at this stage, that patSelectedElectrons and selectedSC have at least one element has already been checked
  // double-check, for protection, to avoid invalid .at(0)
  if(patSelectedElectrons.size()==0 || selectedSC.size()==0) return false;

  // Sort the collections by pT
  std::sort(patSelectedElectrons.begin(),patSelectedElectrons.end(),compareEle()) ; 
  std::sort(selectedSC.begin(),          selectedSC.end(),          compareSC) ; 

  // further selection on mass
  reco::Particle::LorentzVector Z(patSelectedElectrons.at(0).ecalDrivenMomentum());
  math::PtEtaPhiELorentzVectorD theScp4( selectedSC.at(0).energy()/cosh(selectedSC.at(0).eta()), selectedSC.at(0).eta(), selectedSC.at(0).phi(), selectedSC.at(0).energy() );
  Z+=theScp4;

  if ( massMin_ < Z.M() && Z.M() < massMax_ ){
    if(debug_) std::cout << "++debug: mass request satisfied;" << std::endl; }
  else{
    if(debug_) std::cout << "++debug: mass request failed;" << std::endl;
    return false;  }    

  nEvt3++; // count how many pass mass requirement


  if(debug_) std::cout << "++debug: final PASS " << std::endl;
  return true ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
EleFilter::beginJob()
{
  nEvt1 = 0 ;   nEvt2 = 0 ;   nEvt3 = 0 ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EleFilter::endJob() {
  std::cout << "\n\nResults of filtering: " << std::endl ; 
  std::cout << "Saw " << nEvt <<  " events in total"  << std::endl ; 
  std::cout << "    " << nEvt1 << " events with at least required electron(s) (" << minEleNum_ << " with pt > " << minElePt_ << ")" << std::endl ; 
  std::cout << "    " << nEvt2 << " events with additionally the required SC(s) (" << minSCNum_ << " with pt > " << minSCEt_ << " and fabs(eta)<" << maxSCAbsEta_ <<  ")" << std::endl ; 
  std::cout << "    " << nEvt3 << " events with additionally the required mass requirements (" << massMin_ << " < m < " << massMax_ << std::endl ; 
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(EleFilter);
