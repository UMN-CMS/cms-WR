// -*- C++ -*-
//
// Package:    HeavyNuEleTriggerEff
// Class:      HeavyNuEleTriggerEff
// 
/**\class HeavyNuEleTriggerEff HeavyNuEleTriggerEff.cc HeavyNu/HeavyNuEleTriggerEff/src/HeavyNuEleTriggerEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giovanni Franzoni,27 2-013,+41227678347,
//         Created:  Fri May 18 12:18:35 CEST 2012
// $Id: HeavyNuEleTriggerEff.cc,v 1.2 2012/05/23 15:47:12 franzoni Exp $
//
//

// if uncommented, a chunk of code will be ignored
#define momentarilyIngnore 1

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"


#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"


//
// class declaration
//

class HeavyNuEleTriggerEff : public edm::EDAnalyzer {
public:
  explicit HeavyNuEleTriggerEff(const edm::ParameterSet&);
  ~HeavyNuEleTriggerEff();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  bool passOfflineSelection( const edm::Event & , const edm::EventSetup &, std::vector< std::pair<pat::Electron, float> > &, std::vector< std::pair<pat::Jet, float> > & );
  pat::TriggerObjectRefVector findObjectsMatchedToPath( const edm::Handle< pat::TriggerEvent > triggerEvent, const pat::TriggerPathRef seedingAcceptedPath,  const std::vector<std::string>& electronFilters );
  
  
  // ----------member data ---------------------------
  edm::InputTag rhoTag_ ; 
  double        rho_ ;
  edm::InputTag electronTag_;
  double maxAbsEtaOfflEle_;
  double minPtOfflEle_;
  int    heepVersion_;
  double ebScale_;
  double eeScale_;

  edm::InputTag jetTag_;
  int           numOfflJets_;
  double        minPtOfflJets_;
  

  edm::InputTag trigEventTag_;
  edm::InputTag isolatedEmSource_;
  edm::InputTag nonIsolatedEmSource_;
  bool          doDebugMessages_;

  std::vector<std::string> seedHLTForL1Efficiency_;
  std::string              targetL1Algo_;

  std::vector<std::string> seedHLTForHLTEfficiency_;
  double                   minPtOfObjects_;
  std::vector<std::string> targetHLTPaths_;


  std::vector<std::string> electronSeedingFilters_;
  std::vector<std::string> electronTargetFilters_;

  unsigned int counterEvtsAll_;
  unsigned int counterEvtsPassingOfflineSelection_;
  unsigned int counterEvtsWithSeedingFired_;
  unsigned int counterEvtsWithSeedingAbovePt_;
  unsigned int counterEvtsWithTargetValid_;
  unsigned int counterEvtsWithTargetFired_;

  std::string plotFolderName_;
  
  //TFileDirectory* thePlotsDir;
  TH1F*  massDenominator_;
  TH1F*  massNumerator_;
  TH1F*  eventsFate_;

  TH1F*  pTele1Denominator_;
  TH1F*  pTele1Numerator_;
  TH1F*  pTele2Denominator_;
  TH1F*  pTele2Numerator_;


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
HeavyNuEleTriggerEff::HeavyNuEleTriggerEff(const edm::ParameterSet& iConfig) :
  rhoTag_              ( iConfig.getParameter< edm::InputTag > ("rhoForHEEPId") ),
  electronTag_         ( iConfig.getParameter<edm::InputTag>("electronTag") ), 
  maxAbsEtaOfflEle_    ( iConfig.getParameter< double >("maxAbsEtaOfflEle") ), 
  minPtOfflEle_        ( iConfig.getParameter< double >("minPtOfflEle") ), 
  heepVersion_         ( iConfig.getParameter< int >("heepVersion") ), 
  ebScale_             ( iConfig.getParameter< double >("ebScale") ), 
  eeScale_             ( iConfig.getParameter< double >("eeScale") ), 
  jetTag_              ( iConfig.getParameter<edm::InputTag>("jetTag") ), 
  numOfflJets_         ( iConfig.getParameter< int >("numOfflJets") ),
  minPtOfflJets_       ( iConfig.getParameter< double >("minPtOfflJets") ), 
  trigEventTag_        ( iConfig.getParameter< edm::InputTag > ( "trigEventTag" ) ),
  isolatedEmSource_    ( iConfig.getParameter< edm::InputTag > ( "nonIsolatedEmSource" ) ),
  nonIsolatedEmSource_ ( iConfig.getParameter< edm::InputTag > ( "nonIsolatedEmSource" ) ),
  doDebugMessages_     ( iConfig.getParameter< bool > ( "doDebugMessages" )       ),
  seedHLTForL1Efficiency_ ( iConfig.getParameter< std::vector<std::string> > ( "seedHLTForL1Efficiency" ) ),
  targetL1Algo_        ( iConfig.getParameter<std::string> ( "targetL1Algo" ) ),
  seedHLTForHLTEfficiency_ ( iConfig.getParameter< std::vector<std::string> > ( "seedHLTForHLTEfficiency" ) ),
  minPtOfObjects_      ( iConfig.getParameter<double> ( "minPtOfObjects" ) ),
  targetHLTPaths_       ( iConfig.getParameter< std::vector<std::string> > ( "targetHLTPaths" ) ),
  electronSeedingFilters_       ( iConfig.getParameter< std::vector<std::string> > ( "electronSeedingFilters" ) ),
  electronTargetFilters_       ( iConfig.getParameter< std::vector<std::string> > ( "electronTargetFilters" ) ),
  counterEvtsAll_ (0),
  counterEvtsPassingOfflineSelection_ (0),
  counterEvtsWithSeedingFired_ (0),
  counterEvtsWithSeedingAbovePt_(0),
  counterEvtsWithTargetValid_ (0),
  counterEvtsWithTargetFired_ (0),
  plotFolderName_                ( iConfig.getParameter< std::string > ( "plotFolderName" ) )
{

  //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  //TFileDirectory thePlotsDir     =  fs->mkdir(plotFolderName_.c_str());
  massDenominator_ = (TH1F*)   fs->make<TH1F>("massDenominator","massDenominator; m(ee) [GeV]",1000,0.,2000);
  massNumerator_   = (TH1F*)   fs->make<TH1F>("massNumerator","massNumerator; m(ee) [GeV]",1000,0.,2000);
  pTele1Denominator_ = (TH1F*) fs->make<TH1F>("pTele1Denominator","pTele1Denominator; p_{T}(e1) [GeV]",1000,0.,2000.);
  pTele1Numerator_   = (TH1F*) fs->make<TH1F>("pTele1Numerator","pTele1Numerator; p_{T}(e1) [GeV]",1000,0.,2000);
  pTele2Denominator_ = (TH1F*) fs->make<TH1F>("pTele2Denominator","pTele2Denominator; p_{T}(e2) [GeV]",1000,0.,2000.);
  pTele2Numerator_   = (TH1F*) fs->make<TH1F>("pTele2Numerator","pTele2Numerator; p_{T}(e2) [GeV]",1000,0.,2000);
  eventsFate_      = (TH1F*)   fs->make<TH1F>("eventsFate","eventsFate; the events fate",10,0.,10);
  
}


HeavyNuEleTriggerEff::~HeavyNuEleTriggerEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HeavyNuEleTriggerEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  counterEvtsAll_++;
  eventsFate_->Fill(0);

  // require offline selections:
  // - 2 electrons passing HEEP, agove 40 GeV
  // - a numer of jets which is configurable: typically 0 or 2 
  // if failed, don't consider the event

  std::vector< std::pair<pat::Electron, float> > theOfflineCands;
  std::vector< std::pair<pat::Jet,      float> > jetCands;
  if ( passOfflineSelection( iEvent , iSetup, theOfflineCands, jetCands) ){
    if (doDebugMessages_) std::cout << "passOfflineSelection returns: " << passOfflineSelection( iEvent , iSetup, theOfflineCands, jetCands) << std::endl;
    counterEvtsPassingOfflineSelection_ ++;
    eventsFate_->Fill(1);
  }
  else {
    if (doDebugMessages_) std::cout << "(in analyze) size of theOfflineCands: " << theOfflineCands.size() << " and size of jet coll is: " << jetCands.size() << std::endl;
    // move on to the next event if basic offline conditions are not satisfied 
    return;
  }  


  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( trigEventTag_, triggerEvent );
   if ( !triggerEvent.isValid() ) {
     std::cout << "triggerEvent not found " << std::endl;
     assert(0);
   }
   else { if(doDebugMessages_ && 0) std::cout << " triggerEvent was fuond " << std::endl;}
   

   

#ifndef momentarilyIngnore
   if (doDebugMessages_) std::cout << "seedHLTForL1Efficiency has size: " << seedHLTForL1Efficiency_.size() << std::endl;


   for (unsigned int i=0; i<seedHLTForL1Efficiency_.size(); i++) { 

     // std::cout << "Investigating path: " << seedHLTForL1Efficiency_.at(i) << std::endl ; 
     const pat::TriggerPathRef iPath = triggerEvent->pathRef( seedHLTForL1Efficiency_.at(i) ) ; 
     if ( iPath.isNonnull() && iPath->wasAccept() ) { 
       std::cout << "Found path!" << std::endl ; 
       
       std::cout << "Path information: " 
		 << " " << iPath->name() << " " 
		 << " " << iPath->index() << " " 
		 << " pr: " << iPath->prescale() << " " 
		 << " r: " << iPath->wasRun() << " " 
		 << " acc: " << iPath->wasAccept() << " " 
		 << " " << iPath->wasError() << " " 
		 << " last: " << iPath->lastActiveFilterSlot() << " " 
		 << " " << iPath->modules().size() ; 
       std::cout // << iPath->modules().at(iPath->lastActiveFilterSlot()) << " "
	 << iPath->l3Filters() << " " 
	 << iPath->xTrigger() << " " ; 
       // for (unsigned int j=0; j<iPath->filterIndices().size(); j++) { 
       //   std::cout << iPath->modules().at(iPath->filterIndices().at(j)) << " " ; 
       // }
       std::cout << " \n\n\n\n" 
		 << std::endl ; 
       
       std::cout << "Looking for trigger objects involved in this path: " <<  iPath->name() << std::endl ; 
 
       pat::TriggerObjectRefVector objectsInPath = triggerEvent->pathObjects(iPath->name(),true) ;
       pat::TriggerFilterRefVector filtersInPath = triggerEvent->pathFilters(iPath->name(),true) ;
       std::cout  << "size of objectsInPath : " << objectsInPath.size() << std::endl; 

       
       // http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/T>> Entering Package PhysicsTools/PatAlgos
       // // Get a reference to a certain L1 algorithm by name,
       // // NULL, if algorithm is not found
       //    const TriggerAlgorithmRef algorithmRef( const std::string & nameAlgorithm ) const;


       // http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/TriggerEvent.h#308
       // 
       // TriggerObjectRefVector algorithmObjects( const std::string & nameAlgorithm ) const;

       //       http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/TriggerEvent.h#310 
       //       /// Checks, if an object was used in a certain algorithm given by name
       //       bool objectInAlgorithm( const TriggerObjectRef & objectRef, const std::string & nameAlgorithm ) const;


     }
   }

   pat::TriggerAlgorithmRefVector theAlgosActive = triggerEvent->algorithmRefs();
   std::cout << "size of theAlgosActive: " << theAlgosActive.size() << std::endl;

   pat::TriggerAlgorithmRefVector theAlgosPassed = triggerEvent->acceptedAlgorithms();
   std::cout << "size of theAlgosPassed: " << theAlgosPassed.size() << std::endl;

   pat::TriggerAlgorithmRefVector thePhysicsAlgosPassed = triggerEvent->acceptedPhysAlgorithms();
   std::cout << "size of thePhysicsAlgosPassed: " << thePhysicsAlgosPassed.size() << std::endl;

   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTrigger#PATTriggerProducer
   std::cout << "physics algos which passed: " << std::endl;
   for(unsigned int phAl=0; phAl<thePhysicsAlgosPassed.size(); phAl++)
     {
       std::cout << thePhysicsAlgosPassed.at(phAl)->name() << "\n";
       pat::TriggerObjectRefVector objectsInAlgo =  triggerEvent->algorithmObjects( thePhysicsAlgosPassed.at(phAl)->name() );
       // std::cout << "\t\t\t\t size of objectsInAlgo: " << objectsInAlgo.size() << std::endl;
     }
   std::cout << std::endl;



   // Isolated EM particles
   Handle< l1extra::L1EmParticleCollection > isoEmColl ;
   iEvent.getByLabel( isolatedEmSource_, isoEmColl ) ;
   std::cout << "++ HeavyNuEleTriggerEff: Number of isolated EM " << isoEmColl->size() << std::endl ;
   
   for( l1extra::L1EmParticleCollection::const_iterator emItr = isoEmColl->begin() ;
	emItr != isoEmColl->end() ;
	++emItr )
     {
       std::cout << "\t  p4 (" << emItr->px()
	    << ", " << emItr->py()
	    << ", " << emItr->pz()
	    << ", " << emItr->energy()
	    << ") et " << emItr->et()
	    << " eta " << emItr->eta()
	    << " phi " << emItr->phi()
		 << std::endl ;
     }
   
   // Non-isolated EM particles
   Handle< l1extra::L1EmParticleCollection > nonIsoEmColl ;
   iEvent.getByLabel( nonIsolatedEmSource_, nonIsoEmColl ) ;
   std::cout << "++ HeavyNuEleTriggerEff: Number of non-isolated EM " << nonIsoEmColl->size() << std::endl ;
   
   for( l1extra::L1EmParticleCollection::const_iterator emItr = nonIsoEmColl->begin() ;
	emItr != nonIsoEmColl->end() ;
	++emItr )
     {
       std::cout << "\t\t  p4 (" << emItr->px()
	    << ", " << emItr->py()
	    << ", " << emItr->pz()
	    << ", " << emItr->energy()
	    << ") et " << emItr->et()
	    << " eta " << emItr->eta()
	    << " phi " << emItr->phi()
		 << std::endl ;
     }
   

#endif





   ////////////////////////////////////////////////////////////////////////
   // HLT efficiencies from here below ////////////////////////////////////
   // counterEvtsAll_++;  // migrated up at the very beginning of analyze

   //////////////////////////// prepare seeding trigger path /////////////////
   std::vector <std::string> seedHLTForHLTEfficiencyAccepted;
   std::vector <std::string> targetHLTPathsAccepted;
   pat::TriggerPathRefVector theAcceptedPaths =  triggerEvent->acceptedPaths();
   // check if any of the seeding HLT paths have passed
   for ( pat::TriggerPathRefVector::const_iterator iPath = theAcceptedPaths.begin(); iPath != theAcceptedPaths.end(); ++iPath ) {

     for (unsigned int i=0; i<seedHLTForHLTEfficiency_.size(); i++) { 
       if ( seedHLTForHLTEfficiency_.at(i) == (*iPath)->name() ) 	 { seedHLTForHLTEfficiencyAccepted.push_back( (*iPath)->name() ); }
     }// loop over active paths seeking for the seeding path

     for (unsigned int i=0; i<targetHLTPaths_.size(); i++) { 
       if ( targetHLTPaths_.at(i) == (*iPath)->name() ) 	         { targetHLTPathsAccepted.push_back( (*iPath)->name() ); }
     }// loop over active paths seeking for the target path
     
   }// loop over HLT-eff seeding paths
   


   // if the seeding path has not fired, stop execution
   if ( seedHLTForHLTEfficiencyAccepted.size()>0 ) {
     counterEvtsWithSeedingFired_++;
     eventsFate_->Fill(2);
   }
   else                                            return; 




   //////////////////////////// prepare target trigger path /////////////////
   // at least one target HLT path must be provided; if not, bail out
   if (targetHLTPaths_.size()==0) {
     std::cout << "you've not provided targetHLTPaths; bailing out" << std::endl;
     assert(0);
   }
   std::vector <std::string> targetHLTPathsInMenu;
   //check if any of the provided target triggers are in the menu
   pat::TriggerPathRefVector thePathsInMenu =  triggerEvent->pathRefs();
   for (unsigned int i=0; i<targetHLTPaths_.size(); i++) { 
     for ( pat::TriggerPathRefVector::const_iterator iPath = thePathsInMenu.begin(); iPath != thePathsInMenu.end(); ++iPath ) {
       if ( targetHLTPaths_.at(i) == (*iPath)->name() ) 	 { targetHLTPathsInMenu.push_back( (*iPath)->name() ); }
     }// loop over active paths
   }// loop over HLT-eff seeding paths



   // if there's not at least one target path which is validly in the menu, stop execution
   if ( targetHLTPathsInMenu.size()>0 ) {
     counterEvtsWithTargetValid_++;
     eventsFate_->Fill(3);
   }
   else                                 return;


   // fill plots for "BEFORE CUTS" HERE


  
   if (doDebugMessages_) {
     std::cout << "number of active HLT paths: " << thePathsInMenu.size() << std::endl; // why do I care?
     std::cout << "seedHLTForHLTEfficiency has size: " << seedHLTForHLTEfficiency_.size() 
	       << " of which  " << seedHLTForHLTEfficiencyAccepted.size() << " were 'run in the menu' && 'PASSED'" << std::endl;
     std::cout << "targetHLTPaths.at(0) is: " << targetHLTPaths_.at(0) << std::endl;
     std::cout << "targetHLTPaths size: " << targetHLTPaths_.size() 
	       << " of which  " << targetHLTPathsInMenu.size() << " were run in the menu (they may have passed or failed)" << std::endl;
   }
   


   // do the HLT efficnecy study only if there's at least:
   // - a seeding path which was accepted
   // - a target path which was specified in input by the user and is in the menu
   // ==> conditions are guaranteed by the return statements of above
   
   //   if ( seedHLTForHLTEfficiencyAccepted.size()>0 && targetHLTPathsInMenu.size()>0 ) 
   //     {
   
   // one of the seding HLT paths has fired. Chose the first one as a base trigger
   const pat::TriggerPathRef seedingAcceptedPath = triggerEvent->pathRef( seedHLTForHLTEfficiencyAccepted.at(0) ) ; 
   
   if (doDebugMessages_) {
     std::cout << "\n\n\t\tSeeding path information: " 
	       << " " << seedingAcceptedPath->name() << " " 
	       << " " << seedingAcceptedPath->index() << " " 
	       << " pr: " << seedingAcceptedPath->prescale() << " " 
	       << " r: " << seedingAcceptedPath->wasRun() << " " 
	       << " acc: " << seedingAcceptedPath->wasAccept() << " " 
	       << " " << seedingAcceptedPath->wasError() << " " 
	       << " last: " << seedingAcceptedPath->lastActiveFilterSlot() << " " 
	       << " " << seedingAcceptedPath->modules().size() ; 
     std::cout // << seedingAcceptedPath->modules().at(seedingAcceptedPath->lastActiveFilterSlot()) << " "
       << seedingAcceptedPath->l3Filters() << " " 
       << seedingAcceptedPath->xTrigger() << " " ; 
   }
   
   
   // get the HLT objects which match the seeding HLT path
   pat::TriggerObjectRefVector objectsMatchedToSeedingPath = findObjectsMatchedToPath(triggerEvent, seedingAcceptedPath , electronSeedingFilters_);
   if (doDebugMessages_) {  std::cout << "found: "<<objectsMatchedToSeedingPath.size()<<" objects matching the path: "<<seedingAcceptedPath->name()<<"\n" << std::endl; }
   // check how many of those objects have pt>minPtOfObjects_=33 (which is required by target trigger)
   unsigned int counterSeedingObjectsAbovePt(0); 
   for (pat::TriggerObjectRefVector::const_iterator funct = objectsMatchedToSeedingPath.begin(); funct !=objectsMatchedToSeedingPath.end(); funct++){
     if (doDebugMessages_) std::cout << "objectsMatchedToPath SEED = " << (*funct)->pt() << " and eta " << (*funct)->eta() << std::endl ;
     if ( (*funct)->pt() >minPtOfObjects_ ) {
       counterSeedingObjectsAbovePt++;
     }
   }
   
   // the seeding trigger must have at least two of its objects above the PT treshold of the target trigger 
   if ( counterSeedingObjectsAbovePt < 2)
     {
       if (doDebugMessages_) std::cout << "objectsMatchedToPath SEED with pt > " << minPtOfObjects_ << " : " << counterSeedingObjectsAbovePt << std::endl ;
       return;   
     }
   else 
     {
       if (doDebugMessages_) std::cout << "objectsMatchedToPath SEED with pt > " << minPtOfObjects_ << " : " << counterSeedingObjectsAbovePt << std::endl ;
       counterEvtsWithSeedingAbovePt_++;
       eventsFate_->Fill(4);
     }
   
   // how do I chose which of the offline candidates are to be used? For now, chose the first two in the collection.. 
   massDenominator_ -> Fill( ( theOfflineCands.at(0).first.p4() +  theOfflineCands.at(1).first.p4() ).M() ); 
   if ( theOfflineCands.at(0).first.pt() > theOfflineCands.at(1).first.pt()){
     pTele1Denominator_ -> Fill( theOfflineCands.at(0).first.pt() );
   }

   if ( targetHLTPathsAccepted.size()>0) {
     counterEvtsWithTargetFired_++;
     eventsFate_->Fill(5);
     massNumerator_ -> Fill( ( theOfflineCands.at(0).first.p4() +  theOfflineCands.at(1).first.p4() ).M() ); 
   }
   else return;

   // get the HLT objects which match the target HLT path; if more than one target HLT path, chose the first one (should be irrelevant if analysis ORs them all)
   const pat::TriggerPathRef targetHLTPath = triggerEvent->pathRef( targetHLTPathsInMenu.at(0) ) ; 
   pat::TriggerObjectRefVector objectsMatchedToTargetPath = findObjectsMatchedToPath(triggerEvent, targetHLTPath , electronTargetFilters_);
   if (doDebugMessages_) {
     std::cout << "outside found: " << objectsMatchedToTargetPath.size() << " objects matching the path: " << targetHLTPath->name() << "\n" << std::endl;}
   for (pat::TriggerObjectRefVector::const_iterator funct2 = objectsMatchedToTargetPath.begin(); funct2 !=objectsMatchedToTargetPath.end(); funct2++){
     if (doDebugMessages_) std::cout << "objectsMatchedToPath TARGET = " << (*funct2)->pt() << " and eta " << (*funct2)->eta() << std::endl ;
   }
   
   
   
   
   
}// end analyze


bool HeavyNuEleTriggerEff::passOfflineSelection( const edm::Event& iEvent, 
					  const edm::EventSetup& iSetup,
					  std::vector< std::pair<pat::Electron, float> > & eCands,
					  std::vector< std::pair<pat::Jet, float> > & jetCands 
					  )
{

  eCands  .clear();
  jetCands.clear();
  
  // needed by the HEEP 
  edm::Handle<double> electronRhoHandle ; 
  iEvent.getByLabel(rhoTag_, electronRhoHandle) ; 
  if ( electronRhoHandle.isValid() ) 
    {       rho_      = ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.) ; 
      if (doDebugMessages_) std::cout << " valid handle to rho found. rho: " << rho_ << std::endl;     }
  else
    {       std::cout << "NO valid handle to rho found which was expected: " << rhoTag_ << " Bailing out " << std::endl;
      assert (0);      }
  


  // electrons to start with
  edm::Handle<pat::ElectronCollection> patElectronCollection ; 
  iEvent.getByLabel(electronTag_,patElectronCollection) ; 
  if (!patElectronCollection.isValid()) {      
    std::cout << "no pat electron collection: " << electronTag_ << " found; bailing out." << std::endl;         assert(0); } 
  else{
    if (doDebugMessages_) std::cout << " pat electron collection found with size: " << (*patElectronCollection.product()).size() << std::endl;     }


  std::vector< std::pair<pat::Electron, float> >  eList;
  eList =   hnu::getElectronList( patElectronCollection, 
				  maxAbsEtaOfflEle_,
				  minPtOfflEle_, minPtOfflEle_, 
				  heepVersion_,
				  rho_);
  if (eList.size()<2) return false;
  


  // jets to start with
  edm::Handle<pat::JetCollection> pJets;
  iEvent.getByLabel(jetTag_, pJets);
  if (!pJets.isValid()) {      
    std::cout << "no pat pJets collection: " << jetTag_ << " found; bailing out." << std::endl;         assert(0); } 
  else{
    if (doDebugMessages_) std::cout << " pat jets collection found with size: " << (*pJets.product()).size() << std::endl;     }
  // corrections for jets
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord > ().get("AK5Calo", JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *  jecuObj_ = new JetCorrectionUncertainty(JetCorPar);
  
  double maxJetAbsEta  = 4.8; 
  int    applyJECUsign = 0;
  int    jecEta        = 3; 
  bool   isMC          = false;
  int    applyJERsign  = 0 ;
  std::vector< std::pair<pat::Jet, float> > jetCandsBefId =    hnu::getJetList(pJets, 
									  jecuObj_, 
									  minPtOfflJets_, maxJetAbsEta,
									  applyJECUsign, 
									  jecEta, 
									  isMC, 
									  applyJERsign);
  if (jetCandsBefId.size()< static_cast<unsigned int> (numOfflJets_) ) return false;
  // require JET ID
  for( unsigned int Jcand=0; Jcand<jetCandsBefId.size(); Jcand++ ){
    if ( hnu::jetID( jetCandsBefId.at(Jcand).first ) > 1 ) jetCands.push_back( jetCandsBefId.at(Jcand) );
  }



  // check that electrtons are 'far' from any of the jets
  for (unsigned  int theEl=0; theEl< eList.size(); theEl++ )
    {
      bool isFar(true);
      
      for( unsigned int Jcand=0; Jcand<jetCands.size(); Jcand++ ){
	if (  deltaR( eList.at(theEl).first.eta(), eList.at(theEl).first.phi() , 
		      jetCands.at(Jcand).first.eta(), jetCands.at(Jcand).first.phi() ) 
	      < 0.5 ) 	  isFar=false;
      }
      
      if (isFar) eCands.push_back( eList.at(theEl) ) ;
   
    }


  
  if (doDebugMessages_) std::cout << "(inside passOfflineSelection - at least 2 ele and numOfflJets_ raw jets) size of HEEP offline candidates: " << eCands.size() 
				  << " and size of selected jet candidates: " << jetCands.size() 
				  << " (were: " << jetCandsBefId.size() << " before ID) " << (jetCandsBefId.size()-jetCands.size()) <<  std::endl;




  return (
	  eCands.size()   >= 2  &&
	  jetCands.size() >= static_cast<unsigned int> (numOfflJets_)
	  )
    ;
}



// ------------ method which returns references to objects that match a certain trigger  ------------
pat::TriggerObjectRefVector HeavyNuEleTriggerEff::findObjectsMatchedToPath( const edm::Handle< pat::TriggerEvent > triggerEvent , 
								     const pat::TriggerPathRef seedingAcceptedPath,
								     const std::vector<std::string>& electronFilters )
{

  // one of the seding HLT paths has fired. Chose the first one as a base trigger
  //const pat::TriggerPathRef seedingAcceptedPath = triggerEvent->pathRef( seedHLTForHLTEfficiencyAccepted.at(0) ) ; 
  
  if (doDebugMessages_) {
    std::cout << "\n\n\t\tPath information: " 
	      << " " << seedingAcceptedPath->name() << " " 
	      << " " << seedingAcceptedPath->index() << " " 
	      << " pr: " << seedingAcceptedPath->prescale() << " " 
	      << " r: " << seedingAcceptedPath->wasRun() << " " 
	      << " acc: " << seedingAcceptedPath->wasAccept() << " " 
	      << " " << seedingAcceptedPath->wasError() << " " 
	      << " last: " << seedingAcceptedPath->lastActiveFilterSlot() << " " 
	      << " " << seedingAcceptedPath->modules().size() ; 
    std::cout	   << seedingAcceptedPath->l3Filters() << " " 
		   << seedingAcceptedPath->xTrigger() << " " << std::endl;; 
  }

  if (doDebugMessages_) std::cout << "\nLooking for trigger objects involved in this path: " <<  seedingAcceptedPath->name() << std::endl ; 
  
  pat::TriggerObjectRefVector objectsInPath = triggerEvent->pathObjects(seedingAcceptedPath->name(),true) ;
  pat::TriggerFilterRefVector filtersInPath = triggerEvent->pathFilters(seedingAcceptedPath->name(),true) ;
  if (doDebugMessages_) { std::cout << "for path: " << seedingAcceptedPath->name()  << " size of objectsInPath : " << objectsInPath.size() 
				    << " and filtersInPath: " << filtersInPath.size() << std::endl; }
  
  
  pat::TriggerObjectRefVector objectsMatchedToPath;
  
  // find the two objects from the seeding path
  // code largely inspired to Bryan's 
  // =->  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Minnesota/HeavyNu/AnalysisModules/src/HeavyNuTrigger.cc?revision=1.22&view=markup&sortby=date
  for ( pat::TriggerFilterRefVector::const_iterator ifRef = filtersInPath.begin(); ifRef != filtersInPath.end(); ifRef++) { 
    
    pat::TriggerFilterRef filterRef = *ifRef ; 
    if (doDebugMessages_)  std::cout << "filter in path : " << filterRef->label() <<  std::endl; 
    
    if ( filterRef->isFiring() &&
	 ( std::find(electronFilters.begin(), electronFilters.end(), filterRef->label()) != electronFilters.end() ) ) { 
      if (doDebugMessages_)  std::cout << "==> found filter which is among electronFilters: " << filterRef->label() << std::endl;
      
      
      for ( pat::TriggerObjectRefVector::const_iterator iobjRef = objectsInPath.begin(); iobjRef != objectsInPath.end(); ++iobjRef ) {
	pat::TriggerObjectRef objRef = *iobjRef ; 
	if (doDebugMessages_) std::cout << "Found an object with pT = " << objRef->pt() << " and eta " << objRef->eta() << std::endl ; 
	
	if ( triggerEvent->objectInFilter(objRef,filterRef->label()) ) { // Trigger object was used by the filter
	  if (doDebugMessages_)  std::cout << "found candidate requested IN FILTER object with pt: " <<  objRef->pt() << " and eta " << objRef->eta() << std::endl;
	  objectsMatchedToPath.push_back(objRef);
	}
	
      }// if filter found
    }// loop on the _last_ filters of the seed path
  }//loop on  the filters of the path
  
  if (doDebugMessages_) std::cout << "inside found: " << objectsMatchedToPath.size() << " objects matching the path: " 
				  << /*seedHLTForHLTEfficiencyAccepted.at(0)*/ seedingAcceptedPath->name() << std::endl;
  
  return objectsMatchedToPath;
}


// ------------ method called once each job just before starting event loop  ------------
void 
HeavyNuEleTriggerEff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HeavyNuEleTriggerEff::endJob() 
{

  std::cout << "\n\nHeavyNuEleTriggerEff has seen: " << counterEvtsAll_ << " events;" << std::endl;
  
  std::cout << "\t\t ++ of which: " << counterEvtsPassingOfflineSelection_ << " pass the offline selection : 2 HEEP ele pt > "
	    << minPtOfflEle_ << " and" << numOfflJets_ << " jets pt > " << minPtOfflJets_ << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithSeedingFired_ << " with >=1 seeding paths firing" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithTargetValid_  << " with >=1 valid target path" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithSeedingAbovePt_ << " with >=2 seeding objects above pt " << minPtOfObjects_ << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithTargetFired_  << " with >=1 fired target path" << std::endl;
  std::cout << "\n" << std::endl;

}

// ------------ method called when starting to processes a run  ------------
void 
HeavyNuEleTriggerEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HeavyNuEleTriggerEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HeavyNuEleTriggerEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HeavyNuEleTriggerEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HeavyNuEleTriggerEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuEleTriggerEff);
