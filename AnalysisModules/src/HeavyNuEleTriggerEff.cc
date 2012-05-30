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
// $Id: HeavyNuEleTriggerEff.cc,v 1.7 2012/05/29 22:37:40 franzoni Exp $
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
  // the first run for which the run monitoring histogram will monitor
  const unsigned int firstRunForPlotting_;
  bool               verbosityForRunPlotting_;

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
  
  unsigned int runMin_;
  unsigned int runMax_;

  edm::InputTag trigEventTag_;
  edm::InputTag isolatedEmSource_;
  edm::InputTag nonIsolatedEmSource_;
  bool          doDebugMessages_;

  std::vector<std::string> seedHLTForL1Efficiency_;
  std::string              targetL1Algo_;

  std::vector<std::string> seedHLTForHLTEfficiency_;
  double                   minPtOfObjects_;
  std::vector<std::string> targetHLTPaths_;


  std::vector<std::string> electronSeedingPathEndFilter_;
  std::vector<std::string> electronSeedingPathTagFilter_;
  std::vector<std::string> electronTargetFilters_;

  unsigned int counterEvtsAll_;
  unsigned int counterEvtsPassingOfflineSelection_;
  unsigned int counterEvtsWithSeedingFired_;
  unsigned int counterEvtsWithSeedingAbovePt_;
  unsigned int counterEvtsWithTargetValid_;
  unsigned int counterEvtsWithTargetFired_;

  std::string plotFolderName_;

  //TFileDirectory* thePlotsDir;
  TH1F*  monitorRuns_;

  TH1F*  massDenominator_;
  TH1F*  massNumerator_;
  TH1F*  massFail_;
  TH1F*  eventsFate_;

  TH1F*  pTele1Denominator_;
  TH1F*  pTele1Numerator_;
  TH1F*  pTele1Fail_;
  TH1F*  pTele2Denominator_;
  TH1F*  pTele2Numerator_;
  TH1F*  pTele2Fail_;
  TH1F*  pTTagDenominator_;
  TH1F*  pTTagNumerator_;
  TH1F*  pTTagFail_;
  TH1F*  pTProbeDenominator_;
  TH1F*  pTProbeNumerator_;
  TH1F*  pTProbeFail_;

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
  firstRunForPlotting_  (190645),
  verbosityForRunPlotting_(true),
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
  runMin_              ( iConfig.getParameter< int >("runMin") ), 
  runMax_              ( iConfig.getParameter< int >("runMax") ), 
  trigEventTag_        ( iConfig.getParameter< edm::InputTag > ( "trigEventTag" ) ),
  isolatedEmSource_    ( iConfig.getParameter< edm::InputTag > ( "nonIsolatedEmSource" ) ),
  nonIsolatedEmSource_ ( iConfig.getParameter< edm::InputTag > ( "nonIsolatedEmSource" ) ),
  doDebugMessages_     ( iConfig.getParameter< bool > ( "doDebugMessages" )       ),
  seedHLTForL1Efficiency_ ( iConfig.getParameter< std::vector<std::string> > ( "seedHLTForL1Efficiency" ) ),
  targetL1Algo_        ( iConfig.getParameter<std::string> ( "targetL1Algo" ) ),
  seedHLTForHLTEfficiency_ ( iConfig.getParameter< std::vector<std::string> > ( "seedHLTForHLTEfficiency" ) ),
  minPtOfObjects_      ( iConfig.getParameter<double> ( "minPtOfObjects" ) ),
  targetHLTPaths_       ( iConfig.getParameter< std::vector<std::string> > ( "targetHLTPaths" ) ),
  electronSeedingPathEndFilter_       ( iConfig.getParameter< std::vector<std::string> > ( "electronSeedingPathEndFilter" ) ),
  electronSeedingPathTagFilter_       ( iConfig.getParameter< std::vector<std::string> > ( "electronSeedingPathTagFilter" ) ),
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
  monitorRuns_       = (TH1F*)   fs->make<TH1F>("monitorRuns","monitorRuns; run (starting coll 2012)",10000,firstRunForPlotting_,firstRunForPlotting_+20000);
  massDenominator_   = (TH1F*)   fs->make<TH1F>("massDenominator","massDenominator; m(ee) [GeV]",1000,0.,2000);
  massNumerator_     = (TH1F*)   fs->make<TH1F>("massNumerator","massNumerator; m(ee) [GeV]",1000,0.,2000);
  massFail_          = (TH1F*)   fs->make<TH1F>("massFail","massFail; m(ee) [GeV]",1000,0.,2000);
  pTele1Denominator_ = (TH1F*)   fs->make<TH1F>("pTele1Denominator","pTele1Denominator; p_{T}(e1) [GeV]",1000,0.,2000.);
  pTele1Numerator_   = (TH1F*)   fs->make<TH1F>("pTele1Numerator","pTele1Numerator; p_{T}(e1) [GeV]",1000,0.,2000);
  pTele1Fail_        = (TH1F*)   fs->make<TH1F>("pTele1Fail","pTele1Fail; p_{T}(e1) [GeV]",1000,0.,2000);
  pTele2Denominator_ = (TH1F*)   fs->make<TH1F>("pTele2Denominator","pTele2Denominator; p_{T}(e2) [GeV]",1000,0.,2000.);
  pTele2Numerator_   = (TH1F*)   fs->make<TH1F>("pTele2Numerator","pTele2Numerator; p_{T}(e2) [GeV]",1000,0.,2000);
  pTele2Fail_        = (TH1F*)   fs->make<TH1F>("pTele2Fail","pTele2Fail; p_{T}(e1) [GeV]",1000,0.,2000);
  eventsFate_        = (TH1F*)   fs->make<TH1F>("eventsFate","eventsFate; the events fate",10,0.,10);
  
  pTTagDenominator_   = (TH1F*)   fs->make<TH1F>("pTTagDenominator","pTTagDenominator; p_{T}(tag,denom) [GeV]",1000,0.,2000);
  pTTagNumerator_     = (TH1F*)   fs->make<TH1F>("pTTagNumerator","pTTagNumerator; p_{T}(tag,num) [GeV]",1000,0.,2000);
  pTTagFail_          = (TH1F*)   fs->make<TH1F>("pTTagFail","pTTagFail; p_{T}(tag,fail) [GeV]",1000,0.,2000);
  pTProbeDenominator_ = (TH1F*)   fs->make<TH1F>("pTProbeDenominator","pTProbeDenominator; p_{T}(probe,denom) [GeV]",1000,0.,2000);
  pTProbeNumerator_   = (TH1F*)   fs->make<TH1F>("pTProbeNumerator","pTProbeNumerator; p_{T}(probe,num) [GeV]",1000,0.,2000);
  pTProbeFail_        = (TH1F*)   fs->make<TH1F>("pTProbeFail","pTProbeFail; p_{T}(probe,fail) [GeV]",1000,0.,2000);

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

  if (iEvent.run() < runMin_ || runMax_ > iEvent.run()) return;

  if (firstRunForPlotting_ > iEvent.run() && verbosityForRunPlotting_)
    {
      std::cout << "PROBLEM: in HeavyNuEleTriggerEff there's an event from run: " << iEvent.run() 
		<< " while only runs starting from: " << firstRunForPlotting_ 
		<< " are properly monitored; will continue, and not repeat this message" << std::endl;
	verbosityForRunPlotting_ = false; 
    }
  monitorRuns_->Fill( iEvent.run() );
  

  counterEvtsAll_++;
  eventsFate_->Fill(0);


  // require offline selections:
  // - 2 electrons passing HEEP, agove 40 GeV
  // - a numer of jets which is configurable ( typically 0 or 2 ) pass jet ID
  // - the two electrons are 'far enough' from the jets
  // if these requirements fail, don't consider the event

  std::vector< std::pair<pat::Electron, float> > theEleOfflineCands;
  std::vector< std::pair<pat::Jet,      float> > jetCands;

  if ( passOfflineSelection( iEvent , iSetup, theEleOfflineCands, jetCands) ){
    if (doDebugMessages_) std::cout << "passOfflineSelection returns: " << passOfflineSelection( iEvent , iSetup, theEleOfflineCands, jetCands) << std::endl;
    counterEvtsPassingOfflineSelection_ ++;
    eventsFate_->Fill(1);
  }
  else {
    if (doDebugMessages_) { std::cout << "(in analyze) size of theEleOfflineCands: " 
				      << theEleOfflineCands.size() << " and size of jet coll is: " << jetCands.size() 
				      << std::endl; }

    // if basic offline conditions are not satisfied, move on to the next event
    return;
  }  

  // all sorts of information about trigger, from pat
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( trigEventTag_, triggerEvent );
   if ( !triggerEvent.isValid() ) {
     std::cout << "triggerEvent not found " << std::endl;
     assert(0);
   }
   else { if(doDebugMessages_ && 0) std::cout << " triggerEvent was fuond " << std::endl;}
   
   /*********************************************************
 
      there used to be a lot of RELIC L1 code here
      since we're not looking into it, for the moment, 
      I've moved it out of the way (at the end of this file)

   ************************************************************/

   ////////////////////////////////////////////////////////////////////////
   // HLT efficiencies from here below ////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////


   /////////// prepare seeding trigger path /////////////////
   std::vector <std::string> seedHLTForHLTEfficiencyAccepted;
   std::vector <std::string> targetHLTPathsAccepted;
   pat::TriggerPathRefVector theAcceptedPaths =  triggerEvent->acceptedPaths();

   // check if any of the seeding HLT paths have passed
   for ( pat::TriggerPathRefVector::const_iterator iPath = theAcceptedPaths.begin(); iPath != theAcceptedPaths.end(); ++iPath ) {

     for (unsigned int i=0; i<seedHLTForHLTEfficiency_.size(); i++) { 
       if ( seedHLTForHLTEfficiency_.at(i) == (*iPath)->name() ) 	 { seedHLTForHLTEfficiencyAccepted.push_back( (*iPath)->name() ); }
     }// loop over active paths seeking for the __seeding__ path

     for (unsigned int i=0; i<targetHLTPaths_.size(); i++) { 
       if ( targetHLTPaths_.at(i) == (*iPath)->name() ) 	         { targetHLTPathsAccepted.push_back( (*iPath)->name() ); }
     }// loop over active paths seeking for the __target__ path
     
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

  
   if (doDebugMessages_) {
     std::cout << "number of active HLT paths: " << thePathsInMenu.size() << std::endl; // why do I care?
     std::cout << "seedHLTForHLTEfficiency has size: " << seedHLTForHLTEfficiency_.size() 
	       << " of which  " << seedHLTForHLTEfficiencyAccepted.size() << " were 'run in the menu' && 'PASSED'" << std::endl;
     std::cout << "targetHLTPaths.at(0) is: " << targetHLTPaths_.at(0) << std::endl;
     std::cout << "targetHLTPaths size: " << targetHLTPaths_.size() 
	       << " of which  " << targetHLTPathsInMenu.size() << " were run in the menu (they may have passed or failed)" << std::endl;
   }
   


   ////////////////////////////////////////////////////////////////////////////////
   // do the HLT efficnecy study only if there's at least:
   // - a seeding path which was accepted
   // - a target path which was specified in input by the user and is in the menu
   // ==> conditions are guaranteed by the return statements of above
   
   
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
   
   
   // get the HLT objects which match the LAST filter of the seeding HLT path
   pat::TriggerObjectRefVector objectsMatchedToEndOfSeedingPath = findObjectsMatchedToPath(triggerEvent, seedingAcceptedPath , electronSeedingPathEndFilter_);
   if (doDebugMessages_) 
     {  std::cout << "found: "<<objectsMatchedToEndOfSeedingPath.size()<<" objects matching the path: "<<seedingAcceptedPath->name()<<"\n" << std::endl; }

   // check how many of those objects have pt>minPtOfObjects_=33 (which is required by target trigger)
   unsigned int counterSeedingObjectsAbovePt(0); 
   for (pat::TriggerObjectRefVector::const_iterator funct = objectsMatchedToEndOfSeedingPath.begin(); funct !=objectsMatchedToEndOfSeedingPath.end(); funct++){
     if (doDebugMessages_) std::cout << "objectsMatchedToPath SEED = " << (*funct)->pt() << " and eta " << (*funct)->eta() << std::endl ;
     if ( (*funct)->pt() >minPtOfObjects_ ) {
       counterSeedingObjectsAbovePt++;
     }
   }


   // now see if the objects which come from electronSeedingPathEndFilter_ match any of the offline electrons 
   unsigned int counterSeedingPathMatchedObjects(0);
   for (pat::TriggerObjectRefVector::const_iterator seedPathObj = objectsMatchedToEndOfSeedingPath.begin(); 
	seedPathObj !=objectsMatchedToEndOfSeedingPath.end(); seedPathObj++)
     {
       
       // theEleOfflineCands
       for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
	 {
	   
	   if (  deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
			 (*seedPathObj)->eta(), (*seedPathObj)->phi() ) 
		 < 0.1 ) {
	     counterSeedingPathMatchedObjects++;
	     break;
	   }// if there's matching
	   
	 }// loop over offline electrons 
     }//loop over HLT objects from the electronSeedingPathEndFilter_ filter of the seeding path

   // if the online objects are not matched to offline objects , remove this event from denominator
   if(counterSeedingPathMatchedObjects < objectsMatchedToEndOfSeedingPath.size() ) return;
   // if an event passes this point, it has entered the _denominator_
   // i.e. I have two online objects which have  fired the seeding trigger and which are matched to offline objects

   // how do I chose which of the offline candidates are to be used? For now, chose the first two in the collection... REVISIT? 
   massDenominator_ -> Fill( ( theEleOfflineCands.at(0).first.p4() +  theEleOfflineCands.at(1).first.p4() ).M() ); 
   if ( theEleOfflineCands.at(0).first.pt() > theEleOfflineCands.at(1).first.pt()){
     pTele1Denominator_ -> Fill( theEleOfflineCands.at(0).first.pt() );
     pTele2Denominator_ -> Fill( theEleOfflineCands.at(1).first.pt() );
   }


//   // COMMENT FROM HERE
//
//   // get the HLT objects which match the target HLT path;
//   // if more than one target HLT path, chose the first one (should be irrelevant if analysis ORs them all)
//   pat::TriggerPathRef targetHLTPath = triggerEvent->pathRef( targetHLTPathsInMenu.at(0) ) ; 
//   pat::TriggerObjectRefVector objectsMatchedToTargetPath = findObjectsMatchedToPath(triggerEvent, targetHLTPath , electronTargetFilters_);
//   if (doDebugMessages_) {
//     std::cout << "outside found: " << objectsMatchedToTargetPath.size() << " objects matching the path: " << targetHLTPath->name() << "\n" << std::endl;}
//   for (pat::TriggerObjectRefVector::const_iterator funct2 = objectsMatchedToTargetPath.begin(); funct2 !=objectsMatchedToTargetPath.end(); funct2++){
//     if (doDebugMessages_) std::cout << "objectsMatchedToPath TARGET = " << (*funct2)->pt() << " and eta " << (*funct2)->eta() << std::endl ;
//   }
//   
//   
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // now I need to find, among objectsMatchedToEndOfSeedingPath, the TAG object, i.e. the one which has passed the tight requirements 
//   // if an object in objectsMatchedToEndOfSeedingPath matches any of the object in electronSeedingPathTagFilter
//   // => it's tighter than a simple supercluster => it's the TIGHT electron of the ele-SC pair of the seeding path
//   // ==> this one is going to be my TAG 
//   pat::TriggerObjectRefVector objectsMatchedToSeedingPathTagFilter = findObjectsMatchedToPath(triggerEvent, seedingAcceptedPath, electronSeedingPathTagFilter_);
//   if(objectsMatchedToSeedingPathTagFilter.size()==0){
//     std::cout << "problem: found no objects matched to the TAG filter electronSeedingPathTagFilter: " << electronSeedingPathTagFilter_.at(0) 
//	       << " of the seeding path: " << seedingAcceptedPath->name() << ". This is a problem, bailing out. " << std::endl;
//     assert(0);
//   }
//   else
//     {  
//       if (doDebugMessages_) {std::cout << "found " << objectsMatchedToSeedingPathTagFilter.size()
//					<< " HLT objects matched to the TAG filter  " << electronSeedingPathTagFilter_.at(0)
//					<< " of the seeding path: " << seedingAcceptedPath->name() << std::endl ;}
//       eventsFate_->Fill(4);
//     }
//   
//
//   pat::TriggerObjectRefVector SeedingPathTagObjects;
//   for (pat::TriggerObjectRefVector::const_iterator seedPathEndObj = objectsMatchedToEndOfSeedingPath.begin(); 
//	seedPathEndObj !=objectsMatchedToEndOfSeedingPath.end(); seedPathEndObj++)
//     {
//       
//       for (pat::TriggerObjectRefVector::const_iterator seedPathTagObj = objectsMatchedToSeedingPathTagFilter.begin();
//	    seedPathTagObj !=objectsMatchedToTargetPath.end(); seedPathTagObj++)
//	 {
//	   
//	   if (  deltaR( (*seedPathEndObj)->eta(),  (*seedPathEndObj)->phi() , 
//			 (*seedPathTagObj)->eta(),  (*seedPathTagObj)->phi() ) 
//		 < 0.1 ) {
//	     
//	     // If I get here, I've found an objects among those of the seeding path which is tighter than a SC
//	     // => I'll chose this as a tag
//	     SeedingPathTagObjects.push_back( (*seedPathEndObj) );
//	     break;
//	     
//	   }
//	 }//loop over HLT objects from the objectsMatchedToTargetPath_ filter of the seeding path
//     }//loop over HLT objects from the electronSeedingPathEndFilter_ filter of the seeding path
//
//
//   if(SeedingPathTagObjects.size() ==0 ){
//     std::cout << "Problem: found no objects among those at the end of seeding path which match the output of the TAG filter. Bailing out." << std::endl; 
//     assert(0);
//   }
//   else
//     { 
//       if (doDebugMessages_) {
//	 std::cout << "found: " << SeedingPathTagObjects.size() << " objects at the end of seeding path which match the output of the TAG filter " << std::endl; }
//       eventsFate_->Fill(5);
//     }
//   
//   
//   
//   // randomly select a tag among the possible tags
//   unsigned int theSeedingPathTagIndex =  static_cast<unsigned int>(  (( iEvent.orbitNumber()%100 ) *1.0 / 100. ) * SeedingPathTagObjects.size()  ) ;
//   if (doDebugMessages_) {
//     std::cout << "there are " << SeedingPathTagObjects.size() << " HLT objects from the end of the seeding path matched to TAG filter of the seeding path; I've chosen: " << theSeedingPathTagIndex << " -th \t" << ( ( iEvent.orbitNumber()%100 ) *1.0 / 100. )  << std::endl ;
//   }
//   
//
//
//   // determine the PROBE object among the seeding path objects from the end filter
//   int theSeedingPathProbeIndex(-1);
//   for (pat::TriggerObjectRefVector::const_iterator seedPathEndObj = objectsMatchedToEndOfSeedingPath.begin(); 
//	seedPathEndObj !=objectsMatchedToEndOfSeedingPath.end(); seedPathEndObj++)
//     {
//       
//       if (  deltaR( (*seedPathEndObj)->eta(),    (*seedPathEndObj)->phi() , 
//		     SeedingPathTagObjects.at(theSeedingPathTagIndex)->eta(),  
//		     SeedingPathTagObjects.at(theSeedingPathTagIndex)->phi()  
//		     // objectsMatchedToEndOfSeedingPath.at(theSeedingPathTagIndex)->eta(),  // ==> WRONG!
//		     // objectsMatchedToEndOfSeedingPath.at(theSeedingPathTagIndex)->phi()   // ==> WRONG!
//		     ) 
//	     > 0.1 ) {
//	 
//	 
//	 theSeedingPathProbeIndex = ( seedPathEndObj-objectsMatchedToEndOfSeedingPath.begin() );
//	 break;
//	 // stop after you've found the first possible PROBE
//
//       }// if object is distinct from the TAG
//     }// loop over seeding HLT object, in order to determined the probe
//   
//   
//
//   // do the sanity check: there needs to be a PROBE
//   if (theSeedingPathProbeIndex==-1)
//     {
//       // this case should never happen: there needs to be at least a second object from the seeding path distinct from the TAG object
//       std::cout << "no second object from HLT path: " << seedingAcceptedPath->name() 
//		 << " was found to act as a probe. Problem. Bailig out." << std::endl;
//       return;
//     }
//   else
//     {
//       std::cout << "there are " << objectsMatchedToEndOfSeedingPath.size() << " HLT objects from the end of the seeding path;  I've chosen the number: " << theSeedingPathProbeIndex << " to be the probe online object."  << std::endl ;
//       eventsFate_->Fill(5);
//     }
//   
//
//
//   // now that you have the PROBE HLT object from from the end of  the seeding path, match it to one of the offline objects
//   // there must be a matching since this is a requirement earlier, in the selection in order to enter the denominator
//   int theOfflineEleMatchingProbeIndex(-1);
//   for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
//     {
//       
//       if (  deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
//		     objectsMatchedToEndOfSeedingPath.at(theSeedingPathProbeIndex)->eta(),  
//		     objectsMatchedToEndOfSeedingPath.at(theSeedingPathProbeIndex)->phi() ) 
//	     < 0.1 ) {
//	 theOfflineEleMatchingProbeIndex = theEl;
//	 break;
//       }// if there's matching to probe HLT obhect from end of seeding path
//     }//loop over offline objects    
//
//
//   int theOfflineEleMatchingTagIndex(-1);
//   for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
//     {
//       
//       if (  deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
//		     objectsMatchedToEndOfSeedingPath.at(theSeedingPathTagIndex)->eta(),  
//		     objectsMatchedToEndOfSeedingPath.at(theSeedingPathTagIndex)->phi() ) 
//	     < 0.1 ) {
//	 theOfflineEleMatchingTagIndex = theEl;
//	 break;
//       }// if there's matching to tag HLT obhect from end of seeding path
//     }//loop over offline objects    
//   
//   
//
//   // do the sanity check: there need to be an offline object matching both the PROBE and the TAG
//   if (theOfflineEleMatchingProbeIndex==-1 || theOfflineEleMatchingTagIndex==-1)
//     {
//       std::cout << "no offline object fuond which can matche either the probe or the tag (" 
//		 << theOfflineEleMatchingProbeIndex << "\t" << theOfflineEleMatchingTagIndex   
//		 << ") among online objects from the seeding path: " << seedingAcceptedPath->name() 
//		 << " Problem. Bailig out." << std::endl;
//       assert(0);
//     }
//   else
//     {
//       if (doDebugMessages_) { std::cout << "offline electrons found which matches both the online PROBE and the online TAG object for the path: "
//					 << seedingAcceptedPath->name()  << std::endl ; }
//       eventsFate_->Fill(6);
//     }
//   
//   // ALL THE WAY TO HERE

   
   
   //pTTagDenominator_    -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
   //pTProbeDenominator_  -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
   // now make fill the plots for the the numerator! Has the HLT target bit PASSED?
   if ( targetHLTPathsAccepted.size()>0) {
     counterEvtsWithTargetFired_++;
     eventsFate_->Fill(7);
     massNumerator_   -> Fill( ( theEleOfflineCands.at(0).first.p4() +  theEleOfflineCands.at(1).first.p4() ).M() ); 
     pTele1Numerator_ -> Fill( theEleOfflineCands.at(0).first.pt() );
     pTele2Numerator_ -> Fill( theEleOfflineCands.at(1).first.pt() );
     //pTTagNumerator_     -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );  
     //     pTProbeNumerator_   -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
    }
   else {
     // if target PATH has not passed, fill histograms for FAIL
     massFail_   -> Fill( ( theEleOfflineCands.at(0).first.p4() +  theEleOfflineCands.at(1).first.p4() ).M() ); 
     pTele1Fail_ -> Fill( theEleOfflineCands.at(0).first.pt() );
     pTele2Fail_ -> Fill( theEleOfflineCands.at(1).first.pt() );
     // pTTagFail_   -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
     //pTProbeFail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
    return;
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

  if (doDebugMessages_) {
    std::cout << "\nLooking for trigger objects involved in this path: " 
	      <<  seedingAcceptedPath->name() ;
    if(electronFilters.size()>0) std::cout << " pertaining this filter: " << electronFilters.at(0) ;
      std::cout  << std::endl ; 
  }
 
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
      if (doDebugMessages_)  std::cout << "\n==> found filter which is among electronFilters: " << filterRef->label() << std::endl;
      
      
      for ( pat::TriggerObjectRefVector::const_iterator iobjRef = objectsInPath.begin(); iobjRef != objectsInPath.end(); ++iobjRef ) {
	pat::TriggerObjectRef objRef = *iobjRef ; 
	if (doDebugMessages_) std::cout << "Found an object with pT = " << objRef->pt() << " and eta " << objRef->eta() << std::endl ; 
	
	if ( triggerEvent->objectInFilter(objRef,filterRef->label()) ) { // Trigger object was used by the filter
	  if (doDebugMessages_)  {std::cout << "found candidate requested IN FILTER: " << filterRef->label() 
					    << " object with pt: " <<  objRef->pt() 
					    << " and phi " << objRef->phi()
					    << " and eta " << objRef->eta() << std::endl;   }
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






















   float a;

//#ifndef momentarilyIngnore
//   if (doDebugMessages_) std::cout << "seedHLTForL1Efficiency has size: " << seedHLTForL1Efficiency_.size() << std::endl;
//
//
//   for (unsigned int i=0; i<seedHLTForL1Efficiency_.size(); i++) { 
//
//     // std::cout << "Investigating path: " << seedHLTForL1Efficiency_.at(i) << std::endl ; 
//     const pat::TriggerPathRef iPath = triggerEvent->pathRef( seedHLTForL1Efficiency_.at(i) ) ; 
//     if ( iPath.isNonnull() && iPath->wasAccept() ) { 
//       std::cout << "Found path!" << std::endl ; 
//       
//       std::cout << "Path information: " 
//		 << " " << iPath->name() << " " 
//		 << " " << iPath->index() << " " 
//		 << " pr: " << iPath->prescale() << " " 
//		 << " r: " << iPath->wasRun() << " " 
//		 << " acc: " << iPath->wasAccept() << " " 
//		 << " " << iPath->wasError() << " " 
//		 << " last: " << iPath->lastActiveFilterSlot() << " " 
//		 << " " << iPath->modules().size() ; 
//       std::cout // << iPath->modules().at(iPath->lastActiveFilterSlot()) << " "
//	 << iPath->l3Filters() << " " 
//	 << iPath->xTrigger() << " " ; 
//       // for (unsigned int j=0; j<iPath->filterIndices().size(); j++) { 
//       //   std::cout << iPath->modules().at(iPath->filterIndices().at(j)) << " " ; 
//       // }
//       std::cout << " \n\n\n\n" 
//		 << std::endl ; 
//       
//       std::cout << "Looking for trigger objects involved in this path: " <<  iPath->name() << std::endl ; 
// 
//       pat::TriggerObjectRefVector objectsInPath = triggerEvent->pathObjects(iPath->name(),true) ;
//       pat::TriggerFilterRefVector filtersInPath = triggerEvent->pathFilters(iPath->name(),true) ;
//       std::cout  << "size of objectsInPath : " << objectsInPath.size() << std::endl; 
//
//       
//       // http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/T>> Entering Package PhysicsTools/PatAlgos
//       // // Get a reference to a certain L1 algorithm by name,
//       // // NULL, if algorithm is not found
//       //    const TriggerAlgorithmRef algorithmRef( const std::string & nameAlgorithm ) const;
//
//
//       // http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/TriggerEvent.h#308
//       // 
//       // TriggerObjectRefVector algorithmObjects( const std::string & nameAlgorithm ) const;
//
//       //       http://cmslxr.fnal.gov/lxr/source/DataFormats/PatCandidates/interface/TriggerEvent.h#310 
//       //       /// Checks, if an object was used in a certain algorithm given by name
//       //       bool objectInAlgorithm( const TriggerObjectRef & objectRef, const std::string & nameAlgorithm ) const;
//
//
//     }
//   }
//
//   pat::TriggerAlgorithmRefVector theAlgosActive = triggerEvent->algorithmRefs();
//   std::cout << "size of theAlgosActive: " << theAlgosActive.size() << std::endl;
//
//   pat::TriggerAlgorithmRefVector theAlgosPassed = triggerEvent->acceptedAlgorithms();
//   std::cout << "size of theAlgosPassed: " << theAlgosPassed.size() << std::endl;
//
//   pat::TriggerAlgorithmRefVector thePhysicsAlgosPassed = triggerEvent->acceptedPhysAlgorithms();
//   std::cout << "size of thePhysicsAlgosPassed: " << thePhysicsAlgosPassed.size() << std::endl;
//
//   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTrigger#PATTriggerProducer
//   std::cout << "physics algos which passed: " << std::endl;
//   for(unsigned int phAl=0; phAl<thePhysicsAlgosPassed.size(); phAl++)
//     {
//       std::cout << thePhysicsAlgosPassed.at(phAl)->name() << "\n";
//       pat::TriggerObjectRefVector objectsInAlgo =  triggerEvent->algorithmObjects( thePhysicsAlgosPassed.at(phAl)->name() );
//       // std::cout << "\t\t\t\t size of objectsInAlgo: " << objectsInAlgo.size() << std::endl;
//     }
//   std::cout << std::endl;
//
//
//
//   // Isolated EM particles
//   Handle< l1extra::L1EmParticleCollection > isoEmColl ;
//   iEvent.getByLabel( isolatedEmSource_, isoEmColl ) ;
//   std::cout << "++ HeavyNuEleTriggerEff: Number of isolated EM " << isoEmColl->size() << std::endl ;
//   
//   for( l1extra::L1EmParticleCollection::const_iterator emItr = isoEmColl->begin() ;
//	emItr != isoEmColl->end() ;
//	++emItr )
//     {
//       std::cout << "\t  p4 (" << emItr->px()
//	    << ", " << emItr->py()
//	    << ", " << emItr->pz()
//	    << ", " << emItr->energy()
//	    << ") et " << emItr->et()
//	    << " eta " << emItr->eta()
//	    << " phi " << emItr->phi()
//		 << std::endl ;
//     }
//   
//   // Non-isolated EM particles
//   Handle< l1extra::L1EmParticleCollection > nonIsoEmColl ;
//   iEvent.getByLabel( nonIsolatedEmSource_, nonIsoEmColl ) ;
//   std::cout << "++ HeavyNuEleTriggerEff: Number of non-isolated EM " << nonIsoEmColl->size() << std::endl ;
//   
//   for( l1extra::L1EmParticleCollection::const_iterator emItr = nonIsoEmColl->begin() ;
//	emItr != nonIsoEmColl->end() ;
//	++emItr )
//     {
//       std::cout << "\t\t  p4 (" << emItr->px()
//	    << ", " << emItr->py()
//	    << ", " << emItr->pz()
//	    << ", " << emItr->energy()
//	    << ") et " << emItr->et()
//	    << " eta " << emItr->eta()
//	    << " phi " << emItr->phi()
//		 << std::endl ;
//     }
//   
//
//#endif

