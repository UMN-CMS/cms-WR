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
// $Id: HeavyNuEleTriggerEff.cc,v 1.13 2012/06/15 11:48:45 franzoni Exp $
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

#include "TH1F.h"
#include "TH2F.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

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

  pat::TriggerObjectRefVector findObjectsMatchedToPathNew( const edm::Handle< pat::TriggerEvent > triggerEvent, const pat::TriggerPathRef seedingAcceptedPath,  const std::vector<std::string>& electronFilters );
  
  
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

  edm::InputTag trigEventPatTag_;
  edm::InputTag trigEventHLTTag_;
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

  unsigned int counterExecutions_;
  unsigned int counterEvtsAll_;
  unsigned int counterEvtsPassingOfflineSelection_;
  unsigned int counterEvtsWithSeedingFired_;
  unsigned int counterEvtsWithSeedingAbovePt_;
  unsigned int counterEvtsWithTargetValid_;
  unsigned int counterEvtsWithSeedingPathProbesMatchedToOfflineObj_;
  unsigned int counterEvtsWithTagFound_;
  unsigned int counterEvtsWithProbeFound_;
  unsigned int counterEvtsTangProbeMatchedOffl_;
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

  TH2F*  ptProbeVSmassDenom_;
  TH2F*  ptProbeVSmassNum_;
  TH2F*  ptProbeVSmassFail_;
  TH2F*  etaProbeVSmassDenom_;
  TH2F*  etaProbeVSmassNum_;
  TH2F*  etaProbeVSmassFail_;
  

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
  trigEventPatTag_     ( iConfig.getParameter< edm::InputTag > ( "trigEventPatTag" ) ),
  trigEventHLTTag_     ( iConfig.getParameter< edm::InputTag > ( "trigEventHLTTag" ) ),
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
  counterExecutions_(0),
  counterEvtsAll_ (0),
  counterEvtsPassingOfflineSelection_ (0),
  counterEvtsWithSeedingFired_ (0),
  counterEvtsWithSeedingAbovePt_(0),
  counterEvtsWithTargetValid_ (0),
  counterEvtsWithSeedingPathProbesMatchedToOfflineObj_(0),
  counterEvtsWithTagFound_(0),
  counterEvtsWithProbeFound_(0),
  counterEvtsTangProbeMatchedOffl_(0),
  counterEvtsWithTargetFired_ (0),
  plotFolderName_                ( iConfig.getParameter< std::string > ( "plotFolderName" ) )
{

  //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  monitorRuns_       = (TH1F*)   fs->make<TH1F>("monitorRuns","monitorRuns; run (starting coll 2012)",20000,firstRunForPlotting_,firstRunForPlotting_+20000);
  massDenominator_   = (TH1F*)   fs->make<TH1F>("massDenominator","massDenominator; m(ee) [GeV]",2000,0.,2000);
  massNumerator_     = (TH1F*)   fs->make<TH1F>("massNumerator","massNumerator; m(ee) [GeV]",2000,0.,2000);
  massFail_          = (TH1F*)   fs->make<TH1F>("massFail","massFail; m(ee) [GeV]",2000,0.,2000);
  pTele1Denominator_ = (TH1F*)   fs->make<TH1F>("pTele1Denominator","pTele1Denominator; p_{T}(e1) [GeV]",2000,0.,2000.);
  pTele1Numerator_   = (TH1F*)   fs->make<TH1F>("pTele1Numerator","pTele1Numerator; p_{T}(e1) [GeV]",2000,0.,2000);
  pTele1Fail_        = (TH1F*)   fs->make<TH1F>("pTele1Fail","pTele1Fail; p_{T}(e1) [GeV]",2000,0.,2000);
  pTele2Denominator_ = (TH1F*)   fs->make<TH1F>("pTele2Denominator","pTele2Denominator; p_{T}(e2) [GeV]",2000,0.,2000.);
  pTele2Numerator_   = (TH1F*)   fs->make<TH1F>("pTele2Numerator","pTele2Numerator; p_{T}(e2) [GeV]",2000,0.,2000);
  pTele2Fail_        = (TH1F*)   fs->make<TH1F>("pTele2Fail","pTele2Fail; p_{T}(e1) [GeV]",2000,0.,2000);
  eventsFate_        = (TH1F*)   fs->make<TH1F>("eventsFate","eventsFate; the events fate",10,0.,10);
  
  pTTagDenominator_   = (TH1F*)   fs->make<TH1F>("pTTagDenominator","pTTagDenominator; p_{T}(tag,denom) [GeV]",2000,0.,2000);
  pTTagNumerator_     = (TH1F*)   fs->make<TH1F>("pTTagNumerator","pTTagNumerator; p_{T}(tag,num) [GeV]",2000,0.,2000);
  pTTagFail_          = (TH1F*)   fs->make<TH1F>("pTTagFail","pTTagFail; p_{T}(tag,fail) [GeV]",2000,0.,2000);
  pTProbeDenominator_ = (TH1F*)   fs->make<TH1F>("pTProbeDenominator","pTProbeDenominator; p_{T}(probe,denom) [GeV]",2000,0.,2000);
  pTProbeNumerator_   = (TH1F*)   fs->make<TH1F>("pTProbeNumerator","pTProbeNumerator; p_{T}(probe,num) [GeV]",2000,0.,2000);
  pTProbeFail_        = (TH1F*)   fs->make<TH1F>("pTProbeFail","pTProbeFail; p_{T}(probe,fail) [GeV]",2000,0.,2000);
  
  double yAxisPt[8] = {0, 40., 50., 60., 80., 100., 200., 2000};
  ptProbeVSmassDenom_ = (TH2F*)   fs->make<TH2F>("ptProbeVSmassDenom","ptProbeVSmassDenom; m(ee) [GeV]; p_{T}(probe,fail) [GeV]",200,30.,230,7,yAxisPt);
  ptProbeVSmassNum_   = (TH2F*)   fs->make<TH2F>("ptProbeVSmassNum","ptProbeVSmassNum; m(ee) [GeV]; p_{T}(probe,fail) [GeV]",200,30.,230,7,yAxisPt);
  ptProbeVSmassFail_  = (TH2F*)   fs->make<TH2F>("ptProbeVSmassFail","ptProbeVSmassFail; m(ee) [GeV]; p_{T}(probe,fail) [GeV]",200,30.,230,7,yAxisPt);

  double yAxisEta[11] = {-2.5, -2, -1.6, -1.4, -1., 0., 1., 1.4, 1.6, 2., 2.5};
  etaProbeVSmassDenom_ = (TH2F*)   fs->make<TH2F>("etaProbeVSmassDenom","etaProbeVSmassDenom; m(ee) [GeV]; #eta(probe,fail) [GeV]",200,30.,230,10,yAxisEta);
  etaProbeVSmassNum_   = (TH2F*)   fs->make<TH2F>("etaProbeVSmassNum","etaProbeVSmassNum; m(ee) [GeV]; #eta(probe,fail) [GeV]",200,30.,230,10,yAxisEta);
  etaProbeVSmassFail_  = (TH2F*)   fs->make<TH2F>("etaProbeVSmassFail","etaProbeVSmassFail; m(ee) [GeV]; #eta(probe,fail) [GeV]",200,30.,230,10,yAxisEta);
  

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
  
  // count how many times analyzer was executed 
  counterExecutions_++;
  eventsFate_->Fill(0);

  // allow local run selection
  if (iEvent.run() < runMin_ || runMax_ < iEvent.run()) return;


  // warning in case of unusual run range. But not an error.
  if (firstRunForPlotting_ > iEvent.run() && verbosityForRunPlotting_)
    {
      std::cout << "PROBLEM: in HeavyNuEleTriggerEff there's an event from run: " << iEvent.run() 
		<< " while only runs starting from: " << firstRunForPlotting_ 
		<< " are properly monitored; will continue, and not repeat this message" << std::endl;
	verbosityForRunPlotting_ = false; 
    }
  
  monitorRuns_->Fill( iEvent.run() );
  counterEvtsAll_++;
  eventsFate_->Fill(1);


  // require offline selections:
  // - 2 electrons passing HEEP, above 40 GeV
  // - a numer of jets which is configurable ( typically 0 or 2 ) pass jet ID
  // - the two electrons are 'far enough' from the jets
  // if these requirements fail, don't consider the event

  std::vector< std::pair<pat::Electron, float> > theEleOfflineCands;
  std::vector< std::pair<pat::Jet,      float> > jetCands;

  if ( passOfflineSelection( iEvent , iSetup, theEleOfflineCands, jetCands) ){
    if (doDebugMessages_) std::cout << "passOfflineSelection returns: " << passOfflineSelection( iEvent , iSetup, theEleOfflineCands, jetCands) << std::endl;
    counterEvtsPassingOfflineSelection_ ++;
    eventsFate_->Fill(2);
  }
  else {
    if (doDebugMessages_) { std::cout << "(in analyze) size of theEleOfflineCands which pass the offline selections: " 
				      << theEleOfflineCands.size() << " and size of jet coll is: " << jetCands.size() 
				      << std::endl; }

    // if basic offline conditions are not satisfied, move on to the next event
    return;
  }  

  

  // get all sorts of information about trigger, from pat
  edm::Handle< pat::TriggerEvent > triggerPatEvent;
  iEvent.getByLabel( trigEventPatTag_, triggerPatEvent );
   if ( !triggerPatEvent.isValid() ) {
     std::cout << "triggerPatEvent not found " << std::endl;
     assert(0);
   }
   else { if(doDebugMessages_ ) std::cout << " triggerPatEvent was fuond " << std::endl;}
   
   //get trigger HLT event
   edm::Handle<trigger::TriggerEvent> trigEvent; 
   iEvent.getByLabel(trigEventHLTTag_,trigEvent);
   if ( !trigEvent.isValid() ) {
     std::cout << "trigHLTEvent not found. Bailing out. " << std::endl;
     assert(0);
   }
   else { if(doDebugMessages_ ) std::cout << " trigHLTEvent was fuond " << std::endl;}




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
   pat::TriggerPathRefVector theAcceptedPaths =  triggerPatEvent->acceptedPaths();

   // check if any of the seeding / target HLT paths has fired
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
     eventsFate_->Fill(3);
   }
   else                                            return; 





   //////////////////////////// prepare target trigger path /////////////////
   if (targetHLTPaths_.size()==0) {
     std::cout << "you've not provided targetHLTPaths; bailing out" << std::endl;
     assert(0);
   }   // at least one target HLT path must be provided; if not, bail out

   //check if any of the provided target triggers are in the menu
   std::vector <std::string> targetHLTPathsInMenu;
   pat::TriggerPathRefVector thePathsInMenu =  triggerPatEvent->pathRefs();
   for (unsigned int i=0; i<targetHLTPaths_.size(); i++) { 
     for ( pat::TriggerPathRefVector::const_iterator iPath = thePathsInMenu.begin(); iPath != thePathsInMenu.end(); ++iPath ) {
       if ( targetHLTPaths_.at(i) == (*iPath)->name() ) 	 { targetHLTPathsInMenu.push_back( (*iPath)->name() ); }
     }// loop over active paths
   }// loop over HLT-eff seeding paths


   // if there's not at least one target path which is validly in the menu, stop execution
   if ( targetHLTPathsInMenu.size()>0 ) {
     counterEvtsWithTargetValid_++;
     eventsFate_->Fill(4);
   }
   else                                 return;

  
   if (doDebugMessages_) {
     std::cout << "number of active HLT paths: "       << thePathsInMenu.size() << std::endl; // why do I care?
     std::cout << "seedHLTForHLTEfficiency has size: " << seedHLTForHLTEfficiency_.size() 
	       << " of which  "                        << seedHLTForHLTEfficiencyAccepted.size() << " were 'run in the menu' && 'PASSED'" << std::endl;
     std::cout << "targetHLTPaths.at(0) is: "          << targetHLTPaths_.at(0) << std::endl;
     std::cout << "targetHLTPaths size: "              << targetHLTPaths_.size() 
	       << " of which: "                        << targetHLTPathsInMenu.size() 
	       << " were run in the menu (they may have passed or failed)" << std::endl;
   }
   


   ////////////////////////////////////////////////////////////////////////////////
   // in summary, at this point: do the HLT efficnecy study only if there's at least:
   // - a seeding path which was accepted
   // - a target path which was specified in input by the user and is in the menu
   // ==> conditions are guaranteed by the return statements of above
   



   
   // one of the seding HLT paths has fired. Chose the first one as a base trigger to do the efficiency
   const pat::TriggerPathRef seedingAcceptedPath = triggerPatEvent->pathRef( seedHLTForHLTEfficiencyAccepted.at(0) ) ; 
   
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
     std::cout       << seedingAcceptedPath->l3Filters() << " " 
		     << seedingAcceptedPath->xTrigger() << " " ; 
   }

   


   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // get the HLT objects which match the LAST filter of the seeding HLT path  => I'll build the TAG and the PROBE out of these
   std::vector<math::XYZTLorentzVector> trigObjP4sLastFilterSeedPath;
   // process name comes from the trigger event process as they have to be the same...
   trigtools::getP4sOfObsPassingFilter(trigObjP4sLastFilterSeedPath,*trigEvent,electronSeedingPathEndFilter_.at(0),trigEventHLTTag_.process()); 
   if (doDebugMessages_) 
     { 
       std::cout << "\n\nfound: "                                    << trigObjP4sLastFilterSeedPath.size() 
		 <<" (HLT) objects matching the SEED last filter: "<< electronSeedingPathEndFilter_.at(0) <<"\n" << std::endl; 

       for(size_t objNr=0;objNr<trigObjP4sLastFilterSeedPath.size();objNr++){
	 std::cout <<"     et obj "<<objNr<<" : "<<trigObjP4sLastFilterSeedPath[objNr].Et() 
		   << "\t eta: " << trigObjP4sLastFilterSeedPath[objNr].Eta() 
		   << "\t phi: " << trigObjP4sLastFilterSeedPath[objNr].Phi() 
		   <<std::endl;
       }// loop over HLT objects passing last filter of HLT seed path

     }// if debug

   

   // do the objects which fired electronSeedingPathEndFilter match any of the offline electrons ?
   unsigned int counterSeedingPathHLTMatchedObjects(0);
   for(size_t objNr=0; objNr<trigObjP4sLastFilterSeedPath.size(); objNr++)
     {
       for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
	 {
	   
	   if (  deltaR( theEleOfflineCands.at(theEl).first.eta(), 
			 theEleOfflineCands.at(theEl).first.phi() , 
			 trigObjP4sLastFilterSeedPath[objNr].Eta(),
			 trigObjP4sLastFilterSeedPath[objNr].Phi()
			 )  < 0.1 ) 
	     {
	       counterSeedingPathHLTMatchedObjects++;
	       break;
	     }// if there's matching with an offline electron
	   
	 }// loop over offline electrons 
     }//loop over HLT objects from the electronSeedingPathEndFilter_ filter of the seeding path



   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // require that _all_ the online objects which fired the last filter of SEED path be matched to an offline HEEP electron
   // if one of them is not matched an  offline object,
   // reject the event from the sample of the measurement (reduction of background!)
   if(counterSeedingPathHLTMatchedObjects < trigObjP4sLastFilterSeedPath.size() ){
     if (doDebugMessages_) std::cout << "only: " << counterSeedingPathHLTMatchedObjects 
				     << " hlt objects out of: " << trigObjP4sLastFilterSeedPath.size()
				     << " are matched to the : " << theEleOfflineCands.size()
				     << " offline HEEP electrons. Rejecting event from test sample."
				     << std::endl;
     return;
   } 
   else{
     counterEvtsWithSeedingPathProbesMatchedToOfflineObj_++;
     eventsFate_->Fill(5);
   }
   


   

   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // now I need to find, among objectsMatchedToEndOfSeedingPath, the TAG object, 
   // i.e. the one which has passed the tight requirements 
   // if an object in objectsMatchedToEndOfSeedingPath matches any of the object in electronSeedingPathTagFilter
   // => it's tighter than a simple supercluster => it's the TIGHT electron of the ele-SC pair of the seeding path
   // ==> this one is going to be my TAG 

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // get the HLT objects which match the TAG filter of the seeding HLT path
   std::vector<math::XYZTLorentzVector> trigObjP4sTagFilterSeedPath;
   // process name comes from the trigger event process as they have to be the same...
   trigtools::getP4sOfObsPassingFilter(trigObjP4sTagFilterSeedPath,*trigEvent,electronSeedingPathTagFilter_.at(0),trigEventHLTTag_.process()); 
   if (doDebugMessages_) 
     { 
       std::cout << "\nfound: "<< trigObjP4sTagFilterSeedPath.size() <<" (HLT) objects matching the TAG SEED filter: "
		 << electronSeedingPathTagFilter_.at(0) <<"\n" << std::endl; 
       for(size_t objNr=0;objNr<trigObjP4sTagFilterSeedPath.size();objNr++){
	 std::cout <<"     et obj "<<objNr<<" : "<<trigObjP4sTagFilterSeedPath[objNr].Et() 
		   << "\t eta: " << trigObjP4sTagFilterSeedPath[objNr].Eta() 
		   << "\t phi: " << trigObjP4sTagFilterSeedPath[objNr].Phi() 
		   <<std::endl;
       }// loop over HLT objects passing TAG filter of HLT seed path
     }// if debug



   // find one (possible) tag among the objects which passed the end-filter  of the seeding path
   // match objects from 'end-of-path filter' to objects that have passed the 'tag filter';
   std::vector<math::XYZTLorentzVector> SeedingPathTagObjectsHLT;
   for(size_t objNrLast=0; objNrLast< trigObjP4sLastFilterSeedPath.size(); objNrLast++){

     for(size_t objNrTag=0; objNrTag< trigObjP4sTagFilterSeedPath.size(); objNrTag++){


       if (doDebugMessages_) {
	 std::cout << "comparing LastFilterSeedPath : " 
		   << trigObjP4sLastFilterSeedPath[objNrLast].Et() 
		   << "\t" << trigObjP4sLastFilterSeedPath[objNrLast].Eta() << "\t" <<   trigObjP4sLastFilterSeedPath[objNrLast].Phi()
		   << "\nto TagFilterSeedPath      : " 
		   << trigObjP4sTagFilterSeedPath [objNrTag] .Et() << "\t" 
		   << trigObjP4sTagFilterSeedPath [objNrTag] .Eta() << "\t" <<   trigObjP4sTagFilterSeedPath [objNrTag] .Phi()
		   << "\n deltaR: " 
		   << deltaR( trigObjP4sLastFilterSeedPath[objNrLast].Eta(),  trigObjP4sLastFilterSeedPath[objNrLast].Phi() ,
			      trigObjP4sTagFilterSeedPath [objNrTag] .Eta(),  trigObjP4sTagFilterSeedPath [objNrTag] .Phi()  )
		   << std::endl;
       }//debug
       
       if (  deltaR( trigObjP4sLastFilterSeedPath[objNrLast].Eta(),  trigObjP4sLastFilterSeedPath[objNrLast].Phi() , 
		     trigObjP4sTagFilterSeedPath [objNrTag] .Eta(),  trigObjP4sTagFilterSeedPath [objNrTag] .Phi()  ) 
	     < 0.1 ) 
	 {
	   // If I get here, I've found one object among those of the seeding path end which is tighter than a SC
	   // => this is one of the possible tags (I'll randomize later to choose a specific tag, if more than one tag is found)
	   SeedingPathTagObjectsHLT.push_back(  trigObjP4sLastFilterSeedPath[objNrLast]  );
	   break;
	 }// if there's matching
       
     }// loop over HLT objects passing TAG filter of HLT seed path
   }// loop over HLT objects passing LAST filter of HLT seed path
   
   if(SeedingPathTagObjectsHLT.size() ==0 ){
     std::cout << "Problem: found no objects (SeedingPathTagObjectsHLT) among those at the end of seeding path which match the output of the TAG filter."
	       << " Discardig event." << std::endl; 
     //assert(0);
     return;
   }
   else
     { 
       if (doDebugMessages_) {
	 std::cout << "found: " << SeedingPathTagObjectsHLT.size() << " objects at the end of seeding path which match the output of the seeding TAG filter;"
		   << " these are the TAG candidates. " << std::endl; }
       counterEvtsWithTagFound_++;
       eventsFate_->Fill(6);
     }
   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // randomly choose THE tag among the possible tags
   // reproducible random selection based on orbit number
   unsigned int theSeedingPathTagHLTIndex =  static_cast<unsigned int>(  (( iEvent.orbitNumber()%100 ) *1.0 / 100. ) * SeedingPathTagObjectsHLT.size()  ) ;
   if (doDebugMessages_) {
     std::cout << "there are " << SeedingPathTagObjectsHLT.size() << " HLT objects from the end of the seeding path matched to TAG filter of the seeding path; I've chosen: " << theSeedingPathTagHLTIndex << " -th \t" << ( ( iEvent.orbitNumber()%100 ) *1.0 / 100. )  << std::endl ;
     
   }

   


   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // choose the PROBE among the end-of-path objects of the seed path, which is away from the chosen TAG
   int theSeedingPathProbeHLTIndex(-1);
   for(size_t objNrP=0; objNrP< trigObjP4sLastFilterSeedPath.size(); objNrP++){
     
     if (  deltaR( trigObjP4sLastFilterSeedPath[objNrP].Eta(),
		   trigObjP4sLastFilterSeedPath[objNrP].Phi(),
		   SeedingPathTagObjectsHLT[theSeedingPathTagHLTIndex].Eta(),  
		   SeedingPathTagObjectsHLT[theSeedingPathTagHLTIndex].Phi()  
		   ) 
	   > 0.1 ) {
       theSeedingPathProbeHLTIndex = objNrP;
       break;
       // stop after you've found the first possible PROBE   => use that as probe
     }// if probe candidate is separate from chosen TAG
   }// loop over probe candidates
   

   // do the sanity check: there needs to be a PROBE
   if (theSeedingPathProbeHLTIndex==-1)
     {
       // this case should never happen: there needs to be at least a second object from the seeding path last filter distinct from theT AG object
       std::cout << "no second object from HLT path: " << seedingAcceptedPath->name() 
		 << " was found to act as a probe. Problem. Skipping event." << std::endl;
       return;
     }
   else
     {
       if (doDebugMessages_) {
       std::cout << "PROBE (HLT) found; there are " << trigObjP4sLastFilterSeedPath.size() << " HLT objects from the end of the seeding path; "
		 <<" I've chosen the number: " << theSeedingPathProbeHLTIndex << " to be the probe online object."  << std::endl ;

       std::cout << "the PROBE is: " 
		 << trigObjP4sLastFilterSeedPath[theSeedingPathProbeHLTIndex].Et() 
		 << "\t" << trigObjP4sLastFilterSeedPath[theSeedingPathProbeHLTIndex].Eta() 
		 << "\t" << trigObjP4sLastFilterSeedPath[theSeedingPathProbeHLTIndex].Phi()
		 << std::endl;
       }
       counterEvtsWithProbeFound_++;       
       eventsFate_->Fill(7);
     }
   math::XYZTLorentzVector theTag   =  SeedingPathTagObjectsHLT[theSeedingPathTagHLTIndex];
   math::XYZTLorentzVector theProbe =  trigObjP4sLastFilterSeedPath[theSeedingPathProbeHLTIndex];



   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // match tag and probe hlt objects to offline objects
   // the object which fired the last filter of the seeding path have already been matched to the offline objects =>
   // matching hlt objects to offline is guaranteed
   
   int theOfflineEleMatchingProbeIndex(-1);
   for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
     {
       
       if (  
	   deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
		   theProbe.Eta(),  		     theProbe.Phi())  < 0.1  
	   && 
	   deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
		   theTag.Eta(),  		     theTag.Phi())    > 0.1  	     
	   ) 
	 {
	   theOfflineEleMatchingProbeIndex = theEl;
	   break;
	 }// if there's matching to probe HLT obhect from end of seeding path
     }//loop over offline objects    


   int theOfflineEleMatchingTagIndex(-1);
   for (unsigned  int theEl=0; theEl< theEleOfflineCands.size(); theEl++ )
     {
       
       if (  
	   deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
		   theTag.Eta(),  theTag.Phi())        < 0.1  
	   && 
	   deltaR( theEleOfflineCands.at(theEl).first.eta(), theEleOfflineCands.at(theEl).first.phi() , 
		   theProbe.Eta(),  theProbe.Phi())    > 0.1  	     
	   ) 
	 {
	   theOfflineEleMatchingTagIndex = theEl;
	   break;
	 }// if there's matching to tag HLT obhect from end of seeding path
     }//loop over offline objects    


   // do the sanity check: there need to be an offline object matching both the PROBE and the TAG
   if (theOfflineEleMatchingProbeIndex==-1 || theOfflineEleMatchingTagIndex==-1)
     {
       std::cout << "no offline object fuond which can matche either the probe or the tag (" 
		 << theOfflineEleMatchingProbeIndex << "\t" << theOfflineEleMatchingTagIndex   
		 << ") among online objects from the seeding path: " << seedingAcceptedPath->name() 
		 << " Problem. Bailig out." << std::endl;
       assert(0);
     }
   else
     {
       if (doDebugMessages_) { std::cout << "offline electrons found which matches both the online PROBE and the online TAG object for the path: "
					 << seedingAcceptedPath->name()  << std::endl ; }
       counterEvtsTangProbeMatchedOffl_++;
       eventsFate_->Fill(8);
     }



   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // get the HLT objects which match the LAST filter of the TARGET HLT path
   // will be needed later to decide if the probe has passed or not 
   std::vector<math::XYZTLorentzVector> trigObjP4sLastFilterTargetPath;
   // process name comes from the trigger event process as they have to be the same...
   trigtools::getP4sOfObsPassingFilter(trigObjP4sLastFilterTargetPath,*trigEvent,electronTargetFilters_.at(0),trigEventHLTTag_.process()); 
   if (doDebugMessages_) 
     { 
       
       std::cout << "found: "<< trigObjP4sLastFilterTargetPath.size() 
		 <<" (HLT) objects matching the TARGET filter: "<< electronTargetFilters_.at(0) <<"\n" << std::endl; 

       for(size_t objNr=0;objNr<trigObjP4sLastFilterTargetPath.size();objNr++){
	 std::cout <<"     et obj "<<objNr
		   <<" : "<<trigObjP4sLastFilterTargetPath[objNr].Et() 
		   << "\t" << trigObjP4sLastFilterTargetPath[objNr].Eta() 
		   << "\t" << trigObjP4sLastFilterTargetPath[objNr].Phi() 
		   <<std::endl;
       }// loop over HLT objects passing last filter of HLT target path
       
     }// if debug


   
   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // if an event passes this point, it has entered the _denominator_
   // i.e. I have two online objects which have  fired the seeding trigger and which are matched to offline objects
   // how do I chose which of the offline candidates are to be used? For now, chose the first two in the collection... REVISIT? 
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   massDenominator_    -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				.M() ); 
   pTTagDenominator_   -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
   pTProbeDenominator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
   
   ptProbeVSmassDenom_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				.M(),
				theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() ); 
   
   etaProbeVSmassDenom_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				   theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				 .M(),
				 theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.eta() ); 
   
   if ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() > theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt())
     {
       pTele1Denominator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
       pTele2Denominator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
     }
   else 
     {
       pTele1Denominator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
       pTele2Denominator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
     }


   
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   // does the probe match any of the objects from the last filter of the TARGET HLT path?  
   ///////////////////////////////////////////////////////////////////////////////////////////////////
   bool hasProbeFiredTargetPath(false);
   int lastFilterTargetMatchedObjIndex(-1);
   for(size_t objNrP=0; objNrP< trigObjP4sLastFilterTargetPath.size(); objNrP++){
     
     if (  deltaR( trigObjP4sLastFilterTargetPath[objNrP].Eta(),
		   trigObjP4sLastFilterTargetPath[objNrP].Phi(),
		   theProbe.Eta(),  
		   theProbe.Phi()  
		   )  < 0.1 ) 
       {
	 hasProbeFiredTargetPath         =true;
	 counterEvtsWithTargetFired_ ++;
	 lastFilterTargetMatchedObjIndex =objNrP;
	 break;
	 // stop after you've found the first end-of-target-path object matching the probe => the probe has fired the trigger!
       }// if probe is matched to end-of-target-path object 
   }// loop over end-of-target-path object
   

   if (doDebugMessages_) 
     { 
       std::cout << "++ conclusion:\nthe TAG is: et: " 
		 << theTag.Et() 
		 << "\t eta: " << theTag.Eta() 
		 << "\t phi: " << theTag.Phi()
		 << "\nthe PROBE is: et: " 
		 << theProbe.Et() 
		 << "\t eta: " << theProbe.Eta() 
		 << "\t phi: " << theProbe.Phi();
       if(hasProbeFiredTargetPath)
	 {
	   std::cout << "\n PROBE passes and is matched to: et: "
		     << trigObjP4sLastFilterTargetPath[lastFilterTargetMatchedObjIndex].Et() 
		     << "\t eta: " << trigObjP4sLastFilterTargetPath[lastFilterTargetMatchedObjIndex].Eta() 
		     << "\t phi: " << trigObjP4sLastFilterTargetPath[lastFilterTargetMatchedObjIndex].Phi();
	 }
       else 
	 {
	   std::cout << "\n PROBE failed. " << std::endl;
	 }
     }// end debug
   
   
   if(hasProbeFiredTargetPath)
     {
       eventsFate_->Fill(9);
       massNumerator_    -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				    theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				  .M() ); 
       pTTagNumerator_   -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
       pTProbeNumerator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
       ptProbeVSmassNum_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				    theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				  .M(),
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() ); 
       
       etaProbeVSmassNum_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				    theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				  .M(),
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.eta() ); 
       
       if ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() > theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt())
	 {
	   pTele1Numerator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
	   pTele2Numerator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
	 }
       else 
	 {
	   pTele1Numerator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
	   pTele2Numerator_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
	 }
     }// end "if probe passed"
   else
     {// "if probe failed"
       massFail_    -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
			       theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
			     .M() ); 
       pTTagFail_   -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
       pTProbeFail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
       ptProbeVSmassFail_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				    theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				  .M(),
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );        
       etaProbeVSmassFail_ -> Fill( ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.p4() + 
				    theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.p4() )
				  .M(),
				  theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.eta() );
       
       if ( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() > theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt())
	 {
	   pTele1Fail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
	   pTele2Fail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
	 }
       else 
	 {
	   pTele1Fail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingProbeIndex).first.pt() );
	   pTele2Fail_ -> Fill( theEleOfflineCands.at(theOfflineEleMatchingTagIndex).first.pt() );
	 }
     }// end "if probe failed"
   
   
   
}// end analyze



bool HeavyNuEleTriggerEff::passOfflineSelection( const edm::Event& iEvent, 
						 const edm::EventSetup& iSetup,
						 std::vector< std::pair<pat::Electron, float> > & eCands,
						 std::vector< std::pair<pat::Jet, float> > & jetCands 
						 )
{// passOfflineSelection
  
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
	  eCands.size()   >= 2  && // not so good this '2' is hard coded, but leave it for now. FIXME
	  jetCands.size() >= static_cast<unsigned int> (numOfflJets_)
	  )    ;

}// end of passOfflineSelection 



// ------------ method which returns references to objects that match a certain trigger  ------------
pat::TriggerObjectRefVector HeavyNuEleTriggerEff::findObjectsMatchedToPath( const edm::Handle< pat::TriggerEvent > triggerPatEvent , 
								     const pat::TriggerPathRef seedingAcceptedPath,
								     const std::vector<std::string>& electronFilters )
{

  // one of the seding HLT paths has fired. Chose the first one as a base trigger
  //const pat::TriggerPathRef seedingAcceptedPath = triggerPatEvent->pathRef( seedHLTForHLTEfficiencyAccepted.at(0) ) ; 
  
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
 
  pat::TriggerObjectRefVector objectsInPath = triggerPatEvent->pathObjects(seedingAcceptedPath->name(),true) ;
  pat::TriggerFilterRefVector filtersInPath = triggerPatEvent->pathFilters(seedingAcceptedPath->name(),true) ;
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
	
	if ( triggerPatEvent->objectInFilter(objRef,filterRef->label()) ) { // Trigger object was used by the filter
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



// ------------ method which returns references to objects that match a certain trigger  ------------
pat::TriggerObjectRefVector HeavyNuEleTriggerEff::findObjectsMatchedToPathNew( const edm::Handle< pat::TriggerEvent > triggerPatEvent , 
									       const pat::TriggerPathRef seedingAcceptedPath,
									       const std::vector<std::string>& electronFilters )
{

  // one of the seding HLT paths has fired. Chose the first one as a base trigger
  //const pat::TriggerPathRef seedingAcceptedPath = triggerPatEvent->pathRef( seedHLTForHLTEfficiencyAccepted.at(0) ) ; 
  
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
 
  pat::TriggerObjectRefVector objectsInPath = triggerPatEvent->pathObjects(seedingAcceptedPath->name(),true) ;
  pat::TriggerFilterRefVector filtersInPath = triggerPatEvent->pathFilters(seedingAcceptedPath->name(),true) ;
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
	
	if ( triggerPatEvent->objectInFilter(objRef,filterRef->label()) ) { // Trigger object was used by the filter
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

  std::cout << "\n\n\nHeavyNuEleTriggerEff was executed: " << counterExecutions_ << " times " << std::endl;
  std::cout << "\n\nHeavyNuEleTriggerEff and has seen / analyzed: " << counterEvtsAll_ << " events (fitting desired run range);" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsPassingOfflineSelection_ << " pass the offline selection : 2 HEEP ele pt > "
	      << minPtOfflEle_ << " and" << numOfflJets_ << " jets pt > " << minPtOfflJets_ << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithSeedingFired_ << " with >=1 seeding paths firing" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithTargetValid_  << " with >=1 valid target path" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithSeedingAbovePt_ << " with >=2 seeding objects above pt " << minPtOfObjects_ << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithSeedingPathProbesMatchedToOfflineObj_ << " all end-path-objects of seed path are matched to offline ele's " << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithTagFound_ << " tag was found  " << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithProbeFound_ << " probe was found " << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsTangProbeMatchedOffl_ << " tag and probe are matched to offline objs" << std::endl;
  std::cout << "\t\t ++ of which: " << counterEvtsWithTargetFired_  << " the probe passed " << std::endl;
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
