#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"

// #include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
// #include "DataFormats/PatCandidates/interface/TriggerPath.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuTrigger::HeavyNuTrigger(const edm::ParameterSet & iConfig) :
  trigEventTag_( iConfig.getParameter< edm::InputTag >( "trigEventTag" ) ),
  muonTriggers_( iConfig.getParameter< std::vector<std::string> >( "muonTriggers" ) ),
  beginRun_    ( iConfig.getParameter< std::vector<int> >( "firstRun" ) ),
  endRun_      ( iConfig.getParameter< std::vector<int> >( "lastRun" ) ),
  muonMatch_   ( iConfig.getParameter< std::string >  ( "muonMatch"    ) ),
  triggerPt_   ( iConfig.getParameter< double >       ( "triggerPt"    ) ),
  trigEra_     ( iConfig.getParameter< int >          ( "trigEra" ) ),
  johnnyApple_ ( iConfig.getParameter< int >          ( "randomSeed"   ) )
{
  matchingEnabled_ = false;
  if (trigEventTag_.label().size() &&
      muonMatch_.size())
    matchingEnabled_ = true;

  if (!matchingEnabled_) {
    std::cout << "Trigger matching is === DISABLED ===" << std::endl;
    std::cout << "   (Trigger sim. random seed = "<<johnnyApple_<<")"<<std::endl;
    triggerRandom_ = new TRandom(johnnyApple_);
    std::cout << "   (Trigger Era              = "<<trigEra_<<")"<<std::endl;
    std::cout << "   (Trigger pT               = "<<triggerPt_<<")"<<std::endl;
  }

  // Safety: first/last run vectors are not the same size as the trigger list
  if ( muonTriggers_.size() != beginRun_.size() || 
       muonTriggers_.size() != endRun_.size() || 
       beginRun_.size() != endRun_.size() ) { 
    std::cout << "INFO: Size of trigger initialization vectors are not equal." << std::endl ; 
    std::cout << "      Resetting all vectors to size of trigger list and dropping run restriction " << std::endl ; 

    unsigned int maxsize = std::max( muonTriggers_.size(),
				     std::max( beginRun_.size(),endRun_.size() ) ) ; 
    
    beginRun_.clear() ; 
    endRun_.clear() ; 
    for (unsigned int i=0; i<maxsize; i++) { 
      beginRun_.push_back(0) ; 
      endRun_.push_back(999999) ; 
    }
  }
}                                      // HeavyNuTrigger::HeavyNuTrigger

//======================================================================

void
HeavyNuTrigger::book(const TFileDirectory& tdir, trigHistos_t *thist)
{
  if (!matchingEnabled_) return;

  assert(thist);

  thist->trigMatchPtCorrel = tdir.make< TH2D >( "h2d_trigMatchPtCorrel","", 60,0.,300.,60,0.,300. );
  thist->trigMatchPtCorrel->SetTitle ( "#mu_{reco} vs. #mu_{trig} p_{T} (GeV)" );
  thist->trigMatchPtCorrel->SetXTitle( "p_{T}(#mu_{reco}) (GeV)" );
  thist->trigMatchPtCorrel->SetYTitle( "p_{T}(#mu_{trig}) (GeV)" );

  thist->trigMatchDR2      = tdir.make< TH1D >( "h1d_trigMatchDR2", "",50,0.,0.25 );
  thist->trigMatchDR2      ->SetTitle ( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );
  thist->trigMatchDR2      ->SetXTitle( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );

  thist->trigMatchDRDPt    = tdir.make< TH2D >( "h2d_trigMatchDRDPt","",50,0.,0.25,50,-0.5,0.5 );
  thist->trigMatchDRDPt    ->SetTitle ( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2} vs #Delta p_{T}_{rel}(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDRDPt    ->SetXTitle( "(#Delta R(#mu_{reco}, #mu_{trig}))^{2}" );
  thist->trigMatchDRDPt    ->SetYTitle( "#Delta p_{T}_{rel}(#mu_{reco}, #mu_{trig})" );

  thist->trigMatchDetaPhi  = tdir.make< TH2D >( "h2d_trigMatchDetaPhi","",50,-0.25,0.25,50,-0.25,0.25 );
  thist->trigMatchDetaPhi  ->SetTitle ( "#Delta #phi(#mu_{reco}, #mu_{trig}) vs #Delta #eta(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDetaPhi  ->SetXTitle( "#Delta #eta(#mu_{reco}, #mu_{trig})" );
  thist->trigMatchDetaPhi  ->SetYTitle( "#Delta #phi(#mu_{reco}, #mu_{trig})" );

  thist->trigUnmatchedPt   = tdir.make< TH1D >( "h1d_trigUnmatchedPt","", 200,0,2000 );
  thist->trigUnmatchedPt   ->SetTitle( "Unmatched muon p_{T}; #mu p_{T} (GeV)" );

  thist->trigAllCandMuPt   = tdir.make< TH1D >( "h1d_trigAllCandMuPt","", 200,0,2000 );
  thist->trigAllCandMuPt   ->SetTitle( "All trigger-match muon candidates p_{T}; #mu p_{T} (GeV)" );

  thist->trigUnmatchedEtaPhi = tdir.make< TH2D >( "h2d_trigUnmatchEtaPhi","", 50,-2.5,2.5,63,0.,6.3 );
  thist->trigUnmatchedEtaPhi ->SetTitle( "Unmatched muon #eta/#phi; #mu #eta; #mu #phi" );

  thist->trigAllCandMuEtaPhi = tdir.make< TH2D >( "h2d_trigAllCandMuEtaPhi","", 50,-2.5,2.5,63,0.,6.3 );
  thist->trigAllCandMuEtaPhi ->SetTitle( "All trigger-match muon candidates; #mu #eta; #mu #phi" );

}                                                // HeavyNuTrigger::book

//======================================================================

bool
HeavyNuTrigger::isTriggerMatched(const pat::Muon&  m,
				 const edm::Event& iEvent,
				 trigHistos_t *thist)
{
  bool matched=false;
  // bool passTrig=false;
  
  // Only one trigger can be used for matching in a given run
  int run = iEvent.run() ; 
  std::vector<std::string> validHLTpaths ; 
  
  for (unsigned int i=0; i<muonTriggers_.size(); i++) 
    if (run >= beginRun_.at(i) && run <= endRun_.at(i)) validHLTpaths.push_back(muonTriggers_.at(i)) ; 

//   std::cout << "Looking for trigger match for muon with pT " << m.pt() 
//    	    << " and eta " << m.eta() << std::endl ; 

  // muon trigger matching is only allowed within |eta|<2.1
  if ( matchingEnabled_ &&(fabs(m.eta()) < 2.1) ) {
    // PAT trigger information
//     edm::Handle< pat::TriggerEvent > triggerEvent;

//     iEvent.getByLabel( trigEventTag_, triggerEvent );
//     if ( !triggerEvent.isValid() ) {
//       std::cerr << "triggerEvent not found " << std::endl;
//       return false;
//     }

    const pat::TriggerObjectStandAloneCollection muonMatchCollection = m.triggerObjectMatches();
    // std::cout << "Trigger object matches size: " << muonMatchCollection.size() << std::endl ; 

    for (unsigned int i=0; i<muonMatchCollection.size(); i++) { 
      if ( matched ) break ; // Quit as soon as we find a match
      // std::cout << "Trigger object " << i+1 << " of " << muonMatchCollection.size() << std::endl ; 
      pat::TriggerObject muonTrigger = muonMatchCollection.at(i) ; 
      pat::TriggerObjectStandAlone muonTriggerInfo = muonMatchCollection.at(i) ; 
      // Look for a match with one of our paths
      std::vector<std::string> hltPaths = muonTriggerInfo.pathNames(true,false) ; 
      bool hltPathMatch = false ; 
      for (unsigned int j=0; j<hltPaths.size(); j++) { 
	if (hltPathMatch) break ; 
	for (unsigned int k=0; k<validHLTpaths.size(); k++) { 
	  if (hltPaths.at(j) == validHLTpaths.at(k)) { 
            // std::cout << "Found a match to HLT path: " << muonTriggers_.at(k) << std::endl ; 
	    hltPathMatch = true ; 
	    break ; 
	  }
	}
      }
      // Finding a trigger object is not enough.  Need to impose the last filter (pT) 
      // Requirements to see if the trigger would have accepted the event based on this muon
      if ( hltPathMatch && muonTrigger.pt() > triggerPt_ ) { 
	double dr2  = reco::deltaR2 <pat::Muon,pat::TriggerObject>( m,muonTrigger );
	double dpt  = 1.-(muonTrigger.pt()/m.pt());
	double dphi = reco::deltaPhi( m.phi(),muonTrigger.phi() );
	double deta = m.eta() - muonTrigger.eta();

	// One more requirement: make sure that the muon is nearby the trigger object
	if ( sqrt(dr2) < 0.1 ) matched = true ; 

//         std::cout << "pT is " << muonTrigger.pt() << " with dpt = " << dpt << std::endl ; 
//         std::cout << "dR is " << sqrt(dr2) << std::endl ; 

	if (thist) {	  
	  thist->trigMatchPtCorrel->Fill( m.pt(),muonTrigger.pt() );
	  thist->trigMatchDR2     ->Fill( dr2 );
	  thist->trigMatchDRDPt   ->Fill( dr2,dpt );
	  thist->trigMatchDetaPhi ->Fill( deta,dphi );
	}
      }
      if ( !matched ) {
	if ( thist ) {
	  thist->trigUnmatchedPt->Fill( m.pt() );
	  thist->trigUnmatchedEtaPhi->Fill( m.eta(),m.phi() );
	}
      }
    }

    if ( thist ) { 
      thist->trigAllCandMuPt->Fill( m.pt() );
      thist->trigAllCandMuEtaPhi->Fill( m.eta(),m.phi() );
    }
  }
  return ( matched ); 
}                                    // HeavyNuTrigger::isTriggerMatched

//======================================================================

bool
HeavyNuTrigger::simulateForMC(double pt,double eta,int signOfError2apply)
{
  if (matchingEnabled_)
    throw cms::Exception("invalid trigger configuration");

  // Triggers outside |eta| < 2.1 not allowed
  if ( fabs(eta) >= 2.1 ) return false ; 
  // Cannot trigger if you do not meet the minimum pT threshold
  if ( pt < triggerPt_ ) return false ;

  // Trigger studies updated 11 Dec 2011
  // HLT_Mu24 results are used for 30 < pT < 40 GeV
  // All other bins combine Mu24, Mu40, Mu40_eta2p1 efficiencies for run 2011a
  // Run 2011b results use Mu40_eta2p1 only
  const double effslo2011a[]  = {0.868979,0.855499,0.862255,0.874853,0.882812,0.876860,0.823287,0.823287};
  const double effsnom2011a[] = {0.875648,0.864359,0.872848,0.885346,0.897508,0.890611,0.871189,0.871189};
  const double effshi2011a[]  = {0.882316,0.873220,0.883441,0.895840,0.912204,0.904362,0.919090,0.919090};
  const double upedge2011a[]  = {      40,      50,      60,      80,     100,     200,    3500,      -1};

  const double effslo2011b[]  = {0.908256,0.922264,0.944146,0.920008,0.906448,0.917130,0.917130};
  const double effsnom2011b[] = {0.900857,0.913545,0.936445,0.906582,0.893309,0.866667,0.866667};
  const double effshi2011b[]  = {0.893458,0.904825,0.928743,0.893155,0.880171,0.816203,0.816203};
  const double upedge2011b[]  = {      50,      60,      80,     100,     200,    3500,      -1};

  // 2011 A is the default
  const double *effs = ( (trigEra_ == 20111) ? effsnom2011a : effsnom2011b );
  if ( signOfError2apply ) {
    if ( trigEra_ == 20111 ) effs = (signOfError2apply > 0) ? effshi2011a : effslo2011a;
    else                     effs = (signOfError2apply > 0) ? effshi2011b : effslo2011b;
  }

  int i;
  const double *upedge = ( (trigEra_ == 20110) ? upedge2011a : upedge2011b );
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double eff=effs[i];
    
  return (triggerRandom_->Uniform()<eff);
}

//======================================================================

void HeavyNuTrigger::endJob()
{
}
