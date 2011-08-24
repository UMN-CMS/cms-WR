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
  muonMatch_   ( iConfig.getParameter< std::string >  ( "muonMatch"    ) ),
  triggerPt_   ( iConfig.getParameter< double >       ( "triggerPt"    ) ),
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
    std::cout << "   (Trigger pT               = "<<triggerPt_<<")"<<std::endl;
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
  
  // For some muons, we ignore trigger match if pT is too low
  int run = iEvent.run() ; 

  std::cout << "Looking for trigger match for muon with pT " << m.pt() 
	    << " and eta " << m.eta() << std::endl ; 

  // muon trigger matching is only allowed within |eta|<2.1
  if ( matchingEnabled_ &&(fabs(m.eta()) < 2.1) ) {
    // PAT trigger information
    edm::Handle< pat::TriggerEvent > triggerEvent;

    iEvent.getByLabel( trigEventTag_, triggerEvent );
    if ( !triggerEvent.isValid() ) {
      std::cerr << "triggerEvent not found " << std::endl;
      return false;
    }

    const pat::TriggerObjectStandAloneCollection muonMatchCollection = m.triggerObjectMatches();
    std::cout << "Trigger object matches size: " << muonMatchCollection.size() << std::endl ; 

    for (unsigned int i=0; i<muonMatchCollection.size(); i++) { 
      if ( matched ) break ; // Quit as soon as we find a match
      std::cout << "Trigger object " << i+1 << " of " << muonMatchCollection.size() << std::endl ; 
      pat::TriggerObject muonTrigger = muonMatchCollection.at(i) ; 
      if ( muonTrigger.pt() > triggerPt_ ) { 
	matched = true ; 

	double dr2  = reco::deltaR2 <pat::Muon,pat::TriggerObject>( m,muonTrigger );
	double dpt  = 1.-(muonTrigger.pt()/m.pt());
	double dphi = reco::deltaPhi( m.phi(),muonTrigger.phi() );
	double deta = m.eta() - muonTrigger.eta();
	std::cout << "pT is " << muonTrigger.pt() << " with dpt = " << dpt << std::endl ; 
	std::cout << "dR is " << sqrt(dr2) << std::endl ; 

    // PAT trigger helper for trigger matching information
    // const pat::helper::TriggerMatchHelper matchHelper;

    // const pat::TriggerObjectRef
    //   trigRef( matchHelper.triggerMatchObject( reco::CandidateBaseRef(m.originalObject()),
    //    					       muonMatch_,iEvent,*triggerEvent) );
    
    // Finding a trigger object is not enough.  Need to impose the last filter (pT) 
    // Requirements to see if the trigger would have accepted the event based on this muon

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
  
  // determined from HLT_Mu24 studies...updated 23 Aug 2011
  const double effslo24[]  = {0.945,0.939,0.946,0.948,0.935,0.954,0.954};
  const double effsnom24[] = {0.931,0.923,0.928,0.930,0.899,0.924,0.924};
  const double effshi24[]  = {0.957,0.952,0.961,0.963,0.961,0.974,0.974};
  const double upedge24[]  = {   40,   50,   60,   80,  100,  3500,   -1};

  // determined from HLT_Mu40 studies...updated 23 Aug 2011
  const double effslo40[]  = {0.928,0.934,0.942,0.946,0.941,0.941};
  const double effsnom40[] = {0.922,0.925,0.934,0.934,0.930,0.930};
  const double effshi40[]  = {0.934,0.942,0.949,0.956,0.951,0.951};
  const double upedge40[]  = {50,   60,   80,   100,  3500, -1};

  const double *effs = ( (triggerPt_ == 24) ? effsnom24 : effsnom40 );
  if ( signOfError2apply ) {
    if ( triggerPt_ == 24 ) effs = (signOfError2apply > 0) ? effshi24 : effslo24;
    else                    effs = (signOfError2apply > 0) ? effshi40 : effslo40;
  }

  int i;
  const double *upedge = ( (triggerPt_ == 24) ? upedge24 : upedge40 );
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double eff=effs[i];
    
  return (triggerRandom_->Uniform()<eff);
}

//======================================================================

void HeavyNuTrigger::endJob()
{
}
