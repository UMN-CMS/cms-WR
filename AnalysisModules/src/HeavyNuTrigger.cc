#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuTrigger::HeavyNuTrigger(const edm::ParameterSet & iConfig) :
  trigEventTag_( iConfig.getParameter< edm::InputTag >( "trigEventTag" ) ),
  muonMatch_   ( iConfig.getParameter< std::string >  ( "muonMatch"    ) )
{
  matchingEnabled_ = false;
  if (trigEventTag_.label().size() &&
      muonMatch_.size())
    matchingEnabled_ = true;

  if (!matchingEnabled_)
    std::cout << "Trigger matching is === DISABLED ===" << std::endl;

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

}                                                // HeavyNuTrigger::book

//======================================================================

bool
HeavyNuTrigger::isTriggerMatched(const pat::MuonRef& m,
				 const edm::Event& iEvent,
				 trigHistos_t *thist)
{
  if (!matchingEnabled_) return false;

  // muon trigger matching is only allowed within |eta|<2.1
  //
  if (fabs(m->eta()) >= 2.1)
    return false;

  // PAT trigger information
  edm::Handle< pat::TriggerEvent > triggerEvent;

  iEvent.getByLabel( trigEventTag_, triggerEvent );
  if (!triggerEvent.isValid()) {
    std::cerr << "triggerEvent not found " << std::endl;
    return false;
  }

  // PAT trigger helper for trigger matching information
  const pat::helper::TriggerMatchHelper matchHelper;

  const pat::TriggerObjectRef
    trigRef( matchHelper.triggerMatchObject(reco::CandidateBaseRef(m),
					    muonMatch_,iEvent,*triggerEvent) );

  // fill histograms
  if ( trigRef.isAvailable() ) {
    std::cout << "We got one!!!" << std::endl;
    if (thist) {
      double dr2  = reco::deltaR2 <pat::Muon,pat::TriggerObject>(*m,*trigRef);
      double dphi = reco::deltaPhi(m->phi(),trigRef->phi());
      double deta = m->eta() - trigRef->eta();
      double dpt  = 1.-(trigRef->pt()/m->pt());
      
      thist->trigMatchPtCorrel->Fill(m->pt(),trigRef->pt());
      thist->trigMatchDR2     ->Fill(dr2);
      thist->trigMatchDRDPt   ->Fill(dr2,dpt);
      thist->trigMatchDetaPhi ->Fill(deta,dphi);
    }
    return true;
  }

  return false;
}                                    // HeavyNuTrigger::isTriggerMatched

//======================================================================

void HeavyNuTrigger::endJob()
{
}
