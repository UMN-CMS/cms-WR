#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "TMath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//using namespace pat;
//using namespace pat::helper;
//using namespace TMath;

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
}

void HeavyNuTrigger::book(const TFileDirectory& tdir)
{
  if (!matchingEnabled_) return;

  trigMatchPtCorrel = tdir.make< TH2D >( "h2d_trigMatchPtCorrel","", 60,0.,300.,60,0.,300. );
  trigMatchPtCorrel->SetTitle ( "#mu_{1} vs. #mu_{trig} p_{T} (GeV)" );
  trigMatchPtCorrel->SetXTitle( "p_{T}(#mu_{1}) (GeV)" );
  trigMatchPtCorrel->SetYTitle( "p_{T}(#mu_{trig}) (GeV)" );

  trigMatchDR2      = tdir.make< TH1D >( "h1d_trigMatchDR2", "",50,0.,0.25 );
  trigMatchDR2      ->SetTitle ( "(#Delta R(#mu_{1}, #mu_{trig}))^{2}" );
  trigMatchDR2      ->SetXTitle( "(#Delta R(#mu_{1}, #mu_{trig}))^{2}" );

  trigMatchDRDPt    = tdir.make< TH2D >( "h2d_trigMatchDRDPt","",50,0.,0.25,50,0.,0.5 );
  trigMatchDRDPt    ->SetTitle ( "(#Delta R(#mu_{1}, #mu_{trig}))^{2} vs #Delta p_{T}_{rel}(#mu_{1}, #mu_{trig})" );
  trigMatchDRDPt    ->SetXTitle( "(#Delta R(#mu_{1}, #mu_{trig}))^{2}" );
  trigMatchDRDPt    ->SetYTitle( "#Delta p_{T}_{rel}(#mu_{1}, #mu_{trig})" );

  trigMatchDetaPhi  = tdir.make< TH2D >( "h2d_trigMatchDetaPhi","",50,0.,0.25,50,0.,0.25 );
  trigMatchDetaPhi  ->SetTitle ( "#Delta #phi(#mu_{1}, #mu_{trig}) vs #Delta #eta(#mu_{1}, #mu_{trig})" );
  trigMatchDetaPhi  ->SetXTitle( "#Delta #eta(#mu_{1}, #mu_{trig})" );
  trigMatchDetaPhi  ->SetYTitle( "#Delta #phi(#mu_{1}, #mu_{trig})" );
}

bool HeavyNuTrigger::isTriggerMatched(const pat::MuonRef& m,
				      const edm::Event& iEvent)
{
  if (!matchingEnabled_) return false;

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
    double dr2  = reco::deltaR2 <pat::Muon,pat::TriggerObject>(*m,*trigRef);
    double dphi = reco::deltaPhi(m->phi(),trigRef->phi());
    double deta = m->eta() - trigRef->eta();
    double dpt  = 1.-(trigRef->pt()/m->pt());

    trigMatchPtCorrel->Fill(m->pt(),trigRef->pt());
    trigMatchDR2     ->Fill(dr2);
    trigMatchDRDPt   ->Fill(dr2,dpt);
    trigMatchDetaPhi ->Fill(deta,dphi);
    return true;
  }

  return false;
}

void HeavyNuTrigger::endJob()
{
}
