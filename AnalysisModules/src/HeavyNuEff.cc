// system include files
#include <iostream>
#include <algorithm>
#include <vector>

#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuMuHist.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEff.h"

void hnu::studyMuonEff(edm::Event& iEvent, std::vector<pat::Muon>& muCands,
                       HeavyNuEvent& hnuEvent,
                       std::vector< std::pair<pat::Jet, float> >& jetCands,
                       edm::Handle<pat::GenericParticleCollection>& pTracks,
                       edm::Handle<reco::TrackCollection>& gTracks,
                       HeavyNuTrigger *trig_,
                       hnu::CutsStruct& cuts,
                       EffInfo& effinfo,
                       HeavyNu::HistStruct& hists)
{
    std::vector<pat::Muon> tagMuons;
    for (unsigned int i = 0; i < muCands.size(); i++)
    {
        // Muons already tight, in proper detector region, with sufficient pT
        pat::Muon tagCand = muCands.at(i) ;
        if (hnu::muIsolation(tagCand) < cuts.muon_trackiso_limit)
        {
            if ( (trig_->matchingEnabled() && iEvent.isRealData() && trig_->isTriggerMatched(tagCand, iEvent)) ||
                (!iEvent.isRealData() && trig_->simulateForMC(tagCand.pt(), tagCand.eta())) )
            {
                bool jetOverlap = false ;
                for (unsigned int j = 0; j < jetCands.size(); j++)
                {
                    pat::Jet jet = jetCands.at(j).first ;
                    if ( hnu::jetID(jet) > 0 )
                    {
                        double dR = deltaR(tagCand.eta(), tagCand.phi(), jet.eta(), jet.phi()) ;
                        if (dR <= cuts.minimum_muon_jet_dR)
                        {
                            jetOverlap = true ;
                            break ;
                        }
                    }
                }
                // If muon is separated from jets, is isolated, and trigger matched, it is a tag
                if ( !jetOverlap ) tagMuons.push_back( tagCand ) ;
            }
        }
    }

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    if (!iEvent.getByLabel("offlineBeamSpot", beamSpotHandle))
        throw cms::Exception("Trying to do efficiency studies, cannot find beam spot");

    // Make sure that we have at least two jets
    // keeping in mind that we only study efficiency for the nominal case
    std::vector<pat::Jet> validJets ;
    for (unsigned int i = 0; i < jetCands.size(); i++)
    {
        pat::Jet jet = jetCands.at(i).first ;
        if ( hnu::jetID(jet) > 0 && jet.pt() > cuts.minimum_jet_pt ) validJets.push_back( jet ) ;
    }
    // if ( validJets.size() >= 2 && tagMuons.size() > 0 ) {
    if ( tagMuons.size() > 0 )
    { // Allow for checks of ID efficiency for different jet multiplicity
        //--- Create the two lists of probe muons ---//
        std::vector<pat::Muon> tightMuonProbes ;
        std::vector<pat::Muon> cjMuonProbes ;
        std::vector<pat::Muon> trigProbes ;
        for (unsigned int i = 0; i < muCands.size(); i++)
        {
            bool jetOverlap = false ;
            double minDeltaR = 0.8 ;
            for (unsigned int j = 0; j < validJets.size(); j++)
            {
                pat::Jet jet = validJets.at(j) ;
                double dR = deltaR(muCands.at(i).eta(), muCands.at(i).phi(), jet.eta(), jet.phi()) ;
                if (dR < minDeltaR) minDeltaR = dR ;
                if (dR <= cuts.minimum_muon_jet_dR)
                {
                    jetOverlap = true ;
                    break ;
                }
            }
            if ( !jetOverlap )
            {
                tightMuonProbes.push_back( muCands.at(i) ) ;
                if ( minDeltaR < 0.8 ) cjMuonProbes.push_back( muCands.at(i) ) ;
                if ( fabs(muCands.at(i).eta()) < 2.1 &&
                    hnu::muIsolation(muCands.at(i)) < cuts.muon_trackiso_limit )
                    trigProbes.push_back( muCands.at(i) ) ;
            }
        }
        std::vector< std::pair<pat::Muon, pat::Muon> > tagTightProbes =
                hnu::getTagProbePair<pat::Muon > ( tagMuons, tightMuonProbes, effinfo.ZwinMinGeV, effinfo.ZwinMaxGeV,
                ((effinfo.oneTP)?(effinfo.tpRandom->Uniform()):(-1.0)) ) ;
        std::vector< std::pair<pat::Muon, pat::Muon> > tagCJmuonProbes =
                hnu::getTagProbePair<pat::Muon > ( tagMuons, cjMuonProbes, effinfo.ZwinMinGeV, effinfo.ZwinMaxGeV,
                ((effinfo.oneTP)?(effinfo.tpRandom->Uniform()):(-1.0)), false ) ;
        std::vector< std::pair<pat::Muon, pat::Muon> > tagTrigProbes =
                hnu::getTagProbePair<pat::Muon > ( tagMuons, trigProbes, effinfo.ZwinMinGeV, effinfo.ZwinMaxGeV,
                ((effinfo.oneTP)?(effinfo.tpRandom->Uniform()):(-1.0)) ) ;

        std::vector<pat::GenericParticle> trackProbes ;
        pat::GenericParticleCollection trackCands = *(pTracks.product());
        std::sort(trackCands.begin(), trackCands.end(), hnu::pTcompare());
        for (unsigned int i = 0; i < trackCands.size(); i++)
        {
            pat::GenericParticle trackCand = trackCands.at(i) ;
            if ( trackCand.pt() <= cuts.minimum_mu2_pt ) break ; // Sorted collection: quit once below
            if ( fabs(trackCand.eta()) >= cuts.maximum_mu_abseta ) continue ;

            bool jetOverlap  = false ;
            for (unsigned int j = 0; j < validJets.size(); j++)
            {
                pat::Jet jet = validJets.at(j) ;
                double dR = deltaR(trackCand.eta(), trackCand.phi(), jet.eta(), jet.phi()) ;
                if (dR < cuts.minimum_muon_jet_dR)
                {
                    jetOverlap = true ;
                    break ;
                }
            }
            if ( !jetOverlap ) trackProbes.push_back( trackCand ) ;
        }
        std::vector< std::pair<pat::Muon, pat::GenericParticle> > tagTrackProbes =
                hnu::getTagProbePair<pat::GenericParticle > ( tagMuons, trackProbes, effinfo.ZwinMinGeV, effinfo.ZwinMaxGeV,
                ((effinfo.oneTP)?(effinfo.tpRandom->Uniform()):(-1.0)) ) ;

        //if ( tagTightProbes.size() > 0 || tagTrackProbes.size() > 0 ) debuggingEvents = true ;

        // Defined both tag+probe collections
        studyMuonSelectionEff(tagTrackProbes, tagTightProbes, tagCJmuonProbes,
                              muCands, gTracks, beamSpotHandle, hnuEvent.eventWgt, validJets.size(), cuts, hists);

        // Study Trigger Matching efficiency for data only...here 2 or more jets is enforced
        if ( trig_->matchingEnabled() && iEvent.isRealData() && validJets.size() >= 2 )
        {
            for (unsigned int i = 0; i < tagTrigProbes.size(); i++)
            {
                pat::Muon theTag   = tagTrigProbes.at(i).first ;
                pat::Muon theProbe = tagTrigProbes.at(i).second ;
                hists.TightTagTrigProbeInZwin->tapfill( theTag, theProbe, theProbe.trackIso(), hnuEvent.eventWgt ) ;
                if ( trig_->isTriggerMatched(theProbe, iEvent) )
                    hists.TightTagTrigProbePassesInZwin->tapfill( theTag, theProbe, theProbe.trackIso(), hnuEvent.eventWgt ) ;
            }
        }
    }
}


void hnu::studyMuonSelectionEff(const std::vector< std::pair<pat::Muon, pat::GenericParticle> >& trackTP,
                                const std::vector< std::pair<pat::Muon, pat::Muon> >& tightTP,
                                const std::vector< std::pair<pat::Muon, pat::Muon> >& cjTP,
                                const std::vector< pat::Muon >& tightMuons,
                                const edm::Handle<reco::TrackCollection>& gTracks,
                                const edm::Handle<reco::BeamSpot>& beamspot,
                                double wgt, int nJets, hnu::CutsStruct& cuts,
                                HeavyNu::HistStruct& hists)
{
    
    reco::TrackCollection generalTracks = *(gTracks.product()); 

    // First case: Check muon ID
    // Probe is generic track meeting pT, eta requirements
    for (unsigned int i=0; i<trackTP.size(); i++) {
      pat::Muon            theTag   = trackTP.at(i).first ; 
      pat::GenericParticle theProbe = trackTP.at(i).second ; 
      // Calculate the isolation for the probe
      double trkSumPtIsoCone = 0. ; 
      for (unsigned int k=0; k<generalTracks.size(); k++) {
	reco::Track track = generalTracks.at(k) ; 
	double dR = deltaR(track.eta(), track.phi(), theProbe.eta(), theProbe.phi());
	if ( fabs(dR) > 0.3 || fabs(dR) < 0.01 ) continue ;
	if ( fabs(track.dxy(beamspot->position())) > 0.1 ) continue ;
	double dz = fabs( theProbe.vertex().Z() - track.vertex().Z() ) ;
	if ( dz > 0.2 ) continue ;
	trkSumPtIsoCone += track.pt() ;
      }
      if ( nJets >= 0 ) hists.TightTagTrackProbeInZwin0jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTrackProbeInZwin1jet->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTrackProbeInZwin2jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      // Does the probe pass tight ID selection?
      unsigned int muIdx = tightMuons.size() ;
      for (unsigned int k=0; k<tightMuons.size(); k++) {
	pat::Muon tightMuon = tightMuons.at(k) ; 
	if ( deltaR(tightMuon.eta(), tightMuon.phi(), theProbe.eta(), theProbe.phi()) < 0.02 &&
	     fabs(theProbe.pt()-tightMuon.pt())/theProbe.pt() < 0.05 ) {
	  muIdx = k ; break ; 
	}
      }
      if ( muIdx < tightMuons.size() ) {
          if ( nJets >= 0 ) hists.TightTagTrackProbePassesInZwin0jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTrackProbePassesInZwin1jet->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTrackProbePassesInZwin2jets->tapfill( theTag,theProbe,trkSumPtIsoCone,wgt ) ;
      }
    }
        
    // Second case: Check isolation
    // Probe is tight muon separated from jets
    for (unsigned int i=0; i<tightTP.size(); i++) {
      pat::Muon theTag   = tightTP.at(i).first ; 
      pat::Muon theProbe = tightTP.at(i).second ; 
      // Confirmed: we have a valid probe
      if ( nJets >= 0 ) hists.TightTagTightProbeInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTightProbeInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTightProbeInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( hnu::muIsolation(theProbe,1.0) < cuts.muon_trackiso_limit ) { 
          if ( nJets >= 0 ) hists.TightTagTightProbePassesInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTightProbePassesInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTightProbePassesInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      }
    }

    // Third case:
    // Probe is tight muon with 0.5 < dR(mu,j) < 0.8
    for (unsigned int i=0; i<cjTP.size(); i++) {
      pat::Muon theTag   = cjTP.at(i).first ; 
      pat::Muon theProbe = cjTP.at(i).second ; 
      if ( nJets >= 0 ) hists.TightTagTightCJProbeInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 1 ) hists.TightTagTightCJProbeInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( nJets >= 2 ) hists.TightTagTightCJProbeInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      if ( hnu::muIsolation(theProbe,1.0) < cuts.muon_trackiso_limit ) { 
          if ( nJets >= 0 ) hists.TightTagTightCJProbePassesInZwin0jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 1 ) hists.TightTagTightCJProbePassesInZwin1jet->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
          if ( nJets >= 2 ) hists.TightTagTightCJProbePassesInZwin2jets->tapfill( theTag,theProbe,theProbe.trackIso(),wgt ) ;
      }
    }            
}

void hnu::studyElectronEff(edm::Event& iEvent,
                           std::vector< std::pair<pat::Electron, float> >& eCands,
                           edm::Handle<pat::ElectronCollection>& pElecs,
                           HeavyNuEvent& hnuEvent,
                           std::vector< std::pair<pat::Jet, float> >& jetCands,
                           HeavyNuTrigger *trig_,
                           hnu::CutsStruct& cuts,
                           EffInfo& effinfo,
                           HeavyNu::HistStruct& hists)
{
    std::vector<pat::Electron> tagElectrons ;
    for (unsigned int i = 0; i < eCands.size(); i++)
    {
        // Electrons already pass HEEP selection, in proper detector region, with sufficient pT
        pat::Electron tagCand = eCands.at(i).first ;
        if ( (trig_->matchingEnabled() && iEvent.isRealData() && trig_->isTriggerMatched(tagCand, tagCand, iEvent, 1)) ||
            !iEvent.isRealData() )
        { // 1) for data we can check a trigger match...for MC it will happen later
            // 2) in the muon case tags are also selected according to the HLT efficiency;
            //    deemed an overkill + now we're doing HLT efficiencies as a function of mass => ignore this for electrons
            bool jetOverlap = false ;
            int nValidJets = 0 ;
            for (unsigned int j = 0; j < jetCands.size(); j++)
            {
                pat::Jet jet = jetCands.at(j).first ;
                if ( hnu::jetID(jet) > 0 )
                {
                    nValidJets++ ;
                    double dR = deltaR(tagCand.eta(), tagCand.phi(), jet.eta(), jet.phi()) ;
                    if (dR <= cuts.minimum_muon_jet_dR)
                    {
                        jetOverlap = true ;
                        break ;
                    }
                }
            }
            // If electron is separated from jets, is isolated, it is a tag candidate
            // NOTE: Both tag and probe must match to trigger, so full TP pair not settled!
            if ( !jetOverlap ) tagElectrons.push_back( tagCand ) ;
        }
    } // List of all possible electron tag candidates

    // std::cout << "I have a tag list of size " << tagElectrons.size() << std::endl ;

    // Make sure that we have at least two jets
    // keeping in mind that we only study efficiency for the nominal case
    std::vector<pat::Jet> validJets ;
    for (unsigned int i = 0; i < jetCands.size(); i++)
    {
        pat::Jet jet = jetCands.at(i).first ;
        if ( hnu::jetID(jet) > 0 && jet.pt() > cuts.minimum_jet_pt ) validJets.push_back( jet ) ;
    }

    if ( tagElectrons.size() > 0 )
    { // Allow for checks of ID efficiency for different jet multiplicity
        //--- Create the list of probe electrons ---//
        std::vector<pat::Electron> gsfElectronProbes ;
        std::vector<pat::Electron> gsfCands = *(pElecs.product()) ;
        std::sort(gsfCands.begin(), gsfCands.end(), hnu::pTcompare()) ;

        for (unsigned int i = 0; i < gsfCands.size(); i++)
        {
            pat::Electron gsfCand = gsfCands.at(i) ;
            if ( hnu::getElectronEt(gsfCand, false) < cuts.minimum_mu2_pt ) break ; // Sorted collection, quit once below
            bool jetProbeOverlap = false ;
            for (unsigned int j = 0; j < jetCands.size(); j++)
            {
                pat::Jet jet = jetCands.at(j).first ;
                if ( hnu::jetID(jet) > 0 )
                {
                    double dR = deltaR(gsfCand.eta(), gsfCand.phi(), jet.eta(), jet.phi()) ;
                    if (dR <= cuts.minimum_muon_jet_dR)
                    {
                        jetProbeOverlap = true ;
                        break ;
                    }
                }
            }
            if ( jetProbeOverlap )
            {
                continue ;
            }

            bool isEB = ( fabs(hnu::getElectronSCEta(gsfCand)) < 1.4442 ) ;
            bool isEE = ( fabs(hnu::getElectronSCEta(gsfCand)) < 2.5 && fabs(hnu::getElectronSCEta(gsfCand)) > 1.56 ) ;
            if ( !isEB && !isEE ) continue ;
            // Probe must participate in the trigger!
            if ( (trig_->matchingEnabled() && iEvent.isRealData() && trig_->isTriggerMatched(gsfCand, gsfCand, iEvent, 1)) ||
                !iEvent.isRealData() )
            {
                gsfElectronProbes.push_back( gsfCand ) ;
            }
        }

        //
        // Now need to introduce a trigger requirement for MC...something to be done later
        //

        std::vector< std::pair<pat::Electron, pat::Electron> > heepTagGsfProbes =
                hnu::getTagProbePair( tagElectrons, gsfElectronProbes, effinfo.ZwinMinGeV, effinfo.ZwinMaxGeV,
                                     ((effinfo.oneTP)?(effinfo.tpRandom->Uniform()):(-1.0)) ) ;

        // std::cout << "I have " << heepTagGsfProbes.size() << " tag/probe pairs" << std::endl ;

        studyElectronSelectionEff(heepTagGsfProbes, eCands, hnuEvent.eventWgt, validJets.size(), hists);
    }
}

void hnu::studyElectronSelectionEff(const std::vector< std::pair<pat::Electron, pat::Electron> >& TP,
                                    const std::vector< std::pair<pat::Electron, float> >& heepElectrons,
                                    double wgt, int nJets, HeavyNu::HistStruct& hists)
{
    
  // Only case: Check HEEP ID
  // Probe is GSF track that passes pT, eta, trigger requirements
  for (unsigned int i=0; i<TP.size(); i++) {
    pat::Electron theTag   = TP.at(i).first ; 
    pat::Electron theProbe = TP.at(i).second ; 

    if ( nJets >= 0 ) hists.HeepTagGsfProbeInZwin0jets->tapfill( theTag,theProbe,wgt ) ;
    if ( nJets >= 1 ) hists.HeepTagGsfProbeInZwin1jet->tapfill( theTag,theProbe,wgt ) ;
    if ( nJets >= 2 ) hists.HeepTagGsfProbeInZwin2jets->tapfill( theTag,theProbe,wgt ) ;
    
    // std::cout << "Tag electron pT " << hnu::getElectronEt(theTag,false) 
    // 		<< " and probe electron pT " << hnu::getElectronEt(theProbe,false) << std::endl ; 

    // Does the probe pass tight ID selection?
    unsigned int eIdx = heepElectrons.size() ;
    for (unsigned int k=0; k<heepElectrons.size(); k++) {
      pat::Electron heepElectron = heepElectrons.at(k).first ; 
      // std::cout << "Comparing a probe with pT " << hnu::getElectronEt(theProbe,false) 
      // 		<< " and eta " << hnu::getElectronSCEta(theProbe) << " to heep electron with pT " 
      // 		<< hnu::getElectronEt(heepElectron,false) 
      // 		<< " and eta " << hnu::getElectronSCEta(heepElectron) 
      // 		<< std::endl ; 
      // std::cout << "dR: " << deltaR(heepElectron.eta(), heepElectron.phi(), theProbe.eta(), theProbe.phi())
      // 		<< " and dPt: " << fabs(hnu::getElectronEt(theProbe,false)-hnu::getElectronEt(heepElectron,false))/hnu::getElectronEt(theProbe,false)
      // 		<< std::endl ; 
      
      if ( deltaR(heepElectron.eta(), heepElectron.phi(), theProbe.eta(), theProbe.phi()) < 0.02 &&
	   fabs(hnu::getElectronEt(theProbe,false)-hnu::getElectronEt(heepElectron,false))/hnu::getElectronEt(theProbe,false) < 0.05 ) {
	eIdx = k ; break ; 
      }
    }

    // if ( eIdx < heepElectrons.size() ) 
    //   std::cout << "TNP SUCCESS: Found a match in position " << eIdx << std::endl ; 
    // else 
    //   std::cout << "TNP FAILURE: Found no match in the good electrons list!!! " << iEvent.id() << std::endl ; 
    
    if ( eIdx < heepElectrons.size() ) {
      if ( nJets >= 0 ) hists.HeepTagGsfProbePassesInZwin0jets->tapfill( theTag,theProbe,wgt ) ;
      if ( nJets >= 1 ) hists.HeepTagGsfProbePassesInZwin1jet->tapfill( theTag,theProbe,wgt ) ;
      if ( nJets >= 2 ) hists.HeepTagGsfProbePassesInZwin2jets->tapfill( theTag,theProbe,wgt ) ;
    }
  }        
}
