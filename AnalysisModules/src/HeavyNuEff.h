#ifndef HeavyNuEff_h
#define HeavyNuEff_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "HeavyNu/AnalysisModules/src/HeavyNu.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

namespace hnu
{
   void studyMuonEff(edm::Event& iEvent, std::vector<pat::Muon>& muCands,
                     HeavyNuEvent& hnuEvent,
                     std::vector< std::pair<pat::Jet, float> >& jetCands,
                     edm::Handle<pat::GenericParticleCollection>& pTracks,
                     edm::Handle<reco::TrackCollection>& gTracks,
                     HeavyNuTrigger *trig_,
                     hnu::CutsStruct& cuts,
                     EffInfo& effinfo,
                     HeavyNu::HistStruct& hists);

   void studyMuonSelectionEff(const std::vector< std::pair<pat::Muon, pat::GenericParticle> >& trackTP,
                              const std::vector< std::pair<pat::Muon, pat::Muon> >& tightTP,
                              const std::vector< std::pair<pat::Muon, pat::Muon> >& cjTP,
                              const std::vector< pat::Muon >& tightMuons,
                              const edm::Handle<reco::TrackCollection>& gTracks,
                              const edm::Handle<reco::BeamSpot>& beamspot,
                              double wgt, int nJets, hnu::CutsStruct& cuts,
                              HeavyNu::HistStruct& hists);

   void studyElectronEff(edm::Event& iEvent,
                         std::vector< std::pair<pat::Electron, float> >& eCands,
                         edm::Handle<pat::ElectronCollection>& pElecs,
                         HeavyNuEvent& hnuEvent,
                         std::vector< std::pair<pat::Jet, float> >& jetCands,
                         HeavyNuTrigger *trig_,
                         hnu::CutsStruct& cuts,
                         EffInfo& effinfo,
                         HeavyNu::HistStruct& hists);

   void studyElectronSelectionEff(const std::vector< std::pair<pat::Electron, pat::Electron> >& TP,
                                  const std::vector< std::pair<pat::Electron, float> >& heepElectrons,
                                  double wgt, int nJets, HeavyNu::HistStruct& hists);
}

#endif