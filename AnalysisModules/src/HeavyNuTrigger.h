/** -- C++ -- **/
#ifndef HeavyNuTrigger_h_included
#define HeavyNuTrigger_h_included

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"

class HeavyNuTrigger {
 public:

  explicit HeavyNuTrigger(const edm::ParameterSet & iConfig);
  void book(const TFileDirectory& tdir);

  inline bool matchingEnabled() { return matchingEnabled_; }

  bool isTriggerMatched(const pat::MuonRef & muon,
			const edm::Event   & iEvent);

  void endJob();

 private:

  bool matchingEnabled_;

  TH2D *trigMatchPtCorrel;
  TH1D *trigMatchDR2;
  TH2D *trigMatchDRDPt;
  TH2D *trigMatchDetaPhi;

  edm::InputTag trigEventTag_;
  std::string   muonMatch_;
};

#endif
