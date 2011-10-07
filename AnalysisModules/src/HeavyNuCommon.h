#ifndef HEAVY_NU_COMMON_INCLUDED
#define HEAVY_NU_COMMON_INCLUDED 1

#include <vector>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ESHandle.h"
// #include "FWCore/Framework/interface/Event.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
// #include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Math/VectorUtil.h"
#include <math.h>

//#define CMSSW_3XX
#define CMSSW_4XX

namespace hnu {

  class pTcompareRef
  {
  public:

    template <class T> bool operator() (const T& a, const T& b) {
      return a->pt() > b->pt();
    }
  };

  class pTcompare
  {
  public:

    template <class T> bool operator() (const T& a, const T& b) {
      return a.pt() > b.pt();
    }
    template <class T> bool operator() (const T *a, const T *b) {
      return a->pt() > b->pt();
    }
  };

  class scaleCompare
  {
  public:
    
    template <class T> bool operator() (const T& a, const T& b) {
      return (a.first.pt()*a.second) > (b.first.pt()*b.second);
    }
  };

    bool isVBTFloose       ( const pat::Muon& m );
    bool isVBTFtight       ( const pat::Muon& m );
    double muIsolation     ( const pat::Muon& m, const double scale ) ; 
    
    double getElectronEt    ( const pat::Electron& e ) ; 
    double getElectronSCEta ( const pat::Electron& e ) ; 
    bool passesHEEPv31      ( const pat::Electron& e ) ; 

    int  jetID             ( const pat::Jet& j );

    /* bool passesTrigger     (const double mu1pt, const double mu2pt, */
    /* 			    const bool mu1trig, const bool mu2trig, const uint32_t run) ;  */


    int numberOfPrimaryVertices(edm::Handle<reco::VertexCollection> pvHandle) ;

    double avgVertex(const reco::JPTJet &tjet, double maxDeltaVR=1.0);
    double avgVertex(const pat::Jet &tJet, double maxDeltaVR);
    double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR=1.0);

    float jecTotalUncertainty(float jpt, float jeta,
			      JetCorrectionUncertainty *jecUnc,
			      int correctEra, bool isBjet, bool directionIsUp);

    std::vector< std::pair<pat::Jet,float> >      getJetList (edm::Handle<pat::JetCollection>& pJets,
							      JetCorrectionUncertainty* jecUnc,
							      double minPt, double maxAbsEta,
							      int jecSign=0, int jecEra=3);
    std::vector<pat::Muon>                        getMuonList(edm::Handle<pat::MuonCollection>& pMuons,
							      edm::Handle<reco::MuonCollection>& tevMuons,
							      double minPt, double maxAbsEta,
							      double ptScale=1.0);
    std::vector< std::pair<pat::Electron,float> > getElectronList(edm::Handle<pat::ElectronCollection>& pElecs,
								  double maxAbsEta, 
								  double minPtEB, double minPtEE, 
								  float ebScale=1.0, float eeScale=1.0);
    std::vector<float> generate_flat10_mc(const int npt);

    std::vector<float> get_standard_pileup_data(int pileupEra, const int npt); 

    std::pair<float,double> pileupReweighting(edm::Handle< std::vector<PileupSummaryInfo> > pPU, 
					      const edm::LumiReWeighting mcWeight);
}


#endif // HEAVY_NU_COMMON_INCLUDED
