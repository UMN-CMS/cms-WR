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
    double muIsolation     ( const pat::Muon& m, const double scale=1.0 ) ; 
    
    double getElectronEt    ( const pat::Electron& e ) ; 
    double getElectronSCEta ( const pat::Electron& e ) ; 
    bool passesHEEPv31      ( const pat::Electron& e ) ; 

    int  jetID             ( const pat::Jet& j );

    void initPDFSet(int i, std::string name) ; 
    
    double getPDFWeight(float Q, int id1, float x1, int id2, float x2,
                        bool doPDFreweight, int pdfReweightBaseId, int pdfReweightTargetId) ; 
  
    /* bool passesTrigger     (const double mu1pt, const double mu2pt, */
    /* 			    const bool mu1trig, const bool mu2trig, const uint32_t run) ;  */


    int numberOfPrimaryVertices(edm::Handle<reco::VertexCollection> pvHandle) ;

    double avgVertex(const reco::JPTJet &tjet, double maxDeltaVR=1.0);
    double avgVertex(const pat::Jet &tJet, double maxDeltaVR);
    double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR=1.0);

    float jecTotalUncertainty(float jpt, float jeta,
			      JetCorrectionUncertainty *jecUnc,
			      int correctEra, bool isBjet, bool directionIsUp);
    double muScaleLUT(pat::Muon& iM) ; 

    std::vector< std::pair<pat::Jet,float> >      getJetList (edm::Handle<pat::JetCollection>& pJets,
							      JetCorrectionUncertainty* jecUnc,
							      double minPt, double maxAbsEta,
							      int jecSign=0, int jecEra=3, 
							      bool isMC=false, int jerSign=0);
    std::vector<pat::Muon>                        getMuonList(edm::Handle<pat::MuonCollection>& pMuons,
							      edm::Handle<reco::MuonCollection>& tevMuons,
							      double minPt, double maxAbsEta,
							      double mesScale=1.0, bool merUnc=false,
							      bool trackerPt=false);
    std::vector< std::pair<pat::Electron,float> > getElectronList(edm::Handle<pat::ElectronCollection>& pElecs,
								  double maxAbsEta, 
								  double minPtEB, double minPtEE, 
								  float ebScale=1.0, float eeScale=1.0);
    std::vector<float> generate_flat10_mc();

    std::vector<float> get_standard_pileup_data(int pileupEra); 

    std::pair<float,double> pileupReweighting(const edm::Handle< std::vector<PileupSummaryInfo> >& pPU, 
					      edm::LumiReWeighting& mcWeight);


    // Object of this code is to find all tag and probe pairs for a given event
    // The tag will *always* be in the first position, and the probe collection 
    // should always contain all tags.  If combination x contains a two-tag pair, 
    // some other combination y will contain the same pair in reverse order
    template <class T> std::vector< std::pair<pat::Muon,T> > 
    getTagProbePair(const std::vector<pat::Muon>& tags,const std::vector<T>& probes,
		    double minMass, double maxMass,double rval,bool sanityCheck=true) {

        std::vector< std::pair<pat::Muon,T> > tagprobes ; 
    
        // Sanity check to make sure tag collection is subset of probes
        std::vector<unsigned int> tagLocs ; 
        for (unsigned int i=0; i<tags.size(); i++) { 
	  double minDR = 999. ; 
	  for (unsigned int j=0; j<probes.size(); j++) { 
	    double dR = ROOT::Math::VectorUtil::DeltaR(tags.at(i).p4(),probes.at(j).p4()) ; 
	    if ( dR < minDR ) minDR = dR ; 
	    if ( dR < 0.02 ) { // Found the tag in the probe list
	      tagLocs.push_back(j) ;
	      break ; 
	    }
	  }
	  if ( sanityCheck && minDR > 0.02 ) 
	    std::cout << "WARNING: Tag does not match any probe candidate, closest is dR = " << minDR << std::endl ; 
        }
        if ( sanityCheck && tagLocs.size() != tags.size() ) 
            std::cout << "WARNING!!! Expected to find all tags in probe list, but for " 
                      << tags.size() << " tags found " << tagLocs.size() << " in probe list" 
                      << std::endl ; 
        
        for (unsigned int i=0; i<tags.size(); i++) { 
            for (unsigned int j=0; j<probes.size(); j++) { 
	      // By skipping known tags, tag+probe sample is biased towards failing probes
	      // if (std::find(tagLocs.begin(),tagLocs.end(),j) != tagLocs.end()) continue ; 
	      double mass = (tags.at(i).p4()+probes.at(j).p4()).M() ; 
	      if ( mass < minMass || mass > maxMass ) continue ; 
	      std::pair<pat::Muon,T> tp = std::make_pair(tags.at(i),probes.at(j)) ; 
	      tagprobes.push_back( tp ) ; 
            }
        }
    
        // Return all tag+probe combinations
        if ( rval < 0 || tagprobes.size() <= 1 ) return tagprobes ; 
        // Randomly pick one tag+probe combination from those available
        unsigned int idx = (unsigned int)(rval * tagprobes.size()) ; 
        std::vector< std::pair<pat::Muon,T> > singletp ; 
        singletp.push_back( tagprobes.at(idx) ) ; 
        return singletp ; 
    }

}


#endif // HEAVY_NU_COMMON_INCLUDED
