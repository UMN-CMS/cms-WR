#include <TF1.h>

#include "HeavyNuCommon.h"
#include "TVector3.h"
#include "TRandom.h"

#ifdef DO_LHAPDF
#include "LHAPDF/LHAPDF.h"
#endif

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

// const float pileup2010A = 1.2;
// const float pileup2010B = 2.2;
// const float pileup2011A = 5.0;

// Data corrections provided by N. Kypreos.
// MC corrections are inverse of data, and given below
const double muScaleLUTarray100[5] = { 0.379012,0.162107,0.144514,0.0125131,-0.392431 } ; 

namespace hnu {

  bool etCompare::operator() (const std::pair<pat::Electron,float>& a, const std::pair<pat::Electron,float>& b) {
    return (getElectronEt(a.first,false)*a.second) > (getElectronEt(b.first,false)*b.second);
  }

  bool isLooseMuonNoPF(const pat::Muon& m)
  {
    return (m.isGlobalMuon() || m.isTrackerMuon());
  }
    
  bool isLooseMuon(const pat::Muon& m)
  {
    return (m.isPFMuon() && isLooseMuonNoPF(m));
  }

  // This selection takes care of what is common to "Tight" and "High pT Tight" muon selection
  // What is *not* included: PF muon, normalized chi^2, dpT/pT, dxy, dz
  bool isTightMuonCore(const pat::Muon& m)
  {
    reco::TrackRef gt = m.globalTrack();
    reco::TrackRef it = m.innerTrack();
    if (gt.isNull() || it.isNull()) return false;

    return ( (m.isGlobalMuon()) &&
	     (gt->hitPattern().numberOfValidMuonHits() > 0) &&
	     (it->hitPattern().numberOfValidPixelHits() > 0) &&
	     (m.numberOfMatchedStations() > 1) &&
	     (it->hitPattern().trackerLayersWithMeasurement() > 5) ); 
  }
    
  bool isTightMuon(const pat::Muon& m, edm::Handle<reco::VertexCollection> pvHandle)
  {
    if ( !isTightMuonCore(m) ) return false;
    
    if ( m.muonBestTrack().isNull() ) return false;
    if ( !pvHandle.isValid() ) return false;
    reco::Vertex pv = pvHandle->at(0) ; 
   
    return ( m.isPFMuon() && (m.normChi2() < 10) &&
	     (fabs(m.muonBestTrack()->dxy(pv.position())) < 0.2 ) &&
	     (fabs(m.muonBestTrack()->dz(pv.position())) < 0.5) );
  } // HeavyNu::isTightMuon

  // Information taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  // On 19 Nov, 2012
  bool isTightHighPtMuon(const pat::Muon& m, edm::Handle<reco::VertexCollection> pvHandle)
  { 
    if ( !isTightMuonCore(m) ) return false;    

    reco::TrackRef cktTrack = (muon::tevOptimized(m, 200, 40., 17., 0.25)).first;
    if (cktTrack.isNull()) return false;
    if ( !pvHandle.isValid() ) return false;
    reco::Vertex pv = pvHandle->at(0) ; 

    return ( (fabs(cktTrack->dxy(pv.position())) < 0.2 ) &&
	     (fabs(cktTrack->dz(pv.position())) < 0.5) &&
	     (cktTrack->ptError()/cktTrack->pt() < 0.3) ); 
  }

  double getElectronEt(const pat::Electron& e, bool useCorrectedEnergy) { 
    
    double ET     = e.superCluster()->energy() * sin( e.p4().theta() ) ; 
    double corrET = e.caloEnergy() * sin( e.p4().theta() ) ; 

    return ( useCorrectedEnergy ? corrET : ET ) ; 
  } 

  double getElectronSCEta(const pat::Electron& e) { 
    return e.caloPosition().eta() ; 
  }

  double getElectronSCPhi(const pat::Electron& e) {
    return e.caloPosition().phi() ;
  }

  // Updated to HEEP v4.1 ID and new requirements for electron fakes
  bool passesFakeRateMinimum(const pat::Electron& e, edm::Handle<reco::VertexCollection> pvHandle) { 

    if ( !e.ecalDriven() ) return false ; 

    double eEta = getElectronSCEta(e) ; 
    bool isEB   = ( fabs(eEta) < 1.442 ) ; 
    bool isEE   = ( fabs(eEta) < 2.5 && fabs(eEta) > 1.56 ) ; 
    if ( !isEB && !isEE ) return false ; 

    double HoE      = e.hadronicOverEm() ; 
    double sig_iEiE = e.sigmaIetaIeta() ; 
    int nLostHits   = e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ; 
    double absdxy   = (pvHandle.isValid()) ? fabs(e.gsfTrack()->dxy(pvHandle->at(0).position())) : 999.0 ; 

    if ( nLostHits > 1 ) return false ; 
    if ( isEB ) { 
      if ( HoE > 0.15 )       return false ; 
      if ( sig_iEiE > 0.013 ) return false ;
      if ( absdxy > 0.02 )    return false ; 
    } else if ( isEE ) { 
      if ( HoE > 0.10 )       return false ; 
      if ( sig_iEiE > 0.034 ) return false ; 
      if ( absdxy > 0.05 )    return false ; 
    }

    return true ; 
  } 

  double fakeProbability(const pat::Electron& e) { 

    double prob = 0. ; 

    double eEta = getElectronSCEta(e) ; 
    double ePt  = getElectronEt(e,false) ; 

    // Values taken from AN/2012-415 (Z' -> ee), v10
    if ( fabs(eEta) < 1.442 ) { 
      prob =  (( ePt < 189.3 ) ?  ( 0.0179 - 0.000056 * ePt ) : ( 0.0073 )) ; 
    } else if ( fabs(eEta) < 2.0 && fabs(eEta) > 1.56 ) {
      if ( ePt < 96.6 ) prob = exp(-2.31 - 0.011*ePt) ;
      else prob = (( ePt < 178.0 ) ? ( 0.040 - 0.000059 * ePt ) : ( 0.0295 )) ; 
    } else if ( fabs(eEta) < 2.5 && fabs(eEta) > 2.0 ) { 
      prob = (( ePt < 115.4 ) ?  ( 0.099 - 0.00035 * ePt ) : ( 0.0586 )) ; 
    }

    return prob ; 
  }

  // Fake rate calculated for loose muons, within 0.3 of a valid jet, passing ID/iso selection
  double fakeProbability(const pat::Muon& mu) { 

    double eta  = mu.eta() ;
    double pt   = mu.pt() ; 
    double prob = 0. ; 

    if ( abs(eta) < 2.4 ) { // Determined April 22, 2013 using 15/fb of 8 TeV data
      if ( pt >= 100.0 ) {
        if      ( abs(eta) > 2.1 ) prob = 0.194 ; 
        else if ( abs(eta) > 1.2 ) prob = 0.156 ; 
        else                       prob = 0.088 ;
      } else if ( pt >= 40.0 ) {
        if      ( abs(eta) > 2.1 ) prob = 0.199 ; 
        else if ( abs(eta) > 1.2 ) prob = 0.118 ; 
        else                       prob = 0.078 ;
      }
    }

    return prob ; 
  }

  bool passesHEEP(const pat::Electron& e, int heepVersion, double rho, edm::Handle<reco::VertexCollection> pvHandle) { 

    if ( heepVersion != 31 && heepVersion != 32 && abs(heepVersion) != 40  && abs(heepVersion) != 41) return false ; 

    double ePt  = getElectronEt(e,(abs(heepVersion) != 40) && (abs(heepVersion) != 41)) ; 
    // if ( ePt < 35.0 ) std::cout << "Removing low pT electron from consideration" << std::endl ; 
    // if ( ePt < 35.0 ) return false ; 

    // std::cout << "heepVersion: " << heepVersion << " and rho " << rho << std::endl ; 
    
    // All electrons must be ECAL driven
    // if ( !e.ecalDriven() ) std::cout << "FAILURE: electron is not ecal driven!!!" << std::endl ; 
    if ( !e.ecalDriven() ) return false ;

    double eEta = getElectronSCEta(e) ; 
    bool isEB   = ( fabs(eEta) < 1.442 ) ; 
    bool isEE   = ( fabs(eEta) < 2.5 && fabs(eEta) > 1.56 ) ; 
    // std::cout << "Electron pT: " << ePt << std::endl ; 
    // std::cout << "Electron SC eta: " << eEta << std::endl ; 
    // if ( !isEB && !isEE) std::cout << "FAILURE: electron is not in EB or EE: " << eEta << std::endl ; 
    if ( !isEB && !isEE ) return false ; 
    
    double HoE = e.hadronicOverEm() ; 
    // std::cout << "Electron HoE: " << HoE << std::endl ; 
    // if ( HoE > 0.05 ) std::cout << "FAILURE: HOE = " << HoE << std::endl ; 
    if ( HoE > 0.05 ) return false ; 

    double dEtaIn = fabs(e.deltaEtaSuperClusterTrackAtVtx()) ; 
    // std::cout << "Electron dEtaIn: " << dEtaIn << std::endl ; 
    // if ( ( isEB && dEtaIn > 0.005 ) || ( isEE && dEtaIn > 0.007 ) ) std::cout << "FAILURE: Electron dEtaIn: " << dEtaIn << std::endl ; 
    if ( isEB && dEtaIn > 0.005 ) return false ; 
    if ( isEE && dEtaIn > 0.007 ) return false ; 

    double dPhiIn = fabs(e.deltaPhiSuperClusterTrackAtVtx()) ; 
    // std::cout << "Electron dPhiIn: " << dPhiIn << std::endl ; 
    // if ( heepVersion == 40 && dPhiIn > 0.06 ) std::cout << "FAILURE: Electron dPhiIn = " << dPhiIn << std::endl ; 
    if (  heepVersion == 31                                                      && dPhiIn > 0.09 ) return false ; 
    if ( (heepVersion == 32 || abs(heepVersion) == 40 || abs(heepVersion) == 41) && dPhiIn > 0.06 ) return false ; 

    double sig_iEiE = e.sigmaIetaIeta() ; 
    // if ( isEE ) std::cout << "Electron sigiEiE: " << sig_iEiE << std::endl ; 
    // if ( isEE && sig_iEiE > 0.03 ) std::cout << "FAILURE: Electron sigiEiE = " << sig_iEiE << std::endl ; 
    if ( isEE && sig_iEiE > 0.03 ) return false ; 

    double e1x5_5x5 = e.e1x5() / e.e5x5() ; 
    double e2x5_5x5 = e.e2x5Max() / e.e5x5() ; 
    // if ( isEB ) std::cout << "Electron 1x5, 2x5: " << e1x5_5x5 << ", " << e2x5_5x5 << std::endl ; 
    // if ( isEB && (e2x5_5x5 < 0.94 && e1x5_5x5 < 0.83) ) std::cout << "FAILURE: Electron 1x5, 2x5 = " << e1x5_5x5 << ", " << e2x5_5x5 << std::endl ; 
    if ( isEB && (e2x5_5x5 < 0.94 && e1x5_5x5 < 0.83) ) return false ; 

    int nLostHits = e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ; 
    // std::cout << "Electron nLostHits: " << nLostHits << std::endl ; 
    // if ( nLostHits != 0 ) std::cout << "FAILURE: Electron nLostHits = " << nLostHits << std::endl ; 
    if( heepVersion == 31 || heepVersion == 32 || heepVersion == 40) if ( nLostHits != 0 ) return false ; 
    if( heepVersion == 41) if ( nLostHits > 1 ) return false ; 
    
    double ecalHcalIso = e.dr03EcalRecHitSumEt() + e.dr03HcalDepth1TowerSumEt() ; 
    double thresholdEB = 2. + 0.03 * ePt ; 
    double thresholdEE = 2.5 ;
    if ( isEE && ePt > 50 ) thresholdEE += (0.03 * (ePt - 50.)) ; 
    if ( abs(heepVersion) == 40 || abs(heepVersion) == 41 ) { 
      thresholdEB += ( 0.28 * rho ) ; 
      thresholdEE += ( 0.28 * rho ) ; 
    }
    // std::cout << "Electron caloIso: " << ecalHcalIso << " = " << e.dr03EcalRecHitSumEt() << " + " << e.dr03HcalDepth1TowerSumEt() 
    // 	      << std::endl ; 
    // std::cout << "Electron threshold: " << thresholdEB << " (EB) " << thresholdEE << " (EE) " << std::endl ; 
    
    // if ( isEB && ecalHcalIso > thresholdEB ) std::cout << "FAILURE: EB and caloIso" << std::endl ; 
    // if ( isEE && ecalHcalIso > thresholdEE ) std::cout << "FAILURE: EE and caloIso" << std::endl ; 

    if ( isEB && ecalHcalIso > thresholdEB ) return false ; 
    if ( isEE && ecalHcalIso > thresholdEE ) return false ; 

    double hcalIsoDepth2 = e.dr03HcalDepth2TowerSumEt() ; 
    if ( (heepVersion == 31) && isEE && hcalIsoDepth2 > 0.5 ) return false ; 

    double trkIso = e.dr03TkSumPt() ;
    // if ( heepVersion == 40 && trkIso > 5. ) std::cout << "FAILURE: Electron trkIso = " << trkIso << std::endl ; 
    // std::cout << "Electron trkIso: " << trkIso << std::endl ; 
    if ( (heepVersion == 31) && ((isEB && trkIso > 7.5) || (isEE && trkIso > 15.)) ) return false ; 
    if ( (heepVersion == 32 || abs(heepVersion) == 40 || abs(heepVersion) == 41) && (trkIso > 5.) ) return false ; 
    
    if(pvHandle.isValid())
    {
      double absdxy = fabs(e.gsfTrack()->dxy(pvHandle->at(0).position()));
      if(isEB && absdxy > 0.02) return false;
      if(isEE && absdxy > 0.05) return false;
    } else {
      if (abs(heepVersion) == 41) return false ; 
    }
    
    // std::cout << "Congratulations!  The electron passes HEEP selection" << std::endl ; 

    // Passes HEEP selection
    return true ; 
  }


  bool passesHEEPv31(const pat::Electron& e) { 
    if ( !e.ecalDriven() ) return false ; 

    double ePt  = getElectronEt(e,true) ; 
    double eEta = getElectronSCEta(e) ; 

    // common requirements 
    bool HoE    = (e.hadronicOverEm() < 0.05) ; 
    bool dPhiIn = (fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.09) ; 
    if ( !HoE || !dPhiIn ) return false ; 
    int nLostHits = e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ; 
    if ( nLostHits > 0 ) return false ; // conversion rejection

    double ecalIso  = e.dr03EcalRecHitSumEt() ; 
    double hcalIso1 = e.dr03HcalDepth1TowerSumEt() ; 
    double hcalIso2 = e.dr03HcalDepth2TowerSumEt() ; 
    double trkIso   = e.dr03TkSumPt() ; 

    if ( fabs(eEta) < 1.442 ) { // Barrel HEEP electron candidate
      double e15 = e.e1x5() / e.e5x5() ; 
      double e25 = e.e2x5Max() / e.e5x5() ; 
      bool shape = ( e15 > 0.83 ) || ( e25 > 0.94 ) ; 
      bool dEta  = ( fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.005 ) ; 
      bool isolated  = ( (ecalIso + hcalIso1) < (2. + 0.03*ePt) ) && 
	( trkIso < 7.5 ) ; 
      return ( shape && dEta && isolated ) ; 
    } else if ( fabs(eEta) > 1.560 ) { // Endcap HEEP electron candidate 
      bool shape = ( e.sigmaIetaIeta() < 0.03 ) ; 
      bool dEta  = ( fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 ) ; 
      double threshold = ( ePt < 50. ) ? (2.5) : (2.5 + 0.03 * (ePt-50)) ; 
      bool isolated = ((ecalIso + hcalIso1) < threshold) && (hcalIso2 < 0.5) && (trkIso < 15.) ; 
      return ( shape && dEta && isolated ) ; 
    }
    return false ; 
  }

  bool passesHEEPv32(const pat::Electron& e) { 
    if ( !e.ecalDriven() ) return false ; 
    double ePt  = getElectronEt(e,true) ; 
    double eEta = getElectronSCEta(e) ; 

    // common requirements 
    bool HoE    = (e.hadronicOverEm() < 0.05) ; 
    bool dPhiIn = (fabs(e.deltaPhiSuperClusterTrackAtVtx()) < 0.06) ; 
    if ( !HoE || !dPhiIn ) return false ; 
    int nLostHits = e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ; 
    if ( nLostHits > 0 ) return false ; // conversion rejection

    double ecalIso = e.dr03EcalRecHitSumEt() ; 
    double hcalIso = e.dr03HcalDepth1TowerSumEt() ; 
    double trkIso  = e.dr03TkSumPt() ; 
    if ( trkIso > 5.0 ) return false ; 

    if ( fabs(eEta) < 1.442 ) { // Barrel HEEP electron candidate
      double e15 = e.e1x5() / e.e5x5() ; 
      double e25 = e.e2x5Max() / e.e5x5() ; 
      bool shape = ( e15 > 0.83 ) || ( e25 > 0.94 ) ; 
      bool dEta  = ( fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.005 ) ; 
      bool isolated  = ( (ecalIso + hcalIso) < (2. + 0.03*ePt) ) ; 
      return ( shape && dEta && isolated ) ; 
    } else if ( fabs(eEta) > 1.560 ) { // Endcap HEEP electron candidate 
      bool shape = ( e.sigmaIetaIeta() < 0.03 ) ; 
      bool dEta  = ( fabs(e.deltaEtaSuperClusterTrackAtVtx()) < 0.007 ) ; 
      double threshold = ( ePt < 50. ) ? (2.5) : (2.5 + 0.03 * (ePt-50)) ; 
      bool isolated = ((ecalIso + hcalIso) < threshold) ; 
      return ( shape && dEta && isolated ) ; 
    }
    return false ; 
  }

  std::vector< std::pair<pat::Electron,pat::Electron> > 
  getTagProbePair(const std::vector<pat::Electron>& tags,const std::vector<pat::Electron>& probes,
		  double minMass, double maxMass,double rval,bool sanityCheck) {

    std::vector< std::pair<pat::Electron,pat::Electron> > tagprobes ; 
    
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
	// std::cout << "Tag " << (i+1) << " of " << tags.size() << " and probe " << (j+1) 
	// 	  << " of " << probes.size() << " has a mass of " << mass << std::endl ; 
	// if ( mass < minMass || mass > maxMass ) std::cout << "This tag/probe combination is rejected" << std::endl ;
	if ( mass < minMass || mass > maxMass ) continue ; 

	// double dR = ROOT::Math::VectorUtil::DeltaR(tags.at(i).p4(),probes.at(j).p4()) ; 
	// if ( dR < 0.3 ) std::cout << "The tag and probe are too close: " << dR << std::endl ; 

	std::pair<pat::Electron,pat::Electron> tp = std::make_pair(tags.at(i),probes.at(j)) ; 
	tagprobes.push_back( tp ) ; 
      }
    }
    
    // Return all tag+probe combinations
    if ( rval < 0 || tagprobes.size() <= 1 ) return tagprobes ; 
    // Randomly pick one tag+probe combination from those available
    unsigned int idx = (unsigned int)(rval * tagprobes.size()) ; 
    std::vector< std::pair<pat::Electron,pat::Electron> > singletp ; 
    singletp.push_back( tagprobes.at(idx) ) ; 
    return singletp ; 
  }

  double muIsolation(const pat::Muon& m, const double pTscale) {
    double mupt = pTscale * m.pt() ; 
    return (m.trackIso() / mupt) ; 
  }

  int jetID(const pat::Jet& j) {
    int val = -1 ; 

    if (j.isPFJet()) { 
      PFJetIDSelectionFunctor jetIDloose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
      PFJetIDSelectionFunctor jetIDtight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);

      pat::strbitset ret = jetIDloose.getBitTemplate();
      ret.set(false); bool loose  = jetIDloose (j, ret);
      // ret.set(false); bool medium = jetIDmedium(j, ret);
      ret.set(false); bool tight  = jetIDtight (j, ret);
      if (tight)       val = 3 ; 
      // else if (medium) val = 2 ; 
      else if (loose)  val = 1 ;
      else             val = 0 ; 
    } else { // Calo jets (deprecated)
      JetIDSelectionFunctor jetIDloose(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);
      JetIDSelectionFunctor jetIDtight(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT);

      pat::strbitset ret = jetIDloose.getBitTemplate();
      ret.set(false); bool loose = jetIDloose(j, ret);
      ret.set(false); bool tight = jetIDtight(j, ret);
      if ( tight )      val = 3 ; 
      else if ( loose ) val = 1 ; 
      else              val = 0 ; 
    }
    return val ; 
  }

  void initPDFSet(int i, std::string name) {
#ifdef DO_LHAPDF
      LHAPDF::initPDFSet(i,name) ;
#endif
  }
    
  double getPDFWeight(float Q, int id1, float x1, int id2, float x2,
                      bool doPDFreweight, int pdfReweightBaseId, int pdfReweightTargetId) {
  
      if (!doPDFreweight) return 1.0;

      double w = 1.0;

#ifdef DO_LHAPDF
      LHAPDF::usePDFMember(1,pdfReweightBaseId);
      double pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
      double pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  
      LHAPDF::usePDFMember(2,pdfReweightTargetId);
      double newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
      double newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
      
      w=(newpdf1/pdf1*newpdf2/pdf2);
#endif
  
      //  printf("My weight is %f\n",w);
  
      return w;
  }
    
  int numberOfPrimaryVertices(edm::Handle<reco::VertexCollection> pvHandle) { 
    int nvertex = 0 ; 
    
    const reco::VertexCollection& vertices = *pvHandle.product();
    static const int minNDOF = 4;
    static const double maxAbsZ = 15.0;
    static const double maxd0 = 2.0;

    for (reco::VertexCollection::const_iterator vit=vertices.begin(); 
	 vit!=vertices.end(); vit++) {
      if (vit->ndof() > minNDOF && 
	  (fabs(vit->z()) <= maxAbsZ) && 
	  (fabs(vit->position().rho()) <= maxd0)) nvertex++;
    }
    return nvertex ; 
  }

  double avgVertex(const reco::JPTJet &tjet, double maxDeltaVR) {
      double avetz=0, weight=0;
      double thresh=0.0001;
      reco::TrackRefVector::const_iterator itrack;
      for (itrack=tjet.getPionsInVertexInCalo().begin(); itrack!=tjet.getPionsInVertexInCalo().end(); itrack++) {
          if ( ( maxDeltaVR > thresh && sqrt( pow((*itrack)->vx(),2) + pow((*itrack)->vy(),2) ) < maxDeltaVR ) || maxDeltaVR < 0.0 ) {
              avetz+=(*itrack)->vz()/(*itrack)->normalizedChi2();
              weight+=1./(*itrack)->normalizedChi2();
          }
      }
      for (itrack=tjet.getPionsInVertexOutCalo().begin(); itrack!=tjet.getPionsInVertexOutCalo().end(); itrack++) {
          if ( ( maxDeltaVR > thresh && sqrt( pow((*itrack)->vx(),2) + pow((*itrack)->vy(),2) ) < maxDeltaVR ) || maxDeltaVR < 0.0 ) {
              avetz+=(*itrack)->vz()/(*itrack)->normalizedChi2();
              weight+=1./(*itrack)->normalizedChi2();
          }
      }
      for (itrack=tjet.getPionsOutVertexInCalo().begin(); itrack!=tjet.getPionsOutVertexInCalo().end(); itrack++) {
          if ( ( maxDeltaVR > thresh && sqrt( pow((*itrack)->vx(),2) + pow((*itrack)->vy(),2) ) < maxDeltaVR ) || maxDeltaVR < 0.0 ) {
              avetz+=(*itrack)->vz()/(*itrack)->normalizedChi2();
              weight+=1./(*itrack)->normalizedChi2();
          }
      }

      if (weight<0.01) { // throw out weight ~= 0 results
          return -100;
      } else {
          return avetz/weight;
      }
  }

    double avgVertex(const pat::Jet &tJet, double maxDeltaVR)
    {
        double avetz = 0, weight = 0;
        double thresh = 0.0001;
        reco::PFCandidateFwdPtrVector::const_iterator itrack;
        for(itrack = tJet.pfCandidatesFwdPtr().begin(); itrack != tJet.pfCandidatesFwdPtr().end(); itrack++)
        {
            if(((maxDeltaVR > thresh && sqrt(pow((*itrack)->vx(), 2) + pow((*itrack)->vy(), 2)) < maxDeltaVR) || maxDeltaVR < 0.0) && abs((*itrack)->pdgId()) == 211)
            {
                avetz += (*itrack)->vz() / (*itrack)->vertexChi2();
                weight += 1. / (*itrack)->vertexChi2();
            }
        }

        if(weight < 0.01)
        { // throw out weight ~= 0 results
            return -100;
        }
        else
        {
            return avetz / weight;
        }
    }


  double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR){

        // Find best match jptjets
        reco::JPTJetCollection::const_iterator tJet=jptJets.end();
        for (reco::JPTJetCollection::const_iterator i=jptJets.begin(); i != jptJets.end(); i++){
            if (tJet == jptJets.end() || ROOT::Math::VectorUtil::DeltaR(pJet.p4(),i->p4()) < ROOT::Math::VectorUtil::DeltaR(pJet.p4(),tJet->p4())){
                tJet =i;
            }
        }

        // Return Vert
	if ( tJet != jptJets.end() ) return avgVertex(*tJet, maxDeltaVR);
	return -100. ; 
  }

  std::vector<float> get_standard_pileup_mc(int era){
    const double npu_probs_summer12[60] = { 2.344E-05, 2.344E-05, 2.344E-05, 2.344E-05, 4.687E-04, // 0-4
					    4.687E-04, 7.032E-04, 9.414E-04, 1.234E-03, 1.603E-03, // 5-9
					    2.464E-03, 3.250E-03, 5.021E-03, 6.644E-03, 8.502E-03, // 10-14
					    1.121E-02, 1.518E-02, 2.033E-02, 2.608E-02, 3.171E-02, // 15-19
					    3.667E-02, 4.060E-02, 4.338E-02, 4.520E-02, 4.641E-02, // 20-24
					    4.735E-02, 4.816E-02, 4.881E-02, 4.917E-02, 4.909E-02, // 25-29
					    4.842E-02, 4.707E-02, 4.501E-02, 4.228E-02, 3.896E-02, // 30-34
					    3.521E-02, 3.118E-02, 2.702E-02, 2.287E-02, 1.885E-02, // 35-39
					    1.508E-02, 1.166E-02, 8.673E-03, 6.190E-03, 4.222E-03, // 40-44
					    2.746E-03, 1.698E-03, 9.971E-04, 5.549E-04, 2.924E-04, // 45-49
					    1.457E-04, 6.864E-05, 3.054E-05, 1.282E-05, 5.081E-06, // 50-54
					    1.898E-06, 6.688E-07, 2.221E-07, 6.947E-08, 2.047E-08  // 55-59
    }; 
    
    const double Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04,
                                       2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02,
                                       1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02,
                                       4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02,
                                       5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02,
                                       3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02,
                                       2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02,
                                       1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03,
                                       4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03,
                                       1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04,
                                       2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05,
                                       4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06
    };

    const double* npu_probs = npu_probs_summer12 ;
    int npt = sizeof(npu_probs_summer12)/sizeof(double);
    if(era == 20121)
    {
        npu_probs = npu_probs_summer12 ;
        npt = sizeof(npu_probs_summer12)/sizeof(double);
    }
    else if(era == 20122)
    {
        npt = sizeof(Summer2012_S10)/sizeof(double);
        npu_probs = Summer2012_S10 ;
    }
    else std::cout << "!!!INVALID PILEUP ERA!!! " << era << std::endl;
        

    std::vector<float> retval;
    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(npu_probs[i]);
    return retval;
  }


  std::vector<float> get_standard_pileup_data(int era, int systShift) {
    
    const double json_2012ab_run195868[] = { 
         6499.89,    8805.91,    86264.8, 4.14212E06, 5.65696E06, 
       2.1853E06, 1.11201E07, 2.88172E07, 5.71031E07,  8.9966E07, 
      1.21335E08,  1.6028E08, 2.06592E08, 2.50997E08, 2.95318E08, 
      3.25967E08, 3.17726E08, 2.76024E08, 2.29525E08, 1.93205E08, 
      1.66461E08, 1.45222E08, 1.26966E08, 1.10733E08, 9.59231E07, 
      8.21369E07, 6.92822E07,  5.7453E07, 4.67855E07, 3.73841E07, 
      2.92965E07,  2.2509E07,  1.6952E07, 1.25127E07, 9.05136E06, 
      6.41617E06, 4.45669E06, 3.03318E06, 2.02257E06, 1.32128E06, 
          845530,     529977,     325328,     195551,     115083, 
         66298.5,    37383.3,    20628.8,    11138.9,    5884.91, 
         3041.82,    1538.17,    760.925,    368.254,    174.355, 
         80.7658,    36.6066,    16.2357,    7.04718,    2.99396
    } ; 

    const double json_2012ABCD[] = { //Full 2012 "True" pileup distribution (69.3 mb min bias xsec)
        1.226122e+04, 3.285491e+04, 8.971792e+04, 2.624664e+05, 5.464699e+05, 3.037012e+06, 1.750948e+07, 5.160423e+07, 1.219667e+08, 2.468959e+08,
        4.348696e+08, 6.771360e+08, 8.764220e+08, 9.921868e+08, 1.066224e+09, 1.121976e+09, 1.154995e+09, 1.164012e+09, 1.157054e+09, 1.138393e+09,
        1.110421e+09, 1.076601e+09, 1.038365e+09, 9.911932e+08, 9.270334e+08, 8.415625e+08, 7.378647e+08, 6.232810e+08, 5.055952e+08, 3.925056e+08,
        2.912031e+08, 2.065269e+08, 1.399708e+08, 9.045642e+07, 5.556643e+07, 3.237016e+07, 1.788073e+07, 9.391192e+06, 4.715788e+06, 2.282339e+06,
        1.075819e+06, 5.003868e+05, 2.333164e+05, 1.110043e+05, 5.479941e+04, 2.839548e+04, 1.548825e+04, 8.844905e+03, 5.236186e+03, 3.180093e+03,
        1.964044e+03, 1.225146e+03, 7.677785e+02, 4.812791e+02, 3.006442e+02, 1.865584e+02, 1.146868e+02, 6.969379e+01, 4.179288e+01, 2.469788e+01
    };

    const double json_2012ABCD_systDn_5percent[] = { //Full 2012 "True" pileup distribution with 0.95*(69.3 mb min bias xsec)
        1.308767e+04, 4.509325e+04, 1.031247e+05, 3.338053e+05, 7.295486e+05, 5.880321e+06, 2.810640e+07, 7.761540e+07, 1.802559e+08, 3.494281e+08,
        5.971892e+08, 8.504573e+08, 1.010027e+09, 1.101872e+09, 1.168902e+09, 1.211303e+09, 1.225100e+09, 1.218943e+09, 1.198922e+09, 1.167964e+09, 
        1.130255e+09, 1.087150e+09, 1.032406e+09, 9.564408e+08, 8.558476e+08, 7.362674e+08, 6.071279e+08, 4.779073e+08, 3.579670e+08, 2.550159e+08, 
        1.728272e+08, 1.112440e+08, 6.777516e+07, 3.895835e+07, 2.110900e+07, 1.080776e+07, 5.259919e+06, 2.455452e+06, 1.112790e+06, 4.970865e+05, 
        2.229699e+05, 1.025244e+05, 4.924725e+04, 2.499837e+04, 1.341830e+04, 7.552494e+03, 4.403560e+03, 2.629302e+03, 1.593048e+03, 9.726153e+02, 
        5.950851e+02, 3.631980e+02, 2.202578e+02, 1.322851e+02, 7.847219e+01, 4.588065e+01, 2.639711e+01, 1.492736e+01, 8.289706e+00, 4.518162e+00
    };
               
    const double json_2012ABCD_systUp_5percent[] = { //Full 2012 "True" pileup distribution with 1.05*(69.3 mb min bias xsec)
        1.153146e+04, 2.287525e+04, 8.152293e+04, 2.089707e+05, 4.586210e+05, 1.591248e+06, 1.056854e+07, 3.478318e+07, 8.302347e+07, 1.745792e+08, 
        3.155269e+08, 5.165381e+08, 7.319414e+08, 8.819819e+08, 9.691806e+08, 1.031106e+09, 1.077418e+09, 1.102949e+09, 1.108498e+09, 1.101083e+09, 
        1.083675e+09, 1.058277e+09, 1.027777e+09, 9.935883e+08, 9.524635e+08, 8.977388e+08, 8.247911e+08, 7.348082e+08, 6.333137e+08, 5.267692e+08, 
        4.214962e+08, 3.237556e+08, 2.386699e+08, 1.689066e+08, 1.146427e+08, 7.445013e+07, 4.613535e+07, 2.723771e+07, 1.532772e+07, 8.245078e+06, 
        4.261362e+06, 2.131331e+06, 1.041058e+06, 5.022613e+05, 2.426103e+05, 1.191326e+05, 6.035979e+04, 3.190701e+04, 1.767355e+04, 1.022785e+04, 
        6.136158e+03, 3.781529e+03, 2.374045e+03, 1.508239e+03, 9.646131e+02, 6.184683e+02, 3.961197e+02, 2.526689e+02, 1.600844e+02, 1.005201e+02
    };
    
    const double* pileupDist=json_2012ABCD;
    int npt = 60;
    if (era == 20121) 
    {
        npt = sizeof(json_2012ab_run195868)/sizeof(double);
        pileupDist=json_2012ab_run195868;
    }
    else if (era == 20122) 
    {
        if(systShift == 0)
        {
            npt = sizeof(json_2012ABCD)/sizeof(double);
            pileupDist=json_2012ABCD;
        }
        else if(systShift > 0)
        {
            npt = sizeof(json_2012ABCD_systUp_5percent)/sizeof(double);
            pileupDist=json_2012ABCD_systUp_5percent;
        }
        else
        {
            npt = sizeof(json_2012ABCD_systDn_5percent)/sizeof(double);
            pileupDist=json_2012ABCD_systDn_5percent;
        }
    }
    else std::cout << "!!!INVALID PILEUP ERA!!! " << era << std::endl;

    std::vector<float> retval;

    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(pileupDist[i]);
    return retval;
  }

  std::pair<float,double> pileupReweighting(const edm::Handle<std::vector<PileupSummaryInfo> >& pPU,
					    edm::LumiReWeighting& puWeight) { 

    int   nPileup = -1 ; 
    double weight = 1.0 ; 
    // float  avg_nvtx = -1. ;

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for (PVI = pPU->begin(); PVI != pPU->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();

      if (BX == 0) { 
        // nPileup = PVI->getPU_NumInteractions();
        nPileup = PVI->getTrueNumInteractions();
        continue;
      }
    }
    weight = puWeight.weight( nPileup );
    
//     if (pPU.isValid() && pPU->size() > 0) {
//       std::vector<PileupSummaryInfo>::const_iterator puIter ; 

//       float sum_nvtx = 0 ; 
//       for (puIter=pPU->begin(); puIter!=pPU->end(); puIter++) 
// 	sum_nvtx += float( puIter->getPU_NumInteractions() ) ; 

//       avg_nvtx = sum_nvtx / 3. ; 

//       // std::cout << "About to look up weight for " << avg_nvtx << " vertices, initial weight " << weight << std::endl ; 
//       weight = puWeight.weight3BX( avg_nvtx ) ; 
//       // if ( nPileup > int(mcWeight.size()) || nPileup < 0 ) 
//       // 	std::cout << "WARNING: Weight vector is too small, size " << mcWeight.size() << std::endl ; 
//       // weight *= mcWeight[nPileup];
//       // std::cout << "MC weight is now " << weight << std::endl ; 
//     }

    // std::pair<float,double> pileupInfo = std::make_pair(avg_nvtx,weight) ; 
    std::pair<float,double> pileupInfo = std::make_pair(nPileup,weight) ; 
    return pileupInfo ; 
  }

  float jecTotalUncertainty(float jpt, float jeta,
			    JetCorrectionUncertainty *jecUnc,
			    int correctEra,
			    bool isBjet,
			    bool directionIsUp) {
    float offunc; // the "official" eta/pt dependent uncertainty

    jecUnc->setJetPt((float)jpt);
    jecUnc->setJetEta((float)jeta);
    offunc = jecUnc->getUncertainty(directionIsUp);
        // these are no longer necessary, commented out to ensure they never get applied without express intent
        //    float pileupCorrection = (2 * 0.75 * 0.8) / jpt;
        //    if      (correctEra == 1) pileupCorrection *= pileup2010A;
        //    else if (correctEra == 2) pileupCorrection *= pileup2010B;
        //    else if (correctEra == 3) pileupCorrection *= pileup2011A;
        //    else { // Merge 2010A and 2010B corrections
        //      float pileup2010 = ((3.18 * pileup2010A) + (32.96 * pileup2010B)) / 36.14;
        //      pileupCorrection *= pileup2010;
        //    }
    // Calculations taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC, version 23
    float totalUnc = offunc;//sqrt((offunc*offunc) + (pileupCorrection*pileupCorrection) + (0.025*0.025));
    return totalUnc;
  }

  std::vector< std::pair<pat::Jet,float> > getJetList(edm::Handle<pat::JetCollection>& pJets,
						      JetCorrectionUncertainty* jecUnc,
						      double minPt, double maxAbsEta, 
						      int jecSign, int jecEra, 
						      bool isMC, int jerSign) {
    
    std::vector< std::pair<pat::Jet,float> > jetList ; 

    for (unsigned int iJet=0; iJet<pJets->size(); iJet++) {
      pat::Jet iJ = pJets->at(iJet);
      float jpt = iJ.pt();
      float jeta = iJ.eta();
      if (fabs(jeta) > maxAbsEta) continue ; 
      // std::cout << "I have a jet with pt " << jpt 
      //  		<< " and eta " << fabs(iJ.eta()) << std::endl ; 
      int jpdgId = 0;
      if (iJ.genParton()) jpdgId = iJ.genParton()->pdgId();
      bool isBjet = (abs(jpdgId) == 5);
      float jecuscale = 1.0f;
      // Corrections for Jet Energy Resolution differences in data vs. MC
      // Updated for 2012 analysis
      if ( isMC ) { // Known jet energy resolution difference between data and MC
	std::vector<double> etaBin = {   0.5,   1.1,   1.7,   2.3,   5.0 } ; 
	std::vector<double> factor = { 1.052, 1.057, 1.096, 1.134, 1.288 } ; 
	std::vector<double> uncHi  = { 0.063, 0.057, 0.065, 0.094, 0.200 } ; 
	std::vector<double> uncLo  = { 0.062, 0.056, 0.064, 0.092, 0.199 } ; 

	unsigned int ibin = 0 ; 
	for (ibin=0; etaBin.at(ibin) < fabs(iJ.eta()); ibin++) ;
	double corrFactor = factor.at(ibin) ; 
	if      ( jerSign > 0 ) corrFactor += uncHi.at(ibin) ; 
	else if ( jerSign < 0 ) corrFactor -= uncLo.at(ibin) ; 
	// std::cout << "Correction factor (" << jerSign << "): " << corrFactor << std::endl ; 

	const reco::GenJet* iG = iJ.genJet() ; 
	if ( iG ) { 
	  double corr_delta_pt = ( iJ.pt() - iG->pt() ) * corrFactor ; 
	  // std::cout << "Corrected delta pt is: " << corr_delta_pt << " (" << iJ.pt() << "," << iG->pt() 
	  // 	    << "," << factor << ")" << std::endl ; 
	  double jerscale = std::max(0.0,((iG->pt()+corr_delta_pt)/iJ.pt())) ;
	  // std::cout << "jerscale is " << jerscale << std::endl ; 
	  jpt *= jerscale ; 
	  iJ.setP4( iJ.p4()*jerscale ) ; 
	}
      }
      if (jecSign) {
	float jecu = jecTotalUncertainty(jpt, jeta, jecUnc, jecEra, isBjet, (jecSign > 0));
	jecuscale = (1.0 + (float(jecSign) * jecu));
	// std::cout << "jecuscale is " << jecuscale << std::endl ; 
	jpt *= jecuscale;
	iJ.setP4(iJ.p4()*jecuscale);
      }
      // std::cout << "Jet cand with pt " << jpt << " eta " << jeta << " and electron energy fraction " 
      //  		<< iJ.electronEnergyFraction() << "%.  " << iJ.photonEnergy() << " + " 
      //  		<< iJ.electronEnergy() << " attributed to EM photons and electrons (in GeV)" << std::endl ; 
      // std::cout << "Corrected jet pT: " << jpt << std::endl ; 
      if (jpt > minPt) { 
	std::pair<pat::Jet,float> jetCand = std::make_pair(iJ,jecuscale) ;
	jetList.push_back( jetCand ) ; 
      }
    }

    std::sort(jetList.begin(),jetList.end(),scaleCompare()) ; 
    return jetList ; 
  }

  double muScaleLUT(pat::Muon& iM) { 

    const double etastep   = 2.0 * 2.4 / 5 ; 
    const double etavec[6] = {-2.4,(-2.4+1.0*etastep),(-2.4+2.0*etastep),
			      (-2.4+3.0*etastep),(-2.4+4.0*etastep),2.4} ; 
    unsigned int ieta = 0 ; 
    while (ieta < 5 && iM.eta() > etavec[ieta+1]) ieta++ ;  

    // Alignment corrections are given in chg/TeV.
    // To get GeV scale, need to divide by 10^3
    return muScaleLUTarray100[ieta] / 1000. ; 
  } 

  std::vector<pat::Muon> getMuonList(edm::Handle<pat::MuonCollection>& pMuons,
				     edm::Handle<reco::VertexCollection>& pvHandle, 
				     int idEra, double minPt, double maxAbsEta, 
				     double mesScale, bool muBiasUnc, 
				     bool trackerPt) {

    // std::cout << "Please note that the ID era is: " << idEra << std::endl ; 

    double ptScale = mesScale ; 
    std::vector<pat::Muon> muonList ; 
    for (unsigned int iMuon = 0; iMuon < pMuons->size(); iMuon++) {
      pat::Muon iM = pMuons->at(iMuon);
      
      if(iM.globalTrack().isNull() || iM.innerTrack().isNull()) continue;
      
      // For 53 we must recalculate the muon Pt (This needs CMSSW_5_3_6_p1)
      reco::TrackRef cktTrack = (muon::tevOptimized(iM, 200, 40., 17., 0.25)).first;
      if (cktTrack.isNull()) continue;
      reco::Particle::PolarLorentzVector p4(cktTrack->pt(),iM.eta(),iM.phi(),0.1057);
      
      if ( !iM.isGlobalMuon() ) continue ; 
      if ( fabs(iM.eta()) > maxAbsEta ) continue ;

      if ( muBiasUnc ) { 
        int charge = iM.charge() / abs(iM.charge()) ; 
        double k  = 0. ; 
        if ( fabs(iM.eta()) < 1.2 ) k = 0.05 / 1000. ; 
        double pt = (cktTrack.isNull() ? iM.pt() : cktTrack->pt() ) ; 
        ptScale   = ( double(charge) / ( double(charge) + k*pt ) ) ; 
      }

      // Method for getting high pT muons has changed in 53X
      if ( cktTrack.isNonnull() ) {
        iM.setP4( p4 * ptScale ) ; 
      } else { 
        iM.setP4( iM.p4() * ptScale ) ;
      }

      if ( iM.pt() < minPt ) continue ; 
      bool passesID = false ;

      // ID options
      if ( idEra == 0 ) passesID = true ;
      else if ( idEra > 0 ) passesID = isTightHighPtMuon(iM, pvHandle) ; // Default analysis selection
      else if ( idEra < 0 ) passesID = isLooseMuonNoPF(iM) ; // Global or Tracker muon (QCD or special studies only)
      
      if ( !passesID ) continue ;

      // Now take a look at TeV (refit) muons, and see if pT needs adjusting
      if ( trackerPt ) {
        double energy = sqrt( iM.innerTrack()->p()*iM.innerTrack()->p() + 0.1057 * 0.1057 ) ; // Track p^2 + muon_mass^2 
        reco::Particle::LorentzVector trackP4(iM.innerTrack()->px(),iM.innerTrack()->py(),iM.innerTrack()->pz(),energy) ;  
        iM.setP4( trackP4 ) ; 
      }

      muonList.push_back(iM);
    }

    std::sort(muonList.begin(),muonList.end(),pTcompare()) ; 
    return muonList ; 
  }

  double getElectronEscale(double eta1, double eta2) { 
    
    double scaleCorrection = 1.0 ; 
    int nEB = 0 ; 
    int nEE = 0 ; 

    if ( fabs(eta1) < 1.442 ) nEB++ ; 
    if ( fabs(eta2) < 1.442 ) nEB++ ; 
    if ( fabs(eta1) < 2.5 && fabs(eta1) > 1.56 ) nEE++ ; 
    if ( fabs(eta2) < 2.5 && fabs(eta2) > 1.56 ) nEE++ ; 
    
    // All corrections to MC, bringing energy down to match data
    // TAKEN FROM exo-2012_45_V10, PAGE 24 
    double ebebCorr = 1.0 - 0.0057; 
    double ebeeCorr = 1.0 - 0.0068; 
    double eeeeCorr = 1.0 - 0.0056; 

    if ( nEB + nEE < 2 ) { // Special case for top
      if      ( nEB ) scaleCorrection = sqrt( ebebCorr ) ; 
      else if ( nEE ) scaleCorrection = sqrt( eeeeCorr ) ; 
    } else { 
      if      ( nEB == 2 ) scaleCorrection = ebebCorr ; 
      else if ( nEE == 2 ) scaleCorrection = eeeeCorr ; 
      else                 scaleCorrection = ebeeCorr ; 
    }

    return scaleCorrection ; 
  }


  std::vector< std::pair<pat::Electron,float> > getElectronList(edm::Handle<pat::ElectronCollection>& pElecs,
								double maxAbsEta, 
								double minEtEB, double minEtEE, 
								int heepVersion,
                                edm::Handle<reco::VertexCollection> pvHandle,
								double rho,
								float ebScale, float eeScale) {
    
    std::vector< std::pair<pat::Electron,float> > electronList ; 
    // Just in case...
    if(heepVersion > 0)
    { // heepVersion = 0 for no heep selection
        if ( heepVersion != 31 && heepVersion != 32 && heepVersion != 40 && heepVersion != 41 && heepVersion != 0) {
          std::cout << "WARNING Invalid HEEP version: " << heepVersion << std::endl ;
          return electronList ;
        }
    }
    for (unsigned int iElectron=0; iElectron < pElecs->size(); iElectron++) {
      pat::Electron iE = pElecs->at(iElectron); 
      if ( getElectronSCEta(iE) > maxAbsEta ) continue ; 

      // Apply selection if requested
      if ( (heepVersion > 0) && !passesHEEP(iE,heepVersion,rho, pvHandle)) continue ;
      // Apply anti-selection if requested (fake electron sample)
      if ( (heepVersion < 0) && 
	   (!passesFakeRateMinimum(iE,pvHandle) || passesHEEP(iE,abs(heepVersion),rho, pvHandle)) ) continue ;

      float scale      = ( (iE.isEB()) ? ebScale : eeScale ) ;
      float elecEt     = getElectronEt(iE,false) * scale ; 
      bool  passEtCuts = ( (iE.isEB()) ? (elecEt > minEtEB) : (elecEt > minEtEE) ) ;

      if (passEtCuts) {
        iE.setP4(iE.p4() * scale);
        std::pair<pat::Electron,float> electronCand = std::make_pair(iE,scale) ;
        electronList.push_back( electronCand ) ; 
      }
    }
    // std::sort(electronList.begin(),electronList.end(),scaleCompare()) ; 
    std::sort(electronList.begin(),electronList.end(),etCompare()) ; 
    return electronList ; 
    }

    bool addMuon(pat::Muon& m, HeavyNuEvent& hne, double minimum_muon_jet_dR, int nDirtyCands)
    {
        double dRj1 = 999.9, dRj2 = 999.9;
        if(hne.nJets > 0) dRj1 = deltaR(m.eta(), m.phi(), hne.j1.eta(), hne.j1.phi());
        if(hne.nJets > 1) dRj2 = deltaR(m.eta(), m.phi(), hne.j2.eta(), hne.j2.phi());
        double dRm1 = 999. ;
        if ( hne.nLeptons > 0 && nDirtyCands > 0) dRm1 = deltaR(m.eta(), m.phi(), hne.mu1.eta(), hne.mu1.phi());
        if (dRj1 > minimum_muon_jet_dR && dRj2 > minimum_muon_jet_dR && dRm1 > 0.3)
        {
        	hne.nMuons++ ;
            hne.nLeptons++ ;
            if(hne.nMuons == 1)      hne.mu1 = m;
            else if(hne.nMuons == 2) hne.mu2 = m;
            else std::cout << "WARNING: Expected empty muon position" << std::endl;
            return true;
        }
        return false;
    }

    bool addElectron(pat::Electron& e, HeavyNuEvent& hne, double minimum_muon_jet_dR)
    {
        double dRj1 = 999.9, dRj2 = 999.9;
        if(hne.nJets > 0) dRj1 = deltaR(e.eta(), e.phi(), hne.j1.eta(), hne.j1.phi());
        if(hne.nJets > 1) dRj2 = deltaR(e.eta(), e.phi(), hne.j2.eta(), hne.j2.phi());
        if(dRj1 > minimum_muon_jet_dR && dRj2 > minimum_muon_jet_dR)
        {
            hne.nElectrons++ ;
            hne.nLeptons++;
            if(hne.nElectrons == 1) hne.e1 = e;
            else if(hne.nElectrons == 2) 
            {
                if(getElectronEt(hne.e1, false) > getElectronEt(e, false)) hne.e2 = e;
                else
                {
                    hne.e2 = hne.e1;
                    hne.e1 = e;
                }
            }
            else std::cout << "WARNING: Expected empty electron position" << std::endl;
            return true;
        }
        return false;
    }


    const double ELECTRON_MASS = 0.000511;
    const double MUON_MASS     = 0.105;
    const double TAU_MASS      = 1.7768;
    bool isMuon = false;

    TLorentzVector decayTau(const reco::Candidate* tau, bool isMuon)
    {
        const double LEP_MASS = (isMuon)?MUON_MASS:ELECTRON_MASS;
        const double E_MAX    = (TAU_MASS*TAU_MASS + LEP_MASS*LEP_MASS) / (2 * TAU_MASS);
        
        TF1 *energy = new TF1("e_dist", tauDecay, 0.0, E_MAX, 0);

        hnu::isMuon = isMuon;
        TLorentzVector *v = new TLorentzVector(0.0, 0.0, 0.0, 0.0);
        v->SetPtEtaPhiM(energy->GetRandom(), 0.0, 0.0, LEP_MASS);
        
        double rand_theta = acos(gRandom->Uniform(-1, 1)), rand_phi = gRandom->Uniform(0, 2*3.14159);

        //v->RotateX(cos(rand_phi) * sin(rand_theta));
        v->RotateY(rand_theta);
        v->RotateX(rand_phi);
        
        double tauE = sqrt(tau->p()*tau->p()+TAU_MASS*TAU_MASS);
        v->Boost(tau->px()/tauE, tau->py()/tauE, tau->pz()/tauE);

        return *v;
    }

    double tauDecay(double *x, double *par)
    {
        const double LEP_MASS = (isMuon)?MUON_MASS:ELECTRON_MASS;
        const double E_MAX    = (TAU_MASS*TAU_MASS + LEP_MASS*LEP_MASS) / (2 * TAU_MASS);
        return (2/E_MAX)*(3*x[0]*x[0]/(E_MAX*E_MAX) - 2*x[0]*x[0]*x[0]/(E_MAX*E_MAX*E_MAX));
    }
    
    double getBtagSF(double pt)
    {
        // CSVM scale factor 
        // https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_payload_Moriond13.txt
        return 0.726981*((1.+(0.253238*pt))/(1.+(0.188389*pt)));
    }

}





#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

void trigtools::getP4sOfObsPassingFilter(std::vector<math::XYZTLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess)
{
  p4s.clear();

  edm::InputTag filterTag(filterName,"",hltProcess); 
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){ //check that filter is in triggerEvent
    const trigger::Keys& trigKeys = trigEvent.filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      math::XYZTLorentzVector objP4;
      objP4.SetPxPyPzE(obj.px(),obj.py(),obj.pz(),obj.energy());
      p4s.push_back(objP4);
    }//end loop over keys
  }//end check that filter is valid and in trigEvent
}

void trigtools::getP4sOfObsPassingFilter(std::vector<TLorentzVector>& p4s,const trigger::TriggerEvent& trigEvent,const std::string& filterName,const std::string& hltProcess)
{
  p4s.clear();
 
  edm::InputTag filterTag(filterName,"",hltProcess); 
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){ //check that filter is in triggerEvent
    const trigger::Keys& trigKeys = trigEvent.filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      TLorentzVector objP4;
      objP4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
      p4s.push_back(objP4);
    }//end loop over keys
  }//end check that filter is valid and in trigEvent
}



