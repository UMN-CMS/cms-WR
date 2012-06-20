#include "HeavyNuCommon.h"
#include "TVector3.h"

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
  bool isVBTFloose(const pat::Muon& m)
  {
    return m.muonID("AllGlobalMuons")&&(m.numberOfValidHits()>10);
  }
  
  bool isVBTFtight(const pat::Muon& m)
  {
    if( !isVBTFloose(m) ) return false; // this should already have been checked.

    reco::TrackRef gt = m.globalTrack();
    if (gt.isNull()) {
      std::cerr << "Mu global track reference is NULL" << std::endl;
      return false;
    }
    return (m.muonID("AllTrackerMuons") &&
	    (m.dB() < 0.2) &&
	    (m.normChi2() < 10) &&
	    (m.numberOfMatches() > 1) &&
	    (gt->hitPattern().numberOfValidMuonHits()>0) &&
	    (gt->hitPattern().numberOfValidPixelHits()>0) );

  } // HeavyNu::isVBTFtight

  // Information taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  // On 12 May, 2012
  bool is2012MuTight(const pat::Muon& m, edm::Handle<reco::VertexCollection> pvHandle) { 

    reco::TrackRef gt = m.globalTrack();
    reco::TrackRef it = m.innerTrack();
    if (gt.isNull() || it.isNull()) {
      //std::cerr << "Mu global track or inner track reference is NULL" << std::endl;
      return false;
    }
    if ( !pvHandle.isValid() ) { 
      //std::cerr << "No primary vertex available.  Skipping advanced selection." << std::endl;
      return false;
    }      
    reco::Vertex pv = pvHandle->at(0) ; 

    return ( (m.isGlobalMuon()) &&
	     // (m.normChi2() < 10) && --> Ignored for high pT muons
	     (gt->hitPattern().numberOfValidMuonHits() > 0) &&
	     (it->hitPattern().numberOfValidPixelHits() > 0) &&
	     (m.numberOfMatchedStations() > 1) &&
	     (m.dB() < 0.2) &&
	     (fabs(it->dz(pv.position())) < 0.5) &&
	     // (it->hitPattern().trackerLayersWithMeasurement() > 5) ) ; 
	     (it->hitPattern().trackerLayersWithMeasurement() > 8) ) ; // modification for high pT muons
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

  bool passesFakeRateMinimum(const pat::Electron& e) { 

    if ( !e.ecalDriven() ) return false ; 

    double eEta = getElectronSCEta(e) ; 
    bool isEB   = ( fabs(eEta) < 1.442 ) ; 
    bool isEE   = ( fabs(eEta) < 2.5 && fabs(eEta) > 1.56 ) ; 
    if ( !isEB && !isEE ) return false ; 

    double HoE      = e.hadronicOverEm() ; 
    double sig_iEiE = e.sigmaIetaIeta() ; 
    int nLostHits   = e.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ; 

    if ( nLostHits != 0 ) return false ; 
    if ( isEB ) { 
      if ( HoE > 0.15 )       return false ; 
      if ( sig_iEiE > 0.013 ) return false ; 
    } else if ( isEE ) { 
      if ( HoE > 0.10 )       return false ; 
      if ( sig_iEiE > 0.034 ) return false ; 
    }

    return true ; 
  } 

  double fakeProbability(const pat::Electron& e) { 

    double prob = 0. ; 

    double eEta = getElectronSCEta(e) ; 
    double ePt  = getElectronEt(e,false) ; 

    if ( fabs(eEta) < 1.442 ) { 
      prob =  (( ePt < 94.2 ) ?  ( 0.021 - 0.00012 * ePt ) : ( 0.0097 )) ; 
    } else if ( fabs(eEta) < 2.0 && fabs(eEta) > 1.56 ) { 
      prob = (( ePt < 100.3 ) ? ( 0.0711 - 0.00035 * ePt ) : ( 0.036 )) ; 
    } else if ( fabs(eEta) < 2.5 && fabs(eEta) > 2.0 ) { 
      prob = (( ePt < 109.5 ) ?  ( 0.118 - 0.00047 * ePt ) : ( 0.062 + 0.000041 * ePt )) ; 
    }

    return prob ; 
  }

  // Placeholder until real values are known
  double fakeProbability(const pat::Muon& mu) { 

    double eta  = mu.eta() ; 
    double prob = 0. ; 

    if ( abs(eta) < 2.4 ) { // Values determined from data on 20 June 2012
      if      ( eta < -2.1 ) prob = 0.25938 ; 
      else if ( eta < -1.2 ) prob = 0.13414 ; 
      else if ( eta < -0.9 ) prob = 0.105791 ; 
      else if ( eta < 0.0 )  prob = 0.0944799 ; 
      else if ( eta < 0.9 )  prob = 0.095889 ; 
      else if ( eta < 1.2 )  prob = 0.108 ; 
      else if ( eta < 2.1 )  prob = 0.135556 ; 
      else if ( eta < 2.4 )  prob = 0.25242 ; 
    }

    return prob ; 
  }

  bool passesHEEP(const pat::Electron& e, int heepVersion, double rho) { 

    if ( heepVersion != 31 && heepVersion != 32 && abs(heepVersion) != 40 ) return false ; 

    double ePt  = getElectronEt(e,(abs(heepVersion) != 40)) ; 
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
    if (  heepVersion == 31                            && dPhiIn > 0.09 ) return false ; 
    if ( (heepVersion == 32 || abs(heepVersion) == 40) && dPhiIn > 0.06 ) return false ; 

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
    if ( nLostHits != 0 ) return false ; 
    
    double ecalHcalIso = e.dr03EcalRecHitSumEt() + e.dr03HcalDepth1TowerSumEt() ; 
    double thresholdEB = 2. + 0.03 * ePt ; 
    double thresholdEE = 2.5 ;
    if ( isEE && ePt > 50 ) thresholdEE += (0.03 * (ePt - 50.)) ; 
    if ( abs(heepVersion) == 40 ) { 
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
    if ( (heepVersion == 31) &&
	 ((isEB && trkIso > 7.5) || (isEE && trkIso > 15.)) ) return false ; 
    if ( (heepVersion == 32 || abs(heepVersion) == 40) && (trkIso > 5.) ) return false ; 

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

  std::vector<float> generate_flat10_mc(int era){
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    // const double npu_probs[25] = {0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 0-4
    // 				  0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 5-9
    // 				  0.0698146584,0.0630151648,0.0526654164,0.0402754482,0.0292988928,  // 10-14
    // 				  0.0194384503,0.0122016783,0.0072070420,0.0040036370,0.0020278322,  // 15-19
    // 				  0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 // 20-24 
    // };

    // see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities for PU_S4 samples
    // Using the "spike at zero + smearing distribution" as shown on the twiki and recommended for in-time PU corrections
    const double npu_probs_summer11[35] = { 1.45346E-01,6.42802E-02,6.95255E-02,6.96747E-02,6.92955E-02, // 0-4
					    6.84997E-02,6.69528E-02,6.45515E-02,6.09865E-02,5.63323E-02, // 5-9
					    5.07322E-02,4.44681E-02,3.79205E-02,3.15131E-02,2.54220E-02, // 10-14
					    2.00184E-02,1.53776E-02,1.15387E-02,8.47608E-03,6.08715E-03, // 15-19
					    4.28255E-03,2.97185E-03,2.01918E-03,1.34490E-03,8.81587E-04, // 20-24
					    5.69954E-04,3.61493E-04,2.28692E-04,1.40791E-04,8.44606E-05, // 25-29
					    5.10204E-05,3.07802E-05,1.81401E-05,1.00201E-05,5.80004E-06  // 30-34
    }; 

    const double npu_probs_fall11[50] = { 0.003388501, 0.010357558, 0.024724258, 0.042348605, 0.058279812, // 0-4
					  0.068851751, 0.072914824, 0.071579609, 0.066811668, 0.060672356, // 5-9
					  0.054528356,  0.04919354, 0.044886042, 0.041341896,   0.0384679, // 10-14
					  0.035871463,  0.03341952, 0.030915649, 0.028395374, 0.025798107, // 15-19
					  0.023237445, 0.020602754,   0.0180688, 0.015559693, 0.013211063, // 20-24
					  0.010964293, 0.008920993, 0.007080504, 0.005499239, 0.004187022, // 25-29
					  0.003096474, 0.002237361, 0.001566428, 0.001074149, 0.000721755, // 30-34
					  0.000470838,  0.00030268, 0.000184665, 0.000112883, 6.74043E-05, // 35-39
					  3.82178E-05, 2.22847E-05, 1.20933E-05, 6.96173E-06,  3.4689E-06, // 40-44
					  1.96172E-06, 8.49283E-07, 5.02393E-07, 2.15311E-07, 9.56938E-08  // 45-49
    };

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



    bool isSummer11 = ( era == 20111 || era == 20112 ) ; 
    bool isFall11   = ( era == 20113 || era == 20114 ) ; 
    const double* npu_probs = npu_probs_summer12 ;
    int npt = sizeof(npu_probs_summer12)/sizeof(double);
    // const double* npu_probs = npu_probs_summer11 ;
    // int npt = sizeof(npu_probs_summer11)/sizeof(double);
    if ( isSummer11 ) {
        npt = sizeof(npu_probs_summer11)/sizeof(double);
        npu_probs = npu_probs_summer11 ;
    }
    if ( isFall11 ) {
        npt = sizeof(npu_probs_fall11)/sizeof(double);
        npu_probs = npu_probs_fall11 ;
    }

    std::vector<float> retval;
    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(npu_probs[i]);
    return retval;
    // double s = 0.0;
    // for(int npu=0; npu<25; ++npu){
    //   double npu_estimated = dataDist[npu];
    //   result[npu] = npu_estimated / npu_probs[npu];
    //   s += npu_estimated;
    // }
    // // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    // for(int npu=0; npu<25; ++npu){
    //   result[npu] /= s;
    // }
    // return result;
  }


  std::vector<float> get_standard_pileup_data(int era) {
    const double default_pd[] = { 100, 100, 100, 0, 0, 0, 0, 0, 0, 0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
    const double may10_json[]= {
      3920760.81, 6081805.28, 13810357.99, 22505758.94, 28864043.84, 30917427.86, 28721324.56, 23746403.90, 17803439.77, 12274902.61, 
      7868110.47, 4729915.40, 2686011.14, 1449831.56, 747892.03, 370496.38, 177039.19, 81929.35, 36852.78, 16164.45, 
      6932.97, 2914.39, 1202.92, 488.15, 194.93, 76.64, 29.66, 11.30, 4.24, 1.56, 
      0.57, 0.20, 0.07, 0.02, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00 };
    const double dec22_json[]= {
      5525914.08, 9064207.39, 8760233.30, 6230731.01, 3615158.71, 1806273.52, 802707.46, 323868.48, 120245.31, 41458.73, 
      13360.27, 4043.62, 1153.89, 311.48, 79.76, 19.43, 4.51, 1.00, 0.21, 0.04, 
      0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
      0.00 };

    // Pileup histograms assembled from inputs in this directory: 
    // /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp

    // Note: there is a 36th entry with 0.008, but other ingredients have only 35 entries so omitting
    const double json_2011a[] = {
      1.29654E07, 5.58514E07, 1.29329E08, 2.12134E08, 2.76138E08,
      3.03604E08, 2.93258E08, 2.55633E08,  2.0497E08, 1.53264E08,
      1.07936E08, 7.21006E07,  4.5913E07,   2.797E07, 1.63426E07,
      9.17598E06, 4.95861E06, 2.58239E06,  1.2977E06,     629975,
          295784,     134470,    59260.1,    25343.9,    10530.1,
         4255.05,    1673.95,    641.776,    240.022,    87.6504,
          31.281,    10.9195,    3.73146,    1.24923,   0.602368
    } ; 

    const double json_2011b[] = {
          481142, 3.21393E06, 1.15733E07, 2.91676E07, 5.76072E07,
      9.51074E07, 1.36849E08,  1.7665E08,  2.0885E08, 2.29582E08,
      2.37228E08, 2.32243E08, 2.16642E08, 1.93361E08,  1.6564E08,
      1.36514E08, 1.08455E08, 8.31965E07, 6.17147E07, 4.43296E07,
      3.08733E07, 2.08734E07, 1.37166E07, 8.77106E06, 5.46389E06,
      3.31952E06, 1.96896E06,  1.1414E06,     647299,     359460,
          195642,     104449,    54741.4,    28184.3,    28004.9
    } ; 

    // Pileup histograms for Fall11 assembled from inputs in this directory: 
    // /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp
    // Per instructions, the so-called "true" distributions from data are used
    const double json_2011a_44x[] = { 
             0.0,   252573.0, 5.73861E06, 5.05641E07, 2.68197E08,
      5.12977E08, 5.21466E08, 3.75661E08, 2.41467E08,  1.4315E08,
      5.38318E07, 1.18125E07, 1.29073E06,     111569,    6537.93,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0
    } ; 

    const double json_2011b_44x[] = { 
             0.0,    27267.4,      35590,    74493.3,     574590,
      2.90648E06, 3.36311E07, 9.36661E07, 1.38283E08, 1.87624E08,
      2.15647E08,  2.1173E08, 1.87002E08, 1.46693E08, 9.44372E07,
      4.60317E07, 1.69231E07, 5.18161E06, 1.42805E06,     437008,
          102694,     6516.2,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0,
	     0.0,        0.0,        0.0,        0.0,        0.0
    } ; 

    // const double json_2012a_run193557[] = { 
    //        59.3073,      156.1,    6681.18,  2.79236E06, 7.28453E06,
    // 	   66597.5,    60190.1,    48270.5,      155944, 1.19828E06, 
    // 	7.83965E06, 2.16456E07, 3.29565E07,  4.09794E07, 4.95412E07, 
    // 	5.83303E07, 6.37098E07, 6.28154E07,  5.53595E07, 4.38815E07, 
    // 	3.23858E07, 2.31768E07, 1.65636E07,  1.22208E07, 9.41924E06, 
    // 	7.42174E06, 5.82024E06,  4.4711E06,   3.3413E06, 2.42264E06, 
    // 	1.70269E06, 1.15991E06,     766260,      491354,     306205, 
    // 	    185700,     109734,    63248.8,     35581.5,    19541.5, 
    // 	   10475.4,    5477.99,    2792.36,     1386.16,    669.462, 
    // 	   314.265,     143.27,    63.3834,     27.1948,    11.3099, 
    // 	   4.55733,    1.77868,   0.672206,    0.245942,  0.0870985, 
    // 	 0.0298523, 0.00990102, 0.00317743, 0.000986567,     402199
    // } ; 

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

    const double* pileupDist=default_pd;
    int npt = 35;
    if (era == 20121) 
    {
        // npt = sizeof(json_2012a_run193557)/sizeof(double);
        // pileupDist=json_2012a_run193557;
        npt = sizeof(json_2012ab_run195868)/sizeof(double);
        pileupDist=json_2012ab_run195868;
    }
    if (era == 20111) 
    {
        npt = sizeof(json_2011a)/sizeof(double);
        pileupDist=json_2011a;
    }
    if (era == 20112) 
    {
        npt = sizeof(json_2011b)/sizeof(double);
        pileupDist=json_2011b;
    }
    if (era == 20113) 
    {
        npt = sizeof(json_2011a_44x)/sizeof(double);
        pileupDist=json_2011a_44x;
    }
    if (era == 20114) 
    {
        npt = sizeof(json_2011b_44x)/sizeof(double);
        pileupDist=json_2011b_44x;
    }
    if (era == 20110) 
    {
        npt = sizeof(may10_json)/sizeof(double);
        pileupDist=may10_json;
    }
    if (era>=20100 && era<=20109) 
    {
        npt = sizeof(dec22_json)/sizeof(double);
        pileupDist=dec22_json;
    }

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
      // 		<< " and eta " << fabs(iJ.eta()) << std::endl ; 
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
      // 		<< iJ.electronEnergyFraction() << "%.  " << iJ.photonEnergy() << " + " 
      // 		<< iJ.electronEnergy() << " attributed to EM photons and electrons (in GeV)" << std::endl ; 
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
				     edm::Handle<reco::MuonCollection>& tevMuons,
				     edm::Handle<reco::VertexCollection>& pvHandle, 
				     int idEra, double minPt, double maxAbsEta, 
				     double mesScale, bool muBiasUnc, 
				     bool trackerPt) {

    // std::cout << "Please note that the ID era is: " << idEra << std::endl ; 

    double ptScale = mesScale ; 
    std::vector<pat::Muon> muonList ; 
    for (unsigned int iMuon = 0; iMuon < pMuons->size(); iMuon++) {
      pat::Muon iM = pMuons->at(iMuon) ; 
      if ( !iM.isGlobalMuon() ) continue ; 
      // std::cout << "Investigating muon with pt " << iM.pt() << std::endl ; 
      if ( fabs(iM.eta()) > maxAbsEta ) continue ; 
      // std::cout << "Checking tevOptimized" << std::endl ; 
      // std::cout << iM.isGlobalMuon() << std::endl ; 
      // std::cout << iM.globalTrack().isNull() << std::endl ; 
      // std::cout << iM.innerTrack().isNull() << std::endl ; 
      // std::cout << iM.tpfmsTrack().isNull() << std::endl ; 
      // std::cout << iM.pickyTrack().isNull() << std::endl ; 
      // reco::Muon::MuonTrackTypePair muPair = muon::tevOptimized(iM,200.,4.,6.) ; 
      // std::cout << "Got here" << std::endl ; 

      reco::TrackRef muRef = muon::tevOptimized(iM,200.,4.,6.).first ; 
      // reco::TrackRef muRef = muPair.first ; 
      // std::cout << "Is muRef null? " << muRef.isNull() << std::endl ; 

      if ( muBiasUnc ) { 
	int charge = iM.charge() / abs(iM.charge()) ; 
	// std::cout << "Charge is " << charge << std::endl ; 
	// Modification of muon scale uncertainty: 5% per TeV, only for barrel
	// double k   = muScaleLUT(iM) ;
	double k  = 0. ; 
	if ( fabs(iM.eta()) < 1.2 ) k = 0.05 / 1000. ; 
	double pt = (muRef.isNull() ? iM.pt() : muRef->pt() ) ; 
	ptScale   = ( double(charge) / ( double(charge) + k*pt ) ) ; 
      }

      // Method for getting high pT muons has changed in 50X - 52X
      if ( muRef.isNonnull() ) { 
	reco::Particle::PolarLorentzVector muP4( muRef->pt(),muRef->eta(),muRef->phi(),0.1057) ; 
	iM.setP4( muP4 * ptScale ) ; 
      } else { 
	iM.setP4( iM.p4() * ptScale ) ; 
      }

      double mupt = iM.pt() ;
      if ( mupt < minPt ) continue ; 
      // std::cout << "Checking ID " << iM.pt() << std::endl ; 
      bool passesID = ( (idEra == 2012) ? ( is2012MuTight(iM,pvHandle) ) : ( isVBTFtight(iM) ) ) ; 
      // std::cout << "Muon ID: " << passesID << std::endl ; 
      if ( !passesID ) continue ; 
      // if ( !isVBTFtight(iM) ) continue ; 

      // Now take a look at TeV (refit) muons, and see if pT needs adjusting
      if ( trackerPt ) { 
	double energy = sqrt( iM.innerTrack()->p()*iM.innerTrack()->p() + 
			      0.1057 * 0.1057 ) ; // Track p^2 + muon_mass^2 
	reco::Particle::LorentzVector trackP4(iM.innerTrack()->px(),iM.innerTrack()->py(),iM.innerTrack()->pz(),energy) ;  
	iM.setP4( trackP4 ) ; 
      // } else { 
      // 	// Method for getting high pT muons has changed in 50X - 52X
      // 	iM.setP4( (muon::tevOptimized(iM,200.,4.,6.).first->p4()) * ptScale ) ; 

      // 	// for (unsigned int iTeV=0; iTeV<tevMuons->size() ; iTeV++) { 
      // 	//   reco::Muon tevMuon = tevMuons->at(iTeV) ; 
      // 	//   double tevMuPt = tevMuon.pt() * ptScale ; 
      // 	//   if ( tevMuPt < minPt ) continue ;  
      // 	//   double dR = ROOT::Math::VectorUtil::DeltaR(tevMuon.p4(),iM.p4()) ;
      // 	//   if ( dR < 0.001 ) { // Replace muon p4 with refit (TeV) muon p4
      // 	//     iM.setP4( tevMuon.p4() * ptScale ) ;
      // 	//   }
      // 	// }
      }

      muonList.push_back(iM) ; 
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
    // Updated to AN 12/171, v8 results
    double ebebCorr = 1.0 - 0.0088 ; 
    double ebeeCorr = 1.0 - 0.0143 ; 
    double eeeeCorr = 1.0 - 0.0145 ; 

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
								double rho,
								float ebScale, float eeScale) {
    
    std::vector< std::pair<pat::Electron,float> > electronList ; 
    // Just in case...
    if(heepVersion > 0)
    {
        if ( heepVersion < 30 || heepVersion > 40 ) {
          if ( heepVersion <= 2 ) heepVersion += 30 ;
          else                    heepVersion = 40 ;
        }   // the 0 case is meant to provide: "do no HEEP selection"
        if ( heepVersion != 31 && heepVersion != 32 && heepVersion != 40 && heepVersion != 0) {
          std::cout << "WARNING Invalid HEEP version: " << heepVersion << std::endl ;
          return electronList ;
        }
    }
    for (unsigned int iElectron=0; iElectron < pElecs->size(); iElectron++) {
      pat::Electron iE = pElecs->at(iElectron); 
      // std::cout << "Considering electron with pT: " << getElectronEt(iE,(heepVersion != 40)) 
      // 		<< " and eta: " << getElectronSCEta(iE) << std::endl ; 
      // if ( getElectronSCEta(iE) > maxAbsEta ) std::cout << "Event outside boundary, ignoring" << std::endl ; 
      if ( getElectronSCEta(iE) > maxAbsEta ) continue ; 
      // if ( (heepVersion == 31) && !passesHEEPv31(iE) ) continue ; 
      // if ( (heepVersion == 32) && !passesHEEPv32(iE) ) continue ; 
      // if ( (heepVersion == 40) && !passesHEEPv40(iE) ) continue ; 

      // Apply selection if requested
      if ( (heepVersion > 0) && !passesHEEP(iE,heepVersion,rho) ) continue ;
      // Apply anti-selection if requested (fake electron sample)
      if ( (heepVersion < 0) && 
	   (!passesFakeRateMinimum(iE) || passesHEEP(iE,abs(heepVersion),rho)) ) continue ; 

      float scale      = ( (iE.isEB()) ? ebScale : eeScale ) ;
      float elecEt     = getElectronEt(iE,(heepVersion != 40)) * scale ; 
      bool  passEtCuts = ( (iE.isEB()) ? (elecEt > minEtEB) : (elecEt > minEtEE) ) ;

      if (passEtCuts) {
        iE.setP4(iE.p4() * scale);
	std::pair<pat::Electron,float> electronCand = std::make_pair(iE,scale) ;
	electronList.push_back( electronCand ) ; 
      }
    }
    std::sort(electronList.begin(),electronList.end(),scaleCompare()) ; 
    return electronList ; 
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


