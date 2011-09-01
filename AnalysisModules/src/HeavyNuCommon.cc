#include "HeavyNuCommon.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

const float pileup2010A = 1.2;
const float pileup2010B = 2.2;
const float pileup2011A = 5.0;

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

  double getElectronEt(const pat::Electron& e) { 
    reco::Particle::LorentzVector eScaled = e.p4() * (e.caloEnergy() / e.energy()) ;
    return eScaled.Et() ; 
  } 

  double getElectronSCEta(const pat::Electron& e) { 
    return e.superCluster()->eta() ; 
  }

  bool passesHEEPv31(const pat::Electron& e) { 
    if ( !e.ecalDriven() ) return false ; 
    double ePt  = getElectronEt(e) ; 
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


  double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR){

        // Find best match jptjets
        reco::JPTJetCollection::const_iterator tJet=jptJets.end();
        for (reco::JPTJetCollection::const_iterator i=jptJets.begin(); i != jptJets.end(); i++){
            if (tJet == jptJets.end() || ROOT::Math::VectorUtil::DeltaR(pJet.p4(),i->p4()) < ROOT::Math::VectorUtil::DeltaR(pJet.p4(),tJet->p4())){
                tJet =i;
            }
        }

        // Return Vert
        return avgVertex(*tJet, maxDeltaVR);
  }

  std::vector<float> generate_flat10_mc(const int npt){
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    // const double npu_probs[25] = {0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 0-4
    // 				  0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,  // 5-9
    // 				  0.0698146584,0.0630151648,0.0526654164,0.0402754482,0.0292988928,  // 10-14
    // 				  0.0194384503,0.0122016783,0.0072070420,0.0040036370,0.0020278322,  // 15-19
    // 				  0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 // 20-24 
    // };

    // see https://twiki.cern.ch/twiki/bin/view/CMS/PileupMCReweightingUtilities for PU_S4 samples
    // Additional values (25+ vertices) calculated by hand assuming average of 12 pileup collisions
    // Result is not exact but will get rid of problems with n(PU) > 25 in Monte Carlo
    const double npu_probs[51] = { 0.104109000,0.070357300,0.069844500,0.069825400,0.069705400, // 0-4
				   0.069790700,0.069675100,0.069448600,0.068033200,0.065104400, // 5-9
				   0.059803600,0.052739500,0.043951300,0.035220200,0.026671400, // 10-14
				   0.019411000,0.013397400,0.008985360,0.005751600,0.003514930, // 15-19
				   0.002120870,0.001228910,0.000705920,0.000384744,0.000219377, // 20-24
				   0.000105301,4.86004e-05,2.16002e-05,9.25723e-06,3.83058e-06, // 25-29
				   1.53223e-06,5.93122e-07,2.22421e-07,8.08802e-08,2.85460e-08, // 30-34
				   9.78719e-09,3.26240e-09,1.05807e-09,3.34129e-10,1.02809e-10, // 35-39
				   3.08426e-11,9.02712e-12,2.57918e-12,7.19770e-13,1.96301e-13, // 40-44
				   5.23469e-14,1.36557e-14,3.48657e-15,8.71641e-16,2.13463e-16, // 45-49
				   5.12312e-17 }; 

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


  std::vector<float> get_standard_pileup_data(int era, const int npt) {
    const double default_pd[] = { 100, 100, 100, 0, 0, 0, 0, 0, 0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  };
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

    const double* pileupDist=default_pd;

    if (era>=20110 && era<=20112) pileupDist=may10_json;
    if (era>=20100 && era<=20109) pileupDist=dec22_json;

    std::vector<float> retval;
    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(pileupDist[i]);
    return retval;
  }

  std::pair<float,double> pileupReweighting(edm::Handle<std::vector<PileupSummaryInfo> > pPU,
					    edm::LumiReWeighting puWeight) { 

    // int   nPileup = -1 ; 
    double weight   = 1.0 ; 
    float  avg_nvtx = -1. ; 

    if (pPU.isValid() && pPU->size() > 0) {
      std::vector<PileupSummaryInfo>::const_iterator puIter ; 

      float sum_nvtx = 0 ; 
      for (puIter=pPU->begin(); puIter!=pPU->end(); puIter++) 
	sum_nvtx += float( puIter->getPU_NumInteractions() ) ; 

      avg_nvtx = sum_nvtx / 3. ; 

      // std::cout << "About to look up weight for " << avg_nvtx << " vertices, initial weight " << weight << std::endl ; 
      weight = puWeight.weight3BX( avg_nvtx ) ; 
      // if ( nPileup > int(mcWeight.size()) || nPileup < 0 ) 
      // 	std::cout << "WARNING: Weight vector is too small, size " << mcWeight.size() << std::endl ; 
      // weight *= mcWeight[nPileup];
      // std::cout << "MC weight is now " << weight << std::endl ; 
    }

    std::pair<float,double> pileupInfo = std::make_pair(avg_nvtx,weight) ; 
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

    float pileupCorrection = (2 * 0.75 * 0.8) / jpt;
    if      (correctEra == 1) pileupCorrection *= pileup2010A;
    else if (correctEra == 2) pileupCorrection *= pileup2010B;
    else if (correctEra == 3) pileupCorrection *= pileup2011A;
    else { // Merge 2010A and 2010B corrections
      float pileup2010 = ((3.18 * pileup2010A) + (32.96 * pileup2010B)) / 36.14;
      pileupCorrection *= pileup2010;
    }
    // Calculations taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC, version 23
    float totalUnc = sqrt((offunc*offunc) + (pileupCorrection*pileupCorrection) + (0.025*0.025));
    return totalUnc;
  }

  std::vector< std::pair<pat::Jet,float> > getJetList(edm::Handle<pat::JetCollection>& pJets,
						      JetCorrectionUncertainty* jecUnc,
						      double minPt, double maxAbsEta, 
						      int jecSign, int jecEra) {
    
    std::vector< std::pair<pat::Jet,float> > jetList ; 

    for (unsigned int iJet=0; iJet<pJets->size(); iJet++) {
      pat::JetRef iJ = pat::JetRef(pJets, iJet);
      float jpt = (*iJ).pt();
      float jeta = (*iJ).eta();
      if (fabs(jeta) > maxAbsEta) continue ; 
      int jpdgId = 0;
      if (iJ->genParton()) jpdgId = iJ->genParton()->pdgId();
      bool isBjet = (abs(jpdgId) == 5);
      float jecuscale = 1.0f;
      if (jecSign) {
	float jecu = jecTotalUncertainty(jpt, jeta, jecUnc, jecEra, isBjet, (jecSign > 0));
	jecuscale = (1.0 + (float(jecSign) * jecu));
	jpt *= jecuscale;
      }
      if (jpt > minPt) { 
	std::pair<pat::Jet,float> jetCand = std::make_pair((*iJ),jecuscale) ;
	jetList.push_back( jetCand ) ; 
      }
    }

    std::sort(jetList.begin(),jetList.end(),scaleCompare()) ; 
    return jetList ; 
  }

  std::vector<pat::Muon> getMuonList(edm::Handle<pat::MuonCollection>& pMuons,
				     edm::Handle<reco::MuonCollection>& tevMuons,
				     double minPt, double maxAbsEta, 
				     double ptScale) {

    std::vector<pat::Muon> muonList ; 
    for (unsigned int iMuon = 0; iMuon < pMuons->size(); iMuon++) {
      pat::Muon iM = pMuons->at(iMuon) ; 
      double mupt = ptScale * (iM.pt());
      if ( mupt < minPt ) continue ; 
      if ( fabs(iM.eta()) > maxAbsEta ) continue ; 
      if ( !isVBTFtight(iM) ) continue ; 

      // Now take a look at TeV (refit) muons, and see if pT needs adjusting
      for (unsigned int iTeV=0; iTeV<tevMuons->size() ; iTeV++) { 
	reco::Muon tevMuon = tevMuons->at(iTeV) ; 
	double tevMuPt = tevMuon.pt() * ptScale ; 
	if ( tevMuPt < minPt ) continue ;  
	double dR = ROOT::Math::VectorUtil::DeltaR(tevMuon.p4(),iM.p4()) ;
	if ( dR < 0.001 ) { // Replace muon p4 with refit (TeV) muon p4
	  iM.setP4( tevMuon.p4() ) ; 
	}
      }

      muonList.push_back(iM) ; 
    }

    std::sort(muonList.begin(),muonList.end(),pTcompare()) ; 
    return muonList ; 
  }

  std::vector< std::pair<pat::Electron,float> > getElectronList(edm::Handle<pat::ElectronCollection>& pElecs,
								   double maxAbsEta, 
								   double minEtEB, double minEtEE, 
								   float ebScale, float eeScale) {
    
    std::vector< std::pair<pat::Electron,float> > electronList ; 
    for (unsigned int iElectron=0; iElectron < pElecs->size(); iElectron++) {
      pat::ElectronRef iE = pat::ElectronRef(pElecs, iElectron); 
      if ( getElectronSCEta(*iE) > maxAbsEta ) continue ; 
      if ( !passesHEEPv31(*iE) ) continue ; 

      float scale  = ( (iE->isEB()) ? ebScale : eeScale ) ; 
      float elecEt = getElectronEt(*iE) * scale ; 
      bool passEtCuts = ( (iE->isEB()) ? (elecEt > minEtEB) : (elecEt > minEtEE) ) ; 

      if (passEtCuts) { 
	std::pair<pat::Electron,float> electronCand = std::make_pair((*iE),scale) ;
	electronList.push_back( electronCand ) ; 
      }
    }
    std::sort(electronList.begin(),electronList.end(),scaleCompare()) ; 
    return electronList ; 
  }
}
