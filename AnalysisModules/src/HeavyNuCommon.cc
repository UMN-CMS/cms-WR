#include "HeavyNuCommon.h"

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
    
  }                                                // HeavyNu::isVBTFtight



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

  std::vector<double> generate_flat10_weights(const std::vector<double>& dataDist){
    // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
    const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,
				  0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
				  0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
    std::vector<double> result(25);
    double s = 0.0;
    for(int npu=0; npu<25; ++npu){
      double npu_estimated = dataDist[npu];
      result[npu] = npu_estimated / npu_probs[npu];
      s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<25; ++npu){
      result[npu] /= s;
    }
    return result;
  }


  std::vector<double> get_standard_pileup_data(int era) {
    const int npt=50;
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

    const double* pileupDist=default_pd;

    if (era>=20110 && era<=20112) pileupDist=may10_json;

    std::vector<double> retval;
    retval.reserve(npt);
    for (int i=0; i<npt; i++)
      retval.push_back(pileupDist[i]);
    return retval;
  }

 
}
