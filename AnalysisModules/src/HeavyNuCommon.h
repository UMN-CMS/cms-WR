#ifndef HEAVY_NU_COMMON_INCLUDED
#define HEAVY_NU_COMMON_INCLUDED 1

#include <vector>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "Math/VectorUtil.h"
#include <math.h>

//#define CMSSW_3XX
#define CMSSW_4XX

namespace hnu {

    bool isVBTFloose       ( const pat::Muon& m );
    bool isVBTFtight       ( const pat::Muon& m );
    double avgVertex(const reco::JPTJet &tjet, double maxDeltaVR=1.0);
    double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR=1.0);

    std::vector<double> generate_flat10_weights(const std::vector<double>& dataDist);

    std::vector<double> get_standard_pileup_data(int pileupEra); 

}


#endif // HEAVY_NU_COMMON_INCLUDED
