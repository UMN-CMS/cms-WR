#ifndef HEAVY_NU_COMMON_INCLUDED
#define HEAVY_NU_COMMON_INCLUDED 1

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "Math/VectorUtil.h"
#include <math.h>


namespace hnu {

    bool isVBTFloose       ( const pat::Muon& m );
    bool isVBTFtight       ( const pat::Muon& m );
    double avgVertex(const reco::JPTJet &tjet, double maxDeltaVR=1.0);
    double caloJetVertex(const pat::Jet &pJet, const reco::JPTJetCollection &jptJets, double maxDeltaVR=1.0);

}


#endif // HEAVY_NU_COMMON_INCLUDED
