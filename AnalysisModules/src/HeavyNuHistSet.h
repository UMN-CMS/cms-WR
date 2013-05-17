#ifndef HEAVY_NU_HISTSET_INCLUDED
#define HEAVY_NU_HISTSET_INCLUDED 1

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

class HeavyNuHistSet
{
public:
    HeavyNuHistSet() {}
    HeavyNuHistSet(TFileDirectory *td, const std::string& post, int histsToBook = 1);
    //book histogram set w/ common suffix inside the provided TFileDirectory
    virtual void book(TFileDirectory *, const std::string&);
    // fill all histos of the set with the two lepton candidates
    virtual void fill(HeavyNuEvent& hne);

    //book tag and probe
    virtual void tapbook(TFileDirectory *, const std::string&) {}

    // add function to fill tag and probe for muons
    virtual void tapfill(const pat::Muon& theTag,     const pat::Muon& theProbe,            const double probeTrkIso, const double wgt) {}
    virtual void tapfill(const pat::Muon& theTag,     const pat::GenericParticle& theProbe, const double trkIso,      const double wgt) {}
    virtual void tapfill(const pat::Electron& theTag, const pat::Electron& theProbe,                                  const double wgt) {}

private:
    TH1 *evtWeight;
    TH1 *ptL1, *ptL2, *ptJet1, *ptJet2;
    TH1 *etaL1pt30, *etaL1pt40, *etaL2pt30, *etaL2pt40;
    TH1 *phiL1pt30, *phiL1pt40, *phiL2pt30, *phiL2pt40;
    TH2 *ptL1VsPtL2ss, *ptL1VsPtL2os;
    TH1 *etaL1, *etaL2, *etaJet1, *etaJet2;
    TH1 *phiL1, *phiL2, *phiJet1, *phiJet2;
    TH1 *dEtaL, *dPhiL, *dEtaJet, *dPhiJet;
    TH2 *dEtaPhiL, *dEtaPhiJet;
    TH1 *dRminL1jet, *dRminL2jet, *dRminLJet;
    TH1 *hptrelL1, *hptrelL2;
    TH2 *ptrelVsdRminL1jet, *ptrelVsdRminL2jet;
    TH2 *jetID2d;

    TH1 *mindRjet_genjet, *maxdRjet_genjet, *nuLMatchedJets;

    TH1 *mLL, *mLLOS, *mLLSS, *diLCharge, *mLLZoom;
    TH1 *mWR, *mNuR1, *mNuR2, *mJJ, *ptWR, *st, *mLQmin;
    TH2 *mNuR2D, *jetPtvsNum;
    TH1 *L1ptFracWRmass;

    TH2 *mLLvsmWR, *mWRvsmNuR1, *mWRvsmNuR2, *mWRvsNPV, *mLLZoomvsNPV;

    TH3 *mLLZoomvsptL1vsptL2;

    TH1 *btagJet1, *btagJet2;
    TH1 *numBjets, *njets;
    
    TH1 *mWR_1b, *mWR_2b, *mLL_1b, *mLL_2b, *mLLZoom_1b, *mLLZoom_2b;

    TH1* met;

    // Vertex plots
    TH1* vtx_LL, *vtx_jj, *vtx_min_L1j, *vtx_min_L2j, *vtx_min_Lj, *vtx_max_dist;

    // mc type
    TH1* mc_type;

    // Pileup, vertex count information
    TH1 *n_pileup, *n_vertex, *n_vertex_noWgt;
};

#endif
