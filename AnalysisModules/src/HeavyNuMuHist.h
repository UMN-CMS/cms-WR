//Adds specialized histograms for the muon analysis

#ifndef HEAVY_NU_MU_HISTSET_INCLUDED
#define HEAVY_NU_MU_HISTSET_INCLUDED 1

#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuHistSet.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuTrigger.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <string>

class HeavyNuMuHist : public HeavyNuHistSet
{
public:
    HeavyNuMuHist(TFileDirectory *td, const std::string& post, int histsToBook = 1);

    //book histogram set w/ common suffix inside the provided TFileDirectory
    void book(TFileDirectory *, const std::string&);

    //book tag and probe
    void tapbook(TFileDirectory *, const std::string&);

    // fill all histos of the set with the two lepton candidates
    void fill(HeavyNuEvent& hne);

    // add function to fill tag and probe
    void tapfill(const pat::Muon& theTag, const pat::Muon& theProbe,            const double probeTrkIso, const double wgt);
    void tapfill(const pat::Muon& theTag, const pat::GenericParticle& theProbe, const double trkIso,      const double wgt);

    HeavyNuTrigger::trigHistos_t trigHistos;

private:
    //muon specific histograms
    TH1 *qualMu1, *qualMu2;

    TH1 *mu1trackIso, *mu1hcalIso, *mu1ecalIso, *mu1calIso, *mu1dB;
    TH1 *mu2trackIso, *mu2hcalIso, *mu2ecalIso, *mu2calIso, *mu2dB;

    TH1 *mu1trackRelIso, *mu1hcalRelIso, *mu1ecalRelIso, *mu1calRelIso;
    TH1 *mu2trackRelIso, *mu2hcalRelIso, *mu2ecalRelIso, *mu2calRelIso;

    TH1 *dptL1gen, *dptL2gen;
    TH1 *dRL1gen, *dRL2gen;
    TH1 *mLLGenZoom;

    //tag and probe histograms
    TH1 *tpEvtWeight;
    TH1 *ptTag, *etaTag, *phiTag;
    TH1 *ptProbe, *etaProbe, *phiProbe;
    TH1 *etaProbePt30, *etaProbePt40, *phiProbePt30, *phiProbePt40;
    TH1 *ptProbeRiso100, *etaProbeRiso100, *phiProbeRiso100;
    TH1 *etaProbePt30Riso100, *phiProbePt30Riso100;
    TH1 *etaProbePt40Riso100, *phiProbePt40Riso100;
    TH1 *ptProbeRiso50, *etaProbeRiso50, *phiProbeRiso50;
    TH1 *etaProbePt30Riso50, *phiProbePt30Riso50;
    TH1 *etaProbePt40Riso50, *phiProbePt40Riso50;
    TH1 *ptProbeRiso20, *etaProbeRiso20, *phiProbeRiso20;
    TH1 *etaProbePt30Riso20, *phiProbePt30Riso20;
    TH1 *etaProbePt40Riso20, *phiProbePt40Riso20;
    TH1 *ptProbeRiso10, *etaProbeRiso10, *phiProbeRiso10;
    TH1 *etaProbePt30Riso10, *phiProbePt30Riso10;
    TH1 *etaProbePt40Riso10, *phiProbePt40Riso10;
    TH1 *ptProbeRiso5, *etaProbeRiso5, *phiProbeRiso5;
    TH1 *etaProbePt30Riso5, *phiProbePt30Riso5;
    TH1 *etaProbePt40Riso5, *phiProbePt40Riso5;
    TH1 *tagTrackIso, *tagTrackRelIso, *probeTrackIso, *probeTrackRelIso;
    TH1 *probePt30TrackIso, *probePt30TrackRelIso, *probePt40TrackIso, *probePt40TrackRelIso;
    TH1 *mMuMuTP, *mMuMuTPRiso100, *mMuMuTPRiso50, *mMuMuTPRiso20, *mMuMuTPRiso10, *mMuMuTPRiso5;
    TH1 *mMuMuPt30TP, *mMuMuPt30TPRiso100, *mMuMuPt30TPRiso50, *mMuMuPt30TPRiso20, *mMuMuPt30TPRiso10, *mMuMuPt30TPRiso5;
    TH1 *mMuMuPt40TP, *mMuMuPt40TPRiso100, *mMuMuPt40TPRiso50, *mMuMuPt40TPRiso20, *mMuMuPt40TPRiso10, *mMuMuPt40TPRiso5;

    TH2 *ptTag_mass, *etaTag_mass, *phiTag_mass;
    TH2 *ptProbe_mass, *etaProbe_mass, *phiProbe_mass;
    TH2 *etaProbePt30_mass, *etaProbePt40_mass, *phiProbePt30_mass, *phiProbePt40_mass;
    TH2 *ptProbeRiso100_mass, *etaProbeRiso100_mass, *phiProbeRiso100_mass;
    TH2 *etaProbePt30Riso100_mass, *phiProbePt30Riso100_mass;
    TH2 *etaProbePt40Riso100_mass, *phiProbePt40Riso100_mass;
    TH2 *ptProbeRiso50_mass, *etaProbeRiso50_mass, *phiProbeRiso50_mass;
    TH2 *etaProbePt30Riso50_mass, *phiProbePt30Riso50_mass;
    TH2 *etaProbePt40Riso50_mass, *phiProbePt40Riso50_mass;
    TH2 *ptProbeRiso20_mass, *etaProbeRiso20_mass, *phiProbeRiso20_mass;
    TH2 *etaProbePt30Riso20_mass, *phiProbePt30Riso20_mass;
    TH2 *etaProbePt40Riso20_mass, *phiProbePt40Riso20_mass;
    TH2 *ptProbeRiso10_mass, *etaProbeRiso10_mass, *phiProbeRiso10_mass;
    TH2 *etaProbePt30Riso10_mass, *phiProbePt30Riso10_mass;
    TH2 *etaProbePt40Riso10_mass, *phiProbePt40Riso10_mass;
    TH2 *ptProbeRiso5_mass, *etaProbeRiso5_mass, *phiProbeRiso5_mass;
    TH2 *etaProbePt30Riso5_mass, *phiProbePt30Riso5_mass;
    TH2 *etaProbePt40Riso5_mass, *phiProbePt40Riso5_mass;
    TH2 *tagTrackIso_mass, *tagTrackRelIso_mass, *probeTrackIso_mass, *probeTrackRelIso_mass;
    TH2 *probePt30TrackIso_mass, *probePt30TrackRelIso_mass, *probePt40TrackIso_mass, *probePt40TrackRelIso_mass;


};

#endif