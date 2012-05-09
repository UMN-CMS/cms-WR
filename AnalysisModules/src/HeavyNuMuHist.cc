#include "HeavyNu/AnalysisModules/src/HeavyNuMuHist.h"


static const int muonQualityFlags = 4;
static const std::string muonQuality[4] = {"All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"};

HeavyNuMuHist::HeavyNuMuHist(TFileDirectory *td, const std::string& post, int histsToBook)
{
    if(histsToBook & 1) book(td, post);
    if(histsToBook & 2) tapbook(td, post);
}

void HeavyNuMuHist::book(TFileDirectory * td, const std::string& post)
{
    HeavyNuHistSet::book(td, post);

    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    t = "Quality (#mu_{1}) " + post;
    qualMu1 = td->make<TH1F > ("qualMu1", t.c_str(), muonQualityFlags, 0, muonQualityFlags);
    t = "Quality (#mu_{2}) " + post;
    qualMu2 = td->make<TH1F > ("qualMu2", t.c_str(), muonQualityFlags, 0, muonQualityFlags);
    for(int i = 0; i < muonQualityFlags; i++)
    {
        qualMu1->GetXaxis()->SetBinLabel(i + 1, muonQuality[i].c_str());
        qualMu2->GetXaxis()->SetBinLabel(i + 1, muonQuality[i].c_str());
    }

    // isolation

    t = "trackIso(#mu_{1}) " + post;
    mu1trackIso = td->make<TH1F > ("mu1trackIso", t.c_str(), 40, 0., 200.);
    t = "hcalIso(#mu_{1}) " + post;
    mu1hcalIso = td->make<TH1F > ("mu1hcalIso", t.c_str(), 40, 0., 200.);
    t = "ecalIso(#mu_{1}) " + post;
    mu1ecalIso = td->make<TH1F > ("mu1ecalIso", t.c_str(), 40, 0., 200.);
    t = "caloIso(#mu_{1}) " + post;
    mu1calIso = td->make<TH1F > ("mu1caloIso", t.c_str(), 40, 0., 200.);
    t = "Dxy(#mu_{1}) " + post;
    mu1dB = td->make<TH1F > ("mu1dB", t.c_str(), 50, 0., 1.);

    t = "trackIso(#mu_{2}) " + post;
    mu2trackIso = td->make<TH1F > ("mu2trackIso", t.c_str(), 40, 0., 200.);
    t = "hcalIso(#mu_{2}) " + post;
    mu2hcalIso = td->make<TH1F > ("mu2hcalIso", t.c_str(), 40, 0., 200.);
    t = "ecalIso(#mu_{2}) " + post;
    mu2ecalIso = td->make<TH1F > ("mu2ecalIso", t.c_str(), 40, 0., 200.);
    t = "caloIso(#mu_{2}) " + post;
    mu2calIso = td->make<TH1F > ("mu2caloIso", t.c_str(), 40, 0., 200.);
    t = "Dxy(#mu_{2}) " + post;
    mu2dB = td->make<TH1F > ("mu2dB", t.c_str(), 50, 0., 1.);

    t = "trackRelIso(#mu_{1}) " + post;
    mu1trackRelIso = td->make<TH1F > ("mu1trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#mu_{1}) " + post;
    mu1hcalRelIso = td->make<TH1F > ("mu1hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#mu_{1}) " + post;
    mu1ecalRelIso = td->make<TH1F > ("mu1ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#mu_{1}) " + post;
    mu1calRelIso = td->make<TH1F > ("mu1caloRelIso", t.c_str(), 50, 0., 5.);

    t = "trackRelIso(#mu_{2}) " + post;
    mu2trackRelIso = td->make<TH1F > ("mu2trackRelIso", t.c_str(), 50, 0., 5.);
    t = "hcalRelIso(#mu_{2}) " + post;
    mu2hcalRelIso = td->make<TH1F > ("mu2hcalRelIso", t.c_str(), 50, 0., 5.);
    t = "ecalRelIso(#mu_{2}) " + post;
    mu2ecalRelIso = td->make<TH1F > ("mu2ecalRelIso", t.c_str(), 50, 0., 5.);
    t = "caloRelIso(#mu_{2}) " + post;
    mu2calRelIso = td->make<TH1F > ("mu2caloRelIso", t.c_str(), 50, 0., 5.);
    
    t = "#Delta p_{T}(#L_{1},gen) " + post;
    dptL1gen = td->make<TH1F > ("dptL1gen", t.c_str(), 50, -0.50, 0.50);
    t = "#Delta p_{T}(#L_{2},gen) " + post;
    dptL2gen = td->make<TH1F > ("dptL2gen", t.c_str(), 50, -0.50, 0.50);
    t = "#Delta R(#L_{1},gen) " + post;
    dRL1gen = td->make<TH1F > ("dRL1gen", t.c_str(), 50, 0, 0.01);
    t = "#Delta R(#L_{2},gen) " + post;
    dRL2gen = td->make<TH1F > ("dRL2gen", t.c_str(), 50, 0, 0.01);

    t = "M(#L #L)(generated) " + post;
    mLLGenZoom = td->make<TH1F > ("mLLGenZoom", t.c_str(), 400, 0, 400);
}

void HeavyNuMuHist::tapbook(TFileDirectory * td, const std::string& post)
{
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    t = "TP event weight " + post;
    tpEvtWeight = td->make<TH1F > ("tpEvtWeight", t.c_str(), 1000, 0.0, 10.0);

    // ----------  Muon histograms  ----------

    t = "p_{T}(tag) " + post;
    ptTag = td->make<TH1F > ("ptTag", t.c_str(), 100, 0., 1000.);
    t = "#eta(tag) " + post;
    etaTag = td->make<TH1F > ("etaTag", t.c_str(), 50, -2.5, 2.5);
    t = "#phi(tag) " + post;
    phiTag = td->make<TH1F > ("phiTag", t.c_str(), 30, -3.14159, 3.14159);

    t = "p_{T}(probe) " + post;
    ptProbe = td->make<TH1F > ("ptProbe", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 100% probe relIso " + post;
    ptProbeRiso100 = td->make<TH1F > ("ptProbeRiso100", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 50% probe relIso " + post;
    ptProbeRiso50 = td->make<TH1F > ("ptProbeRiso50", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 20% probe relIso " + post;
    ptProbeRiso20 = td->make<TH1F > ("ptProbeRiso20", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 10% probe relIso " + post;
    ptProbeRiso10 = td->make<TH1F > ("ptProbeRiso10", t.c_str(), 100, 0., 1000.);
    t = "p_{T}(probe) 5% probe relIso " + post;
    ptProbeRiso5 = td->make<TH1F > ("ptProbeRiso5", t.c_str(), 100, 0., 1000.);

    t = "#eta(probe) " + post;
    etaProbe = td->make<TH1F > ("etaProbe", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) " + post;
    etaProbePt30 = td->make<TH1F > ("etaProbePt30", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) " + post;
    etaProbePt40 = td->make<TH1F > ("etaProbePt40", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe) 100% probe relIso " + post;
    etaProbeRiso100 = td->make<TH1F > ("etaProbeRiso100", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) 100% probe relIso " + post;
    etaProbePt30Riso100 = td->make<TH1F > ("etaProbePt30Riso100", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) 100% probe relIso " + post;
    etaProbePt40Riso100 = td->make<TH1F > ("etaProbePt40Riso100", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe) 50% probe relIso " + post;
    etaProbeRiso50 = td->make<TH1F > ("etaProbeRiso50", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) 50% probe relIso " + post;
    etaProbePt30Riso50 = td->make<TH1F > ("etaProbePt30Riso50", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) 50% probe relIso " + post;
    etaProbePt40Riso50 = td->make<TH1F > ("etaProbePt40Riso50", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe) 20% probe relIso " + post;
    etaProbeRiso20 = td->make<TH1F > ("etaProbeRiso20", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) 20% probe relIso " + post;
    etaProbePt30Riso20 = td->make<TH1F > ("etaProbePt30Riso20", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) 20% probe relIso " + post;
    etaProbePt40Riso20 = td->make<TH1F > ("etaProbePt40Riso20", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe) 10% probe relIso " + post;
    etaProbeRiso10 = td->make<TH1F > ("etaProbeRiso10", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) 10% probe relIso " + post;
    etaProbePt30Riso10 = td->make<TH1F > ("etaProbePt30Riso10", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) 10% probe relIso " + post;
    etaProbePt40Riso10 = td->make<TH1F > ("etaProbePt40Riso10", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe) 5% probe relIso " + post;
    etaProbeRiso5 = td->make<TH1F > ("etaProbeRiso5", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 30 GeV) 5% probe relIso " + post;
    etaProbePt30Riso5 = td->make<TH1F > ("etaProbePt30Riso5", t.c_str(), 50, -2.5, 2.5);
    t = "#eta(probe p_{T} > 40 GeV) 5% probe relIso " + post;
    etaProbePt40Riso5 = td->make<TH1F > ("etaProbePt40Riso5", t.c_str(), 50, -2.5, 2.5);

    t = "#phi(probe) " + post;
    phiProbe = td->make<TH1F > ("phiProbe", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) " + post;
    phiProbePt30 = td->make<TH1F > ("phiProbePt30", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) " + post;
    phiProbePt40 = td->make<TH1F > ("phiProbePt40", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 100% probe relIso " + post;
    phiProbeRiso100 = td->make<TH1F > ("phiProbeRiso100", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) 100% probe relIso " + post;
    phiProbePt30Riso100 = td->make<TH1F > ("phiProbePt30Riso100", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) 100% probe relIso " + post;
    phiProbePt40Riso100 = td->make<TH1F > ("phiProbePt40Riso100", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 50% probe relIso " + post;
    phiProbeRiso50 = td->make<TH1F > ("phiProbeRiso50", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) 50% probe relIso " + post;
    phiProbePt30Riso50 = td->make<TH1F > ("phiProbePt30Riso50", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) 50% probe relIso " + post;
    phiProbePt40Riso50 = td->make<TH1F > ("phiProbePt40Riso50", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 20% probe relIso " + post;
    phiProbeRiso20 = td->make<TH1F > ("phiProbeRiso20", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) 20% probe relIso " + post;
    phiProbePt30Riso20 = td->make<TH1F > ("phiProbePt30Riso20", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) 20% probe relIso " + post;
    phiProbePt40Riso20 = td->make<TH1F > ("phiProbePt40Riso20", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 10% probe relIso " + post;
    phiProbeRiso10 = td->make<TH1F > ("phiProbeRiso10", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) 10% probe relIso " + post;
    phiProbePt30Riso10 = td->make<TH1F > ("phiProbePt30Riso10", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) 10% probe relIso " + post;
    phiProbePt40Riso10 = td->make<TH1F > ("phiProbePt40Riso10", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe) 5% probe relIso " + post;
    phiProbeRiso5 = td->make<TH1F > ("phiProbeRiso5", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 30 GeV) 5% probe relIso " + post;
    phiProbePt30Riso5 = td->make<TH1F > ("phiProbePt30Riso5", t.c_str(), 30, -3.14159, 3.14159);
    t = "#phi(probe p_{T} > 40 GeV) 5% probe relIso " + post;
    phiProbePt40Riso5 = td->make<TH1F > ("phiProbePt40Riso5", t.c_str(), 30, -3.14159, 3.14159);

    // isolation
    t = "trackIso(tag) " + post;
    tagTrackIso = td->make<TH1F > ("tagTrackIso", t.c_str(), 200, 0., 10.);
    t = "trackIso(probe) " + post;
    probeTrackIso = td->make<TH1F > ("probeTrackIso", t.c_str(), 200, 0., 10.);
    t = "trackIso(probe p_{T} > 30 GeV) " + post;
    probePt30TrackIso = td->make<TH1F > ("probePt30TrackIso", t.c_str(), 200, 0., 10.);
    t = "trackIso(probe p_{T} > 40 GeV) " + post;
    probePt40TrackIso = td->make<TH1F > ("probePt40TrackIso", t.c_str(), 200, 0., 10.);
    t = "trackRelIso(tag) " + post;
    tagTrackRelIso = td->make<TH1F > ("tagTrackRelIso", t.c_str(), 100, 0., 2.);
    t = "trackRelIso(probe) " + post;
    probeTrackRelIso = td->make<TH1F > ("probeTrackRelIso", t.c_str(), 100, 0., 2.);
    t = "trackRelIso(probe p_{T} > 30 GeV) " + post;
    probePt30TrackRelIso = td->make<TH1F > ("probePt30TrackRelIso", t.c_str(), 100, 0., 2.);
    t = "trackRelIso(probe p_{T} > 40 GeV) " + post;
    probePt40TrackRelIso = td->make<TH1F > ("probePt40TrackRelIso", t.c_str(), 100, 0., 2.);

    // Two dimensional histograms: quantity vs. dimuon mass

    t = "p_{T}(tag) vs. #mu#mu mass " + post;
    ptTag_mass = td->make<TH2F > ("ptTag_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "#eta(tag) vs. #mu#mu mass " + post;
    etaTag_mass = td->make<TH2F > ("etaTag_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#phi(tag) vs. #mu#mu mass " + post;
    phiTag_mass = td->make<TH2F > ("phiTag_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);

    t = "p_{T}(probe) vs. #mu#mu mass " + post;
    ptProbe_mass = td->make<TH2F > ("ptProbe_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "p_{T}(probe) 100% probe relIso vs. #mu#mu mass " + post;
    ptProbeRiso100_mass = td->make<TH2F > ("ptProbeRiso100_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "p_{T}(probe) 50% probe relIso vs. #mu#mu mass " + post;
    ptProbeRiso50_mass = td->make<TH2F > ("ptProbeRiso50_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "p_{T}(probe) 20% probe relIso vs. #mu#mu mass " + post;
    ptProbeRiso20_mass = td->make<TH2F > ("ptProbeRiso20_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "p_{T}(probe) 10% probe relIso vs. #mu#mu mass " + post;
    ptProbeRiso10_mass = td->make<TH2F > ("ptProbeRiso10_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);
    t = "p_{T}(probe) 5% probe relIso vs. #mu#mu mass " + post;
    ptProbeRiso5_mass = td->make<TH2F > ("ptProbeRiso5_mass", t.c_str(), 100, 0., 1000., 60, 60., 120.);

    t = "#eta(probe) vs. #mu#mu mass " + post;
    etaProbe_mass = td->make<TH2F > ("etaProbe_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) vs. #mu#mu mass " + post;
    etaProbePt30_mass = td->make<TH2F > ("etaProbePt30_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) vs. #mu#mu mass " + post;
    etaProbePt40_mass = td->make<TH2F > ("etaProbePt40_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe) 100% probe relIso vs. #mu#mu mass " + post;
    etaProbeRiso100_mass = td->make<TH2F > ("etaProbeRiso100_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) 100% probe relIso vs. #mu#mu mass " + post;
    etaProbePt30Riso100_mass = td->make<TH2F > ("etaProbePt30Riso100_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) 100% probe relIso vs. #mu#mu mass " + post;
    etaProbePt40Riso100_mass = td->make<TH2F > ("etaProbePt40Riso100_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe) 50% probe relIso vs. #mu#mu mass " + post;
    etaProbeRiso50_mass = td->make<TH2F > ("etaProbeRiso50_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) 50% probe relIso vs. #mu#mu mass " + post;
    etaProbePt30Riso50_mass = td->make<TH2F > ("etaProbePt30Riso50_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) 50% probe relIso vs. #mu#mu mass " + post;
    etaProbePt40Riso50_mass = td->make<TH2F > ("etaProbePt40Riso50_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe) 20% probe relIso vs. #mu#mu mass " + post;
    etaProbeRiso20_mass = td->make<TH2F > ("etaProbeRiso20_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) 20% probe relIso vs. #mu#mu mass " + post;
    etaProbePt30Riso20_mass = td->make<TH2F > ("etaProbePt30Riso20_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) 20% probe relIso vs. #mu#mu mass " + post;
    etaProbePt40Riso20_mass = td->make<TH2F > ("etaProbePt40Riso20_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe) 10% probe relIso vs. #mu#mu mass " + post;
    etaProbeRiso10_mass = td->make<TH2F > ("etaProbeRiso10_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) 10% probe relIso vs. #mu#mu mass " + post;
    etaProbePt30Riso10_mass = td->make<TH2F > ("etaProbePt30Riso10_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) 10% probe relIso vs. #mu#mu mass " + post;
    etaProbePt40Riso10_mass = td->make<TH2F > ("etaProbePt40Riso10_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe) 5% probe relIso vs. #mu#mu mass " + post;
    etaProbeRiso5_mass = td->make<TH2F > ("etaProbeRiso5_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 30 GeV) 5% probe relIso vs. #mu#mu mass " + post;
    etaProbePt30Riso5_mass = td->make<TH2F > ("etaProbePt30Riso5_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);
    t = "#eta(probe p_{T} > 40 GeV) 5% probe relIso vs. #mu#mu mass " + post;
    etaProbePt40Riso5_mass = td->make<TH2F > ("etaProbePt40Riso5_mass", t.c_str(), 50, -2.5, 2.5, 60, 60., 120.);

    t = "#phi(probe) vs. #mu#mu mass " + post;
    phiProbe_mass = td->make<TH2F > ("phiProbe_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) vs. #mu#mu mass " + post;
    phiProbePt30_mass = td->make<TH2F > ("phiProbePt30_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) vs. #mu#mu mass " + post;
    phiProbePt40_mass = td->make<TH2F > ("phiProbePt40_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe) 100% probe relIso vs. #mu#mu mass " + post;
    phiProbeRiso100_mass = td->make<TH2F > ("phiProbeRiso100_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) 100% probe relIso vs. #mu#mu mass " + post;
    phiProbePt30Riso100_mass = td->make<TH2F > ("phiProbePt30Riso100_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) 100% probe relIso vs. #mu#mu mass " + post;
    phiProbePt40Riso100_mass = td->make<TH2F > ("phiProbePt40Riso100_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe) 50% probe relIso vs. #mu#mu mass " + post;
    phiProbeRiso50_mass = td->make<TH2F > ("phiProbeRiso50_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) 50% probe relIso vs. #mu#mu mass " + post;
    phiProbePt30Riso50_mass = td->make<TH2F > ("phiProbePt30Riso50_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) 50% probe relIso vs. #mu#mu mass " + post;
    phiProbePt40Riso50_mass = td->make<TH2F > ("phiProbePt40Riso50_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe) 20% probe relIso vs. #mu#mu mass " + post;
    phiProbeRiso20_mass = td->make<TH2F > ("phiProbeRiso20_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) 20% probe relIso vs. #mu#mu mass " + post;
    phiProbePt30Riso20_mass = td->make<TH2F > ("phiProbePt30Riso20_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) 20% probe relIso vs. #mu#mu mass " + post;
    phiProbePt40Riso20_mass = td->make<TH2F > ("phiProbePt40Riso20_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe) 10% probe relIso vs. #mu#mu mass " + post;
    phiProbeRiso10_mass = td->make<TH2F > ("phiProbeRiso10_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) 10% probe relIso vs. #mu#mu mass " + post;
    phiProbePt30Riso10_mass = td->make<TH2F > ("phiProbePt30Riso10_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) 10% probe relIso vs. #mu#mu mass " + post;
    phiProbePt40Riso10_mass = td->make<TH2F > ("phiProbePt40Riso10_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe) 5% probe relIso vs. #mu#mu mass " + post;
    phiProbeRiso5_mass = td->make<TH2F > ("phiProbeRiso5_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 30 GeV) 5% probe relIso vs. #mu#mu mass " + post;
    phiProbePt30Riso5_mass = td->make<TH2F > ("phiProbePt30Riso5_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);
    t = "#phi(probe p_{T} > 40 GeV) 5% probe relIso vs. #mu#mu mass " + post;
    phiProbePt40Riso5_mass = td->make<TH2F > ("phiProbePt40Riso5_mass", t.c_str(), 30, -3.14159, 3.14159, 60, 60., 120.);


    // isolation
    t = "trackIso(tag) vs. #mu#mu mass " + post;
    tagTrackIso_mass = td->make<TH2F > ("tagTrackIso_mass", t.c_str(), 200, 0., 10., 60, 60., 120.);
    t = "trackIso(probe) vs. #mu#mu mass " + post;
    probeTrackIso_mass = td->make<TH2F > ("probeTrackIso_mass", t.c_str(), 200, 0., 10., 60, 60., 120.);
    t = "trackIso(probe p_{T} > 30 GeV) vs. #mu#mu mass " + post;
    probePt30TrackIso_mass = td->make<TH2F > ("probePt30TrackIso_mass", t.c_str(), 200, 0., 10., 60, 60., 120.);
    t = "trackIso(probe p_{T} > 40 GeV) vs. #mu#mu mass " + post;
    probePt40TrackIso_mass = td->make<TH2F > ("probePt40TrackIso_mass", t.c_str(), 200, 0., 10., 60, 60., 120.);
    t = "trackRelIso(tag) vs. #mu#mu mass " + post;
    tagTrackRelIso_mass = td->make<TH2F > ("tagTrackRelIso_mass", t.c_str(), 100, 0., 2., 60, 60., 120.);
    t = "trackRelIso(probe) vs. #mu#mu mass " + post;
    probeTrackRelIso_mass = td->make<TH2F > ("probeTrackRelIso_mass", t.c_str(), 100, 0., 2., 60, 60., 120.);
    t = "trackRelIso(probe p_{T} > 30 GeV) vs. #mu#mu mass " + post;
    probePt30TrackRelIso_mass = td->make<TH2F > ("probePt30TrackRelIso_mass", t.c_str(), 100, 0., 2., 60, 60., 120.);
    t = "trackRelIso(probe p_{T} > 40 GeV) vs. #mu#mu mass " + post;
    probePt40TrackRelIso_mass = td->make<TH2F > ("probePt40TrackRelIso_mass", t.c_str(), 100, 0., 2., 60, 60., 120.);


    // ----------  Composite histograms  ----------
    t = "M(#mu #mu) " + post;
    mMuMuTP = td->make<TH1F > ("mMuMuTP", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) (probe p_{T} > 30 GeV)" + post;
    mMuMuPt30TP = td->make<TH1F > ("mMuMuPt30TP", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) (probe p_{T} > 40 GeV)" + post;
    mMuMuPt40TP = td->make<TH1F > ("mMuMuPt40TP", t.c_str(), 60, 60, 120);

    t = "M(#mu #mu) 5% probe relIso " + post;
    mMuMuTPRiso5 = td->make<TH1F > ("mMuMuTPRiso5", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 5% probe (p_{T} > 30 GeV) relIso " + post;
    mMuMuPt30TPRiso5 = td->make<TH1F > ("mMuMuPt30TPRiso5", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 5% probe (p_{T} > 40 GeV) relIso " + post;
    mMuMuPt40TPRiso5 = td->make<TH1F > ("mMuMuPt40TPRiso5", t.c_str(), 60, 60, 120);

    t = "M(#mu #mu) 10% probe relIso " + post;
    mMuMuTPRiso10 = td->make<TH1F > ("mMuMuTPRiso10", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 10% probe (p_{T} > 30 GeV) relIso " + post;
    mMuMuPt30TPRiso10 = td->make<TH1F > ("mMuMuPt30TPRiso10", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 10% probe (p_{T} > 40 GeV) relIso " + post;
    mMuMuPt40TPRiso10 = td->make<TH1F > ("mMuMuPt40TPRiso10", t.c_str(), 60, 60, 120);

    t = "M(#mu #mu) 20% probe relIso " + post;
    mMuMuTPRiso20 = td->make<TH1F > ("mMuMuTPRiso20", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 20% probe (p_{T} > 30 GeV) relIso " + post;
    mMuMuPt30TPRiso20 = td->make<TH1F > ("mMuMuPt30TPRiso20", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 20% probe (p_{T} > 40 GeV) relIso " + post;
    mMuMuPt40TPRiso20 = td->make<TH1F > ("mMuMuPt40TPRiso20", t.c_str(), 60, 60, 120);

    t = "M(#mu #mu) 50% probe relIso " + post;
    mMuMuTPRiso50 = td->make<TH1F > ("mMuMuTPRiso50", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 50% probe (p_{T} > 30 GeV) relIso " + post;
    mMuMuPt30TPRiso50 = td->make<TH1F > ("mMuMuPt30TPRiso50", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 50% probe (p_{T} > 40 GeV) relIso " + post;
    mMuMuPt40TPRiso50 = td->make<TH1F > ("mMuMuPt40TPRiso50", t.c_str(), 60, 60, 120);

    t = "M(#mu #mu) 100% probe relIso " + post;
    mMuMuTPRiso100 = td->make<TH1F > ("mMuMuTPRiso100", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 100% probe (p_{T} > 30 GeV) relIso " + post;
    mMuMuPt30TPRiso100 = td->make<TH1F > ("mMuMuPt30TPRiso100", t.c_str(), 60, 60, 120);
    t = "M(#mu #mu) 100% probe (p_{T} > 40 GeV) relIso " + post;
    mMuMuPt40TPRiso100 = td->make<TH1F > ("mMuMuPt40TPRiso100", t.c_str(), 60, 60, 120);
}

void HeavyNuMuHist::fill(HeavyNuEvent& hne)
{
    HeavyNuHistSet::fill(hne);

    double wgt = hne.eventWgt;

    if(hne.nLeptons >= 1)
    {
        mu1trackIso->Fill(hne.mu1.trackIso(), wgt);
        mu1hcalIso->Fill(hne.mu1.hcalIso(), wgt);
        mu1ecalIso->Fill(hne.mu1.ecalIso(), wgt);
        mu1calIso->Fill(hne.mu1.caloIso(), wgt);
        mu1dB->Fill(hne.mu1.dB(), wgt);

        mu1trackRelIso->Fill(hne.mu1.trackIso() / hne.mu1.pt(), wgt);
        mu1hcalRelIso->Fill(hne.mu1.hcalIso() / hne.mu1.pt(), wgt);
        mu1ecalRelIso->Fill(hne.mu1.ecalIso() / hne.mu1.pt(), wgt);
        mu1calRelIso->Fill(hne.mu1.caloIso() / hne.mu1.pt(), wgt);
    }

    if(hne.nLeptons >= 2)
    {
        mu2trackIso->Fill(hne.mu2.trackIso(), wgt);
        mu2hcalIso->Fill(hne.mu2.hcalIso(), wgt);
        mu2ecalIso->Fill(hne.mu2.ecalIso(), wgt);
        mu2calIso->Fill(hne.mu2.caloIso(), wgt);
        mu2dB->Fill(hne.mu2.dB(), wgt);

        mu2trackRelIso->Fill(hne.mu2.trackIso() / hne.mu2.pt(), wgt);
        mu2hcalRelIso->Fill(hne.mu2.hcalIso() / hne.mu2.pt(), wgt);
        mu2ecalRelIso->Fill(hne.mu2.ecalIso() / hne.mu2.pt(), wgt);
        mu2calRelIso->Fill(hne.mu2.caloIso() / hne.mu2.pt(), wgt);
    }

    for(int i = 0; i < muonQualityFlags; i++)
    {
        if(hne.nLeptons >= 1 && hne.mu1.muonID(muonQuality[i])) qualMu1->Fill(i, wgt);
        if(hne.nLeptons >= 2 && hne.mu2.muonID(muonQuality[i])) qualMu2->Fill(i, wgt);
    }

    if(hne.isMC)
    {
        if(hne.nLeptons >= 1 && hne.mu1.genLepton() != 0)
        {
            float dpt = hne.mu1.pt() - hne.mu1.genLepton()->pt();
            float dR = deltaR(hne.mu1.eta(), hne.mu1.phi(), hne.mu1.genLepton()->eta(), hne.mu1.genLepton()->phi());
            dptL1gen->Fill(dpt / hne.mu1.genLepton()->pt());
            dRL1gen->Fill(dR);
        }
        if(hne.nLeptons >= 2 && hne.mu2.genLepton() != 0)
        {
            float dpt = hne.mu2.pt() - hne.mu2.genLepton()->pt();
            float dR = deltaR(hne.mu2.eta(), hne.mu2.phi(), hne.mu2.genLepton()->eta(), hne.mu2.genLepton()->phi());
            dptL2gen->Fill(dpt / hne.mu2.genLepton()->pt());
            dRL2gen->Fill(dR);
        }
        if(hne.nLeptons >= 2 && (hne.mu1.genLepton() != 0) && (hne.mu2.genLepton() != 0))
        {
            reco::Particle::LorentzVector L1gp4 = hne.mu1.genLepton()->p4();
            reco::Particle::LorentzVector L2gp4 = hne.mu2.genLepton()->p4();
            mLLGenZoom->Fill((L1gp4 + L2gp4).M());
        }
    }

}

void HeavyNuMuHist::tapfill(const pat::Muon& theTag, const pat::Muon& theProbe, const double probeTrkIso, const double wgt)
{
    tpEvtWeight->Fill(wgt) ;
    reco::Particle::LorentzVector mumu = theTag.p4() + theProbe.p4();

    ptTag->Fill(theTag.pt(), wgt);
    etaTag->Fill(theTag.eta(), wgt);
    phiTag->Fill(theTag.phi(), wgt);
    ptTag_mass->Fill(theTag.pt(), mumu.M(), wgt);
    etaTag_mass->Fill(theTag.eta(), mumu.M(), wgt);
    phiTag_mass->Fill(theTag.phi(), mumu.M(), wgt);

    ptProbe->Fill(theProbe.pt(), wgt);
    etaProbe->Fill(theProbe.eta(), wgt);
    phiProbe->Fill(theProbe.phi(), wgt);
    ptProbe_mass->Fill(theProbe.pt(), mumu.M(), wgt);
    etaProbe_mass->Fill(theProbe.eta(), mumu.M(), wgt);
    phiProbe_mass->Fill(theProbe.phi(), mumu.M(), wgt);

    tagTrackIso ->Fill(theTag.trackIso(), wgt);
    probeTrackIso ->Fill(probeTrkIso, wgt);
    tagTrackIso_mass ->Fill(theTag.trackIso(), mumu.M(), wgt);
    probeTrackIso_mass ->Fill(probeTrkIso, mumu.M(), wgt);

    tagTrackRelIso->Fill(theTag.trackIso() / theTag.pt(), wgt);
    probeTrackRelIso->Fill(probeTrkIso / theProbe.pt(), wgt);
    tagTrackRelIso_mass->Fill(theTag.trackIso() / theTag.pt(), mumu.M(), wgt);
    probeTrackRelIso_mass->Fill(probeTrkIso / theProbe.pt(), mumu.M(), wgt);

    mMuMuTP->Fill(mumu.M(), wgt);

    // Special histograms for extra checking
    if ( theProbe.pt() > 30. )
    {
        etaProbePt30->Fill(theProbe.eta(), wgt) ;
        phiProbePt30->Fill(theProbe.phi(), wgt) ;
        etaProbePt30_mass->Fill(theProbe.eta(), mumu.M(), wgt) ;
        phiProbePt30_mass->Fill(theProbe.phi(), mumu.M(), wgt) ;

        probePt30TrackIso ->Fill(probeTrkIso, wgt);
        probePt30TrackIso_mass ->Fill(probeTrkIso, mumu.M(), wgt);

        probePt30TrackRelIso->Fill(probeTrkIso / theProbe.pt(), wgt);
        probePt30TrackRelIso_mass->Fill(probeTrkIso / theProbe.pt(), mumu.M(), wgt);

        mMuMuPt30TP->Fill(mumu.M(), wgt);

        if ( theProbe.pt() > 40. )
        {
            etaProbePt40->Fill(theProbe.eta(), wgt) ;
            phiProbePt40->Fill(theProbe.phi(), wgt) ;
            etaProbePt40_mass->Fill(theProbe.eta(), mumu.M(), wgt) ;
            phiProbePt40_mass->Fill(theProbe.phi(), mumu.M(), wgt) ;

            probePt40TrackIso ->Fill(probeTrkIso, wgt);
            probePt40TrackIso_mass ->Fill(probeTrkIso, mumu.M(), wgt);

            probePt40TrackRelIso->Fill(probeTrkIso / theProbe.pt(), wgt);
            probePt40TrackRelIso_mass->Fill(probeTrkIso / theProbe.pt(), mumu.M(), wgt);

            mMuMuPt40TP->Fill(mumu.M(), wgt);
        }
    }

    double probe_relIso = probeTrkIso / theProbe.pt() ;
    if ( probe_relIso < 1.0 )
    {
        ptProbeRiso100->Fill(theProbe.pt(), wgt);
        etaProbeRiso100->Fill(theProbe.eta(), wgt);
        phiProbeRiso100->Fill(theProbe.phi(), wgt);
        ptProbeRiso100_mass->Fill(theProbe.pt(), mumu.M(), wgt);
        etaProbeRiso100_mass->Fill(theProbe.eta(), mumu.M(), wgt);
        phiProbeRiso100_mass->Fill(theProbe.phi(), mumu.M(), wgt);
        mMuMuTPRiso100->Fill(mumu.M(), wgt);
        if ( probe_relIso < 0.5 )
        {
            ptProbeRiso50->Fill(theProbe.pt(), wgt);
            etaProbeRiso50->Fill(theProbe.eta(), wgt);
            phiProbeRiso50->Fill(theProbe.phi(), wgt);
            ptProbeRiso50_mass->Fill(theProbe.pt(), mumu.M(), wgt);
            etaProbeRiso50_mass->Fill(theProbe.eta(), mumu.M(), wgt);
            phiProbeRiso50_mass->Fill(theProbe.phi(), mumu.M(), wgt);
            mMuMuTPRiso50->Fill(mumu.M(), wgt);
            if ( probe_relIso < 0.2 )
            {
                ptProbeRiso20->Fill(theProbe.pt(), wgt);
                etaProbeRiso20->Fill(theProbe.eta(), wgt);
                phiProbeRiso20->Fill(theProbe.phi(), wgt);
                ptProbeRiso20_mass->Fill(theProbe.pt(), mumu.M(), wgt);
                etaProbeRiso20_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                phiProbeRiso20_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                mMuMuTPRiso20->Fill(mumu.M(), wgt);
                if ( probe_relIso < 0.1 )
                {
                    ptProbeRiso10->Fill(theProbe.pt(), wgt);
                    etaProbeRiso10->Fill(theProbe.eta(), wgt);
                    phiProbeRiso10->Fill(theProbe.phi(), wgt);
                    ptProbeRiso10_mass->Fill(theProbe.pt(), mumu.M(), wgt);
                    etaProbeRiso10_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                    phiProbeRiso10_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                    mMuMuTPRiso10->Fill(mumu.M(), wgt);
                    if ( probe_relIso < 0.05 )
                    {
                        ptProbeRiso5->Fill(theProbe.pt(), wgt);
                        etaProbeRiso5->Fill(theProbe.eta(), wgt);
                        phiProbeRiso5->Fill(theProbe.phi(), wgt);
                        ptProbeRiso5_mass->Fill(theProbe.pt(), mumu.M(), wgt);
                        etaProbeRiso5_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                        phiProbeRiso5_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                        mMuMuTPRiso5->Fill(mumu.M(), wgt);
                    }
                }
            }
        }
    }

    if ( theProbe.pt() > 30. )
    {
        if ( probe_relIso < 1.0 )
        {
            etaProbePt30Riso100->Fill(theProbe.eta(), wgt);
            phiProbePt30Riso100->Fill(theProbe.phi(), wgt);
            etaProbePt30Riso100_mass->Fill(theProbe.eta(), mumu.M(), wgt);
            phiProbePt30Riso100_mass->Fill(theProbe.phi(), mumu.M(), wgt);
            mMuMuPt30TPRiso100->Fill(mumu.M(), wgt);
            if ( probe_relIso < 0.5 )
            {
                etaProbePt30Riso50->Fill(theProbe.eta(), wgt);
                phiProbePt30Riso50->Fill(theProbe.phi(), wgt);
                etaProbePt30Riso50_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                phiProbePt30Riso50_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                mMuMuPt30TPRiso50->Fill(mumu.M(), wgt);
                if ( probe_relIso < 0.2 )
                {
                    etaProbePt30Riso20->Fill(theProbe.eta(), wgt);
                    phiProbePt30Riso20->Fill(theProbe.phi(), wgt);
                    etaProbePt30Riso20_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                    phiProbePt30Riso20_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                    mMuMuPt30TPRiso20->Fill(mumu.M(), wgt);
                    if ( probe_relIso < 0.1 )
                    {
                        etaProbePt30Riso10->Fill(theProbe.eta(), wgt);
                        phiProbePt30Riso10->Fill(theProbe.phi(), wgt);
                        etaProbePt30Riso10_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                        phiProbePt30Riso10_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                        mMuMuPt30TPRiso10->Fill(mumu.M(), wgt);
                        if ( probe_relIso < 0.05 )
                        {
                            etaProbePt30Riso5->Fill(theProbe.eta(), wgt);
                            phiProbePt30Riso5->Fill(theProbe.phi(), wgt);
                            etaProbePt30Riso5_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                            phiProbePt30Riso5_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                            mMuMuPt30TPRiso5->Fill(mumu.M(), wgt);
                        }
                    }
                }
            }
        }
        if ( theProbe.pt() > 40. )
        {
            if ( probe_relIso < 1.0 )
            {
                etaProbePt40Riso100->Fill(theProbe.eta(), wgt);
                phiProbePt40Riso100->Fill(theProbe.phi(), wgt);
                etaProbePt40Riso100_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                phiProbePt40Riso100_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                mMuMuPt40TPRiso100->Fill(mumu.M(), wgt);
                if ( probe_relIso < 0.5 )
                {
                    etaProbePt40Riso50->Fill(theProbe.eta(), wgt);
                    phiProbePt40Riso50->Fill(theProbe.phi(), wgt);
                    etaProbePt40Riso50_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                    phiProbePt40Riso50_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                    mMuMuPt40TPRiso50->Fill(mumu.M(), wgt);
                    if ( probe_relIso < 0.2 )
                    {
                        etaProbePt40Riso20->Fill(theProbe.eta(), wgt);
                        phiProbePt40Riso20->Fill(theProbe.phi(), wgt);
                        etaProbePt40Riso20_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                        phiProbePt40Riso20_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                        mMuMuPt40TPRiso20->Fill(mumu.M(), wgt);
                        if ( probe_relIso < 0.1 )
                        {
                            etaProbePt40Riso10->Fill(theProbe.eta(), wgt);
                            phiProbePt40Riso10->Fill(theProbe.phi(), wgt);
                            etaProbePt40Riso10_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                            phiProbePt40Riso10_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                            mMuMuPt40TPRiso10->Fill(mumu.M(), wgt);
                            if ( probe_relIso < 0.05 )
                            {
                                etaProbePt40Riso5->Fill(theProbe.eta(), wgt);
                                phiProbePt40Riso5->Fill(theProbe.phi(), wgt);
                                etaProbePt40Riso5_mass->Fill(theProbe.eta(), mumu.M(), wgt);
                                phiProbePt40Riso5_mass->Fill(theProbe.phi(), mumu.M(), wgt);
                                mMuMuPt40TPRiso5->Fill(mumu.M(), wgt);
                            }
                        }
                    }
                }
            }
        }
    }
}

void HeavyNuMuHist::tapfill(const pat::Muon& theTag, const pat::GenericParticle& theProbe, const double trkIso, const double wgt)
{
    pat::Muon muonProbe;
    muonProbe.setP4( theProbe.p4() );
    HeavyNuMuHist::tapfill( theTag, muonProbe, trkIso, wgt );
}
