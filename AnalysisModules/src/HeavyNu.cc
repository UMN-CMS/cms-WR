// -*- C++ -*-
//
// Package:    HeavyNu
// Class:      HeavyNu
// 
/**\class HeavyNu HeavyNu.cc HeavyNu/AnalyzerModules/src/HeavyNu.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Jeremy M Mans
//         Created:  Mon May 31 07:00:26 CDT 2010
// $Id: HeavyNu.cc,v 1.123 2013/04/22 21:42:45 bdahmes Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
#include <TRandom.h>

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID 
// Order valid for 38X only.  Can be moved after Frameworkfwd.h in 39X
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
// Needed for 39X
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/FileBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TLorentzVector.h"

#include "HeavyNu/AnalysisModules/src/HeavyNu.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEff.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "HeavyNu.h"

template <class T>
inline std::string int2str(T i)
{
    std::ostringstream ss;
    ss << i;
    return ss.str();
}

template <class T> void outputCandidate(const T& p)
{
    std::cout << "pt=" << p.pt() << " GeV, eta=" << p.eta() << ", phi=" << p.phi();
}

static std::string btagName;
static double minBtagDiscVal; // for discriminating B-tagged jets.

const int muonQualityFlags = 4;
const std::string muonQuality[] = {
    "All", "AllGlobalMuons", "AllStandAloneMuons", "AllTrackerMuons"
};

inline void labelMuonQualAxis(TAxis *ax){
    //for(int i = 0; i < muonQualityFlags; i++)
    //{
    //    ax->SetBinLabel(i + 1, muonQuality[i].c_str());
    //    ax->SetBinLabel(i + 1, muonQuality[i].c_str());
    //}
}

bool HeavyNu::isWrDaughter(const reco::Candidate* mother, int pdgid)
{
    for(size_t i = 0; i < mother->numberOfMothers(); i++)
    {
        if((mother->mother(i)->pdgId() == pdgid || isWrDaughter(mother->mother(i), pdgid))) return true;
    }
    return false;
}

reco::Candidate* HeavyNu::findGenTau(edm::Handle<reco::GenParticleCollection> gpc, bool primary)
{
    for(reco::GenParticleCollection::const_iterator igp = gpc->begin(); igp != gpc->end(); ++igp)
    {
        if(abs(igp->pdgId()) == 15 && igp->numberOfDaughters() >= 1)
        {
            if( primary && isWrDaughter((reco::Candidate*)&(*igp), 9900024) && !isWrDaughter((reco::Candidate*)&(*igp), 9900016)) return (reco::Candidate*)&(*igp);
            if(!primary && isWrDaughter((reco::Candidate*)&(*igp), 9900024) &&  isWrDaughter((reco::Candidate*)&(*igp), 9900016)) return (reco::Candidate*)&(*igp);
        }
    }
    return NULL;
}

void HeavyNu::fill(pat::MuonCollection muons,
                   pat::ElectronCollection electrons,
                   pat::JetCollection jets,
                   const HeavyNuEvent& hnuEvent,
                   const bool goodJets,
                   const bool goodLeps,
                   HeavyNuHistSet *hnmh)
{
    if(analysisMode_ == HeavyNuEvent::TAUX) return;

    HeavyNuEvent hne(hnuEvent);

    std::sort(muons.begin(), muons.end(), hnu::pTcompare());
    std::sort(electrons.begin(), electrons.end(), hnu::pTcompare());
    std::sort(jets.begin(), jets.end(), hnu::pTcompare());

    if(!goodLeps)
    {
        hne.nLeptons = 0;
        if((analysisMode_ == HeavyNuEvent::HNUMU || analysisMode_ == HeavyNuEvent::TOP) && muons.size() > 0)
        {
            hne.mu1 = muons[0];
            hne.nLeptons++;
        }
        if((analysisMode_ == HeavyNuEvent::HNUE || analysisMode_ == HeavyNuEvent::TOP) && electrons.size() > 0)
        {
            hne.e1 = electrons[0];
            hne.nLeptons++;
        }

        if(analysisMode_ == HeavyNuEvent::HNUMU && muons.size() > 1)
        {
            hne.mu2 = muons[1];
            hne.nLeptons++;
        }
        if(analysisMode_ == HeavyNuEvent::HNUE && electrons.size() > 1)
        {
            hne.e2 = electrons[1];
            hne.nLeptons++;
        }
    }

    hne.btagName = btagName;
    if(!goodJets)
    {
        if(jets.size() > 0)
        {
            hne.j1 = jets[0];
            double j1bdisc = hne.j1.bDiscriminator(btagName);
            hne.isBJet1 = j1bdisc >= minBtagDiscVal;
            hne.nJets = 1;
        }

        if(jets.size() > 1)
        {
            hne.j2 = jets[1];
            double j1bdisc = hne.j1.bDiscriminator(btagName);
            hne.isBJet2 = j1bdisc >= minBtagDiscVal;
            hne.nJets++;
        }
    }
    //hne.met1 = metc[0];

    //hne.calculateLL(correctEscale_);
    hne.calculate(correctEscale_);

    hnmh->fill(hne);
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HeavyNu::HeavyNu(const edm::ParameterSet& iConfig)
{
    // ==================== Get all parameters ====================
    //
    dolog_ = iConfig.getParameter<bool>("DoLog");

    muonTag_ = iConfig.getParameter< edm::InputTag > ("muonTag");
    jetTag_ = iConfig.getParameter< edm::InputTag > ("jetTag");
    metTag_ = iConfig.getParameter< edm::InputTag > ("metTag");
    elecTag_ = iConfig.getParameter< edm::InputTag > ("electronTag");
    trackTag_ = iConfig.getParameter< edm::InputTag > ("trackTag");

    elecRhoTag_ = iConfig.getParameter< edm::InputTag > ("electronRho");

    btagName = iConfig.getParameter<std::string > ("BtagName");
    minBtagDiscVal = iConfig.getParameter<double>("minBtagDiscr");

    cuts.minimum_mu1_pt = iConfig.getParameter<double>("minMu1pt");
    cuts.minimum_mu2_pt = iConfig.getParameter<double>("minMu2pt");
    cuts.minimum_jet_pt = iConfig.getParameter<double>("minJetPt");
    cuts.maximum_mu_abseta = iConfig.getParameter<double>("maxMuAbsEta");
    cuts.maximum_elec_abseta = iConfig.getParameter<double>("maxElecAbsEta");
    cuts.maximum_jet_abseta = iConfig.getParameter<double>("maxJetAbsEta");
    cuts.minimum_mumu_mass = iConfig.getParameter<double>("minMuMuMass");
    cuts.minimum_mWR_mass = iConfig.getParameter<double>("min4objMass");
    cuts.minimum_muon_jet_dR = iConfig.getParameter<double>("minMuonJetdR");
    cuts.qcd_jet_muon_dRoverlap = iConfig.getParameter<double>("qcdMuJetOverlap");
    cuts.muon_trackiso_limit = iConfig.getParameter<double>("muonTrackRelIsoLimit");
    cuts.maxVertexZsep = iConfig.getParameter<double>("maxVertexZsepCM");
    cuts.maxJetVZsepCM = iConfig.getParameter<double>("maxJetVZsepCM");

    effinfo_.ZwinMinGeV = iConfig.getParameter<double>("ZmassWinMinGeV");
    effinfo_.ZwinMaxGeV = iConfig.getParameter<double>("ZmassWinMaxGeV");

    jecVal_ = iConfig.getParameter<int>("jecEra");

    pileupEra_ = iConfig.getParameter<int>("pileupEra");
    puShift_   = iConfig.getParameter<int>("systPileupShift") ;
    disableTriggerCorrection_ = iConfig.getParameter<bool>("DisableTriggerCorrection");

    applyJECUsign_ = iConfig.getParameter<int>("applyJECUsign");
    if(applyJECUsign_) applyJECUsign_ /= abs(applyJECUsign_); // ensure -1,0,+1

    applyJERsign_ = iConfig.getParameter<int>("applyJERsign");
    if(applyJERsign_) applyJERsign_ /= abs(applyJERsign_); // ensure -1,0,+1

    applyMESfactor_ = iConfig.getParameter<double>("applyMESfactor");
    merUncertainty_ = iConfig.getParameter<bool>("checkMERUnc");
    correctEscale_  = iConfig.getParameter<bool>("correctEscale");

    applyMuIDCorrections_ = iConfig.getParameter<bool>("applyMuIDEffcorr");
    applyMuIDEffsign_ = iConfig.getParameter<int>("applyMuIDEffsign");
    if(applyMuIDEffsign_) applyMuIDEffsign_ /= abs(applyMuIDEffsign_); // ensure -1,0,+1

    highestPtTriggerOnly_ = iConfig.getParameter<bool>("highestPtTriggerOnly");
    applyTrigEffsign_ = iConfig.getParameter<int>("applyTrigEffsign");
    if(applyTrigEffsign_) applyTrigEffsign_ /= abs(applyTrigEffsign_); // ensure -1,0,+1

    studyMuonSelectionEff_ = iConfig.getParameter<bool>("studyMuSelectEff");
    effinfo_.oneTP = iConfig.getParameter<bool>("oneTPcand");
    if (studyMuonSelectionEff_ && effinfo_.oneTP)
    {
        tpSeed_   = iConfig.getParameter<int>("tpRandomSeed");
        effinfo_.tpRandom = new TRandom(tpSeed_);
    }

    studyScaleFactorEvolution_ = iConfig.getParameter<bool>("studyScaleFactor");
    studyAlternativeSelection_ = iConfig.getParameter<bool>("alternativeSelections");

    isPFJets_ = iConfig.getParameter<bool>("isPFJets");
    studyRatePerRun_ = iConfig.getParameter<bool>("studyRatePerRun");

    addSlopeTree_ = iConfig.getUntrackedParameter<bool>("addSlopeTree");

    randseed_ =  iConfig.getUntrackedParameter<unsigned int>("randseed");

    // Default HEEP version is 4.0 (2012 selection)
    heepVersion_ = iConfig.getUntrackedParameter<int>("heepVersion", 41);
    if ( heepVersion_ < 30 || heepVersion_ > 41 )
    {
        std::cout << "!!!!!!!!INVALID HEEP VERSION: " << heepVersion_ << " (setting HEEP 4.1)!!!!!!!!" << std::endl;
        heepVersion_ = 41;
    }


    std::string am = iConfig.getUntrackedParameter<std::string>("analysisMode");
    if(!am.compare("HNUMU")) analysisMode_ = HeavyNuEvent::HNUMU;
    else if(!am.compare("HNUE")) analysisMode_ = HeavyNuEvent::HNUE;
    else if(!am.compare("TOP")) analysisMode_ = HeavyNuEvent::TOP;
    else if(!am.compare("QCD")) analysisMode_ = HeavyNuEvent::QCD;
    else if(!am.compare("CLO")) analysisMode_ = HeavyNuEvent::CLO;
    else if(!am.compare("TAUX")) analysisMode_ = HeavyNuEvent::TAUX;
    else std::cout << "!!!!!!!!INVALID ANALYSIS MODE : " << am << " !!!!!!!!\noptions are: HNUMU, HNUE, TOP, QCD, CLO" << std::endl;

    switch(analysisMode_)
    {
        case HeavyNuEvent::HNUMU:
            pdgid_ = 9900014;
            break;
        case HeavyNuEvent::HNUE:
            pdgid_ = 9900012;
            break;
        case HeavyNuEvent::TAUX:
            pdgid_ = 9900016;
            break;
        case HeavyNuEvent::TOP:
        case HeavyNuEvent::QCD:
        case HeavyNuEvent::CLO:
            break;
    }

    nDirtyCands_ = iConfig.getUntrackedParameter<int>("nFakeLeptons", 0);
    muonsInJets_ = iConfig.getUntrackedParameter<bool>("qcdMuonsInJets", false);

    doPDFreweight_ = iConfig.getUntrackedParameter<bool>("doPDFReweight", false);
    if (doPDFreweight_)
    {
        pdfReweightBaseName = iConfig.getUntrackedParameter<std::string>("pdfReweightBaseName");
        pdfReweightTargetName = iConfig.getUntrackedParameter<std::string>("pdfReweightTargetName");
        pdfReweightBaseId = iConfig.getUntrackedParameter<int>("pdfReweightBaseId", 0);
        pdfReweightTargetId = iConfig.getUntrackedParameter<int>("pdfReweightTargetId", 0);
        std::cout << "PDF Reweighting from " << pdfReweightBaseName << ":" << pdfReweightBaseId
                << " to " << pdfReweightTargetName << ":" << pdfReweightTargetId
                << std::endl;
    }

    // ==================== Init other members ====================
    //

    //    nnif_ = new HeavyNu_NNIF(iConfig);
    trig_ = new HeavyNuTrigger(iConfig.getParameter<edm::ParameterSet > ("trigMatchPset"));
    muid_ = new HeavyNuID(iConfig.getParameter<edm::ParameterSet > ("muIDPset"));
    // ==================== Book the histos ====================
    //
    edm::Service<TFileService> fs;

    hists.mc_type = fs->make<TH1D > ("mc_type", "MC Type Code", 100, -0.5, 99.5);

    //this must be after at least 1 call of fs->make<...>(...) in order for the directory to have been created.
    if(addSlopeTree_) hnuTree_ = new HeavyNuTree(*fs->getBareDirectory(), true);

    if ( studyRatePerRun_ )
    {
        hists.z2jetPerRun = fs->make<TH1I > ("z2jetPerRun", "M3 Z #to ee/#mu#mu events per run", 20000, 190000, 210000) ;
    }
    hists.nelec = fs->make<TH1D > ("nelec", "N(e^{#pm})", 10, -0.5, 9.5);
    hists.nmuAll = fs->make<TH1D > ("nmuAll", "N(#mu^{#pm})", 10, -0.5, 9.5);
    hists.nmuLoose = fs->make<TH1D > ("nmuLoose", "N(#mu^{#pm}) passes Loose", 10, -0.5, 9.5);
    hists.nmuTight = fs->make<TH1D > ("nmuTight", "N(#mu^{#pm}) passes Tight", 10, -0.5, 9.5);
    hists.cutlevel = fs->make<TH1D > ("cutlevel", "Cut Level", 11, -1.5, 9.5);
    hists.cutlevel->GetXaxis()->SetBinLabel(1, "Raw");
    hists.cutlevel->GetXaxis()->SetBinLabel(2, "No cuts");
    hists.cutlevel->GetXaxis()->SetBinLabel(3, "mmjj p_{T}");
    hists.cutlevel->GetXaxis()->SetBinLabel(4, "M1 (Trigger)");
    hists.cutlevel->GetXaxis()->SetBinLabel(5, "M2 (Vertex)");
    hists.cutlevel->GetXaxis()->SetBinLabel(6, "M3 (high p_{T})");
    hists.cutlevel->GetXaxis()->SetBinLabel(7, "M4 (M_{#mu#mu})");
    hists.cutlevel->GetXaxis()->SetBinLabel(8, "M5 (M(W_{R})");
    hists.weights = fs->make<TH1D > ("weights", "event weights", 1000, 0.0, 10.0);
    hists.njet = fs->make<TH1D > ("njet", "N(Jet)", 50, -0.5, 49.5);
    hists.nmet = fs->make<TH1D > ("nmet", "N(MET)", 50, -0.5, 49.5);
    hists.met = fs->make<TH1D > ("met", "MET distribution", 100, 0, 2000);
    hists.muPt = fs->make<TH1D > ("muPt", "#mu p_{T} distribution", 100, 0, 2000);
    hists.muEta = fs->make<TH1D > ("muEta", "#mu #eta distribution", 50, -2.5, 2.5);
    hists.muPhi = fs->make<TH1D > ("muPhi", "#mu #phi distribution", 60, -3.14159, 3.14159);
    hists.jetPt = fs->make<TH1D > ("jetPt", "jet p_{T} distribution", 100, 0, 2000);
    hists.jetEta = fs->make<TH1D > ("jetEta", "jet #eta distribution", 50, -5, 5);
    hists.jetPhi = fs->make<TH1D > ("jetPhi", "jet #phi distribution", 60, -3.14159, 3.14159);
    hists.jetID = fs->make<TH1I > ("jetID", "Jet ID", 3, 0, 3);
    hists.jetPtvsNum = fs->make<TH2D > ("jetPtvsNum", "Jet P_{T} vs. Jet # ", 11, -0.5, 10.5, 200, 0., 2000.);
    hists.dVzMuJets = fs->make<TH2D > ("dVzMuJets", "vertex dz Mu Jets vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);
    hists.dVzMuMus = fs->make<TH2D > ("dVzMuMus", "vertex dz Mu Mus vs npue", 500, -5.0, 5.0, 50, -0.5, 49.5);

    hists.jetID->GetXaxis()->SetBinLabel(1, "Neither");
    hists.jetID->GetXaxis()->SetBinLabel(2, "Loose");
    hists.jetID->GetXaxis()->SetBinLabel(3, "Tight");

    hists.trkIsoStudy = fs->make<TH1D > ("trkIsoStudy", "Tracker Relative Isolation", 100, 0., 1.);
    hists.closejetMu1tagMu2probeInZwin = fs->make<TH1D > ("closejetMu1tagMu2probeInZwin", "Probe (#mu_{2}) p_{T}", 100, 0., 1000.);
    hists.closejetMu1tagMu2passInZwin = fs->make<TH1D > ("closejetMu1tagMu2passInZwin", "Passing probe (#mu_{2}) p_{T}", 100, 0., 1000.);
    hists.closejetMu2tagMu1probeInZwin = fs->make<TH1D > ("closejetMu2tagMu1probeInZwin", "Probe (#mu_{1}) p_{T}", 100, 0., 1000.);
    hists.closejetMu2tagMu1passInZwin = fs->make<TH1D > ("closejetMu2tagMu1passInZwin", "Passing probe (#mu_{1}) p_{T}", 100, 0., 1000.);

    if(applyMESfactor_ == 1.0)
    { // otherwise don't bother

        // Loose/Tight vs Pt
        hists.looseMuPt = fs->make<TH1D > ("looseMuPt", "#mu p_{T}, passes Loose", 100, 0, 2000);
        hists.tightMuPt = fs->make<TH1D > ("tightMuPt", "#mu p_{T}, passes Tight", 100, 0, 2000);

        //  Muon quality variables vs muon p_T
        //
        hists.muNvalidHitsVsPt = fs->make<TH2D > ("muNvalidHitsVsPt",
                "#mu # Valid hits vs p_{T}; p_{T}(#mu) (GeV); # Tracker+Pixel Hits",
                150, 0, 3000, 50, 0, 50);
        hists.mudBvsPt = fs->make<TH2D > ("mudBvsPt",
                "#mu dXY vs p_{T}; p_{T}(#mu) (GeV); dXY(#mu)",
                150, 0, 3000, 40, 0., 20.);
        hists.muNormChi2vsPt = fs->make<TH2D > ("muNormChi2vsPt",
                "#mu Norm #chi^{2} vs p_{T}; p_{T}(#mu) (GeV); norm #chi^{2}",
                150, 0, 3000, 25, 0, 50);
        hists.muQualVsPt = fs->make<TH2D > ("muIDvsPt",
                "Qual(#mu) vs p_{T}(#mu); p_{T}(#mu) (GeV)",
                150, 0, 3000, 4, 0, 4);
        hists.muNmatchesVsPt = fs->make<TH2D > ("muNmatchesVsPt",
                "#mu # Matches vs p_{T}; p_{T}(#mu) (GeV); # Matches",
                150, 0, 3000, 10, 0, 10);
        hists.muNvalidMuonHitsVsPt = fs->make<TH2D > ("muNvalidMuonHitsVsPt",
                "# Valid Muon hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Muon Hits",
                150, 0, 3000, 50, 0, 50);
        hists.muNvalidPixelHitsVsPt = fs->make<TH2D > ("muNvalidPixelHitsVsPt",
                "# Valid Pixel hits vs p_{T}; p_{T}(#mu) (GeV); # Valid Pixel Hits",
                150, 0, 3000, 10, 0, 10);

        labelMuonQualAxis(hists.muQualVsPt->GetXaxis());

        hists.muTrckIsoVsPt = fs->make<TH2D > ("muTrckIsoVsPt", "trackIso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);trackIso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muHcalIsoVsPt = fs->make<TH2D > ("muHcalIsoVsPt", "HCAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);HCAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muEcalIsoVsPt = fs->make<TH2D > ("muEcalIsoVsPt", "ECAL Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);ECAL Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
        hists.muCaloIsoVsPt = fs->make<TH2D > ("muCaloIsoVsPt", "Calo Iso(#mu) vs p_{T}(#mu);p_{T}(#mu)(GeV);Calo Iso (GeV)", 150, 0., 3000., 30, 0., 300.);
    }

    // Histos per cut:
    //
    // hists.twoL->book(new TFileDirectory(fs->mkdir("cutm1_LL")), "(two lepton:-1)");
    if(analysisMode_ == HeavyNuEvent::HNUMU)
    {
        hists.noCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.JJptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0a_JJpt")), "(2 jets with ptcuts:0a)");
        hists.LLptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0b_LLpt")), "(2 leptons with ptcuts:0b)");
        hists.LLJJptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(L1 High pt cut:4)");
        //hists.Mu1HighPtCut_1bjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4a_L1HighPt_1b")), "(L1 High pt cut:4a)");
        //hists.Mu1HighPtCut_2bjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4b_L1HighPt_2b")), "(L1 High pt cut:4b)");
        hists.ZRegionCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4c_ZPeak")), "(Z-Peak:4c)");
        hists.diLmassCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(mumu mass cut:5)");
        hists.mWRmassCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(mumujj mass cut:6)");
        hists.OneJetTrigHighPt = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4clo_TrigHighPt")), "(3objects with trigger and pt cuts:4clo)");
        hists.OneJetdiLmassCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5clo_diLmass")), "(mumu mass cut:5clo)");

        if ( studyAlternativeSelection_ )
        {
            hists.AlternativeElecChanPt = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu1Pt80Mu2Pt40")), "(Alternative: Electron pT selection)");
            hists.AlternativeMu1Pt40 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu1Pt40")), "(Alternative: mu1 pT 40 GeV)");
            hists.AlternativeMu2Pt40 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu2Pt40")), "(Alternative: mu2 pT 40 GeV)");
            hists.AlternativeMu2Pt60 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltMu2Pt60")), "(Alternative: mu2 pT 60 GeV)");
            hists.AlternativeJetPt60 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltJetPt60")), "(Alternative: jet pT 60 GeV)");
            hists.AlternativeBarrelLoose = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltBarrelLoose")), "(Alternative: one barrel muon, loose)");
            hists.AlternativeBarrelTight = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltBarrelTight")), "(Alternative: one barrel muon, tight)");
            hists.AlternativeAtLeastOneBjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltOneBjet")), "(Alternative: at least one b-jet)");
            hists.AlternativeTwoBjets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltTwoBjets")), "(Alternative: two b-jets)");
            hists.AlternativeDimuonMass120 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("AltDimuon120")), "(Alternative: dimuon mass 120 GeV)");
        }

        if(studyScaleFactorEvolution_)
        {
            hists.Mu1Pt40GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt40GeV")), "(Mu1 40 GeV pt cut)");
            hists.Mu1Pt50GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt50GeV")), "(Mu1 50 GeV pt cut)");
            hists.Mu1Pt60GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt60GeV")), "(Mu1 60 GeV pt cut)");
            hists.Mu1Pt80GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt80GeV")), "(Mu1 80 GeV pt cut)");
            hists.Mu1Pt100GeVCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Pt100GeV")), "(Mu1 100 GeV pt cut)");

            hists.Mu1HighPtCutVtxEq1 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtxEq1")), "(Mu1 60 GeV pt cut, 1 vtx)");
            hists.Mu1HighPtCutVtx2to5 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtx2to5")), "(Mu1 60 GeV pt cut, 2-5 vtx)");
            hists.Mu1HighPtCutVtxGt5 = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtVtxGt5")), "(Mu1 60 GeV pt cut, 6+ vtx)");

            hists.Mu1HighPtCutNoJets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPtNoJets")), "(Mu1 60 GeV pt cut, no jets)");
            hists.Mu1HighPtCut1Jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1HighPt1Jet")), "(Mu1 60 GeV pt cut, 1 jet)");
        }

        if(trig_->matchingEnabled())
        {
            TFileDirectory *tdirm1 = new TFileDirectory(fs->mkdir("Mu1TrigMatchesInZwin"));
            TFileDirectory *tdirm2 = new TFileDirectory(fs->mkdir("Mu2TrigMatchesInZwin"));
            hists.Mu1TrigMatchesInZwin = new HeavyNuMuHist(tdirm1, "(#mu1 trigger match in Z mass Window)");
            hists.Mu2TrigMatchesInZwin = new HeavyNuMuHist(tdirm2 , "(#mu2 Trigger match in Z mass Window)");
            // hists.Mu1Mu2TrigMatchesInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("Mu1Mu2TrigMatchesInZwin")), "(#mu1,#mu2 Trigger match in Z mass Window)");
            trig_->book(*tdirm1, &(((HeavyNuMuHist*)hists.Mu1TrigMatchesInZwin)->trigHistos));
            trig_->book(*tdirm2, &(((HeavyNuMuHist*)hists.Mu2TrigMatchesInZwin)->trigHistos));
        }

        if(studyMuonSelectionEff_)
        {
            hists.TightTagTrackProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin0jets")), "(probe_{0}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin1jet")), "(probe_{1}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbeInZwin2jets")), "(probe_{2}, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin0jets")), "(probe_{0} passes, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin1jet")), "(probe_{1} passes, ID in Z mass Window)", 2);
            hists.TightTagTrackProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrackProbePassesInZwin2jets")), "(probe_{2} passes, ID in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin0jets")), "(probe_{0}, iso in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin1jet")), "(probe_{1}, iso in Z mass Window)", 2);
            hists.TightTagTightProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbeInZwin2jets")), "(probe_{2}, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin0jets")), "(probe_{0} passes, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin1jet")), "(probe_{1} passes, iso in Z mass Window)", 2);
            hists.TightTagTightProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightProbePassesInZwin2jets")), "(probe_{2} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin0jets")), "(CJ probe_{0}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin1jet")), "(CJ probe_{1}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbeInZwin2jets")), "(CJ probe_{2}, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin0jets")), "(CJ probe_{0} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin1jet")), "(CJ probe_{1} passes, iso in Z mass Window)", 2);
            hists.TightTagTightCJProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTightCJProbePassesInZwin2jets")), "(CJ probe_{2} passes, iso in Z mass Window)", 2);
            hists.TightTagTrigProbeInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrigProbeInZwin")), "(Trig probe, in Z mass Window)", 2);
            hists.TightTagTrigProbePassesInZwin = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("TightTagTrigProbePassesInZwin")), "(Trig probe passes, in Z mass Window)", 2);
        }

    }
    else if(analysisMode_ == HeavyNuEvent::HNUE)
    {
        hists.noCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.JJptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0a_JJpt")), "(2 jets with ptcuts:0a)");
        hists.LLptCuts = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut0b_LLpt")), "(2 leptons with ptcuts:0b)");
        hists.LLJJptCuts =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches =  new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(L1 High pt cut:4)");
        //hists.Mu1HighPtCut_1bjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4a_L1HighPt_1b")), "(L1 High pt cut:4a)");
        //hists.Mu1HighPtCut_2bjet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4b_L1HighPt_2b")), "(L1 High pt cut:4b)");
        hists.ZRegionCut = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("cut4c_ZPeak")), "(Z-Peak:4c)");
        hists.diLmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(ee mass cut:5)");
        hists.mWRmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(eejj mass cut:6)");
        hists.OneJetTrigHighPt = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4clo_TrigHighPt")), "(3objects with trigger and pt cuts:4clo)");
        hists.OneJetdiLmassCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5clo_diLmass")), "(ee mass cut:5clo)");

        if(studyMuonSelectionEff_)
        {
            hists.HeepTagGsfProbeInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbeInZwin0jets")), "(probe_{0}, ID in Z mass Window)", 2);
            hists.HeepTagGsfProbeInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbeInZwin1jet")), "(probe_{1}, ID in Z mass Window)", 2);
            hists.HeepTagGsfProbeInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbeInZwin2jets")), "(probe_{2}, ID in Z mass Window)", 2);
            hists.HeepTagGsfProbePassesInZwin0jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbePassesInZwin0jets")), "(probe_{0} passes, ID in Z mass Window)", 2);
            hists.HeepTagGsfProbePassesInZwin1jet = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbePassesInZwin1jet")), "(probe_{1} passes, ID in Z mass Window)", 2);
            hists.HeepTagGsfProbePassesInZwin2jets = new HeavyNuMuHist(new TFileDirectory(fs->mkdir("HeepTagGsfProbePassesInZwin2jets")), "(probe_{2} passes, ID in Z mass Window)", 2);
        }
    }
    else if(analysisMode_ == HeavyNuEvent::TOP)
    {
        hists.noCuts = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.JJptCuts = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut0a_JJpt")), "(2 jets with ptcuts:0a)");
        hists.LLptCuts = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut0b_LLpt")), "(2 leptons with ptcuts:0b)");
        hists.LLJJptCuts =   new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches =  new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(L1 High pt cut:4)");
        //hists.Mu1HighPtCut_1bjet = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut4a_L1HighPt_1b")), "(L1 High pt cut:4a)");
        //hists.Mu1HighPtCut_2bjet = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut4b_L1HighPt_2b")), "(L1 High pt cut:4b)");
        hists.ZRegionCut = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut4c_ZPeak")), "(Z-Peak:4c)");
        hists.diLmassCut =   new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(emu mass cut:5)");
        hists.mWRmassCut =   new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(emujj mass cut:6)");
        hists.OneJetTrigHighPt = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut4clo_TrigHighPt")), "(3objects with trigger and pt cuts:4clo)");
        hists.OneJetdiLmassCut = new HeavyNuTopHist(new TFileDirectory(fs->mkdir("cut5clo_diLmass")), "(emu mass cut:5clo)");
    }
    else if(analysisMode_ == HeavyNuEvent::TAUX)
    {
        hists.noCuts = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)");
        hists.JJptCuts = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut0a_JJpt")), "(2 jets with ptcuts:0a)");
        hists.LLptCuts = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut0b_LLpt")), "(2 leptons with ptcuts:0b)");
        hists.LLJJptCuts =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)");
        hists.TrigMatches =  new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut2_TrigMatches")), "(Trigger match:2)");
        hists.Mu1HighPtCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4_L1HighPt")), "(L1 High pt cut:4)");
        //hists.Mu1HighPtCut_1bjet = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4a_L1HighPt_1b")), "(L1 High pt cut:4a)");
        //hists.Mu1HighPtCut_2bjet = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4b_L1HighPt_2b")), "(L1 High pt cut:4b)");
        hists.ZRegionCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4c_ZPeak")), "(Z-Peak:4c)");
        hists.diLmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(ee mass cut:5)");
        hists.mWRmassCut =   new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(eejj mass cut:6)");
        hists.OneJetTrigHighPt = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut4clo_TrigHighPt")), "(3objects with trigger and pt cuts:4clo)");
        hists.OneJetdiLmassCut = new HeavyNuHistSet(new TFileDirectory(fs->mkdir("cut5clo_diLmass")), "(ll mass cut:5clo)");
    }

    hists.rundir = new TFileDirectory(fs->mkdir("RunDir"));

    init_ = false;

    if(applyJECUsign_ > 0)
    {
        hists.jecUncHi = fs->make<TH1F > ("jecUncHi", "JEC Uncertainty (high); Uncertainty (%)", 100, 0, 10);
        hists.jecUncHiVsEtaPt = fs->make<TProfile2D > ("jecUncHiVsEtaPt",
                "JEC Uncertainty (high)(%);Jet #eta;Jet p_{T} (GeV)",
                50, -2.5, 2.5, 100, 0, 1000);
    }
    else if(applyJECUsign_ < 0)
    {
        hists.jecUncLo = fs->make<TH1F > ("jecUncLo", "JEC Uncertainty (low);  Uncertainty (%)", 100, 0, 10);
        hists.jecUncLoVsEtaPt = fs->make<TProfile2D > ("jecUncLoVsEtaPt",
                "JEC Uncertainty (low)(%);Jet #eta;Jet p_{T} (GeV)",
                50, -2.5, 2.5, 100, 0, 1000);
    }

    MCweightByVertex_ = edm::LumiReWeighting(hnu::get_standard_pileup_mc(pileupEra_),
                                             hnu::get_standard_pileup_data(pileupEra_, puShift_));

    // For the record...
    std::cout << "Configurable cut values applied:" << std::endl;
    std::cout << "Analysis Mode     = " << analysisMode_ << std::endl;
    std::cout << "muonTag           = " << muonTag_ << std::endl;
    std::cout << "jetTag            = " << jetTag_ << std::endl;
    std::cout << "metTag            = " << metTag_ << std::endl;
    std::cout << "electronTag       = " << elecTag_ << std::endl;
    std::cout << "heepVersion       = " << heepVersion_ << std::endl;
    std::cout << "trackTag          = " << trackTag_ << std::endl;
    std::cout << "btagName          = " << btagName << std::endl;
    std::cout << "minBtagDiscr      = " << minBtagDiscVal << std::endl;
    std::cout << "ZmassWinMinGeV    = " << effinfo_.ZwinMinGeV << " GeV" << std::endl;
    std::cout << "ZmassWinMaxGeV    = " << effinfo_.ZwinMaxGeV << " GeV" << std::endl;
    std::cout << "minMu1pt          = " << cuts.minimum_mu1_pt << " GeV" << std::endl;
    std::cout << "minMu2pt          = " << cuts.minimum_mu2_pt << " GeV" << std::endl;
    std::cout << "minJetPt          = " << cuts.minimum_jet_pt << " GeV" << std::endl;
    std::cout << "maxElecAbsEta     = " << cuts.maximum_elec_abseta << std::endl;
    std::cout << "maxMuAbsEta       = " << cuts.maximum_mu_abseta << std::endl;
    std::cout << "maxJetAbsEta      = " << cuts.maximum_jet_abseta << std::endl;
    std::cout << "minMuonJetdR      = " << cuts.minimum_muon_jet_dR << std::endl;
    std::cout << "muonTrackRelIso   = " << cuts.muon_trackiso_limit << std::endl;
    std::cout << "minMuMuMass       = " << cuts.minimum_mumu_mass << " GeV" << std::endl;
    std::cout << "min4objMass       = " << cuts.minimum_mWR_mass << " GeV" << std::endl;
    std::cout << "jecEra            = " << jecVal_ << std::endl;
    std::cout << "applyJECUsign     = " << applyJECUsign_ << std::endl;
    std::cout << "applyMESfactor    = " << applyMESfactor_ << std::endl;
    std::cout << "Mu Scale Unc.     = " << merUncertainty_ << std::endl;
    std::cout << "Correct E Scale   = " << correctEscale_ << std::endl;
    std::cout << "applyMuIDEffcorr  = " << applyMuIDCorrections_ << std::endl;
    std::cout << "applyMuIDEffsign  = " << applyMuIDEffsign_ << std::endl;
    std::cout << "applyTrigEffsign  = " << applyTrigEffsign_ << std::endl;
    std::cout << "studyMuSelectEff  = " << studyMuonSelectionEff_ << std::endl;
    std::cout << "studyScaleFactor  = " << studyScaleFactorEvolution_ << std::endl;
    std::cout << "pileup era        = " << pileupEra_ << std::endl;
    std::cout << "DisableTriggerCorrection = " << disableTriggerCorrection_ << std::endl;
    std::cout << "isPFJets          = " << isPFJets_ << std::endl;
    std::cout << "addSlopeTree      = " << addSlopeTree_ << std::endl;
    if ( studyRatePerRun_ )
    {
        std::cout << "Study rate of Z+2 jets events for run range " << hists.z2jetPerRun->GetXaxis()->GetBinLowEdge(1)
                << " to " << hists.z2jetPerRun->GetXaxis()->GetBinUpEdge(hists.z2jetPerRun->GetNbinsX()) << std::endl ;
    }

}

HeavyNu::~HeavyNu(){

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

void HeavyNu::fillBasicMuHistos(const pat::Muon& m, edm::Handle<reco::VertexCollection> pvHandle)
{
    double mupt = m.pt();
    hists.muPt->Fill(applyMESfactor_ * mupt);
    hists.muEta->Fill(m.eta());
    hists.muPhi->Fill(m.phi());

    if(applyMESfactor_ == 1.0)
    {
        hists.mudBvsPt->Fill(mupt, m.dB());

        if(hnu::isLooseMuon(m) && m.globalTrack().isNonnull())
        {
            hists.looseMuPt ->Fill(mupt);
            hists.muNvalidHitsVsPt ->Fill(mupt, m.numberOfValidHits());
            hists.muNormChi2vsPt ->Fill(mupt, m.normChi2());
            hists.muNmatchesVsPt ->Fill(mupt, m.numberOfMatches());

            reco::TrackRef gt = m.globalTrack();
            hists.muNvalidMuonHitsVsPt ->Fill(mupt, gt->hitPattern().numberOfValidMuonHits());
            hists.muNvalidPixelHitsVsPt->Fill(mupt, gt->hitPattern().numberOfValidPixelHits());

            if(hnu::isTightHighPtMuon(m, pvHandle)) hists.tightMuPt->Fill(mupt);
        }
        hists.muQualVsPt->Fill(mupt, 0);
        for(int i = 1; i < muonQualityFlags; i++)
            if(m.muonID(muonQuality[i]))
                hists.muQualVsPt->Fill(mupt, i);

        hists.muTrckIsoVsPt->Fill(mupt, m.trackIso());
        hists.muHcalIsoVsPt->Fill(mupt, m.hcalIso());
        hists.muEcalIsoVsPt->Fill(mupt, m.ecalIso());
        hists.muCaloIsoVsPt->Fill(mupt, m.caloIso());
    }
} // HeavyNu::fillBasicMuHistos

void HeavyNu::fillBasicJetHistos(const pat::Jet& j, int jetnum)
{
    double jpt = j.pt(), jeta = j.eta();
    float totalunc = 0.0f;

    int jpdgId = 0;
    if(j.genParton()) jpdgId = j.genParton()->pdgId();
    bool isBjet = (abs(jpdgId) == 5);

    float jecuHi = 0.;
    float jecuLo = 0.;

    if(applyJECUsign_)
    {
        jecuHi = hnu::jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, true);
        jecuLo = hnu::jecTotalUncertainty(jpt, jeta, jecuObj_, jecVal_, isBjet, false);
        totalunc = (applyJECUsign_ > 0)?jecuHi:jecuLo;
        jpt *= (1.0 + (applyJECUsign_ * totalunc));
    }

    hists.jetPt ->Fill(jpt);
    hists.jetEta->Fill(jeta);
    hists.jetPhi->Fill(j.phi());
    if(isPFJets_) hists.jetID ->Fill(hnu::jetID(j));
    hists.jetPtvsNum->Fill(jetnum, jpt);

    jecuHi *= 100.f;
    jecuLo *= 100.f;
    if(jecuHi != jecuLo) std::cout << jeta << "\t" << jpt << "\t" << jecuHi << "\t" << jecuLo << std::endl;

    if(applyJECUsign_ > 0)
    {
        hists.jecUncHi->Fill(jecuHi);
        hists.jecUncHiVsEtaPt->Fill(jeta, j.pt(), (double)jecuHi);
    }
    else if(applyJECUsign_ < 0)
    {
        hists.jecUncLo->Fill(jecuLo);
        hists.jecUncLoVsEtaPt->Fill(jeta, j.pt(), (double)jecuLo);
    }
} // HeavyNu::fillBasicJetHistos

TH1 * HeavyNu::bookRunHisto(uint32_t runNumber)
{
    std::string runstr = int2str<uint32_t > (runNumber);
    return hists.rundir->make <TH1I > (runstr.c_str(), runstr.c_str(), 1, 1, 2);
}

void HeavyNu::studyJetVertex(edm::Handle<pat::JetCollection>& pJets,
                             edm::Handle<reco::JPTJetCollection>& jptJets,
                             edm::Handle<pat::MuonCollection>& pMuons, int npue)
{
    if(pMuons->size() < 2) return;

    const pat::Muon& m1 = pMuons->at(0);
    const pat::Muon& m2 = pMuons->at(1);

    double mLL = (m1.p4() + m2.p4()).M();
    double dzll = m1.vz() - m2.vz();

    hists.dVzMuMus->Fill(dzll, npue);

    if(mLL < 71 || mLL > 111 || fabs(dzll) > 0.1) return;

    for(size_t iJet = 0; iJet < pJets->size(); iJet++)
    {
        pat::JetRef iJ = pat::JetRef(pJets, iJet);

        if(hnu::jetID(*iJ) > 1 && iJ->pt() > 40)
        {
            double jvz = hnu::caloJetVertex(*iJ, *jptJets);
            hists.dVzMuJets->Fill(m1.vz() - jvz, npue);
        }
    }
}

void HeavyNu::studyJetMatching(HeavyNuEvent& hnuEvent, edm::Handle<std::vector<reco::GenJet> > genjets)
{
    if(hnuEvent.nJets < 2) return;  //safty against to few jets.  

    //here we match the reco jets to gen jets
    double mindR1 = 100, mindR2 = 100;
    for(std::vector<reco::GenJet>::const_iterator igj = genjets->begin(); igj != genjets->end(); ++igj)
    {
        double dR1 = deltaR(hnuEvent.j1.p4(), igj->p4());
        double dR2 = deltaR(hnuEvent.j2.p4(), igj->p4());

        if(dR1 < mindR1)
        {
            mindR1 = dR1;
            hnuEvent.gj1 = *igj;
        }
        if(dR2 < mindR2)
        {
            mindR2 = dR2;
            hnuEvent.gj2 = *igj;
        }
    }

    //after finding the matching gen jets we track their parentage and try to match them to a Nu_L
    bool gmj1 = false, gmj2 = false;
    std::vector<const reco::GenParticle*> mothers = hnuEvent.gj1.getGenConstituents();
    for(std::vector<const reco::GenParticle*>::const_iterator iM = mothers.begin(); iM != mothers.end(); ++iM)
    {
        if((gmj1 |= isWrDaughter(*iM, pdgid_))) break;
    }
    mothers = hnuEvent.gj2.getGenConstituents();
    for(std::vector<const reco::GenParticle*>::const_iterator iM = mothers.begin(); iM != mothers.end(); ++iM)
    {
        if((gmj2 |= isWrDaughter(*iM, pdgid_))) break;
    }
    hnuEvent.numNuLJetsMatched = (int)gmj1 + (int)gmj2;
}

void HeavyNu::selectMuons(std::vector<pat::Muon>& muCands,
                          std::vector< std::pair<pat::Jet, float> >& jetCands,
                          edm::Handle<reco::VertexCollection> pvHandle,
                          HeavyNuEvent& hnuEvent)
{
    unsigned int dirtyPosition = muCands.size() + 1 ;
    unsigned int cleanPosition = muCands.size() + 1 ;
    for(unsigned int i = 0; i < muCands.size(); i++)
    {
        if(hnuEvent.nLeptons == 2) break;
        pat::Muon iM = muCands.at(i);
        if ( nDirtyCands_ > 0 )
        {
            double caloIso = (iM.ecalIso() + iM.hcalIso()) / iM.pt() ;
            bool isDirty = (hnu::muIsolation(iM) > cuts.muon_trackiso_limit) ;
            if (muid_->idEra() < 0) isDirty = isDirty && !hnu::isTightHighPtMuon(iM, pvHandle) ;
            if ( muonsInJets_ && isDirty )
            {
                bool overlap = false ;
                for (unsigned int j = 0; j < jetCands.size(); j++)
                {
                    pat::Jet iJ = jetCands.at(j).first ;
                    double dRjm = deltaR(iM.eta(), iM.phi(), iJ.eta(), iJ.phi());
                    if ( dRjm < cuts.qcd_jet_muon_dRoverlap )
                    {
                        overlap = true ;
                        break ;
                    }
                }
                isDirty = isDirty && overlap ;
            }
            else
            {
                isDirty = isDirty && ( caloIso > 0.10 ) ;
            }
            double jetSeparation = (muonsInJets_ ? -1.0 : cuts.minimum_muon_jet_dR) ;

            if ( nDirtyCands_ == 2 )
            {
                if ( isDirty )
                {
                    hnu::addMuon(iM, hnuEvent, jetSeparation, nDirtyCands_);
                }
            }
                // Muon candidate list already sorted, so dirty/clean filled in proper order
            else if ( nDirtyCands_ == 1 )
            {
                if ( isDirty && dirtyPosition > muCands.size() )
                {
                    if(hnu::addMuon(iM, hnuEvent, jetSeparation, nDirtyCands_)) dirtyPosition = i;
                }
                else if ( !isDirty && cleanPosition > muCands.size() && hnu::muIsolation(iM) < cuts.muon_trackiso_limit)
                {
                    if (muid_->idEra() > 0 || (hnu::isTightHighPtMuon(iM, pvHandle)))
                        if (hnu::addMuon(iM, hnuEvent, jetSeparation, 0)) cleanPosition = i;
                }
            }
        }
        else if(hnu::muIsolation(iM) < cuts.muon_trackiso_limit)
        {
            if (muid_->idEra() > 0 || (hnu::isTightHighPtMuon(iM, pvHandle)))
                hnu::addMuon(iM, hnuEvent, cuts.minimum_muon_jet_dR, 0);
        }
    }

    if ( nDirtyCands_ > 0 && hnuEvent.nLeptons > 1 )
    {
        // "Fake" muons must fail nominal isolation requirements...correct for this
        if ( nDirtyCands_ == 2 )
        {
            hnuEvent.eventWgt *= (hnu::fakeProbability(hnuEvent.mu1) / (1.0 - hnu::fakeProbability(hnuEvent.mu1))) ;
            hnuEvent.eventWgt *= (hnu::fakeProbability(hnuEvent.mu2) / (1.0 - hnu::fakeProbability(hnuEvent.mu2))) ;
        }
        if ( nDirtyCands_ == 1 )
            hnuEvent.eventWgt *= (hnu::fakeProbability(muCands.at(dirtyPosition)) / (1.0 - hnu::fakeProbability(muCands.at(dirtyPosition)))) ;
    }

    if(applyMuIDCorrections_ && hnuEvent.isMC)
    {
        //double mu1wgt = (hnuEvent.nLeptons > 0)? (muid_->weightForMC((hnuEvent.mu1.pt()), applyMuIDEffsign_)):1.0;
        //double mu2wgt = (hnuEvent.nLeptons > 1)? (muid_->weightForMC((hnuEvent.mu2.pt()), applyMuIDEffsign_)):1.0;

        double mu1wgt = (hnuEvent.nLeptons > 0)? (muid_->weightForMCbyEta((hnuEvent.mu1.eta()), applyMuIDEffsign_)):1.0;
        double mu2wgt = (hnuEvent.nLeptons > 1)? (muid_->weightForMCbyEta((hnuEvent.mu2.eta()), applyMuIDEffsign_)):1.0;

        hnuEvent.eventWgt *= (mu1wgt * mu2wgt);
    }
}

void HeavyNu::selectElectrons(std::vector< std::pair<pat::Electron, float> >& eCands, edm::Handle<pat::ElectronCollection>& pElecs, edm::Handle<reco::VertexCollection> pvHandle, HeavyNuEvent& hnuEvent)
{
    // std::cout << "DEBUG: selectElectrons with " << eCands.size() << " candidates, need " << nDirtyCands_ << " dirty" << std::endl ; 
    unsigned int dirtyPosition = eCands.size() + 1 ;
    unsigned int cleanPosition = eCands.size() + 1 ;
    std::vector< std::pair<pat::Electron, float> > eDirtyCands = ( nDirtyCands_ > 0 ) ?
            hnu::getElectronList(pElecs, cuts.maximum_elec_abseta, cuts.minimum_mu2_pt, cuts.minimum_mu2_pt,
                                 (-1 * heepVersion_), pvHandle, elecRho_) : eCands ;
    if ( nDirtyCands_ == 0 )
    {
        for(unsigned int i = 0; i < eCands.size(); i++)
        {
            if(hnuEvent.nLeptons == 2) break;
            pat::Electron iE = eCands.at(i).first;
            hnu::addElectron(iE, hnuEvent, cuts.minimum_muon_jet_dR);
        }
    }
    else
    {
        // std::cout << "DEBUG: Looking for " << nDirtyCands_ << " dirty electrons" << std::endl ; 
        // std::cout << "DEBUG: There are " << eDirtyCands.size() << " dirty electrons" << std::endl ; 
        dirtyPosition = eDirtyCands.size() + 1 ;
        if ( nDirtyCands_ == 2 )
        {
            for(unsigned int i = 0; i < eDirtyCands.size(); i++)
            {
                // std::cout << "DEBUG: Looking at dirty electron " << i+1 << " of " << eDirtyCands.size() << std::endl ; 
                if(hnuEvent.nLeptons == 2) break;
                pat::Electron iE = eDirtyCands.at(i).first;
                // std::cout << "DEBUG: electron has pt " << hnu::getElectronEt(iE,false) << std::endl ; 
                hnu::addElectron(iE, hnuEvent, cuts.minimum_muon_jet_dR);
                // std::cout << "DEBUG: Now have " << hnuEvent.nLeptons << " leptons" << std::endl ; 
            }
            if ( hnuEvent.nLeptons == 2 )
            {
                // Need to correct for the fact that "fake" electrons are required to fail nominal ID
                hnuEvent.eventWgt *= ( hnu::fakeProbability(hnuEvent.e1) / ( 1.0 - hnu::fakeProbability(hnuEvent.e1) ) ) ;
                hnuEvent.eventWgt *= ( hnu::fakeProbability(hnuEvent.e2) / ( 1.0 - hnu::fakeProbability(hnuEvent.e2) ) ) ;
            }
        }
        else if ( nDirtyCands_ == 1 )
        {
            // Dirty electrons are in a separate collection than clean electrons
            // addElectron does the sorting automatically, and takes care of counting electrons
            for (unsigned int i = 0; i < eDirtyCands.size(); i++)
            {
                // std::cout << "DEBUG: Looking at dirty electron " << i+1 << " of " << eDirtyCands.size() << std::endl ; 
                pat::Electron iE = eDirtyCands.at(i).first;
                // std::cout << "DEBUG: electron has pt " << hnu::getElectronEt(iE,false) << std::endl ; 
                if (hnu::addElectron(iE, hnuEvent, cuts.minimum_muon_jet_dR))
                {
                    // std::cout << "DEBUG: Setting dirtyPosition to " << i << std::endl ; 
                    dirtyPosition = i ;
                    // std::cout << "DEBUG: dirtyPosition is " << dirtyPosition << std::endl ; 
                    break ;
                }
            }
            for (unsigned int i = 0; i < eCands.size(); i++)
            {
                // std::cout << "DEBUG: Looking at clean electron " << i+1 << " of " << eCands.size() << std::endl ; 
                pat::Electron iE = eCands.at(i).first;
                // std::cout << "DEBUG: electron has pt " << hnu::getElectronEt(iE,false) << std::endl ; 
                if (hnu::addElectron(iE, hnuEvent, cuts.minimum_muon_jet_dR))
                {
                    // std::cout << "DEBUG: Setting cleanPosition to " << i << std::endl ; 
                    cleanPosition = i ;
                    // std::cout << "DEBUG: cleanPosition is " << cleanPosition << std::endl ; 
                    break ;
                }
            }
            if ( hnuEvent.nLeptons == 2 )
                hnuEvent.eventWgt *= (hnu::fakeProbability(eDirtyCands.at(dirtyPosition).first) / (1.0 - hnu::fakeProbability(eDirtyCands.at(dirtyPosition).first))) ;
        }
    }

    if(applyMuIDCorrections_ && hnuEvent.isMC) // Only care about applying weights if you have two candidates
    {
        if ( hnuEvent.nLeptons > 1 )
        { // Only care about applying weights if you have two candidates
            double e1wgt = muid_->weightElectronsForMC(hnu::getElectronSCEta(hnuEvent.e1), applyMuIDEffsign_) ;
            double e2wgt = muid_->weightElectronsForMC(hnu::getElectronSCEta(hnuEvent.e2), applyMuIDEffsign_) ;
            hnuEvent.eventWgt *= (e1wgt * e2wgt);
        }
    }
}

void HeavyNu::reselectJets(std::vector< std::pair<pat::Jet, float> >& jetCands, HeavyNuEvent& hnuEvent)
{

    double dRj1m1 = 999.0, dRj1m2 = 999.0,  dRj2m1 = 999.0, dRj2m2 = 999.0 ;
    if (hnuEvent.nLeptons > 1)
    {
        if (hnuEvent.nJets > 0)
        {
            dRj1m1 = deltaR(hnuEvent.mu1.eta(), hnuEvent.mu1.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
            dRj1m2 = deltaR(hnuEvent.mu2.eta(), hnuEvent.mu2.phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
        }
        if (hnuEvent.nJets > 1)
        {
            dRj2m1 = deltaR(hnuEvent.mu1.eta(), hnuEvent.mu1.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
            dRj2m2 = deltaR(hnuEvent.mu2.eta(), hnuEvent.mu2.phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
        }

        if ( dRj1m1 < cuts.minimum_muon_jet_dR || dRj1m2 < cuts.minimum_muon_jet_dR )
        { // Must replace 1st jet
            hnuEvent.nJets-- ;
            if ( dRj2m1 < cuts.minimum_muon_jet_dR || dRj2m2 < cuts.minimum_muon_jet_dR )
            { // Must also replace 2nd jet
                hnuEvent.nJets-- ;
                if ( jetCands.size() >= 3 )
                {
                    hnuEvent.j1 = jetCands.at(2).first ;
                    hnuEvent.j1scale = jetCands.at(2).second ;
                    hnuEvent.isBJet1 = hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal;
                    if ( jetCands.size() >= 4 )
                    {
                        hnuEvent.j2 = jetCands.at(3).first ;
                        hnuEvent.j2scale = jetCands.at(3).second ;
                        hnuEvent.isBJet2 = hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal;
                    }
                }
            }
            else
            { // Promote second jet
                hnuEvent.j1 = hnuEvent.j2 ;
                hnuEvent.j1scale = hnuEvent.j2scale ;
                hnuEvent.isBJet1 = hnuEvent.isBJet2 ;
                if ( jetCands.size() >= 3 )
                {
                    hnuEvent.j2 = jetCands.at(2).first ;
                    hnuEvent.j2scale = jetCands.at(2).second ;
                    hnuEvent.isBJet2 = hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal;
                }
            }
        }
        else if ( dRj2m1 < cuts.minimum_muon_jet_dR || dRj2m2 < cuts.minimum_muon_jet_dR )
        { // Must only replace 2nd jet
            hnuEvent.nJets-- ;
            if ( jetCands.size() >= 3 )
            {
                hnuEvent.j2 = jetCands.at(2).first ;
                hnuEvent.j2scale = jetCands.at(2).second ;
                hnuEvent.isBJet2 = hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal;
            }
        }
    }
}

void HeavyNu::selectTop(std::vector< std::pair<pat::Electron, float> >& eCands, std::vector<pat::Muon>& muCands, edm::Handle<reco::VertexCollection> pvHandle, HeavyNuEvent& hnuEvent)
{
    for(unsigned int i = 0; i < muCands.size(); i++)
    {
        if(hnuEvent.nMuons >= 1) break;
        pat::Muon iM = muCands.at(i);
        if(hnu::muIsolation(iM) < cuts.muon_trackiso_limit)
        {
            hnu::addMuon(iM, hnuEvent, cuts.minimum_muon_jet_dR, 0);
        }
    }

    if(applyMuIDCorrections_ && hnuEvent.isMC)
    {
        double mu1wgt = (hnuEvent.nMuons >= 1)? (muid_->weightForMCbyEta((hnuEvent.mu1.eta()), applyMuIDEffsign_)):1.0;

        hnuEvent.eventWgt *= mu1wgt;
    }

    for(unsigned int i = 0; i < eCands.size(); i++)
    {
        if(hnuEvent.nElectrons >= 1) break;
        pat::Electron iE = eCands.at(i).first;
        hnu::addElectron(iE, hnuEvent, cuts.minimum_muon_jet_dR);
    }

    if(applyMuIDCorrections_ && hnuEvent.isMC) // Only care about applying weights if you have two candidates
    {
        if ( hnuEvent.nElectrons >= 1 )
        { // Only care about applying weights if you have two candidates
            double e1wgt = muid_->weightElectronsForMC(hnu::getElectronSCEta(hnuEvent.e1), applyMuIDEffsign_) ;
            hnuEvent.eventWgt *= e1wgt;
        }
    }
}

void HeavyNu::alternativeSelection(HeavyNuEvent& hnuEvent)
{
    if ( studyAlternativeSelection_ && hnuEvent.nJets >= 2)
    {
        double l1pt = hnuEvent.l1pt ;
        double l2pt = hnuEvent.l2pt ;
        //double j1pt  = hnuEvent.j1scale * hnuEvent.j1.pt() ;
        double j2pt  = hnuEvent.j2.pt() ;
        if ( hnuEvent.mLL >= cuts.minimum_mumu_mass )
        { // Standard dimuon requirement
            if ( l1pt >= 40. ) hists.AlternativeMu1Pt40->fill(hnuEvent) ;
            if ( l1pt >=  cuts.minimum_mu1_pt )
            { // Standard mu1 pT requirement
                if ( l2pt >= 40. )
                {
                    hists.AlternativeMu2Pt40->fill(hnuEvent) ;
                    if ( l1pt >= 80. ) hists.AlternativeElecChanPt->fill(hnuEvent) ;
                }
                if ( l2pt >= 60. ) hists.AlternativeMu2Pt60->fill(hnuEvent) ;
                if ( j2pt >= 60. )  hists.AlternativeJetPt60->fill(hnuEvent) ;
                //--- Requirement for at least one muon to be barrel ---//
                bool atLeastOneBarrelMuonLoose = ( fabs(hnuEvent.l1->eta()) < 1.2 || fabs(hnuEvent.l2->eta()) < 1.2 ) ;
                bool atLeastOneBarrelMuonTight = ( fabs(hnuEvent.l1->eta()) < 0.8 || fabs(hnuEvent.l2->eta()) < 0.8 ) ;
                if ( atLeastOneBarrelMuonLoose )
                {
                    hists.AlternativeBarrelLoose->fill(hnuEvent) ;
                    if ( atLeastOneBarrelMuonTight ) hists.AlternativeBarrelTight->fill(hnuEvent) ;
                }
                //--- Checking for b-tagged jets using TCHE Loose ---//
                int nBtags = 0 ;
                if ( hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal ) nBtags++ ;
                if ( hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal ) nBtags++ ;
                if ( nBtags > 0 )  hists.AlternativeAtLeastOneBjet->fill(hnuEvent) ;
                if ( nBtags == 2 ) hists.AlternativeTwoBjets->fill(hnuEvent) ;
            }
        }
        if ( hnuEvent.mLL >= 120. )
        {
            if ( l1pt >=  cuts.minimum_mu1_pt )
            { // Standard mu1 pT requirement
                hists.AlternativeDimuonMass120->fill(hnuEvent) ;
            }
        }
    }

    if(studyScaleFactorEvolution_)
    {
        double l1pt = hnuEvent.l1pt;
        if(l1pt > 40.) hists.Mu1Pt40GeVCut->fill(hnuEvent);
        if(l1pt > 50.) hists.Mu1Pt50GeVCut->fill(hnuEvent);
        if(l1pt > 60.) hists.Mu1Pt60GeVCut->fill(hnuEvent);
        if(l1pt > 80.) hists.Mu1Pt80GeVCut->fill(hnuEvent);
        if(l1pt > 100.) hists.Mu1Pt100GeVCut->fill(hnuEvent);

        if(hnuEvent.l1pt < cuts.minimum_mu1_pt)
        {
            if(hnuEvent.n_primary_vertex == 1)
                hists.Mu1HighPtCutVtxEq1->fill(hnuEvent);
            else if(hnuEvent.n_primary_vertex <= 5)
                hists.Mu1HighPtCutVtx2to5->fill(hnuEvent);
            else if(hnuEvent.n_primary_vertex > 5)
                hists.Mu1HighPtCutVtxGt5->fill(hnuEvent);
        }
    }
}

// ------------ method called to for each event  ------------

bool HeavyNu::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    HeavyNuEvent hnuEvent(analysisMode_);
    if(addSlopeTree_) hnuTree_->clear();

    evtCounter++;
    // std::cout << "Event number: " << evtCounter << std::endl ; 

    hnuEvent.isMC = !iEvent.isRealData();
    hnuEvent.pfJets = isPFJets_;
    hnuEvent.scaleMuE(applyMESfactor_);

    if(iEvent.isRealData())
    {
        if((applyMESfactor_ != 1.0) || (applyJECUsign_ != 0.0))
            throw cms::Exception("Energy scale studies not allowed on data currently");

        uint32_t runn = iEvent.id().run();
        std::map<uint32_t, TH1 *>::const_iterator it = m_runHistos_.find(runn);
        TH1 *runh;
        if(it == m_runHistos_.end())
        {
            runh = bookRunHisto(runn);
            m_runHistos_[runn] = runh;
        }
        else
            runh = it->second;
        runh->Fill(1);
    }

    edm::Handle<reco::JPTJetCollection> jptJets;
    iEvent.getByLabel("JetPlusTrackZSPCorJetAntiKt5", jptJets); // Some day we should make this readable from the cfg file

    edm::Handle<pat::MuonCollection> pMuons;
    iEvent.getByLabel(muonTag_, pMuons);

    edm::Handle<pat::ElectronCollection> pElecs;
    iEvent.getByLabel(elecTag_, pElecs);

    edm::Handle<double> electronRhoHandle ;
    iEvent.getByLabel(elecRhoTag_, electronRhoHandle) ;
    elecRho_ = ((electronRhoHandle.isValid()) ? (*(electronRhoHandle.product())) : 0.) ;

    edm::Handle<pat::JetCollection> pJets;
    iEvent.getByLabel(jetTag_, pJets);

    edm::Handle<reco::TrackCollection> gTracks ;
    iEvent.getByLabel("generalTracks", gTracks) ;

    edm::Handle<pat::GenericParticleCollection> pTracks;
    iEvent.getByLabel(trackTag_, pTracks);

    edm::Handle<pat::METCollection> pMET;
    iEvent.getByLabel(metTag_, pMET);

    //Shirpa reweighting info
    edm::Handle<GenEventInfoProduct> geneventinfo;
    iEvent.getByLabel("generator", geneventinfo);

    //Gen jets for gen matching
    edm::Handle<std::vector<reco::GenJet> > genjets;

    if(hnuEvent.isMC)
    {
        edm::Handle<std::vector<PileupSummaryInfo> > pPU;
        iEvent.getByLabel("addPileupInfo", pPU);
        if(pileupEra_ < 20100) hnuEvent.eventWgt = 1.0;
        else
        {
            std::pair<float, double> pileup = hnu::pileupReweighting(pPU, MCweightByVertex_);
            hnuEvent.n_pue = pileup.first ; // Will only be used for studies, thus no syst. correction necessary
            hnuEvent.eventWgt *= pileup.second;
            //if ( fabs(puShift_) > 0.001 ) hnuEvent.eventWgt *= poissonNvtxShifter_.ShiftWeight( pileup.first ) ; //only for 2011A
        }
        //Shirpa reweighting
        hnuEvent.eventWgt *= geneventinfo->weight();
        // PDF reweighting
#ifdef DO_LHAPDF
        if (doPDFreweight_)
        {
            edm::Handle<GenEventInfoProduct> geip;
            iEvent.getByLabel("generator", geip);

            float Q = geip->pdf()->scalePDF;
            int id1 = geip->pdf()->id.first;
            int id2 = geip->pdf()->id.second;
            float x1 = geip->pdf()->x.first;
            float x2 = geip->pdf()->x.second;

            hnuEvent.eventWgt *= hnu::getPDFWeight(Q, id1, x1, id2, x2,
                                                   doPDFreweight_,
                                                   pdfReweightBaseId, pdfReweightTargetId);
        }
#endif

        hists.weights->Fill(hnuEvent.eventWgt);
        // generator information
        edm::Handle<reco::GenParticleCollection> genInfo;
        if(iEvent.getByLabel("genParticles", genInfo))
        {
            hnuEvent.decayID(*genInfo);
        }
        iEvent.getByLabel("ak5GenJets", genjets);
    }
    
    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
    hnuEvent.n_primary_vertex = hnu::numberOfPrimaryVertices(pvHandle);

    if(!pElecs.isValid() || !pMuons.isValid() || !pJets.isValid() || !(pMET.isValid() && (pMET->size() > 0)))
    {
        std::cout << "Exiting as valid PAT objects not found" << std::endl;
        std::cout << "Electrons: " << pElecs.isValid() << std::endl;
        std::cout << "Muons:     " << pMuons.isValid() << std::endl;
        std::cout << "Jets:      " << pJets.isValid() << std::endl;
        std::cout << "MET:       " << pMET.isValid() << std::endl;
        return false;
    }

    if(firstEvent_)
    {
        // handle the jet corrector parameters collection,
        // get the jet corrector parameters collection from the global tag
        //
        edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
        if ( isPFJets_ ) iSetup.get<JetCorrectionsRecord > ().get("AK5PFchs", JetCorParColl) ;
        else             iSetup.get<JetCorrectionsRecord > ().get("AK5Calo", JetCorParColl) ;

        // get the uncertainty parameters from the collection,
        // instantiate the jec uncertainty object
        if(applyJECUsign_)
        {
            JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
            jecuObj_ = new JetCorrectionUncertainty(JetCorPar);
        }

        // Some sanity checks inserted
        // If running on Monte Carlo, expect
        //   - pileup configuration to make sense
        //   - to apply trigger and ID corrections
        // If running on data
        //   - no corrections are applied
        // Of course, this changes completely when doing systematics checks
        if(studyMuonSelectionEff_) std::cout << "Histograms for studying muon reco/ID efficiency will be created" << std::endl;
        if(studyScaleFactorEvolution_) std::cout << "Histograms for Z scale factor cross checks will be created" << std::endl;
        if(applyJECUsign_)
        {
            std::cout << "Studies will be used to estimate JEC uncertainty" << std::endl;
            if(studyMuonSelectionEff_ || studyScaleFactorEvolution_)
            {
                // std::cout << "WARNING: You are performing studies with modified jets.  This is most likely wrong!" << std::endl;
                std::cout << "WARNING: Tag-and-probe and/or scale factor studies with modified jets.  Disabling!" << std::endl;
                studyMuonSelectionEff_     = false ;
                studyScaleFactorEvolution_ = false ;
            }
        }
        else std::cout << "Nominal Jet corrections applied" << std::endl;
        if(applyMESfactor_ != 1.0)
        {
            std::cout << "Studies will be used to estimate MES uncertainty: " << applyMESfactor_ << std::endl;
            if(studyMuonSelectionEff_ || studyScaleFactorEvolution_)
            {
                std::cout << "WARNING: Tag-and-probe and/or scale factor studies with modified muons.  Disabling!" << std::endl;
                studyMuonSelectionEff_     = false ;
                studyScaleFactorEvolution_ = false ;
            }
        }
        else std::cout << "No MES corrections applied" << std::endl;
        if(!disableTriggerCorrection_)
        {
            if(applyTrigEffsign_) std::cout << "Studies will be used to estimate trigger efficiency uncertainty" << std::endl;
            else std::cout << "Nominal trigger corrections applied" << std::endl;
        }

        if(hnuEvent.isMC)
        {
            int pileupYear = pileupEra_ ;
            int idYear = muid_->idEra();

            if ( abs(idYear) > 1 || (idYear < 0 && nDirtyCands_ == 0) )
                std::cout << "WARNING: Using bad(?) choice for ID: " << idYear << std::endl ;

            // bool allErasMatch = true;
            if(applyMuIDCorrections_)
            {
                if(applyMuIDEffsign_) std::cout << "Studies will be used to estimate Mu ID uncertainty: " << applyMuIDEffsign_
                        << std::endl;
                else std::cout << "Nominal Mu ID corrections applied" << std::endl;
            }
            else std::cout << "You have disabled Mu ID corrections.  Are you sure?" << std::endl;
            if(disableTriggerCorrection_)
            {
                std::cout << "You have disabled the trigger correction.  Are you sure?" << std::endl;
                // if(pileupYear != idYear) allErasMatch = false;
            }
                //             else
                //             {
                //                 allErasMatch = (pileupYear == idYear);
                //             }
                //             if(!allErasMatch)
                //             {
                //                 std::cout << "WARNING: You do not appear to have consistent corrections applied!" << std::endl;
                //                 std::cout << "         pileup era is " << pileupEra_ << ", year for mu ID is " << idYear
                //                         << std::endl;
                //             }
            else
            {
                std::cout << "Looking at corrections, I assume you are running with the " << pileupYear << " year settings" << std::endl;
            }
        }
        else
        {
            if(disableTriggerCorrection_)
                std::cout << "WARNING: You have disabled the trigger correction in data.  What are you doing?!?" << std::endl;
        }
        std::cout << "===============================================" << std::endl;

        firstEvent_ = false;
    }

    hists.mc_type->Fill(hnuEvent.mc_class, hnuEvent.eventWgt);
    hists.nelec->Fill(pElecs->size());
    hists.nmuAll->Fill(pMuons->size());
    hists.njet->Fill(pJets->size());
    //hists.nmet->Fill(pMET->size());

    studyJetVertex(pJets, jptJets, pMuons, hnuEvent.n_primary_vertex);

    //if(pMET->size())
    //    hists.met->Fill(pMET->at(0).pt());
    //else
    //    hists.met->Fill(0);

    for(size_t iJet = 0; iJet < pJets->size(); iJet++)
    {
        pat::JetRef iJ = pat::JetRef(pJets, iJet);
        fillBasicJetHistos(*iJ, iJet + 1);
    }

    int nloose = 0, ntight = 0;
    for(size_t iMuon = 0; iMuon < pMuons->size(); iMuon++)
    {
        pat::MuonRef iM = pat::MuonRef(pMuons, iMuon);
        if(!iM.isAvailable()) continue;
        if ( iM->globalTrack().isNonnull() )
        {
            reco::TrackRef cktTrack = (muon::tevOptimized(*iM, 200, 40., 17., 0.25)).first;
            fillBasicMuHistos(*iM, pvHandle);
            if(hnu::isLooseMuonNoPF(*iM))
            {
                nloose++;
                if ( !cktTrack.isNull() && hnu::isTightHighPtMuon(*iM, pvHandle) )
                    ntight++;
            }
        }
    }
    hists.nmuLoose->Fill(nloose);
    hists.nmuTight->Fill(ntight);

    hists.cutlevel->Fill(-1.0, hnuEvent.eventWgt);

    // Look for valid jets and put them in the event

    std::vector< std::pair<pat::Jet, float> > jetCands =
            hnu::getJetList(pJets, jecuObj_, cuts.minimum_jet_pt, cuts.maximum_jet_abseta, applyJECUsign_, jecVal_, hnuEvent.isMC, applyJERsign_);
    hnuEvent.nJets = jetCands.size();

    if(hnuEvent.nJets >= 2)
    {
        hists.cutlevel->Fill(0.0, hnuEvent.eventWgt);
    }
    fill(*pMuons, *pElecs, *pJets, hnuEvent, false, false, hists.noCuts);

    // Look for valid muons
    std::vector<pat::Muon> muCands =
            hnu::getMuonList(pMuons, pvHandle, muid_->idEra(),
                             cuts.minimum_mu2_pt, cuts.maximum_mu_abseta, applyMESfactor_, merUncertainty_, false);

    // Look for valid electrons
    std::vector< std::pair<pat::Electron, float> > eCands =
            hnu::getElectronList(pElecs, cuts.maximum_elec_abseta, cuts.minimum_mu2_pt, cuts.minimum_mu2_pt,
                                 heepVersion_, pvHandle, elecRho_);

    // std::cout << "DEBUG: Got basic objects" << std::endl ; 

    // In order to avoid bias, it is necessary to perform muon studies
    // immediately after first creating the muon list
    //bool debuggingEvents = false ;
    if (studyMuonSelectionEff_ && analysisMode_ == HeavyNuEvent::HNUMU)
    {
        hnu::studyMuonEff(iEvent, muCands, hnuEvent, jetCands, pTracks, gTracks, trig_, cuts, effinfo_, hists);
    }

    if (studyMuonSelectionEff_ && analysisMode_ == HeavyNuEvent::HNUE)
    {
        hnu::studyElectronEff(iEvent, eCands, pElecs, hnuEvent, jetCands, trig_, cuts, effinfo_, hists);
    }

    //if(hnuEvent.nJets < 2) return false;

    if(hnuEvent.nJets >= 1)hnuEvent.j1 = jetCands.at(0).first;
    if(hnuEvent.nJets >= 2)hnuEvent.j2 = jetCands.at(1).first;
    if(hnuEvent.nJets >= 1)hnuEvent.j1scale = jetCands.at(0).second;
    if(hnuEvent.nJets >= 2)hnuEvent.j2scale = jetCands.at(1).second;

    hnuEvent.btagName = btagName;
    if(hnuEvent.nJets >= 1) hnuEvent.isBJet1 = hnuEvent.j1.bDiscriminator(btagName) >= minBtagDiscVal;
    if(hnuEvent.nJets >= 2) hnuEvent.isBJet2 = hnuEvent.j2.bDiscriminator(btagName) >= minBtagDiscVal;

    //conduct gen study on jets
    if(hnuEvent.isMC)
    {
        studyJetMatching(hnuEvent, genjets);
    }

    if(hnuEvent.nJets >= 1) hnuEvent.tjV1 = hnu::caloJetVertex(hnuEvent.j1, *jptJets);
    if(hnuEvent.nJets >= 2) hnuEvent.tjV2 = hnu::caloJetVertex(hnuEvent.j2, *jptJets);

    bool l1trig = false;
    bool l2trig = false;


    if(analysisMode_ == HeavyNuEvent::HNUMU)
    {
        selectMuons(muCands, jetCands, pvHandle, hnuEvent);
        if (muonsInJets_) reselectJets(jetCands, hnuEvent) ;

        //--- Trigger Matching needed for efficiency studies ---//
        if(trig_->matchingEnabled() && iEvent.isRealData())
        {
            l1trig = (hnuEvent.nLeptons > 0) &&
                    trig_->isTriggerMatched(hnuEvent.mu1, iEvent, &(((HeavyNuMuHist*)hists.Mu1TrigMatchesInZwin)->trigHistos));
            l2trig = (hnuEvent.nLeptons > 1) &&
                    trig_->isTriggerMatched(hnuEvent.mu2, iEvent, &(((HeavyNuMuHist*)hists.Mu2TrigMatchesInZwin)->trigHistos));
        }
        else if(!iEvent.isRealData())
        {
            if(disableTriggerCorrection_)
            {
                l1trig = true;
                l2trig = true;
            }
            else
            {
                l1trig = (hnuEvent.nLeptons > 0) && trig_->simulateForMC(hnuEvent.mu1.pt(), hnuEvent.mu1.eta(), applyTrigEffsign_) && hnuEvent.mu1.pt() > 40;
                l2trig = (hnuEvent.nLeptons > 1) && trig_->simulateForMC(hnuEvent.mu2.pt(), hnuEvent.mu2.eta(), applyTrigEffsign_) && hnuEvent.mu2.pt() > 40;
            }
        }
    }// HeavyNuEvent::HNUMU
    else if(analysisMode_ == HeavyNuEvent::HNUE)
    {
        selectElectrons(eCands, pElecs, pvHandle, hnuEvent);

        // at this stage we need to assume we have >=2 electrons
        // use the first two in the list eCands to build a di-lepton mass
        double diEleMass(0.);
        double dummyEta(0.);
        if( hnuEvent.nLeptons >= 2 )
        {
            diEleMass = ( hnuEvent.e1.p4() +  hnuEvent.e2.p4() ).M() ;
        }
        bool l12trig = false ;
        if ( iEvent.isRealData() )
            l12trig = (hnuEvent.nLeptons > 1) && trig_->isTriggerMatched(hnuEvent.e1, hnuEvent.e2, iEvent) ;
        else
            l12trig = (hnuEvent.nLeptons > 1) && trig_->simulateForMCdiEleMass(diEleMass, dummyEta, applyTrigEffsign_) ;

        l1trig = l2trig = l12trig ;
    }
    else if(analysisMode_ == HeavyNuEvent::TOP)
    {
        selectTop(eCands, muCands, pvHandle, hnuEvent);

        //--- Trigger Matching needed for efficiency studies ---//
        if(trig_->matchingEnabled() && iEvent.isRealData())
        {
            l1trig = (hnuEvent.nMuons > 0) && trig_->isTriggerMatched(hnuEvent.mu1, iEvent, &(((HeavyNuMuHist*)hists.Mu1TrigMatchesInZwin)->trigHistos));
        }
        else if(!iEvent.isRealData())
        {
            if(disableTriggerCorrection_)
            {
                l1trig = true;
            }
            else
            {
                l1trig = (hnuEvent.nMuons > 0) && trig_->simulateForMC(hnuEvent.mu1.pt(), hnuEvent.mu1.eta(), applyTrigEffsign_) && hnuEvent.mu1.pt() > 40;
            }
        }

        l2trig = false; //the electron never triggers.
    }
    else if(analysisMode_ == HeavyNuEvent::TAUX && !iEvent.isRealData())
    {
        edm::Handle<reco::GenParticleCollection> genInfo;
        iEvent.getByLabel("genParticles", genInfo);

        reco::Candidate* ptau = findGenTau(genInfo, true);
        reco::Candidate* stau = findGenTau(genInfo, false);

        if(!ptau || !stau) return false;

        bool firstIsMuon = gRandom->Uniform(0, 1) > 0.5;
        hnuEvent.tl1 = hnu::decayTau(ptau,  firstIsMuon);
        hnuEvent.tl2 = hnu::decayTau(stau, !firstIsMuon);

        double dRj1 = 999.9, dRj2 = 999.9;
        if(hnuEvent.nJets > 0) dRj1 = deltaR(hnuEvent.tl1.Eta(), hnuEvent.tl1.Phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
        if(hnuEvent.nJets > 1) dRj2 = deltaR(hnuEvent.tl1.Eta(), hnuEvent.tl1.Phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
        if (dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR) hnuEvent.nLeptons++ ;
        dRj1 = 999.9, dRj2 = 999.9;
        if(hnuEvent.nJets > 0) dRj1 = deltaR(hnuEvent.tl2.Eta(), hnuEvent.tl2.Phi(), hnuEvent.j1.eta(), hnuEvent.j1.phi());
        if(hnuEvent.nJets > 1) dRj2 = deltaR(hnuEvent.tl2.Eta(), hnuEvent.tl2.Phi(), hnuEvent.j2.eta(), hnuEvent.j2.phi());
        if (dRj1 > cuts.minimum_muon_jet_dR && dRj2 > cuts.minimum_muon_jet_dR) hnuEvent.nLeptons++ ;

        if(hnuEvent.nLeptons < 2 || hnuEvent.tl1.Pt() < 40.0 || hnuEvent.tl2.Pt() < 40) return false;

        //--- Trigger Matching needed for efficiency studies ---//
        if(disableTriggerCorrection_)
        {
            l1trig = true;
        }
        else if(firstIsMuon)
        {
            l1trig = trig_->simulateForMC(hnuEvent.tl1.Pt(), hnuEvent.tl1.Eta(), 0) && hnuEvent.tl1.Pt() > 40;
        }
        else if(!firstIsMuon)
        {
            l1trig = trig_->simulateForMC(hnuEvent.tl2.Pt(), hnuEvent.tl2.Eta(), 0) && hnuEvent.tl2.Pt() > 40;
        }

        l2trig = false; //the electron never triggers.
    }

    // std::cout << "DEBUG: Checked Trigger Conditions" << std::endl ; 
    // std::cout << "DEBUG: Values --> " << hnuEvent.nJets << " " << jetCands.size() << ", " << hnuEvent.nLeptons << " " << muCands.size() << " " << eCands.size() << ", " ;
    // std::cout << "DEBUG: TrigValues --> " << l1trig << " " << l2trig << std::endl ; 

    hnuEvent.scaleMuE(applyMESfactor_);
    // std::cout << "DEBUG: Applied scaling" << std::endl ; 
    hnuEvent.calculate(correctEscale_); // calculate various details
    // std::cout << "DEBUG: Calculated...ready for basic plots" << std::endl ; 

    // JJ and LL cutlevel must be made here before definite number cuts are applied
    if(hnuEvent.nJets >= 2) fill(*pMuons, *pElecs, *pJets, hnuEvent, true, false, hists.JJptCuts);
    // std::cout << "DEBUG: Filled JJptCuts" << std::endl ; 
    if(hnuEvent.nLeptons >= 2) fill(*pMuons, *pElecs, *pJets, hnuEvent, false, true, hists.LLptCuts);
    // std::cout << "DEBUG: Filled LLptCuts" << std::endl ; 

    if (hnuEvent.nLeptons < 2) return false ;

    // QCD closure test: Exactly one jet
    if (hnuEvent.nJets == 1 && hnu::jetID(hnuEvent.j1) > 0)
    {
        if ( (l1trig && l2trig) &&
            (hnuEvent.l1pt > cuts.minimum_mu1_pt) )
        {
            // std::cout << "DEBUG: Filling OneJetTrigHighPt" << std::endl ; 
            hists.OneJetTrigHighPt->fill(hnuEvent) ;
            // std::cout << "DEBUG: Successful filling of OneJetTrigHighPt" << std::endl ; 

            if ( hnuEvent.mLL > cuts.minimum_mumu_mass )
                hists.OneJetdiLmassCut->fill(hnuEvent) ;
        }
    }

    if(hnuEvent.nJets < 2) return false;

    if((hnuEvent.nJets >= 1 && hnu::jetID(hnuEvent.j1) < 1) || (hnuEvent.nJets >= 2 && hnu::jetID(hnuEvent.j2) < 1)) return false;

    // std::cout << "DEBUG: Jet requirements" << std::endl ; 

    hists.cutlevel->Fill(1.0, hnuEvent.eventWgt); // Two highest pT muons that are isolated, separated from chosen jets
    hists.LLJJptCuts->fill(hnuEvent);
    if(pMET->size()) hnuEvent.met1 = pMET->at(0);

    // Fill slope fit tuple here
    if(addSlopeTree_)
    {
        hnuTree_->event_.mll = hnuEvent.mLL;
        hnuTree_->event_.mlljj = hnuEvent.mWR;
        hnuTree_->event_.l1pt = hnuEvent.l1pt;
        hnuTree_->event_.l1eta = hnuEvent.l1eta;
        hnuTree_->event_.l1phi = hnuEvent.l1phi;
        hnuTree_->event_.l2pt = hnuEvent.l2pt;
        hnuTree_->event_.l2eta = hnuEvent.l2eta;
        hnuTree_->event_.l2phi = hnuEvent.l2phi;
        hnuTree_->event_.l1jdR = 999.9;
        hnuTree_->event_.l2jdR = 999.9;
        for(size_t iJet = 0; iJet < pJets->size(); iJet++)
        {
            pat::JetRef iJ = pat::JetRef(pJets, iJet);
            if(iJ->pt() < 40) continue;
            hnuTree_->event_.l1jdR = std::min((double)hnuTree_->event_.l1jdR, deltaR(hnuEvent.l1eta, hnuEvent.l1phi, iJ->eta(), iJ->phi()));
            hnuTree_->event_.l2jdR = std::min((double)hnuTree_->event_.l2jdR, deltaR(hnuEvent.l2eta, hnuEvent.l2phi, iJ->eta(), iJ->phi()));
        }
        if ( analysisMode_ == HeavyNuEvent::HNUMU )
        {
            short id = 0;
            if(hnu::isLooseMuonNoPF(hnuEvent.mu1)) id |= 0x01;
            if(hnu::isLooseMuon(hnuEvent.mu1)) id |= 0x02;
            if(hnu::isTightMuonCore(hnuEvent.mu1)) id |= 0x04;
            if(hnu::isTightMuon(hnuEvent.mu1, pvHandle)) id |= 0x08;
            if(hnu::isTightHighPtMuon(hnuEvent.mu1, pvHandle)) id |= 0x10;
            hnuTree_->event_.l1id = id;
            id = 0;
            if(hnu::isLooseMuonNoPF(hnuEvent.mu2)) id |= 0x01;
            if(hnu::isLooseMuon(hnuEvent.mu2)) id |= 0x02;
            if(hnu::isTightMuonCore(hnuEvent.mu2)) id |= 0x04;
            if(hnu::isTightMuon(hnuEvent.mu2, pvHandle)) id |= 0x08;
            if(hnu::isTightHighPtMuon(hnuEvent.mu2, pvHandle)) id |= 0x10;
            hnuTree_->event_.l2id = id;
        }
        else if ( analysisMode_ == HeavyNuEvent::HNUE )
        {
            hnuTree_->event_.l1id = (hnu::passesHEEP(hnuEvent.e1,heepVersion_,elecRho_, pvHandle))?1:0;
            hnuTree_->event_.l2id = (hnu::passesHEEP(hnuEvent.e2,heepVersion_,elecRho_, pvHandle))?1:0;
        }
        else
        {
            hnuTree_->event_.l1id = 0;
            hnuTree_->event_.l2id = 0;
        }
        hnuTree_->event_.j1pt = hnuEvent.j1.pt();
        hnuTree_->event_.j1eta = hnuEvent.j1.eta();
        hnuTree_->event_.j1phi = hnuEvent.j1.phi();
        hnuTree_->event_.j2pt = hnuEvent.j2.pt();
        hnuTree_->event_.j2eta = hnuEvent.j2.eta();
        hnuTree_->event_.j2phi = hnuEvent.j2.phi();
        hnuTree_->event_.weight = hnuEvent.eventWgt;
        hnuTree_->event_.flavor = hnuEvent.mode;
        hnuTree_->event_.n_pileup = hnuEvent.n_pue;
        hnuTree_->event_.n_primaryVertex = hnuEvent.n_primary_vertex;
        hnuTree_->event_.cutlevel = 1;
    }

    //--- Trigger code needs to be updated...placeholder for now ---//  is this still a palce holder? -Joe 9/24/12
    if(!l1trig && !l2trig)
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false;
    }
    hists.cutlevel->Fill(2.0, hnuEvent.eventWgt); // Event meets trigger requirements
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 2;
    hists.TrigMatches->fill(hnuEvent);

    // hists.LLJJptCuts->fill(hnuEvent);

    //--- Impose vertex requirement here ---//
    float deltaVzJ1J2 = 0.0, deltaVzJ1M1 = 0.0, deltaVzJ2M2 = 0.0, deltaVzJ1M2 = 0.0, deltaVzJ2M1 = 0.0;
    if(hnuEvent.nJets >= 1 && analysisMode_ != HeavyNuEvent::TAUX)
    {
        deltaVzJ1M1 = fabs(hnuEvent.tjV1 - hnuEvent.l1->vertex().Z());
        deltaVzJ1M2 = fabs(hnuEvent.tjV1 - hnuEvent.l2->vertex().Z());
    }
    if(hnuEvent.nJets >= 2 && analysisMode_ != HeavyNuEvent::TAUX)
    {
        deltaVzJ2M2 = fabs(hnuEvent.tjV2 - hnuEvent.l2->vertex().Z());
        deltaVzJ2M1 = fabs(hnuEvent.tjV2 - hnuEvent.l1->vertex().Z());
        deltaVzJ1J2 = fabs(hnuEvent.tjV1 - hnuEvent.tjV2);
    }
    if( (cuts.maxJetVZsepCM > 0) &&
       ((deltaVzJ1J2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M1 >= cuts.maxJetVZsepCM) ||
       (deltaVzJ2M2 >= cuts.maxJetVZsepCM) || (deltaVzJ1M2 >= cuts.maxJetVZsepCM) ||
       (deltaVzJ2M1 >= cuts.maxJetVZsepCM)) )
        return false;
    float deltaVzM1M2 = (analysisMode_ != HeavyNuEvent::TAUX)?(fabs(hnuEvent.l1->vertex().Z() - hnuEvent.l2->vertex().Z())):(0.0);
    if(cuts.maxVertexZsep > 0 && deltaVzM1M2 >= cuts.maxVertexZsep)
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false;
    }

    hists.cutlevel->Fill(3.0, hnuEvent.eventWgt); // Event meets vertex requirements
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 3;
    //hists.VertexCuts->fill(hnuEvent);

    //--- The "basic" object, trigger, and (possibly) vertex requirements should be done ---//
    //--- Consider alternative selection requirements ---//
    alternativeSelection(hnuEvent);

    if(hnuEvent.l1pt < cuts.minimum_mu1_pt)
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false;
    }
    // std::cout << "DEBUG: Min pT requirement" << std::endl ; 

    hists.cutlevel->Fill(4.0, hnuEvent.eventWgt); // Event meets high muon pT requirements
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 4;
    hists.Mu1HighPtCut->fill(hnuEvent);
    //if(hnuEvent.isBJet1 || hnuEvent.isBJet2) hists.Mu1HighPtCut_1bjet->fill(hnuEvent);
    //if(hnuEvent.isBJet1 && hnuEvent.isBJet2) hists.Mu1HighPtCut_2bjet->fill(hnuEvent);

    if ( studyRatePerRun_ && inZmassWindow(hnuEvent.mLL) )
        hists.z2jetPerRun->Fill( iEvent.id().run() ) ;

    //Z-peak cuts
    if(hnuEvent.mLL > 71 && hnuEvent.mLL < 111)
    {
        hists.ZRegionCut->fill(hnuEvent);
    }

    if(hnuEvent.mLL < cuts.minimum_mumu_mass)
    {
        if(addSlopeTree_) hnuTree_->fill();
        return false; // dimuon mass cut
    }
    hists.cutlevel->Fill(5.0, hnuEvent.eventWgt); // Event meets dimuon mass requirements
    if(addSlopeTree_) hnuTree_->event_.cutlevel = 5;
    hists.diLmassCut->fill(hnuEvent);
    // std::cout << "DEBUG: Satisfied DiL mass requirement" << std::endl ; 

    // Change the final logic of the filter based on LQ meeting discussion:
    // Interest in seeing events that pass the dilepton mass requirement
    // if ( hnuEvent.mWR<cuts.minimum_mWR_mass ) return false;  // 4-object mass cut
    if(hnuEvent.mWR >= cuts.minimum_mWR_mass)
    {
        hists.cutlevel->Fill(6.0, hnuEvent.eventWgt); // Event meets W_R mass requirements
        if(addSlopeTree_) hnuTree_->event_.cutlevel = 6;
        hists.mWRmassCut->fill(hnuEvent);
    }
    if(addSlopeTree_) hnuTree_->fill();

    if(iEvent.isRealData())
    {
        bool l1posChg = (hnuEvent.l1->charge() > 0) ;
        bool l2posChg = (hnuEvent.l2->charge() > 0) ;

        std::string mode = "UNKNOWN" ;
        if ( analysisMode_ == HeavyNuEvent::HNUMU ) mode = "MUON" ;
        if ( analysisMode_ == HeavyNuEvent::HNUE )  mode = "ELECTRON" ;
        if ( analysisMode_ == HeavyNuEvent::TOP )   mode = "TOP" ;

        if ( hnuEvent.l1->charge() == 0 || hnuEvent.l2->charge() == 0 )
            std::cout << "WARNING: found lepton with zero charge"
                << hnuEvent.l1->charge() << " " << hnuEvent.l2->charge() << std::endl ;

        if ( hnuEvent.mWR > 2000.0 ) std::cout << " (HIGH MASS) " ;
        if ( hnuEvent.mWR > 800.0 && hnuEvent.mWR < 1400.0 ) std::cout << " (1 TeV BUMP) " ;
        if ( hnuEvent.mWR > 1600.0 && hnuEvent.mWR < 2200.0 && analysisMode_ == HeavyNuEvent::HNUE )
            std::cout << " (2 TeV SPIKE) " ;
        std::cout << "\t" << mode << " " << iEvent.id() << std::endl ;
        std::cout << "\tM(W_R)  = " << hnuEvent.mWR << " GeV";
        std::cout << ", M(NuR1) = " << hnuEvent.mNuR1 << " GeV";
        std::cout << ", M(NuR2) = " << hnuEvent.mNuR2 << " GeV" << std::endl;
        std::cout << "\tM(mumu) = " << hnuEvent.mLL << " GeV";
        std::cout << ", M(JJ) = " << hnuEvent.mJJ << " GeV" << std::endl;
        std::cout << "\tJets:   j1 ";
        outputCandidate(hnuEvent.j1);
        std::cout << ", j2 ";
        outputCandidate(hnuEvent.j2);
        std::cout << std::endl;
        std::cout << "\tLeptons: l1, l" << (l1posChg ? "+":"-");
        outputCandidate(*(hnuEvent.l1));
        std::cout << ", l2, l" << (l2posChg ? "+":"-");
        outputCandidate(*(hnuEvent.l2));
        std::cout << std::endl;
    }

    return true;
}

// #include "LHAPDF/LHAPDF.h"

// double HeavyNu::getPDFWeight(float Q, int id1, float x1, int id2, float x2) {


//   if (!doPDFreweight_) return 1.0;

//   LHAPDF::usePDFMember(1,pdfReweightBaseId);
//   double pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
//   double pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;

//   LHAPDF::usePDFMember(2,pdfReweightTargetId);
//   double newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
//   double newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;

//   double w=(newpdf1/pdf1*newpdf2/pdf2);

//   //  printf("My weight is %f\n",w);

//   return w;

// }


// ------------ method called once each job just before starting event loop  ------------

void HeavyNu::beginJob()
{
    //    nnif_->beginJob();
    firstEvent_ = true;
    evtCounter = 0;

    if (doPDFreweight_)
    {
        hnu::initPDFSet(1, pdfReweightBaseName);
        hnu::initPDFSet(2, pdfReweightTargetName);
    }

}

// ------------ method called once each job just after ending the event loop  ------------

void HeavyNu::endJob()
{
    //    nnif_->endJob();
    trig_->endJob();
    muid_->endJob();
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNu);



