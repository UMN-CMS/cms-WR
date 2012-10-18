/** -- C++ -- **/
#ifndef HeavyNuEvent_h_included
#define HeavyNuEvent_h_included

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TH1F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

/** The purpose of this class is contain the key items for
    a HeavyNuEvent and provide a simple way to pass this information
    between sections of the code.
 */
typedef std::pair<int, int> hNuMassHypothesis;

class HeavyNuEvent
{
public:

    enum anal_type
    {
        HNUMU, HNUE, TOP, QCD, CLO
    };
    anal_type mode;

    HeavyNuEvent(anal_type theMode = HNUMU);
    HeavyNuEvent(const HeavyNuEvent& hnue);

    // void initialize(int mode);

    void scaleMuE(double mufactor = 1.0, double efactor = 1.0);
    void calculateLL(bool correctEscale = false);
    void calculate(bool correctEscale = false);
    void decayID(const reco::GenParticleCollection& gpc);

    // These generic particles are used to pass features common to all leptons.
    // Only options common to all leptons should included here.
    pat::GenericParticle *l1, *l2;

    double l1pt, l2pt, l1eta, l2eta, l1phi, l2phi;

    // This function sets the generic leptons
    void regularize();

    bool isMC, pfJets;
    int numNuLJetsMatched;
    // mc_class=0 (something else), 1=ee, 2=mm, 3=tau tau
    int mc_class;

    pat::Muon mu1, mu2;
    pat::Jet j1, j2;
    pat::Electron e1, e2;

    reco::GenParticle gm1, gm2;
    reco::GenJet gj1, gj2;

    bool isBJet1, isBJet2;

    int nLeptons, nJets, nMuons, nElectrons, numBJets;

    double tjV1, tjV2;
    int n_primary_vertex, n_pue;

    pat::MET met1;

    // separately stored for JEC Uncertainty studies
    // (saves space and time not copying whole jet objects,
    //  particularly during the jet selection)
    //
    float j1scale, j2scale;

    float MuScale, ElecScale;

    int cutlevel;
    double eventWgt;

    reco::Particle::LorentzVector vLL;
    reco::Particle::LorentzVector vJJ;
    reco::Particle::LorentzVector lv_evt;
    reco::Particle::LorentzVector WR;

    double dRminL1jet, dRminL2jet;
    double ptrelL1, ptrelL2;

    double mWR, mJJ, mLL, mNuR1, mNuR2;

    std::string btagName;
};

#endif
