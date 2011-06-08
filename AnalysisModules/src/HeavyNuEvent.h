/** -- C++ -- **/
#ifndef HeavyNuEvent_h_included
#define HeavyNuEvent_h_included

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

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


/** The purpose of this class is contain the key items for
    a HeavyNuEvent and provide a simple way to pass this information
    between sections of the code.
*/
typedef std::pair<int,int> hNuMassHypothesis;

class HeavyNuEvent {
public:
  HeavyNuEvent() { eventWgt = 1.0 ; }
  void regularize();
  void calculate(int nMu);
  void calculate() { return calculate(2) ; }
  void calculateMuMu(double muptfactor);
  void calculateMuE(double muptfactor,double elefactor);
  void decayID(const HepMC::GenEvent& genE);

  bool isMC;
  // mc_class=0 (something else), 1=ee, 2=mm, 3=tau tau
  int mc_class;

  pat::MuonRef mu1, mu2, mu[2];
  pat::JetRef  j1,  j2,  j[2];
  pat::ElectronRef e1, e2, e[2]; 

  double tjV1, tjV2, tjV[2];
  int n_primary_vertex, n_pue;

  pat::METRef met1;

  // separately stored for JEC Uncertainty studies
  // (saves space and time not copying whole jet objects,
  //  particularly during the jet selection)
  //
  float j1scale, j2scale;

  float MESscale;
  float EEScale;

  double eventWgt ; 

  reco::Particle::LorentzVector vMuMu;
  reco::Particle::LorentzVector vJJ;
  reco::Particle::LorentzVector lv_evt;
  reco::Particle::LorentzVector WR;

  double ctheta_mu1_jj, cthetaz_mu1_jj;
  double ctheta_mu2_jj, cthetaz_mu2_jj;

  double czeta_mumu; // cosine of 3D angle between the muon 3-mom vectors

  double ctheta_mumu, cthetaz_mumu;
  double ctheta_jj,   cthetaz_jj;

  double dRminMu1jet, dRminMu2jet;
  double ptrelMu1,    ptrelMu2;

  double mWR, mJJ, mMuMu, mNuR1, mNuR2;

  std::vector<float> nnoutputs;
};

#endif
