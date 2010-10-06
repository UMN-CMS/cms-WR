/** -- C++ -- **/
#ifndef HeavyNuEvent_h_included
#define HeavyNuEvent_h_included

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

/** The purpose of this class is contain the key items for
    a HeavyNuEvent and provide a simple way to pass this information
    between sections of the code.
*/
class HeavyNuEvent {
public:
  HeavyNuEvent() { reset(); }
  void reset();
  void regularize();
  void calculate();
  void calculateMuMu();

  bool isMC;

  const pat::Muon *mu1, *mu2, *mu[2];
  const pat::Jet  *j1, *j2, *j[2];

  reco::Particle::LorentzVector vMuMu;
  reco::Particle::LorentzVector vJJ;
  reco::Particle::LorentzVector lv_evt;
  reco::Particle::LorentzVector WR;

  double ctheta_mu1_jj, cthetaz_mu1_jj;
  double ctheta_mu2_jj, cthetaz_mu2_jj;

  double ctheta_mumu, cthetaz_mumu;
  double ctheta_jj,   cthetaz_jj;

  double dRminMu1jet, dRminMu2jet;
  double ptrelMu1, ptrelMu2;

  double mWR, mJJ, mMuMu, mNuR1, mNuR2;

};

#endif
