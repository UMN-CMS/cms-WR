#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "TVector3.h"

void HeavyNuEvent::regularize() {
  mu[0]=mu1;
  mu[1]=mu2;
  j[0]=j1;
  j[1]=j2;
}

static double planeCosAngle(const reco::Particle::Vector& plane1,
			    const reco::Particle::Vector& plane2,
			    const reco::Particle::Vector& vect) {
  reco::Particle::Vector cp=plane1.unit().Cross(plane2.unit()).unit();
  double rv=fabs(vect.unit().Dot(cp));
  //  std::cout << plane1 << " " << plane2 << " " << vect << " " << rv << std::endl;
  return rv;
}

void HeavyNuEvent::calculateMuMu() {
  vMuMu  = mu1->p4()+mu2->p4();
  czeta_mumu   =mu1->momentum().unit().Dot(mu2->momentum().unit());
  ctheta_mumu  =planeCosAngle(mu1->momentum(),mu2->momentum(),reco::Particle::Vector(0,0,1));
  mMuMu = vMuMu.M();
}

void HeavyNuEvent::calculate() {

  vMuMu  = mu1->p4()+mu2->p4();
  vJJ    = j1->p4()+j2->p4();
  lv_evt = vMuMu+vJJ;

  czeta_mumu   =fabs(mu1->momentum().unit().Dot(mu2->momentum().unit()));

  ctheta_mumu  =planeCosAngle(mu1->momentum(),mu2->momentum(),reco::Particle::Vector(0,0,1));
  ctheta_jj    =planeCosAngle(j1->momentum(),j2->momentum(),reco::Particle::Vector(0,0,1));
  ctheta_mu1_jj=planeCosAngle(j1->momentum(),j2->momentum(),mu1->momentum());
  ctheta_mu2_jj=planeCosAngle(j1->momentum(),j2->momentum(),mu2->momentum());

  // LorentzVector of just the Z deboost.
  reco::Particle::LorentzVector deboostz(0,0,-lv_evt.pz(),lv_evt.pz());
  
  reco::Particle::LorentzVector mu1z=mu1->p4()+deboostz;
  reco::Particle::LorentzVector mu2z=mu2->p4()+deboostz;
  reco::Particle::LorentzVector j1z=j1->p4()+deboostz;
  reco::Particle::LorentzVector j2z=j2->p4()+deboostz;
  
  cthetaz_mumu   = planeCosAngle(mu1z.Vect(),mu2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_jj     = planeCosAngle(j1z.Vect(),j2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_mu1_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),mu1z.Vect());
  cthetaz_mu2_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),mu2z.Vect());

  float dRmu1jet1 = deltaR(mu1->eta(), mu1->phi(), j1->eta(), j1->phi()) ; 
  float dRmu1jet2 = deltaR(mu1->eta(), mu1->phi(), j2->eta(), j2->phi()) ; 
  float dRmu2jet1 = deltaR(mu2->eta(), mu2->phi(), j1->eta(), j1->phi()) ; 
  float dRmu2jet2 = deltaR(mu2->eta(), mu2->phi(), j2->eta(), j2->phi()) ; 

  // find the closest jets
  dRminMu1jet = std::min(dRmu1jet1,dRmu1jet2);
  dRminMu2jet = std::min(dRmu2jet1,dRmu2jet2);

  // what are the muon transverse momenta relative to the closest jets?
  pat::JetRef  j4mu1, j4mu2;
  j4mu1 = (dRminMu1jet == dRmu1jet1) ? j1 : j2;
  j4mu2 = (dRminMu2jet == dRmu2jet1) ? j1 : j2;

  TVector3 mu1vec(mu1->momentum().X(), mu1->momentum().Y(), mu1->momentum().Z());
  TVector3 mu2vec(mu2->momentum().X(), mu2->momentum().Y(), mu2->momentum().Z());

  TVector3 jt1vec(j4mu1->p4().Vect().X(), j4mu1->p4().Vect().Y(), j4mu1->p4().Vect().Z() );
  TVector3 jt2vec(j4mu2->p4().Vect().X(), j4mu2->p4().Vect().Y(), j4mu2->p4().Vect().Z() );

  ptrelMu1 = mu1vec.Perp(jt1vec);
  ptrelMu2 = mu2vec.Perp(jt2vec);

  // Composite objects
  mJJ   = vJJ.M();
  mMuMu = vMuMu.M();

  mWR   = lv_evt.M();

  mNuR1 = (vJJ + mu1->p4()).M();
  mNuR2 = (vJJ + mu2->p4()).M();
}
