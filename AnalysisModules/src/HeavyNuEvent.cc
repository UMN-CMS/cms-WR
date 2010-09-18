#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

void HeavyNuEvent::reset() {
  mu1=0;
  mu2=0;
  mu[0]=0;
  mu[1]=0;
  j1=0;
  j2=0;
  j[0]=0;
  j[1]=0;
}


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


void HeavyNuEvent::calculate() {

  mumu=mu1->p4()+mu2->p4();
  jj=j1->p4()+j2->p4();
  lv_evt=mumu+jj;

  
  ctheta_mumu=planeCosAngle(mu1->momentum(),mu2->momentum(),reco::Particle::Vector(0,0,1));
  ctheta_jj=planeCosAngle(j1->momentum(),j2->momentum(),reco::Particle::Vector(0,0,1));
  ctheta_mu1_jj=planeCosAngle(j1->momentum(),j2->momentum(),mu1->momentum());
  ctheta_mu2_jj=planeCosAngle(j1->momentum(),j2->momentum(),mu2->momentum());

  // LorentzVector of just the Z deboost.
  reco::Particle::LorentzVector deboostz(0,0,-lv_evt.pz(),lv_evt.pz());
  
  reco::Particle::LorentzVector mu1z=mu1->p4()+deboostz;
  reco::Particle::LorentzVector mu2z=mu2->p4()+deboostz;
  reco::Particle::LorentzVector j1z=j1->p4()+deboostz;
  reco::Particle::LorentzVector j2z=j2->p4()+deboostz;
  
  cthetaz_mumu=planeCosAngle(mu1z.Vect(),mu2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_jj=planeCosAngle(j1z.Vect(),j2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_mu1_jj=planeCosAngle(j1z.Vect(),j2z.Vect(),mu1z.Vect());
  cthetaz_mu2_jj=planeCosAngle(j1z.Vect(),j2z.Vect(),mu2z.Vect());



}
