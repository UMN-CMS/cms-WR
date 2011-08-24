#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "TVector3.h"

HeavyNuEvent::HeavyNuEvent(anal_type theMode) { 
  mode = theMode ; 
  eventWgt = 1.0 ; 
  cutlevel = -1 ; 
  nMuons = 0 ; 
  nJets = 0 ; 
  nElectrons = 0 ; 
}

// Fills the HNE object with relevant information
// void HeavyNuEvent::initialize(int mode, 
// 			      std::vector<pat::MuonRef> muons, 
// 			      std::vector< std::pair<pat::JetRef,float> > jets,
// 			      std::vector<pat::ElectronRef> electrons,
// 			      double muPtScale, double muIsoLimit, double muJdR) {

//   cutlevel = 0 ; // We do not pass any cuts to start

//   if ( mode == HNU ) { // Nominal running

//     if ( jets.size() < 2 ) return ; 
//     j1 = jets.at(0) ; j2 = jets.at(1) ; 
//     if ( (hnu::jetID(*j1) < 1) || (hnu::jetID(*j2) < 1) ) return ; 
//     cutlevel = 1 ; // Two highest pT jets in the event pass ID requirements

//     for (unsigned int i=0; i<muons.size(); i++) { 
//       if ( !mu2.isNull() ) break ; 
//       pat::MuonRef iM = muons.at(i) ; 
//       if ( hnu::muIsolation((*iM),muPtScale) < muIsoLimit ) {
// 	double dRj1 = deltaR(iM->eta(), iM->phi(), j1->eta(), j1->phi()) ; 
// 	double dRj2 = deltaR(iM->eta(), iM->phi(), j2->eta(), j2->phi()) ; 
// 	if (dRj1 > muJdR && dRj2 > muJdR) { 
// 	  if      ( mu1.isNull() ) mu1 = iM ; 
// 	  else if ( mu2.isNull() ) mu2 = iM ; 
// 	  else    std::cout << "WARNING: Expected empty muon position" << std::endl ; 
// 	}
//       }
//     }
//     if ( mu2.isNull() ) return ;
//     cutlevel = 2 ; // Two highest pT muons that are isolated, separated from chosen jets
    
//         if(trig_->matchingEnabled() &&
//             iEvent.isRealData())
//     {
//         if(inZmassWindow(hnuEvent.mMuMu))
//             hists.MuTightInZwin.fill(hnuEvent, v_null);

//         // require that one muon be BOTH tight and trigger-matched
//         //
//         mu1trig = mu1trig &&
//                 trig_->isTriggerMatched(hnuEvent.mu1, iEvent,
//                 &(hists.Mu1TrigMatchesInZwin.trigHistos));

//         mu2trig = mu2trig &&
//                 trig_->isTriggerMatched(hnuEvent.mu2, iEvent,
//                 &(hists.Mu2TrigMatchesInZwin.trigHistos));

//     }
//     else if(!iEvent.isRealData() && !disableTriggerCorrection_)
//     {
//         mu1trig = mu1trig && trig_->simulateForMC(applyMESfactor_ * hnuEvent.mu1->pt(), applyTrigEffsign_);
//         mu2trig = mu2trig && trig_->simulateForMC(applyMESfactor_ * hnuEvent.mu2->pt(), applyTrigEffsign_);
//     }
//   } else if ( mode == TOP ) { 
//   } else if ( mode == QCD ) { 
//   }
// }

void HeavyNuEvent::regularize() {
  mu[0]=mu1;
  if ( mode == HNU || mode == QCD ) mu[1] = mu2 ; 
  if ( mode == TOP ) e[1] = e1 ; 

  j[0]=j1;
  j[1]=j2;
  tjV[0]=tjV1;
  tjV[1]=tjV2;
}

void HeavyNuEvent::scaleMuE(double mufactor, double efactor) {
  MuScale   = mufactor ; 
  ElecScale = efactor ; 
}

static double planeCosAngle(const reco::Particle::Vector& plane1,
			    const reco::Particle::Vector& plane2,
			    const reco::Particle::Vector& vect) {
  reco::Particle::Vector cp=plane1.unit().Cross(plane2.unit()).unit();
  double rv=fabs(vect.unit().Dot(cp));
  //  std::cout << plane1 << " " << plane2 << " " << vect << " " << rv << std::endl;
  return rv;
}

static double triangleArea(const reco::Particle::Vector& v1,
			   const reco::Particle::Vector& v2,
			   const reco::Particle::Vector& v3) { 

  // All actions assume/enforce unit sphere
  double cos_theta12 = v1.unit().Dot(v2.unit()) ; 
  double cos_theta13 = v1.unit().Dot(v3.unit()) ; 
  double cos_theta23 = v2.unit().Dot(v3.unit()) ; 

  // No triangle: vectors are either on top of each other or back-to-back
  if ( fabs(cos_theta12) == 1 || fabs(cos_theta13) == 1 || fabs(cos_theta23) == 1 ) return -1. ; 

  // Now calculate angles on sphere surface using identities
  double sin_theta12 = sqrt( 1.0 - cos_theta12 * cos_theta12 ) ; 
  double sin_theta13 = sqrt( 1.0 - cos_theta13 * cos_theta13 ) ; 
  double sin_theta23 = sqrt( 1.0 - cos_theta23 * cos_theta23 ) ; 

  double angle12_sphere = (cos_theta12 - (cos_theta13*cos_theta23)) / (sin_theta13*sin_theta23) ; 
  double angle13_sphere = (cos_theta13 - (cos_theta12*cos_theta23)) / (sin_theta12*sin_theta23) ; 
  double angle23_sphere = (cos_theta23 - (cos_theta12*cos_theta13)) / (sin_theta12*sin_theta13) ;

  double pi = 3.14159265 ; 
  double area = ( angle12_sphere + angle13_sphere + angle23_sphere - pi ) ; 

  return area ; 
}

void HeavyNuEvent::calculateMuMu() {
  vMuMu    = MuScale*(mu1.p4()+mu2.p4());
  mMuMu    = vMuMu.M();
  czeta_mumu   =mu1.momentum().unit().Dot(mu2.momentum().unit());
  ctheta_mumu  =planeCosAngle(mu1.momentum(),mu2.momentum(),reco::Particle::Vector(0,0,1));
}

void HeavyNuEvent::calculateMuE() {
  vMuMu    = MuScale*mu1.p4() + ElecScale*e1.p4();
  mMuMu    = vMuMu.M();
}

void HeavyNuEvent::calculate() {

  reco::Particle::LorentzVector j1p4 = j1.p4();
  reco::Particle::LorentzVector j2p4 = j2.p4();

  reco::Particle::LorentzVector lep1p4 = mu1.p4();
  reco::Particle::LorentzVector lep2p4 = (mode == CLO) ? mu1.p4() : ((mode == TOP) ? e1.p4() : mu2.p4()) ;  

  // if doing JECU studies, apply scaling factor here
  //
  if( j1scale != 1.0 ) j1p4 *= j1scale;
  if( j2scale != 1.0 ) j2p4 *= j2scale;

  // std::cout << j1scale << "+++" << j2scale << std::endl ; 
 
  reco::Particle::Vector j1mom = j1p4.Vect();
  reco::Particle::Vector j2mom = j2p4.Vect();

  // if doing MES studies, apply scaling factor here
  //
  if ( mode == TOP ) {
    if ( MuScale != 1.0 )   lep1p4 *= MuScale;  
    if ( ElecScale != 1.0 ) lep2p4 *= ElecScale; 
  } else { 
    if ( MuScale != 1.0 ) { lep1p4 *= MuScale; lep2p4 *= MuScale; }
  }

  reco::Particle::Vector lep1mom = lep1p4.Vect();
  reco::Particle::Vector lep2mom = lep2p4.Vect();

  vMuMu  = lep1p4+lep2p4;
  vJJ    = j1p4+j2p4;
  lv_evt = vMuMu+vJJ;

  // Calculate opening angle of m1/m2 + jj on a sphere
  area_1jj = triangleArea(lep1mom,j1mom,j2mom) ; 
  area_2jj = triangleArea(lep2mom,j1mom,j2mom) ; 

  czeta_mumu   =fabs(lep1mom.unit().Dot(lep2mom.unit()));

  ctheta_mumu  =planeCosAngle(lep1mom,lep2mom,reco::Particle::Vector(0,0,1));
  ctheta_jj    =planeCosAngle(j1mom,j2mom,reco::Particle::Vector(0,0,1));
  ctheta_mu1_jj=planeCosAngle(j1mom,j2mom,lep1mom);
  ctheta_mu2_jj=planeCosAngle(j1mom,j2mom,lep2mom);

  // LorentzVector of just the Z deboost.
  reco::Particle::LorentzVector deboostz(0,0,-lv_evt.pz(),lv_evt.pz());
  
  reco::Particle::LorentzVector lep1z=lep1p4+deboostz;
  reco::Particle::LorentzVector lep2z=lep2p4+deboostz;
  reco::Particle::LorentzVector j1z=j1p4+deboostz;
  reco::Particle::LorentzVector j2z=j2p4+deboostz;
  
  cthetaz_mumu   = planeCosAngle(lep1z.Vect(),lep2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_jj     = planeCosAngle(j1z.Vect(),j2z.Vect(),reco::Particle::Vector(0,0,1));
  cthetaz_mu1_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),lep1z.Vect());
  cthetaz_mu2_jj = planeCosAngle(j1z.Vect(),j2z.Vect(),lep2z.Vect());

  float dRlep1jet1 = deltaR( mu1.eta(), mu1.phi(), j1.eta(), j1.phi() ) ; 
  float dRlep1jet2 = deltaR( mu1.eta(), mu1.phi(), j2.eta(), j2.phi() ) ; 
  float dRlep2jet1 = (( mode == TOP ) ? 
		      (deltaR( e1.eta(), e1.phi(), j1.eta(), j1.phi() )) :  
		      (( mode == CLO ) ? 
		       dRlep1jet1 : 
		       (deltaR( mu2.eta(), mu2.phi(), j1.eta(), j1.phi()))) ) ;  
  float dRlep2jet2 = (( mode == TOP ) ? 
		      (deltaR( e1.eta(), e1.phi(), j2.eta(), j2.phi() )) : 
		      (( mode == CLO ) ? 
		       dRlep1jet2 : 
		       (deltaR( mu2.eta(), mu2.phi(), j2.eta(), j2.phi()))) ) ; 
  
  // find the closest jets
  dRminMu1jet = std::min( dRlep1jet1,dRlep1jet2 );
  dRminMu2jet = std::min( dRlep2jet1,dRlep2jet2 );

  // what are the muon transverse momenta relative to the closest jets?
  reco::Particle::Vector jmom4lep1 = (dRminMu1jet == dRlep1jet1) ? j1mom : j2mom;
  reco::Particle::Vector jmom4lep2 = (dRminMu2jet == dRlep2jet1) ? j1mom : j2mom;

  TVector3 lep1vec( lep1mom.X(), lep1mom.Y(), lep1mom.Z() );
  TVector3 lep2vec( lep2mom.X(), lep2mom.Y(), lep2mom.Z() );

  TVector3 jt4lep1vec( jmom4lep1.X(), jmom4lep1.Y(), jmom4lep1.Z() );
  TVector3 jt4lep2vec( jmom4lep2.X(), jmom4lep2.Y(), jmom4lep2.Z() );

  ptrelMu1 = lep1vec.Perp( jt4lep1vec );
  ptrelMu2 = lep2vec.Perp( jt4lep2vec );

  // Composite objects
  mJJ   = vJJ.M();
  mMuMu = vMuMu.M();

  mWR   = lv_evt.M();

  mNuR1 = (vJJ + lep1p4).M();
  mNuR2 = (vJJ + lep2p4).M();
}
/*
void HeavyNuEvent::decayID(const HepMC::GenEvent& genE) {

  HepMC::GenEvent::vertex_const_iterator vtex;
  HepMC::GenVertex::particles_out_const_iterator Pout;

  mc_class=0;

  for (vtex=genE.vertices_begin();vtex!=genE.vertices_end();vtex++){
    if(((*vtex)->particles_in_size())==1){ // We want a decay, not collision
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==9900024){ // Is a WR
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 1;
            break;
	    
          case 13:	// mu
            mc_class = 2;
            break;
	    
          case 15:  // tau
            mc_class = 3;
            break;
	    
          default:	// else
	    break;
	  }
	  if (mc_class!=0) break;
	
	}
	break; //end of id loop
      }
      if(abs((*((*vtex)->particles_in_const_begin()))->pdg_id())==23){ // Is a Z
	for(Pout=(*vtex)->particles_out_const_begin(); Pout!=(*vtex)->particles_out_const_end(); Pout++){
	  int pdf_id = abs((*Pout)->pdg_id());
	  switch ( pdf_id ) {
          case 11:	// e
            mc_class = 11;
            break;
	    
          case 13:	// mu
            mc_class = 12;
            break;
	    
          case 15:  // tau
            mc_class = 13;
            break;
	    
          default:	// else
            break;
	  }
	  if (mc_class!=0) break;	
	  
	}
	break; //end of id loop
      }
    }

  }
}
*/
void HeavyNuEvent::decayID(const reco::GenParticleCollection& gpp) {
  mc_class=0;

  reco::GenParticleCollection::const_iterator i;
  for (i=gpp.begin(); i!=gpp.end(); i++) {
    if (abs(i->pdgId())==9900024 && i->numberOfDaughters()>=2) {
      int apid=abs(i->daughter(0)->pdgId());
      if (apid==11) mc_class=1;
      if (apid==13) mc_class=2;
      if (apid==15) mc_class=3;
      break;
    }
    if (abs(i->pdgId())==23 && i->numberOfDaughters()>=2) {
      int apid=abs(i->daughter(0)->pdgId());
      if (apid==11) mc_class=11;
      if (apid==13) mc_class=12;
      if (apid==15) mc_class=13;
      break;
      /*
      std::cout << i->pdgId() << " " << i->numberOfDaughters() << "  ";
      for (unsigned int j=0; j<i->numberOfDaughters(); j++)
	std::cout  << i->daughter(j)->pdgId() << " " ;
      std::cout << std::endl;
      */
    }
  }
}
