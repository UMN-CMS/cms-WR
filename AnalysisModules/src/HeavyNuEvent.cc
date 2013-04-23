#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"
#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "TVector3.h"

HeavyNuEvent::HeavyNuEvent(anal_type theMode)
{
	//anal_type theMode = HNUMU;
    mode = theMode ;
    eventWgt = 1.0 ;
    cutlevel = -1 ;
    nLeptons = 0 ;
    nJets = 0 ;
    nElectrons = 0 ;
    nMuons = 0;
    numNuLJetsMatched = 0;
    isBJet1 = isBJet2 = false;

    // Protection: set all scale factors to 1 by default
    j1scale = 1.0 ;
    j2scale = 1.0 ;
    MuScale = 1.0 ;
    ElecScale = 1.0 ;

    l1 = l2 = NULL;
}

HeavyNuEvent::HeavyNuEvent(const HeavyNuEvent& hne)
{
    mode = hne.mode;
    eventWgt = hne.eventWgt;
    cutlevel = hne.cutlevel;
    nLeptons = hne.nLeptons;
    nJets = hne.nJets;
    nMuons = hne.nMuons;
    nElectrons = hne.nElectrons;
    numNuLJetsMatched = hne.numNuLJetsMatched;
    isBJet1 = hne.isBJet1;
    isBJet2 = hne.isBJet2;
    numBJets = hne.numBJets;

    j1scale = hne.j1scale;
    j2scale = hne.j2scale;
    MuScale = hne.MuScale;
    ElecScale = hne.ElecScale;
    
    mu1 = hne.mu1;
    mu2 = hne.mu2;
    j1 = hne.j1;
    j2 = hne.j2;
    e1 = hne.e1;
    e2 = hne.e2;

    gm1 = hne.gm1;
    gm2 = hne.gm2;
    gj1 = hne.gj1;
    gj2 = hne.gj2;

    tjV1 = hne.tjV1;
    tjV2 = hne.tjV2;
    
    n_primary_vertex = hne.n_primary_vertex;
    n_pue = hne.n_pue;

    met1 = hne.met1;
}

void HeavyNuEvent::scaleMuE(double mufactor, double efactor)
{
    MuScale   = mufactor ;
    ElecScale = efactor ;
}

void HeavyNuEvent::calculateLL(bool correctEscale)
{
	std::cout << "calculateLL is depricated, stop using it, NOW!!!!!!!!!!!!!!!!!!" << std::endl;
    //double scale = 1.0 ; 
    //switch(mode)
    //{
    //    case HNUE:
    //if ( isMC && correctEscale ) 
    //  scale = hnu::getElectronEscale( hnu::getElectronSCEta(e1),hnu::getElectronSCEta(e2) ) ;  
    //vLL = scale * ( e1.p4() + e2.p4() ) ;
    //break;
    //    case TOP:
    //if ( isMC && correctEscale ) 
    //  scale = hnu::getElectronEscale( hnu::getElectronSCEta(e1) ) ;  
    //        vLL = (scale * e1.p4()) + mu1.p4();
    //        break;
    //    case HNUMU:
    //        vLL = mu1.p4() + mu2.p4();
    //        break;
    //    case QCD:
    //    case CLO:
    //        break;
    //}
    //mLL = vLL.M();
}

void HeavyNuEvent::calculate(bool correctEscale)
{
    regularize();

    reco::Particle::LorentzVector j1p4 = j1.p4();
    reco::Particle::LorentzVector j2p4 = j2.p4();

    reco::Particle::LorentzVector lep1p4, lep2p4;

    switch(mode)
    {
        case HNUE:
        case TAUX:
            lep1p4 = reco::Particle::PolarLorentzVector(l1pt, l1eta, l1phi, 0);
            lep2p4 = reco::Particle::PolarLorentzVector(l2pt, l2eta, l2phi, 0);
            break;
        case TOP:
            lep1p4 = reco::Particle::PolarLorentzVector(l1pt, l1eta, l1phi, (nLeptons >= 2 && mu1.pt() > hnu::getElectronEt(e1, false))?0.105:0.0);
            lep2p4 = reco::Particle::PolarLorentzVector(l2pt, l2eta, l2phi, (nLeptons >= 2 && mu1.pt() > hnu::getElectronEt(e1, false))?0.0:0.105);
            break;
        case HNUMU:
            lep1p4 = mu1.p4();
            lep2p4 = mu2.p4();
            break;
        case QCD:
            lep1p4 = mu1.p4();
            lep2p4 = mu2.p4();
            break;
        case CLO:
            break;
    }

    // if doing JECU studies, apply scaling factor here
    //
    if( j1scale != 1.0 ) j1p4 *= j1scale;
    if( j2scale != 1.0 ) j2p4 *= j2scale;

    // std::cout << j1scale << "+++" << j2scale << std::endl ;

    reco::Particle::Vector j1mom = j1p4.Vect();
    reco::Particle::Vector j2mom = j2p4.Vect();

    // if doing MES studies, apply scaling factor here
    //
    //if ( mode == TOP ) {
    //  if ( MuScale != 1.0 )   lep1p4 *= MuScale;
    //  if ( ElecScale != 1.0 ) lep2p4 *= ElecScale;
    //}
    //else {
    //  if ( MuScale != 1.0 ) { lep1p4 *= MuScale; lep2p4 *= MuScale; }
    //}

    reco::Particle::Vector lep1mom = lep1p4.Vect();
    reco::Particle::Vector lep2mom = lep2p4.Vect();

    double scale = 1.0 ; 
    switch(mode)
    {
        case HNUE:
	    if(isMC && correctEscale) scale = hnu::getElectronEscale(l1eta, l2eta); 
	    vLL = scale * ( lep1p4 + lep2p4 ) ;
	    break;
        case TOP:
	    if(isMC && correctEscale) scale = hnu::getElectronEscale(l2eta); 
            vLL = lep1p4 + (scale * lep2p4) ; 
            break;
        default:
            vLL = lep1p4 + lep2p4 ;
            break;
    }

    // vLL  = lep1p4 + lep2p4;
    vJJ    = j1p4 + j2p4;
    lv_evt = vLL + vJJ;

    // LorentzVector of just the Z deboost.
    //reco::Particle::LorentzVector deboostz(0, 0, -lv_evt.pz(), lv_evt.pz());

    //reco::Particle::LorentzVector lep1z=lep1p4+deboostz;
    //reco::Particle::LorentzVector lep2z=lep2p4+deboostz;
    //reco::Particle::LorentzVector j1z=j1p4+deboostz;
    //reco::Particle::LorentzVector j2z=j2p4+deboostz;

    float dRlep1jet1 = deltaR( lep1p4.eta(), lep1p4.phi(), j1.eta(), j1.phi() ) ;
    float dRlep1jet2 = deltaR( lep1p4.eta(), lep1p4.phi(), j2.eta(), j2.phi() ) ;
    float dRlep2jet1 = deltaR( lep2p4.eta(), lep2p4.phi(), j1.eta(), j1.phi() ) ;
    float dRlep2jet2 = deltaR( lep2p4.eta(), lep2p4.phi(), j2.eta(), j2.phi() ) ;

    // find the closest jets
    dRminL1jet = std::min( dRlep1jet1, dRlep1jet2 );
    dRminL2jet = std::min( dRlep2jet1, dRlep2jet2 );

    // what are the muon transverse momenta relative to the closest jets?
    reco::Particle::Vector jmom4lep1 = (dRminL1jet == dRlep1jet1) ? j1mom : j2mom;
    reco::Particle::Vector jmom4lep2 = (dRminL2jet == dRlep2jet1) ? j1mom : j2mom;

    TVector3 lep1vec( lep1mom.X(), lep1mom.Y(), lep1mom.Z() );
    TVector3 lep2vec( lep2mom.X(), lep2mom.Y(), lep2mom.Z() );

    TVector3 jt4lep1vec( jmom4lep1.X(), jmom4lep1.Y(), jmom4lep1.Z() );
    TVector3 jt4lep2vec( jmom4lep2.X(), jmom4lep2.Y(), jmom4lep2.Z() );

    ptrelL1 = lep1vec.Perp( jt4lep1vec );
    ptrelL2 = lep2vec.Perp( jt4lep2vec );

    // Composite objects
    mJJ = vJJ.M();
    mLL = vLL.M();
    
    WR = vJJ + vLL;

    mWR   = lv_evt.M();

    mNuR1 = (vJJ + lep1p4).M();
    mNuR2 = (vJJ + lep2p4).M();
}

void HeavyNuEvent::decayID(const reco::GenParticleCollection& gpp)
{
    mc_class = 0;

    reco::GenParticleCollection::const_iterator i;
    for (i = gpp.begin(); i != gpp.end(); i++)
    {
        if (abs(i->pdgId()) == 9900024 && i->numberOfDaughters() >= 2)
        {
            int apid = abs(i->daughter(0)->pdgId());
            if (apid == 11) mc_class = 1;
            if (apid == 13) mc_class = 2;
            if (apid == 15) mc_class = 3;
            break;
        }
        if (abs(i->pdgId()) == 23 && i->numberOfDaughters() >= 2)
        {
            int apid = abs(i->daughter(0)->pdgId());
            if (apid == 11) mc_class = 11;
            if (apid == 13) mc_class = 12;
            if (apid == 15) mc_class = 13;
            break;

        }
    }
}

void HeavyNuEvent::regularize()
{
    switch(mode)
    {
        case HNUE:
            if(nLeptons >= 1)
            {
                l1 = new pat::GenericParticle(e1);
                l1pt = hnu::getElectronEt(e1, false);
                l1eta = hnu::getElectronSCEta(e1);
                l1phi = hnu::getElectronSCPhi(e1);
            }
            if(nLeptons >= 2)
            {
                l2 = new pat::GenericParticle(e2);
                l2pt = hnu::getElectronEt(e2, false);
                l2eta = hnu::getElectronSCEta(e2);
                l2phi = hnu::getElectronSCPhi(e2);
            }
            break;
        case TOP:
            if(nLeptons >= 2)
            {
                if(mu1.pt() > hnu::getElectronEt(e1, false))
                {
                    l1 = new pat::GenericParticle(mu1);
                    l1pt = mu1.pt();
                    l1eta = mu1.eta();
                    l1phi = mu1.phi();

                    l2 = new pat::GenericParticle(e1);
                    l2pt = hnu::getElectronEt(e1, false);
                    l2eta = hnu::getElectronSCEta(e1);
                    l2phi = hnu::getElectronSCPhi(e1);
                }
                else
                {
                    l2 = new pat::GenericParticle(mu1);
                    l2pt = mu1.pt();
                    l2eta = mu1.eta();
                    l2phi = mu1.phi();

                    l1 = new pat::GenericParticle(e1);
                    l1pt = hnu::getElectronEt(e1, false);
                    l1eta = hnu::getElectronSCEta(e1);
                    l1phi = hnu::getElectronSCPhi(e1);
                }
            }
            break;
        case HNUMU:
            if(nLeptons >= 1)
            {
                l1 = new pat::GenericParticle(mu1);
                l1pt = mu1.pt();
                l1eta = mu1.eta();
                l1phi = mu1.phi();
            }
            if(nLeptons >= 2)
            {
                l2 = new pat::GenericParticle(mu2);
                l2pt = mu2.pt();
                l2eta = mu2.eta();
                l2phi = mu2.phi();
            }
            break;
        case TAUX:
            if(nLeptons >= 1)
            {
                l1 = NULL;
                l1pt = tl1.Pt();
                l1eta = tl1.Eta();
                l1phi = tl1.Phi();
            }
            if(nLeptons >= 2)
            {
                l2 = NULL;
                l2pt = tl2.Pt();
                l2eta = tl2.Eta();
                l2phi = tl2.Phi();
            }
            break;
        case QCD:
        case CLO:
            break;
    }
}
