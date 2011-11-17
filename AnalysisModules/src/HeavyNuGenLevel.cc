// -*- C++ -*-
//
// Package:    HeavyNuGenLevel
// Class:      HeavyNuGenLevel
// 
/**\class HeavyNuGenLevel HeavyNuGenLevel.cc JetCheck/HeavyNuGenLevel/src/HeavyNuGenLevel.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Alexander Gude
//         Created:  Thu May 12 11:15:22 CDT 2011
// $Id: HeavyNuGenLevel.cc,v 1.4 2011/11/15 20:15:54 mansj Exp $
//
//


// system include files
#include <memory>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HeavyNu/AnalysisModules/src/HeavyNuCommon.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "LHAPDF/LHAPDF.h"


#include "Math/VectorUtil.h"

#include <iostream>

//
// class declaration
//

class HeavyNuGenLevel : public edm::EDAnalyzer
{
public:
    explicit HeavyNuGenLevel(const edm::ParameterSet&);
    ~HeavyNuGenLevel();

private:
    virtual void beginJob();
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();

  double getWeight(float Q, int id1, float x1, int id2, float x2,
		   const reco::Particle::LorentzVector& l1, const reco::Particle::LorentzVector& l2);


    struct CutStruct
    {
        double minJetPT;
        double minL1PT;
        double minL2PT;
        double minLJdR;
        double minLLMass;
        double min4objMass;
    } cuts;

    struct HistStruct
    {
      void book(TFileDirectory *, const std::string&, float* p_weight);
      void fill(const reco::Particle::LorentzVector& l1, const reco::Particle::LorentzVector& l2,
                reco::GenJetCollection::const_iterator j1, reco::GenJetCollection::const_iterator j2);
      
      TH1 *cutProgress;
      TH1 *ptl1, *ptl2;
      TH1 *ptj1, *ptj2;
      TH1 *ljdR, *min_ljdR;
      TH1 *mll, *m4obj;
      float* p_weight;
    } hists, cut0, cut1, cut2, cut3, cut4, cut5;

  float evt_weight;
  bool doPDFreweight_, pdfReweightAddZmass_;
  std::string pdfReweightBaseName, pdfReweightTargetName;
  int pdfReweightBaseId, pdfReweightTargetId;

    // ----------member data ---------------------------
};

void HeavyNuGenLevel::HistStruct::book(TFileDirectory *td, const std::string& post, float* pw)
{
  p_weight=pw;
    std::string t; // histogram title string;

    TH1::SetDefaultSumw2();

    //t = "Progress through cut series : " + post;
    //cutProgress = td->make<TH1F > ("cutProgress", t.c_str(), 20, -0.5, 19.5);
    t = "Lepton 1 PT : " + post;
    ptl1 = td->make<TH1F > ("ptl1", t.c_str(), 50, 0, 500);
    t = "Lepton 2 PT : " + post;
    ptl2 = td->make<TH1F > ("ptl2", t.c_str(), 50, 0, 400);
    t = "Jet 1 PT : " + post;
    ptj1 = td->make<TH1F > ("ptj1", t.c_str(), 50, 0, 500);
    t = "Jet 2 PT : " + post;
    ptj2 = td->make<TH1F > ("ptj2", t.c_str(), 50, 0, 500);

    t = "Lepton-Jet dR : " + post;
    ljdR = td->make<TH1F > ("ljdR", t.c_str(), 50, 0, 5);
    t = "Minimum Lepton-Jet dR : " + post;
    min_ljdR = td->make<TH1F > ("min_ljdR", t.c_str(), 50, 0, 5);

    t = "MLL : " + post;
    mll = td->make<TH1F > ("mll", t.c_str(), 150, 0, 1500);
    t = "MWR : " + post;
    m4obj = td->make<TH1F > ("m4obj", t.c_str(), 250, 0, 2500);
}

void HeavyNuGenLevel::HistStruct::fill(const reco::Particle::LorentzVector& l1, const reco::Particle::LorentzVector& l2,
        reco::GenJetCollection::const_iterator j1, reco::GenJetCollection::const_iterator j2)
{
  ptl1->Fill(l1.pt(),*p_weight);
  ptl2->Fill(l2.pt(),*p_weight);
  ptj1->Fill(j1->pt(),*p_weight);
  ptj2->Fill(j2->pt(),*p_weight);
  
  double dR11 = deltaR(l1.eta(), l1.phi(),
		       j1->eta(), j1->phi());
  double dR12 = deltaR(l1.eta(), l1.phi(),
		       j2->eta(), j2->phi());
  double dR21 = deltaR(l2.eta(), l2.phi(),
		       j1->eta(), j1->phi());
  double dR22 = deltaR(l2.eta(), l2.phi(),
		       j2->eta(), j2->phi());
  
  ljdR->Fill(dR11,*p_weight);
  ljdR->Fill(dR12,*p_weight);
  ljdR->Fill(dR21,*p_weight);
  ljdR->Fill(dR22,*p_weight);
  min_ljdR->Fill(std::min(std::min(dR11, dR12), std::min(dR21, dR22)),*p_weight);
  
  
    reco::Particle::LorentzVector ll = l1 + l2;

    mll->Fill(ll.M(),*p_weight);

    reco::Particle::LorentzVector jj = j1->p4() + j2->p4();
    reco::Particle::LorentzVector wr = jj + ll;

    m4obj->Fill(wr.M(),*p_weight);
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

HeavyNuGenLevel::HeavyNuGenLevel(const edm::ParameterSet& iConfig) :
  doPDFreweight_(iConfig.getUntrackedParameter<bool>("doPDFReweight",false))
{

  if (doPDFreweight_) {
    pdfReweightBaseName=iConfig.getUntrackedParameter<std::string>("pdfReweightBaseName");
    pdfReweightTargetName=iConfig.getUntrackedParameter<std::string>("pdfReweightTargetName");
    pdfReweightBaseId=iConfig.getUntrackedParameter<int>("pdfReweightBaseId",0);
    pdfReweightTargetId=iConfig.getUntrackedParameter<int>("pdfReweightTargetId",0);
    pdfReweightAddZmass_=iConfig.getUntrackedParameter<bool>("pdfReweightAddZMass",true); // fix POWHEG bug 
    std::cout << "PDF Reweighting from " << pdfReweightBaseName << ":" << pdfReweightBaseId
	      << " to " << pdfReweightTargetName << ":" << pdfReweightTargetId
	      << " AddZMass = " << pdfReweightAddZmass_ << std::endl;
  }

    //now do what ever initialization is needed
    cuts.minJetPT = iConfig.getParameter<double>("minJetPt");
    cuts.minL1PT = iConfig.getParameter<double>("minMu1pt");
    cuts.minL2PT = iConfig.getParameter<double>("minMu2pt");
    cuts.minLJdR = iConfig.getParameter<double>("minMuonJetdR");
    cuts.minLLMass = iConfig.getParameter<double>("minMuMuMass");
    cuts.min4objMass = iConfig.getParameter<double>("min4objMass");

    evt_weight=1.0;

    edm::Service<TFileService> fs;

    hists.cutProgress = fs->make<TH1F > ("cutProgress", "Progress through cut series", 21, -1.5, 19.5);
    hists.ptl1 = fs->make<TH1F > ("ptl1", "Lepton 1 PT", 50, 0, 500);
    hists.ptl2 = fs->make<TH1F > ("ptl2", "Lepton 2 PT", 50, 0, 400);
    hists.ptj1 = fs->make<TH1F > ("ptj1", "Jet 1 PT", 50, 0, 500);
    hists.ptj2 = fs->make<TH1F > ("ptj2", "Jet 2 PT", 50, 0, 500);

    hists.ljdR = fs->make<TH1F > ("ljdR", "Lepton-Jet dR", 50, 0, 5);
    hists.min_ljdR = fs->make<TH1F > ("min_ljdR", "Minimum Lepton-Jet dR", 50, 0, 5);

    hists.mll = fs->make<TH1F > ("mll", "MLL", 50, 0, 700);
    hists.m4obj = fs->make<TH1F > ("m4obj", "MWR", 50, 0, 2500);

    cut0.book(new TFileDirectory(fs->mkdir("cut0_none")), "(no cuts)",&evt_weight);
    //cut1.book(new TFileDirectory(fs->mkdir("cutX_LLpt")), "(dileptons with ptcuts:X)",&evt_weight);
    cut2.book(new TFileDirectory(fs->mkdir("cut1_LLJJpt")), "(4objects with ptcuts:1)",&evt_weight);
    cut3.book(new TFileDirectory(fs->mkdir("cut4_Mu1HighPt")), "(Mu1 High pt cut:4)",&evt_weight);
    cut4.book(new TFileDirectory(fs->mkdir("cut5_diLmass")), "(mumu mass cut:5a)",&evt_weight);
    cut5.book(new TFileDirectory(fs->mkdir("cut6_mWRmass")), "(mumujj mass cut:6)",&evt_weight);

}

HeavyNuGenLevel::~HeavyNuGenLevel()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // Need recoCaloJets_ak5CaloJets_RECO
    // recoJPTJets_JetPlusTrackZSPCorJetAntiKt5_RECO
}


//
// member functions
//

// ------------ method called to for each event  ------------

void HeavyNuGenLevel::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //using namespace edm;
    //using namespace ROOT::Math::VectorUtil; // DeltaR

    edm::Handle<edm::HepMCProduct> hepMCEvt;
    edm::Handle<reco::GenJetCollection> genJets;
    hists.cutProgress->Fill(-1);

    iEvent.getByLabel("ak5GenJetsNoMuNoNu", genJets);


    reco::Particle::LorentzVector l1,l2;
    reco::GenJetCollection::const_iterator ji, jet1 = genJets->end(), jet2 = genJets->end();
    float Q, x1, x2;
    int id1, id2;

    if (iEvent.getByLabel("generator", hepMCEvt)) {

      const HepMC::GenEvent* genEvt = hepMCEvt->GetEvent();
      const HepMC::GenEvent* genE = genEvt;

      HepMC::GenEvent::vertex_const_iterator vtex;
      HepMC::GenVertex::particles_out_const_iterator Pout;
      HepMC::GenParticle* theLep = 0;
      HepMC::GenParticle* theLep2 = 0;

      for(vtex = genE->vertices_begin(); vtex != genE->vertices_end(); vtex++)
	{
	  for(Pout = (*vtex)->particles_out_const_begin(); Pout != (*vtex)->particles_out_const_end(); Pout++)
	    {
	      if(abs((*Pout)->pdg_id()) == 13 && (*Pout)->status() == 1)
		{
		  if(theLep == 0 || theLep->momentum().perp()<(*Pout)->momentum().perp())
		    {
		      theLep2 = theLep;
		      theLep = *Pout;
		    }
		  else if(theLep2 == 0 || theLep2->momentum().perp()<(*Pout)->momentum().perp())
		    {
		      theLep2 = *Pout;
		    }
		}
	    }
	}
      if(theLep2 == 0)
	{
	  //std::cout << "Got less than two!\n";
	  return;
	}
      l1=reco::Particle::LorentzVector(theLep->momentum().px(), theLep->momentum().py(), theLep->momentum().pz(), theLep->momentum().e());
      l2=reco::Particle::LorentzVector(theLep2->momentum().px(), theLep2->momentum().py(), theLep2->momentum().pz(), theLep2->momentum().e());

      const HepMC::PdfInfo* pdfstuff = genEvt->pdf_info(); 
      Q = pdfstuff->scalePDF();

      id1 = pdfstuff->id1();
      x1 = pdfstuff->x1();
      id2 = pdfstuff->id2();
      x2 = pdfstuff->x2();
      
    } else {
      edm::Handle<reco::GenParticleCollection> gpp;
      iEvent.getByLabel("genParticles",gpp); 


      reco::GenParticleCollection::const_iterator i;
      int nfound=0;
      for (i=gpp->begin(); i!=gpp->end(); i++) {
	if (i->status()!=1 || abs(i->pdgId())!=13) continue; // only final-state muons need apply
	if (i->pt()>l1.pt()) {
	  l2=l1;
	  l1=i->p4();
	  nfound++;
	} else if (i->pt()>l2.pt()) {
	  l2=i->p4();
	  nfound++;
	}
      }
      if (nfound<2) return;

      edm::Handle<GenEventInfoProduct> geip;
      iEvent.getByLabel("generator",geip);
      
      Q=geip->pdf()->scalePDF;
      id1=geip->pdf()->id.first;
      id2=geip->pdf()->id.second;
      x1=geip->pdf()->x.first;
      x2=geip->pdf()->x.second;
    }
    evt_weight=getWeight(Q,id1,x1,id2,x2,l1,l2);



    hists.cutProgress->Fill(0.0,evt_weight);    

    for(ji = genJets->begin(); ji != genJets->end(); ji++)
    {
        if(fabs(ji->eta()) > 2.5) continue;
        if(jet1 == genJets->end() || jet1->pt() < ji->pt())
        {
            jet2 = jet1;
            jet1 = ji;
        }
        else if(jet2 == genJets->end() || jet2->pt() < ji->pt())
        {
            jet2 = ji;
        }
    }


    if(jet1 != genJets->end() && jet2 != genJets->end()) cut0.fill(l1, l2, jet1, jet2);

    hists.cutProgress->Fill(1,evt_weight);

    hists.ptl1->Fill(l1.pt(),evt_weight);
    hists.ptl2->Fill(l2.pt(),evt_weight);

    if(l1.pt() > cuts.minL2PT &&
            l2.pt() > cuts.minL2PT &&
            fabs(l1.eta()) < 2.4 &&
            fabs(l2.eta()) < 2.4 &&
            (fabs(l1.eta()) < 2.1 || fabs(l2.eta()) < 2.1)
            )
    {

        hists.cutProgress->Fill(2,evt_weight);
        //cut1.fill(l1, l2, NULL, NULL);

        if(jet1 != genJets->end() && jet2 != genJets->end())
        {
            hists.ptj1->Fill(jet1->pt(),evt_weight);
            hists.ptj2->Fill(jet2->pt(),evt_weight);

            if(jet1->pt() > cuts.minJetPT && jet2->pt() > cuts.minJetPT)
            {
                hists.cutProgress->Fill(3,evt_weight);
                cut2.fill(l1, l2, jet1, jet2);

                if(l1.pt() < cuts.minL1PT) return;

                double dR11 = deltaR(l1.eta(), l1.phi(),
                        jet1->eta(), jet1->phi());
                double dR12 = deltaR(l1.eta(), l1.phi(),
                        jet2->eta(), jet2->phi());
                double dR21 = deltaR(l2.eta(), l2.phi(),
                        jet1->eta(), jet1->phi());
                double dR22 = deltaR(l2.eta(), l2.phi(),
                        jet2->eta(), jet2->phi());

                hists.ljdR->Fill(dR11,evt_weight);
                hists.ljdR->Fill(dR12,evt_weight);
                hists.ljdR->Fill(dR21,evt_weight);
                hists.ljdR->Fill(dR22,evt_weight);
                //hists.ljdR->Fill(dR11);
                hists.min_ljdR->Fill(std::min(std::min(dR11, dR12), std::min(dR21, dR22)),evt_weight);

                if(dR11 > cuts.minLJdR && dR12 > cuts.minLJdR &&
                        dR21 > cuts.minLJdR && dR22 > cuts.minLJdR)
                {
                    hists.cutProgress->Fill(4);
                    cut3.fill(l1, l2, jet1, jet2);

                    reco::Particle::LorentzVector ll = l1+l2;

                    hists.mll->Fill(ll.M(),evt_weight);

                    if(ll.M() > cuts.minLLMass)
                    {
                        hists.cutProgress->Fill(5);
                        cut4.fill(l1, l2, jet1, jet2);

                        reco::Particle::LorentzVector jj = jet1->p4() + jet2->p4();
                        reco::Particle::LorentzVector wr = jj + ll;

                        hists.m4obj->Fill(wr.M(),evt_weight);

                        if(wr.M() > cuts.min4objMass)
                        {
                            hists.cutProgress->Fill(6);
                            cut5.fill(l1, l2, jet1, jet2);
                        }

                    }


                }
            }

        }
    }

}


// ------------ method called once each job just before starting event loop  ------------

void
HeavyNuGenLevel::beginJob()
{
  if (doPDFreweight_) {
    LHAPDF::initPDFSet(1,pdfReweightBaseName);
    LHAPDF::initPDFSet(2,pdfReweightTargetName);
  }
}

// ------------ method called once each job just after ending the event loop  ------------

void
HeavyNuGenLevel::endJob()
{
}


double HeavyNuGenLevel::getWeight(float Q, int id1, float x1, int id2, float x2,
				  const reco::Particle::LorentzVector& l1, const reco::Particle::LorentzVector& l2) {

  if (!doPDFreweight_) return 1.0;

  if (pdfReweightAddZmass_) {
    reco::Particle::LorentzVector ll = l1 + l2;

    Q=sqrt(Q*Q+ll.M2());
  }

  LHAPDF::usePDFMember(1,pdfReweightBaseId);
  double pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
  double pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  
  LHAPDF::usePDFMember(2,pdfReweightTargetId);
  double newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
  double newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
  
  double w=(newpdf1/pdf1*newpdf2/pdf2);

  //  printf("My weight is %f\n",w);

  return w;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(HeavyNuGenLevel);
