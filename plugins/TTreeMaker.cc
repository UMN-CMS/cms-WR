#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

class TTreeMaker : public edm::EDAnalyzer {
public:
  explicit TTreeMaker(const edm::ParameterSet&);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  const edm::InputTag muons_src;
  const edm::InputTag jets_src;
  const edm::InputTag genparticles_src;
  const edm::InputTag genjets_src;
  const edm::InputTag primary_vertex_src;

  const float Mlljj_cut;
  const float Mll_cut;
  const float isolation_dR;
  const float leading_lepton_pt_cut;
  const float lepton_eta_cut;
  const float subleading_lepton_pt_cut;
  const float leading_jet_pt_cut;
  const float jet_eta_cut;
  const float subleading_jet_pt_cut;

  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned event;
    
    float leading_mu_pt;
    float subleading_mu_pt;
    float leading_jet_pt;
    float subleading_jet_pt;

    float leading_mu_eta;
    float subleading_mu_eta;
    float leading_jet_eta;
    float subleading_jet_eta;

    float leading_mu_phi;
    float subleading_mu_phi;
    float leading_jet_phi;
    float subleading_jet_phi;

    float dilepton_mass;
    float Mlljj;
    float dR_leadLepton_leadJet;
    float dR_leadLepton_subleadJet;
    float dR_subleadLepton_leadJet;
    float dR_subleadLepton_subleadJet;


    // Generator level quantities
    // Leading lepton corresponds to the WR daughter
    // Subleading muon corresponds to the NR daughter

    float gen_leading_mu_pt;
    float gen_subleading_mu_pt;
    float gen_leading_jet_pt;
    float gen_subleading_jet_pt;

    float gen_leading_mu_eta;
    float gen_subleading_mu_eta;
    float gen_leading_jet_eta;
    float gen_subleading_jet_eta;

    float gen_leading_mu_phi;
    float gen_subleading_mu_phi;
    float gen_leading_jet_phi;
    float gen_subleading_jet_phi;

    float gen_dilepton_mass;
    float gen_Mlljj;
    float gen_dR_leadLepton_leadJet;
    float gen_dR_leadLepton_subleadJet;
    float gen_dR_subleadLepton_leadJet;
    float gen_dR_subleadLepton_subleadJet;
    
    //std::vector<float> pv_x;
    
    tree_t() { clear(); }

    void clear() {
      run = lumi = event = 0;
      leading_mu_pt = leading_mu_eta = leading_mu_phi = subleading_mu_pt = subleading_mu_eta = subleading_mu_phi = -999;
      leading_jet_pt = leading_jet_eta = leading_jet_phi = subleading_jet_pt = subleading_jet_eta = subleading_jet_phi = -999;
      gen_leading_mu_pt = gen_leading_mu_eta = gen_leading_mu_phi = gen_subleading_mu_pt = gen_subleading_mu_eta = gen_subleading_mu_phi = -999;
      gen_leading_jet_pt = gen_leading_jet_eta = gen_leading_jet_phi = gen_subleading_jet_pt = gen_subleading_jet_eta = gen_subleading_jet_phi = -999;
      dilepton_mass = gen_dilepton_mass = Mlljj = gen_Mlljj = -999;
      dR_leadLepton_leadJet = dR_leadLepton_subleadJet = dR_subleadLepton_leadJet = dR_subleadLepton_subleadJet = 9;
      gen_dR_leadLepton_leadJet = gen_dR_leadLepton_subleadJet = gen_dR_subleadLepton_leadJet = gen_dR_subleadLepton_subleadJet = 9;
      //pv_x.clear();
      
    }
  };

  TTree* tree;
  tree_t nt;
};

TTreeMaker::TTreeMaker(const edm::ParameterSet& cfg)
  : muons_src(cfg.getParameter<edm::InputTag>("muons_src")),
    jets_src(cfg.getParameter<edm::InputTag>("jets_src")),
    genparticles_src(cfg.getParameter<edm::InputTag>("genparticles_src")),
    genjets_src(cfg.getParameter<edm::InputTag>("genjets_src")),
    primary_vertex_src(cfg.getParameter<edm::InputTag>("primary_vertex_src")),
    Mlljj_cut(cfg.getParameter<double>("Mlljj_cut")),
    Mll_cut(cfg.getParameter<double>("Mll_cut")),
    isolation_dR(cfg.getParameter<double>("isolation_dR")),
    leading_lepton_pt_cut(cfg.getParameter<double>("leading_lepton_pt_cut")),
    lepton_eta_cut(cfg.getParameter<double>("lepton_eta_cut")),
    subleading_lepton_pt_cut(cfg.getParameter<double>("subleading_lepton_pt_cut")),
    leading_jet_pt_cut(cfg.getParameter<double>("leading_jet_pt_cut")),
    jet_eta_cut(cfg.getParameter<double>("jet_eta_cut")),
    subleading_jet_pt_cut(cfg.getParameter<double>("subleading_jet_pt_cut"))
{
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &nt.run, "run/i");
  tree->Branch("lumi", &nt.lumi, "lumi/i");
  tree->Branch("event", &nt.event, "event/i");

  tree->Branch("leading_mu_pt", &nt.leading_mu_pt, "leading_mu_pt/F");
  tree->Branch("leading_mu_eta", &nt.leading_mu_eta, "leading_mu_eta/F");
  tree->Branch("leading_mu_phi", &nt.leading_mu_phi, "leading_mu_phi/F");
  tree->Branch("subleading_mu_pt", &nt.subleading_mu_pt, "subleading_mu_pt/F");
  tree->Branch("subleading_mu_eta", &nt.subleading_mu_eta, "subleading_mu_eta/F");
  tree->Branch("subleading_mu_phi", &nt.subleading_mu_phi, "subleading_mu_phi/F");
  tree->Branch("leading_jet_pt", &nt.leading_jet_pt, "leading_jet_pt/F");
  tree->Branch("leading_jet_eta", &nt.leading_jet_eta, "leading_jet_eta/F");
  tree->Branch("leading_jet_phi", &nt.leading_jet_phi, "leading_jet_phi/F");
  tree->Branch("subleading_jet_pt", &nt.subleading_jet_pt, "subleading_jet_pt/F");
  tree->Branch("subleading_jet_eta", &nt.subleading_jet_eta, "subleading_jet_eta/F");
  tree->Branch("subleading_jet_phi", &nt.subleading_jet_phi, "subleading_jet_phi/F");

  tree->Branch("gen_leading_mu_pt", &nt.gen_leading_mu_pt, "gen_leading_mu_pt/F");
  tree->Branch("gen_leading_mu_eta", &nt.gen_leading_mu_eta, "gen_leading_mu_eta/F");
  tree->Branch("gen_leading_mu_phi", &nt.gen_leading_mu_phi, "gen_leading_mu_phi/F");
  tree->Branch("gen_subleading_mu_pt", &nt.gen_subleading_mu_pt, "gen_subleading_mu_pt/F");
  tree->Branch("gen_subleading_mu_eta", &nt.gen_subleading_mu_eta, "gen_subleading_mu_eta/F");
  tree->Branch("gen_subleading_mu_phi", &nt.gen_subleading_mu_phi, "gen_subleading_mu_phi/F");
  tree->Branch("gen_leading_jet_pt", &nt.gen_leading_jet_pt, "gen_leading_jet_pt/F");
  tree->Branch("gen_leading_jet_eta", &nt.gen_leading_jet_eta, "gen_leading_jet_eta/F");
  tree->Branch("gen_leading_jet_phi", &nt.gen_leading_jet_phi, "gen_leading_jet_phi/F");
  tree->Branch("gen_subleading_jet_pt", &nt.gen_subleading_jet_pt, "gen_subleading_jet_pt/F");
  tree->Branch("gen_subleading_jet_eta", &nt.gen_subleading_jet_eta, "gen_subleading_jet_eta/F");
  tree->Branch("gen_subleading_jet_phi", &nt.gen_subleading_jet_phi, "gen_subleading_jet_phi/F");
  
  tree->Branch("Mlljj", &nt.Mlljj, "Mlljj/F");
  tree->Branch("gen_Mlljj", &nt.gen_Mlljj, "gen_Mlljj/F");
  tree->Branch("dilepton_mass", &nt.dilepton_mass, "dilepton_mass/F");
  tree->Branch("gen_dilepton_mass", &nt.gen_dilepton_mass, "gen_dilepton_mass/F");
  tree->Branch("dR_leadLepton_leadJet", &nt.dR_leadLepton_leadJet, "dR_leadLepton_leadJet/F");
  tree->Branch("dR_leadLepton_subleadJet", &nt.dR_leadLepton_subleadJet, "dR_leadLepton_subleadJet/F");
  tree->Branch("dR_subleadLepton_leadJet", &nt.dR_subleadLepton_leadJet, "dR_subleadLepton_leadJet/F");
  tree->Branch("dR_subleadLepton_subleadJet", &nt.dR_subleadLepton_subleadJet, "dR_subleadLepton_subleadJet/F");
  tree->Branch("gen_dR_leadLepton_leadJet", &nt.gen_dR_leadLepton_leadJet, "gen_dR_leadLepton_leadJet/F");
  tree->Branch("gen_dR_leadLepton_subleadJet", &nt.gen_dR_leadLepton_subleadJet, "gen_dR_leadLepton_subleadJet/F");
  tree->Branch("gen_dR_subleadLepton_leadJet", &nt.gen_dR_subleadLepton_leadJet, "gen_dR_subleadLepton_leadJet/F");
  tree->Branch("gen_dR_subleadLepton_subleadJet", &nt.gen_dR_subleadLepton_subleadJet, "gen_dR_subleadLepton_subleadJet/F");

}

void TTreeMaker::analyze(const edm::Event& event, const edm::EventSetup&) {
  nt.clear();
  nt.run = event.id().run();
  nt.lumi = event.luminosityBlock();
  nt.event = event.id().event();

  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muons_src, muons);
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jets_src, jets);
  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(genparticles_src, gen_particles);
  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel(genjets_src, gen_jets);
  
  std::vector<reco::GenParticle> gmuons(2);
  std::vector<reco::GenParticle> gen_quarks;
  std::vector<reco::GenJet> gjets(2);
  
  // Skip events with a tau NR or an electron NR

  bool double_nu = false;
  for(auto gen : *gen_particles){
    if(abs(gen.pdgId()) == 9900012 || abs(gen.pdgId()) == 9900016)
      double_nu = true;        
  }
  
  std::vector<pat::Muon> mus;
  std::vector<pat::Jet> js;
  edm::Handle<reco::VertexCollection> primary_vertex;
  event.getByLabel("offlineSlimmedPrimaryVertices", primary_vertex);
  
  for(auto mu : *muons){         
    if(mu.pt() > leading_lepton_pt_cut && fabs(mu.eta()) < lepton_eta_cut && mu.isHighPtMuon(primary_vertex->at(0))){
      mus.push_back(mu);	     
    }    
  }
  for(auto j : *jets){  
    if(j.pt() > leading_jet_pt_cut && fabs(j.eta()) < jet_eta_cut){
      bool isolated = true;
      for(auto m : mus){
	if(deltaR(m,j) < isolation_dR) isolated = false;
      }    
      if(isolated) js.push_back(j);  
    }
  }
  
  if(mus.size() > 0){
    nt.leading_mu_pt = mus[0].pt();
    nt.leading_mu_eta = mus[0].eta();
    nt.leading_mu_phi = mus[0].phi();
    if(js.size() > 0) nt.dR_leadLepton_leadJet = deltaR(mus[0],js[0]);
    if(js.size() > 1) nt.dR_leadLepton_subleadJet = deltaR(mus[0],js[1]);
    if(mus.size() > 1){
      nt.subleading_mu_pt = mus[1].pt();
      nt.subleading_mu_eta = mus[1].eta();
      nt.subleading_mu_phi = mus[1].phi();
      nt.dilepton_mass = (mus[0].p4() + mus[1].p4()).M();
      if(js.size() > 0) nt.dR_subleadLepton_leadJet = deltaR(mus[1],js[0]);
      if(js.size() > 1){
	nt.dR_subleadLepton_subleadJet = deltaR(mus[1],js[1]);
	nt.Mlljj = (mus[0].p4() + mus[1].p4() + js[0].p4() + js[1].p4()).M();
      }
    }
  }
  if(js.size() > 0){
    nt.leading_jet_pt = js[0].pt();
    nt.leading_jet_eta = js[0].eta();
    nt.leading_jet_phi = js[0].phi();
    if(js.size() > 1){
      nt.subleading_jet_pt = js[1].pt();
      nt.subleading_jet_eta = js[1].eta();
      nt.subleading_jet_phi = js[1].phi();
    }
  }

  if(!double_nu)
    tree->Fill();
}

DEFINE_FWK_MODULE(TTreeMaker);
