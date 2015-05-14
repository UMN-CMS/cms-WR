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

class FlatTTree : public edm::EDAnalyzer {
public:
  explicit FlatTTree(const edm::ParameterSet&);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  const edm::InputTag muons_src;
  const edm::InputTag jets_src;
  const edm::InputTag genparticles_src;
  const edm::InputTag genjets_src;
  const edm::InputTag primary_vertex_src;

  const float Mlljj_cut;
  const float Mll_cut;
  const float matching_dR;
  const bool matching;

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
      leading_mu_pt = leading_mu_eta = leading_mu_phi = subleading_mu_pt = subleading_mu_eta = subleading_mu_phi = 0;
      leading_jet_pt = leading_jet_eta = leading_jet_phi = subleading_jet_pt = subleading_jet_eta = subleading_jet_phi = 0;
      gen_leading_mu_pt = gen_leading_mu_eta = gen_leading_mu_phi = gen_subleading_mu_pt = gen_subleading_mu_eta = gen_subleading_mu_phi = 0;
      gen_leading_jet_pt = gen_leading_jet_eta = gen_leading_jet_phi = gen_subleading_jet_pt = gen_subleading_jet_eta = gen_subleading_jet_phi = 0;
      dilepton_mass = gen_dilepton_mass = Mlljj = gen_Mlljj = 0;
      dR_leadLepton_leadJet = dR_leadLepton_subleadJet = dR_subleadLepton_leadJet = dR_subleadLepton_subleadJet = 0;
      gen_dR_leadLepton_leadJet = gen_dR_leadLepton_subleadJet = gen_dR_subleadLepton_leadJet = gen_dR_subleadLepton_subleadJet = 0;
      //pv_x.clear();
      
    }
  };

  TTree* tree;
  tree_t nt;
};

FlatTTree::FlatTTree(const edm::ParameterSet& cfg)
  : muons_src(cfg.getParameter<edm::InputTag>("muons_src")),
    jets_src(cfg.getParameter<edm::InputTag>("jets_src")),
    genparticles_src(cfg.getParameter<edm::InputTag>("genparticles_src")),
    genjets_src(cfg.getParameter<edm::InputTag>("genjets_src")),
    primary_vertex_src(cfg.getParameter<edm::InputTag>("primary_vertex_src")),
    Mlljj_cut(cfg.getParameter<double>("Mlljj_cut")),
    Mll_cut(cfg.getParameter<double>("Mll_cut")),
    matching_dR(cfg.getParameter<double>("matching_dR")),
    matching(cfg.getParameter<bool>("matching"))
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

void FlatTTree::analyze(const edm::Event& event, const edm::EventSetup&) {
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

  // Find the two muons from the WR and the neutrino
  // and the two quarks from the neutrino decay
  for(auto gen : *gen_particles){
    if(abs(gen.pdgId()) == 13 && gen.mother() != 0){
      if(abs(gen.mother()->pdgId()) == 9900024 ) gmuons[0] = gen;
      else if(abs(gen.mother()->pdgId()) == 9900014 ) gmuons[1] = gen;
    }
    if(abs(gen.pdgId()) < 7 && gen.mother() != 0){
      if(abs(gen.mother()->pdgId()) == 9900014)	gen_quarks.push_back(gen);	
    }
  }
  // Match the quarks to genjets with a dR cut
  for(auto gj: *gen_jets){
    if(deltaR(gj,gen_quarks[0]) < 0.1) gjets[0] = gj;
    else if(deltaR(gj,gen_quarks[1]) < 0.1) gjets[1] = gj;
  }
  
  if(gmuons.size() > 0){
    nt.gen_leading_mu_pt = gmuons[0].pt();
    nt.gen_leading_mu_eta = gmuons[0].eta();
    nt.gen_leading_mu_phi = gmuons[0].phi();
    if(gjets.size() > 0) nt.gen_dR_leadLepton_leadJet = deltaR(gmuons[0],gjets[0]);
    if(gjets.size() > 1) nt.gen_dR_leadLepton_subleadJet = deltaR(gmuons[0],gjets[1]);
    if(gmuons.size() > 1){
      nt.gen_subleading_mu_pt = gmuons[1].pt();
      nt.gen_subleading_mu_eta = gmuons[1].eta();
      nt.gen_subleading_mu_phi = gmuons[1].phi();
      nt.gen_dilepton_mass = (gmuons[0].p4() + gmuons[1].p4()).M();
      if(gjets.size() > 0) nt.gen_dR_subleadLepton_leadJet = deltaR(gmuons[1],gjets[0]);
      if(gjets.size() > 1){
	nt.gen_dR_subleadLepton_subleadJet = deltaR(gmuons[1],gjets[1]);
	nt.gen_Mlljj = (gmuons[0].p4() + gmuons[1].p4() + gjets[0].p4() + gjets[1].p4()).M();
      }
    }
  }
  if(gjets.size() > 0){
    nt.gen_leading_jet_pt = gjets[0].pt();
    nt.gen_leading_jet_eta = gjets[0].eta();
    nt.gen_leading_jet_phi = gjets[0].phi();
    if(gjets.size() > 1){
      nt.gen_subleading_jet_pt = gjets[1].pt();
      nt.gen_subleading_jet_eta = gjets[1].eta();
      nt.gen_subleading_jet_phi = gjets[1].phi();
    }
  }
  
  std::vector<pat::Muon> mus(2);
  std::vector<pat::Jet> js(2);
  
  // Find the reco muons, either by matching to the gen muons
  // or just the first two muons in the collection
  if(matching){
    for(auto mu : *muons){      
      if(gmuons.size() > 0){
	if(deltaR(gmuons[0],mu) < matching_dR) mus[0] = mu;      
      }
      if(gmuons.size() > 1){
	if(deltaR(gmuons[1],mu) < matching_dR) mus[1] = mu;      
      }
    }
  }
  else {
    if(muons->size() > 0) mus[0] = muons->at(0);
    if(muons->size() > 1) mus[1] = muons->at(1);
  }

  // Find the reco jets, either by matching to the gen jets
  // or just the first two jets in the collection
  if(matching){
    for(auto j : *jets){
      if(gjets.size() > 0){
	if(deltaR(gjets[0],j) < matching_dR) js[0] = j;      
      }
      if(gjets.size() > 1){
	if(deltaR(gjets[1],j) < matching_dR) js[1] = j;      
      }
    }
  }
  else {
    if(jets->size() > 0) js[0] = jets->at(0);
    if(jets->size() > 1) js[1] = jets->at(1);
  }

  if(mus.size() > 0){
    nt.leading_mu_pt = mus[0].pt();
    nt.leading_mu_eta = mus[0].eta();
    nt.leading_mu_phi = mus[0].phi();
    if(gjets.size() > 0) nt.dR_leadLepton_leadJet = deltaR(mus[0],js[0]);
    if(gjets.size() > 1) nt.dR_leadLepton_subleadJet = deltaR(mus[0],js[1]);
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

  tree->Fill();
}

DEFINE_FWK_MODULE(FlatTTree);
