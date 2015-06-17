#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"


class Matching : public edm::EDAnalyzer {
public:
  explicit Matching(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&);

private:
  const edm::InputTag gen_src;
  const edm::InputTag vertex_src;
  const edm::InputTag gen_jet_src;
  const edm::InputTag electron_src;
  const edm::InputTag muon_src;
  const edm::InputTag jet_src;

  TH1F* mu_pre0;
  

};

Matching::Matching(const edm::ParameterSet& cfg)
  : gen_src(cfg.getParameter<edm::InputTag>("gen_src")),
    vertex_src(cfg.getParameter<edm::InputTag>("vertex_src")),
    gen_jet_src(cfg.getParameter<edm::InputTag>("gen_jet_src")),
    electron_src(cfg.getParameter<edm::InputTag>("electron_src")),
    muon_src(cfg.getParameter<edm::InputTag>("muon_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src"))
{
  edm::Service<TFileService> fs;
  mu_pre0 = fs->make<TH1F>("mu_pre0", "", 4, 0, 3);

    
}


void Matching::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<pat::ElectronCollection> electrons;
  event.getByLabel(electron_src, electrons);
  edm::Handle<pat::MuonCollection> muons;
  event.getByLabel(muon_src, muons);
  edm::Handle<reco::VertexCollection> primary_vertex;
  event.getByLabel(vertex_src, primary_vertex);
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jet_src, jets);
  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(gen_src, gen_particles);
  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel(gen_jet_src, gen_jets);

  std::vector<pat::Electron> eles;
  std::vector<pat::Muon> mus(2);
  std::vector<reco::GenParticle> WR_daughters;
  std::vector<reco::GenParticle> gen_mus;
  std::vector<reco::GenParticle> gjets;
  std::vector<reco::GenParticle> gmus_W;  
  std::vector<pat::Jet> js(2);

  edm::Handle<reco::VertexCollection> vertices;
  event.getByLabel("offlineSlimmedPrimaryVertices", vertices);
  cout<<"event"<<endl;
  for(auto g : *gen_particles){
    if(abs(g.pdgId()) == 9900012 || abs(g.pdgId()) == 9900016)
      cout<<"flavor flav"<<endl;
    if(abs(g.pdgId()) < 15 && (abs(g.mother()->pdgId()) == 9900024|| abs(g.mother()->pdgId()) == 9900014 ) ){
      //cout<<g.pdgId()<<" "<<g.mother()->pdgId()<<endl;
      WR_daughters.push_back(g);
    }    
  }

  if(WR_daughters.size() == 4){
    if(abs(WR_daughters[0].pdgId()) != 13)
      cout<<"BREAK"<<endl;

    double tmp_dR_0 = 0.4;
    double tmp_dR_1 = 0.4;
    for(auto m: *muons){
      if(deltaR(WR_daughters[0],m) < tmp_dR_0){
	mus[0] = m;
	tmp_dR_0 = deltaR(WR_daughters[0],m);
      }
      if(deltaR(WR_daughters[1],m) < tmp_dR_1){
	mus[1] = m;
	tmp_dR_1 = deltaR(WR_daughters[1],m);
      }
    }
    double tmp_dR_2 = 0.4;
    double tmp_dR_3 = 0.4;
    for(auto j: *jets){
      if(deltaR(WR_daughters[2],j) < tmp_dR_2){
	js[0] = j;
	tmp_dR_2 = deltaR(WR_daughters[2],j);
      }
      if(deltaR(WR_daughters[3],j) < tmp_dR_3){
	js[1] = j;
	tmp_dR_3 = deltaR(WR_daughters[3],j);
      }
    }
    cout<<mus[0].pt()<<" "<<js[0].pt()<<" "<<mus[1].pt()<<" "<<js[1].pt()<<endl;
  }

  if(mus[0].pt() == 0)
    mu_pre0->Fill(0);
  if(mus[1].pt() == 0)
    mu_pre0->Fill(1);
  if(js[0].pt() == 0)
    mu_pre0->Fill(2);
  if(js[1].pt() == 0)
    mu_pre0->Fill(3);

}

DEFINE_FWK_MODULE(Matching);
