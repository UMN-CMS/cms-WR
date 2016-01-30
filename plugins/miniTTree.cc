#include "TTree.h"
#include "TLorentzVector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/ValueMap.h"

typedef double JECUnc_t;
typedef edm::ValueMap<JECUnc_t> JECUnc_Map;

class miniTTree : public edm::EDAnalyzer {
public:
  explicit miniTTree(const edm::ParameterSet&);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  edm::EDGetToken electronsMiniAODToken_;
  edm::EDGetToken muonsMiniAODToken_;
  edm::EDGetToken jetsMiniAODToken_;

  const std::string jec_unc_src;

  const bool is_mc;

  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned long long event;

    std::vector<TLorentzVector> electrons_p4;
    std::vector<TLorentzVector> muons_p4;
    std::vector<TLorentzVector> jets_p4;
    
    std::vector<float> jet_uncertainty;
    std::vector<float> electron_scale;
    std::vector<float> electron_smearing;

    float PU;
    float weight;
        
    tree_t() { clear(); }

    void clear() {
      run = lumi = event = 0;
    
      electrons_p4.clear();
      muons_p4.clear();
      jets_p4.clear();

      jet_uncertainty.clear();
      electron_scale.clear();
      electron_smearing.clear();

      PU = -999.;
      weight = 0.0;
      
    }
  };

  TTree* tree;
  tree_t nt;
};

miniTTree::miniTTree(const edm::ParameterSet& cfg)
  : jec_unc_src(cfg.getParameter<std::string>("jec_unc_src")),
    is_mc(cfg.getParameter<bool>("is_mc"))
{
  electronsMiniAODToken_   = mayConsume<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons_src"));
  muonsMiniAODToken_ = mayConsume<edm::View<pat::Muon> >(cfg.getParameter<edm::InputTag>("muons_src"));
  jetsMiniAODToken_ = mayConsume<edm::View<pat::Jet> >(cfg.getParameter<edm::InputTag>("jets_src"));

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &nt.run);
  tree->Branch("lumi", &nt.lumi);
  tree->Branch("event", &nt.event);

  tree->Branch("electrons_p4", &nt.electrons_p4,32000,-1);
  tree->Branch("muons_p4", &nt.muons_p4,32000,-1);
  tree->Branch("jets_p4", &nt.jets_p4,32000,-1);

  tree->Branch("jet_uncertainty",&nt.jet_uncertainty);
  tree->Branch("electron_scale",&nt.electron_scale);
  tree->Branch("electron_smearing",&nt.electron_smearing);

  tree->Branch("PU", &nt.PU, "PU/F");
  tree->Branch("weight",&nt.weight,"weight/F");

}

void miniTTree::analyze(const edm::Event& event, const edm::EventSetup&) {
  nt.clear();
  nt.run = event.id().run();
  nt.lumi = event.luminosityBlock();
  nt.event = event.id().event();

  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  event.getByToken(electronsMiniAODToken_,electrons);
  edm::Handle<edm::View<pat::Muon> > muons;
  event.getByToken(muonsMiniAODToken_,muons);
  edm::Handle<edm::View<pat::Jet> > jets;
  event.getByToken(jetsMiniAODToken_,jets);
  edm::Handle<GenEventInfoProduct> evinfo;

  edm::Handle< edm::ValueMap<double> > jec_unc;
  event.getByLabel(jec_unc_src, "JECUncertainty", jec_unc);

  if(is_mc) 
    {
      event.getByLabel("generator", evinfo);
      nt.weight = evinfo->weight();
    }

  nt.PU = 1.0;
  
  for (size_t i = 0; i < electrons->size(); ++i){
    const auto ele = electrons->ptrAt(i);
    TLorentzVector p4;
    p4.SetPtEtaPhiM(ele->pt(),ele->eta(),ele->phi(),ele->mass());
    nt.electrons_p4.push_back(p4);
    nt.electron_scale.push_back(1.0);
    nt.electron_smearing.push_back(1.0);
  }
  for (size_t i = 0; i < muons->size(); ++i){
    const auto mu = muons->ptrAt(i);
    TLorentzVector p4;
    p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
    nt.muons_p4.push_back(p4);
  }
  for (size_t i = 0; i < jets->size(); ++i){
    const auto jet = jets->ptrAt(i);
    TLorentzVector p4;
    p4.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());
    nt.jets_p4.push_back(p4);
    nt.jet_uncertainty.push_back((*jec_unc)[jet]);
  }



  
  tree->Fill();
}

DEFINE_FWK_MODULE(miniTTree);
