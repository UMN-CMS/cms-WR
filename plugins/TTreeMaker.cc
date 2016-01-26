#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/Common/interface/ValueMap.h"


class TTreeMaker : public edm::EDAnalyzer {
public:
  explicit TTreeMaker(const edm::ParameterSet&);

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  static bool wayToSort(pat::Muon m1, pat::Muon m2);

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::InputTag muons_src;
  edm::EDGetToken electronsMiniAODToken_;
  const edm::InputTag jets_src;
  const edm::InputTag met_src;
  const edm::InputTag genparticles_src;
  const edm::InputTag genjets_src;
  const edm::InputTag primary_vertex_src;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleHEEPIdFullInfoMapToken_;
  bool verboseIdFlag_;

  const std::vector<std::string> bDiscriminators;

  const float Mlljj_cut;
  const float Mll_cut;
  const float isolation_dR;
  const float leading_lepton_pt_cut;
  const float lepton_eta_cut;
  const float subleading_lepton_pt_cut;
  const float leading_jet_pt_cut;
  const float jet_eta_cut;
  const float subleading_jet_pt_cut;
  const bool muon_mode;
  const bool electron_mode;
  const bool is_mc;

  struct tree_t {
    unsigned run;
    unsigned lumi;
    unsigned long long event;
    
    float leading_lepton_pt;
    float subleading_lepton_pt;
    float leading_jet_pt;
    float subleading_jet_pt;

    float leading_lepton_eta;
    float subleading_lepton_eta;
    float leading_jet_eta;
    float subleading_jet_eta;

    float leading_lepton_phi;
    float subleading_lepton_phi;
    float leading_jet_phi;
    float subleading_jet_phi;

    int leading_lepton_charge;
    int subleading_lepton_charge;

    float dilepton_mass;
    float Mlljj;
    float dR_leadLepton_leadJet;
    float dR_leadLepton_subleadJet;
    float dR_subleadLepton_leadJet;
    float dR_subleadLepton_subleadJet;
    
    float dR_leadLepton_subleadLepton;
    float dR_leadJet_subleadJet;
    
    int nleptons;
    int njets;
    int nvertices;

    float lepton_pt;
    float jet_pt;

    float weight;
    float angle3D;

    // Muon ID
    std::vector<bool> isGlobal;
    std::vector<int> numberOfValidMuonHits;
    std::vector<int> numberOfMatchedStations;
    std::vector<float> sigmapt;
    std::vector<float> dxy;
    std::vector<float> dz;
    std::vector<int> numberOfValidPixelHits;
    std::vector<int> trackerLayersWithMeasurement;

    // Btagging
    std::vector<float> leading_bTags;
    std::vector<float> subleading_bTags;

    // MET
    float met_pt;
    float met_phi;

    // Generator level quantities
    // Leading lepton corresponds to the WR daughter
    // Subleading lepton corresponds to the NR daughter

    float gen_leading_lepton_pt;
    float gen_subleading_lepton_pt;
    float gen_leading_jet_pt;
    float gen_subleading_jet_pt;

    float gen_leading_lepton_eta;
    float gen_subleading_lepton_eta;
    float gen_leading_jet_eta;
    float gen_subleading_jet_eta;

    float gen_leading_lepton_phi;
    float gen_subleading_lepton_phi;
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
      leading_lepton_pt = leading_lepton_eta = leading_lepton_phi = subleading_lepton_pt = subleading_lepton_eta = subleading_lepton_phi = -999;
      leading_lepton_charge = subleading_lepton_charge = -999;
      leading_jet_pt = leading_jet_eta = leading_jet_phi = subleading_jet_pt = subleading_jet_eta = subleading_jet_phi = -999;
      gen_leading_lepton_pt = gen_leading_lepton_eta = gen_leading_lepton_phi = gen_subleading_lepton_pt = gen_subleading_lepton_eta = gen_subleading_lepton_phi = -999;
      gen_leading_jet_pt = gen_leading_jet_eta = gen_leading_jet_phi = gen_subleading_jet_pt = gen_subleading_jet_eta = gen_subleading_jet_phi = -999;
      dilepton_mass = gen_dilepton_mass = Mlljj = gen_Mlljj = -999;
      dR_leadLepton_leadJet = dR_leadLepton_subleadJet = dR_subleadLepton_leadJet = dR_subleadLepton_subleadJet = 9;
      gen_dR_leadLepton_leadJet = gen_dR_leadLepton_subleadJet = gen_dR_subleadLepton_leadJet = gen_dR_subleadLepton_subleadJet = 9;

      dR_leadLepton_subleadLepton = dR_leadJet_subleadJet = 9;
      isGlobal.clear();
      numberOfValidMuonHits.clear();
      numberOfMatchedStations.clear();
      sigmapt.clear();
      dxy.clear();
      dz.clear();
      numberOfValidPixelHits.clear();
      trackerLayersWithMeasurement.clear();

      leading_bTags.clear();
      subleading_bTags.clear();

      met_pt = met_phi = -999;
    
      nleptons = njets = nvertices = -1;
      lepton_pt = jet_pt = -999;

      angle3D = -999;

      //pv_x.clear();
      
    }
  };

  TTree* tree;
  tree_t nt;
};

TTreeMaker::TTreeMaker(const edm::ParameterSet& cfg)
  : muons_src(cfg.getParameter<edm::InputTag>("muons_src")),
    jets_src(cfg.getParameter<edm::InputTag>("jets_src")),
    met_src(cfg.getParameter<edm::InputTag>("met_src")),
    genparticles_src(cfg.getParameter<edm::InputTag>("genparticles_src")),
    genjets_src(cfg.getParameter<edm::InputTag>("genjets_src")),
    primary_vertex_src(cfg.getParameter<edm::InputTag>("primary_vertex_src")),

    eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eleVetoIdMap"))),
    eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eleLooseIdMap"))),
    eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eleMediumIdMap"))),
    eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eleTightIdMap"))),
    eleHEEPIdMapToken_(consumes<edm::ValueMap<bool> >(cfg.getParameter<edm::InputTag>("eleHEEPIdMap"))),
    eleHEEPIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >
			       (cfg.getParameter<edm::InputTag>("eleHEEPIdFullInfoMap"))),
    verboseIdFlag_(cfg.getParameter<bool>("eleIdVerbose")),
    
    bDiscriminators(cfg.getParameter<std::vector<std::string> >("bDiscriminators")),

    Mlljj_cut(cfg.getParameter<double>("Mlljj_cut")),
    Mll_cut(cfg.getParameter<double>("Mll_cut")),
    isolation_dR(cfg.getParameter<double>("isolation_dR")),
    leading_lepton_pt_cut(cfg.getParameter<double>("leading_lepton_pt_cut")),
    lepton_eta_cut(cfg.getParameter<double>("lepton_eta_cut")),
    subleading_lepton_pt_cut(cfg.getParameter<double>("subleading_lepton_pt_cut")),
    leading_jet_pt_cut(cfg.getParameter<double>("leading_jet_pt_cut")),
    jet_eta_cut(cfg.getParameter<double>("jet_eta_cut")),
    subleading_jet_pt_cut(cfg.getParameter<double>("subleading_jet_pt_cut")),
    muon_mode(cfg.getParameter<bool>("muon_mode")),
    electron_mode(cfg.getParameter<bool>("electron_mode")),
    is_mc(cfg.getParameter<bool>("is_mc"))
{
  beamSpotToken_    = consumes<reco::BeamSpot>(cfg.getParameter <edm::InputTag>("beamSpot"));
  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >(cfg.getParameter<edm::InputTag>("conversionsMiniAOD"));
  electronsMiniAODToken_   = mayConsume<edm::View<reco::GsfElectron> >(cfg.getParameter<edm::InputTag>("electrons_src"));
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &nt.run);
  tree->Branch("lumi", &nt.lumi);
  tree->Branch("event", &nt.event);

  tree->Branch("leading_lepton_pt", &nt.leading_lepton_pt, "leading_lepton_pt/F");
  tree->Branch("leading_lepton_eta", &nt.leading_lepton_eta, "leading_lepton_eta/F");
  tree->Branch("leading_lepton_phi", &nt.leading_lepton_phi, "leading_lepton_phi/F");
  tree->Branch("subleading_lepton_pt", &nt.subleading_lepton_pt, "subleading_lepton_pt/F");
  tree->Branch("subleading_lepton_eta", &nt.subleading_lepton_eta, "subleading_lepton_eta/F");
  tree->Branch("subleading_lepton_phi", &nt.subleading_lepton_phi, "subleading_lepton_phi/F");
  tree->Branch("leading_jet_pt", &nt.leading_jet_pt, "leading_jet_pt/F");
  tree->Branch("leading_jet_eta", &nt.leading_jet_eta, "leading_jet_eta/F");
  tree->Branch("leading_jet_phi", &nt.leading_jet_phi, "leading_jet_phi/F");
  tree->Branch("subleading_jet_pt", &nt.subleading_jet_pt, "subleading_jet_pt/F");
  tree->Branch("subleading_jet_eta", &nt.subleading_jet_eta, "subleading_jet_eta/F");
  tree->Branch("subleading_jet_phi", &nt.subleading_jet_phi, "subleading_jet_phi/F");

  tree->Branch("leading_lepton_charge", &nt.leading_lepton_charge, "leading_lepton_charge/I");
  tree->Branch("subleading_lepton_charge", &nt.subleading_lepton_charge, "subleading_lepton_charge/I");

  tree->Branch("gen_leading_lepton_pt", &nt.gen_leading_lepton_pt, "gen_leading_lepton_pt/F");
  tree->Branch("gen_leading_lepton_eta", &nt.gen_leading_lepton_eta, "gen_leading_lepton_eta/F");
  tree->Branch("gen_leading_lepton_phi", &nt.gen_leading_lepton_phi, "gen_leading_lepton_phi/F");
  tree->Branch("gen_subleading_lepton_pt", &nt.gen_subleading_lepton_pt, "gen_subleading_lepton_pt/F");
  tree->Branch("gen_subleading_lepton_eta", &nt.gen_subleading_lepton_eta, "gen_subleading_lepton_eta/F");
  tree->Branch("gen_subleading_lepton_phi", &nt.gen_subleading_lepton_phi, "gen_subleading_lepton_phi/F");
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
  tree->Branch("dR_leadLepton_subleadLepton", &nt.dR_leadLepton_subleadLepton, "dR_leadLepton_subleadLepton/F");
  tree->Branch("dR_leadJet_subleadJet", &nt.dR_leadJet_subleadJet, "dR_leadJet_subleadJet/F");

  tree->Branch("gen_dR_leadLepton_leadJet", &nt.gen_dR_leadLepton_leadJet, "gen_dR_leadLepton_leadJet/F");
  tree->Branch("gen_dR_leadLepton_subleadJet", &nt.gen_dR_leadLepton_subleadJet, "gen_dR_leadLepton_subleadJet/F");
  tree->Branch("gen_dR_subleadLepton_leadJet", &nt.gen_dR_subleadLepton_leadJet, "gen_dR_subleadLepton_leadJet/F");
  tree->Branch("gen_dR_subleadLepton_subleadJet", &nt.gen_dR_subleadLepton_subleadJet, "gen_dR_subleadLepton_subleadJet/F");

  tree->Branch("isGlobal",&nt.isGlobal);
  tree->Branch("numberOfValidMuonHits",&nt.numberOfValidMuonHits);
  tree->Branch("numberOfMatchedStations",&nt.numberOfMatchedStations);
  tree->Branch("sigmapt",&nt.sigmapt);
  tree->Branch("dxy",&nt.dxy);
  tree->Branch("dz",&nt.dz);
  tree->Branch("numberOfValidPixelHits",&nt.numberOfValidPixelHits);
  tree->Branch("trackerLayersWithMeasurement",&nt.trackerLayersWithMeasurement);

  tree->Branch("leading_bTags",&nt.leading_bTags);
  tree->Branch("subleading_bTags",&nt.subleading_bTags);

  tree->Branch("met_pt",&nt.met_pt);
  tree->Branch("met_phi",&nt.met_phi);

  tree->Branch("nleptons",&nt.nleptons,"nleptons/i");
  tree->Branch("njets",&nt.njets,"njets/i");
  tree->Branch("nvertices",&nt.nvertices,"nvertices/i");
  tree->Branch("lepton_pt",&nt.lepton_pt,"lepton_pt/F");
  tree->Branch("jet_pt",&nt.jet_pt,"jet_pt/F");

  tree->Branch("weight",&nt.weight,"weight/F");
  tree->Branch("angle3D",&nt.angle3D,"angle3D/F");

}

bool TTreeMaker::wayToSort(pat::Muon m1, pat::Muon m2) { return m1.pt() > m2.pt(); }

void TTreeMaker::analyze(const edm::Event& event, const edm::EventSetup&) {
  nt.clear();
  nt.run = event.id().run();
  nt.lumi = event.luminosityBlock();
  nt.event = event.id().event();

  edm::Handle<pat::MuonRefVector> muons;
  event.getByLabel(muons_src, muons);
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  event.getByToken(electronsMiniAODToken_,electrons);
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jets_src, jets);
  edm::Handle<pat::METCollection> mets;
  event.getByLabel(met_src, mets);
  edm::Handle<reco::GenParticleCollection> gen_particles;
  event.getByLabel(genparticles_src, gen_particles);
  edm::Handle<reco::GenJetCollection> gen_jets;
  event.getByLabel(genjets_src, gen_jets);
  
  std::vector<reco::GenParticle> gmuons(2);
  std::vector<reco::GenParticle> geles(2);
  std::vector<reco::GenParticle> gen_quarks;
  std::vector<reco::GenJet> gjets(2);
    
  std::vector<pat::Muon> mus;
  std::vector<reco::GsfElectron> eles;
  std::vector<pat::Jet> js;
  edm::Handle<reco::VertexCollection> primary_vertex;
  event.getByLabel("offlineSlimmedPrimaryVertices", primary_vertex);
  edm::Handle<GenEventInfoProduct> evinfo;

  if(is_mc) 
    {
      event.getByLabel("generator", evinfo);
      nt.weight = evinfo->weight();
    }

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  event.getByToken(beamSpotToken_,theBeamSpot);  

  if (primary_vertex->empty()) return; // skip the event if no PV found
  
  // Find the first vertex in the collection that passes
  // good quality criteria
  reco::VertexCollection::const_iterator firstGoodVertex = primary_vertex->end();
  int firstGoodVertexIdx = 0;
  for (reco::VertexCollection::const_iterator vtx = primary_vertex->begin(); 
       vtx != primary_vertex->end(); ++vtx, ++firstGoodVertexIdx) {
    // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake
	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==primary_vertex->end() )
    return; // skip event if there are no good PVs

  // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  event.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  event.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
  event.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  event.getByToken(eleTightIdMapToken_,tight_id_decisions);
  event.getByToken(eleHEEPIdMapToken_ ,heep_id_decisions);
  // Full cut flow info for one of the working points:
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_cutflow_data;
  event.getByToken(eleHEEPIdFullInfoMapToken_,heep_id_cutflow_data);

  nt.met_pt = mets->at(0).pt();
  nt.met_phi = mets->at(0).phi();

  nt.nvertices = primary_vertex->size();
  
  for(auto j : *jets){  
    nt.jet_pt = j.pt();
    if(j.pt() > leading_jet_pt_cut && fabs(j.eta()) < jet_eta_cut){
      js.push_back(j);  
    }
  }
  nt.njets = js.size();

  // Muon Mode
  if(muon_mode && !electron_mode){
    for(auto muon : *muons){
		pat::Muon mu = *muon;
      if(mu.pt() > leading_lepton_pt_cut && fabs(mu.eta()) < lepton_eta_cut && mu.isHighPtMuon(primary_vertex->at(0)) && (mu.isolationR03().sumPt/mu.pt()) < 0.1 ){
	if(js.size()>1)
	  if(deltaR(js[0],mu)>isolation_dR && deltaR(js[1],mu)>isolation_dR)
	    mus.push_back(mu);	     
      }    
    }  

    sort(mus.begin(), mus.end(), wayToSort);    
    nt.nleptons = mus.size();

    if(mus.size() > 0){
      nt.isGlobal.push_back(mus[0].isGlobalMuon());
      nt.numberOfValidMuonHits.push_back(mus[0].numberOfValidHits());
      nt.numberOfMatchedStations.push_back(mus[0].numberOfMatchedStations());
      nt.sigmapt.push_back(mus[0].tunePMuonBestTrack()->ptError()/mus[0].pt());
      nt.dxy.push_back(mus[0].dB());
      nt.dz.push_back(mus[0].tunePMuonBestTrack()->dz(primary_vertex->at(0).position()));
      if(mus[0].innerTrack().isAvailable()){
	nt.numberOfValidPixelHits.push_back(mus[0].innerTrack()->hitPattern().numberOfValidPixelHits());
	nt.trackerLayersWithMeasurement.push_back(mus[0].innerTrack()->hitPattern().trackerLayersWithMeasurement());
      }
      nt.leading_lepton_pt = mus[0].pt();
      nt.leading_lepton_eta = mus[0].eta();
      nt.leading_lepton_phi = mus[0].phi();
      nt.leading_lepton_charge = mus[0].charge();
      if(js.size() > 0) nt.dR_leadLepton_leadJet = deltaR(mus[0],js[0]);
      if(js.size() > 1) nt.dR_leadLepton_subleadJet = deltaR(mus[0],js[1]);
      if(mus.size() > 1){
	nt.isGlobal.push_back(mus[1].isGlobalMuon());
	nt.numberOfValidMuonHits.push_back(mus[1].numberOfValidHits());
	nt.numberOfMatchedStations.push_back(mus[1].numberOfMatchedStations());
	nt.sigmapt.push_back(mus[1].tunePMuonBestTrack()->ptError()/mus[1].pt());
	nt.dxy.push_back(mus[1].dB());
	nt.dz.push_back(mus[1].tunePMuonBestTrack()->dz(primary_vertex->at(0).position()));
	if(mus[1].innerTrack().isAvailable()){
	  nt.numberOfValidPixelHits.push_back(mus[1].innerTrack()->hitPattern().numberOfValidPixelHits());
	  nt.trackerLayersWithMeasurement.push_back(mus[1].innerTrack()->hitPattern().trackerLayersWithMeasurement());
	}
	nt.angle3D = angle(mus[0].momentum().x(),mus[0].momentum().y(),mus[0].momentum().z(),mus[1].momentum().x(),mus[1].momentum().y(),mus[1].momentum().z());
	nt.subleading_lepton_pt = mus[1].pt();
	nt.subleading_lepton_eta = mus[1].eta();
	nt.subleading_lepton_phi = mus[1].phi();
	nt.subleading_lepton_charge = mus[1].charge();
	nt.dilepton_mass = (mus[0].p4() + mus[1].p4()).M();
	nt.dR_leadLepton_subleadLepton = deltaR(mus[0],mus[1]);
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
      for(auto b:bDiscriminators)
	nt.leading_bTags.push_back(js[0].bDiscriminator(b));
      if(js.size() > 1){
	nt.subleading_jet_pt = js[1].pt();
	nt.subleading_jet_eta = js[1].eta();
	nt.subleading_jet_phi = js[1].phi();
	nt.dR_leadJet_subleadJet = deltaR(js[0],js[1]);
	for(auto b:bDiscriminators)
	  nt.subleading_bTags.push_back(js[1].bDiscriminator(b));
      }
    }
  }

  // Muon Electron Mode
  // Muon = leading, electron = subleading
  if(muon_mode && electron_mode){
    for(auto muon : *muons){    
		pat::Muon mu = *muon;
      if(mu.pt() > leading_lepton_pt_cut && fabs(mu.eta()) < lepton_eta_cut && mu.isHighPtMuon(primary_vertex->at(0)) && (mu.isolationR03().sumPt/mu.pt()) < 0.1 ){
	if(js.size()>1)
	  if(deltaR(js[0],mu)>isolation_dR && deltaR(js[1],mu)>isolation_dR)
	    mus.push_back(mu);	     
      }    
    }
    for (size_t i = 0; i < electrons->size(); ++i){
      const auto ele = electrons->ptrAt(i);
      bool isPassHEEP = (*heep_id_decisions)[ele];
      if(ele->pt() > leading_lepton_pt_cut && fabs(ele->eta()) < lepton_eta_cut && isPassHEEP){
	if(js.size()>1)
	  if(deltaR(js[0],*ele)>isolation_dR && deltaR(js[1],*ele)>isolation_dR)
	    eles.push_back(*ele);	     
      }
    }

    sort(mus.begin(), mus.end(), wayToSort);
    nt.nleptons = mus.size() + eles.size();

    if(mus.size() > 0){
      nt.isGlobal.push_back(mus[0].isGlobalMuon());
      nt.numberOfValidMuonHits.push_back(mus[0].numberOfValidHits());
      nt.numberOfMatchedStations.push_back(mus[0].numberOfMatchedStations());
      nt.sigmapt.push_back(mus[0].tunePMuonBestTrack()->ptError()/mus[0].pt());
      nt.dxy.push_back(mus[0].dB());
      nt.dz.push_back(mus[0].tunePMuonBestTrack()->dz(primary_vertex->at(0).position()));
      if(mus[0].innerTrack().isAvailable()){
	nt.numberOfValidPixelHits.push_back(mus[0].innerTrack()->hitPattern().numberOfValidPixelHits());
	nt.trackerLayersWithMeasurement.push_back(mus[0].innerTrack()->hitPattern().trackerLayersWithMeasurement());
      }      nt.leading_lepton_pt = mus[0].pt();
      nt.leading_lepton_eta = mus[0].eta();
      nt.leading_lepton_phi = mus[0].phi();
      nt.leading_lepton_charge = mus[0].charge();
      if(js.size() > 0) nt.dR_leadLepton_leadJet = deltaR(mus[0],js[0]);
      if(js.size() > 1) nt.dR_leadLepton_subleadJet = deltaR(mus[0],js[1]);
    }

    if(eles.size() > 0){
      //std::cout<<"ele"<<std::endl;
      nt.subleading_lepton_pt = eles[0].pt();
      nt.subleading_lepton_eta = eles[0].eta();
      nt.subleading_lepton_phi = eles[0].phi();
      if(js.size() > 0) nt.dR_subleadLepton_leadJet = deltaR(eles[0],js[0]);
      if(js.size() > 1) nt.dR_subleadLepton_subleadJet = deltaR(eles[0],js[1]);
    }

    if(mus.size() > 0 && eles.size() > 0){
      //std::cout<<"mu and ele"<<std::endl;
      nt.dilepton_mass = (mus[0].p4() + eles[0].p4()).M();
      nt.dR_leadLepton_subleadLepton = deltaR(mus[0],eles[0]);
      if(js.size() > 1) nt.Mlljj = (mus[0].p4() + eles[0].p4() + js[0].p4() + js[1].p4()).M();
    }

    if(js.size() > 0){
      //std::cout<<"jet"<<std::endl;
      nt.leading_jet_pt = js[0].pt();
      nt.leading_jet_eta = js[0].eta();
      nt.leading_jet_phi = js[0].phi();
      for(auto b:bDiscriminators)
	nt.leading_bTags.push_back(js[0].bDiscriminator(b));
      if(js.size() > 1){
	nt.subleading_jet_pt = js[1].pt();
	nt.subleading_jet_eta = js[1].eta();
	nt.subleading_jet_phi = js[1].phi();
	nt.dR_leadJet_subleadJet = deltaR(js[0],js[1]);
	for(auto b:bDiscriminators)
	  nt.subleading_bTags.push_back(js[1].bDiscriminator(b));
      }
    }
  }
  // Electron Mode
  else if(!muon_mode && electron_mode){
    for (size_t i = 0; i < electrons->size(); ++i){
      const auto ele = electrons->ptrAt(i);
      bool isPassHEEP = (*heep_id_decisions)[ele];
      if(ele->pt() > leading_lepton_pt_cut && fabs(ele->eta()) < lepton_eta_cut && isPassHEEP){
	if(js.size()>1)
	  if(deltaR(js[0],*ele)>isolation_dR && deltaR(js[1],*ele)>isolation_dR)
	    eles.push_back(*ele);	     
      }
    }
    
    nt.nleptons = eles.size();

    if(eles.size() > 0){
      nt.leading_lepton_pt = eles[0].pt();
      nt.leading_lepton_eta = eles[0].eta();
      nt.leading_lepton_phi = eles[0].phi();
      if(js.size() > 0) nt.dR_leadLepton_leadJet = deltaR(eles[0],js[0]);
      if(js.size() > 1) nt.dR_leadLepton_subleadJet = deltaR(eles[0],js[1]);
      if(eles.size() > 1){
	nt.angle3D = angle(eles[0].momentum().x(),eles[0].momentum().y(),eles[0].momentum().z(),eles[1].momentum().x(),eles[1].momentum().y(),eles[1].momentum().z());
	nt.subleading_lepton_pt = eles[1].pt();
	nt.subleading_lepton_eta = eles[1].eta();
	nt.subleading_lepton_phi = eles[1].phi();
	nt.dilepton_mass = (eles[0].p4() + eles[1].p4()).M();
	nt.dR_leadLepton_subleadLepton = deltaR(eles[0],eles[1]);
	if(js.size() > 0) nt.dR_subleadLepton_leadJet = deltaR(eles[1],js[0]);
	if(js.size() > 1){
	  nt.dR_subleadLepton_subleadJet = deltaR(eles[1],js[1]);
	  nt.Mlljj = (eles[0].p4() + eles[1].p4() + js[0].p4() + js[1].p4()).M();
	}
      }
    }
    if(js.size() > 0){
      nt.leading_jet_pt = js[0].pt();
      nt.leading_jet_eta = js[0].eta();
      nt.leading_jet_phi = js[0].phi();
      for(auto b:bDiscriminators)
	nt.leading_bTags.push_back(js[0].bDiscriminator(b));
      if(js.size() > 1){
	nt.subleading_jet_pt = js[1].pt();
	nt.subleading_jet_eta = js[1].eta();
	nt.subleading_jet_phi = js[1].phi();
	nt.dR_leadJet_subleadJet = deltaR(js[0],js[1]);
	for(auto b:bDiscriminators)
	  nt.subleading_bTags.push_back(js[1].bDiscriminator(b));
      }
    }
  }
  
  tree->Fill();
}

DEFINE_FWK_MODULE(TTreeMaker);
