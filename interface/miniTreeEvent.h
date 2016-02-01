class miniTreeEvent{
 public:
  Int_t run;
  Int_t lumi;
  Long64_t event;

  // RunTime 
  // UInt_t run_time;
 
  std::vector<TLorentzVector> * electrons_p4;
  std::vector<TLorentzVector> * muons_p4;
  std::vector<TLorentzVector> * jets_p4;
    
  std::vector<Float_t> * jec_uncertainty;
  std::vector<Float_t> * electron_scale;
  std::vector<Float_t> * electron_smearing;
  std::vector<Int_t> * electron_charge;
  std::vector<Int_t> * muon_charge;
  std::vector<Float_t> * muon_IDSF_central;
  std::vector<Float_t> * muon_IsoSF_central;
  std::vector<Float_t> * muon_IDSF_error;
  std::vector<Float_t> * muon_IsoSF_error;

  Float_t nPU;
  Int_t nPV;
  Float_t weight;
  Float_t PU_reweight;
  
  miniTreeEvent();

  void clear();
  void SetBranches(TTree* tree);
  void SetBranchAddresses(TTree* tree, miniTreeEvent& event);

};

miniTreeEvent::miniTreeEvent():
electrons_p4(new std::vector<TLorentzVector>),
  muons_p4(new std::vector<TLorentzVector>),
  jets_p4(new std::vector<TLorentzVector>),
  jec_uncertainty(new std::vector<Float_t>),
  electron_scale(new std::vector<Float_t>),
  electron_smearing(new std::vector<Float_t>),
  electron_charge(new std::vector<Int_t>),
  muon_charge(new std::vector<Int_t>),
  muon_IDSF_central(new std::vector<Float_t>),
  muon_IsoSF_central(new std::vector<Float_t>),
  muon_IDSF_error(new std::vector<Float_t>),
  muon_IsoSF_error(new std::vector<Float_t>)
{
  clear();
};


void miniTreeEvent::clear() {
  run = lumi = event = 0;
    
  electrons_p4->clear();
  muons_p4->clear();
  jets_p4->clear();

  jec_uncertainty->clear();
  electron_scale->clear();
  electron_smearing->clear();
  electron_charge->clear();
  muon_charge->clear();
  muon_IDSF_central->clear();
  muon_IsoSF_central->clear();
  muon_IDSF_error->clear();
  muon_IsoSF_error->clear();

  nPU = -999.;
  nPV = 0.;
  weight = 0.0;
  PU_reweight = 0.0;

}

void miniTreeEvent::SetBranches(TTree* tree) {
  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("event", &event);

  tree->Branch("electrons_p4", electrons_p4,32000,-1);
  tree->Branch("muons_p4", muons_p4,32000,-1);
  tree->Branch("jets_p4", jets_p4,32000,-1);

  tree->Branch("jec_uncertainty",jec_uncertainty);
  tree->Branch("electron_scale",electron_scale);
  tree->Branch("electron_smearing",electron_smearing);
  tree->Branch("electron_charge",electron_charge);
  tree->Branch("muon_charge",muon_charge);
  tree->Branch("muon_IDSF_central",muon_IDSF_central);
  tree->Branch("muon_IsoSF_central",muon_IsoSF_central);
  tree->Branch("muon_IDSF_error",muon_IDSF_error);
  tree->Branch("muon_IsoSF_error",muon_IsoSF_error);

  tree->Branch("nPU", &nPU);
  tree->Branch("nPV", &nPV);
  tree->Branch("weight",&weight);
  tree->Branch("PU_reweight",&PU_reweight);

}

void miniTreeEvent::SetBranchAddresses(TTree* tree, miniTreeEvent& event) {

  delete event.electrons_p4;
  delete event.muons_p4;
  delete event.jets_p4;

  delete event.jec_uncertainty;
  delete event.electron_scale;
  delete event.electron_smearing;
  delete event.electron_charge;
  delete event.muon_charge;
  delete event.muon_IDSF_central;
  delete event.muon_IsoSF_central;
  delete event.muon_IDSF_error;
  delete event.muon_IsoSF_error;

  event.electrons_p4 = 0;
  event.muons_p4 = 0;
  event.jets_p4 = 0;

  event.jec_uncertainty = 0;
  event.electron_scale = 0;
  event.electron_smearing = 0;
  event.electron_charge = 0;
  event.muon_charge = 0;
  event.muon_IDSF_central = 0;
  event.muon_IsoSF_central = 0;
  event.muon_IDSF_error = 0;
  event.muon_IsoSF_error = 0;

  tree->SetBranchAddress("run",&event.run);
  tree->SetBranchAddress("lumi", &event.lumi);
  tree->SetBranchAddress("event", &event.event);

  tree->SetBranchAddress("electrons_p4", &event.electrons_p4);
  tree->SetBranchAddress("muons_p4", &event.muons_p4);
  tree->SetBranchAddress("jets_p4", &event.jets_p4);

  tree->SetBranchAddress("jec_uncertainty",&event.jec_uncertainty);
  tree->SetBranchAddress("electron_scale",&event.electron_scale);
  tree->SetBranchAddress("electron_smearing",&event.electron_smearing);
  tree->SetBranchAddress("electron_charge",&event.electron_charge);
  tree->SetBranchAddress("muon_charge",&event.muon_charge);
  tree->SetBranchAddress("muon_IDSF_central",&event.muon_IDSF_central);
  tree->SetBranchAddress("muon_IsoSF_central",&event.muon_IsoSF_central);
  tree->SetBranchAddress("muon_IDSF_error",&event.muon_IDSF_error);
  tree->SetBranchAddress("muon_IsoSF_error",&event.muon_IsoSF_error);

  tree->SetBranchAddress("nPU", &event.nPU);
  tree->SetBranchAddress("nPV", &event.nPV);
  tree->SetBranchAddress("weight",&event.weight);
  tree->SetBranchAddress("PU_reweight",&event.PU_reweight);

}
