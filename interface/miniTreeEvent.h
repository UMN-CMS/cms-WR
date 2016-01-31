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
    
  std::vector<float> * jec_uncertainty;
  std::vector<float> * electron_scale;
  std::vector<float> * electron_smearing;

  Int_t nPU;
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
  jec_uncertainty(new std::vector<float>),
  electron_scale(new std::vector<float>),
  electron_smearing(new std::vector<float>)
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

  tree->Branch("nPU", &nPU);
  tree->Branch("nPV", &nPV);
  tree->Branch("weight",&weight);
  tree->Branch("PU_reweight",&PU_reweight);

}

void miniTreeEvent::SetBranchAddresses(TTree* tree, miniTreeEvent& event) {
  
  event.electrons_p4 = 0;
  event.muons_p4 = 0;
  event.jets_p4 = 0;

  event.jec_uncertainty = 0;
  event.electron_scale = 0;
  event.electron_smearing = 0;

  tree->SetBranchAddress("run",&event.run);
  tree->SetBranchAddress("lumi", &event.lumi);
  tree->SetBranchAddress("event", &event.event);

  tree->SetBranchAddress("electrons_p4", &event.electrons_p4);
  tree->SetBranchAddress("muons_p4", &event.muons_p4);
  tree->SetBranchAddress("jets_p4", &event.jets_p4);

  tree->SetBranchAddress("jec_uncertainty",&event.jec_uncertainty);
  tree->SetBranchAddress("electron_scale",&event.electron_scale);
  tree->SetBranchAddress("electron_smearing",&event.electron_smearing);

  tree->SetBranchAddress("nPU", &event.nPU);
  tree->SetBranchAddress("nPV", &event.nPV);
  tree->SetBranchAddress("weight",&event.weight);
  tree->SetBranchAddress("PU_reweight",&event.PU_reweight);

}
