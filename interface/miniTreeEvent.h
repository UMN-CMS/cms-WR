class miniTreeEvent{
 public:
  Int_t run;
  Int_t lumi;
  Long64_t event;

  // RunTime 
  // UInt_t run_time;
 
  std::vector<TLorentzVector> electrons_p4;
  std::vector<TLorentzVector> muons_p4;
  std::vector<TLorentzVector> jets_p4;
    
  std::vector<float> jet_uncertainty;
  std::vector<float> electron_scale;
  std::vector<float> electron_smearing;

  Int_t nPU;
  Int_t nPV;
  Float_t weight;

  miniTreeEvent(void){ clear(); };
  void clear();
  void SetBranches(TTree* tree);
  void SetBranchAddresses(TTree* tree, miniTreeEvent& event);

  std::vector<TLorentzVector> * electrons_p4_ptr;
  std::vector<TLorentzVector> * muons_p4_ptr;
  std::vector<TLorentzVector> * jets_p4_ptr;

  std::vector<float> * jet_uncertainty_ptr;
  std::vector<float> * electron_scale_ptr;
  std::vector<float> * electron_smearing_ptr;

};


void miniTreeEvent::clear() {
  run = lumi = event = 0;
    
  electrons_p4.clear();
  muons_p4.clear();
  jets_p4.clear();

  jet_uncertainty.clear();
  electron_scale.clear();
  electron_smearing.clear();

  nPU = -999.;
  nPV = 0.;
  weight = 0.0;
      
}

void miniTreeEvent::SetBranches(TTree* tree) {
  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("event", &event);

  tree->Branch("electrons_p4", &electrons_p4,32000,-1);
  tree->Branch("muons_p4", &muons_p4,32000,-1);
  tree->Branch("jets_p4", &jets_p4,32000,-1);

  tree->Branch("jet_uncertainty",&jet_uncertainty);
  tree->Branch("electron_scale",&electron_scale);
  tree->Branch("electron_smearing",&electron_smearing);

  tree->Branch("nPU", &nPU);
  tree->Branch("nPV", &nPV);
  tree->Branch("weight",&weight);

}

void miniTreeEvent::SetBranchAddresses(TTree* tree, miniTreeEvent& event) {
  
  event.electrons_p4_ptr = 0;
  event.muons_p4_ptr = 0;
  event.jets_p4_ptr = 0;

  event.jet_uncertainty_ptr = 0;
  event.electron_scale_ptr = 0;
  event.electron_smearing_ptr = 0;

  tree->SetBranchAddress("run",&event.run);
  tree->SetBranchAddress("lumi", &event.lumi);
  tree->SetBranchAddress("event", &event.event);

  tree->SetBranchAddress("electrons_p4", &event.electrons_p4_ptr);
  tree->SetBranchAddress("muons_p4", &event.muons_p4_ptr);
  tree->SetBranchAddress("jets_p4", &event.jets_p4_ptr);

  tree->SetBranchAddress("jet_uncertainty",&event.jet_uncertainty_ptr);
  tree->SetBranchAddress("electron_scale",&event.electron_scale_ptr);
  tree->SetBranchAddress("electron_smearing",&event.electron_smearing_ptr);

  tree->SetBranchAddress("nPU", &event.nPU);
  tree->SetBranchAddress("nPV", &event.nPV);
  tree->SetBranchAddress("weight",&event.weight);

}
