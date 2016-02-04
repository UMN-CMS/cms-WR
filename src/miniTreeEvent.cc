// #include "ExoAnalysis/cmsWR/interface/miniTreeEvent.h"
#include "../interface/miniTreeEvent.h"

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
  muon_IsoSF_error(new std::vector<Float_t>),
  _owningMembers(true)
{
	
  clear();
};


miniTreeEvent::miniTreeEvent(const miniTreeEvent& otherEvent):
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
  muon_IsoSF_error(new std::vector<Float_t>),
  _owningMembers(true)
{
  clear();
  *electrons_p4 = *(otherEvent.electrons_p4);
  *muons_p4 = *(otherEvent.muons_p4);
  *jets_p4 = *(otherEvent.jets_p4);
  *jec_uncertainty = *(otherEvent.jec_uncertainty);
  *electron_scale = *(otherEvent.electron_scale);
  *electron_smearing = *(otherEvent.electron_smearing);
  *electron_charge = *(otherEvent.electron_charge);
  *muon_charge = *(otherEvent.muon_charge);
  *muon_IDSF_central = *(otherEvent.muon_IDSF_central);
  *muon_IsoSF_central = *(otherEvent.muon_IsoSF_central);
  *muon_IDSF_error = *(otherEvent.muon_IDSF_error);
  *muon_IsoSF_error = *(otherEvent.muon_IsoSF_error);
  _owningMembers = false;

  run = otherEvent.run;
  lumi = otherEvent.lumi;
  event = otherEvent.event;

  nPU = otherEvent.nPU;
  nPV = otherEvent.nPV;
  weight = otherEvent.weight;
  PU_reweight = otherEvent.PU_reweight;
  

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

void miniTreeEvent::SetBranchAddresses(TChain* tree) {

  delete electrons_p4;
  delete muons_p4;
  delete jets_p4;

  delete jec_uncertainty;
  delete electron_scale;
  delete electron_smearing;
  delete electron_charge;
  delete muon_charge;
  delete muon_IDSF_central;
  delete muon_IsoSF_central;
  delete muon_IDSF_error;
  delete muon_IsoSF_error;

  _owningMembers=false;

  electrons_p4 = 0;
  muons_p4 = 0;
  jets_p4 = 0;

  jec_uncertainty = 0;
  electron_scale = 0;
  electron_smearing = 0;
  electron_charge = 0;
  muon_charge = 0;
  muon_IDSF_central = 0;
  muon_IsoSF_central = 0;
  muon_IDSF_error = 0;
  muon_IsoSF_error = 0;

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("lumi", &lumi);
  tree->SetBranchAddress("event", &event);

  // tree->SetBranchAddress("electrons_p4", &electrons_p4);
  // tree->SetBranchAddress("muons_p4", &muons_p4);
  // tree->SetBranchAddress("jets_p4", &jets_p4);

  // tree->SetBranchAddress("jec_uncertainty",&jec_uncertainty);
  // tree->SetBranchAddress("electron_scale",&electron_scale);
  // tree->SetBranchAddress("electron_smearing",&electron_smearing);
  // tree->SetBranchAddress("electron_charge",&electron_charge);
  // tree->SetBranchAddress("muon_charge",&muon_charge);
  // tree->SetBranchAddress("muon_IDSF_central",&muon_IDSF_central);
  // tree->SetBranchAddress("muon_IsoSF_central",&muon_IsoSF_central);
  // tree->SetBranchAddress("muon_IDSF_error",&muon_IDSF_error);
  // tree->SetBranchAddress("muon_IsoSF_error",&muon_IsoSF_error);

  // tree->SetBranchAddress("nPU", &nPU);
  // tree->SetBranchAddress("nPV", &nPV);
  // tree->SetBranchAddress("weight",&weight);
  // tree->SetBranchAddress("PU_reweight",&PU_reweight);

}
