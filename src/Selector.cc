#include <iostream>
#include <algorithm>
#include <vector>
#include "ExoAnalysis/cmsWR/interface/Selector.h"
#include "DataFormats/Math/interface/deltaR.h"

float dR_TLV(TLorentzVector t1,TLorentzVector t2) {return deltaR(t1.Eta(),t1.Phi(),t2.Eta(),t2.Phi()); };
void goodJets(myJetCollection *evJets, myJetCollection *selJets){
  for(auto j:*evJets){
    if(j.p4.Pt() > 40 && fabs(j.p4.Eta()) < 2.4)
      selJets->push_back(j);
  }
};
void goodEles(myElectronCollection *evEles, myElectronCollection *selEles){
  for(auto e:*evEles){
    if(e.p4.Pt() > 40 && fabs(e.p4.Eta()) < 2.4 && fabs(e.p4.Eta()) < 1.4222 && fabs(e.p4.Eta()) > 1.566)
      selEles->push_back(e);
  }
};
void goodMuons(myMuonCollection *evMuons, myMuonCollection *selMuons){
  for(auto m:*evMuons){
    if(m.p4.Pt() > 40 && fabs(m.p4.Eta()) < 2.4)
      selMuons->push_back(m);
  }
};

Selector::Selector(const miniTreeEvent& myEvent) : 
	WR_mass(-1),
	dilepton_mass(-1)
{
  int nele = myEvent.electrons_p4->size();
  for(int i = 0; i < nele; i++){
    myElectron ele;
    ele.p4 = myEvent.electrons_p4->at(i);
    ele.scale = myEvent.electron_scale->at(i);
    ele.smearing = myEvent.electron_smearing->at(i);
    ele.charge = myEvent.electron_charge->at(i);
    electrons.push_back(ele);
  }
  int nmu = myEvent.muons_p4->size();
  for(int i = 0; i < nmu; i++){
    myMuon mu;
    mu.p4 = myEvent.muons_p4->at(i);
    mu.IDSF = myEvent.muon_IDSF_central->at(i);
    mu.IsoSF = myEvent.muon_IsoSF_central->at(i);
    mu.IDSF_error = myEvent.muon_IDSF_error->at(i);
    mu.IsoSF_error = myEvent.muon_IsoSF_error->at(i);
    mu.charge = myEvent.muon_charge->at(i);
    muons.push_back(mu);
  }
  int njet = myEvent.jets_p4->size();
  for(int i = 0; i < njet; i++){
    myJet jet;
    jet.p4 = myEvent.jets_p4->at(i);
    jet.jec_uncertainty = myEvent.jec_uncertainty->at(i);
    jets.push_back(jet);
  }

  Clear();
};

Selector::Selector(){
  Clear();
};

void Selector::Clear() {
  WR_mass = dilepton_mass = weight = 0.0;
};


bool Selector::isPassing(Int_t tag){

  myJetCollection gJets;
  myElectronCollection gEles;
  myMuonCollection gMuons;

  // Basic Kinematic cuts
  goodJets(&jets, &gJets);
  goodEles(&electrons, &gEles);
  goodMuons(&muons, &gMuons);

  std::sort(gJets.begin(), gJets.end(),
	    [](myJet const &a, myJet const &b) {
	      return a.p4.Pt() < b.p4.Pt(); 
	    });
  std::sort(gEles.begin(), gEles.end(),
	    [](myElectron const &a, myElectron const &b) {
	      return a.p4.Pt() < b.p4.Pt(); 
	    });
  std::sort(gMuons.begin(), gMuons.end(),
	    [](myMuon const &a, myMuon const &b) {
	      return a.p4.Pt() < b.p4.Pt(); 
	    });

  jets = gJets;
  electrons = gEles;
  muons = gMuons;

  // Assert at least 2 good jets
  if(gJets.size() < 2){
    return false;
  }
  
#ifdef DEBUG
  if(gJets[0].p4.Pt() != (jets_p4)->at(0).p4.Pt()){
    cout<<"UNSORTED JETS"<<endl;
    cout<<gJets[0].Pt()<< " " <<(jets->at(0)).p4.Pt()<<endl;
    cout<<gJets[1].Pt()<< " " <<(jets_p4)->at(1).Pt()<<endl;
    cout<<gJets[0].Eta()<< " " <<(jets_p4)->at(0).Eta()<<endl;
    cout<<gJets[1].Eta()<< " " <<(jets_p4)->at(1).Eta()<<endl;
  }
#endif
  
  if(tag == 0) // EEJJ Channel
    {
      // Assert at least 2 good leptons
      if(gEles.size() < 2){	
	return false;
      }
      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      WR_mass = (gEles[0].p4 + gEles[1].p4 + gJets[0].p4 + gJets[1].p4).M();
      dilepton_mass = (gEles[0].p4 + gEles[1].p4).M();
      // Apply dR cuts and leading lepton cut
      if(gEles[0].p4.Pt() < 60)
	return false;
      if(dR_TLV(gEles[0].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gEles[0].p4,gJets[1].p4) < 0.4)
	return false;
      if(dR_TLV(gEles[1].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gEles[1].p4,gJets[1].p4) < 0.4)
	return false;
      return true;      
    }

  if(tag == 1) // MuMuJJ Channel
    {
      // Assert at least 2 good leptons
      if(gMuons.size() < 2){
	return false;
      }
#ifdef DEBUG
      if(gMuons[0].Pt() != (muons_p4)->at(0).Pt()){
	cout<<"UNSORTED MUONS"<<endl;
	cout<<gMuons[0].p4.Pt()<< " " <<(muons.p4)->at(0).Pt()<<endl;
	cout<<gMuons[1].p4.Pt()<< " " <<(muons.p4)->at(1).Pt()<<endl;
      }
#endif
      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      WR_mass = (gMuons[0].p4 + gMuons[1].p4 + gJets[0].p4 + gJets[1].p4).M();
      dilepton_mass = (gMuons[0].p4 + gMuons[1].p4).M();
      // Apply dR cuts and leading lepton cut
      if(gMuons[0].p4.Pt() < 60)
	return false;
      if(dR_TLV(gMuons[0].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gMuons[0].p4,gJets[1].p4) < 0.4)
	return false;
      if(dR_TLV(gMuons[1].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gMuons[1].p4,gJets[1].p4) < 0.4)
	return false;
      return true;
    }
  if(tag == 2) // EMuJJ Channel
    {
      // Assert at least 2 good leptons
      if(gEles.size() < 1 || gMuons.size() < 1){
	return false;
      }
      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      WR_mass = (gEles[0].p4 + gMuons[0].p4 + gJets[0].p4 + gJets[1].p4).M();
      dilepton_mass = (gEles[0].p4 + gMuons[0].p4).M();
      // Apply dR cuts and leading lepton cut
      if(gEles[0].p4.Pt() < 60 && gMuons[0].p4.Pt() < 60)
	return false;
      if(dR_TLV(gEles[0].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gEles[0].p4,gJets[1].p4) < 0.4)
	return false;
      if(dR_TLV(gMuons[0].p4,gJets[0].p4) < 0.4)
	return false;
      if(dR_TLV(gMuons[0].p4,gJets[1].p4) < 0.4)
	return false;
      return true;
    }
  return false;  
}

void Selector::SetBranches(TTree* tree) {

  Float_t lead_lepton_pt;
  Float_t sublead_lepton_pt;
  Float_t lead_lepton_eta;
  Float_t sublead_lepton_eta;
  Float_t lead_lepton_phi;
  Float_t sublead_lepton_phi;
  Float_t lead_jet_pt;
  Float_t sublead_jet_pt;
  Float_t lead_jet_eta;
  Float_t sublead_jet_eta;
  Float_t lead_jet_phi;
  Float_t sublead_jet_phi;

  Float_t weight;    

  tree->Branch("lead_lepton_pt", &lead_lepton_pt);
  tree->Branch("sublead_lepton_pt", &sublead_lepton_pt);
  tree->Branch("lead_lepton_eta", &lead_lepton_eta);
  tree->Branch("sublead_lepton_eta", &sublead_lepton_eta);
  tree->Branch("lead_lepton_phi", &lead_lepton_phi);
  tree->Branch("sublead_lepton_phi", &sublead_lepton_phi);
  tree->Branch("lead_jet_pt", &lead_jet_pt);
  tree->Branch("sublead_jet_pt", &sublead_jet_pt);
  tree->Branch("lead_jet_eta", &lead_jet_eta);
  tree->Branch("sublead_jet_eta", &sublead_jet_eta);
  tree->Branch("lead_jet_phi", &lead_jet_phi);
  tree->Branch("sublead_jet_phi", &sublead_jet_phi);

  tree->Branch("weight", &weight);
  tree->Branch("WR_mass", &WR_mass);
  tree->Branch("dilepton_mass", &dilepton_mass);
  
}

void Selector::SetBranchAddresses(TTree* tree) {
  
  Float_t lead_lepton_pt;
  Float_t sublead_lepton_pt;
  Float_t lead_lepton_eta;
  Float_t sublead_lepton_eta;
  Float_t lead_lepton_phi;
  Float_t sublead_lepton_phi;
  Float_t lead_jet_pt;
  Float_t sublead_jet_pt;
  Float_t lead_jet_eta;
  Float_t sublead_jet_eta;
  Float_t lead_jet_phi;
  Float_t sublead_jet_phi;

  Float_t weight;  

  tree->SetBranchAddress("leading_lepton_pt",&lead_lepton_pt);
  /*
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
  */
}

