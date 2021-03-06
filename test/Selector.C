#include "../interface/Selector.h"
#include "DataFormats/Math/interface/deltaR.h"
//#define DEBUG

float dR_TLV(TLorentzVector t1,TLorentzVector t2) {return deltaR(t1.Eta(),t1.Phi(),t2.Eta(),t2.Phi()); }
void goodJets(std::vector<TLorentzVector> *evJets, std::vector<TLorentzVector> *selJets){
  for(auto j:*evJets){
    if(j.Pt() > 40 && fabs(j.Eta()) < 2.4)
      selJets->push_back(j);
  }
}
void goodEles(std::vector<TLorentzVector> *evEles, std::vector<TLorentzVector> *selEles){
  for(auto e:*evEles){
    if(e.Pt() > 40 && fabs(e.Eta()) < 2.4)
      selEles->push_back(e);
  }
}
void goodMuons(std::vector<TLorentzVector> *evMuons, std::vector<TLorentzVector> *selMuons){
  for(auto m:*evMuons){
    if(m.Pt() > 40 && fabs(m.Eta()) < 2.4)
      selMuons->push_back(m);
  }
}

int SelectEvent(Selector *selEvent, Int_t tag){

  std::vector<TLorentzVector> gJets;
  std::vector<TLorentzVector> gEles;
  std::vector<TLorentzVector> gMuons;

  // Basic Kinematic cuts
  goodJets(selEvent->jets_p4, &gJets);
  goodEles(selEvent->electrons_p4, &gEles);
  goodMuons(selEvent->muons_p4, &gMuons);
  
  selEvent->pass_selection = true;

  // Assert at least 2 good jets
  if(gJets.size() < 2){
    selEvent->pass_selection = false;
    return 0;
  }

#ifdef DEBUG
  if(gJets[0].Pt() != (selEvent->jets_p4)->at(0).Pt()){
    cout<<"UNSORTED JETS"<<endl;
    cout<<gJets[0].Pt()<< " " <<(selEvent->jets_p4)->at(0).Pt()<<endl;
    cout<<gJets[1].Pt()<< " " <<(selEvent->jets_p4)->at(1).Pt()<<endl;
    cout<<gJets[0].Eta()<< " " <<(selEvent->jets_p4)->at(0).Eta()<<endl;
    cout<<gJets[1].Eta()<< " " <<(selEvent->jets_p4)->at(1).Eta()<<endl;
  }
#endif

  if(tag == 0) // EEJJ Channel
    {
      // Assert at least 2 good leptons
      if(gEles.size() < 2){
	selEvent->pass_selection = false;
	return 0;
      }

      if(gEles[0].Pt() != (selEvent->electrons_p4)->at(0).Pt())
	cout<<"UNSORTED ELECTRONS"<<endl;

      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      selEvent->WR_mass = (gEles[0] + gEles[1] + gJets[0] + gJets[1]).M();
      selEvent->dilepton_mass = (gEles[0] + gEles[1]).M();
      // Apply dR cuts and leading lepton cut
      if(gEles[0].Pt() < 60 || dR_TLV(gEles[0],gJets[0]) < 0.4 || dR_TLV(gEles[0],gJets[1]) < 0.4 || dR_TLV(gEles[1],gJets[0]) < 0.4 || dR_TLV(gEles[1],gJets[1]) < 0.4)
	selEvent->pass_selection = false;
    }

  if(tag == 1) // MuMuJJ Channel
    {
      // Assert at least 2 good leptons
      if(gMuons.size() < 2){
	selEvent->pass_selection = false;
	return 0;
      }
#ifdef DEBUG
      if(gMuons[0].Pt() != (selEvent->muons_p4)->at(0).Pt()){
	cout<<"UNSORTED MUONS"<<endl;
	cout<<gMuons[0].Pt()<< " " <<(selEvent->muons_p4)->at(0).Pt()<<endl;
	cout<<gMuons[1].Pt()<< " " <<(selEvent->muons_p4)->at(1).Pt()<<endl;
      }
#endif
      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      selEvent->WR_mass = (gMuons[0] + gMuons[1] + gJets[0] + gJets[1]).M();
      selEvent->dilepton_mass = (gMuons[0] + gMuons[1]).M();
      // Apply dR cuts and leading lepton cut
      if(gMuons[0].Pt() < 60 || dR_TLV(gMuons[0],gJets[0]) < 0.4 || dR_TLV(gMuons[0],gJets[1]) < 0.4 || dR_TLV(gMuons[1],gJets[0]) < 0.4 || dR_TLV(gMuons[1],gJets[1]) < 0.4)
	selEvent->pass_selection = false;
    }
  if(tag == 2) // EMuJJ Channel
    {
      // Assert at least 2 good leptons
      if(gEles.size() < 1 || gMuons.size() < 1){
	selEvent->pass_selection = false;
	return 0;
      }
      // Build the WR mass and dilepton mass with the 2 highest pT jets and 2 highest pT leptons
      selEvent->WR_mass = (gEles[0] + gMuons[0] + gJets[0] + gJets[1]).M();
      selEvent->dilepton_mass = (gEles[0] + gMuons[0]).M();
      // Apply dR cuts and leading lepton cut
      if((gEles[0].Pt() < 60 && gMuons[0].Pt() < 60) || dR_TLV(gEles[0],gJets[0]) < 0.4 || dR_TLV(gEles[0],gJets[1]) < 0.4 || dR_TLV(gMuons[0],gJets[0]) < 0.4 || dR_TLV(gMuons[0],gJets[1]) < 0.4)
	selEvent->pass_selection = false;
    }

  return 0;

}
