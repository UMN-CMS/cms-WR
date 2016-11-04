#include<fstream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TFile.h"

#include "ExoAnalysis/cmsWR/interface/JetResolution.h"

#include <iostream>


void JetResolution(miniTreeEvent *myEvent, TRandom3 rand, bool isData)
{
	int Iterator = 0;
	if(isData)
	  return;
	for(auto jets : *(myEvent->jets_p4)) {
	  float smeared_pt = 0.0;
	  float reco_pt = (*(myEvent->jets_p4))[Iterator].Pt();
	  float gen_pt = (*(myEvent->genJetPt))[Iterator];
	  float JER = (*(myEvent->jetResolution))[Iterator];
	  float sf = (*(myEvent->JER_sf))[Iterator];
	  
	  if((*(myEvent->genJetMatch))[Iterator])
	    smeared_pt = TMath::Max(float(0.),gen_pt+(sf*(reco_pt-gen_pt)));
	  else
	    smeared_pt = rand.Gaus(reco_pt,TMath::Sqrt(sf*sf-1)*JER);

	  (*(myEvent->jets_p4))[Iterator].SetPtEtaPhiM(smeared_pt, (*(myEvent->jets_p4))[Iterator].Eta(), (*(myEvent->jets_p4))[Iterator].Phi(), (*(myEvent->jets_p4))[Iterator].M());

	  Iterator++;
	}


}

