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
#include "ExoAnalysis/cmsWR/interface/RoccoR.h"
#include "ExoAnalysis/cmsWR/interface/rochcor2015.h"
#include "ExoAnalysis/cmsWR/interface/ToyThrower.h"

#ifdef DEBUG
#include <iostream>
#endif

void ToyThrower(miniTreeEvent *myEvent,  float rand_smear[], float rand_up_down[], int random_seed, std::vector<std::string> list, bool isData)
{

	float qter = 1.0;
	int Iterator = 0;
	int Iterator_Up_Down = 0;
	int Flag_Smear_Muon_Scale = 0, Flag_Smear_Muon_ID_Iso = 0, Flag_Smear_Muon_Resolution = 0, Flag_Smear_Jet_Scale = 0;
	double Smear_ID = rand_up_down[Iterator_Up_Down++];
	double Smear_ISO = rand_up_down[Iterator_Up_Down++];
	double Smear_Res = rand_up_down[Iterator_Up_Down++];
	double Smear_Jet_Scale = rand_up_down[Iterator_Up_Down++];

	for(unsigned int iii = 0; iii < list.size(); iii++) {
		if(list[iii] == "Smear_Muon_Scale")     Flag_Smear_Muon_Scale = 1;
		else if(list[iii] == "Smear_Muon_ID_Iso")    Flag_Smear_Muon_ID_Iso = 1;
		else if(list[iii] == "Smear_Jet_Scale")      Flag_Smear_Jet_Scale = 1;
		else if(list[iii] == "Smear_Muon_Resolution")      Flag_Smear_Muon_Resolution = 1;
	}

	rochcor2015 *rmcor = new rochcor2015(random_seed);

	for(auto muons : * (myEvent->muons_p4)) {
#ifdef DEBUG
		std::cout << std::endl << " Muon number= " << Iterator << " Muon Pt Before = " << (*(myEvent->muons_p4))[Iterator].Pt() << " Muon Eta Before = " << (*(myEvent->muons_p4))[Iterator].Eta() << std::endl;
#endif

		if(Flag_Smear_Muon_ID_Iso && !isData ) {
			if(Smear_ID >= 0.) (*(myEvent->muon_IDSF_central))[Iterator] += Smear_ID * (*(myEvent->muon_IDSF_error))[Iterator];
			else            (*(myEvent->muon_IDSF_central))[Iterator] -= Smear_ID * (*(myEvent->muon_IDSF_error))[Iterator];

			if(Smear_ISO >= 0.) (*(myEvent->muon_IsoSF_central))[Iterator] += Smear_ISO * (*(myEvent->muon_IsoSF_error))[Iterator];
			else            (*(myEvent->muon_IsoSF_central))[Iterator] -= Smear_ISO * (*(myEvent->muon_IsoSF_error))[Iterator];
		}

		if(Flag_Smear_Muon_Scale) {
			if(isData)
				rmcor->momcor_data((*(myEvent->muons_p4))[Iterator], (*(myEvent->muon_charge))[Iterator], 0, qter);
			else rmcor->momcor_mc((*(myEvent->muons_p4))[Iterator], (*(myEvent->muon_charge))[Iterator], 0, qter);

		}

		if(Flag_Smear_Muon_Resolution) {
		        if(Smear_Res >= 0.) *(myEvent->muons_p4).Pt() += Smear_Res * 0.1;
			else *(myEvent->muons_p4).Pt() -= Smear_Res * 0.1;

		}

#ifdef DEBUG
		std::cout << std::endl << " Muon number= " << Iterator << " Muon Pt After = " << (*(myEvent->muons_p4))[Iterator].Pt() << " Muon Eta After = " << (*(myEvent->muons_p4))[Iterator].Eta() << std::endl;
#endif
		Iterator++;
	}

	Iterator = 0;
	for(auto jets : * (myEvent->jets_p4)) {
#ifdef DEBUG
		std::cout << std::endl << " Jet number= " << Iterator << " " << (*(myEvent->jec_uncertainty))[Iterator] << " " << Smear << " Jet Pt Before = " << (*(myEvent->jets_p4))[Iterator].Pt() << " Jet Eta Before = " << (*(myEvent->jets_p4))[Iterator].Eta() << std::endl;
#endif

		if(Flag_Smear_Jet_Scale)
			(*(myEvent->jets_p4))[Iterator] = (1 + (Smear_Jet_Scale) * (*(myEvent->jec_uncertainty))[Iterator]) * (*(myEvent->jets_p4))[Iterator];
#ifdef DEBUG
		std::cout << std::endl << " Jet number= " << Iterator << " Jet Pt Before = " << (*(myEvent->jets_p4))[Iterator].Pt() << " Jet Eta After = " << (*(myEvent->jets_p4))[Iterator].Eta() << std::endl;
#endif

		Iterator++;
	}

	delete rmcor;
}

