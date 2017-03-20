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
#include "ExoAnalysis/cmsWR/interface/muresolution_run2.h"
#include "ExoAnalysis/cmsWR/interface/rochcor2015.h"
#include "ExoAnalysis/cmsWR/interface/ToyThrower.h"
#include "math.h"

#ifdef DEBUG
#include <iostream>
#endif

void ToyThrower(miniTreeEvent *myEvent,  float rand_smear[], float rand_up_down[], int random_seed, std::vector<std::string> list, bool isData, bool firstLoop)
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
		else if(list[iii] == "Smear_Muon_Resolution")      Flag_Smear_Muon_Resolution = 1;
		else if(list[iii] == "Smear_Jet_Scale")      Flag_Smear_Jet_Scale = 1;
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

		if(Flag_Smear_Muon_Resolution && !isData) {
			float smeared_pt = (*(myEvent->muons_p4))[Iterator].Pt();
			if(Smear_Res >= 0.)
				smeared_pt += Smear_Res * 0.04 * smeared_pt;
			else smeared_pt -= Smear_Res * 0.04 * smeared_pt;

			(*(myEvent->muons_p4))[Iterator].SetPtEtaPhiM(smeared_pt, (*(myEvent->muons_p4))[Iterator].Eta(), (*(myEvent->muons_p4))[Iterator].Phi(), (*(myEvent->muons_p4))[Iterator].M());

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

		if(Flag_Smear_Jet_Scale){
			//Smear_Jet_Scale is a random number pulled from a standard normal distribution with mean 0.0 and variance 1.0
			if(isData || firstLoop){
				//smear resolution of jets in data based on jet energy correction uncertainty
				(*(myEvent->jets_p4))[Iterator] = (1 + (Smear_Jet_Scale) * (*(myEvent->jec_uncertainty))[Iterator]) * (*(myEvent->jets_p4))[Iterator];
			}
			else{
				//skip jets with very low pT and very high pT, as instructed on Twiki titled
				//"Jet Energy Resolution: Official Software Tools for retrieving JER and scale factors"
				if((*(myEvent->jets_p4))[Iterator].Pt() < 15.0 || std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 2.40){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 2.30 && (*(myEvent->jets_p4))[Iterator].Pt() > 806.4){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 2.10 && (*(myEvent->jets_p4))[Iterator].Pt() > 831.5){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.90 && (*(myEvent->jets_p4))[Iterator].Pt() > 1090.){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.70 && (*(myEvent->jets_p4))[Iterator].Pt() > 1126.){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.30 && (*(myEvent->jets_p4))[Iterator].Pt() > 1648.){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.10 && (*(myEvent->jets_p4))[Iterator].Pt() > 2175.){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.80 && (*(myEvent->jets_p4))[Iterator].Pt() > 2651.){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.50 && (*(myEvent->jets_p4))[Iterator].Pt() > 3155){
					Iterator++;
					continue;
				}
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.0 && (*(myEvent->jets_p4))[Iterator].Pt() > 3197.){
					Iterator++;
					continue;
				}

				//////////////////////////////////////////
				//smear resolution of jets in MC based on data/MC scale SF and jet energy resolution vs eta
				//////////////////////////////////////////
#ifdef DEBUG
				std::cout<<"\t"<< std::endl;
#endif
				//first determine the PT scale factor which was already applied to MC, and the uncertainty on this SF
				//this SF accounts for different jet energy resolutions in data and MC (worse in data)

#ifdef DEBUG	
				std::cout<<"before smearing jet pt = "<< ((*(myEvent->jets_p4))[Iterator].Pt()) << " and jet eta = "<< ((*(myEvent->jets_p4))[Iterator].Eta()) << std::endl;
#endif
				double ptSF = 1.00, ptSFUnc = 0.0;
				double pt = (*(myEvent->jets_p4))[Iterator].Pt();
				if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) <= 0.80) ptSF = 1.061, ptSFUnc = 0.023;
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) <= 1.30) ptSF = 1.088, ptSFUnc = 0.029;
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) <= 1.90) ptSF = 1.106, ptSFUnc = 0.030;
				else if(std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) <= 2.50) ptSF = 1.126, ptSFUnc = 0.094;
				ptSF += ptSFUnc*Smear_Jet_Scale;
				double ptSFSqd = std::fmax((ptSF*ptSF - 1.0)/pt, 0.0);
#ifdef DEBUG
				std::cout<<"ptSFSqd = "<< ptSFSqd << std::endl;
#endif
				
				//now determine the jet energy resolution uncertainty for the jet based on its eta and rho
				double jerSigma = 1.0;
				if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 2.30 ){
					jerSigma = std::sqrt((std::pow(4.341,2)/(pt*pt)) + (std::pow(0.8304,2)/(std::pow(pt,-1.043))) + std::pow(0.04588,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 2.10 ){
					jerSigma = std::sqrt((std::pow(4.05,2)/(pt*pt)) + (std::pow(0.7653,2)/(std::pow(pt,-1.016))) + std::pow(0.04598,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.90 ){
					jerSigma = std::sqrt((std::pow(3.402,2)/(pt*pt)) + (std::pow(0.9552,2)/(std::pow(pt,-1.062))) + std::pow(0.04044,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.70 ){
					jerSigma = std::sqrt((std::pow(2.093,2)/(pt*pt)) + (std::pow(1.419,2)/(std::pow(pt,-1.178))) + std::pow(0.04109,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.30 ){
					jerSigma = std::sqrt((std::pow(2.831,2)/(pt*pt)) + (std::pow(1.089,2)/(std::pow(pt,-1.01))) + std::pow(0.04586,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 1.10 ){
					jerSigma = std::sqrt((std::pow(3.065,2)/(pt*pt)) + (std::pow(0.6603,2)/(std::pow(pt,-0.8076))) + std::pow(0.04368,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.80 ){
					jerSigma = std::sqrt((std::pow(3.298,2)/(pt*pt)) + (std::pow(0.415,2)/(std::pow(pt,-0.6451))) + std::pow(0.02981,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.50 ){
					jerSigma = std::sqrt((std::pow(3.101,2)/(pt*pt)) + (std::pow(0.3634,2)/(std::pow(pt,-0.6054))) + std::pow(0.02263,2) );
				}
				else if( std::fabs((*(myEvent->jets_p4))[Iterator].Eta()) >= 0.0 ){
					jerSigma = std::sqrt((std::pow(3.096,2)/(pt*pt)) + (std::pow(0.3558,2)/(std::pow(pt,-0.5949))) + std::pow(0.01943,2) );
				}

#ifdef DEBUG
				std::cout<<"jerSigma = "<< jerSigma << std::endl;
#endif

				//now compute the correction factor C to apply to the jet pT
				TRandom3 Rand;
				double C = Rand.Gaus(0.0, jerSigma)*sqrt(ptSFSqd);

#ifdef DEBUG
				std::cout<<"correction to jet pT = "<< C << std::endl;
#endif

				//smear the jet pt
				double smearedJetPt = (1+C)*((*(myEvent->jets_p4))[Iterator].Pt());
				(*(myEvent->jets_p4))[Iterator].SetPtEtaPhiM(smearedJetPt, (*(myEvent->jets_p4))[Iterator].Eta(), (*(myEvent->jets_p4))[Iterator].Phi(), (*(myEvent->jets_p4))[Iterator].M());


			}
		}


#ifdef DEBUG
		std::cout << std::endl << " Jet number= " << Iterator << " Jet Pt Before = " << (*(myEvent->jets_p4))[Iterator].Pt() << " Jet Eta After = " << (*(myEvent->jets_p4))[Iterator].Eta() << std::endl;
#endif

		Iterator++;
	}

	delete rmcor;
}

