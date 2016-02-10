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

#ifdef DEBUG
#include <iostream>
#endif

void ToyThrower(miniTreeEvent myEvent,  float rand_smear[], float rand_up_down[], int random_seed, std::vector<std::string> list, int isData){

      float qter = 1.0;
      double Smear = 0.0;
      int Iterator = 0;
      int Iterator_Up_Down = 0;
      int Iterator_Smear = 0;  
      int Flag_Smear_Muon_Scale=0,Flag_Smear_Muon_ID_Iso=0,Flag_Smear_Electron_Scale=0,Flag_Smear_Jet_Scale=0;

      for(unsigned int iii=0;iii<list.size();iii++){
              if(list[iii] == "Smear_Muon_Scale")     Flag_Smear_Muon_Scale=1;
         else if(list[iii] == "Smear_Muon_ID_Iso")    Flag_Smear_Muon_ID_Iso=1;
         else if(list[iii] == "Smear_Electron_Scale") Flag_Smear_Electron_Scale=1;
         else if(list[iii] == "Smear_Jet_Scale")      Flag_Smear_Jet_Scale=1;
        }

      rochcor2015 *rmcor = new rochcor2015(random_seed);

      for(auto muons:*(myEvent.muons_p4)){
#ifdef DEBUG
		  std::cout<<std::endl<<" Muon number= "<<Iterator<<" Muon Pt Before = "<<(*(myEvent.muons_p4))[Iterator].Pt()<<" Muon Eta Before = "<<(*(myEvent.muons_p4))[Iterator].Eta()<<std::endl;
#endif

          if(Flag_Smear_Muon_ID_Iso && !isData ){
             Smear = rand_up_down[Iterator_Up_Down++];
             if(Smear >= 0.) (*(myEvent.muon_IDSF_central))[Iterator] += Smear*(*(myEvent.muon_IDSF_error))[Iterator];  
             else            (*(myEvent.muon_IDSF_central))[Iterator] -= Smear*(*(myEvent.muon_IDSF_error))[Iterator];

             Smear=rand_up_down[Iterator_Up_Down++];
             if(Smear >= 0.) (*(myEvent.muon_IsoSF_central))[Iterator] += Smear*(*(myEvent.muon_IsoSF_error))[Iterator];                            
             else            (*(myEvent.muon_IsoSF_central))[Iterator] -= Smear*(*(myEvent.muon_IsoSF_error))[Iterator];
            }

          if(Flag_Smear_Muon_Scale){
             if(isData) 
                  rmcor->momcor_data((*(myEvent.muons_p4))[Iterator], (*(myEvent.muon_charge))[Iterator], 0, qter);
             else rmcor->momcor_mc((*(myEvent.muons_p4))[Iterator], (*(myEvent.muon_charge))[Iterator], 0, qter);

            
            }

#ifdef DEBUG
		  std::cout<<std::endl<<" Muon number= "<<Iterator<<" Muon Pt After = "<<(*(myEvent.muons_p4))[Iterator].Pt()<<" Muon Eta After = "<<(*(myEvent.muons_p4))[Iterator].Eta()<<std::endl;
#endif
          Iterator++;
        }

      Iterator = 0;

      for(auto electrons:*(myEvent.electrons_p4)){
#ifdef DEBUG
		  std::cout<<std::endl<<" Electron number= "<<Iterator<<" Electron Pt Before = "<<(*(myEvent.electrons_p4))[Iterator].Pt()<<" Electron Eta Before = "<<(*(myEvent.electrons_p4))[Iterator].Eta()<<std::endl;
#endif

          if(Flag_Smear_Electron_Scale){
            if(isData){
              Smear=rand_up_down[Iterator_Up_Down++];
              (*(myEvent.electrons_p4))[Iterator] = (1 + (Smear)*(*(myEvent.electron_scale))[Iterator])*(*(myEvent.electrons_p4))[Iterator];
             }
            else{
                 Smear = rand_smear[Iterator_Smear++];  
                 (*(myEvent.electrons_p4))[Iterator] = (1 + (Smear)*(*(myEvent.electron_smearing))[Iterator])*(*(myEvent.electrons_p4))[Iterator];
                  } 
      }
#ifdef DEBuG
		  std::cout << std::endl<<" Electron number= "<<Iterator<<" Electron Pt After = "<<(*(myEvent.electrons_p4))[Iterator].Pt()<<" Electron Eta After = "<<(*(myEvent.electrons_p4))[Iterator].Eta()<<std::endl;
#endif
          Iterator++;  
        }

      Iterator = 0;

      for(auto jets:*(myEvent.jets_p4)){
          Smear = rand_smear[Iterator_Up_Down++];
#ifdef DEBUG
		  std::cout<<std::endl<<" Jet number= "<<Iterator<<" Jet Pt Before = "<<(*(myEvent.jets_p4))[Iterator].Pt()<<" Jet Eta Before = "<<(*(myEvent.jets_p4))[Iterator].Eta()<<std::endl;
#endif
          if(Flag_Smear_Jet_Scale)
            (*(myEvent.jets_p4))[Iterator] = (1 + (Smear)*(*(myEvent.jec_uncertainty))[Iterator])*(*(myEvent.jets_p4))[Iterator];
#ifdef DEBUG
		  std::cout<<std::endl<<" Jet number= "<<Iterator<<" Jet Pt Before = "<<(*(myEvent.jets_p4))[Iterator].Pt()<<" Jet Eta After = "<<(*(myEvent.jets_p4))[Iterator].Eta()<<std::endl;
#endif

          Iterator++;
        }

  delete rmcor;
}
