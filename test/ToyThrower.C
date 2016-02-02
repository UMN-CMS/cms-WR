#include<fstream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "../plugins/rochcor2015.h"
#include "../plugins/muresolution_run2.h"
using namespace std;

void ToyThrower(miniTreeEvent myEvent, int RandomSeed){
      int DEBUG = 0;
      float qter = 1.0;
      TRandom3* Rand=new TRandom3(); 
      Rand->SetSeed(RandomSeed);
      double Smear = 0.0;
      int Iterator = 0;
//      rochcor2015 *rmcor = new rochcor2015(RandomSeed);
      for(auto muons:*(myEvent.muons_p4)){
          Smear=Rand->Gaus(0.0,1);
          if(DEBUG) cout<<endl<<" Muon number= "<<Iterator<<" Muon Pt Before = "<<(*(myEvent.muons_p4))[Iterator].Pt()<<" Muon Eta Before = "<<(*(myEvent.muons_p4))[Iterator].Eta()<<endl;
          if(Smear >= 0.) (*(myEvent.muon_IDSF_central))[Iterator] += Smear*(*(myEvent.muon_IDSF_error))[Iterator];  
          else            (*(myEvent.muon_IDSF_central))[Iterator] -= Smear*(*(myEvent.muon_IDSF_error))[Iterator];

          Smear=Rand->Gaus(0.0,1);
          if(Smear >= 0.) (*(myEvent.muon_IsoSF_central))[Iterator] += Smear*(*(myEvent.muon_IsoSF_error))[Iterator];                            
          else            (*(myEvent.muon_IsoSF_central))[Iterator] -= Smear*(*(myEvent.muon_IsoSF_error))[Iterator];

//          rmcor->momcor_mc((*(myEvent.muons_p4))[Iterator], (*(myEvent.muon_charge))[Iterator], 0, qter);
          if(DEBUG) cout<<endl<<" Muon number= "<<Iterator<<" Muon Pt After = "<<(*(myEvent.muons_p4))[Iterator].Pt()<<" Muon Eta After = "<<(*(myEvent.muons_p4))[Iterator].Eta()<<endl;
          Iterator++;
        }

      Iterator = 0;

      for(auto electrons:*(myEvent.electrons_p4)){
          Smear=Rand->Gaus(0.0,1);
          if(DEBUG) cout<<endl<<" Electron number= "<<Iterator<<" Electron Pt Before = "<<(*(myEvent.electrons_p4))[Iterator].Pt()<<" Electron Eta Before = "<<(*(myEvent.electrons_p4))[Iterator].Eta()<<endl;
          (*(myEvent.electrons_p4))[Iterator] = (1 + (Smear)*(*(myEvent.electron_smearing))[Iterator])*(*(myEvent.electrons_p4))[Iterator]; 
          if(DEBUG) cout<<endl<<" Electron number= "<<Iterator<<" Electron Pt After = "<<(*(myEvent.electrons_p4))[Iterator].Pt()<<" Electron Eta After = "<<(*(myEvent.electrons_p4))[Iterator].Eta()<<endl;
          Iterator++;  
        }

      Iterator = 0;

      for(auto jets:*(myEvent.jets_p4)){
          Smear=Rand->Gaus(0.0,1);
          if(DEBUG) cout<<endl<<" Jet number= "<<Iterator<<" Jet Pt Before = "<<(*(myEvent.jets_p4))[Iterator].Pt()<<" Jet Eta Before = "<<(*(myEvent.jets_p4))[Iterator].Eta()<<endl;
          (*(myEvent.jets_p4))[Iterator] = (1 + (Smear)*(*(myEvent.jec_uncertainty))[Iterator])*(*(myEvent.jets_p4))[Iterator];
          if(DEBUG) cout<<endl<<" Jet number= "<<Iterator<<" Jet Pt Before = "<<(*(myEvent.jets_p4))[Iterator].Pt()<<" Jet Eta After = "<<(*(myEvent.jets_p4))[Iterator].Eta()<<endl;
          Iterator++;
        }

//      delete rmcor;
}

