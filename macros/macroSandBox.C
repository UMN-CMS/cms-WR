#include <TFile.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TObjArray.h>
#include "TCollection.h"
#include <TCut.h>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <TEventList.h>
#include <TEntryList.h>
#include <TEntryListArray.h>
#include "dumpTreePlots.C"

using namespace std;

#define DEBUG


///use this macro to develop and run new macros which don't have a central theme, other than being useful to the WR analysis
///the first use of this macro was to plot reco pT/gen pT for reco jets and leptons matched to GEN counterparts

void macroSandBox(){

}///end macroSandBox()

