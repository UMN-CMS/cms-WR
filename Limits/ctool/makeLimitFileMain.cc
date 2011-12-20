#include "makeLimitFile.hh"
#include "systematics.h"
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"

int main(int argc, char* argv[]) {
  if (argc<9) {
    printf("Usage: makeLimitFile [lumi] [mwr] [mn] [xsec] [data file] [signal file] [outfilename] [systematicsdb] [optional special syst mode]\n");
    return 1;
  }
  float lumi=atof(argv[1]);
  MassPoint pt;
  pt.mwr=atoi(argv[2]);
  pt.mnr=atoi(argv[3]);
  float xsec=atof(argv[4]);
  TFile df(argv[5]);
  TFile sf(argv[6]);

  SystematicsDB syst;
  syst.load(argv[8]);
  syst.standardSystematics(); // define the standard systematics
  if (argc==10)
    int special_syst_mode=atoi(argv[9]);

  

  makeLimitFile(lumi,xsec,pt,&df,&sf,argv[7],syst);

  return 0;
}
