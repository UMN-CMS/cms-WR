#include "makeLimitFile.hh"
#include "systematics.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<6) {
    printf("Usage: makeLimitFile2011 [lumi] [mwr] [data file] [signal file] [outfilename] [optional special syst mode]\n");
    return 1;
  }
  float lumi=atof(argv[1]);
  int mwr=atoi(argv[2]);
  TFile df(argv[3]);
  TFile sf(argv[4]);

  if (argc==7)
    special_syst_mode=atoi(argv[6]);

  makeLimitFile2011(lumi,mwr,&df,&sf,argv[5]);

  return 0;
}
