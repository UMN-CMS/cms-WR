#include "makeLimitFile.hh"
#include "systematics.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<7) {
    printf("Usage: makeLimitFile2011 [lumi] [xsec] [mwr] [data file] [signal file] [outfilename] [optional special syst mode]\n");
    return 1;
  }
  float lumi=atof(argv[1]);
  float xsec=atof(argv[2]);
  int mwr=atoi(argv[3]);
  TFile df(argv[4]);
  TFile sf(argv[5]);

  if (argc==8)
    special_syst_mode=atoi(argv[7]);

  makeLimitFile2011(lumi,xsec,mwr,&df,&sf,argv[6]);

  return 0;
}
