#include "makeLimitFile.hh"
#include <stdio.h>
#include "systematics.h"
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<8) {
    printf("Usage: makeLimitFile2Y [mwr] [lumi 2011] [data file 2011] [signal file 2011] [data file 2010] [signal file 2010] [outfilename] [optional special syst mode]\n");
    return 1;
  }
  float lumi=atof(argv[2]);
  int mwr=atoi(argv[1]);
  TFile df11(argv[3]);
  TFile sf11(argv[4]);
  TFile df10(argv[5]);
  TFile sf10(argv[6]);

  if (argc==9)
    special_syst_mode=atoi(argv[8]);


  makeLimitFileTwoYear(lumi,mwr,&df11,&sf11,&df10,&sf10,argv[7]);

  return 0;
}
