#include "makeLimitFile.hh"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<7) {
    printf("Usage: makeLimitFile2Y [lumi 2011] [data file 2011] [signal file 2011] [data file 2010] [signal file 2010] [outfilename]\n");
    return 1;
  }
  float lumi=atof(argv[1]);
  TFile df11(argv[2]);
  TFile sf11(argv[3]);
  TFile df10(argv[4]);
  TFile sf10(argv[5]);

  makeLimitFileTwoYear(lumi,&df11,&sf11,&df10,&sf10,argv[6]);

  return 0;
}
