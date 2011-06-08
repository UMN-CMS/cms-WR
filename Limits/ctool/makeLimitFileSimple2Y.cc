#include "makeLimitFileSimple.hh"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<8) {
    printf("Usage: makeLimitFileSimple2Y [mwr] [minbin] [lumi 2011] [data file 2011] [signal file 2011] [data file 2010] [signal file 2010] [outfilename]\n");
    return 1;
  }
  int mwr=atoi(argv[1]);
  int minbin=atoi(argv[2]);
  float lumi=atof(argv[3]);
  TFile df11(argv[4]);
  TFile sf11(argv[5]);
  TFile df10(argv[6]);
  TFile sf10(argv[7]);

  makeLimitFileTwoYear(lumi,minbin,mwr,&df11,&sf11,&df10,&sf10,argv[8]);

  return 0;
}
