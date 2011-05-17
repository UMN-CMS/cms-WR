#include "makeLimitFile.hh"
#include <stdio.h>

int main(int argc, char* argv[]) {
  if (argc<4) {
    printf("Usage: makeLimitFile2010 [data file] [signal file] [outfilename]\n");
    return 1;
  }
  TFile df(argv[1]);
  TFile sf(argv[2]);

  makeLimitFile2010(&df,&sf,argv[3]);

  return 0;
}
