#include "makeLimitFile.hh"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
  if (argc<5) {
    printf("Usage: makeLimitFile2010 [mwr] [data file] [signal file] [outfilename]\n");
    return 1;
  }
  int mwr=atoi(argv[1]);
  TFile df(argv[2]);
  TFile sf(argv[3]);

  makeLimitFile2010(mwr, &df,&sf,argv[4]);

  return 0;
}
