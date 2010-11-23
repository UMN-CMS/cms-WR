/*
Author : Jeremy Mans
Purpose : Analyze a results file and return #back in bins for a signal efficiency of <given>
*/
#define __USE_XOPEN2K8
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>

#define MAXSIGNAL 100000
#define MAXBACKGROUND 400000

int compare_fn(const void* a, const void* b) {
   return (*((float*)a)>(*((float*)b)))?1:-1;
}

int main(int argc, char* argv[]) {
  float eff=0.5f, intended, was, req=0.0;
  FILE* f;
  char buffer[4096];
  unsigned int i,nb,ns;

  float* signal,*background;

  signal=new float[MAXSIGNAL];
  background=new float[MAXBACKGROUND];

  if (argc==1 || argc>4) {
    printf("Usage of nnbackeff:\n");
    printf("nnbackeff [resultsfile.res] {eff}\n");
    printf("  The efficiency value is optional and defaults to 0.5\n");
    return 1;
  }
  if (argc>=3) eff=atof(argv[2]);
  if (argc>=4) req=atof(argv[3]);

  if ((f=fopen(argv[1],"rt"))==NULL) {
    printf("Cannot open file....\n");
    return 2;
  }


  nb=0;
  ns=0;

  while (!feof(f)) {
    fgets(buffer,4096,f);
    buffer[strlen(buffer)-1]=0; // trim off the \n
    if (buffer[0]=='#') { // beginning of a set
      fgets(buffer,4096,f);
      intended=atof(buffer);
      fgets(buffer,4096,f);
      was=atof(buffer);
      if (intended<0.25) background[nb++]=was;
      else if (intended>0.75) signal[ns++]=was;
    }
    buffer[0]=0;
  }

  qsort(background,nb,sizeof(float),compare_fn);
  qsort(signal,ns,sizeof(float),compare_fn);

  float cpv=signal[int(float(ns)*(1.0f-eff))];
  for (i=0; i<nb && background[i]<cpv; i++);

  if (signal[9*ns/10]<req) i=0;

  printf("%d background events in the %f efficiency region\n",nb-i,eff);
  return 0;
}
