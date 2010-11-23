/*
Author : Jeremy Mans
Purpose : Analyze a results file and return the difference in LLRs for the signal and background
*/

#define __USE_XOPEN2K8
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NBINS 100

const float left=0.0;
const float right=1.0;

#define MAX(a,b) ((a)>(b)?(a):(b))

double LLRcalculate(float* back, float* sig, double& blike, double& slike) {

/*            
*     Computes Log(Likelihood Ratio) b/(s+b) estimator from binned histograms.
*     This is essentially the Delphi estimator and should be fast.
*
*     Input:   
      real back(*)            ! Array of expected backround events per bin.
      real sig(*)             ! Array of expected signal events per bin.
*     Output:  
      double precision slike  ! (Log Likelihood ratio) that thrown s+b vs s+b expectation.
      double precision blike  ! (Log Likelihood) that thrown b fits b expectation.
*     The likelihood ratio is used in clfast to order s+b, b, d...         
*/
	double ztests,npreds,ztestb,npredb,ls,lb;
	int ib;

	slike=0;
	blike=0;
	for(ib=0; ib<NBINS; ib++) {
	       if (back[ib]<1e-8) continue;
		ztests=sig[ib]+back[ib];
		npreds=MAX(ztests,1.0e-6);
		ztestb=back[ib];
		npredb=MAX(back[ib],1.0e-6);
		ls=-npreds+ztests*log(npreds);
		lb=-npredb+ztests*log(npredb);
		slike+=lb-ls;
		ls=-npreds+ztestb*log(npreds);
		lb=-npredb+ztestb*log(npredb);
		blike+=lb-ls;
	}
	return slike-blike;
}

void Cleanup(float* b, float* s) {
  int i=0,ilast=-1;
  float cs=0;
  while (b[i]==0 && i<NBINS) {
    cs+=s[i];
    s[i]=0;
    i++;    
  }
  
  for (; i<NBINS; i++) {
    if (b[i]>0) {
      s[i]+=cs;
      cs=0;
      ilast=i;
    } else {
      cs+=s[i];
      s[i]=0;
    }
  }
  if (ilast!=NBINS-1 && ilast>=0) {
    s[ilast]+=cs;
  }
}

int main(int argc, char* argv[]) {
  float intended, was;
  FILE* f;
  char buffer[4096];
  unsigned int i,nb,ns;
  float sweight=1.0, bweight=1.0;

  float* signal,*background;

  signal=new float[NBINS];
  background=new float[NBINS];

  if (argc==1 || argc>4) {
    printf("Usage of nnllr:\n");
    printf("nnllr.exe [resultsfile.res] {signal weight} {background weight}\n");
    printf("  The efficiency value is optional and defaults to 0.5\n");
    return 1;
  }
  if (argc>=3) sweight=atof(argv[2]);
  if (argc==4) bweight=atof(argv[3]);


  if ((f=fopen(argv[1],"rt"))==NULL) {
    printf("Cannot open file....\n");
    return 2;
  }

  for (i=0; i<NBINS; i++) {
    signal[i]=0;
    background[i]=0;
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
      if (intended<0.25) {
	if (was<left) background[0]+=bweight;
	else if (was>=right) background[NBINS-1]+=bweight;
	else background[int((was-left)/(right-left)*float(NBINS))]+=bweight;
      }
      else if (intended>0.75) {
	if (was<left) signal[0]+=sweight;
	else if (was>=right) signal[NBINS-1]+=sweight;
	else signal[int((was-left)/(right-left)*float(NBINS))]+=sweight;
      }
    }
    buffer[0]=0;
  }

  Cleanup(background,signal);

  double llrsb,llrb;
  LLRcalculate(background, signal, llrb, llrsb);

  printf("%f delta LLR (%f and %f), weight = s=%f,b=%f\n",llrb-llrsb,llrb,llrsb,sweight, bweight);
  return 0;
}
