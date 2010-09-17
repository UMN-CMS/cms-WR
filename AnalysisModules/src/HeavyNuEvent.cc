#include "HeavyNu/AnalysisModules/src/HeavyNuEvent.h"

void HeavyNuEvent::reset() {
  mu1=0;
  mu2=0;
  mu[0]=0;
  mu[1]=0;
  j1=0;
  j2=0;
  j[0]=0;
  j[1]=0;
}


void HeavyNuEvent::regularize() {
  mu[0]=mu1;
  mu[1]=mu2;
  j[0]=j1;
  j[1]=j2;
}
