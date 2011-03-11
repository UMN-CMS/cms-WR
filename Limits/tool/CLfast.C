#include <math.h>
#include <TRandom3.h>
#include <algorithm>

template <class T>
inline const T& MAXof(const T& a, const T& b) { return (a>b)?a:b; }

struct CLpoint {
  
  /// observed confidence level for signal (from data)
  double clobs;
  /// background confidence level
  double clb;
  /// median expected confidence level in the presence of signal
  double clmed;
  /// +1 sigma deviation in the median expected confidence level in the presence of signal
  double clmed_p1s;
  /// +2 sigma deviation in the median expected confidence level in the presence of signal
  double clmed_p2s;
  /// -1 sigma deviation in the median expected confidence level in the presence of signal
  double clmed_m1s;
  /// -2 sigma deviation in the median expected confidence level in the presence of signal
  double clmed_m2s;

  /// observed log-likelihood ratio (from data)
  double llrobs;
  /// expected log-likelihood ratio in the absence of signal
  double llrb;
  /// +1 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_p1s;
  /// +2 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_p2s;
  /// -1 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_m1s;
  /// -2 sigma deviation in the expected log-likelihood ratio in the absence of signal
  double llrb_m2s;

  /// expected log-likelihood ratio in the presence of signal
  double llrsb;
  /// +1 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_p1s;
  /// +2 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_p2s;
  /// -1 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_m1s;
  /// -2 sigma deviation in the expected log-likelihood ratio in the presence of signal
  double llrsb_m2s;

  /// confidence level sought for in CrossSectionLimit code
  double xsec_cl;
  /// factor which the signal must be scaled by to obtain an observed limit of [xsec_cl]
  double xsec_obsfactor;
  /// factor which the signal must be scaled by to obtain a median expected limit of [xsec_cl]
  double xsec_medfactor;
};

double calculateLLR(double signal, double bkgd, double data) {
  double llr;
  if (bkgd<1e-6) llr=-signal+data*log(1+signal*1e6);
  else llr=-signal+data*log(1+signal/bkgd);

  double ztests,npreds,npredb,ls,lb;

  ztests=data;
  npreds=std::max(signal+bkgd,1e-6);
  npredb=std::max(bkgd,1e-6);

  ls=-npreds+ztests*log(npreds);
  lb=-npredb+ztests*log(npredb);
  llr=lb-ls;

  return llr;
}

void CLfast(double signal, double bkgd, double data, CLpoint& CLs, int its=100000, TRandom* aRand=0) {
  static const int MaxKeep = (1024*1024)/4/2;
  static float btrial[MaxKeep];
  static float sbtrial[MaxKeep];

  TRandom* m_rand=aRand;

  if (aRand==0) m_rand=gRandom;
  
  double llr_d, llr_b_expected, llr_sb_expected, llr_sb, llr_b;
  double sblike, blike;
  
  int numer=0;
  int denom=0;
  
  llr_d=calculateLLR(signal,bkgd,data);

  sblike=signal+bkgd;
  blike=bkgd;

  llr_sb_expected=calculateLLR(signal,bkgd,sblike);
  llr_b_expected=calculateLLR(signal,bkgd,blike);

  CLs.llrobs=llr_d;
  CLs.llrb=llr_b_expected;
  CLs.llrsb=llr_sb_expected;
  
  int i;
  double zback,zsig,zback2;
  for (i=0; i<5*its && denom<its; i++) {
    if (bkgd>0) {
      zback=m_rand->Poisson(bkgd);
      zback2=zback;
    } else { zback=0; zback2=0; }
    if (signal>0) {
      zsig=m_rand->Poisson(signal);
    } else zsig=0;
    sblike=zback+zsig;
    blike=zback2;

    llr_sb=calculateLLR(signal,bkgd,sblike);
    llr_b=calculateLLR(signal,bkgd,blike);

    //    printf("%d %.0f %.0f %.0f %.0f %4.2f %4.2f %4.2f\n",i,zback,zsig,sblike,blike,llr_sb, llr_b, llr_d);
    

    
    if (llr_d<=llr_b) { // background "less than" data
      denom++;
      if (llr_d<=llr_sb) numer++; // s+b "less than" data
    }
    if (i<MaxKeep) {
      btrial[i]=llr_b;
      sbtrial[i]=llr_sb;
    }
  }
  
  int nmcdone=i;

  //printf("%d %d \n",numer,denom);
  
  CLs.clb=(double(denom)/double(nmcdone));
  double clobs=1.0;
  if (denom>0) clobs=1.0-double(numer)/double(denom); // fraction of trials where data < s+b
  if (signal<0.001) clobs=0.0;
  CLs.clobs=clobs;
  
  int lim=(nmcdone<MaxKeep)?nmcdone:MaxKeep;
  
  std::sort(btrial,btrial+lim);
  std::sort(sbtrial,sbtrial+lim);
  
  CLs.clmed=-1.0;
  CLs.clmed_m1s=-1.0;
  CLs.clmed_m2s=-1.0;
  CLs.clmed_p1s=-1.0;
  CLs.clmed_p2s=-1.0;

  int m1s_pos = int(lim*(1.0-0.3173/2.0)+0.5)+1;
  int p1s_pos = int(lim*(    0.3173/2.0)+0.5)+1;
  int m2s_pos = int(lim*(1.0-0.0455/2.0)+0.5)+1;
  int p2s_pos = int(lim*(    0.0455/2.0)+0.5)+1;

  CLs.llrsb_m2s=sbtrial[m2s_pos];
  CLs.llrsb_m1s=sbtrial[m1s_pos];
  CLs.llrsb_p1s=sbtrial[p1s_pos];
  CLs.llrsb_p2s=sbtrial[p2s_pos];

  CLs.llrb_m2s=btrial[m2s_pos];
  CLs.llrb_m1s=btrial[m1s_pos];
  CLs.llrb_p1s=btrial[p1s_pos];
  CLs.llrb_p2s=btrial[p2s_pos];
  
  float val=btrial[lim/2];
  for (i=0; i<lim; i++) {
    if (sbtrial[i]>=val && CLs.clmed<0) CLs.clmed=MAXof(0.0,float(i)/(lim/2)-1.0);
    if (sbtrial[i]>=btrial[m1s_pos] && CLs.clmed_m1s<0) CLs.clmed_m1s=MAXof(0.0,float(i)/(lim/2)-1.0);
    if (sbtrial[i]>=btrial[m2s_pos] && CLs.clmed_m2s<0) CLs.clmed_m2s=MAXof(0.0,float(i)/(lim/2)-1.0);
    if (sbtrial[i]>=btrial[p1s_pos] && CLs.clmed_p1s<0) CLs.clmed_p1s=MAXof(0.0,float(i)/(lim/2)-1.0);
    if (sbtrial[i]>=btrial[p2s_pos] && CLs.clmed_p2s<0) CLs.clmed_p2s=MAXof(0.0,float(i)/(lim/2)-1.0);
  }
  if (CLs.clmed==-1.0) CLs.clmed=1.0;
  if (CLs.clmed_m1s==-1.0) CLs.clmed_m1s=1.0;
  if (CLs.clmed_m2s==-1.0) CLs.clmed_m2s=1.0;
  if (CLs.clmed_p1s==-1.0) CLs.clmed_p1s=1.0;
  if (CLs.clmed_p2s==-1.0) CLs.clmed_p2s=1.0;
  
  if (CLs.clmed>1.0 && CLs.clmed<1.001) CLs.clmed=1.0;
  if (CLs.clmed<0.0) CLs.clmed=0.0;
}

static double func(double sig, double back) {
  static CLpoint pt;
  CLfast(sig,back,0,pt,10000);
  return pt.clmed-0.95;
}

static void bracket(double& x1, double& x2, double sig, double back) {

  for (int i=0; i<100; i++) {
    double value=func(sig*x2,back);
    //    printf("%d %f\n",i,value);
    if (x2>1.0) {
      if (value>0.0) break;
      x1=x2;
      x2*=2.0;
    } else {
      if (value<0.0) break;
      x1=x2;
      x2/=2.0;
    }
  }
  if (x1>x2) {
    x1*=1.5;
    x2/=1.5;
  } else {
    x1/=1.5;
    x2*=1.5;
  }
}

#define SIGN(a,b) ((b)>=0.0 ? ((a)<0.0?-(a):(a)):((a)<0.0?(a):-(a)))


static float zriddr(float x1, float x2, double signal, double bkgd, float acc) {
  int j;
  float ans, fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  
  fl=func(x1*signal,bkgd);
  fh=func(x2*signal,bkgd);
  if ((fl>0.0f && fh<0.0f) || (fl< 0.0f && fh >0.0f)) {
    xl=x1;
    xh=x2;
    xm=0.5*(xl+xh);
    ans=xm;
    for (j=0; j<60; j++) {
      xm=0.5*(xl+xh);
      fm=func(xm*signal,bkgd);
      s=sqrt(fm*fm-fl*fh);
      //      printf("   CBR: got %.2f%% at %.4e (%f)\n",fm*100.0,xm,s);
      if (s==0.0) return ans;
      xnew=xm+(xm-xl)*((fl>=fh?1.0:-1.0)*fm/s);
      //                      if (fabs(xnew-ans)<=xacc) return ans;
      if (fabs(fm)<=acc) return ans;
      ans=xnew;
      fnew=func(ans*signal,bkgd);
      if (fnew==0.0) return ans;
      if (SIGN(fm,fnew)!=fm) {
	xl=xm;
	fl=fm;
	xh=ans;
	fh=fnew;
      } else if (SIGN(fl,fnew)!=fl) {
	xh=ans;
	fh=fnew;
      } else if (SIGN(fh,fnew)!=fh) {
	xl=ans;
	fl=fnew;
      } else printf("huh?");
      //                      if (fabs(xh-xl)<=xacc) return ans;
      if (fabs(fnew)<=acc) return ans;
    }
  }
  return -1.0;
}


double CLfast_goalSeek(double signal, double background) {
  double high_factor=2.0;
  double low_factor=0.5;

  
  while (func(signal*low_factor,background)<0)  {    
    //    printf("LF %f %f\n",low_factor,func(signal*low_factor,background));
    low_factor*=2.0;
  }

  while (func(signal*high_factor,background)>0)  {    
    //    printf("HF %f %f\n",high_factor,func(signal*high_factor,background));
    high_factor/=2.0;
  }


  bracket(low_factor,high_factor,signal,background);
  /*
  printf("brackets: %f %f %f %f\n",low_factor,high_factor,
	 func(signal*low_factor,background),
	 func(signal*high_factor,background)
	 );
  */
  double factor=zriddr(low_factor,high_factor,signal,background,0.01);
  //  printf("%f %f\n",factor,func(signal*factor,background));

  return factor;
}
