#include "../interface/RooDoubleCB.h"
#include "TMath.h"

///constructor from user defined args
///two crystal ball fxns with the same mean but different sigma, alpha, and power
RooDoubleCB::RooDoubleCB(const char * name, const char * title,
		RooAbsReal& _x, RooAbsReal& _mean,
		RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		RooAbsReal& _alphaL, RooAbsReal& _alphaR,
		RooAbsReal& _powerL, RooAbsReal& _powerR
		):
	RooAbsPdf(name,title),
	x("x","Observable", this, _x),
	mean("mean","Mean", this, _mean),
	sigmaL("sigmaL","sigmaL", this, _sigmaL),
	sigmaR("sigmaR","sigmaR", this, _sigmaR),
	alphaL("alphaL","alphaL", this, _alphaL),
	alphaR("alphaR","alphaR", this, _alphaR),
	powerL("powerL","powerL", this, _powerL),
	powerR("powerR","powerR", this, _powerR){

	}


///constructor using existing RooDoubleCB object - needed for clone() method
RooDoubleCB::RooDoubleCB(const RooDoubleCB& other, const char* name):
	RooAbsPdf(other,name), 
	x("x", this, other.x),
	mean("mean", this, other.mean),
	sigmaL("sigmaL", this, other.sigmaL),
	sigmaR("sigmaR", this, other.sigmaR),
	alphaL("alphaL", this, other.alphaL),
	alphaR("alphaR", this, other.alphaR),
	powerL("powerL", this, other.powerL),
	powerR("powerR", this, other.powerR){

	}

///fxn which is called automatically by RooFit methods like fitTo
Double_t RooDoubleCB::evaluate() const{
	Double_t absAlphaL = TMath::Abs((Double_t) alphaL);
	Double_t absAlphaR = TMath::Abs((Double_t) alphaR);
	Double_t a, b;	///<used in power law tail

	Double_t s = (x-mean), t, alpha, power;
	if(s<0){//to the left of the mean
		t = fabs(s/sigmaL), alpha = absAlphaL, power = powerL;
	}
	else{//to the right of the mean
		t = fabs(s/sigmaR), alpha = absAlphaR, power = powerR;
	}///end to the right of the mean
	
	if(t > alpha){//in the power law tail to the left of the mean
		a = pow(power/alpha,power)*exp(-0.5*alpha*alpha);
		b = (power/alpha) - alpha;
		return a/pow(b+t,power);
	}
	else{//in the gaussian region to the left of the mean
		return exp(-0.5*t*t);
	}
}///end evaluate()

