#include "../interface/RooCruijff.h"

///constructor from user defined args
///two Gaussians with the same mean and different sigmas, which far from the mean
///transform into two exponential decays with different decay powers alphaL and alphaR
RooCruijff::RooCruijff(const char * name, const char * title,
		RooAbsReal& _x, RooAbsReal& _mean,
		RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		RooAbsReal& _alphaL, RooAbsReal& _alphaR
		):
	RooAbsPdf(name,title),
	x("x","Observable", this, _x),
	mean("mean","Mean", this, _mean),
	sigmaL("sigmaL","sigmaL", this, _sigmaL),
	sigmaR("sigmaR","sigmaR", this, _sigmaR),
	alphaL("alphaL","alphaL", this, _alphaL),
	alphaR("alphaR","alphaR", this, _alphaR){

	}

///constructor using existing RooCruijff object - needed for clone() method
RooCruijff::RooCruijff(const RooCruijff& other, const char* name):
	RooAbsPdf(other,name), 
	x("x", this, other.x),
	mean("mean", this, other.mean),
	sigmaL("sigmaL", this, other.sigmaL),
	sigmaR("sigmaR", this, other.sigmaR),
	alphaL("alphaL", this, other.alphaL),
	alphaR("alphaR", this, other.alphaR){

	}

///fxn which is called automatically by RooFit methods like fitTo
Double_t RooCruijff::evaluate() const{

	Double_t sigma, alpha;
	if(x < mean) {
		sigma=sigmaL;
		alpha=alphaL;
	} else {
		sigma = sigmaR;
		alpha = alphaR;
	}

	Double_t delta_x = x - mean;
	Double_t delta_x2 = delta_x * delta_x;

	return exp(- delta_x2 / (2 * sigma * sigma + alpha * delta_x2));
}

