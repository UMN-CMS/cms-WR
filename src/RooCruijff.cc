
#include "../interface/RooCruijff.hh"

//ClassImp(RooCruijff)

  RooCruijff::RooCruijff(const char *name, const char *title,
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

RooCruijff::RooCruijff(const RooCruijff& other, const char *name):
  RooAbsPdf(other,name), 
    x("x", this, other.x),
    mean("mean", this, other.mean),
    sigmaL("sigmaL", this, other.sigmaL),
    sigmaR("sigmaR", this, other.sigmaR),
    alphaL("alphaL", this, other.alphaL),
    alphaR("alphaR", this, other.alphaR){

}


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


