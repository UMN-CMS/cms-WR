
//#ifndef __CINT__
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
//#endif

/**/
class RooDoubleCB : public RooAbsPdf {
	public:

		RooDoubleCB(const char *name, const char *title,
				RooAbsReal& _x, RooAbsReal& _mean, 
				RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
				RooAbsReal& _alphaL, RooAbsReal& _alphaR,
				RooAbsReal& _powerL, RooAbsReal& _powerR
				);
		RooDoubleCB(const RooDoubleCB& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooDoubleCB(*this,newname); }
		inline virtual ~RooDoubleCB() { }

		//  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
		//  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

		//  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
		//  void generateEvent(Int_t code);

	protected:
		
		RooRealProxy x;
		RooRealProxy mean;
		RooRealProxy sigmaL;
		RooRealProxy sigmaR;
		RooRealProxy alphaL;
		RooRealProxy alphaR;
		RooRealProxy powerL;
		RooRealProxy powerR;

		Double_t evaluate() const;


};
/**/
