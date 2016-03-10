#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"


namespace Fits
{
///consts needed for TTBar or MuonEG fit rescaling
Float_t nEMuDataEvts = 464.;
Float_t intLumi = 2488.245, ttBarXsxn = 56.7, nTTBarEvts = 24521141.f, emuInter = 1.165;
//Float_t ttBarMcEEtoEMuRatioSlope = -0.000028, ttBarMcEEtoEMuRatioIntercept = 0.44, ttBarMcMuMutoEMuRatioSlope = 0.000012, ttBarMcMuMutoEMuRatioIntercept = 0.64;
Float_t ttBarMcEEtoEMuRatioIntercept = 0.42, ttBarMcMuMutoEMuRatioIntercept = 0.65;


///RooFit functions, and prerequisite RooRealVars and RooFormulaVars, used in bkgndFitRooMacro.C and other RooFit macros
RooRealVar massWR("fourObjectMass", "fourObjectMass", 600, 6500);
RooRealVar genEvtWeights("evWeightSign", "evWeightSign", -2, 2);
RooArgSet vars(massWR, genEvtWeights);

RooRealVar expPowerForLinearExp("expPowerForLinearExp", "", -1., -0.00005);

RooRealVar expPower("expPower", "expPower", -1., -0.00005);
RooExponential expPdf("rescaledExpPdf", "", massWR, expPower);
RooAbsPdf * expPdfRooAbsPdf = &expPdf;


RooExponential mumuExpPdf("rescaledExpPdf", "", massWR, expPower);
RooAbsPdf * mumuExpPdfRooAbsPdf = &mumuExpPdf;


RooRealVar mcExpPower("mcExpPower", "", -0.05, -0.00005);
RooExponential mcExpPdf("mcExpPdf", "mcExpPdf", massWR, mcExpPower);
RooRealVar pwrOne("pwrOne", "", 0., 100.);
RooRealVar pwrTwo("pwrTwo", "", 0., 10.);
RooRealVar dOne("dOne", "dOne", 12.0, 14.);
RooRealVar dTwo("dTwo", "dTwo", 4.1, 8.5);
RooRealVar dThree("dThree", "dThree", 0.3, 1.);

RooRealVar systExpPower("systExpPower", "", -0.05, -0.000005);
RooExponential systExpPdf("systExpPdf", "systExpPdf", massWR, systExpPower);

//sum of two exponentials
RooRealVar expPowerForSumTwoExp("expPowerForSumTwoExp", "expPowerForSumTwoExp", -1., -0.00005);
RooRealVar offsetScdExpPow("offsetScdExpPow", "", -0.1, -0.00001);
RooFormulaVar scdExpPow("scdExpPow", "@0+@1", RooArgSet(expPowerForSumTwoExp, offsetScdExpPow));
RooExponential scdExpPdf("scdExpPdf", "", massWR, scdExpPow);
RooRealVar coef("coef", "", 0., 1.);

RooRealVar linearSlope("linearSlope", "", 0., 0.001); ///<slope of linear fit to TTBar LLJJ mass/EMuJJ mass ratio
RooPolynomial linearPdf("linearPdf", "", massWR, RooArgList(linearSlope));
RooRealVar systLinearSlope("systLinearSlope", "", 0., 0.001); ///<slope of systLinear fit to TTBar LLJJ mass/EMuJJ mass ratio
RooPolynomial systLinearPdf("systLinearPdf", "", massWR, RooArgList(systLinearSlope));


//modified exponential
RooRealVar expPowerForModExp("expPowerForModExp", "expPowerForModExp", -1., -0.00005);
RooRealVar modExpPow("modExpPow", "", 0.0, 1.5);

//quad exponential
RooRealVar expPowerForQuadExp("expPowerForQuadExp", "expPowerForQuadExp", -1., -0.00005);
RooRealVar multExpPowFactor("multExpPowFactor", "", 0.0, 0.0001);
RooFormulaVar expPowTwo("expPowTwo", "(@0)*(@1)", RooArgSet(expPowerForQuadExp, multExpPowFactor));


RooProdPdf linearWithExpPdf("linearWithExpPdf", "", expPdf, linearPdf);
RooProdPdf systematicLinearWithExpPdf("systematicLinearWithExpPdf", "", systExpPdf, systLinearPdf);



RooGenericPdf dataFitPowerLawNatLog("powerLawNatLog", "std::pow((1-(@0/13000)),@1)/std::pow((@0/13000),(@2 + (@3)*(TMath::Log(@0/13000)) ))", RooArgSet(massWR, dOne, dTwo, dThree));
RooAbsPdf * pwrLawNatLogRooAbsPdf = &dataFitPowerLawNatLog;


RooGenericPdf dataFitPowerLaw("rescaledPowerLawSystematicPdf", "std::pow((1-(@0/13000)),@1)/std::pow((@0/13000),@2)", RooArgSet(massWR, pwrOne, pwrTwo));
RooAbsPdf * pwrLawRooAbsPdf = &dataFitPowerLaw;

RooProdPdf dataFitLinearExp("linearExp", "", systExpPdf, systLinearPdf);
RooAbsPdf * linearExpRooAbsPdf = &dataFitLinearExp;

RooGenericPdf dataFitModExp("modExp", "TMath::Exp((@1)*(TMath::Power(@0,@2)))", RooArgSet(massWR, expPowerForModExp, modExpPow));
RooAbsPdf * modExpRooAbsPdf = &dataFitModExp;

RooGenericPdf dataFitQuadraticExp("quadraticExp", "TMath::Exp((@1)*(@0) + (@2)*(TMath::Power(@0,2)))", RooArgSet(massWR, expPowerForQuadExp, expPowTwo));
RooAbsPdf * quadExpRooAbsPdf = &dataFitQuadraticExp;

RooAddPdf dataFitSumTwoExp("sumTwoExp", "", expPdf, scdExpPdf, coef);
RooAbsPdf * sumTwoExpRooAbsPdf = &dataFitSumTwoExp;


///USE THESE LINES OF CODE for multiplying a 1st degree RooPolynomial with a fit fxn defined above
//RooPolynomial systLinearPdf("systLinearPdf","", massWR, RooArgList(systLinearSlope));
//RooPolynomial mumuSystLinearPdf("mumuSystLinearPdf","", massWR, RooArgList(mumuSystLinearSlope));
//systLinearSlope.setVal(ttBarMcEEtoEMuRatioSlope/ttBarMcEEtoEMuRatioIntercept);
//mumuSystLinearSlope.setVal(ttBarMcMuMutoEMuRatioSlope/ttBarMcMuMutoEMuRatioIntercept);

//RooProdPdf rescaledEEPowerLawSystematicPdf("rescaledPowerLawSystematicPdf","",dataFitPowerLaw,systLinearPdf);
//RooAbsPdf * rescaledEEPwrLawSystPdf = &rescaledEEPowerLawSystematicPdf;

//RooProdPdf rescaledMuMuPowerLawSystematicPdf("rescaledPowerLawSystematicPdf","",dataFitPowerLaw,mumuSystLinearPdf);
//RooAbsPdf * rescaledMuMuPwrLawSystPdf = &rescaledMuMuPowerLawSystematicPdf;

//RooProdPdf rescaledEEExp("rescaledExpPdf","",expPdf,systLinearPdf);
//RooAbsPdf * rescaledEEExpPdf = &rescaledEEExp;

//RooProdPdf rescaledMuMuExp("rescaledExpPdf","",expPdf,mumuSystLinearPdf);
//RooAbsPdf * rescaledMuMuExpPdf = &rescaledMuMuExp;

//RooProdPdf rescaledEESumTwoExp("rescaledSumTwoExpPdf","", dataFitSumTwoExp,systLinearPdf);
//RooAbsPdf * rescaledEESumTwoExpPdf = &rescaledEESumTwoExp;

//RooProdPdf rescaledEEModExp("rescaledModExpPdf","", dataFitModExp, systLinearPdf);
//RooAbsPdf * rescaledEEModExpPdf = &rescaledEEModExp;

//RooProdPdf rescaledMuMuModExp("rescaledModExpPdf","", dataFitModExp, mumuSystLinearPdf);
//RooAbsPdf * rescaledMuMuModExpPdf = &rescaledMuMuModExp;

//RooProdPdf rescaledEEQuadExp("rescaledQuadExpPdf","", dataFitQuadraticExp, systLinearPdf);
//RooAbsPdf * rescaledEEQuadExpPdf = &rescaledEEQuadExp;

//RooProdPdf rescaledMuMuQuadExp("rescaledQuadExpPdf","", dataFitQuadraticExp, mumuSystLinearPdf);
//RooAbsPdf * rescaledMuMuQuadExpPdf = &rescaledMuMuQuadExp;

}
