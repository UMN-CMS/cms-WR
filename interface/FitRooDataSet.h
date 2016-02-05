#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"

void fitRooDataSet(RooFitResult * ftRslt, RooDataSet * dataSetIn, RooAbsPdf * pdf)
{

	ftRslt = pdf->fitTo(*dataSetIn, RooFit::SumW2Error(kTRUE), RooFit::Save(kTRUE));

}
