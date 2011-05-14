#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuID::HeavyNuID()
{
}                                      // HeavyNuID::HeavyNuID

//======================================================================

double
HeavyNuID::weightForMC(double pt,int signOfError2apply)
{
  // determined from 2011 data...corrections for loose+iso efficiency only
  const double scalelo[]  = {0.9935,0.9919,0.9901,0.9808,0.9862,0.9649,0.8927,0.8927};
  const double scalenom[] = {0.9951,0.9930,0.9925,0.9857,0.9957,0.9756,0.9492,0.9492};
  const double scalehi[]  = {0.9967,0.9941,0.9949,0.9906,1.0000,0.9863,1.0000,1.0000};
  const double upedge[]= {   40,   50,   60,   70, 80, 150, 3500,   -1};

  const double *scale = scalenom;
  if ( signOfError2apply )
    scale = (signOfError2apply > 0) ? scalehi : scalelo;

  int i;
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double factor=scale[i];
    
  return factor ; 
}

//======================================================================

void HeavyNuID::endJob()
{
}
