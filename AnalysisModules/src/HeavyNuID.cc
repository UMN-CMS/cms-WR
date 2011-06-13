#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuID::HeavyNuID(const edm::ParameterSet & iConfig):
  idEra_ ( iConfig.getParameter< int > ( "eraForId" ) )
{
  std::cout << "Correction era for ID: " << idEra_ << std::endl ; 
}                                      // HeavyNuID::HeavyNuID

//======================================================================

double
HeavyNuID::weightForMC(double pt,int signOfError2apply)
{
  // determined from 2011 42x data and Spring11 MC
  const double scalelo2011[]  = {1.00088,0.983438,0.988129,0.986851,0.979319,0.953469,0.953469} ; 
  const double scalenom2011[] = {1.00329,0.984218,0.989802,0.99039,0.985222,0.966939,0.966939} ; 
  const double scalehi2011[]  = {1.00571,0.984998,0.991475,0.993928,0.991126,0.980408,0.980408} ; 
  const double upedge2011[]   = {   40,   50,   60,   70, 100, 3500,   -1};

  const double scalelo2010[]  = {0.990575,0.975468,0.984249,0.985711,0.982309,0.982309} ; 
  const double scalenom2010[] = {0.999596,0.977869,0.98594,0.988714,0.990021,0.990021} ; 
  const double scalehi2010[]  = {1.00862,0.98027,0.987632,0.991717,0.997733,0.997733} ; 
  const double upedge2010[]   = {   30,   40,   50,   70, 3500,   -1};

  const double *scale  = ( (idEra_ == 2010) ? scalenom2010 : scalenom2011 );
  const double *upedge = ( (idEra_ == 2010) ? upedge2010 : upedge2011 );

  if ( signOfError2apply ) {
    if ( idEra_ == 2010 ) scale = (signOfError2apply > 0) ? scalehi2010 : scalelo2010 ;
    if ( idEra_ == 2011 ) scale = (signOfError2apply > 0) ? scalehi2011 : scalelo2011 ;
  }

  int i;
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double factor=scale[i];
    
  return factor ; 
}

//======================================================================

void HeavyNuID::endJob()
{
}
