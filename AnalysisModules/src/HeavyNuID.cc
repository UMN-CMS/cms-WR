#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuID::HeavyNuID(const edm::ParameterSet & iConfig):
  idEra_ ( iConfig.getParameter< int > ( "eraForId" ) )
{
    if ( idEra_ == 20111 )      std::cout << "Applying run 2011A muon ID corrections: " << std::endl ; 
    else if ( idEra_ == 20112 ) std::cout << "Applying run 2011B muon ID corrections: " << std::endl ; 
    else {
        std::cout << "WARNING: Unknown correction.  Applying 2011A" << std::endl ;
        idEra_ = 20111 ;
    }
}                                      // HeavyNuID::HeavyNuID

//======================================================================

double
HeavyNuID::weightForMC(double pt,int signOfError2apply)
{
  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
  const double scalelo2011A[]  = {1.00139,0.99304,1.00710,0.968835,0.938102,0.945127,0.774347,0.774347} ; 
  const double scalenom2011A[] = {1.00713,1.00382,1.02449,0.988440,0.962890,0.968289,0.899621,0.899621} ; 
  const double scalehi2011A[]  = {1.01286,1.01459,1.04188,1.008040,0.987678,0.991451,1.024900,1.024900} ; 
  const double upedge2011A[]   = {      40,      50,      60,      80,     100,     200,    3500,      -1} ;

  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
  const double scalelo2011B[]  = {0.995869,0.982769,0.998694,0.988053,0.987864,0.996661,0.894987,0.894987} ; 
  const double scalenom2011B[] = {1.000590,0.991830,1.015460,1.005710,1.011420,1.018780,1.000000,1.000000} ; 
  const double scalehi2011B[]  = {1.005300,1.000890,1.032230,1.023370,1.034980,1.040900,1.105010,1.105010} ; 
  const double upedge2011B[]   = {      40,      50,      60,      80,     100,     200,    3500,      -1} ;

  const double *scale  = ( (idEra_ == 20111) ? scalenom2011A : scalenom2011B );
  const double *upedge = ( (idEra_ == 20111) ? upedge2011A : upedge2011B );

  if ( signOfError2apply ) {
    if ( idEra_ == 20111 ) scale = (signOfError2apply > 0) ? scalehi2011A : scalelo2011A ;
    if ( idEra_ == 20112 ) scale = (signOfError2apply > 0) ? scalehi2011B : scalelo2011B ;
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
