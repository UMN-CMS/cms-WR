#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuID::HeavyNuID(const edm::ParameterSet & iConfig):
  idEra_ ( iConfig.getParameter< int > ( "eraForId" ) )
{
    if ( idEra_ == 20110 )      std::cout << "Applying run 2011A muon ID corrections: " << std::endl ; 
    else if ( idEra_ == 20111 ) std::cout << "Applying run 2011B muon ID corrections: " << std::endl ; 
    else {
        std::cout << "WARNING: Unknown correction.  Applying 2011A" << std::endl ;
        idEra_ = 20110 ;
    }
}                                      // HeavyNuID::HeavyNuID

//======================================================================

double
HeavyNuID::weightForMC(double pt,int signOfError2apply)
{
  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, October 27, 2011
  const double scalelo2011A[]  = {0.994,0.999,0.999} ; 
  const double scalenom2011A[] = {0.996,1.000,1.000} ; 
  const double scalehi2011A[]  = {0.998,1.001,1.001} ; 
  const double upedge2011A[]   = {   40, 3500,   -1} ;

  // Not yet determined...defaults taken as 2011A for safety
  const double scalelo2011B[]  = {0.994,0.999,0.999} ; 
  const double scalenom2011B[] = {0.996,1.000,1.000} ; 
  const double scalehi2011B[]  = {0.998,1.001,1.001} ; 
  const double upedge2011B[]   = {   40, 3500,   -1} ;

  const double *scale  = ( (idEra_ == 20110) ? scalenom2011A : scalenom2011B );
  const double *upedge = ( (idEra_ == 20110) ? upedge2011A : upedge2011B );

  if ( signOfError2apply ) {
    if ( idEra_ == 20110 ) scale = (signOfError2apply > 0) ? scalehi2011A : scalelo2011A ;
    if ( idEra_ == 20111 ) scale = (signOfError2apply > 0) ? scalehi2011B : scalelo2011B ;
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
