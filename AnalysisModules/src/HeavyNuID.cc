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
  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
  const double scalelo2011A[]  = {0.978303,0.992186,0.988027,0.957773,1.005050,0.958547,0.990894,0.990894} ; 
  const double scalenom2011A[] = {0.983719,1.001520,1.001770,0.973712,1.029180,0.982626,1.107540,1.107540} ; 
  const double scalehi2011A[]  = {0.989135,1.010850,1.015520,0.989651,1.053310,1.006700,1.224180,1.224180} ; 
  const double upedge2011A[]   = {      40,      50,      60,      80,     100,     200,    3500,      -1} ;

  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
  const double scalelo2011B[]  = {0.963595,0.929206,0.938946,1.010270,1.031300,1.016630,0.744268,0.744268} ; 
  const double scalenom2011B[] = {0.968602,0.937838,0.952025,1.024040,1.053860,1.037920,0.843972,0.843972} ; 
  const double scalehi2011B[]  = {0.973609,0.946470,0.965104,1.037820,1.076420,1.059220,0.943676,0.943676} ; 
  const double upedge2011B[]   = {      40,      50,      60,      80,     100,     200,    3500,      -1} ;

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
