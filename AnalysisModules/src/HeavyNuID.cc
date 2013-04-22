#include "HeavyNu/AnalysisModules/src/HeavyNuID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//======================================================================

HeavyNuID::HeavyNuID(const edm::ParameterSet & iConfig):
  idEra_ ( iConfig.getParameter< int > ( "eraForId" ) )
{
  if ( idEra_ != 0 ) idEra_ = idEra_ / abs(idEra_) ;    
  std::cout << "ID definition for high pT muons: " ;
  if ( idEra_ == 0 ) std::cout << "NONE" << std::endl ; 
  else               std::cout << (idEra_ == -1 ? "LOOSE" : "TIGHT" ) << std::endl ;
}                                      // HeavyNuID::HeavyNuID

//======================================================================

double
HeavyNuID::weightForMC(double pt,int signOfError2apply)
{
  // determined from 2011A 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
//   const double scalelo2011A[]  = {1.00139,0.99304,1.00710,0.968835,0.938102,0.945127,0.774347,0.774347} ; 
//   const double scalenom2011A[] = {1.00713,1.00382,1.02449,0.988440,0.962890,0.968289,0.899621,0.899621} ; 
//   const double scalehi2011A[]  = {1.01286,1.01459,1.04188,1.008040,0.987678,0.991451,1.024900,1.024900} ; 

  // determined from 2011B 42x data and Summer11 sherpa Z+jets MC, December 11, 2011
//   const double scalelo2011B[]  = {0.995869,0.982769,0.998694,0.988053,0.987864,0.996661,0.894987,0.894987} ; 
//   const double scalenom2011B[] = {1.000590,0.991830,1.015460,1.005710,1.011420,1.018780,1.000000,1.000000} ; 
//   const double scalehi2011B[]  = {1.005300,1.000890,1.032230,1.023370,1.034980,1.040900,1.105010,1.105010} ; 

  // No corrections for 2012, assign 3% uncertainty across the board
  const double scalelo2012[]  = {0.970000,0.970000,0.970000,0.970000,0.970000,0.970000,0.970000} ; 
  const double scalenom2012[] = {1.000000,1.000000,1.000000,1.000000,1.000000,1.000000,1.000000} ; 
  const double scalehi2012[]  = {1.030000,1.030000,1.030000,1.030000,1.030000,1.030000,1.030000} ; 

//  const double upedge2011[]   = {      40,      50,      60,      80,     100,     200,    3500,      -1} ;
  const double upedge2012[]   = {      50,      60,      80,     100,     200,    4000,      -1} ;

  const double *scale  = scalenom2012 ; 
//   if (idEra_ == 20111) scale = scalenom2011A ; 
//   if (idEra_ == 20112) scale = scalenom2011B ; 
  // const double *upedge = ( (idEra_ == 20121) ? upedge2012 : upedge2011 );
  const double *upedge = upedge2012;

  if ( signOfError2apply ) {
//     if ( idEra_ == 20111 ) scale = (signOfError2apply > 0) ? scalehi2011A : scalelo2011A ;
//     if ( idEra_ == 20112 ) scale = (signOfError2apply > 0) ? scalehi2011B : scalelo2011B ;
    scale = (signOfError2apply > 0) ? scalehi2012 : scalelo2012 ;
  }

  int i;
  for (i=0; upedge[i]>0 && upedge[i]<pt; i++);
  double factor=scale[i];
    
  return factor ; 
}

double HeavyNuID::weightForMCbyEta(double eta,int signOfError2apply)
{
  // Latest full 19 ifb 45+ GeV ID efficiencies
  // https://twiki.cern.ch/twiki/pub/CMS/TWikiEXO-MUO/muon_eff_NewHighPt_Run_ABCD.pdf
  // systematics here - may change for 100 GeV+ muons
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonTagAndProbe
  //const double scalelo2012[]  = {0.970000,0.970000,0.970000,0.970000}; 
  const double scalenom2012[] = {0.989900,0.988400,0.995800,0.991000}; 
  //const double scalehi2012[]  = {1.030000,1.030000,1.030000,1.030000}; 

  const double upedge2012[]   = {     0.9,     1.2,     2.1,     2.4, -1.0};

  const double systematic = 0.005;
  double scale = 1.0;

  if ( signOfError2apply ) {
    scale = (signOfError2apply > 0)?(1 + systematic):(1 - systematic);
  }

  int i;
  for (i=0; upedge2012[i]>0 && upedge2012[i]<eta; i++);
  double factor=scalenom2012[i]*scale;
    
  return factor ; 
}

double
HeavyNuID::weightElectronsForMC(double eta, int signOfError2apply) {

  // EB/EE efficiency scale factors
  const double ebScale    = 1.0000; // Bias of probe  
  const double ebScaleErr = 0.0200; // Deviation from unity
  // const double ebScaleErr = 0.004742 ; // Stat error on scale factor

  const double eeScale    = 1.0000; 
  const double eeScaleErr = 0.0400; // Deviation from unity
  // const double eeScaleErr = 0.009480 ; // Stat error on scale factor

  double factor  = 1.0 ;
  int systUpDown = 0.0 ;
  if ( signOfError2apply) systUpDown = signOfError2apply / abs(signOfError2apply) ; 
  bool isEB = ( fabs(eta) < 1.442 ) ; 
  bool isEE = ( fabs(eta) < 2.5 && fabs(eta) > 1.56 ) ; 
  if ( isEB ) factor = ebScale + (ebScaleErr * double(systUpDown)) ; 
  if ( isEE ) factor = eeScale + (eeScaleErr * double(systUpDown)) ; 

  return factor ; 
}

//======================================================================

void HeavyNuID::endJob()
{
}
