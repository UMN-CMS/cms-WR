// #include "ExoAnalysis/cmsWR/interface/AnalysisResult.h"
#include "../interface/AnalysisResult.h"

AnalysisResult::AnalysisResult() : nparam(0), nmasses(0)
{
}


void AnalysisResult::SetBranches(TTree* tree)
{
	tree->Branch("Normalization", &normalization);
	tree->Branch("nparam", &nparam);
	tree->Branch("FitParameters", &fit_parameters, "FitParameters[nparam]/F");
	tree->Branch("FitParameterErrors", &fit_parameter_errors, "FitParameterErrors[nparam]/F");
	tree->Branch("nmasses", &nmasses);
	tree->Branch("NEventsInRange", &events_in_range, "NEventsInRange[nmasses]/F");
	tree->Branch("UnweightedNEventsInRange", &unweighted_events_in_range, "UnweightedNEventsInRange[nmasses]/i");
	tree->Branch("ErrorEventsInRange", &error_in_range, "ErrorEventsInRange[nmasses]/F");
	tree->Branch("FitIntegralInRange", &fit_integral_in_range, "FitIntegralInRange[nmasses]/F");
}

//void AnalysisResult::SetBranchAddress(TTree* tree)
//{
//	tree->SetBranchAddress("Normalization", &normalization);
//	tree->SetBranchAddress("nparam", &nparam);
//	tree->SetBranchAddress("FitParameters", &fit_parameters);
//	tree->SetBranchAddress("FitParameterErrors", &fit_parameter_errors);
//	tree->SetBranchAddress("nmasses", &nmasses);
//	tree->SetBranchAddress("NEventsInRange", &events_in_range);
//	tree->SetBranchAddress("FitIntegralInRange", &fit_integral_in_range);
//}
