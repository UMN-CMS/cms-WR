#ifndef AnalysisResult_h
#define AnalysisResult_h
#include <TLorentzVector.h>
#include <TTree.h>
#include <TChain.h>

/** \class AnalysisResult AnalysisResult.h AnalysisResult.cc
 * This class defines the content of one fit
 * There are also methods to read and to write from a TTree the event
 */
class AnalysisResult
{
public:

// public members to be filled by your program
	Float_t normalization;
	UInt_t nparam;
	UInt_t nmasses;
	Float_t fit_parameters[16], fit_parameter_errors[16];
	Float_t events_in_range[64];
	Float_t fit_integral_in_range[64];

	AnalysisResult(); ///< default contructor (empty)

	void SetBranches(TTree* tree);
	//void SetBranchAddresses(TTree* tree);

};


#endif
