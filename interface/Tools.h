#ifndef ExoAnalysis_cmsWR_Tools
#define ExoAnalysis_cmsWR_Tools

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

std::vector<const reco::Candidate*> ancestry_muon(const pat::Muon& muon)
{
	std::vector<const reco::Candidate*> result;
	const reco::Candidate* tmp = muon.genParticle();
	while(tmp->mother() != 0) {
		result.push_back(tmp);
		tmp = tmp->mother();
	}
	return result;
}
std::vector<const reco::Candidate*> ancestry_electron(const pat::Electron& electron)
{
	std::vector<const reco::Candidate*> result;
	const reco::Candidate* tmp = electron.genParticle();
	while(tmp->mother() != 0) {
		result.push_back(tmp);
		tmp = tmp->mother();
	}
	return result;
}


#endif
