#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TKey.h>
#include <TTree.h>
#include <TChain.h>

#include <string>
#include <iostream>

class producePileupWeight : public edm::EDProducer
{

public:
	producePileupWeight(const edm::ParameterSet&);

	virtual void produce(edm::Event&, const edm::EventSetup&);

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
	std::string outputCollName;
	std::string PileupMCFilename;
	std::string PileupDataFilename;
	//Tokens
	edm::EDGetToken pileUpInfoToken_;
	// typedef
	typedef std::map<float, double> PUweights_t;
	PUweights_t pu_weights;

	static const int MAX_PU_REWEIGHT = 59;

};

producePileupWeight::producePileupWeight(const edm::ParameterSet& cfg):
	outputCollName(cfg.getParameter<std::string>("outputCollectionName")),
	PileupMCFilename(cfg.getParameter<std::string>("PileupMCFilename")),
	PileupDataFilename(cfg.getParameter<std::string>("PileupDataFilename")),
	pileUpInfoToken_(mayConsume<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo")))
{
	produces <float>(outputCollName);

	TFile data(PileupDataFilename.c_str());
	TFile mc(PileupMCFilename.c_str());

	if(!data.IsOpen() || !mc.IsOpen()) {
		std::cerr << "[ERROR] data or mc PU file not opened" << std::endl;
		return;
	}

	TH1D* puMC_hist = (TH1D*) mc.Get("pileup");
	if(puMC_hist == NULL) {
		std::cerr << "[ERROR] mc pileup histogram is NULL" << std::endl;
		return;
	}

	TH1D *puData_hist = (TH1D *) data.Get("pileup");
	if(puData_hist == NULL) {
		std::cerr << "[ERROR] data pileup histogram is NULL" << std::endl;
		exit(1);
	}

	// they should have bin width = 1
	if(puMC_hist->GetBinWidth(2) != 1) {
		std::cerr << "[ERROR] Bin width for mc pileup distribution different from 1" << std::endl;
		return;
	}

	if(puData_hist->GetBinWidth(2) != 1) {
		std::cerr << "[ERROR] Bin width for data pileup distribution different from 1" << std::endl;
		return;
	}

	// they should have the same binning otherwise PU weights until the max between them
	int nBins = MAX_PU_REWEIGHT;
	if(puData_hist->GetNbinsX() != puMC_hist->GetNbinsX()) {
		std::cerr << "[WARNING] Different binning between data and mc pileup distributions:" << puData_hist->GetNbinsX() << "\t" << puMC_hist->GetNbinsX() << std::endl;
		nBins = std::min(puData_hist->GetNbinsX(), puMC_hist->GetNbinsX());
		//std::cerr << "Pileup reweighted until nPU max = " << nBins;
	}

	nBins = std::min(MAX_PU_REWEIGHT, nBins);

	// Normalize Histograms
	float puMC_int = puMC_hist->Integral(0, nBins);
	float puData_int = puData_hist->Integral(0, nBins);
	puMC_hist->Scale(1. / puMC_int);
	puData_hist->Scale(1. / puData_int);

	double puMCweight_int = 0;
	for (int i = 0; i < nBins; i++) {
		double binContent = puMC_hist->GetBinContent(i + 1);
		if(binContent == 0 && puData_hist->GetBinContent(i + 1) != 0) {
			if(puData_hist->GetBinContent(i + 1) < 1e-4 || i<4) {
				std::cerr << "[WARNING] mc bin empty while data not: iBin = " << i + 1 << std::endl;
				std::cerr << "          data bin = " << puData_hist->GetBinContent(i + 1) << std::endl;
			} else {
				std::cerr << "[ERROR] mc bin empty while data not: iBin = " << i + 1 << std::endl;
				std::cerr << "        data bin = " << puData_hist->GetBinContent(i + 1) << std::endl;
				exit(1);
			}
		}
		double weight = binContent > 0 ? puData_hist->GetBinContent(i + 1) / binContent : 0;
		// the MC normalization should not change, so calculate the
		// integral of the MC and rescale the weights by that
		puMCweight_int += weight * binContent;
		pu_weights[i] = weight;
	}

	for (int i = 0; i < nBins; i++) {
		pu_weights[i] /= puMCweight_int;
	}
}

void producePileupWeight::produce(edm::Event& event, const edm::EventSetup&)
{
	float nPU = 0;
	float PU_weight = 1.0;
	if(!event.isRealData()) {
		edm::Handle<edm::View<PileupSummaryInfo> > PU_Info;
		event.getByToken(pileUpInfoToken_, PU_Info);
		for(auto p : *PU_Info) {
			int BX = p.getBunchCrossing();
			if(BX == 0)
				nPU = p.getTrueNumInteractions();
		}
		PU_weight = pu_weights[nPU];
	}
	std::auto_ptr<float> PU_weight_ptr(new float(PU_weight));
	event.put(PU_weight_ptr, outputCollName);
}

void
producePileupWeight::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
	edm::ParameterSetDescription desc;
	desc.add<std::string>("outputCollectionName", "PileupWeights");
	desc.add<std::string>("PileupMCFilename");
	desc.add<std::string>("PileupDataFilename");
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(producePileupWeight);
