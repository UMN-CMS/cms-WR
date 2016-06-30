#include "../interface/SelectorHist.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
#include <fstream>

SelectorHist sel::hists;

hist_p SelectorHist::operator()(std::string name, int nbins, float low, float hi)
{
	if(histograms.find(name) == histograms.end()) {
		histograms[name] = std::make_shared<TH1F>(name.c_str(), name.c_str(), nbins, low, hi);
		histogram_names.push_back(name);
	}
	return histograms[name];
}

void SelectorHist::PrintEntries(std::string folder, std::string tag)
{
	std::ofstream o(folder + "/" + "selhists_" + tag + ".txt");

	for(auto name : histogram_names) {
		o << name << "\t" << histograms[name]->GetEntries() << std::endl;
	}
}
void SelectorHist::Draw(std::string folder, std::string tag)
{
	//ourmen
	gStyle->SetOptStat(110010);
	TCanvas c("c", "c", 600, 600);
	for(auto name : histogram_names) {
		if( name.find("cut") != std::string::npos) continue;
		histograms[name]->Draw();
		if (histograms.find(name + "_cut") != histograms.end()) {
			histograms[name + "_cut"]->SetLineColor(kRed);
			histograms[name + "_cut"]->Draw("same");
		}
		c.SaveAs((folder + "/selhists_" + tag + "_" + name + ".png").c_str());
	}
}

void SelectorHist::Clear()
{
	histograms.clear();
	histogram_names.clear();
}
