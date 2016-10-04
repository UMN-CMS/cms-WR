#include "../interface/SelectorHist.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <algorithm>

SelectorHist sel::hists;

hist_p SelectorHist::operator()(std::string name, int nbins, float low, float hi)
{
	if(histograms.find(name) == histograms.end()) {
		histograms[name] = std::make_shared<TH1F>(name.c_str(), name.c_str(), nbins, low, hi);
		histogram_names.push_back(name);
	}
	return histograms[name];
}

void SelectorHist::operator()(std::string name, int nbins, float low, float hi, float value, bool pass)
{
	if(nminusones.find(name) == nminusones.end()) {
		nminusones[name] = std::make_shared<TH1F>(name.c_str(), name.c_str(), nbins, low, hi);
		pass_all[name] = std::make_shared<TH1F>((name + "_pass").c_str(), (name + "_pass").c_str(), nbins, low, hi);
	}

	cut c;
	c.name = name;
	c.value = value;
	c.pass = pass;

	cuts.push_back(c);
}


void SelectorHist::FillNMinusOnes()
{
	unsigned int npass = 0;
	auto nminus_ic = cuts.end();
	for (auto ic = cuts.begin(); ic != cuts.end(); ++ic) {
		npass += int(ic->pass);
		if(!ic->pass) nminus_ic = ic;
	}

	if(npass == cuts.size()) {
		for(auto c : cuts) {
			nminusones[c.name]->Fill(c.value);
			pass_all[c.name]->Fill(c.value);
		}
	} else if(npass == cuts.size() - 1) {
		nminusones[nminus_ic->name]->Fill(nminus_ic->value);
	}

	cuts.clear();
}
void SelectorHist::PrintEntries(std::string folder, std::string tag)
{
	std::ofstream o(folder + "/" + "selhists_" + tag + ".txt");

	for(auto name : histogram_names) {
		o << name << "\t" << histograms[name]->GetEntries() << std::endl;
	}

	for(auto i : nminusones) {
		std::string name = i.first;
		o << "nminusones" << '\t' << name << "\t" << i.second->GetEntries() << '\t' << pass_all[name]->GetEntries() << std::endl;
	}
}
void SelectorHist::Draw(std::string folder, std::string tag)
{
	gStyle->SetOptStat(0);
	TCanvas c("c", "c", 600, 600);
	TLatex tex;
	tex.SetNDC(true);
	tex.SetTextSize(0.03);
	for(auto name : histogram_names) {
		if( name.find("cut") != std::string::npos) continue;
		histograms[name]->Draw();
		if (histograms.find(name + "_cut") != histograms.end()) {
			histograms[name + "_cut"]->SetLineColor(kRed);
			histograms[name + "_cut"]->Draw("same");
			tex.DrawLatex(.11, .85, ("before: " + std::to_string((int)histograms[name]->GetEntries())).c_str());
			tex.DrawLatex(.11, .8, ("after: " + std::to_string((int)histograms[name + "_cut"]->GetEntries())).c_str());
		}
		c.SaveAs((folder + "/selhists_" + tag + "_" + name + ".png").c_str());
	}
	for(auto i : nminusones) {
		i.second->Draw();
		pass_all[i.first]->SetLineColor(kRed);
		pass_all[i.first]->Draw("same");
		tex.DrawLatex(.11, .85, ("before: " + std::to_string((int)i.second->GetEntries())).c_str());
		tex.DrawLatex(.11, .8, ("after: " + std::to_string((int)pass_all[i.first]->GetEntries())).c_str());
		c.SaveAs((folder + "/nminusones_" + tag + "_" + i.first + ".png").c_str());
	}
}

void SelectorHist::Clear()
{
	histograms.clear();
	histogram_names.clear();
}
