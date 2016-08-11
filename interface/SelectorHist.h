
#ifndef selectorhist_h
#define selectorhist_h

#include "TH1F.h"
#include <map>
#include <memory>
#include <string>

typedef std::shared_ptr<TH1F> hist_p;
class SelectorHist
{
public:
	SelectorHist() {}
	~SelectorHist()
	{
		histograms.clear();
	}

	struct cut {
		std::string name;
		float value;
		bool pass;
	};

	hist_p operator()(std::string name, int nbins, float low, float hi);
	void operator()(std::string name, int nbins, float low, float hi, float value, bool pass);
	void FillNMinusOnes();
	void PrintEntries(std::string folder, std::string tag);
	void Draw(std::string folder, std::string tag);
	void Clear();

private:
	std::map<std::string, hist_p > histograms;
	std::vector<std::string> histogram_names;
	std::map<std::string, hist_p > nminusones;
	std::map<std::string, hist_p > pass_all;
	std::vector<cut> cuts;

};

namespace sel
{
extern SelectorHist hists;
}

#endif
