
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

	hist_p operator()(std::string name, int nbins, float low, float hi);
	void PrintEntries(std::string folder, std::string tag);
	void Draw(std::string folder, std::string tag);
	void Clear();

private:
	std::map<std::string, hist_p > histograms;
	std::vector<std::string> histogram_names;

};

namespace sel
{
extern SelectorHist hists;
}

#endif
