import ExoAnalysis.cmsWR.combineTools as combineTools

import argparse


def getRange(hist, perc):
	integral = hist.Integral()
	min_range = (hist.GetNbinsX(), 1, hist.GetNbinsX(), 1)
	for i in range(1, hist.GetNbinsX()):
		good = False
		for j in range(i + 1, hist.GetNbinsX()):
			range_hist = hist.Integral(i, j)
			if perc <= range_hist/integral:
				good = True
				break
		if good:
			min_range = min( min_range, (j - i, i, j, range_hist/integral) )
	size, low, hi, p = min_range
	return (hist.GetBinLowEdge(low), 
			hist.GetBinLowEdge(hi+1), p)

parser = argparse.ArgumentParser(description='Make datacards')
parser.add_argument('-d', '--dir', dest='basedir',
		default="./",
		help='base dir for analysis tree files')
parser.add_argument('-t', '--tag', dest='tag',
		default="",
		help='tag name for analysis tree files')

args = parser.parse_args()

minitrees = combineTools.miniTreeInterface(
			base=args.basedir,
			tag =args.tag,
			)

for channel in ["ee", "mumu"]:
	for mass in combineTools.mass_cut:
		try:
			h = minitrees.getMassHisto(mass, channel, "signal")
			for p in [.7, .8, .9]:
				low, hi, a_p = getRange(h, p)
				print channel, p, '%04d' % mass, low, hi, a_p
		except:
			pass


