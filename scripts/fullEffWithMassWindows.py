import ROOT
import numpy as np
from ExoAnalysis.cmsWR.PlotUtils import customROOTstyle
import ExoAnalysis.cmsWR.combineTools as combineTools

customROOTstyle()
ROOT.gROOT.SetBatch(True)   #disable visual rendering of plots

#declare numpy arrays which hold WR mass points, and the efficiency of WR events at each mass point to pass the entire online and offline selection including mass window cuts
binsize = 50
nGeneratedEvents = 50000.
mass = []   #holds WR mass points as strings without decimal point
floatMass = np.ndarray(25)   #holds WR mass points as floats, with decimal point
eleLowBound = np.ndarray(25)
eleUpBound = np.ndarray(25)
ele_results = np.ndarray(25)
muLowBound = np.ndarray(25)
muUpBound = np.ndarray(25)
mu_results = np.ndarray(25)
linei=0

#fill the arrays named eleLowBound and eleUpBound with the lower and upper mass window bounds
#fill the array named mass with the MWR values
#read this info from mass_cuts.txt in configs dir
with open("configs/mass_cuts.txt","r") as c:
	for line in c:
		if line[0] == '#': continue
		line = map(str,line.split())
		#now line is split into a 4 element array, and each element is a string NOT a float or int
		if line[0] == 'EE':
			if line[1] == '0': continue
			mass.append(line[1])
			floatMass[linei] = float(line[1])
			eleLowBound[linei] = float(line[2])
			eleUpBound[linei] = float(line[3])
			linei += 1
		if line[0] == 'MuMu':
			if line[1] == '0':
				linei = 0
				continue
			muLowBound[linei] = float(line[2])
			muUpBound[linei] = float(line[3])
			linei += 1

#end loop over lines in mass_cuts.txt

#fill ele_results and mu_results using the array named mass (converted into an int array) to read the appropriate root files
#since mass is a 1D array, m is a pair
#ele_results and mu_results will be filled with the efficiency of all online + offline + mass window cuts with respect to the
#total number of events generated
for m in enumerate(mass):
	tempMuNumerHist = combineTools.getOtherMassHisto("dilepton_mass>200. && WR_mass>"+str(muLowBound[m[0]])+" && WR_mass<"+str(muUpBound[m[0]]),"Tree_Iter0", ("WeightedWithWindowWR_mumujj_" + str(int(m[1])), "analysisCppOutputRootFiles/selected_tree_WRtoMuMuJJ_%i_%i_signal_%s.root" % (int(m[1]), (int(m[1])/2), "mumuMuMu") ), binsize=binsize)
	mu_results[m[0]] = float(tempMuNumerHist.Integral(1,tempMuNumerHist.GetNbinsX())/nGeneratedEvents)
	
	tempEleNumerHist = combineTools.getOtherMassHisto("dilepton_mass>200. && WR_mass>"+str(eleLowBound[m[0]])+" && WR_mass<"+str(eleUpBound[m[0]]),"Tree_Iter0", ("WeightedWithWindowWR_eejj_" + str(int(m[1])), "analysisCppOutputRootFiles/selected_tree_WRtoEEJJ_%i_%i_signal_%s.root" % (int(m[1]), (int(m[1])/2), "eeEE") ), binsize=binsize)
	ele_results[m[0]] = float(tempEleNumerHist.Integral(1,tempEleNumerHist.GetNbinsX())/nGeneratedEvents)


##now that the mass window cut efficiencies in both channels have been obtained
##make two graphs (one for each lepton channel, to be consistent with the AN) showing the mass window cut efficiency
##versus WR mass of the centrally produced WR signal samples
h = ROOT.TH1F("h","",1,700, 6200)


h.SetMinimum(0)
h.SetMaximum(1)

##electrons
c = ROOT.TCanvas("c","c",800,800)
h.Draw()
h.SetXTitle("M_{EEJJ} [GeV]")
h.SetYTitle("Selection+Mass Window Eff*Accept")

leg = ROOT.TLegend(.3, .8, .85, .9)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs = []
graphs.append( ROOT.TGraph(len(floatMass), floatMass, ele_results))
g = graphs[-1]
g.SetLineColor(ROOT.kBlack)
g.SetLineWidth(3)
g.Draw("L")
leg.AddEntry(g, "Global", "L")

leg.Draw("same")
c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.png")
c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.pdf")
c.SaveAs("plots/eleChnlWrSignalEffDiagonalMNuWithMassWindowCuts.C")


h.Draw()
h.SetXTitle("M_{MuMuJJ} [GeV]")

##muons
leg = ROOT.TLegend(.3, .8, .85, .9)
leg.SetTextFont(42)
leg.SetTextSize(0.032)
leg.SetBorderSize(0)

graphs.append( ROOT.TGraph(len(floatMass), floatMass, mu_results))
g = graphs[-1]
g.SetLineColor(ROOT.kBlack)
g.SetLineWidth(3)
g.Draw("L")
leg.AddEntry(g, "Global", "L")


leg.Draw("same")
c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.png")
c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.pdf")
c.SaveAs("plots/muChnlWrSignalEffDiagonalMNuWithMassWindowCuts.C")
