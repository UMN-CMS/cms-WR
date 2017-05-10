import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np
import math

#run this script from one directory above scripts/, and pipe the output to a .tex file using the '>' operator

#change fileDir to point to the directory from which expected bkgnd and signal, and observed files should be read
#then everything else will run correctly
treeName="syst_tree"
fileDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysAllSystAprilTwentyThree/"

expTopMuMuFile = "selected_tree_data_flavoursidebandEMuMuMu.root"
expDyMuMuFile = "selected_tree_DYAMC_signal_mumuMuMu.root"
obsMuMuFile = "selected_tree_data_signal_mumuMuMu.root"
expTopEEFile = "selected_tree_data_flavoursidebandEMuEE.root"
expDyEEFile = "selected_tree_DYAMC_signal_eeEE.root"
obsEEFile = "selected_tree_data_signal_eeEE.root"


##now make a LaTex table with the following format (copied twice, once for each lepton channel):
## WR mass hypothesis | exp WR evts +/- unc | exp DY evts +/- unc | exp top evts +/- unc | exp bkgnd +/- unc | data |
tex= "\\centering\n\\begin{{tabular}}{{{col}}}\n{table}\\end{{tabular}}"
header = "  \\multirow{3}{2cm}{\\boldsymbol{\\MWR} \\bfseries(GeV)} & \\multirow{3}{*}{} & \\multicolumn{5}{c|}{Electron channel} & \\multicolumn{5}{c|}{Muon channel}  \\\\\n & \\multirow{2}{2cm}{Signal events (stat + syst)} & \\multirow{2}{2cm}{$Z/\\gamma^*$ events (stat + syst)} & \\multirow{2}{2cm}{Top quark events (stat + syst)} & \\multirow{2}{2cm}{Total bkgnd events (stat + syst)} & \\multirow{2}{2cm}{Data events} & \\multirow{2}{2cm}{Signal events (stat + syst)} & \\multirow{2}{2cm}{$Z/\\gamma^*$ events (stat + syst)} & \\multirow{2}{2cm}{Top quark events (stat + syst)} & \\multirow{2}{2cm}{Total bkgnd events (stat + syst)} & \\multirow{2}{2cm}{Data events} \\\\\n"
table = {}  #dictionary named table
masses = set()	#unordered sequence which cannot be indexed, and does not record order of insertion or element position
index = 0   #mass window index used to read correct element in an array branch from root file
with open("configs/mass_cuts.txt") as f:
	for line in f:
		if line[0] == "#": continue
		
		#now split the iterator named line into four strings, two of which (low and hi) will not be used
		ch, mass, low, hi = line.split()

		#initialize floating point variables to hold the number of mean expected and observed evts, and uncertainties in a mass window
		numWREvts = 0.0
		wrStatUnc = 0.0
		wrSystUnc = 0.0
		numDYEvts = 0.0
		dyStatUnc = 0.0
		dySystUnc = 0.0
		numTopEvts = 0.0
		topStatUnc = 0.0
		topSystUnc = 0.0
		numExpEvts = 0.0
		expEvtStatUnc = 0.0
		expEvtSystUnc = 0.0
		numObsEvts = 0.0
		
		#based on the string named ch, fill the two num Evts vars declared immediately above
		#the mass_cuts.txt file is ordered from lowest mass to highest mass, so a simple index which increments after each
		#iteration of this loop will allow the correct array element to be read from the input root files
		if ch == 'EE':
			numDYEvts = combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
			dyStatUnc = combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "ErrorEventsInRange["+str(index)+"]")
			dySystUnc = math.sqrt(pow(0.41*combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchStdDev(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]"), 2) )
			numTopEvts = (0.4315)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
			topStatUnc = (0.4315)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
			topSystUnc = math.sqrt( pow(0.1*(0.4315)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]")), 2)  + pow((0.4315)*(combineTools.getBranchStdDev(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]")), 2) )
			if numTopEvts < 0.1 : topStatUnc = 0.4315
			numExpEvts = numDYEvts + numTopEvts
			expEvtStatUnc = math.sqrt(pow(dyStatUnc,2) + pow(topStatUnc,2) )
			expEvtSystUnc = math.sqrt(pow(dySystUnc,2) + pow(topSystUnc,2) )
			numObsEvts = combineTools.getBranchMean(fileDir+obsEEFile, treeName, "NEventsInRange["+str(index)+"]")
		if ch == 'MuMu':
			numDYEvts = combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			dyStatUnc = combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]")
			dySystUnc = math.sqrt(pow(0.41*combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchStdDev(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]"), 2) )
			numTopEvts = (0.6592)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
			topStatUnc = (0.6592)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
			topSystUnc = math.sqrt( pow(0.1*(0.6592)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]")), 2)  + pow((0.6592)*(combineTools.getBranchStdDev(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]")), 2) )
			if numTopEvts < 0.1 : topStatUnc = 0.6592
			numExpEvts = numDYEvts + numTopEvts
			expEvtStatUnc = math.sqrt(pow(dyStatUnc,2) + pow(topStatUnc,2) )
			expEvtSystUnc = math.sqrt(pow(dySystUnc,2) + pow(topSystUnc,2) )
			numObsEvts = combineTools.getBranchMean(fileDir+obsMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			

		index += 1
		
		#mass_cuts.txt begins with EE channel, then changes to MuMu channel
		#reset index to 1 when line[1] equals 6000
		if mass == "6000": index = 0
		
		#now add a new key and value to the dictionary. the key is a pair containing the mass and lepton channel, and the value is a string containing the number of expected and observed events
		#the numbers 6 and 7 control how many digits and white space are shown in each dict value
		#round controls the number of decimal points shown
		table[(mass , ch)] = "& {sig:<6} & {dy:<6} & {top:<6} & {total:<6} & {obs:<7}".format(sig = (str(round(numWREvts,2))+" +/- "+str(round(wrStatUnc,2))+" +/- "+str(round(wrSystUnc,2)) ), dy = (str(round(numDYEvts,2))+" +/- "+str(round(dyStatUnc,2))+" +/- "+str(round(dySystUnc,2))), top = (str(round(numTopEvts,2))+" +/- "+str(round(topStatUnc,2))+" +/- "+str(round(topSystUnc,2))), total = (str(round(numExpEvts,2))+" +/- "+str(round(expEvtStatUnc,2))+" +/- "+str(round(expEvtSystUnc,2))), obs=str(numObsEvts))
		masses.add(mass)
#end reading mass_cuts.txt

#define the table structure (number of columns and column titles)
tex_table = header
for mass in sorted(masses):
	if mass == "0": continue
	tex_table  += "{mass:4d} {ee} {mumu} \\\\ \\hline\n".format(mass=int(mass), ee=table[(mass, "EE")], mumu = table[(mass, "MuMu")])


print tex.format(table=tex_table,   col="c|c|c|c|c|c|c|c|c|c|c|")
