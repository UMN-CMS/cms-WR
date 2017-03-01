import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np
import math

#change fileDir to point to the directory from which expected bkgnd and signal, and observed files should be read
treeName="syst_tree"
fileDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysAllSystSmoothedWindowsFebrTwentyOne/"

expTopMuMuFile = "selected_tree_data_flavoursidebandEMuMuMu.root"
expDyMuMuFile = "selected_tree_DYAMC_signal_mumuMuMu.root"
obsMuMuFile = "selected_tree_data_signal_mumuMuMu.root"
expTopEEFile = "selected_tree_data_flavoursidebandEMuEE.root"
expDyEEFile = "selected_tree_DYAMC_signal_eeEE.root"
obsEEFile = "selected_tree_data_signal_eeEE.root"


##now make a LaTex table with one column showing the WR mass hypothesis, two columns showing the number of expected and obs events in the ele channel
##and two columns showing the number of expected and obs events in the muon channel
tex= "\n\\begin{{tabular}}{{{col}}}\\hline\n{table}\\end{{tabular}}"
header = " & \\multicolumn{2}{c|}{Electrons}  & \\multicolumn{2}{c|}{Muons}  \\\\  \\hline\nMass(GeV) & Expected(Evts) & Observed(Evts) & Expected(Evts) & Observed(Evts) \\\\  \\hline\n"
table = {}  #dictionary named table
masses = set()	#unordered sequence which cannot be indexed, and does not record order of insertion or element position
index = 0   #mass window index used to read correct element in an array branch from root file
with open("configs/mass_cuts.txt") as f:
	for line in f:
		if line[0] == "#": continue
		
		#now split the iterator named line into four strings, two of which (low and hi) will not be used
		ch, mass, low, hi = line.split()
		
		#initialize two floating point variables which will hold the number of expected and observed evts in one mass window
		numObsEvts = 0.0
		numExpEvts = 0.0
		expEvtUnc = 0.0
		
		#based on the string named ch, fill the two num Evts vars declared immediately above
		#the mass_cuts.txt file is ordered from lowest mass to highest mass, so a simple index which increments after each
		#iteration of this loop will allow the correct array element to be read from the input root files
		if ch == 'EE':
			numObsEvts = combineTools.getBranchMean(fileDir+obsEEFile, treeName, "NEventsInRange["+str(index)+"]")
			numExpEvts += combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
			numExpEvts += (0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
			expEvtUnc += math.sqrt(pow(0.4*combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "ErrorEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchStdDev(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(0.05*(0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]")), 2)  + pow((0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "ErrorEventsInRange["+str(index)+"]")), 2) + pow((0.4194)*(combineTools.getBranchStdDev(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]")), 2) )
		if ch == 'MuMu':
			numObsEvts = combineTools.getBranchMean(fileDir+obsMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			numExpEvts += combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			numExpEvts += (0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
			expEvtUnc += math.sqrt(pow(0.4*combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]"), 2) + pow(combineTools.getBranchStdDev(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]"), 2) + pow(0.05*(0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]")), 2)  + pow((0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]")), 2) + pow((0.6563)*(combineTools.getBranchStdDev(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]")), 2) )
	

		index += 1
		
		#mass_cuts.txt begins with EE channel, then changes to MuMu channel
		#reset index to 1 when line[1] equals 6000
		if mass == "6000": index = 0
		
		#now add a new key and value to the dictionary. the key is a pair containing the mass and lepton channel, and the value is a string containing the number of expected and observed events
		#the numbers 6 and 7 control how many digits and white space are shown in each dict value
		table[(mass , ch)] = "& {low:<6} & {hi:<7}".format(low = (str(round(numExpEvts,2))+" +/- "+str(round(expEvtUnc,2))), hi=str(numObsEvts))
		masses.add(mass)
#end reading mass_cuts.txt

#define the table structure (number of columns and column titles)
tex_table = header
for mass in sorted(masses):
	if mass == "0": continue
	tex_table  += "{mass:4d} {ee} {mumu} \\\\ \\hline\n".format(mass=int(mass), ee=table[(mass, "EE")], mumu = table[(mass, "MuMu")])



print tex.format(table=tex_table,   col="|c|c|c|c|c|c|c|")
