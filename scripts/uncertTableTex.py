import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np
import math

##################
#this script makes a table showing the expected number of bkgnd and signal events, and the errors on these expectations due to
#due to uncertainties in the MLLJJ normalization and shape
#this script prints the table to stdout, and this should be directed to a file named uncertainty_table.tex
##################

#change fileDir to point to the directory from which expected bkgnd and signal, and observed files should be read
treeName="syst_tree"

#these three directories are needed to estimate the uncertainty from all shape systematics (jets and leptons), only jet syst, and only lept syst
fileDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysAllSystSmoothedWindowsFebrTwentyOne/"
jetShapeSystDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysJetSystNewMassWindows/"
leptShapeSystDir="/afs/cern.ch/work/s/skalafut/public/WR_starting2015/processedWithAnalysisCpp/3200toysAllLeptSystNewMassWindows/"

expTopMuMuFile = "selected_tree_data_flavoursidebandEMuMuMu.root"
expDyMuMuFile = "selected_tree_DYMadInclAndHT_signal_mumuMuMu_withMllWeight.root"
expTopEEFile = "selected_tree_data_flavoursidebandEMuEE.root"
expDyEEFile = "selected_tree_DYMadInclAndHT_signal_eeEE_withMllWeight.root"


#final version tex= "\\centering\n\\begin{{tabular}}{{{col}}}\n{table}\\end{{tabular}}"
#temp
tex= "\\centering\n\\begin{{tabular}}{{{col}}}\n{table}\n\\hline\n\\end{{tabular}}"

#in header:
#  multirow{3}{1cm}{\MWR (GeV)} makes 3 rows each with one entry which is 1 centimeter tall, and filled with the string \MWR (Gev)
#  multirow{3}{*}{} appends 3 new elements to the right of the 3 existing rows (so now each row has 2 elements), but does not fill the new elements with anything
#  multicolumn{5}{c||}{Electron channel } appends a double verticle bar followed by 5 new elements to the right side of only the top row (so now the top row has 7 elements, while the lower rows have 2).
#    Written across all five elements is the string Electron channel, and the lower rows have the structure to hold 5 columns of elements, but nothing yet is shown in those elements 
#  multicolumn{5}{c||}{Muon channel} makes a similar addition only to the top row, and displays the string Muon channel. this is the last addition to the top row.
#  NEXT LINE
#  multirow{2}{1cm}{Exp. events} makes 2 rows which are 1 centimeter tall, and sets the top element to Exp., and the bottom element to events. this structure is used such that different
#    systematic uncertainties can be shown under the umbrella of "Syst. Unc" written in the top row
#  multirow{2}{1cm}{Stat. Unc.} appends 2 new elements, each 1 centimeter tall, on the right end of the two new rows. now each row has 2 elements
#  multicolumn{3}{c||}{Syst. Unc} appends a double verticle bar followed by 3 new elements on the right side of the top row. Now the top row has 5 elements labeled "Exp. Events | Stat. Unc. |   Syst. Unc  |"
#    underneath the Electron channel header
#  the next two multirow and multicolumn commands make an identical 5 column structure under the Muon channel header
#  NEXT LINE
#  a single row is added with several columns, most of which are empty. this row identifies the columns under the Syst. Unc header which represent all shape syst unc, shape syst unc due to jets, and shape syst unc due to leptons
header = "\\multirow{3}{1cm}{\MWR (GeV)} & \\multirow{3}{*}{} & \\multicolumn{5}{c||}{Electron channel } & \\multicolumn{5}{c||}{Muon channel}                                                                                          \\\\\n                              &                          & \\multirow{2}{1cm}{Exp. events}          & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{3}{c||}{Syst. Unc} & \\multirow{2}{1cm}{Exp. events} & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{3}{c||}{Syst. Unc.}  \\\\  \n                              &                          &                                         &            & all                             & jets                        & leptons       &  &  & all & jets & leptons       \\\\ \n\\hline"

table = {}  #dictionary named table
masses = set()	#unordered sequence which cannot be indexed, and does not record order of insertion or element position
index = 0   #mass window index used to read correct element in an array branch from root file

#read mass_cuts.txt file and only write contents to the table when specific WR mass hypotheses are read from mass_cuts.txt
with open("configs/mass_cuts.txt") as f:
	for line in f:
		if line[0] == "#": continue
		
		#now split the iterator named line into four strings, two of which (low and hi) will not be used
		ch, mass, low, hi = line.split()
		
		#initialize floating point variables which will hold the number of expected evts and their uncertainties in one mass window
		dyEvts = 0.0
		dyStatUnc = 0.0
		dyShapeUncFromAll = 0.0
		dyShapeUncFromJets = 0.0
		dyShapeUncFromLepts = 0.0
		topEvts = 0.0
		topStatUnc = 0.0
		topShapeUncFromAll = 0.0
		topShapeUncFromJets = 0.0
		topShapeUncFromLepts = 0.0
		wrEvts = 0.0
		wrStatUnc = 0.0
		wrShapeUncFromAll = 0.0
		wrShapeUncFromJets = 0.0
		wrShapeUncFromLepts = 0.0
		
		#the mass_cuts.txt file is ordered from lowest mass to highest mass, so a simple index which increments after each
		#iteration of this loop will allow the correct array element to be read from the input root files
		if mass == "1600" or mass == "2200" or mass == "2800":
			if ch == 'EE':
				dyEvts += combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
				dyStatUnc += combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "ErrorEventsInRange["+str(index)+"]")
				dyShapeUncFromAll += combineTools.getBranchStdDev(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
				dyShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
				dyShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
				
				topEvts += (0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
				topStatUnc += (0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
				topShapeUncFromAll += (0.4194)*(combineTools.getBranchStdDev(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
				topShapeUncFromJets += (0.4194)*(combineTools.getBranchStdDev(jetShapeSystDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
				topShapeUncFromLepts += (0.4194)*(combineTools.getBranchStdDev(leptShapeSystDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
				
				wrEvts +=
				
			
			if ch == 'MuMu':
				dyEvts += combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
				dyStatUnc += combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]")
				dyShapeUncFromAll += combineTools.getBranchStdDev(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
				dyShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
				dyShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
				
				topEvts += (0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
				topStatUnc += (0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
				topShapeUncFromAll += (0.6563)*(combineTools.getBranchStdDev(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
				topShapeUncFromJets += (0.6563)*(combineTools.getBranchStdDev(jetShapeSystDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
				topShapeUncFromLepts += (0.6563)*(combineTools.getBranchStdDev(leptShapeSystDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
				
	
		#end if mass == 1600, 2200 or 2800 selection
	

		#increment the index even if nothing is added to the table. this index is only used to access the correct mass window in various branches of the input tree
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



print tex.format(table=tex_table,   col="c|c||c|c|>{\\bfseries}c*{2}{c}||c|c|>{\\bfseries}c*{2}{c}||")
