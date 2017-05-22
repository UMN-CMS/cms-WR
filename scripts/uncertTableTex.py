import subprocess
import ExoAnalysis.cmsWR.combineTools as combineTools
import numpy as np
import math

##################
#THIS SCRIPT IS NO LONGER NEEDED
#the table produced by this script has been replaced by the table produced by expEvtsWithUncsAndObs_tex.py
#
#
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

#####################
#in header for original style table:
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
#header for original style table, which showed the shape syst unc in three categories (all, jet, lepton) for both lepton channels
#header = "\\multirow{3}{1cm}{\MWR (GeV)} & \\multirow{3}{*}{} & \\multicolumn{5}{c||}{Electron channel } & \\multicolumn{5}{c||}{Muon channel}                                                                                          \\\\\n                              &                          & \\multirow{2}{1cm}{Exp. events}          & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{3}{c||}{Syst. Unc} & \\multirow{2}{1cm}{Exp. events} & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{3}{c||}{Syst. Unc.}  \\\\  \n                              &                          &                                         &            & all                             & jets                        & leptons       &  &  & all & jets & leptons       \\\\ \n\\hline\n"
####################

#####################
#in header for new style table with only total shape syst unc and normalization unc (one less column per lepton channel than the original table defined above):
header = "\\multirow{3}{1cm}{\MWR (GeV)} & \\multirow{3}{*}{} & \\multicolumn{4}{c||}{Electron channel } & \\multicolumn{4}{c||}{Muon channel}                                                                                          \\\\\n                              &                          & \\multirow{2}{1cm}{Exp. events}          & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{2}{c||}{Syst. Unc} & \\multirow{2}{1cm}{Exp. events} & \\multirow{2}{1cm}{Stat. Unc.} & \\multicolumn{2}{c||}{Syst. Unc.}  \\\\  \n                              &                          &                                         &            & shape                             & norm                        &   &  & shape & norm      \\\\ \n\\hline\n"
####################

table = {}  #dictionary named table
masses = set()	#unordered sequence which cannot be indexed, and does not record order of insertion or element position
index = 0   #mass window index used to read correct element in an array branch from root file

#read mass_cuts.txt file and only write contents to the table when specific WR mass hypotheses are read from mass_cuts.txt
#because of how the table header is defined, electron and muon results should be saved to temporary float vars at the same time
with open("configs/mass_cuts.txt") as f:
	for line in f:
		if line[0] == "#": continue
		
		#now split the iterator named line into four strings, two of which (low and hi) will not be used
		ch, mass, low, hi = line.split()
		
		#initialize floating point variables which will hold the number of expected evts and their uncertainties in one mass window
		#electron channel
		dyEEEvts = 0.0
		dyEEStatUnc = 0.0
		dyEEShapeUncFromAll = 0.0
		dyEENormUnc = 0.0
		#dyEEShapeUncFromJets = 0.0
		#dyEEShapeUncFromLepts = 0.0
		topEEEvts = 0.0
		topEEStatUnc = 0.0
		topEEShapeUncFromAll = 0.0
		topEENormUnc = 0.0
		#topEEShapeUncFromJets = 0.0
		#topEEShapeUncFromLepts = 0.0
		wrEEEvts = 0.0
		wrEEStatUnc = 0.0
		wrEEShapeUncFromAll = 0.0
		wrEENormUnc = 0.0
		#wrEEShapeUncFromJets = 0.0
		#wrEEShapeUncFromLepts = 0.0
		
		#muon channel
		dyMuMuEvts = 0.0
		dyMuMuStatUnc = 0.0
		dyMuMuShapeUncFromAll = 0.0
		dyMuMuNormUnc = 0.0
		#dyMuMuShapeUncFromJets = 0.0
		#dyMuMuShapeUncFromLepts = 0.0
		topMuMuEvts = 0.0
		topMuMuStatUnc = 0.0
		topMuMuShapeUncFromAll = 0.0
		topMuMuNormUnc = 0.0
		#topMuMuShapeUncFromJets = 0.0
		#topMuMuShapeUncFromLepts = 0.0
		wrMuMuEvts = 0.0
		wrMuMuStatUnc = 0.0
		wrMuMuShapeUncFromAll = 0.0
		wrMuMuNormUnc = 0.0
		#wrMuMuShapeUncFromJets = 0.0
		#wrMuMuShapeUncFromLepts = 0.0
		
		#the mass_cuts.txt file is ordered from lowest mass to highest mass, so a simple index which increments after each
		#iteration of this loop will allow the correct array element to be read from the input root files
		if mass == "1600" or mass == "2200" or mass == "2800":
			##electron channel
			dyEEEvts += combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
			dyEEStatUnc += combineTools.getBranchMean(fileDir+expDyEEFile, treeName, "ErrorEventsInRange["+str(index)+"]")
			dyEEShapeUncFromAll += combineTools.getBranchStdDev(fileDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
			dyEENormUnc += 0.4*dyEEEvts
			#dyEEShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")
			#dyEEShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+expDyEEFile, treeName, "NEventsInRange["+str(index)+"]")

			topEEEvts += (0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
			topEEStatUnc += (0.4194)*(combineTools.getBranchMean(fileDir+expTopEEFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
			topEEShapeUncFromAll += (0.4194)*(combineTools.getBranchStdDev(fileDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
			topEENormUnc += topEEEvts*math.sqrt(0.027*0.027 + 0.05*0.05)
			#topEEShapeUncFromJets += (0.4194)*(combineTools.getBranchStdDev(jetShapeSystDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))
			#topEEShapeUncFromLepts += (0.4194)*(combineTools.getBranchStdDev(leptShapeSystDir+expTopEEFile, treeName, "NEventsInRange["+str(index)+"]"))

			wrEEEvts += combineTools.getBranchMean(fileDir+"selected_tree_WRtoEEJJ_%i_%i_signal_eeEE.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			wrEEStatUnc += combineTools.getBranchMean(fileDir+"selected_tree_WRtoEEJJ_%i_%i_signal_eeEE.root" % (int(mass), (int(mass)/2)), treeName, "ErrorEventsInRange["+str(index)+"]")
			wrEEShapeUncFromAll += combineTools.getBranchStdDev(fileDir+"selected_tree_WRtoEEJJ_%i_%i_signal_eeEE.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			wrEENormUnc += wrEEEvts*0.027
			#wrEEShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+"selected_tree_WRtoEEJJ_%i_%i_signal_eeEE.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			#wrEEShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+"selected_tree_WRtoEEJJ_%i_%i_signal_eeEE.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			
			##muon channel
			dyMuMuEvts += combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			dyMuMuStatUnc += combineTools.getBranchMean(fileDir+expDyMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]")
			dyMuMuShapeUncFromAll += combineTools.getBranchStdDev(fileDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			dyMuMuNormUnc += 0.4*dyMuMuEvts
			#dyMuMuShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
			#dyMuMuShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+expDyMuMuFile, treeName, "NEventsInRange["+str(index)+"]")
				
			topMuMuEvts += (0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
			topMuMuStatUnc += (0.6563)*(combineTools.getBranchMean(fileDir+expTopMuMuFile, treeName, "ErrorEventsInRange["+str(index)+"]"))
			topMuMuShapeUncFromAll += (0.6563)*(combineTools.getBranchStdDev(fileDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
			topMuMuNormUnc += topMuMuEvts*math.sqrt(0.027*0.027 + 0.05*0.05)
			#topMuMuShapeUncFromJets += (0.6563)*(combineTools.getBranchStdDev(jetShapeSystDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
			#topMuMuShapeUncFromLepts += (0.6563)*(combineTools.getBranchStdDev(leptShapeSystDir+expTopMuMuFile, treeName, "NEventsInRange["+str(index)+"]"))
				
			wrMuMuEvts += combineTools.getBranchMean(fileDir+"selected_tree_WRtoMuMuJJ_%i_%i_signal_mumuMuMu.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			wrMuMuStatUnc += combineTools.getBranchMean(fileDir+"selected_tree_WRtoMuMuJJ_%i_%i_signal_mumuMuMu.root" % (int(mass), (int(mass)/2)), treeName, "ErrorEventsInRange["+str(index)+"]")
			wrMuMuShapeUncFromAll += combineTools.getBranchStdDev(fileDir+"selected_tree_WRtoMuMuJJ_%i_%i_signal_mumuMuMu.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			wrMuMuNormUnc += wrMuMuEvts*0.027
			#wrMuMuShapeUncFromJets += combineTools.getBranchStdDev(jetShapeSystDir+"selected_tree_WRtoMuMuJJ_%i_%i_signal_mumuMuMu.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
			#wrMuMuShapeUncFromLepts += combineTools.getBranchStdDev(leptShapeSystDir+"selected_tree_WRtoMuMuJJ_%i_%i_signal_mumuMuMu.root" % (int(mass), (int(mass)/2)), treeName, "NEventsInRange["+str(index)+"]")
	
			##now add a new key and value to the dictionary. the key is the WR mass hypothesis. the value is everything shown in one row - expected events and uncertainties for both lepton channels
			multirowStr = "\\multirow{3}{*}{"+str(mass)+"}"
			
			#original row entry for table showing 3 shape uncertainties - all, jet only, lepton only
			#rowEntryStr = " & \\WR  & {expcEEWR} & {statEEWR} & {allShpEEWR} & {jetShpEEWR} & {leptShpEEWR} & {expcMuMuWR} & {statMuMuWR} & {allShpMuMuWR} & {jetShpMuMuWR} & {leptShpMuMuWR} \\\\ \n & DY  & {expcEEDY} & {statEEDY} & {allShpEEDY} & {jetShpEEDY} & {leptShpEEDY} & {expcMuMuDY} & {statMuMuDY} & {allShpMuMuDY} & {jetShpMuMuDY} & {leptShpMuMuDY} \\\\ \n & top  & {expcEEtop} & {statEEtop} & {allShpEEtop} & {jetShpEEtop} & {leptShpEEtop} & {expcMuMutop} & {statMuMutop} & {allShpMuMutop} & {jetShpMuMutop} & {leptShpMuMutop}".format(expcEEWR = str(round(wrEEEvts,2)), statEEWR = str(round(wrEEStatUnc,2)), allShpEEWR = str(round(wrEEShapeUncFromAll,2)), jetShpEEWR = str(round(wrEEShapeUncFromJets,2)), leptShpEEWR = str(round(wrEEShapeUncFromLepts,2)), expcMuMuWR = str(round(wrMuMuEvts,2)), statMuMuWR = str(round(wrMuMuStatUnc,2)), allShpMuMuWR = str(round(wrMuMuShapeUncFromAll,2)), jetShpMuMuWR = str(round(wrMuMuShapeUncFromJets,2)), leptShpMuMuWR = str(round(wrMuMuShapeUncFromLepts,2)), expcEEDY = str(round(dyEEEvts,2)), statEEDY = str(round(dyEEStatUnc,2)), allShpEEDY = str(round(dyEEShapeUncFromAll,2)), jetShpEEDY = str(round(dyEEShapeUncFromJets,2)), leptShpEEDY = str(round(dyEEShapeUncFromLepts,2)), expcMuMuDY = str(round(dyMuMuEvts,2)), statMuMuDY = str(round(dyMuMuStatUnc,2)), allShpMuMuDY = str(round(dyMuMuShapeUncFromAll,2)), jetShpMuMuDY = str(round(dyMuMuShapeUncFromJets,2)), leptShpMuMuDY = str(round(dyMuMuShapeUncFromLepts,2)), expcEEtop = str(round(topEEEvts,2)), statEEtop = str(round(topEEStatUnc,2)), allShpEEtop = str(round(topEEShapeUncFromAll,2)), jetShpEEtop = str(round(topEEShapeUncFromJets,2)), leptShpEEtop = str(round(topEEShapeUncFromLepts,2)), expcMuMutop = str(round(topMuMuEvts,2)), statMuMutop = str(round(topMuMuStatUnc,2)), allShpMuMutop = str(round(topMuMuShapeUncFromAll,2)), jetShpMuMutop = str(round(topMuMuShapeUncFromJets,2)), leptShpMuMutop = str(round(topMuMuShapeUncFromLepts,2)) )
			
			#row entry for table showing 2 syst uncertainties - lept+jet shape unc, and normalization unc
			rowEntryStr = " & \\WR  & {expcEEWR} & {statEEWR} & {allShpEEWR} & {normEEWR} & {expcMuMuWR} & {statMuMuWR} & {allShpMuMuWR} & {normMuMuWR} \\\\ \n & DY  & {expcEEDY} & {statEEDY} & {allShpEEDY} & {normEEDY} & {expcMuMuDY} & {statMuMuDY} & {allShpMuMuDY} & {normMuMuDY} \\\\ \n & top  & {expcEEtop} & {statEEtop} & {allShpEEtop} & {normEEtop} & {expcMuMutop} & {statMuMutop} & {allShpMuMutop} & {normMuMutop}".format(expcEEWR = str(round(wrEEEvts,2)), statEEWR = str(round(wrEEStatUnc,2)), allShpEEWR = str(round(wrEEShapeUncFromAll,2)), normEEWR = str(round(wrEENormUnc,2)), expcMuMuWR = str(round(wrMuMuEvts,2)), statMuMuWR = str(round(wrMuMuStatUnc,2)), allShpMuMuWR = str(round(wrMuMuShapeUncFromAll,2)), normMuMuWR = str(round(wrMuMuNormUnc,2)), expcEEDY = str(round(dyEEEvts,2)), statEEDY = str(round(dyEEStatUnc,2)), allShpEEDY = str(round(dyEEShapeUncFromAll,2)), normEEDY = str(round(dyEENormUnc,2)), expcMuMuDY = str(round(dyMuMuEvts,2)), statMuMuDY = str(round(dyMuMuStatUnc,2)), allShpMuMuDY = str(round(dyMuMuShapeUncFromAll,2)), normMuMuDY = str(round(dyMuMuNormUnc,2)), expcEEtop = str(round(topEEEvts,2)), statEEtop = str(round(topEEStatUnc,2)), allShpEEtop = str(round(topEEShapeUncFromAll,2)), normEEtop = str(round(topEENormUnc,2)), expcMuMutop = str(round(topMuMuEvts,2)), statMuMutop = str(round(topMuMuStatUnc,2)), allShpMuMutop = str(round(topMuMuShapeUncFromAll,2)), normMuMutop = str(round(topMuMuNormUnc,2)) )
			
			table[(mass)] = multirowStr + rowEntryStr
			masses.add(mass)

		##end if mass == 1600, 2200 or 2800
	
		#mass_cuts.txt begins with EE channel, then changes to MuMu channel
		#leave this loop after the last EE channel line in mass_cuts.txt is read
		if mass == "6000": break
	
		##increment the index even if nothing is added to the table. this index is only used to access the correct mass window in various branches of the input tree
		index += 1
		
##end reading mass_cuts.txt

##define the table structure (number of columns and column titles)
tex_table = header
for mass in sorted(masses):
	#below the header defined earlier, add a row showing the expected number of events, stat unc, and syst unc in both lepton channels for all signal and bkgnd processes
	#tex_table  += "{mass:4d} {ee} {mumu} \\\\ \n\\hline\n".format(mass=int(mass), ee=table[(mass, "EE")], mumu = table[(mass, "MuMu")])
	tex_table  += "{entry} \\\\ \n\\hline\n".format(entry=table[(mass)])


#this col string, setup for the original table with 3 shape uncertainties, creates a double verticle line after the 2nd column, makes the 5th
#column bold, and adds a double verticle line 2 columns to the right of the bolded column
#print tex.format(table=tex_table,   col="c|c||c|c|>{\\bfseries}c*{2}{c}||c|c|>{\\bfseries}c*{2}{c}||")

#col string for the new table, which shows one lept+jet shape syst unc, and one norm syst unc. no bolded columns
print tex.format(table=tex_table,   col="c|c||c|c|c|c||c|c|c|c||")
