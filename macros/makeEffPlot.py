#!/usr/bin/env python

#run this macro via ./makeEffPlot configEffPlot.json
#where configEffPlot.json defines the variables used to make the plot

import argparse
import sys
import os
import json

"""
Setup argument parser
"""

parser = argparse.ArgumentParser(description="This program takes an input JSON config and extracts plots from Tag-and-Probe ROOT files. The output consists of a plot with superimposed graphs from multiple TnP files and the according fit canvases.")
parser.add_argument("inputJsonConfig", help="Path to the input JSON config file")
parser.add_argument("-f", "--fast", default=0, action="count", help="Skip fetching and saving the fit canvases for each plot")
parser.add_argument("-v", "--verbosity", default=1, help="Increase or decrease output verbosity")
args = parser.parse_args()

"""
Parse JSON file
"""

with open(args.inputJsonConfig, 'r') as f:
	data = json.loads(f.read())

"""
Go through plots defined in config JSON
"""

from ROOT import * # import this here, otherwise it overwrites the argparse stuff
gROOT.SetBatch(True) # set ROOT to batch mode, this suppresses printing canvases
gROOT.ProcessLine("gErrorIgnoreLevel = 1001;") # suppress stdout pollution of canvas.Print(...)

for keyPlot in data:
	if args.verbosity==1:
		print('Processing plot config: {}'.format(keyPlot))
		print('Config comment: {}'.format(data[keyPlot]['comment']))

	# Get result plots from fit canvases (in dir fit_eff_plots)

	# Get input parameters
	inputFilenames = []
	inputPaths = []
	inputBinnedVariables = []
	inputLegendEntries = []
	inputLabels = []
	for keyInputs in data[keyPlot]['inputs']:
		inputFilenames.append(data[keyPlot]['inputs'][keyInputs]['filename'])
		inputPaths.append(os.path.join(data[keyPlot]['inputs'][keyInputs]['directory'],'fit_eff_plots'))
		inputBinnedVariables.append(data[keyPlot]['inputs'][keyInputs]['binnedVariable'])
		inputLegendEntries.append(data[keyPlot]['inputs'][keyInputs]['legendEntry'])
		inputLabels.append(data[keyPlot]['inputs'][keyInputs]['label'])

	# Get input graphs from files
	inputGraphs = []
	for iGraph in range(len(inputFilenames)):
		inputFile = TFile.Open(inputFilenames[iGraph])
		if not inputFile:
			print('[ERROR] File not found: {}'.format(inputFilenames[iGraph]))
			sys.exit()
		inputDir = inputFile.GetDirectory(inputPaths[iGraph])
		inputName = None
		for keys in inputDir.GetListOfKeys():
			if keys.GetName()[0:len(inputBinnedVariables[iGraph])]==inputBinnedVariables[iGraph]:
				inputName = keys.GetName()
		if args.verbosity==1:
			print('Load plot \'{}\': {}'.format(inputLabels[iGraph], inputName))
		inputGraphs.append(inputDir.Get(inputName).GetPrimitive('hxy_fit_eff'))

	# Set line color and marker style for each graph using given maps from config
	colorMap = data[keyPlot]['plot']['colorMap']
	markerMap = data[keyPlot]['plot']['markerMap']
	if args.verbosity==1:
		print('Using colormap: {}'.format(colorMap))
		print('Using markermap: {}'.format(markerMap))
	if len(colorMap)<len(inputGraphs):
		print('[ERROR] The defined colormap has not enough entries for the number of defined input graphs.')
		sys.exit()
	if len(markerMap)<len(inputGraphs):
		print('[ERROR] The defined markermap has not enough entries for the number of defined input graphs.')
		sys.exit()

	for iGraph in range(len(inputGraphs)):
		inputGraphs[iGraph].SetLineColor(colorMap[iGraph])
		inputGraphs[iGraph].SetMarkerStyle(markerMap[iGraph])
		inputGraphs[iGraph].SetMarkerColor(colorMap[iGraph])

	# Generate superimposed graph using TMultiGraph
	mg = TMultiGraph()
	for graph in inputGraphs:
		mg.Add(graph)

	# Setup canvas with all elements
	canvas = TCanvas('canvas', 'canvas', 800, 800)

	pad = TPad('pad', 'pad', 0.01, 0.00, 1.00, 1.00)
	pad.Draw()
	pad.cd()

	plotX = data[keyPlot]['plot']['x']
	plotY = data[keyPlot]['plot']['y']
	mg.Draw('AP')
	mg.GetXaxis().SetRangeUser(plotX[0], plotX[1])
	mg.GetXaxis().SetTitle(plotX[2])
	mg.GetXaxis().SetLabelSize(22)
	mg.GetXaxis().SetTitleFont(63)
	mg.GetXaxis().SetLabelFont(43)
	mg.GetXaxis().SetTitleSize(22)
	mg.GetXaxis().SetLabelSize(20)
	mg.GetXaxis().SetTitleOffset(1.2)
	mg.GetYaxis().SetRangeUser(plotY[0], plotY[1])
	mg.GetYaxis().SetTitle(plotY[2])
	mg.GetYaxis().SetLabelSize(22)
	mg.GetYaxis().SetTitleFont(63)
	mg.GetYaxis().SetLabelFont(43)
	mg.GetYaxis().SetTitleSize(22)
	mg.GetYaxis().SetLabelSize(20)
	mg.GetYaxis().SetTitleOffset(1.5)

	canvas.cd()
	leg = TLegend(0.37, 0.76, 0.75, 0.88)
	leg.SetHeader(data[keyPlot]['plot']['legendTitle'])
	header = leg.GetListOfPrimitives().First()
	header.SetTextColor(1)
	header.SetTextFont(43)
	header.SetTextSize(20)
	for iGraph in range(len(inputGraphs)):
		leg.AddEntry(inputGraphs[iGraph], inputLegendEntries[iGraph], 'LP')
	leg.SetBorderSize(0)
	leg.SetTextFont(43)
	leg.SetTextSize(20)
	leg.Draw()

	canvas.cd()
	latex = TLatex()
	latex.SetNDC()
	latex.SetTextFont(61)
	latex.SetTextSize(0.06)
	latex.DrawLatex(0.16, 0.82, data[keyPlot]['plot']['logo'][0])
	latex.SetTextFont(52)
	latex.SetTextSize(0.04)
	latex.SetTextAlign(11);
	latex.DrawLatex(0.16, 0.77, data[keyPlot]['plot']['logo'][1])
	latex.SetTextFont(42)
	latex.SetTextSize(0.038)
	latex.SetTextAlign(31);
	latex.DrawLatex(0.90, 0.91, data[keyPlot]['plot']['caption'])
	canvas.Update()

	# Save plot
	outputDirectory = data[keyPlot]['output']['directory']
	if args.verbosity==1:
		print('Output directory: {}'.format(outputDirectory))
	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)
	for fileType in data[keyPlot]['output']['fileType']:
		canvas.SaveAs(os.path.join(outputDirectory,data[keyPlot]['output']['filenamePlot']+'.'+fileType))

	# Skip fetching fit canvases if flag is set
	if args.fast != 0:
		if args.verbosity==1:
			print('')
		continue

	# Make directories for fit canvases
	if args.verbosity==1:
		print('Output directories for fit canvases: {}'.format(inputLabels))
	outputDirectoryInputs = []
	for label in inputLabels:
		outputDirectoryInputs.append(os.path.join(outputDirectory, label))
		if not os.path.exists(outputDirectoryInputs[-1]):
			os.makedirs(outputDirectoryInputs[-1])

	# Get fit canvases and store them to output sub-directory
	inputPathsFit = []
	for keyInputs in data[keyPlot]['inputs']:
		inputPathsFit.append(data[keyPlot]['inputs'][keyInputs]['directory'])
	for iFile in range(len(inputFilenames)):
		inputFile = TFile(inputFilenames[iFile])
		if not inputFile:
			print('[ERROR] Cannot find file: {}'.format(inputFilenames[iFile]))
			sys.exit()
		inputDir = inputFile.GetDirectory(inputPathsFit[iFile])
		if not inputDir:
			print('[ERROR] Cannot find directory (file, directory): {}, {}'.format(inputFilenames[iFile], inputPathsFit[iFile]))
			sys.exit()
		for keys in inputDir.GetListOfKeys():
			if TString(keys.GetName()).Contains(inputBinnedVariables[iFile]):
				pathDirBin = os.path.join(inputPathsFit[iFile], keys.GetName())
				dirBin = inputFile.GetDirectory(pathDirBin)
				for fileType in data[keyPlot]['output']['fileType']:
					fitCanvas = dirBin.Get('fit_canvas') # NOTE you have to put this here, toherwise the loop won't work
					pathOutput = os.path.join(outputDirectoryInputs[iFile], keys.GetName()+'.'+fileType)
					if not fitCanvas:
						print('[WARNING] Found missing fit canvas (file, directory, directory bin): {}, {}, {}'.format(inputFilenames[iFile],inputPathsFit[iFile],pathDirBin))
					else:
						fitCanvas.SaveAs(pathOutput)

	if args.verbosity==1:
		print('')
