#additional VarParsing options which are used in multiple python cfg files in test/
import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('standard')

options.register('GT',
		"",
		VarParsing.VarParsing.multiplicity.singleton,
		VarParsing.VarParsing.varType.string,
		"global tag name")

options.maxEvents = -1
options.parseArguments()

#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v1', '') for Run2015C data
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v2', '') for Run2015D data
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v2', '') for re-MINIAOD Spring15MC
#process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '') for WR signal Spring15MC



