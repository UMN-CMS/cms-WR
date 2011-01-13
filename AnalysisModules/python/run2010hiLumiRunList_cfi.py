import FWCore.ParameterSet.Config as cms

# Reading configuration file from  runreg_7TeVCollisions10StreamExpress.cfg
# You asked for the runreg info in the run range:132440-149442
# for dataset: 
# /StreamExpress/%Express%:132440-149181
# /Run2010B-PromptReco-v2/DQM:149182-149442
# with the following quality flags:
# Muon GOOD
# Egam GOOD
# Jmet GOOD
# Lumi GOOD
# Ecal GOOD
# L1t GOOD
# Hlt GOOD
# Csc GOOD
# Pix GOOD
# Track GOOD
# Rpc GOOD
# Strip GOOD
# Hcal GOOD
# Dt GOOD
# Es GOOD
# and with the following DCS status:
# Bpix
# Fpix
# Tibtid
# TecM
# TecP
# Tob
# Ebminus
# Ebplus
# EeMinus
# EePlus
# EsMinus
# EsPlus
# HbheA
# HbheB
# HbheC
# Ho
# Hf
# Dtminus
# Dtplus
# Dt0
# CscMinus
# CscPlus
# Rpc
# Manual bad LS in comment column: True
# Beam energy requested: 3500
# Cross-check requested agains DBS in PDs: ['/MinimumBias/Commissioning10-v4/RAW', ' /MinimumBias/Run2010A-v1/RAW', ' /MinimumBias/Run2010B-v1/RAW']
# RunRegistry from:  http://pccmsdqm04.cern.ch/runregistry/xmlrpc
 
# Accessing run registry....
 
# Building up beam energy per run info... please wait it can be long
# Retrieving energy from the RR run table and making some check....
# WARNING: Something wrong with energies in run 147757
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 147048
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 143665
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 142558
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 142076
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 140401
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 140359
# WARNING: Getting: " from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 133082
# WARNING: Getting: 3143877 from RR.Using default value of:3500.0
# WARNING: Something wrong with energies in run 133031
# WARNING: Getting: 7864 from RR.Using default value of:3500.0
# Accessing DBS to get lumi mapping: this can be VERY long....
 
#-------------------------------------------
#Json file:  Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.txt  written.
#-------------------------------------------
 
# CFG snippet to select:
lumisToProcess = cms.untracked.VLuminosityBlockRange(
	'147196:1-147196:90',
	'147214:1-147214:79',
	'147216:1-147216:63',
	'147217:1-147217:193',
	'147218:1-147218:45',
	'147219:1-147219:293',
	'147219:309-147219:320',
	'147222:1-147222:444',
	'147284:32-147284:306',
	'147390:1-147390:478',
	'147390:480-147390:839',
	'147450:80-147450:166',
	'147451:1-147451:116',
	'147451:118-147451:129',
	'147452:1-147452:44',
	'147453:1-147453:146',
	'147454:1-147454:97',
	'147754:1-147754:167',
	'147754:170-147754:377',
	'147755:81-147755:231',
	'147757:1-147757:363',
	'147926:77-147926:548',
	'147927:1-147927:152',
	'147929:1-147929:266',
	'147929:272-147929:618',
	'147929:620-147929:643',
	'148002:92-148002:203',
	'148029:50-148029:483',
	'148029:485-148029:569',
	'148029:571-148029:571',
	'148031:1-148031:341',
	'148031:472-148031:757',
	'148031:759-148031:855',
	'148032:1-148032:199',
	'148058:1-148058:97',
	'148822:1-148822:446',
	'148829:1-148829:73',
	'148829:75-148829:240',
	'148829:244-148829:303',
	'148860:1-148860:39',
	'148862:1-148862:18',
	'148862:20-148862:108',
	'148862:110-148862:149',
	'148862:151-148862:165',
	'148862:224-148862:258',
	'148862:262-148862:297',
	'148862:299-148862:366',
	'148862:368-148862:504',
	'148862:512-148862:679',
	'148864:1-148864:31',
	'148864:33-148864:141',
	'148864:224-148864:236',
	'148864:238-148864:476',
	'148864:478-148864:680',
	'148952:70-148952:257',
	'148953:1-148953:100',
	'149003:84-149003:238',
	'149011:1-149011:341',
	'149011:343-149011:706',
	'149058:1-149058:65',
	'149063:1-149063:102',
	'149181:229-149181:1840',
	'149181:1844-149181:1920',
	'149182:1-149182:16',
	'149182:18-149182:62',
	'149182:84-149182:169',
	'149182:171-149182:444',
	'149291:79-149291:79',
	'149291:82-149291:786',
	'149291:788-149291:788',
	'149291:790-149291:790',
	'149291:794-149291:794',
	'149294:1-149294:171',
)
