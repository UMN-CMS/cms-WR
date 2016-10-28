import FWCore.ParameterSet.Config as cms
import pickle
import math

Additional_ID_Systematics = 0.01
Additional_ISO_Systematics = 0.01

f = open('python/MuonID_Z_RunBCD_prompt80X_7p65.pkl', 'r')
g = open('python/MuonIso_Z_RunBCD_prompt80X_7p65.pkl', 'r')
ID_results = pickle.load(f)
Iso_results = pickle.load(g)

ii = 0
SF_ID_C = []
SF_ID_E = []
SF_ISO_C = []
SF_ISO_E = []

for i in xrange(0,13):
    SF_ID_C.append(0)
    SF_ID_E.append(0)
    SF_ISO_C.append(0)
    SF_ISO_E.append(0)

for key, result in sorted(ID_results["MC_NUM_HighPtIDPt20andIPCut_DEN_genTracks_PAR_eta"]["eta_ratio"].iteritems()) : 
    SF_ID_C[ii] = result["value"]
    SF_ID_E[ii] = math.sqrt(pow(result["error"],2) + pow(Additional_ID_Systematics,2))
    ii += 1

ii=0
for key, result in sorted(Iso_results["MC_NUM_LooseRelTkIso_DEN_HighPtID_PAR_eta"]["eta_ratio"].iteritems()) :
    SF_ISO_C[ii] = result["value"]
    SF_ISO_E[ii] = math.sqrt(pow(result["error"],2) + pow(Additional_ISO_Systematics,2))
    ii += 1

    

# make a collection of TuneP muons which pass isHighPt ID
muonIdIsoSF = cms.EDProducer("produceLepIdIsoScaleFactors",
		src = cms.InputTag("wRminiTreeMuon"),
                OutputCollectionName1 = cms.string("MuonSFIdCentral"),
                OutputCollectionName2 = cms.string("MuonSFIdError"),
                OutputCollectionName3 = cms.string("MuonSFIsoCentral"),
                OutputCollectionName4 = cms.string("MuonSFIsoError"),
                Scale_Factor_ID_Central = cms.vdouble(SF_ID_C),
                Scale_Factor_ID_Error = cms.vdouble(SF_ID_E),
                Scale_Factor_ISO_Central = cms.vdouble(SF_ISO_C),
                Scale_Factor_ISO_Error = cms.vdouble(SF_ISO_E)
)

#MuonIdIsoSFProdSequence = cms.Sequence(MuonIdIsoSFProd)
