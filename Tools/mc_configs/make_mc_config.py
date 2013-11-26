CONFIG_TEXT = """
## WR -> Nu + l -> j + j + l + l


import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.),
    crossSection = cms.untracked.double(0.0343),
    ExternalDecays = cms.PSet(
        Tauola = cms.untracked.PSet(
            TauolaPolar,
            TauolaDefaultInputCards
        ),
        parameterSets = cms.vstring('Tauola')
    ),
    UseExternalGenerators = cms.untracked.bool(True),        
    comEnergy = cms.double(8000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring(
            'MSEL = 0      ! User defined process',
            'MSUB(354) = 1 ! Heavy neutrino',
            'MSTJ(41) = 2  ! Pythia QED bremsshtrahlung',
            'PMAS(C9900024,1) = %(WRMass)i. ! WR Mass',
            'PMAS(C9900012,1) = %(NuMass)i. ! nu_Re',
            'PMAS(C9900014,1) = %(NuMass)i. ! nu_Rmu',
            'PMAS(C9900016,1) = %(NuMass)i. ! nu_Rtau',
            '9900024:alloff',
            '9900024:onifany = 9900012 9900014 9900016'),
        parameterSets = cms.vstring('pythiaUESettings', 'processParameters')
    )
)

ProductionFilterSequence = cms.Sequence(generator)
"""

WR_MASSES = (1000, 1500, 2000, 3000)
FACTORS = (1/16., 1/8., 3/16., 1/4., 5/6.)
for wr_mass in WR_MASSES:
    for factor in FACTORS:
        nu_mass = wr_mass * factor
        str_dict = {'WRMass':wr_mass, 'NuMass':nu_mass}
        file_name = "WRToNuLeptonToLLJJ_MWR_%(WRMass)i_MNu_%(NuMass)i_TuneZ2star_8TeV_pythia6_tauola_cff.py" % str_dict
        text = CONFIG_TEXT % str_dict
        f = open(file_name, 'w')
        f.write(text)
        f.close()
