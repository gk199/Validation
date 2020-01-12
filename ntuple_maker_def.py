# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=ntuple_maker_def.py -n 10 --no_output --era=Run3 --data --conditions=106X_upgrade2021_realistic_v4 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/A9F2F183-B6CF-334E-B11B-9555E62B06C6.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RAW2DIGI',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/1C14A4C8-A1F9-EE43-93EB-EAD32B0EDF0F.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/A9F2F183-B6CF-334E-B11B-9555E62B06C6.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/DC6047E8-D25F-8445-B716-8919C8851935.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/81D8B636-6C78-D44A-837A-8EFCC0C0FD4D.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/DD24FA22-51E6-2C45-B9C3-DFA904CB470B.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/F0C8B255-A341-654D-855E-703954F3D2C4.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/D5D7F509-64B0-084A-A72A-43C3C1B8409C.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/F630853D-23E7-9342-A97C-B769337A5A07.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/780573E2-A564-144D-9041-F0CC02736B94.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/AAB27338-2C04-7246-A061-1891589E70A2.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/D1ABAAF2-91BB-E146-BA8E-A1A0C9E130E4.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/51A22145-1CA4-3F44-A8DE-9399A2BAF2ED.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/EF7E9E51-A8C2-F84A-9806-06ED18AFF56E.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/3973AA2D-8B74-9C46-B40E-264AAE8A6BD6.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/C2065E51-152D-6748-A61A-9DA36E270F17.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/1200E14F-B3A1-594B-9600-5F1F07E12687.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/34FE9804-9E8B-A448-9D3F-4AADD2204A2F.root',
                                                              '/store/relval/CMSSW_10_6_1_patch1/RelValNuGun/GEN-SIM-DIGI-RAW/PU_106X_mcRun3_2021_realistic_v3_rsb-v2/20000/ECC183E8-FD9C-984E-99A8-2E50DB284414.root'
),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2021_realistic_v4', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAWsimHcalTP 

#call to customisation function L1TReEmulFromRAWsimHcalTP imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAWsimHcalTP(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleRAWEMU 

#call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process)

# End of customisation functions

# Customisation from command line

process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
