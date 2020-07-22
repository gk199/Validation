# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: l1Ntuple -s RAW2DIGI --python_filename=ntuple_maker_def.py -n 10 --no_output --era=Run3 --data --conditions=106X_upgrade2021_realistic_v4 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimHcalTP --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleGEN --customise_commands=process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False) --filein=file:/eos/cms/store/user/lowang/mh1000_pl10000_step1.root --no_exec
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
# files from DAS with TDC - June 2020
# Higgs -> XX -> bbbb, mh=125, mff=50, ctau=3000mm
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/1324D824-A31D-A240-BA0E-63AC08699698.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/692C5437-303C-CF45-9E1C-9BB93D7437FB.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/D5097967-94E0-3246-97F2-CC45B1A6D282.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/F39A8599-916F-FA44-A3EE-370A623E3226.root'), # only non-TAPE files
# Higgs -> XX -> bbbb, mh=125, mff=50, ctau=30000mm TAPE (all)
# Higgs -> XX -> bbbb, mh=350, mff=160, ctau=1000mm TAPE (mostly)
# Higgs -> XX -> bbbb, mh=250, mff=120, ctau=1000mm
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/F4C824A6-CC41-B24A-AE34-3752F607F909.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/EC854807-A260-4940-8860-DA7BE9113291.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/E9ECD2BD-FE90-3942-B0FA-A170D92E7AD7.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-250_MFF-120_CTau-1000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/DDB82C98-2E22-D743-8C8F-B0BF6DED54A0.root'),
# Higgs -> XX -> bbbb, mh=1000, mff=450, ctau=10000mm
                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/CB5F0BDD-1F67-744B-917C-627523A7AB37.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/C61600E0-BD7B-F04F-A6BC-21C6C7255978.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/A5F6039D-BA76-DA45-A535-DE9D2C458B1A.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/230000/9C995599-5833-1A43-9249-1D277EF03D09.root'),
# QCD, flat pt 15-3000
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240002/8A122C34-2DA8-9B45-8442-59A906FC2D33.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240002/89DA7E24-BB5B-4A4A-9DDD-7963A4ACBB14.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240002/89C323F3-88FC-CA4C-80E2-B9C4CC69E717.root','root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240002/8878ADE4-057D-F74A-AA21-DD5A196C13F3.root'),
# neutrino gun TAPE
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/Neutrino_Pt-2to20_gun_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240000/FFBB7939-1DFC-4948-9BEB-4159B15714C3.root', 'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/Neutrino_Pt-2to20_gun_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240000/FD24971A-F204-AF40-864C-3C9E769EF9A4.root', 'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/Neutrino_Pt-2to20_gun_HCAL/GEN-SIM-DIGI-RAW/packHBTDC_110X_mcRun3_2021_realistic_v6-v1/240000/FC81B680-E617-104C-9F32-9BA7C5A3C756.root'),
# min bias, tune CP5
#                            fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20GS/MinBias_TuneCP5_14TeV-pythia8/GEN-SIM/110X_mcRun3_2021_realistic_v6-v1/00009/E66C34AE-ECB4-5A4F-9A27-8CEE7404C55A.root'),
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
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')

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
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleGEN

#call to customisation function L1NtupleRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleRAWEMU(process)
process = L1NtupleGEN(process) 

# End of customisation functions

# Customisation from command line

process.HcalTPGCoderULUT.LUTGenerationMode=cms.bool(False)
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
