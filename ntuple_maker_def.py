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
#   fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-1000_MFF-450_CTau-500mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh1000_pl500_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-1000_MFF-450_CTau-1000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh1000_pl1000_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-1000_MFF-450_CTau-10000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh1000_pl10000_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/mh1000_mh450_pl500_step1.root'), # this is the no pileup sample!!      
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-350_MFF-180_CTau-500mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh350_ml160_pl500_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-350_MFF-180_CTau-1000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh350_ml160_pl1000_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-350_MFF-180_CTau-10000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh350_ml160_pl10000_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-250_MFF-120_CTau-500mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh250_ml120_pl500_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-250_MFF-120_CTau-1000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh250_ml120_pl1000_step1.root'),     
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-250_MFF-120_CTau-10000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh250_ml120_pl10000_step1.root'),     
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-125_MFF-50_CTau-500mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh125_pl500_step1.root'),
#fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-125_MFF-50_CTau-1000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh125_pl1000_step1.root'),
fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/LLP_htobbbb/step1/MH-125_MFF-50_CTau-10000mm_step1.root'), #'file:/eos/cms/store/user/lowang/mh125_pl10000_step1.root'),
#                            fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_1.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_3.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_4.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_5.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_6.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_7.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_8.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_9.root','file:/eos/cms/store/user/lowang/RelValNuGun/RelValNuGun_RAW_TDC/200129_165224/0000/RelValNuGun_PU_step1_10.root'),
#    fileNames = cms.untracked.vstring('file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_1.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_2.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_3.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_4.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_5.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_6.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_7.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_8.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_9.root','file:/eos/cms/store/user/lowang/SMP-PhaseIITDRFall17DR-00002_step1_10.root'), # QCD files from Long
# files from DAS with TDC - March 2020
# mh 350 GeV mX 160 GeV ct 500mm = 0.5m
#fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/50000/99864D7A-E892-1549-8968-17C448835ECC.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/50000/8F1C77A9-60C7-864A-B8A5-76E14C6FB5CB.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/FC790DA9-A472-E54D-8F50-298BA27C2075.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/DDCC69F6-F53B-4B43-BE3D-9CFD20D0DCBC.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/D8C33D7B-6AE7-2D45-A521-D8BA8FA10F19.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/CE74CFCC-24F4-DD46-A0D2-18BD723B5BD6.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/C259CC1E-8357-7F49-8629-541C80B60DA1.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/B97E0570-E868-A742-B18E-5DD9439A35CB.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/B14E8173-AEF0-444C-940D-961D72FFFC67.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/99F09F37-DFF3-B24A-AC21-D21C85F064A1.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/992A6601-AD30-B647-B037-4D283EB226DB.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/94B2E145-DF30-6F49-B6F0-CC18F838F4B5.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/86FCA72A-047C-5840-B2AD-E191D5BDDE72.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/85EC4665-1513-FA47-B1C2-53EAD918CFA3.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/23501A15-31CA-6C4C-B115-B232C686B96A.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/186FF856-106E-E341-944F-BCDA00840E27.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/16233FE0-F658-D348-8799-A30724352E73.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/15980B93-8FA0-F34E-858D-A20C1EF57B02.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/09C1CC74-385D-794F-98F7-A49C3EB43ACA.root',
#                                  'root://cms-xrd-global.cern.ch///store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-350_MFF-160_CTau-500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/10000/046D0D55-0920-EF4A-84BD-7640FABCE8E6.root'),
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
