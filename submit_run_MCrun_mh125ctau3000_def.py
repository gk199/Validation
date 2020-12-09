from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'mh125_ctau3000_TDC'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuple_maker_def.py'

config.Data.inputDataset = '/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm_TuneCP5_14TeV_pythia8_HCAL/Run3Winter20DRPremixMiniAOD-packHBTDC_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_mh125_ctau3000_TDC_MC'

# Where the output files will be transmitted to                                                                                                                                                   
config.Data.outLFNDirBase = '/store/group/dpg_hcal/comm_hcal/gillian/LLP_Run3/TimingTrigger/Dec110X/'
config.Site.storageSite = 'T2_CH_CERN'



# info from template file
#NEWCONDITIONS = False
#OUTPUTSITE = 'T2_CH_CERN'
#DATASET = '/Neutrino_Pt-2to20_gun_HCAL/Run3Winter20DRPremixMiniAOD-packHBTDC_110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-DIGI-RAW'

