from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'OUTFILENAME'
config.General.workArea = 'resultsAna_JOBTAG/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/templateData_Prompt2023Cv4_JECdB.py'
#config.JobType.outputFiles = ['OUTFILENAME.root']
#config.JobType.scriptExe = 'submitFileCrab.sh'
config.JobType.inputFiles  = ['/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/Summer23Prompt23_RunCv4_V1_DATA.db', '/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/Summer23Prompt23_RunCv4_JRV1_DATA.db', '/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/Summer23Prompt23_V1_MC_UncertaintySources_AK4PFPuppi.txt']
config.Data.inputDBS = 'global'
config.Data.inputDataset = 'DATASETNAME'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.lumiMask = '/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/Utilities/crab/Cert_Collisions2023_366442_370790_Golden.json'
config.JobType.numCores = 2
config.Data.publication = True
# This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'QCDpreEE_PT120to170_ntuple'

config.Site.storageSite = 'T2_IT_Bari'
#config.JobType.allowUndistributedCMSSW = True
#config.Site.ignoreGlobalBlacklist  = True
