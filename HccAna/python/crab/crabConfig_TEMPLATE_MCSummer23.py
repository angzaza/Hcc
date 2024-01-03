from CRABClient.UserUtilities import config, getUsername
config = config()

config.General.requestName = 'OUTFILENAME'
config.General.workArea = 'resultsAna_JOBTAG/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/templateMC_QCDZqq_Summer23.py'
#config.JobType.outputFiles = ['OUTFILENAME.root']
#config.JobType.scriptExe = 'submitFileCrab.sh'
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = 'DATASETNAME'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.JobType.numCores = 2
config.Data.publication = True
# This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'QCDpreEE_PT120to170_ntuple'

config.Site.storageSite = 'T2_IT_Bari'
#config.JobType.allowUndistributedCMSSW = True
#config.Site.ignoreGlobalBlacklist  = True
