import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("HccAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.categories.append('HccAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
#process.GlobalTag.globaltag='102X_upgrade2018_realistic_v15'
#process.GlobalTag.globaltag='102X_upgrade2018_realistic_v18'
process.GlobalTag.globaltag='130X_dataRun3_Prompt_v3' 


### uncomment if you apply JECS ###########################
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'Summer23Prompt23_RunCv123_V1_DATA','')
## se vuoi prendere le JEC dalla GlobalTag, devi utilizzare una GT che le contenga
## altrimenti puoi usare una GT che non le contiene e prendere le JEC dal dB file con PoolDBESSource
## in questo caso, commenta le due righe sopra legate alla GT e usa la forma sopra process.GlobalTag.globaltag='126X_mcRun3_2023_forPU65_v4'


process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.options = cms.untracked.PSet(
        numberOfThreads = cms.untracked.uint32(2),
				#SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.options.numberOfConcurrentLuminosityBlocks = 1

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(False),
  printOnlyHardInteraction = cms.untracked.bool(False), # Print only status=3 particles. This will not work for Pythia8, which does not have any such particles.
  src = cms.InputTag("genParticles")
)


myfilelist = cms.untracked.vstring(
#''/store/mc/Run3Summer22MiniAODv3/QCD_PT-120to170_TuneCP5_13p6TeV_pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v2/80000/06ab5fe5-59c7-4dad-9dc2-fe0b75af4379.root',',
'/store/data/Run2023C/JetMET0/MINIAOD/PromptReco-v3/000/367/661/00000/0dacd30d-c8b3-4056-8d17-c6a8fd52d96f.root',
#'/store/data/Run2023C/JetMET0/MINIAOD/PromptReco-v2/000/367/516/00000/0a71927c-d393-467e-85b7-160f0fc95d5c.root',
)


process.source = cms.Source("PoolSource",fileNames = myfilelist,
           #                 duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
           #eventsToProcess = cms.untracked.VEventRange('1:12:66017')
                            )

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("prova.root',")
                                   fileName = cms.string("out_JetMET_Prompt2023Cv123.root")
)

# jet energy corrections

import os
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
                           CondDBSetup,
                           #for execution in local
                           connect = cms.string("sqlite_file:/afs/cern.ch/work/a/azaza/HccAna/CMSSW_13_0_13/src/Hcc/HccAna/python/Summer23Prompt23_RunCv123_V1_DATA.db"),
                           #for ccrab
                           #connect = cms.string("sqlite_file:src/Hcc/HccAna/python/Summer23Prompt23_RunCv123_V1_DATA.db"),
                           toGet =  cms.VPSet(
                              cms.PSet(
                                 record = cms.string("JetCorrectionsRecord"),
                                 tag = cms.string("JetCorrectorParametersCollection_Summer23Prompt23_RunCv123_V1_DATA_AK4PFPuppi"),
                                 label= cms.untracked.string("AK4PFPuppi")
                              ),
              )
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsPuppi'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
)

## la collezione di jet con le correzioni applicate si chiama "updatedPatJetsUpdatedJEC"

process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)


#read JER-SF from a db file
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jer = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        toGet = cms.VPSet(
            # Resolution
            cms.PSet(
                record = cms.string('JetResolutionRcd'),
                tag    = cms.string('JR_Summer23Prompt23_RunCv123_JRV1_DATA_PtResolution_AK4PFPuppi'),
                label  = cms.untracked.string('AK4PFPuppi_pt')
                ),

            # Scale factors
            cms.PSet(
                record = cms.string('JetResolutionScaleFactorRcd'),
                tag    = cms.string('JR_Summer23Prompt23_RunCv123_JRV1_DATA_SF_AK4PFPuppi'),
                label  = cms.untracked.string('AK4PFPuppi')
                ),
            ),
        connect = cms.string('sqlite:Summer23Prompt23_RunCv123_JRV1_DATA.db')
        #for crab
        #connect = cms.string('sqlite:src/Hcc/HccAna/python/Summer23Prompt23_RunCv123_JRV1_DATA.db')

        )

process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

#smearing the jet resolution
process.slimmedJetsSmeared = cms.EDProducer('SmearedPATJetProducer',
       src = cms.InputTag('updatedPatJetsUpdatedJEC'),
       enabled = cms.bool(True),
       rho = cms.InputTag("fixedGridRhoFastjetAll"),
       algo = cms.string('AK4PFPuppi'),
       algopt = cms.string('AK4PFPuppi_pt'),
       #resolutionFile = cms.FileInPath('Autumn18_V7_MC_PtResolution_AK4PFchs.txt'),
       #scaleFactorFile = cms.FileInPath('combined_SFs_uncertSources.txt'),

       genJets = cms.InputTag('slimmedGenJets'),
       dRMax = cms.double(0.2),
       dPtMaxFactor = cms.double(3),

       debug = cms.untracked.bool(False),
   # Systematic variation
   # 0: Nominal
   # -1: -1 sigma (down variation)
   # 1: +1 sigma (up variation)
   variation = cms.int32(0),  # If not specified, default to 0
   uncertaintySource = cms.string(""), # If not specified, default to Total
       )
#process.p=cms.Path(process.slimmedJetsSmeared)
process.jerSequence=cms.Sequence(process.slimmedJetsSmeared)

# example for computing total uncertainties on jet resolution
process.slimmedJetsSmearedDown=process.slimmedJetsSmeared.clone(variation=cms.int32(-1))
#process.p+=process.slimmedJetsSmearedDown
process.slimmedJetsSmearedUp=process.slimmedJetsSmeared.clone(variation=cms.int32(1))
#process.p+=process.slimmedJetsSmearedUp


###########################   PuJetID-Update
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId
process.pileupJetIdUpdated = pileupJetId.clone(
  jets=cms.InputTag("updatedPatJetsUpdatedJEC"),
  inputIsCorrected=True,
  applyJec=False,
  vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
)
#patAlgosToolsTask.add(process.pileupJetIdUpdated)

updateJetCollection(
   process,
   labelName = 'PileupJetID',
   jetSource = cms.InputTag('slimmedJetsSmeared'),
)

process.updatedPatJetsPileupJetID.userData.userInts.src = ['pileupJetIdUpdated:fullId']
process.updatedPatJetsPileupJetID.userData.userFloats.src = ['pileupJetIdUpdated:fullDiscriminant']

###########################   QGTagger
process.load('RecoJets.JetProducers.QGTagger_cfi')
#process.QGTagger.srcJets=cms.InputTag("selectedUpdatedPatJetsPileupJetID")
#process.QGTagger.srcJets=cms.InputTag("updatedPatJetsPileupJetID")
process.QGTagger.srcJets=cms.InputTag("updatedPatJetsUpdatedJEC")
process.QGTagger.srcVertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices")

updateJetCollection(
   process,
   labelName = 'QGT',
   #jetSource = cms.InputTag('updatedPatJetsPileupJetID'),
   jetSource = cms.InputTag('updatedPatJetsUpdatedJEC'),
)
process.updatedPatJetsQGT.userData.userFloats.src = ['QGTagger:qgLikelihood']




# Recompute MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD(process,
            isData=True,
            )


process.unpackedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
    patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
    triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
    unpackFilterLabels = cms.bool(True)
)

# Analyzer
process.Ana = cms.EDAnalyzer('HccAna',
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              #electronUnSSrc  = cms.untracked.InputTag("electronsMVA"),
                              electronUnSSrc  = cms.untracked.InputTag("selectedElectrons"),
                              #electronUnSSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              #electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("slimmedMuons"),
                              #muonSrc      = cms.untracked.InputTag("boostedMuons"),
                              #tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              #jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              #AK4PuppiJetSrc       = cms.InputTag("updatedPatJetsUpdatedJEC"),
                              #AK4PuppiJetSrc       = cms.InputTag("slimmedJetsPuppi"),
                              AK4PuppiJetSrc       = cms.InputTag("updatedPatJetsQGT"),
                              AK4PuppiJetSmearSrc       = cms.InputTag("slimmedJetsSmeared"),
                              AK4PuppiJetSmearUpSrc       = cms.InputTag("slimmedJetsSmearedUp"),
                              AK4PuppiJetSmearDownSrc       = cms.InputTag("slimmedJetsSmearedDown"),
                              AK4PuppiJetUncorrSrc       = cms.InputTag("slimmedJetsPuppi"),
			      AK8PuppiJetSrc       = cms.untracked.InputTag("slimmedJetsAK8"),
                              #hltPFJetForBtagSrc  = cms.InputTag("hltPFJetForBtag", "", "HLT"),
                              #hltAK4PFJetsCorrectedSrc  = cms.InputTag("hltAK4PFJetsCorrected", "", "HLT"),
                              hltAK4CaloJetsCorrectedSrc  = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed", "", "HLT"),
                              hltAK4PFJetsCorrectedSrc  = cms.InputTag("hltAK4PFJetsLooseIDCorrected", "", "HLT"),
                              #pfJetTagCollectionParticleNetprobcSrc = cms.InputTag("hltParticleNetONNXJetTags","probc","HLT"),
                              #pfJetTagCollectionParticleNetprobbSrc = cms.InputTag("hltParticleNetONNXJetTags","probb","HLT"),
                              #pfJetTagCollectionParticleNetprobudsSrc = cms.InputTag("hltParticleNetONNXJetTags","probuds","HLT"),
                              #pfJetTagCollectionParticleNetprobgSrc = cms.InputTag("hltParticleNetONNXJetTags","probg","HLT"),
                              #pfJetTagCollectionParticleNetprobtauhSrc = cms.InputTag("hltParticleNetONNXJetTags","probtauh","HLT"),
                              #jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              #mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              bxvCaloJetSrc =  cms.InputTag("caloStage2Digis","Jet"),
                              bxvCaloMuonSrc =  cms.InputTag("gmtStage2Digis","Muon"),
                              bxvCaloHTSrc =  cms.InputTag("caloStage2Digis","EtSum"),
                              #mergedjetSrc = cms.untracked.InputTag("slimmedJets"),
                              metSrc       = cms.untracked.InputTag("slimmedMETsPuppi"),
                              #metSrc       = cms.untracked.InputTag("slimmedMETs","","Hcc"),
                              #metSrc       = cms.untracked.InputTag("slimmedMETs"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
                              conversionSrc  = cms.untracked.InputTag("reducedEgamma","reducedConversions"),
                              isMC         = cms.untracked.bool(False),
                              isHcc         = cms.untracked.bool(True),
                              isZqq         = cms.untracked.bool(False),
                              isZcc         = cms.untracked.bool(False),
                              isZbb         = cms.untracked.bool(False),
                              isQCD         = cms.untracked.bool(False),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              CrossSection = cms.untracked.double(1),#DUMMYCROSSSECTION),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              pileupSrc     = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              pfCandsSrc   = cms.untracked.InputTag("packedPFCandidates"),
                              fsrPhotonsSrc = cms.untracked.InputTag("boostedFsrPhotons"),
                              prunedgenParticlesSrc = cms.untracked.InputTag("prunedGenParticles"),
                              packedgenParticlesSrc = cms.untracked.InputTag("packedGenParticles"),
                              genJetsSrc = cms.untracked.InputTag("slimmedGenJets"),
                              generatorSrc = cms.untracked.InputTag("generator"),
                              lheInfoSrc = cms.untracked.InputTag("externalLHEProducer"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.InputTag("TriggerResults","","HLT"),
                              #triggerObjects = cms.InputTag("slimmedPatTrigger"),
                              triggerObjects = cms.InputTag("unpackedPatTrigger"),
                              doJER = cms.untracked.bool(True),
                              doJEC = cms.untracked.bool(True),
                              algInputTag = cms.InputTag("gtStage2Digis"),
                              doTriggerMatching = cms.untracked.bool(False),
                              triggerList = cms.untracked.vstring(
                                #VBFHToCC
                                  'HLT_QuadPFJet100_88_70_30_PNet1CvsAll0p5_VBF3Tight_v',
                                  'HLT_QuadPFJet100_88_70_30_v',
                                  'HLT_QuadPFJet105_88_75_30_PNet1CvsAll0p5_VBF3Tight_v',
                                  'HLT_QuadPFJet105_88_75_30_v',
                                  'HLT_QuadPFJet111_90_80_30_PNet1CvsAll0p6_VBF3Tight_v',
                                  'HLT_QuadPFJet111_90_80_30_v',
                                  'HLT_PFJet80_v',
                                  'HLT_PFJet60_v',
                                # Control path VBFHToBB 
                                  'HLT_QuadPFJet103_88_75_15_v',
                                  'HLT_QuadPFJet105_88_76_15_v',
                                  'HLT_QuadPFJet111_90_80_15_v',
                                  'HLT PFJet80_v',
                                  # Single Jet VBFHToBB path:
                                  'HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2_v',
                                  'HLT_QuadPFJet105_88_76_15_PFBTagDeepJet_1p3_VBF2_v',
                                  'HLT_QuadPFJet111_90_80_15_PFBTagDeepJet_1p3_VBF2_v',
                                  # DoubleJet VBFHToBB path
                                  'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v',
                                  'HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v',
                                  'HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepJet_1p3_7p7_VBF1_v',
                                  #AK8 for Zqq
                                  'PFHT1050_v',
                                  'PFJet500_v',
                                  'AK8PFJet500_v',
                                  'AK8PFJet400_TrimMass30_v',
                                  'AK8PFJet420_TrimMass30_v',
                                  'AK8PFHT800_TrimMass50_v',

                              ),
                              verbose = cms.untracked.bool(False),              
                              skimLooseLeptons = cms.untracked.int32(0),              
                              skimTightLeptons = cms.untracked.int32(0),              
                              #bestCandMela = cms.untracked.bool(False),
                              year = cms.untracked.int32(2023),####for year put 2022 use 20220 for preEE and 20221 for postEE
                              isCode4l = cms.untracked.bool(True), 
                              JECUncFileAK4Src = cms.string("Summer23Prompt23_RunCv123_V1_DATA_UncertaintySources_AK4PFPuppi.txt"),
payload = cms.string("AK4PFPuppi"),


                             )


process.p = cms.Path(process.jecSequence* 
                     #process.slimmedJetsSmeared*  --uncomment-- comm
                     #process.pileupJetIdUpdated*  --uncomment-- comm 
                     #process.updatedPatJetsPileupJetID*  --uncomment-- comm
                     process.jerSequence*
                     process.slimmedJetsSmearedDown*
                     process.slimmedJetsSmearedUp*
                     process.QGTagger*
                     process.updatedPatJetsQGT*
                     process.unpackedPatTrigger* 
                     process.Ana
                     )

