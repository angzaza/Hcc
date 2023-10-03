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

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
'/store/data/Run2023C/JetMET0/MINIAOD/PromptReco-v2/000/367/516/00000/0a43c8b3-0528-4157-b6eb-e93f3f45da2b.root',
'/store/data/Run2023C/JetMET0/MINIAOD/PromptReco-v2/000/367/516/00000/0a71927c-d393-467e-85b7-160f0fc95d5c.root',
)


process.source = cms.Source("PoolSource",fileNames = myfilelist,
           #                 duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
           #eventsToProcess = cms.untracked.VEventRange('1:12:66017')
                            )

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("prova.root',")
                                   fileName = cms.string("out_Data2023.root")
)

# jet energy corrections

#import os
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jec = cms.ESSource("PoolDBESSource",
#                           CondDBSetup,
#                           connect = cms.string("sqlite_file:/afs/cern.ch/work/a/azaza/HccAna/prova/CMSSW_13_1_0_pre2/src/Hcc/HccAna/python/Winter23Prompt23_RunA_V1_DATA.db"),
#                           toGet =  cms.VPSet(
#                              cms.PSet(
#                                 record = cms.string("JetCorrectionsRecord"),
#                                 tag = cms.string("JetCorrectorParametersCollection_Winter23Prompt23_RunA_V1_DATA_AK4PFPuppi"),
#                                 label= cms.untracked.string("AK4PFPuppi")
#                              ),
#              )
#)



#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#updateJetCollection(
#   process,
#   jetSource = cms.InputTag('slimmedJetsPuppi'),
#   labelName = 'UpdatedJEC',
#   jetCorrections = ('AK4PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)

#)
#from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask

#process.jecSequence = cms.Sequence(process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        #initialSeed = cms.untracked.uint32(SEED), # for HPC
        initialSeed = cms.untracked.uint32(123456), # for crab
        engineName = cms.untracked.string('TRandom3')
    )
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
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJets"),
                              #AK4PuppiJetSrc       = cms.InputTag("updatedPatJetsUpdatedJEC"),
                              AK4PuppiJetSrc       = cms.InputTag("slimmedJetsPuppi"),
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
                              mergedjetSrc = cms.untracked.InputTag("slimmedJets"),
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
                              #elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              #muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              #rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
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
                                  'HLT_DiPFJetAve80_v',
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
                              year = cms.untracked.int32(2018),####for year put 2016,2017, or 2018 to select correct setting

payload = cms.string("AK4PFchs"),


                             )


process.p = cms.Path(#process.fsrPhotonSequence*
                     #process.boostedMuons*
                     #process.jecSequence*
                     #process.calibratedMuons*
                     #process.regressionApplication*
                     #process.selectedElectrons*
                     process.unpackedPatTrigger* 
                     #process.calibratedPatElectrons*
                     #process.egmGsfElectronIDSequence*
                     #process.electronMVAValueMapProducer*
                     #process.electronsMVA*
                     #process.egmPhotonIDSequence*
                     #process.egammaPostRecoSeq*
 	             #process.calibratedPatElectrons*
                     #process.jetCorrFactors*
                     #process.pileupJetIdUpdated*
                     #process.slimmedJetsJEC*
                     #process.QGTagger*
                     #process.AK8PFJetCorrFactors*
                     #process.slimmedJetsAK8JEC*
                     #process.fullPatMetSequence*
                     #process.corrJets*
                     #process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS*#process.rivetProducerHZZFid*
		     #process.prefiringweight *
                   #  process.printTree*
                     process.Ana
                     
               )

