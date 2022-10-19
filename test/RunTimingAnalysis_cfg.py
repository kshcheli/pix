
import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('PPSTiming2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
                        "/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/357/482/00000/c0f85b9d-9108-404a-8031-a1236afed164.root"
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/807/00000/05c6bf1e-0ad0-4750-aefe-db84ab2cc55d.root"
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/087e7349-282f-4372-8182-490a61c3931e.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/10a7ac22-984a-440a-b8bb-fb0262f0ae88.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/6a28ffe0-9d45-4e61-b37b-3f60ea36a6f9.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/6da027ca-dc4d-4016-b273-507ada786ab1.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/8c4496ea-2836-46cd-8bdd-294b2731a301.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/8d8381cb-62fa-4479-98d1-94873083115f.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/94297045-edc5-4fa6-ad69-844639d47583.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/d589a70f-48a3-4e2e-9479-57c4c951e85e.root",
			#"/store/data/Run2022C/ZeroBias/AOD/PromptReco-v1/000/355/933/00000/faa1abb5-0c17-4e43-b511-88227fc9daf0.root"
                         )
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import GlobalTag
process.GlobalTag.globaltag = '124X_dataRun3_Prompt_frozen_v4'
process.GlobalTag.toGet = cms.VPSet()

# JH: recipe from C. Misan to pick up new timing calibrations                                                                                                                           
process.GlobalTag.toGet=cms.VPSet(
  cms.PSet(record = cms.string("PPSTimingCalibrationRcd"),
           tag = cms.string("PPSDiamondTimingCalibration_Run3_recovered_v1"),
           label = cms.untracked.string('PPSTestCalibration'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
          )
)

# JH - rerun reco sequence with new timing conditions                                                                                                                   
process.load("RecoPPS.Configuration.recoCTPPS_cff")
process.ctppsDiamondRecHits.timingCalibrationTag=cms.string("GlobalTag:PPSTestCalibration")

process.ctppsDiamondLocalTracks.recHitsTag = cms.InputTag("ctppsDiamondRecHits","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.tagDiamondTrack = cms.InputTag("ctppsDiamondLocalTracks","","PPSTiming2")
process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","PPSTiming2")
process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(True)
process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)


process.mydiamonds = cms.EDAnalyzer(
    'PPSTimingAnalyzer',
    #    lhcInfoLabel = cms.string(''),
    verticesTag = cms.InputTag('offlinePrimaryVertices'),
    tracksTag = cms.InputTag('generalTracks'),
    #
    # Take PPS information from the existing Prompt RECO AOD
    #
    #    tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits"),
    #    tagTrackLites = cms.InputTag( "ctppsLocalTrackLiteProducer"),
    #    ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP"),
    #    ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP"),
    #    #
    # Alternatively, uncomment these lines to take PPS information from on-the-fly re-RECO
    #
    tagDiamondRecHits = cms.InputTag("ctppsDiamondRecHits",""),                                               
    tagTrackLites = cms.InputTag( "ctppsLocalTrackLiteProducer", ""), 
    tagRPixRecHit = cms.InputTag("ctppsPixelRecHits", ""), 
    ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP"),
    ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP"),
    maxVertices = cms.uint32(1000),
    outfilename = cms.untracked.string( "output_ZeroBias.root" )
)

# Trigger                                                                                                                  
#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#from HLTrigger.HLTfilters.HLTLevel1Seed_cfi import *

#process.load('L1Trigger.Skimmer.l1Filter_cfi')
#process.l1Filter.algorithms = cms.vstring("L1_FirstCollisionInOrbit")

process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('L1_FirstCollisionInOrbit')


#process.hltFilter = copy.deepcopy(HLTLevel1Seed)
#process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilter.HLTPaths = ['L1_FirstCollisionInOrbit']

process.ALL = cms.Path(
#    process.hltLevel1GTSeed *
  # Uncomment these lines, to re-run the PPS local+proton timing reconstruction starting from AOD
 #   process.ctppsDiamondRecHits *
 #   process.ctppsDiamondLocalTracks *
 #   process.ctppsPixelLocalTracks * 
 #   process.ctppsLocalTrackLiteProducer *
 #   process.ctppsProtons *
   process.mydiamonds 
                       )

process.schedule = cms.Schedule(process.ALL)

#print(process.dumpPython())
