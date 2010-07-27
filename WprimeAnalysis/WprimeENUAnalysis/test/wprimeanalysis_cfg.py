import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START3X_V26A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
   '/store/relval/CMSSW_3_5_7/RelValWE/GEN-SIM-RECO/START3X_V26-v1/0012/EC9F278B-6949-DF11-99F2-003048678B0E.root'
                                                )
)

#Define PAT sequence
# Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff");


# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')


# get the jet corrections
##from PhysicsTools.PatAlgos.tools.jetTools import *
##switchJECSet( process, "Summer09_7TeV_ReReco332")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
## uncomment this line to run on an 35X input sample
run36xOn35xInput(process)

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
   eleLabel   = cms.InputTag("patElectrons"),
   barrelCuts = cms.PSet(heepBarrelCuts),
   endcapCuts = cms.PSet(heepEndcapCuts)
)


#Analysis
process.myanalysis = cms.EDAnalyzer('WprimeTree',
                                    
   recHitCollection_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
   recHitCollection_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
   superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
   superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

   electronTag         = cms.InputTag("heepPatElectrons"),
   jetTag              = cms.InputTag("patJets"),
   calometTag          = cms.InputTag("patMETs"),
   tcmetTag            = cms.InputTag("patMETsTC"),
   pfmetTag            = cms.InputTag("patMETsPF"),
   muonTag             = cms.InputTag("patMuons"),
   electronID          = cms.untracked.string("eidRobustHighEnergy"),
   btagAlgo            = cms.untracked.string("jetBProbabilityBJetTags"),
   HLTInputTag         = cms.InputTag("TriggerResults::HLT"),
   L1InputTag          = cms.InputTag("gtDigis"),

   runOnMC             = cms.bool(True),                     
   storePDFWeights     = cms.bool(False),
   pdfWeightsTag       = cms.InputTag("pdfWeights:cteq65")                            
)


process.TFileService = cms.Service("TFileService",
   fileName = cms.string("test.root")
)


process.p = cms.Path(
   process.patDefaultSequence *
   process.heepPatElectrons* #heepifies the pat electrons (resets energy to ecal energy and adds heep id)
   process.myanalysis
   ) 



## process.load("Configuration.EventContent.EventContent_cff")
## process.out = cms.OutputModule("PoolOutputModule",
##             process.FEVTSIMEventContent,
##             fileName = cms.untracked.string('file:/tmp/malberti/testPat.root')
## )

## process.out.outputCommands.append('drop *_*_*_*')
## process.out.outputCommands.append('keep *_*_*_myprocess')
## process.out.outputCommands.append('keep *_*_*_HLT')

## process.outpath = cms.EndPath(process.wprimePATSequence *  process.out)

