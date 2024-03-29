
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('START3X_V26A::All')  #35x
process.GlobalTag.globaltag = cms.string('START36_V10::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(106)
    )

process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring(
        #'file:/media/amassiro/deguio/Datasets/Summer10_Wenu_GEN-SIM-RECO_START36_V9_S09-v1.root'
        'file:/media/amassiro/deguio/Datasets/Fall10_WToENu_TuneZ2_7TeV-pythia6_AODSIM_START38_V12-v1.root'
        #'file:/media/amassiro/deguio/Datasets/Spring10_WprimeToENu_M-800_7TeV-pythia6_GEN-SIM-RECO_START3X_V26-v1.root'

       
       )
)

#filter to divide overlapped pthat for QCD MC
process.genFilter = cms.EDFilter("MCProcessFilter",
                                 MinPthat = cms.untracked.vdouble(30),
                                 MaxPthat = cms.untracked.vdouble(80)
                                 )

#PAT configuration
process.load("PhysicsTools.PatAlgos.patSequences_cff");

from PhysicsTools.PatAlgos.tools.coreTools import *

# turn off MC matching for the process
#removeMCMatching(process, ['All'])


# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *

## uncomment this line to run on an 35X input sample
#run36xOn35xInput(process)

#Analysis

# FILTERS 
process.highetele = cms.EDFilter("GsfElectronSelector",
        src = cms.InputTag("gsfElectrons"),
        cut = cms.string("superCluster().get().energy()*sin(theta())> 20")
)


process.highetFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("highetele"),
                                    minNumber = cms.uint32(1),
                                    )

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
   btagAlgoHighEff     = cms.untracked.string("simpleSecondaryVertexHighEffBJetTags"),
   btagAlgoHighPur     = cms.untracked.string("simpleSecondaryVertexHighPurBJetTags"),
   HLTInputTag         = cms.InputTag("TriggerResults::REDIGI36X"),
   L1InputTag          = cms.InputTag("gtDigis"),

   runOnMC             = cms.bool(True),
   storePDFWeights     = cms.bool(False),
   pdfWeightsTag       = cms.InputTag("pdfWeights:cteq65")

)

#save HLT infos
from PhysicsTools.NtupleUtils.HLTrigResultsDumper_cfi import *
process.TriggerResults = HLTrigResultsDumper.clone()
process.TriggerResults.HLTriggerResults = cms.InputTag("TriggerResults::REDIGI36X")
process.TriggerResults.HLTPaths = cms.vstring('HLT_Photon10_L1R','HLT_Ele10_LW_L1R','HLT_Ele15_LW_L1R','HLT_Ele15_SW_L1R','HLT_Ele15_SW_CaloEleId_L1R','HLT_Ele17_SW_CaloEleId_L1R','HLT_Ele22_SW_CaloEleId_L1R')   # provide list of HLT paths (or patterns) you want

# filter on primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
   minimumNDOF = cms.uint32(4) ,
   maxAbsZ = cms.double(24),
   maxd0 = cms.double(2)
)


process.eventsCounterTotal = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterTotal"))
process.eventsCounterGoodEvt = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterGoodEvt"))
process.eventsCounterHighEtEle = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterHighEtEle"))
process.eventsCounterPatElectronSequence = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterPatElectronSequence"))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("WPrimeAnalysisTree_MC.root")
)


# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)



process.p = cms.Path(

    process.eventsCounterTotal *  #<<---
    
#    process.genFilter *
    process.noscraping *
    process.primaryVertexFilter *



    process.eventsCounterGoodEvt *  #<<---

    process.highetele *
    process.highetFilter *

    process.eventsCounterHighEtEle *  #<<---

    process.patDefaultSequence *
    process.heepPatElectrons *

    process.eventsCounterPatElectronSequence *  #<<---

    process.TriggerResults *
    
    process.myanalysis
    ) 



