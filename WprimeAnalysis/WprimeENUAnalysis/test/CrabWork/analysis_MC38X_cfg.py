import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START38_V12::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1006)
    )

process.source = cms.Source("PoolSource",
       #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),    
       fileNames = cms.untracked.vstring(
        #'file:/media/amassiro/deguio/Datasets/Summer10_Wenu_GEN-SIM-RECO_START36_V9_S09-v1.root'
        'file:/media/amassiro/deguio/Datasets/Fall10_WToENu_TuneZ2_7TeV-pythia6_AODSIM_START38_V12-v1.root'
        #'file:/media/amassiro/deguio/Datasets/Spring10_WprimeToENu_M-800_7TeV-pythia6_GEN-SIM-RECO_START3X_V26-v1.root'
        #'file:/media/amassiro/deguio/Datasets/WJets_36XReReco.root'  #madgraph
       )
)


#PAT configuration
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patJets.addTagInfos     = cms.bool(False)

from PhysicsTools.PatAlgos.tools.coreTools import *

#turn off MC matching for the process
#removeMCMatching(process, ['All'])

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')
process.patMETsAK5Calo = process.patMETs

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *


# FILTERS 
process.genFilter = cms.EDFilter("MCProcessFilter",
                                 MinPthat = cms.untracked.vdouble(30),
                                 MaxPthat = cms.untracked.vdouble(80)
                                 )

process.highetele = cms.EDFilter("GsfElectronSelector",
                                 src = cms.InputTag("gsfElectrons"),
                                 cut = cms.string("superCluster().get().energy()*sin(theta())> 25")
                                 )

process.highetFilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("highetele"),
                                    minNumber = cms.uint32(1)
                                    )

from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel   = cms.InputTag("selectedPatElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts)
                                          )

# filter on primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
   minimumNDOF = cms.uint32(4) ,
   maxAbsZ = cms.double(24),
   maxd0 = cms.double(2)
)

# FilterOutScraping
process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)



#Analysis
process.myanalysis = cms.EDAnalyzer('WprimeTreePAT',
   
   recHitCollection_EB = cms.InputTag("reducedEcalRecHitsEB"),
   recHitCollection_EE = cms.InputTag("reducedEcalRecHitsEE"),
   superClusterCollection_EB = cms.InputTag("correctedHybridSuperClusters"),
   superClusterCollection_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),

   PVTag               = cms.InputTag("offlinePrimaryVerticesWithBS"),
   electronTag         = cms.InputTag("heepPatElectrons"),
   jetTag              = cms.InputTag("selectedPatJets"),
   calometTag          = cms.InputTag("patMETsAK5Calo"),
   tcmetTag            = cms.InputTag("patMETsTC"),
   pfmetTag            = cms.InputTag("patMETsPF"),
   muonTag             = cms.InputTag("selectedPatMuons"),
   electronID          = cms.untracked.string("eidRobustHighEnergy"),
   btagAlgoHighEff     = cms.untracked.string("simpleSecondaryVertexHighEffBJetTags"),
   btagAlgoHighPur     = cms.untracked.string("simpleSecondaryVertexHighPurBJetTags"),
   HLTInputTag         = cms.InputTag("TriggerResults::REDIGI38X"),
   L1InputTag          = cms.InputTag("gtDigis"),

   runOnMC             = cms.bool(True),
   storePDFWeights     = cms.bool(False),
   pdfWeightsTag       = cms.InputTag("pdfWeights:cteq65")

)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("WPrimeAnalysisTree_MC.root")
)

#save HLT infos
from PhysicsTools.NtupleUtils.HLTrigResultsDumper_cfi import *
process.TriggerResults = HLTrigResultsDumper.clone()
process.TriggerResults.HLTriggerResults = cms.InputTag("TriggerResults::REDIGI38X")
process.TriggerResults.HLTPaths = cms.vstring('HLT_Photon10_L1R','HLT_Ele10_LW_L1R','HLT_Ele15_LW_L1R','HLT_Ele15_SW_L1R','HLT_Ele15_SW_CaloEleId_L1R','HLT_Ele17_SW_CaloEleId_L1R','HLT_Ele22_SW_CaloEleId_L1R','HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1')   # provide list of HLT paths (or patterns) you want


#save event number after the selections
process.eventsCounterTotal = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterTotal"))
process.eventsCounterGoodEvt = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterGoodEvt"))
process.eventsCounterHighEtEle = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterHighEtEle"))
process.eventsCounterPatElectronSequence = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterPatElectronSequence"))



process.p = cms.Path(

    process.eventsCounterTotal *  #<<--- keep track of the number of event it run on
    
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



