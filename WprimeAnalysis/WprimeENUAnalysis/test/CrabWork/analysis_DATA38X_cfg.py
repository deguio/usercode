import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)


process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_38X_V13::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
    )

process.source = cms.Source("PoolSource",
       #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),    
       fileNames = cms.untracked.vstring(
        #'file:/tmp/malberti/reco_7TeV_380_pat_keepAll.root'
        #'file:/tmp/deguio/reco_7TeV_380_pat.root'
        '/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/146/431/B83E5CC1-C0C6-DF11-9834-0030487A3232.root'
       )
)


#PAT configuration
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patJets.addTagInfos     = cms.bool(False)

from PhysicsTools.PatAlgos.tools.coreTools import *

#turn off MC matching for the process
removeMCMatching(process, ['All'])

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')
process.patMETsAK5Calo = process.patMETs

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *


# FILTERS 
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
                                          eleLabel = cms.InputTag("selectedPatElectrons"),
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

# filter on PhysDeclared bit
process.skimming = cms.EDFilter("PhysDecl",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
                                )

# filter on bit = and (40 || 41) and !(bit36 || bit37 || bit38 || bit39)
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')



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
   HLTInputTag         = cms.InputTag("TriggerResults::HLT"),
   L1InputTag          = cms.InputTag("gtDigis"),

   runOnMC             = cms.bool(False),
   storePDFWeights     = cms.bool(False),
   pdfWeightsTag       = cms.InputTag("pdfWeights:cteq65")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("WPrimeAnalysisTree_DATA.root")
)

#save HLT infos
from PhysicsTools.NtupleUtils.HLTrigResultsDumper_cfi import *
process.TriggerResults = HLTrigResultsDumper.clone()
process.TriggerResults.HLTriggerResults = cms.InputTag("TriggerResults::HLT")
process.TriggerResults.HLTPaths = cms.vstring('HLT_Ele10_LW_L1R',
                                              'HLT_Ele15_LW_L1R',
                                              'HLT_Ele15_SW_L1R',
                                              'HLT_Ele15_SW_CaloEleId_L1R',
                                              'HLT_Ele17_SW_CaloEleId_L1R',
                                              'HLT_Ele22_SW_CaloEleId_L1R',
                                              'HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1',
                                              'HLT_Ele22_SW_TighterEleId_L1R_v2',
                                              'HLT_Ele22_SW_TighterEleId_L1R_v3',
                                              'HLT_Photon10_L1R',
                                              'HLT_Photon17_Isol_SC17HE_L1R_v1',
                                              'HLT_Photon20_Cleaned_L1R',
                                              'HLT_Photon22_SC22HE_L1R_v1',
                                              'HLT_Photon25_Cleaned_L1R',
                                              'HLT_Photon30_Cleaned_L1R',
                                              'HLT_Photon35_Isol_Cleaned_L1R',
                                              'HLT_Photon40_CaloId_Cleaned_L1R_v1',
                                              'HLT_Photon40_Isol_Cleaned_L1R_v1',
                                              'HLT_Photon50_Cleaned_L1R',
                                              'HLT_Photon70_Cleaned_L1R_v1',
                                              )   # provide list of HLT paths (or patterns) you want

#save event number after the selections
process.eventsCounterTotal = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterTotal"))
process.eventsCounterGoodEvt = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterGoodEvt"))
process.eventsCounterHighEtEle = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterHighEtEle"))
process.eventsCounterPatElectronSequence = cms.EDFilter("eventsCounter", histoName = cms.string("eventsCounterPatElectronSequence"))


process.p = cms.Path(

    
    process.eventsCounterTotal *  #<<--- keep track of the number of event it run on

    #process.skimming *            #deprecated phys decl
    #process.hltLevel1GTSeed *    #deprecated
    process.noscraping *
    process.primaryVertexFilter *

    process.eventsCounterGoodEvt *  #<<---

    process.highetele *
    process.highetFilter *

    process.eventsCounterHighEtEle *  #<<---

    process.patDefaultSequence *
    process.heepPatElectrons *

    process.eventsCounterPatElectronSequence *  #<<---

    process.TriggerResults *  #save HLT trigger infos
    process.myanalysis
    )



