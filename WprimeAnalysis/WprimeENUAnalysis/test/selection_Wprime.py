import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.selections = cms.PSet(
    
   nEvents   =  cms.untracked.int32(-1),
   MCpresent = cms.untracked.bool(False),
   nonIso = cms.untracked.bool(False),
    
   maxJetPt  = cms.untracked.double(100.),
    
   minEleEt  = cms.untracked.double(30.),
   eleId     =  cms.untracked.int32(0),  #hardcoded!!
    
   minEtOverMet = cms.untracked.double(0.4),
   maxEtOverMet = cms.untracked.double(1.5),

   eleMetPhiMin =  cms.untracked.double(2.5),

   #CAMBIARE ANCHE MCPRESENT
    evtList = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_ALL/WPrimeEventList"),
    outFile = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_ALL/WPrimeAnalysisTree.root"),

    #evtList = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_2210/WPrimeEventList"),
    #outFile = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_2210/WPrimeAnalysisTree.root"),
    #evtList = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_1510/WPrimeEventList"),
    #outFile = cms.string("/media/amassiro/deguio/Wprime/PROVA/DATA_1510/WPrimeAnalysisTree.root"),

#   evtList = cms.string("PROVA/DATA_PROMPTRECO/WPrimeEventList"),
#   outFile = cms.string("PROVA/DATA_PROMPTRECO/WPrimeAnalysisTree.root"),

#   evtList = cms.string("PROVA/Wenu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_noScraping/WPrimeEventList"),
#   outFile = cms.string("PROVA/Wenu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_noScraping/WPrimeAnalysisTree.root"),

   #CAMBIARE ANCHE MCPRESENT
#   evtList = cms.string("PROVA/WToENu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1_AODSIM/WPrimeEventList"),
#   outFile = cms.string("PROVA/WToENu_TuneZ2_7TeV-pythia6_Fall10-START38_V12-v1_AODSIM/WPrimeAnalysisTree.root"),

#   evtList = cms.string("PROVA/WJets-madgraph_Summer10-START37_V5_S09-v1_GEN-SIM-RECO/WPrimeEventList"),
#   outFile = cms.string("PROVA/WJets-madgraph_Summer10-START37_V5_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),




        )


process.inputNtuples = cms.PSet(
    inputFiles = cms.vstring(

        '/media/amassiro/Wprime/DATA_08122010/EG_Run2010A-Nov4ReReco_v1_RECO_merged.root',
        '/media/amassiro/Wprime/DATA_08122010/Electron_Run2010B-Nov4ReReco_v1_RECO_merged.root',

#        '/media/amassiro/Wprime/DATA_11112010/EG_Run2010A-Sep17ReReco_v2_RECO_merged.root',
#        '/media/amassiro/Wprime/DATA_11112010/Electron_Run2010B-PromptReco-v2_RECO_merged.root',
        
#        '/media/amassiro/Wprime/DATA_31102010/Photon_Run2010B-PromptReco-v2_RECO_merged.root',

        
    )
)
