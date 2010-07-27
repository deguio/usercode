import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.selections = cms.PSet(

   nEvents   =  cms.untracked.int32(-1),
   MCpresent = cms.untracked.bool(True),

   maxJetPt  = cms.untracked.double(100.),

   minEleEt  = cms.untracked.double(30.),
   eleId     =  cms.untracked.int32(0),  #hardcoded!!

   minEtOverMet = cms.untracked.double(0.6),
   maxEtOverMet = cms.untracked.double(1.5),

   eleMetPhiMin =  cms.untracked.double(2.5),

   #outFile = cms.string("PROVA/DATA_ALL/WPrimeAnalysisTree.root"),

   #outFile = cms.string("PROVA/Wenu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   #outFile = cms.string("PROVA/Wtaunu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   #outFile = cms.string("PROVA/Zee_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   #outFile = cms.string("PROVA/TTbar_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_EMEnriched_Pt80to170_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_EMEnriched_Pt30to80_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_EMEnriched_Pt20to30_Summer10-START36_V9_S09-v2_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_BCtoE_Pt20to30_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_BCtoE_Pt30to80_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   ##outFile = cms.string("PROVA/QCD_BCtoE_Pt80to170_Summer10-START36_V9_S09-v1_GEN-SIM-RECO/WPrimeAnalysisTree.root"),
   #outFile = cms.string("PROVA/QCD_Pt15_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG/WPrimeAnalysisTree.root"),
   #outFile = cms.string("PROVA/QCD_Pt30_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG/WPrimeAnalysisTree.root"),
   outFile = cms.string("PROVA/QCD_Pt80_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG/WPrimeAnalysisTree.root"),

   #test
   #outFile = cms.string("test.root"),
        )


process.inputNtuples = cms.PSet(
            inputFiles = cms.vstring(

        #'/media/amassiro/deguio/Wprime/DATA_17072010/EG_Run2010A-Jul16thReReco-v2_RECO_17072010_merged.root'

        #'/media/amassiro/deguio/Wprime/MC_17072010/Wenu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        #'/media/amassiro/deguio/Wprime/MC_17072010/Wtaunu_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        #'/media/amassiro/deguio/Wprime/MC_17072010/Zee_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        #'/media/amassiro/deguio/Wprime/MC_17072010/TTbar_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_EMEnriched_Pt80to170_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_EMEnriched_Pt30to80_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_EMEnriched_Pt20to30_Summer10-START36_V9_S09-v2_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_BCtoE_Pt20to30_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_BCtoE_Pt30to80_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        ##'/media/amassiro/deguio/Wprime/MC_17072010/QCD_BCtoE_Pt80to170_Summer10-START36_V9_S09-v1_GEN-SIM-RECO_17072010_merged.root'
        #'/media/amassiro/deguio/Wprime/MC_17072010/QCD_Pt15_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG_17072010_merged.root'
        #'/media/amassiro/deguio/Wprime/MC_17072010/QCD_Pt30_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG_17072010_merged.root'
        '/media/amassiro/deguio/Wprime/MC_17072010/QCD_Pt80_Summer10-START36_V9_S09-v1_GEN-SIM-RECODEBUG_17072010_merged.root'

        #test
        #'/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_6_3/src/WprimeAnalysis/WprimeENUAnalysis/test/CrabWork/WPrimeAnalysisTree_MC.root',
        #'/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_6_3/src/WprimeAnalysis/WprimeENUAnalysis/test/CrabWork/WPrimeAnalysisTree_MC2.root'
    )
)
