import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.selections = cms.PSet(

   nEvents   =  cms.untracked.int32(-1),
   MCpresent = cms.untracked.bool(True),

   maxJetPt  = cms.untracked.double(100.),

   minEleEt  = cms.untracked.double(30.),
   eleId     =  cms.untracked.int32(0),  #hardcoded!!

   minEtOverMet = cms.untracked.double(0.4),
   maxEtOverMet = cms.untracked.double(1.5),

   eleMetPhiMin =  cms.untracked.double(2.5),

   outFile = cms.string("PROVA/DATA_ALL/WPrimeAnalysisTree.root"),

 
        )


process.inputNtuples = cms.PSet(
    inputFiles = cms.vstring(

#        '/media/amassiro/deguio/Wprime/DATA_15092010_correctTag/EG_Run2010A-Jul16thReReco-v2_RECO_15092010_correctGlobalTag_merged.root',
#        '/media/amassiro/deguio/Wprime/DATA_15092010_correctTag/EG_Run2010A-PromptReco-v4_15092010_correctGlobalTag_merged.root'

        '/media/amassiro/deguio/Wprime/DATA_30092010/EG_Run2010A-Sep17ReReco_v2_RECO_merged.root',
        '/media/amassiro/deguio/Wprime/DATA_30092010/EG_Run2010A-Sep17ReReco_v2_Jul16_merged.root'

    )
)
