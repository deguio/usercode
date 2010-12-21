import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.selections = cms.PSet(

   nEvents   =  cms.untracked.int32(-1),
   MCpresent = cms.untracked.bool(True),
   nonIso = cms.untracked.bool(False),

   maxJetPt  = cms.untracked.double(100.),

   minEleEt  = cms.untracked.double(30.),
   eleId     =  cms.untracked.int32(0),  #hardcoded!!

   minEtOverMet = cms.untracked.double(0.4),
   maxEtOverMet = cms.untracked.double(1.5),

   eleMetPhiMin =  cms.untracked.double(2.5),

   evtList = cms.string("EVTLIST"),
   outFile = cms.string("OUTFILE"),
        )


process.inputNtuples = cms.PSet(
    inputFiles = cms.vstring(

        'LISTOFFILES'
    )
)
