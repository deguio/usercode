import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)


process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

    #prodotto io
    #'file:/data/deguio/Dalitz/PYTHIA6_QCDpt_80_120_10TeV_cff_py_GEN.root'

    #/QCDDiJetPt80to120/Summer08_IDEAL_V9_v1/GEN-SIM-RECO
    'file:/data/deguio/Dalitz/QCDDiJetPt80to120/009AC3E3-BF97-DD11-93B5-00093D13BB43.root',
    'file:/data/deguio/Dalitz/QCDDiJetPt80to120/0E2CC0F1-BF97-DD11-B4F6-0015C5EC47A2.root',
    'file:/data/deguio/Dalitz/QCDDiJetPt80to120/122F43EF-BF97-DD11-822B-0015C5E5B9C5.root',
    'file:/data/deguio/Dalitz/QCDDiJetPt80to120/2A56CC78-2397-DD11-A133-0015C5E5B9C5.root'

    )
)


process.myanalysis = cms.EDAnalyzer("QCDAnalyzer",
        HepMCLabel     = cms.untracked.string("source")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("testQCD_80_120.root")
)

process.p = cms.Path(process.myanalysis)
