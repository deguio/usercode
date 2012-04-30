#ifndef HggVertexAnalysis_h
#define HggVertexAnalysis_h

///==== include ====
#include "h2gglobeEventReader.h"

#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include <time.h>
#include <sstream>
#include "MyTest.h"

#include "VertexAlgoParameters.h"
#include "HggVertexAnalyzer.h"
#include "HggVertexFromConversions.h"
#include "PhotonInfo.h"

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"
#include "Math/Vector4D.h"

class HggVertexAnalysis {
 public :

  HggVertexAnalysis(h2gglobeEventReader*  reader=0);
  ~HggVertexAnalysis();

  void     bookHistos();
  void     saveHistos(TFile *fout);
  void     analyze(int nentries, int isData, int useWeights, TH1F* h, int useKfactors, TH1F* hkfact, int doVBFselection);
  
  int      DiphotonCategory (float eta1, float eta2, float r9_1, float r9_2);

  //void     findMCHiggsPhotons( TClonesArray* gp_p4, Short_t* gp_status, Short_t* gp_pdgid, Short_t* gp_mother, int gp_n, 
  //			       TClonesArray* pho_calopos, TClonesArray* pho_p4, int pho_n,
  //			       int& passSelection, int& mc1, int& mc2 , int& i1, int& i2  );

  void     findMCHiggsPhotons( h2gglobeEventReader* ev_,
  			       int& passSelection, int& mc1, int& mc2 , int& i1, int& i2  );


  double   KfactorsWeight(TH1F* hkfact, TClonesArray* gp_p4, int mc1, int mc2);
 
  std::pair<int, int> Select2HighestPtJets(h2gglobeEventReader* ev_, TLorentzVector& leadpho, TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt);
  bool DijetTag(std::pair<int,int> highestPtJets, TLorentzVector& dipho);

  // histograms 
  
  //histos
  TH1F* PtGen;

  TH1F* PtAll_Baseline ;
  TH1F* PtGood_Baseline;
  TH1F* PtAll_BDT ;
  TH1F* PtGood_BDT;
 
  TH1F* EtaAll_Baseline;
  TH1F* EtaGood_Baseline;
  TH1F* EtaAll_BDT;
  TH1F* EtaGood_BDT;
  
  TH1F* NvtAll_Baseline;
  TH1F* NvtGood_Baseline;
  TH1F* NvtAll_BDT;
  TH1F* NvtGood_BDT;
 
  TH1F* NpuAll_Baseline;
  TH1F* NpuGood_Baseline;
  TH1F* NpuAll_BDT;
  TH1F* NpuGood_BDT;

  TH1F* InvMassAll_Baseline;
  TH1F* InvMassGood_Baseline;
  TH1F* InvMassAll_BDT;
  TH1F* InvMassGood_BDT;

  TH1F* InvMassAll_Baseline_cat[4];
  TH1F* InvMassGood_Baseline_cat[4];
  TH1F* InvMassAll_BDT_cat[4];
  TH1F* InvMassGood_BDT_cat[4];


  TH1F *hsumpt2_sig;
  TH1F *hsumpt2_bkg;

  TH1F *hptasym_sig;
  TH1F *hptasym_bkg;

  TH1F *hptbal_sig;
  TH1F *hptbal_bkg;

  
 private :
  
  h2gglobeEventReader* ev_;
 
  

};

#endif
