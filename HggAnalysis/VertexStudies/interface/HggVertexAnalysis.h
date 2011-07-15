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
  void     analyze(int nentries, int isData, int useWeights, TH1F* h);
  
  int      DiphotonCategory (float eta1, float eta2, float r9_1, float r9_2);

  void     findMCHiggsPhotons( TClonesArray* gp_p4, Short_t* gp_status, Short_t* gp_pdgid, Short_t* gp_mother, int gp_n, 
			       TClonesArray* pho_calopos, TClonesArray* pho_p4, int pho_n,
			       int& passSelection, int& i1, int& i2  );



  // histograms 
  
  //histos
  TH1F* PtAll_sumpt2 ;
  TH1F* PtGood_sumpt2;
  TH1F* PtAll_rank ;
  TH1F* PtGood_rank;
 
  TH1F* EtaAll_sumpt2;
  TH1F* EtaGood_sumpt2;
  TH1F* EtaAll_rank;
  TH1F* EtaGood_rank;
  
  TH1F* NvtAll_sumpt2;
  TH1F* NvtGood_sumpt2;
  TH1F* NvtAll_rank;
  TH1F* NvtGood_rank;
 
  TH1F* NpuAll_sumpt2;
  TH1F* NpuGood_sumpt2;
  TH1F* NpuAll_rank;
  TH1F* NpuGood_rank;

  TH1F* InvMassAll_sumpt2;
  TH1F* InvMassGood_sumpt2;
  TH1F* InvMassAll_rank;
  TH1F* InvMassGood_rank;

  TH1F* InvMassAll_sumpt2_cat[4];
  TH1F* InvMassGood_sumpt2_cat[4];
  TH1F* InvMassAll_rank_cat[4];
  TH1F* InvMassGood_rank_cat[4];

  
 private :
  
  h2gglobeEventReader* ev_;
 
  

};

#endif
