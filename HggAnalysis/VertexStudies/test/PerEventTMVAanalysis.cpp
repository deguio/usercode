
///==== include ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "readJSONFile.h"

#include "TH1F.h"
#include "TH2F.h"
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
#include "../src/selection.cc"
#include "../src/eleId95.cc"

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif



#define R_ECAL    129
#define Z_ENDCAP  317

#define etaEB   1.4442 
#define etaEE   1.566


using namespace std;
 

// --------- MAIN -------------------

int main(int argc, char** argv)
{ 
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> VertexStudiesAnalysis::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
  
  // Parse the config file
  parseConfigFile (argv[1]) ;

  std::string baseDir   = gConfigParser -> readStringOption("Input::baseDir");
  std::string inputFile = gConfigParser -> readStringOption("Input::inputFile");
  
  std::string treeName  = gConfigParser -> readStringOption("Input::treeName");

  std::string tmvaMethod    = gConfigParser -> readStringOption("Input::tmvaMethod");
  std::string tmvaWeights   = gConfigParser -> readStringOption("Input::tmvaWeights");

  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  
  
  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  std::cout << std::endl;
  std::cout << std::endl;
  
  // Chain
  TChain* tree = new TChain(treeName.c_str());
  tree->Add((baseDir+inputFile).c_str());
  
  // read tree
  int evtNumber, nVertices, npu;
  float trueVtx_z;
  int nConv;
  int isSig[5];
  float ptbal[5], ptasym[5], logsumpt2[5], diphopt[5], pv_z[5], mva[5];
  float  photonSC_eta[2], photonSC_phi[2], photonSC_E[2];   

  tree -> SetBranchAddress("evtNumber", &evtNumber);
  tree -> SetBranchAddress("nVertices", &nVertices);
  tree -> SetBranchAddress("npu", &npu);
  tree -> SetBranchAddress("photonSC_eta",         photonSC_eta);
  tree -> SetBranchAddress("photonSC_phi",         photonSC_phi);
  tree -> SetBranchAddress("photonSC_E",           photonSC_E);
  tree -> SetBranchAddress("logsumpt2",   logsumpt2);
  tree -> SetBranchAddress("ptbal",    ptbal);
  tree -> SetBranchAddress("ptasym",   ptasym);
  tree -> SetBranchAddress("diphopt",  diphopt);
  tree -> SetBranchAddress("mva",  mva);
  tree -> SetBranchAddress("pv_z", pv_z);
  tree -> SetBranchAddress("isSig", isSig);
  tree -> SetBranchAddress("trueVtx_z", &trueVtx_z);
  tree -> SetBranchAddress("nConv", &nConv);
 
  float dZ[3];
  

  // TMVA per event reader
  float diphoPt0, nVert, MVA0, MVA1, MVA2, dZ1, dZ2, nconv;

  TMVA::Reader *tmvaPerEvtReader_ = new TMVA::Reader( "!Color:!Silent" );
  tmvaPerEvtReader_->AddVariable( "diphoPt0", &diphoPt0 );
  tmvaPerEvtReader_->AddVariable( "nVert"   , &nVert   );
  tmvaPerEvtReader_->AddVariable( "MVA0"    , &MVA0	   );
  tmvaPerEvtReader_->AddVariable( "MVA1"    , &MVA1	   );
  tmvaPerEvtReader_->AddVariable( "dZ1"     , &dZ1	   );
  tmvaPerEvtReader_->AddVariable( "MVA2"    , &MVA2	   ); 
  tmvaPerEvtReader_->AddVariable( "dZ2"     , &dZ2	   );
  tmvaPerEvtReader_->AddVariable( "nConv"   , &nconv	   );
  tmvaPerEvtReader_->BookMVA(  tmvaMethod.c_str() , tmvaWeights.c_str() );


  //histograms 

  TH1F *tmvaPerEventOutput       = new TH1F ("tmvaPerEventOutput","tmvaPerEventOutput",400,-2.,2.);
  TH1F *tmvaPerEventOutput_RightVtx = new TH1F("tmvaPerEventOutput_RightVtx" ,"tmvaPerEventOutput_RightVtx" ,200,-1,1);
  TH1F *tmvaPerEventOutput_WrongVtx = new TH1F("tmvaPerEventOutput_WrongVtx" ,"tmvaPerEventOutput_WrongVtx" ,200,-1,1);

  TH2F *tmvaPerEvent_vs_dzSelVtx = new TH2F ("tmvaPerEvent_vs_dzSelVtx","tmvaPerEvent_vs_dzSelVtx",100,-1.,1., 100, -10,10);
  TH2F *tmvaPerEvent_vs_diphopt  = new TH2F ("tmvaPerEvent_vs_diphopt","tmvaPerEvent_vs_diphopt",100,-1.,1., 200, 0,200);
 
  //start loop over entries
  for (int u = 0; u < tree->GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      tree->GetEntry(u);
      
      diphoPt0 = diphopt[0];
      nVert = (float)nVertices;
      MVA0  = mva[0];
      MVA1  = mva[1];
      MVA2  = mva[2];
      dZ1   = pv_z[1]-pv_z[0];
      dZ2   = pv_z[2]-pv_z[0];
      nconv = (float)nConv;

      float evtmva = tmvaPerEvtReader_->EvaluateMVA(tmvaMethod);

      // fill histos
      tmvaPerEventOutput->Fill(evtmva);
      tmvaPerEvent_vs_dzSelVtx -> Fill(evtmva, pv_z[0]- trueVtx_z);
      tmvaPerEvent_vs_diphopt -> Fill(evtmva, diphoPt0);

      if (fabs(pv_z[0]- trueVtx_z) < 1.) 
	tmvaPerEventOutput_RightVtx -> Fill(evtmva);
      else
	tmvaPerEventOutput_WrongVtx -> Fill(evtmva);


    } // end loop over entries
  
   
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
  
  tmvaPerEventOutput->Write();
  tmvaPerEvent_vs_dzSelVtx -> Write();
  tmvaPerEvent_vs_diphopt -> Write();
  tmvaPerEventOutput_RightVtx -> Write();
  tmvaPerEventOutput_WrongVtx -> Write();

	
  ff.Close();
  return 0; 
  
}











