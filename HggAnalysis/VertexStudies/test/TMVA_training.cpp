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

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"
#include "TMVA/Factory.h"

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
    
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  

  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  std::cout << std::endl;
  std::cout << std::endl;
  

  // Chain
  TChain* tree = new TChain(treeName.c_str());
  tree->Add((baseDir+inputFile).c_str());
  

  // Declaration of leaf types
  Int_t           evtNumber;
  Int_t           nVertices;
  Float_t         diphopt;
  Float_t         logsumpt2;
  Int_t           nch;
  Int_t           nchthr;
  Float_t         ptmax;
  Float_t         ptmax3;
  Float_t         ptbal;
  Float_t         ptasym;
  Float_t         sumtwd;
  Float_t         dk;
  Int_t           isSig;
  
  tree->SetBranchAddress("diphopt", &diphopt);
  tree->SetBranchAddress("logsumpt2", &logsumpt2);
  tree->SetBranchAddress("nch", &nch);
  tree->SetBranchAddress("nchthr", &nchthr);
  tree->SetBranchAddress("ptmax", &ptmax);
  tree->SetBranchAddress("ptmax3", &ptmax3);
  tree->SetBranchAddress("ptbal", &ptbal);
  tree->SetBranchAddress("ptasym", &ptasym);
  tree->SetBranchAddress("sumtwd", &sumtwd);
  tree->SetBranchAddress("dk", &dk);
  tree->SetBranchAddress("isSig", &isSig);

  

  // Create a new root output file.
  TString outfileName( (outputRootFilePath+outputRootFileName).c_str());
  outfileName = outfileName+".root";
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
  TMVA::Factory *factory = new TMVA::Factory( outputRootFileName.c_str(), outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  factory->AddVariable( "logsumpt2" );
  factory->AddVariable( "ptbal" );
  factory->AddVariable( "ptasym" );
  factory->AddVariable( "nch" );
//   factory->AddVariable( "nchthr" );
//   factory->AddVariable( "ptmax" );
//   factory->AddVariable( "ptmax3" );
//   factory->AddVariable( "sumtwd" );
//   factory->AddVariable( "dk" );
  
  
  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  
  // ====== register trees ====================================================
  //
  // the following method is the prefered one:
  // you can add an arbitrary number of signal or background trees
  factory->AddSignalTree    ( tree,     signalWeight     );
  factory->AddBackgroundTree( tree,     backgroundWeight );
   

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "isSig == 1"; 
   TCut mycutb = "isSig == 0"; 

   // tell the factory to use all remaining events in the trees after training for testing:
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "SplitMode=Random:NormMode=NumEvents:!V" );


   // Boosted Decision Trees: use BDTG ( Gradient Boost )
   factory->BookMethod( TMVA::Types::kBDT, "BDTG",
			"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
}
