///==== include ====
#include "HggVertexAnalysis.h"
#include <iostream>

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
  std::string weights   = gConfigParser -> readStringOption("Input::weights");
  
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  
  
  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  int isData   = gConfigParser -> readIntOption("Options::isData");

  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  int doVBFselection = gConfigParser -> readDoubleOption("Options::doVBFselection");

  int useWeights      = gConfigParser -> readIntOption("Options::useWeights");
  int poissonWeights  = gConfigParser -> readIntOption("Options::poissonWeights");
  int nAvePU          = gConfigParser -> readIntOption("Options::nAvePU");
  std::string puweightsFileName    = gConfigParser -> readStringOption("Options::puweightsFileName");  
  
  int useKfactors                  = gConfigParser -> readIntOption("Options::useKfactors");
  std::string kfactorsFileName     = gConfigParser -> readStringOption("Options::kfactorsFileName");  
  std::string kfactorsHistoName    = gConfigParser -> readStringOption("Options::kfactorsHistoName");  

  std::cout << std::endl;
  std::cout << std::endl;
  
  //------ use PU weights for MC
  TH1F* hweights;
  if (useWeights){
    TFile weightsFile(puweightsFileName.c_str(),"READ");  
    if ( poissonWeights ){
      std::cout << "<N_PU> for reweighting : " << nAvePU << std::endl;
      char hname[100];
      sprintf(hname,"hwmc%d",nAvePU);
      hweights = (TH1F*)weightsFile.Get(hname)->Clone("hweights");
      hweights->SetDirectory(0);
      if (!hweights) {
	cout << " >>> File for weights not found !!! " << endl;
	return 0;
      }
    }
    else {
      hweights = (TH1F*)weightsFile.Get("hweights")->Clone("hweights");
      hweights->SetDirectory(0);
      cout << hweights ->GetMaximum() << endl;
    }
    weightsFile.Close();
  }


  //------ use K -factors (*** only for Glu-Glu ***)
  TH1F* hkfactors;
  if (useKfactors){
    TFile kfactorsFile(kfactorsFileName.c_str(),"READ");  
    hkfactors = (TH1F*)kfactorsFile.Get(kfactorsHistoName.c_str())->Clone("hkfactors");
    hkfactors->SetDirectory(0);
    cout << hkfactors ->GetMaximum() << endl;
    kfactorsFile.Close();
  }


  
  // Chain
  TChain* chain = new TChain(treeName.c_str());
  chain->Add((baseDir+inputFile).c_str());
   
  h2gglobeEventReader *reader = new h2gglobeEventReader((TTree*)(chain));
  std::cout<<"Found "<< reader->GetEntries() <<" entries"<<std::endl;

  int nentries = entryMAX ; 
  if ( entryMAX < 0 ) nentries = reader->GetEntries();
  
  TFile *fout = new TFile((outputRootFilePath+outputRootFileName).c_str(),"RECREATE");  
  
  HggVertexAnalysis hggAna(reader);
  hggAna.bookHistos();
  //  hggAna.analyze(nentries, isData, useWeights, hweights);
  hggAna.analyze(nentries, isData, useWeights, hweights, useKfactors, hkfactors, doVBFselection);
  hggAna.saveHistos(fout);
  
  std::cout << "CIAOOOOOOO !!!" << std::endl;
  
}

  





