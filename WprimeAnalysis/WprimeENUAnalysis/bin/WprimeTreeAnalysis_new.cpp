
//---- CMSSW includes ----
#include "PhysicsTools/NtupleUtils/interface/hFactory.h"
#include "PhysicsTools/NtupleUtils/interface/h2Factory.h"
#include "PhysicsTools/NtupleUtils/interface/hChain.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//---- tree content ----
#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

//---- std include ----
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

//---- root includes ----
#include "TH1.h"
#include "TH2.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define PI 3.14159265
#define TWOPI 6.28318530


//<<<<------------------------------------------------------------------------------------------>>>>
//<<<<------------------------------------------------------------------------------------------>>>>


//==============================
//==== get input files list ====
//==============================
double deltaPhi(double phi1,double phi2) 
{
  double deltaphi = fabs(phi1-phi2);
  if (deltaphi > TWOPI) deltaphi -= TWOPI;
  if (deltaphi > PI ) deltaphi = TWOPI - deltaphi;
  return deltaphi;
}


//<<<<------------------------------------------------------------------------------------------>>>>
//<<<<------------------------------------------------------------------------------------------>>>>



//======================
//==== main program ====
//======================
int main (int argc, char ** argv)
{



  //==============================
  //==== get input files list ====
  //==============================
  std::string fileName(argv[1]);
  boost::shared_ptr<edm::ProcessDesc> processDesc = edm::readConfigFile(fileName);
  boost::shared_ptr<edm::ParameterSet> parameterSet = processDesc->getProcessPSet();
  edm::ParameterSet subPSetSelections =  parameterSet -> getParameter<edm::ParameterSet>("selections");
  int nEvents_ = subPSetSelections.getUntrackedParameter<int>("nEvents", -1);
  bool MCpresent_ =  subPSetSelections.getUntrackedParameter<bool>("MCpresent", false);
  double maxJetPt_ = subPSetSelections.getUntrackedParameter<double>("maxJetPt", 100.);
  double minEleEt_= subPSetSelections.getUntrackedParameter<double>("minEleEt", 30.);
  int eleId_ = subPSetSelections.getUntrackedParameter<int>("eleId", 0);
  double minEtOverMet_ = subPSetSelections.getUntrackedParameter<double>("minEtOverMet", 0.4);
  double maxEtOverMet_ = subPSetSelections.getUntrackedParameter<double>("maxEtOverMet", 1.5);
  double eleMetPhiMin_ = subPSetSelections.getUntrackedParameter<double>("eleMetPhiMin", 2.5);

  std::string outFile_ = subPSetSelections.getParameter<std::string> ("outFile");

  edm::ParameterSet subPSetInput = parameterSet -> getParameter<edm::ParameterSet>("inputNtuples");
  std::vector<std::string> inputFiles = subPSetInput.getParameter<std::vector<std::string> > ("inputFiles");

  

  //==========================
  //==== global variables ====
  //==========================
  
  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;
 


  //=======================
  // === load the tree ====
  //=======================
  TChain *chain = new TChain ("myanalysis/WprimeAnalysisTree") ;
  WprimeTreeContent treeVars ;
  setBranchAddresses (chain, treeVars) ;


  // Input files
  for(std::vector<std::string>::const_iterator listIt = inputFiles.begin();
      listIt != inputFiles.end(); 
      ++listIt)
    {
      chain -> Add (listIt -> c_str());
      std::cout << *listIt << std::endl;
    }

  int nevents = 0;
  if (nEvents_ == -1) nevents = chain->GetEntries ();
  else nevents = nEvents_;

  std::cout << "FOUND " << nevents << " ENTRIES\n" ;




  //========================================
  //==== define histo entries and steps ====
  //========================================
  int nSteps = 9;
  TH1F* events = new TH1F("events", "events", nSteps, 0., 1.*nSteps);
  std::map<int, int> stepEvents;
  std::map<int, std::string> stepName;

  

  
  //=========================== 
  //==== define hfactory ====== 
  //=========================== 
  hFactory*  histograms  = new hFactory(outFile_);
  h2Factory* histograms2 = new h2Factory(outFile_);
  //EB
  histograms -> add_h1("hEt_EB",     "",  3000,   0.,  1500.,   nSteps+1);
  histograms -> add_h1("hMt_EB",     "",  3000,   0.,  3000.,   nSteps+1);
  histograms -> add_h1("hElePhi_EB",     "",  2000,   -3.5,  3.5,   nSteps+1);
  histograms -> add_h1("hDPhiEleMet_EB",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hEtOverMet_EB",  "",  1000, 0.,  10., nSteps+1);
  //EE
  histograms -> add_h1("hEt_EE",     "",  3000,   0.,  1500.,   nSteps+1);
  histograms -> add_h1("hMt_EE",     "",  3000,   0.,  3000.,   nSteps+1);
  histograms -> add_h1("hElePhi_EE",     "",  2000,   -3.5,  3.5,   nSteps+1);
  histograms -> add_h1("hDPhiEleMet_EE",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hEtOverMet_EE",  "",  1000, 0.,  10., nSteps+1);
  //EB + EE
  histograms -> add_h1("hMet",       "",  3000, 0.,  1500.,  nSteps+1);
  histograms -> add_h1("hEleEta",     "",  2000,   -3.,  3.,   nSteps+1);
  histograms -> add_h1("hNEle",        "",  10,    0.,  10., nSteps+1);
  histograms -> add_h1("hNJets",     "",  50,   0.,  50.,    nSteps+1);

  //fake studies
  histograms -> add_h1("hSCEt",     "",  3000,   0.,  1500.,    nSteps+1);
  
  //======================================
  //==== define additional histograms ====
  //======================================

  
  
  //============================
  //==== preselection steps ====
  //============================
  std::map<int,int> eventsCounter;
  for(std::vector<std::string>::const_iterator listIt = inputFiles.begin();
      listIt != inputFiles.end(); 
      ++listIt)
    {
      TFile* fin = new TFile(listIt -> c_str(), "READ");
      eventsCounter[1]    +=   ((TH1F*)(fin -> Get("eventsCounterTotal/eventsCounterTotal")) )                             -> GetBinContent(1);
      eventsCounter[2]    +=   ((TH1F*)(fin -> Get("eventsCounterGoodEvt/eventsCounterGoodEvt")) )                         -> GetBinContent(1);
      eventsCounter[3]    +=   ((TH1F*)(fin -> Get("eventsCounterHighEtEle/eventsCounterHighEtEle")) )                     -> GetBinContent(1);
      eventsCounter[4]    +=   ((TH1F*)(fin -> Get("eventsCounterPatElectronSequence/eventsCounterPatElectronSequence")) ) -> GetBinContent(1);
      
      delete fin;
    }
  
  int step = 1;
  stepEvents[step] = eventsCounter[step];
  stepName[step] = "1) total events";
  
  step = 2;
  stepEvents[step] = eventsCounter[step];
  stepName[step] = "2) good events";
  
  step = 3;
  stepEvents[step] = eventsCounter[step];
  stepName[step] = "3) high Electrons";
  
  
  
  
  //===========================
  //==== loop over entries ====
  //===========================
  for (int entry = 0; entry < nevents; ++entry) 
    {
      //===================
      //==== get entry ====
      //===================
      chain->GetEntry (entry) ;
      if(entry%100000 == 0) std::cout << "event n. " << entry << std::endl;
      
      //=================================
      //==== step 4: All the events  ====
      //=================================
      step = 4;
      stepEvents[step] ++;
      stepName[step] = "4) after preselections";

      double met = treeVars.pfMet;
      double metPhi = treeVars.pfMetPhi;
      double mex = treeVars.pfMex;
      double mey = treeVars.pfMey;

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hNEle",             step,    treeVars.nElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);
      


      //================================
      //==== step 5: HLT selection  ====
      //================================
      //if (treeVars.HLT_Photon10_L1R == 0) continue;  //prescaled!
      if (treeVars.HLT_Ele15_LW_L1R == 0) continue;
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "5) HLT Photon10 L1R";
      

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hNEle",             step,    treeVars.nElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);




      //================================
      //==== loop over SuperCluster ====
      //================================
      for (int ii = 0; ii < treeVars.nSuperClusters; ++ii)
	{




	}




      //============================
      //==== loop over ele cand ====
      //============================
      int nGoodElectrons = 0;
      int chosenEle = 0;
      for (int i = 0; i < treeVars.nElectrons; ++i)
	{
	  float et  = treeVars.eleEt[i] ;
	  float eta = treeVars.eleEta[i] ;
	  

	  //--------------------------------------------
	  // keep only electrons in ECAL fiducial volume
	  //--------------------------------------------
	  if ( fabs(eta) > EtaCutEB && fabs(eta) < EtaCutEE) continue;
	  

	  //--------------------------------------------
	  // keep only electrons with ET > 30 GeV
	  //--------------------------------------------
	  if ( et < minEleEt_) continue ;
	  

	  //--------------------------------------------
	  // electron ID
	  //--------------------------------------------
	  //instructions here: 
	  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h?revision=1.5&view=markup&pathrev=HEAD
	  //https://twiki.cern.ch/twiki/bin/view/CMS/HEEPSelector

	  //if ( treeVars.eleId[i] != eleId_ ) continue;
	  if ( (treeVars.eleId[i] != 0x0) && ((treeVars.eleId[i] &~ 0x0010) != 0x0) ) continue;  //FIXME!! distinguere barrel ed endcap


	  /*
	  //--------------------------------------------
	  // electron Isolation
	  //--------------------------------------------
	    float combIsoCut  = 0.;
	    float combinedIso = 0.;
	    
	    // isolation
	    if  ( fabs(eta) < EtaCutEB &&  treeVars.eleTrkIso[i] > 7.5)  continue;
	    if  ( fabs(eta) > EtaCutEE &&  treeVars.eleTrkIso[i] > 15.)  continue;
	    
	    combinedIso = treeVars.eleEcalIso[i] + treeVars.eleHcalIsoD1[i] ;
	    if ( fabs(eta) < EtaCutEB ) combIsoCut = 2 + 0.03*et;
	    if ( fabs(eta) > EtaCutEE && et < 50.) combIsoCut = 2.5;
	    if ( fabs(eta) > EtaCutEE && et > 50.) combIsoCut = 2.5 + 0.03*(et-50);
	    if ( combinedIso > combIsoCut ) continue;
	    
	    if ( treeVars.eleHcalIsoD2[i] > 0.5 && fabs(eta) > EtaCutEE) continue;
	  */
	  
	  nGoodElectrons++;
	  chosenEle = i;
	  
	}// end loop over ele cand
      
      
      //================================
      //======== step 6: eleId =========
      //================================
      if( nGoodElectrons == 0 ) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "6) eleId";

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);



      //============================
      //=== step 7: only one ele ===
      //============================
      if ( nGoodElectrons != 1 ) continue;
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "7) single electron";

      double eleEt  = treeVars.eleEt[chosenEle];
      double elePx  = treeVars.elePx[chosenEle];
      double elePy  = treeVars.elePy[chosenEle];
      double elePhi = treeVars.elePhi[chosenEle];
      double eleEta = treeVars.eleEta[chosenEle];
      double elePt  = sqrt(elePx*elePx + elePy*elePy);      
      double cphi   = (elePx*mex + elePy*mey ) / (met*elePt);
      double mt     = sqrt(2*eleEt*met*(1 - cphi));
      double dPhiEleMet = deltaPhi(metPhi,elePhi);
      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(fabs (eleEta) < 1.479)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);
	}
      else
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);
	}




      //=============================
      //=== step 8: met selection ===
      //=============================
      if (eleEt/met < minEtOverMet_ || eleEt/met > maxEtOverMet_) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "8) met selection";
     

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(fabs (eleEta) < 1.479)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);
	}
      else
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);
	}





      //==============================
      //=== step 9: ele - met btob ===
      //==============================
      if (dPhiEleMet < eleMetPhiMin_ ) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "9) ele - met btob";


      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(fabs (eleEta) < 1.479)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);
	}
      else
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);
	}


      
    } // end loop over entries
  



  //======================
  //==== write histos ====
  //======================

  TFile *fout = new TFile(outFile_.c_str(),"RECREATE");
  for(int step = 1; step <= nSteps; ++step)
    {
      events -> SetBinContent(step, stepEvents[step]);
      events -> GetXaxis() -> SetBinLabel(step, stepName[step].c_str());
    }

  events -> Write();

  fout -> Close();

  delete histograms;
  delete histograms2;

  return 0;
}
