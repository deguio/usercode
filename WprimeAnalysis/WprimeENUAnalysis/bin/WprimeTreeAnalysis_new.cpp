//---- CMSSW includes ----
#include "PhysicsTools/NtupleUtils/interface/hFactory.h"
#include "PhysicsTools/NtupleUtils/interface/h2Factory.h"
#include "PhysicsTools/NtupleUtils/interface/hChain.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEleSelector.h"

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
  bool nonIso_ =  subPSetSelections.getUntrackedParameter<bool>("nonIso", false);
  double maxJetPt_ = subPSetSelections.getUntrackedParameter<double>("maxJetPt", 100.);
  double minEleEt_= subPSetSelections.getUntrackedParameter<double>("minEleEt", 25.);
  int eleId_ = subPSetSelections.getUntrackedParameter<int>("eleId", 0);
  double minEtOverMet_ = subPSetSelections.getUntrackedParameter<double>("minEtOverMet", 0.4);
  double maxEtOverMet_ = subPSetSelections.getUntrackedParameter<double>("maxEtOverMet", 1.5);
  double eleMetPhiMin_ = subPSetSelections.getUntrackedParameter<double>("eleMetPhiMin", 2.5);

  std::string outFile_ = subPSetSelections.getParameter<std::string> ("outFile");
  std::string evtList_ = subPSetSelections.getParameter<std::string> ("evtList");

  edm::ParameterSet subPSetInput = parameterSet -> getParameter<edm::ParameterSet>("inputNtuples");
  std::vector<std::string> inputFiles = subPSetInput.getParameter<std::vector<std::string> > ("inputFiles");

  

  //==========================
  //==== global variables ====
  //==========================
  
  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;
 
  std::ofstream* outFile100;
  std::ofstream* outFile150;
  std::ofstream* outFile200;

  if (MCpresent_ == false)
    {
      outFile100 = new std::ofstream((evtList_+"_100.txt").c_str(), std::ios::out);
      outFile150 = new std::ofstream((evtList_+"_150.txt").c_str(), std::ios::out);
      outFile200 = new std::ofstream((evtList_+"_200.txt").c_str(), std::ios::out);
    }


  //=======================
  // === load the tree ====
  //=======================
  TChain *chain = new TChain ("myanalysis/WprimeAnalysisTree") ;
  WprimeTreeContent treeVars ;
  setBranchAddresses (chain, treeVars) ;

  //===========================
  // === load the HLT tree ====
  //===========================
  int HLTnpaths;
  int HLTwasrun[20];
  int HLTaccept[20];
  int HLTerror[20];

  TChain *chainHLT = new TChain ("TriggerResults/HLTree") ;
  chainHLT -> SetBranchAddress("HLTnpaths", &HLTnpaths);
  chainHLT -> SetBranchAddress("HLTwasrun", HLTwasrun);
  chainHLT -> SetBranchAddress("HLTaccept", HLTaccept);
  chainHLT -> SetBranchAddress("HLTerror" , HLTerror);
  

  // Input files
  for(std::vector<std::string>::const_iterator listIt = inputFiles.begin();
      listIt != inputFiles.end(); 
      ++listIt)
    {
      chain -> Add (listIt -> c_str());
      chainHLT -> Add (listIt -> c_str());
      std::cout << *listIt << std::endl;
    }

  int nevents = 0;
  if (nEvents_ == -1) nevents = chain->GetEntries ();
  else nevents = nEvents_;

  std::cout << "FOUND " << nevents << " ENTRIES\n" ;




  //========================================
  //==== define histo entries and steps ====
  //========================================
  int nSteps = 14;
  TH1F* events = new TH1F("events", "events", nSteps, 0., 1.*nSteps);
  std::map<int, int> stepEvents;
  std::map<int, std::string> stepName;

  

  
  //=========================== 
  //==== define hfactory ====== 
  //=========================== 
  h2Factory* histograms2 = new h2Factory(outFile_);
  hFactory*  histograms  = new hFactory(outFile_);
  //EB
  histograms -> add_h1("hEt_EB",     "",  3000,   0.,  1500.,   nSteps+1);
  histograms -> add_h1("hMt_EB",     "",  3000,   0.,  3000.,   nSteps+1);
  histograms -> add_h1("hElePhi_EB",     "",  2000,   -3.5,  3.5,   nSteps+1);
  histograms -> add_h1("hDPhiEleMet_EB",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hEtOverMet_EB",  "",  1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEtOverMet_vsE_EB",  "", 3000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEoverP_vsE_EB",  "", 1000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("heleSeedTime_vsEleSeedEne_EB", "", 3000, 0., 1500., 400, -20., 20., nSteps+1);
  histograms -> add_h1("hRecoFlag_EB",     "",  20,   0.,  20.,   nSteps+1);
  //EE
  histograms -> add_h1("hEt_EE",     "",  3000,   0.,  1500.,   nSteps+1);
  histograms -> add_h1("hMt_EE",     "",  3000,   0.,  3000.,   nSteps+1);
  histograms -> add_h1("hElePhi_EE",     "",  2000,   -3.5,  3.5,   nSteps+1);
  histograms -> add_h1("hDPhiEleMet_EE",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hEtOverMet_EE",  "",  1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEtOverMet_vsE_EE",  "",  3000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEoverP_vsE_EE",  "",  1000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("heleSeedTime_vsEleSeedEne_EE", "", 3000, 0., 1500., 400, -20., 20., nSteps+1);
  histograms -> add_h1("hRecoFlag_EE",     "",  20,   0.,  20.,   nSteps+1);
  //EB + EE
  histograms -> add_h1("hEt",     "",  3000,   0.,  1500.,   nSteps+1);
  histograms -> add_h1("hMt",     "",  3000,   0.,  3000.,   nSteps+1);
  histograms -> add_h1("hElePhi",     "",  2000,   -3.5,  3.5,   nSteps+1);
  histograms -> add_h1("hDPhiEleMet",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hEtOverMet",  "",  1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEtOverMet_vsE",  "",  3000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("hEoverP_vsE",  "",  1000, 0., 1500., 1000, 0.,  50., nSteps+1);
  histograms2 -> add_h2("heleSeedTime_vsEleSeedEne", "", 3000, 0., 1500., 400, -20., 20., nSteps+1);
  histograms -> add_h1("hSquareSelection_eleEt", "", 3000, 0., 1500., nSteps+1);
  histograms2 -> add_h2("hSquareSelection_eleSeedTime_vsEleSeedEne", "", 3000, 0., 1500., 400, -20., 20., nSteps+1);
  histograms -> add_h1("hCorrected_DPhiEleMet",  "",  1000, 0.,  4., nSteps+1);
  histograms -> add_h1("hCorrected_EtOverMet",  "",  1000, 0.,  50., nSteps+1);
  histograms -> add_h1("hRecoFlag",     "",  20,   0.,  20.,   nSteps+1);

  histograms -> add_h1("hMet",       "",  3000, 0.,  1500.,  nSteps+1);
  histograms -> add_h1("hEleEta",     "",  2000,   -3.,  3.,   nSteps+1);
  histograms -> add_h1("hNEle",        "",  10,    0.,  10., nSteps+1);
  histograms -> add_h1("hNJets",     "",  50,   0.,  50.,    nSteps+1);

  histograms -> add_h1("hMex",        "", 6000, -1500., 1500., nSteps+1);
  histograms -> add_h1("hMey",        "", 6000, -1500., 1500., nSteps+1);

  //fake studies
  histograms -> add_h1("hSCEt",     "",  3000,   0.,  1500.,    nSteps+1);
  



  //===============================================
  //==== define additional histograms and tree ====
  //===============================================

  TH1F* histoEt = new TH1F("histoEt","histoEt", 1500, 0., 1500.);
  TH1F* histoMt = new TH1F("histoMt","histoMt", 3000, 0., 3000.);
  TH1F* hEt_cumulative = new TH1F("h_0_hEt_cumulative","h_0_hEt_cumulative", 1500, 0., 1500.);
  TH1F* hMt_cumulative = new TH1F("h_0_hMt_cumulative","h_0_hMt_cumulative", 3000, 0., 3000.);

  TH1F* histoEt_EB = new TH1F("histoEt_EB","histoEt_EB", 1500, 0., 1500.);
  TH1F* histoMt_EB = new TH1F("histoMt_EB","histoMt_EB", 3000, 0., 3000.);
  TH1F* hEt_cumulative_EB = new TH1F("h_0_hEt_cumulative_EB","h_0_hEt_cumulative_EB", 1500, 0., 1500.);
  TH1F* hMt_cumulative_EB = new TH1F("h_0_hMt_cumulative_EB","h_0_hMt_cumulative_EB", 3000, 0., 3000.);

  TH1F* histoEt_EE = new TH1F("histoEt_EE","histoEt_EE", 1500, 0., 1500.);
  TH1F* histoMt_EE = new TH1F("histoMt_EE","histoMt_EE", 3000, 0., 3000.);
  TH1F* hEt_cumulative_EE = new TH1F("h_0_hEt_cumulative_EE","h_0_hEt_cumulative_EE", 1500, 0., 1500.);
  TH1F* hMt_cumulative_EE = new TH1F("h_0_hMt_cumulative_EE","h_0_hMt_cumulative_EE", 3000, 0., 3000.);

  //tree variables
  //double eleEt;
  //double mt;
  
  //TTree* tree = new TTree("tree","tree");
  //tree -> Branch("eleEt",       &eleEt,             "eleEt/D");
  //tree -> Branch("mt",          &mt,                "mt/D");
  
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
  stepName[step] = "3) high Et electrons";
  
  
  
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

      double met = treeVars.tcMet;
      double metPhi = treeVars.tcMetPhi;
      double mex = treeVars.tcMex;
      double mey = treeVars.tcMey;

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hNEle",             step,    treeVars.nElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);
      
      
      //==============================
      //==== step 5: MC selection ====
      //==============================
      //      if (MCpresent_ == true && fabs(treeVars.pdgId[0]) != 11) continue;  //FIXME dimensione vettore!!!
      //      if (treeVars.nGenParticles != 1) std::cout << "PROBLEM: more than one lepton @gen level" << std::endl;

      step ++;
      stepEvents[step] ++;
      stepName[step] = "5) GEN level selection (only for MC)";
      



      //================================
      //==== step 6: HLT selection  ====
      //================================

      //FIXME hardcoded
      //'0 -> HLT_Photon10_L1R',
      //'1 -> HLT_Ele10_LW_L1R',
      //'2 -> HLT_Ele15_LW_L1R',
      //'3 -> HLT_Ele15_SW_L1R',
      //'4 -> HLT_Ele15_SW_CaloEleId_L1R'
      //'5 -> HLT_Ele17_SW_CaloEleId_L1R'
      //'6 -> HLT_Ele22_SW_CaloEleId_L1R'
      //'7 -> HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1'

      chainHLT->GetEntry(entry);

      int accept = 0;
      if(MCpresent_ == true)  //darren's reccomendation
	{
	  for (int ii = 0; ii < HLTnpaths; ++ ii)    //fixme!! naive approach
	    accept += HLTaccept[ii];
	}

      else
	{
	  if ( treeVars.runId >= 135059 &&  treeVars.runId <= 140041  ) accept = HLTaccept[1];
	  if ( treeVars.runId >= 140042 &&  treeVars.runId <= 141900  ) accept = HLTaccept[3];
	  if ( treeVars.runId >= 141901 &&  treeVars.runId <= 146427  ) accept = HLTaccept[4];
	  if ( treeVars.runId >= 146428 &&  treeVars.runId <= 147116  ) accept = HLTaccept[5];	  
	  if ( treeVars.runId >= 147117 &&  treeVars.runId <= 999999  ) accept = HLTaccept[7];	  
	}

      if (accept == 0) continue;



      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "6) HLT";
      

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hNEle",             step,    treeVars.nElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);



      //=============================
      //==== step 7: HBHE noise  ====
      //=============================
      if (MCpresent_ == false && treeVars.hcalnoiseLoose == 0) continue; //fixme (fintanto che non e' simulato)

      step ++;
      stepEvents[step] ++;
      stepName[step] = "7) HBHE noise";

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hNEle",             step,    treeVars.nElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);




      //============================
      //==== loop over ele cand ====
      //============================
      int nGoodElectrons = 0;
      int chosenEle = 0;
      for (int i = 0; i < treeVars.nElectrons; ++i)
	{

	  double eleEt  = treeVars.eleEt[i];
	  double elePx  = treeVars.elePx[i];
	  double elePy  = treeVars.elePy[i];
	  double elePhi = treeVars.elePhi[i];
	  double eleEta = treeVars.eleEta[i];
	  double elePt  = sqrt(elePx*elePx + elePy*elePy);
	  double eleId  = treeVars.eleId[i];
	  double cphi   = (elePx*mex + elePy*mey ) / (met*elePt);
	  double mt     = sqrt(2*eleEt*met*(1 - cphi));
	  double dPhiEleMet = deltaPhi(metPhi,elePhi);
	  bool eleIsEB = treeVars.eleIsEB[i];
	  bool eleIsEE = treeVars.eleIsEE[i];


	  //--------------------------------------------
	  // keep only electrons in ECAL fiducial volume
	  //--------------------------------------------
	  if ( treeVars.eleIsGap[i] == 1 ) continue;

	  //--------------------------------------------
	  // keep only electrons with ET > XX GeV
	  //--------------------------------------------
	  if ( eleEt < minEleEt_) continue ;
	  
	  //--------------------------------------------
	  // swiss cross e timing (FIXME hardcoded)
	  //--------------------------------------------
	  //if ( treeVars.eleSeedSwissCross[i] > 0.95 ) continue; //included in 38x
	  //if (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) continue; //solo se posso fidarmi della calibrazione del timing

	  //--------------------------------------------
	  // electron ID (FIXME hardcoded)
	  //--------------------------------------------
	  //instructions here: 
	  //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h?revision=1.5&view=markup&pathrev=HEAD
	  //https://twiki.cern.ch/twiki/bin/view/CMS/HEEPSelector

	  // barrel + endcap (same selections in 38X)
	  if (nonIso_ == false)
	    {
	      if ( !(treeVars.eleId[i] == 0x0) ) continue;
	    }
	  else
	    {
	      if ( eleIsEB == 1 && 
		   !(heep::CutCodes::passCuts(eleId,"ecalDriven:et:detEta:dEtaIn:dPhiIn:e2x5Over5x5:hadem:isolPtTrks") && 
		     !heep::CutCodes::passCuts(eleId,"isolEmHadDepth1") ) ) continue;
	      if ( eleIsEE == 1 && 
		   !(heep::CutCodes::passCuts(eleId,"ecalDriven:et:detEta:dEtaIn:dPhiIn:sigmaIEtaIEta:hadem:isolPtTrks") && 
		     !(heep::CutCodes::passCuts(eleId,"isolEmHadDepth1") && heep::CutCodes::passCuts(eleId,"isolHadDepth2")) ) ) continue;
	    }
	  
	  nGoodElectrons++;
	  chosenEle = i;
	  
	}// end loop over ele cand
      
      
      //================================
      //======== step 8: eleId =========
      //================================
      if( nGoodElectrons == 0 ) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "8) eleId";

      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);



      //============================
      //=== step 9: only one ele ===
      //============================
      if ( nGoodElectrons != 1 ) continue;
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "9) single electron";

      double eleE  = treeVars.eleE[chosenEle];
      double eleEmax  = treeVars.eleSeedEnergy[chosenEle];
      double eleTime  = treeVars.eleSeedTime[chosenEle];
      double elePx  = treeVars.elePx[chosenEle];
      double elePy  = treeVars.elePy[chosenEle];
      double elePz  = treeVars.elePz[chosenEle];
      double eleP = sqrt(elePx*elePx + elePy+elePy + elePz*elePz);
      double elePhi = treeVars.elePhi[chosenEle];
      double eleEta = treeVars.eleEta[chosenEle];
      double elePt  = sqrt(elePx*elePx + elePy*elePy);      
      double cphi   = (elePx*mex + elePy*mey ) / (met*elePt);
      double dPhiEleMet = deltaPhi(metPhi,elePhi);
      bool eleIsEB = treeVars.eleIsEB[chosenEle];
      bool eleIsEE = treeVars.eleIsEE[chosenEle];
      //already defined for TTree
      double eleEt  = treeVars.eleEt[chosenEle];
      double mt     = sqrt(2*eleEt*met*(1 - cphi));

      if ( (eleIsEB == true) &&
	   (MCpresent_ == false) &&
	   (treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
	   (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
	   (treeVars.eleSeedEnergy[chosenEle] > 150.) &&
	   (treeVars.eleSeedTime[chosenEle] > -4.)
	   )
	{
	  histograms -> Fill("hSquareSelection_eleEt",  step,    eleEt);
	  histograms2 -> Fill("hSquareSelection_eleSeedTime_vsEleSeedEne", step,  eleEmax, eleTime);
	}


      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms2 -> Fill("hEtOverMet_vsE",   step,    eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",      step,    eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne", step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(eleIsEB == true)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EB",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EB",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EB",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EB",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);
	}
      else if(eleIsEE == true)
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EE",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EE",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EE",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EE",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);
	}



      //==============================================
      //==== computing corrected MET if necessary ====
      //==============================================
      // if  ( (eleIsEB == true) &&
      // 	    (MCpresent_ == false) &&
      // 	    (treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
      // 	    (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
      // 	    (treeVars.eleSeedEnergy[chosenEle] > 150.) &&
      // 	    (treeVars.eleSeedTime[chosenEle] > -4.)
      // 	    )
      // 	{
      // 	  mex = mex - (eleEmax*elePx/eleE);
      // 	  mey = mey - (eleEmax*elePy/eleE);
      // 	  met = sqrt(mex*mex + mey*mey);
      // 	  cphi = (elePx*mex + elePy*mey ) / (met*elePt);
      // 	  mt   = sqrt(2*eleEt*met*(1 - cphi));
      // 	  dPhiEleMet = acos (cphi);
      // 	  histograms -> Fill("hCorrected_EtOverMet",        step,    eleEt/met);
      // 	  histograms -> Fill("hCorrected_DPhiEleMet",       step,    dPhiEleMet);
      // 	}
      

      //===============================
      //=== step 10: ele - met btob ===
      //===============================
      
      if (dPhiEleMet < eleMetPhiMin_) continue; //original
      // bool conditionDeltaphiOK = ( (dPhiEleMet > eleMetPhiMin_) || 
      // 				   ( (eleIsEB == true) &&
      // 				     (MCpresent_ == false) &&
      // 				     (treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
      // 				     (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
      // 				     (treeVars.eleSeedEnergy[chosenEle] > 150.) &&
      // 				     (treeVars.eleSeedTime[chosenEle] > -4.)
      // 				     )
      // 				   );
      // if (!conditionDeltaphiOK) continue;
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "10) ele - met btob";

      if ( (eleIsEB == true) &&
	   (MCpresent_ == false) &&
	   (treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
	   (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
	   (treeVars.eleSeedEnergy[chosenEle] > 150.) &&
	   (treeVars.eleSeedTime[chosenEle] > -4.)
	   )
	{
	  histograms -> Fill("hSquareSelection_eleEt",  step,    eleEt);
	  histograms2 -> Fill("hSquareSelection_eleSeedTime_vsEleSeedEne", step,  eleEmax, eleTime);
	}
      
      
      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms2 -> Fill("hEtOverMet_vsE",    step,   eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",        step,   eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne",  step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(eleIsEB == true)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EB",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EB",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EB",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EB",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);

	}
      else if(eleIsEE == true)
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EE",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EE",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EE",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EE",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);

	}

      
      //==============================
      //=== step 11: met selection ===
      //==============================

      if (eleEt/met < minEtOverMet_ || eleEt/met > maxEtOverMet_) continue;  //original
      // bool conditionMetOK = ( (eleEt/met > minEtOverMet_ && eleEt/met < maxEtOverMet_) || 
      // 			      ( (eleIsEB == true) &&
      // 				(MCpresent_ == false) &&
      // 				(treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
      // 				(treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
      // 				(treeVars.eleSeedEnergy[chosenEle] > 150.) &&
      // 				(treeVars.eleSeedTime[chosenEle] > -4.)
      // 				)
      // 			      );
      // if (!conditionMetOK) continue;  //fixme hardcoded
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "11) met selection";

      
      if ( (eleIsEB == true) &&
	   (MCpresent_ == false) &&
	   (treeVars.eleSeedSwissCross[chosenEle] < 0.95) &&
	   (treeVars.ecalRecHitRecoFlag[chosenEle] == 2) && //koutoftime
	   (treeVars.eleSeedEnergy[chosenEle] > 150.) &&
	   (treeVars.eleSeedTime[chosenEle] > -4.)
	   )
	{
	  histograms -> Fill("hSquareSelection_eleEt",  step,    eleEt);
	  histograms2 -> Fill("hSquareSelection_eleSeedTime_vsEleSeedEne", step,  eleEmax, eleTime);
	}
      
     
      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms2 -> Fill("hEtOverMet_vsE",    step,   eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",        step,   eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne",  step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if(eleIsEB == true)
	{
	  histograms -> Fill("hEt_EB",               step,    eleEt);
	  histograms -> Fill("hMt_EB",               step,    mt);
	  histograms -> Fill("hEtOverMet_EB",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EB",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EB",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EB",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EB",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EB",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EB",           step,    elePhi);
	}
      else if(eleIsEE == true)
	{
	  histograms -> Fill("hEt_EE",               step,    eleEt);
	  histograms -> Fill("hMt_EE",               step,    mt);
	  histograms -> Fill("hEtOverMet_EE",        step,    eleEt/met);
	  histograms2 -> Fill("hEtOverMet_vsE_EE",    step,   eleE, eleEt/met);
	  histograms2 -> Fill("hEoverP_vsE_EE",        step,   eleE, eleE/eleP);
	  histograms2 -> Fill("heleSeedTime_vsEleSeedEne_EE",  step,  eleEmax, eleTime);
	  histograms -> Fill("hRecoFlag_EE",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
	  histograms -> Fill("hDPhiEleMet_EE",       step,    dPhiEleMet);
	  histograms -> Fill("hElePhi_EE",           step,    elePhi);
	}


      //==================================
      // ==== END ANALYSIS SELECTIONS ====
      //==================================



      histoEt -> Fill(eleEt);
      histoMt -> Fill(mt);
      
      if(eleIsEB == true)
	{
	  histoEt_EB -> Fill(eleEt);
	  histoMt_EB -> Fill(mt);
	}
      else if(eleIsEE == true)
	{
	  histoEt_EE -> Fill(eleEt);
	  histoMt_EE -> Fill(mt);
	}

      //tree -> Fill();

      //========================================================
      //=== step 12: saving interesting events above 100 GeV ===
      //========================================================
      if(mt < 100.) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "12) mT > 100 GeV";
     
      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms2 -> Fill("hEtOverMet_vsE",    step,   eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",        step,   eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne",  step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);
            
      if (MCpresent_ == false)
	(*outFile100) << "Run = " << treeVars.runId << "; LS = " << treeVars.lumiId << "; evt = " << treeVars.eventId
		      << "; mT = " << mt << "; Et = " << eleEt << "; tcMET = " << met << "\n";
      
      //========================================================
      //=== step 13: saving interesting events above 150 GeV ===
      //========================================================
      if(mt < 150.) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "13) mT > 150 GeV";

      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms2 -> Fill("hEtOverMet_vsE",    step,   eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",        step,   eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne",  step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);
      
      if (MCpresent_ == false)
	(*outFile150) << "Run = " << treeVars.runId << "; LS = " << treeVars.lumiId << "; evt = " << treeVars.eventId
		      << "; mT = " << mt << "; Et = " << eleEt << "; tcMET = " << met << "\n";

      //========================================================
      //=== step 14: saving interesting events above 200 GeV ===
      //========================================================
      if(mt < 200.) continue;

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "14) mT > 200 GeV";

      histograms -> Fill("hEt",               step,    eleEt);
      histograms -> Fill("hMt",               step,    mt);
      histograms -> Fill("hEtOverMet",        step,    eleEt/met);
      histograms -> Fill("hDPhiEleMet",       step,    dPhiEleMet);
      histograms2 -> Fill("hEtOverMet_vsE",    step,   eleE, eleEt/met);
      histograms2 -> Fill("hEoverP_vsE",        step,   eleE, eleE/eleP);
      histograms2 -> Fill("heleSeedTime_vsEleSeedEne",  step,  eleEmax, eleTime);
      histograms -> Fill("hRecoFlag",         step,    treeVars.ecalRecHitRecoFlag[chosenEle]);
      histograms -> Fill("hElePhi",           step,    elePhi);      
      histograms -> Fill("hMet",              step,    met);
      histograms -> Fill("hMex",              step,    mex);
      histograms -> Fill("hMey",              step,    mey);
      histograms -> Fill("hEleEta",           step,    eleEta);
      histograms -> Fill("hNEle",             step,    nGoodElectrons);
      histograms -> Fill("hNJets",            step,    treeVars.nJets);

      if (MCpresent_ == false)
	(*outFile200) << "Run = " << treeVars.runId << "; LS = " << treeVars.lumiId << "; evt = " << treeVars.eventId
		      << "; mT = " << mt << "; Et = " << eleEt << "; tcMET = " << met << "\n";


    } // end loop over entries
  
  //====================================
  //==== evalutate other quantities ====
  //====================================

  for (int kk = 1; kk < histoEt->GetNbinsX()+1; ++kk)
    hEt_cumulative -> SetBinContent(kk, histoEt->Integral(kk,histoEt->GetNbinsX()));
  for (int kk = 1; kk < histoMt->GetNbinsX()+1; ++kk)
    hMt_cumulative -> SetBinContent(kk, histoMt->Integral(kk,histoMt->GetNbinsX()));

  for (int kk = 1; kk < histoEt_EB->GetNbinsX()+1; ++kk)
    hEt_cumulative_EB -> SetBinContent(kk, histoEt_EB->Integral(kk,histoEt_EB->GetNbinsX()));
  for (int kk = 1; kk < histoMt_EB->GetNbinsX()+1; ++kk)
    hMt_cumulative_EB -> SetBinContent(kk, histoMt_EB->Integral(kk,histoMt_EB->GetNbinsX()));

  for (int kk = 1; kk < histoEt_EE->GetNbinsX()+1; ++kk)
    hEt_cumulative_EE -> SetBinContent(kk, histoEt_EE->Integral(kk,histoEt_EE->GetNbinsX()));
  for (int kk = 1; kk < histoMt_EE->GetNbinsX()+1; ++kk)
    hMt_cumulative_EE -> SetBinContent(kk, histoMt_EE->Integral(kk,histoMt_EE->GetNbinsX()));

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
  //tree -> Write();


  fout -> mkdir("hEt_cumulative");
  fout -> cd("hEt_cumulative");
  hEt_cumulative -> Write();
  fout -> cd();

  fout -> mkdir("hMt_cumulative");
  fout -> cd("hMt_cumulative");
  hMt_cumulative -> Write();
  fout -> cd();

  fout -> mkdir("hEt_cumulative_EB");
  fout -> cd("hEt_cumulative_EB");
  hEt_cumulative_EB -> Write();
  fout -> cd();

  fout -> mkdir("hMt_cumulative_EB");
  fout -> cd("hMt_cumulative_EB");
  hMt_cumulative_EB -> Write();
  fout -> cd();

  fout -> mkdir("hEt_cumulative_EE");
  fout -> cd("hEt_cumulative_EE");
  hEt_cumulative_EE -> Write();
  fout -> cd();

  fout -> mkdir("hMt_cumulative_EE");
  fout -> cd("hMt_cumulative_EE");
  hMt_cumulative_EE -> Write();
  fout -> cd();

  fout -> Close();
  
  if (MCpresent_ == false)
    {
      outFile100 -> close();
      outFile150 -> close();
      outFile200 -> close();
    }

  delete histograms;
  delete histograms2;
  //delete tree;

  return 0;
}
