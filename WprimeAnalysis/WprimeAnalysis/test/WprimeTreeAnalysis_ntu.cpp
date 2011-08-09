//---- CMSSW includes ----
#include "treeReader.h"
#include "WprimeTreeVariables.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "hChain.h"

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

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

#define EtaCutEB 1.442
#define EtaCutEEmin 1.560
#define EtaCutEEmax 2.5

#define PI 3.141592654


//======================
//==== main program ====
//======================
int main (int argc, char ** argv)
{

  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> WprimeAnalysis usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }

  //==============================
  //==== get input files list ====
  //==============================

  // Parse the config file
  parseConfigFile (argv[1]) ;

  std::string inputFileList_ = gConfigParser -> readStringOption("Input::InputFileList");
  std::string treeName_      = gConfigParser -> readStringOption("Input::TreeName");
  std::string JSON_          = gConfigParser -> readStringOption("Input::jsonFileName");
  std::string WMass_         = gConfigParser -> readStringOption("Input::WMass");
  
  int         nEvents_       = gConfigParser ->  readIntOption("Options::nEvents");
  int         dataFlag_      = gConfigParser ->  readIntOption("Options::dataFlag");
  int         nonIso_        = gConfigParser ->  readIntOption("Options::nonIso");
  int         doVBTFeleId_        = gConfigParser ->  readIntOption("Options::doVBTFeleId");
  float crossSection_        = gConfigParser -> readFloatOption("Options::crossSection");

  std::string outFile_       = gConfigParser -> readStringOption("Output::OutputFile");  

  float minEleEt_fir_            = gConfigParser -> readFloatOption("Cuts::minEleEtFir");  
  float minEleEt_sec_            = gConfigParser -> readFloatOption("Cuts::minEleEtSec");  
  float minMET_            = gConfigParser -> readFloatOption("Cuts::minMET");  
  float minEtOverMet_        = gConfigParser -> readFloatOption("Cuts::minEtOverMet");  
  float maxEtOverMet_        = gConfigParser -> readFloatOption("Cuts::maxEtOverMet");  
  float minEleMetDPhi_       = gConfigParser -> readFloatOption("Cuts::minEleMetDPhi");  

  //---------------------------
  std::map<int, int> totalEvents           = GetTotalEvents("AllEvents/totalEvents", inputFileList_.c_str());
  std::map<int, int> preselectionsEvents   = GetTotalEvents("ElectronsFilterEvents/totalEvents", inputFileList_.c_str());
  //std::map<int, int> preselectionsEvents   = GetTotalEvents("ElectronsPFlowFilterEvents/totalEvents", inputFileList_.c_str());
  //std::map<int, int> preselectionsEvents   = GetTotalEvents("PhotonsFilterEvents/totalEvents", inputFileList_.c_str());


  std::cout << std::endl;  
  std::cout << std::endl;  

  //==========================
  //==== global variables ====
  //==========================
  //define the map with the events
  std::map<std::pair<int,std::pair<int,int> >,int> eventsMap;  

  //=======================
  // === load the tree ====
  //=======================
  // Open tree
  TChain* chain = new TChain(treeName_.c_str());
  if(!FillChain(*chain, inputFileList_.c_str())) return 1;

  treeReader reader((TTree*)(chain));

  std::cout << std::endl;
  std::cout << std::endl;

  //====================
  //==== parse JSON ====
  //====================

  // Get run/LS map from JSON file
  std::cout << ">>> Get run/LS map from JSON file" << std::endl;
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(JSON_);

  

  int nevents = 0;
  if (nEvents_ == -1) nevents = chain->GetEntries ();
  else nevents = nEvents_;


  std::cout<<"running over "<< nevents <<" entries"<<std::endl;

  //========================================
  //==== define histo entries and steps ====
  //========================================
  int nSteps = 10;
  TH1F* events = new TH1F("events", "events", nSteps, 0., 1.*nSteps);
  std::map<int, int> stepEvents;
  std::map<int, std::string> stepName;


  //===================================
  //==== define variable container ====
  //===================================
  WprimeVariables vars;
  InitializeTree(vars, outFile_);

  // define clone trees
  int firstStep = 8;
  int lastStep = 10;
  std::map<int, TTree*> cloneTrees;
  for(int step = firstStep; step <= lastStep; ++step)
    {
      char treeName[50];
      sprintf(treeName, "ntu_%d", step);
      cloneTrees[step] = CloneTree(vars);
      cloneTrees[step] -> SetName(treeName); 
    }  
  

  //==========================
  //==== define HLT paths ====
  //==========================

  std::vector<std::string> HLTPathNamesMC;

  //trigger Spring11
  HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v1");
  HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v2");
  HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v3");

  //trigger Fall10
  HLTPathNamesMC.push_back("HLT_Ele12_SW_TightEleId_L1R");



  //------------------------------------
  std::vector<std::string> HLTPathNamesDATA;

  if(nonIso_ == 0)
    {
      //   HLTPathNamesDATA.push_back("HLT_Ele45_CaloIdVT_TrkIdT_v1");
      //   HLTPathNamesDATA.push_back("HLT_Ele45_CaloIdVT_TrkIdT_v2");
      //   HLTPathNamesDATA.push_back("HLT_Ele45_CaloIdVT_TrkIdT_v3");
      
      //ele27 e ele32 dovrebbero essere esclusivi
      HLTPathNamesDATA.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1");
      HLTPathNamesDATA.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
      HLTPathNamesDATA.push_back("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3");
      
      HLTPathNamesDATA.push_back("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1");
      HLTPathNamesDATA.push_back("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
      HLTPathNamesDATA.push_back("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3");
      
      //1E33 menu
      HLTPathNamesDATA.push_back("HLT_Ele25_WP80_PFMT40_v1");
      HLTPathNamesDATA.push_back("HLT_Ele27_WP80_PFMT50_v1");

      HLTPathNamesDATA.push_back("HLT_Ele52_CaloIdVT_TrkIdT_v1");
      HLTPathNamesDATA.push_back("HLT_Ele52_CaloIdVT_TrkIdT_v2");
      HLTPathNamesDATA.push_back("HLT_Ele52_CaloIdVT_TrkIdT_v3");
      HLTPathNamesDATA.push_back("HLT_Ele52_CaloIdVT_TrkIdT_v4");

      HLTPathNamesDATA.push_back("HLT_Ele65_CaloIdVT_TrkIdT_v1");
      HLTPathNamesDATA.push_back("HLT_Ele65_CaloIdVT_TrkIdT_v2");
      HLTPathNamesDATA.push_back("HLT_Ele65_CaloIdVT_TrkIdT_v3");
      HLTPathNamesDATA.push_back("HLT_Ele65_CaloIdVT_TrkIdT_v4");
    }
  else
    {
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v1");
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v2");
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v3");  //presc 500
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v4");  //presc 1000
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v5");  
      HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v6");  
    }


  //============================
  //==== preselection steps ====
  //============================
  
  int step = 1;
  stepEvents[step] = totalEvents[1];
  stepName[step] = "1) total events";
  
  step = 2;
  stepEvents[step] = preselectionsEvents[1];
  stepName[step] = "2) preselections";
  
  
  
  //===========================
  //==== loop over entries ====
  //===========================
  for (int entry = 0; entry < nevents; ++entry) 
    {
      std::vector<std::string>* HLT_Names = reader.GetString("HLT_Names");
      std::vector<float>* HLT_WasRun = reader.GetFloat("HLT_WasRun");
      std::vector<float>* HLT_Accept = reader.GetFloat("HLT_Accept");
      std::vector<float>* HLT_Error = reader.GetFloat("HLT_Error");
      std::vector<int>* HLT_Prescale = reader.GetInt("HLT_Prescale");
      std::vector<float>* mc_ptHat = reader.GetFloat("mc_ptHat");
      std::vector<int>* mc_PU_NumInteractions = reader.GetInt("mc_PUit_NumInteractions");

      std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
      std::vector<float>* electrons_charge = reader.GetFloat("electrons_charge");
      std::vector<float>* electrons_z = reader.GetFloat("electrons_z");
      std::vector<float>* electrons_dB = reader.GetFloat("electrons_dB");
      std::vector<float>* electrons_edB = reader.GetFloat("electrons_edB");
      std::vector<float>* electrons_dxy_BS = reader.GetFloat("electrons_dxy_BS");
      std::vector<float>* electrons_dz_BS = reader.GetFloat("electrons_dz_BS");
      std::vector<float>* electrons_dxy_PV = reader.GetFloat("electrons_dxy_PV");
      std::vector<float>* electrons_dz_PV = reader.GetFloat("electrons_dz_PV");
      std::vector<float>* electrons_tkIsoR03 = reader.GetFloat("electrons_tkIsoR03");
      std::vector<float>* electrons_tkIsoR04 = reader.GetFloat("electrons_tkIsoR04");
      std::vector<float>* electrons_emIsoR03 = reader.GetFloat("electrons_emIsoR03");
      std::vector<float>* electrons_emIsoR04 = reader.GetFloat("electrons_emIsoR04");
      std::vector<float>* electrons_hadIsoR03_depth1 = reader.GetFloat("electrons_hadIsoR03_depth1");
      std::vector<float>* electrons_hadIsoR03_depth2 = reader.GetFloat("electrons_hadIsoR03_depth2");
      std::vector<float>* electrons_hadIsoR04_depth1 = reader.GetFloat("electrons_hadIsoR04_depth1");
      std::vector<float>* electrons_hadIsoR04_depth2 = reader.GetFloat("electrons_hadIsoR04_depth2");
      std::vector<int>* electrons_isEB = reader.GetInt("electrons_isEB");
      std::vector<int>* electrons_ecalDrivenSeed = reader.GetInt("electrons_ecalDrivenSeed");
      std::vector<int>* electrons_trackerDrivenSeed = reader.GetInt("electrons_trackerDrivenSeed");
      std::vector<float>* electrons_mva = reader.GetFloat("electrons_mva");
      std::vector<float>* electrons_eSC = reader.GetFloat("electrons_eSC");
      std::vector<float>* electrons_pin = reader.GetFloat("electrons_pin");
      std::vector<float>* electrons_pout = reader.GetFloat("electrons_pout");
      std::vector<float>* electrons_pcalo = reader.GetFloat("electrons_pcalo");
      std::vector<float>* electrons_eSCOverP = reader.GetFloat("electrons_eSCOverP");
      std::vector<float>* electrons_eSeedOverP = reader.GetFloat("electrons_eSeedOverP");
      std::vector<int>* electrons_classification = reader.GetInt("electrons_classification");
      std::vector<float>* electrons_fbrem = reader.GetFloat("electrons_fbrem");
      std::vector<float>* electrons_hOverE = reader.GetFloat("electrons_hOverE");
      std::vector<float>* electrons_deltaPhiIn = reader.GetFloat("electrons_deltaPhiIn");
      std::vector<float>* electrons_deltaEtaIn = reader.GetFloat("electrons_deltaEtaIn");
      std::vector<float>* electrons_sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta");
      std::vector<float>* electrons_e1x5 = reader.GetFloat("electrons_e1x5");
      std::vector<float>* electrons_e2x5Max = reader.GetFloat("electrons_e2x5Max");
      std::vector<float>* electrons_e5x5 = reader.GetFloat("electrons_e5x5");
      std::vector<int>* electrons_mishits = reader.GetInt("electrons_mishits");
      std::vector<int>* electrons_nAmbiguousGsfTracks = reader.GetInt("electrons_nAmbiguousGsfTracks");

      //===================
      //==== get entry ====
      //===================
      reader.GetEntry (entry) ;
      ClearWprimeVariables(vars);
      if(entry%100000 == 0) std::cout << "event n. " << entry << std::endl;

      //===================================
      //==== REMOVE DUPLICATES IN DATA ====
      //===================================
      if( dataFlag_ == 1 )
	{
	  std::pair<int,int> eventLSandID(reader.GetInt("lumiId")->at(0), reader.GetInt("eventId")->at(0));
	  std::pair<int,std::pair<int,int> > eventRUNandLSandID(reader.GetInt("runId")->at(0), eventLSandID);
	  
	  if( eventsMap[eventRUNandLSandID] == 1 ) continue;
	  else eventsMap[eventRUNandLSandID] = 1;
	}
      
      //=============================
      //==== Pt hat cut for Wenu ====
      //=============================

      //if (mc_ptHat->at(0) > 100. ) continue;

      //=================================
      //==== step 3: All the events  ====
      //=================================

      vars.mW = atof(WMass_.c_str());
      vars.totEvents = totalEvents[1];
      vars.crossSection = crossSection_;
      vars.dataFlag = dataFlag_;
      vars.runId   = reader.GetInt("runId")->at(0);
      vars.lumiId  = reader.GetInt("lumiId")->at(0);
      vars.eventId = reader.GetInt("eventId")->at(0);

      if (dataFlag_ == 0) SetPUVariables(vars, reader);
      SetPVVariables(vars, reader);
      SetMetVariables(vars, reader);

      step = 3 ;
      stepEvents[step] ++ ;
      stepName[step] = "3) after preselections";

      double met    = vars.met.pt();
      double metPhi = vars.met.phi();
      double mex    = vars.met.px();
      double mey    = vars.met.py();


      //========================
      //==== JSON selection ====
      //========================
      
      bool skipEvent = false;
      if( dataFlag_ == 1 )
	{
	  int runId  = vars.runId;
	  int lumiId = vars.lumiId;
	  if(AcceptEventByRunAndLumiSection(runId, lumiId, jsonMap) == false) skipEvent = true;      
	}
      
      if( skipEvent == true ) continue;
      
      step ++;
      stepEvents[step] ++;
      stepName[step] = "4) JSON";
      
      //================================
      //==== step 6: HLT selection  ====
      //================================

      int prescale = 0;
      if(dataFlag_ == 0)
	{
	  for(unsigned int HLTIt = 0; HLTIt < HLTPathNamesMC.size(); ++HLTIt)
	    {
	      if ( GetPrescale(reader, HLTPathNamesMC.at(HLTIt)) > 0)
		prescale = GetPrescale(reader, HLTPathNamesMC.at(HLTIt));
	    }
	  // prescale = 1;
   	}
      else
	{
	  for(unsigned int HLTIt = 0; HLTIt < HLTPathNamesDATA.size(); ++HLTIt)
	    {
	      if ( GetPrescale(reader, HLTPathNamesDATA.at(HLTIt)) > 0)
		prescale = GetPrescale(reader, HLTPathNamesDATA.at(HLTIt));
	    }
	}
      
      if( prescale == 0 ) continue;
      vars.hltPrescale = prescale;


      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "5) HLT";
      

      //=============================
      //==== step 7: HBHE noise  ====
      //=============================
      if (dataFlag_ == 1 && reader.GetInt("HCAL_noise")->at(0) == 0 ) continue;
      
      step ++;
      stepEvents[step] ++;
      stepName[step] = "7) HBHE noise";
     

      //============================
      //==== loop over ele cand ====
      //============================
      int nGoodElectrons = 0;
      int chosenEle = 0;
      for (unsigned int i = 0; i < reader.Get4V("electrons")->size(); ++i)
  	{
	  
  	  double eleEt  = electrons->at(i).pt(); //FIXME: VOGLIO ESC*COS(THETATRACCIA)
  	  double elePt  = electrons->at(i).pt(); //FIXME: VOGLIO ESC*COS(THETATRACCIA)
	  
  	  double elePhi = electrons->at(i).phi(); //FIXME: VOGLIO ETA E PHI DELLA TRACCIA
  	  double eleEta = electrons->at(i).eta(); //FIXME: VOGLIO ETA E PHI DELLA TRACCIA
	  	  
  	  bool eleIsEB = electrons_isEB->at(i);
	  
	  float e2x5_overE5x5 = electrons_e2x5Max->at(i)/electrons_e5x5->at(i);
	  float e1x5_overE5x5 = electrons_e1x5->at(i)/electrons_e5x5->at(i);
	  
	  float tkIso  = electrons_tkIsoR03->at(i); 
	  float emIso  = electrons_emIsoR03->at(i);
	  float hadIso = electrons_hadIsoR03_depth1->at(i) + electrons_hadIsoR03_depth2->at(i);
	  
	  float sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta")->at(i);
	  float DphiIn = fabs(reader.GetFloat("electrons_deltaPhiIn")->at(i));
	  float DetaIn = fabs(reader.GetFloat("electrons_deltaEtaIn")->at(i));
	  float HOverE = reader.GetFloat("electrons_hOverE")->at(i);
	  	  
	  
  	  //--------------------------------------------
  	  // keep only electrons in ECAL fiducial volume
  	  //--------------------------------------------
  	  //if ( treeVars.eleIsGap[i] == 1 ) continue;
  	  if ( fabs(eleEta) > EtaCutEB && fabs(eleEta) < EtaCutEEmin ) continue;
	  if ( fabs(eleEta) > EtaCutEEmax) continue;
  	  //--------------------------------------------
  	  // keep only electrons with ET > XX GeV
  	  //--------------------------------------------
  	  if ( eleEt < minEleEt_sec_) continue ; //at least 25 as in 2010
	  
	  //================================================================================================================================================
	  //--------------------
	  //---- HEEP eleId ----
	  //--------------------
	  
	  if (eleIsEB == true && doVBTFeleId_ == 0) //EB HEEP
	    {
	      if ( electrons_ecalDrivenSeed->at(i) != 1 ) continue;
	  
	      if ( fabs(electrons_deltaEtaIn->at(i)) > 0.005 ) continue;
	      if ( fabs(electrons_deltaPhiIn->at(i)) > 0.09 ) continue;
	      if ( electrons_hOverE->at(i) > 0.05 ) continue;
	      
	      if ( e2x5_overE5x5 < 0.94 && e1x5_overE5x5 < 0.83) continue;
	      
	      //iso variables
	      bool iso_EM_HAD = ( (electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i)) < (2.+0.03*eleEt) );
	      bool iso_tk     = ( electrons_tkIsoR03->at(i) < 7.5 );
	      if (nonIso_ == 0)
		{
		  //if ( ( electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i) ) > (2.+0.03*eleEt) ) continue;
		  //if ( electrons_tkIsoR03->at(i) > 7.5 ) continue;

		  if (!iso_EM_HAD) continue;
		  if (!iso_tk) continue;
		}
	      else
		{
		  //not ISO in at least one of the isolation def
		  if (iso_EM_HAD) continue;
                  if (iso_tk) continue;
		}
	    }
	  if (eleIsEB == false && doVBTFeleId_ == 0) //EE HEEP
	    {
	      if ( electrons_ecalDrivenSeed->at(i) != 1 ) continue;

	      if ( fabs(electrons_deltaEtaIn->at(i))  > 0.007 ) continue;
	      if ( fabs(electrons_deltaPhiIn->at(i))  > 0.09 ) continue;
	      if ( electrons_hOverE->at(i)      > 0.05 ) continue;
	      
	      if ( electrons_sigmaIetaIeta->at(i) > 0.03 ) continue;
	      
	      //iso variables
	      bool iso_EM_HAD = ( (eleEt < 50.)  && ( (electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i) ) < 2.5 ) ) ||
		                ( (eleEt >= 50.) && ( (electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i) ) < ( 2.+0.03*(eleEt-50.) ) ) );
	      bool iso_HAD    = ( electrons_hadIsoR03_depth2->at(i) < 0.5 );
	      bool iso_tk     = ( electrons_tkIsoR03->at(i) < 15. );

	      if (nonIso_ == 0)
		{
		  //if ( eleEt <  50. && (electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i) ) > 2.5 ) continue;
		  //if ( eleEt >= 50. && (electrons_emIsoR03->at(i) + electrons_hadIsoR03_depth1->at(i) ) >  ( 2.+0.03*(eleEt-50.) ) )  continue;
 		  //if ( electrons_hadIsoR03_depth2->at(i) > 0.5 ) continue;
 		  //if ( electrons_tkIsoR03->at(i) > 15. ) continue;
		  
		  if (!iso_EM_HAD) continue;
		  if (!iso_HAD) continue;
		  if (!iso_tk) continue;
		}
	      else
		{
		  //not ISO in at least one of the isolation def
		  if (iso_EM_HAD) continue;
                  if (iso_HAD) continue;
                  if (iso_tk) continue;
		}
	    }
	  //================================================================================================================================================
	  //--------------------
	  //---- VBTF eleId ----
	  //--------------------

	  if (eleIsEB == true && doVBTFeleId_ == 1) //EB VBTF
	    {
	      //SOTTRARRE IL PU COME QUI: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011#ID_and_Conversion_Rejection
	      float combIso = ( tkIso + std::max(0., emIso - 1.) + hadIso) / elePt ;

	      //if(  combIso      > 0.07 ) continue;

	      //riproduco i tagli @ HLT level
	      if( tkIso / elePt    > 0.09 ) continue;
	      if( emIso / elePt    > 0.07 ) continue;
	      if( hadIso/ elePt    > 0.10 ) continue;

	      // eleId VBTF 2011 80% 
	      if( sigmaIetaIeta > 0.01  ) continue;
	      if(        DphiIn > 0.06 ) continue;
	      if(        DetaIn > 0.004 ) continue;
	      if(        HOverE > 0.04 ) continue;
	      
	    }
	  if (eleIsEB == false && doVBTFeleId_ == 1) //EB VBTF
	    {
	      //SOTTRARRE IL PU COME QUI: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011#ID_and_Conversion_Rejection
	      float combIso = ( tkIso + emIso + hadIso) / elePt ;

	      //if(  combIso      > 0.06 ) continue;

	      //riproduco i tagli @ HLT level
	      if( tkIso / elePt    > 0.04 ) continue;
	      if( emIso / elePt    > 0.05 ) continue;
	      if( hadIso/ elePt    > 0.025 ) continue;

	      // eleId VBTF 80%
	      if( sigmaIetaIeta > 0.03  ) continue;
	      if(        DphiIn > 0.03 ) continue;
	      if(        DetaIn > 0.007 ) continue;
	      if(        HOverE > 0.025 ) continue;
	    }
	  
	  //================================================================================================================================================
	  
	  
	  nGoodElectrons++;
	  chosenEle = i;
	  
	}// end loop over ele cand
      
      
      
      //================================
      //======== step 8: eleId =========
      //================================
      if( nGoodElectrons == 0 ) continue;
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "6) eleId";
      

      //============================
      //=== step 9: only one ele ===
      //============================
      if ( nGoodElectrons != 1 ) continue;
      if (electrons->at(chosenEle).pt() < minEleEt_fir_) continue;

      vars.selectIt_ele = chosenEle;
      SetElectronVariables(vars, reader);
      SetMetVariables(vars, reader);

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "7) single electron";


      cloneTrees[step] -> Fill();        

      //===============================
      //=== step 10: ele - met btob ===
      //===============================
      
      if (vars.eleMet_Dphi < minEleMetDPhi_) continue; //original
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "8) ele - met btob";

      cloneTrees[step] -> Fill();        
      
      //==============================
      //=== step 12: met selection ===
      //==============================

      if (vars.ele.pt()/vars.met.pt() < minEtOverMet_ || vars.ele.pt()/vars.met.pt() > maxEtOverMet_) continue;  //original
      //if (met < minMET_) continue;  //original
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "9) met selection";

      cloneTrees[step] -> Fill();    

    } // end loop over entries
  

  for(step = firstStep; step <= lastStep; ++step)
  {
    cloneTrees[step] -> Write(); //Write the tree
  }

  DeleteWprimeVariables(vars); //Close the TFile

  //======================
  //==== write histos ====
  //======================
  TFile fout(outFile_.c_str(), "UPDATE"); //riapro lo stesso TFile in UPDATE
  fout.cd();

  for(int step = 1; step <= nSteps; ++step)
    {
      events -> SetBinContent(step, stepEvents[step]);
      events -> GetXaxis() -> SetBinLabel(step, stepName[step].c_str());

    }

  events -> Write();

  //================
  //==== saving ====
  //================


  fout.Close(); //richiudo il file

  return 0;
}
