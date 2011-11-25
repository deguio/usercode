//---- CMSSW includes ----
#include "treeReader.h"
#include "WprimeTreeVariables.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "hChain.h"

#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "Math/GenVector/VectorUtil.h"
#include "TRandom3.h"
#include <time.h>
#include <sstream>

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
  int         evalFakeRate_  = gConfigParser ->  readIntOption("Options::evalFakeRate");
  int         useFakeRate_   = gConfigParser ->  readIntOption("Options::useFakeRate");
  int         doVBTFeleId_   = gConfigParser ->  readIntOption("Options::doVBTFeleId");
  int         doRecoilCorrection_   = gConfigParser ->  readIntOption("Options::doRecoilCorrection");
  float       crossSection_  = gConfigParser -> readFloatOption("Options::crossSection");

  std::string outFile_       = gConfigParser -> readStringOption("Output::OutputFile");  

  float minEleEt_fir_            = gConfigParser -> readFloatOption("Cuts::minEleEtFir");  
  float minEleEt_sec_            = gConfigParser -> readFloatOption("Cuts::minEleEtSec");  
  float minEtOverMet_        = gConfigParser -> readFloatOption("Cuts::minEtOverMet");  
  float maxEtOverMet_        = gConfigParser -> readFloatOption("Cuts::maxEtOverMet");  
  float minEleMetDPhi_       = gConfigParser -> readFloatOption("Cuts::minEleMetDPhi");  

  //---------------------------
  std::map<int, int> totalEvents           = GetTotalEvents("AllEvents/totalEvents", inputFileList_.c_str()); 
  //std::map<int, int> preselectionsEvents   = GetTotalEvents("PhotonsFilterEvents/totalEvents", inputFileList_.c_str()); 
  //std::map<int, int> preselectionsEvents   = GetTotalEvents("ElectronsFilterEvents/totalEvents", inputFileList_.c_str()); 
  

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

  //===============================================
  //==== define additional histograms and tree ====
  //===============================================

  TH1F* numerator_EB = new TH1F("numerator_EB","numerator_EB", 150, 0., 1500.);
  TH1F* numerator_EE = new TH1F("numerator_EE","numerator_EE", 150, 0., 1500.);

  TH1F* denominator_EB = new TH1F("denominator_EB","denominator_EB", 150, 0., 1500.);
  TH1F* denominator_EE = new TH1F("denominator_EE","denominator_EE", 150, 0., 1500.);

  // TGraphAsymmErrors* fakeRate_EB = new TGraphAsymmErrors();
  // TGraphAsymmErrors* fakeRate_EE = new TGraphAsymmErrors();


  //fake rate fits HEEP
  //TF1* func_fakeRateEB = new TF1("fakeRateEB","0.005633 + 0.00007476*x");
  //TF1* func_fakeRateEE = new TF1("fakeRateEE","0.02251  + 0.00002687*x");

  //fake rate fits WP80 May10
  //TF1* func_fakeRateEB = new TF1("fakeRateEB","0.002077 + 0.0001818*x");
  //TF1* func_fakeRateEE = new TF1("fakeRateEE","0.003678 + 0.0002549*x");

  //fake rate fits WP80 05Jul Photon20 HLT
  TF1* func_fakeRateEB = new TF1("fakeRateEB","0.002943 + 0.0001625*x");
  TF1* func_fakeRateEE = new TF1("fakeRateEE","0.006 + 0.000189*x");



  //===================================
  //==== define variable container ====
  //===================================
  WprimeVariables vars;
  InitializeTree(vars, outFile_);

  // define clone trees
  int firstStep = 6;
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

  //no HLT requirement on the MC. No loss of efficiency due to the tracking

  std::vector<std::string> HLTPathNamesMC;

  // //trigger Spring11
  // HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v1");
  // HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v2");
  // HLTPathNamesMC.push_back("HLT_Ele22_SW_TighterEleId_L1R_v3");

  // //trigger Fall10
  // HLTPathNamesMC.push_back("HLT_Ele12_SW_TightEleId_L1R");



  //------------------------------------
  //preoccuparsi di trigger e prescales

  std::vector<std::string> HLTPathNamesDATA;
  
  //devono essere esclusivi se non voglio complicarmi la vita con i prescales
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v1");
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v2");
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v3");  //presc 500
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v4");  //presc 1000
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v5");  //presc 2000
  HLTPathNamesDATA.push_back("HLT_Photon30_CaloIdVL_v6");  //presc 2000 coperto fino a fine 05Jul

  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v1");
  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v2");
  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v3");
  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v4");
  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v5");
  HLTPathNamesDATA.push_back("HLT_Photon90_CaloIdVL_v6");

  //============================
  //==== preselection steps ====
  //============================
  
  int step = 1;
  stepEvents[step] = totalEvents[1];
  stepName[step] = "1) total events";
  
  // step = 2;
  // stepEvents[step] = preselectionsEvents[1];
  // stepName[step] = "2) preselections";
  
  
  
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
      std::vector<ROOT::Math::XYZVector>* electrons_p_atVtx = reader.Get3V("electrons_p_atVtx");
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

      std::vector<ROOT::Math::XYZTVector>* photons = reader.Get4V("photons");

      std::vector<ROOT::Math::XYZTVector>* PFMet = reader.Get4V("PFMet");

      //===================
      //==== get entry ====
      //===================
      reader.GetEntry (entry) ;
      ClearWprimeVariables(vars);
      if(entry%100000 == 0) std::cout << "event n. " << entry << std::endl;

      step = 2 ;
      stepEvents[step] ++ ;
      stepName[step] = "after preselections";


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

      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "remove duplication";
      

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

      if (dataFlag_ == 0) SetPUVariables(vars, reader);   //per i dati mi va bene l'inizializzazione a -1
      SetPVVariables(vars, reader);
      SetMetVariables(vars, reader, doRecoilCorrection_);


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
      stepName[step] = "JSON";
      
      //================================
      //==== step 6: HLT selection  ====
      //================================

      //std::cout << "---> New event" << std::endl;

      int prescale = 0;
      if(dataFlag_ == 0)
	{
	  // for(unsigned int HLTIt = 0; HLTIt < HLTPathNamesMC.size(); ++HLTIt)
	  //   {
	  //     if ( GetPrescale(reader, HLTPathNamesMC.at(HLTIt)) > 0)
	  // 	prescale = GetPrescale(reader, HLTPathNamesMC.at(HLTIt));
	  //   }
	  prescale = 1;
   	}
      else
	{
	  for(unsigned int HLTIt = 0; HLTIt < HLTPathNamesDATA.size(); ++HLTIt)
	    {
	      if ( GetPrescale(reader, HLTPathNamesDATA.at(HLTIt)) > 0)
		prescale = GetPrescale(reader, HLTPathNamesDATA.at(HLTIt));
	    }
	}
      
      //std::cout << "prescale = " << prescale << std::endl;
      if( prescale == 0 ) continue;
      vars.hltPrescale = prescale;
      
      step ++ ;
      stepEvents[step] ++ ;
      stepName[step] = "HLT";
      

      //=============================
      //==== step 7: HBHE noise  ====
      //=============================
      if (dataFlag_ == 1 && reader.GetInt("HCAL_noise")->at(0) == 0 ) continue;
      
      step ++;
      stepEvents[step] ++;
      stepName[step] = "HBHE noise";

      cloneTrees[step] -> Fill();
     
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      //============================
      //==== loop over pho cand ====
      //============================
      int nGoodPhotons = 0;
      int chosenPho = 0;

      //APPLICARE UNA MINIMA DI ISOLAMENTO E H/E (v. definizione triggers)
      //v. qui per CaloIdVL @ HLT level: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaWorkingPointsv3
      for (unsigned int i = 0; i < photons->size(); ++i)
  	{
  	  double phoEt  = photons->at(i).pt();
	  if (phoEt < 25.) continue;            //minimum pt requirement on the photons

	  ++nGoodPhotons;
	  chosenPho = i;
	} // end loop over photons

      //=======================
      //==== exactly 1 pho ====
      //=======================
      if (nGoodPhotons != 1) continue;  //eventi a singolo photone
      if (photons->at(chosenPho).pt() < 35.) continue;

      double phoEt  = photons->at(chosenPho).pt();
      double phoEta = photons->at(chosenPho).eta();
      double phoPhi = photons->at(chosenPho).phi();

      vars.selectIt_pho = chosenPho;
      SetPhotonVariables(vars, reader);
      SetMetVariables(vars, reader, doRecoilCorrection_);  //second time to set photonBranches  (deltaPhi, Mt, ...)

      step ++ ;  //step 7
      stepEvents[step] ++ ;
      stepName[step] = "exactly 1 good pho";

      cloneTrees[step] -> Fill();


      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      if (evalFakeRate_ == 1)
	{
	  //=======================
	  //==== MET selection ====
	  //=======================
	  
	  if ( PFMet->at(0).pt() > 20.) continue;  //reduce the contamination from W (ma non negli endcap)
	  
	  //fill denominator
	  if (fabs(phoEta) < EtaCutEB) denominator_EB->Fill(phoEt);
	  else denominator_EE->Fill(phoEt);

	  step ++ ;
	  stepEvents[step] ++ ;
	  stepName[step] = "MET < 20 GeV";
	  
	  cloneTrees[step] -> Fill();


	  //========================
	  //==== match with ele ====
	  //========================
	  
	  double deltaR_min = 10000;
	  int eleIndex = -1;
	  for (unsigned int i = 0; i < electrons->size(); ++i)
	    {
	      double eleEta = electrons_p_atVtx->at(i).eta();
	      double elePhi = electrons_p_atVtx->at(i).phi();
	      
	      double deltaR_ele_pho = deltaR(eleEta, elePhi, phoEta, phoPhi);
	      
	      if (deltaR_ele_pho < deltaR_min) 
		{
		  deltaR_min = deltaR_ele_pho;
		  eleIndex = i;
		}
	    }
	  
	  if (deltaR_min > 0.1) continue;
	  
	  
	  //ho scelto l'elettrone matchato al photon
	  // double eleEt  = electrons->at(eleIndex).pt(); //FIXME: VOGLIO ESC*COS(THETATRACCIA)
	  // double elePt  = electrons->at(eleIndex).pt(); //FIXME: VOGLIO ESC*COS(THETATRACCIA)

	  // double elePhi = electrons->at(eleIndex).phi(); //FIXME: VOGLIO ETA E PHI DELLA TRACCIA
	  // double eleEta = electrons->at(eleIndex).eta(); //FIXME: VOGLIO ETA E PHI DELLA TRACCIA

  	  double eleEt  = electrons_p_atVtx->at(eleIndex).Rho()/electrons_p_atVtx->at(eleIndex).R()*electrons_eSC->at(eleIndex);
  	  double elePhi = electrons_p_atVtx->at(eleIndex).phi();
  	  double eleEta = electrons_p_atVtx->at(eleIndex).eta();
	  
	  bool eleIsEB = electrons_isEB->at(eleIndex);
	  
	  float e2x5_overE5x5 = electrons_e2x5Max->at(eleIndex)/electrons_e5x5->at(eleIndex);
	  float e1x5_overE5x5 = electrons_e1x5->at(eleIndex)/electrons_e5x5->at(eleIndex);
	  
	  float tkIso  = electrons_tkIsoR03->at(eleIndex); 
	  float emIso  = electrons_emIsoR03->at(eleIndex);
	  float hadIso = electrons_hadIsoR03_depth1->at(eleIndex) + electrons_hadIsoR03_depth2->at(eleIndex);
	  
	  float sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta")->at(eleIndex);
	  float DphiIn = fabs(reader.GetFloat("electrons_deltaPhiIn")->at(eleIndex));
	  float DetaIn = fabs(reader.GetFloat("electrons_deltaEtaIn")->at(eleIndex));
	  float HOverE = reader.GetFloat("electrons_hOverE")->at(eleIndex);
	  
	  //salvo le info dell'ele matchato con pho
	  vars.selectIt_ele = eleIndex;
	  SetElectronVariables(vars, reader);
	  SetMetVariables(vars, reader, doRecoilCorrection_); //third time to set electron deltaPhi and Mt
	  
	  step ++ ;
	  stepEvents[step] ++ ;
	  stepName[step] = "ele-pho matching";
	  
	  cloneTrees[step] -> Fill();
	  
	  //================================================================================================================================================
	  //--------------------
	  //---- HEEP eleID ----
	  //--------------------
	  if (eleIsEB == true && doVBTFeleId_ == 0) //EB HEEP
	    {
	      if ( electrons_ecalDrivenSeed->at(eleIndex) != 1 ) continue;
	      
	      if ( fabs(electrons_deltaEtaIn->at(eleIndex)) > 0.005 ) continue;
	      if ( fabs(electrons_deltaPhiIn->at(eleIndex)) > 0.09 ) continue;
	      if ( electrons_hOverE->at(eleIndex) > 0.05 ) continue;
	      
	      if ( e2x5_overE5x5 < 0.94 && e1x5_overE5x5 < 0.83) continue;
	      
	      //iso variables
	      bool iso_EM_HAD = ( (electrons_emIsoR03->at(eleIndex) + electrons_hadIsoR03_depth1->at(eleIndex)) < (2.+0.03*eleEt) );
	      bool iso_tk     = ( electrons_tkIsoR03->at(eleIndex) < 7.5 );
	      
	      if (!iso_EM_HAD) continue;
	      if (!iso_tk) continue;
	    }
	  if (eleIsEB == false && doVBTFeleId_ == 0) //EE HEEP
	    {
	      if ( electrons_ecalDrivenSeed->at(eleIndex) != 1 ) continue;

	      if ( fabs(electrons_deltaEtaIn->at(eleIndex))  > 0.007 ) continue;
	      if ( fabs(electrons_deltaPhiIn->at(eleIndex))  > 0.09 ) continue;
	      if ( electrons_hOverE->at(eleIndex)      > 0.05 ) continue;
	      
	      if ( electrons_sigmaIetaIeta->at(eleIndex) > 0.03 ) continue;
	      
	      //iso variables
	      bool iso_EM_HAD = ( (eleEt < 50.)  && ( (electrons_emIsoR03->at(eleIndex) + electrons_hadIsoR03_depth1->at(eleIndex) ) < 2.5 ) ) ||
		                ( (eleEt >= 50.) && ( (electrons_emIsoR03->at(eleIndex) + electrons_hadIsoR03_depth1->at(eleIndex) ) < ( 2.+0.03*(eleEt-50.) ) ) );
	      bool iso_HAD    = ( electrons_hadIsoR03_depth2->at(eleIndex) < 0.5 );
	      bool iso_tk     = ( electrons_tkIsoR03->at(eleIndex) < 15. );
	      
	      if (!iso_EM_HAD) continue;
	      if (!iso_HAD) continue;
	      if (!iso_tk) continue;
	      
	    }
	  
	  //================================================================================================================================================
	  //--------------------
	  //---- VBTF eleId ----
	  //--------------------

	  if (eleIsEB == true && doVBTFeleId_ == 1) //EB VBTF
	    {
	      //SOTTRARRE IL PU COME QUI: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011#ID_and_Conversion_Rejection
	      float combIso = ( tkIso + std::max(0., emIso - 1.) + hadIso) / eleEt ;

	      //if(  combIso      > 0.07 ) continue;

	      //riproduco i tagli @ HLT level
	      if( tkIso / eleEt    > 0.09 ) continue;
	      if( emIso / eleEt    > 0.07 ) continue;
	      if( hadIso/ eleEt    > 0.10 ) continue;

	      // eleId VBTF 2011 80% 
	      if( sigmaIetaIeta > 0.01  ) continue;
	      if(        DphiIn > 0.06 ) continue;
	      if(        DetaIn > 0.004 ) continue;
	      if(        HOverE > 0.04 ) continue;
	      
	    }
	  if (eleIsEB == false && doVBTFeleId_ == 1) //EB VBTF
	    {
	      //SOTTRARRE IL PU COME QUI: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SimpleCutBasedEleID2011#ID_and_Conversion_Rejection
	      float combIso = ( tkIso + emIso + hadIso) / eleEt ;

	      //if(  combIso      > 0.06 ) continue;

	      //riproduco i tagli @ HLT level
	      if( tkIso / eleEt    > 0.04 ) continue;
	      if( emIso / eleEt    > 0.05 ) continue;
	      if( hadIso/ eleEt    > 0.025 ) continue;

	      // eleId VBTF 80%
	      if( sigmaIetaIeta > 0.03  ) continue;
	      if(        DphiIn > 0.03 ) continue;
	      if(        DetaIn > 0.007 ) continue;
	      if(        HOverE > 0.025 ) continue;
	    }
	  
	  //================================================================================================================================================
	  
	  
	  step ++ ;
	  stepEvents[step] ++ ;
	  stepName[step] = "pass the eleId";
	  
	  cloneTrees[step] -> Fill();

	  //fill numerator
	  if (fabs(phoEta) < EtaCutEB) numerator_EB->Fill(phoEt);
	  else numerator_EE->Fill(phoEt);	
	}//end eval fake rate
      
      //-----------------------------------------------------------------------------------------------------------------------------------------------------
      if (useFakeRate_ == 1)
	{//step 8 = step 7

	  //FIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXME
	  if (fabs(phoEta) < EtaCutEB) vars.pho_weight = func_fakeRateEB->Eval( vars.pho.Et() );
	  else vars.pho_weight = func_fakeRateEE->Eval( vars.pho.Et() );
	  
	  //set ele variables con pho info per PrintPlot  
	  vars.ele = vars.pho;
	  vars.p_ele = vars.p_pho;
	  vars.eleMet_Dphi = vars.phoMet_Dphi;
	  vars.eleMet_mt = vars.phoMet_mt;
	  //FIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXMEFIXME

	  step = 8 ;
	  stepEvents[step] ++ ;
	  stepName[step] = "exactly one good pho";
	  
	  cloneTrees[step] -> Fill();

	  //===============================
	  //=== step 9: ele - met btob ===
	  //===============================	  
	  if ( deltaPhi(vars.pho.phi(), vars.met.phi()) < minEleMetDPhi_ ) continue; //original
	  
	  step = 9 ;
	  stepEvents[step] ++ ;
	  stepName[step] = "deltaPhi ele-MET";
	  
	  cloneTrees[step] -> Fill();
	  
	  //=================================
	  //=== step 10: et/met selection ===
	  //=================================
	  if (vars.pho.pt()/vars.met.pt() < minEtOverMet_ || vars.pho.pt()/vars.met.pt() > maxEtOverMet_) continue;  //original

	  step = 10 ;
	  stepEvents[step] ++ ;
	  stepName[step] = "et/MET";

	  cloneTrees[step] -> Fill();

	}//end use fake rate
      


      
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

  if ( evalFakeRate_ ==1 )
    {
      numerator_EB->Write();
      numerator_EE->Write();
      
      denominator_EB->Write();
      denominator_EE->Write();
      
      // fakeRate_EB->BayesDivide(numerator_EB, denominator_EB);
      // fakeRate_EE->BayesDivide(numerator_EE, denominator_EE);
      
      // fakeRate_EB->Write("EB");
      // fakeRate_EE->Write("EE");
    }

  //================
  //==== saving ====
  //================


  fout.Close(); //richiudo il file

  return 0;
}


int AcceptHLTPath(treeReader& reader, const std::string& HLTPathName)
{
  int prescale = 0;
  
  std::vector<std::string> HLT_names = *(reader.GetString("HLT_Names"));
  for(unsigned int HLTIt = 0; HLTIt < HLT_names.size(); ++HLTIt)
    if( (reader.GetString("HLT_Names")->at(HLTIt) == HLTPathName) && (reader.GetFloat("HLT_Accept")->at(HLTIt) == 1) )
      {
	prescale = reader.GetInt("HLT_Prescale")->at(HLTIt) ;
	// std::cout << "chiamata" << std::endl;
	//std::cout << HLTPathName << " "<<  reader.GetInt("HLT_Prescale")->at(HLTIt) << std::endl;
	// std::cout << reader.GetFloat("HLT_Accept")->at(HLTIt) << std::endl;
      }
  
  return prescale;
}
