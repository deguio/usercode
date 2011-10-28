
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

  std::string inputFileList = gConfigParser -> readStringOption("Input::inputFileList");

  std::string treeName  = gConfigParser -> readStringOption("Input::treeName");
  
  std::string tmvaMethod    = gConfigParser -> readStringOption("Input::tmvaMethod");
  std::string tmvaWeights   = gConfigParser -> readStringOption("Input::tmvaWeights");
  
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  

  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  int isData   = gConfigParser -> readIntOption("Options::isData");

  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  int addConversionToMva = gConfigParser -> readIntOption("Options::addConversionToMva");
  int useMvaRanking      = gConfigParser -> readIntOption("Options::useMvaRanking");
  

  int useWeights      = gConfigParser -> readIntOption("Options::useWeights");
  int poissonWeights  = gConfigParser -> readIntOption("Options::poissonWeights");
  int nAvePU          = gConfigParser -> readIntOption("Options::nAvePU");
  
  std::string puweightsFileName = gConfigParser -> readStringOption("Options::puweightsFileName");  

  int useJSON      = gConfigParser -> readIntOption("Options::useJSON");
  std::string jsonFileName = gConfigParser -> readStringOption("Options::jsonFileName");  

  //******* Get run/LS map from JSON file *******
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  if ( isData && useJSON ) {
    std::cout << ">>> Getting GOOD  run/LS from JSON file" << std::endl;
    jsonMap = readJSONFile(jsonFileName);
    std::cout << std::endl;
    std::cout << std::endl;
  }

  
  //****** Get weights for MC ****** 
  float nmax;
  float w[50];
  TRandom *gRandom = new TRandom();
  
  if (useWeights){
    
    TFile weightsFile(puweightsFileName.c_str(),"READ");  
    TH1F* hweights;
     
    if ( poissonWeights ){
      std::cout << "N ave PU for Poisson PU reweighting : " << nAvePU << std::endl;
      char hname[100];
      sprintf(hname,"hwmc%d",nAvePU);
      hweights = (TH1F*)weightsFile.Get(hname);
    }
    else {
      hweights = (TH1F*)weightsFile.Get("hweights");
    }
    nmax = hweights ->GetMaximum();
    std::cout << " Max weight " << nmax << std::endl;
    
    for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
      w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
    }
    weightsFile.Close();
  }
  
  
  //****** Parameters for vertex finding algo ****** 
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 0.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;
  vtxAlgoParams_.removeTracksInCone = 1;
  vtxAlgoParams_.coneSize           = 0.05;
  
  //--- sumpt2 ordering
  vector<string> ranksumpt2_;
  ranksumpt2_.push_back("logsumpt2");  

  //--- ranking product variables
  vector<string> rankVariables_;
  // variables order matters to resolve ties
  rankVariables_.push_back("ptbal"), rankVariables_.push_back("ptasym"),rankVariables_.push_back("logsumpt2");  
  
  //--- tmva variables
  vector<string> tmvaPerVtxVariables_;
  tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"),tmvaPerVtxVariables_.push_back("logsumpt2");  
  if( addConversionToMva ) {
    tmvaPerVtxVariables_.push_back("limPullToConv");
    tmvaPerVtxVariables_.push_back("nConv");
  } 
  
  //--- book TMVA
  TMVA::Reader *tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
  HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaPerVtxVariables_ );
  tmvaReader_->BookMVA( tmvaMethod.c_str(), tmvaWeights.c_str() );


  //****** BOOK OUTPUT HISTOGRAMS ******

  TH1F PtAll("PtAll","Pt of boson all",80,0,400);
  TH1F PtGood("PtGood","Pt of boson good",80,0,400);
  TH1F PtGood_BDT("PtGood_BDT","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK("PtGood_RANK","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EBEB("PtAll_EBEB","Pt of boson all",80,0,400);
  TH1F PtGood_EBEB("PtGood_EBEB","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EBEB("PtGood_BDT_EBEB","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EBEB("PtGood_RANK_EBEB","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EBEE("PtAll_EBEE","Pt of boson all",80,0,400);
  TH1F PtGood_EBEE("PtGood_EBEE","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EBEE("PtGood_BDT_EBEE","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EBEE("PtGood_RANK_EBEE","Pt of boson good (RANKING)",80,0,400);

  TH1F PtAll_EEEE("PtAll_EEEE","Pt of boson all",80,0,400);
  TH1F PtGood_EEEE("PtGood_EEEE","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_EEEE("PtGood_BDT_EEEE","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_EEEE("PtGood_RANK_EEEE","Pt of boson good (RANKING)",80,0,400);

  TH1F EtaAll("EtaAll","Eta of max SC",50,-5,5);
  TH1F EtaGood("EtaGood","Eta of max SC good",50,-5,5);
  TH1F EtaGood_BDT("EtaGood_BDT","Eta of max SC good (BDT)",50,-5,5);
  TH1F EtaGood_RANK("EtaGood_RANK","Eta of max SC good (RANKING)",50,-5,5);
  
  TH1F NvtAll("NvtAll","number of PV all",50,0,50);
  TH1F NvtGood("NvtGood","number of PV good",50,0,50);
  TH1F NvtGood_BDT("NvtGood_BDT","number of PV good (BDT)",50,0,50);
  TH1F NvtGood_RANK("NvtGood_RANK","number of PV good (RANKING)",50,0,50);

  TH1F NpuAll("NpuAll","number of PV all",50,0,50);
  TH1F NpuGood("NpuGood","number of PV good",50,0,50);
  TH1F NpuGood_BDT("NpuGood_BDT","number of PV good (BDT)",50,0,50);
  TH1F NpuGood_RANK("NpuGood_RANK","number of PV good (RANKING)",50,0,50);

  TH1F PtGood_matchedClosest("PtGood_matchedClosest","Pt of boson good",80,0,400);
  TH1F PtGood_BDT_matchedClosest("PtGood_BDT_matchedClosest","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK_matchedClosest("PtGood_RANK_matchedClosest","Pt of boson good (RANKING)",80,0,400);

  TH1F EtaGood_matchedClosest("EtaGood_matchedClosest","Eta of max SC good",50,-5,5);
  TH1F EtaGood_BDT_matchedClosest("EtaGood_BDT_matchedClosest","Eta of max SC good (BDT)",50,-5,5);
  TH1F EtaGood_RANK_matchedClosest("EtaGood_RANK_matchedClosest","Eta of max SC good (RANKING)",50,-5,5);
  
  TH1F NvtGood_matchedClosest("NvtGood_matchedClosest","number of PV good",50,0,50);
  TH1F NvtGood_BDT_matchedClosest("NvtGood_BDT_matchedClosest","number of PV good (BDT)",50,0,50);
  TH1F NvtGood_RANK_matchedClosest("NvtGood_RANK_matchedClosest","number of PV good (RANKING)",50,0,50);

  TH1F NpuGood_matchedClosest("NpuGood_matchedClosest","number of PV good",50,0,50);
  TH1F NpuGood_BDT_matchedClosest("NpuGood_BDT_matchedClosest","number of PV good (BDT)",50,0,50);
  TH1F NpuGood_RANK_matchedClosest("NpuGood_RANK_matchedClosest","number of PV good (RANKING)",50,0,50);

  TH2F hdist("hdist"," hdist",80,0,200,400,-10,10);
  TH2F hdiff_dZ_muons("hdiff_dZ_muons","hdiff_dZ_muons",80,0,200,400,-10,10);
  TH2F hdiff_dZ_electrons("hdiff_dZ_electrons","hdiff_dZ_electrons",80,0,200,400,-10,10);

  TH1F BDToutput("BDToutput","BDT output",500,-1,1);
  TH1F BDToutput_sig("BDToutput_sig","BDT output - signal vertices",500,-1,1);
  TH1F BDToutput_bkg("BDToutput_bkg","BDT output - background",500,-1,1);
  

  TH1F pt2h("pt2h","pt2 H",500,0,500);
  TH1F pt2bkg("pt2bkg","pt2 bkg",500,0,500);


  TH2F * hAcceptedLumis = new TH2F("hAcceptedLumis","hAcceptedLumis",20000, 160000, 180000, 10000, 0, 10000);


  float ww = 1;
  float r9cut = 0.93;
  
  //****** LOAD TREE ******
  TChain* chain = new TChain(treeName.c_str());
  FillChain(*chain, inputFileList.c_str());
  treeReader reader((TTree*)(chain));
  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  //****** Start loop over entries ******
  int runId, lumiId;

  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      //*** filter bad runs/lumis
      runId = reader.GetInt("runId")->at(0);
      lumiId = reader.GetInt("lumiId")->at(0);
      
      bool skipEvent = false;
      if( isData && useJSON ){
	if(AcceptEventByRunAndLumiSection(runId,lumiId,jsonMap) == false)
	  skipEvent = true;
      }
      if( skipEvent == true ) continue;
      hAcceptedLumis -> Fill(runId, lumiId);

      //*** pu weights
      std::vector<float>*PU_z ;
      std::vector<int>* mc_PUit_NumInteractions; 
      int npu ;

      if ( !isData ){
      	mc_PUit_NumInteractions  = reader.GetInt("mc_PUit_NumInteractions");
	npu = mc_PUit_NumInteractions->at(0);
      
	//--- use weights 
	if (useWeights){
	  float myrnd = gRandom->Uniform(0,nmax);
	  if (myrnd > w[npu]) continue;
	}
      }

     
      //*** setup common branches ***
      std::vector<int>* PV_nTracks;
      std::vector<float>* PV_z;
      std::vector<float>* PV_d0;
      std::vector<ROOT::Math::XYZVector>* PVtracks;
      std::vector<int>* PVtracks_PVindex;
      std::vector<int>* tracks_PVindex;
      std::vector<float>* tracks_dz ; //?
      std::vector<float>* tracks_dz_PV ; //?
      std::vector<float>* tracks_dxy_PV ; //?
      std::vector<ROOT::Math::XYZTVector>* sc; // supercluster
      
      int accept = 0;
      int indpho1 = -100;
      int indpho2 = -100;
	
      ROOT::Math::XYZTVector sum2pho;
      float etaMaxSC ;
      float TrueVertex_Z;

      //*** selections for Hgg ***
      if (isHiggs){
	std::vector<ROOT::Math::XYZVector>* mc_H_vertex = reader.Get3V("mc_H_vertex");
	std::vector<ROOT::Math::XYZTVector>* mcV1 = reader.Get4V("mcV1");
	std::vector<ROOT::Math::XYZTVector>* mcV2 = reader.Get4V("mcV2");
	std::vector<ROOT::Math::XYZTVector>* photons = reader.Get4V("photons");
	std::vector<float>* photons_r9 = reader.GetFloat("photons_r9");
	
	if (mc_H_vertex->size() != 1) continue;
	
	PV_nTracks       = reader.GetInt("PV_nTracks");
	PV_z             = reader.GetFloat("PV_z");
	PV_d0            = reader.GetFloat("PV_d0");
	PVtracks         = reader.Get3V("PVtracks");
	PVtracks_PVindex = reader.GetInt("PVtracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	sc               = reader.Get4V("photons_SC");

	hggSelection(mcV1, mcV2, photons, sc, photons_r9, accept, indpho1, indpho2);
	if (!accept) continue;

	if ( photons_r9->at(indpho1) < r9cut ) continue;
	if ( photons_r9->at(indpho2) < r9cut ) continue;
	
	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	TrueVertex_Z = mc_H_vertex->at(0).Z();
      }// end Hgg selection


      //*** selections for Zee ***
      if (isZee){
	std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
	// std::vector<float>* eleid = reader.GetFloat("simpleEleId95cIso");
	// no eleID95 available for data 2011 --> compute it by hand
	std::vector<float>* eleid = new std::vector<float>;
	eleid->clear();
	if ( electrons->size() < 2) continue; 
	
	for (unsigned int iele = 0; iele < electrons->size(); iele++){
	  
	  float pt = electrons->at(iele).pt();
	  float tkIso   = reader.GetFloat("electrons_tkIsoR03")->at(iele);
	  float emIso   = reader.GetFloat("electrons_emIsoR03")->at(iele);
	  float hadIso  = reader.GetFloat("electrons_hadIsoR03_depth1")->at(iele) + reader.GetFloat("electrons_hadIsoR03_depth2")->at(iele);
	  float combIso = tkIso + emIso + hadIso;
	
	  int isEB = reader.GetInt("electrons_isEB")->at(iele);
	  float sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta")->at(iele);
	  float DetaIn        = reader.GetFloat("electrons_deltaEtaIn")->at(iele);
	  float DphiIn        = reader.GetFloat("electrons_deltaPhiIn")->at(iele);
	  float HOverE        = reader.GetFloat("electrons_hOverE")->at(iele);
	  int mishits         = reader.GetInt("electrons_mishits")->at(iele);

	  float id = eleId95 ( pt, tkIso, emIso, hadIso, combIso, isEB, sigmaIetaIeta, DetaIn, DphiIn, HOverE, mishits);	
	  id *= 7.; // to emulate simpleEleId95cIso

	  eleid->push_back( id );
	}

	std::vector<float>* electrons_dz_PV_noEle = reader.GetFloat("electrons_dz_PV_noEle");

	PV_nTracks       = reader.GetInt("PV_noEle_nTracks");
	PV_z             = reader.GetFloat("PV_noEle_z");
	PV_d0            = reader.GetFloat("PV_noEle_d0");
	PVtracks         = reader.Get3V("PVEleLessTracks");
	PVtracks_PVindex = reader.GetInt("PVEleLessTracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	//sc               = reader.Get4V("electrons_SC");
	sc               = reader.Get4V("electrons");
	
	zeeSelection(electrons, eleid, accept, indpho1, indpho2);
	if (!accept) continue;

	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = electrons->at(indpho1)+ electrons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (electrons_dz_PV_noEle->at(indpho1) + electrons_dz_PV_noEle->at(indpho2))/2.;

	hdiff_dZ_electrons.Fill( sum2pho.pt(), electrons_dz_PV_noEle->at(indpho1) - electrons_dz_PV_noEle->at(indpho2) );
      }

      //*** selections for Zmumu
      if ( isZmumu ){
	std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
	std::vector<int>* muons_global = reader.GetInt("muons_global");
	std::vector<int>* muons_tracker = reader.GetInt("muons_tracker");
	std::vector<float>* muons_dz_PV_noMuon = reader.GetFloat("muons_dz_PV_noMuon");
	
	PV_nTracks       = reader.GetInt("PV_noMuon_nTracks");
	PV_z             = reader.GetFloat("PV_noMuon_z");
	PV_d0            = reader.GetFloat("PV_noMuon_d0");
	PVtracks         = reader.Get3V("PVMuonLessTracks");
	PVtracks_PVindex = reader.GetInt("PVMuonLessTracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz        = reader.GetFloat("tracks_dz");
	sc               = reader.Get4V("muons");  // use muon info for SC
	
	zmumuSelection(muons,muons_global,muons_tracker, accept, indpho1, indpho2);
	
	if (!accept) continue;
	
	etaMaxSC = muons->at(indpho1).eta();
	sum2pho  = muons->at(indpho1)+ muons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indpho1) + muons_dz_PV_noMuon->at(indpho2))/2.;
	
	hdiff_dZ_muons.Fill( sum2pho.pt(), muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2) );
	
	if ( fabs(muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2))> 0.5) continue;
      }//Zmumu end
      
      // branches buffers
      int nvtx_;
      float  vtxx_[1000], vtxy_[1000], vtxz_[1000];
      int ntracks_;
      float tkpx_[1000], tkpy_[1000], tkpz_[1000], tkPtErr_[1000], tkWeight_[1000], 
	tkd0_[1000], tkd0Err_[1000], tkdz_[1000], tkdzErr_[1000];
      int tkVtxId_[1000];
      bool tkIsHighPurity_[1000];
      
      float phocalox_[100], phocaloy_[100], phocaloz_[100], phoen_[100];
      
      // set variables
      
      // vertices 
      nvtx_    = (int) PV_z->size();
      for ( int iv = 0; iv < nvtx_; iv++){
	vtxx_[iv] =  0;
	vtxy_[iv] =  0;
	vtxz_[iv] =  PV_z->at(iv) ;
      }
      
      // tracks
      ntracks_ = PVtracks->size();
      for (int itrk = 0; itrk <ntracks_; itrk++ ){
	tkpx_[itrk]    = PVtracks->at(itrk).X();
	tkpy_[itrk]    = PVtracks->at(itrk).Y();
	tkpz_[itrk]    = PVtracks->at(itrk).Z();
	tkPtErr_[itrk] = 0;
	tkVtxId_[itrk] = PVtracks_PVindex->at(itrk);
	
	tkWeight_[itrk]= 1.;
	tkd0_[itrk]    = 0;
	tkd0Err_[itrk] = 0;
	tkdz_[itrk]    = 0;
	tkdzErr_[itrk] = 0;
	tkIsHighPurity_[itrk]= 1.;
      }
      
      // photons
      for (int ipho = 0 ; ipho < sc->size(); ipho++){
	float px = sc->at(ipho).X();
	float py = sc->at(ipho).Y();
	float pz = sc->at(ipho).Z();
	float pt = sqrt ( px*px+py*py ); 
	float theta = 2*atan(exp(-sc->at(ipho).eta()) );
	float tantheta = tan(theta);

	if ( fabs(sc->at(ipho).eta()) < etaEB ) {
	  phocalox_[ipho] = R_ECAL*px/pt;
	  phocaloy_[ipho] = R_ECAL*py/pt;
	  phocaloz_[ipho] = R_ECAL/tantheta;
	} 

	if ( fabs(sc->at(ipho).eta()) > etaEE ) {
	  float r_endcap  = fabs(Z_ENDCAP * tantheta);
	  phocalox_[ipho] = r_endcap * px/pt;
	  phocaloy_[ipho] = r_endcap * py/pt;
	  if (pz > 0) phocaloz_[ipho] = Z_ENDCAP;
	  else phocaloz_[ipho] = -Z_ENDCAP;
	} 
	
	phoen_[ipho] = sc->at(ipho).E();
	
      }


      float eta1 = sc->at(indpho1).eta();
      float eta2 = sc->at(indpho2).eta();

      if ( (fabs(eta1) > etaEB && fabs(eta1) < etaEE) || fabs(eta1) > 2.5) continue;
      if ( (fabs(eta2) > etaEB && fabs(eta2) < etaEE) || fabs(eta2) > 2.5) continue;
   
      //*** set vertex info
      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
           
      //*** set photon info
      PhotonInfo pho1(indpho1, TVector3(phocalox_[indpho1],phocaloy_[indpho1],phocaloz_[indpho1]),phoen_[indpho1]); 
      PhotonInfo pho2(indpho2, TVector3(phocalox_[indpho2],phocaloy_[indpho2],phocaloz_[indpho2]),phoen_[indpho2]); 
     
      //*** vertex analyzer
      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);
      vAna.analyze(vinfo,pho1,pho2);
            
      //*** preselect vertices 
      std::vector<int> presel;
      for(int i=0; i<nvtx_; i++) {
	presel.push_back(i); 
      }
      vAna.preselection(presel);
      
      
      //*** Look if the H vertex matches one of the PV vertices
      float dmin = 10000;
      int iClosest = -1;
      for ( int uu = 0; uu < nvtx_; uu++){
	float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	if ( distt < dmin)   { dmin = distt; iClosest = uu; }
      }

       


      //*** NOW FILL HISTOGRAMS
      
      PtAll.Fill( sum2pho.pt(),ww );
      if (fabs(eta1) < etaEB && fabs(eta2) < etaEB)  PtAll_EBEB.Fill( sum2pho.pt(),ww );
      if (fabs(eta1) > etaEE && fabs(eta2) > etaEE)  PtAll_EEEE.Fill( sum2pho.pt(),ww );
      if ( (fabs(eta1) < etaEB && fabs(eta2) > etaEE) || (fabs(eta2) < etaEB && fabs(eta1) > etaEE))  PtAll_EBEE.Fill( sum2pho.pt(),ww );

      EtaAll.Fill( etaMaxSC ,ww);
      NvtAll.Fill( nvtx_,ww );
      if (!isData) NpuAll.Fill(npu,ww);

      //*** SUMPT2 CRITERION
      vector<int> ranksumpt2 = vAna.rankprod(ranksumpt2_);
      //-- matching 1cm
      if ( fabs( TrueVertex_Z - PV_z->at(ranksumpt2[0]) ) < 1.) {
	PtGood.Fill( sum2pho.pt(),ww );
	EtaGood.Fill( etaMaxSC ,ww);
	NvtGood.Fill( nvtx_ ,ww);
	if (!isData) NpuGood.Fill(npu,ww);
      }
      //-- matching closest vtx
      if ( iClosest == ranksumpt2[0]){
	PtGood_matchedClosest.Fill( sum2pho.pt(),ww );
	EtaGood_matchedClosest.Fill( etaMaxSC ,ww);
	NvtGood_matchedClosest.Fill( nvtx_ ,ww);
	if (!isData) NpuGood_matchedClosest.Fill(npu,ww);
      }        

   
      //*** RANKING PRODUCT
      vector<int> rankprod = vAna.rankprod(rankVariables_);
      //-- matching 1cm
      if ( fabs( TrueVertex_Z - PV_z->at(rankprod[0]) ) < 1.) {
	PtGood_RANK.Fill( sum2pho.pt(),ww );
	EtaGood_RANK.Fill( etaMaxSC ,ww);
	NvtGood_RANK.Fill( nvtx_,ww );
	if (!isData) NpuGood_RANK.Fill(npu,ww);
      }
      //-- matching closest vtx
      if ( iClosest == rankprod[0]){
	PtGood_RANK_matchedClosest.Fill( sum2pho.pt(),ww );
	EtaGood_RANK_matchedClosest.Fill( etaMaxSC ,ww);
	NvtGood_RANK_matchedClosest.Fill( nvtx_,ww );
	if (!isData) NpuGood_RANK_matchedClosest.Fill(npu,ww);
      }
      

      //*** BDT 
      vector<int> ranktmva = vAna.rank(*tmvaReader_,tmvaMethod);
      float vtxmva;
      for (int iv=0;iv<ranktmva.size();iv++) {
	vtxmva = vAna.mva(ranktmva[iv]);
	BDToutput.Fill( vtxmva );
	if ( iClosest == ranktmva[iv] ) BDToutput_sig.Fill( vtxmva );
	else BDToutput_bkg.Fill( vtxmva );
      }
      //-- matching 1cm
      if ( fabs( TrueVertex_Z - PV_z->at(ranktmva[0]) ) < 1.) {
	PtGood_BDT.Fill( sum2pho.pt(),ww );
	EtaGood_BDT.Fill( etaMaxSC ,ww);
	NvtGood_BDT.Fill( nvtx_ ,ww);
	if (!isData) NpuGood_BDT.Fill(npu,ww);
      }
      //-- matching closest vtx
      if ( iClosest == ranktmva[0]){
	PtGood_BDT_matchedClosest.Fill( sum2pho.pt(),ww );
	EtaGood_BDT_matchedClosest.Fill( etaMaxSC ,ww);
	NvtGood_BDT_matchedClosest.Fill( nvtx_,ww );
	if (!isData) NpuGood_BDT_matchedClosest.Fill(npu,ww);
      }
      
      
      
    }// end loop over entries

  std::cout << "END LOOP OVER ENTRIES" << std::endl;
  std::cout << "Saving histos on file ..." << std::endl;
  
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");

  hAcceptedLumis -> Write();
  
  PtAll.Write();
  PtGood.Write();
  PtGood_BDT.Write();
  PtGood_RANK.Write();
  
  EtaAll.Write();
  EtaGood.Write();
  EtaGood_BDT.Write();
  EtaGood_RANK.Write();
  
  NvtAll.Write();
  NvtGood.Write();
  NvtGood_BDT.Write(); 
  NvtGood_RANK.Write();

  NpuAll.Write();
  NpuGood.Write();
  NpuGood_BDT.Write(); 
  NpuGood_RANK.Write();
    
  PtGood_matchedClosest.Write();
  PtGood_BDT_matchedClosest.Write();
  PtGood_RANK_matchedClosest.Write();
  
  EtaGood_matchedClosest.Write();
  EtaGood_BDT_matchedClosest.Write();
  EtaGood_RANK_matchedClosest.Write();
  
  NvtGood_matchedClosest.Write();
  NvtGood_BDT_matchedClosest.Write(); 
  NvtGood_RANK_matchedClosest.Write();

  NpuGood_matchedClosest.Write();
  NpuGood_BDT_matchedClosest.Write(); 
  NpuGood_RANK_matchedClosest.Write();

  hdist.Write();
  hdiff_dZ_muons.Write();
  hdiff_dZ_electrons.Write();
  
  BDToutput.Write();
  BDToutput_sig.Write();
  BDToutput_bkg.Write();
   
  ff.Close();
  
  std::cout << "BYE BYE !!!! " << std::endl;

  return 0; 
  

}











