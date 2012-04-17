
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
    
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  

  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
 
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  int isData   = gConfigParser -> readIntOption("Options::isData");
  
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
      sprintf(hname,"hwpoisson_%d",nAvePU);
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
  
  //****** open TH1 with tracks pt for kolmogorov test
  TH1F *hTrackPt;
  TFile ptFile("hTracksPt_PU25_upTo10GeV.root","READ"); 
  hTrackPt = (TH1F*)ptFile.Get("hTracksPt")->Clone("hTrackPt");
  hTrackPt->SetDirectory(0);
  ptFile.Close();
  std::cout << hTrackPt->GetEntries() << std::endl;

  TH1F *htemp = (TH1F*)hTrackPt->Clone("htemp");
  htemp->Reset();
 
    
  //****** Parameters for vertex finding algo ****** 
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 1.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;
  vtxAlgoParams_.removeTracksInCone = 1;
  vtxAlgoParams_.coneSize           = 0.05;
  
      
  //****** BOOK TREE ******
  int isSig, evtNumber, nVertices;
  float diphopt; 
  float logsumpt2;
  float ptbal, ptasym, ptmax, ptmax3, sumpt, nchthr, nch, sumtwd, pzasym;
  float photons_eta[2], photons_phi[2], photons_E[2], photons_pt[2], photonsSC_eta[2], photonsSC_phi[2], photonsSC_E[2], photonsMC_eta[2], photonsMC_phi[2], photonsMC_E[2];   
  
  // no ptratio : 100% con ptbal

  float trackPt[1000];
  float dk;

  TTree* outTree = new TTree("TMVA_vertexTree","TMVA_vertexTree");
  outTree -> Branch("evtNumber", &evtNumber, "evtNumber/I");
  outTree -> Branch("nVertices", &nVertices, "nVertices/I");

  outTree -> Branch("diphopt",   &diphopt,   "diphopt/F");
  outTree -> Branch("logsumpt2", &logsumpt2, "logsumpt2/F");
  outTree -> Branch("nch",       &nch,       "nch/F");
  outTree -> Branch("nchthr",    &nchthr,    "nchthr/F");
  outTree -> Branch("ptmax",     &ptmax,     "ptmax/F");
  outTree -> Branch("ptmax3",    &ptmax3,    "ptmax3/F");
  outTree -> Branch("ptbal",     &ptbal,     "ptbal/F");
  outTree -> Branch("ptasym",    &ptasym,    "ptasym/F");
  outTree -> Branch("sumtwd",    &sumtwd,    "sumtwd/F");
  outTree -> Branch("dk",        &dk,        "dk/F");

  outTree -> Branch("photons_eta",           photons_eta,              "photons_eta[2]/F");
  outTree -> Branch("photons_phi",           photons_phi,              "photons_phi[2]/F");
  outTree -> Branch("photons_E",             photons_E,                "photons_E[2]/F");
  outTree -> Branch("photons_pt",            photons_pt,               "photons_pt[2]/F");
  outTree -> Branch("photonsSC_eta",         photonsSC_eta,            "photonsSC_eta[2]/F");
  outTree -> Branch("photonsSC_phi",         photonsSC_phi,            "photonsSC_phi[2]/F");
  outTree -> Branch("photonsSC_E",           photonsSC_E,              "photonsSC_E[2]/F");
  outTree -> Branch("photonsMC_eta",         photonsMC_eta,            "photonsMC_eta[2]/F");
  outTree -> Branch("photonsMC_phi",         photonsMC_phi,            "photonsMC_phi[2]/F");
  outTree -> Branch("photonsMC_E",           photonsMC_E,              "photonsMC_E[2]/F");

  outTree -> Branch("isSig", &isSig, "isSig/I");
    

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
      std::vector<ROOT::Math::XYZTVector>* photons_SC; // supercluster
      
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
	photons_SC       = reader.Get4V("photons_SC");

	hggSelection(mcV1, mcV2, photons, photons_SC, photons_r9, accept, indpho1, indpho2);
	if (!accept) continue;

	if ( photons_r9->at(indpho1) < r9cut ) continue;
	if ( photons_r9->at(indpho2) < r9cut ) continue;
	
	etaMaxSC = photons_SC->at(indpho1).eta();
	sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	TrueVertex_Z = mc_H_vertex->at(0).Z();
	
	if (photons->at(indpho1).E() > photons->at(indpho2).E())
	  {
	    photons_eta[0] = photons->at(indpho1).eta();
	    photons_eta[1] = photons->at(indpho2).eta();
	    photons_phi[0] = photons->at(indpho1).phi();
	    photons_phi[1] = photons->at(indpho2).phi();
	    photons_E[0]   = photons->at(indpho1).E();
	    photons_E[1]   = photons->at(indpho2).E();
	    photons_pt[0]  = photons->at(indpho1).pt();
	    photons_pt[1]  = photons->at(indpho2).pt();
	    
	    photonsSC_eta[0] = photons_SC->at(indpho1).eta();
	    photonsSC_eta[1] = photons_SC->at(indpho2).eta();
	    photonsSC_phi[0] = photons_SC->at(indpho1).phi();
	    photonsSC_phi[1] = photons_SC->at(indpho2).phi();
	    photonsSC_E[0] = photons_SC->at(indpho1).E();
	    photonsSC_E[1] = photons_SC->at(indpho2).E();
	  }
	else
	  {
	    photons_eta[1] = photons->at(indpho1).eta();
	    photons_eta[0] = photons->at(indpho2).eta();
	    photons_phi[1] = photons->at(indpho1).phi();
	    photons_phi[0] = photons->at(indpho2).phi();
	    photons_E[1]   = photons->at(indpho1).E();
	    photons_E[0]   = photons->at(indpho2).E();
	    photons_pt[1]  = photons->at(indpho1).pt();
	    photons_pt[0]  = photons->at(indpho2).pt();
	    
	    photonsSC_eta[1] = photons_SC->at(indpho1).eta();
	    photonsSC_eta[0] = photons_SC->at(indpho2).eta();
	    photonsSC_phi[1] = photons_SC->at(indpho1).phi();
	    photonsSC_phi[0] = photons_SC->at(indpho2).phi();
	    photonsSC_E[1] = photons_SC->at(indpho1).E();
	    photonsSC_E[0] = photons_SC->at(indpho2).E();
	  }
	
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
	photons_SC       = reader.Get4V("electrons");
	
	zeeSelection(electrons, eleid, accept, indpho1, indpho2);
	if (!accept) continue;

	etaMaxSC = photons_SC->at(indpho1).eta();
	sum2pho  = electrons->at(indpho1)+ electrons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (electrons_dz_PV_noEle->at(indpho1) + electrons_dz_PV_noEle->at(indpho2))/2.;
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
	photons_SC       = reader.Get4V("muons");  // use muon info for SC
	
	zmumuSelection(muons,muons_global,muons_tracker, accept, indpho1, indpho2);
	
	if (!accept) continue;
	
	etaMaxSC = muons->at(indpho1).eta();
	sum2pho  = muons->at(indpho1)+ muons->at(indpho2);
	TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indpho1) + muons_dz_PV_noMuon->at(indpho2))/2.;
	
	if (muons->at(indpho1).E() > muons->at(indpho2).E())
	  {
	    photons_eta[0] = muons->at(indpho1).eta();
	    photons_eta[1] = muons->at(indpho2).eta();
	    photons_phi[0] = muons->at(indpho1).phi();
	    photons_phi[1] = muons->at(indpho2).phi();
	    photons_E[0]   = muons->at(indpho1).E();
	    photons_E[1]   = muons->at(indpho2).E();
	    photons_pt[0]  = muons->at(indpho1).pt();
	    photons_pt[1]  = muons->at(indpho2).pt();
	    
	    photonsSC_eta[0] = muons->at(indpho1).eta();
	    photonsSC_eta[1] = muons->at(indpho2).eta();
	    photonsSC_phi[0] = muons->at(indpho1).phi();
	    photonsSC_phi[1] = muons->at(indpho2).phi();
	    photonsSC_E[0]   = muons->at(indpho1).E();
	    photonsSC_E[1]   = muons->at(indpho2).E();
	  }
	else
	  {
	    photons_eta[1] = muons->at(indpho1).eta();
	    photons_eta[0] = muons->at(indpho2).eta();
	    photons_phi[1] = muons->at(indpho1).phi();
	    photons_phi[0] = muons->at(indpho2).phi();
	    photons_E[1]   = muons->at(indpho1).E();
	    photons_E[0]   = muons->at(indpho2).E();
	    photons_pt[1]  = muons->at(indpho1).pt();
	    photons_pt[0]  = muons->at(indpho2).pt();
	    
	    photonsSC_eta[1] = muons->at(indpho1).eta();
	    photonsSC_eta[0] = muons->at(indpho2).eta();
	    photonsSC_phi[1] = muons->at(indpho1).phi();
	    photonsSC_phi[0] = muons->at(indpho2).phi();
	    photonsSC_E[1]   = muons->at(indpho1).E();
	    photonsSC_E[0]   = muons->at(indpho2).E();
	  }
	
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
      for (int ipho = 0 ; ipho < photons_SC->size(); ipho++){
	float px = photons_SC->at(ipho).X();
	float py = photons_SC->at(ipho).Y();
	float pz = photons_SC->at(ipho).Z();
	float pt = sqrt ( px*px+py*py ); 
	float theta = 2*atan(exp(-photons_SC->at(ipho).eta()) );
	float tantheta = tan(theta);

	if ( fabs(photons_SC->at(ipho).eta()) < etaEB ) {
	  phocalox_[ipho] = R_ECAL*px/pt;
	  phocaloy_[ipho] = R_ECAL*py/pt;
	  phocaloz_[ipho] = R_ECAL/tantheta;
	} 

	if ( fabs(photons_SC->at(ipho).eta()) > etaEE ) {
	  float r_endcap  = fabs(Z_ENDCAP * tantheta);
	  phocalox_[ipho] = r_endcap * px/pt;
	  phocaloy_[ipho] = r_endcap * py/pt;
	  if (pz > 0) phocaloz_[ipho] = Z_ENDCAP;
	  else phocaloz_[ipho] = -Z_ENDCAP;
	} 
	
	phoen_[ipho] = photons_SC->at(ipho).E();
	
      }


      float eta1 = photons_SC->at(indpho1).eta();
      float eta2 = photons_SC->at(indpho2).eta();

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
      
      //*** NOW FILL THE TREE
      	  
      for ( int uu = 0; uu < nvtx_; uu++){
			
	diphopt   = vAna.diphopt(uu);
	logsumpt2 = vAna.logsumpt2(uu);
	nch       = vAna.nch(uu);
	nchthr    = vAna.nchthr(uu);
	ptbal     = vAna.ptbal(uu);
	ptasym    = vAna.ptasym(uu);
	ptmax     = vAna.ptmax(uu);
	ptmax3    = vAna.ptmax3(uu);
	sumtwd    = vAna.sumtwd(uu);

	htemp->Reset();
	for (unsigned int kk = 0; kk < PVtracks->size(); ++kk){
	  if (PVtracks_PVindex->at(kk) == uu){
	    float tkpt = sqrt(PVtracks->at(kk).perp2()) ;
	    htemp ->Fill( (tkpt < 10. ? tkpt : 10.) );
	  }
	}

	if ( htemp->GetEntries()!=0 ) 
	  dk = htemp->KolmogorovTest(hTrackPt);

	
	isSig = 0;
	if ( fabs( PV_z->at(uu) - TrueVertex_Z ) < 1. )
	  isSig = 1;
	
	evtNumber = u;
	nVertices = nvtx_;
      
	
	outTree -> Fill();
	

      }// end vertex loop
      

      
      
      
    }// end loop over entries

  std::cout << "END LOOP OVER ENTRIES" << std::endl;
  std::cout << "Saving histos on file ..." << std::endl;
  
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");

  hAcceptedLumis -> Write();

  outTree -> Write();   
  htemp->Write();

  ff.Close();
  
  std::cout << "BYE BYE !!!! " << std::endl;

  return 0; 
  

}











