
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
  std::string treeName      = gConfigParser -> readStringOption("Input::treeName");

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

  // Get run/LS map from JSON file
  std::cout << ">>> Getting GOOD  run/LS from JSON file" << std::endl;
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(jsonFileName);


  std::cout << std::endl;
  std::cout << std::endl;
  
  //------ use weights for MC
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

    for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
      // trick to skip low lumi runs that have very large weights
      if ( ibin < 6 ) hweights->SetBinContent(ibin,0.);
      w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
    }
    nmax = hweights ->GetMaximum();
    std::cout << " Max weight " << nmax << std::endl;
    weightsFile.Close();
  }
  
  // Chain
  TChain* chain = new TChain(treeName.c_str());
  FillChain(*chain, inputFileList.c_str());
  treeReader reader((TTree*)(chain));
  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  // prepare stuff for vertex analysis

  // algo parameters
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 0.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;
  vtxAlgoParams_.removeTracksInCone = 1;
  vtxAlgoParams_.coneSize           = 0.05;
  
  //variables for RANKING AND mva based on sumpt2, ptbal, ptasym to rank vertices 
  vector<string> rankVariables_;
  rankVariables_.push_back("ptbal"), rankVariables_.push_back("ptasym"),rankVariables_.push_back("logsumpt2");  
  
  vector<string> tmvaPerVtxVariables_;
  tmvaPerVtxVariables_.push_back("ptbal"), tmvaPerVtxVariables_.push_back("ptasym"),tmvaPerVtxVariables_.push_back("logsumpt2");  
  if( addConversionToMva ) {
    tmvaPerVtxVariables_.push_back("limPullToConv");
    tmvaPerVtxVariables_.push_back("nConv");
  } 

  TMVA::Reader *tmvaReader_ = new TMVA::Reader( "!Color:!Silent" );
  HggVertexAnalyzer::bookVariables( *tmvaReader_, tmvaPerVtxVariables_ );
  tmvaReader_->BookMVA( tmvaMethod.c_str(), tmvaWeights.c_str() );



  //histograms 

  TH1F *tmvaPerVertexOutput[5];
  char hname[100];
  for (int i = 0; i < 5; i++){
    sprintf(hname,"tmvaPerVertexOutput_%d",i);
    tmvaPerVertexOutput[i] = new TH1F(hname,hname,200,-1,1);
  }
  
  TH1F *tmvaPerVertexOutput_RightVtx = new TH1F("tmvaPerVertexOutput_RightVtx" ,"tmvaPerVertexOutput_RightVtx" ,200,-1,1);
  TH1F *tmvaPerVertexOutput_WrongVtx = new TH1F("tmvaPerVertexOutput_WrongVtx" ,"tmvaPerVertexOutput_WrongVtx" ,200,-1,1);

  TH2F *ptbal_vs_ptasym = new TH2F("ptbal_vs_ptasym","ptbal_vs_ptasym",400, -1, 1 , 400,-200,200);
  TH2F *ptbal_vs_sumpt2 = new TH2F("ptbal_vs_sumpt2","ptbal_vs_sumpt2",300, -30, 30 , 400,-200,200);
  TH2F *ptasym_vs_sumpt2 = new TH2F("ptasym_vs_sumpt2","ptasym_vs_sumpt2",300, -30, 30 , 400,-1,1);


  TH1F * hdmin = new TH1F("hdmin","hdmin",500, -2,2);

  // output tree
  int evtNumber, nVertices, npu;
  float trueVtx_z;
  int isSig[5];
  float ptbal[5], ptasym[5], logsumpt2[5], diphopt[5], pv_z[5], mva[5];
  float  photonSC_eta[2], photonSC_phi[2], photonSC_E[2];   

  int nConv;
  float pullToConv[5];

  TTree* outTree = new TTree("TMVAtree","TMVAtree");
  outTree -> Branch("evtNumber", &evtNumber, "evtNumber/I");
  outTree -> Branch("nVertices", &nVertices, "nVertices/I");
  outTree -> Branch("npu", &npu, "npu/I");
  outTree -> Branch("photonSC_eta",         photonSC_eta,            "photonSC_eta[2]/F");
  outTree -> Branch("photonSC_phi",         photonSC_phi,            "photonSC_phi[2]/F");
  outTree -> Branch("photonSC_E",           photonSC_E,              "photonSC_E[2]/F");
  outTree -> Branch("logsumpt2",   logsumpt2,    "logsumpt2[5]/F");
  outTree -> Branch("ptbal",    ptbal,     "ptbal[5]/F");
  outTree -> Branch("ptasym",   ptasym,    "ptasym[5]/F");
  outTree -> Branch("diphopt",  diphopt,   "diphopt[5]/F");
  outTree -> Branch("mva",  mva,   "mva[5]/F");
  outTree -> Branch("pv_z", pv_z, "pv_z[5]/F");
  outTree -> Branch("isSig", isSig, "isSig[5]/I");
  outTree -> Branch("trueVtx_z", &trueVtx_z, "trueVtx_z/F");
  outTree -> Branch("nConv", &nConv, "nConv/I");
  outTree -> Branch("pullToConv", pullToConv,"pullToConv[5]/F");


  float ww = 1;
  float r9cut = 0.93;
 
  //start loop over entries
  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      // filter bad runs/lumis
      int runId = reader.GetInt("runId")->at(0);
      int lumiId = reader.GetInt("lumiId")->at(0);

      bool skipEvent = false;
      if( isData && useJSON ){
	if(AcceptEventByRunAndLumiSection(runId,lumiId,jsonMap) == false)
	  skipEvent = true;
      }

      if( skipEvent == true ) continue;

      //*** pu weights
      std::vector<float>*PU_z ;
      std::vector<int>* mc_PUit_NumInteractions; 
      std::vector<float>* mc_PUit_TrueNumInteractions; 
      int npu ;
      float npuTrue;
      
      if ( !isData ){
	mc_PUit_NumInteractions  = reader.GetInt("mc_PUit_NumInteractions");
	npu = mc_PUit_NumInteractions->at(0);

      	mc_PUit_TrueNumInteractions  = reader.GetFloat("mc_PUit_TrueNumInteractions"); // needed for 2012 PU reweighting
	npuTrue = mc_PUit_TrueNumInteractions->at(0);

	//--- use weights 
	if (useWeights){
	  float myrnd = gRandom->Uniform(0,nmax);
	  //if (myrnd > w[npu]) continue; // used in 2011
	  if (myrnd > w[int(npuTrue)]) continue; // for 2012
	  //ww = w[int(npuTrue)];
	}
      }

      //setup common branches
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

      // selections for Hgg
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


      //selections for Zee
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


      }

      // selections for Zmumu
      if ( isZmumu )
	 {
	   std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
	   std::vector<int>* muons_global = reader.GetInt("muons_global");
	   std::vector<int>* muons_tracker = reader.GetInt("muons_tracker");
	   std::vector<float>* muons_tkIsoR03 = reader.GetFloat("muons_tkIsoR03");
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
	   
	   zmumuSelection(muons,muons_global,muons_tracker, muons_tkIsoR03,accept, indpho1, indpho2);

	   if (!accept) continue;
	   
	   etaMaxSC = muons->at(indpho1).eta();
	   sum2pho  = muons->at(indpho1)+ muons->at(indpho2);
	   TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indpho1) + muons_dz_PV_noMuon->at(indpho2))/2.;
	   
	   //if ( fabs(muons_dz_PV_noMuon->at(indpho1) - muons_dz_PV_noMuon->at(indpho2))> 0.5) continue;

   
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
      
      //photons
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
   
      
      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
      

     
      PhotonInfo pho1(indpho1, TVector3(phocalox_[indpho1],phocaloy_[indpho1],phocaloz_[indpho1]),phoen_[indpho1]); 
      PhotonInfo pho2(indpho2, TVector3(phocalox_[indpho2],phocaloy_[indpho2],phocaloz_[indpho2]),phoen_[indpho2]); 
     
      //look if the H vertex matches one of the PV vertices
      float dmin = 10000;
      int iClosest = -1;
      for ( int uu = 0; uu < nvtx_; uu++)
	{
	  float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	  if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	}
      hdmin->Fill(dmin);

      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);
      vAna.analyze(vinfo,pho1,pho2);
      
      // if Zmumu : setNconv to zero
      vAna.setNConv(0);


      // order vertices using rankprod algo
      vector<int> rankprod = vAna.rankprod(rankVariables_);
      // order vertices using mva 3 variables 
      if ( useMvaRanking )
	rankprod = vAna.rank(*tmvaReader_,tmvaMethod);


      // preselect first 3 vertices
      vAna.preselection(rankprod);
      vAna.evaluate(*tmvaReader_,tmvaMethod);

      for (int iv=0;iv<rankprod.size();iv++) {
	float vtxmva = vAna.mva(rankprod[iv]);

	ptbal_vs_ptasym->Fill(vAna.ptasym(rankprod[iv]),vAna.ptbal(rankprod[iv]));
	ptbal_vs_sumpt2->Fill(vAna.logsumpt2(rankprod[iv]), vAna.ptbal(rankprod[iv]));
	ptasym_vs_sumpt2->Fill(vAna.logsumpt2(rankprod[iv]),vAna.ptasym(rankprod[iv]));

	if (rankprod[iv] == iClosest) {
	  tmvaPerVertexOutput_RightVtx->Fill(vtxmva);
	}
	else  
	  tmvaPerVertexOutput_WrongVtx->Fill(vtxmva);
     
 	// only first 5 ranked vertices
	if ( iv < 5 ){
	  tmvaPerVertexOutput[iv]->Fill(vtxmva);
 	  isSig[iv] = 0;  
 	  //if ( fabs( TrueVertex_Z - PV_z->at( rankprod[iv] ) ) < 1.) // matching within 1 cm
	  if (rankprod[iv] == iClosest)  // closest vtx
 	    isSig[iv] = 1;
 	  logsumpt2[iv] = vAna.logsumpt2(rankprod[iv]) ;
 	  ptbal [iv]    = vAna.ptbal(rankprod[iv])  ;
	  ptasym[iv]    = vAna.ptasym(rankprod[iv]) ;
	  diphopt[iv]   = vAna.diphopt(rankprod[iv]);
	  pullToConv[iv]= vAna.limpulltoconv(rankprod[iv]);
	  pv_z[iv]      = PV_z->at(rankprod[iv])    ;
	  mva[iv]       = vtxmva;
 	}

      } // end loop over vertices

      
    
      // filling output mini-tree
	  
	  
      evtNumber = u;
      nVertices = nvtx_;
      trueVtx_z = TrueVertex_Z;
      nConv     = vAna.nconv(vAna.pairID(indpho1,indpho2));
      
      photonSC_eta[0] = pho1.caloPosition().Eta(); 
      photonSC_phi[0] = pho1.caloPosition().Phi(); 
      photonSC_E[0]   = pho1.energy();

      photonSC_eta[1] = pho2.caloPosition().Eta(); 
      photonSC_phi[1] = pho2.caloPosition().Phi(); 
      photonSC_E[1]   = pho2.energy();
      
      outTree -> Fill();	  


    } // end loop over entries
  
   
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
  
  for (int i = 0; i < 5; i++){
    tmvaPerVertexOutput[i]->Write();
  }

  tmvaPerVertexOutput_RightVtx->Write();
  tmvaPerVertexOutput_WrongVtx->Write();

  ptbal_vs_ptasym->Write();
  ptbal_vs_sumpt2->Write();
  ptasym_vs_sumpt2->Write();

  hdmin->Write();

  outTree->Write();

  ff.Close();
  return 0; 
  
}











