
///==== test program ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

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

#include "../src/selection.cc"
#include "../src/eleId95.cc"

#define etaEB 1.444
#define etaEE 1.560

#define PI 3.141592654


bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta);


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
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  int isData   = gConfigParser -> readIntOption("Options::isData");
  
  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  int useWeights      = gConfigParser -> readIntOption("Options::useWeights");
  int poissonWeights  = gConfigParser -> readIntOption("Options::poissonWeights");
  int nAvePU          = gConfigParser -> readIntOption("Options::nAvePU");
  
  std::string puweightsFileName = gConfigParser -> readStringOption("Options::puweightsFileName");  



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
    
    nmax = hweights ->GetMaximum();
    std::cout << nmax << std::endl;
    
    for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
      w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
    }
    weightsFile.Close();
  }
  

  std::cout << std::endl;
  std::cout << std::endl;
  

  //Filling the ntuple
  int isSig, evtNumber, nVertices, nTracks, nTracksPt05, nTracksPt4;
  float sumPt2, modSumPt, modSumPtInCone_30, modSumPtInCone_45, deltaPhi_HSumPt, sum2PhoPt, normalizedChi2, sumModPt;
  float ptbal, ptasym;
  float photons_eta[2], photons_phi[2], photons_E[2], photons_pt[2], photonsSC_eta[2], photonsSC_phi[2], photonsSC_E[2], photonsMC_eta[2], photonsMC_phi[2], photonsMC_E[2];   
  
  float trackPx[1000], trackPy[1000], trackPz[1000];

  TTree* outTree = new TTree("TMVA_vertexTree","TMVA_vertexTree");
  outTree -> Branch("evtNumber", &evtNumber, "evtNumber/I");
  outTree -> Branch("nVertices", &nVertices, "nVertices/I");

  outTree -> Branch("sumPt2",   &sumPt2,    "sumPt2/F");
  outTree -> Branch("modSumPt", &modSumPt,  "modSumPt/F");
  outTree -> Branch("sumModPt", &sumModPt,  "sumModPt/F");
  outTree -> Branch("modSumPtInCone_30", &modSumPtInCone_30,  "modSumPtInCone_30/F");
  outTree -> Branch("modSumPtInCone_45", &modSumPtInCone_45,  "modSumPtInCone_45/F");
  outTree -> Branch("deltaPhi_HSumPt", &deltaPhi_HSumPt,  "deltaPhi_HSumPt/F");
  outTree -> Branch("sum2PhoPt",  &sum2PhoPt,   "sum2PhoPt/F");
  outTree -> Branch("nTracks",&nTracks, "nTracks/I");
  outTree -> Branch("nTracksPt05",&nTracksPt05, "nTracksPt05/I");
  outTree -> Branch("nTracksPt4",&nTracksPt4, "nTracksPt4/I");
  outTree -> Branch("ptbal",  &ptbal,   "ptbal/F");
  outTree -> Branch("ptasym",  &ptasym,   "ptasym/F");
  outTree -> Branch("normalizedChi2",&normalizedChi2, "normalizedChi2/F");
  
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

  outTree -> Branch("trackPx",trackPx, "trackPx[nTracks]/F");
  outTree -> Branch("trackPy",trackPy, "trackPy[nTracks]/F");
  outTree -> Branch("trackPz",trackPz, "trackPz[nTracks]/F");

  outTree -> Branch("isSig", &isSig, "isSig/I");

  
  //Chain
  TChain* chain = new TChain(treeName.c_str());
  chain->Add((baseDir+inputFile).c_str());
  treeReader reader((TTree*)(chain));
  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;
  
  
  //start loop over entries
  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      //setup common branches
      #include "../includeMe_forTMVA_check.h"
      
      std::vector<float>*PU_z ;
      std::vector<int>* mc_PUit_NumInteractions; 
      int npu ;
      

      //--- use weights
      
      if (useWeights){
	
	//PU_z  = reader.GetFloat("PU_z");
	//npu = PU_z->size();
	
	mc_PUit_NumInteractions  = reader.GetInt("mc_PUit_NumInteractions");
	npu = mc_PUit_NumInteractions->at(0);
	
	float myrnd = gRandom->Uniform(0,nmax);
       	if (myrnd > w[npu]) continue;
      }
      
      
      //global variables
      float TrueVertex_Z;
      
      std::vector<float>* PV_normalizedChi2;
      std::vector<int>* PV_ndof;
      std::vector<int>* PV_nTracks;
      std::vector<float>* PV_z;
      std::vector<float>* PV_d0;
      std::vector<float>* PV_SumPt2;
      std::vector<float>* PV_SumPt;
      
      std::vector<ROOT::Math::XYZVector>* PVtracks;
      std::vector<int>* PVtracks_PVindex;
      std::vector<int>* PVtracks_numberOfValidHits;
      std::vector<float>* PVtracks_normalizedChi2;
      
      
      ROOT::Math::XYZTVector sum2pho;
      //---------------------------
      //--- set up Hgg branches ---
      //---------------------------
      if ( isHiggs )
	{
	  
	  //taking H variables
	  std::vector<ROOT::Math::XYZVector>* TrueVertex = reader.Get3V("mc_H_vertex");
	  std::vector<ROOT::Math::XYZTVector>* mcH = reader.Get4V("mc_H");
	  
	  std::vector<ROOT::Math::XYZTVector>* mcV1 = reader.Get4V("mcV1");
	  std::vector<ROOT::Math::XYZTVector>* mcV2 = reader.Get4V("mcV2");
	  
	  if ( mcV1->size() != 1 ||  mcV2->size() != 1) continue;
	  if (mcV1->at(0).E() > mcV2->at(0).E())
	    {
	      photonsMC_eta[0] = mcV1->at(0).Eta();
	      photonsMC_eta[1] = mcV2->at(0).Eta();
	      photonsMC_phi[0] = mcV1->at(0).Phi();
	      photonsMC_phi[1] = mcV2->at(0).Phi();
	      photonsMC_E[0]   = mcV1->at(0).E();
	      photonsMC_E[1]   = mcV2->at(0).E();
	    }
	  else
	    {
	      photonsMC_eta[1] = mcV1->at(0).Eta();
	      photonsMC_eta[0] = mcV2->at(0).Eta();
	      photonsMC_phi[1] = mcV1->at(0).Phi();
	      photonsMC_phi[0] = mcV2->at(0).Phi();
	      photonsMC_E[1]   = mcV1->at(0).E();
	      photonsMC_E[0]   = mcV2->at(0).E();
	    }
	  
	  
	  if (TrueVertex->size() != 1) continue;
	  TrueVertex_Z = TrueVertex->at(0).Z();
	  
	  PV_normalizedChi2 = reader.GetFloat("PV_normalizedChi2");
	  PV_ndof = reader.GetInt("PV_ndof");
	  PV_nTracks = reader.GetInt("PV_nTracks");
	  PV_z = reader.GetFloat("PV_z");
	  PV_d0 = reader.GetFloat("PV_d0");
	  PV_SumPt2 = reader.GetFloat("PV_SumPt2");
	  PV_SumPt = reader.GetFloat("PV_SumPt");
	  
	  PVtracks = reader.Get3V("PVtracks");
	  PVtracks_PVindex = reader.GetInt("PVtracks_PVindex");
	  PVtracks_numberOfValidHits = reader.GetInt("PVtracks_numberOfValidHits");
	  PVtracks_normalizedChi2 = reader.GetFloat("PVtracks_normalizedChi2");
	  
	  int indpho1 = -100, indpho2 = -100;
	  int ngood = 0;
	  
	  double dR_1_min = 10000.;
	  double dR_2_min = 10000.;
	  
	  for(unsigned int u=0; u < photons->size(); u++)
	    {
	      double dR_1 = deltaR( photonsMC_eta[0], photonsMC_phi[0], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
	      if (dR_1 < dR_1_min) { dR_1_min = dR_1; indpho1 = u; }
	      
	      double dR_2 = deltaR( photonsMC_eta[1], photonsMC_phi[1], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
	      if (dR_2 < dR_2_min) { dR_2_min = dR_2; indpho2 = u; }
	    }
	  
	  if(dR_1_min > 0.15 || dR_2_min > 0.15) continue;
	  if(photons_r9->at(indpho1) < 0.93 || photons_r9->at(indpho2) < 0.93) continue;
	  
	  sum2pho = photons->at(indpho1)+ photons->at(indpho2);
	  
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
	  
	}//Hgg
      //---------------------------
      //--- set up Zee branches ---
      //---------------------------
      if ( isZee )
	{
	  std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
	  std::vector<ROOT::Math::XYZTVector>* electrons_SC = reader.Get4V("electrons");
	  
	  //taking Zee variables
	  //	  std::vector<float>* eleid = reader.GetFloat("simpleEleId95cIso");
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
	  PV_normalizedChi2 = reader.GetFloat("PV_noEle_normalizedChi2");
	  PV_ndof = reader.GetInt("PV_noEle_ndof");
	  PV_nTracks = reader.GetInt("PV_noEle_nTracks");
	  PV_z = reader.GetFloat("PV_noEle_z");
	  PV_d0 = reader.GetFloat("PV_noEle_d0");
	  PV_SumPt2 = reader.GetFloat("PV_noEle_SumPt2");
	  PV_SumPt = reader.GetFloat("PV_noEle_SumPt");
	  
	  PVtracks = reader.Get3V("PVEleLessTracks");
	  PVtracks_PVindex = reader.GetInt("PVEleLessTracks_PVindex");
	  PVtracks_numberOfValidHits = reader.GetInt("PVEleLessTracks_numberOfValidHits");
	  PVtracks_normalizedChi2 = reader.GetFloat("PVEleLessTracks_normalizedChi2");
	  
	  int indele1 = -100, indele2 = -100;
	  int ngood = 0;
	  for(unsigned int uu = 0; uu < eleid->size(); uu++)
	    {
	      if (  eleid->at(uu) == 7 && indele1 < 0 )      {indele1 = uu; ngood++;}
	      else if ( eleid->at(uu) == 7 && indele2 < 0 )  {indele2 = uu; ngood++;}
	      else if( eleid->at(uu) == 7 )                  { ngood++;}
	    }
	  if ( ngood != 2){continue;}
	  
	  sum2pho = electrons->at(indele1) + electrons->at(indele2);
	  
	  if (electrons->at(indele1).E() > electrons->at(indele2).E())
	    {
	      photons_eta[0] = electrons->at(indele1).eta();
	      photons_eta[1] = electrons->at(indele2).eta();
	      photons_phi[0] = electrons->at(indele1).phi();
	      photons_phi[1] = electrons->at(indele2).phi();
	      photons_E[0]   = electrons->at(indele1).E();
	      photons_E[1]   = electrons->at(indele2).E();
	      photons_pt[0]  = electrons->at(indele1).pt();
	      photons_pt[1]  = electrons->at(indele2).pt();
	      
	      photonsSC_eta[0] = electrons_SC->at(indele1).eta();
	      photonsSC_eta[1] = electrons_SC->at(indele2).eta();
	      photonsSC_phi[0] = electrons_SC->at(indele1).phi();
	      photonsSC_phi[1] = electrons_SC->at(indele2).phi();
	      photonsSC_E[0]   = electrons_SC->at(indele1).E();
	      photonsSC_E[1]   = electrons_SC->at(indele2).E();
	    }
	  else
	    {
	      photons_eta[1] = electrons->at(indele1).eta();
	      photons_eta[0] = electrons->at(indele2).eta();
	      photons_phi[1] = electrons->at(indele1).phi();
	      photons_phi[0] = electrons->at(indele2).phi();
	      photons_E[1] = electrons->at(indele1).E();
	      photons_E[0] = electrons->at(indele2).E();
	      photons_pt[1] = electrons->at(indele1).pt();
	      photons_pt[0] = electrons->at(indele2).pt();
	      
	      photonsSC_eta[1] = electrons_SC->at(indele1).eta();
	      photonsSC_eta[0] = electrons_SC->at(indele2).eta();
	      photonsSC_phi[1] = electrons_SC->at(indele1).phi();
	      photonsSC_phi[0] = electrons_SC->at(indele2).phi();
	      photonsSC_E[1] = electrons_SC->at(indele1).E();
	      photonsSC_E[0] = electrons_SC->at(indele2).E();
	    }
	  
	  
	  if ( fabs(sum2pho.M() - 91) > 8 ) continue;
	  
	  TrueVertex_Z = PV_z->at(0) + (electrons_dz_PV_noEle->at(indele1) + electrons_dz_PV_noEle->at(indele2))/2.;
	  
	}//Zee end
      //---------------------------
      //--- set up Zee branches ---
      //---------------------------
      if ( isZmumu )
	{
	  
	  //taking Zmumu variables
	  PV_normalizedChi2 = reader.GetFloat("PV_noMuon_normalizedChi2");
	  PV_ndof = reader.GetInt("PV_noMuon_ndof");
	  PV_nTracks = reader.GetInt("PV_noMuon_nTracks");
	  PV_z = reader.GetFloat("PV_noMuon_z");
	  PV_d0 = reader.GetFloat("PV_noMuon_d0");
	  PV_SumPt2 = reader.GetFloat("PV_noMuon_SumPt2");
	  PV_SumPt = reader.GetFloat("PV_noMuon_SumPt");
	  
	  PVtracks = reader.Get3V("PVMuonLessTracks");
	  PVtracks_PVindex = reader.GetInt("PVMuonLessTracks_PVindex");
	  PVtracks_numberOfValidHits = reader.GetInt("PVMuonLessTracks_numberOfValidHits");
	  PVtracks_normalizedChi2 = reader.GetFloat("PVMuonLessTracks_normalizedChi2");
	  
	  int indmu1 = -100, indmu2 = -100;
	  int ngood = 0;
	  for( int uu = 0; uu < muons->size(); uu++)
	    {
	      if ( muons->at(uu).pt() > 10 && indmu1 < 0 )        { indmu1 = uu; ngood++; }
	      else if ( muons->at(uu).pt() > 10 && indmu2 < 0 )   { indmu2 = uu; ngood++; }
	      else if ( muons->at(uu).pt() > 10 )                 { ngood++;}
	    }
	  
	  if ( ngood != 2){continue;}
	  
	  sum2pho = muons->at(indmu1) + muons->at(indmu2);
	  
	  if (muons->at(indmu1).E() > muons->at(indmu2).E())
	    {
	      photons_eta[0] = muons->at(indmu1).eta();
	      photons_eta[1] = muons->at(indmu2).eta();
	      photons_phi[0] = muons->at(indmu1).phi();
	      photons_phi[1] = muons->at(indmu2).phi();
	      photons_E[0]   = muons->at(indmu1).E();
	      photons_E[1]   = muons->at(indmu2).E();
	      photons_pt[0]  = muons->at(indmu1).pt();
	      photons_pt[1]  = muons->at(indmu2).pt();
	      
	      photonsSC_eta[0] = muons->at(indmu1).eta();
	      photonsSC_eta[1] = muons->at(indmu2).eta();
	      photonsSC_phi[0] = muons->at(indmu1).phi();
	      photonsSC_phi[1] = muons->at(indmu2).phi();
	      photonsSC_E[0]   = muons->at(indmu1).E();
	      photonsSC_E[1]   = muons->at(indmu2).E();
	    }
	  else
	    {
	      photons_eta[1] = muons->at(indmu1).eta();
	      photons_eta[0] = muons->at(indmu2).eta();
	      photons_phi[1] = muons->at(indmu1).phi();
	      photons_phi[0] = muons->at(indmu2).phi();
	      photons_E[1]   = muons->at(indmu1).E();
	      photons_E[0]   = muons->at(indmu2).E();
	      photons_pt[1]  = muons->at(indmu1).pt();
	      photons_pt[0]  = muons->at(indmu2).pt();
	      
	      photonsSC_eta[1] = muons->at(indmu1).eta();
	      photonsSC_eta[0] = muons->at(indmu2).eta();
	      photonsSC_phi[1] = muons->at(indmu1).phi();
	      photonsSC_phi[0] = muons->at(indmu2).phi();
	      photonsSC_E[1]   = muons->at(indmu1).E();
	      photonsSC_E[0]   = muons->at(indmu2).E();
	    }
	  
	  
	  if ( fabs(sum2pho.M() - 91) > 8 ) continue;
	  TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indmu1) + muons_dz_PV_noMuon->at(indmu2))/2.;
	  
	}//Zmumu end
      
      
      
      
      //vettore degli indici dei PV buoni (no splitting)
      std::vector<int> goodIndex;
      for (unsigned int vItr = 0; vItr < PV_z->size(); ++vItr)
	{
	  //if (PV_nTracks->at(vItr) > 3 || PV_SumPt2->at(vItr) > 10.) goodIndex.push_back(vItr);
	  goodIndex.push_back(vItr);
	}
      
      int npv = goodIndex.size();
      if (npv == 0 ){continue;}
      
      //look if the H vertex matches one of the PV vertices
      float dmin = 10000;
      int iClosest = -1;
      int Vmatched = 0;
      for ( int uu = 0; uu < npv; uu++)
	{
	  float distt = fabs( PV_z->at(goodIndex[uu]) - TrueVertex_Z );
	  if ( distt < 0.3 )   { Vmatched++; }
	  if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	}
      
      //pt balance for vertexId
      for ( int uu = 0; uu < npv; uu++)
	{
	  nTracks     = 0;
	  nTracksPt05 = 0;
	  nTracksPt4  = 0;
	  
	  ROOT::Math::XYZVector sumTracks;
	  ROOT::Math::XYZVector sumTracksInCone_30;
	  ROOT::Math::XYZVector sumTracksInCone_45;
	  float vtxPx = 0;
	  float vtxPy = 0;
	  float vtxPz = 0;
	  float vtxPt = 0;
	  
	  sumModPt = 0;
	  sumPt2 = 0;
	  
	  for (unsigned int kk = 0; kk < PVtracks->size(); ++kk)
	    {
	      if (PVtracks_PVindex->at(kk) == goodIndex[uu])
		{
		  if ( PVtracks->at(kk).perp2() < trackThr*trackThr ) continue;
		  
// 		  double dr1 = deltaR( photonsSC_eta[0], photonsSC_phi[0], PVtracks->at(kk).eta(), PVtracks->at(kk).phi());
// 		  double dr2 = deltaR( photonsSC_eta[1], photonsSC_phi[1], PVtracks->at(kk).eta(), PVtracks->at(kk).phi());
// 		  if(dr1 < 0.05 || dr2 < 0.05) continue;
		  
		  vtxPx += PVtracks->at(kk).X();
		  vtxPy += PVtracks->at(kk).Y();
		  vtxPz += PVtracks->at(kk).Z();
		  
		  sumTracks += PVtracks->at(kk);
		  sumModPt += sqrt((PVtracks->at(kk)).perp2());
		  sumPt2 += PVtracks->at(kk).perp2();
		  if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*11./12. ) sumTracksInCone_30 += PVtracks->at(kk);
		  if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*7./8. ) sumTracksInCone_45 += PVtracks->at(kk);
		  
		  trackPx[nTracks] = PVtracks->at(kk).X();
		  trackPy[nTracks] = PVtracks->at(kk).Y();
		  trackPz[nTracks] = PVtracks->at(kk).Z();
		  
		  nTracks ++;
		  if ( PVtracks->at(kk).perp2() > 0.25 ) nTracksPt05++;
		  if ( PVtracks->at(kk).perp2() > 16. )  nTracksPt4++;
		  
		}
	    }//tracks loop
	  
	  deltaPhi_HSumPt = deltaPhi( sumTracks.phi(),sum2pho.phi());
	  modSumPt = sqrt( sumTracks.perp2() );
	  modSumPtInCone_30 = sqrt( sumTracksInCone_30.perp2() );
	  modSumPtInCone_45 = sqrt( sumTracksInCone_45.perp2() );
	  sum2PhoPt = sum2pho.pt();	   
	  normalizedChi2 = PV_normalizedChi2->at(goodIndex[uu]);
	  
	  ptbal  = -(vtxPx*sum2pho.X() + vtxPy*sum2pho.Y())/sum2pho.pt(); 
	  vtxPt  = sqrt(vtxPx*vtxPx+vtxPy*vtxPy);
	  ptasym = (vtxPt - sum2PhoPt)/ (vtxPt + sum2PhoPt);
	  
	  
	  isSig = 0.;
	  if (uu == iClosest)
	    isSig = 1;
	  
	  
	  
	  evtNumber = u;
	  nVertices = npv;
	  
	  outTree -> Fill();	   
	}//vertex loop
    } //evt loop
  
  
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
  
  outTree -> Write();   
  
  ff.Close();
  return 0;
  
  
}


bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta ) 
{
  
  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );
  
  return (iso && shape && HoE<0.05);

}
