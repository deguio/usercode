
///==== include ====

#include "treeReader.h"
#include "hFactory.h"
#include "hFunctions.h"
#include "stdHisto.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
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

#include <iostream>

#include "TClonesArray.h"
#include "TMatrix.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif



#define R_ECAL    129
#define Z_ENDCAP  317

#define etaEB   1.442 //?
#define etaEE   1.560


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
  
  std::string baseDir = gConfigParser -> readStringOption("Input::baseDir");
  std::string inputFile = gConfigParser -> readStringOption("Input::inputFile");
  std::string treeName = gConfigParser -> readStringOption("Input::treeName");
  std::string weights = gConfigParser -> readStringOption("Input::weights");
  
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  
  
  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  int useWeights  = gConfigParser -> readIntOption("Options::useWeights");

  std::cout << std::endl;
  std::cout << std::endl;
  
  //------ use weights for MC
  float nmax;
  float w[50];
  TRandom *gRandom = new TRandom();
  
  if (useWeights){
  
    //TFile weightsFile("PUweightsFromReco.root","READ");
    //TH1F* hweights = (TH1F*)weightsFile.Get("hWeightsReco");

    //TFile weightsFile("testPUweights.root","READ");
    //TH1F* hweights = (TH1F*)weightsFile.Get("hwdata");
    
    TFile weightsFile("PUweightsFromUnfolding.root","READ");
    TH1F* hweights = (TH1F*)weightsFile.Get("hWeightsUnfolded");

    nmax = hweights ->GetMaximum();
    std::cout << nmax << std::endl;
    
    for (int ibin = 1; ibin < hweights->GetNbinsX()+1; ibin++){
      w[ibin-1] = hweights->GetBinContent(ibin);  // bin 1 --> nvtx = 0 
    }
    weightsFile.Close();
  }
  
  
  // Chain
  TChain* chain = new TChain(treeName.c_str());
  chain->Add((baseDir+inputFile).c_str());
  treeReader reader((TTree*)(chain));


  //histos  - NB 3 categories

  TH1F *NvtAll[3];
  TH1F *PtAll[3];
  TH1F *EtaAll[3];
  TH1F *R9All[3];

  TH1F *NvtGood[3];
  TH1F *PtGood[3];
  TH1F *EtaGood[3];
  TH1F *R9Good[3];
  
  char hname[100];

  for (int i=0; i< 3; i++){
    sprintf(hname,"NvtAll_cat%d",i+1);
    NvtAll[i] = new TH1F(hname, hname, 50,0,50);

    sprintf(hname,"NvtGood_cat%d",i+1);
    NvtGood[i] = new TH1F(hname, hname, 50,0,50);

    sprintf(hname,"PtAll_cat%d",i+1);
    PtAll[i] = new TH1F(hname, hname, 80,0,400);

    sprintf(hname,"PtGood_cat%d",i+1);
    PtGood[i] = new TH1F(hname, hname, 80,0,400);

    sprintf(hname,"EtaAll_cat%d",i+1);
    EtaAll[i] = new TH1F(hname, hname, 50,-5,5);

    sprintf(hname,"EtaGood_cat%d",i+1);
    EtaGood[i] = new TH1F(hname, hname, 50,-5,5);

    sprintf(hname,"R9All_cat%d",i+1);
    R9All[i] = new TH1F(hname, hname, 100,0,1);

    sprintf(hname,"R9Good_cat%d",i+1);
    R9Good[i] = new TH1F(hname, hname, 100,0,1);
    
  }


  TH1F *hrankConv = new TH1F("hrankConv","hrankConv",50,0,50);

  TH1F* hdz1 = new TH1F("hdz1","#sigma_{z} from conversion",200,-10,10);
  TH1F* hdz2 = new TH1F("hdz2","#sigma_{z} from conversion",200,-10,10);

  TH1F* hdeltaz1FromTrueVertex = new TH1F("hdeltaz1FromTrueVertex","#Delta(z_{conv} - z_{true})",1000,-500,500);
  TH1F* hdeltaz2FromTrueVertex = new TH1F("hdeltaz2FromTrueVertex","#Delta(z_{conv} - z_{true})",1000,-500,500);
 
  TH1F* hdzconv = new TH1F("hdzconv"," #Delta(z_{conv1} - z_{conv2})",200,-10,10);

  TH1F* hNvtxInDeltaZ1 = new TH1F("hNvtxInDeltaZ1","Number of vertices within deltaZ ",20,0,20);
  TH1F* hNvtxInDeltaZ2 = new TH1F("hNvtxInDeltaZ2","Number of vertices within deltaZ ",20,0,20);

  TH1F* hpho1_r9greater093_isConverted = new TH1F("hpho1_r9greater093_isConverted","hpho1_r9greater093_isConverted",5,-0.5,4.5);
  TH1F* hpho2_r9greater093_isConverted = new TH1F("hpho2_r9greater093_isConverted","hpho2_r9greater093_isConverted",5,-0.5,4.5);
  TH1F* hpho1_r9lower093_isConverted = new TH1F("hpho1_r9lower093_isConverted","hpho1_r9lower093_isConverted",5,-0.5,4.5);
  TH1F* hpho2_r9lower093_isConverted = new TH1F("hpho2_r9lower093_isConverted","hpho2_r9lower093_isConverted",5,-0.5,4.5);
  
  TH1F* hpho1_r9lower093_hasConversionCandidate = new TH1F("hpho1_r9lower093_hasConversionCandidate","hpho1_r9lower093_hasConversionCandidate",5,-0.5,4.5);
  TH1F* hpho2_r9lower093_hasConversionCandidate = new TH1F("hpho2_r9lower093_hasConversionCandidate","hpho2_r9lower093_hasConversionCandidate",5,-0.5,4.5);



  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  // algo parameters
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 0.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;
  vtxAlgoParams_.removeTracksInCone = 1;
  vtxAlgoParams_.coneSize           = 0.05;

  vtxAlgoParams_.sigmaPix  = 0.06;
  vtxAlgoParams_.sigmaTib  = 0.67;
  vtxAlgoParams_.sigmaTob  = 2.04;
  vtxAlgoParams_.sigmaFwd1 = 0.18;
  vtxAlgoParams_.sigmaFwd2 = 0.61;
  vtxAlgoParams_.sigmaFwd3 = 0.99;

 

  vector<string> ranksumpt2_;
  // variables order matters to resolve ties
  ranksumpt2_.push_back("logsumpt2");  

  vector<string> rankVariables_;
  // variables order matters to resolve ties
  rankVariables_.push_back("ptbal"), rankVariables_.push_back("ptasym"),rankVariables_.push_back("logsumpt2");  

  float ww = 1;

  

  //start loop over entries
  for (int u = 0; u < reader.GetEntries(); u++ )
    {
      if(u == entryMAX) break;
      if(u < entryMIN)  continue;
      if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
      reader.GetEntry(u);
      
      //--- use weights
      if (useWeights){
	std::vector<float>* PU_z       = reader.GetFloat("PU_z");
	//std::vector<float>* PU_z       = reader.GetFloat("PV_noMuon_z");
	int npu = PU_z->size(); // number of sim PU vtx
	float myrnd = gRandom->Uniform(0,nmax);
	//if (myrnd > w[npu]) continue;
	if (myrnd > w[npu+1]) continue;
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

      
      std::vector<float>* BS_x0 = reader.GetFloat("BS_x0");
      std::vector<float>* BS_y0 = reader.GetFloat("BS_y0");
      std::vector<float>* BS_z0 = reader.GetFloat("BS_z0");

      std::vector<ROOT::Math::XYZTVector>* photons = reader.Get4V("photons");
      std::vector<float>* photons_r9 = reader.GetFloat("photons_r9");
      std::vector<ROOT::Math::XYZTVector>* sc = reader.Get4V("photons_SC");
      std::vector<ROOT::Math::XYZVector>*  photons_SCpos    = reader.Get3V("photons_SCpos");
      
      std::vector<ROOT::Math::XYZVector>* photons_convVtx = reader.Get3V("photons_convVtx");
      std::vector<int>*  photons_convNtracks = reader.GetInt("photons_convNtracks");
      std::vector<int>* photons_convVtxIsValid = reader.GetInt("photons_convVtxIsValid");
      std::vector<float>* photons_convVtxChi2 = reader.GetFloat("photons_convVtxChi2");
      std::vector<float>* photons_convVtxNDOF = reader.GetFloat("photons_convVtxNDOF");
      std::vector<float>* photons_convEoverP  = reader.GetFloat("photons_convEoverP");
      

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

	hggSelection(mcV1, mcV2, photons, sc, photons_r9, accept, indpho1, indpho2);
	if (!accept) continue;
	
	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	TrueVertex_Z = mc_H_vertex->at(0).Z();

      }// end Hgg selection

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
	vtxx_[iv] =  0; //FIXME
	vtxy_[iv] =  0; //FIXME
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
      for (int ipho = 0 ; ipho < photons_SCpos->size(); ipho++){
      	phocalox_[ipho] = photons_SCpos->at(ipho).x();
	phocaloy_[ipho] = photons_SCpos->at(ipho).y();
	phocaloz_[ipho] = photons_SCpos->at(ipho).z();
	phoen_[ipho]    = sc->at(ipho).E();
      } 

      if ( (fabs(sc->at(indpho1).eta()) > etaEB && fabs(sc->at(indpho1).eta()) < etaEE) || fabs(sc->at(indpho1).eta()) > 2.5) continue;
      if ( (fabs(sc->at(indpho2).eta()) > etaEB && fabs(sc->at(indpho2).eta()) < etaEE) || fabs(sc->at(indpho2).eta()) > 2.5) continue;
   



      // set vtx info
      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);


      // set all variables needed for conversions
      const TVector3 bs(BS_x0->at(0),BS_y0->at(0),BS_z0->at(0));
      
      const TVector3 convVtx1(photons_convVtx->at(indpho1).x(),photons_convVtx->at(indpho1).y(),photons_convVtx->at(indpho1).z());
      const TVector3 convVtx2(photons_convVtx->at(indpho2).x(),photons_convVtx->at(indpho2).y(),photons_convVtx->at(indpho2).z());
      
      int iDet1 = 0;
      int iDet2 = 0;
      if ( fabs(sc->at(indpho1).eta()) < etaEB ) iDet1 = 1;
      if ( fabs(sc->at(indpho2).eta()) < etaEB ) iDet2 = 1;
       
      double convProb1 = TMath::Prob(double(photons_convVtxChi2->at(indpho1)),int(photons_convVtxNDOF->at(indpho1)));
      double convProb2 = TMath::Prob(double(photons_convVtxChi2->at(indpho2)),int(photons_convVtxNDOF->at(indpho2)));
      

      PhotonInfo pho1(TVector3(phocalox_[indpho1],phocaloy_[indpho1],phocaloz_[indpho1]),
		      bs,
		      convVtx1,
		      phoen_[indpho1],
		      iDet1,
		      photons_convNtracks->at(indpho1),
		      photons_convVtxIsValid->at(indpho1),
		      convProb1,
		      photons_convEoverP->at(indpho1)		      
		      ); 
            
      PhotonInfo pho2(TVector3(phocalox_[indpho2],phocaloy_[indpho2],phocaloz_[indpho2]),
		      bs,
		      convVtx2,
		      phoen_[indpho2],
		      iDet2,
		      photons_convNtracks->at(indpho2),
		      photons_convVtxIsValid->at(indpho2),
		      convProb2,
		      photons_convEoverP->at(indpho2)		      
		      ); 
  

      //------------------------------------------------------------------    
      // control histograms for conversions
      HggVertexFromConversions vAnaFromConv(vtxAlgoParams_);
      
      double z1 = 0;
      double dz1 = 0;
      double z2 = 0;
      double dz2 = 0;
 
      if (pho1.isAConversion()){
	z1  = vAnaFromConv.vtxZ(pho1);
	dz1 = vAnaFromConv.vtxdZ(pho1);
      }
      
      if (pho2.isAConversion()){
	z2  = vAnaFromConv.vtxZ(pho2);
	dz2 = vAnaFromConv.vtxdZ(pho2);
		
	if ( pho1.isAConversion()) 
	  hdzconv->Fill(z1 - z2);
      }
      
      
      if ( photons_r9->at(indpho1) > 0.93 ) hpho1_r9greater093_isConverted->Fill( pho1.isAConversion() );
      if ( photons_r9->at(indpho2) > 0.93 ) hpho2_r9greater093_isConverted->Fill( pho2.isAConversion() );

      if ( photons_r9->at(indpho1) < 0.93 ) {
	hpho1_r9lower093_isConverted->Fill( pho1.isAConversion() );
	
	if ( photons_convVtx->at(indpho1).x()!=-999 ){
	  hpho1_r9lower093_hasConversionCandidate->Fill(1);}
	else hpho1_r9lower093_hasConversionCandidate->Fill(0);
      }
      


      if ( photons_r9->at(indpho2) < 0.93 ) {
	hpho2_r9lower093_isConverted->Fill( pho2.isAConversion() );

	if ( photons_convVtx->at(indpho2).x()!= -999 ){
	  hpho2_r9lower093_hasConversionCandidate->Fill(1);}
	  else hpho2_r9lower093_hasConversionCandidate->Fill(0);
      }
      

      //--------------------------------------------------------------------------------------
      

      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);
      vAna.analyze(vinfo,pho1,pho2);
      
      
      // VERTEX SELECTION

      // find the PV closest to the MC Higgs vertex
      float dmin = 10000;
      int iClosest = -1;
      for ( int uu = 0; uu < nvtx_; uu++)
	{
	  float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	  if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	}


      //---- category 1 : R9 > 0.93 for both photons
      
      if (photons_r9->at(indpho1) >0.93 && photons_r9->at(indpho2) >0.93) {
	
	PtAll[0]->Fill( sum2pho.pt(),ww );
	EtaAll[0]->Fill( etaMaxSC ,ww);
	NvtAll[0]->Fill( nvtx_,ww );
	R9All[0]->Fill( photons_r9->at(indpho1) ,ww );
	
	// preselect vertexes 
	std::vector<int> presel;
	for(int i=0; i<nvtx_; i++) {
	  presel.push_back(i); 
	}
	vAna.preselection(presel);
	
       	// rankprod criterion
	vector<int> rankprod = vAna.rankprod(rankVariables_);
	if ( iClosest == rankprod[0]){
	  PtGood[0]->Fill( sum2pho.pt(),ww );
	  EtaGood[0]->Fill( etaMaxSC ,ww);
	  NvtGood[0]->Fill( nvtx_ ,ww);
	  R9Good[0]->Fill( photons_r9->at(indpho1) ,ww );
	}
      }	

      
      // category 2 : R9 < 0.93 for at least one photon, but no good conversions
      if ( (photons_r9->at(indpho1) <0.93 || photons_r9->at(indpho2) <0.93) &&  
	   !pho1.isAConversion() && !pho2.isAConversion())  
	{

	  PtAll[1]->Fill( sum2pho.pt(),ww );
	  EtaAll[1]->Fill( etaMaxSC ,ww);
	  NvtAll[1]->Fill( nvtx_,ww );
	  R9All[1]->Fill( photons_r9->at(indpho1) ,ww );
	  
	  // preselect vertexes 
	  std::vector<int> presel;
	  for(int i=0; i<nvtx_; i++) {
	    presel.push_back(i); 
	  }
	  vAna.preselection(presel);
	  
	  // rankprod criterion
	  vector<int> rankprod = vAna.rankprod(rankVariables_);
	  if ( iClosest == rankprod[0]){
	    PtGood[1]->Fill( sum2pho.pt(),ww );
	    EtaGood[1]->Fill( etaMaxSC ,ww);
	    NvtGood[1]->Fill( nvtx_ ,ww);
	    R9Good[1]->Fill( photons_r9->at(indpho1) ,ww );
	  }
	  
	}
      

      

      // category 3 : R9 < 0.93 for at least one photon, and good conversion
      if ( (photons_r9->at(indpho1) <0.93 || photons_r9->at(indpho2) <0.93) &&  
	   (pho1.isAConversion() || pho2.isAConversion()) )  
	{
	  
	  PtAll[2]->Fill( sum2pho.pt(),ww );
	  EtaAll[2]->Fill( etaMaxSC ,ww);
	  NvtAll[2]->Fill( nvtx_,ww );
	  R9All[2]->Fill( photons_r9->at(indpho1) ,ww );

	  float zconv  = 0;
	  float dzconv = 0;
	  
	  // 2 good conversions  ????
	  if ( pho1.isAConversion() && pho2.isAConversion() ){
	    zconv  = (z1 + z2)/2; // FIXME : usare la media pesata
	    dzconv = dz1;
	  }
	  
	  // 1 good conversion 
	  if ( pho1.isAConversion() && !pho2.isAConversion() ){
	    
	    zconv  = z1;
	    dzconv = dz1;
	    
	    hdz1->Fill(dz1);
	    hdeltaz1FromTrueVertex->Fill(z1 - TrueVertex_Z);
	    
	    int nVtxCloseToConvVtx1 = 0;
 	    for ( int uu = 0; uu < nvtx_; uu++){
	      float distt1 = fabs( PV_z->at(uu) - z1 );
	      if ( distt1 < dz1)   { nVtxCloseToConvVtx1++;}
	    }
	    hNvtxInDeltaZ1->Fill(nVtxCloseToConvVtx1);
	  
	  }
	  
	  if ( pho2.isAConversion() && !pho1.isAConversion() ){
	    zconv  = z2;
	    dzconv = dz2;
	    
	    hdz2->Fill(dz2);
	    hdeltaz2FromTrueVertex->Fill(z2 - TrueVertex_Z);
	    
	    int nVtxCloseToConvVtx2 = 0;      
	    for ( int uu = 0; uu < nvtx_; uu++){
	      float distt2 = fabs( PV_z->at(uu) - z2 );
	      if ( distt2 < dz2)   { nVtxCloseToConvVtx2++;}
	    }
	    
	    hNvtxInDeltaZ2->Fill(nVtxCloseToConvVtx2 );
	    
	  }

	  
	  // preselect vertexes 
	  std::vector<int> presel;
	  for(int i=0; i<nvtx_; i++) {
	    if ( fabs(zconv - PV_z->at(i)) < dzconv ) 
	      presel.push_back(i); 
	  }
	  vAna.preselection(presel);
	  
	  if (presel.size()!=0) {
	  
	  // rankprod criterion
	  vector<int> rankprod = vAna.rankprod(rankVariables_);
	  //std::cout << rankprod .size() << std::endl;
	  if ( iClosest == rankprod[0]){
	    PtGood[2]->Fill( sum2pho.pt(),ww );
	    EtaGood[2]->Fill( etaMaxSC ,ww);
	    NvtGood[2]->Fill( nvtx_ ,ww);
	    R9Good[2]->Fill( photons_r9->at(indpho1) ,ww );
	  }
	  }
	}
      
      
    }// end loop over entries



  // compute efficiencies
  TGraphAsymmErrors *EffVsNvtx[3];
  TGraphAsymmErrors *EffVsPt[3];
  TGraphAsymmErrors *EffVsEta[3];
  TGraphAsymmErrors *EffVsR9[3];
  
  for (int i=0; i< 3; i++){
    NvtAll[i]->Sumw2();
    NvtGood[i]->Sumw2();
    PtAll[i]->Sumw2();
    PtGood[i]->Sumw2();
    EtaAll[i]->Sumw2();
    EtaGood[i]->Sumw2();
    R9All[i]->Sumw2();
    R9Good[i]->Sumw2();

    EffVsNvtx[i] = new TGraphAsymmErrors();
    EffVsPt[i] = new TGraphAsymmErrors();
    EffVsEta[i] = new TGraphAsymmErrors();
    EffVsR9[i] = new TGraphAsymmErrors();

    EffVsNvtx[i]-> BayesDivide(NvtGood[i], NvtAll[i], "cp");
    EffVsPt[i]  -> BayesDivide(PtGood[i], PtAll[i], "cp");
    EffVsEta[i] -> BayesDivide(EtaGood[i], EtaAll[i], "cp"); 
    EffVsR9[i]  -> BayesDivide(R9Good[i], R9All[i], "cp");

  }




  /*
  EffVsNvtx_cat1.BayesDivide(NvtGood[0], NvtAll[0], "cp"); 

  TGraphAsymmErrors EffVsPt_cat1;
  EffVsPt_cat1.BayesDivide(PtGood[0], PtAll[0], "cp"); 
    
  TGraphAsymmErrors EffVsEta_cat1;
  
   
  TGraphAsymmErrors EffVsR9_cat1;
  EffVsR9_cat1.BayesDivide(R9Good[0], R9All[0], "cp"); 


  TGraphAsymmErrors EffVsNvtx_cat2;
  EffVsNvtx_cat2.BayesDivide(NvtGood[0], NvtAll[0], "cp"); 

  TGraphAsymmErrors EffVsPt_cat2;
  EffVsPt_cat2.BayesDivide(PtGood[0], PtAll[0], "cp"); 
    
  TGraphAsymmErrors EffVsEta_cat2;
  EffVsEta_cat2.BayesDivide(EtaGood[0], EtaAll[0], "cp"); 
   
  TGraphAsymmErrors EffVsR9_cat2;
  EffVsR9_cat2.BayesDivide(R9Good[0], R9All[0], "cp"); 
  */


  
  
   
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
  
  for (int i=0; i< 3; i++){
   
    NvtAll[i]->Write();
    NvtGood[i]->Write();

    PtAll[i]->Write();
    PtGood[i]->Write();

    EtaAll[i]->Write();
    EtaGood[i]->Write();

    R9All[i]->Write();
    R9Good[i]->Write();
  
    sprintf(hname,"EffVsNvtx_cat%d",i+1);
    EffVsNvtx[i]->Write(hname);

    sprintf(hname,"EffVsPt_cat%d",i+1);
    EffVsPt[i]->Write(hname);

    sprintf(hname,"EffVsEta_cat%d",i+1);
    EffVsEta[i]->Write(hname);

    sprintf(hname,"EffVsR9_cat%d",i+1);
    EffVsR9[i]->Write(hname);
    
  }
  
  
  hdz1->Write();
  hdz2->Write();
  hdeltaz1FromTrueVertex->Write();
  hdeltaz2FromTrueVertex->Write();
  hdzconv->Write();
  hNvtxInDeltaZ1->Write();
  hNvtxInDeltaZ2->Write();


  hpho1_r9greater093_isConverted->Write();
  hpho2_r9greater093_isConverted->Write();

  hpho1_r9lower093_hasConversionCandidate->Write();
  hpho2_r9lower093_hasConversionCandidate->Write();

  hpho1_r9lower093_isConverted->Write();
  hpho2_r9lower093_isConverted->Write();
  


  ff.Close();
  return 0; 
  
}













