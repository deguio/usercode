
///==== include ====

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


#include "HggVertexAnalyzer.h"
#include "../src/selection.cc"

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
    TFile weightsFile("testPUweights.root","READ");
    TH1F* hweights = (TH1F*)weightsFile.Get("hwdata");
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
  

  //histos

  TH1F PtAll("PtAll","Pt of boson all",80,0,400);
  TH1F PtGood("PtGood","Pt of boson good",80,0,400);
  TH1F PtGood_BDT("PtGood_BDT","Pt of boson good (BDT)",80,0,400);
  TH1F PtGood_RANK("PtGood_RANK","Pt of boson good (RANKING)",80,0,400);

  TH1F EtaAll("EtaAll","Eta of max SC",50,-5,5);
  TH1F EtaGood("EtaGood","Eta of max SC good",50,-5,5);
  TH1F EtaGood_BDT("EtaGood_BDT","Eta of max SC good (BDT)",50,-5,5);
  TH1F EtaGood_RANK("EtaGood_RANK","Eta of max SC good (RANKING)",50,-5,5);
  
  TH1F NvtAll("NvtAll","number of PV all",50,0,50);
  TH1F NvtGood("NvtGood","number of PV good",50,0,50);
  TH1F NvtGood_BDT("NvtGood_BDT","number of PV good (BDT)",50,0,50);
  TH1F NvtGood_RANK("NvtGood_RANK","number of PV good (RANKING)",50,0,50);

  TH1F bdth("bdtH"," bdt H",500,-1,1);
  TH1F bdtbkg("bdtBkg"," bdt bkg",500,-1,1);

  TH1F pt2h("pt2h","pt2 H",500,0,500);
  TH1F pt2bkg("pt2bkg","pt2 bkg",500,0,500);

  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;

  // algo parameters
  VertexAlgoParameters vtxAlgoParams_;
  vtxAlgoParams_.rescaleTkPtByError = false;
  vtxAlgoParams_.trackCountThr      = 0.;
  vtxAlgoParams_.highPurityOnly     = false;
  vtxAlgoParams_.maxD0Signif        = 9999999.;
  vtxAlgoParams_.maxDzSignif        = 9999999.;


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
	int npu = PU_z->size(); // number of sim PU vtx
	float myrnd = gRandom->Uniform(0,nmax);
	if (myrnd > w[npu]) continue;
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
	
	etaMaxSC = sc->at(indpho1).eta();
	sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	TrueVertex_Z = mc_H_vertex->at(0).Z();
      }// end Hgg selection


      // selections for Zee
      if (isZee){
	
	std::vector<float>* eleid = reader.GetFloat("simpleEleId95cIso");
	std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
	std::vector<float>* electrons_dz_PV_noEle = reader.GetFloat("electrons_dz_PV_noEle");

	PV_nTracks       = reader.GetInt("PV_noEle_nTracks");
	PV_z             = reader.GetFloat("PV_noEle_z");
	PV_d0            = reader.GetFloat("PV_noEle_d0");
	PVtracks         = reader.Get3V("PVEleLessTracks");
	PVtracks_PVindex = reader.GetInt("PVEleLessTracks_PVindex");
	tracks_PVindex   = reader.GetInt("tracks_PVindex");
	tracks_dxy_PV    = reader.GetFloat("tracks_dxy_PV");
	tracks_dz_PV     = reader.GetFloat("tracks_dz_PV");
	tracks_dz     = reader.GetFloat("tracks_dz");
	sc               = reader.Get4V("electrons_SC");
	
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
	   
	   zmumuSelection(muons, accept, indpho1, indpho2);

	   if (!accept) continue;
	   
	   etaMaxSC = muons->at(indpho1).eta();
	   sum2pho  = muons->at(indpho1)+ muons->at(indpho2);
	   TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indpho1) + muons_dz_PV_noMuon->at(indpho2))/2.;
	  	   
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

// 	for (int ii = 0 ; ii < tracks_dz_PV->size(); ii++){
// 	  if ( tracks_PVindex->at(ii) == PVtracks_PVindex->at(itrk)){
// 	    tkWeight_[itrk]= 1.;
// 	    tkd0_[itrk]    = 1;
// 	    tkd0Err_[itrk] = 1;
// 	    tkdz_[itrk]    = 1;
// 	    tkdzErr_[itrk] = 1;
// 	    tkIsHighPurity_[itrk]= 1.;
// 	  }
// 	}
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


      //if ( (fabs(sc->at(indpho1).eta()) > etaEB && fabs(sc->at(indpho1).eta()) < etaEE) || fabs(sc->at(indpho1).eta()) > 2.5) continue;
      //if ( (fabs(sc->at(indpho2).eta()) > etaEB && fabs(sc->at(indpho2).eta()) < etaEE) || fabs(sc->at(indpho2).eta()) > 2.5) continue;
   

      TupleVertexInfo vinfo( nvtx_, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
      
      PhotonInfo pho1(TVector3(phocalox_[indpho1],phocaloy_[indpho1],phocaloz_[indpho1]),phoen_[indpho1]); 
      PhotonInfo pho2(TVector3(phocalox_[indpho2],phocaloy_[indpho2],phocaloz_[indpho2]),phoen_[indpho2]); 
          
      HggVertexAnalyzer vAna(vtxAlgoParams_,nvtx_);

      vAna.analyze(vinfo,pho1,pho2);
      
      // preselect vertexes 
      std::vector<int> presel;
      for(int i=0; i<nvtx_; i++) {
	presel.push_back(i); 
      }
      
      vAna.preselection(presel);
      

      
      //look if the H vertex matches one of the PV vertices
      float dmin = 10000;
      int iClosest = -1;
      for ( int uu = 0; uu < nvtx_; uu++)
	{
	  float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	  if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	}

      
      //fill histos
      PtAll.Fill( sum2pho.pt(),ww );
      EtaAll.Fill( etaMaxSC ,ww);
      NvtAll.Fill( nvtx_,ww );
       

      // sumpt2 criterion
      if( iClosest == 0 ) {
	PtGood.Fill( sum2pho.pt(),ww );
	EtaGood.Fill( etaMaxSC ,ww);
	NvtGood.Fill( nvtx_ ,ww);
      }
              

      // ranking 
      vector<int> rankprod = vAna.rankprod(rankVariables_);
      if ( iClosest == rankprod[0]){
	PtGood_RANK.Fill( sum2pho.pt(),ww );
	EtaGood_RANK.Fill( etaMaxSC ,ww);
	NvtGood_RANK.Fill( nvtx_,ww );
      }

      // BDT - FIXME : this is a placeholder
      //vector<int> ranktmva = vAna.rankprod(*tmvaReader_,tmvaMethod_);//
      if ( iClosest == rankprod[0]){
	PtGood_BDT.Fill( sum2pho.pt(),ww );
	EtaGood_BDT.Fill( etaMaxSC ,ww);
	NvtGood_BDT.Fill( nvtx_ ,ww);
      }
      

   
    }// end loop over entries
  
  TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
  
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
    
  bdtbkg.Write();
  bdth.Write();
  
  pt2bkg.Write();
  pt2h.Write();
   


  ff.Close();
  return 0; 
  
}











