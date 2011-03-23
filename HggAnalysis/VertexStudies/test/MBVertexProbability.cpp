
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

#define etaEB 1.444
#define etaEE 1.560

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta);

double fact (int n);


int main(int argc, char** argv)
{ 
   
  TChain* chain = new TChain("MiBiCommonNTTwoPhotons_005_8/SimpleNtuple");
  //chain->Add("/media/amassiro/VertexStudies/RelValMinBias_CMSSW_3_9_7-START39_V8-v1_GEN-SIM-RECO.root");
  chain->Add("/data2/VertexStudies/RelValMinBias_CMSSW_3_9_7-START39_V8-v1_GEN-SIM-RECO.root");
  
  treeReader reader((TTree*)(chain));
  
  int entryMAX = reader.GetEntries();
  
 
  // book histograms
  
  TH1F *hNvtx = new TH1F("hNvtx","number of vertices",75,0,25);
  hNvtx -> GetXaxis() -> SetTitle("Number of vertices}");

  TH1F *hTrkPt = new TH1F("hTrkPt","track pt",10000,0,100);
  hTrkPt -> GetXaxis() -> SetTitle("track p_{T} (GeV)");

  TH1F *hLog10TrkPt = new TH1F("hLog10TrkPt","Log10 track pt",20000,-2,2);
  hLog10TrkPt -> GetXaxis() -> SetTitle("log_{10}(p_{T})");

  TH1F *hPMBtrack = new TH1F("hPMBtrack","hPMBtrack",100,0,1);

  TH1F *hPMBvertex = new TH1F("hPMBvertex","hPMBvertex",100,0,1);
  TH2F *hPMB_vs_nTracks = new TH2F("hPMB_vs_nTracks","hPMB_vs_nTracks",1000,0,100,100,0,1);


  std::cout<<"found "<< entryMAX<<" entries"<<std::endl;
  
  float pTthr = 0.5;


  for (int u=0; u < entryMAX; u++ ){
    if(u%20000 == 0){std::cout<<"reading event "<< u <<std::endl;}
    reader.GetEntry(u);
    
    #include "../includeMe.h"
        
    int npv = PV_z->size();
    if ( npv > 0 && PV_nTracks->at(0) < 1) {npv = 0;}
    hNvtx->Fill( npv );
    if (npv == 0 ){continue;}

    for (int itrk=0;  itrk<PVtracks->size(); itrk++ ) { // loop over tracks in the vertices
      float pt = sqrt(PVtracks->at(itrk).perp2());
      hTrkPt -> Fill(pt);
      hLog10TrkPt -> Fill(log10(pt));
    }
    
  }// end loop over entries
  


  //----------




  float ntot = hTrkPt -> GetEntries();
  hTrkPt ->Scale(1./ntot);
  
  ntot = hLog10TrkPt -> GetEntries();
  hLog10TrkPt ->Scale(1./ntot);
  
  int nbins = hTrkPt -> GetNbinsX();
  int nbins10 = hLog10TrkPt -> GetNbinsX();

  TH1F *hprob = new TH1F("hprob","hprob",0,1,10000);

  for (int u=0; u < entryMAX; u++ ){
    if(u%20000 == 0){std::cout<<"reading event "<< u <<std::endl;}
    
    reader.GetEntry(u);
    
    #include "../includeMe.h"

    int npv = PV_z->size();
    if ( npv > 0 && PV_nTracks->at(0) < 1) {npv = 0;}
    hNvtx->Fill( npv );
    if (npv == 0 ){continue;}

    for (int iv = 0; iv < npv; iv++){ // loop over vertices
      
      float p = 1;
      for (int itrk=0;  itrk < PVtracks->size(); itrk++ ) { // loop over tracks in the vertex
	if (PVtracks_PVindex->at(itrk) == iv){
          float pt = sqrt(PVtracks->at(itrk).perp2());
	  if (pt > pTthr) {
	    int bin  = hLog10TrkPt -> FindBin(log10(pt));
	    int binThr = hLog10TrkPt -> FindBin(log10(pTthr));
	    float prob = hLog10TrkPt -> Integral(bin,nbins10)/hLog10TrkPt -> Integral(binThr,nbins10);
	    hPMBtrack -> Fill ( prob );
	    p = p * prob; // ??? perche'????
	  }
	}  
      }
      
      
      double mbp = 0;
      int nGoodTracks = 0;

      for (int itrk=0;  itrk < PVtracks->size(); itrk++ ) { // loop over tracks in the vertex
	if (PVtracks_PVindex->at(itrk) == iv){
	  float pt = sqrt(PVtracks->at(itrk).perp2());
	  if (pt > pTthr){
	    mbp += TMath::Power(-TMath::Log(p),nGoodTracks)/TMath::Gamma(nGoodTracks+1);  // ??? perche'????
	    nGoodTracks++;

	  }  
	}
      }
    
 
      mbp = mbp*p; // normalized product of probabilities
      //mbp = p; // simple product
      
      

      if (nGoodTracks!=0){
	hPMB_vs_nTracks -> Fill ( nGoodTracks,mbp);
	hPMBvertex-> Fill (mbp);
      }
    }
      
    
    
  }// end loop over entries
  




  
  TFile ff("MBVertexProbability.root","recreate");

  hNvtx->Write();
  hTrkPt->Write();
  hLog10TrkPt->Write();
  hPMB_vs_nTracks -> Write();
  hPMBtrack->Write();
  hPMBvertex->Write();

  ff.Close();
  return 0;
  

}

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta )

{

  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );

  return (iso && shape && HoE<0.05);

}



double fact(int n)

{
  int result = 1;
  for (int i = n; i>0; i--){
    result = result * i; 
  }

  
  return (result);
}
