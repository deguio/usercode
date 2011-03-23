
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

void hggSelection (std::vector<ROOT::Math::XYZTVector>* mcV1, 
		   std::vector<ROOT::Math::XYZTVector>* mcV2, 
		   std::vector<ROOT::Math::XYZTVector>* photons,
		   std::vector<ROOT::Math::XYZTVector>* photons_SC,
		   std::vector<float>* photons_r9,
		   int& passSelection, int& i1, int& i2);


double fact (int n);


//// ------- MAIN --------- 

int main(int argc, char** argv)
{ 
   
  TChain* chain = new TChain("MiBiCommonNTTwoPhotons_005_8/SimpleNtuple");
  chain->Add("/data2/VertexStudies/MC_23022011/GluGluToHToGG_M-120_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1_AODSIM_merged.root");
  
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

  //for (int u=0; u < entryMAX; u++ ){
  for (int u=0; u < 50000; u++ ){
    if(u%20000 == 0){std::cout<<"reading event "<< u <<std::endl;}
    reader.GetEntry(u);
    

    #include "../includeMe.h"

    //std::cout << " ciao "<<std::endl;        
    int npv = PV_z->size();
    if ( npv > 0 && PV_nTracks->at(0) < 1) {npv = 0;}
    hNvtx->Fill( npv );
    if (npv == 0 ){continue;}

    int accept;
    int indpho1 = -100;
    int indpho2 = -100;
    
    hggSelection(mcV1, mcV2, photons, photons_SC, photons_r9, accept, indpho1, indpho2);
    if (!accept) continue;

    //look if the H vertex matches one of the PV vertices
    if (mc_H_vertex->size() != 1) continue;
    float    TrueVertex_Z = mc_H_vertex->at(0).Z();
    
    float dmin = 10000;
    int iClosest = -1;
    for ( int uu = 0; uu < npv; uu++) {
      float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
      if ( distt < dmin)   { dmin = distt; iClosest = uu; }
    }
    

    for (int itrk=0;  itrk<PVtracks->size(); itrk++ ) { // loop over tracks in the vertices
      if (PVtracks_PVindex->at(itrk) != iClosest){
	float pt = sqrt(PVtracks->at(itrk).perp2());
	hTrkPt -> Fill(pt);
	hLog10TrkPt -> Fill(log10(pt));
      }
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

    //selections
    int accept;
    int indpho1 = -100;
    int indpho2 = -100;
    
    hggSelection(mcV1, mcV2, photons, photons_SC, photons_r9, accept, indpho1, indpho2);
    if (!accept) continue;

    
    //look if the H vertex matches one of the PV vertices
    if (mc_H_vertex->size() != 1) continue;
    float TrueVertex_Z = mc_H_vertex->at(0).Z();

    
    float dmin = 10000;
    int iClosest = -1;
    for ( int uu = 0; uu < npv; uu++)
      {
	float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	if ( distt < dmin)   { dmin = distt; iClosest = uu; }
      }
    
    
    for (int iv = 0; iv < npv; iv++){ // loop over vertices
      
      float p = 1;
      for (int itrk=0;  itrk < PVtracks->size(); itrk++ ) { // loop over tracks in the vertex
	if (PVtracks_PVindex->at(itrk) == iv && PVtracks_PVindex->at(itrk) != iClosest){
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
	if (PVtracks_PVindex->at(itrk) == iv && PVtracks_PVindex->at(itrk) != iClosest){
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
  




  
  TFile ff("MBVertexProbability_MBfromHgg.root","recreate");

  hNvtx->Write();
  hTrkPt->Write();
  hLog10TrkPt->Write();
  hPMB_vs_nTracks -> Write();
  hPMBtrack->Write();
  hPMBvertex->Write();

  ff.Close();
  return 0;
  

}


//-------------------------------------------------------------------------------------------------
bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta )
{

  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );

  return (iso && shape && HoE<0.05);

}
//-------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------
double fact(int n)
{
  int result = 1;
  for (int i = n; i>0; i--){
    result = result * i; 
  }
  return (result);
}
//-------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------
void hggSelection (std::vector<ROOT::Math::XYZTVector>* mcV1, 
		   std::vector<ROOT::Math::XYZTVector>* mcV2, 
		   std::vector<ROOT::Math::XYZTVector>* photons ,
		   std::vector<ROOT::Math::XYZTVector>* photons_SC,
		   std::vector<float>* photons_r9,
		   int& passSelection, int& i1, int& i2
		   )
{
  float photonsMC_eta[2], photonsMC_phi[2];
  
  if ( mcV1->size() != 1 ||  mcV2->size() != 1 ||  photons->size() == 0){
    passSelection = 0;
    i1 = -100;
    i2=-100;
    return;
  }
  
  
  // MC matching 
  if ( mcV1->at(0).E() > mcV2->at(0).E()) {
    photonsMC_eta[0] = mcV1->at(0).Eta();
    photonsMC_eta[1] = mcV2->at(0).Eta();
    photonsMC_phi[0] = mcV1->at(0).Phi();
    photonsMC_phi[1] = mcV2->at(0).Phi();
  }
  else      {
    photonsMC_eta[1] = mcV1->at(0).Eta();
    photonsMC_eta[0] = mcV2->at(0).Eta();
    photonsMC_phi[1] = mcV1->at(0).Phi();
    photonsMC_phi[0] = mcV2->at(0).Phi();
  }
  
  int index1 = -100, index2 = -100;
  double dR_1_min = 10000.;
  double dR_2_min = 10000.;
  for(unsigned int u=0; u < photons->size(); u++)
    {
      double dR_1 = deltaR( photonsMC_eta[0], photonsMC_phi[0], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
      if (dR_1 < dR_1_min) { dR_1_min = dR_1; index1 = u; }
      
      double dR_2 = deltaR( photonsMC_eta[1], photonsMC_phi[1], photons_SC->at(u).eta(), photons_SC->at(u).phi() );
      if (dR_2 < dR_2_min) { dR_2_min = dR_2; index2 = u; }
      
    }

  
  if (photons->at(index1).E() > photons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  


  bool is_mcmatched   = (dR_1_min < 0.15) && (dR_2_min < 0.15);
  bool is_unconverted = (photons_r9->at(index1) > 0.93) &&  (photons_r9->at(index2) > 0.93);
  bool pass_kincuts   = (photons->at(index1).pt()> 40.) && ( photons->at(index2).pt() > 30.) ;
  
  if( is_mcmatched && is_unconverted && pass_kincuts ) {
    passSelection = 1;    
  }


}
//-------------------------------------------------------------------------------------------------
