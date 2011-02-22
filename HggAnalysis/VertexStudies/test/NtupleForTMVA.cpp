
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

#define PI 3.141592654

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta);

int main(int argc, char** argv)
{ 
  int isSig;
  float sumPt2, sumPt, sumPtMod, sumPtModInCone_30, sumPtModInCone_45, deltaPhi_HSumPt, sum2PhoPt, tracksNum, normalizedChi2;
  TTree* outTree = new TTree("TMVA_vertexTree","TMVA_vertexTree");
  outTree -> Branch("sumPt2",   &sumPt2,    "sumPt2/F");
  outTree -> Branch("sumPt",    &sumPt,     "sumPt/F");
  outTree -> Branch("sumPtMod", &sumPtMod,  "sumPtMod/F");
  outTree -> Branch("sumPtModInCone_30", &sumPtModInCone_30,  "sumPtModInCone_30/F");
  outTree -> Branch("sumPtModInCone_45", &sumPtModInCone_45,  "sumPtModInCone_45/F");
  outTree -> Branch("deltaPhi_HSumPt", &deltaPhi_HSumPt,  "deltaPhi_HSumPt/F");
  outTree -> Branch("sum2PhoPt",  &sum2PhoPt,   "sum2PhoPt/F");
  outTree -> Branch("tracksNum",&tracksNum, "tracksNum/F");
  outTree -> Branch("normalizedChi2",&normalizedChi2, "normalizedChi2/F");

  outTree -> Branch("isSig", &isSig, "isSig/I");


  TChain* chain = new TChain("MiBiCommonNTTwoPhotons_005_8/SimpleNtuple");
  chain->Add("/data2/VertexStudies/GluGluToHToGG_M-120_7TeV-powheg-pythia6_Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1_AODSIM_merged.root");
  treeReader reader((TTree*)(chain));
  
  int entryMAX = reader.GetEntries();
  //int entryMAX = 10000;

   std::cout<<"found "<< entryMAX<<" entries"<<std::endl;
   


   for (int u=0; u < entryMAX; u++ )
     {
       if(u%20000 == 0){std::cout<<"reading event "<< u <<std::endl;}
       reader.GetEntry(u);
       #include "../includeMe.h"
       std::vector<ROOT::Math::XYZVector>* MCvertex = reader.Get3V("mc_H_vertex");
       TClonesArray* phoRechit_E = reader.GetTClonesArray("photons_rechitE");  
       
       //bkg selection
       int indpho1 = -100, indpho2 = -100;
       int ngood = 0;
       for(unsigned int u=0; u < photons->size(); u++)
	 {
	   bool id =  PhotonId( photons->at(u).pt(), photons_SC->at(u).eta() , photons_ecalIso->at(u), photons_hcalIso->at(u), photons_hadronicOverEm->at(u), photons_trkSumPtHollowConeDR04->at(u), photons_sigmaIetaIeta->at(u) );
	   if( id && photons->at(u).pt() > 10 && indpho1 < 0 ) 
	     {
	       ngood++;
	       indpho1 = u;
	     }
	   else if(id && photons->at(u).pt() > 10 && indpho2 < 0)
	     {
	       ngood++;
	       indpho2 = u;
	     }
	   else if(id && photons->at(u).pt() > 10 )
	     {	     
	       ngood++;
	     }
	 }
       if ( ngood != 2) continue;

       ROOT::Math::XYZTVector sum2pho = photons->at(indpho1)+ photons->at(indpho2);

       int npv = PV_z->size();
       if ( npv > 0 && PV_nTracks->at(0) < 1) {npv = 0;}
       if (npv == 0 ) continue;
    
       //look if the H vertex matches one of the PV vertices
       float dmin = 10000;
       int iClosest = -1;
       for ( int uu = 0; uu < npv; uu++)
	 {
	   float distt = fabs( PV_z->at(uu) - MCvertex->at(0).Z() );
	   if ( distt < dmin)
	     {
	       dmin = distt;
	       iClosest = uu;
	     }
	 }
       
       if (dmin > 0.3) continue;
       
       for ( int uu = 0; uu < npv; uu++)
	 {
	   //BKG TTree fill

	     ROOT::Math::XYZVector sumTracks;
	     ROOT::Math::XYZVector sumTracksInCone_30;
	     ROOT::Math::XYZVector sumTracksInCone_45;
	     for (unsigned int kk = 0; kk < PVtracks->size(); ++kk)
	       if (PVtracks_PVindex->at(kk) == uu)
		 {
		   sumTracks += PVtracks->at(kk);
		   if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*11./12. ) sumTracksInCone_30 += PVtracks->at(kk);
		   if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*7./8. ) sumTracksInCone_45 += PVtracks->at(kk);
		 }

	     //variable filling
	     sumPt2 = PV_SumPt2->at(uu);
	     sumPt  = PV_SumPt->at(uu);
	     deltaPhi_HSumPt = deltaPhi( sumTracks.phi(),sum2pho.phi());
	     sumPtMod = sqrt( sumTracks.perp2() );
	     sumPtModInCone_30 = sqrt( sumTracksInCone_30.perp2() );
	     sumPtModInCone_45 = sqrt( sumTracksInCone_45.perp2() );
	     sum2PhoPt = sum2pho.pt();
	     tracksNum = PV_nTracks->at(uu);
	     normalizedChi2 = PV_normalizedChi2->at(uu);	 

	     isSig = 0;
	     if (uu == iClosest)
	       isSig = 1;

	     outTree -> Fill();	   
	 }
     }
   

   TFile fout("NtupleForTMVA.root","RECREATE");
   
   outTree -> Write();

   fout.Close();
   return 0;
  

}

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta )

{

  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );

  return (iso && shape && HoE<0.05);

}

