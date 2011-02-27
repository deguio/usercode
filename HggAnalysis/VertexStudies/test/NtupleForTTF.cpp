
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
  
  std::string outputRootFilePath = gConfigParser -> readStringOption("Output::outputRootFilePath");
  std::string outputRootFileName = gConfigParser -> readStringOption("Output::outputRootFileName");  
  
  int entryMIN = gConfigParser -> readIntOption("Options::entryMIN");
  int entryMAX = gConfigParser -> readIntOption("Options::entryMAX");
  int isZee    = gConfigParser -> readIntOption("Options::isZee");
  int isZmumu  = gConfigParser -> readIntOption("Options::isZmumu");
  int isHiggs  = gConfigParser -> readIntOption("Options::isHiggs");

  double trackThr = gConfigParser -> readDoubleOption("Options::trackThr");

  std::cout << std::endl;
  std::cout << std::endl;


  //Filling the ntuple
  int isSig[100];

  int nVertices, nPVTracks[100];
  float vertexZ[100], sumPt2[100], modSumPt[100], modSumPtInCone_30[100], modSumPtInCone_45[100], deltaPhi_HSumPt[100], sum2PhoPt, normalizedChi2[100], sumModPt[100];
  
  int nTracks, trackPVIndex[1000];
  float trackPx[1000], trackPy[1000], trackPz[1000], trackDz[1000];

  TTree* outTree = new TTree("TMVA_vertexTree","TMVA_vertexTree");
  outTree -> Branch("nVertices",            &nVertices,             "nVertices/I");
  outTree -> Branch("vertexZ",              vertexZ,                "vertexZ[nVertices]/F");
  outTree -> Branch("sumPt2",               sumPt2,                 "sumPt2[nVertices]/F");
  outTree -> Branch("modSumPt",             modSumPt,               "modSumPt[nVertices]/F");
  outTree -> Branch("sumModPt",             sumModPt,               "sumModPt[nVertices]/F");
  outTree -> Branch("modSumPtInCone_30",    modSumPtInCone_30,      "modSumPtInCone_30[nVertices]/F");
  outTree -> Branch("modSumPtInCone_45",    modSumPtInCone_45,      "modSumPtInCone_45[nVertices]/F");
  outTree -> Branch("deltaPhi_HSumPt",      deltaPhi_HSumPt,        "deltaPhi_HSumPt[nVertices]/F");
  outTree -> Branch("sum2PhoPt",            &sum2PhoPt,              "sum2PhoPt/F");
  outTree -> Branch("nPVTracks",            nPVTracks,              "nPVTracks[nVertices]/I");


  outTree -> Branch("nTracks",              &nTracks,               "nTracks/I");
  outTree -> Branch("trackPx",              trackPx,                "trackPx[nTracks]/F");
  outTree -> Branch("trackPy",              trackPy,                "trackPy[nTracks]/F");
  outTree -> Branch("trackPz",              trackPz,                "trackPz[nTracks]/F");
  outTree -> Branch("trackDz",              trackDz,                "trackDz[nTracks]/F");
  outTree -> Branch("trackPVIndex",         trackPVIndex,           "trackPVIndex[nTracks]/I");

  outTree -> Branch("isSig", &isSig, "isSig[nVertices]/I");


  //Chain + histos
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

       //global variables
       float TrueVertex_Z;

       std::vector<float>* PV_normalizedChi2;
       std::vector<int>*   PV_ndof;
       std::vector<int>*   PV_nTracks;
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
	   for(unsigned int u=0; u < photons->size(); u++)
	     {
	       bool id =  PhotonId( photons->at(u).pt(), photons_SC->at(u).eta() , photons_ecalIso->at(u), photons_hcalIso->at(u), photons_hadronicOverEm->at(u), photons_trkSumPtHollowConeDR04->at(u), photons_sigmaIetaIeta->at(u) );
	       if( id && photons->at(u).pt() > 10 && indpho1 < 0 )
		 {ngood++; indpho1 = u;}
	       else if(id && photons->at(u).pt() > 10 && indpho2 < 0)
		 {ngood++; indpho2 = u;}
	       else if(id && photons->at(u).pt() > 10 )
		 {ngood++;}
	     }
	   if ( ngood != 2) continue;
	   
	   sum2pho = photons->at(indpho1)+ photons->at(indpho2);
	   if ( fabs(sum2pho.M() - 120) > 4 ) continue;

	 }//Hgg
       //---------------------------
       //--- set up Zee branches ---
       //---------------------------
       if ( isZee )
	 {
	   //taking Zee variables
	   std::vector<float>* eleid = reader.GetFloat("simpleEleId95cIso");

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
	   for( unsigned int uu = 0; uu < muons->size(); uu++)
	     {
	       if ( muons->at(uu).pt() > 10 && indmu1 < 0 )        { indmu1 = uu; ngood++; }
	       else if ( muons->at(uu).pt() > 10 && indmu2 < 0 )   { indmu2 = uu; ngood++; }
	       else if ( muons->at(uu).pt() > 10 )                 { ngood++;}
	     }
	   if ( ngood != 2){continue;}


	   sum2pho = muons->at(indmu1) + muons->at(indmu2);
	   if ( fabs(sum2pho.M() - 91) > 8 ) continue;

	   TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indmu1) + muons_dz_PV_noMuon->at(indmu2))/2.;
	   
	 }//Zmumu end
       
       // //vettore degli indici dei PV buoni (no splitting)
       // std::vector<int> goodIndex;
       // for (unsigned int vItr = 0; vItr < PV_z->size(); ++vItr)
       // 	 {
       // 	   if (PV_nTracks->at(vItr) > 3 || PV_SumPt2->at(vItr) > 10.) goodIndex.push_back(vItr);
       // 	 }

       // int npv = goodIndex.size();
       // if (npv == 0 ){continue;}
              
       
       //look if the H vertex matches one of the PV vertices
       float dmin = 10000;
       int iClosest = -1;
       int Vmatched = 0;
       for ( unsigned int uu = 0; uu < PV_z->size(); uu++)
	 {
	   float distt = fabs( PV_z->at(uu) - TrueVertex_Z );
	   if ( distt < 0.3 )   { Vmatched++; }
	   if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	 }
       
       //vertices loop
       for ( unsigned int uu = 0; uu < PV_z->size(); uu++)
	 {
	   nPVTracks[uu]     = 0;

	   ROOT::Math::XYZVector sumTracks;
	   ROOT::Math::XYZVector sumTracksInCone_30;
	   ROOT::Math::XYZVector sumTracksInCone_45;

	   sumModPt[uu] = 0;
	   sumPt2[uu]   = 0;
	   
	   for (unsigned int kk = 0; kk < PVtracks->size(); ++kk)
	     {
	       if ( PVtracks->at(kk).perp2() < trackThr*trackThr ) continue;
	       sumTracks     += PVtracks->at(kk);
	       sumModPt[uu]  += sqrt((PVtracks->at(kk)).perp2());
	       sumPt2[uu]    += PVtracks->at(kk).perp2();
	       if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*11./12. ) sumTracksInCone_30 += PVtracks->at(kk);
	       if ( deltaPhi(PVtracks->at(kk).phi(), sum2pho.phi()) > PI*7./8. )   sumTracksInCone_45   += PVtracks->at(kk);
	       
	       nPVTracks[uu] ++;
	       
	     }//PV tracks loop
	   
	   vertexZ[uu]            = PV_z->at(uu);
	   deltaPhi_HSumPt[uu]    = deltaPhi( sumTracks.phi(),sum2pho.phi());
	   modSumPt[uu]           = sqrt( sumTracks.perp2() );
	   modSumPtInCone_30[uu]  = sqrt( sumTracksInCone_30.perp2() );
	   modSumPtInCone_45[uu]  = sqrt( sumTracksInCone_45.perp2() );
	   normalizedChi2[uu]     = PV_normalizedChi2->at(uu);
	   

	   isSig[uu]   = 0;
	   if (uu == iClosest)
	     isSig[uu] = 1;
	   
	 }//vertices loop

       nVertices = PV_z->size();
       sum2PhoPt          = sum2pho.pt();	   



       //general tracks loop
       nTracks = tracks->size();
       for(unsigned int ii = 0; ii < tracks->size(); ++ii)
	 {
	   trackPx[ii] = tracks->at(ii).X();
	   trackPy[ii] = tracks->at(ii).Y();
	   trackPz[ii] = tracks->at(ii).Z();
	   trackDz[ii] = tracks_dz->at(ii);
	   trackPVIndex[ii] = (int)tracks_PVindex->at(ii);
	 }


       //one entry per event
       outTree -> Fill();	   
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

