
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

void zeeSelection (std::vector<ROOT::Math::XYZTVector>* electrons ,
		   std::vector<float>* eleid,
		   int& passSelection, int& i1, int& i2
		   );

void zmumuSelection (std::vector<ROOT::Math::XYZTVector>* muons ,
		     int& passSelection, int& i1, int& i2
		     );

//------MAIN------

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

  std::cout << std::endl;
  std::cout << std::endl;


  //Chain + histos
  TChain* chain = new TChain(treeName.c_str());
  chain->Add((baseDir+inputFile).c_str());

  treeReader reader((TTree*)(chain));
  
  TH1F PtAll("PtAll","Pt of boson all",40,0,200);
  TH1F PtGood("PtGood","Pt of boson good",40,0,200);
  TH1F PtGood_BDT("PtGood_BDT","Pt of boson good BDT",40,0,200);

  TH1F EtaAll("EtaAll","max photon Eta ",50,-5,5);
  TH1F EtaGood("EtaGood","max photon eta good",50,-5,5);
  TH1F EtaGood_BDT("EtaGood_BDT","max photon Eta good BDT",50,-5,5);
   
  TH1F NvtAll("NvtAll","number of PV all",20,0,20);
  TH1F NvtGood("NvtGood","number of PV good",20,0,20);
  TH1F NvtGood_BDT("NvtGood_BDT","number of PV good BDT",20,0,20);
 
  TH2F PtAll_vsNvtAll   ("PtAll_vsNvtAll","PtAll_vsNvtAll",20,0,20,40,0,200);
  TH2F PtGood_vsNvtGood ("PtGood_vsNvtGood","PtGood_vsNvtGood",20,0,20,40,0,200);
  TH2F PtBDT_vsNvtBDT   ("PtBDT_vsNvtBDT","PtBDT_vsNvtBDT",20,0,20,40,0,200);
  TH1F bdth("bdtH"," bdt H",500,-1,1);
  TH1F bdtbkg("bdtBkg"," bdt bkg",500,-1,1);

  TH1F pt2h("pt2h","pt2 H",500,0,500);
  TH1F pt2bkg("pt2bkg","pt2 bkg",500,0,500);

  TH1F nmatched("nmatched","nm",100,0,10);


  TH1F *hMBTrackProbPV  = new TH1F("hMBTrackProbPV","MB track probability for PV tracks",100,0,1);
  TH1F *hMBTrackProbMB  = new TH1F("hMBTrackProbMB","MB track probability for MB tracks",100,0,1);
  TH1F *hMBVertexProbPV = new TH1F("hMBVertexProbPV","MB vertex probability for PV tracks",100,0,1);
  TH1F *hMBVertexProbMB = new TH1F("hMBVertexProbMB","MB vertex probability for MB tracks",100,0,1);

  std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;
  
  // read prob. density ofr MB tracks from file 
  TFile *fMB = TFile::Open("MBVertexProbability_MBfromHgg.root");
  TH1F *hLog10TrkPt = (TH1F*)fMB->Get("hLog10TrkPt"); 
  int nbins10 = hLog10TrkPt -> GetNbinsX();
  int binThr = hLog10TrkPt -> FindBin(log10(trackThr));

   //start loop over entries
   for (int u = 0; u < reader.GetEntries(); u++ )
     {
       if(u == entryMAX) break;
       if(u < entryMIN)  continue;
       if(u%10000 == 0) std::cout<<"reading event "<< u <<std::endl;
       reader.GetEntry(u);

       //setup common branches
       #include "../includeMe_forTMVA_check.h"
       
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
       

       int accept  = 0;

       ROOT::Math::XYZTVector sum2pho;
       float etaMaxSC ;
       float TrueVertex_Z;

       //---------------------------
       //--- set up Hgg branches ---
       //---------------------------
       if ( isHiggs )
	 {
	   //taking H variables
	   if (mc_H_vertex->size() != 1) continue;

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

	   int indpho1 = -100;
	   int indpho2 = -100;

	   hggSelection(mcV1, mcV2, photons, photons_SC, photons_r9, accept, indpho1, indpho2);

	   if (!accept) continue;
	   
	   etaMaxSC = photons_SC->at(indpho1).eta();
	   sum2pho  = photons->at(indpho1)+ photons->at(indpho2);
	   TrueVertex_Z = mc_H_vertex->at(0).Z();	   

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
	   
	   int indele1 = -100;
	   int indele2 = -100;
       
	   zeeSelection(electrons, eleid, accept, indele1, indele2);

	   if (!accept) continue;
	   
	   etaMaxSC = electrons_SC->at(indele1).eta();
	   sum2pho  = electrons->at(indele1)+ electrons->at(indele2);
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
	   

	   int indmu1 = -100;
	   int indmu2 = -100;

	   zmumuSelection(muons, accept, indmu1, indmu2);

	   if (!accept) continue;
	   
	   etaMaxSC = muons->at(indmu1).eta();
	   sum2pho  = muons->at(indmu1)+ muons->at(indmu2);
	   TrueVertex_Z = PV_z->at(0) + (muons_dz_PV_noMuon->at(indmu1) + muons_dz_PV_noMuon->at(indmu2))/2.;
	   
	 }//Zmumu end
       
       //vettore degli indici dei PV buoni (not splitted)
       std::vector<int> goodIndex;
       for (unsigned int vItr = 0; vItr < PV_z->size(); ++vItr)
	 {
	   //	   if (PV_nTracks->at(vItr) > 3 || PV_SumPt2->at(vItr) > 10.) 
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
	   if ( distt < 0.6 )   { Vmatched++; }
	   if ( distt < dmin)   { dmin = distt; iClosest = uu; }
	 }
       nmatched.Fill(Vmatched);
       

       PtAll.Fill( sum2pho.pt() );
       EtaAll.Fill( etaMaxSC );
       NvtAll.Fill(npv);
       PtAll_vsNvtAll.Fill(npv,sum2pho.pt());
       
       //is sumPt2 the good choice?
       if( iClosest == 0 ) {
	 PtGood.Fill( sum2pho.pt() );
	 EtaGood.Fill( etaMaxSC );
	 NvtGood.Fill(npv);
	 PtGood_vsNvtGood.Fill(npv,sum2pho.pt());
       }
              
       // Use Minimum Bias Probability

       int MBPind = -1;
       float MBPmin = 999999.;
       
       for  (int iv = 0; iv < npv; iv++){ // loop over vertices
	 
	 float p  = 1;
	 
	 for (int itrk=0;  itrk<PVtracks->size(); itrk++ ) { // loop over tracks in the vertex
	   if (PVtracks_PVindex->at(itrk) == goodIndex[iv]){
	     float pt = sqrt(PVtracks->at(itrk).perp2());
	     if (pt > trackThr){  
	       int bin    = hLog10TrkPt -> FindBin(log10(pt));
	       float prob = hLog10TrkPt -> Integral(bin,nbins10)/hLog10TrkPt -> Integral(binThr,nbins10);
	       
	       if ( goodIndex[iv] == iClosest ){
		 hMBTrackProbPV->Fill(prob);
	       }
	       else {
		 hMBTrackProbMB->Fill(prob);
	       }
	       p = p * prob;
	     }
	   }
	 }

	 double mbp = 0;
	 int nGoodTracks = 0;
	 bool goodvtxfound = false;
	 	 
	 for (int itrk=0;  itrk < PVtracks->size(); itrk++ ) { // loop over tracks in the vertex
	   if (PVtracks_PVindex->at(itrk) == goodIndex[iv]){
	     goodvtxfound = true;
	     float pt = sqrt(PVtracks->at(itrk).perp2());
	     if (pt > trackThr){
	       mbp += TMath::Power(-TMath::Log(p),nGoodTracks)/TMath::Gamma(nGoodTracks+1);  // ??? perche'????
	       nGoodTracks++;	 
	     }  
	   }
	 }// 
	   
	   
	 if (nGoodTracks!=0){
	   if ( goodIndex[iv] == iClosest ){
	     hMBVertexProbPV->Fill(p*mbp);
	   }
	   else {
	     hMBVertexProbMB->Fill(p*mbp);
	   }
	 }
	 
	 //mbp = p; // simple product of probabilities
	 mbp = mbp*p; // normalized product
	 	 
	 if (mbp < MBPmin && goodvtxfound && nGoodTracks!=0){ // quando nGoodTracks e' zero lo considero MB
 	   MBPmin = mbp;
 	   MBPind = iv;
 	 }
	 
       } // end loop over vertices
     
       if (MBPind == -1 ) continue;

       //if( fabs(PV_z->at(goodIndex[MBPind]) - TrueVertex_Z ) < 0.6 )
       if( iClosest == goodIndex[MBPind] )
	 {
	   PtGood_BDT.Fill( sum2pho.pt() );
	   EtaGood_BDT.Fill( etaMaxSC );
	   NvtGood_BDT.Fill( npv );
 	   PtBDT_vsNvtBDT.Fill( npv,sum2pho.pt() );
 	 }
       
     } //evt loop
   
   
   TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
   
   PtAll.Write();
   PtGood.Write();
   PtGood_BDT.Write();

   EtaAll.Write();
   EtaGood.Write();
   EtaGood_BDT.Write();
  
   NvtAll.Write();
   NvtGood.Write();
   NvtGood_BDT.Write(); 

   PtAll_vsNvtAll.Write();   
   PtGood_vsNvtGood.Write();
   PtBDT_vsNvtBDT.Write();

   hMBTrackProbPV  -> Write();
   hMBTrackProbMB  -> Write();
   hMBVertexProbPV -> Write();
   hMBVertexProbMB -> Write();
    


   bdtbkg.Write();
   bdth.Write();
   
   pt2bkg.Write();
   pt2h.Write();
   
   nmatched.Write();

   ff.Close();
   return 0;
   
   
}

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta ) 
{

  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );

  return (iso && shape && HoE<0.05);
}




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





void zeeSelection (std::vector<ROOT::Math::XYZTVector>* electrons ,
		   std::vector<float>* eleid,
		   int& passSelection, int& i1, int& i2
		   )

{
  passSelection = 0;
  int index1 = -100, index2 = -100;
  int ngood = 0;
  for(unsigned int uu = 0; uu < eleid->size(); uu++)
    {
      if (  eleid->at(uu) == 7 && index1 < 0 )      {index1 = uu; ngood++;}
      else if ( eleid->at(uu) == 7 && index2 < 0 )  {index2 = uu; ngood++;}
      else if ( eleid->at(uu) == 7 )                { ngood++;}
    }
  
  if ( ngood != 2){
    passSelection = 0;
    i1 = -100;
    i2 = -100;
    return;
  }
  
  if (electrons->at(index1).E() > electrons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  

  ROOT::Math::XYZTVector v = electrons->at(index1) + electrons->at(index2);
  bool pass_mcut      = (fabs( v.M()  - 91. ) > 8. ) ;
  bool pass_kincuts   = (electrons->at(index1).pt()> 40.) && ( electrons->at(index2).pt() > 30.) ;
  
  if( pass_mcut && pass_kincuts ) {
    passSelection = 1;    
  }

}


//------------------------------------------------------------------------

void zmumuSelection (std::vector<ROOT::Math::XYZTVector>* muons ,
		     int& passSelection, int& i1, int& i2
		     )

{
  passSelection = 0;
  int index1 = -100, index2 = -100;
  int ngood = 0;

  for( int uu = 0; uu < muons->size(); uu++)
    {
      if ( muons->at(uu).pt() > 10 && index1 < 0 )        { index1 = uu; ngood++; }
      else if ( muons->at(uu).pt() > 10 && index2 < 0 )   { index2 = uu; ngood++; }
      else if ( muons->at(uu).pt() > 10 )                 { ngood++;}
    }
    
  if ( ngood != 2){
    passSelection = 0;
    i1 = -100;
    i2 = -100;
    return;
  }
  
  if (muons->at(index1).E() > muons->at(index2).E())
    {
      i1 = index1;
      i2 = index2;
    }
  else
    {
      i1 = index2;
      i2 = index1;
    }
  

  ROOT::Math::XYZTVector v = muons->at(index1) + muons->at(index2);
  bool pass_mcut      = (fabs( v.M() - 91. ) > 8. ) ;
  bool pass_kincuts   = (muons->at(index1).pt()> 40.) && ( muons->at(index2).pt() > 30.) ;
  
  if( pass_mcut && pass_kincuts ) {
    passSelection = 1;    
  }

}
