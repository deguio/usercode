
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
  
  TH1F NvtAll("NvtAll","number of PV all",20,0,20);
  TH1F NvtGood("NvtGood","number of PV good",20,0,20);
  TH1F NvtGood_BDT("NvtGood_BDT","number of PV good BDT",20,0,20);

  TH1F EtaAll("EtaAll","Eta of boson all",100,-5,5);
  TH1F EtaGood("EtaGood","Eta of boson good",100,-5,5);
  TH1F EtaGood_BDT("EtaGood_BDT","Eta of boson good BDT",100,-5,5);

  TH2F PtAll_vsNvtAll   ("PtAll_vsNvtAll","PtAll_vsNvtAll",20,0,20,40,0,200);
  TH2F PtGood_vsNvtGood ("PtGood_vsNvtGood","PtGood_vsNvtGood",20,0,20,40,0,200);
  TH2F PtBDT_vsNvtBDT   ("PtBDT_vsNvtBDT","PtBDT_vsNvtBDT",20,0,20,40,0,200);

  TH1F bdth("bdtH"," bdt H",500,-1,1);
  TH1F bdtbkg("bdtBkg"," bdt bkg",500,-1,1);

  TH1F pt2h("pt2h","pt2 H",500,0,500);
  TH1F pt2bkg("pt2bkg","pt2 bkg",500,0,500);

  TH1F nmatched("nmatched","nm",100,0,10);

  TH1F sum2Pho_overMC("sum2Pho_overMC","Sum2Pho/mcH",180,-3.,3.);

   std::cout<<"found "<< reader.GetEntries() <<" entries"<<std::endl;
   
    //
   // create the Reader object
   //
   TMVA::Reader *TMVAreader = new TMVA::Reader( "!Color:!Silent" );    
   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use

   Float_t sumPt2_TMVA, tracksNum_TMVA, deltaPhi_HSumPt_TMVA, sumPtModOversum2PhoPt_TMVA, sum2PhoPt_TMVA;
   TMVAreader->AddVariable( "sumPt2"          , &sumPt2_TMVA);
   TMVAreader->AddVariable( "nTracks"       , &tracksNum_TMVA);
   TMVAreader->AddVariable( "deltaPhi_HSumPt" , &deltaPhi_HSumPt_TMVA);
   TMVAreader->AddVariable( "sum2PhoPt"       , &sum2PhoPt_TMVA);
   TMVAreader->AddVariable( "modSumPt_sum2PhoPt" , &sumPtModOversum2PhoPt_TMVA);


   TMVAreader->BookMVA( "BDTmethod", weights.c_str() ); 

   

   // correction function for Zee
   // TFile *fr = TFile::Open("hRatio.root");
   //TH1F *hRatio = (TH1F*)fr->Get("hRatio");  
   //float w = 1;
   
   TFile *fr = TFile::Open("correction.root");
   TH2F *SumPt2_vs_BosonPt_Hgg = (TH2F*)fr->Get("SumPt2_vs_BosonPt_Hgg");  
   TH2F *SumPt2_vs_BosonPt_Zee = (TH2F*)fr->Get("SumPt2_vs_BosonPt_Zee");  
   float w = 1;
   
   float lowEdge[10]  = {0. , 10., 20., 30., 40., 50., 70. , 100., 150., 200.};
   float highEdge[10] = {10., 20., 30., 40., 50., 70., 100., 150., 200., 1000.};

   TH1D* hHggSumPt2_bin[10];
   TH1D* hZeeSumPt2_bin[10];
   char hname[100];
   for (int i = 0; i < 10; i++){
     sprintf(hname,"hHggSumPt2_bin_%d",i);
     hHggSumPt2_bin[i] = SumPt2_vs_BosonPt_Hgg->ProjectionY(hname,lowEdge[i],highEdge[i]);
     hHggSumPt2_bin[i]->Scale(1./hHggSumPt2_bin[i]->Integral());
     sprintf(hname,"hZeeSumPt2_bin_%d",i);
     hZeeSumPt2_bin[i] = SumPt2_vs_BosonPt_Zee->ProjectionY(hname,lowEdge[i],highEdge[i]);
     hZeeSumPt2_bin[i]->Scale(1./hZeeSumPt2_bin[i]->Integral());
   }

  

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

	   if (mcH->at(0).Pt() > 5.) sum2Pho_overMC.Fill( sum2pho.Pt()/mcH->at(0).Pt() );

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
	   for( int uu = 0; uu < muons->size(); uu++)
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
       
       //vettore degli indici dei PV buoni (not splitted)
       std::vector<int> goodIndex;
       for (unsigned int vItr = 0; vItr < PV_z->size(); ++vItr)
	 {
	   if (PV_nTracks->at(vItr) > 3 || PV_SumPt2->at(vItr) > 10.) goodIndex.push_back(vItr);
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


       // compute logsumpt2
       float mysumpt2 = 0;
       for (unsigned int kk = 0; kk < PVtracks->size(); ++kk)
	 {
	   if (PVtracks_PVindex->at(kk) == goodIndex[iClosest])
	     {
	       if ( PVtracks->at(kk).perp2() < trackThr*trackThr ) continue;
	       mysumpt2 += PVtracks->at(kk).perp2();
	     }
	 }//tracks loop
       
       
       int ptindex = 0; 
       for (int jj = 0; jj < 10; jj++){
	 if ( sum2pho.pt()>lowEdge[jj] && sum2pho.pt()<highEdge[jj] ) ptindex = jj;
	   }

       w = hHggSumPt2_bin[ptindex]->Interpolate(mysumpt2) / hZeeSumPt2_bin[ptindex]->Interpolate(mysumpt2);
       
       PtAll.Fill( sum2pho.pt(),w );
       EtaAll.Fill( sum2pho.eta(),w );
       NvtAll.Fill(npv,w);
       PtAll_vsNvtAll.Fill(npv,sum2pho.pt(),w);



       //is sumPt2 the good choice?
       if( fabs(PV_z->at(goodIndex[0]) - TrueVertex_Z ) < 0.6 )
	 {
	   PtGood.Fill( sum2pho.pt(),w );
	   EtaGood.Fill( sum2pho.eta(),w );
	   NvtGood.Fill(npv,w);
	   PtGood_vsNvtGood.Fill(npv,sum2pho.pt(),w);
	 }
       
       
       


       // ------- TMVA ---------------

       
       int TMVAind = -1;
       float TMVAmax = -1000.;
           
       //pt balance for vertexId
       for ( int uu = 0; uu < npv; uu++)
	 {
	   int ngoodTracks = 0;
	   ROOT::Math::XYZVector sumTracks;
	   float summodpt = 0;
	   float sumpt2 = 0;
	   //  ROOT::Math::XYZVector sumTracksWrong;
	   
	   for (unsigned int kk = 0; kk < PVtracks->size(); ++kk)
	     {
	       if (PVtracks_PVindex->at(kk) == goodIndex[uu])
		 {
		   if ( PVtracks->at(kk).perp2() < trackThr*trackThr ) continue;
		   sumTracks += PVtracks->at(kk);
		   summodpt += sqrt((PVtracks->at(kk)).perp2());
		   sumpt2 += PVtracks->at(kk).perp2();
		   ngoodTracks ++;
		 }
	     }//tracks loop
	   
	   
	   //TMVA variables
	   sumPt2_TMVA = sumpt2;
	   tracksNum_TMVA = ngoodTracks;
	   deltaPhi_HSumPt_TMVA = deltaPhi( sumTracks.phi(),sum2pho.phi() );
	   sum2PhoPt_TMVA = sum2pho.pt();
	   sumPtModOversum2PhoPt_TMVA = sqrt( sumTracks.perp2() ) / sum2pho.pt();
	   
	   //Evaluate TMVA
	   Double_t mva = TMVAreader->EvaluateMVA( "BDTmethod" ); 
	   
	   if( uu == iClosest ) 
	     {
	       bdth.Fill(mva,w);
	       pt2h.Fill(sumpt2,w); 
	     }
	   else
	     {
	       bdtbkg.Fill(mva,w);
	       pt2bkg.Fill(sumpt2,w);
	     }
	   
	   //solo il meglio
	   if (mva > TMVAmax) 
	     {
	       TMVAmax = mva;
	       TMVAind = uu;
	     }
	   
	   
	 }//vertex loop
       
       if( fabs(PV_z->at(goodIndex[TMVAind]) - TrueVertex_Z ) < 0.6 )
	 {
	   PtGood_BDT.Fill( sum2pho.pt(), w );
	   EtaGood_BDT.Fill( sum2pho.eta(), w );
	   NvtGood_BDT.Fill( npv, w);
	   PtBDT_vsNvtBDT.Fill( npv,sum2pho.pt(), w );
	 }
       
       
     } //evt loop
   
   
   TFile ff( (outputRootFilePath+outputRootFileName).c_str(),"recreate");
   
   PtAll.Write();
   PtGood.Write();
   NvtAll.Write();
   NvtGood.Write();
   EtaAll.Write();
   EtaGood.Write();

   PtAll_vsNvtAll.Write();
   PtGood_vsNvtGood.Write();
   
   PtGood_BDT.Write();
   EtaGood_BDT.Write();
   NvtGood_BDT.Write();
   PtBDT_vsNvtBDT.Write();
   
   bdtbkg.Write();
   bdth.Write();
   
   pt2bkg.Write();
   pt2h.Write();
   
   nmatched.Write();

   sum2Pho_overMC.Write();


   for (int i = 0; i < 10 ; i ++){
     hHggSumPt2_bin[i]->Write();
     hZeeSumPt2_bin[i]->Write();

   }


   ff.Close();
   return 0;
   
   
}

bool PhotonId( float et, float eta, float Eiso, float Hiso, float HoE, float Tiso, float setaeta ) 
{

  bool iso = Eiso < (4.2 + 0.006*et ) && Hiso < (2.2 + 0.0025*et) && Tiso < (3.5 + 0.001*et);
  bool shape  = ( fabs (eta) < 1.45 && setaeta < 0.013 ) || ( fabs (eta) > 1.47 && setaeta < 0.030 );

  return (iso && shape && HoE<0.05);

}

