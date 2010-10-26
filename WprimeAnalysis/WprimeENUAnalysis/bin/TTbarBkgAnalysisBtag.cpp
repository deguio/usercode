#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//---- tree content ----
#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"


//---- CMSSW includes
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEleSelector.h"


//---- root includes
#include "TH1.h"
#include "TH2.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>


#define PI 3.14159265
#define TWOPI 6.28318530

using namespace std;

double deltaPhi(double phi1,double phi2) {
 double deltaphi = fabs(phi1-phi2);  
 if (deltaphi > TWOPI) deltaphi -= TWOPI;  
 if (deltaphi > PI ) deltaphi = TWOPI - deltaphi;  
 return deltaphi; 
}

//! main program
int main (int argc, char** argv)
{

  // Reading input parameters
  
  //---- Input Variables ----
  /*std::string FileIn;
  std::string FileOut;
  float Xsec;
  float GenEff;
  float Lumi;
  float IsData;
    
  std::string inputName = argv[1];
  char buffer[1000];
  ifstream file(inputName.c_str());
  
  std::string variableNameFileIn  = "FileIn";
  std::string variableNameFileOut = "FileOut";
  std::string variableNameXsec    = "Xsec";
  std::string variableNameGenEff  = "GenEff";
  std::string variableNameLumi    = "Lumi";
  std::string variableNameIsData  = "IsData";
  std::cerr << " Reading  " << inputName << " ... " << std::endl;

  file >>  variableNameFileIn  >> FileIn 
       >>  variableNameFileOut >> FileOut 
       >>  variableNameXsec    >> Xsec
       >>  variableNameGenEff  >> GenEff
       >>  variableNameLumi    >> Lumi
       >>  variableNameIsData  >> IsData ;

  cout <<  variableNameFileIn  << "  "  << FileIn  << "  "  
       <<  variableNameFileOut << "  "  << FileOut << "  "  
       <<  variableNameXsec    << "  "  << Xsec    << "  "  
       <<  variableNameGenEff  << "  "  << GenEff  << "  "  
       <<  variableNameLumi    << "  "  << Lumi    << "  "  
       <<  variableNameIsData  << "  "  << IsData  << "  "  
       <<  endl;
  */

  // read input parameters
  std::string InFileName  = argv[1];
  std::string OutFileName = argv[2] ;
  float  Xsec    = atof(argv[3]);
  float  IsData  = atof(argv[4]);

  std::string FileIn  ;
  if (IsData) FileIn = "/media/amassiro/Wprime/DATA_20102010/" + InFileName;
  else FileIn = "/media/amassiro/Wprime/MC_20102010/" + InFileName;

  std::string FileOut =  "outputHistos_" + OutFileName + "_7TeV.root";;

  cout << FileIn << endl;
  cout << FileOut << endl;



  //================================
  // === load the analysis tree ====
  //================================
  TChain *chain = new TChain ("myanalysis/WprimeAnalysisTree") ;
  WprimeTreeContent treeVars ;
  setBranchAddresses (chain, treeVars) ;
  chain->Add(FileIn.c_str());
  int nEntries = chain->GetEntries();
  std::cout << "------> Number of events "<< nEntries << " <------\n" ;
  

  //===========================
  // === load the HLT tree ====
  //===========================
  int HLTnpaths;
  int HLTwasrun[20];
  int HLTaccept[20];
  int HLTerror[20];

  TChain *chainHLT = new TChain ("TriggerResults/HLTree") ;
  chainHLT -> SetBranchAddress("HLTnpaths", &HLTnpaths);
  chainHLT -> SetBranchAddress("HLTwasrun", HLTwasrun);
  chainHLT -> SetBranchAddress("HLTaccept", HLTaccept);
  chainHLT -> SetBranchAddress("HLTerror" , HLTerror);
  chainHLT -> Add(FileIn.c_str());


  //================================
  // calculate event weight for MC =
  //================================
  float w;
  float Lumi = 1.; // norm all MC to 1 pb^-1
  if (IsData) 
    w = 1;
  else{
    TFile* fin = new TFile(FileIn.c_str(), "READ");
    TH1F *hevents = (TH1F*)fin ->Get("eventsCounterTotal/eventsCounterTotal");
    w = Xsec*Lumi/(float)hevents->GetBinContent(1);
    cout <<"xsec = " << Xsec   << "     w =" << w << endl;
  }
  


  // file to save histos
  
  TFile *fout = new TFile(FileOut.c_str(),"recreate");

  //INITIALIZING HISTOGRAMS

  TH1F *het       = new TH1F("het","het",1500,0,1500);
  TH1F *hmet      = new TH1F("hmet","hmet",1500,0,1500);
  TH1F *hmt       = new TH1F("hmt","hmt",3000,0,3000);
  TH1F *hetOmet   = new TH1F("hetOmet","hetOmet",150,0,15);
  TH1F *hdphi     = new TH1F("hdphi","hdphi",320,0,3.2);

  TH1F *hetEB     = new TH1F("hetEB","hetEB",1500,0,1500);
  TH1F *hmetEB    = new TH1F("hmetEB","hmetEB",1500,0,1500);
  TH1F *hmtEB     = new TH1F("hmtEB","hmtEB",3000,0,3000);
  TH1F *hetOmetEB = new TH1F("hetOmetEB","hetOmetEB",150,0,15);
  TH1F *hdphiEB   = new TH1F("hdphiEB","hdphiEB",320,0,3.2);

  TH1F *hetEE     = new TH1F("hetEE","hetEE",1500,0,1500);
  TH1F *hmetEE    = new TH1F("hmetEE","hmetEE",1500,0,1500);
  TH1F *hmtEE     = new TH1F("hmtEE","hmtEE",3000,0,3000);
  TH1F *hetOmetEE = new TH1F("hetOmetEE","hetOmetEE",150,0,15);
  TH1F *hdphiEE   = new TH1F("hdphiEE","hdphiEE",320,0,3.2);

  TH1F *het1b = new TH1F("het1b","het1b",1500,0,1500);
  TH1F *het2b = new TH1F("het2b","het2b",1500,0,1500);

  TH1F *hmet1b = new TH1F("hmet1b","hmet1b",1500,0,1500);
  TH1F *hmet2b = new TH1F("hmet2b","hmet2b",1500,0,1500);

  TH1F *hmt1b = new TH1F("hmt1b","hmt1b",3000,0,3000);
  TH1F *hmt2b = new TH1F("hmt2b","hmt2b",3000,0,3000);

  TH1F *hNjets = new TH1F("hNjets","hNjets",50,0,50);
  TH1F *hNjets1b = new TH1F("hNjets1b","hNjets1b",50,0,50);
  TH1F *hNjets2b = new TH1F("hNjets2b","hNjets2b",50,0,50);

  TH1F *hBdisc = new TH1F("hBdisc","hBdisc",100,0,10);
  TH2F *hBdisc_vs_mt = new TH2F("hBdisc_vs_mt","hBdisc_vs_mt",500,0,500,100,0,10);

  TH1F *het1bNonIso = new TH1F("het1bNonIso","het1bNonIso",1500,0,1500);
  TH1F *het2bNonIso = new TH1F("het2bNonIso","het2bNonIso",1500,0,1500);

  TH1F *hmet1bNonIso = new TH1F("hmet1bNonIso","hmet1bNonIso",1500,0,1500);
  TH1F *hmet2bNonIso = new TH1F("hmet2bNonIso","hmet2bNonIso",1500,0,1500);

  TH1F *hmt1bNonIso = new TH1F("hmt1bNonIso","hmt1bNonIso",3000,0,3000);
  TH1F *hmt2bNonIso = new TH1F("hmt2bNonIso","hmt2bNonIso",3000,0,3000);
 

  float EtaCutEB    = 1.4442;
  float EtaCutEE    = 1.560;

  float n1b = 0;
  float n2b = 0;
  float bdiscThr = 3.0;

  int nGoodElectrons = 0;
  int chosenEle = 0;

  float counter = 0;

  // START LOOP OVER ENTRIES
  for (int entry = 0 ; entry < nEntries ; ++entry) {

    chain->GetEntry(entry) ;
    if(entry%100000==0) std::cout << "------> reading entry " << entry << " <------\n" ;

    //if ( (!IsData) == true  &&  fabs(treeVars.pdgId[0]) != 11 ) continue; 


    //================================
    //==== HLT selection          ====
    //================================
    //FIXME hardcoded
    //'0 -> HLT_Photon10_L1R',
    //'1 -> HLT_Ele10_LW_L1R',
    //'2 -> HLT_Ele15_LW_L1R',
    //'3 -> HLT_Ele15_SW_L1R',
    //'4 -> HLT_Ele15_SW_CaloEleId_L1R'
    //'5 -> HLT_Ele17_SW_CaloEleId_L1R'
    //'6 -> HLT_Ele22_SW_CaloEleId_L1R'
    //'7 -> HLT_Ele27_SW_CaloEleIdTrack_L1R_v1'
    
    chainHLT->GetEntry(entry);

    int accept = 0;

       

    if (IsData){
      if ( treeVars.runId >= 135059 &&  treeVars.runId <= 140041  ) accept = HLTaccept[1];
      if ( treeVars.runId >= 140042 &&  treeVars.runId <= 141900  ) accept = HLTaccept[3];
      if ( treeVars.runId >= 141901 &&  treeVars.runId <= 146427  ) accept = HLTaccept[4];
      if ( treeVars.runId >= 146428 &&  treeVars.runId <= 147116  ) accept = HLTaccept[5];  
      if ( treeVars.runId >= 147117 &&  treeVars.runId <= 999999  ) accept = HLTaccept[7];   
    }

    if (!IsData){
      //accept = HLTaccept[2] + HLTaccept[3] ; // cmssw36X
      for (int ii = 0; ii < 8; ii++){
	accept+=HLTaccept[ii];
      }
    }

    if (accept == 0) continue;
    


    double met    = treeVars.tcMet;
    double metPhi = treeVars.tcMetPhi;
    double mex    = treeVars.tcMex;
    double mey    = treeVars.tcMey;

    nGoodElectrons = 0;
    chosenEle = 0;

    // loop over electron candidates
    for (int i = 0; i < treeVars.nElectrons; ++i)
      {

	double eleEt  = treeVars.eleEt[i];
	double elePx  = treeVars.elePx[i];
	double elePy  = treeVars.elePy[i];
	double elePhi = treeVars.elePhi[i];
	double eleEta = treeVars.eleEta[i];
	double elePt  = sqrt(elePx*elePx + elePy*elePy);      
	double eleId  = treeVars.eleId[i];
	double cphi   = (elePx*mex + elePy*mey ) / (met*elePt);
	double mt     = sqrt(2*eleEt*met*(1 - cphi));
	double dPhiEleMet = deltaPhi(metPhi,elePhi);
	bool eleIsEB = treeVars.eleIsEB[i];
	bool eleIsEE = treeVars.eleIsEE[i];
	



	//--------------------------------------------
	// QCD data driven
	//--------------------------------------------
	// check all selection but invert isolation

	// btagging
	int nbtagNonIso = 0;
	int njetsNonIso = 0;
	double deta, dphi, dR;
	    

	if ( eleEt> 30 && dPhiEleMet > 2.5 &&  eleEt/met > 0.4 && eleEt/met < 1.5  ) {

	  if ( (eleIsEB == 1 && 
		heep::CutCodes::passCuts(eleId,"ecalDriven:et:detEta:dEtaIn:dPhiIn:e2x5Over5x5:hadem:isolPtTrks") && 
	        !heep::CutCodes::passCuts(eleId,"isolEmHadDepth1")) ||
	       (eleIsEE == 1 && 
		heep::CutCodes::passCuts(eleId,"ecalDriven:et:detEta:dPhiIn:sigmaIEtaIEta:hadem:isolPtTrks") && 
		!(heep::CutCodes::passCuts(eleId,"isolEmHadDepth1") && 
		  heep::CutCodes::passCuts(eleId,"isolHadDepth2")))
	       )
	    
	    {
	      
	      for (int j=0; j<treeVars.nJets; j++){
		deta = treeVars.jetEta[j] -  eleEta;
		dphi = deltaPhi(treeVars.jetPhi[j],elePhi);
		dR = sqrt(deta*deta+dphi*dphi);
		if (dR<0.5) continue;
		if (fabs(treeVars.jetEta[j])<2.4 && treeVars.jetPt[j] > 20.) {
		  float bdisc = treeVars.jetBdiscHighPur[j];
		  if (bdisc > bdiscThr) nbtagNonIso++; // high purity btag
		  njetsNonIso++;
		}
	      }// end loop over jets
	      
	      if (nbtagNonIso==1) {
		het1bNonIso->Fill(eleEt,w);
		hmet1bNonIso->Fill(met,w);
		hmt1bNonIso->Fill(mt,w);
	      }
	      
	      if (nbtagNonIso==2) {
		het2bNonIso->Fill(eleEt,w);
		hmet2bNonIso->Fill(met,w);
		hmt2bNonIso->Fill(mt,w);
	      }
	      
	    }
	}
	               	    
     
	//--------------------------------------------
	// keep only electrons in ECAL fiducial volume
	//--------------------------------------------
	
	if ( treeVars.eleIsGap[i] == 1 ) continue;

	//--------------------------------------------
	// keep only electrons with ET > 30 GeV
	//--------------------------------------------
	if ( eleEt < 30 ) continue ;
	
	//--------------------------------------------
	// electron ID
	//--------------------------------------------
	//instructions here: 
	//http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h?revision=1.5&view=markup&pathrev=HEAD
	//https://twiki.cern.ch/twiki/bin/view/CMS/HEEPSelector

	//if ( treeVars.eleId[i] != eleId_ ) continue;
	  
	if ( treeVars.eleId[i] != 0x0 ) continue; 
	
	// barrel
	//if ( fabs(eta) < EtaCutEB && (treeVars.eleId[i] != 0x0) ) continue; 
	//endcap
	//if ( fabs(eta) > EtaCutEE && (treeVars.eleId[i] != 0x0) && ((treeVars.eleId[i] &~ 0x0010) != 0x0) ) continue;// remove eta cut in EE for 36X reco  
	  
	  nGoodElectrons++;
	  chosenEle = i;
	  
      }// end loop over ele cand
    
  

    if ( nGoodElectrons != 1 ) continue;
   
    double eleEt  = treeVars.eleEt[chosenEle];
    double elePx  = treeVars.elePx[chosenEle];
    double elePy  = treeVars.elePy[chosenEle];
    double elePhi = treeVars.elePhi[chosenEle];
    double eleEta = treeVars.eleEta[chosenEle];
    double elePt  = sqrt(elePx*elePx + elePy*elePy);      
    double cphi   = (elePx*mex + elePy*mey ) / (met*elePt);
    double mt     = sqrt(2*eleEt*met*(1 - cphi));
    double dPhiEleMet = deltaPhi(metPhi,elePhi);


    // EB only
    if (fabs(eleEta)<EtaCutEB){
      if (dPhiEleMet > 2.5 ) hetOmetEB->Fill(eleEt/met,w);
      if (eleEt/met > 0.4 && eleEt/met < 1.5) hdphiEB ->Fill(dPhiEleMet,w);
    }

    // EE only
    if (fabs(eleEta)>EtaCutEB){
      if (dPhiEleMet > 2.5 ) hetOmetEE->Fill(eleEt/met,w);
      if (eleEt/met > 0.4 && eleEt/met < 1.5) hdphiEE ->Fill(dPhiEleMet,w);
    }

    // EB+EE
    hetOmet -> Fill(eleEt/met,w);
    hdphi   -> Fill(dPhiEleMet,w);
    
    if (eleEt/met < 0.4 || eleEt/met > 1.5) continue;
    if (dPhiEleMet < 2.5 ) continue;
    
    het    -> Fill(eleEt,w);
    hmet   -> Fill(met,w);
    hmt    -> Fill(mt,w);

    // EB only
    if (fabs(eleEta)<EtaCutEB){
      hetEB    -> Fill(eleEt,w);
      hmtEB    -> Fill(mt,w);
      hmetEB   -> Fill(met,w);
    }
    // EE only
    if (fabs(eleEta)>EtaCutEB){
      hetEE    -> Fill(eleEt,w);
      hmtEE    -> Fill(mt,w);
      hmetEE   -> Fill(met,w);
    }

    if (mt > 100 ) counter+=w;

    // btagging
    int nbtag = 0;
    int njets = 0;
    double deta, dphi, dR;

    for (int j=0; j<treeVars.nJets; j++){
      deta = treeVars.jetEta[j] -  eleEta;
      dphi = deltaPhi(treeVars.jetPhi[j],elePhi);
      dR = sqrt(deta*deta+dphi*dphi);
      if (dR<0.5) continue;
      if (fabs(treeVars.jetEta[j])<2.4 && treeVars.jetPt[j] > 20.) {
	float bdisc = treeVars.jetBdiscHighPur[j];
	hBdisc->Fill(bdisc,w);
	hBdisc_vs_mt ->Fill(mt,bdisc,w);
	if (bdisc > bdiscThr) nbtag++; // high purity btag
	njets++;
      }
    }// end loop over jets

    hNjets -> Fill(njets,w);
 
    
    if (nbtag==1) {
      n1b+=w;
      het1b->Fill(eleEt,w);
      hmet1b->Fill(met,w);
      hmt1b->Fill(mt,w);
      hNjets1b->Fill(njets,w);
    }
    
    if (nbtag==2) {
      n2b+=w;
      het2b->Fill(eleEt,w);
      hmet2b->Fill(met,w);
      hmt2b->Fill(mt,w);
      hNjets2b->Fill(njets,w);
    }
    
  }// END LOOP OVER ENTRIES
 
  het->Sumw2();
  het1b->Sumw2();
  het2b->Sumw2();

  hmet->Sumw2();
  hmet1b->Sumw2();
  hmet2b->Sumw2();

  hmt->Sumw2();
  hmt1b->Sumw2();
  hmt2b->Sumw2();

  hNjets->Sumw2();
  hNjets1b->Sumw2();
  hNjets2b->Sumw2();

  hBdisc ->Sumw2();

  cout << "Numer of events mt > 100 : " << counter << endl; 

  // compute N(TT) and b-tagging efficiency

  std::cout << "Number of events with 1 b-jet = " << n1b << std::endl;
  std::cout << "Number of events with 2 b-jet = " << n2b << std::endl;

  float A1 = 0.118; // sunghyun
  float A2 = 0.848; // sunghyun
  float eA1 = 0.01;
  float eA2 = 0.01;

  float eb = (A1/A2+2)/(n1b/n2b+2);
  float eeb = (pow(eA1/A1,2) + pow(eA2/A2,2))/(pow(A1/A2+2,2)) + (1./n1b + 1./n2b)/(pow(n1b/n2b+2,2));
  eeb = eb * sqrt(eeb);
  cout<<"b-tagging efficiency = "<< eb <<" +/- "<<eeb<<endl;

  float nexp = n2b/(eb*eb*A2);
  float errn = nexp*sqrt(1./n2b + eA2/A2*eA2/A2 + pow(2*eeb/eb,2));
  cout<<"True # of ttbar events = "<<het->GetSumOfWeights()<<endl;
  cout<<"Estimated N(ttbar) from b-tag method = "<< nexp <<" +/- "<< errn <<endl;
  cout<<endl;



  // write histos 
  fout->Write();
  fout->Close();


}
