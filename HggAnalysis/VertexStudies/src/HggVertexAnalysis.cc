///==== include ====
#include "HggVertexAnalysis.h"

#define etaEB   1.4442 
#define etaEE   1.566


using namespace std;


//------------------------------------------------------------------------------------------------------------------------
HggVertexAnalysis::HggVertexAnalysis(h2gglobeEventReader *e)
{
  if (e == 0) {
    cout << "No input tree!!!" << endl;
  }
  ev_ = e;
}
//------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------
HggVertexAnalysis::~HggVertexAnalysis()
{
}
//------------------------------------------------------------------------------------------------------------------------


//--- RECO-MC PHOTONS MATCHING --------------------------------------------------------------------------------------------------------
void HggVertexAnalysis::findMCHiggsPhotons(h2gglobeEventReader *ev_, int& passSelection, int& mc1, int& mc2, int& i1, int& i2  )
{

  TClonesArray* gp_p4 = ev_->gp_p4;
  Short_t* gp_status  = ev_->gp_status;
  Short_t* gp_pdgid   = ev_->gp_pdgid;
  Short_t* gp_mother  = ev_->gp_mother;
  int gp_n            = ev_->gp_n;
  
  TClonesArray* pho_p4      = ev_->pho_p4;
 
  
  for ( int i = 0; i< gp_n ; i++ ){
    int pid        = gp_pdgid[i];
    int status     = gp_status[i];

    if (pid != 22 || status!=1) continue;
    
    // uncomment this is mother info available
    //     int mompid  = gp_pdgid[gp_mother[i]];
    //     int gmompid = gp_pdgid[gp_mother[gp_mother[i]]];
    //     //std::cout << i << "  " <<pid << "  " << mompid <<"  "<< gmompid<< std::endl;     
    //     if ( (mompid ==25 || gmompid ==25) && mc1 < 0 ) { mc1  = i; }
    //     else if ( (mompid ==25 || gmompid ==25) && mc2 < 0 ) { mc2  = i; }
    
    if ( mc1 < 0 ) { mc1  = i; }
    else if ( mc2 < 0 ) { mc2  = i; }
    
  }
  
  if (mc1 < 0 || mc2 < 0 ) return;
  
  TLorentzVector *mcV1 = (TLorentzVector*)gp_p4->At(mc1);  
  TLorentzVector *mcV2 = (TLorentzVector*)gp_p4->At(mc2);  

  int index1 = -100, index2 = -100;
  double dr1min = 10000.;
  double dr2min = 10000.;

  for (int i = 0; i < ev_->pho_n; i++){

    TVector3 *scpos = (TVector3*) ev_->pho_calopos -> At(i);
    
    double dr1 = scpos->DeltaR(mcV1->Vect());
    if (dr1 < dr1min) { dr1min = dr1; index1 = i; }

    double dr2 = scpos->DeltaR(mcV2->Vect());
    if (dr2 < dr2min) { dr2min = dr2; index2 = i; }
  }
 
  TLorentzVector *photon1 = (TLorentzVector*)pho_p4->At(index1);
  TLorentzVector *photon2 = (TLorentzVector*)pho_p4->At(index2);

  // order photons by pt
  if (photon1->Pt() > photon2->Pt()) {
    i1 = index1;
    i2 = index2;
  }
  else{
    i1 = index2;
    i2 = index1;
  }
  
  bool is_mcmatched   = (dr1min < 0.15) && (dr2min< 0.15);
  
  if( is_mcmatched ) {
    passSelection = 1;
    return;
  }
}



//---DI-PHOTON CATEGORIES --------------------------------------------------
//-- just for test
int HggVertexAnalysis::DiphotonCategory (float eta1, float eta2, float r9_1, float r9_2){
  
  int category = -1;
  
  if ( fabs(eta1) < etaEB && fabs(eta2) < etaEB && min(r9_1, r9_2) > 0.94 )  category = 0;
  if ( fabs(eta1) < etaEB && fabs(eta2) < etaEB && min(r9_1, r9_2) < 0.94 )  category = 1;
  
  if ( ( fabs(eta1) > etaEB || fabs(eta2) > etaEB) && min(r9_1, r9_2) > 0.94 )  category = 2;
  if ( ( fabs(eta1) > etaEB || fabs(eta2) > etaEB) && min(r9_1, r9_2) < 0.94 )  category = 3;
  return category ;
  
}

//---WEIGHT for KFACTORS----------------------------------------------------------------------
double HggVertexAnalysis::KfactorsWeight(TH1F* hkfact, TClonesArray* gp_p4, int mc1, int mc2){
  
  TLorentzVector *mcV1 = (TLorentzVector*)gp_p4->At(mc1);  
  TLorentzVector *mcV2 = (TLorentzVector*)gp_p4->At(mc2);  

  TLorentzVector mcH  = (*mcV1+*mcV2);

  float genpt = mcH.Pt();
  float bin   = hkfact->FindBin(genpt);
  float w     = hkfact-> GetBinContent(bin);
  return (w);

}



//---BOOK HISTOGRAMS-------------------------------------------------
void HggVertexAnalysis::bookHistos()
{
  //histos
  PtGen = new TH1F("PtGen","PtGen",80,0,400);

  PtAll_Baseline   = new TH1F("PtAll_Baseline","PtAll_Baseline",80,0,400);
  PtGood_Baseline  = new TH1F("PtGood_Baseline","PtGood_Baseline",80,0,400);
   
  EtaAll_Baseline  = new TH1F("EtaAll_Baseline","EtaAll_Baseline",50,-5,5);
  EtaGood_Baseline = new TH1F("EtaGood_Baseline","EtaGood_Baseline",50,-5,5);
    
  NvtAll_Baseline  = new TH1F("NvtAll_Baseline","number of PV all",50,0,50);
  NvtGood_Baseline = new TH1F("NvtGood_Baseline","number of PV good",50,0,50);
   
  NpuAll_Baseline  = new TH1F("NpuAll_Baseline","number of PU all",50,0,50);
  NpuGood_Baseline = new TH1F("NpuGood_Baseline","number of PU good",50,0,50);
  
  InvMassAll_Baseline  = new TH1F("InvMassAll_Baseline","Invariant mass all",4000,0,200);
  InvMassGood_Baseline = new TH1F("InvMassGood_Baseline","Invariant mass good",4000,0,200);

  PtAll_BDT   = new TH1F("PtAll_BDT","PtAll_BDT",80,0,400);
  PtGood_BDT  = new TH1F("PtGood_BDT","PtGood_BDT",80,0,400);
   
  EtaAll_BDT  = new TH1F("EtaAll_BDT","EtaAll_BDT",50,-5,5);
  EtaGood_BDT = new TH1F("EtaGood_BDT","EtaGood_BDT",50,-5,5);
    
  NvtAll_BDT  = new TH1F("NvtAll_BDT","number of PV all",50,0,50);
  NvtGood_BDT = new TH1F("NvtGood_BDT","number of PV good",50,0,50);
   
  NpuAll_BDT  = new TH1F("NpuAll_BDT","number of PU all",50,0,50);
  NpuGood_BDT = new TH1F("NpuGood_BDT","number of PU good",50,0,50);
  
  InvMassAll_BDT  = new TH1F("InvMassAll_BDT","Invariant mass all",4000,0,200);
  InvMassGood_BDT = new TH1F("InvMassGood_BDT","Invariant mass good",4000,0,200);

  char hname[100];
  char htitle[100];
  for (int i = 0; i < 4; i++){
    if (i==0)  sprintf(htitle,"#gamma#gamma Invariant mass, both photons in EB, min(R9)>0.94");
    if (i==2)  sprintf(htitle,"#gamma#gamma Invariant mass, both photons in EB, min(R9)<0.94");
    if (i==2)  sprintf(htitle,"#gamma#gamma Invariant mass, at least one photon in EE , min(R9)>0.94");
    if (i==3)  sprintf(htitle,"#gamma#gamma Invariant mass, at least one photon in EE , min(R9)<0.94");

    sprintf(hname,"InvMassAll_Baseline_cat%02d",i);
    InvMassAll_Baseline_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassGood_Baseline_cat%02d",i);
    InvMassGood_Baseline_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassAll_BDT_cat%02d",i);
    InvMassAll_BDT_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassGood_BDT_cat%02d",i);
    InvMassGood_BDT_cat[i] = new TH1F(hname,htitle,400,0,200); 

  }

  hsumpt2_sig = new TH1F("hsumpt2_sig","hsumpt2_sig",50,-5,15);
  hsumpt2_sig->SetLineColor(kBlue);
  hsumpt2_bkg = new TH1F("hsumpt2_bkg","hsumpt2_bkg",50,-5,15);
  hsumpt2_bkg->SetLineColor(kRed);

  hptasym_sig = new TH1F("hptasym_sig","hptasym_sig",50,-1,1);
  hptasym_sig->SetLineColor(kBlue);
  hptasym_bkg = new TH1F("hptasym_bkg","hptasym_bkg",50,-1,1);
  hptasym_bkg->SetLineColor(kRed);

  hptbal_sig  = new TH1F("hptbal_sig","hptbal_sig",100,-100,150);
  hptbal_sig->SetLineColor(kBlue);
  hptbal_bkg  = new TH1F("hptbal_bkg","hptbal_bkg",100,-100,150);
  hptbal_bkg->SetLineColor(kRed);
}

//---ANALYZE---------------------------------------------------------------------------------------------------------------------------

void HggVertexAnalysis::analyze(int nentries, int isData, int useWeights, TH1F* hpu, int useKfactors, TH1F* hkfact, int doVBFselection )
{

  if ( ev_== 0) return;

  cout << "Start analyzing ..." << endl;
  
  //--- preparing weights for PU
  TRandom *gRandom = new TRandom();
  float nmax = 0; 
  float w[50];
  if (!isData && useWeights){
    nmax = hpu -> GetMaximum();
    std::cout << "Maximum weight = " << nmax << std::endl;
    for (int ibin = 1; ibin < hpu->GetNbinsX()+1; ibin++){
      w[ibin-1] = hpu->GetBinContent(ibin);  // bin 1 --> npu = 0 
    }
  }

  //--- preparing weights for k-factors (glu-glu)
  float wkfact[500];
  float nmaxk = 0;
  if (!isData && useKfactors){
    nmaxk = hkfact -> GetMaximum();
    std::cout << "Maximum weight from k-factors = " << nmaxk << std::endl;
    for (int ibin = 1; ibin < hkfact->GetNbinsX()+1; ibin++){
      wkfact[ibin-1] = hkfact->GetBinContent(ibin);  // bin 1 --> npu = 0 
    }
  }
  
  //**** start loop over entries
  float myrnd;
  float ww = 1.;

  for (int jentry = 0; jentry < nentries ; jentry++ ){
    if(jentry%1000 == 0) std::cout<<"reading event "<< jentry << "\r" <<std::flush;
  
    ev_->GetEntry(jentry);
    
    //--- selections for Hgg: find reco photons matched to MC
    int accept = 0;
    int indpho1 = -100;
    int indpho2 = -100;
    int mc1     = -1;
    int mc2     = -1;
    findMCHiggsPhotons(ev_, accept, mc1, mc2, indpho1, indpho2);
    
    if (!accept) continue;
    if (mc1 < 0 || mc2 < 0) {
      cout << "No MC matching found for photons !!!! " << endl;
      continue;  // mc matching not found
    }    

    double TrueVertex_Z = ((TVector3*)ev_->gv_pos->At(0)) -> Z();      
    
    //     TLorentzVector *sc1 = (TLorentzVector*) ev_->sc_p4->At(ev_->pho_scind[indpho1]);
    //     TLorentzVector *sc2 = (TLorentzVector*) ev_->sc_p4->At(ev_->pho_scind[indpho2]);
    
    TLorentzVector *sc1 = (TLorentzVector*) ev_->pho_p4->At(indpho1);
    TLorentzVector *sc2 = (TLorentzVector*) ev_->pho_p4->At(indpho2);
       
    TVector3 *pho1_calopos = (TVector3*) ev_->pho_calopos->At(indpho1);
    TVector3 *pho2_calopos = (TVector3*) ev_->pho_calopos->At(indpho2);
    
    float eta1 = pho1_calopos -> Eta();
    float eta2 = pho2_calopos -> Eta();

    float pho1_r9 = ev_->pho_r9[indpho1];
    float pho2_r9 = ev_->pho_r9[indpho2];
    int  cat      = DiphotonCategory( eta1, eta2, pho1_r9, pho2_r9);
    
    if ( (fabs(eta1) > etaEB && fabs(eta1) < etaEE) || fabs(eta1) > 2.5) continue;
    if ( (fabs(eta2) > etaEB && fabs(eta2) < etaEE) || fabs(eta2) > 2.5) continue;
    

    //--- photons
    PhotonInfo pho1(indpho1, *pho1_calopos,  sc1->E()); 
    PhotonInfo pho2(indpho2, *pho2_calopos,  sc2->E()); 
    
    int nvtx = ev_->vtx_std_n;
    
    //--- event weights
    int npu = ev_->pu_n ;
    if (useWeights){
      if (useKfactors) 
	ww = w[npu]*KfactorsWeight(hkfact, ev_->gp_p4, mc1, mc2);
      else 
	ww = w[npu];
      //myrnd = gRandom->Uniform(0,nmax*nmaxk);
      //if (myrnd > ww) continue;
    }

          
    
    //--- gen Pt(gamma-gamma spectrum)
    TLorentzVector *mcV1 = (TLorentzVector*)ev_->gp_p4->At(mc1);  
    TLorentzVector *mcV2 = (TLorentzVector*)ev_->gp_p4->At(mc2);  
    TLorentzVector mcH   = (*mcV1+*mcV2);
     

    //************** vtx baseline selection (sumpt2 criterion)
    // recompute di-photon pt wrt the selected vtx      
    TVector3 * vtxpos= (TVector3*) ev_->vtx_std_xyz->At(0);
    TLorentzVector p1, p2, dipho;
    p1    = pho1.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    p2    = pho2.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    dipho = (p1+p2);

    //---- di-jet tagging (for VBF category)
    bool isVBF;
    if (doVBFselection) {
      std::pair<int,int> highestPtJets(-1,-1);
      highestPtJets = Select2HighestPtJets(ev_, p1, p2, -999.,-999.  );
      isVBF = DijetTag(highestPtJets, dipho);
    }

    if ((doVBFselection && isVBF) || (!doVBFselection)){
      //if ( p1.Pt() > (dipho.M()/3) && p2.Pt() > (dipho.M()/4) ){ // new 
      if ( (p1.Pt()/dipho.M())>(55./120.) && p2.Pt() > 25 ){ // vbf cuts 
	
	PtGen->Fill(mcH.Pt(),ww);

	PtAll_Baseline->Fill( dipho.Pt(),ww );
	InvMassAll_Baseline->Fill( dipho.M(),ww );
	InvMassAll_Baseline_cat[cat]-> Fill( dipho.M(),ww );
	
	EtaAll_Baseline->Fill( p1.Eta(),ww );
	NvtAll_Baseline->Fill( nvtx,ww );
	
	if (!isData) NpuAll_Baseline->Fill(npu,ww);
	
	// matching within 1cm
	if ( fabs( TrueVertex_Z - vtxpos->Z() ) < 1.) {
	  PtGood_Baseline->Fill( dipho.Pt() , ww);
	  EtaGood_Baseline->Fill( p1.Eta(), ww );
	  NvtGood_Baseline->Fill( nvtx,ww );
	  InvMassGood_Baseline->Fill( dipho.M(), ww );
	  InvMassGood_Baseline_cat[cat]-> Fill( dipho.M() , ww);
	  if (!isData) NpuGood_Baseline->Fill(npu, ww);
	}
      }
    }
    //********
    

    //************** vtx combined selection 
    //recompute di-photon pt wrt the selected vtx      
    int selvtx = ev_ -> vtx_std_sel;
    vtxpos= (TVector3*) ev_->vtx_std_xyz->At(selvtx);
    
    p1    = pho1.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    p2    = pho2.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    dipho = (p1+p2);
    
    //---- di-jet tagging (for VBF category)
    if (doVBFselection) {
      std::pair<int,int> highestPtJets(-1,-1);
      highestPtJets = Select2HighestPtJets(ev_, p1, p2, -999.,-999.  );
      isVBF = DijetTag(highestPtJets, dipho);
    }

    if ((doVBFselection && isVBF) || (!doVBFselection)){
      //if ( p1.Pt() > (dipho.M()/3) && p2.Pt() > (dipho.M()/4)){ // new 
      if ( (p1.Pt()/dipho.M())>(55./120.) && p2.Pt() > 25 ){ // vbf cuts 
   	PtAll_BDT->Fill( dipho.Pt(), ww );
	InvMassAll_BDT->Fill( dipho.M(), ww );
	InvMassAll_BDT_cat[cat]-> Fill( dipho.M(), ww );
	EtaAll_BDT->Fill( p1.Eta(), ww );
	NvtAll_BDT->Fill( nvtx , ww);
	if (!isData) NpuAll_BDT->Fill(npu, ww);
	
	// matching within 1cm
	if ( fabs( TrueVertex_Z - vtxpos->Z() ) < 1.) {
	  PtGood_BDT->Fill( dipho.Pt(), ww );
	  EtaGood_BDT->Fill( p1.Eta(), ww );
	  NvtGood_BDT->Fill( nvtx, ww );
	  if (!isData) NpuGood_BDT->Fill(npu, ww);
	  InvMassGood_BDT->Fill( dipho.M(), ww );
	  InvMassGood_BDT_cat[cat]-> Fill( dipho.M(), ww );
	}
      }
    }


    // --- vtx id variables
    if (doVBFselection && isVBF){
      for (int idipho = 0; idipho < ev_->dipho_n; idipho++){
	if ( ev_->dipho_leadind[idipho]!= indpho1 || ev_->dipho_subleadind[idipho]!= indpho2) continue;
	for (int iv=0; iv<nvtx; iv++){
	  TVector3 * vtxpos= (TVector3*) ev_->vtx_std_xyz->At(iv);
	  std::vector<std::vector<float> > * vtx_std_sumpt2 =  ev_->vtx_std_sumpt2;
	  std::vector<std::vector<float> > * vtx_std_ptasym =  ev_->vtx_std_ptasym;
	  std::vector<std::vector<float> > * vtx_std_ptbal  =  ev_->vtx_std_ptbal;
	  float sumpt2 =  (*vtx_std_sumpt2)[idipho][iv];
	  float ptasym =  (*vtx_std_ptasym)[idipho][iv];
	  float ptbal  =  (*vtx_std_ptbal)[idipho][iv];
	  if ( fabs( TrueVertex_Z - vtxpos->Z() ) < 1.) {
	    hsumpt2_sig -> Fill(log(sumpt2),ww);
	    hptasym_sig -> Fill(ptasym,ww);
	    hptbal_sig  -> Fill(ptbal,ww);
	  }
	  else{
	    hsumpt2_bkg -> Fill(log(sumpt2),ww);
	    hptasym_bkg -> Fill(ptasym,ww);
	    hptbal_bkg  -> Fill(ptbal,ww);
	  }
	}
      }
    }



    
  }// end loop over entries
  
}

//---SAVING HISTOGRAMS---------------------------------------------------------------------------------------

void HggVertexAnalysis::saveHistos(TFile * fout)
{
  
  fout->cd();
  
  PtGen->Write();

  PtAll_Baseline->Write();
  PtGood_Baseline->Write();
  InvMassAll_Baseline->Write();
  InvMassGood_Baseline->Write();
  EtaAll_Baseline->Write();
  EtaGood_Baseline->Write();
  NvtAll_Baseline->Write();
  NvtGood_Baseline->Write();
  NpuAll_Baseline->Write();
  NpuGood_Baseline->Write();
  
  PtAll_BDT->Write();
  PtGood_BDT->Write();
  InvMassAll_BDT->Write();
  InvMassGood_BDT->Write();
  EtaAll_BDT->Write();
  EtaGood_BDT->Write();
  NvtAll_BDT->Write();
  NvtGood_BDT->Write();
  NpuAll_BDT->Write();
  NpuGood_BDT->Write();


  for (int i = 0; i < 4; i++ ){
    InvMassAll_Baseline_cat[i]-> Write();
    InvMassGood_Baseline_cat[i]-> Write();
    InvMassAll_BDT_cat[i]-> Write();
    InvMassGood_BDT_cat[i]-> Write();
  }
  
  hsumpt2_sig -> Write();
  hsumpt2_bkg -> Write();
  hptasym_sig -> Write();
  hptasym_bkg -> Write();
  hptbal_sig -> Write();
  hptbal_bkg -> Write();

  fout->Close();
  
  return;
}

  
//--- SELECT 2 HIGHEST PT JETS (taken from h2gglobe) --------------------------------------------------------------------------------------

std::pair<int, int> HggVertexAnalysis::Select2HighestPtJets(h2gglobeEventReader* ev_, TLorentzVector& leadpho, 
							    TLorentzVector& subleadpho, float jtLMinPt, float jtTMinPt)
{
  std::pair<int, int> myJets(-1,-1);
  std::pair<int, int> myJetsnew(-1,-1);
  std::pair<float, float> myJetspt(-1.,-1.);

  float dr2pho = 0.5;
  float dr2jet = 0.5;
  
  TLorentzVector* j1p4;
  TLorentzVector* j2p4;
  float j1pt=-1;
  float j2pt=-1;

  // select highest pt jets
  // veto jets close to photons or each other
  for(int j1_i=0; j1_i<ev_->jet_algoPF1_n; j1_i++){

    j1p4 = (TLorentzVector*) ev_->jet_algoPF1_p4->At(j1_i);

    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();
    if(j1pt<jtTMinPt) continue;
    for(int j2_i=j1_i+1; j2_i<ev_->jet_algoPF1_n; j2_i++){
      j2p4 = (TLorentzVector*) ev_->jet_algoPF1_p4->At(j2_i);
      if(fabs(j2p4->Eta()) > 4.7) continue;
      if(j2p4->DeltaR(leadpho) < dr2pho) continue;
      if(j2p4->DeltaR(subleadpho) < dr2pho) continue;
      if(j2p4->DeltaR(*j1p4) < dr2jet) continue;
      j2pt=j2p4->Pt();
      
      if(j2pt<jtTMinPt) continue;
      if(std::max(j1pt,j2pt)<jtLMinPt) continue;

      if(j1pt>j2pt){
        jtLMinPt=j1pt;
        jtTMinPt=j2pt;
        myJets.first = j1_i;
        myJets.second = j2_i;
      } 
      else {
        jtLMinPt=j2pt;
        jtTMinPt=j1pt;
        myJets.first = j2_i;
        myJets.second = j1_i;
      }
    }
  }

  for(int j1_i=0; j1_i<ev_->jet_algoPF1_n; j1_i++){
    j1p4 = (TLorentzVector*) ev_->jet_algoPF1_p4->At(j1_i);
    if(fabs(j1p4->Eta()) > 4.7) continue;
    if(j1p4->DeltaR(leadpho) < dr2pho) continue;
    if(j1p4->DeltaR(subleadpho) < dr2pho) continue;
    j1pt=j1p4->Pt();

    if(j1pt>myJetspt.first) {
      myJetsnew.second=myJetsnew.first;
      myJetspt.second=myJetspt.first;
      myJetspt.first=j1pt;
      myJetsnew.first=j1_i;
    }
    else if(j1pt>myJetspt.second) {
      myJetspt.second=j1pt;
      myJetsnew.second=j1_i;
    }
  }

  
  return myJetsnew;

}
//------------------------------------------------------------------------------------------------------------------------------


// --- DIJET TAGGING FOR VBF ---------------------------------------------------------------------------------------------------
bool HggVertexAnalysis::DijetTag(std::pair<int,int> highestPtJets, TLorentzVector &dipho){
  
  float myVBFLeadJPt;
  float myVBFSubJPt;
  float myVBF_Mjj;
  float myVBFdEta;
  float myVBFZep;
  float myVBFdPhi;
  bool  isDijetTagged = false;

  bool VBFpresel = (highestPtJets.first!=-1)&&(highestPtJets.second!=-1);
  
  if(VBFpresel){
    
    TLorentzVector* jet1 = (TLorentzVector*)ev_->jet_algoPF1_p4->At(highestPtJets.first);
    TLorentzVector* jet2 = (TLorentzVector*)ev_->jet_algoPF1_p4->At(highestPtJets.second);
    TLorentzVector dijet = (*jet1) + (*jet2);
    myVBFLeadJPt = jet1->Pt();
    myVBFSubJPt  = jet2->Pt();
    myVBF_Mjj    = dijet.M();
    myVBFdEta    = fabs(jet1->Eta() - jet2->Eta());
    myVBFZep     = fabs(dipho.Eta() - 0.5*(jet1->Eta() + jet2->Eta()));
    myVBFdPhi    = fabs(dipho.DeltaPhi(dijet));
    
    if (
	myVBFLeadJPt > 30  && 
	myVBFSubJPt  > 20  && 
	myVBFdEta    > 3.5 &&
	myVBFZep     < 2.5 && 
	myVBFdPhi    > 2.6 &&
	myVBF_Mjj    > 350.
	)                       
      isDijetTagged = true;
  }
  
  return isDijetTagged;
}
//---------------------------------------------------------------------------------------------------------------------------
