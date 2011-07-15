///==== include ====
#include "HggVertexAnalysis.h"

#define etaEB   1.4442 
#define etaEE   1.566


using namespace std;


//------------------------------------------------------------------------------------------------------------------------
HggVertexAnalysis::HggVertexAnalysis(h2gglobeEventReader *ev_)
{
  if (ev_ == 0) {
    cout << "No input tree!!!" << endl;
  }
}
//------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------
HggVertexAnalysis::~HggVertexAnalysis()
{
}
//------------------------------------------------------------------------------------------------------------------------


//-- RECO-MC PHOTONS MATCHING --------------------------------------------------------------------------------------------------------
void HggVertexAnalysis::findMCHiggsPhotons( TClonesArray* gp_p4, Short_t* gp_status, Short_t* gp_pdgid, Short_t* gp_mother, int gp_n, 
 					    TClonesArray* pho_calopos, TClonesArray* pho_p4, int pho_n,
					    int& passSelection, int& i1, int& i2  )
{
  int mc1 = -1;
  int mc2 = -1;
  
  for ( int i = 0; i< gp_n ; i++ ){
    int pid        = gp_pdgid[i];
    int status     = gp_status[i];
    
    if (pid != 22 || status!=1) continue;
    int mompid  = gp_pdgid[gp_mother[i]];
    int gmompid = gp_pdgid[gp_mother[gp_mother[i]]];
    
    if ( (mompid ==25 || gmompid ==25) && mc1 < 0 ) { mc1  = i; }
    else if ( (mompid ==25 || gmompid ==25) && mc2 < 0 ) { mc2  = i; }
  }
    
  if (mc1 < 0 || mc2 < 0 ) return;

  TLorentzVector *mcV1 = (TLorentzVector*)gp_p4->At(mc1);  
  TLorentzVector *mcV2 = (TLorentzVector*)gp_p4->At(mc2);  

  int index1 = -100, index2 = -100;
  double dr1min = 10000.;
  double dr2min = 10000.;

  for (int i = 0; i < pho_n; i++){

    TVector3 *scpos = (TVector3*) pho_calopos -> At(i);
    
    double dr1 = scpos->DeltaR(mcV1->Vect());
    if (dr1 < dr1min) { dr1min = dr1; index1 = i; }

    double dr2 = scpos->DeltaR(mcV2->Vect());
    if (dr2 < dr2min) { dr2min = dr2; index2 = i; }

  }
 
  TLorentzVector *photon1 = (TLorentzVector*)pho_p4->At(index1);
  TLorentzVector *photon2 = (TLorentzVector*)pho_p4->At(index2);
  
  if (photon1->E() > photon2->E()) {
    i1 = index1;
    i2 = index2;
  }
  else{
    i1 = index2;
    i2 = index1;
  }
    
  bool is_mcmatched   = (dr1min < 0.15) && (dr2min< 0.15);
  //bool is_unconverted = (photons_r9->at(index1) > 0.93) &&  (photons_r9->at(index2) > 0.93);
  bool pass_kincuts   = (photon1->Pt()> 40.) && ( photon2->Pt() > 30.) ;
 
  //if( is_mcmatched && is_unconverted && pass_kincuts ) {
  if( is_mcmatched && pass_kincuts ) {
    passSelection = 1;
    return;
  }
}



//---DI-PHOTON CATEGORIES --------------------------------------------------

int HggVertexAnalysis::DiphotonCategory (float eta1, float eta2, float r9_1, float r9_2){
  
  int category = -1;
  
  if ( fabs(eta1) < etaEB && fabs(eta2) < etaEB && min(r9_1, r9_2) > 0.94 )  category = 0;
  if ( fabs(eta1) < etaEB && fabs(eta2) < etaEB && min(r9_1, r9_2) < 0.94 )  category = 1;
  
  if ( ( fabs(eta1) > etaEB || fabs(eta2) > etaEB) && min(r9_1, r9_2) > 0.94 )  category = 2;
  if ( ( fabs(eta1) > etaEB || fabs(eta2) > etaEB) && min(r9_1, r9_2) < 0.94 )  category = 3;
  return category ;
  
}


//----------------------------------------------------

void HggVertexAnalysis::bookHistos()
{
  //histos
  PtAll_sumpt2  = new TH1F("PtAll_sumpt2","PtAll_sumpt2",80,0,400);
  PtGood_sumpt2 = new TH1F("PtGood_sumpt2","PtGood_sumpt2",80,0,400);
  PtAll_rank    = new TH1F("PtAll_rank","PtAll_rank",80,0,400);
  PtGood_rank   = new TH1F("PtGood_rank","PtGood_rank",80,0,400);
 
  EtaAll_sumpt2  = new TH1F("EtaAll_sumpt2","EtaAll_sumpt2",50,-5,5);
  EtaGood_sumpt2 = new TH1F("EtaGood_sumpt2","EtaGood_sumpt2",50,-5,5);
  EtaAll_rank    = new TH1F("EtaAll_rank","EtaAll_rank",50,-5,5);
  EtaGood_rank   = new TH1F("EtaGood_rank","EtaGood_rank",50,-5,5);
  
  NvtAll_sumpt2  = new TH1F("NvtAll_sumpt2","number of PV all",50,0,50);
  NvtGood_sumpt2 = new TH1F("NvtGood_sumpt2","number of PV good",50,0,50);
  NvtAll_rank    = new TH1F("NvtAll_rank","number of PV all",50,0,50);
  NvtGood_rank   = new TH1F("NvtGood_rank","number of PV good",50,0,50);
 
  NpuAll_sumpt2  = new TH1F("NpuAll_sumpt2","number of PU all",50,0,50);
  NpuGood_sumpt2 = new TH1F("NpuGood_sumpt2","number of PU good",50,0,50);
  NpuAll_rank    = new TH1F("NpuAll_rank","number of PU all",50,0,50);
  NpuGood_rank   = new TH1F("NpuGood_rank","number of PU good",50,0,50);

  InvMassAll_sumpt2  = new TH1F("InvMassAll_sumpt2","Invariant mass all",400,0,200);
  InvMassGood_sumpt2 = new TH1F("InvMassGood_sumpt2","Invariant mass good",400,0,200);
  InvMassAll_rank    = new TH1F("InvMassAll_rank","Invariant mass all",400,0,200);
  InvMassGood_rank   = new TH1F("InvMassGood_rank","Invariant mass good",400,0,200);

  char hname[100];
  char htitle[100];
  for (int i = 0; i < 4; i++){
    if (i==0)  sprintf(htitle,"#gamma#gamma Invariant mass, both photons in EB, min(R9)>0.94");
    if (i==2)  sprintf(htitle,"#gamma#gamma Invariant mass, both photons in EB, min(R9)<0.94");
    if (i==2)  sprintf(htitle,"#gamma#gamma Invariant mass, at least one photon in EE , min(R9)>0.94");
    if (i==3)  sprintf(htitle,"#gamma#gamma Invariant mass, at least one photon in EE , min(R9)<0.94");

    sprintf(hname,"InvMassAll_sumpt2_cat%02d",i);
    InvMassAll_sumpt2_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassGood_sumpt2_cat%02d",i);
    InvMassGood_sumpt2_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassAll_rank_cat%02d",i);
    InvMassAll_rank_cat[i] = new TH1F(hname,htitle,400,0,200); 

    sprintf(hname,"InvMassGood_rank_cat%02d",i);
    InvMassGood_rank_cat[i] = new TH1F(hname,htitle,400,0,200); 

  }

}

//---------------------------------------------------------------------------------------------------------------

void HggVertexAnalysis::analyze(int nentries, int isData, int useWeights, TH1F* h)
{

  if ( ev_==0) return;

  cout << "Start analyzing ..." << endl;
  
  TRandom *gRandom = new TRandom();
  int nmax = 0; 
  float w[50];
  if (!isData && useWeights){
    nmax = h ->GetMaximum();
    std::cout << "Maximum weight = " << nmax << std::endl;
    for (int ibin = 1; ibin < h->GetNbinsX()+1; ibin++){
      w[ibin-1] = h->GetBinContent(ibin);  // bin 1 --> npu = 0 
    }
  }
  
   cout << " entries in tree " << ev_->GetEntries() << endl;    

//   //start loop over entries
//   for (int jentry = 0; jentry < nentries ; jentry++ ){
//     if(jentry%1000 == 0) std::cout<<"reading event "<< jentry <<std::endl;
    
//     Long64_t ientry = ev_->LoadTree(jentry);
    //if (ientry < 0) break;
    //ev_->GetEntry(jentry);
    
    /*
    int npu = 0;
    //--- use weights 
    if (useWeights){
      npu = ev_->pu_n ;
      float myrnd = gRandom->Uniform(0,nmax);
      if (myrnd > w[npu]) continue;
    }

    
    // selections for Hgg
    if ( ev_ -> gv_n != 1) continue;
    
    TClonesArray* gp_p4 = ev_->gp_p4;
    Short_t* gp_status  = ev_->gp_status;
    Short_t* gp_pdgid   = ev_->gp_pdgid;
    Short_t* gp_mother  = ev_->gp_mother;
    int gp_n            = ev_->gp_n;
    
    TClonesArray* sc_p4       = ev_->sc_p4;
    TClonesArray* pho_p4      = ev_->pho_p4;
    TClonesArray* pho_calopos = ev_->pho_calopos;  
    
    int accept = 0;
    int indpho1 = -100;
    int indpho2 = -100;
    findMCHiggsPhotons(gp_p4, gp_status, gp_pdgid, gp_mother, gp_n, pho_calopos, pho_p4, ev_->pho_n, accept, indpho1, indpho2);
    if (!accept) continue;
    
    double TrueVertex_Z = ((TVector3*)ev_->gv_pos->At(0)) -> Z();      

    TLorentzVector *sc1 = (TLorentzVector*) ev_->sc_p4->At(ev_->pho_scind[indpho1]);
    TLorentzVector *sc2 = (TLorentzVector*) ev_->sc_p4->At(ev_->pho_scind[indpho2]);
    
    TVector3 *pho1_calopos = (TVector3*) pho_calopos->At(indpho1);
    TVector3 *pho2_calopos = (TVector3*) pho_calopos->At(indpho2);
    
    float eta1 = pho1_calopos -> Eta();
    float eta2 = pho2_calopos -> Eta();

    float pho1_r9 = ev_->pho_r9[indpho1];
    float pho2_r9 = ev_->pho_r9[indpho2];
    int  cat      = DiphotonCategory( eta1, eta2, pho1_r9, pho2_r9);
    
    if ( (fabs(eta1) > etaEB && fabs(eta1) < etaEE) || fabs(eta1) > 2.5) continue;
    if ( (fabs(eta2) > etaEB && fabs(eta2) < etaEE) || fabs(eta2) > 2.5) continue;
    
    PhotonInfo pho1( *pho1_calopos,  sc1->E()); 
    PhotonInfo pho2( *pho2_calopos,  sc2->E()); 
    
    int nvtx = ev_->vtx_std_n;

    // if ( runvtxid ) {
    
    //TupleVertexInfo vinfo( nvtx, vtxx_ , vtxy_, vtxz_, ntracks_, tkpx_, tkpy_, tkpz_, tkPtErr_, tkVtxId_, tkWeight_, tkd0_, tkd0Err_,tkdz_, tkdzErr_ , tkIsHighPurity_);
    
    //}

    //************** vtx baseline selection (sumpt2 criterion)
    // recompute di-photon pt wrt the selected vtx      
    TVector3 * vtxpos= (TVector3*) ev_->vtx_std_xyz->At(0);
    TLorentzVector p1, p2, dipho;
    p1    = pho1.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    p2    = pho2.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    dipho = (p1+p2);
    
    PtAll_sumpt2->Fill( dipho.Pt() );
    InvMassAll_sumpt2->Fill( dipho.M() );
    InvMassAll_sumpt2_cat[cat]-> Fill( dipho.M() );
        
    EtaAll_sumpt2->Fill( p1.Eta() );
    NvtAll_sumpt2->Fill( nvtx );
    if (!isData) NpuAll_sumpt2->Fill(npu);
    

    // matching within 1cm
    if ( fabs( TrueVertex_Z - vtxpos->Z() ) < 1.) {
      PtGood_sumpt2->Fill( dipho.Pt() );
      EtaGood_sumpt2->Fill( p1.Eta() );
      NvtGood_sumpt2->Fill( nvtx );
      InvMassGood_sumpt2->Fill( dipho.M() );
      InvMassGood_sumpt2_cat[cat]-> Fill( dipho.M() );
      if (!isData) NpuGood_sumpt2->Fill(npu);
    }
    //********
    

    //************** vtx combined selection 
    // recompute di-photon pt wrt the selected vtx      
    int selvtx = ev_ -> vtx_std_sel;
    vtxpos= (TVector3*) ev_->vtx_std_xyz->At(selvtx);
    
    p1    = pho1.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    p2    = pho2.p4(vtxpos->X(), vtxpos->Y(), vtxpos->Z());
    dipho = (p1+p2);
    
    PtAll_rank->Fill( dipho.Pt() );
    InvMassAll_rank->Fill( dipho.M() );
    InvMassAll_rank_cat[cat]-> Fill( dipho.M() );
    EtaAll_rank->Fill( p1.Eta() );
    NvtAll_rank->Fill( nvtx );
    if (!isData) NpuAll_rank->Fill(npu);

    // matching within 1cm
    if ( fabs( TrueVertex_Z - vtxpos->Z() ) < 1.) {
      PtGood_rank->Fill( dipho.Pt() );
      EtaGood_rank->Fill( p1.Eta() );
      NvtGood_rank->Fill( nvtx );
      if (!isData) NpuGood_rank->Fill(npu);
      InvMassGood_rank->Fill( dipho.M() );
      InvMassGood_rank_cat[cat]-> Fill( dipho.M() );
    }


    */
  //}// end loop over entries
  
}

//------------------------------------------------------------------------------------------

void HggVertexAnalysis::saveHistos(TFile * fout)
{
  
  fout->cd();
  
  PtAll_sumpt2->Write();
  PtGood_sumpt2->Write();
  InvMassAll_sumpt2->Write();
  InvMassGood_sumpt2->Write();
  EtaAll_sumpt2->Write();
  EtaGood_sumpt2->Write();
  NvtAll_sumpt2->Write();
  NvtGood_sumpt2->Write();
  NpuAll_sumpt2->Write();
  NpuGood_sumpt2->Write();
  
  PtAll_rank->Write();
  PtGood_rank->Write();
  InvMassAll_rank->Write();
  InvMassGood_rank->Write();
  EtaAll_rank->Write();
  EtaGood_rank->Write();
  NvtAll_rank->Write();
  NvtGood_rank->Write();
  NpuAll_rank->Write();
  NpuGood_rank->Write();


  for (int i = 0; i < 4; i++ ){
    InvMassAll_sumpt2_cat[i]-> Write();
    InvMassGood_sumpt2_cat[i]-> Write();
    InvMassAll_rank_cat[i]-> Write();
    InvMassGood_rank_cat[i]-> Write();
  }
  

  
  fout->Close();
  
  return;
}

  





