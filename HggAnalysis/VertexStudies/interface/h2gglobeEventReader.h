//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 11:34:55 2011 by ROOT version 5.27/06b
// from TTree event/Reduced tree
// found on file: /data2/VertexStudies/globe_V09_00_pm_11_07_01_02_red/GluGlu_M-120_FlatPU35.root
//////////////////////////////////////////////////////////

#ifndef h2gglobeEventReader_h
#define h2gglobeEventReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TClonesArray.h"

const int kMaxl1_labels = 1;

class h2gglobeEventReader {
 public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           lumis;
   Int_t           run;
   Int_t           bx;
   TClonesArray    *sc_p4;
   TClonesArray    *sc_xyz;
   Int_t           pho_n;
   Float_t         pho_feta[9][5];   //[pho_n]
   Float_t         pho_crackcorr[9];   //[pho_n]
   Float_t         pho_localcorr[9];   //[pho_n]
   Int_t           pho_isEB[9];   //[pho_n]
   Int_t           pho_isEE[9];   //[pho_n]
   Float_t         pho_see[9];   //[pho_n]
   Float_t         pho_sieie[9];   //[pho_n]
   Float_t         pho_sipip[9];   //[pho_n]
   Float_t         pho_sieip[9];   //[pho_n]
   Float_t         pho_e1x5[9];   //[pho_n]
   Float_t         pho_e2x5[9];   //[pho_n]
   Float_t         pho_e3x3[9];   //[pho_n]
   Float_t         pho_e5x5[9];   //[pho_n]
   Float_t         pho_emaxxtal[9];   //[pho_n]
   Float_t         pho_hoe[9];   //[pho_n]
   Float_t         pho_r1x5[9];   //[pho_n]
   Float_t         pho_r2x5[9];   //[pho_n]
   Float_t         pho_r9[9];   //[pho_n]
   Int_t           pho_isEBGap[9];   //[pho_n]
   Int_t           pho_isEEGap[9];   //[pho_n]
   Int_t           pho_isEBEEGap[9];   //[pho_n]
   Float_t         pho_zernike20[9];   //[pho_n]
   Float_t         pho_zernike42[9];   //[pho_n]
   Float_t         pho_e2nd[9];   //[pho_n]
   Float_t         pho_e2x5right[9];   //[pho_n]
   Float_t         pho_e2x5left[9];   //[pho_n]
   Float_t         pho_e2x5Top[9];   //[pho_n]
   Float_t         pho_e2x5bottom[9];   //[pho_n]
   Float_t         pho_eright[9];   //[pho_n]
   Float_t         pho_eleft[9];   //[pho_n]
   Float_t         pho_etop[9];   //[pho_n]
   Float_t         pho_ebottom[9];   //[pho_n]
   Float_t         pho_seed_time[9];   //[pho_n]
   Float_t         pho_seed_outoftimechi2[9];   //[pho_n]
   Float_t         pho_seed_chi2[9];   //[pho_n]
   Float_t         pho_seed_recoflag[9];   //[pho_n]
   Float_t         pho_seed_severity[9];   //[pho_n]
   Float_t         pho_ecalsumetconedr04[9];   //[pho_n]
   Float_t         pho_hcalsumetconedr04[9];   //[pho_n]
   Float_t         pho_trksumptsolidconedr04[9];   //[pho_n]
   Float_t         pho_trksumpthollowconedr04[9];   //[pho_n]
   Float_t         pho_ntrksolidconedr04[9];   //[pho_n]
   Float_t         pho_ntrkhollowconedr04[9];   //[pho_n]
   Float_t         pho_ecalsumetconedr03[9];   //[pho_n]
   Float_t         pho_hcalsumetconedr03[9];   //[pho_n]
   Float_t         pho_trksumptsolidconedr03[9];   //[pho_n]
   Float_t         pho_trksumpthollowconedr03[9];   //[pho_n]
   Float_t         pho_ntrksolidconedr03[9];   //[pho_n]
   Float_t         pho_ntrkhollowconedr03[9];   //[pho_n]
   Int_t           pho_barrel[9];   //[pho_n]
   Int_t           pho_haspixseed[9];   //[pho_n]
   Int_t           pho_hasconvtks[9];   //[pho_n]
   Int_t           pho_nconv[9];   //[pho_n]
   Int_t           pho_conv_ntracks[9];   //[pho_n]
   Float_t         pho_conv_pairinvmass[9];   //[pho_n]
   Float_t         pho_conv_paircotthetasep[9];   //[pho_n]
   Float_t         pho_conv_eoverp[9];   //[pho_n]
   Float_t         pho_conv_zofprimvtxfromtrks[9];   //[pho_n]
   Float_t         pho_conv_distofminapproach[9];   //[pho_n]
   Float_t         pho_conv_dphitrksatvtx[9];   //[pho_n]
   Float_t         pho_conv_dphitrksatecal[9];   //[pho_n]
   Float_t         pho_conv_detatrksatecal[9];   //[pho_n]
   Float_t         pho_conv_tk1_d0[9];   //[pho_n]
   Float_t         pho_conv_tk1_pout[9];   //[pho_n]
   Float_t         pho_conv_tk1_pin[9];   //[pho_n]
   Float_t         pho_conv_tk2_d0[9];   //[pho_n]
   Float_t         pho_conv_tk2_pout[9];   //[pho_n]
   Float_t         pho_conv_tk2_pin[9];   //[pho_n]
   Float_t         pho_conv_tk1_dz[9];   //[pho_n]
   Float_t         pho_conv_tk1_dzerr[9];   //[pho_n]
   Int_t           pho_conv_tk1_nh[9];   //[pho_n]
   Float_t         pho_conv_tk2_dz[9];   //[pho_n]
   Float_t         pho_conv_tk2_dzerr[9];   //[pho_n]
   Int_t           pho_conv_tk2_nh[9];   //[pho_n]
   Int_t           pho_conv_ch1ch2[9];   //[pho_n]
   Float_t         pho_conv_chi2[9];   //[pho_n]
   Float_t         pho_conv_chi2_probability[9];   //[pho_n]
   Int_t           pho_conv_validvtx[9];   //[pho_n]
   Int_t           pho_conv_MVALikelihood[9];   //[pho_n]
   TClonesArray    *pho_p4;
   TClonesArray    *pho_calopos;
   TClonesArray    *pho_conv_vtx;
   TClonesArray    *pho_conv_pair_momentum;
   TClonesArray    *pho_conv_refitted_momentum;
   TClonesArray    *pho_conv_vertexcorrected_p4;
   Int_t           pho_scind[9];   //[pho_n]
   Int_t           conv_n;
   TClonesArray    *conv_p4;
   Int_t           conv_ntracks[151];   //[conv_n]
   Float_t         conv_pairinvmass[151];   //[conv_n]
   Float_t         conv_paircotthetasep[151];   //[conv_n]
   Float_t         conv_eoverp[151];   //[conv_n]
   Float_t         conv_distofminapproach[151];   //[conv_n]
   Float_t         conv_dphitrksatvtx[151];   //[conv_n]
   Float_t         conv_dphitrksatecal[151];   //[conv_n]
   Float_t         conv_detatrksatecal[151];   //[conv_n]
   Float_t         conv_dxy[151];   //[conv_n]
   Float_t         conv_dz[151];   //[conv_n]
   Float_t         conv_lxy[151];   //[conv_n]
   Float_t         conv_lz[151];   //[conv_n]
   Float_t         conv_zofprimvtxfromtrks[151];   //[conv_n]
   std::vector<std::vector<unsigned short> > *conv_nHitsBeforeVtx;
   Int_t           conv_nSharedHits[151];   //[conv_n]
   Int_t           conv_validvtx[151];   //[conv_n]
   Int_t           conv_MVALikelihood[151];   //[conv_n]
   Float_t         conv_chi2[151];   //[conv_n]
   Float_t         conv_chi2_probability[151];   //[conv_n]
   Float_t         conv_vtx_xErr[151];   //[conv_n]
   Float_t         conv_vtx_yErr[151];   //[conv_n]
   Float_t         conv_vtx_zErr[151];   //[conv_n]
   Float_t         conv_tk1_dz[151];   //[conv_n]
   Float_t         conv_tk2_dz[151];   //[conv_n]
   Float_t         conv_tk1_dzerr[151];   //[conv_n]
   Float_t         conv_tk2_dzerr[151];   //[conv_n]
   Int_t           conv_tk1_nh[151];   //[conv_n]
   Int_t           conv_tk2_nh[151];   //[conv_n]
   Int_t           conv_ch1ch2[151];   //[conv_n]
   Float_t         conv_tk1_d0[151];   //[conv_n]
   Float_t         conv_tk1_pout[151];   //[conv_n]
   Float_t         conv_tk1_pin[151];   //[conv_n]
   Float_t         conv_tk2_d0[151];   //[conv_n]
   Float_t         conv_tk2_pout[151];   //[conv_n]
   Float_t         conv_tk2_pin[151];   //[conv_n]
   TClonesArray    *conv_vtx;
   TClonesArray    *conv_pair_momentum;
   TClonesArray    *conv_refitted_momentum;
   Int_t           process_id;
   Float_t         weight;
   Float_t         pthat;
   Int_t           gp_n;
   TClonesArray    *gp_p4;
   Short_t         gp_status[1667];   //[gp_n]
   Short_t         gp_pdgid[1667];   //[gp_n]
   Short_t         gp_mother[1667];   //[gp_n]
   Int_t           gv_n;
   TClonesArray    *gv_pos;
   Int_t           pu_n;
   std::vector<float>   *pu_zpos;
   std::vector<float>   *pu_sumpt_lowpt;
   std::vector<float>   *pu_sumpt_highpt;
   std::vector<int>     *pu_ntrks_lowpt;
   std::vector<int>     *pu_ntrks_highpt;
   Int_t           vtx_std_n;
   Int_t           vtx_std_ntks[45];   //[vtx_std_n]
   Float_t         vtx_std_x2dof[45];   //[vtx_std_n]
   TClonesArray    *vtx_std_xyz;
   TClonesArray    *vtx_std_dxdydz;
   Float_t         vtx_std_ndof[45];   //[vtx_std_n]
   Float_t         rho;
   TClonesArray    *bs_xyz;
   Float_t         bs_sigmaZ;
   Float_t         bs_x0Error;
   Float_t         bs_y0Error;
   Float_t         bs_z0Error;
   Float_t         bs_sigmaZ0Error;
   Float_t         met_tcmet;
   Float_t         met_phi_tcmet;
   std::vector<unsigned short> *hlt1_bit;
   Int_t           hlt_n;
   std::vector<std::vector<unsigned short> > *hlt_candpath;
   std::vector<std::string>  *hlt_path_names_HLT1;
   TClonesArray    *hlt_p4;
   Int_t           l1emiso_n;
   Float_t         l1emiso_et[4];   //[l1emiso_n]
   Float_t         l1emiso_eta[4];   //[l1emiso_n]
   Float_t         l1emiso_phi[4];   //[l1emiso_n]
   Int_t           l1emnoniso_n;
   Float_t         l1emnoniso_et[4];   //[l1emnoniso_n]
   Float_t         l1emnoniso_eta[4];   //[l1emnoniso_n]
   Float_t         l1emnoniso_phi[4];   //[l1emnoniso_n]
   std::vector<int>     *l1bits_phy;
   Int_t           l1_labels_;
   std::string     l1_labels_first[kMaxl1_labels];
   Int_t           l1_labels_second[kMaxl1_labels];   //[l1_labels_]
   std::vector<std::vector<float> > *vtx_std_diphopt;
   std::vector<std::vector<float> > *vtx_std_nch;
   std::vector<std::vector<float> > *vtx_std_ptmax;
   std::vector<std::vector<float> > *vtx_std_sumpt;
   std::vector<std::vector<float> > *vtx_std_ptvtx;
   std::vector<std::vector<float> > *vtx_std_acosA;
   std::vector<std::vector<float> > *vtx_std_ptasym;
   std::vector<std::vector<float> > *vtx_std_ptbal;
   std::vector<std::vector<float> > *vtx_std_nchthr;
   std::vector<std::vector<float> > *vtx_std_ptmax3;
   std::vector<std::vector<float> > *vtx_std_thrust;
   std::vector<std::vector<float> > *vtx_std_sumweight;
   std::vector<std::vector<float> > *vtx_std_sumpt2;
   std::vector<std::vector<float> > *vtx_std_ptratio;
   std::vector<std::vector<float> > *vtx_std_pzasym;
   std::vector<std::vector<float> > *vtx_std_spher;
   std::vector<std::vector<float> > *vtx_std_aplan;
   std::vector<std::vector<float> > *vtx_std_sumpr;
   std::vector<std::vector<float> > *vtx_std_sumawy;
   std::vector<std::vector<float> > *vtx_std_sumtrv;
   std::vector<std::vector<float> > *vtx_std_sumtwd;
   std::vector<std::vector<float> > *vtx_std_awytwdasym;
   Int_t           vtx_std_pho1;
   Int_t           vtx_std_pho2;
   std::vector<int>     *pho_matchingConv;
   std::vector<std::vector<int> > *vtx_std_ranked_list;
   Int_t           vtx_std_sel;
   std::vector<std::vector<float> > *pho_tkiso_recvtx_030_002_0000_10_01;
   Float_t         pho_tkiso_badvtx_040_002_0000_10_01[9];   //[pho_n]
   Int_t           pho_tkiso_badvtx_id[9];   //[pho_n]
   Float_t         pho_drtotk_25_99[9];   //[pho_n]
   Int_t           dipho_n;
   Int_t           dipho_leadind[6];   //[dipho_n]
   Int_t           dipho_subleadind[6];   //[dipho_n]
   Int_t           dipho_vtxind[6];   //[dipho_n]
   Float_t         dipho_sumpt[6];   //[dipho_n]
   std::vector<std::vector<short> > *pho_cic6cutlevel_lead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic6passcuts_lead;
   std::vector<std::vector<short> > *pho_cic6cutlevel_sublead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic6passcuts_sublead;
   std::vector<std::vector<short> > *pho_cic4cutlevel_lead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4passcuts_lead;
   std::vector<std::vector<short> > *pho_cic4cutlevel_sublead;
   std::vector<std::vector<std::vector<unsigned int> > > *pho_cic4passcuts_sublead;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_run;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_sc_p4;   //!
   TBranch        *b_sc_xyz;   //!
   TBranch        *b_pho_n;   //!
   TBranch        *b_pho_feta;   //!
   TBranch        *b_pho_crackcorr;   //!
   TBranch        *b_pho_localcorr;   //!
   TBranch        *b_pho_isEB;   //!
   TBranch        *b_pho_isEE;   //!
   TBranch        *b_pho_see;   //!
   TBranch        *b_pho_sieie;   //!
   TBranch        *b_pho_sipip;   //!
   TBranch        *b_pho_sieip;   //!
   TBranch        *b_pho_e1x5;   //!
   TBranch        *b_pho_e2x5;   //!
   TBranch        *b_pho_e3x3;   //!
   TBranch        *b_pho_e5x5;   //!
   TBranch        *b_pho_emaxxtal;   //!
   TBranch        *b_pho_hoe;   //!
   TBranch        *b_pho_r1x5;   //!
   TBranch        *b_pho_r2x5;   //!
   TBranch        *b_pho_r9;   //!
   TBranch        *b_pho_isEBGap;   //!
   TBranch        *b_pho_isEEGap;   //!
   TBranch        *b_pho_isEBEEGap;   //!
   TBranch        *b_pho_zernike20;   //!
   TBranch        *b_pho_zernike42;   //!
   TBranch        *b_pho_e2nd;   //!
   TBranch        *b_pho_e2x5right;   //!
   TBranch        *b_pho_e2x5left;   //!
   TBranch        *b_pho_e2x5Top;   //!
   TBranch        *b_pho_e2x5bottom;   //!
   TBranch        *b_pho_eright;   //!
   TBranch        *b_pho_eleft;   //!
   TBranch        *b_pho_etop;   //!
   TBranch        *b_pho_ebottom;   //!
   TBranch        *b_pho_seed_time;   //!
   TBranch        *b_pho_seed_outoftimechi2;   //!
   TBranch        *b_pho_seed_chi2;   //!
   TBranch        *b_pho_seed_recoflag;   //!
   TBranch        *b_pho_seed_severity;   //!
   TBranch        *b_pho_ecalsumetconedr04;   //!
   TBranch        *b_pho_hcalsumetconedr04;   //!
   TBranch        *b_pho_trksumptsolidconedr04;   //!
   TBranch        *b_pho_trksumpthollowconedr04;   //!
   TBranch        *b_pho_ntrksolidconedr04;   //!
   TBranch        *b_pho_ntrkhollowconedr04;   //!
   TBranch        *b_pho_ecalsumetconedr03;   //!
   TBranch        *b_pho_hcalsumetconedr03;   //!
   TBranch        *b_pho_trksumptsolidconedr03;   //!
   TBranch        *b_pho_trksumpthollowconedr03;   //!
   TBranch        *b_pho_ntrksolidconedr03;   //!
   TBranch        *b_pho_ntrkhollowconedr03;   //!
   TBranch        *b_pho_barrel;   //!
   TBranch        *b_pho_haspixseed;   //!
   TBranch        *b_pho_hasconvtks;   //!
   TBranch        *b_pho_nconv;   //!
   TBranch        *b_pho_conv_ntracks;   //!
   TBranch        *b_pho_conv_pairinvmass;   //!
   TBranch        *b_pho_conv_paircotthetasep;   //!
   TBranch        *b_pho_conv_eoverp;   //!
   TBranch        *b_pho_conv_zofprimvtxfromtrks;   //!
   TBranch        *b_pho_conv_distofminapproach;   //!
   TBranch        *b_pho_conv_dphitrksatvtx;   //!
   TBranch        *b_pho_conv_dphitrksatecal;   //!
   TBranch        *b_pho_conv_detatrksatecal;   //!
   TBranch        *b_pho_conv_tk1_d0;   //!
   TBranch        *b_pho_conv_tk1_pout;   //!
   TBranch        *b_pho_conv_tk1_pin;   //!
   TBranch        *b_pho_conv_tk2_d0;   //!
   TBranch        *b_pho_conv_tk2_pout;   //!
   TBranch        *b_pho_conv_tk2_pin;   //!
   TBranch        *b_pho_conv_tk1_dz;   //!
   TBranch        *b_pho_conv_tk1_dzerr;   //!
   TBranch        *b_pho_conv_tk1_nh;   //!
   TBranch        *b_pho_conv_tk2_dz;   //!
   TBranch        *b_pho_conv_tk2_dzerr;   //!
   TBranch        *b_pho_conv_tk2_nh;   //!
   TBranch        *b_pho_conv_ch1ch2;   //!
   TBranch        *b_pho_conv_chi2;   //!
   TBranch        *b_pho_conv_chi2_probability;   //!
   TBranch        *b_pho_conv_validvtx;   //!
   TBranch        *b_pho_conv_MVALikelihood;   //!
   TBranch        *b_pho_p4;   //!
   TBranch        *b_pho_calopos;   //!
   TBranch        *b_pho_conv_vtx;   //!
   TBranch        *b_pho_conv_pair_momentum;   //!
   TBranch        *b_pho_conv_refitted_momentum;   //!
   TBranch        *b_pho_conv_vertexcorrected_p4;   //!
   TBranch        *b_pho_scind;   //!
   TBranch        *b_conv_n;   //!
   TBranch        *b_conv_p4;   //!
   TBranch        *b_conv_ntracks;   //!
   TBranch        *b_conv_pairinvmass;   //!
   TBranch        *b_conv_paircotthetasep;   //!
   TBranch        *b_conv_eoverp;   //!
   TBranch        *b_conv_distofminapproach;   //!
   TBranch        *b_conv_dphitrksatvtx;   //!
   TBranch        *b_conv_dphitrksatecal;   //!
   TBranch        *b_conv_detatrksatecal;   //!
   TBranch        *b_conv_dxy;   //!
   TBranch        *b_conv_dz;   //!
   TBranch        *b_conv_lxy;   //!
   TBranch        *b_conv_lz;   //!
   TBranch        *b_conv_zofprimvtxfromtrks;   //!
   TBranch        *b_conv_nHitsBeforeVtx;   //!
   TBranch        *b_conv_nSharedHits;   //!
   TBranch        *b_conv_validvtx;   //!
   TBranch        *b_conv_MVALikelihood;   //!
   TBranch        *b_conv_chi2;   //!
   TBranch        *b_conv_chi2_probability;   //!
   TBranch        *b_conv_vtx_xErr;   //!
   TBranch        *b_conv_vtx_yErr;   //!
   TBranch        *b_conv_vtx_zErr;   //!
   TBranch        *b_conv_tk1_dz;   //!
   TBranch        *b_conv_tk2_dz;   //!
   TBranch        *b_conv_tk1_dzerr;   //!
   TBranch        *b_conv_tk2_dzerr;   //!
   TBranch        *b_conv_tk1_nh;   //!
   TBranch        *b_conv_tk2_nh;   //!
   TBranch        *b_conv_ch1ch2;   //!
   TBranch        *b_conv_tk1_d0;   //!
   TBranch        *b_conv_tk1_pout;   //!
   TBranch        *b_conv_tk1_pin;   //!
   TBranch        *b_conv_tk2_d0;   //!
   TBranch        *b_conv_tk2_pout;   //!
   TBranch        *b_conv_tk2_pin;   //!
   TBranch        *b_conv_vtx;   //!
   TBranch        *b_conv_pair_momentum;   //!
   TBranch        *b_conv_refitted_momentum;   //!
   TBranch        *b_process_id;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_gp_n;   //!
   TBranch        *b_gp_p4;   //!
   TBranch        *b_gp_status;   //!
   TBranch        *b_gp_pdgid;   //!
   TBranch        *b_gp_mother;   //!
   TBranch        *b_gv_n;   //!
   TBranch        *b_gv_pos;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_zpos;   //!
   TBranch        *b_pu_sumpt_lowpt;   //!
   TBranch        *b_pu_sumpt_highpt;   //!
   TBranch        *b_pu_ntrks_lowpt;   //!
   TBranch        *b_pu_ntrks_highpt;   //!
   TBranch        *b_vtx_std_n;   //!
   TBranch        *b_vtx_std_ntks;   //!
   TBranch        *b_vtx_std_x2dof;   //!
   TBranch        *b_vtx_std_xyz;   //!
   TBranch        *b_vtx_std_dxdydz;   //!
   TBranch        *b_vtx_std_ndof;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_bs_xyz;   //!
   TBranch        *b_bs_sigmaZ;   //!
   TBranch        *b_bs_x0Error;   //!
   TBranch        *b_bs_y0Error;   //!
   TBranch        *b_bs_z0Error;   //!
   TBranch        *b_bs_sigmaZ0Error;   //!
   TBranch        *b_met_tcmet;   //!
   TBranch        *b_met_phi_tcmet;   //!
   TBranch        *b_hlt1_bit;   //!
   TBranch        *b_hlt_n;   //!
   TBranch        *b_hlt_candpath;   //!
   TBranch        *b_hlt_path_names_HLT1;   //!
   TBranch        *b_hlt_p4;   //!
   TBranch        *b_l1emiso_n;   //!
   TBranch        *b_l1emiso_et;   //!
   TBranch        *b_l1emiso_eta;   //!
   TBranch        *b_l1emiso_phi;   //!
   TBranch        *b_l1emnoniso_n;   //!
   TBranch        *b_l1emnoniso_et;   //!
   TBranch        *b_l1emnoniso_eta;   //!
   TBranch        *b_l1emnoniso_phi;   //!
   TBranch        *b_l1bits_phy;   //!
   TBranch        *b_l1_labels_;   //!
   TBranch        *b_l1_labels_first;   //!
   TBranch        *b_l1_labels_second;   //!
   TBranch        *b_vtx_std_diphopt;   //!
   TBranch        *b_vtx_std_nch;   //!
   TBranch        *b_vtx_std_ptmax;   //!
   TBranch        *b_vtx_std_sumpt;   //!
   TBranch        *b_vtx_std_ptvtx;   //!
   TBranch        *b_vtx_std_acosA;   //!
   TBranch        *b_vtx_std_ptasym;   //!
   TBranch        *b_vtx_std_ptbal;   //!
   TBranch        *b_vtx_std_nchthr;   //!
   TBranch        *b_vtx_std_ptmax3;   //!
   TBranch        *b_vtx_std_thrust;   //!
   TBranch        *b_vtx_std_sumweight;   //!
   TBranch        *b_vtx_std_sumpt2;   //!
   TBranch        *b_vtx_std_ptratio;   //!
   TBranch        *b_vtx_std_pzasym;   //!
   TBranch        *b_vtx_std_spher;   //!
   TBranch        *b_vtx_std_aplan;   //!
   TBranch        *b_vtx_std_sumpr;   //!
   TBranch        *b_vtx_std_sumawy;   //!
   TBranch        *b_vtx_std_sumtrv;   //!
   TBranch        *b_vtx_std_sumtwd;   //!
   TBranch        *b_vtx_std_awytwdasym;   //!
   TBranch        *b_vtx_std_pho1;   //!
   TBranch        *b_vtx_std_pho2;   //!
   TBranch        *b_pho_matchingConv;   //!
   TBranch        *b_vtx_std_ranked_list;   //!
   TBranch        *b_vtx_std_sel;   //!
   TBranch        *b_pho_tkiso_recvtx_030_002_0000_10_01;   //!
   TBranch        *b_pho_tkiso_badvtx_040_002_0000_10_01;   //!
   TBranch        *b_pho_tkiso_badvtx_id;   //!
   TBranch        *b_pho_drtotk_25_99;   //!
   TBranch        *b_dipho_n;   //!
   TBranch        *b_dipho_leadind;   //!
   TBranch        *b_dipho_subleadind;   //!
   TBranch        *b_dipho_vtxind;   //!
   TBranch        *b_dipho_sumpt;   //!
   TBranch        *b_pho_cic6cutlevel_lead;   //!
   TBranch        *b_pho_cic6passcuts_lead;   //!
   TBranch        *b_pho_cic6cutlevel_sublead;   //!
   TBranch        *b_pho_cic6passcuts_sublead;   //!
   TBranch        *b_pho_cic4cutlevel_lead;   //!
   TBranch        *b_pho_cic4passcuts_lead;   //!
   TBranch        *b_pho_cic4cutlevel_sublead;   //!
   TBranch        *b_pho_cic4passcuts_sublead;   //!

   h2gglobeEventReader(TTree *tree=0);
   virtual ~h2gglobeEventReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


};

#endif





