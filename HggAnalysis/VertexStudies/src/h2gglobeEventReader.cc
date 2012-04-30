//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 11:34:55 2011 by ROOT version 5.27/06b
// from TTree event/Reduced tree
// found on file: /data2/VertexStudies/globe_V09_00_pm_11_07_01_02_red/GluGlu_M-120_FlatPU35.root
//////////////////////////////////////////////////////////

#include "h2gglobeEventReader.h"

using namespace std;

//#ifdef h2gglobeEventReader_cc

h2gglobeEventReader::h2gglobeEventReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
 //      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data2/VertexStudies/globe_V09_00_pm_11_07_01_02_red/GluGlu_M-120_FlatPU35.root");
//       if (!f) {
//          f = new TFile("/data2/VertexStudies/globe_V09_00_pm_11_07_01_02_red/GluGlu_M-120_FlatPU35.root");
     //}
      tree = (TTree*)gDirectory->Get("event");
   }
   
   
   Init(tree);
}

h2gglobeEventReader::~h2gglobeEventReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t h2gglobeEventReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Int_t h2gglobeEventReader::GetEntries()
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntries();
}

Long64_t h2gglobeEventReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void h2gglobeEventReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   sc_p4 = 0;
   sc_xyz = 0;
   pho_p4 = 0;
   pho_calopos = 0;
   pho_conv_vtx = 0;
   pho_conv_pair_momentum = 0;
   pho_conv_refitted_momentum = 0;
   pho_conv_vertexcorrected_p4 = 0;
   conv_p4 = 0;
   conv_nHitsBeforeVtx = 0;
   conv_vtx = 0;
   conv_pair_momentum = 0;
   conv_refitted_momentum = 0;
   gp_p4 = 0;
   gv_pos = 0;
   pu_zpos = 0;
   pu_sumpt_lowpt = 0;
   pu_sumpt_highpt = 0;
   pu_ntrks_lowpt = 0;
   pu_ntrks_highpt = 0;
   vtx_std_xyz = 0;
   vtx_std_dxdydz = 0;
   bs_xyz = 0;
   hlt1_bit = 0;
   hlt_candpath = 0;
   hlt_path_names_HLT1 = 0;
   hlt_p4 = 0;
   l1bits_phy = 0;
   jet_algoPF1_p4 = 0;
   vtx_std_diphopt = 0;
   vtx_std_nch = 0;
   vtx_std_ptmax = 0;
   vtx_std_sumpt = 0;
   vtx_std_ptvtx = 0;
   vtx_std_acosA = 0;
   vtx_std_ptasym = 0;
   vtx_std_ptbal = 0;
   vtx_std_nchthr = 0;
   vtx_std_ptmax3 = 0;
   vtx_std_thrust = 0;
   vtx_std_sumweight = 0;
   vtx_std_sumpt2 = 0;
   vtx_std_ptratio = 0;
   vtx_std_pzasym = 0;
   vtx_std_spher = 0;
   vtx_std_aplan = 0;
   vtx_std_sumpr = 0;
   vtx_std_sumawy = 0;
   vtx_std_sumtrv = 0;
   vtx_std_sumtwd = 0;
   vtx_std_awytwdasym = 0;
   pho_matchingConv = 0;
   vtx_std_ranked_list = 0;
   pho_tkiso_recvtx_030_002_0000_10_01 = 0;
   pho_cic6cutlevel_lead = 0;
   pho_cic6passcuts_lead = 0;
   pho_cic6cutlevel_sublead = 0;
   pho_cic6passcuts_sublead = 0;
   pho_cic4cutlevel_lead = 0;
   pho_cic4passcuts_lead = 0;
   pho_cic4cutlevel_sublead = 0;
   pho_cic4passcuts_sublead = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("sc_p4", &sc_p4, &b_sc_p4);
   fChain->SetBranchAddress("sc_xyz", &sc_xyz, &b_sc_xyz);
   fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
   fChain->SetBranchAddress("pho_feta", pho_feta, &b_pho_feta);
   fChain->SetBranchAddress("pho_crackcorr", pho_crackcorr, &b_pho_crackcorr);
   fChain->SetBranchAddress("pho_localcorr", pho_localcorr, &b_pho_localcorr);
   fChain->SetBranchAddress("pho_isEB", pho_isEB, &b_pho_isEB);
   fChain->SetBranchAddress("pho_isEE", pho_isEE, &b_pho_isEE);
   fChain->SetBranchAddress("pho_see", pho_see, &b_pho_see);
   fChain->SetBranchAddress("pho_sieie", pho_sieie, &b_pho_sieie);
   fChain->SetBranchAddress("pho_sipip", pho_sipip, &b_pho_sipip);
   fChain->SetBranchAddress("pho_sieip", pho_sieip, &b_pho_sieip);
   fChain->SetBranchAddress("pho_e1x5", pho_e1x5, &b_pho_e1x5);
   fChain->SetBranchAddress("pho_e2x5", pho_e2x5, &b_pho_e2x5);
   fChain->SetBranchAddress("pho_e3x3", pho_e3x3, &b_pho_e3x3);
   fChain->SetBranchAddress("pho_e5x5", pho_e5x5, &b_pho_e5x5);
   fChain->SetBranchAddress("pho_emaxxtal", pho_emaxxtal, &b_pho_emaxxtal);
   fChain->SetBranchAddress("pho_hoe", pho_hoe, &b_pho_hoe);
   fChain->SetBranchAddress("pho_r1x5", pho_r1x5, &b_pho_r1x5);
   fChain->SetBranchAddress("pho_r2x5", pho_r2x5, &b_pho_r2x5);
   fChain->SetBranchAddress("pho_r9", pho_r9, &b_pho_r9);
   //fChain->SetBranchAddress("pho_isEBGap", pho_isEBGap, &b_pho_isEBGap);
   //fChain->SetBranchAddress("pho_isEEGap", pho_isEEGap, &b_pho_isEEGap);
   //fChain->SetBranchAddress("pho_isEBEEGap", pho_isEBEEGap, &b_pho_isEBEEGap);
   fChain->SetBranchAddress("pho_zernike20", pho_zernike20, &b_pho_zernike20);
   fChain->SetBranchAddress("pho_zernike42", pho_zernike42, &b_pho_zernike42);
   fChain->SetBranchAddress("pho_e2nd", pho_e2nd, &b_pho_e2nd);
   fChain->SetBranchAddress("pho_e2x5right", pho_e2x5right, &b_pho_e2x5right);
   fChain->SetBranchAddress("pho_e2x5left", pho_e2x5left, &b_pho_e2x5left);
   fChain->SetBranchAddress("pho_e2x5Top", pho_e2x5Top, &b_pho_e2x5Top);
   fChain->SetBranchAddress("pho_e2x5bottom", pho_e2x5bottom, &b_pho_e2x5bottom);
   fChain->SetBranchAddress("pho_eright", pho_eright, &b_pho_eright);
   fChain->SetBranchAddress("pho_eleft", pho_eleft, &b_pho_eleft);
   fChain->SetBranchAddress("pho_etop", pho_etop, &b_pho_etop);
   fChain->SetBranchAddress("pho_ebottom", pho_ebottom, &b_pho_ebottom);
   fChain->SetBranchAddress("pho_seed_time", pho_seed_time, &b_pho_seed_time);
   fChain->SetBranchAddress("pho_seed_outoftimechi2", pho_seed_outoftimechi2, &b_pho_seed_outoftimechi2);
   fChain->SetBranchAddress("pho_seed_chi2", pho_seed_chi2, &b_pho_seed_chi2);
   fChain->SetBranchAddress("pho_seed_recoflag", pho_seed_recoflag, &b_pho_seed_recoflag);
   fChain->SetBranchAddress("pho_seed_severity", pho_seed_severity, &b_pho_seed_severity);
   fChain->SetBranchAddress("pho_ecalsumetconedr04", pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
   fChain->SetBranchAddress("pho_hcalsumetconedr04", pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
   fChain->SetBranchAddress("pho_trksumptsolidconedr04", pho_trksumptsolidconedr04, &b_pho_trksumptsolidconedr04);
   fChain->SetBranchAddress("pho_trksumpthollowconedr04", pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
   fChain->SetBranchAddress("pho_ntrksolidconedr04", pho_ntrksolidconedr04, &b_pho_ntrksolidconedr04);
   fChain->SetBranchAddress("pho_ntrkhollowconedr04", pho_ntrkhollowconedr04, &b_pho_ntrkhollowconedr04);
   fChain->SetBranchAddress("pho_ecalsumetconedr03", pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
   fChain->SetBranchAddress("pho_hcalsumetconedr03", pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
   fChain->SetBranchAddress("pho_trksumptsolidconedr03", pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
   fChain->SetBranchAddress("pho_trksumpthollowconedr03", pho_trksumpthollowconedr03, &b_pho_trksumpthollowconedr03);
   fChain->SetBranchAddress("pho_ntrksolidconedr03", pho_ntrksolidconedr03, &b_pho_ntrksolidconedr03);
   fChain->SetBranchAddress("pho_ntrkhollowconedr03", pho_ntrkhollowconedr03, &b_pho_ntrkhollowconedr03);
   fChain->SetBranchAddress("pho_barrel", pho_barrel, &b_pho_barrel);
   fChain->SetBranchAddress("pho_haspixseed", pho_haspixseed, &b_pho_haspixseed);
   fChain->SetBranchAddress("pho_hasconvtks", pho_hasconvtks, &b_pho_hasconvtks);
   fChain->SetBranchAddress("pho_nconv", pho_nconv, &b_pho_nconv);
   fChain->SetBranchAddress("pho_conv_ntracks", pho_conv_ntracks, &b_pho_conv_ntracks);
   fChain->SetBranchAddress("pho_conv_pairinvmass", pho_conv_pairinvmass, &b_pho_conv_pairinvmass);
   fChain->SetBranchAddress("pho_conv_paircotthetasep", pho_conv_paircotthetasep, &b_pho_conv_paircotthetasep);
   fChain->SetBranchAddress("pho_conv_eoverp", pho_conv_eoverp, &b_pho_conv_eoverp);
   fChain->SetBranchAddress("pho_conv_zofprimvtxfromtrks", pho_conv_zofprimvtxfromtrks, &b_pho_conv_zofprimvtxfromtrks);
   fChain->SetBranchAddress("pho_conv_distofminapproach", pho_conv_distofminapproach, &b_pho_conv_distofminapproach);
   fChain->SetBranchAddress("pho_conv_dphitrksatvtx", pho_conv_dphitrksatvtx, &b_pho_conv_dphitrksatvtx);
   fChain->SetBranchAddress("pho_conv_dphitrksatecal", pho_conv_dphitrksatecal, &b_pho_conv_dphitrksatecal);
   fChain->SetBranchAddress("pho_conv_detatrksatecal", pho_conv_detatrksatecal, &b_pho_conv_detatrksatecal);
   fChain->SetBranchAddress("pho_conv_tk1_d0", pho_conv_tk1_d0, &b_pho_conv_tk1_d0);
   fChain->SetBranchAddress("pho_conv_tk1_pout", pho_conv_tk1_pout, &b_pho_conv_tk1_pout);
   fChain->SetBranchAddress("pho_conv_tk1_pin", pho_conv_tk1_pin, &b_pho_conv_tk1_pin);
   fChain->SetBranchAddress("pho_conv_tk2_d0", pho_conv_tk2_d0, &b_pho_conv_tk2_d0);
   fChain->SetBranchAddress("pho_conv_tk2_pout", pho_conv_tk2_pout, &b_pho_conv_tk2_pout);
   fChain->SetBranchAddress("pho_conv_tk2_pin", pho_conv_tk2_pin, &b_pho_conv_tk2_pin);
   fChain->SetBranchAddress("pho_conv_tk1_dz", pho_conv_tk1_dz, &b_pho_conv_tk1_dz);
   fChain->SetBranchAddress("pho_conv_tk1_dzerr", pho_conv_tk1_dzerr, &b_pho_conv_tk1_dzerr);
   fChain->SetBranchAddress("pho_conv_tk1_nh", pho_conv_tk1_nh, &b_pho_conv_tk1_nh);
   fChain->SetBranchAddress("pho_conv_tk2_dz", pho_conv_tk2_dz, &b_pho_conv_tk2_dz);
   fChain->SetBranchAddress("pho_conv_tk2_dzerr", pho_conv_tk2_dzerr, &b_pho_conv_tk2_dzerr);
   fChain->SetBranchAddress("pho_conv_tk2_nh", pho_conv_tk2_nh, &b_pho_conv_tk2_nh);
   fChain->SetBranchAddress("pho_conv_ch1ch2", pho_conv_ch1ch2, &b_pho_conv_ch1ch2);
   fChain->SetBranchAddress("pho_conv_chi2", pho_conv_chi2, &b_pho_conv_chi2);
   fChain->SetBranchAddress("pho_conv_chi2_probability", pho_conv_chi2_probability, &b_pho_conv_chi2_probability);
   fChain->SetBranchAddress("pho_conv_validvtx", pho_conv_validvtx, &b_pho_conv_validvtx);
   fChain->SetBranchAddress("pho_conv_MVALikelihood", pho_conv_MVALikelihood, &b_pho_conv_MVALikelihood);
   fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
   fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
   fChain->SetBranchAddress("pho_conv_vtx", &pho_conv_vtx, &b_pho_conv_vtx);
   fChain->SetBranchAddress("pho_conv_pair_momentum", &pho_conv_pair_momentum, &b_pho_conv_pair_momentum);
   fChain->SetBranchAddress("pho_conv_refitted_momentum", &pho_conv_refitted_momentum, &b_pho_conv_refitted_momentum);
   fChain->SetBranchAddress("pho_conv_vertexcorrected_p4", &pho_conv_vertexcorrected_p4, &b_pho_conv_vertexcorrected_p4);
   fChain->SetBranchAddress("pho_scind", pho_scind, &b_pho_scind);
   fChain->SetBranchAddress("conv_n", &conv_n, &b_conv_n);
   fChain->SetBranchAddress("conv_p4", &conv_p4, &b_conv_p4);
   fChain->SetBranchAddress("conv_ntracks", conv_ntracks, &b_conv_ntracks);
   fChain->SetBranchAddress("conv_pairinvmass", conv_pairinvmass, &b_conv_pairinvmass);
   fChain->SetBranchAddress("conv_paircotthetasep", conv_paircotthetasep, &b_conv_paircotthetasep);
   fChain->SetBranchAddress("conv_eoverp", conv_eoverp, &b_conv_eoverp);
   fChain->SetBranchAddress("conv_distofminapproach", conv_distofminapproach, &b_conv_distofminapproach);
   fChain->SetBranchAddress("conv_dphitrksatvtx", conv_dphitrksatvtx, &b_conv_dphitrksatvtx);
   fChain->SetBranchAddress("conv_dphitrksatecal", conv_dphitrksatecal, &b_conv_dphitrksatecal);
   fChain->SetBranchAddress("conv_detatrksatecal", conv_detatrksatecal, &b_conv_detatrksatecal);
   fChain->SetBranchAddress("conv_dxy", conv_dxy, &b_conv_dxy);
   fChain->SetBranchAddress("conv_dz", conv_dz, &b_conv_dz);
   fChain->SetBranchAddress("conv_lxy", conv_lxy, &b_conv_lxy);
   fChain->SetBranchAddress("conv_lz", conv_lz, &b_conv_lz);
   fChain->SetBranchAddress("conv_zofprimvtxfromtrks", conv_zofprimvtxfromtrks, &b_conv_zofprimvtxfromtrks);
   fChain->SetBranchAddress("conv_nHitsBeforeVtx", &conv_nHitsBeforeVtx, &b_conv_nHitsBeforeVtx);
   fChain->SetBranchAddress("conv_nSharedHits", conv_nSharedHits, &b_conv_nSharedHits);
   fChain->SetBranchAddress("conv_validvtx", conv_validvtx, &b_conv_validvtx);
   fChain->SetBranchAddress("conv_MVALikelihood", conv_MVALikelihood, &b_conv_MVALikelihood);
   fChain->SetBranchAddress("conv_chi2", conv_chi2, &b_conv_chi2);
   fChain->SetBranchAddress("conv_chi2_probability", conv_chi2_probability, &b_conv_chi2_probability);
   fChain->SetBranchAddress("conv_vtx_xErr", conv_vtx_xErr, &b_conv_vtx_xErr);
   fChain->SetBranchAddress("conv_vtx_yErr", conv_vtx_yErr, &b_conv_vtx_yErr);
   fChain->SetBranchAddress("conv_vtx_zErr", conv_vtx_zErr, &b_conv_vtx_zErr);
   fChain->SetBranchAddress("conv_tk1_dz", conv_tk1_dz, &b_conv_tk1_dz);
   fChain->SetBranchAddress("conv_tk2_dz", conv_tk2_dz, &b_conv_tk2_dz);
   fChain->SetBranchAddress("conv_tk1_dzerr", conv_tk1_dzerr, &b_conv_tk1_dzerr);
   fChain->SetBranchAddress("conv_tk2_dzerr", conv_tk2_dzerr, &b_conv_tk2_dzerr);
   fChain->SetBranchAddress("conv_tk1_nh", conv_tk1_nh, &b_conv_tk1_nh);
   fChain->SetBranchAddress("conv_tk2_nh", conv_tk2_nh, &b_conv_tk2_nh);
   fChain->SetBranchAddress("conv_ch1ch2", conv_ch1ch2, &b_conv_ch1ch2);
   fChain->SetBranchAddress("conv_tk1_d0", conv_tk1_d0, &b_conv_tk1_d0);
   fChain->SetBranchAddress("conv_tk1_pout", conv_tk1_pout, &b_conv_tk1_pout);
   fChain->SetBranchAddress("conv_tk1_pin", conv_tk1_pin, &b_conv_tk1_pin);
   fChain->SetBranchAddress("conv_tk2_d0", conv_tk2_d0, &b_conv_tk2_d0);
   fChain->SetBranchAddress("conv_tk2_pout", conv_tk2_pout, &b_conv_tk2_pout);
   fChain->SetBranchAddress("conv_tk2_pin", conv_tk2_pin, &b_conv_tk2_pin);
   fChain->SetBranchAddress("conv_vtx", &conv_vtx, &b_conv_vtx);
   fChain->SetBranchAddress("conv_pair_momentum", &conv_pair_momentum, &b_conv_pair_momentum);
   fChain->SetBranchAddress("conv_refitted_momentum", &conv_refitted_momentum, &b_conv_refitted_momentum);
   fChain->SetBranchAddress("process_id", &process_id, &b_process_id);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("gp_n", &gp_n, &b_gp_n);
   fChain->SetBranchAddress("gp_p4", &gp_p4, &b_gp_p4);
   fChain->SetBranchAddress("gp_status", gp_status, &b_gp_status);
   fChain->SetBranchAddress("gp_pdgid", gp_pdgid, &b_gp_pdgid);
   fChain->SetBranchAddress("gp_mother", gp_mother, &b_gp_mother);
   fChain->SetBranchAddress("gv_n", &gv_n, &b_gv_n);
   fChain->SetBranchAddress("gv_pos", &gv_pos, &b_gv_pos);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("pu_zpos", &pu_zpos, &b_pu_zpos);
   fChain->SetBranchAddress("pu_sumpt_lowpt", &pu_sumpt_lowpt, &b_pu_sumpt_lowpt);
   fChain->SetBranchAddress("pu_sumpt_highpt", &pu_sumpt_highpt, &b_pu_sumpt_highpt);
   fChain->SetBranchAddress("pu_ntrks_lowpt", &pu_ntrks_lowpt, &b_pu_ntrks_lowpt);
   fChain->SetBranchAddress("pu_ntrks_highpt", &pu_ntrks_highpt, &b_pu_ntrks_highpt);
   fChain->SetBranchAddress("vtx_std_n", &vtx_std_n, &b_vtx_std_n);
   fChain->SetBranchAddress("vtx_std_ntks", vtx_std_ntks, &b_vtx_std_ntks);
   fChain->SetBranchAddress("vtx_std_x2dof", vtx_std_x2dof, &b_vtx_std_x2dof);
   fChain->SetBranchAddress("vtx_std_xyz", &vtx_std_xyz, &b_vtx_std_xyz);
   fChain->SetBranchAddress("vtx_std_dxdydz", &vtx_std_dxdydz, &b_vtx_std_dxdydz);
   fChain->SetBranchAddress("vtx_std_ndof", vtx_std_ndof, &b_vtx_std_ndof);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("bs_xyz", &bs_xyz, &b_bs_xyz);
   fChain->SetBranchAddress("bs_sigmaZ", &bs_sigmaZ, &b_bs_sigmaZ);
   fChain->SetBranchAddress("bs_x0Error", &bs_x0Error, &b_bs_x0Error);
   fChain->SetBranchAddress("bs_y0Error", &bs_y0Error, &b_bs_y0Error);
   fChain->SetBranchAddress("bs_z0Error", &bs_z0Error, &b_bs_z0Error);
   fChain->SetBranchAddress("bs_sigmaZ0Error", &bs_sigmaZ0Error, &b_bs_sigmaZ0Error);
   fChain->SetBranchAddress("met_tcmet", &met_tcmet, &b_met_tcmet);
   fChain->SetBranchAddress("met_phi_tcmet", &met_phi_tcmet, &b_met_phi_tcmet);
   fChain->SetBranchAddress("hlt1_bit", &hlt1_bit, &b_hlt1_bit);
   fChain->SetBranchAddress("hlt_n", &hlt_n, &b_hlt_n);
   fChain->SetBranchAddress("hlt_candpath", &hlt_candpath, &b_hlt_candpath);
   fChain->SetBranchAddress("hlt_path_names_HLT1", &hlt_path_names_HLT1, &b_hlt_path_names_HLT1);
   fChain->SetBranchAddress("hlt_p4", &hlt_p4, &b_hlt_p4);
   fChain->SetBranchAddress("l1emiso_n", &l1emiso_n, &b_l1emiso_n);
   fChain->SetBranchAddress("l1emiso_et", l1emiso_et, &b_l1emiso_et);
   fChain->SetBranchAddress("l1emiso_eta", l1emiso_eta, &b_l1emiso_eta);
   fChain->SetBranchAddress("l1emiso_phi", l1emiso_phi, &b_l1emiso_phi);
   fChain->SetBranchAddress("l1emnoniso_n", &l1emnoniso_n, &b_l1emnoniso_n);
   fChain->SetBranchAddress("l1emnoniso_et", l1emnoniso_et, &b_l1emnoniso_et);
   fChain->SetBranchAddress("l1emnoniso_eta", l1emnoniso_eta, &b_l1emnoniso_eta);
   fChain->SetBranchAddress("l1emnoniso_phi", l1emnoniso_phi, &b_l1emnoniso_phi);
   fChain->SetBranchAddress("l1bits_phy", &l1bits_phy, &b_l1bits_phy);
   fChain->SetBranchAddress("l1_labels", &l1_labels_, &b_l1_labels_);
   fChain->SetBranchAddress("l1_labels.first", &l1_labels_first, &b_l1_labels_first);
   fChain->SetBranchAddress("l1_labels.second", &l1_labels_second, &b_l1_labels_second);
   fChain->SetBranchAddress("jet_algoPF1_n", &jet_algoPF1_n, &b_jet_algoPF1_n);
   fChain->SetBranchAddress("jet_algoPF1_erescale", jet_algoPF1_erescale, &b_jet_algoPF1_erescale);
   fChain->SetBranchAddress("jet_algoPF1_p4", &jet_algoPF1_p4, &b_jet_algoPF1_p4);
   fChain->SetBranchAddress("vtx_std_diphopt", &vtx_std_diphopt, &b_vtx_std_diphopt);
   fChain->SetBranchAddress("vtx_std_nch", &vtx_std_nch, &b_vtx_std_nch);
   fChain->SetBranchAddress("vtx_std_ptmax", &vtx_std_ptmax, &b_vtx_std_ptmax);
   fChain->SetBranchAddress("vtx_std_sumpt", &vtx_std_sumpt, &b_vtx_std_sumpt);
   fChain->SetBranchAddress("vtx_std_ptvtx", &vtx_std_ptvtx, &b_vtx_std_ptvtx);
   fChain->SetBranchAddress("vtx_std_acosA", &vtx_std_acosA, &b_vtx_std_acosA);
   fChain->SetBranchAddress("vtx_std_ptasym", &vtx_std_ptasym, &b_vtx_std_ptasym);
   fChain->SetBranchAddress("vtx_std_ptbal", &vtx_std_ptbal, &b_vtx_std_ptbal);
   fChain->SetBranchAddress("vtx_std_nchthr", &vtx_std_nchthr, &b_vtx_std_nchthr);
   fChain->SetBranchAddress("vtx_std_ptmax3", &vtx_std_ptmax3, &b_vtx_std_ptmax3);
   fChain->SetBranchAddress("vtx_std_thrust", &vtx_std_thrust, &b_vtx_std_thrust);
   fChain->SetBranchAddress("vtx_std_sumweight", &vtx_std_sumweight, &b_vtx_std_sumweight);
   fChain->SetBranchAddress("vtx_std_sumpt2", &vtx_std_sumpt2, &b_vtx_std_sumpt2);
   fChain->SetBranchAddress("vtx_std_ptratio", &vtx_std_ptratio, &b_vtx_std_ptratio);
   fChain->SetBranchAddress("vtx_std_pzasym", &vtx_std_pzasym, &b_vtx_std_pzasym);
   fChain->SetBranchAddress("vtx_std_spher", &vtx_std_spher, &b_vtx_std_spher);
   fChain->SetBranchAddress("vtx_std_aplan", &vtx_std_aplan, &b_vtx_std_aplan);
   fChain->SetBranchAddress("vtx_std_sumpr", &vtx_std_sumpr, &b_vtx_std_sumpr);
   fChain->SetBranchAddress("vtx_std_sumawy", &vtx_std_sumawy, &b_vtx_std_sumawy);
   fChain->SetBranchAddress("vtx_std_sumtrv", &vtx_std_sumtrv, &b_vtx_std_sumtrv);
   fChain->SetBranchAddress("vtx_std_sumtwd", &vtx_std_sumtwd, &b_vtx_std_sumtwd);
   fChain->SetBranchAddress("vtx_std_awytwdasym", &vtx_std_awytwdasym, &b_vtx_std_awytwdasym);
   fChain->SetBranchAddress("vtx_std_pho1", &vtx_std_pho1, &b_vtx_std_pho1);
   fChain->SetBranchAddress("vtx_std_pho2", &vtx_std_pho2, &b_vtx_std_pho2);
   fChain->SetBranchAddress("pho_matchingConv", &pho_matchingConv, &b_pho_matchingConv);
   fChain->SetBranchAddress("vtx_std_ranked_list", &vtx_std_ranked_list, &b_vtx_std_ranked_list);
   fChain->SetBranchAddress("vtx_std_sel", &vtx_std_sel, &b_vtx_std_sel);
   fChain->SetBranchAddress("pho_tkiso_recvtx_030_002_0000_10_01", &pho_tkiso_recvtx_030_002_0000_10_01, &b_pho_tkiso_recvtx_030_002_0000_10_01);
   fChain->SetBranchAddress("pho_tkiso_badvtx_040_002_0000_10_01", pho_tkiso_badvtx_040_002_0000_10_01, &b_pho_tkiso_badvtx_040_002_0000_10_01);
   fChain->SetBranchAddress("pho_tkiso_badvtx_id", pho_tkiso_badvtx_id, &b_pho_tkiso_badvtx_id);
   fChain->SetBranchAddress("pho_drtotk_25_99", pho_drtotk_25_99, &b_pho_drtotk_25_99);
   fChain->SetBranchAddress("dipho_n", &dipho_n, &b_dipho_n);
   fChain->SetBranchAddress("dipho_leadind", dipho_leadind, &b_dipho_leadind);
   fChain->SetBranchAddress("dipho_subleadind", dipho_subleadind, &b_dipho_subleadind);
   fChain->SetBranchAddress("dipho_vtxind", dipho_vtxind, &b_dipho_vtxind);
   fChain->SetBranchAddress("dipho_sumpt", dipho_sumpt, &b_dipho_sumpt);
   fChain->SetBranchAddress("pho_cic6cutlevel_lead", &pho_cic6cutlevel_lead, &b_pho_cic6cutlevel_lead);
   fChain->SetBranchAddress("pho_cic6passcuts_lead", &pho_cic6passcuts_lead, &b_pho_cic6passcuts_lead);
   fChain->SetBranchAddress("pho_cic6cutlevel_sublead", &pho_cic6cutlevel_sublead, &b_pho_cic6cutlevel_sublead);
   fChain->SetBranchAddress("pho_cic6passcuts_sublead", &pho_cic6passcuts_sublead, &b_pho_cic6passcuts_sublead);
   fChain->SetBranchAddress("pho_cic4cutlevel_lead", &pho_cic4cutlevel_lead, &b_pho_cic4cutlevel_lead);
   fChain->SetBranchAddress("pho_cic4passcuts_lead", &pho_cic4passcuts_lead, &b_pho_cic4passcuts_lead);
   fChain->SetBranchAddress("pho_cic4cutlevel_sublead", &pho_cic4cutlevel_sublead, &b_pho_cic4cutlevel_sublead);
   fChain->SetBranchAddress("pho_cic4passcuts_sublead", &pho_cic4passcuts_sublead, &b_pho_cic4passcuts_sublead);
   Notify();
}

Bool_t h2gglobeEventReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void h2gglobeEventReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t h2gglobeEventReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef h2gglobeEventReader_cc
