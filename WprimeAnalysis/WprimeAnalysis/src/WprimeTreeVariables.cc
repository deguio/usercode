#include "WprimeTreeVariables.h"



void InitializeTree(WprimeVariables& vars, const std::string& outputRootFileName)
{
  //-------------
  // Reduced tree
  //-------------
  
  vars.m_outputRootFile = new TFile(outputRootFileName.c_str(), "RECREATE");  
  
  vars.m_reducedTree = new TTree("ntu", "ntu");
  vars.m_reducedTree -> SetDirectory(vars.m_outputRootFile);
  
  vars.m_reducedTree -> Branch("mW",           &vars.mW,                     "mW/F");
  vars.m_reducedTree -> Branch("totEvents",    &vars.totEvents,       "totEvents/I");
  vars.m_reducedTree -> Branch("crossSection", &vars.crossSection, "crossSection/F");
  vars.m_reducedTree -> Branch("dataFlag",     &vars.dataFlag,         "dataFlag/I");
  vars.m_reducedTree -> Branch("runId",        &vars.runId,               "runId/I");
  vars.m_reducedTree -> Branch("lumiId",       &vars.lumiId,             "lumiId/I");
  vars.m_reducedTree -> Branch("eventId",      &vars.eventId,           "eventId/I");
  vars.m_reducedTree -> Branch("hltPrescale",  &vars.hltPrescale,       "hltPrescale/I");
  
  //PU variables
  vars.m_reducedTree -> Branch("mc_PU_NumInteractions",             &vars.mc_PU_NumInteractions,                         "mc_PU_NumInteractions/I");  

  // PV variables
  vars.m_reducedTree -> Branch("PV_d0",             &vars.PV_d0,                         "PV_d0/F");
  vars.m_reducedTree -> Branch("PV_z",              &vars.PV_z,                           "PV_z/F");
  vars.m_reducedTree -> Branch("PV_nTracks",        &vars.PV_nTracks,               "PV_nTracks/I");
  vars.m_reducedTree -> Branch("PV_ndof",           &vars.PV_ndof,                     "PV_ndof/I");
  vars.m_reducedTree -> Branch("PV_normalizedChi2", &vars.PV_normalizedChi2, "PV_normalizedChi2/F");
  
  
  // lepton variables
  vars.m_reducedTree -> Branch("ele", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &vars.p_ele);
  vars.m_reducedTree -> Branch("ele_eSC",  &vars.ele_eSC,   "ele_eSC/F");
  vars.m_reducedTree -> Branch("ele_eSeed",  &vars.ele_eSeed,   "ele_eSeed/F");
  vars.m_reducedTree -> Branch("ele_timeSeed",  &vars.ele_timeSeed,   "ele_timeSeed/F");
  vars.m_reducedTree -> Branch("ele_flagSeed",  &vars.ele_flagSeed,   "ele_flagSeed/I");
  vars.m_reducedTree -> Branch("ele_swissCrossSeed",  &vars.ele_swissCrossSeed,   "ele_swissCrossSeed/F");
  vars.m_reducedTree -> Branch("ele_e1x5",  &vars.ele_e1x5,   "ele_ex5/F");
  vars.m_reducedTree -> Branch("ele_e2x5",  &vars.ele_e2x5,   "ele_e2x5/F");
  vars.m_reducedTree -> Branch("ele_e5x5",  &vars.ele_e5x5,   "ele_e5x5/F");
  vars.m_reducedTree -> Branch("ele_charge",  &vars.ele_charge,   "ele_charge/F");
  vars.m_reducedTree -> Branch("ele_dxy",     &vars.ele_dxy,         "ele_dxy/F");
  vars.m_reducedTree -> Branch("ele_dz",      &vars.ele_dz,           "ele_dz/F");
  vars.m_reducedTree -> Branch("ele_tkIso",   &vars.ele_tkIso,     "ele_tkIso/F");
  vars.m_reducedTree -> Branch("ele_emIso",   &vars.ele_emIso,     "ele_emIso/F");
  vars.m_reducedTree -> Branch("ele_hadIso",  &vars.ele_hadIso,   "ele_hadIso/F");
  vars.m_reducedTree -> Branch("ele_hadIso_d1",  &vars.ele_hadIso_d1,   "ele_hadIso_d1/F");
  vars.m_reducedTree -> Branch("ele_hadIso_d2",  &vars.ele_hadIso_d2,   "ele_hadIso_d2/F");
  vars.m_reducedTree -> Branch("ele_isEB",          &vars.ele_isEB,                   "ele_isEB/I");
  vars.m_reducedTree -> Branch("ele_isEcalDriven",          &vars.ele_isEcalDriven,       "ele_isEcalDriven/I");
  vars.m_reducedTree -> Branch("ele_sigmaIetaIeta", &vars.ele_sigmaIetaIeta, "ele_sigmaIetaIeta/F");
  vars.m_reducedTree -> Branch("ele_DphiIn",        &vars.ele_DphiIn,               "ele_DphiIn/F");
  vars.m_reducedTree -> Branch("ele_DetaIn",        &vars.ele_DetaIn,               "ele_DetaIn/F");
  vars.m_reducedTree -> Branch("ele_HOverE",        &vars.ele_HOverE,               "ele_HOverE/F");
  vars.m_reducedTree -> Branch("ele_fbrem",         &vars.ele_fbrem,                 "ele_fbrem/F");
  vars.m_reducedTree -> Branch("ele_EOverP",        &vars.ele_EOverP,               "ele_EOverP/F");
  
  //photons variables
  vars.m_reducedTree -> Branch("pho", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &vars.p_pho);
  vars.m_reducedTree -> Branch("pho_weight",  &vars.pho_weight,   "pho_weight/F");

  // met variables
  vars.m_reducedTree -> Branch("met", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &vars.p_met);
  vars.m_reducedTree -> Branch("eleMet_mt",    &vars.eleMet_mt,       "eleMet_mt/F");
  vars.m_reducedTree -> Branch("eleMet_Dphi",  &vars.eleMet_Dphi,   "eleMet_Dphi/F");  
  vars.m_reducedTree -> Branch("phoMet_mt",    &vars.phoMet_mt,       "phoMet_mt/F");
  vars.m_reducedTree -> Branch("phoMet_Dphi",  &vars.phoMet_Dphi,   "phoMet_Dphi/F");  
}

void SetBranchAddresses(WprimeVariables& vars, TTree* tree)
{
  tree -> SetBranchAddress("runId",        &vars.runId);
  tree -> SetBranchAddress("lumiId",       &vars.lumiId);
  tree -> SetBranchAddress("eventId",      &vars.eventId);
  tree -> SetBranchAddress("eleMet_mt",    &vars.eleMet_mt);

  vars.p_ele = &vars.ele;
  tree -> SetBranchAddress("ele",          &vars.p_ele);
  tree -> SetBranchAddress("ele_eSC",      &vars.ele_eSC);
  tree -> SetBranchAddress("ele_eSeed",    &vars.ele_eSeed);
  tree -> SetBranchAddress("ele_timeSeed", &vars.ele_timeSeed);
  tree -> SetBranchAddress("ele_flagSeed", &vars.ele_flagSeed);
  tree -> SetBranchAddress("ele_swissCrossSeed", &vars.ele_swissCrossSeed);
  tree -> SetBranchAddress("ele_e1x5",     &vars.ele_e1x5);
  tree -> SetBranchAddress("ele_e2x5",     &vars.ele_e2x5);
  tree -> SetBranchAddress("ele_e5x5",     &vars.ele_e5x5);
  tree -> SetBranchAddress("ele_charge",   &vars.ele_charge);
  tree -> SetBranchAddress("ele_dxy",      &vars.ele_dxy);
  tree -> SetBranchAddress("ele_dz",       &vars.ele_dz);
  tree -> SetBranchAddress("ele_tkIso",    &vars.ele_tkIso);
  tree -> SetBranchAddress("ele_emIso",    &vars.ele_emIso);
  tree -> SetBranchAddress("ele_hadIso",   &vars.ele_hadIso);
  tree -> SetBranchAddress("ele_hadIso_d1",  &vars.ele_hadIso_d1);
  tree -> SetBranchAddress("ele_hadIso_d2",  &vars.ele_hadIso_d2);
  tree -> SetBranchAddress("ele_isEB",          &vars.ele_isEB);
  tree -> SetBranchAddress("ele_isEcalDriven",  &vars.ele_isEcalDriven);
  tree -> SetBranchAddress("ele_sigmaIetaIeta", &vars.ele_sigmaIetaIeta);
  tree -> SetBranchAddress("ele_DphiIn",        &vars.ele_DphiIn);
  tree -> SetBranchAddress("ele_DetaIn",        &vars.ele_DetaIn);
  tree -> SetBranchAddress("ele_HOverE",        &vars.ele_HOverE);
  tree -> SetBranchAddress("ele_fbrem",         &vars.ele_fbrem);
  tree -> SetBranchAddress("ele_EOverP",        &vars.ele_EOverP);
  
  //photons variables
  vars.p_pho = &vars.pho;
  tree -> SetBranchAddress("pho",               &vars.p_pho);
  tree -> SetBranchAddress("pho_weight",        &vars.pho_weight);

  // met variables
  vars.p_met = &vars.met;
  tree -> SetBranchAddress("met",          &vars.p_met);
  tree -> SetBranchAddress("eleMet_mt",    &vars.eleMet_mt);
  tree -> SetBranchAddress("eleMet_Dphi",  &vars.eleMet_Dphi);  
  tree -> SetBranchAddress("phoMet_mt",    &vars.phoMet_mt);
  tree -> SetBranchAddress("phoMet_Dphi",  &vars.phoMet_Dphi);  
}


void ClearWprimeVariables(WprimeVariables& vars)
{
  //run variables
  vars.runId = -1;
  vars.lumiId = -1;
  vars.eventId = -1;
  vars.hltPrescale = -1;

  //PU variables
  vars.mc_PU_NumInteractions = -1;

  // PV variables
  vars.PV_d0 = -1.;
  vars.PV_z = -99.;  
  vars.PV_nTracks = -1;
  vars.PV_ndof = -1;
  vars.PV_normalizedChi2 = -1.;
  
  
  // eleton variables
  vars.selectIt_ele = -1;

  vars.ele = ROOT::Math::XYZTVector(0., 0., 0., 0.);
  vars.p_ele = NULL;
  
  vars.pho_weight = 1.;
  vars.ele_eSC = -1.;
  vars.ele_eSeed = -1.;
  vars.ele_timeSeed = -99;
  vars.ele_flagSeed = -99;
  vars.ele_swissCrossSeed = -99;
  vars.ele_e1x5 = -1.;
  vars.ele_e2x5 = -1.;
  vars.ele_e5x5 = -1.;
  vars.ele_charge = -1.;
  vars.ele_dxy = -1.;
  vars.ele_dz = -99.;
  vars.ele_tkIso = -1.;
  vars.ele_emIso = -1.;
  vars.ele_hadIso = -1.;
  vars.ele_hadIso_d1 = -1.;
  vars.ele_hadIso_d2 = -1.;
  vars.ele_isEB = -1;
  vars.ele_isEcalDriven = -1;
  vars.ele_sigmaIetaIeta = -1.;
  vars.ele_DphiIn = -99.;
  vars.ele_DetaIn = -99.;
  vars.ele_HOverE = -1.;
  vars.ele_fbrem = -1.;
  vars.ele_EOverP = -1.;

  // photon variables
  vars.selectIt_pho = -1;

  vars.pho = ROOT::Math::XYZTVector(0., 0., 0., 0.);
  vars.p_pho = NULL;

  // met variables 
  vars.met = ROOT::Math::XYZTVector(0., 0., 0., 0.);
  vars.p_met = NULL;
  
  vars.eleMet = ROOT::Math::XYZTVector(0., 0., 0., 0.);
  
  vars.eleMet_mt = -1.;
  vars.eleMet_Dphi = -1.;
  vars.phoMet_mt = -1.;
  vars.phoMet_Dphi = -1.;
  
 
}



void DeleteWprimeVariables(WprimeVariables& vars)
{
  // save tree
  vars.m_outputRootFile -> cd();
  //vars.m_reducedTree -> Write();
  vars.m_outputRootFile -> Close();
}

TTree* CloneTree(WprimeVariables& vars)
{
  //clone tree
  return vars.m_reducedTree -> CloneTree(0);
}

void SetPUVariables(WprimeVariables& vars, treeReader& reader)
{
  vars.mc_PU_NumInteractions = reader.GetInt("mc_PUit_NumInteractions")->at(0);
}

void SetPVVariables(WprimeVariables& vars, treeReader& reader)
{
  vars.PV_d0 = reader.GetFloat("PV_d0")->at(0);
  vars.PV_nTracks = reader.GetInt("PV_nTracks")->at(0);
  vars.PV_ndof = reader.GetInt("PV_ndof")->at(0);
  vars.PV_normalizedChi2 = reader.GetFloat("PV_normalizedChi2")->at(0);
  vars.PV_z = reader.GetFloat("PV_z")->at(0);
}


void SetElectronVariables(WprimeVariables& vars, treeReader& reader)
{
  vars.ele = reader.Get4V("electrons")->at(vars.selectIt_ele);
  vars.p_ele = &vars.ele;

  vars.ele_eSC  = reader.GetFloat("electrons_eSC")->at(vars.selectIt_ele);
  vars.ele_eSeed  = reader.GetFloat("electrons_eSeed")->at(vars.selectIt_ele);
  vars.ele_timeSeed  = reader.GetFloat("electrons_timeSeed")->at(vars.selectIt_ele);
  vars.ele_flagSeed  = reader.GetInt("electrons_flagSeed")->at(vars.selectIt_ele);
  vars.ele_swissCrossSeed  = reader.GetFloat("electrons_swissCrossSeed")->at(vars.selectIt_ele);
  vars.ele_e1x5  = reader.GetFloat("electrons_e1x5")->at(vars.selectIt_ele);
  vars.ele_e2x5  = reader.GetFloat("electrons_e2x5Max")->at(vars.selectIt_ele);
  vars.ele_e5x5  = reader.GetFloat("electrons_e5x5")->at(vars.selectIt_ele);
  vars.ele_charge  = reader.GetFloat("electrons_charge")->at(vars.selectIt_ele);
  vars.ele_dxy     = reader.GetFloat("electrons_dxy_PV")->at(vars.selectIt_ele);
  vars.ele_dz      = reader.GetFloat("electrons_dz_PV")->at(vars.selectIt_ele);
  vars.ele_tkIso   = reader.GetFloat("electrons_tkIsoR03")->at(vars.selectIt_ele);
  vars.ele_emIso   = reader.GetFloat("electrons_emIsoR03")->at(vars.selectIt_ele);
  vars.ele_hadIso  = reader.GetFloat("electrons_hadIsoR03_depth1")->at(vars.selectIt_ele)+reader.GetFloat("electrons_hadIsoR03_depth2")->at(vars.selectIt_ele);
  vars.ele_hadIso_d1  = reader.GetFloat("electrons_hadIsoR03_depth1")->at(vars.selectIt_ele);
  vars.ele_hadIso_d2  = reader.GetFloat("electrons_hadIsoR03_depth2")->at(vars.selectIt_ele);
  vars.ele_isEB = reader.GetInt("electrons_isEB")->at(vars.selectIt_ele);
  vars.ele_isEcalDriven = reader.GetInt("electrons_ecalDrivenSeed")->at(vars.selectIt_ele);
  vars.ele_sigmaIetaIeta = reader.GetFloat("electrons_sigmaIetaIeta")->at(vars.selectIt_ele);
  vars.ele_DphiIn = reader.GetFloat("electrons_deltaPhiIn")->at(vars.selectIt_ele);
  vars.ele_DetaIn = reader.GetFloat("electrons_deltaEtaIn")->at(vars.selectIt_ele);
  vars.ele_HOverE = reader.GetFloat("electrons_hOverE")->at(vars.selectIt_ele);
  vars.ele_fbrem  = reader.GetFloat("electrons_fbrem")->at(vars.selectIt_ele);
  vars.ele_EOverP = reader.GetFloat("electrons_eSCOverP")->at(vars.selectIt_ele);
 
}

void SetPhotonVariables(WprimeVariables& vars, treeReader& reader)
{
  vars.pho = reader.Get4V("photons")->at(vars.selectIt_pho);
  vars.p_pho = &vars.pho;
}


void SetMetVariables(WprimeVariables& vars, treeReader& reader)
{
  vars.met = reader.Get4V("PFMet")->at(0);
  vars.p_met = &vars.met;
  
  vars.eleMet = vars.ele + vars.met;
  
  vars.eleMet_Dphi = deltaPhi(vars.ele.phi(), vars.met.phi());
  vars.eleMet_mt = sqrt( 2. * vars.ele.pt() * vars.met.pt() * ( 1 - cos(vars.eleMet_Dphi) ) ) ;

  vars.phoMet_Dphi = deltaPhi(vars.pho.phi(), vars.met.phi());
  vars.phoMet_mt = sqrt( 2. * vars.pho.pt() * vars.met.pt() * ( 1 - cos(vars.phoMet_Dphi) ) ) ;
}
