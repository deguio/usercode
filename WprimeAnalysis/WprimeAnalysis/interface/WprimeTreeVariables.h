#include "treeReader.h"
#include "ConfigParser.h"
#include "ntpleUtils.h"
#include "stdHisto.h"

#include "TH1F.h"
#include "TProfile.h"
#include "TObject.h"
#include "TTree.h"



struct WprimeVariables
{
  // tree definition
  TFile* m_outputRootFile;
  TTree* m_reducedTree;
  
  
  // input parameters
  float mW;
  int totEvents;
  float crossSection;
  int dataFlag; 
  int runId; 
  int lumiId; 
  int eventId; 
  int hltPrescale;
  
  //PU
  int mc_PU_NumInteractions;

  
  // PV variables
  float PV_d0;
  float PV_z;
  int PV_nTracks;
  int PV_ndof;
  float PV_normalizedChi2;

  int selectIt_ele;
  
  ROOT::Math::XYZTVector ele;
  ROOT::Math::XYZTVector* p_ele;
  
  float pho_weight;
  float ele_eSeed;
  float ele_timeSeed;
  int ele_flagSeed;
  float ele_swissCrossSeed;
  float ele_eSC;
  float ele_e1x5;
  float ele_e2x5;
  float ele_e5x5;
  float ele_charge;
  float ele_dxy;
  float ele_dz;
  float ele_tkIso;
  float ele_emIso;
  float ele_hadIso;
  float ele_hadIso_d1;
  float ele_hadIso_d2;
  int ele_isEB;
  int ele_isEcalDriven;
  float ele_sigmaIetaIeta;
  float ele_DphiIn;
  float ele_DetaIn;
  float ele_HOverE;
  float ele_fbrem;
  float ele_EOverP;

  //pho variables
  int selectIt_pho;
  
  ROOT::Math::XYZTVector pho;
  ROOT::Math::XYZTVector* p_pho;
  
  
  // met variables
  ROOT::Math::XYZTVector met;
  ROOT::Math::XYZTVector* p_met;
  
  ROOT::Math::XYZTVector eleMet;
  float eleMet_mt;
  float eleMet_Dphi;
  float phoMet_mt;
  float phoMet_Dphi;
  
  
};



void InitializeTree(WprimeVariables&, const std::string&);
void SetBranchAddresses(WprimeVariables& vars, TTree* tree);

void ClearWprimeVariables(WprimeVariables&);
void DeleteWprimeVariables(WprimeVariables&);
TTree* CloneTree(WprimeVariables& vars);

void SetPUVariables(WprimeVariables& vars, treeReader& reader);
void SetPVVariables(WprimeVariables& vars, treeReader& reader);
void SetElectronVariables(WprimeVariables& vars, treeReader& reader);
void SetPhotonVariables(WprimeVariables& vars, treeReader& reader);
void SetMetVariables(WprimeVariables& vars, treeReader& reader);

