#ifndef WprimeTreeContent_h
#define WprimeTreeContent_h

#include "TChain.h" 

#define MAXVERTICES 25
#define MAXGENPARTICLES 10
#define MAXELECTRONS 10
#define MAXSUPERCLUSTERS 100
#define MAXMUONS 10
#define MAXJETS 100
#define MAXTECHL1BITS 64
#define MAXALGOL1BITS 128


struct WprimeTreeContent
{
  // Flags
  static bool vertexVariables;
  static bool genVariables;
  static bool electronVariables;
  static bool superclusterVariables;
  static bool metVariables;
  static bool jetVariables;
  static bool muonVariables;
  static bool HLTrigVariables;
  static bool L1TrigVariables;
  
  unsigned int BX;
  unsigned int lumiId;
  unsigned int runId;
  unsigned int eventId;
  unsigned int eventNaiveId;
  unsigned int hcalnoiseLoose;
  unsigned int hcalnoiseTight;
 

  // electron variables
  int nGenParticles;
  float mc_V_pdgId[MAXGENPARTICLES];
  float mc_V_charge[MAXGENPARTICLES];
  float mc_V_E[MAXGENPARTICLES];
  float mc_V_px[MAXGENPARTICLES];
  float mc_V_py[MAXGENPARTICLES];
  float mc_V_pz[MAXGENPARTICLES];

  float mcF1_fromV_pdgId[MAXGENPARTICLES];
  float mcF1_fromV_charge[MAXGENPARTICLES];
  float mcF1_fromV_E[MAXGENPARTICLES];
  float mcF1_fromV_px[MAXGENPARTICLES];
  float mcF1_fromV_py[MAXGENPARTICLES];
  float mcF1_fromV_pz[MAXGENPARTICLES];

  float mcF2_fromV_pdgId[MAXGENPARTICLES];
  float mcF2_fromV_charge[MAXGENPARTICLES];
  float mcF2_fromV_E[MAXGENPARTICLES];
  float mcF2_fromV_px[MAXGENPARTICLES];
  float mcF2_fromV_py[MAXGENPARTICLES];
  float mcF2_fromV_pz[MAXGENPARTICLES];

  int nVertices;
  int nTracksVertex[MAXVERTICES];

  int nElectrons;
  float elePx[MAXELECTRONS];
  float elePy[MAXELECTRONS];
  float elePz[MAXELECTRONS];
  float eleE[MAXELECTRONS];
  float eleEt[MAXELECTRONS];
  float eleEta[MAXELECTRONS];
  float elePhi[MAXELECTRONS];
  int   eleId[MAXELECTRONS];
  float eleSigmaIEtaIEta[MAXELECTRONS];
  float eleE1x5[MAXELECTRONS];
  float eleE2x5[MAXELECTRONS];
  float eleE5x5[MAXELECTRONS];
  float eleSeedSwissCross[MAXELECTRONS];
  int   eleCharge[MAXELECTRONS];

  float eleHoE[MAXELECTRONS];
  float eleTrkIso[MAXELECTRONS];
  float eleEcalIso[MAXELECTRONS];
  float eleHcalIsoD1[MAXELECTRONS];
  float eleHcalIsoD2[MAXELECTRONS];
  int eleIsEB[MAXELECTRONS];
  int eleIsEE[MAXELECTRONS];
  int eleIsGap[MAXELECTRONS];
  float eleSeedEnergy[MAXELECTRONS];
  float eleSeedTime[MAXELECTRONS];
  int ecalRecHitRecoFlag[MAXELECTRONS];

  float genelePt[MAXELECTRONS];
  float geneleEta[MAXELECTRONS];
  float genelePhi[MAXELECTRONS];

  //SUPERCLUSTER VARIABLES 
  int nSuperClusters;
  float scE[MAXSUPERCLUSTERS];
  float scEt[MAXSUPERCLUSTERS];
  float scEta[MAXSUPERCLUSTERS];
  float scPhi[MAXSUPERCLUSTERS];
  float scHoE[MAXSUPERCLUSTERS];


  // MET VARIABLES
  float caloMet;
  float caloMex;
  float caloMey;
  float caloMetPhi;

  float tcMet;
  float tcMex;
  float tcMey;
  float tcMetPhi;

  float pfMet;
  float pfMex;
  float pfMey;
  float pfMetPhi;
 

  // JET VARIABLES
  int nJets;
  float jetPx[MAXJETS];
  float jetPy[MAXJETS];
  float jetPz[MAXJETS];
  float jetPt[MAXJETS];
  float jetEta[MAXJETS];
  float jetPhi[MAXJETS];
  float jetBdiscHighEff[MAXJETS];
  float jetBdiscHighPur[MAXJETS];

  float genjetPt[MAXJETS];
  float genjetEta[MAXJETS];
  float genjetPhi[MAXJETS];


  // MUON VARIABLES
  int nMuons;
  float muonPx[MAXMUONS];
  float muonPy[MAXMUONS];
  float muonPz[MAXMUONS];
  float muonPt[MAXMUONS];
  float muonEta[MAXMUONS];
  float muonPhi[MAXMUONS];
  

  // HLT VARIABLES
  int HLT_Ele15_LW_L1R;
  int HLT_Photon10_L1R;
  int HLT_Photon15_L1R;
  int HLT_Photon20_L1R;

  // L1 trigger variables
  int techL1Bit[MAXTECHL1BITS];
  int algoL1Bit[MAXALGOL1BITS];


};







// ------------------------------------------------------------------------
//! branch addresses settings

void setBranchAddresses(TTree* chain, WprimeTreeContent& treeVars);






// ------------------------------------------------------------------------
//! create branches for a tree

void setBranches(TTree* chain, WprimeTreeContent& treeVars);






// ------------------------------------------------------------------------
//! initialize branches

void initializeBranches(TTree* chain, WprimeTreeContent& treeVars);



#endif
