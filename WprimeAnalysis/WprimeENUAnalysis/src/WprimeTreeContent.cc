#include "WprimeAnalysis/WprimeENUAnalysis/interface/WprimeTreeContent.h"

bool WprimeTreeContent::genVariables           = true;
bool WprimeTreeContent::vertexVariables        = true;
bool WprimeTreeContent::electronVariables      = true;
bool WprimeTreeContent::superclusterVariables  = false;
bool WprimeTreeContent::metVariables           = true;
bool WprimeTreeContent::jetVariables           = true;
bool WprimeTreeContent::muonVariables          = true;
bool WprimeTreeContent::HLTrigVariables        = false;
bool WprimeTreeContent::L1TrigVariables        = false;
  
void setBranchAddresses(TTree* chain, WprimeTreeContent& treeVars)
{
  chain -> SetBranchAddress("BX",            &treeVars.BX);
  chain -> SetBranchAddress("lumiId",        &treeVars.lumiId);
  chain -> SetBranchAddress("runId",         &treeVars.runId);
  chain -> SetBranchAddress("eventId",       &treeVars.eventId);
  chain -> SetBranchAddress("eventNaiveId",  &treeVars.eventNaiveId);
  chain -> SetBranchAddress("hcalnoiseLoose",&treeVars.hcalnoiseLoose);
  chain -> SetBranchAddress("hcalnoiseTight",&treeVars.hcalnoiseTight);

  //GEN VARIABLES
  if(WprimeTreeContent::genVariables) { 
    chain -> SetBranchAddress("nGenParticles",             &treeVars.nGenParticles);

    chain -> SetBranchAddress("mc_V_pdgId",                 treeVars.mc_V_pdgId);
    chain -> SetBranchAddress("mc_V_E",                     treeVars.mc_V_E);
    chain -> SetBranchAddress("mc_V_px",                    treeVars.mc_V_px);
    chain -> SetBranchAddress("mc_V_py",                    treeVars.mc_V_py);
    chain -> SetBranchAddress("mc_V_pz",                    treeVars.mc_V_pz);
    chain -> SetBranchAddress("mc_V_charge",                treeVars.mc_V_charge);

    chain -> SetBranchAddress("mcF1_fromV_E",               treeVars.mcF1_fromV_E);
    chain -> SetBranchAddress("mcF1_fromV_px",              treeVars.mcF1_fromV_px);
    chain -> SetBranchAddress("mcF1_fromV_py",              treeVars.mcF1_fromV_py);
    chain -> SetBranchAddress("mcF1_fromV_pz",              treeVars.mcF1_fromV_pz);
    chain -> SetBranchAddress("mcF1_fromV_charge",          treeVars.mcF1_fromV_charge);
    chain -> SetBranchAddress("mcF1_fromV_pdgId",           treeVars.mcF1_fromV_pdgId);

    chain -> SetBranchAddress("mcF2_fromV_E",               treeVars.mcF2_fromV_E);
    chain -> SetBranchAddress("mcF2_fromV_px",              treeVars.mcF2_fromV_px);
    chain -> SetBranchAddress("mcF2_fromV_py",              treeVars.mcF2_fromV_py);
    chain -> SetBranchAddress("mcF2_fromV_pz",              treeVars.mcF2_fromV_pz);
    chain -> SetBranchAddress("mcF2_fromV_charge",          treeVars.mcF2_fromV_charge);
    chain -> SetBranchAddress("mcF2_fromV_pdgId",           treeVars.mcF2_fromV_pdgId);
  }

  //VERTEX VARIABLES
  if(WprimeTreeContent::vertexVariables) {
    
    chain -> SetBranchAddress("nVertices",             &treeVars.nVertices);
    chain -> SetBranchAddress("nTracksVertex",          treeVars.nTracksVertex);
  }



  // ELECTRON VARIABLES  
  if(WprimeTreeContent::electronVariables) {
  
    chain -> SetBranchAddress("nElectrons",             &treeVars.nElectrons);
    chain -> SetBranchAddress("elePx",                   treeVars.elePx);
    chain -> SetBranchAddress("elePy",                   treeVars.elePy);
    chain -> SetBranchAddress("elePz",                   treeVars.elePz);
    chain -> SetBranchAddress("eleE",                    treeVars.eleE);
    chain -> SetBranchAddress("eleEt",                   treeVars.eleEt);
    chain -> SetBranchAddress("eleEta",                  treeVars.eleEta);
    chain -> SetBranchAddress("elePhi",                  treeVars.elePhi);
    chain -> SetBranchAddress("eleId",                   treeVars.eleId);
    chain -> SetBranchAddress("eleSigmaIEtaIEta",        treeVars.eleSigmaIEtaIEta);
    chain -> SetBranchAddress("eleE1x5",                 treeVars.eleE1x5);
    chain -> SetBranchAddress("eleE2x5",                 treeVars.eleE2x5);
    chain -> SetBranchAddress("eleE5x5",                 treeVars.eleE5x5);
    chain -> SetBranchAddress("eleSeedSwissCross",       treeVars.eleSeedSwissCross);
    chain -> SetBranchAddress("eleCharge",               treeVars.eleCharge);
    chain -> SetBranchAddress("eleHoE",                  treeVars.eleHoE);
    chain -> SetBranchAddress("eleTrkIso",               treeVars.eleTrkIso);
    chain -> SetBranchAddress("eleEcalIso",              treeVars.eleEcalIso);
    chain -> SetBranchAddress("eleHcalIsoD1",            treeVars.eleHcalIsoD1);
    chain -> SetBranchAddress("eleHcalIsoD2",            treeVars.eleHcalIsoD2);
    chain -> SetBranchAddress("eleIsEB",                 treeVars.eleIsEB);
    chain -> SetBranchAddress("eleIsEE",                 treeVars.eleIsEE);
    chain -> SetBranchAddress("eleIsGap",                treeVars.eleIsGap);
    chain -> SetBranchAddress("eleSeedEnergy",           treeVars.eleSeedEnergy);
    chain -> SetBranchAddress("eleSeedTime",             treeVars.eleSeedTime);
    chain -> SetBranchAddress("ecalRecHitRecoFlag",      treeVars.ecalRecHitRecoFlag);
    
    chain -> SetBranchAddress("genelePt",                treeVars.genelePt); 
    chain -> SetBranchAddress("geneleEta",               treeVars.geneleEta); 
    chain -> SetBranchAddress("genelePhi",               treeVars.genelePhi); 
    
  } // ELECTRON VARIABLES
  
 
  // SUPERCLUSTER VARIABLES  
  if(WprimeTreeContent::superclusterVariables) {
    
    chain -> SetBranchAddress("nSuperClusters",  &treeVars.nSuperClusters);
    chain -> SetBranchAddress("scE",              treeVars.scE);
    chain -> SetBranchAddress("scEt",             treeVars.scEt);
    chain -> SetBranchAddress("scEta",            treeVars.scEta);
    chain -> SetBranchAddress("scPhi",            treeVars.scPhi);
    chain -> SetBranchAddress("scHoE",            treeVars.scHoE);
  
  }


  // MET VARIABLES  
  if(WprimeTreeContent::metVariables)  {  

    chain -> SetBranchAddress("caloMet",                   &treeVars.caloMet);
    chain -> SetBranchAddress("caloMex",                   &treeVars.caloMex);
    chain -> SetBranchAddress("caloMey",                   &treeVars.caloMey);
    chain -> SetBranchAddress("caloMetPhi",                &treeVars.caloMetPhi);
    
    chain -> SetBranchAddress("tcMet",                   &treeVars.tcMet);
    chain -> SetBranchAddress("tcMex",                   &treeVars.tcMex);
    chain -> SetBranchAddress("tcMey",                   &treeVars.tcMey);
    chain -> SetBranchAddress("tcMetPhi",                &treeVars.tcMetPhi);
    
    chain -> SetBranchAddress("pfMet",                   &treeVars.pfMet);
    chain -> SetBranchAddress("pfMex",                   &treeVars.pfMex);
    chain -> SetBranchAddress("pfMey",                   &treeVars.pfMey);
    chain -> SetBranchAddress("pfMetPhi",                &treeVars.pfMetPhi);
    
  } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables)  {  
    
    chain -> SetBranchAddress("nJets",             &treeVars.nJets);
    chain -> SetBranchAddress("jetPx",              treeVars.jetPx);
    chain -> SetBranchAddress("jetPy",              treeVars.jetPy);
    chain -> SetBranchAddress("jetPz",              treeVars.jetPz);
    chain -> SetBranchAddress("jetPt",              treeVars.jetPt);
    chain -> SetBranchAddress("jetEta",             treeVars.jetEta);
    chain -> SetBranchAddress("jetPhi",             treeVars.jetPhi);
    chain -> SetBranchAddress("jetBdiscHighEff",    treeVars.jetBdiscHighEff);
    chain -> SetBranchAddress("jetBdiscHighPur",    treeVars.jetBdiscHighPur);
    
    chain -> SetBranchAddress("genjetPt",           treeVars.genjetPt);
    chain -> SetBranchAddress("genjetEta",          treeVars.genjetEta);
    chain -> SetBranchAddress("genjetPhi",          treeVars.genjetPhi);
    
  } // JET VARIABLES
  

  
  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables)  {  
  
    chain -> SetBranchAddress("nMuons",             &treeVars.nMuons);
    chain -> SetBranchAddress("muonPx",              treeVars.muonPx);
    chain -> SetBranchAddress("muonPy",              treeVars.muonPy);
    chain -> SetBranchAddress("muonPz",              treeVars.muonPz);
    chain -> SetBranchAddress("muonPt",              treeVars.muonPt);
    chain -> SetBranchAddress("muonEta",             treeVars.muonEta);
    chain -> SetBranchAddress("muonPhi",             treeVars.muonPhi);
  
  } // MUON VARIABLES
  


  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables){
    
    chain -> SetBranchAddress("techL1Bit",     treeVars.techL1Bit);
    chain -> SetBranchAddress("algoL1Bit",     treeVars.algoL1Bit);
  
  }  //L1 VARIABLES
  
  
  //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)  {

    chain -> SetBranchAddress("HLT_Ele15_LW_L1R",           &treeVars.HLT_Ele15_LW_L1R);
    chain -> SetBranchAddress("HLT_Photon10_L1R",           &treeVars.HLT_Photon10_L1R);
    chain -> SetBranchAddress("HLT_Photon15_L1R",           &treeVars.HLT_Photon15_L1R);
    chain -> SetBranchAddress("HLT_Photon20_L1R",           &treeVars.HLT_Photon20_L1R);
    
  }//HLT VARIABLES
  
}


 



void setBranches(TTree* chain, WprimeTreeContent& treeVars)
{
  chain -> Branch("BX",            &treeVars.BX,                       "BX/i");
  chain -> Branch("lumiId",        &treeVars.lumiId,               "lumiId/i");
  chain -> Branch("runId",         &treeVars.runId,                 "runId/i");
  chain -> Branch("eventId",       &treeVars.eventId,             "eventId/i");
  chain -> Branch("eventNaiveId",  &treeVars.eventNaiveId,   "eventNaiveId/i");
  chain -> Branch("hcalnoiseLoose",  &treeVars.hcalnoiseLoose,   "hcalnoiseLoose/i");
  chain -> Branch("hcalnoiseTight",  &treeVars.hcalnoiseTight,   "hcalnoiseTight/i");
  
  // GENPARTICLES  VARIABLES  
  if(WprimeTreeContent::genVariables)  {
    
    chain -> Branch("nGenParticles",             &treeVars.nGenParticles,      "nGenParticles/I");
    chain -> Branch("mc_V_pdgId",                 treeVars.mc_V_pdgId,         "mc_V_pdgId[nGenParticles]/F");
    chain -> Branch("mc_V_E",                     treeVars.mc_V_E,             "mc_V_E[nGenParticles]/F");
    chain -> Branch("mc_V_px",                    treeVars.mc_V_px,            "mc_V_px[nGenParticles]/F");
    chain -> Branch("mc_V_py",                    treeVars.mc_V_py,            "mc_V_py[nGenParticles]/F");
    chain -> Branch("mc_V_pz",                    treeVars.mc_V_pz,            "mc_V_pz[nGenParticles]/F");
    chain -> Branch("mc_V_charge",                treeVars.mc_V_charge,        "mc_V_charge[nGenParticles]/F");
    
    chain -> Branch("mcF1_fromV_E",               treeVars.mcF1_fromV_E,       "mcF1_fromV_E[nGenParticles]/F");
    chain -> Branch("mcF1_fromV_px",              treeVars.mcF1_fromV_px,      "mcF1_fromV_px[nGenParticles]/F");
    chain -> Branch("mcF1_fromV_py",              treeVars.mcF1_fromV_py,      "mcF1_fromV_py[nGenParticles]/F");
    chain -> Branch("mcF1_fromV_pz",              treeVars.mcF1_fromV_pz,      "mcF1_fromV_pz[nGenParticles]/F");
    chain -> Branch("mcF1_fromV_charge",          treeVars.mcF1_fromV_charge,  "mcF1_fromV_charge[nGenParticles]/F");
    chain -> Branch("mcF1_fromV_pdgId",           treeVars.mcF1_fromV_pdgId,   "mcF1_fromV_pdgId[nGenParticles]/F");
    
    chain -> Branch("mcF2_fromV_E",               treeVars.mcF2_fromV_E,       "mcF2_fromV_E[nGenParticles]/F");
    chain -> Branch("mcF2_fromV_px",              treeVars.mcF2_fromV_px,      "mcF2_fromV_px[nGenParticles]/F");
    chain -> Branch("mcF2_fromV_py",              treeVars.mcF2_fromV_py,      "mcF2_fromV_py[nGenParticles]/F");
    chain -> Branch("mcF2_fromV_pz",              treeVars.mcF2_fromV_pz,      "mcF2_fromV_pz[nGenParticles]/F");
    chain -> Branch("mcF2_fromV_charge",          treeVars.mcF2_fromV_charge,  "mcF2_fromV_charge[nGenParticles]/F");
    chain -> Branch("mcF2_fromV_pdgId",           treeVars.mcF2_fromV_pdgId,   "mcF2_fromV_pdgId[nGenParticles]/F");
		    
}

  // VERTEX  VARIABLES  
  if(WprimeTreeContent::vertexVariables)  {
    
    chain -> Branch("nVertices",         &treeVars.nVertices,       "nVertices/I");
    chain -> Branch("nTracksVertex",     treeVars.nTracksVertex,    "nTracksVertex[nVertices]/I");
  
  }



  // ELECTRON  VARIABLES  
  if(WprimeTreeContent::electronVariables)  {
    
    chain -> Branch("nElectrons",        &treeVars.nElectrons,       "nElectrons/I");
    chain -> Branch("elePx",              treeVars.elePx,            "elePx[nElectrons]/F");
    chain -> Branch("elePy",              treeVars.elePy,            "elePy[nElectrons]/F");
    chain -> Branch("elePz",              treeVars.elePz,            "elePz[nElectrons]/F");
    chain -> Branch("eleE",               treeVars.eleE,             "eleE[nElectrons]/F");
    chain -> Branch("eleEt",              treeVars.eleEt,            "eleEt[nElectrons]/F");
    chain -> Branch("eleEta",             treeVars.eleEta,           "eleEta[nElectrons]/F");
    chain -> Branch("elePhi",             treeVars.elePhi,           "elePhi[nElectrons]/F");
    chain -> Branch("eleId",              treeVars.eleId,            "eleId[nElectrons]/I");
    chain -> Branch("eleSigmaIEtaIEta",   treeVars.eleSigmaIEtaIEta, "eleSigmaIEtaIEta[nElectrons]/F");
    chain -> Branch("eleE1x5",            treeVars.eleE1x5,          "eleE1xE5[nElectrons]/F");
    chain -> Branch("eleE2x5",            treeVars.eleE2x5,          "eleE2xE5[nElectrons]/F");
    chain -> Branch("eleE5x5",            treeVars.eleE5x5,          "eleE5xE5[nElectrons]/F");
    chain -> Branch("eleSeedSwissCross",  treeVars.eleSeedSwissCross,"eleSeedSwisscross[nElectrons]/F");
    chain -> Branch("eleCharge",          treeVars.eleCharge,       "eleCharge[nElectrons]/I");
    chain -> Branch("eleHoE",             treeVars.eleHoE,          "eleHoE[nElectrons]/F");
    chain -> Branch("eleTrkIso",          treeVars.eleTrkIso,       "eleTrkIso[nElectrons]/F");
    chain -> Branch("eleEcalIso",         treeVars.eleEcalIso,      "eleEcalIso[nElectrons]/F");
    chain -> Branch("eleHcalIsoD1",       treeVars.eleHcalIsoD1,    "eleHcalIsoD1[nElectrons]/F");
    chain -> Branch("eleHcalIsoD2",       treeVars.eleHcalIsoD2,    "eleHcalIsoD2[nElectrons]/F");
    chain -> Branch("eleIsEB",            treeVars.eleIsEB,         "eleIsEB[nElectrons]/I" );
    chain -> Branch("eleIsEE",            treeVars.eleIsEE,         "eleIsEE[nElectrons]/I" );
    chain -> Branch("eleIsGap",           treeVars.eleIsGap,        "eleIsGap[nElectrons]/I" );
    chain -> Branch("eleSeedEnergy",      treeVars.eleSeedEnergy,   "eleSeedEnergy[nElectrons]/F" );
    chain -> Branch("eleSeedTime",        treeVars.eleSeedTime,     "eleSeedTime[nElectrons]/F" );
    chain -> Branch("ecalRecHitRecoFlag", treeVars.ecalRecHitRecoFlag,"ecalRecHitRecoFlag[nElectrons]/I");

    chain -> Branch("genelePt",           treeVars.genelePt,        "genelePt[nElectrons]/F");
    chain -> Branch("geneleEta",          treeVars.geneleEta,       "geneleEta[nElectrons]/F");
    chain -> Branch("genelePhi",          treeVars.genelePhi,       "genelePhi[nElectrons]/F");
    
  }
  
 
  // SUPERCLUSTER  VARIABLES  
  if(WprimeTreeContent::superclusterVariables)  {

    chain -> Branch("nSuperClusters",    &treeVars.nSuperClusters,       "nSuperClusters/I");
    chain -> Branch("scE",               treeVars.scE,             "scE[nSuperClusters]/F");
    chain -> Branch("scEt",              treeVars.scEt,            "scEt[nSuperClusters]/F");
    chain -> Branch("scEta",             treeVars.scEta,           "scEta[nSuperClusters]/F");
    chain -> Branch("scPhi",             treeVars.scPhi,           "scPhi[nSuperClusters]/F");
    chain -> Branch("scHoE",             treeVars.scHoE,           "scHoE[nSuperClusters]/F");
  

  }

  // MET VARIABLES  
  if(WprimeTreeContent::metVariables)  {
    
    chain -> Branch("caloMet",    &treeVars.caloMet,   "caloMet/F");
    chain -> Branch("caloMex",    &treeVars.caloMex,   "caloMex/F"); 
    chain -> Branch("caloMey",    &treeVars.caloMey,   "caloMey/F");
    chain -> Branch("caloMetPhi", &treeVars.caloMetPhi,"caloMetPhi/F");
    
    chain -> Branch("tcMet",      &treeVars.tcMet,     "tcMet/F");
    chain -> Branch("tcMex",      &treeVars.tcMex,     "tcMex/F"); 
    chain -> Branch("tcMey",      &treeVars.tcMey,     "tcMey/F");
    chain -> Branch("tcMetPhi",   &treeVars.tcMetPhi,  "tcMetPhi/F");
    
    chain -> Branch("pfMet",      &treeVars.pfMet,     "pfMet/F");
    chain -> Branch("pfMex",      &treeVars.pfMex,     "pfMex/F"); 
    chain -> Branch("pfMey",      &treeVars.pfMey,     "pfMey/F");
    chain -> Branch("pfMetPhi",   &treeVars.pfMetPhi,  "pfMetPhi/F");
    
  } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables) {  
      
    chain -> Branch("nJets",        &treeVars.nJets,       "nJets/I");
    chain -> Branch("jetPx",        treeVars.jetPx,        "jetPx[nJets]/F");
    chain -> Branch("jetPy",        treeVars.jetPy,        "jetPy[nJets]/F");
    chain -> Branch("jetPz",        treeVars.jetPz,        "jetPz[nJets]/F");
    chain -> Branch("jetPt",        treeVars.jetPt,        "jetPt[nJets]/F");
    chain -> Branch("jetEta",       treeVars.jetEta,       "jetEta[nJets]/F");
    chain -> Branch("jetPhi",       treeVars.jetPhi,       "jetPhi[nJets]/F");
    chain -> Branch("jetBdiscHighEff",     treeVars.jetBdiscHighEff,     "jetBdiscHighEff[nJets]/F");
    chain -> Branch("jetBdiscHighPur",     treeVars.jetBdiscHighPur,     "jetBdiscHighPur[nJets]/F");
    
    chain -> Branch("genjetPt",     treeVars.genjetPt,     "genjetPt[nJets]/F");
    chain -> Branch("genjetEta",    treeVars.genjetEta,    "genjetEta[nJets]/F");
    chain -> Branch("genjetPhi",    treeVars.genjetPhi,    "genjetPhi[nJets]/F");
    
  } // JET VARIABLES
  


  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables){
  
    chain -> Branch("nMuons",     &treeVars.nMuons,         "nMuons/I");
    chain -> Branch("muonPx",      treeVars.muonPx,         "muonPx[nMuons]/F");
    chain -> Branch("muonPy",      treeVars.muonPy,         "muonPy[nMuons]/F");
    chain -> Branch("muonPz",      treeVars.muonPz,         "muonPz[nMuons]/F");
    chain -> Branch("muonPt",      treeVars.muonPt,         "muonPt[nMuons]/F");
    chain -> Branch("muonEta",     treeVars.muonEta,        "muonEta[nMuons]/F");
    chain -> Branch("muonPhi",     treeVars.muonPhi,        "muonPhi[nMuons]/F");
    
  } // MUON VARIABLES



 

  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables)  {

    chain -> Branch("techL1Bit",     treeVars.techL1Bit,    "techL1Bit[64]/I");
    chain -> Branch("algoL1Bit",     treeVars.algoL1Bit,    "algoL1Bit[128]/I");
    
  }  //L1 VARIABLES
  
  


  //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)  {

    chain -> Branch("HLT_Ele15_LW_L1R", &treeVars.HLT_Ele15_LW_L1R, "HLT_Ele15_LW_L1R/I");
    chain -> Branch("HLT_Photon10_L1R", &treeVars.HLT_Photon10_L1R, "HLT_Photon10_L1R/I");
    chain -> Branch("HLT_Photon15_L1R", &treeVars.HLT_Photon15_L1R, "HLT_Photon15_L1R/I");
    chain -> Branch("HLT_Photon20_L1R", &treeVars.HLT_Photon20_L1R, "HLT_Photon20_L1R/I");
    
  }//HLT VARIABLES
  
}



void initializeBranches(TTree* chain, WprimeTreeContent& treeVars)
{
  treeVars.BX = 0;
  treeVars.lumiId = 0;
  treeVars.runId = 0;
  treeVars.eventId = 0; 
  treeVars.eventNaiveId = 0;
  treeVars.hcalnoiseLoose = 0;
  treeVars.hcalnoiseTight = 0;
  
  // VERTICES VARIABLES  
  if(WprimeTreeContent::vertexVariables) {    
    for(int i = 0; i < MAXVERTICES; ++i){
      
      treeVars.nTracksVertex[i] = 0;
    }
    treeVars.nVertices = 0;
  }
    
  // GENPARTICLES VARIABLES  
  if(WprimeTreeContent::genVariables) {    
    for(int i = 0; i < MAXGENPARTICLES; ++i){
      
      treeVars.mc_V_pdgId[i] = -9999;
      treeVars.mc_V_charge[i] = -9999;
      treeVars.mc_V_E[i] = -9999;
      treeVars.mc_V_px[i] = -9999;
      treeVars.mc_V_py[i] = -9999;
      treeVars.mc_V_pz[i] = -9999;

      treeVars.mcF1_fromV_pdgId[i] = -9999;
      treeVars.mcF1_fromV_charge[i] = -9999;
      treeVars.mcF1_fromV_E[i] = -9999;
      treeVars.mcF1_fromV_px[i] = -9999;
      treeVars.mcF1_fromV_py[i] = -9999;
      treeVars.mcF1_fromV_pz[i] = -9999;

      treeVars.mcF2_fromV_pdgId[i] = -9999;
      treeVars.mcF2_fromV_charge[i] = -9999;
      treeVars.mcF2_fromV_E[i] = -9999;
      treeVars.mcF2_fromV_px[i] = -9999;
      treeVars.mcF2_fromV_py[i] = -9999;
      treeVars.mcF2_fromV_pz[i] = -9999;

    }
    treeVars.nGenParticles = 0;
  }
    






  // ELECTRONS VARIABLES  
  if(WprimeTreeContent::electronVariables) {    
    for(int i = 0; i < MAXELECTRONS; ++i){
     
      treeVars.elePx[i] = -9999;
      treeVars.elePy[i] = -9999;
      treeVars.elePz[i] = -9999;
      treeVars.eleE[i] = -9999;
      treeVars.eleEt[i] = -9999;
      treeVars.eleEta[i] = -9999;
      treeVars.elePhi[i] = -9999;
      treeVars.eleId[i] = -9999;
      treeVars.eleSigmaIEtaIEta[i] = -9999;
      treeVars.eleE1x5[i] = -9999;
      treeVars.eleE2x5[i] = -9999;
      treeVars.eleE5x5[i] = -9999;
      treeVars.eleSeedSwissCross[i] = -9999;
      treeVars.eleCharge[i] = -9999;
      treeVars.eleHoE[i] = -9999;
      treeVars.eleTrkIso[i] = -9999;
      treeVars.eleEcalIso[i] = -9999;
      treeVars.eleHcalIsoD1[i] = -9999;
      treeVars.eleHcalIsoD2[i] = -9999;
      treeVars.eleIsEB[i] = -9999;
      treeVars.eleIsEE[i] = -9999;
      treeVars.eleIsGap[i] = -9999;
      treeVars.eleSeedEnergy[i] = -9999;
      treeVars.eleSeedTime[i] = -9999;
      treeVars.ecalRecHitRecoFlag[i] = -9999;
      
      treeVars.genelePt[i] = -9999;
      treeVars.geneleEta[i] = -9999;
      treeVars.genelePhi[i] = -9999;
    }
    
    treeVars.nElectrons = 0;
  } // ELECTRONS VARIABLES
  

// ELECTRONS VARIABLES  
  if(WprimeTreeContent::electronVariables) {    
    for(int i = 0; i < MAXSUPERCLUSTERS; ++i){
     
      treeVars.scE[i]   = -9999;
      treeVars.scEt[i]  = -9999;
      treeVars.scEta[i] = -9999;
      treeVars.scPhi[i] = -9999;
      treeVars.scHoE[i] = -9999;
   
    }
    treeVars.nSuperClusters = 0;
  }

 // MET VARIABLES  
  if(WprimeTreeContent::metVariables) {
  
    treeVars.caloMet = -9999;
    treeVars.caloMex = -9999;
    treeVars.caloMey = -9999;
    treeVars.caloMetPhi = -9999;
    
    treeVars.tcMet = -9999;
    treeVars.tcMex = -9999;
    treeVars.tcMey = -9999;
    treeVars.tcMetPhi = -9999;
    
    treeVars.pfMet = -9999;
    treeVars.pfMex = -9999;
    treeVars.pfMey = -9999;
    treeVars.pfMetPhi = -9999;
    
  } // MET VARIABLES
  

  // JET VARIABLES  
  if(WprimeTreeContent::jetVariables) {  
    for(int i = 0; i < MAXJETS; ++i) {

      treeVars.jetPx[i] = -9999;
      treeVars.jetPy[i] = -9999;
      treeVars.jetPz[i] = -9999;
      treeVars.jetPt[i] = -9999;
      treeVars.jetEta[i] = -9999;
      treeVars.jetPhi[i] = -9999;
      treeVars.jetBdiscHighEff[i] = -9999;
      treeVars.jetBdiscHighPur[i] = -9999;
      
      treeVars.genjetPt[i]  = -9999;
      treeVars.genjetEta[i] = -9999;
      treeVars.genjetPhi[i] = -9999;
    }
    
    treeVars.nJets = 0;
    
  } // JET VARIABLES
  


  // MUON VARIABLES  
  if(WprimeTreeContent::muonVariables) {  
    for(int i = 0; i < MAXMUONS; ++i) {

      treeVars.muonPx[i] = -9999;
      treeVars.muonPy[i] = -9999;
      treeVars.muonPz[i] = -9999;
      treeVars.muonPt[i] = -9999;
      treeVars.muonEta[i] = -9999;
      treeVars.muonPhi[i] = -9999;
    }
    
    treeVars.nMuons = 0;
    
    
  } // MUON VARIABLES
  

  
  //L1 VARIABLES
  if(WprimeTreeContent::L1TrigVariables) {
    for (int i = 0; i < 64 ; i++){
      treeVars.techL1Bit[i] = -9999;
    }
    
    for (int i = 0; i < 128 ; i++){
      treeVars.algoL1Bit[i] = -9999;
    }
  }  //L1 VARIABLES
  
  
  //HLT VARIABLES
  if(WprimeTreeContent::HLTrigVariables)  {
  
    treeVars.HLT_Ele15_LW_L1R = -9999;
    treeVars.HLT_Photon10_L1R = -9999; 
    treeVars.HLT_Photon15_L1R = -9999; 
    treeVars.HLT_Photon20_L1R = -9999; 
  
  }//HLT VARIABLES
  

}
