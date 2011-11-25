#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "plotUtils_ntu.h"

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <vector>
#include <iostream>
#include <string>
#include <fstream>


int main(int argc, char** argv)
{
  setTDRStyle();
    

  //std::string inputDir("/gwteraz/users/deguio/PLOTS/NTUPLES/");
  //std::string outputDir("/gwteraz/users/deguio/PLOTS/FIGURES_hadRecoil/");
  //std::string outputDir("/gwteraz/users/deguio/PLOTS/FIGURES_pfElectrons/");
  //std::string outputDir("/gwteraz/users/deguio/PLOTS/FIGURES_useFakeRate/");
  std::string outputDir("/gwteraz/users/deguio/PLOTS/FIGURES_test/");
  std::string inputDir("/gwteraz/users/deguio/PLOTS/TEST/");

  int step = 10;
  //float lumi = 204.16;//May10 ReReco
  //float lumi = 264.25;//27052011
  //float lumi = 348.95;//03062011
  //float lumi = 498.61;//10062011
  //float lumi = 715.04;//17062011
  //float lumi = 880.82;//24062011
  //float lumi = 976.20;//01072011
  //float lumi = 1091.40;//06072011
  //float lumi = 67.25;//PRv5 29072011
  //float lumi = 1466;//05Jul + 05Aug
  //float lumi = 2129;//  Cert_160404-173692_7TeV_PromptReco_Collisions11_JSON.txt
  //float lumi = 1466;//  05Jul_bis+05Aug_bis
  
  float lumi = 1132;//05Jul
  //float lumi = 1000;

  std::string drawMode = "eventsScaled";
  //std::string drawMode = "sameAreaStack";
  //std::string drawMode = "sameAreaNoStack";
  
  
  // draw plots
  drawTStack_ntu* stack = new drawTStack_ntu(inputDir, "config/crossSections_wPrime_Summer11_ntu.txt", "ntu_WprimeTreeAnalysis", outputDir);
  //drawTStack_ntu* stack = new drawTStack_ntu(inputDir, "config/crossSections_wPrime_Summer11_ntu_test.txt", "ntu_WprimeTreeAnalysis", outputDir);
    

  //=============
  //==== Eff ====  
  //=============  

  //ERRORE NEL NUMERO DI EVENTI SE USO CORREZIONI COME CUTS (fondi DD, prescale o altro)  -->> adattare plotUtils

  stack -> DrawEvents("events", lumi, step, true);
  stack -> DrawEvents("eventsScaled", lumi, step, true);
  stack -> DrawEvents("eventsScaledStack", lumi, step, true);
  stack -> DrawEvents("efficiencies", lumi, step, true);
  stack -> DrawEvents("efficienciesRelative", lumi, step, true);
  stack -> DrawEvents("significance", lumi, step, false);


  std::vector<std::string> variableNames;
  variableNames.push_back("");
  std::vector<std::string> variableNames2;
  variableNames2.push_back("");
  variableNames2.push_back("");
  
  std::vector<std::string>* cuts = new std::vector<std::string>;
  cuts->push_back("");
  std::vector<std::string>* cuts2 = new std::vector<std::string>;
  cuts2->push_back("");  
  cuts2->push_back("");
  
  std::string histoName;
  
  //===========================
  //==== Numbers for limit ====  
  //===========================
  variableNames.at(0) = "eleMet_mt";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> numbersForLimit(variableNames,lumi,step,cuts,25.,0.,2500.);

  //============================
  //==== interesting events ====
  //============================
  stack -> printMtAboveThr(step,400.);
   
//   //============
//   //==== PV ====
//   //============

//   variableNames.at(0) = "PV_d0";  
//   histoName    = "PV_d0";
//   stack -> SetXaxisRange(0., 0.2);
//   stack -> SetXaxisTitle("d0_{PV}");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);
  
//   variableNames.at(0) = "PV_z";
//   histoName    = "PV_z";
//   stack -> SetXaxisRange(-30., 30.);
//   stack -> SetXaxisTitle("z_{PV}");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);
  
//   variableNames.at(0) = "PV_ndof";
//   histoName    = "PV_ndof";
//   stack -> SetXaxisRange(0., 250);
//   stack -> SetXaxisTitle("ndof_{PV}");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 125, false);
  
//   variableNames.at(0) = "PV_normalizedChi2";
//   histoName    = "PV_normalizedChi2";
//   stack -> SetXaxisRange(0., 2.);
//   stack -> SetXaxisTitle("#chi^{2}/ndof_{PV}");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);



  //=============
  //==== MET ====  
  //=============  

//   variableNames.at(0) = "met.Et()";
//   histoName    = "met";
//   cuts->at(0) = "pho_weight*hltPrescale";
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("PFMet (GeV)");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  variableNames.at(0) = "met.Et()";
  histoName    = "met";
  cuts->at(0) = "pho_weight*hltPrescale*((ele_corr.Et()/met.Et()) > 0.4 && (ele_corr.Et()/met.Et()) < 1.5)";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetXaxisTitle("PFMet (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  variableNames.at(0) = "calomet.Et()";
  histoName    = "calomet";
  cuts->at(0) = "pho_weight*hltPrescale*((ele_corr.Et()/calomet.Et()) > 0.4 && (ele_corr.Et()/calomet.Et()) < 1.5)";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetXaxisTitle("CaloMet (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  variableNames.at(0) = "tcmet.Et()";
  histoName    = "tcmet";
  cuts->at(0) = "pho_weight*hltPrescale*((ele_corr.Et()/tcmet.Et()) > 0.4 && (ele_corr.Et()/tcmet.Et()) < 1.5)";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetXaxisTitle("TcMet (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  // variableNames.at(0) = "met.px()";
  // histoName    = "mex";
  // cuts->at(0) = "PURescaleFactor((PUit_n+PUoot_n)/3.)*pho_weight*hltPrescale";
  // stack -> SetXaxisRange(-500., 500.);
  // stack -> SetXaxisTitle("PFMex (GeV)");
  // stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  // variableNames.at(0) = "met.py()";
  // histoName    = "mey";
  // cuts->at(0) = "PURescaleFactor((PUit_n+PUoot_n)/3.)*pho_weight*hltPrescale";
  // stack -> SetXaxisRange(-500., 500.);
  // stack -> SetXaxisTitle("PFMey (GeV)");
  // stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

//   //==============
//   //==== nEle ====
//   //==============
//   histoNames.at(0) = "hNEle";
//   stack -> SetXaxisRange(0., 10.);
//   stack -> SetXaxisTitle("nEle");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 1., false);


//   //===============
//   //==== nJets ====  
//   //===============
//   histoNames.at(0) = "hNJets";
//   stack -> SetXaxisRange(0., 50.);
//   stack -> SetXaxisTitle("nJets");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 1., false);


  
  //================
  //==== EleEta ====
  //================

  variableNames.at(0) = "ele_corr.eta()";
  histoName    = "eleEta";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(-3., 3.);
  stack -> SetXaxisTitle("#eta^{ele}");
  //stack -> SetYaxisRange(0.001, 10000.);
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);


  //================
  //==== ElePhi ====
  //================

  variableNames.at(0) = "ele_corr.phi()";
  histoName    = "elePhi";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("#phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);

  variableNames.at(0) = "ele_corr.phi()";
  histoName    = "elePhi_EB";
  cuts->at(0) = "pho_weight*hltPrescale*(ele_isEB == 1)";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EB #phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);

  variableNames.at(0) = "ele_corr.phi()";
  histoName    = "elePhi_EE";
  cuts->at(0) = "pho_weight*hltPrescale*(ele_isEB == 0)";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EE #phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);



  //============
  //==== Et ====
  //============
  variableNames.at(0) = "ele_corr.Et()";
  histoName    = "eleEt";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetYaxisRange(0.003, 1000000.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  variableNames.at(0) = "ele_corr.Et()";
  histoName    = "eleEt_cumulative";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);


  //===================
  //==== EtOverMet ====  
  //===================
  //N-1 selection
  variableNames.at(0) = "ele_corr.Et()/met.Et()";
  histoName    = "etOverMet";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 10.);
  stack -> SetXaxisTitle("E_{T}^{ele}/E_{T}^{miss}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step-1, 200, false, cuts, false);

//   histoNames.at(0) = "hEtOverMetNonIso";
//   //stack -> SetXaxisRange(0., 10.);
//   stack -> SetXaxisRange(0., 50.);
//   stack -> SetXaxisTitle("Et over Met non iso");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 20., false);

//   histoNames.at(0) = "hEtOverMetNonIso_EB";
//   //stack -> SetXaxisRange(0., 10.);
//   stack -> SetXaxisRange(0., 50.);
//   stack -> SetXaxisTitle("Et over Met non iso EB");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 20., false);

//   histoNames.at(0) = "hEtOverMetNonIso_EE";
//   //stack -> SetXaxisRange(0., 10.);
//   stack -> SetXaxisRange(0., 50.);
//   stack -> SetXaxisTitle("Et over Met non iso EE");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 20., false);


  //====================
  //==== DPhiEleMet ====  
  //====================
  //N-1 selection
  variableNames.at(0) = "eleMet_Dphi";
  histoName    = "eleMet_Dphi";
  cuts->at(0) = "pho_weight*hltPrescale * ((ele_corr.Et()/met.Et()) > 0.4) * ((ele_corr.Et()/met.Et()) < 1.5)";
  stack -> SetXaxisRange(0., 3.15);
  stack -> SetXaxisTitle("#Delta#phi_{e E_{T}^{miss}}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step-2, 200, false, cuts, false);

//   histoNames.at(0) = "hDPhiEleMet";
//   stack -> SetXaxisRange(0., 3.15);
//   stack -> SetXaxisTitle("#Delta#phi_{e E_{T}^{miss}}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 10., false);

//   histoNames.at(0) = "hDPhiEleMet_EB";
//   stack -> SetXaxisRange(0., 3.15);
//   stack -> SetXaxisTitle("EB #Delta#phi_{e E_{T}^{miss}}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 10., false);

//   histoNames.at(0) = "hDPhiEleMet_EE";
//   stack -> SetXaxisRange(0., 3.15);
//   stack -> SetXaxisTitle("EE #Delta#phi_{e E_{T}^{miss}}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 10., false);




  //==============
  //==== hMt =====  FIXME AGGIUNGERE *hltPrescale NEL CASO DI FAKERATE
  //==============

//   variableNames.at(0) = "eleMet_mt";
//   histoName    = "mT_outoftimeRescale";
//   cuts->at(0) = "pho_weight*hltPrescale";
//   stack -> SetXaxisRange(0., 1500.);
//   stack -> SetYaxisRange(0.001, 1000000.);
//   stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

//   variableNames.at(0) = "eleMet_mt";
//   histoName    = "mT_outoftimeRescale_cumulative";
//   cuts->at(0) = "pho_weight*hltPrescale";
//   stack -> SetXaxisRange(0., 1500.);
//   stack -> SetYaxisRange(0.001, 1000000.);
//   stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
//   stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_intimeRescale";
  cuts->at(0) = "PURescaleFactor(PUit_n)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_intimeRescale_cumulative";
  cuts->at(0) = "PURescaleFactor(PUit_n)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_noReweight";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_noReweight_cumulative";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);


  //==============
  //==== hEoP =====
  //==============

  variableNames.at(0) = "ele_EOverP";
  histoName    = "EoP";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 100.);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("E/P");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts, false);

  variableNames.at(0) = "ele_EOverP";
  histoName    = "EoP_cumulative";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 100.);
  stack -> SetXaxisTitle("E/P");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts, true);

  //==============
  //==== hHoE ====
  //==============

  variableNames.at(0) = "ele_HOverE";
  histoName    = "HoE";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 0.1);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("H/E");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts, false);


  //==============
  //==== hDPhiIn =====
  //==============

  variableNames.at(0) = "ele_DphiIn";
  histoName    = "DphiIn";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(-0.1, 0.1);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("DphiIn");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 40, true, cuts, false);

  //==============
  //==== hDEtaIn =====
  //==============

  variableNames.at(0) = "ele_DetaIn";
  histoName    = "DetaIn";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(-0.01, 0.01);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("DetaIn");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 40, true, cuts, false);

  //=================
  //==== iso var ====
  //=================

  variableNames.at(0) = "ele_tkIso";
  histoName    = "tkIso";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 20.);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("tkIso");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 40, true, cuts, false);

  variableNames.at(0) = "ele_hadIso";
  histoName    = "hadIso";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 20.);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("hadIso");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 40, true, cuts, false);

  variableNames.at(0) = "ele_emIso";
  histoName    = "emIso";
  cuts->at(0) = "pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 20.);
  //stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("emIso");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 40, true, cuts, false);


  
 }
