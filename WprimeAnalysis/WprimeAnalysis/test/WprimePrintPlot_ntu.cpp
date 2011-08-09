#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "plotUtils_ntu.h"
#include "plotUtils.h"

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
    

  std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/NTUPLES_06072011/");
  //std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/NTUPLES_May10/"); //HEEP
  //std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/hltStudies/");
  //std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/NTUPLES_27052011/");   //HEEP
  //std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/MC_ntu/");
  //std::string outputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/PLOTS_fakeRate/");
  std::string outputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/PLOTS/");

  int step = 10;
  //float lumi = 204.16;//May10 ReReco
  //float lumi = 264.25;//27052011
  //float lumi = 348.95;//03062011
  //float lumi = 498.61;//10062011
  //float lumi = 715.04;//17062011
  //float lumi = 880.82;//24062011
  //float lumi = 976.20;//01072011
  float lumi = 1091.40;//06072011
  //float lumi = 67.25;//PRv5 29072011

  std::string drawMode = "eventsScaled";
  //std::string drawMode = "sameAreaStack";
  //std::string drawMode = "sameAreaNoStack";
  
  
  // draw plots
  drawTStack_ntu* stack = new drawTStack_ntu(inputDir, "config/crossSections_wPrime_Spring11_ntu.txt", "ntu_WprimeTreeAnalysis", outputDir);
    

  //=============
  //==== Eff ====  
  //=============  

  //ERRORE NEL NUMERO DI EVENTI SE USO CORREZIONI COME CUTS (fondi DD, prescale o altro)  -->> adattare plotUtils

  stack -> DrawEvents("events", lumi, step, true);
  stack -> DrawEvents("eventsScaled", lumi, step, true);
  stack -> DrawEvents("eventsScaledStack", lumi, step, true);
  stack -> DrawEvents("efficiencies", lumi, step, true);
  stack -> DrawEvents("efficienciesRelative", lumi, step, true);
  //stack -> DrawEvents("significance", lumi, step, false);


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
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
  stack -> numbersForLimit(variableNames,lumi,step,cuts,25.,0.,1500.);

  //============================
  //==== interesting events ====
  //============================
  stack -> printMtAboveThr(step,400.);
   
  //============
  //==== PV ====
  //============

  variableNames.at(0) = "PV_d0";  
  histoName    = "PV_d0";
  stack -> SetXaxisRange(0., 0.2);
  stack -> SetXaxisTitle("d0_{PV}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);
  
  variableNames.at(0) = "PV_z";
  histoName    = "PV_z";
  stack -> SetXaxisRange(-30., 30.);
  stack -> SetXaxisTitle("z_{PV}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);
  
  variableNames.at(0) = "PV_ndof";
  histoName    = "PV_ndof";
  stack -> SetXaxisRange(0., 250);
  stack -> SetXaxisTitle("ndof_{PV}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 125, false);
  
  variableNames.at(0) = "PV_normalizedChi2";
  histoName    = "PV_normalizedChi2";
  stack -> SetXaxisRange(0., 2.);
  stack -> SetXaxisTitle("#chi^{2}/ndof_{PV}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false);



  //=============
  //==== MET ====  
  //=============  

  variableNames.at(0) = "met.Et()";
  histoName    = "met";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetXaxisTitle("PFMet (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  variableNames.at(0) = "met.px()";
  histoName    = "mex";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight";
  stack -> SetXaxisRange(-500., 500.);
  stack -> SetXaxisTitle("PFMex (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

  variableNames.at(0) = "met.py()";
  histoName    = "mey";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight";
  stack -> SetXaxisRange(-500., 500.);
  stack -> SetXaxisTitle("PFMey (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, true, cuts);

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

  variableNames.at(0) = "ele.eta()";
  histoName    = "eleEta";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(-3., 3.);
  stack -> SetXaxisTitle("#eta^{ele}");
  stack -> SetYaxisRange(0.001, 10000.);
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);


  //================
  //==== ElePhi ====
  //================

  variableNames.at(0) = "ele.phi()";
  histoName    = "elePhi";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("#phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);

  variableNames.at(0) = "ele.phi()";
  histoName    = "elePhi_EB";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*(ele_isEB == 1)*pho_weight";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EB #phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);

  variableNames.at(0) = "ele.phi()";
  histoName    = "elePhi_EE";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*(ele_isEB == 0)*pho_weight";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EE #phi^{ele}");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 80, false, cuts);




  //------------------ old -------------------------

//   histoNames.at(0) = "hElePhi";
//   stack -> SetXaxisRange(-3.5, 3.5);
//   stack -> SetXaxisTitle("#phi^{ele}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 25., false);

//   histoNames.at(0) = "hElePhi_EB";
//   stack -> SetXaxisRange(-3.5, 3.5);
//   stack -> SetXaxisTitle("EB #phi^{ele}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 25., false);

//   histoNames.at(0) = "hElePhi_EE";
//   stack -> SetXaxisRange(-3.5, 3.5);
//   stack -> SetXaxisTitle("EE #phi^{ele}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 25., false);


  //============
  //==== Et ====
  //============
  variableNames.at(0) = "ele.Et()";
  histoName    = "eleEt";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetYaxisRange(0.003, 1000000.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  variableNames.at(0) = "ele.Et()";
  histoName    = "eleEt_cumulative";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 500.);
  stack -> SetYaxisRange(0.003, 1000000.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, true);


  // histoNames.at(0) = "hEt";
  // //stack -> SetXaxisRange(0., 200.);
  // stack -> SetXaxisRange(0., 500.);
  // stack -> SetXaxisTitle("ele Et (GeV)");
  // stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);

//   histoNames.at(0) = "hEt_EB";
//   //stack -> SetXaxisRange(0., 200.);
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("EB ele Et (GeV)");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);

//   histoNames.at(0) = "hEt_EE";
//   //stack -> SetXaxisRange(0., 200.);
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("EE ele Et (GeV)");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);


//   //===================
//   //==== EtOverMet ====  
//   //===================
//   histoNames.at(0) = "hEtOverMet";
//   stack -> SetXaxisRange(0., 10.);
//   //stack -> SetXaxisRange(0., 2.);
//   stack -> SetXaxisTitle("E_{T}^{ele}/E_{T}^{miss}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., false);

//   histoNames.at(0) = "hEtOverMet_EB";
//   stack -> SetXaxisRange(0., 10.);
//   //stack -> SetXaxisRange(0., 2.);
//   stack -> SetXaxisTitle("EB E_{T}^{ele}/E_{T}^{miss}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., false);

//   histoNames.at(0) = "hEtOverMet_EE";
//   stack -> SetXaxisRange(0., 10.);
//   //stack -> SetXaxisRange(0., 2.);
//   stack -> SetXaxisTitle("EE E_{T}^{ele}/E_{T}^{miss}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., false);


  //===================
  //==== EtOverMet ====  
  //===================
  //N-1 selection
  variableNames.at(0) = "ele.Et()/met.Et()";
  histoName    = "etOverMet";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
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
  cuts->at(0) = "(PURescaleFactor(mc_PU_NumInteractions) && (ele.Et()/met.Et() > 0.4) && (ele.Et()/met.Et() < 1.5))*pho_weight*hltPrescale";
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
  //==== hMt =====
  //==============

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_cumulative";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*pho_weight*hltPrescale";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetYaxisRange(0.001, 1000000.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_cumulative_noReweight";
  cuts->at(0) = "pho_weight";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, true);


//   histoNames.at(0) = "hMt";
//   stack -> SetXaxisRange(0., 1500.);
//   stack -> SetYaxisRange(0.001, 500000.);
//   stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);

//   histoNames.at(0) = "hMt_EB";
//   stack -> SetXaxisRange(0., 1000.);
//   stack -> SetXaxisTitle("EB M_{T} (GeV/c^{2})");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);

//   histoNames.at(0) = "hMt_EE";
//   stack -> SetXaxisRange(0., 1000.);
//   stack -> SetXaxisTitle("EE M_{T} (GeV/c^{2})");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, step, 5., true);

// //   histoNames.at(0) = "hMtNonIso";
// //   stack -> SetXaxisRange(0., 250.);
// //   stack -> SetXaxisTitle("transverse Mass non iso (GeV)");
// //   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 4., false);

// //   histoNames.at(0) = "hMtNonIso_EB";
// //   stack -> SetXaxisRange(0., 250.);
// //   stack -> SetXaxisTitle("transverse Mass non iso EB (GeV)");
// //   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 4., false);

// //   histoNames.at(0) = "hMtNonIso_EE";
// //   stack -> SetXaxisRange(0., 250.);
// //   stack -> SetXaxisTitle("transverse Mass non iso EE (GeV)");
// //   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 4., false);


//   //===========================
//   //==== cumulative plots =====
//   //===========================
//   histoNames.at(0) = "hMt_cumulative";
//   stack -> SetXaxisRange(0., 1000.);
//   stack -> SetXaxisTitle("M_{T} cumulative (GeV/c^{2})");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

//   histoNames.at(0) = "hEt_cumulative";
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("Et cumulative (GeV)");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

//   histoNames.at(0) = "hMt_cumulative_EB";
//   stack -> SetXaxisRange(0., 1000.);
//   stack -> SetXaxisTitle("M_{T} cumulative EB (GeV/c^{2}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

//   histoNames.at(0) = "hMt_cumulative_EE";
//   stack -> SetXaxisRange(0., 1000.);
//   stack -> SetXaxisTitle("M_{T} cumulative EE (GeV/c^{2}");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

//   histoNames.at(0) = "hEt_cumulative_EE";
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("Et cumulative EE (GeV)");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

//   histoNames.at(0) = "hEt_cumulative_EB";
//   stack -> SetXaxisRange(0., 500.);
//   stack -> SetXaxisTitle("Et cumulative EB (GeV)");
//   stack -> Draw(histoNames, drawMode.c_str(), lumi, 0, 1., true);

  
 }
