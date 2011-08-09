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
    

  std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/FAKERATE/evalFakeRate/");
  //std::string inputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/FAKERATE/useFakeRate/");
  std::string outputDir("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/PLOTS_fakeRate/");

  //int step = 10;//numerator & useFakeRate & contamination(senza applicare FR)
  int step = 8;//denominator
  float lumi = 201.38; //May10
  

  std::string drawMode = "eventsScaled";
  //std::string drawMode = "sameAreaStack";
  
  
  // draw plots
  drawTStack_ntu* stack = new drawTStack_ntu(inputDir, "config/crossSections_wPrime_Spring11_fakeRate.txt", "ntu_WprimeTreeAnalysis", outputDir);
    
  //=============
  //==== Eff ====  
  //=============  
  stack -> DrawEvents("events", lumi, step, true);
  stack -> DrawEvents("eventsScaled", lumi, step, true);
  stack -> DrawEvents("eventsScaledStack", lumi, step, true);
  stack -> DrawEvents("efficiencies", lumi, step, true);
  stack -> DrawEvents("efficienciesRelative", lumi, step, true);

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

  //=============
  //==== MET ====  
  //=============  
  variableNames.at(0) = "met.Et()";
  histoName    = "met";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(0., 100.);
  stack -> SetXaxisTitle("PFMet (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false, cuts);

  variableNames.at(0) = "met.px()";
  histoName    = "mex";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(-500., 500.);
  stack -> SetXaxisTitle("PFMex (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false, cuts);

  variableNames.at(0) = "met.py()";
  histoName    = "mey";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(-500., 500.);
  stack -> SetXaxisTitle("PFMey (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 100, false, cuts);

  //================
  //==== Pho Et ====
  //================
  variableNames.at(0) = "pho.Et()";
  histoName    = "phoEt";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*hltPrescale";
  stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisTitle("pho Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, false);

  variableNames.at(0) = "pho.Et()";
  histoName    = "phoEt_cumulative";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)*hltPrescale";
  stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisTitle("pho Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, true);

  //================
  //==== Ele Et ====  
  //================
  variableNames.at(0) = "ele.Et()";
  histoName    = "eleEt";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, false);

  variableNames.at(0) = "ele.Et()";
  histoName    = "eleEt_cumulative";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, true);

  //==============
  //==== hMt =====
  //==============
  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(0., 400.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, false);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_cumulative";
  cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, true);

  variableNames.at(0) = "eleMet_mt";
  histoName    = "mT_cumulative_noReweight";
  cuts->at(0) = "";
  stack -> SetXaxisRange(0., 1500.);
  stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, true);
  
  // //========================
  // //==== hMt fake rate =====
  // //========================
  // variableNames.at(0) = "sqrt( 2. * pho.pt() * met.pt() * ( 1 - cos(eleMet_Dphi) ) )";
  // histoName    = "mT";
  // cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  // stack -> SetXaxisRange(0., 1500.);
  // stack -> SetXaxisTitle("M_{T} (GeV/c^{2})");
  // stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, true, cuts, false);

  // //==========================
  // //==== Ele Et fake rate ====  
  // //==========================
  // variableNames.at(0) = "ele.Et()";
  // histoName    = "eleEt";
  // cuts->at(0) = "PURescaleFactor(mc_PU_NumInteractions)";
  // stack -> SetXaxisRange(0., 150.);
  // stack -> SetXaxisTitle("ele Et (GeV)");
  // stack -> Draw(variableNames, histoName, drawMode, lumi, step, 200, false, cuts, false);

}
