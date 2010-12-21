#include "PhysicsTools/NtupleUtils/interface/setTDRStyle.h"
//#include "ntpleUtils.h"
#include "PhysicsTools/NtupleUtils/interface/plotUtils.h"

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
    

  std::string inputDir("/media/amassiro/deguio/Wprime/PROVA/");
  std::string outputDir("/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_8_4_patch2/src/WprimeAnalysis/WprimeENUAnalysis/bin/PROVA/");

  int step = 12;  
  //float lumi = 14.63;
  
  //float lumi = 2.855837;
  //float lumi = 4.293942;
  //float lumi = 7.307259;
  
  //float lumi = 4.954393;
  float lumi = 36.1;

  // draw plots
  drawTStack* stack = new drawTStack(inputDir, "/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_8_4_patch2/src/WprimeAnalysis/WprimeENUAnalysis/bin/crossSections_wPrime_fall10.txt", "WPrimeAnalysisTree", outputDir);
  
  std::string histoName;
  std::vector<std::string> histoNames;
  histoNames.push_back("");
  std::vector<std::string> histoNames2;
  histoNames2.push_back("");
  histoNames2.push_back("");
  

  //==============
  //==== hMt =====
  //==============
  histoNames.at(0) = "hMt";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("transverse Mass (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  //==============
  //==== hEt =====
  //==============
  histoNames.at(0) = "hEt";
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  //=============
  //==== MET ====  
  //=============  
  histoNames.at(0) = "hMet";
  //stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("TcMet (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  //================
  //==== EleEta ====
  //================
  histoNames.at(0) = "hEleEta";
  stack -> SetXaxisRange(-3., 3.);
  stack -> SetXaxisTitle("ele Eta");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 50., false);

  //================
  //==== ElePhi ====
  //================
  histoNames.at(0) = "hElePhi";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("ElePhi");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 50., false);

  //===========================
  //==== cumulative plots =====
  //===========================
  histoNames.at(0) = "hMt_cumulative";
  stack -> SetXaxisRange(0., 400.);
  stack -> SetXaxisTitle("transverse Mass cumulative (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 1., true);

  
}
