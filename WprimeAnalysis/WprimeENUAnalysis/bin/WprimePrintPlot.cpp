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
    

  std::string inputDir("/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_6_3/src/WprimeAnalysis/WprimeENUAnalysis/bin/PROVA/");
  std::string outputDir("/afs/cern.ch/user/d/deguio/scratch0/Wprime/CMSSW_3_6_3/src/WprimeAnalysis/WprimeENUAnalysis/bin/PROVA/");

  int step = 9;
  float lumi = 3.00269;  // pb-1 -> eventScaled
  
  
  
  
  // draw plots
  drawTStack* stack = new drawTStack(inputDir, "crossSections_wPrime.txt", "WPrimeAnalysisTree", outputDir);
  
  std::string histoName;
  std::vector<std::string> histoNames;
  histoNames.push_back("");
  std::vector<std::string> histoNames2;
  histoNames2.push_back("");
  histoNames2.push_back("");
  
  //=============
  //==== Eff ====  
  //=============  
  stack -> DrawEvents("events", lumi, step, true);
  stack -> DrawEvents("eventsScaled", lumi, step, true);
  stack -> DrawEvents("eventsScaledStack", lumi, step, true);
  stack -> DrawEvents("efficiencies", lumi, step, true);
  stack -> DrawEvents("efficienciesRelative", lumi, step, true);

  
  //=============
  //==== MET ====  
  //=============  
  histoNames.at(0) = "hMet";
  //stack -> SetXaxisRange(0., 150.);
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("TcMet (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., false);
  

  //==============
  //==== nEle ====
  //==============
  histoNames.at(0) = "hNEle";
  stack -> SetXaxisRange(0., 10.);
  stack -> SetXaxisTitle("nEle");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 1., false);


  //===============
  //==== nJets ====  
  //===============
  histoNames.at(0) = "hNJets";
  stack -> SetXaxisRange(0., 50.);
  stack -> SetXaxisTitle("nJets");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 1., false);


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

  histoNames.at(0) = "hElePhi_EB";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EB ElePhi");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 50., false);

  histoNames.at(0) = "hElePhi_EE";
  stack -> SetXaxisRange(-3.5, 3.5);
  stack -> SetXaxisTitle("EE ElePhi");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 50., false);


  //============
  //==== Et ====  
  //============
  histoNames.at(0) = "hEt";
  //stack -> SetXaxisRange(0., 200.);
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("ele Et (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., false);

  histoNames.at(0) = "hEt_EB";
  //stack -> SetXaxisRange(0., 200.);
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("EB ele Et (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., false);

  histoNames.at(0) = "hEt_EE";
  //stack -> SetXaxisRange(0., 200.);
  stack -> SetXaxisRange(0., 120.);
  stack -> SetXaxisTitle("EE ele Et (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., false);


  //===================
  //==== EtOverMet ====  
  //===================
  histoNames.at(0) = "hEtOverMet";
  stack -> SetXaxisRange(0., 10.);
  //stack -> SetXaxisRange(0., 2.);
  stack -> SetXaxisTitle("Et over Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 5., false);

  histoNames.at(0) = "hEtOverMet_EB";
  stack -> SetXaxisRange(0., 10.);
  //stack -> SetXaxisRange(0., 2.);
  stack -> SetXaxisTitle("EB Et over Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 5., false);

  histoNames.at(0) = "hEtOverMet_EE";
  stack -> SetXaxisRange(0., 10.);
  //stack -> SetXaxisRange(0., 2.);
  stack -> SetXaxisTitle("EE Et over Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 5., false);


  //==========================
  //==== EtOverMet NonIso ====  
  //==========================
  histoNames.at(0) = "hEtOverMetNonIso";
  //stack -> SetXaxisRange(0., 10.);
  stack -> SetXaxisRange(0., 50.);
  stack -> SetXaxisTitle("Et over Met non iso");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 20., false);

  histoNames.at(0) = "hEtOverMetNonIso_EB";
  //stack -> SetXaxisRange(0., 10.);
  stack -> SetXaxisRange(0., 50.);
  stack -> SetXaxisTitle("Et over Met non iso EB");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 20., false);

  histoNames.at(0) = "hEtOverMetNonIso_EE";
  //stack -> SetXaxisRange(0., 10.);
  stack -> SetXaxisRange(0., 50.);
  stack -> SetXaxisTitle("Et over Met non iso EE");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 20., false);


  //====================
  //==== DPhiEleMet ====  
  //====================
  histoNames.at(0) = "hDPhiEleMet";
  stack -> SetXaxisRange(2., 4.);
  stack -> SetXaxisTitle("DPhi Ele Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 10., false);

  histoNames.at(0) = "hDPhiEleMet_EB";
  stack -> SetXaxisRange(2., 4.);
  stack -> SetXaxisTitle("EB DPhi Ele Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 10., false);

  histoNames.at(0) = "hDPhiEleMet_EE";
  stack -> SetXaxisRange(2., 4.);
  stack -> SetXaxisTitle("EE DPhi Ele Met");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 10., false);




  //==============
  //==== hMt =====
  //==============
  histoNames.at(0) = "hMt";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("transverse Mass (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  histoNames.at(0) = "hMt_EB";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("EB transverse Mass (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  histoNames.at(0) = "hMt_EE";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("EE transverse Mass (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, step, 4., true);

  histoNames.at(0) = "hMtNonIso";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("transverse Mass non iso (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 4., false);

  histoNames.at(0) = "hMtNonIso_EB";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("transverse Mass non iso EB (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 4., false);

  histoNames.at(0) = "hMtNonIso_EE";
  stack -> SetXaxisRange(0., 250.);
  stack -> SetXaxisTitle("transverse Mass non iso EE (GeV)");
  stack -> Draw(histoNames, "eventsScaled", lumi, 0, 4., false);

  
}
