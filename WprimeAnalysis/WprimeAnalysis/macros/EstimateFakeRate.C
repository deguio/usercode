#include <iostream>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TArrow.h"
#include "TCanvas.h"
#include <stdio.h>
#include <TLatex.h>


void EstimateFakeRate()
{

  // int numeratorStep = 10;

  // char treeName[50];
  // sprintf(treeName, "ntu_%d", numeratorStep);

  // TChain* dataChain = new TChain(treeName);
  // dataChain -> Add("/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/output/FAKERATE/estimateFakeRate/Photon_Run2011A-May10ReReco-v1_AOD/ntu_WprimeTreeAnalysis.root");

  TH1F* numerator_EB   = (TH1F*)_file0->Get("numerator_EB");
  TH1F* numerator_EE   = (TH1F*)_file0->Get("numerator_EE");
  TH1F* denominator_EB = (TH1F*)_file0->Get("denominator_EB");
  TH1F* denominator_EE = (TH1F*)_file0->Get("denominator_EE");

  TGraphAsymmErrors* fakeRate_EB = new TGraphAsymmErrors;
  TGraphAsymmErrors* fakeRate_EE = new TGraphAsymmErrors;

  fakeRate_EB->BayesDivide(numerator_EB, denominator_EB);
  fakeRate_EE->BayesDivide(numerator_EE, denominator_EE);

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetGridx();
  c1->SetGridy();
  fakeRate_EB->GetXaxis()->SetTitle("photon Et [GeV] EB");
  fakeRate_EB->GetYaxis()->SetTitle("fake rate EB");
  fakeRate_EB->Draw("AP");

  TCanvas* c2 = new TCanvas("c2","c2");
  c2->SetGridx();
  c2->SetGridy();
  fakeRate_EE->GetXaxis()->SetTitle("photon Et [GeV] EE");
  fakeRate_EE->GetYaxis()->SetTitle("fake rate EE");
  fakeRate_EE->Draw("AP");
    
}
