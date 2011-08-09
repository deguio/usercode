// compilo con 
// c++ -o QCD_data_driven_roofit_def `root-config --cflags --ldflags --glibs` -lRooFit -lRooFitCore -lMinuit QCD_data_driven_roofit_def.C setTDRStyle.cc

#include "setTDRStyle.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"

#include <cstdlib>
#include <string>
#include <map>
#include <utility>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/stat.h> 

using namespace RooFit;

//main
int main(int argc, char** argv)
{

  setTDRStyle();
  TApplication *theApp = new TApplication( "app", &argc, argv );
  
  //parameters
  std::ifstream listFile( "/data2/deguio/Wprime/WprimeAnalysis/WprimeAnalysis/config/crossSections_wPrime_Spring11.txt" );
  float lumi = 43.40;

  //int step = 12; // FIXME
  int rebin = 4; // FIXME

  //-----------------------
  //---- opening files ----
  //-----------------------

  //data
  TFile* fin_DATA       = new TFile("output/DATA_isoHLT/histo_WprimeTreeAnalysis.root","READ");
  TH1F* etOverMet_DATA  = (TH1F*)fin_DATA->Get("hEtOverMet/h_8_hEtOverMet"); //FIXME step -1 et_over_met rimane l'ultimo?
  etOverMet_DATA -> Rebin(rebin);
  etOverMet_DATA -> Sumw2();
  TH1F* mt_DATA         = (TH1F*)fin_DATA->Get("hMt/h_9_hMt"); //FIXME
  mt_DATA -> SetMarkerStyle(20);
  mt_DATA -> Rebin(rebin);
  //mt_DATA -> Sumw2();
  TH1F* hMt_cumulative_DATA      = (TH1F*)fin_DATA->Get("hMt_cumulative/h_0_hMt_cumulative");
  hMt_cumulative_DATA -> SetMarkerStyle(20);
  hMt_cumulative_DATA -> Rebin(rebin);
  //hMt_cumulative_DATA -> Sumw2();

  //file qcd template from data non iso
  TFile* fin_DATA_NONISO       = new TFile("output/DATA_nonIso/histo_WprimeTreeAnalysis.root","READ");
  TH1F* etOverMet_QCD   = (TH1F*)fin_DATA_NONISO->Get("hEtOverMet/h_8_hEtOverMet");
  etOverMet_QCD -> Rebin(rebin);
  etOverMet_QCD -> Sumw2();
  TH1F* mt_QCD          = (TH1F*)fin_DATA_NONISO->Get("hMt/h_9_hMt");  
  mt_QCD -> SetLineColor(kGray+1);
  mt_QCD -> SetFillColor(kGray+1);
  mt_QCD -> SetFillStyle(3003);
  mt_QCD -> SetLineWidth(2);
  mt_QCD -> Rebin(rebin);
  //mt_QCD -> Sumw2();
  TH1F* hMt_cumulative_QCD      = (TH1F*)fin_DATA_NONISO->Get("hMt_cumulative/h_0_hMt_cumulative");
  //hMt_cumulative_QCD -> SetLineColor(kGray+1);
  //hMt_cumulative_QCD -> SetFillColor(kGray+1);
  //hMt_cumulative_QCD -> SetFillStyle(3003);
  //hMt_cumulative_QCD -> SetLineWidth(2);
  hMt_cumulative_QCD -> Rebin(rebin);
  //hMt_cumulative_QCD -> Sumw2();

  //file W template from MC
  TFile* fin_WenuMC = new TFile("output/MC/WToENu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_AODSIM/histo_WprimeTreeAnalysis.root","READ");
  TH1F* etOverMet_Wenu = (TH1F*)fin_WenuMC->Get("hEtOverMet/h_8_hEtOverMet");
  etOverMet_Wenu -> Rebin(rebin);
  etOverMet_Wenu -> Sumw2();
  TH1F* mt_Wenu = (TH1F*)fin_WenuMC->Get("hMt/h_9_hMt");
  mt_Wenu -> SetLineColor(kAzure+1);
  mt_Wenu -> SetFillColor(kAzure+1);
  mt_Wenu -> SetFillStyle(3003);
  mt_Wenu -> SetLineWidth(2);
  mt_Wenu -> Rebin(rebin);
  //mt_Wenu -> Sumw2();
  TH1F* hMt_cumulative_Wenu      = (TH1F*)fin_WenuMC->Get("hMt_cumulative/h_0_hMt_cumulative");
  //hMt_cumulative_Wenu -> SetLineColor(kAzure+1);
  //hMt_cumulative_Wenu -> SetFillColor(kAzure+1);
  //hMt_cumulative_Wenu -> SetFillStyle(3003);
  //hMt_cumulative_Wenu -> SetLineWidth(2);
  hMt_cumulative_Wenu -> Rebin(rebin);
  //hMt_cumulative_Wenu -> Sumw2();
  
  //file other MC with fixed normalization
  TH1F* etOverMet_OTHER = NULL;
  TH1F* mt_OTHER = NULL;
  TH1F* mt_cumulative_OTHER = NULL;
  bool isFirst = true;
  std::map<std::string, TH1F*> mt_other_map;
  std::map<std::string, TH1F*> mt_cumulative_other_map;

  while(!listFile.eof())
    {
      std::string sample;
      std::string sumName;
      int dataFlag;
      int isDD;
      double mH;
      double crossSection;
      int color;
      
      listFile >> sample >> sumName >> dataFlag >> isDD >> mH >> crossSection >> color;
      
      if(sample.size() == 0)
	continue;
      if(sample.at(0) == '#')
	continue;
      if(isDD != 0)
	continue;
      if(dataFlag != 0)
	continue;
      
      std::cout << sample << "\t"
		<< sumName << "\t"
		<< dataFlag << "\t"
		<< isDD << "\t"
		<< mH << "\t"
		<< crossSection << "\t" 
		<< color << "\t" 
		<< std::endl;
      


      //apro il file per ogni fondo
      TFile* fin_dummy      = new TFile(("output/MC/"+sample+"/histo_WprimeTreeAnalysis.root").c_str(),"READ");

      TH1F* events          = (TH1F*)fin_dummy->Get("events");

      TH1F* dummy_etOverMET = (TH1F*)fin_dummy->Get("hEtOverMet/h_8_hEtOverMet");
      dummy_etOverMET -> Rebin(rebin);
      dummy_etOverMET -> Sumw2();
      dummy_etOverMET -> Scale(lumi*crossSection/events -> GetBinContent(1));

      TH1F* dummy_mt        = (TH1F*)fin_dummy->Get("hMt/h_9_hMt");
      dummy_mt        -> Rebin(rebin);
      //dummy_mt        -> Sumw2();
      dummy_mt        -> Scale(lumi*crossSection/events -> GetBinContent(1));
      dummy_mt    -> SetLineColor(color);
      dummy_mt    -> SetFillColor(color);
      dummy_mt    -> SetFillStyle(3003);
      dummy_mt    -> SetLineWidth(2);
      
      TH1F* dummy_mt_cumulative      = (TH1F*)fin_dummy->Get("hMt_cumulative/h_0_hMt_cumulative");
      dummy_mt_cumulative -> Rebin(rebin);
      //dummy_mt_cumulative        -> Sumw2();
      dummy_mt_cumulative -> Scale(lumi*crossSection/events -> GetBinContent(1));
      //dummy_mt_cumulative -> SetLineColor(color);
      //dummy_mt_cumulative -> SetFillColor(color);
      //dummy_mt_cumulative -> SetFillStyle(3003);
      //dummy_mt_cumulative -> SetLineWidth(2);
      
      std::map<std::string, TH1F*>::const_iterator mapItr = mt_other_map.find(sumName);
      if (mapItr == mt_other_map.end())
	{
	  mt_other_map[sumName] = (TH1F*)dummy_mt -> Clone();
	  mt_cumulative_other_map[sumName] = (TH1F*)dummy_mt_cumulative -> Clone();
	}
      else
	{
	  mt_other_map[sumName] -> Add(dummy_mt);
	  mt_cumulative_other_map[sumName] -> Add(dummy_mt_cumulative);	  
	}
      
      if (isFirst == true)
	{
	  etOverMet_OTHER = (TH1F*)dummy_etOverMET -> Clone();
	  mt_OTHER        = (TH1F*)dummy_mt -> Clone();
	  mt_cumulative_OTHER = (TH1F*)dummy_mt_cumulative -> Clone();
	  isFirst         = false;
	}
      else
	{
	  etOverMet_OTHER -> Add(dummy_etOverMET);
	  mt_OTHER        -> Add(dummy_mt);
	  mt_cumulative_OTHER -> Add(dummy_mt_cumulative);
	}
    }
  
  //----------------------------
  //---- creating templates ----
  //----------------------------

  // define the RooRealVars
  RooRealVar x("x", "Et over Met", 0., 500.);
  
  
  // define the RooDataHist
  RooDataHist rooDataHisto("rooDataHisto",   "data",        RooArgList(x), etOverMet_DATA);
  RooDataHist rooQCDHisto("rooQCDHisto",     "QCD non iso", RooArgList(x), etOverMet_QCD);
  RooDataHist rooWenuHisto("rooWenuHisto",   "Wenu",        RooArgList(x), etOverMet_Wenu);
  RooDataHist rooOtherHisto("rooOtherHisto", "other bkg",   RooArgList(x), etOverMet_OTHER);
  
  
  // define the signal/bkg shapes from histograms
  RooHistPdf QCDPdf("QCDPdf",     "QCDPdf",   x, rooQCDHisto,   0) ;
  RooHistPdf WenuPdf("WenuPdf",   "WenuPdf",  x, rooWenuHisto,  0) ;
  RooHistPdf OtherPdf("OtherPdf", "OtherPdf", x, rooOtherHisto, 0) ;
  
  
  // define the total shape
  RooRealVar NQCD(  "NQCD",   "NQCD",   150., 0., 10000000.);
  RooRealVar NWenu( "NWenu",  "NWenu",  40.,  0., 10000000.);
  //RooRealVar NOther( "NOther",  "NOther",  10., 0., 10000000.);
  RooRealVar NOther("NOther", "NOther", etOverMet_OTHER -> Integral(), etOverMet_OTHER -> Integral()-1, etOverMet_OTHER -> Integral()+1);  //FIXME -> come fisso il parametro?
  
  RooAddPdf totPdf("totPdf", "tot", RooArgList(QCDPdf,WenuPdf,OtherPdf), RooArgList(NQCD,NWenu,NOther));
  
  // fit the histo
  totPdf.fitTo(rooDataHisto, RooFit::Extended(1), RooFit::Minos(1), RooFit::SumW2Error(kTRUE));


  // plot  Et / MET
  RooPlot* rooPlot = x.frame();
  rooDataHisto.plotOn(rooPlot, DataError(RooAbsData::SumW2));
  totPdf.plotOn(rooPlot, LineColor(kBlue));
  totPdf.plotOn(rooPlot, Components(QCDPdf), LineStyle(kDashed), LineColor(kGray+2));
  totPdf.plotOn(rooPlot, Components(WenuPdf), LineStyle(kDashed), LineColor(kAzure+1));
  totPdf.plotOn(rooPlot, Components(OtherPdf), LineStyle(kDashed), LineColor(kGreen));
  TCanvas c1;
  rooPlot->Draw();
  
  char QCDString[50];
  sprintf(QCDString, "NQCD = %.3f^{+%.1f}_{-%.1f}", NQCD.getVal(), NQCD.getAsymErrorHi(), NQCD.getAsymErrorLo());
  
  char WenuString[50];
  sprintf(WenuString, "NWenu = %.3f^{+%.1f}_{-%.1f}", NWenu.getVal(), NWenu.getAsymErrorHi(), NWenu.getAsymErrorLo());
  
  char OtherString[50];
  sprintf(OtherString, "NOther = %.3f^{+%.1f}_{-%.1f}", NOther.getVal(), NOther.getAsymErrorHi(), NOther.getAsymErrorLo());

  TLatex QCDLatex(0.65, 0.85, QCDString);
  QCDLatex.SetNDC();
  QCDLatex.SetTextFont(42);
  QCDLatex.SetTextSize(0.03);
  QCDLatex.Draw("same");
  
  TLatex WenuLatex(0.65, 0.75, WenuString);
  WenuLatex.SetNDC();
  WenuLatex.SetTextFont(42);
  WenuLatex.SetTextSize(0.03);
  WenuLatex.Draw("same");
  
  TLatex OtherLatex(0.65, 0.65, OtherString);
  OtherLatex.SetNDC();
  OtherLatex.SetTextFont(42);
  OtherLatex.SetTextSize(0.03);
  OtherLatex.Draw("same");


  //------------Integral in the range -> scale to mT-----------------
  //                ---non lo so fare con rooFit---
  etOverMet_QCD -> Scale(NQCD.getVal() / etOverMet_QCD->Integral());  //assegno l'area ottenuta col fit
  double mT_QCD_scaleFactor = etOverMet_QCD -> Integral(etOverMet_QCD -> FindBin(0.4), etOverMet_QCD -> FindBin(1.5)); //integro nel range di interesse
  
  etOverMet_Wenu -> Scale(NWenu.getVal() / etOverMet_Wenu->Integral());
  double mT_Wenu_scaleFactor = etOverMet_Wenu -> Integral(etOverMet_Wenu -> FindBin(0.4), etOverMet_Wenu -> FindBin(1.5));
  
  std::cout << std::endl;
  std::cout << "0.4 < NWenu < 1.5 : " << mT_Wenu_scaleFactor << std::endl;
  std::cout << "Wenu scale factor = " << mT_Wenu_scaleFactor / mt_Wenu -> Integral() << std::endl;
  std::cout << "0.4 < NQCD < 1.5 : " << mT_QCD_scaleFactor << std::endl;
  std::cout << "QCD scale factor = " << mT_QCD_scaleFactor / mt_QCD -> Integral() << std::endl;
  std::cout << std::endl;

  //--- scaling ---
  mt_QCD -> Scale(mT_QCD_scaleFactor / mt_QCD -> Integral()); //scalo mt rispettivamente
  mt_Wenu -> Scale(mT_Wenu_scaleFactor / mt_Wenu -> Integral());

  //-------- plots histos --------
  //--- common latex ---

  std::string cmsinfo = "CMS Preliminary 2011";
 TLatex CMSInfoLatex(0.38, 0.88, cmsinfo.c_str());
 CMSInfoLatex.SetNDC();
 CMSInfoLatex.SetTextFont(42);
 CMSInfoLatex.SetTextSize(0.04);

 std::string lumiinfo = "#int L dt = 43.40 pb^{-1}";
 TLatex LUMIInfoLatex(0.45, 0.82, lumiinfo.c_str());
 LUMIInfoLatex.SetNDC();
 LUMIInfoLatex.SetTextFont(42);
 LUMIInfoLatex.SetTextSize(0.03);

 std::string comeinfo = "#sqrt{s} = 7 TeV";
 TLatex COMEInfoLatex(0.48, 0.75, comeinfo.c_str());
 COMEInfoLatex.SetNDC();
 COMEInfoLatex.SetTextFont(42);
 COMEInfoLatex.SetTextSize(0.03);

  //--- mt ---
  TLegend legend(0.68, 0.78, 0.99, 0.99);
  legend.SetFillColor(kWhite);

  THStack* stack_mt = new THStack("stack_mt","stack_mt");
  for (std::map<std::string, TH1F*>::const_iterator itr = mt_other_map.begin();
       itr != mt_other_map.end();
       ++itr)
    {
      stack_mt -> Add(itr -> second);
      legend.AddEntry(itr -> second, (itr -> first).c_str(), "F");
    }

  TCanvas c2;
  c2.SetLogy();
  c2.cd();

  stack_mt -> Add(mt_QCD);
  legend.AddEntry(mt_QCD, "08_QCD", "F");
  stack_mt -> Add(mt_Wenu);
  legend.AddEntry(mt_Wenu, "10_Wenu", "F");
  legend.AddEntry(mt_DATA, "DATA", "P");
      
  stack_mt -> Draw();
  mt_DATA->Draw("E same");
  legend.Draw("same");
  CMSInfoLatex.Draw("same");
  LUMIInfoLatex.Draw("same");
  COMEInfoLatex.Draw("same");

  stack_mt -> GetXaxis() -> SetTitle("transverse mass (GeV)");
  stack_mt -> GetXaxis() -> SetRangeUser(0.,500.);
  stack_mt -> SetMaximum(10000.);
  stack_mt -> SetMinimum(0.1);

  //--- cumulative ---  
  TLegend legend2(0.74,0.78,0.92,0.92);
 legend2.SetFillColor(kWhite);
 legend2.SetTextFont(42);
 legend2.SetTextSize(0.03);

 THStack* stack_cumul = new THStack("stack_cumul","stack_cumul");
 for (std::map<std::string, TH1F*>::const_iterator itr = mt_cumulative_other_map.begin();
      itr != mt_cumulative_other_map.end();
      ++itr)
   {
     stack_cumul -> Add(itr -> second);
   }
 
 TCanvas c3;
 c3.SetLogy();
 c3.cd();
 c3.SetGridx();
 c3.SetGridy();

 mt_Wenu->Rebin(rebin);
 hMt_cumulative_Wenu->Rebin(rebin);
 mt_QCD->Rebin(rebin);
 hMt_cumulative_QCD->Rebin(rebin);
 mt_OTHER->Rebin(rebin);
 mt_cumulative_OTHER->Rebin(rebin);
 mt_DATA->Rebin(rebin);
 hMt_cumulative_DATA->Rebin(rebin);


 for (int kk = 1; kk < mt_Wenu->GetNbinsX()+1; ++kk)
   hMt_cumulative_Wenu -> SetBinContent(kk, mt_Wenu->Integral(kk,mt_Wenu->GetNbinsX()));
 for (int kk = 1; kk < mt_QCD->GetNbinsX()+1; ++kk)
   hMt_cumulative_QCD -> SetBinContent(kk, mt_QCD->Integral(kk,mt_QCD->GetNbinsX()));
 for (int kk = 1; kk < mt_OTHER->GetNbinsX()+1; ++kk)
   mt_cumulative_OTHER -> SetBinContent(kk, mt_OTHER->Integral(kk,mt_OTHER->GetNbinsX()));
 
  mt_cumulative_OTHER -> Add(hMt_cumulative_QCD);
  mt_cumulative_OTHER -> Add(hMt_cumulative_Wenu);
  mt_cumulative_OTHER -> SetFillColor(kAzure+1);
  //mt_cumulative_OTHER -> SetLineColor(kAzure+1);
  //mt_cumulative_OTHER -> SetFillStyle(0);
  mt_cumulative_OTHER -> GetXaxis() -> SetTitle("M_{T}");
  mt_cumulative_OTHER -> GetYaxis() -> SetTitle("Events > x");
  mt_cumulative_OTHER -> GetXaxis() -> SetRangeUser(0.,500.);
  mt_cumulative_OTHER -> SetMaximum(100000.);
  mt_cumulative_OTHER -> SetMinimum(0.1);
  
  mt_cumulative_OTHER -> Draw();
  legend2.AddEntry(hMt_cumulative_DATA, "Data");
  legend2.AddEntry(mt_cumulative_OTHER, "Total Bkg");

  for (int kk = 1; kk < mt_QCD->GetNbinsX()+1; ++kk)
    hMt_cumulative_DATA -> SetBinContent(kk, mt_DATA->Integral(kk,mt_DATA->GetNbinsX()));

  hMt_cumulative_DATA->Draw("E same");

  legend2.Draw("same");
  CMSInfoLatex.Draw("same");
  LUMIInfoLatex.Draw("same");
  COMEInfoLatex.Draw("same");

 
  TFile* out = new TFile("mt_QCD_DD.root","RECREATE");
  mt_QCD -> Write();
  out -> Close();

  
  theApp -> Run();
  return 0;
}
