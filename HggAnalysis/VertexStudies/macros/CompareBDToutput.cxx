{

  
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  int saveScaleFactors = 0;

  TFile *f[2];
  
  //RUN2011A
  f[0] = TFile::Open("outputEff/test_Zmumu_Run2011A.root");  
  f[1] = TFile::Open("outputEff/test_DYJetsToLL_Fall11_S6_PUweights_2011A_0100_73500.root");

  //RUN2011B
  //f[0] = TFile::Open("outputEff/test_Zmumu_Run2011B-PromptReco-v1.root");  
  //f[1] = TFile::Open("outputEff/test_DYJetsToLL_Fall11_S6_PUweights_2011B_0100_73500.root");

  // FULL 2011
  // f[0] = TFile::Open("outputEff/test_Zmumu_2011.root");  
  //   f[1] = TFile::Open("outputEff/test_DYJetsToLL_Fall11_S6_PUweights_2011_0100_73500.root");

  TLatex *latex = new TLatex(0.55,0.85,"#splitline{          CMS preliminary}{#sqrt{s} = 7 TeV  Run2011A}");
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);

  
  TH1F* BDToutput_sig[2];
  TH1F* BDToutput_bkg[2];
  int nRe = 10;

  for (int i = 0; i < 2 ; i++){
    BDToutput_sig[i]= (TH1F*)f[i]->Get("BDToutput_sig");
    BDToutput_bkg[i]= (TH1F*)f[i]->Get("BDToutput_bkg");
    BDToutput_sig[i]->Sumw2();
    BDToutput_bkg[i]->Sumw2();
    
    BDToutput_sig[i]->Rebin(nRe);
    BDToutput_bkg[i]->Rebin(nRe);

    if (i==0){
      BDToutput_bkg[i]-> SetMarkerColor(kRed+2);
      BDToutput_bkg[i]-> SetMarkerStyle(20);
      BDToutput_sig[i]-> SetMarkerColor(kGreen+2);
      BDToutput_sig[i]-> SetMarkerStyle(20);
    }
    else{

      BDToutput_bkg[1]-> SetFillColor(kRed+2);
      BDToutput_bkg[1]-> SetFillStyle(3005);  

      BDToutput_sig[1]-> SetFillColor(kGreen+2);
      BDToutput_sig[1]-> SetFillStyle(3004);  
    
    }

  }

  
  TCanvas *c = new TCanvas("c","c",500,500);
  BDToutput_bkg[1]->GetXaxis()->SetTitle("BDT output");
  BDToutput_bkg[1]->DrawNormalized("histo");
  BDToutput_bkg[0]->DrawNormalized("esame"); 
  BDToutput_sig[1]->DrawNormalized("histo same");
  BDToutput_sig[0]->DrawNormalized("esame");
  latex->Draw("same");

  TLegend leg2 (0.6, 0.6,0.89,0.75);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(BDToutput_sig[1],"right vertex MC","F");
  leg2->AddEntry(BDToutput_bkg[1],"wrong vertex MC","F");
  leg2->AddEntry(BDToutput_sig[0],"right vertex DATA","LP");
  leg2->AddEntry(BDToutput_bkg[0],"wrong vertex DATA","LP");
  leg2->Draw("same");

  
}
