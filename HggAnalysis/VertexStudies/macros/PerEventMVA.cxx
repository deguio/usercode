
{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  TFile *f[2];
  f[0] = TFile::Open("outputTMVA/PerEventMVA_DYJetsToLL_S9_minBiasXsec68300_190456-194479.root");  
  f[1] = TFile::Open("outputTMVA/PerEventMVA_DoubleMu_190456-194479.root");
  
  
  // histograms
  TH1F* tmvaPerEventOutput_RightVtx[2];
  TH1F* tmvaPerEventOutput_WrongVtx[2];

  TH1F* tmvaPerEventOutput[2];

  int nre = 4;

  for (int i = 0; i < 2; i++){
    
    tmvaPerEventOutput_RightVtx[i] = (TH1F*)f[i]->Get("tmvaPerEventOutput_RightVtx");
    tmvaPerEventOutput_WrongVtx[i] = (TH1F*)f[i]->Get("tmvaPerEventOutput_WrongVtx");
    tmvaPerEventOutput_RightVtx[i] -> Sumw2();
    tmvaPerEventOutput_WrongVtx[i] -> Sumw2();
    tmvaPerEventOutput_RightVtx[i] -> Rebin(nre);
    tmvaPerEventOutput_WrongVtx[i] -> Rebin(nre);


    tmvaPerEventOutput[i] = (TH1F*) f[i]->Get("tmvaPerEventOutput");

    tmvaPerEventOutput[i] -> Sumw2();
    tmvaPerEventOutput[i] -> Rebin(nre);
      
    if (i==0){
      tmvaPerEventOutput[i] -> SetFillColor(kAzure+1);
      tmvaPerEventOutput_RightVtx[i] -> SetFillColor(kGreen+1);
      tmvaPerEventOutput_WrongVtx[i] -> SetFillColor(kRed+1);
      tmvaPerEventOutput_RightVtx[i] -> SetFillStyle(3002);
      tmvaPerEventOutput_WrongVtx[i] -> SetFillStyle(3005);
    }
      
    if (i==1){
      tmvaPerEventOutput[i] -> SetMarkerStyle(20);
      tmvaPerEventOutput[i] -> SetMarkerSize(0.8);
      
      tmvaPerEventOutput_RightVtx[i] -> SetMarkerStyle(20);
      tmvaPerEventOutput_WrongVtx[i] -> SetMarkerStyle(20);
      tmvaPerEventOutput_RightVtx[i] -> SetMarkerSize(0.8);
      tmvaPerEventOutput_WrongVtx[i] -> SetMarkerSize(0.8);
      tmvaPerEventOutput_RightVtx[i] -> SetMarkerColor(kGreen+2);
      tmvaPerEventOutput_WrongVtx[i] -> SetMarkerColor(kRed+2);
    }
    
  }


  TLegend leg (0.6, 0.6,0.89,0.75);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(tmvaPerEventOutput[0],"MC Z#rightarrow#mu#mu","F");
  leg->AddEntry(tmvaPerEventOutput[1],"Data Z#rightarrow#mu#mu","LP");

  TLatex *latex = new TLatex(0.55,0.85,"#splitline{          CMS preliminary}{#sqrt{s} = 8 TeV L = 1.92 fb^{-1}}");
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);


  TCanvas *c0 = new TCanvas("c0","c0",500,500);
  tmvaPerEventOutput[0] -> GetXaxis() -> SetTitle("MVA_{event}"); 
  tmvaPerEventOutput[0] -> GetXaxis() -> SetRangeUser(-1,1); 
  tmvaPerEventOutput[0] -> DrawNormalized("histo");
  tmvaPerEventOutput[1] ->DrawNormalized("esame");
  leg->Draw("same");
  latex->Draw("same");
  
  TLegend leg2 (0.6, 0.6,0.89,0.75);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(tmvaPerEventOutput_RightVtx[0],"right vertex MC","F");
  leg2->AddEntry(tmvaPerEventOutput_WrongVtx[0],"wrong vertex MC","F");
  leg2->AddEntry(tmvaPerEventOutput_RightVtx[1],"right vertex DATA","LP");
  leg2->AddEntry(tmvaPerEventOutput_WrongVtx[1],"wrong vertex DATA","LP");

  TCanvas *cright = new TCanvas("cright","cright",500,500);
  tmvaPerEventOutput_RightVtx[0] -> GetXaxis() -> SetTitle("MVA_{event}");
  tmvaPerEventOutput_RightVtx[0] -> DrawNormalized("histo");
  tmvaPerEventOutput_WrongVtx[0] -> DrawNormalized("histo same");
  tmvaPerEventOutput_RightVtx[1] -> DrawNormalized("esame");
  tmvaPerEventOutput_WrongVtx[1] -> DrawNormalized("esame");
  leg2->Draw("same");
  latex->Draw("same");

 

}
