{

  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  TFile *f[4];

  TFile *f[0] = TFile::Open("output/vtxCombined_Hgg-M120_Spring11-PU_S1_1sigmaz.root");
  TFile *f[1] = TFile::Open("output/vtxCombined_Hgg-M120_Spring11-PU_S1_2sigmaz.root");
  TFile *f[2] = TFile::Open("output/vtxCombined_Hgg-M120_Spring11-PU_S1_3sigmaz.root");
  TFile *f[3] = TFile::Open("output/vtxCombined_Hgg-M120_Spring11-PU_S1_5sigmaz.root");

  TH1F *hNvtxInDeltaZ[4];

  
  TGraphAsymmErrors *EffVsNvtx_cat3[4];
  TGraphAsymmErrors *EffVsPt_cat3[4];
  TGraphAsymmErrors *EffVsR9_cat3[4];

  for ( int i = 0 ; i < 4 ; i++){
    hNvtxInDeltaZ[i] = (TH1F*)f[i] ->Get("hNvtxInDeltaZ1");
    hNvtxInDeltaZ[i] ->SetLineColor(i+1);
    hNvtxInDeltaZ[i] ->SetLineWidth(2);
    hNvtxInDeltaZ[i] ->GetXaxis()->SetTitle("N_{vtx}");
    hNvtxInDeltaZ[i] ->GetYaxis()->SetTitle("a.u.");

    EffVsNvtx_cat3[i] = (TGraphAsymmErrors*)f[i] ->Get("EffVsNvtx_cat3");
    EffVsPt_cat3[i]   = (TGraphAsymmErrors*)f[i] ->Get("EffVsPt_cat3");
    EffVsR9_cat3[i]   = (TGraphAsymmErrors*)f[i] ->Get("EffVsR9_cat3");
        
    EffVsNvtx_cat3[i] ->SetLineColor(i+1);
    EffVsNvtx_cat3[i] ->SetMarkerColor(i+1);
    EffVsNvtx_cat3[i] ->GetXaxis()->SetTitle("number of reco PV");
    EffVsNvtx_cat3[i] ->GetYaxis()->SetTitle("efficiency");

    EffVsPt_cat3[i] ->SetLineColor(i+1);
    EffVsPt_cat3[i] ->SetMarkerColor(i+1);
    EffVsPt_cat3[i] ->GetXaxis()->SetTitle("boson p^{T} (GeV)");
    EffVsPt_cat3[i] ->GetYaxis()->SetTitle("efficiency");

    EffVsR9_cat3[i] ->SetLineColor(i+1);
    EffVsR9_cat3[i] ->SetMarkerColor(i+1);
    EffVsR9_cat3[i] ->GetXaxis()->SetTitle("max photon R9");
    EffVsR9_cat3[i] ->GetYaxis()->SetTitle("efficiency");


  }

  TLegend leg(0.68, 0.78, 0.99, 0.99);
  leg.SetFillColor(kWhite);
  leg.SetBorderSize(1);
  leg.AddEntry(hNvtxInDeltaZ[0],"1#sigma_{z}","L");
  leg.AddEntry(hNvtxInDeltaZ[1],"2#sigma_{z}","L");
  leg.AddEntry(hNvtxInDeltaZ[2],"3#sigma_{z}","L");
  leg.AddEntry(hNvtxInDeltaZ[3],"5#sigma_{z}","L");
 
  TCanvas *c1 = new TCanvas("c1","c1");
  hNvtxInDeltaZ[0]->GetXaxis()->SetRangeUser(0,10);
  hNvtxInDeltaZ[0]->DrawNormalized("");
  hNvtxInDeltaZ[1]->DrawNormalized("same");
  hNvtxInDeltaZ[2]->DrawNormalized("same");
  hNvtxInDeltaZ[3]->DrawNormalized("same");
  leg.Draw("same");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->SetGridx();
  c2->SetGridy();
  EffVsNvtx_cat3[0]->GetXaxis()->SetRangeUser(0,15);
  EffVsNvtx_cat3[0]->GetYaxis()->SetRangeUser(0.4,1.1);
  EffVsNvtx_cat3[0]->Draw("ap");
  EffVsNvtx_cat3[1]->Draw("psame");
  EffVsNvtx_cat3[2]->Draw("psame");
  EffVsNvtx_cat3[3]->Draw("psame");
  leg.Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetGridx();
  c3->SetGridy();
  EffVsPt_cat3[0]->GetXaxis()->SetRangeUser(0,200);
  EffVsPt_cat3[0]->GetYaxis()->SetRangeUser(0.4,1.1);
  EffVsPt_cat3[0]->Draw("ap");
  EffVsPt_cat3[1]->Draw("psame");
  EffVsPt_cat3[2]->Draw("psame");
  EffVsPt_cat3[3]->Draw("psame");
  leg.Draw("same");

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->SetGridx();
  c4->SetGridy();
  EffVsR9_cat3[0]->GetXaxis()->SetRangeUser(0,1);
  EffVsR9_cat3[0]->GetYaxis()->SetRangeUser(0.4,1.1);
  EffVsR9_cat3[0]->Draw("ap");
  EffVsR9_cat3[1]->Draw("psame");
  EffVsR9_cat3[2]->Draw("psame");
  EffVsR9_cat3[3]->Draw("psame");
  leg.Draw("same");  




}
