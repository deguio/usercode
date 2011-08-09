
{
  TH1F* histoEt_VBTF_All = (TH1F*)_file0->Get("histoEt_VBTF_All");
  TH1F* histoEt_VBTF_hlt = (TH1F*)_file0->Get("histoEt_VBTF_hlt");
  TH1F* histoEt_VBTF_wp80 = (TH1F*)_file0->Get("histoEt_VBTF_wp80");

  histoEt_VBTF_All->Sumw2();
  histoEt_VBTF_hlt->Sumw2();
  histoEt_VBTF_wp80->Sumw2();

  histoEt_VBTF_All->Rebin(10);
  histoEt_VBTF_hlt->Rebin(10);
  histoEt_VBTF_wp80->Rebin(10);

  TGraphAsymmErrors* ratio_hlt = new TGraphAsymmErrors;
  ratio_hlt->BayesDivide(histoEt_VBTF_hlt,histoEt_VBTF_All);

  TGraphAsymmErrors* ratio_wp80 = new TGraphAsymmErrors;
  ratio_wp80->BayesDivide(histoEt_VBTF_wp80,histoEt_VBTF_All);

  TCanvas* pp = new TCanvas("pp","pp");
  pp->SetGridx();
  pp->SetGridy();

  ratio_wp80->Draw("AP");

  ratio_hlt->SetMarkerColor(kRed);
  ratio_hlt->SetLineColor(kRed);
  //ratio_hlt->Draw("P,same");

  
}
