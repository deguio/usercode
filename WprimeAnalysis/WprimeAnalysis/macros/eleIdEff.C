
{
  TH1F* histoEt_HEEP = (TH1F*)_file0->Get("histoEt_HEEP");
  TH1F* histoEt_VBTF = (TH1F*)_file0->Get("histoEt_VBTF");


  histoEt_HEEP->Sumw2();
  histoEt_VBTF->Sumw2();

  histoEt_HEEP->Rebin(30);
  histoEt_VBTF->Rebin(30);

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors;
  ratio->BayesDivide(histoEt_VBTF,histoEt_HEEP);

  TCanvas* pp = new TCanvas("pp","pp");
  pp->SetGridx();
  pp->SetGridy();

  ratio->Draw("AP");


  histoEt_HEEP->SetMarkerColor(kRed);
  histoEt_HEEP->SetLineColor(kRed);
  histoEt_VBTF->SetMarkerColor(kBlack);
  histoEt_VBTF->SetLineColor(kBlack);

  TCanvas* ppp = new TCanvas("ppp","ppp");

  ppp->SetGridx();
  ppp->SetGridy();

  histoEt_HEEP->Draw();
  histoEt_VBTF->Draw("sames");
  
}
