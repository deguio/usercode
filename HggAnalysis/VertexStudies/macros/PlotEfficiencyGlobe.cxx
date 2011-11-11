{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  bool writeEffFile = false;
  bool useVariableBinning = false;

  TFile *f = TFile::Open("outputEff/test_globe_Glu121_S6_PUweights_2011A_0100_73500.root");

  int nRePt = 1;

  TH1F *NvtAll;
  TH1F *NvtGood;
  
  TH1F *PtAll;
  TH1F *PtGood;

  //double xbins[13] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,50.,70.,110.,250.};
  double xbins[9] = {0.,10.,20.,30.,40.,50.,70.,110.,250.};

  TGraphAsymmErrors *EffVsNvtx = new TGraphAsymmErrors();
  TGraphAsymmErrors *EffVsPt   = new TGraphAsymmErrors();

  NvtAll  = (TH1F*)f->Get("NpuAll_BDT");
  NvtGood = (TH1F*)f->Get("NpuGood_BDT");
  
  PtAll  = (TH1F*)f->Get("PtAll_BDT");
  PtGood = (TH1F*)f->Get("PtGood_BDT");

  if (useVariableBinning){
    PtAll->Rebin(8,"PtAllRebinned",xbins);
    PtGood->Rebin(8,"PtGoodRebinned",xbins);
  }

  float ptlow = 0.;
  float pthigh = 400.;
  int bin1 = PtAll->FindBin(ptlow);
  int bin2 = PtAll->FindBin(pthigh);

  double eff;
  double err;

  double eff = PtGood->Integral(bin1,bin2)/ PtAll->Integral(bin1,bin2);
  double err = sqrt(eff*(1-eff)/PtAll->Integral(bin1,bin2));
  
  cout << "Efficiency integrated in boson pt : [" << ptlow << ","<< pthigh<<"] GeV"<< endl;
  cout << "Total (no categorories): eff = " << eff << " +/- " << err<< endl;  
 
  
  EffVsNvtx -> Divide(NvtGood, NvtAll, "cp");
  if (useVariableBinning)
    EffVsPt   -> Divide(PtGoodRebinned, PtAllRebinned, "cp");
  else
    EffVsPt   -> Divide(PtGood, PtAll, "cp");


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGridx();
  c1->SetGridy();
  EffVsNvtx->GetXaxis()->SetRangeUser(1,15);
  EffVsNvtx->GetYaxis()->SetRangeUser(0.0,1.1);
  EffVsNvtx->GetXaxis()->SetTitle("number of PV"); 
  EffVsNvtx->GetYaxis()->SetTitle("efficiency"); 
  EffVsNvtx->Draw("ap");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->SetGridx();
  c2->SetGridy();
  EffVsPt->GetXaxis()->SetRangeUser(1,250);
  EffVsPt->GetYaxis()->SetRangeUser(0.0,1.1);
  EffVsPt->GetXaxis()->SetTitle("p_{T}(H) (GeV/c)"); 
  EffVsPt->GetYaxis()->SetTitle("Fraction of events"); 
  EffVsPt->Draw("ap");

  TPaveText *pt = new TPaveText(15,0.05,90,0.3);
  pt->SetFillStyle(1001);
  pt->SetFillColor(0);
  //pt->SetBorderSize(1);
  pt->AddText("#splitline{#splitline{CMS preliminary}{Simulation}}{#LTn_{PU}#GT=6.6}");
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->Draw("same");

  TLegend *leg = new TLegend(0.6,0.65,0.89,0.7);
  //leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg ->AddEntry(EffVsPt,"z_{reco}-z_{true} < 10 mm","LP");
  leg->Draw("same");

  if (writeEffFile){
    TFile *fileout = new TFile("vtxIdEff_vs_bosonPt_globe_Glu121_S6_PUweights_2011B_0100_73500.root","recreate");
    EffVsPt->SetTitle("efficiency");
    EffVsPt->Write("efficiency");
  }
}


