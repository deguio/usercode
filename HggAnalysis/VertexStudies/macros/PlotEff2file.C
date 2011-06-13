{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  int saveScaleFactors = 1;

  TFile *f[2];
  f[0] = TFile::Open("outputEff/test_Zmumu_DATA2011_160404-165542.root"); 
  f[1] = TFile::Open("outputEff/test_ZmumuMC_Summer2011_PUweights_160404-165542.root");
 
//   string legtitle1 = "H#rightarrow#gamma#gamma <PU>=6 (sumpt2)";
//   string legtitle1bdt = "H#rightarrow#gamma#gamma <PU>=6 (Ranking)";

//   string legtitle2 = "H#rightarrow#gamma#gamma <PU>=10 (sumpt2)";
//   string legtitle2bdt = "H#rightarrow#gamma#gamma <PU>=10 (Ranking)";
  
  // legend
  string legtitle1 = "DATA Z#rightarrow#mu#mu(sumpt2)";
  string legtitle1bdt = "DATA Z#rightarrow#mu#mu(Ranking)";
  string legtitle2 = "MC Z#rightarrow#mu#mu(sumpt2)";
  string legtitle2bdt = "MC Z#rightarrow#mu#mu(Ranking)";

  // rebinning
  int nRePt = 2;
  // if variable size bins
  double xbins[13] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,50.,70.,110.,400.};

  // histograms
  TH1F* NvtxAll[2];
  TH1F* NvtxSumpt2[2];
  TH1F* NvtxBDT[2];

  TH1F* PtAll[2];
  TH1F* PtSumpt2[2];
  TH1F* PtBDT[2];

  TH1* PtAllRebinned[2];
  TH1* PtSumpt2Rebinned[2];
  TH1* PtBDTRebinned[2];

  TH1F* EtaAll[2];
  TH1F* EtaSumpt2[2];
  TH1F* EtaBDT[2];

  TGraphAsymmErrors *effVsNvtx_BDT[2];
  TGraphAsymmErrors *effVsNvtx_sumpt2[2];

  TGraphAsymmErrors *effVsPt_BDT[2];
  TGraphAsymmErrors *effVsPt_sumpt2[2];

  TGraphAsymmErrors *effVsEta_BDT[2];
  TGraphAsymmErrors *effVsEta_sumpt2[2];

  int mycolor, mystyle;
  char hname[100];

  for (int i = 0; i < 2; i++){
    
    NvtxAll[i] = (TH1F*) f[i]->Get("NvtAll");
    NvtxSumpt2[i] = (TH1F*) f[i]->Get("NvtGood");
    NvtxBDT[i] = (TH1F*) f[i]->Get("NvtGood_RANK");

    PtAll[i] = (TH1F*) f[i]->Get("PtAll");
    PtSumpt2[i] = (TH1F*) f[i]->Get("PtGood");
    PtBDT[i] = (TH1F*) f[i]->Get("PtGood_RANK");

    EtaAll[i] = (TH1F*) f[i]->Get("EtaAll");
    EtaSumpt2[i] = (TH1F*) f[i]->Get("EtaGood");
    EtaBDT[i] = (TH1F*) f[i]->Get("EtaGood_RANK");
  
    sprintf(hname,"PtAllRebinned_%d",i);
    PtAllRebinned[i] = PtAll[i]   ->Rebin(12,hname,xbins);
    sprintf(hname,"PtSumpt2Rebinned_%d",i);
    PtSumpt2Rebinned[i] = PtSumpt2[i]->Rebin(12,hname,xbins);
    sprintf(hname,"PtBDTRebinned_%d",i);
    PtBDTRebinned[i] = PtBDT[i]   ->Rebin(12,"PtBDTRebinned",xbins);

    effVsNvtx_BDT[i]     = new TGraphAsymmErrors();
    effVsNvtx_sumpt2[i]  = new TGraphAsymmErrors();
    effVsPt_BDT[i]     = new TGraphAsymmErrors();
    effVsPt_sumpt2[i]  = new TGraphAsymmErrors();
    effVsEta_BDT[i]    = new TGraphAsymmErrors();
    effVsEta_sumpt2[i] = new TGraphAsymmErrors();

    if (i==0) {
      mycolor = kGreen+1;
    }
    else {mycolor = kRed;
    }
    
      effVsNvtx_BDT[i]->SetLineColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerStyle(24);

    effVsPt_BDT[i]->SetLineColor(mycolor);
    effVsPt_BDT[i]->SetMarkerColor(mycolor);
    effVsPt_BDT[i]->SetMarkerStyle(24);

    effVsEta_BDT[i]->SetLineColor(mycolor);
    effVsEta_BDT[i]->SetMarkerColor(mycolor);
    effVsEta_BDT[i]->SetMarkerStyle(24);

    effVsNvtx_sumpt2[i]->SetLineColor(mycolor);
    effVsNvtx_sumpt2[i]->SetMarkerColor(mycolor);
    effVsNvtx_sumpt2[i]->SetMarkerStyle(20);

    effVsPt_sumpt2[i]->SetLineColor(mycolor);
    effVsPt_sumpt2[i]->SetMarkerColor(mycolor);
    effVsPt_sumpt2[i]->SetMarkerStyle(20);

    effVsEta_sumpt2[i]->SetLineColor(mycolor);
    effVsEta_sumpt2[i]->SetMarkerColor(mycolor);
    effVsEta_sumpt2[i]->SetMarkerStyle(20);

  }

  
  for (int i = 0; i<2; i++){

    effVsNvtx_sumpt2[i]->BayesDivide(NvtxSumpt2[i],NvtxAll[i], "cp"); 
    effVsNvtx_BDT[i]   ->BayesDivide(NvtxBDT[i],NvtxAll[i], "cp");
 
   //  effVsPt_sumpt2[i]->BayesDivide(PtSumpt2[i],PtAll[i], "cp"); 
//     effVsPt_BDT[i]   ->BayesDivide(PtBDT[i],PtAll[i], "cp"); 
    effVsPt_sumpt2[i]->BayesDivide(PtSumpt2Rebinned[i],PtAllRebinned[i], "cp"); 
    effVsPt_BDT[i]   ->BayesDivide(PtBDTRebinned[i],PtAllRebinned[i], "cp"); 

    effVsEta_sumpt2[i]->BayesDivide(EtaSumpt2[i],EtaAll[i], "cp"); 
    effVsEta_BDT[i]   ->BayesDivide(EtaBDT[i],EtaAll[i], "cp"); 
    
  }


  
  int nMBMax =16;
  
  TLegend legend1(0.68, 0.78, 0.99, 0.99);
  legend1.SetFillColor(kWhite);
  legend1.SetBorderSize(1);
  legend1.AddEntry(effVsNvtx_sumpt2[0],legtitle1.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_sumpt2[1],legtitle2.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[0],legtitle1bdt.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[1],legtitle2bdt.c_str(),"LP");



  //*** EFF vs NVTX
  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.,1.1);
  cc.SetStats(0); 
  // cc.GetXaxis()->SetTitle("number of reconstructed vertices"); 
  cc.GetXaxis()->SetTitle("number of simulated PU vertices"); 
  cc.GetYaxis()->SetTitle("efficiency");
  cc.Draw();
  for (int i = 0; i< 2; i++){
    effVsNvtx_sumpt2[i]->Draw("p,same");
    effVsNvtx_BDT[i]->Draw("p,same");
  }
  legend1.Draw("same");


  //*** EFF vs BOSON PT  
  TCanvas c2;
  c2.SetGridx();
  c2.SetGridy();
  TH2F dd("dd","",500,0,500,1000,0.0,1.1);
  dd.SetStats(0); 
  dd.GetXaxis()->SetTitle("boson p^{T} (GeV) "); 
  dd.GetYaxis()->SetTitle("efficiency");
  dd.Draw();
  for (int i = 0; i< 2; i++){
    effVsPt_sumpt2[i]->Draw("p,same");
    effVsPt_BDT[i]->Draw("p,same");
  }
  legend1.Draw("same");
  
  
  float ptlow = 0.;
  float pthigh = 400.;

  int bin1 = PtSumpt2[0]->FindBin(ptlow);
  int bin2 = PtSumpt2[0]->FindBin(pthigh);

   
  double eff1 = PtSumpt2[0]->Integral(bin1,bin2)/ PtAll[0]->Integral(bin1,bin2);
  double eff2 = PtSumpt2[1]->Integral(bin1,bin2)/ PtAll[1]->Integral(bin1,bin2);
  
  double eff1bdt = PtBDT[0]->Integral(bin1,bin2)/ PtAll[0]->Integral(bin1,bin2);
  double eff2bdt = PtBDT[1]->Integral(bin1,bin2)/ PtAll[1]->Integral(bin1,bin2);
  
  cout << "Efficiency integrated in boson pt : [" << ptlow << ","<< pthigh<<"] GeV"<< endl;
  cout << legtitle1.c_str() << " --> eff = " << eff1 << " +/- " 
       << sqrt(eff1*(1-eff1)/PtAll[0]->Integral(bin1,bin2))<< endl;  
  cout << legtitle2.c_str() << " --> eff = " << eff2 << " +/- " 
       << sqrt(eff2*(1-eff2)/PtAll[1]->Integral(bin1,bin2))<<endl;  
  cout << legtitle1bdt.c_str() << " --> eff = " << eff1bdt << " +/- " 
       << sqrt(eff1bdt*(1-eff1bdt)/PtAll[0]->Integral(bin1,bin2))<<endl;  
  cout << legtitle2bdt.c_str() << " --> eff = " << eff2bdt << " +/- " 
       << sqrt(eff2bdt*(1-eff2bdt)/PtAll[1]->Integral(bin1,bin2))<<endl;  
  
  //RATIOs

  //------ RATIO EFF vs PT -------------------------------------------
  TGraphErrors * ratioEffVsPt_sumpt2 = new TGraphErrors();
  ratioEffVsPt_sumpt2 ->SetMarkerStyle(20);
  
  int nbins = effVsPt_sumpt2[0]->GetHistogram()->GetNbinsX();
  
  double pt, err1, err2, ratio , err;
  
  for (int i = 0; i < nbins; i++){
    effVsPt_sumpt2[0]->GetPoint(i,pt,eff1);
    effVsPt_sumpt2[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsPt_sumpt2[0]->GetErrorY(i);
    err2 =  effVsPt_sumpt2[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsPt_sumpt2->SetPoint(i,pt,ratio);
    ratioEffVsPt_sumpt2->SetPointError(i, effVsPt_sumpt2[0]->GetErrorX(i) ,err);
  }
  
  TGraphErrors * ratioEffVsPt_BDT = new TGraphErrors();
  ratioEffVsPt_BDT ->SetMarkerStyle(24);

  for (int i = 0; i < nbins; i++){
    effVsPt_BDT[0]->GetPoint(i,pt,eff1);
    effVsPt_BDT[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsPt_BDT[0]->GetErrorY(i);
    err2 =  effVsPt_BDT[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsPt_BDT->SetPoint(i,pt,ratio);
    ratioEffVsPt_BDT->SetPointError(i,effVsPt_BDT[1]->GetErrorX(i),err);
  }
  
    
  
  TCanvas cRatioPt("cRatioPt","cRatioPt",600,300);
  cRatioPt.SetGridx();
  cRatioPt.SetGridy();

  ratioEffVsPt_sumpt2->GetHistogram()->GetXaxis()->SetTitle("boson p^{T} (GeV)");
  ratioEffVsPt_sumpt2->GetHistogram()->GetYaxis()->SetTitle("eff(data)/eff(MC)");
  ratioEffVsPt_sumpt2->GetHistogram()->GetXaxis()->SetRangeUser(0,200);
  ratioEffVsPt_sumpt2->GetHistogram()->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsPt_sumpt2->Draw("ap");
  ratioEffVsPt_BDT->Draw("psame");

  TLegend legend4(0.15, 0.2, 0.45, 0.4);
  legend4.SetFillColor(kWhite);
  legend4.SetBorderSize(1);
  legend4.AddEntry(ratioEffVsPt_sumpt2,"sumpt2","LP");
  legend4.AddEntry(ratioEffVsPt_BDT,"ranking","LP");
  legend4.Draw("same");



  //------ RATIO EFF vs NVTX -------------------------------------------
  TGraphErrors * ratioEffVsNvtx_sumpt2 = new TGraphErrors();
  ratioEffVsNvtx_sumpt2 ->SetMarkerStyle(20);
  
  int nbins = effVsNvtx_sumpt2[0]->GetHistogram()->GetNbinsX();
  
  double pt, err1, err2, ratio , err;
  
  for (int i = 0; i < nbins; i++){
    effVsNvtx_sumpt2[0]->GetPoint(i,pt,eff1);
    effVsNvtx_sumpt2[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsNvtx_sumpt2[0]->GetErrorY(i);
    err2 =  effVsNvtx_sumpt2[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsNvtx_sumpt2->SetPoint(i,pt,ratio);
    ratioEffVsNvtx_sumpt2->SetPointError(i,effVsNvtx_sumpt2[0]->GetErrorX(i),err);
  }
  
  TGraphErrors * ratioEffVsNvtx_BDT = new TGraphErrors();
  ratioEffVsNvtx_BDT ->SetMarkerStyle(24);

  for (int i = 0; i < nbins; i++){
    effVsNvtx_BDT[0]->GetPoint(i,pt,eff1);
    effVsNvtx_BDT[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsNvtx_BDT[0]->GetErrorY(i);
    err2 =  effVsNvtx_BDT[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsNvtx_BDT->SetPoint(i,pt,ratio);
    ratioEffVsNvtx_BDT->SetPointError(i,effVsNvtx_BDT[0]->GetErrorX(i),err);
  }

  TCanvas cRatioNvtx("cRatioNvtx","cRatioNvtx",600,300);
  cRatioNvtx.SetGridx();
  cRatioNvtx.SetGridy();

  ratioEffVsNvtx_sumpt2->GetHistogram()->GetXaxis()->SetTitle("number of reconstructed vertices");
  ratioEffVsNvtx_sumpt2->GetHistogram()->GetYaxis()->SetTitle("eff(data)/eff(MC)");
  ratioEffVsNvtx_sumpt2->GetHistogram()->GetXaxis()->SetRangeUser(0,16);
  ratioEffVsNvtx_sumpt2->GetHistogram()->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsNvtx_sumpt2->Draw("ap");
  ratioEffVsNvtx_BDT->Draw("psame");


  TLegend legend5(0.15, 0.2, 0.45, 0.4);
  legend5.SetFillColor(kWhite);
  legend5.SetBorderSize(1);
  legend5.AddEntry(ratioEffVsNvtx_sumpt2,"sumpt2","LP");
  legend5.AddEntry(ratioEffVsNvtx_BDT,"ranking","LP");

  legend5.Draw("same");

  if (saveScaleFactors){
    TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_Summer11_test.root","recreate");
    ratioEffVsPt_BDT->SetTitle("scaleFactor");
    ratioEffVsPt_BDT->Write("scaleFactor");
  }

  
}
