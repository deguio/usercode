{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  int saveScaleFactors = 0;

  TFile *f[2];
  f[0] = TFile::Open("outputEff/test_Zmumu_Run2011B-PromptReco-v1_AOD_upto178677.root");  
  f[1] = TFile::Open("outputEff/test_Zmumu_Run2011B-PromptReco-v1_AOD_upto178677.root");
  
  // legend
  string legtitle1    = "DATA Z#rightarrow#mu#mu (BDT)";
  string legtitle1bdt = "DATA Z#rightarrow#mu#mu (RANK)";
  string legtitle2    = "MC Z#rightarrow#mu#mu (BDT)";
  string legtitle2bdt = "MC Z#rightarrow#mu#mu (RANK)";

  // rebinning
  int nRePt = 1;

  // histograms
  TH1F* NvtxAll[2];
  TH1F* NvtxBaseline[2];
  TH1F* NvtxBDT[2];

  TH1F* PtAll[2];
  TH1F* PtBaseline[2];
  TH1F* PtBDT[2];

  TH1* PtAllRebinned[2];
  TH1* PtBaselineRebinned[2];
  TH1* PtBDTRebinned[2];

  TGraphAsymmErrors *effVsNvtx_BDT[2];
  TGraphAsymmErrors *effVsNvtx_Baseline[2];

  TGraphAsymmErrors *effVsPt_BDT[2];
  TGraphAsymmErrors *effVsPt_Baseline[2];

  TGraphAsymmErrors *effVsEta_BDT[2];
  TGraphAsymmErrors *effVsEta_Baseline[2];

  int mycolor, mystyle;
  char hname[100];

  for (int i = 0; i < 2; i++){
    
    NvtxAll[i]    = (TH1F*) f[i]->Get("NvtAll");
    NvtxBaseline[i] = (TH1F*) f[i]->Get("NvtGood_RANK");
    NvtxBDT[i]    = (TH1F*) f[i]->Get("NvtGood_BDT");

    PtAll[i]      = (TH1F*) f[i]->Get("PtAll");
    PtBaseline[i]   = (TH1F*) f[i]->Get("PtGood_RANK");
    PtBDT[i]      = (TH1F*) f[i]->Get("PtGood_BDT");

    PtAll[i]      -> Rebin(nRePt);
    PtBaseline[i]   -> Rebin(nRePt);
    PtBDT[i]      -> Rebin(nRePt);

    effVsNvtx_BDT[i]   = new TGraphAsymmErrors();
    effVsNvtx_Baseline[i]= new TGraphAsymmErrors();
    effVsPt_BDT[i]     = new TGraphAsymmErrors();
    effVsPt_Baseline[i]  = new TGraphAsymmErrors();
  
    if (i==0) {mycolor = kGreen+1;}
    else {mycolor = kRed;}
    
    effVsNvtx_BDT[i]->SetLineColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerStyle(24);

    effVsPt_BDT[i]->SetLineColor(mycolor);
    effVsPt_BDT[i]->SetMarkerColor(mycolor);
    effVsPt_BDT[i]->SetMarkerStyle(24);

    effVsNvtx_Baseline[i]->SetLineColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerStyle(20);

    effVsPt_Baseline[i]->SetLineColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerStyle(20);

  }

  
  for (int i = 0; i<2; i++){

    effVsNvtx_Baseline[i]->BayesDivide(NvtxBaseline[i],NvtxAll[i], "cp"); 
    effVsNvtx_BDT[i]   ->BayesDivide(NvtxBDT[i],NvtxAll[i], "cp");
    
    effVsPt_Baseline[i]->BayesDivide(PtBaseline[i],PtAll[i], "cp"); 
    effVsPt_BDT[i]   ->BayesDivide(PtBDT[i],PtAll[i], "cp"); 
        
  }

  int nMBMax = 30;
  
  TLegend legend1(0.68, 0.78, 0.99, 0.99);
  legend1.SetFillColor(kWhite);
  legend1.SetBorderSize(1);
  legend1.AddEntry(effVsNvtx_Baseline[0],legtitle1.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_Baseline[1],legtitle2.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[0],legtitle1bdt.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[1],legtitle2bdt.c_str(),"LP");

  //*** EFF vs NVTX
  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.,1.1);
  cc.SetStats(0); 
  cc.GetXaxis()->SetTitle("number of reconstructed vertices"); 
  //cc.GetXaxis()->SetTitle("number of simulated PU vertices"); 
  cc.GetYaxis()->SetTitle("fraction of events");
  cc.Draw();
  for (int i = 0; i< 2; i++){
    effVsNvtx_Baseline[i]->Draw("p,same");
    effVsNvtx_BDT[i]->Draw("p,same");
  }
  legend1.Draw("same");


  //*** EFF vs BOSON PT  
  TCanvas c2;
  c2.SetGridx();
  c2.SetGridy();
  TH2F dd("dd","",200,0,200,1000,0.0,1.1);
  dd.SetStats(0); 
  dd.GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)"); 
  dd.GetYaxis()->SetTitle("fraction of events");
  dd.Draw();
  for (int i = 0; i< 2; i++){
    effVsPt_Baseline[i]->Draw("p,same");
    effVsPt_BDT[i]->Draw("p,same");
  }
  legend1.Draw("same");
    
  float ptlow = 0.;
  float pthigh = 400.;

  int bin1 = PtBaseline[0]->FindBin(ptlow);
  int bin2 = PtBaseline[0]->FindBin(pthigh);

   
  double eff1 = PtBaseline[0]->Integral(bin1,bin2)/ PtAll[0]->Integral(bin1,bin2);
  double eff2 = PtBaseline[1]->Integral(bin1,bin2)/ PtAll[1]->Integral(bin1,bin2);
  
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
  
  //EFFICIENCY RATIOs

  //------ RATIO EFF vs PT -------------------------------------------
  TGraphErrors * ratioEffVsPt_Baseline = new TGraphErrors();
  ratioEffVsPt_Baseline ->SetMarkerStyle(20);
  
  int nbins = effVsPt_Baseline[0]->GetHistogram()->GetNbinsX();
  
  double pt, err1, err2, ratio , err;
  
  for (int i = 0; i < nbins; i++){
    effVsPt_Baseline[0]->GetPoint(i,pt,eff1);
    effVsPt_Baseline[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsPt_Baseline[0]->GetErrorY(i);
    err2 =  effVsPt_Baseline[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsPt_Baseline->SetPoint(i,pt,ratio);
    ratioEffVsPt_Baseline->SetPointError(i, effVsPt_Baseline[0]->GetErrorX(i) ,err);
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

  ratioEffVsPt_Baseline->GetHistogram()->GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)");
  ratioEffVsPt_Baseline->GetHistogram()->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsPt_Baseline->GetHistogram()->GetXaxis()->SetRangeUser(0,200);
  ratioEffVsPt_Baseline->GetHistogram()->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsPt_Baseline->Draw("ap");
  ratioEffVsPt_BDT->Draw("psame");

  TLegend legend4(0.15, 0.2, 0.45, 0.4);
  legend4.SetFillColor(kWhite);
  legend4.SetBorderSize(1);
  legend4.AddEntry(ratioEffVsPt_Baseline,"RANK","LP");
  legend4.AddEntry(ratioEffVsPt_BDT,"BDT","LP");
  legend4.Draw("same");



  //------ RATIO EFF vs NVTX -------------------------------------------
  TGraphErrors * ratioEffVsNvtx_Baseline = new TGraphErrors();
  ratioEffVsNvtx_Baseline ->SetMarkerStyle(20);
  
  int nbins = effVsNvtx_Baseline[0]->GetHistogram()->GetNbinsX();
  
  double pt, err1, err2, ratio , err;
  
  for (int i = 0; i < nbins; i++){
    effVsNvtx_Baseline[0]->GetPoint(i,pt,eff1);
    effVsNvtx_Baseline[1]->GetPoint(i,pt,eff2);
    
    err1 =  effVsNvtx_Baseline[0]->GetErrorY(i);
    err2 =  effVsNvtx_Baseline[1]->GetErrorY(i);
    
    if (eff2!=0 && eff1!=0)   {
      ratio = eff1/eff2;
      err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); }
    else {
      ratio = 0;
      err = 0;
    }
    
    ratioEffVsNvtx_Baseline->SetPoint(i,pt,ratio);
    ratioEffVsNvtx_Baseline->SetPointError(i,effVsNvtx_Baseline[0]->GetErrorX(i),err);
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

  ratioEffVsNvtx_Baseline->GetHistogram()->GetXaxis()->SetTitle("number of reconstructed vertices");
  ratioEffVsNvtx_Baseline->GetHistogram()->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsNvtx_Baseline->GetHistogram()->GetXaxis()->SetRangeUser(0,16);
  ratioEffVsNvtx_Baseline->GetHistogram()->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsNvtx_Baseline->Draw("ap");
  ratioEffVsNvtx_BDT->Draw("psame");


  TLegend legend5(0.15, 0.2, 0.45, 0.4);
  legend5.SetFillColor(kWhite);
  legend5.SetBorderSize(1);
  legend5.AddEntry(ratioEffVsNvtx_Baseline,"RANK","LP");
  legend5.AddEntry(ratioEffVsNvtx_BDT,"BDT","LP");

  legend5.Draw("same");

  if (saveScaleFactors){
    //    TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_Summer11_Puweights160404-167151.root","recreate");
    //TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_Summer11_Puweights_lp.root","recreate");
    //TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_Summer11_Puweights_eps.root","recreate");
    TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_Summer11_Puweights_after_eps.root","recreate");
    ratioEffVsPt_BDT->SetTitle("scaleFactor");
    ratioEffVsPt_BDT->Write("scaleFactor");
  }

  
}
