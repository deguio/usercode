{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();
  gStyle->SetErrorX(0.5);

  int saveScaleFactors    = 0;
  bool useVariableBinning = false;

  // 0 : data
  // 1 : mc

  TFile *f[2];

  //2012ABC
  f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DoubleMu_Run2012ABC_dz01/TMVA_check_DoubleMu_Run2012ABC.root");
  f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABC_dz01/testEfficiency.root");
  //f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABC_dz01_BSrw1/testEfficiency.root");


  TLatex *latex = new TLatex(0.55,0.85,"#splitline{CMS preliminary}{#sqrt{s} = 8 TeV L = 12 fb^{-1}}");
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);

  // legend
  string legtitle1    = "DATA Z#rightarrow#mu#mu (RANK)";
  string legtitle1bdt = "DATA Z#rightarrow#mu#mu ";
  string legtitle2    = "MC Z#rightarrow#mu#mu (RANK)";
  string legtitle2bdt = "MC Z#rightarrow#mu#mu ";

  // rebinning
  int nRePt = 2;

  // if variable size bins
  double xbins[10] = {0.,10.,20.,30.,40.,50.,70.,110.,250.,400.};

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

  TH1F *effVsNvtx_BDT[2];
  TH1F *effVsNvtx_Baseline[2];

  TH1F *effVsPt_BDT[2];
  TH1F *effVsPt_Baseline[2];

  TH1F *effVsEta_BDT[2];
  TH1F *effVsEta_Baseline[2];

  TH1F *BDToutput[2];
  TH1F *BDToutput_sig[2];
  TH1F *BDToutput_bkg[2];

  TH1F *evtBDToutput[2];
  TH1F *evtBDToutput_sig[2];
  TH1F *evtBDToutput_bkg[2];

  int mycolor, mystyle;
  char hname[100];
  std::string evtType[2] = {"data","mc"};
  
  for (int i = 0; i < 2; i++){

    BDToutput_sig[i] = (TH1F*) f[i]->Get("BDToutput_sig");
    BDToutput_bkg[i] = (TH1F*) f[i]->Get("BDToutput_bkg");
    evtBDToutput[i]     = (TH1F*) f[i]->Get("perEventBDToutput");
    evtBDToutput_sig[i] = (TH1F*) f[i]->Get("perEventBDToutput_sig");
    evtBDToutput_bkg[i] = (TH1F*) f[i]->Get("perEventBDToutput_bkg");

    BDToutput_sig[i]->Sumw2();
    BDToutput_bkg[i]->Sumw2();
    evtBDToutput[i]->Sumw2();
    evtBDToutput_sig[i]->Sumw2();
    evtBDToutput_bkg[i]->Sumw2();
    
    BDToutput_sig[i]->Rebin(10);
    BDToutput_bkg[i]->Rebin(10);
    evtBDToutput[i]->Rebin(10);
    evtBDToutput_sig[i]->Rebin(10);
    evtBDToutput_bkg[i]->Rebin(10);

    NvtxAll[i]      = (TH1F*) f[i]->Get("NvtAll");
    NvtxBaseline[i] = (TH1F*) f[i]->Get("NvtGood_RANK");
    NvtxBDT[i]      = (TH1F*) f[i]->Get("NvtGood_BDT");

    PtAll[i]      = (TH1F*) f[i]->Get("PtAll");
    PtBaseline[i] = (TH1F*) f[i]->Get("PtGood_RANK");
    PtBDT[i]      = (TH1F*) f[i]->Get("PtGood_BDT");

    NvtxAll[i]     -> Sumw2();
    NvtxBaseline[i]-> Sumw2();
    NvtxBDT[i]     -> Sumw2();

    PtAll[i]       -> Sumw2();
    PtBaseline[i]  -> Sumw2();
    PtBDT[i]       -> Sumw2();

    sprintf(hname,"PtAllRebinned_%d",i);
    PtAllRebinned[i] = PtAll[i]   ->Rebin(9,hname,xbins);
    sprintf(hname,"PtBaselineRebinned_%d",i);
    PtBaselineRebinned[i] = PtBaseline[i]->Rebin(9,hname,xbins);
    sprintf(hname,"PtBDTRebinned_%d",i);
    PtBDTRebinned[i] = PtBDT[i]   ->Rebin(9,"PtBDTRebinned",xbins);

    sprintf(hname, "effVsNvtx_Baseline_%s",evtType[i].c_str());
    effVsNvtx_Baseline[i] = (TH1F*)NvtxBaseline[i]->Clone(hname);
    sprintf(hname, "effVsNvtx_BDT_%s",evtType[i].c_str());
    effVsNvtx_BDT[i]      = (TH1F*)NvtxBDT[i]->Clone(hname);

    if (useVariableBinning){
      sprintf(hname, "effVsPt_Baseline_%s",evtType[i].c_str());
      effVsPt_Baseline[i]   = (TH1F*)PtBaselineRebinned[i]->Clone(hname);
      sprintf(hname, "effVsPt_BDT_%s",evtType[i].c_str());
      effVsPt_BDT[i]        = (TH1F*)PtBDTRebinned[i]->Clone(hname); 
    }
    else {
      sprintf(hname, "effVsPt_Baseline_%s",evtType[i].c_str());
      effVsPt_Baseline[i]   = (TH1F*)PtBaseline[i]->Clone(hname);
      sprintf(hname, "effVsPt_BDT_%s",evtType[i].c_str());
      effVsPt_BDT[i]        = (TH1F*)PtBDT[i]->Clone(hname);
    }

    if (i==0) {mycolor = kBlue+1;}
    else {mycolor = kRed;}
    
    effVsNvtx_BDT[i]->SetLineColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerColor(mycolor);
    effVsNvtx_BDT[i]->SetMarkerStyle(20);

    effVsPt_BDT[i]->SetLineColor(mycolor);
    effVsPt_BDT[i]->SetMarkerColor(mycolor);
    effVsPt_BDT[i]->SetMarkerStyle(20);

    effVsNvtx_Baseline[i]->SetLineColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerColor(mycolor);
    effVsNvtx_Baseline[i]->SetMarkerStyle(24);

    effVsPt_Baseline[i]->SetLineColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerColor(mycolor);
    effVsPt_Baseline[i]->SetMarkerStyle(24);

  }

  
  for (int i = 0; i<2; i++){
    effVsNvtx_Baseline[i] -> Divide(NvtxBaseline[i],NvtxAll[i],1,1 ,"B"); 
    effVsNvtx_BDT[i]      -> Divide(NvtxBDT[i],NvtxAll[i],1,1, "B");
   
    if (!useVariableBinning){ 
      effVsPt_Baseline[i]   -> Divide(PtBaseline[i],PtAll[i],1,1, "B"); 
      effVsPt_BDT[i]        -> Divide(PtBDT[i],PtAll[i],1,1, "B");  
    }    
    else {
      effVsPt_Baseline[i]   -> Divide(PtBaselineRebinned[i],PtAllRebinned[i],1,1, "B"); 
      effVsPt_BDT[i]        -> Divide(PtBDTRebinned[i],PtAllRebinned[i],1,1, "B");  
    }
  }

  int nMBMax = 40;
  
  TLegend legend1(0.68, 0.18, 0.92, 0.38);
  legend1.SetFillColor(kWhite);
  legend1.SetBorderSize(1);
  //legend1.AddEntry(effVsNvtx_Baseline[0],legtitle1.c_str(),"LP");
  //legend1.AddEntry(effVsNvtx_Baseline[1],legtitle2.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[0],legtitle1bdt.c_str(),"LP");
  legend1.AddEntry(effVsNvtx_BDT[1],legtitle2bdt.c_str(),"LP");

  //*** EFF vs NVTX
  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.,1.1);
  cc.SetStats(0); 
  cc.GetXaxis()->SetTitle("number of reconstructed vertices"); 
  cc.GetYaxis()->SetTitle("fraction of events");
  cc.Draw();
  for (int i = 0; i< 2; i++){
    //    effVsNvtx_Baseline[i]->Draw("e1,same");
    effVsNvtx_BDT[i]->Draw("e1,same");
  }
  legend1.Draw("same");
  latex->Draw("same");

  //*** EFF vs BOSON PT  
  TCanvas c2;
  c2.SetGridx();
  c2.SetGridy();
  TH2F dd("dd","",250,0,250,1000,0.0,1.1);
  dd.SetStats(0); 
  dd.GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)"); 
  dd.GetYaxis()->SetTitle("fraction of events");
  dd.Draw();
  for (int i = 0; i< 2; i++){
    //effVsPt_Baseline[i]->Draw("e1,same");
    effVsPt_BDT[i]->Draw("e1,same");
  }
  legend1.Draw("same");
  latex->Draw("same");
    
  float ptlow = 0.;
  float pthigh = 250.;

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
  TH1F *ratioEffVsPt_Baseline = (TH1F*)effVsPt_Baseline[0]->Clone("ratioEffVsPt_Baseline");
  ratioEffVsPt_Baseline ->Divide(effVsPt_Baseline[1]);

  TH1F *ratioEffVsPt_BDT = (TH1F*)effVsPt_BDT[0]->Clone("ratioEffVsPt_BDT");
  ratioEffVsPt_BDT ->Divide(effVsPt_BDT[1]);
  
  TCanvas cRatioPt("cRatioPt","cRatioPt",600,300);
  cRatioPt.SetGridx();
  cRatioPt.SetGridy();

  ratioEffVsPt_BDT->SetMarkerColor(1);
  ratioEffVsPt_BDT->SetLineColor(1);
  ratioEffVsPt_BDT->GetXaxis()->SetTitle("p_{T}(Z) (GeV/c)");
  ratioEffVsPt_BDT->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsPt_BDT->GetXaxis()->SetRangeUser(0,400);
  ratioEffVsPt_BDT->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsPt_BDT->Draw("e1");
  //ratioEffVsPt_Baseline->Draw("e1same");

  //-- save also in TGraphErrors format
  TGraphErrors *gratioEffVsPt_BDT = new TGraphErrors();
  for (int ibin = 0; ibin < ratioEffVsPt_BDT->GetNbinsX();  ibin++){
    float x  = ratioEffVsPt_BDT->GetBinCenter(ibin+1);
    float ex = ratioEffVsPt_BDT->GetBinWidth(ibin+1)/2;
    float y  = ratioEffVsPt_BDT->GetBinContent(ibin+1);
    float ey = ratioEffVsPt_BDT->GetBinError(ibin+1);
    gratioEffVsPt_BDT->SetPoint(ibin,x,y);
    gratioEffVsPt_BDT->SetPointError(ibin,ex,ey);
  }

  TLegend legend4(0.15, 0.2, 0.45, 0.4);
  legend4.SetFillColor(kWhite);
  legend4.SetBorderSize(1);
  legend4.AddEntry(ratioEffVsPt_Baseline,"RANK","LP");
  legend4.AddEntry(ratioEffVsPt_BDT,"BDT","LP");
  //legend4.Draw("same");



  //------ RATIO EFF vs NVTX -------------------------------------------
  TH1F *ratioEffVsNvtx_Baseline = (TH1F*)effVsNvtx_Baseline[0]->Clone("ratioEffVsNvtx_Baseline");
  ratioEffVsNvtx_Baseline ->Divide(effVsNvtx_Baseline[1]);

  TH1F *ratioEffVsNvtx_BDT = (TH1F*)effVsNvtx_BDT[0]->Clone("ratioEffVsNvtx_BDT");
  ratioEffVsNvtx_BDT ->Divide(effVsNvtx_BDT[1]);

  TCanvas cRatioNvtx("cRatioNvtx","cRatioNvtx",600,300);
  cRatioNvtx.SetGridx();
  cRatioNvtx.SetGridy();

  ratioEffVsNvtx_BDT->GetXaxis()->SetTitle("number of reconstructed vertices");
  ratioEffVsNvtx_BDT->GetYaxis()->SetTitle("#epsilon(data)/#epsilon(MC)");
  ratioEffVsNvtx_BDT->GetXaxis()->SetRangeUser(0,30);
  ratioEffVsNvtx_BDT->GetYaxis()->SetRangeUser(0.7,1.3);
  ratioEffVsNvtx_BDT->Draw("e1");
  //ratioEffVsNvtx_Baseline->Draw("e1");


  TLegend legend5(0.15, 0.2, 0.45, 0.4);
  legend5.SetFillColor(kWhite);
  legend5.SetBorderSize(1);
  legend5.AddEntry(ratioEffVsNvtx_Baseline,"RANK","LP");
  legend5.AddEntry(ratioEffVsNvtx_BDT,"BDT","LP");
  //legend5.Draw("same");



  // NVTX control plot
  TCanvas *cNvtx = new TCanvas("cNvtx","cNvtx",500,500);
  NvtxAll[1]->SetFillColor(kRed);
  NvtxAll[1]->SetFillStyle(3004);
  NvtxAll[1]->GetXaxis()->SetTitle("Number of vertices");
  NvtxAll[1]->DrawNormalized("histo");
  NvtxAll[0]->DrawNormalized("esame");
  latex->Draw("same");
  float p = NvtxAll[0]->KolmogorovTest(NvtAll[1],"");
  cout << "DATA < N_vtx > = " << NvtxAll[0]->GetMean() << endl;
  cout << "MC   < N_vtx > = " << NvtxAll[1]->GetMean() << endl;
  cout << "NVTX Kolmogorov : p = " << p << endl;
  TLegend leg (0.6, 0.6,0.89,0.75);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(NvtxAll[1],"MC Z#rightarrow#mu#mu","F");
  leg->AddEntry(NvtxAll[0],"Data Z#rightarrow#mu#mu","LP");
  leg->Draw("same");



  //----------------- MVA control plots
  // per vertex mva
  BDToutput_sig[0] -> SetMarkerStyle(20);
  BDToutput_sig[0] -> SetMarkerSize(0.8);
  BDToutput_sig[0] -> SetMarkerColor(kGreen+2);
  BDToutput_bkg[0] -> SetMarkerStyle(20);
  BDToutput_bkg[0] -> SetMarkerSize(0.8);
  BDToutput_bkg[0] -> SetMarkerColor(kRed+2);
  BDToutput_sig[1] -> SetFillColor(kGreen+1);
  BDToutput_sig[1] -> SetFillStyle(3002);
  BDToutput_bkg[1] -> SetFillColor(kRed+1);
  BDToutput_bkg[1] -> SetFillStyle(3005);

  TLegend legVtxMva (0.6, 0.6,0.89,0.75);
  legVtxMva->SetFillColor(0);
  legVtxMva->SetBorderSize(0);
  legVtxMva->AddEntry(evtBDToutput_sig[1],"right vertex MC","F");
  legVtxMva->AddEntry(evtBDToutput_bkg[1],"wrong vertex MC","F");
  legVtxMva->AddEntry(evtBDToutput_sig[0],"right vertex DATA","LP");
  legVtxMva->AddEntry(evtBDToutput_bkg[0],"wrong vertex DATA","LP");

  TCanvas *cVertexMva = new TCanvas("cVertexMva","cVertexMva",500,500);
  BDToutput_bkg[1] -> GetXaxis() -> SetTitle("MVA_{vtx}");
  BDToutput_bkg[1] -> DrawNormalized("histo");
  BDToutput_sig[1] -> DrawNormalized("histo same");
  BDToutput_bkg[0] -> DrawNormalized("esame");
  BDToutput_sig[0] -> DrawNormalized("esame");
  legVtxMva->Draw("same");
  latex->Draw("same");

  // per event mva
  evtBDToutput_sig[0] -> SetMarkerStyle(20);
  evtBDToutput_sig[0] -> SetMarkerSize(0.8);
  evtBDToutput_sig[0] -> SetMarkerColor(kGreen+2);
  evtBDToutput_bkg[0] -> SetMarkerStyle(20);
  evtBDToutput_bkg[0] -> SetMarkerSize(0.8);
  evtBDToutput_bkg[0] -> SetMarkerColor(kRed+2);
  evtBDToutput_sig[1] -> SetFillColor(kGreen+1);
  evtBDToutput_sig[1] -> SetFillStyle(3002);
  evtBDToutput_bkg[1] -> SetFillColor(kRed+1);
  evtBDToutput_bkg[1] -> SetFillStyle(3005);

  TLegend legEvtMva (0.6, 0.6,0.89,0.75);
  legEvtMva->SetFillColor(0);
  legEvtMva->SetBorderSize(0);
  legEvtMva->AddEntry(evtBDToutput_sig[1],"right vertex MC","F");
  legEvtMva->AddEntry(evtBDToutput_bkg[1],"wrong vertex MC","F");
  legEvtMva->AddEntry(evtBDToutput_sig[0],"right vertex DATA","LP");
  legEvtMva->AddEntry(evtBDToutput_bkg[0],"wrong vertex DATA","LP");

  TCanvas *cEventMva = new TCanvas("cEventMva","cEventMva",500,500);
  evtBDToutput_sig[1] -> GetXaxis() -> SetTitle("MVA_{event}");
  evtBDToutput_sig[1] -> DrawNormalized("histo");
  evtBDToutput_bkg[1] -> DrawNormalized("histo same");
  evtBDToutput_sig[0] -> DrawNormalized("esame");
  evtBDToutput_bkg[0] -> DrawNormalized("esame");
  legEvtMva->Draw("same");
  latex->Draw("same");

  // vertex probability
  TH1F *hVertexProbability[2];
  hVertexProbability[0]= (TH1F*)evtBDToutput_sig[0]->Clone("hVertexProbability");
  hVertexProbability[1]= (TH1F*)evtBDToutput_sig[1]->Clone("hVertexProbability");
  hVertexProbability[0]->Divide(evtBDToutput[0]);
  hVertexProbability[1]->Divide(evtBDToutput[1]);
  hVertexProbability[0]->SetMarkerColor(kBlue+1);
  hVertexProbability[1]->SetMarkerColor(kRed+1);
  TF1 *fprob = new TF1("fprob","1+[0]*(x+1)",-1,1);
  fprob->SetParameter(0,-0.49);
  fprob->SetLineColor(kGray+1);
  TCanvas *cProbability = new TCanvas("cProbability","cProbability",500,300);
  cProbability ->SetGridx();
  cProbability ->SetGridy();
  hVertexProbability[0]->GetXaxis()->SetTitle("MVA");
  hVertexProbability[0]->GetYaxis()->SetTitle("probability");
  hVertexProbability[0]->Draw();
  hVertexProbability[1]->Draw("same");
  fprob->Draw("same");
  legend1.Draw("same");

  if (saveScaleFactors){
    TFile *fileout = new TFile("vtxIdScaleFactorFromZmumu_dz01_BSrw1.root","recreate");
    ratioEffVsPt_BDT->SetTitle("hscaleFactor");
    ratioEffVsPt_BDT->Write("hscaleFactor");
    gratioEffVsPt_BDT->SetTitle("scaleFactor");
    gratioEffVsPt_BDT->Write("scaleFactor");
    fileout->Close();
  }

  cout << "ciao ciao " << endl;

  
}
