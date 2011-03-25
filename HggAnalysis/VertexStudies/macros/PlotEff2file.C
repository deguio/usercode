{
  gROOT->SetStyle("Plain");
  gStyle->SetMarkerSize(1.0);

  string legtitle1 = "Hgg (MC matching + R9)";
  string legtitle1bdt = "Hgg ";
  
  string legtitle2 = "Hgg (MC matching + R9 + mass cut)";
  string legtitle2bdt = "Zee PU2010";



  TH1F* NvtAll_1 = (TH1F*) _file0->Get("NvtAll");
  TH1F* NvtPt2_1 = (TH1F*) _file0->Get("NvtGood");
  TH1F* NvtBDT_1 = (TH1F*) _file0->Get("NvtGood_BDT");

  TH1F* PtAll_1 = (TH1F*) _file0->Get("PtAll");
  TH1F* PtPt2_1 = (TH1F*) _file0->Get("PtGood");
  TH1F* PtBDT_1 = (TH1F*) _file0->Get("PtGood_BDT");

//   TH1F* EtaAll_1 = (TH1F*) _file0->Get("EtaAll");
//   TH1F* EtaPt2_1 = (TH1F*) _file0->Get("EtaGood");
//   TH1F* EtaBDT_1 = (TH1F*) _file0->Get("EtaGood_BDT");

  TH1F* EtaAll_1 = (TH1F*) _file0->Get("PtAll");
  TH1F* EtaPt2_1 = (TH1F*) _file0->Get("PtGood");
  TH1F* EtaBDT_1 = (TH1F*) _file0->Get("PtGood_BDT");
  
  TH2F* PtAll_vsNvtAll_1 = (TH2F*)_file0->Get("PtAll_vsNvtAll");
  TH2F* PtPt2_vsNvtPt2_1 = (TH2F*)_file0->Get("PtGood_vsNvtGood");
  TH2F* PtBDT_vsNvtBDT_1 = (TH2F*)_file0->Get("PtBDT_vsNvtBDT");

  TGraphAsymmErrors effVsNv_BDT_1;
  TGraphAsymmErrors effVsNv_Pt_1;

  NvtPt2_1->Sumw2();  
  NvtBDT_1->Sumw2(); 
  NvtAll_1->Sumw2();
  
  effVsNv_BDT_1.BayesDivide(NvtBDT_1,NvtAll_1, "cp"); 
  effVsNv_BDT_1.SetLineColor(kGreen+1); 
  effVsNv_BDT_1.SetMarkerColor(kGreen+1); 
  effVsNv_BDT_1.SetMarkerStyle(24); 
  
  effVsNv_Pt_1.BayesDivide(NvtPt2_1,NvtAll_1, "cp"); 
  effVsNv_Pt_1.SetLineColor(kGreen+1); 
  effVsNv_Pt_1.SetMarkerColor(kGreen+1); 
  effVsNv_Pt_1.SetMarkerStyle(20); 

  TGraphAsymmErrors effVsPt_Pt_1;
  TGraphAsymmErrors effVsPt_BDT_1;

  PtPt2_1->Sumw2();  
  PtBDT_1->Sumw2(); 
  PtAll_1->Sumw2();

  effVsPt_BDT_1.BayesDivide(PtBDT_1,PtAll_1,"cp"); 
  
  effVsPt_BDT_1.SetLineColor(kGreen+1); 
  effVsPt_BDT_1.SetMarkerColor(kGreen+1); 
  effVsPt_BDT_1.SetMarkerStyle(24); 
  
  effVsPt_Pt_1.BayesDivide(PtPt2_1,PtAll_1,"cp");
  effVsPt_Pt_1.SetLineColor(kGreen+1); 
  effVsPt_Pt_1.SetMarkerColor(kGreen+1); 
  effVsPt_Pt_1.SetMarkerStyle(20); 

  
  TGraphAsymmErrors effVsEta_Pt_1;
  TGraphAsymmErrors effVsEta_BDT_1;

  EtaPt2_1->Sumw2();  
  EtaBDT_1->Sumw2(); 
  EtaAll_1->Sumw2();

  effVsEta_BDT_1.BayesDivide(EtaBDT_1,EtaAll_1,"cp"); 
  effVsEta_BDT_1.SetLineColor(kGreen+1); 
  effVsEta_BDT_1.SetMarkerColor(kGreen+1); 
  effVsEta_BDT_1.SetMarkerStyle(24); 

  effVsEta_Pt_1.BayesDivide(EtaPt2_1,EtaAll_1,"cp");
  effVsEta_Pt_1.SetLineColor(kGreen+1); 
  effVsEta_Pt_1.SetMarkerColor(kGreen+1); 
  effVsEta_Pt_1.SetMarkerStyle(20); 

  TH2F* eff2D_Pt_1 = new TH2F("eff2D_Pt_1","eff2D_Pt_1",
			      PtPt2_vsNvtPt2_1->GetNbinsX(),
			      PtPt2_vsNvtPt2_1->GetXaxis()->GetXmin(),
			      PtPt2_vsNvtPt2_1->GetXaxis()->GetXmax(),
			      PtPt2_vsNvtPt2_1->GetNbinsY(),
			      PtPt2_vsNvtPt2_1->GetYaxis()->GetXmin(),
			      PtPt2_vsNvtPt2_1->GetYaxis()->GetXmax() );
  eff2D_Pt_1 -> Divide(PtPt2_vsNvtPt2_1,PtAll_vsNvtAll_1);
  
  TH2F* eff2D_BDT_1 = new TH2F("eff2D_BDT_1","eff2D_BDT_1",
			       PtBDT_vsNvtBDT_1->GetNbinsX(),
			       PtBDT_vsNvtBDT_1->GetXaxis()->GetXmin(),
			       PtBDT_vsNvtBDT_1->GetXaxis()->GetXmax(),
			       PtBDT_vsNvtBDT_1->GetNbinsY(),
			       PtBDT_vsNvtBDT_1->GetYaxis()->GetXmin(),
			       PtBDT_vsNvtBDT_1->GetYaxis()->GetXmax() );
  eff2D_BDT_1 -> Divide(PtBDT_vsNvtBDT_1,PtAll_vsNvtAll_1);


  // second sample--------------

  TH1F* NvtAll_2 = (TH1F*) _file1->Get("NvtAll");
  TH1F* NvtPt2_2 = (TH1F*) _file1->Get("NvtGood");
  TH1F* NvtBDT_2 = (TH1F*) _file1->Get("NvtGood_BDT");
  
  TH1F* PtAll_2 = (TH1F*) _file1->Get("PtAll");
  TH1F* PtPt2_2 = (TH1F*) _file1->Get("PtGood");
  TH1F* PtBDT_2 = (TH1F*) _file1->Get("PtGood_BDT");

//   TH1F* EtaAll_2 = (TH1F*) _file1->Get("EtaAll");
//   TH1F* EtaPt2_2 = (TH1F*) _file1->Get("EtaGood");
//   TH1F* EtaBDT_2 = (TH1F*) _file1->Get("EtaGood_BDT");
  TH1F* EtaAll_2 = (TH1F*) _file1->Get("PtAll");
  TH1F* EtaPt2_2 = (TH1F*) _file1->Get("PtGood");
  TH1F* EtaBDT_2 = (TH1F*) _file1->Get("PtGood_BDT");

  TH2F* PtAll_vsNvtAll_2 = (TH2F*)_file1->Get("PtAll_vsNvtAll");
  TH2F* PtPt2_vsNvtPt2_2 = (TH2F*)_file1->Get("PtGood_vsNvtGood");
  TH2F* PtBDT_vsNvtBDT_2 = (TH2F*)_file1->Get("PtBDT_vsNvtBDT");

  TGraphAsymmErrors effVsNv_Pt_2;
  TGraphAsymmErrors effVsNv_BDT_2;
  
  NvtPt2_2->Sumw2();  
  NvtBDT_2->Sumw2(); 
  NvtAll_2->Sumw2();
  
  effVsNv_BDT_2.BayesDivide(NvtBDT_2,NvtAll_2, "cp"); 
  effVsNv_BDT_2.SetLineColor(kRed); 
  effVsNv_BDT_2.SetMarkerColor(kRed);
  effVsNv_BDT_2.SetMarkerStyle(24);
   
  effVsNv_Pt_2.BayesDivide(NvtPt2_2,NvtAll_2, "cp"); 
  effVsNv_Pt_2.SetLineColor(kRed); 
  effVsNv_Pt_2.SetMarkerColor(kRed); 
  effVsNv_Pt_2.SetMarkerStyle(20);
  
  TGraphAsymmErrors effVsPt_Pt_2;
  TGraphAsymmErrors effVsPt_BDT_2;
  
  PtPt2_2->Sumw2();  
  PtBDT_2->Sumw2(); 
  PtAll_2->Sumw2();
  
  effVsPt_BDT_2.BayesDivide(PtBDT_2,PtAll_2, "cp"); 
  effVsPt_BDT_2.SetLineColor(kRed); 
  effVsPt_BDT_2.SetMarkerColor(kRed); 
  effVsPt_BDT_2.SetMarkerStyle(24); 
  
  effVsPt_Pt_2.BayesDivide(PtPt2_2,PtAll_2, "cp"); 
  effVsPt_Pt_2.SetLineColor(kRed); 
  effVsPt_Pt_2.SetMarkerColor(kRed); 
  effVsPt_Pt_2.SetMarkerStyle(20); 

  TGraphAsymmErrors effVsEta_Pt_2;
  TGraphAsymmErrors effVsEta_BDT_2;
  
  EtaPt2_2->Sumw2();  
  EtaBDT_2->Sumw2(); 
  EtaAll_2->Sumw2();
  
  effVsEta_BDT_2.BayesDivide(EtaBDT_2,EtaAll_2, "cp"); 
  effVsEta_BDT_2.SetLineColor(kRed); 
  effVsEta_BDT_2.SetMarkerColor(kRed); 
  effVsEta_BDT_2.SetMarkerStyle(24); 
  
  effVsEta_Pt_2.BayesDivide(EtaPt2_2,EtaAll_2, "cp"); 
  effVsEta_Pt_2.SetLineColor(kRed); 
  effVsEta_Pt_2.SetMarkerColor(kRed); 
  effVsEta_Pt_2.SetMarkerStyle(20); 

  TH2F* eff2D_Pt_2 = new TH2F("eff2D_Pt_2","eff2D_Pt_2",
			      PtPt2_vsNvtPt2_2->GetNbinsX(),
			      PtPt2_vsNvtPt2_2->GetXaxis()->GetXmin(),
			      PtPt2_vsNvtPt2_2->GetXaxis()->GetXmax(),
			      PtPt2_vsNvtPt2_2->GetNbinsY(),
			      PtPt2_vsNvtPt2_2->GetYaxis()->GetXmin(),
			      PtPt2_vsNvtPt2_2->GetYaxis()->GetXmax() );
  eff2D_Pt_2 -> Divide(PtPt2_vsNvtPt2_2,PtAll_vsNvtAll_2);

  TH2F* eff2D_BDT_2 = new TH2F("eff2D_BDT_2","eff2D_BDT_1",
			       PtBDT_vsNvtBDT_2->GetNbinsX(),
			       PtBDT_vsNvtBDT_2->GetXaxis()->GetXmin(),
			       PtBDT_vsNvtBDT_2->GetXaxis()->GetXmax(),
			       PtBDT_vsNvtBDT_2->GetNbinsY(),
			       PtBDT_vsNvtBDT_2->GetYaxis()->GetXmin(),
			       PtBDT_vsNvtBDT_2->GetYaxis()->GetXmax() );
  eff2D_BDT_2 -> Divide(PtBDT_vsNvtBDT_2,PtAll_vsNvtAll_2);



  int nMBMax =25;
  
  TLegend legend1(0.68, 0.78, 0.99, 0.99);
  legend1.SetFillColor(kWhite);
  
  TCanvas c1;
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.0,1.1);
  cc.SetStats(0); cc.GetXaxis()->SetTitle("number of PV"); cc.GetYaxis()->SetTitle("efficiency");
  cc.Draw();
  
  effVsNv_Pt_1.Draw("p,same");
  legend1.AddEntry(&effVsNv_Pt_1,legtitle1.c_str(),"LP");
  
  effVsNv_BDT_1.Draw("p,same");
  legend1.AddEntry(&effVsNv_BDT_1, legtitle1bdt.c_str(),"LP");
  
  effVsNv_Pt_2.Draw("p,same");
  legend1.AddEntry(&effVsNv_Pt_2, legtitle2.c_str() ,"LP");
  
  effVsNv_BDT_2.Draw("P,same");
  legend1.AddEntry(&effVsNv_BDT_2, legtitle2bdt.c_str(),"LP"); 
  
  legend1.Draw();
  c1.SetGridx();
  c1.SetGridy();

  TLegend legend2(0.68, 0.78, 0.99, 0.99);
  legend2.SetFillColor(kWhite);
  
  TCanvas c2;
  TH2F dd("dd","",200,0,200,1000,0.0,1.1);
  dd.SetStats(0); dd.GetXaxis()->SetTitle("Pt of the boson "); dd.GetYaxis()->SetTitle("efficiency");
  dd.Draw();

   
  effVsPt_Pt_1.Draw("P,same");
  legend2.AddEntry(&effVsPt_Pt_1,legtitle1.c_str(),"LP");
  
  //   effVsPt_BDT_1.Draw("P,same");
  //    legend2.AddEntry(&effVsPt_BDT_1,legtitle1bdt.c_str(),"LP");
  
  effVsPt_Pt_2.Draw("P,same");
  legend2.AddEntry(&effVsPt_Pt_2, legtitle2.c_str() ,"LP");
  
  //    effVsPt_BDT_2.Draw("P,same");
  //    legend2.AddEntry(&effVsPt_BDT_2,legtitle2bdt.c_str() ,"LP");
  
  legend2.Draw();
  c2.SetGridx();
  c2.SetGridy();
  
  

  TLegend legend3(0.68, 0.78, 0.99, 0.99);
  legend3.SetFillColor(kWhite);
  
  TCanvas c3;
  TH2F ddd("ddd","",60,-6,6,1000,0.0,1.1);
  ddd.SetStats(0); ddd.GetXaxis()->SetTitle("Eta max SC"); ddd.GetYaxis()->SetTitle("efficiency");
  ddd.Draw();
  
  
  effVsEta_Pt_1.Draw("P,same");
  legend3.AddEntry(&effVsEta_Pt_1,legtitle1.c_str(),"LP");
  
  effVsEta_BDT_1.Draw("P,same");
  legend3.AddEntry(&effVsEta_BDT_1,legtitle1bdt.c_str(),"LP");
  
  effVsEta_Pt_2.Draw("P,same");
  legend3.AddEntry(&effVsEta_Pt_2,legtitle2.c_str(),"LP");
  
  effVsEta_BDT_2.Draw("P,same");
  legend3.AddEntry(&effVsEta_BDT_2,legtitle2bdt.c_str(),"LP");
  
  legend3.Draw();
  c3.SetGridx();
  c3.SetGridy();
  
  
  double eff1 = PtPt2_1->GetEntries()/ PtAll_1->GetEntries();
  double eff2 = PtPt2_2->GetEntries()/ PtAll_2->GetEntries();
  
  double eff1bdt = PtBDT_1->GetEntries()/ PtAll_1->GetEntries();
  double eff2bdt = PtBDT_2->GetEntries()/ PtAll_2->GetEntries();
  
  cout << legtitle1.c_str() << " --> eff = " << eff1 << endl;  
  cout << legtitle2.c_str() << " --> eff = " << eff2 << endl;  
  cout << legtitle1bdt.c_str() << " --> eff = " << eff1bdt << endl;  
  cout << legtitle2bdt.c_str() << " --> eff = " << eff2bdt << endl;  
  
  

  TGraphErrors * ratioEff_vs_Pt = new TGraphErrors();
  ratioEff_vs_Pt->SetMarkerStyle(20);
  ratioEff_vs_Pt->SetMarkerSize(0.5);

  int nbins = effVsPt_Pt_1->GetHistogram()->GetNbinsX();
  double pt, err1, err2, ratio , err;
  for (int i = 0; i < nbins; i++){
    effVsPt_Pt_1->GetPoint(i,pt,eff1);
    effVsPt_Pt_2->GetPoint(i,pt,eff2);
 
    err1 =  effVsPt_Pt_1->GetErrorY(i);
    err2 =  effVsPt_Pt_2->GetErrorY(i);
    
    ratio = eff1/eff2;
    err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); 

    ratioEff_vs_Pt->SetPoint(i,pt,ratio);
    ratioEff_vs_Pt->SetPointError(i,0,err);
  }
  

  TGraphErrors * ratioEffBDTBDT_vs_Pt = new TGraphErrors();
  ratioEffBDT_vs_Pt->SetMarkerStyle(20);
  ratioEffBDT_vs_Pt->SetMarkerSize(0.5);

  int nbins = effVsPt_Pt_1->GetHistogram()->GetNbinsX();
  for (int i = 0; i < nbins; i++){
    effVsPt_BDT_1->GetPoint(i,pt,eff1);
    effVsPt_BDT_2->GetPoint(i,pt,eff2);
 
    err1 =  effVsPt_BDT_1->GetErrorY(i);
    err2 =  effVsPt_BDT_2->GetErrorY(i);
    
    ratio = eff1/eff2;
    err   = ratio * sqrt( pow( err1/eff1 ,2) + pow( err2/eff2,2) ); 

    ratioEffBDT_vs_Pt->SetPoint(i,pt,ratio);
    ratioEffBDT_vs_Pt->SetPointError(i,0,err);


  }
    
  
  TCanvas cRatioBDT;
  cRatioBDT.SetGridx();
  cRatioBDT.SetGridy();

  ratioEffBDT_vs_Pt->GetHistogram()->GetXaxis()->SetTitle("boson Pt");
  ratioEffBDT_vs_Pt->GetHistogram()->GetYaxis()->SetTitle("Efficiency (Hgg) / Efficiency (Zee)");
  ratioEffBDT_vs_Pt->Draw("ap");
  

  
  
  //    TCanvas c33;
  //    eff2D_Pt_1 -> SetTitle(legtitle1.c_str());
  //    eff2D_Pt_1 -> Draw("colz");
  //    TCanvas c4;
  //    eff2D_BDT_1 -> SetTitle(legtitle1bdt.c_str());
  //    eff2D_BDT_1 -> Draw("colz");
  //    TCanvas c5;
  //    eff2D_Pt_2 -> SetTitle(legtitle2.c_str() );   
  //    eff2D_Pt_2 -> Draw("colz");
  //    TCanvas c6;
  //    eff2D_BDT_2 -> SetTitle(legtitle2bdt.c_str() );
  //    eff2D_BDT_2 -> Draw("colz");
}
