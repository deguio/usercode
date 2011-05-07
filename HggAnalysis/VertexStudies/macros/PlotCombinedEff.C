{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  //TFile *f = TFile::Open("output/vtxCombined_Hgg-M120_PU_AVE6.root");
  TFile *f = TFile::Open("output/vtxCombined_Hgg-M120_Spring11-PU_S1_1sigmaz.root");
    
  string legtitle1 = "AVE PU 6";
  string legtitle2 = "AVE PU 12";


  TH1F *NvtAll[4];
  TH1F *NvtGood[4];

  TH1F *PtAll[4];
  TH1F *PtGood[4];
  
  TH1F *R9All[4];
  TH1F *R9Good[4];
  

  TGraphAsymmErrors *EffVsNvtx[4];
  TGraphAsymmErrors *EffVsPt[4];
  TGraphAsymmErrors *EffVsR9[4];

  char hname[100];

  for (int i = 0; i< 4; i++){
   
    sprintf(hname,"NvtAll_cat%d",i);
    NvtAll[i] = (TH1F*)f->Get(hname);
  
    sprintf(hname,"NvtGood_cat%d",i);
    NvtGood[i] = (TH1F*)f->Get(hname);

    sprintf(hname,"PtAll_cat%d",i);
    PtAll[i] = (TH1F*)f->Get(hname);
    sprintf(hname,"PtGood_cat%d",i);
    PtGood[i] = (TH1F*)f->Get(hname);

    sprintf(hname,"R9All_cat%d",i);
    R9All[i] = (TH1F*)f->Get(hname);
    sprintf(hname,"R9Good_cat%d",i);
    R9Good[i] = (TH1F*)f->Get(hname);

    EffVsNvtx[i] = new TGraphAsymmErrors();
    EffVsPt[i] = new TGraphAsymmErrors();
    EffVsR9[i] = new TGraphAsymmErrors();

    EffVsNvtx[i]-> BayesDivide(NvtGood[i], NvtAll[i], "cp");
    EffVsPt[i]  -> BayesDivide(PtGood[i], PtAll[i], "cp");
    EffVsR9[i]  -> BayesDivide(R9Good[i], R9All[i], "cp");
 
    EffVsNvtx[i]->SetLineColor(i+1);
    EffVsNvtx[i]->SetMarkerColor(i+1);

    EffVsPt[i]->SetLineColor(i+1);
    EffVsPt[i]->SetMarkerColor(i+1);

    EffVsR9[i]->SetLineColor(i+1);
    EffVsR9[i]->SetMarkerColor(i+1);
      
  }
    
  
  float ptlow = 0.;
  float pthigh = 400.;

  int bin1 = PtAll[0]->FindBin(ptlow);
  int bin2 = PtAll[0]->FindBin(pthigh);

  double eff[4];
  double err[4];


  for (int i = 0; i < 4 ; i++){
    
    eff[i] = PtGood[i]->Integral(bin1,bin2)/ PtAll[i]->Integral(bin1,bin2);
    err[i] = sqrt(eff[i]*(1-eff[i])/PtAll[i]->Integral(bin1,bin2));
  }
  
  cout << "Efficiency integrated in boson pt : [" << ptlow << ","<< pthigh<<"] GeV"<< endl;
  cout << "Global (no categorories): eff = " << eff[0] << " +/- " << err[0]<< endl;  
  cout << "Category 1              : eff = " << eff[1] << " +/- " << err[1]<< endl;  
  cout << "Category 2              : eff = " << eff[2] << " +/- " << err[2]<< endl;  
  cout << "Category 3              : eff = " << eff[3] << " +/- " << err[3]<< endl;  
  
  
   
  TLegend leg(0.68, 0.78, 0.99, 0.99);
  leg.SetFillColor(kWhite);
  leg.SetBorderSize(1);
  leg.AddEntry(EffVsNvtx[1],"cat1","LP");
  leg.AddEntry(EffVsNvtx[2],"cat2","LP");
  leg.AddEntry(EffVsNvtx[3],"cat3","LP");

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetGridx();
  c1->SetGridy();
  EffVsNvtx[1]->GetXaxis()->SetRangeUser(1,15);
  EffVsNvtx[1]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsNvtx[1]->GetXaxis()->SetTitle("number of PV"); 
  EffVsNvtx[1]->GetYaxis()->SetTitle("efficiency"); 
  EffVsNvtx[1]->Draw("ap");
  EffVsNvtx[2]->Draw("psame");
  EffVsNvtx[3]->Draw("psame");
  leg.Draw("same");


  TCanvas *c2 = new TCanvas("c2","c2");
  c2->SetGridx();
  c2->SetGridy();
  EffVsPt[1]->GetXaxis()->SetRangeUser(0,200);
  EffVsPt[1]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsPt[1]->GetXaxis()->SetTitle("boson p^{T} (GeV)"); 
  EffVsPt[1]->GetYaxis()->SetTitle("efficiency"); 
  EffVsPt[1]->Draw("ap");
  EffVsPt[2]->Draw("psame");
  EffVsPt[3]->Draw("psame");
  leg.Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetGridx();
  c3->SetGridy();
  EffVsR9[2]->GetXaxis()->SetRangeUser(0.,1.);
  EffVsR9[2]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsR9[2]->GetXaxis()->SetTitle("max photon R9"); 
  EffVsR9[2]->GetYaxis()->SetTitle("efficiency"); 
  EffVsR9[2]->Draw("ap");
  EffVsR9[1]->Draw("psame");
  EffVsR9[3]->Draw("psame");
  leg.Draw("same");

    
  TCanvas *c1all = new TCanvas("c1all","c1all");
  c1all->SetGridx();
  c1all->SetGridy();
  EffVsNvtx[0]->GetXaxis()->SetRangeUser(1,15);
  EffVsNvtx[0]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsNvtx[0]->GetXaxis()->SetTitle("number of PV"); 
  EffVsNvtx[0]->GetYaxis()->SetTitle("efficiency"); 
  EffVsNvtx[0]->Draw("ap");
  
  TCanvas *c2all = new TCanvas("c2all","c2all");
  c2all->SetGridx();
  c2all->SetGridy();
  EffVsPt[0]->GetXaxis()->SetRangeUser(0,200);
  EffVsPt[0]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsPt[0]->GetXaxis()->SetTitle("boson p^{T} (GeV)"); 
  EffVsPt[0]->GetYaxis()->SetTitle("efficiency"); 
  EffVsPt[0]->Draw("ap");
  
  TCanvas *c3all = new TCanvas("c3all","c3all");
  c3all->SetGridx();
  c3all->SetGridy();
  EffVsR9[0]->GetXaxis()->SetRangeUser(0,1);
  EffVsR9[0]->GetYaxis()->SetRangeUser(0.3,1.1);
  EffVsR9[0]->GetXaxis()->SetTitle("max photon R9"); 
  EffVsR9[0]->GetYaxis()->SetTitle("efficiency"); 
  EffVsR9[0]->Draw("ap");
}
