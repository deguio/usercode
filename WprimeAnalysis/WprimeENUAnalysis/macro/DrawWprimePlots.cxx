
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleBorderSize(0);

   float lumi = 4.375;
  
  //float lumi = 1.372;  //2010B
  //float lumi = 0.119; // jul16
  //float lumi = 2.884; // 

  int nRe = 2;

  TFile *fdata = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Data_7TeV.root");
  //TFile *fdata = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Data_Run2010B_7TeV.root");
  //TFile *fdata = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Data_Run137437-144114_7TeV.root");
  //TFile *fdata = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Data_Run139559-140159_7TeV.root");
  hetData = (TH1F*)fdata->Get("het");
  hmetData= (TH1F*)fdata->Get("hmet");
  hmtData = (TH1F*)fdata->Get("hmt");
  hetData->Rebin(nRe);
  hmetData->Rebin(nRe);
  hmtData->Rebin(nRe);

  hetData ->SetMarkerStyle(20);
  hmetData ->SetMarkerStyle(20);
  hmtData ->SetMarkerStyle(20);




  TFile *f[10];
  int N = 9;
  f[0] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_WprimeM800_7TeV.root");
  f[1] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_WprimeM1000_7TeV.root");
  f[2] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_WprimeM1500_7TeV.root");
  f[3] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_QCD_7TeV.root");
  f[4] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Zee_7TeV.root");
  f[5] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Dibosons_7TeV.root");
  f[6] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Wtaunu_7TeV.root");
  f[7] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_TTbar_7TeV.root");
  f[8] = TFile::Open("HistosForTTbarBkg_HighPurBdisc/outputHistos_Wenu_7TeV.root");
   

  TH1F *het[10];
  TH1F *hmet[10];
  TH1F *hmt[10];
 
  TH1F *hetOmet[10];
  TH1F *hdphi[10];
  
  

  int mycolor = 1;

  for (int i = 0; i < N; i++){
    het[i] = (TH1F*)f[i]->Get("het");
    hmet[i]= (TH1F*)f[i]->Get("hmet");
    hmt[i] = (TH1F*)f[i]->Get("hmt");

    het[i] ->Scale(lumi);
    hmet[i]->Scale(lumi);
    hmt[i] ->Scale(lumi);
   

    het[i] ->Rebin(nRe);
    hmet[i]->Rebin(nRe);
    hmt[i] ->Rebin(nRe);
   
    het[i]->SetLineWidth(2);
    hmet[i]->SetLineWidth(2);
    hmt[i]->SetLineWidth(2);
   
    het[i] ->SetFillStyle(3001);
    hmt[i] ->SetFillStyle(3001);
    hmet[i] ->SetFillStyle(3001);
   
    if (i==0) mycolor = 2; //W'
    if (i==1) mycolor = 3; //W'
    if (i==2) mycolor = 4; //W'
    if (i==3) mycolor = 1; // QCD
    if (i==4) mycolor = 5; // Zee
    if (i==5) mycolor = 8; // dibosons 
    if (i==6) mycolor = 7; // Wtaunu
    if (i==7) mycolor = 6; // ttbar
    if (i==8) mycolor = 9; // Wenu
 
    het[i] ->SetFillColor(mycolor);
    hmet[i]->SetFillColor(mycolor);
    hmt[i] ->SetFillColor(mycolor);
    het[i] ->SetLineColor(mycolor);
    hmet[i]->SetLineColor(mycolor);
    hmt[i] ->SetLineColor(mycolor);
      
//     hetOmet[i] = (TH1F*)f[i]->Get("hetOmet");
//     //hetOmet[i] -> Rebin(5);
//     hetOmet[i] -> SetLineWidth(2);
//     hetOmet[i] -> SetLineColor(mycolor);

//     hdphi[i] = (TH1F*)f[i]->Get("hdphi");
//     //hdphi[i] -> Rebin(5);
//     hdphi[i] -> SetLineWidth(2);
//     hdphi[i] -> SetLineColor(mycolor);
  }

   
  int bin1 = het[0]->FindBin(200);
  int bin2 = het[0]->FindBin(1500);
  cout << "Number of events with ET > 100: "<< endl;
  cout << "W' M= 800 : " << het[0]->Integral(bin1,bin2) << endl;
  cout << "W' M=1000 : " << het[1]->Integral(bin1,bin2) << endl;
  cout << "W' M=1500 : " << het[2]->Integral(bin1,bin2) << endl;
  cout << "W         : " << het[8]->Integral(bin1,bin2) << endl;
  cout << "TTbar     : " << het[7]->Integral(bin1,bin2) << endl;
  cout << "Zee       : " << het[4]->Integral(bin1,bin2) << endl;
  cout << "Wtaunu    : " << het[6]->Integral(bin1,bin2) << endl;
  cout << "QCD       : " << het[3]->Integral(bin1,bin2) << endl;
  cout << "dibosons  : " << het[5]->Integral(bin1,bin2) << endl;

  bin1 = hmt[0]->FindBin(150);
  bin2 = hmt[0]->FindBin(3000);
  cout << "Number of events with MT > 150: "<< endl;
  cout << "W' M= 800 : " << hmt[0]->Integral(bin1,bin2) << endl;
  cout << "W' M=1000 : " << hmt[1]->Integral(bin1,bin2) << endl;
  cout << "W' M=1500 : " << hmt[2]->Integral(bin1,bin2) << endl;
  cout << "W         : " << hmt[8]->Integral(bin1,bin2) << endl;
  cout << "TTbar     : " << hmt[7]->Integral(bin1,bin2) << endl;
  cout << "Zee       : " << hmt[4]->Integral(bin1,bin2) << endl;
  cout << "Wtaunu    : " << hmt[6]->Integral(bin1,bin2) << endl;
  cout << "QCD       : " << hmt[3]->Integral(bin1,bin2) << endl;
  cout << "dibosons  : " << hmt[5]->Integral(bin1,bin2) << endl;
  cout << "**** DATA : " << hmtData->Integral(bin1,bin2) << endl;

  TLegend *leg = new TLegend (0.4,0.6,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(het[0],"W', M =  800 GeV","F");
  //leg->AddEntry(het[1],"W', M = 1000 GeV","F");
  //leg->AddEntry(het[2],"W', M = 1500 GeV","F");
  leg->AddEntry(het[3],"QCD ","F");
  leg->AddEntry(het[4],"Z #rightarrow ee","F");
  leg->AddEntry(het[5],"Dibosons","F");
  leg->AddEntry(het[6],"W #rightarrow #tau #nu","F");
  leg->AddEntry(het[7],"t#bar{t}","F");
  leg->AddEntry(het[8],"W #rightarrow e #nu","F");


  TLatex latex;
  latex.SetTextSize(0.035);

  THStack *hetall = new THStack("hetall","hetall");
  hetall->Add(het[5]);
  hetall->Add(het[7]);
  hetall->Add(het[4]);
  hetall->Add(het[6]);
  hetall->Add(het[3]);
  hetall->Add(het[8]);

  
  TCanvas *cet = new TCanvas("cet","cet",600,600);
  cet->SetLogy();
  hetall->Draw("a histo");
  hetall->SetMinimum(0.01);
  hetall->SetMaximum(10000);
  hetall->GetHistogram()->GetXaxis()->SetRangeUser(0.,500);
  hetall->GetHistogram()->GetXaxis()->SetTitle("electron transverse energy (GeV)");
  hetall->GetHistogram()->GetYaxis()->SetTitle("events/5 GeV");
  het[0] ->SetFillStyle(0);
  het[0]->Draw("histo same");
  hetData->Draw("esame");
  leg->Draw("same");

  THStack *hmetall = new THStack("hmetall","hmetall");
  hmetall->Add(hmet[5]);
  hmetall->Add(hmet[7]);
  hmetall->Add(hmet[4]);
  hmetall->Add(hmet[6]);
  hmetall->Add(hmet[3]);
  hmetall->Add(hmet[8]);

  TCanvas *cmet = new TCanvas("cmet","cmet",600,600);
  cmet->SetLogy();
  hmetall->Draw("histo");
  hmetall->SetMinimum(0.01);
  hmetall->SetMaximum(10000);
  hmetall->GetHistogram()->GetXaxis()->SetRangeUser(0.,500);
  hmetall->GetHistogram()->GetXaxis()->SetTitle("missing transverse energy (GeV)");
  hmetall->GetHistogram()->GetYaxis()->SetTitle("events/1 GeV");
  hmet[0] ->SetFillStyle(0);
  hmet[0]->Draw("histo same");
  hmetData->Draw("esame");
  leg->Draw("histo same");

  THStack *hmtall = new THStack("hmtall","hmtall");
  hmtall->Add(hmt[5]);
  hmtall->Add(hmt[7]);
  hmtall->Add(hmt[4]);
  hmtall->Add(hmt[6]);
  hmtall->Add(hmt[3]);
  hmtall->Add(hmt[8]);
  
  TCanvas *cmt = new TCanvas("cmt","cmt",600,600);
  cmt->SetLogy();
  hmtall->Draw("histo");
  hmtall->SetMinimum(0.01);
  hmtall->SetMaximum(10000);
  hmtall->GetHistogram()->GetXaxis()->SetRangeUser(0.,1000);
  hmtall->GetHistogram()->GetXaxis()->SetTitle("transverse mass (GeV)");
  hmtall->GetHistogram()->GetYaxis()->SetTitle("events/1 GeV");
  hmt[0] ->SetFillStyle(0);
  hmt[0]->Draw("histo same");
  hmtData->Draw("esame");
  leg->Draw("same");


  /*
  THStack *hratioall = new THStack("hratioall","hratioall");
  hratioall->Add(hetOmet[5]);
  hratioall->Add(hetOmet[7]);
  hratioall->Add(hetOmet[4]);
  hratioall->Add(hetOmet[6]);
  hratioall->Add(hetOmet[3]);
  hratioall->Add(hetOmet[8]);
  
  TCanvas *cmt = new TCanvas("cmt","cmt",600,600);
  cmt->SetLogy();
  hratioall->Draw("histo");
  hratioall->SetMinimum(0.01);
  hratioall->SetMaximum(10000);
  hratioall->GetHistogram()->GetXaxis()->SetRangeUser(0.,1000);
  hratioall->GetHistogram()->GetXaxis()->SetTitle("Et/MET");
  hratioall->GetHistogram()->GetYaxis()->SetTitle("events");
  hetOmet[0] ->SetFillStyle(0);
  hetOmet[0]->Draw("histo same");
  //hetOmetData->Draw("esame");
  leg->Draw("same");
  */
}
