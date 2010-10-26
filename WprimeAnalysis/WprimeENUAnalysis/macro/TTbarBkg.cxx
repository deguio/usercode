

{
  gROOT->LoadMacro("setTDRStyle.C");
  setTDRStyle();
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  float lumi = 15.0416; //all
      
  int N = 10;
  
  int nRe  = 10;
  int nReB = 1;
  

  // DATA histos
  TFile *fData =  TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_Data_7TeV.root");
  TH1F * hmtData    = (TH1F*)fData->Get("hmt");
  TH1F * hmt1bData  = (TH1F*)fData->Get("hmt1b");
  TH1F * hmt2bData  = (TH1F*)fData->Get("hmt2b");
  TH1F * hBdiscData = (TH1F*)fData->Get("hBdisc");
  TH1F * hNjetsData = (TH1F*)fData->Get("hNjets");
  hmtData    -> SetMarkerStyle(20);
  hmt1bData  -> SetMarkerStyle(20);
  hmt2bData  -> SetMarkerStyle(20);
  hBdiscData -> SetMarkerStyle(20);
  hNjetsData -> SetMarkerStyle(20);

  hmtData   -> Rebin(nRe);
  hmt1bData -> Rebin(nRe);
  hmt2bData -> Rebin(nRe);
  hBdiscData ->Rebin(nReB);

  // MC histos
  TFile *f[10];
  f[0] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_WprimeM800_7TeV.root");
  f[1] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_WprimeM1000_7TeV.root");
  f[2] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_WprimeM1500_7TeV.root");
  f[3] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_QCD_7TeV.root");
  f[4] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_DY_7TeV.root");
  f[5] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_Dibosons_7TeV.root");
  f[6] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_Wtaunu_7TeV.root");
  f[7] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_TTbar_7TeV.root");
  f[8] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_Wenu_7TeV.root");
  f[9] = TFile::Open("../bin/HistosTTbar_Btag3.0_Eta2.4_Pt20/outputHistos_GammaPlusJets_7TeV.root");

  int mycolor;

  TH1F *hmt[10];
  TH1F *hmt1b[10];
  TH1F *hmt2b[10];
  TH1F *hBdisc[10];
  TH1F *hNjets[10];

  for (int i = 0; i< N; i++){
    hmt[i]    = (TH1F*)f[i]->Get("hmt");
    hmt1b[i]  = (TH1F*)f[i]->Get("hmt1b");
    hmt2b[i]  = (TH1F*)f[i]->Get("hmt2b");
    hBdisc[i] = (TH1F*)f[i]->Get("hBdisc");
    hNjets[i] = (TH1F*)f[i]->Get("hNjets");

    hmt[i]->Scale(lumi);
    hmt1b[i]->Scale(lumi);
    hmt2b[i]->Scale(lumi);
    hBdisc[i]->Scale(lumi);
    hNjets[i]->Scale(lumi);

    hmt[i]->Rebin(nRe);
    hmt1b[i]->Rebin(nRe);
    hmt2b[i]->Rebin(nRe);
    hBdisc[i]->Rebin(nReB);

    if (i==0) mycolor = kRed; //W'
    if (i==1) mycolor = kRed; //W'
    if (i==2) mycolor = kRed; //W'
    if (i==3) mycolor = kGray+2; // QCD
    if (i==4) mycolor = kOrange+1; // Zee
    if (i==5) mycolor = kGreen-3; // dibosons 
    if (i==6) mycolor = kYellow; // Wtaunu
    if (i==7) mycolor = kMagenta; // ttbar
    if (i==8) mycolor = kAzure+1; // Wenu
    if (i==9) mycolor = kPink-2; // Wenu
    
    hmt[i]   -> SetFillColor(mycolor);
    hmt1b[i] -> SetFillColor(mycolor);
    hmt2b[i] -> SetFillColor(mycolor);
    hBdisc[i]-> SetFillColor(mycolor);
    hNjets[i]-> SetFillColor(mycolor);

  }

 
  TLegend *leg = new TLegend (0.6,0.6,0.89,0.89);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(hmt[0],"W', M =  800 GeV","F");
  //leg->AddEntry(hmt[1],"W', M = 1000 GeV","F");
  //leg->AddEntry(hmt[2],"W', M = 1500 GeV","F");
  leg->AddEntry(hmt[8],"W #rightarrow e #nu","F");
  leg->AddEntry(hmt[4],"DY #rightarrow ll","F");
  leg->AddEntry(hmt[5],"Dibosons","F");
  leg->AddEntry(hmt[6],"W #rightarrow #tau #nu","F");
  leg->AddEntry(hmt[7],"t#bar{t}","F");
  leg->AddEntry(hmt[3],"QCD ","F");
  leg->AddEntry(hmt[9],"#gamma+jets","F");
  leg->AddEntry(hmtData,"data","PL");

  THStack *hNjetsAll = new THStack("hNjetsAll","hNjetsAll");
  hNjetsAll->Add(hNjets[3]);
  hNjetsAll->Add(hNjets[9]);
  hNjetsAll->Add(hNjets[5]);
  hNjetsAll->Add(hNjets[4]);
  hNjetsAll->Add(hNjets[6]);
  hNjetsAll->Add(hNjets[8]);
  hNjetsAll->Add(hNjets[7]);

  THStack *hBdiscAll = new THStack("hBdiscAll","hBdiscAll");
  hBdiscAll->Add(hBdisc[3]);
  hBdiscAll->Add(hBdisc[9]);
  hBdiscAll->Add(hBdisc[5]);
  hBdiscAll->Add(hBdisc[4]);
  hBdiscAll->Add(hBdisc[6]);
  hBdiscAll->Add(hBdisc[8]);
  hBdiscAll->Add(hBdisc[7]);

  THStack *hmt1bAll = new THStack("hmt1bAll","hmt1bAll");
  hmt1bAll->Add(hmt1b[3]);
  hmt1bAll->Add(hmt1b[9]);
  hmt1bAll->Add(hmt1b[5]);
  hmt1bAll->Add(hmt1b[4]);
  hmt1bAll->Add(hmt1b[6]);
  hmt1bAll->Add(hmt1b[8]);
  hmt1bAll->Add(hmt1b[7]);

  THStack *hmt2bAll = new THStack("hmt2bAll","hmt2bAll");
  hmt2bAll->Add(hmt2b[3]);
  hmt2bAll->Add(hmt2b[9]);
  hmt2bAll->Add(hmt2b[5]);
  hmt2bAll->Add(hmt2b[4]);
  hmt2bAll->Add(hmt2b[6]);
  hmt2bAll->Add(hmt2b[8]);
  hmt2bAll->Add(hmt2b[7]);

  TCanvas * cnjets = new TCanvas("cnjets","cnjets",500,500);
  cnjets->SetLogy();
  hNjetsAll->Draw("a histo");
  hNjetsAll->SetMinimum(0.01);
  hNjetsAll->SetMaximum(1000000);
  hNjetsAll->GetHistogram()->GetXaxis()->SetRangeUser(0.,6);
  hNjetsAll->GetHistogram()->GetXaxis()->SetTitle("number of jets");
  hNjetsAll->GetHistogram()->GetYaxis()->SetTitle("events");
  hNjetsAll->GetHistogram()->GetXaxis()->SetRangeUser(0,6);
  hNjetsAll->Draw("histo");
  hNjetsData->Draw("esame");
  leg->Draw("same");
  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.5,500000,"CMS 2010 Preliminary");
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,100000,"#sqrt{s} = 7 TeV");
  latex.DrawLatex(0.5,50000,"#intLdt = 15.041 pb^{-1}");

  TCanvas * cbdisc = new TCanvas("cbdisc","cbdisc",500,500);
  cbdisc->SetLogy();
  hBdiscAll->Draw("a histo");
  hBdiscAll->SetMinimum(0.01);
  hBdiscAll->SetMaximum(200);
  hBdiscAll->GetHistogram()->GetXaxis()->SetRangeUser(0.,6);
  hBdiscAll->GetHistogram()->GetXaxis()->SetTitle("B discriminant");
  hBdiscAll->GetHistogram()->GetYaxis()->SetTitle("events");
  hBdiscAll->GetHistogram()->GetXaxis()->SetRangeUser(0,6);
  hBdiscAll->Draw("histo");
  hBdiscData->Draw("esame");
  leg->Draw("same");
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.5,150,"CMS 2010 Preliminary");
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.5,80,"#sqrt{s} = 7 TeV");
  latex.DrawLatex(0.5,30,"#intLdt = 15.041 pb^{-1}");


  TCanvas * c1b = new TCanvas("c1b","c1b",500,500);
  c1b->SetLogy();
  hmt1bAll->Draw("a histo");
  hmt1bAll->SetMinimum(0.01);
  hmt1bAll->SetMaximum(100);
  hmt1bAll->GetHistogram()->GetXaxis()->SetRangeUser(0.,500);
  hmt1bAll->GetHistogram()->GetXaxis()->SetTitle("M^{T} (GeV)");
  hmt1bAll->GetHistogram()->GetYaxis()->SetTitle("events/10 GeV");
  hmt1bAll->GetHistogram()->GetXaxis()->SetRangeUser(0,500);
  hmt1bAll->Draw("histo");
  hmt1bData->Draw("esame");
  leg->Draw("same");
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(20,90,"CMS 2010 Preliminary");
  latex.SetTextSize(0.035);
  latex.DrawLatex(20,42,"#sqrt{s} = 7 TeV");
  latex.DrawLatex(300,0.5,"#intLdt = 15.041 pb^{-1}");

  TCanvas * c2b = new TCanvas("c2b","c2b",500,500);
  c2b->SetLogy();
  hmt2bAll->Draw("a histo");
  hmt2bAll->SetMinimum(0.01);
  hmt2bAll->SetMaximum(100.0);
  hmt2bAll->GetHistogram()->GetXaxis()->SetRangeUser(0.,500);
  hmt2bAll->GetHistogram()->GetXaxis()->SetTitle("M^{T} (GeV)");
  hmt2bAll->GetHistogram()->GetYaxis()->SetTitle("events/10 GeV");
  hmt2bAll->GetHistogram()->GetXaxis()->SetRangeUser(0,500);
  hmt2bAll->Draw("histo");
  hmt2bData->Draw("esame");
  leg->Draw("same");
  latex.SetTextFont(42);
  latex.SetTextSize(0.04);
  latex.DrawLatex(20,90,"CMS 2010 Preliminary");
  latex.SetTextSize(0.035);
  latex.DrawLatex(20,42,"#sqrt{s} = 7 TeV");
  latex.DrawLatex(300,0.5,"#intLdt = 15.041 pb^{-1}");

  int bin1  = hmt1bAll->GetHistogram()->FindBin(100.);
  int bin2  = hmt1bAll->GetHistogram()->GetNbinsX()+1;
  float n1b = hmt1b[8]->Integral(bin1,bin2)+hmt1b[7]->Integral(bin1,bin2)+
              hmt1b[6]->Integral(bin1,bin2)+hmt1b[5]->Integral(bin1,bin2)+ 
              hmt1b[4]->Integral(bin1,bin2)+hmt1b[3]->Integral(bin1,bin2)+
              hmt1b[9]->Integral(bin1,bin2);
  float n2b = hmt2b[8]->Integral(bin1,bin2)+hmt2b[7]->Integral(bin1,bin2)+
              hmt2b[6]->Integral(bin1,bin2)+hmt2b[5]->Integral(bin1,bin2)+ 
              hmt2b[4]->Integral(bin1,bin2)+hmt2b[3]->Integral(bin1,bin2)+
              hmt2b[9]->Integral(bin1,bin2);
  //    float n1b = hmt1b[7]->Integral(bin1,bin2);
  //   float n2b = hmt2b[7]->Integral(bin1,bin2);

  cout << "******* MC " << endl;
  cout << "Total number of events with 1 b-jet (MC) = " << n1b << std::endl;
  cout << "Total number of events with 2 b-jet (MC) = " << n2b << std::endl;
  cout << "---- n1b in each sample " << endl; 
  cout << "W' M= 800 : " << hmt1b[0]->Integral(bin1,bin2) << endl;
  cout << "W' M=1000 : " << hmt1b[1]->Integral(bin1,bin2) << endl;
  cout << "W' M=1500 : " << hmt1b[2]->Integral(bin1,bin2) << endl;
  cout << "Wenu      : " << hmt1b[8]->Integral(bin1,bin2) << endl;
  cout << "TTbar     : " << hmt1b[7]->Integral(bin1,bin2) << endl;
  cout << "Zee       : " << hmt1b[4]->Integral(bin1,bin2) << endl;
  cout << "Wtaunu    : " << hmt1b[6]->Integral(bin1,bin2) << endl;
  cout << "QCD       : " << hmt1b[3]->Integral(bin1,bin2) << endl;
  cout << "dibosons  : " << hmt1b[5]->Integral(bin1,bin2) << endl;
  cout << "---- n2b in each sample " << endl; 
  cout << "W' M= 800 : " << hmt2b[0]->Integral(bin1,bin2) << endl;
  cout << "W' M=1000 : " << hmt2b[1]->Integral(bin1,bin2) << endl;
  cout << "W' M=1500 : " << hmt2b[2]->Integral(bin1,bin2) << endl;
  cout << "Wenu      : " << hmt2b[8]->Integral(bin1,bin2) << endl;
  cout << "TTbar     : " << hmt2b[7]->Integral(bin1,bin2) << endl;
  cout << "Zee       : " << hmt2b[4]->Integral(bin1,bin2) << endl;
  cout << "Wtaunu    : " << hmt2b[6]->Integral(bin1,bin2) << endl;
  cout << "QCD       : " << hmt2b[3]->Integral(bin1,bin2) << endl;
  cout << "dibosons  : " << hmt2b[5]->Integral(bin1,bin2) << endl;
  

  float A1 = 0.118; 
  float A2 = 0.848; 
  float eA1 = 0.01;
  float eA2 = 0.01;

  float eb = (A1/A2+2)/(n1b/n2b+2);
  float eeb = (pow(eA1/A1,2) + pow(eA2/A2,2))/(pow(A1/A2+2,2)) + (1./n1b + 1./n2b)/(pow(n1b/n2b+2,2));
  eeb = eb * sqrt(eeb);
  cout<<"b-tagging efficiency = "<< eb <<" +/- "<<eeb<<endl;

  float nexp = n2b/(eb*eb*A2);
  float errn = nexp*sqrt(1./n2b + eA2/A2*eA2/A2 + pow(2*eeb/eb,2));
  cout<<"True # of ttbar events = "<<hmt[7]->Integral(bin1,bin2)<<endl;
  cout<<"Estimated N(ttbar) from b-tag method = "<< nexp <<" +/- "<< errn <<endl;
  cout<<endl;
 


  cout << "**** DATA " << endl;
  float n1bData = hmt1bData->Integral(bin1,bin2);
  float n2bData = hmt2bData->Integral(bin1,bin2);
  
  cout << "Number of events with 1 b-jet (DATA) = " << n1bData << std::endl;
  cout << "Number of events with 2 b-jet (DATA) = " << n2bData << std::endl;

  
  n2bData = hmt2bData->Integral(0,bin1)*hmt2b[7]->Integral(bin1,bin2)/hmt2b[7]->Integral(0,bin1);
  cout << "estimated number of 2b in data above 100 GeV (PORCATA!) = " << n2bData << endl; 

  
  if (n1bData!=0 && n2bData!=0 )
    {
      eb = (A1/A2+2)/(n1bData/n2bData+2);
      eeb = (pow(eA1/A1,2) + pow(eA2/A2,2))/(pow(A1/A2+2,2)) + (1./n1bData + 1./n2bData)/(pow(n1bData/n2bData+2,2));
      eeb = eb * sqrt(eeb);
      cout<<"b-tagging efficiency = "<< eb <<" +/- "<<eeb<<endl;
      
      float nexpData = n2bData/(eb*eb*A2);
      errn = nexpData*sqrt(1./n2bData + eA2/A2*eA2/A2 + pow(2*eeb/eb,2));
      cout<<"MC # of ttbar events = "<<hmt[7]->Integral(bin1,bin2)<<endl;
      cout<<"Estimated N(ttbar) from b-tag method (DATA) = "<< nexpData <<" +/- "<< errn <<endl;
      cout<<endl;
    }
}
