{

  std::string cutString = "isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 2000 && photons_pt[0]>0 && photons_pt[1]>0";

  TTree* ntu1 = (TTree*)_file0->Get("TMVA_vertexTree");
  TTree* ntu2 = (TTree*)_file1->Get("TMVA_vertexTree");
  
  TH1F logSumPt2_1("logSumPt2_1","log(sumPt2)",200,-10.,30.);
  logSumPt2_1.GetXaxis()->SetTitle("log(SumPt2)");
  logSumPt2_1.GetYaxis()->SetTitle("a.u.");
  logSumPt2_1.SetLineColor(2);
  logSumPt2_1.SetFillColor(2);  
  logSumPt2_1.SetFillStyle(3005);

  TH1F logSumPt2_2("logSumPt2_2","log(sumPt2)",200,-10.,30.);
  logSumPt2_2.GetXaxis()->SetTitle("log(SumPt2)");
  logSumPt2_2.GetYaxis()->SetTitle("a.u.");
  logSumPt2_2.SetLineColor(kGreen+2);
  logSumPt2_2.SetFillColor(kGreen+2);  
  logSumPt2_2.SetFillStyle(3004);

  TH1F tracksNum1("tracksNum1","Number of tracks",100,0.,200.);
  tracksNum1.GetXaxis()->SetTitle("number of tracks");
  tracksNum1.GetYaxis()->SetTitle("a.u.");
  tracksNum1.SetLineColor(2);
  tracksNum1.SetFillColor(2);
  tracksNum1.SetFillStyle(3005);

  TH1F tracksNum2("tracksNum2","Number of tracks",100,0.,200.);
  tracksNum2.GetXaxis()->SetTitle("number of tracks");
  tracksNum2.GetYaxis()->SetTitle("a.u.");
  tracksNum2.SetLineColor(kGreen+2);
  tracksNum2.SetFillColor(kGreen+2);
  tracksNum2.SetFillStyle(3004);

  TH1F deltaPhi1("deltaPhi1","DeltaPhi",100,0.,3.15);
  deltaPhi1.GetXaxis()->SetTitle("#Delta#phi");
  deltaPhi1.GetYaxis()->SetTitle("a.u.");
  deltaPhi1.SetLineColor(2);
  deltaPhi1.SetFillColor(2);
  deltaPhi1.SetFillStyle(3005);

  TH1F deltaPhi2("deltaPhi2","DeltaPhi",100,0.,3.15);
  deltaPhi2.GetXaxis()->SetTitle("#Delta#phi");
  deltaPhi2.GetYaxis()->SetTitle("a.u.");
  deltaPhi2.SetLineColor(kGreen+2);
  deltaPhi2.SetFillColor(kGreen+2);
  deltaPhi2.SetFillStyle(3004);

  TH1F ptRatio1("ptRatio1","|SumPt| / bosonPt",100,0.,5.);
  ptRatio1.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptRatio1.GetYaxis()->SetTitle("a.u.");
  ptRatio1.SetLineColor(2);  
  ptRatio1.SetFillColor(2);  
  ptRatio1.SetFillStyle(3005);

  TH1F ptRatio2("ptRatio2","|SumPt| / bosonPt",100,0.,5.);
  ptRatio2.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptRatio2.GetYaxis()->SetTitle("a.u.");
  ptRatio2.SetLineColor(kGreen+2);  
  ptRatio2.SetFillColor(kGreen+2);  
  ptRatio2.SetFillStyle(3004);


  TH1F ptbal1("ptbal1","ptbal1",400,-200.,200.);
  ptbal1.GetXaxis()->SetTitle("ptbal");
  ptbal1.GetYaxis()->SetTitle("a.u.");
  ptbal1.SetLineColor(2);  
  ptbal1.SetFillColor(2);  
  ptbal1.SetFillStyle(3005);
  
  TH1F ptbal2("ptbal2","ptbal",400,-200.,200.);
  ptbal2.GetXaxis()->SetTitle("ptbal");
  ptbal2.GetYaxis()->SetTitle("a.u.");
  ptbal2.SetLineColor(kGreen+2);
  ptbal2.SetFillColor(kGreen+2);
  ptbal2.SetFillStyle(3004);

  TH1F ptasymm1("ptasymm1","ptasymm",100,-1.,1.);
  ptasymm1.GetXaxis()->SetTitle("ptasymm");
  ptasymm1.GetYaxis()->SetTitle("a.u.");
  ptasymm1.SetLineColor(2);  
  ptasymm1.SetFillColor(2);  
  ptasymm1.SetFillStyle(3005);
  
  TH1F ptasymm2("ptasymm2","ptasymm",100,-1.,1.);
  ptasymm2.GetXaxis()->SetTitle("ptasymm");
  ptasymm2.GetYaxis()->SetTitle("a.u.");
  ptasymm2.SetLineColor(kGreen+2);
  ptasymm2.SetFillColor(kGreen+2);
  ptasymm2.SetFillStyle(3004);


  TH1F SumPt2_1("SumPt2_1","sumPt2",200,0.,2000.);
  SumPt2_1.GetXaxis()->SetTitle("SumPt2");
  SumPt2_1.GetYaxis()->SetTitle("a.u.");
  SumPt2_1.SetLineColor(2);
  SumPt2_1.SetFillColor(2);  
  SumPt2_1.SetFillStyle(3005);

  TH1F SumPt2_2("SumPt2_2","sumPt2",200,0.,2000.);
  SumPt2_2.GetXaxis()->SetTitle("SumPt2");
  SumPt2_2.GetYaxis()->SetTitle("a.u.");
  SumPt2_2.SetLineColor(kGreen+2);
  SumPt2_2.SetFillColor(kGreen+2);  
  SumPt2_2.SetFillStyle(3004);

  TH2F SumPt2_vs_BosonPt_1("SumPt2_vs_BosonPt_1","sumPt2 vs boson Pt",500,0,500,200,0.,2000.);
  SumPt2_vs_BosonPt_1.GetXaxis()->SetTitle("boson pT");
  SumPt2_vs_BosonPt_1.GetYaxis()->SetTitle("SumPt2");

  TH2F SumPt2_vs_BosonPt_2("SumPt2_vs_BosonPt_2","sumPt2 vs boson Pt",500,0,500,200,0.,2000.);
  SumPt2_vs_BosonPt_2.GetXaxis()->SetTitle("boson pT");
  SumPt2_vs_BosonPt_2.GetYaxis()->SetTitle("SumPt2");


  ntu1->Draw("TMath::Log(sumPt2) >> logSumPt2_1",cutString.c_str(),"goff");
  ntu2->Draw("TMath::Log(sumPt2) >> logSumPt2_2",cutString.c_str(),"goff");

  ntu1->Draw("nTracks >> tracksNum1",cutString.c_str(),"goff");
  ntu2->Draw("nTracks >> tracksNum2",cutString.c_str(),"goff");

  ntu1->Draw("deltaPhi_HSumPt >> deltaPhi1",cutString.c_str(),"goff");
  ntu2->Draw("deltaPhi_HSumPt >> deltaPhi2",cutString.c_str(),"goff");

  ntu1->Draw("modSumPt / sum2PhoPt >> ptRatio1",cutString.c_str(),"goff");
  ntu2->Draw("modSumPt / sum2PhoPt >> ptRatio2",cutString.c_str(),"goff");

  ntu1->Draw("ptbal >> ptbal1",cutString.c_str(),"goff");
  ntu2->Draw("ptbal >> ptbal2",cutString.c_str(),"goff");

  ntu1->Draw("ptasym >> ptasymm1",cutString.c_str(),"goff");
  ntu2->Draw("ptasym >> ptasymm2",cutString.c_str(),"goff");
  
  ntu1->Draw("sumPt2 >> SumPt2_1",cutString.c_str(),"goff");
  ntu2->Draw("sumPt2 >> SumPt2_2",cutString.c_str(),"goff");

  ntu1->Draw("sumPt2:sum2PhoPt >> SumPt2_vs_BosonPt_1",cutString.c_str(),"goff");
  ntu2->Draw("sumPt2:sum2PhoPt >> SumPt2_vs_BosonPt_2",cutString.c_str(),"goff");



  TLegend leg (0.5, 0.7,0.8,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(&tracksNum1,"MC Hgg","F");
  leg->AddEntry(&tracksNum2,"MC Zee","F");

  TCanvas c1("c1","c1",500,500);
  tracksNum2.DrawNormalized();
  tracksNum1.DrawNormalized("sames");
  leg ->Draw("same");

  TCanvas c2("c2","c2",500,500);
  deltaPhi1.DrawNormalized();
  deltaPhi2.DrawNormalized("sames");
  leg ->Draw("same");

  TCanvas c3("c3","c3",500,500);
  logSumPt2_2.DrawNormalized();
  logSumPt2_1.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c4("c4","c4",500,500);
  ptRatio2.DrawNormalized();
  ptRatio1.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c5("c5","c5",500,500);
  ptbal2.DrawNormalized();
  ptbal1.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c6("c6","c6",500,500);
  ptasymm1.DrawNormalized();
  ptasymm2.DrawNormalized("sames");
  leg->Draw("same");  


 //  SumPt2_1.Scale(1./SumPt2_1.GetEntries());
//   SumPt2_2.Scale(1./SumPt2_2.GetEntries());

//   new TCanvas();
//   TH1F *hRatio = (TH1F*)SumPt2_1->Clone("hRatio");
//   hRatio->Divide(hRatio,&SumPt2_2);
//   hRatio->Draw();
 
//   TFile *fout = new TFile("correction.root","create");
 
//   new TCanvas();
//   SumPt2_vs_BosonPt_1.Draw("colz");

//   new TCanvas();
//   SumPt2_vs_BosonPt_2.Draw("colz");

 
//   SumPt2_vs_BosonPt_1.Write("SumPt2_vs_BosonPt_Hgg");
//   SumPt2_vs_BosonPt_2.Write("SumPt2_vs_BosonPt_Zee");
  
}
