{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


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


  ntu1->Draw("TMath::Log(sumPt2) >> logSumPt2_1","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10 ","goff");
  ntu2->Draw("TMath::Log(sumPt2) >> logSumPt2_2","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10 ","goff");
  
  ntu1->Draw("nTracks >> tracksNum1","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10","goff");
  ntu2->Draw("nTracks >> tracksNum2","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10","goff");

  ntu1->Draw("deltaPhi_HSumPt >> deltaPhi1","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10 ","goff");
  ntu2->Draw("deltaPhi_HSumPt >> deltaPhi2","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10","goff");

  ntu1->Draw("modSumPt / sum2PhoPt >> ptRatio1","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10","goff");
  ntu2->Draw("modSumPt / sum2PhoPt >> ptRatio2","isSig == 1 && sum2PhoPt > 0 && sum2PhoPt < 10","goff");

  TLegend leg (0.5, 0.7,0.8,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(&tracksNum1,"MC Hgg","F");
  leg->AddEntry(&tracksNum2,"MC Zee","F");

  TCanvas c1("c1","c1",500,500);
  tracksNum2.DrawNormalized();
  tracksNum1.DrawNormalized("sames");
  leg ->Draw("same");

  TCanvas c2("c2","c2",500,500);
  deltaPhi2.DrawNormalized();
  deltaPhi1.DrawNormalized("sames");
  leg ->Draw("same");

  TCanvas c3("c3","c3",500,500);
  logSumPt2_2.DrawNormalized();
  logSumPt2_1.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c4("c4","c4",500,500);
  ptRatio2.DrawNormalized();
  ptRatio1.DrawNormalized("sames");
  leg->Draw("same");  
}
