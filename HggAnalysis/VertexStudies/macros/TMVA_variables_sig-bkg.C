{
  gROOT->SetStyle("Plain");

  TTree* ntu = (TTree*)_file0->Get("TMVA_vertexTree");
  
  TH1F logSumPt2_sig("logSumPt2_sig","log(sumPt2)",400,-10.,30.);
  logSumPt2_sig.GetXaxis()->SetTitle("SumPt2");
  logSumPt2_sig.GetYaxis()->SetTitle("a.u.");
  logSumPt2_sig.SetLineColor(2);
  logSumPt2_sig.SetFillColor(2);  
  logSumPt2_sig.SetFillStyle(3002);
  
  TH1F logSumPt2_bkg("logSumPt2_bkg","log(sumPt2)",400,-10.,30.);
  logSumPt2_bkg.GetXaxis()->SetTitle("SumPt2");
  logSumPt2_bkg.GetYaxis()->SetTitle("a.u.");
  logSumPt2_bkg.SetLineColor(4);
  logSumPt2_bkg.SetFillColor(4);
  logSumPt2_bkg.SetFillStyle(3002);

  TH1F tracksNum_sig("tracksNum_sig","Number of tracks",200,0.,200.);
  tracksNum_sig.GetXaxis()->SetTitle("number of tracks");
  tracksNum_sig.GetYaxis()->SetTitle("a.u.");
  tracksNum_sig.SetLineColor(2);
  tracksNum_sig.SetFillColor(2);
  tracksNum_sig.SetFillStyle(3002);
  
  TH1F tracksNum_bkg("tracksNum_bkg","Number of tracks",200,0.,200.);
  tracksNum_bkg.GetXaxis()->SetTitle("number of tracks");
  tracksNum_bkg.GetYaxis()->SetTitle("a.u.");
  tracksNum_bkg.SetLineColor(4);
  tracksNum_bkg.SetFillColor(4);
  tracksNum_bkg.SetFillStyle(3002);

  TH1F deltaPhi_sig("deltaPhi_sig","deltaPhi",100,0.,3.15);
  deltaPhi_sig.GetXaxis()->SetTitle("#Delta#phi");
  deltaPhi_sig.GetYaxis()->SetTitle("a.u.");
  deltaPhi_sig.SetLineColor(2);
  deltaPhi_sig.SetFillColor(2);
  deltaPhi_sig.SetFillStyle(3002);

  TH1F deltaPhi_bkg("deltaPhi_bkg","deltaPhi",100,0.,3.15);
  deltaPhi_bkg.GetXaxis()->SetTitle("#Delta#phi");
  deltaPhi_bkg.GetYaxis()->SetTitle("a.u.");
  deltaPhi_bkg.SetLineColor(4);
  deltaPhi_bkg.SetFillColor(4);
  deltaPhi_bkg.SetFillStyle(3002);

  TH1F ptRatio_sig("ptRatio_sig","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_sig.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptRatio_sig.GetYaxis()->SetTitle("a.u.");
  ptRatio_sig.SetLineColor(2);  
  ptRatio_sig.SetFillColor(2);  
  ptRatio_sig.SetFillStyle(3002);
  
  TH1F ptRatio_bkg("ptRatio_bkg","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_bkg.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptRatio_bkg.GetYaxis()->SetTitle("a.u.");
  ptRatio_bkg.SetLineColor(4);
  ptRatio_bkg.SetFillColor(4);
  ptRatio_bkg.SetFillStyle(3002);

  TH1F ptbal_sig("ptbal_sig","ptbal",100,-100.,200.);
  ptbal_sig.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptbal_sig.GetYaxis()->SetTitle("a.u.");
  ptbal_sig.SetLineColor(2);  
  ptbal_sig.SetFillColor(2);  
  ptbal_sig.SetFillStyle(3002);
  
  TH1F ptbal_bkg("ptbal_bkg","ptbal",100,-100.,200.);
  ptbal_bkg.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptbal_bkg.GetYaxis()->SetTitle("a.u.");
  ptbal_bkg.SetLineColor(4);
  ptbal_bkg.SetFillColor(4);
  ptbal_bkg.SetFillStyle(3002);

  TH1F ptasymm_sig("ptasymm_sig","ptasymm",100,-1.,1.);
  ptasymm_sig.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptasymm_sig.GetYaxis()->SetTitle("a.u.");
  ptasymm_sig.SetLineColor(2);  
  ptasymm_sig.SetFillColor(2);  
  ptasymm_sig.SetFillStyle(3002);
  
  TH1F ptasymm_bkg("ptasymm_bkg","ptasymm",100,-1.,1.);
  ptasymm_bkg.GetXaxis()->SetTitle("|SumPt|/p^{T}_{H}");
  ptasymm_bkg.GetYaxis()->SetTitle("a.u.");
  ptasymm_bkg.SetLineColor(4);
  ptasymm_bkg.SetFillColor(4);
  ptasymm_bkg.SetFillStyle(3002);

  ntu->Draw("TMath::Log(sumPt2) >> logSumPt2_sig","isSig == 1 && sumPt2>0","goff");
  ntu->Draw("TMath::Log(sumPt2) >> logSumPt2_bkg","isSig == 0 && sumPt2>0","goff");

  ntu->Draw("nTracks >> tracksNum_sig","isSig == 1 ","goff");
  ntu->Draw("nTracks >> tracksNum_bkg","isSig == 0 ","goff");

  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_sig","isSig == 1 ","goff");
  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_bkg","isSig == 0 ","goff");

  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_sig","isSig == 1 ","goff");
  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_bkg","isSig == 0 ","goff");

  ntu->Draw("ptbal >> ptbal_sig","isSig == 1 ","goff");
  ntu->Draw("ptbal >> ptbal_bkg","isSig == 0 ","goff");

  ntu->Draw("ptasym >> ptasymm_sig","isSig == 1 ","goff");
  ntu->Draw("ptasym >> ptasymm_bkg","isSig == 0 ","goff");

  TLegend leg (0.5, 0.7,0.8,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(&tracksNum_bkg,"Background","F");
  leg->AddEntry(&tracksNum_sig,"Signal","F");


  TCanvas c1("c1","c1",500,500);
  tracksNum_bkg.DrawNormalized();
  tracksNum_sig.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c2("c2","c2",500,500);
  deltaPhi_sig.DrawNormalized();
  deltaPhi_bkg.DrawNormalized("sames");
  leg->Draw("same");  
 
  TCanvas c3("c3","c3",500,500);
  logSumPt2_sig.DrawNormalized();
  logSumPt2_bkg.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c4("c4","c4",500,500);
  ptRatio_bkg.DrawNormalized();
  ptRatio_sig.DrawNormalized("sames");
  leg->Draw("same");  


  TCanvas c5("c5","c5",500,500);
  ptbal_bkg.DrawNormalized();
  ptbal_sig.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c6("c6","c6",500,500);
  ptasymm_bkg.DrawNormalized();
  ptasymm_sig.DrawNormalized("sames");
  leg->Draw("same");  

}
