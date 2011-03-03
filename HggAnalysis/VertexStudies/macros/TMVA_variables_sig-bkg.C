{
  gROOT->SetStyle("Plain");

  TTree* ntu = (TTree*)_file0->Get("TMVA_vertexTree");
  
  TH1F logSumPt2_sig("logSumPt2_sig","log(sumPt2)",400,-10.,30.);
  logSumPt2_sig.SetLineColor(2);
  logSumPt2_sig.SetFillColor(2);  
  logSumPt2_sig.SetFillStyle(3002);
  
  TH1F logSumPt2_bkg("logSumPt2_bkg","log(sumPt2)",400,-10.,30.);
  logSumPt2_bkg.SetLineColor(4);
  logSumPt2_bkg.SetFillColor(4);
  logSumPt2_bkg.SetFillStyle(3002);

  TH1F tracksNum_sig("tracksNum_sig","tracksNum_sig",200,0.,200.);
  tracksNum_sig.SetLineColor(2);
  tracksNum_sig.SetFillColor(2);
  tracksNum_sig.SetFillStyle(3002);
  
  TH1F tracksNum_bkg("tracksNum_bkg","tracksNum_bkg",200,0.,200.);
  tracksNum_bkg.SetLineColor(4);
  tracksNum_bkg.SetFillColor(4);
  tracksNum_bkg.SetFillStyle(3002);

  TH1F deltaPhi_sig("deltaPhi_sig","deltaPhi_sig",100,0.,3.15);
  deltaPhi_sig.SetLineColor(2);
  deltaPhi_sig.SetFillColor(2);
  deltaPhi_sig.SetFillStyle(3002);

  TH1F deltaPhi_bkg("deltaPhi_bkg","deltaPhi_bkg",100,0.,3.15);
  deltaPhi_bkg.SetLineColor(4);
  deltaPhi_bkg.SetFillColor(4);
  deltaPhi_bkg.SetFillStyle(3002);

  TH1F ptRatio_sig("ptRatio_sig","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_sig.SetLineColor(2);  
  ptRatio_sig.SetFillColor(2);  
  ptRatio_sig.SetFillStyle(3002);
  
  TH1F ptRatio_bkg("ptRatio_bkg","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_bkg.SetLineColor(4);
  ptRatio_bkg.SetFillColor(4);
  ptRatio_bkg.SetFillStyle(3002);


  ntu->Draw("log(sumPt2) >> logSumPt2_sig","isSig == 1 ","goff");
  ntu->Draw("log(sumPt2) >> logSumPt2_bkg","isSig == 0 ","goff");

  ntu->Draw("nTracks >> tracksNum_sig","isSig == 1 ","goff");
  ntu->Draw("nTracks >> tracksNum_bkg","isSig == 0 ","goff");

  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_sig","isSig == 1 ","goff");
  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_bkg","isSig == 0 ","goff");

  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_sig","isSig == 1 ","goff");
  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_bkg","isSig == 0 ","goff");

  TLegend leg (0.6, 0.8,0.8,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(&tracksNum_bkg,"Background","F");
  leg->AddEntry(&tracksNum_sig,"Signal","F");


  TCanvas c1;
  tracksNum_bkg.DrawNormalized();
  tracksNum_sig.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c2;
  deltaPhi_sig.DrawNormalized();
  deltaPhi_bkg.DrawNormalized("sames");
  leg->Draw("same");  
 
  TCanvas c3;
  logSumPt2_sig.DrawNormalized();
  logSumPt2_bkg.DrawNormalized("sames");
  leg->Draw("same");  

  TCanvas c4;
  ptRatio_bkg.DrawNormalized();
  ptRatio_sig.DrawNormalized("sames");
  leg->Draw("same");  

}
