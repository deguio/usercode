{
  TTree* ntu = (TTree*)_file0->Get("TMVA_vertexTree");
  
  TH1F logSumPt2_sig("logSumPt2_sig","log(sumPt2)",400,-10.,30.);
  logSumPt2_sig.SetLineColor(2);
  TH1F logSumPt2_bkg("logSumPt2_bkg","log(sumPt2)",400,-10.,30.);
  logSumPt2_bkg.SetLineColor(4);

  TH1F tracksNum_sig("tracksNum_sig","tracksNum_sig",200,0.,200.);
  tracksNum_sig.SetLineColor(2);
  TH1F tracksNum_bkg("tracksNum_bkg","tracksNum_bkg",200,0.,200.);
  tracksNum_bkg.SetLineColor(4);

  TH1F deltaPhi_sig("deltaPhi_sig","deltaPhi_sig",100,0.,3.15);
  deltaPhi_sig.SetLineColor(2);
  TH1F deltaPhi_bkg("deltaPhi_bkg","deltaPhi_bkg",100,0.,3.15);
  deltaPhi_bkg.SetLineColor(4);

  TH1F ptRatio_sig("ptRatio_sig","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_sig.SetLineColor(2);
  TH1F ptRatio_bkg("ptRatio_bkg","|SumPt| / bosonPt",100,0.,5.);
  ptRatio_bkg.SetLineColor(4);

  TH1F ptRatio_sig("ptRatio_sig","|SumPt| / bosonPt",100,0.,10.);
  ptRatio_sig.SetLineColor(2);




  ntu->Draw("log(sumPt2) >> logSumPt2_sig","isSig == 1 ","goff");
  ntu->Draw("log(sumPt2) >> logSumPt2_bkg","isSig == 0 ","goff");

  ntu->Draw("tracksNum >> tracksNum_sig","isSig == 1 ","goff");
  ntu->Draw("tracksNum >> tracksNum_bkg","isSig == 0 ","goff");

  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_sig","isSig == 1 ","goff");
  ntu->Draw("deltaPhi_HSumPt >> deltaPhi_bkg","isSig == 0 ","goff");

  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_sig","isSig == 1 ","goff");
  ntu->Draw("modSumPt / sum2PhoPt >> ptRatio_bkg","isSig == 0 ","goff");

  TCanvas c1;
  tracksNum_bkg.DrawNormalized();
  tracksNum_sig.DrawNormalized("sames");
  
  TCanvas c2;
  deltaPhi_sig.DrawNormalized();
  deltaPhi_bkg.DrawNormalized("sames");

  TCanvas c3;
  logSumPt2_sig.DrawNormalized();
  logSumPt2_bkg.DrawNormalized("sames");

  TCanvas c4;
  ptRatio_bkg.DrawNormalized();
  ptRatio_sig.DrawNormalized("sames");


}
