{
  TTree* ntu1 = (TTree*)_file0->Get("TMVA_vertexTree");
  TTree* ntu2 = (TTree*)_file1->Get("TMVA_vertexTree");
  
  TH1F tracksNum1("tracksNum1","tracksNum1",200,0.,200.);
  tracksNum1.SetLineColor(2);
  TH1F tracksNum2("tracksNum2","tracksNum2",200,0.,200.);
  tracksNum2.SetLineColor(4);

  TH1F deltaPhi1("deltaPhi1","deltaPhi1",150,0.,3.15);
  deltaPhi1.SetLineColor(2);
  TH1F deltaPhi2("deltaPhi2","deltaPhi2",150,0.,3.15);
  deltaPhi2.SetLineColor(4);

  ntu1->Draw("tracksNumPt05 >> tracksNum1","isSig == 1 && sum2PhoPt > 10 && sum2PhoPt < 15 ","goff");
  ntu2->Draw("tracksNumPt05 >> tracksNum2","isSig == 1 && sum2PhoPt > 10 && sum2PhoPt < 15","goff");

  ntu1->Draw("deltaPhi_HSumPt >> deltaPhi1","isSig == 1 && sum2PhoPt > 10 && sum2PhoPt < 15 ","goff");
  ntu2->Draw("deltaPhi_HSumPt >> deltaPhi2","isSig == 1 && sum2PhoPt > 10 && sum2PhoPt < 15","goff");

  TCanvas c1;
  tracksNum2.DrawNormalized();
  tracksNum1.DrawNormalized("sames");

  TCanvas c2;
  deltaPhi2.DrawNormalized();
  deltaPhi1.DrawNormalized("sames");


}
