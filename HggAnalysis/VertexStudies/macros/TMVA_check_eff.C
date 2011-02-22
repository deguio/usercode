{
  TFile fmva("/Users/deguio/Documents/Universita/Lavoro/VertexStudies/Analysis/MyTest/output/TMVA_check.root");
  TH1F* mva_s = (TH1F*) fmva.Get("bdtH");
  TH1F* mva_b = (TH1F*) fmva.Get("bdtBkg");
  TH1F* sigPt2 = (TH1F*) fmva.Get("pt2h");
  TH1F* bkgPt2 = (TH1F*) fmva.Get("pt2bkg");

  TH1F* fullAll = (TH1F*) fmva.Get("NvtAll");
  TH1F* fullPt2 = (TH1F*) fmva.Get("NvtGood");
  TH1F* fullBDT = (TH1F*) fmva.Get("NvtGood_BDT");


  const int nMBMax = 20;
  TH1F goodPt("goodPtt","apt",nMBMax,1,nMBMax+1);
  TH1F goodMva("goodmvat","amva",nMBMax,1,nMBMax+1);
  TH1F all("all","al",nMBMax,1,nMBMax+1);


  for (int i=0; i< 20000;i++){
    float s =   sigPt2->GetRandom();
    if ( gRandom->Uniform() < sigPt2->GetBinContent(sigPt2->GetNbinsX()+1)/sigPt2->GetEntries() ) s = 1000;  //tengo conto degli overflow
    float b[nMBMax];
    float bmva[nMBMax];

    float bmax = -10;
    for (int j=0; j<nMBMax;j++){
      float bt = bkgPt2->GetRandom();
      if (bt > bmax){ bmax = bt;}
      b[j]=bmax;
    }

    float smva =   mva_s->GetRandom();
    bmax = -100000;
    for (int j=0; j<nMBMax;j++){
      float bt = mva_b->GetRandom();
      if (bt > bmax){ bmax = bt;}
      bmva[j]=bmax;
    }
    
    for (int j=0; j<nMBMax;j++){
      all.Fill(j+2);
      if(s > b[j]) goodPt.Fill(j+2);
      if(smva > bmva[j]) goodMva.Fill(j+2);
    }
      
  }

  TLegend legend(0.68, 0.78, 0.99, 0.99);
  legend.SetFillColor(kWhite);



  //all.Draw();
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.5,1.1);
  cc.Draw();
  
   TH1F effVsNvToyPt ("ToyPt","eff sumpt2 toy mc", goodPt.GetNbinsX(), goodPt.GetXaxis()->GetXmin(), goodPt.GetXaxis()->GetXmax()) ;
   TH1F effVsNvToyBDT ("ToyBDT","eff bdt toy mc", goodPt.GetNbinsX(), goodPt.GetXaxis()->GetXmin(), goodPt.GetXaxis()->GetXmax()) ;

   TH1F effVsNvPt ("effPt","eff sumpt2", fullAll->GetNbinsX(), fullAll->GetXaxis()->GetXmin(), fullAll->GetXaxis()->GetXmax()) ;
   TH1F effVsNvBDT ("effBDT","eff BDT", fullAll->GetNbinsX(), fullAll->GetXaxis()->GetXmin(), fullAll->GetXaxis()->GetXmax()) ;

   effVsNvToyPt.Sumw2();effVsNvToyBDT.Sumw2();effVsNvPt.Sumw2();effVsNvBDT.Sumw2();
   goodPt.Sumw2();goodMva.Sumw2(); fullPt2->Sumw2();  fullBDT->Sumw2();

   effVsNvToyPt.Divide(&goodPt,&all,1,1,"B"); 
   effVsNvToyBDT.Divide(&goodMva,&all,1,1,"B"); effVsNvToyBDT.SetLineColor(4); effVsNvToyBDT.SetMarkerColor(4);

   effVsNvBDT.Divide(fullBDT,fullAll,1,1,"B"); effVsNvBDT.SetLineColor(kGreen+1); effVsNvBDT.SetMarkerColor(kGreen+1); 
   effVsNvPt.Divide(fullPt2,fullAll,1,1,"B"); effVsNvPt.SetLineColor(2); effVsNvPt.SetMarkerColor(2); 

   effVsNvToyPt.Draw("esame");
   effVsNvToyBDT.Draw("esame");
   effVsNvPt.Draw("esame");
   effVsNvBDT.Draw("esame");


   legend.AddEntry(&effVsNvToyPt   , "toy MC sumPt2");
   legend.AddEntry(&effVsNvToyBDT  , "toy MC BDT");
   legend.AddEntry(&effVsNvPt      , "sumPt2");
   legend.AddEntry(&effVsNvBDT     , "BDT");

   legend.Draw("same");


}
