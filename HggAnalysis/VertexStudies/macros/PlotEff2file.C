{
  gROOT->Reset();
  //TFile f1("TMVA_check.root");
  TFile f1("../output/TMVA_check_Zee.root");
 
  TH1F* NvtAll_1 = (TH1F*) f1.Get("NvtAll");
  TH1F* NvtPt2_1 = (TH1F*) f1.Get("NvtGood");
  TH1F* NvtBDT_1 = (TH1F*) f1.Get("NvtGood_BDT");

  TH1F* PtAll_1 = (TH1F*) f1.Get("PtAll");
  TH1F* PtPt2_1 = (TH1F*) f1.Get("PtGood");
  TH1F* PtBDT_1 = (TH1F*) f1.Get("PtGood_BDT");

  TFile f2("TMVA_check_goodVtx_Zmumu.root");

  TH1F* NvtAll_2 = (TH1F*) f2.Get("NvtAll");
  TH1F* NvtPt2_2 = (TH1F*) f2.Get("NvtGood");
  TH1F* NvtBDT_2 = (TH1F*) f2.Get("NvtGood_BDT");
  
  TH1F* PtAll_2 = (TH1F*) f2.Get("PtAll");
  TH1F* PtPt2_2 = (TH1F*) f2.Get("PtGood");
  TH1F* PtBDT_2 = (TH1F*) f2.Get("PtGood_BDT");

  

  int nMBMax =20;
  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.5,1.1);
  cc.SetStats(0); cc.GetXaxis()->SetTitle("number of PV"); cc.GetYaxis()->SetTitle("efficiency");
  cc.Draw();
  
   TH1F effVsNvPt_1 ("effPt_1","eff sumpt2", NvtAll_1->GetNbinsX(), NvtAll_1->GetXaxis()->GetXmin(), NvtAll_1->GetXaxis()->GetXmax()) ;
   TH1F effVsNvBDT_1 ("effBDT_1","eff BDT", NvtAll_1->GetNbinsX(), NvtAll_1->GetXaxis()->GetXmin(), NvtAll_1->GetXaxis()->GetXmax()) ;

   effVsNvPt_1.Sumw2();effVsNvBDT_1.Sumw2();
   NvtPt2_1->Sumw2();  NvtBDT_1->Sumw2(); NvtAll_1->Sumw2();

   effVsNvBDT_1.Divide(NvtBDT_1,NvtAll_1, 1,1,"B"); effVsNvBDT_1.SetLineColor(kGreen+1); effVsNvBDT_1.SetMarkerColor(kGreen+1); 
   effVsNvPt_1.Divide(NvtPt2_1,NvtAll_1, 1,1,"B"); effVsNvPt_1.SetLineColor(2); effVsNvPt_1.SetMarkerColor(2); 


   TH1F effVsNvPt_2 ("effPt_2","eff sumpt2", NvtAll_2->GetNbinsX(), NvtAll_2->GetXaxis()->GetXmin(), NvtAll_2->GetXaxis()->GetXmax()) ;
   TH1F effVsNvBDT_2 ("effBDT_2","eff BDT", NvtAll_2->GetNbinsX(), NvtAll_2->GetXaxis()->GetXmin(), NvtAll_2->GetXaxis()->GetXmax()) ;

   effVsNvPt_2.Sumw2();effVsNvBDT_2.Sumw2();
   NvtPt2_2->Sumw2();  NvtBDT_2->Sumw2(); NvtAll_2->Sumw2();

   effVsNvBDT_2.Divide(NvtBDT_2,NvtAll_2, 1,1,"B"); effVsNvBDT_2.SetLineColor(6); effVsNvBDT_2.SetMarkerColor(6); 
   effVsNvPt_2.Divide(NvtPt2_2,NvtAll_2, 1,1,"B"); effVsNvPt_2.SetLineColor(4); effVsNvPt_2.SetMarkerColor(4); 


   effVsNvPt_1.Draw("esame");
   effVsNvBDT_1.Draw("esame");

   effVsNvPt_2.Draw("esame");
   effVsNvBDT_2.Draw("esame");


   TCanvas c2;
   TH2F dd("dd","",200,0,200,1000,0.5,1.1);
   dd.SetStats(0); dd.GetXaxis()->SetTitle("Pt of the boson "); dd.GetYaxis()->SetTitle("efficiency");
   dd.Draw();


   TH1F effVsPt_Pt_1 ("effPt_1","eff sumpt2", PtAll_1->GetNbinsX(), PtAll_1->GetXaxis()->GetXmin(), PtAll_1->GetXaxis()->GetXmax()) ;
   TH1F effVsPt_BDT_1 ("effBDT_1","eff BDT", PtAll_1->GetNbinsX(), PtAll_1->GetXaxis()->GetXmin(), PtAll_1->GetXaxis()->GetXmax()) ;

   effVsPt_Pt_1.Sumw2();effVsPt_BDT_1.Sumw2();
   PtPt2_1->Sumw2();  PtBDT_1->Sumw2(); PtAll_1->Sumw2();

   effVsPt_BDT_1.Divide(PtBDT_1,PtAll_1, 1,1,"B"); effVsPt_BDT_1.SetLineColor(kGreen+1); effVsPt_BDT_1.SetMarkerColor(kGreen+1); 
   effVsPt_Pt_1.Divide(PtPt2_1,PtAll_1, 1,1,"B"); effVsPt_Pt_1.SetLineColor(2); effVsPt_Pt_1.SetMarkerColor(2); 


   TH1F effVsPt_Pt_2 ("effPt_2","eff sumpt2", PtAll_2->GetNbinsX(), PtAll_2->GetXaxis()->GetXmin(), PtAll_2->GetXaxis()->GetXmax()) ;
   TH1F effVsPt_BDT_2 ("effBDT_2","eff BDT", PtAll_2->GetNbinsX(), PtAll_2->GetXaxis()->GetXmin(), PtAll_2->GetXaxis()->GetXmax()) ;

   effVsPt_Pt_2.Sumw2();effVsPt_BDT_2.Sumw2();
   PtPt2_2->Sumw2();  PtBDT_2->Sumw2(); PtAll_2->Sumw2();

   effVsPt_BDT_2.Divide(PtBDT_2,PtAll_2, 1,1,"B"); effVsPt_BDT_2.SetLineColor(6); effVsPt_BDT_2.SetMarkerColor(6); 
   effVsPt_Pt_2.Divide(PtPt2_2,PtAll_2, 1,1,"B"); effVsPt_Pt_2.SetLineColor(4); effVsPt_Pt_2.SetMarkerColor(4); 


   effVsPt_Pt_1.Draw("esame");
   effVsPt_BDT_1.Draw("esame");

   effVsPt_Pt_2.Draw("esame");
   effVsPt_BDT_2.Draw("esame");
}
