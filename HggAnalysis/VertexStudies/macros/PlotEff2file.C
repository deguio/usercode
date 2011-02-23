{
  //gROOT->Reset();
  //TFile f1("TMVA_check.root");
  //TFile f1("./output/TMVA_check_ZeeMC.root");
 
  TH1F* NvtAll_1 = (TH1F*) _file0->Get("NvtAll");
  TH1F* NvtPt2_1 = (TH1F*) _file0->Get("NvtGood");
  TH1F* NvtBDT_1 = (TH1F*) _file0->Get("NvtGood_BDT");

  TH1F* PtAll_1 = (TH1F*) _file0->Get("PtAll");
  TH1F* PtPt2_1 = (TH1F*) _file0->Get("PtGood");
  TH1F* PtBDT_1 = (TH1F*) _file0->Get("PtGood_BDT");

  TGraphAsymmErrors effVsNv_BDT_1;
  TGraphAsymmErrors effVsNv_Pt_1;

  NvtPt2_1->Sumw2();  NvtBDT_1->Sumw2(); NvtAll_1->Sumw2();
  
  effVsNv_BDT_1.BayesDivide(NvtBDT_1,NvtAll_1, "cp"); effVsNv_BDT_1.SetLineColor(kGreen+1); effVsNv_BDT_1.SetMarkerColor(kGreen+1); 
  effVsNv_Pt_1.BayesDivide(NvtPt2_1,NvtAll_1, "cp"); effVsNv_Pt_1.SetLineColor(2); effVsNv_Pt_1.SetMarkerColor(2); 
  
  
  TGraphAsymmErrors effVsPt_Pt_1;
  TGraphAsymmErrors effVsPt_BDT_1;

  PtPt2_1->Sumw2();  PtBDT_1->Sumw2(); PtAll_1->Sumw2();

  effVsPt_BDT_1.BayesDivide(PtBDT_1,PtAll_1,"cp"); effVsPt_BDT_1.SetLineColor(kGreen+1); effVsPt_BDT_1.SetMarkerColor(kGreen+1); 
  effVsPt_Pt_1.BayesDivide(PtPt2_1,PtAll_1,"cp"); effVsPt_Pt_1.SetLineColor(2); effVsPt_Pt_1.SetMarkerColor(2); 



   //TFile f2("./output/TMVA_check_ZmumuMC.root");

  TH1F* NvtAll_2 = (TH1F*) _file1->Get("NvtAll");
  TH1F* NvtPt2_2 = (TH1F*) _file1->Get("NvtGood");
  TH1F* NvtBDT_2 = (TH1F*) _file1->Get("NvtGood_BDT");
  
  TH1F* PtAll_2 = (TH1F*) _file1->Get("PtAll");
  TH1F* PtPt2_2 = (TH1F*) _file1->Get("PtGood");
  TH1F* PtBDT_2 = (TH1F*) _file1->Get("PtGood_BDT");

  TGraphAsymmErrors effVsNv_Pt_2;
  TGraphAsymmErrors effVsNv_BDT_2;
  
  NvtPt2_2->Sumw2();  NvtBDT_2->Sumw2(); NvtAll_2->Sumw2();
  
  effVsNv_BDT_2.BayesDivide(NvtBDT_2,NvtAll_2, "cp"); effVsNv_BDT_2.SetLineColor(6); effVsNv_BDT_2.SetMarkerColor(6); 
  effVsNv_Pt_2.BayesDivide(NvtPt2_2,NvtAll_2, "cp"); effVsNv_Pt_2.SetLineColor(4); effVsNv_Pt_2.SetMarkerColor(4); 
  

  TGraphAsymmErrors effVsPt_Pt_2;
  TGraphAsymmErrors effVsPt_BDT_2;

   PtPt2_2->Sumw2();  PtBDT_2->Sumw2(); PtAll_2->Sumw2();

   effVsPt_BDT_2.BayesDivide(PtBDT_2,PtAll_2, "cp"); effVsPt_BDT_2.SetLineColor(6); effVsPt_BDT_2.SetMarkerColor(6); 
   effVsPt_Pt_2.BayesDivide(PtPt2_2,PtAll_2, "cp"); effVsPt_Pt_2.SetLineColor(4); effVsPt_Pt_2.SetMarkerColor(4); 

  
   int nMBMax =12;

   TLegend legend1(0.68, 0.78, 0.99, 0.99);
   legend1.SetFillColor(kWhite);

   TCanvas c1;
   TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.2,1.1);
   cc.SetStats(0); cc.GetXaxis()->SetTitle("number of PV"); cc.GetYaxis()->SetTitle("efficiency");
   cc.Draw();
   
   effVsNv_Pt_1.Draw("p,same");
   legend1.AddEntry(&effVsNv_Pt_1,"Pt_1","LP");

   effVsNv_BDT_1.Draw("p,same");
   legend1.AddEntry(&effVsNv_BDT_1,"BDT_1","LP");

   effVsNv_Pt_2.Draw("p,same");
   legend1.AddEntry(&effVsNv_Pt_2,"Pt_2","LP");

   effVsNv_BDT_2.Draw("P,same");
   legend1.AddEntry(&effVsNv_BDT_2,"BDT_2","LP");

   legend1.Draw();
   c1.SetGridx();
   c1.SetGridy();

   TLegend legend2(0.68, 0.78, 0.99, 0.99);
   legend2.SetFillColor(kWhite);

   TCanvas c2;
   TH2F dd("dd","",200,0,200,1000,0.2,1.1);
   dd.SetStats(0); dd.GetXaxis()->SetTitle("Pt of the boson "); dd.GetYaxis()->SetTitle("efficiency");
   dd.Draw();

   
   effVsPt_Pt_1.Draw("P,same");
   legend2.AddEntry(&effVsPt_Pt_1,"Pt_1","LP");

   effVsPt_BDT_1.Draw("P,same");
   legend2.AddEntry(&effVsPt_BDT_1,"BDT_1","LP");

   effVsPt_Pt_2.Draw("P,same");
   legend2.AddEntry(&effVsPt_Pt_2,"Pt_2","LP");

   effVsPt_BDT_2.Draw("P,same");
   legend2.AddEntry(&effVsPt_BDT_2,"BDT_2","LP");
   
   legend2.Draw();
   c2.SetGridx();
   c2.SetGridy();
}
