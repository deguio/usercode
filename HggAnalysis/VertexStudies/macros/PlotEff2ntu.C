{
  //-----------------------
  //--- firstFile (Hgg) ---
  //-----------------------
  TFile* fileHgg = new TFile("output/NtupleForTTF_Hgg.root","READ");
  TTree* ntu_1 = (TTree*)fileHgg->Get("TMVA_vertexTree");

  //eff VS bosonPt
  TH2F* SumPt2_vsBosonPt_1 = new TH2F("SumPt2_vsBosonPt_1","SumPt2_vsBosonPt_1", 50, 0., 200., 500, 0., 5000.);
  ntu_1->Draw("sumPt2[0]:sum2PhoPt >> SumPt2_vsBosonPt_1","isSig[0] == 1","goff"); //possible bias if I ask signal as first vtx

  TH1F* PtAll_1 = new TH1F("PtAll_1","PtAll_1", 50, 0., 200.);
  ntu_1->Draw("sum2PhoPt >> PtAll_1","","goff");

  TH1F* PtPt2_1 = new TH1F("PtPt2_1","PtPt2_1", 50, 0., 200.);
  ntu_1->Draw("sum2PhoPt >> PtPt2_1","isSig[0] == 1","goff");

  TGraphAsymmErrors effVsPt_Pt_1;
  PtPt2_1->Sumw2();  PtAll_1->Sumw2();

  //------------------------
  //--- secondFile (Zee) ---
  //------------------------
  TFile* fileZee = new TFile("output/NtupleForTTF_ZeeMC.root","READ");
  TTree* ntu_2 = (TTree*)fileZee->Get("TMVA_vertexTree");

  //eff VS bosonPt
  TH2F* SumPt2_vsBosonPt_2 = new TH2F("SumPt2_vsBosonPt_2","SumPt2_vsBosonPt_2", 50, 0., 200., 500, 0., 5000.);
  ntu_2->Draw("sumPt2[0]:sum2PhoPt >> SumPt2_vsBosonPt_2","isSig[0] == 1","goff");

  TH1F* PtAll_2 = new TH1F("PtAll_2","PtAll_2", 50, 0., 200.);
  //ntu_2->Draw("sum2PhoPt >> PtAll_2","","goff");

  TH1F* PtPt2_2 = new TH1F("PtPt2_2","PtPt2_2", 50, 0., 200.); 
  //ntu_2->Draw("sum2PhoPt >> PtPt2_2","isSig == 1","goff");
 
  TGraphAsymmErrors effVsPt_Pt_2;
  PtPt2_2->Sumw2();  PtAll_2->Sumw2();

  //----------------------------------------
  //--- setup branches and loop over Zee ---
  //----------------------------------------
  int isSig[100], nVertices;
  float sumPt2[100], sum2PhoPt;

  ntu_2->SetBranchAddress("sum2PhoPt", &sum2PhoPt);
  ntu_2->SetBranchAddress("nVertices", &nVertices);
  
  ntu_2->SetBranchAddress("sumPt2"   , sumPt2);
  ntu_2->SetBranchAddress("isSig"    , isSig);

  int nBosonPtBins = 1;
  std::cout << "Entries: " << ntu_2->GetEntries() << std::endl;

  for(int ii = 0; ii < ntu_2->GetEntries(); ++ii)
  //for(int ii = 0; ii < 50000; ++ii)
    {
      if (ii%10000 == 0) std::cout << "event n. " << ii << std::endl;
      
      ntu_2->GetEntry(ii);
      
      double binWidth = SumPt2_vsBosonPt_2->GetNbinsX() / nBosonPtBins;
      int bin = (int)(SumPt2_vsBosonPt_2->GetXaxis()->FindBin(sum2PhoPt) / binWidth);

      //TH1D* SumPt2_bin_1 = SumPt2_vsBosonPt_1->ProjectionY("SumPt2_bin_1", bin*binWidth, (bin+1)*binWidth);
      //TH1D* SumPt2_bin_2 = SumPt2_vsBosonPt_2->ProjectionY("SumPt2_bin_2", bin*binWidth, (bin+1)*binWidth);
      TH1D* SumPt2_bin_1 = SumPt2_vsBosonPt_1->ProjectionY("SumPt2_bin_1", 0, 50);
      TH1D* SumPt2_bin_2 = SumPt2_vsBosonPt_2->ProjectionY("SumPt2_bin_2", 0, 50);

      SumPt2_bin_1 -> Scale( 1./SumPt2_bin_1->Integral() );      
      SumPt2_bin_2 -> Scale( 1./SumPt2_bin_2->Integral() );      

      //--------------------------
      //--- loop over vertices ---
      //--------------------------
      int signalIndex = -1;
      for (int kk = 0; kk < nVertices; ++kk)
	if (isSig[kk] == 1) signalIndex = kk;

      //std::cout << "signalIndex = " << signalIndex << std::endl;


      double ratio = SumPt2_bin_1->Interpolate(sumPt2[signalIndex]) / SumPt2_bin_2->Interpolate(sumPt2[signalIndex]);
      
      PtAll_2 -> Fill(sum2PhoPt,ratio);
      if (signalIndex == 0) PtPt2_2 -> Fill(sum2PhoPt,ratio);


      // std::cout << "sum2PhoPt = " << sum2PhoPt << std:: endl;
      // std::cout << " SumPt2_vsBosonPt_2->GetNbinsX() = " <<  SumPt2_vsBosonPt_2->GetNbinsX() << std::endl;
      // std::cout << "SumPt2_vsBosonPt_2->GetXaxis()->FindBin(sum2PhoPt) = " << SumPt2_vsBosonPt_2->GetXaxis()->FindBin(sum2PhoPt) << std::endl;
      // std::cout << "binWidth = "     << binWidth << std::endl;
      // std::cout << "bin*binWidth = " << bin*binWidth << std::endl;
      // std::cout << "(bin+1)*binWidth = " << (bin+1)*binWidth << std::endl;
      // std::cout << std::endl;

    }
  
  TCanvas c1;
  c1.cd();
  //Hgg
  effVsPt_Pt_1.BayesDivide(PtPt2_1,PtAll_1,"cp");
  effVsPt_Pt_1.SetLineColor(4);
  effVsPt_Pt_1.SetMarkerColor(4);
  //Zee
  effVsPt_Pt_2.BayesDivide(PtPt2_2,PtAll_2,"cp");
  effVsPt_Pt_2.SetLineColor(2);
  effVsPt_Pt_2.SetMarkerColor(2);

  effVsPt_Pt_2.Draw("AP");
  effVsPt_Pt_1.Draw("P,sames");

  TCanvas c2;
  c2.cd();
  SumPt2_bin_2->SetLineColor(2);
  SumPt2_bin_2->Draw();
  SumPt2_bin_1->SetLineColor(4);
  SumPt2_bin_1->Draw("sames");
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------






  // //Draw common
  //  int nMBMax =12;

  //  TLegend legend1(0.68, 0.78, 0.99, 0.99);
  //  legend1.SetFillColor(kWhite);

  //  TCanvas c1;
  //  TH2F cc("cc","",nMBMax+1,0,nMBMax+1,1000,0.2,1.1);
  //  cc.SetStats(0); cc.GetXaxis()->SetTitle("number of PV"); cc.GetYaxis()->SetTitle("efficiency");
  //  cc.Draw();
   
  //  effVsNv_Pt_1.Draw("p,same");
  //  legend1.AddEntry(&effVsNv_Pt_1,"SumPt2 Zee MC","LP");

  //  effVsNv_BDT_1.Draw("p,same");
  //  legend1.AddEntry(&effVsNv_BDT_1,"BDT Zee MC","LP");
   
  //  effVsNv_Pt_2.Draw("p,same");
  //  legend1.AddEntry(&effVsNv_Pt_2,"SumPt2 Hgg","LP");

  //  effVsNv_BDT_2.Draw("P,same");
  //  legend1.AddEntry(&effVsNv_BDT_2,"BDT Hgg","LP");
   
  //  legend1.Draw();
  //  c1.SetGridx();
  //  c1.SetGridy();

  //  TLegend legend2(0.68, 0.78, 0.99, 0.99);
  //  legend2.SetFillColor(kWhite);

  //  TCanvas c2;
  //  TH2F dd("dd","",200,0,200,1000,0.2,1.1);
  //  dd.SetStats(0); dd.GetXaxis()->SetTitle("Pt of the boson "); dd.GetYaxis()->SetTitle("efficiency");
  //  dd.Draw();

   
  //  effVsPt_Pt_1.Draw("P,same");
  //  legend2.AddEntry(&effVsPt_Pt_1,"SumPt2 Zee MC","LP");

  //  effVsPt_BDT_1.Draw("P,same");
  //  legend2.AddEntry(&effVsPt_BDT_1,"BDT Zee MC","LP");

  //  effVsPt_Pt_2.Draw("P,same");
  //  legend2.AddEntry(&effVsPt_Pt_2,"SumPt2 Hgg","LP");

  //  effVsPt_BDT_2.Draw("P,same");
  //  legend2.AddEntry(&effVsPt_BDT_2,"BDT Hgg","LP");
   
  //  legend2.Draw();
  //  c2.SetGridx();
  //  c2.SetGridy();

  //  TCanvas c3;
  //  eff2D_Pt_1 -> SetTitle("SumPt2 Zee MC");
  //  eff2D_Pt_1 -> Draw("colz");
  //  TCanvas c4;
  //  eff2D_BDT_1 -> SetTitle("BDT Zee MC");
  //  eff2D_BDT_1 -> Draw("colz");
  //  TCanvas c5;
  //  eff2D_Pt_2 -> SetTitle("SumPt2 Hgg");   
  //  eff2D_Pt_2 -> Draw("colz");
  //  TCanvas c6;
  //  eff2D_BDT_2 -> SetTitle("BDT Hgg");
  //  eff2D_BDT_2 -> Draw("colz");
}
