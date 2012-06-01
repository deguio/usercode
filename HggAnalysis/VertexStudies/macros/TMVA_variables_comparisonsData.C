{

  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  // 1 : MC
  // 2 : std

  std::string cutString1 = "isSig == 1 && diphopt > 0 && diphopt < 1000";
  std::string cutString2 = "isSig == 1 && diphopt > 0 && diphopt < 1000";

  int mycolor = 2;

  TTree* ntu1 = (TTree*)_file0->Get("TMVA_vertexTree");
  TTree* ntu2 = (TTree*)_file1->Get("TMVA_vertexTree");
  
  TH1F logSumPt2_1("logSumPt2_1","log(sumPt2)",50,-5.,15.);
  logSumPt2_1.GetXaxis()->SetTitle("log(SumPt2)");
  logSumPt2_1.GetYaxis()->SetTitle("a.u.");
  logSumPt2_1.SetLineColor(mycolor);
  logSumPt2_1.SetFillColor(mycolor);  
  logSumPt2_1.SetFillStyle(3004);

  TH1F logSumPt2_2("logSumPt2_2","log(sumPt2)",50,-5.,15.);
  logSumPt2_2.GetXaxis()->SetTitle("log(SumPt2)");
  logSumPt2_2.GetYaxis()->SetTitle("a.u.");
  logSumPt2_2.SetLineColor(kBlack);

  TH1F tracksNum1("tracksNum1","Number of tracks",50,0.,200.);
  tracksNum1.GetXaxis()->SetTitle("number of tracks");
  tracksNum1.GetYaxis()->SetTitle("a.u.");
  tracksNum1.SetLineColor(mycolor);
  tracksNum1.SetFillColor(mycolor);
  tracksNum1.SetFillStyle(3004);

  TH1F tracksNum2("tracksNum2","Number of tracks",50,0.,200.);
  tracksNum2.GetXaxis()->SetTitle("number of tracks");
  tracksNum2.GetYaxis()->SetTitle("a.u.");
  tracksNum2.SetLineColor(kBlack);
  tracksNum2.SetFillColor(kBlack);
  tracksNum2.SetFillStyle(3004);


  TH1F ptbal1("ptbal1","ptbal1",100,-100.,150.);
  ptbal1.GetXaxis()->SetTitle("ptbal");
  ptbal1.GetYaxis()->SetTitle("a.u.");
  ptbal1.SetLineColor(mycolor);  
  ptbal1.SetFillColor(mycolor);  
  ptbal1.SetFillStyle(3004);
  
  TH1F ptbal2("ptbal2","ptbal",100,-100.,150.);
  ptbal2.GetXaxis()->SetTitle("ptbal");
  ptbal2.GetYaxis()->SetTitle("a.u.");
  ptbal2.SetLineColor(kBlack);
  ptbal2.SetFillColor(kBlack);
  ptbal2.SetFillStyle(3004);

  TH1F ptasymm1("ptasymm1","ptasymm",50,-1.,1.);
  ptasymm1.GetXaxis()->SetTitle("ptasymm");
  ptasymm1.GetYaxis()->SetTitle("a.u.");
  ptasymm1.SetLineColor(mycolor);  
  ptasymm1.SetFillColor(mycolor);  
  ptasymm1.SetFillStyle(3004);
  
  TH1F ptasymm2("ptasymm2","ptasymm",50,-1.,1.);
  ptasymm2.GetXaxis()->SetTitle("ptasymm");
  ptasymm2.GetYaxis()->SetTitle("a.u.");
  ptasymm2.SetLineColor(kBlack);
  ptasymm2.SetFillColor(kBlack);
  ptasymm2.SetFillStyle(3004);

  TH1F nvtx1("nvtx1","number of vertexes",50,0.,50.);
  nvtx1.GetXaxis()->SetTitle("number of vertices");
  nvtx1.GetYaxis()->SetTitle("a.u.");
  nvtx1.SetLineColor(mycolor);
  nvtx1.SetFillColor(mycolor);  
  nvtx1.SetFillStyle(3004);

  TH1F nvtx2("nvtx2","number of vertexes",50,0.,50.);
  nvtx2.GetXaxis()->SetTitle("nvtx");
  nvtx2.GetYaxis()->SetTitle("a.u.");
  nvtx2.SetLineColor(kBlack);
  nvtx2.SetFillColor(kBlack);  
  nvtx2.SetFillStyle(3004);


  TH1F bosonPt1("bosonPt1","boson Pt",100,0.,200.);
  bosonPt1.GetXaxis()->SetTitle("p_{T}(GeV)");
  bosonPt1.GetYaxis()->SetTitle("a.u.");
  bosonPt1.SetLineColor(mycolor);
  bosonPt1.SetFillColor(mycolor);  
  bosonPt1.SetFillStyle(3004);

  TH1F bosonPt2("bosonPt2","boson Pt",100,0.,200.);
  bosonPt2.GetXaxis()->SetTitle("p_{T} (GeV)");
  bosonPt2.GetYaxis()->SetTitle("a.u.");
  bosonPt2.SetLineColor(kBlack);
  bosonPt2.SetFillColor(kBlack);  
  bosonPt2.SetFillStyle(3004);

  TH1F invariantMass1("invariantMass1","boson Pt",100,70.,110.);
  invariantMass1.GetXaxis()->SetTitle("invariant mass (GeV)");
  invariantMass1.GetYaxis()->SetTitle("a.u.");
  invariantMass1.SetLineColor(mycolor);
  invariantMass1.SetFillColor(mycolor);  
  invariantMass1.SetFillStyle(3004);

  TH1F invariantMass2("invariantMass2","boson Pt",100,70.,110.);
  invariantMass2.GetXaxis()->SetTitle("invariant mass (GeV)");
  invariantMass2.GetYaxis()->SetTitle("a.u.");
  invariantMass2.SetLineColor(kBlack);
  invariantMass2.SetFillColor(kBlack);  
  invariantMass2.SetFillStyle(3004);

  ntu1->Draw("logsumpt2 >> logSumPt2_1",cutString1.c_str(),"goff");
  ntu2->Draw("logsumpt2 >> logSumPt2_2",cutString2.c_str(),"goff");

  ntu1->Draw("nch >> tracksNum1",cutString1.c_str(),"goff");
  ntu2->Draw("nch >> tracksNum2",cutString2.c_str(),"goff");

  ntu1->Draw("ptbal >> ptbal1",cutString1.c_str(),"goff");
  ntu2->Draw("ptbal >> ptbal2",cutString2.c_str(),"goff");

  ntu1->Draw("ptasym >> ptasymm1",cutString1.c_str(),"goff");
  ntu2->Draw("ptasym >> ptasymm2",cutString2.c_str(),"goff");
 
  ntu1->Draw("diphopt >> bosonPt1",cutString1.c_str(),"goff");
  ntu2->Draw("diphopt >> bosonPt2",cutString2.c_str(),"goff");

  ntu1->Draw("(diphoM*1.0) >> invariantMass1",cutString1.c_str(),"goff");
  ntu2->Draw("diphoM >> invariantMass2",cutString2.c_str(),"goff");

  ntu1->Draw("nVertices >> nvtx1",cutString1.c_str(),"goff");
  ntu2->Draw("nVertices >> nvtx2",cutString2.c_str(),"goff");


  int n1 = tracksNum1->GetEntries();
  int n2 = tracksNum2->GetEntries();

  tracksNum1.SetNormFactor(1.);
  logSumPt2_1.SetNormFactor(1);
  ptbal1.SetNormFactor(1);
  ptasymm1.SetNormFactor(1);
  nvtx1.SetNormFactor(1);
  bosonPt1.SetNormFactor(1);
  invariantMass1.SetNormFactor(1);

  tracksNum2.SetNormFactor(1);
  logSumPt2_2.SetNormFactor(1);
  ptbal2.SetNormFactor(1);
  ptasymm2.SetNormFactor(1);
  nvtx2.SetNormFactor(1);
  bosonPt2.SetNormFactor(1);
  invariantMass2.SetNormFactor(1);

  TLegend leg (0.6, 0.6,0.89,0.75);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(&tracksNum1,"MC Z#rightarrow#mu#mu","F");
  leg->AddEntry(&tracksNum2,"Data Z#rightarrow#mu#mu","LP");

  TLatex *latex = new TLatex(0.55,0.85,"#splitline{CMS preliminary}{#sqrt{s}=8TeV  L=1.921 fb^{-1}}");
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);

  
  TCanvas c1("c1","c1",500,500);
  tracksNum1.Draw("");
  tracksNum2.Draw("esame");
  leg ->Draw("same");
  latex->Draw("same");

  TCanvas c3("c3","c3",500,500);
  logSumPt2_1.Draw("");
  logSumPt2_2.Draw("esame");
  leg->Draw("same");  
  latex->Draw("same");

  TCanvas c5("c5","c5",500,500);
  ptbal1.GetXaxis()->SetRangeUser(-50.,100);
  ptbal1.Draw("");
  ptbal2.Draw("esame");
  leg->Draw("same");  
  latex->Draw("same");

  TCanvas c6("c6","c6",500,500);
  ptasymm1.Draw("");
  ptasymm2.Draw("esame");
  leg->Draw("same");  
  latex->Draw("same");
  
  TCanvas c7("c7","c7",500,500);
  //  c7->SetLogy();
  nvtx1.Draw("");
  nvtx2.Draw("esame");
  leg->Draw("same");  
  latex->Draw("same");
  latex->Draw("same");
  float p = nvtx1.KolmogorovTest(&nvtx2,"");
  cout << "DATA < N_vtx > = " << nvtx2.GetMean() << endl;
  cout << "MC   < N_vtx > = " << nvtx1.GetMean() << endl;
  cout << "NVTX Kolmogorov : p = " << p << endl;
  
  TCanvas c8("c8","c8",500,500);
  bosonPt1.Draw("");
  bosonPt2.Draw("esame");
  leg->Draw("same");    
  latex->Draw("same");


  TCanvas c9("c9","c9",500,500);
  invariantMass1.Draw("");
  invariantMass2.Draw("esame");
  leg->Draw("same");    
  latex->Draw("same");
  


  //  TFile *fout = new TFile("NvtxZmumu_mc_69400.root","create");
  //  nvtx1->Write("nvtx");


}
