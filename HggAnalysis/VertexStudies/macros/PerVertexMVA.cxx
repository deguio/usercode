
{
  gROOT->LoadMacro("~/setTDRStyle.C");
  setTDRStyle();

  TFile *f[2];
  f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/PerVertexMVA_DYJetsToLL_S9_minBiasXsec68300_190456-194479/PerVertexMVA_DYJetsToLL_S9_minBiasXsec68300_190456-194479.root");  
  f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/PerVertexMVA_DoubleMu_190456-194479/PerVertexMVA_DoubleMu_190456-194479.root");
  
  
  // histograms
  TH1F* tmvaPerVertexOutput_RightVtx[2];
  TH1F* tmvaPerVertexOutput_WrongVtx[2];

  TH1F* tmvaPerVertexOutput[5][2];
  TH2F *ptbal_vs_ptasym[2];
  TH2F *ptbal_vs_sumpt2[2];
  TH2F *ptasym_vs_sumpt2[2];

  int nre = 4;

  for (int i = 0; i < 2; i++){

    ptbal_vs_ptasym[i]  = (TH2F*)f[i]->Get("ptbal_vs_ptasym");
    ptbal_vs_sumpt2[i]  = (TH2F*)f[i]->Get("ptbal_vs_sumpt2");
    ptasym_vs_sumpt2[i] = (TH2F*)f[i]->Get("ptasym_vs_sumpt2");
    
    tmvaPerVertexOutput_RightVtx[i] = (TH1F*)f[i]->Get("tmvaPerVertexOutput_RightVtx");
    tmvaPerVertexOutput_WrongVtx[i] = (TH1F*)f[i]->Get("tmvaPerVertexOutput_WrongVtx");
    tmvaPerVertexOutput_RightVtx[i] -> Sumw2();
    tmvaPerVertexOutput_WrongVtx[i] -> Sumw2();
    tmvaPerVertexOutput_RightVtx[i] -> Rebin(nre);
    tmvaPerVertexOutput_WrongVtx[i] -> Rebin(nre);


    tmvaPerVertexOutput[0][i] = (TH1F*) f[i]->Get("tmvaPerVertexOutput_0");
    tmvaPerVertexOutput[1][i] = (TH1F*) f[i]->Get("tmvaPerVertexOutput_1");
    tmvaPerVertexOutput[2][i] = (TH1F*) f[i]->Get("tmvaPerVertexOutput_2");
    tmvaPerVertexOutput[3][i] = (TH1F*) f[i]->Get("tmvaPerVertexOutput_3");
    tmvaPerVertexOutput[4][i] = (TH1F*) f[i]->Get("tmvaPerVertexOutput_4");

    for (int j = 0; j < 5 ; j++){
      tmvaPerVertexOutput[j][i] -> Sumw2();
      tmvaPerVertexOutput[j][i] -> Rebin(nre);
      
      if (i==0){
	tmvaPerVertexOutput[j][i] -> SetFillColor(kAzure+1);
	tmvaPerVertexOutput_RightVtx[i] -> SetFillColor(kGreen+1);
	tmvaPerVertexOutput_WrongVtx[i] -> SetFillColor(kRed+1);
	tmvaPerVertexOutput_RightVtx[i] -> SetFillStyle(3002);
	tmvaPerVertexOutput_WrongVtx[i] -> SetFillStyle(3005);
      }
      
      if (i==1){
	tmvaPerVertexOutput[j][i] -> SetMarkerStyle(20);
	tmvaPerVertexOutput[j][i] -> SetMarkerSize(0.8);

	tmvaPerVertexOutput_RightVtx[i] -> SetMarkerStyle(20);
	tmvaPerVertexOutput_WrongVtx[i] -> SetMarkerStyle(20);
	tmvaPerVertexOutput_RightVtx[i] -> SetMarkerSize(0.8);
	tmvaPerVertexOutput_WrongVtx[i] -> SetMarkerSize(0.8);
	tmvaPerVertexOutput_RightVtx[i] -> SetMarkerColor(kGreen+2);
	tmvaPerVertexOutput_WrongVtx[i] -> SetMarkerColor(kRed+2);
      }
    }
  }


  TLegend leg (0.6, 0.6,0.89,0.75);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(tmvaPerVertexOutput[0][0],"MC Z#rightarrow#mu#mu","F");
  leg->AddEntry(tmvaPerVertexOutput[0][1],"Data Z#rightarrow#mu#mu","LP");

  TLatex *latex = new TLatex(0.55,0.85,"#splitline{          CMS preliminary}{#sqrt{s} = 8 TeV L = 1.92 fb^{-1}}");
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.04);


  TCanvas *c0 = new TCanvas("c0","c0",500,500);
  tmvaPerVertexOutput[0][0] -> GetXaxis() -> SetTitle("MVA_{vtx}^{0}"); 
  tmvaPerVertexOutput[0][0] -> DrawNormalized("histo");
  tmvaPerVertexOutput[0][1] ->DrawNormalized("esame");
  leg->Draw("same");
  latex->Draw("same");
  
  TCanvas *c1 = new TCanvas("c1","c1",500,500);
  tmvaPerVertexOutput[1][0] -> GetXaxis() -> SetTitle("MVA_{vtx}^{1}"); 
  tmvaPerVertexOutput[1][0] ->DrawNormalized("histo");
  tmvaPerVertexOutput[1][1] ->DrawNormalized("esame");
  leg->Draw("same");
  latex->Draw("same");
  
  TCanvas *c2 = new TCanvas("c2","c2",500,500);
  tmvaPerVertexOutput[2][0] -> GetXaxis() -> SetTitle("MVA_{vtx}^{2}"); 
  tmvaPerVertexOutput[2][0] ->DrawNormalized("histo");
  tmvaPerVertexOutput[2][1] ->DrawNormalized("esame");
  leg->Draw("same");
  latex->Draw("same");


  TLegend leg2 (0.6, 0.6,0.89,0.75);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(tmvaPerVertexOutput_RightVtx[0],"right vertex MC","F");
  leg2->AddEntry(tmvaPerVertexOutput_WrongVtx[0],"wrong vertex MC","F");
  leg2->AddEntry(tmvaPerVertexOutput_RightVtx[1],"right vertex DATA","LP");
  leg2->AddEntry(tmvaPerVertexOutput_WrongVtx[1],"wrong vertex DATA","LP");

  TCanvas *cright = new TCanvas("cright","cright",500,500);
  tmvaPerVertexOutput_RightVtx[0] -> GetXaxis() -> SetTitle("MVA_{vtx}");
  tmvaPerVertexOutput_RightVtx[0] -> DrawNormalized("histo");
  tmvaPerVertexOutput_WrongVtx[0] -> DrawNormalized("histo same");
  tmvaPerVertexOutput_RightVtx[1] -> DrawNormalized("esame");
  tmvaPerVertexOutput_WrongVtx[1] -> DrawNormalized("esame");
  leg2->Draw("same");
  latex->Draw("same");


  // correlation plots
  TTree* ntu1 = (TTree*)f[0]->Get("TMVAtree");
  TTree* ntu2 = (TTree*)f[1]->Get("TMVAtree");

  TH2F *ptbal_vs_ptasym_mc   = new TH2F("ptbal_vs_ptasym_mc","",50,-1,1,200,-100,100);
  TH2F *ptbal_vs_ptasym_data = new TH2F("ptbal_vs_ptasym_data","",50,-1,1,200,-100,100);

  TH2F *diphopt_vs_mva0_mc   = new TH2F("diphopt_vs_mva0_mc","diphopt_vs_mva",50,-1,1,200,0,200);
  TH2F *diphopt_vs_mva0_data = new TH2F("diphopt_vs_mva0_data","diphopt_vs_data",50,-1,1,200,0,200);

 
  ntu1->Draw("ptbal[0]:ptasym[0] >> ptbal_vs_ptasym_mc","","goff");
  ntu2->Draw("ptbal[0]:ptasym[0] >> ptbal_vs_ptasym_data","","goff");

  ntu1->Draw("diphopt[0]:mva[0] >> diphopt_vs_mva0_mc","","goff");
  ntu2->Draw("diphopt[0]:mva[0] >> diphopt_vs_mva0_data","","goff");

  TCanvas *cdiphopt_vs_mva0 = new TCanvas("cdiphopt_vs_mva0","cdiphopt_vs_mva",500,500);
  cdiphopt_vs_mva0->SetRightMargin(4.5);
  cdiphopt_vs_mva0->SetLogz();
  diphopt_vs_mva0_data -> GetXaxis()->SetTitle("MVA_{vtx}^{0}");
  diphopt_vs_mva0_data -> GetYaxis()->SetTitle("boson p_{T} (GeV)");
  diphopt_vs_mva0_data -> Draw("colz");


}
