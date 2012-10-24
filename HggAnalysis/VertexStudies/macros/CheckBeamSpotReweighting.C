
void CheckBeamSpotReweighting()

{
  TFile *f[3];

  f[0] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DoubleMu_Run2012ABC_dz1/TMVA_check_DoubleMu_Run2012ABC.root");
  f[1] = TFile::Open("/afs/cern.ch/work/m/malberti/private/test/testEfficiency.root");
  f[2] = TFile::Open("/afs/cern.ch/work/m/malberti/private/Eff_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABC_dz1/testEfficiency.root");
  
  
  TH1F *h[3];
  TH1F *hh[3];
 
  for (int i=0; i<3 ; i++){
    h[i]  = (TH1F*)f[i]->Get("ChosenVertexDz_BDT");
    h[i]->SetMarkerStyle(20);
    h[i]->SetMarkerSize(0.2);
  }
  
  hh[0]=(TH1F*)h[0]->Clone("hdata");
  hh[1]=(TH1F*)h[1]->Clone("hmc");
  hh[2]=(TH1F*)h[2]->Clone("hmc_noRW");
   

  for (int i=0; i<3 ; i++){
    hh[i]->Rebin(100);
    
  }
  
  int b1 =   hh[0]->FindBin(-0.1);
  int b2 =   hh[0]->FindBin(0.1);

//   for (int ibin=1; ibin<hh[0]->GetNbinsX()+1;ibin++){
//     if (ibin<b1 || ibin > b2){
//       hh[0]->SetBinContent(ibin,0);
//       hh[1]->SetBinContent(ibin,0);
//     }
//   }
   for (int ibin=b1; ibin<b2+1;ibin++){
      hh[0]->SetBinContent(ibin,0);
      hh[1]->SetBinContent(ibin,0);
      hh[2]->SetBinContent(ibin,0);
   }

  
  float newBSmean2  = 1.8902e-01;
  float newBSnorm2  = 4.1813e+03;
  float newBSsigma2 = 7.0811e+00;

  float oldBSmean2  = 4.9986e-01;
  float oldBSnorm2  = 4.0258e+01;
  float oldBSsigma2 = 8.5356e+00;

  TF1 *fda = new TF1("fda","gaus",-20,20);
  fda->SetParameters(newBSnorm2,newBSmean2,newBSsigma2);
  fda->SetLineColor(1);
  TF1 *fmc = new TF1("fmc","gaus",-20,20);
  fmc->SetParameters(oldBSmean2,oldBSnorm2,oldBSsigma2);

  TF1 *ff = new TF1("ff","fda/fmc",-20,20);
  
  
  float nda  = hh[0]->GetSumOfWeights();
  float nmcr = hh[1]->GetSumOfWeights();
  float nmc  = hh[2]->GetSumOfWeights();
 
  cout <<  nda  << "   " << nmcr << "   " <<  nmc  << endl;
  
  hh[0]->Sumw2();
  hh[1]->Sumw2();
  hh[2]->Sumw2();

  hh[0]->Scale(1./nda);
  hh[1]->Scale(1./nmcr);
  hh[2]->Scale(1./nmc);

//   cout << "Fitting data" << std::endl;
//   hh[0]->Fit("fda","RS+");

//   cout << "Fitting MC" << std::endl;
//   hh[1]->Fit("fmc","RS+");

  
  TCanvas *c = new TCanvas("c","c",800,600);
  hh[0]->GetXaxis()->SetRangeUser(-30,30);
  hh[0]->GetXaxis()->SetTitle("z_{chosen}-z_{true} (cm)");
  hh[0]->GetYaxis()->SetTitle("a.u.");
  hh[0]->SetLineColor(1);
  hh[0]->SetMarkerColor(1);
  hh[0]->Draw("es");
  hh[1]->SetLineColor(kGreen+1);
  hh[1]->SetMarkerColor(kGreen+1);
  hh[1]->Draw("esames");
  hh[2]->SetLineColor(kRed);
  hh[2]->SetMarkerColor(kRed);
  hh[2]->Draw("esames");


  TLegend *l = new TLegend(0.12,0.7,0.4,0.89);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(hh[0], "data", "PL" );
  l->AddEntry(hh[2], "MC w/o beam spot re-weighting", "PL" );
  l->AddEntry(hh[1], "MC w/  beam spot re-weighting", "PL" );
  l->Draw("same");

}
