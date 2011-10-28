//****** simple macro to compute PU weights ******
{
  //*** mc file
  TFile *fmc = TFile::Open("HistosForPUReweighting_DYJetsToLL_Summer11_S4.root");
  TH1F *hmc  = (TH1F*)fmc->Get("hnpumc");
 
  //*** data file 
  TFile *fda = TFile::Open("./lp.json.root");
  TH1F *pileup = (TH1F*)fda->Get("pileup");

  TH1F *hdata = new TH1F("hdata","hdata",50,0,50);
  for (int ibin = 1; ibin < 51; ibin++){
    hdata->SetBinContent(ibin, pileup->GetBinContent(ibin));
  }
 
  TCanvas *c1 = new TCanvas("c1","c1");
  hmc  ->DrawNormalized();
  hdata->DrawNormalized("same");
  
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  TCanvas *c2 = new TCanvas("c2","c2");
  hweights->Draw("");

  TFile *fout = new TFile("./PUweights_lp_DYJetsToLL_Summer11_S4.root","create");
  hweights->Write("hweights");

}
