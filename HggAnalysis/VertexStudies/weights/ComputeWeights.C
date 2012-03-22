//****** simple macro to compute PU weights ******
{
  //*** mc file
  TFile *fmc = TFile::Open("HistosForPUReweighting_DYJetsToLL_Fall11_S6.root");
  TH1F *hmc  = (TH1F*)fmc->Get("hnpumc");
 
  //*** data file 
  //TFile *fda = TFile::Open("/afs/cern.ch/user/a/adavidzh/public/json/111105/pileup/2011A_0100_73500.pileup.root");
  TFile *fda = TFile::Open("/afs/cern.ch/user/a/adavidzh/public/json/111105/pileup/2011B_0100_73500.pileup.root");
  //  TFile *fda = TFile::Open("/afs/cern.ch/user/a/adavidzh/public/json/111105/pileup/2011_0100_73500.pileup.root");

  TH1F *pileup = (TH1F*)fda->Get("pileup");

  TH1F *hdata = new TH1F("hdata","hdata",60,0,60);
  for (int ibin = 1; ibin < 61; ibin++){
    hdata->SetBinContent(ibin, pileup->GetBinContent(ibin));
  }
 

  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  //TFile *fout = new TFile("./PUweights_2011A_0100_73500_DYJetsToLL_Fall11_S6.root","create");
  TFile *fout = new TFile("./PUweights_2011B_0100_73500_DYJetsToLL_Fall11_S6.root","create");
  //TFile *fout = new TFile("./PUweights_2011_0100_73500_DYJetsToLL_Fall11_S6.root","create");
  hweights->Write("hweights");


}
