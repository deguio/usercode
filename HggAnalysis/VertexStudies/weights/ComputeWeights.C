//****** simple macro to compute PU weights ******
{
  //*** mc file
  TFile *fmc = TFile::Open("../weights/HistosForPUReweighting_DYJetsToLL_Summer12_S9.root");
  TH1F *hmc  = (TH1F*)fmc->Get("htruenpumc");
 
  //*** data file 
  TFile *fda = TFile::Open("../pileup/pileup_190389-191859_minBiasXsec62460.root");
  TH1F *hdata = (TH1F*)fda->Get("pileup");
 
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  TFile *fout = new TFile("./PUweights_2012_DYJetsToLL_Summer12_S9_minBiasXsec62460.root","create");
  hweights->Write("hweights");


}
