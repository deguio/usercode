//****** simple macro to compute PU weights ******
{
  
  

  //*** mc pileup file
  //TFile *fmc = TFile::Open("../pileup/nPU-Summer12_DD3.190456-200245-29Jun_Prompt.observed.root"); // mc DY 53 -- from Shervin
  //TH1F *hmc  = (TH1F*)fmc->Get("hnpumc");
  //cout << hmc->GetNbinsX()<< endl;

  //*** mc ntuple 
  TChain* chain = new TChain("MiBiCommonNTTwoPhotons/SimpleNtuple");
  chain->Add("root://eoscms//eos/cms/store/cmst3/user/malberti/HIGGS/VERTEX/2012/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM/MiBiCommonNT_*.root");
  TH1F *hmc = new TH1F("hmc","hmc",60,0,60);
  cout << hmc->GetNbinsX()<< endl;
  chain->Draw("mc_PUit_NumInteractions>>hmc");

  //*** data file 
  TFile *fda = TFile::Open("../pileup/pileup_190456-203002_minBiasXsec69400_corr_observed.root");
  //TFile *fda = TFile::Open("../pileup/pileup_198041-201229_minBiasXsec69400_corr_observed.root");
  //TFile *fda = TFile::Open("../pileup/pileup_190456-196531_minBiasXsec69400_corr_observed.root");
  TH1F *hdata = (TH1F*)fda->Get("pileup");
  cout << hdata->GetNbinsX()<< endl;
  
  //*** compute weights
  TH1F *hweights = (TH1F*)hdata->Clone("hweights");
  hweights->Divide(hdata,hmc,1./hdata->GetSumOfWeights(),1./hmc->GetSumOfWeights());

  TFile *fout = new TFile("./PUweights_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012ABC.root","create");
  //TFile *fout = new TFile("./PUweights_DYJetsToLL_Summer12_DR53X-PU_S10_minBiasXsec69400_corr_observed_Run2012AB.root","create");
  hweights->Write("hweights");
  hdata->Write("hdata");
  hmc->Write("hmc");


}
