{
  // gROOT->LoadMacro("~/setTDRStyle.C");
  // setTDRStyle();


  TFile *f = TFile::Open("outputEff/test_globe_Glu121_S6_PUweights_2011_0100_73500.root");
 
  int nRePt = 1;
  
  TH1F *PtAll;
  TH1F *PtGood;
  PtAll  = (TH1F*)f->Get("PtAll_BDT");
  PtGood = (TH1F*)f->Get("PtGood_BDT");

  float ptlow = 0.;
  float pthigh = 400;
  int bin1 = PtAll->FindBin(ptlow);
  int bin2 = PtAll->FindBin(pthigh);

  double eff = PtGood->Integral(bin1,bin2)/ PtAll->Integral(bin1,bin2);
  double err = sqrt(eff*(1-eff)/PtAll->Integral(bin1,bin2));
  
  cout << "Efficiency integrated in boson pt : [" << ptlow << ","<< pthigh<<"] GeV"<< endl;
  cout << "Hgg MC (POWHEG) efficiency = " << eff << " +/- " << err<< endl;  
  cout << endl;

  float cont;
  float pt ;
  //-----------------------------------------------------------------------------------------------------
  //variations for K-factors
  //-----------------------------------------------------------------------------------------------------

  // suppose to have a TGraphError with k-factors
  // new pt spectrum = PtAll x kfact
  // new pt specrtum after vtx id choice = PtAll x kfact x eff_vtxId = PtGood x kfact
  // new efficiency = ratio of the two above  

  TFile *fcorrK = TFile::Open("/afs/cern.ch/user/n/nckw/public/KFactors_AllScales_interpolated.root");
  TH1F * gcorrK[3];
  gcorrK[0] = (TH1F*) fcorrK->Get("kfact120_0");
  gcorrK[1] = (TH1F*) fcorrK->Get("kfact120_1");
  gcorrK[2] = (TH1F*) fcorrK->Get("kfact120_6");

  //TFile *fcorrK = TFile::Open("kfactors/Kfactors_120_AllScales.root");
  //   TH1F * gcorrK[3];
  //   gcorrK[0] = (TH1F*) fcorrK->Get("kfactors/kfact_mh120_ren120_fac120");
  //   gcorrK[1] = (TH1F*) fcorrK->Get("kfactors/kfact_mh120_ren60_fac60");
  //   gcorrK[2] = (TH1F*) fcorrK->Get("kfactors/kfact_mh120_ren240_fac240");

  gcorrK[0]->Rebin(5);
  gcorrK[1]->Rebin(5);
  gcorrK[2]->Rebin(5);

  TH1F *PtGoodK[3] ;
  TH1F *PtAllK[3]; 

  PtGoodK[0] = new TH1F("PtGoodK0", "PtGoodK0",50,0,250);
  PtAllK[0] = new TH1F("PtAllK0", "PtAllK0",50,0,250);

  PtGoodK[1] = new TH1F("PtGoodK1", "PtGoodK1",50,0,250);
  PtAllK[1] = new TH1F("PtAllK1", "PtAllK1",50,0,250);

  PtGoodK[2] = new TH1F("PtGoodK6", "PtGoodK6",50,0,250);
  PtAllK[2] = new TH1F("PtAllK6", "PtAllK6",50,0,250);
  
  float kfact[3];

  for (int jscale=0 ; jscale<3; jscale++){ 
    for (int ibin = 1; ibin < PtAll->GetNbinsX()+1; ibin++ ){
      pt   = PtGood->GetBinCenter(ibin);
      int mybinK = gcorrK[jscale]->FindBin(pt);
      kfact[jscale] = gcorrK[jscale]->GetBinContent(mybinK);
      PtGoodK[jscale]->SetBinContent(ibin, kfact[jscale]*PtGood->GetBinContent(ibin));
      PtAllK[jscale] ->SetBinContent(ibin, kfact[jscale]*PtAll->GetBinContent(ibin));
    }
  } 


  float effK[3];
  float errK[3];

  for (int jscale=0 ; jscale<3; jscale++){ 
    effK[jscale] = PtGoodK[jscale]->Integral(bin1,bin2)/ PtAllK[jscale]->Integral(bin1,bin2);
    errK[jscale] = sqrt(effK[jscale]*(1-effK[jscale])/ PtAllK[jscale]->Integral(bin1,bin2));
  
  }

  cout << "**** Variations for NLO-NNLO k-factors:" << endl;
  cout << "Scaled eff = " << effK[0] << " +/- " << err*effK[0]/eff<< endl;  
  
  cout << "**** Uncertainty from k-factors (varying mu_R, mu_F):" << endl;
  cout << "mu x 2 --> Scaled eff = " << effK[1] << " +/- " << err*effK[1]/eff<< endl;  
  cout << "  -->> Vtx ID eff variation :" << effK[0]-effK[1] << endl;


  cout << "**** Uncertainty from k-factors (varying mu_R, mu_F):" << endl;
  cout << "mu x 0.5 --> Scaled eff = " << effK[2] << " +/- " << err*effK[2]/eff<< endl;  
  cout << "   -->> Vtx ID eff variation :" << effK[0]-effK[2] << endl;
  cout << endl;

  //-----------------------------------------------------------------------------------------------------
  // --- variations for measured vtxId efficiency
  //-----------------------------------------------------------------------------------------------------
  TFile *fcorr = TFile::Open("scaleFactors/BDT_vtxIdScaleFactorFromZmumu_DYJetsToLL_Fall11_S6_Run2011all.root"); // 

  TGraphErrors *gcorr = (TGraphErrors*) fcorr->Get("scaleFactor");   
  
  TH1F *Pt[3];
  char hname[100];
  for (int i = 0 ; i < 3; i++){
    sprintf(hname,"PtGood_%d",i);
    Pt[i] = new TH1F(hname, hname,50,0,250);
  }

  TH1F *PtAllZ = new TH1F("PtAllZ","PtAllZ",50,0,250);

  float corr[3];

   for (int ibin = 1; ibin < PtAll->GetNbinsX()+1; ibin++ ){
    cont = PtGoodK[0]->GetBinContent(ibin);
    pt   = PtGoodK[0]->GetBinCenter(ibin);
    int mybin = gcorr->GetHistogram()->FindBin(pt);
    corr[0] = gcorr->GetY()[ mybin ];
    corr[1] = corr[0]+gcorr->GetErrorY( mybin );
    corr[2] = corr[0]-gcorr->GetErrorY( mybin );
    Pt[0]->SetBinContent(ibin, cont * corr[0] );
    Pt[1]->SetBinContent(ibin, cont * corr[1] );
    Pt[2]->SetBinContent(ibin, cont * corr[2] );
    PtAllZ->SetBinContent(ibin, PtAllK[0]->GetBinContent(ibin));
  }



  double effCorr[3];
  double errCorr[3];
  
  cout << "**** Uncertainties for measured vtxId eff in DATA:" << endl;
  for (int i=0; i < 3; i++){
    effCorr[i] = Pt[i]->Integral(bin1,bin2)/ PtAllZ->Integral(bin1,bin2);
    errCorr[i] = sqrt(effCorr[i]*(1-effCorr[i])/ PtAllZ->Integral(bin1,bin2));
 }

  cout << "Scaled (0sigma)  eff = " << effCorr[0] << " +/- " << err*effCorr[0]/eff<< endl;  
  cout << "Scaled (+1sigma) eff = " << effCorr[1] << " +/- " << err*effCorr[1]/eff<< endl;  
  cout << "Scaled (-1sigma) eff = " << effCorr[2] << " +/- " << err*effCorr[2]/eff<< endl;  

  cout << "  -->> Vtx ID eff uncertainty :" << effCorr[0]-effCorr[1] << "  "  << effCorr[0]-effCorr[2] << endl;
  cout << endl;


  cout << "Total efficiency (after corrections) = " << effCorr[0] 
       << "+/-" << err*effCorr[0]/eff << " (stat.)  +/-" 
       << sqrt(pow((effK[1]-effK[2])/2,2) +  pow( (effCorr[1]-effCorr[2])/2   ,2)) << "(syst.)" << endl;

}
