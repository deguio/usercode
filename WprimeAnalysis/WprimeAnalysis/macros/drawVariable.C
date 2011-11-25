//simple macro to draw the Wenu - multijet comparison for eleID variables
{

  //Get File
  TFile* Wenu_file = new TFile("/gwteraz/users/deguio/MiBiCommonPAT/Summer11/MC_16062011/WToENu_TuneZ2_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2_AODSIM/MiBiCommonNT.root","READ");
  TFile* QCD_30_80_file = new TFile("/gwteraz/users/deguio/MiBiCommonPAT/Summer11/MC_16062011/QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia_Summer11-PU_S4_START42_V11-v1_AODSIM/MiBiCommonNT.root","READ");

  //Get Tree
  TTree* Wenu      = (TTree*)Wenu_file      ->Get("MiBiCommonNTOneElectron/SimpleNtuple");
  TTree* QCD_30_80 = (TTree*)QCD_30_80_file ->Get("MiBiCommonNTOneElectron/SimpleNtuple");

//   //pt
//   TH1F* pt_Wenu      = new TH1F("pt_Wenu","pt_Wenu",1000,0.,1000);
//   pt_Wenu->GetXaxis()->SetTitle("p_{T} [GeV]");
//   pt_Wenu->GetYaxis()->SetTitle("a.u.");
//   pt_Wenu->SetLineColor(2);
//   pt_Wenu->SetLineWidth(2);
  

//   TH1F* pt_QCD_30_80 = new TH1F("pt_QCD_30_80","pt_QCD_30_80",1000,0.,1000);
//   pt_QCD_30_80->GetXaxis()->SetTitle("p_{T} [GeV]");
//   pt_QCD_30_80->GetYaxis()->SetTitle("a.u.");
//   pt_QCD_30_80->SetLineColor(4);
//   pt_QCD_30_80->SetLineWidth(2);

//   Wenu     ->Draw("electrons[0]->pt() >> pt_Wenu","electrons_isEB == 1","goff");
//   QCD_30_80->Draw("electrons[0]->pt() >> pt_QCD_30_80","electrons_isEB == 1","goff");
  
//   TCanvas pt;
//   pt_Wenu->DrawNormalized();
//   pt_QCD_30_80->DrawNormalized("sames");
//   pt.Print("pt.pdf");

//   //deltaEtaIn
//   TH1F* deltaEtaIn_Wenu      = new TH1F("deltaEtaIn_Wenu","deltaEtaIn_Wenu",600,-0.02,0.02);
//   deltaEtaIn_Wenu->GetXaxis()->SetTitle("#Delta#eta_{in}");
//   deltaEtaIn_Wenu->GetYaxis()->SetTitle("a.u.");
//   deltaEtaIn_Wenu->SetLineColor(2);
//   deltaEtaIn_Wenu->SetLineWidth(2);;

//   TH1F* deltaEtaIn_QCD_30_80 = new TH1F("deltaEtaIn_QCD_30_80","deltaEtaIn_QCD_30_80",600,-0.03,0.03);
//   deltaEtaIn_QCD_30_80->GetXaxis()->SetTitle("#Delta#eta_{in}");
//   deltaEtaIn_QCD_30_80->GetYaxis()->SetTitle("a.u.");
//   deltaEtaIn_QCD_30_80->SetLineColor(4);
//   deltaEtaIn_QCD_30_80->SetLineWidth(2);

//   Wenu     ->Draw("electrons_deltaEtaIn[0] >> deltaEtaIn_Wenu","electrons_isEB == 1","goff");
//   QCD_30_80->Draw("electrons_deltaEtaIn[0] >> deltaEtaIn_QCD_30_80","electrons_isEB == 1","goff");
  
//   TCanvas deltaEtaIn;
//   deltaEtaIn_Wenu->DrawNormalized();
//   deltaEtaIn_QCD_30_80->DrawNormalized("sames");
//   deltaEtaIn.Print("deltaEtaIn.pdf");


//   //deltaPhiIn
//   TH1F* deltaPhiIn_Wenu      = new TH1F("deltaPhiIn_Wenu","deltaPhiIn_Wenu",600,-0.15,0.15);
//   deltaPhiIn_Wenu->GetXaxis()->SetTitle("#Delta#phi_{in}");
//   deltaPhiIn_Wenu->GetYaxis()->SetTitle("a.u.");
//   deltaPhiIn_Wenu->SetLineColor(2);
//   deltaPhiIn_Wenu->SetLineWidth(2);

//   TH1F* deltaPhiIn_QCD_30_80 = new TH1F("deltaPhiIn_QCD_30_80","deltaPhiIn_QCD_30_80",600,-0.15,0.15);
//   deltaPhiIn_QCD_30_80->GetXaxis()->SetTitle("#Delta#phi_{in}");
//   deltaPhiIn_QCD_30_80->GetYaxis()->SetTitle("a.u.");
//   deltaPhiIn_QCD_30_80->SetLineColor(4);
//   deltaPhiIn_QCD_30_80->SetLineWidth(2);

//   Wenu     ->Draw("electrons_deltaPhiIn[0] >> deltaPhiIn_Wenu","electrons_isEB == 1","goff");
//   QCD_30_80->Draw("electrons_deltaPhiIn[0] >> deltaPhiIn_QCD_30_80","electrons_isEB == 1","goff");
  
//   TCanvas deltaPhiIn;
//   deltaPhiIn_Wenu->DrawNormalized();
//   deltaPhiIn_QCD_30_80->DrawNormalized("sames");
//   deltaPhiIn.Print("deltaPhiIn.pdf");

//   //sigmaIEtaIEta
//   TH1F* sigmaIEtaIEta_Wenu      = new TH1F("sigmaIEtaIEta_Wenu","sigmaIEtaIEta_Wenu",200,0.,0.02);
//   sigmaIEtaIEta_Wenu->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
//   sigmaIEtaIEta_Wenu->GetYaxis()->SetTitle("a.u.");
//   sigmaIEtaIEta_Wenu->SetLineColor(2);
//   sigmaIEtaIEta_Wenu->SetLineWidth(2);

//   TH1F* sigmaIEtaIEta_QCD_30_80 = new TH1F("sigmaIEtaIEta_QCD_30_80","sigmaIEtaIEta_QCD_30_80",200,0.,0.02);
//   sigmaIEtaIEta_QCD_30_80->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");
//   sigmaIEtaIEta_QCD_30_80->GetYaxis()->SetTitle("a.u.");
//   sigmaIEtaIEta_QCD_30_80->SetLineColor(4);
//   sigmaIEtaIEta_QCD_30_80->SetLineWidth(2);

//   Wenu     ->Draw("electrons_sigmaIetaIeta[0] >> sigmaIEtaIEta_Wenu","electrons_isEB == 1","goff");
//   QCD_30_80->Draw("electrons_sigmaIetaIeta[0] >> sigmaIEtaIEta_QCD_30_80","electrons_isEB == 1","goff");

//   TCanvas sigmaIEtaIEta;
//   sigmaIEtaIEta_Wenu->DrawNormalized();
//   sigmaIEtaIEta_QCD_30_80->DrawNormalized("sames");
//   sigmaIEtaIEta.Print("sigmaIEtaIEta.pdf");


  //HoverE
  TH1F* HoverE_Wenu      = new TH1F("HoverE_Wenu","HoverE_Wenu",300,0.,0.15);
  HoverE_Wenu->GetXaxis()->SetTitle("H/E");
  HoverE_Wenu->GetYaxis()->SetTitle("a.u.");
  HoverE_Wenu->SetLineColor(2);
  HoverE_Wenu->SetLineWidth(2);

  TH1F* HoverE_QCD_30_80 = new TH1F("HoverE_QCD_30_80","HoverE_QCD_30_80",300,0.,0.15);
  HoverE_QCD_30_80->GetXaxis()->SetTitle("H/E");
  HoverE_QCD_30_80->GetYaxis()->SetTitle("a.u.");
  HoverE_QCD_30_80->SetLineColor(4);
  HoverE_QCD_30_80->SetLineWidth(2);

  Wenu     ->Draw("electrons_hOverE[0] >> HoverE_Wenu","electrons_isEB == 1","goff");
  QCD_30_80->Draw("electrons_hOverE[0] >> HoverE_QCD_30_80","electrons_isEB == 1","goff");

  TCanvas HoverE;
  HoverE_Wenu->DrawNormalized();
  HoverE_QCD_30_80->DrawNormalized("sames");
  HoverE.SetLogy();
  HoverE.Print("HoverE.pdf");


  //trackIso
  TH1F* trackIso_Wenu      = new TH1F("trackIso_Wenu","trackIso_Wenu",500,0.,10.);
  trackIso_Wenu->GetXaxis()->SetTitle("track Iso");
  trackIso_Wenu->GetYaxis()->SetTitle("a.u.");
  trackIso_Wenu->SetLineColor(2);
  trackIso_Wenu->SetLineWidth(2);

  TH1F* trackIso_QCD_30_80 = new TH1F("trackIso_QCD_30_80","trackIso_QCD_30_80",500,0.,10.);
  trackIso_QCD_30_80->GetXaxis()->SetTitle("track Iso");
  trackIso_QCD_30_80->GetYaxis()->SetTitle("a.u.");
  trackIso_QCD_30_80->SetLineColor(4);
  trackIso_QCD_30_80->SetLineWidth(2);

  Wenu     ->Draw("electrons_tkIsoR03[0] >> trackIso_Wenu","electrons_isEB == 1","goff");
  QCD_30_80->Draw("electrons_tkIsoR03[0] >> trackIso_QCD_30_80","electrons_isEB == 1","goff");

  TCanvas trackIso;
  trackIso_Wenu->DrawNormalized();
  trackIso_QCD_30_80->DrawNormalized("sames");
  trackIso.SetLogy();
  trackIso.Print("trackIso.pdf");

//   //ecalIso
//   TH1F* ecalIso_Wenu      = new TH1F("ecalIso_Wenu","ecalIso_Wenu",500,0.,10.);
//   ecalIso_Wenu->GetXaxis()->SetTitle("ECAL Iso");
//   ecalIso_Wenu->GetYaxis()->SetTitle("a.u.");
//   ecalIso_Wenu->SetLineColor(2);
//   ecalIso_Wenu->SetLineWidth(2);

//   TH1F* ecalIso_QCD_30_80 = new TH1F("ecalIso_QCD_30_80","ecalIso_QCD_30_80",500,0.,10.);
//   ecalIso_QCD_30_80->GetXaxis()->SetTitle("ecal Iso");
//   ecalIso_QCD_30_80->GetYaxis()->SetTitle("a.u.");
//   ecalIso_QCD_30_80->SetLineColor(4);
//   ecalIso_QCD_30_80->SetLineWidth(2);

//   Wenu     ->Draw("electrons_emIsoR03[0] >> ecalIso_Wenu","electrons_isEB == 1","goff");
//   QCD_30_80->Draw("electrons_emIsoR03[0] >> ecalIso_QCD_30_80","electrons_isEB == 1","goff");

//   TCanvas ecalIso;
//   ecalIso_Wenu->DrawNormalized();
//   ecalIso_QCD_30_80->DrawNormalized("sames");
//   ecalIso.Print("ecalIso.pdf");

}

