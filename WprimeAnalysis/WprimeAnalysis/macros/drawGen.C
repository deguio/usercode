//simple macro to draw the Wprime main variables at gen level
{
  //Get File
  TFile* Wprime_1500_f = new TFile("/gwteraz/users/deguio/MiBiCommonPAT/Summer11/MC_09112011/WprimeToENu_M-1500_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM/MiBiCommonNT.root","READ");
  TFile* Wprime_2000_f = new TFile("/gwteraz/users/deguio/MiBiCommonPAT/Summer11/MC_09112011/WprimeToENu_M-2000_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM/MiBiCommonNT.root","READ");
  TFile* Wprime_2500_f = new TFile("/gwteraz/users/deguio/MiBiCommonPAT/Summer11/MC_09112011/WprimeToENu_M-2500_7TeV-pythia6_Summer11-PU_S4_START42_V11-v1_AODSIM/MiBiCommonNT.root","READ");

  //Get Tree
  TTree* Wprime_1500  = (TTree*)Wprime_1500_f  ->Get("MiBiCommonNTOneElectron/SimpleNtuple");
  TTree* Wprime_2000  = (TTree*)Wprime_2000_f  ->Get("MiBiCommonNTOneElectron/SimpleNtuple");
  TTree* Wprime_2500  = (TTree*)Wprime_2500_f  ->Get("MiBiCommonNTOneElectron/SimpleNtuple");

  //TLegend
  TLegend legend(0.68, 0.78, 0.99, 0.99);
  legend.SetFillColor(kWhite);
  TLegend legendenu(0.68, 0.78, 0.99, 0.99);
  legendenu.SetFillColor(kWhite);

  //pt
  //ele
  TH1F* pt_1500      = new TH1F("pt_1500","pt_1500",80,0.,2000);
  pt_1500->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_1500->GetYaxis()->SetTitle("a.u.");
  pt_1500->SetLineColor(2);
  pt_1500->SetLineWidth(2);
  
  TH1F* pt_2000      = new TH1F("pt_2000","pt_2000",80,0.,2000);
  pt_2000->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_2000->GetYaxis()->SetTitle("a.u.");
  pt_2000->SetLineColor(3);
  pt_2000->SetLineWidth(2);
  
  TH1F* pt_2500      = new TH1F("pt_2500","pt_2500",80,0.,2000);
  pt_2500->GetXaxis()->SetTitle("p_{T} [GeV]");
  pt_2500->GetYaxis()->SetTitle("a.u.");
  pt_2500->SetLineColor(4);
  pt_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF1_fromV[0]->pt() >> pt_1500","","goff");
  Wprime_2000     ->Draw("mcF1_fromV[0]->pt() >> pt_2000","","goff");
  Wprime_2500     ->Draw("mcF1_fromV[0]->pt() >> pt_2500","","goff");

  //nu
  TH1F* ptnu_1500      = new TH1F("ptnu_1500","ptnu_1500",80,0.,2000);
  ptnu_1500->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptnu_1500->GetYaxis()->SetTitle("a.u.");
  ptnu_1500->SetLineColor(95);
  ptnu_1500->SetLineWidth(2);
  
  TH1F* ptnu_2000      = new TH1F("ptnu_2000","ptnu_2000",80,0.,2000);
  ptnu_2000->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptnu_2000->GetYaxis()->SetTitle("a.u.");
  ptnu_2000->SetLineColor(86);
  ptnu_2000->SetLineWidth(2);
  
  TH1F* ptnu_2500      = new TH1F("ptnu_2500","ptnu_2500",80,0.,2000);
  ptnu_2500->GetXaxis()->SetTitle("p_{T} [GeV]");
  ptnu_2500->GetYaxis()->SetTitle("a.u.");
  ptnu_2500->SetLineColor(7);
  ptnu_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF2_fromV[0]->pt() >> ptnu_1500","","goff");
  Wprime_2000     ->Draw("mcF2_fromV[0]->pt() >> ptnu_2000","","goff");
  Wprime_2500     ->Draw("mcF2_fromV[0]->pt() >> ptnu_2500","","goff");
  
  legend.AddEntry(pt_1500,"Wprime 1500","F");
  legend.AddEntry(pt_2000,"Wprime 2000","F");
  legend.AddEntry(pt_2500,"Wprime 2500","F");
  legendenu.AddEntry(pt_1500,"Wprime 1500: ele","F");
  legendenu.AddEntry(pt_2000,"Wprime 2000: ele","F");
  legendenu.AddEntry(pt_2500,"Wprime 2500: ele","F");
  legendenu.AddEntry(ptnu_1500,"Wprime 1500: nu","F");
  legendenu.AddEntry(ptnu_2000,"Wprime 2000: nu","F");
  legendenu.AddEntry(ptnu_2500,"Wprime 2500: nu","F");

  TCanvas pt;
  pt_1500->DrawNormalized();
  pt_2000->DrawNormalized("sames");
  pt_2500->DrawNormalized("sames");
  ptnu_1500->DrawNormalized("sames");
  ptnu_2000->DrawNormalized("sames");
  ptnu_2500->DrawNormalized("sames");
  legendenu.Draw("sames");
  pt.Print("pt.pdf");

  //eta
  //e
  TH1F* eta_1500      = new TH1F("eta_1500","eta_1500",100,-3.,3.);
  eta_1500->GetXaxis()->SetTitle("#eta");
  eta_1500->GetYaxis()->SetTitle("a.u.");
  eta_1500->SetLineColor(2);
  eta_1500->SetLineWidth(2);
  
  TH1F* eta_2000      = new TH1F("eta_2000","eta_2000",100,-3.,3.);
  eta_2000->GetXaxis()->SetTitle("#eta");
  eta_2000->GetYaxis()->SetTitle("a.u.");
  eta_2000->SetLineColor(3);
  eta_2000->SetLineWidth(2);
  
  TH1F* eta_2500      = new TH1F("eta_2500","eta_2500",100,-3.,3.);
  eta_2500->GetXaxis()->SetTitle("#eta");
  eta_2500->GetYaxis()->SetTitle("a.u.");
  eta_2500->SetLineColor(4);
  eta_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF1_fromV[0]->eta() >> eta_1500","","goff");
  Wprime_2000     ->Draw("mcF1_fromV[0]->eta() >> eta_2000","","goff");
  Wprime_2500     ->Draw("mcF1_fromV[0]->eta() >> eta_2500","","goff");
   
  //nu
  TH1F* etanu_1500      = new TH1F("etanu_1500","etanu_1500",100,-3.,3.);
  etanu_1500->GetXaxis()->SetTitle("#eta");
  etanu_1500->GetYaxis()->SetTitle("a.u.");
  etanu_1500->SetLineColor(95);
  etanu_1500->SetLineWidth(2);
  
  TH1F* etanu_2000      = new TH1F("etanu_2000","etanu_2000",100,-3.,3.);
  etanu_2000->GetXaxis()->SetTitle("#eta");
  etanu_2000->GetYaxis()->SetTitle("a.u.");
  etanu_2000->SetLineColor(86);
  etanu_2000->SetLineWidth(2);
  
  TH1F* etanu_2500      = new TH1F("etanu_2500","etanu_2500",100,-3.,3.);
  etanu_2500->GetXaxis()->SetTitle("#eta");
  etanu_2500->GetYaxis()->SetTitle("a.u.");
  etanu_2500->SetLineColor(7);
  etanu_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF2_fromV[0]->eta() >> etanu_1500","","goff");
  Wprime_2000     ->Draw("mcF2_fromV[0]->eta() >> etanu_2000","","goff");
  Wprime_2500     ->Draw("mcF2_fromV[0]->eta() >> etanu_2500","","goff");

  TCanvas eta;
  eta_1500->DrawNormalized();
  eta_2000->DrawNormalized("sames");
  eta_2500->DrawNormalized("sames");
  etanu_1500->DrawNormalized("sames");
  etanu_2000->DrawNormalized("sames");
  etanu_2500->DrawNormalized("sames");
  legendenu.Draw("sames");
  eta.Print("eta.pdf");

  //phi
  //e
  TH1F* phi_1500      = new TH1F("phi_1500","phi_1500",50,-4.,4.);
  phi_1500->GetXaxis()->SetTitle("#phi");
  phi_1500->GetYaxis()->SetTitle("a.u.");
  phi_1500->SetLineColor(95);
  phi_1500->SetLineWidth(2);
  
  TH1F* phi_2000      = new TH1F("phi_2000","phi_2000",50,-4.,4.);
  phi_2000->GetXaxis()->SetTitle("#phi");
  phi_2000->GetYaxis()->SetTitle("a.u.");
  phi_2000->SetLineColor(86);
  phi_2000->SetLineWidth(2);
  
  TH1F* phi_2500      = new TH1F("phi_2500","phi_2500",50,-4.,4.);
  phi_2500->GetXaxis()->SetTitle("#phi");
  phi_2500->GetYaxis()->SetTitle("a.u.");
  phi_2500->SetLineColor(7);
  phi_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF1_fromV[0]->phi() >> phi_1500","","goff");
  Wprime_2000     ->Draw("mcF1_fromV[0]->phi() >> phi_2000","","goff");
  Wprime_2500     ->Draw("mcF1_fromV[0]->phi() >> phi_2500","","goff");
  
  //nu
  TH1F* phinu_1500      = new TH1F("phinu_1500","phinu_1500",50,-4.,4.);
  phinu_1500->GetXaxis()->SetTitle("#phi");
  phinu_1500->GetYaxis()->SetTitle("a.u.");
  phinu_1500->SetLineColor(2);
  phinu_1500->SetLineWidth(2);
  
  TH1F* phinu_2000      = new TH1F("phinu_2000","phinu_2000",50,-4.,4.);
  phinu_2000->GetXaxis()->SetTitle("#phi");
  phinu_2000->GetYaxis()->SetTitle("a.u.");
  phinu_2000->SetLineColor(3);
  phinu_2000->SetLineWidth(2);
  
  TH1F* phinu_2500      = new TH1F("phinu_2500","phinu_2500",50,-4.,4.);
  phinu_2500->GetXaxis()->SetTitle("#phi");
  phinu_2500->GetYaxis()->SetTitle("a.u.");
  phinu_2500->SetLineColor(4);
  phinu_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mcF2_fromV[0]->phi() >> phinu_1500","","goff");
  Wprime_2000     ->Draw("mcF2_fromV[0]->phi() >> phinu_2000","","goff");
  Wprime_2500     ->Draw("mcF2_fromV[0]->phi() >> phinu_2500","","goff");


  TCanvas phi;
  phi_1500->DrawNormalized();
  phi_2000->DrawNormalized("sames");
  phi_2500->DrawNormalized("sames");
  phinu_1500->DrawNormalized("sames");
  phinu_2000->DrawNormalized("sames");
  phinu_2500->DrawNormalized("sames");
  legendenu.Draw("sames");
  phi.Print("phi.pdf");

  //Weta
  TH1F* Weta_1500      = new TH1F("Weta_1500","Weta_1500",100,-3.,3.);
  Weta_1500->GetXaxis()->SetTitle("#eta");
  Weta_1500->GetYaxis()->SetTitle("a.u.");
  Weta_1500->SetLineColor(2);
  Weta_1500->SetLineWidth(2);
  
  TH1F* Weta_2000      = new TH1F("Weta_2000","Weta_2000",100,-3.,3.);
  Weta_2000->GetXaxis()->SetTitle("#eta");
  Weta_2000->GetYaxis()->SetTitle("a.u.");
  Weta_2000->SetLineColor(3);
  Weta_2000->SetLineWidth(2);
  
  TH1F* Weta_2500      = new TH1F("Weta_2500","Weta_2500",100,-3.,3.);
  Weta_2500->GetXaxis()->SetTitle("#eta");
  Weta_2500->GetYaxis()->SetTitle("a.u.");
  Weta_2500->SetLineColor(4);
  Weta_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mc_V[0]->eta() >> Weta_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->eta() >> Weta_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->eta() >> Weta_2500","","goff");
  
  TCanvas Weta;
  Weta_1500->DrawNormalized();
  Weta_2000->DrawNormalized("sames");
  Weta_2500->DrawNormalized("sames");
  legend.Draw("sames");
  Weta.Print("Weta.pdf");

  //Wphi
  TH1F* Wphi_1500      = new TH1F("Wphi_1500","Wphi_1500",50,-4.,4.);
  Wphi_1500->GetXaxis()->SetTitle("#phi");
  Wphi_1500->GetYaxis()->SetTitle("a.u.");
  Wphi_1500->SetLineColor(2);
  Wphi_1500->SetLineWidth(2);
  
  TH1F* Wphi_2000      = new TH1F("Wphi_2000","Wphi_2000",50,-4.,4.);
  Wphi_2000->GetXaxis()->SetTitle("#phi");
  Wphi_2000->GetYaxis()->SetTitle("a.u.");
  Wphi_2000->SetLineColor(3);
  Wphi_2000->SetLineWidth(2);
  
  TH1F* Wphi_2500      = new TH1F("Wphi_2500","Wphi_2500",50,-4.,4.);
  Wphi_2500->GetXaxis()->SetTitle("#phi");
  Wphi_2500->GetYaxis()->SetTitle("a.u.");
  Wphi_2500->SetLineColor(4);
  Wphi_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mc_V[0]->phi() >> Wphi_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->phi() >> Wphi_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->phi() >> Wphi_2500","","goff");
  
  TCanvas Wphi;
  Wphi_1500->DrawNormalized();
  Wphi_2000->DrawNormalized("sames");
  Wphi_2500->DrawNormalized("sames");
  legend.Draw("sames");
  Wphi.Print("Wphi.pdf");

  //WPT
  TH1F* WPT_1500      = new TH1F("WPT_1500","WPT_1500",100,0.,300.);
  WPT_1500->GetXaxis()->SetTitle("p_{T} [GeV]");
  WPT_1500->GetYaxis()->SetTitle("a.u.");
  WPT_1500->SetLineColor(2);
  WPT_1500->SetLineWidth(2);
  
  TH1F* WPT_2000      = new TH1F("WPT_2000","WPT_2000",100,0.,300.);
  WPT_2000->GetXaxis()->SetTitle("p_{T} [GeV]");
  WPT_2000->GetYaxis()->SetTitle("a.u.");
  WPT_2000->SetLineColor(3);
  WPT_2000->SetLineWidth(2);
  
  TH1F* WPT_2500      = new TH1F("WPT_2500","WPT_2500",100,0.,300.);
  WPT_2500->GetXaxis()->SetTitle("p_{T} [GeV]");
  WPT_2500->GetYaxis()->SetTitle("a.u.");
  WPT_2500->SetLineColor(4);
  WPT_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mc_V[0]->Pt() >> WPT_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->Pt() >> WPT_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->Pt() >> WPT_2500","","goff");
  
  TCanvas WPT;
  WPT.SetLogy();
  WPT_2500->DrawNormalized();
  WPT_2000->DrawNormalized("sames");
  WPT_1500->DrawNormalized("sames");
  legend.Draw("sames");
  WPT.Print("WPT.pdf");

  //WPL
  TH1F* WPL_1500      = new TH1F("WPL_1500","WPL_1500",100,-4000.,4000.);
  WPL_1500->GetXaxis()->SetTitle("p_{L} [GeV]");
  WPL_1500->GetYaxis()->SetTitle("a.u.");
  WPL_1500->SetLineColor(2);
  WPL_1500->SetLineWidth(2);
  
  TH1F* WPL_2000      = new TH1F("WPL_2000","WPL_2000",100,-4000.,4000.);
  WPL_2000->GetXaxis()->SetTitle("p_{L} [GeV]");
  WPL_2000->GetYaxis()->SetTitle("a.u.");
  WPL_2000->SetLineColor(3);
  WPL_2000->SetLineWidth(2);
  
  TH1F* WPL_2500      = new TH1F("WPL_2500","WPL_2500",100,-4000.,4000.);
  WPL_2500->GetXaxis()->SetTitle("p_{L} [GeV]");
  WPL_2500->GetYaxis()->SetTitle("a.u.");
  WPL_2500->SetLineColor(4);
  WPL_2500->SetLineWidth(2);
  
  Wprime_1500     ->Draw("mc_V[0]->pz() >> WPL_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->pz() >> WPL_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->pz() >> WPL_2500","","goff");
  
  TCanvas WPL;
  WPL_2500->DrawNormalized();
  WPL_2000->DrawNormalized("sames");
  WPL_1500->DrawNormalized("sames");
  legend.Draw("sames");
  WPL.Print("WPL.pdf");
  
  //mt
  TH1F* mt_1500      = new TH1F("mt_1500","mt_1500",200,0.,5000);
  mt_1500->GetXaxis()->SetTitle("m_{T} [GeV]");
  mt_1500->GetYaxis()->SetTitle("a.u.");
  mt_1500->SetLineColor(2);
  mt_1500->SetLineWidth(2);
  
  TH1F* mt_2000      = new TH1F("mt_2000","mt_2000",200,0.,5000);
  mt_2000->GetXaxis()->SetTitle("m_{T} [GeV]");
  mt_2000->GetYaxis()->SetTitle("a.u.");
  mt_2000->SetLineColor(3);
  mt_2000->SetLineWidth(2);
  
  TH1F* mt_2500      = new TH1F("mt_2500","mt_2500",200,0.,5000);
  mt_2500->GetXaxis()->SetTitle("m_{T} [GeV]");
  mt_2500->GetYaxis()->SetTitle("a.u.");
  mt_2500->SetLineColor(4);
  mt_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mc_V[0]->Mt() >> mt_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->Mt() >> mt_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->Mt() >> mt_2500","","goff");
  
  TCanvas mt;
  mt_1500->DrawNormalized();
  mt_2000->DrawNormalized("sames");
  mt_2500->DrawNormalized("sames");
  legend.Draw("sames");
  mt.Print("mt.pdf");

  //minv
  TH1F* minv_1500      = new TH1F("minv_1500","minv_1500",200,0.,5000);
  minv_1500->GetXaxis()->SetTitle("m_{inv} [GeV]");
  minv_1500->GetYaxis()->SetTitle("a.u.");
  minv_1500->SetLineColor(2);
  minv_1500->SetLineWidth(2);
  
  TH1F* minv_2000      = new TH1F("minv_2000","minv_2000",200,0.,5000);
  minv_2000->GetXaxis()->SetTitle("m_{inv} [GeV]");
  minv_2000->GetYaxis()->SetTitle("a.u.");
  minv_2000->SetLineColor(3);
  minv_2000->SetLineWidth(2);
  
  TH1F* minv_2500      = new TH1F("minv_2500","minv_2500",200,0.,5000);
  minv_2500->GetXaxis()->SetTitle("m_{inv} [GeV]");
  minv_2500->GetYaxis()->SetTitle("a.u.");
  minv_2500->SetLineColor(4);
  minv_2500->SetLineWidth(2);

  Wprime_1500     ->Draw("mc_V[0]->M() >> minv_1500","","goff");
  Wprime_2000     ->Draw("mc_V[0]->M() >> minv_2000","","goff");
  Wprime_2500     ->Draw("mc_V[0]->M() >> minv_2500","","goff");
  
  TCanvas minv;
  minv_1500->DrawNormalized();
  minv_2000->DrawNormalized("sames");
  minv_2500->DrawNormalized("sames");
  legend.Draw("sames");
  minv.Print("minv.pdf");

}
