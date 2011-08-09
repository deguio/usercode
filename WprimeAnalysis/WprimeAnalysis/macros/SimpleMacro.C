{
  TH1F hmt("hmt","hmt",120,0,1200);
  TFile fin("/data2/ghezzi/CMSSW_4_2_3/src/merged_ntu_WprimeTreeAnalysis.root");
  TTree * tr = (TTree*) fin.Get("ntu_10");

  float mt;
  
  tr->SetBranchAddress("eleMet_mt", &mt);

   for(int ientry = 0; ientry < tr->GetEntries(); ++ientry)
    {
      if ( ientry%3 != 2 ) continue; 
      
      if( (ientry%100000 == 0) ) std::cout << "reading MC entry " << ientry << std::endl;
      tr->GetEntry(ientry);
      hmt.Fill(mt);

    }

   hmt.Draw();
   //TFile fout("myh.root","recreate");
   //hmt.Write();
   //return;

   // TF1* myf = new TF1("myf","[0]/(x*x + [1]*x + [2])^[3]",0,1500);
   //myf->SetParameters(2000000,-300,30000,10.5);

    TF1* myf = new TF1("myf","[0]/(x + [1])^[2]",0,1500);
    myf->SetParameters(1000000,50,15.5);

   for(int u=0; u < 20; u++){
   hmt.Fit(myf,"L","",170,500);
   cout<<myf->Integral(700,900)<<endl;
   }
}

