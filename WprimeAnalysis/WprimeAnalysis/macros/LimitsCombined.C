#include <iostream>
#include <stdlib.h>
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TArrow.h"
#include "TCanvas.h"
#include <stdio.h>
#include <TLatex.h>
using namespace std;


	
void LimitsCombined()
{
  // Number of massPoints ele channel
  const int nEle=18;
  // Number of different search windows
  const int numEle=37;
  
  //Luminosity
  TString SlumiEle="1091.400";
  double lumiEle=1132;
  double slumiEle=lumiEle*0.045;
  
  // Calculate Limits?
  double doEle=true;
  // Optimize the searchwindow?
  bool optEle=true;
  
  // Searchwindow lower threshold
  int bin[numEle];
  
  // Ele Arrays
  int dataEle[numEle];
  double totalbackgroundEle[numEle];
  double stotalbackgroundEle[numEle];
  double sigeffEle[nEle][numEle];
  double sigefferrEle[nEle][numEle];
  
  //read Ele		
  // ELECTRON Signal eff. + background + data info (format threshold + #dataevents + background + wprime sig effs
  ifstream in("numbersForLimit.txt");
  if (!in)
    cout << "No File in this folder!"<<endl;
  else 
    for (int i = 0; i < numEle; i++) 
      {
	in >> bin[i] >> dataEle[i] >> totalbackgroundEle[i] >> stotalbackgroundEle[i];
	for (int j = 0; j < nEle; j++)
	  {
	    in >> sigeffEle[j][i] >> sigefferrEle[j][i];	
	    //std::cout << sigeffEle[j][i] << std::endl;
	  }
      }
	
  in.close();
  
  // Defining search windows
  Double_t massEle[nEle] = {1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,3000,3500,4000};
  
  // Theoretical cross section and k-factors
  Double_t crossneu[18]={0.34608,0.22210,0.14402,0.09485,0.06333,0.04237,0.02851,0.01940,0.01346,0.00937,0.00661,0.00472,0.00340,0.00248,0.00143,0.00071,0.00030,0.00015};
  Double_t kfactorneu[18]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
  Double_t massneu[18]={1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2700,3000,3500,4000};
  Double_t pdfuncertainty[18]={0.06516,0.07349,0.07760,0.08471,0.09211,0.09887,0.1044,0.1074,0.1107,0.1030,0.1029,0.10340,0.09862,0.09850,0.09846,0.06221,0.04385,0.03000};
  
  
  
  Double_t xsneu[18];
  Double_t sigmaxsneu[18];
  for (int i=0; i<18; i++)
    {
      xsneu[i]=crossneu[i]*kfactorneu[i];
      sigmaxsneu[i]=xsneu[i]*pdfuncertainty[i];
    }

  // Arrays for saving limit values
  Double_t limit_Ele[nEle][numEle];
  Double_t low_Ele[nEle];
  Int_t low_i_Ele[nEle];
  
  //optimize ele
  for (int i = 0; i < nEle; i++)
    {
      low_Ele[i]=500.;
      low_i_Ele[i]=99.;
      for (int j = 0; j < numEle; j++)
	{
	  cout << endl << "Step i: " << i << " j " << j << endl <<endl;
	  limit_Ele[i][j]=roostats_cla(lumiEle, slumiEle, sigeffEle[i][j], sigefferrEle[i][j], totalbackgroundEle[j], stotalbackgroundEle[j],1);
	  if (limit_Ele[i][j] < low_Ele[i])
	    {	
	      low_Ele[i]=limit_Ele[i][j];
	      low_i_Ele[i]=j;
	      cout << limit_Ele[i][j] << " Lowest so far: " << low_Ele[i]   <<  endl;
	    }
	}
    }
  FILE *file2;
  file2 = fopen("limitwindowsEle.txt","w");
  for(int i = 0; i < nEle; i++)
    {
      fprintf(file2,"%d %d   \n",i,low_i_Ele[i]);
    }
  fclose(file2);
  
  //Some dummy value for scaling the y axis
  double crossw = 1.;
  
  double limitEle[nEle];
  double limitEleExpected[nEle];
  double limitEleExpected1S[2*nEle];
  double limitEleExpected2S[2*nEle];
  
  FILE *file3;
  // Number of iterations for the errorband calculation
  Int_t N=100;//fede (era 100)
  
  // do the actual calculation
  for (int i=0; i<nEle; i++)
    {
      cout << "ELE STEP " << i << endl;
      file3 = fopen("limitEle.txt","a");
      limitEle[i]=roostats_cl95(lumiEle, slumiEle, sigeffEle[i][low_i_Ele[i]], sigefferrEle[i][low_i_Ele[i]], totalbackgroundEle[low_i_Ele[i]], stotalbackgroundEle[low_i_Ele[i]], dataEle[low_i_Ele[i]], false, 1);
      
      std::cout << "------>>>>> lumi = " << lumiEle << " slumi = " << slumiEle << " sigEff = " << sigeffEle[i][low_i_Ele[i]] << " sigEffErr = " << sigefferrEle[i][low_i_Ele[i]]
		<< " totBkg = " << totalbackgroundEle[low_i_Ele[i]] << " totBkgErr = " << stotalbackgroundEle[low_i_Ele[i]] << " data = " << dataEle[low_i_Ele[i]] << std::endl;
      std::cout << "Limit = " << limitEle[i] << std::endl;

      LimitResult result=roostats_clm(lumiEle, slumiEle, sigeffEle[i][low_i_Ele[i]], sigefferrEle[i][low_i_Ele[i]], totalbackgroundEle[low_i_Ele[i]], stotalbackgroundEle[low_i_Ele[i]], N, 1);
      cout << i << " " << 2*nEle-i-1 << endl;
      limitEleExpected[i]=result.GetExpectedLimit();
      limitEleExpected1S[i]=result.GetOneSigmaLowRange();
      limitEleExpected1S[2*nEle-i-1]=result.GetOneSigmaHighRange();
      limitEleExpected2S[i]=result.GetTwoSigmaLowRange();
      limitEleExpected2S[2*nEle-i-1]=result.GetTwoSigmaHighRange();
      limitEleExpected[i]=roostats_cla(lumiEle, slumiEle, sigeffEle[i][low_i_Ele[i]], sigefferrEle[i][low_i_Ele[i]], totalbackgroundEle[low_i_Ele[i]], stotalbackgroundEle[low_i_Ele[i]],1);
      fprintf(file3,"%i & %i & %4.3f PM %4.3f & %4.3f PM %4.3f & %i & %4.3f & %4.3f kkk\n",massEle[i], bin[low_i_Ele[i]], sigeffEle[i][low_i_Ele[i]], sigefferrEle[i][low_i_Ele[i]],totalbackgroundEle[low_i_Ele[i]], stotalbackgroundEle[low_i_Ele[i]],dataEle[low_i_Ele[i]],limitEleExpected[i],limitEle[i]);
      fclose(file3);
    }
  
  file3 = fopen("limitExpectedBands.txt","a");
  fprintf(file3,"limitEleExpected1S[2*nEle]={");
  for (int i=0; i<2*nEle; i++)
    {
      if (i<(2*nEle-1)) fprintf(file3,"%g,",limitEleExpected1S[i]);
      else fprintf(file3,"%g};\n",limitEleExpected1S[i]);
    }
  fprintf(file3,"limitEleExpected2S[2*nEle]={");
  for (int i=0; i<2*nEle; i++){
    if (i<(2*nEle-1)) fprintf(file3,"%g,",limitEleExpected2S[i]);
    else fprintf(file3,"%g};\n",limitEleExpected2S[i]);
  }
  
  fclose(file3);  
  
  
  
  
  Double_t massbandEle[2*nEle];
  for (int i=0; i<nEle; i++)
    {
      massbandEle[i]=massEle[i];
      massbandEle[2*nEle-i-1]=massEle[i];
    }
  
  
  //~  Dummy for text in canvas
  //~
  
  TText *cmspre;
  TLatex intlumi;
  TLatex intlumie;
  TLatex ecm;
  TCanvas *c1= new TCanvas("c1","c1");
  c1->SetGrid();
  
  
  c1->SetLogy();
  
  TGraph* GraphEle=new TGraph(nEle,massEle,limitEle);
  TGraph* GraphEleExpected=new TGraph(nEle,massEle,limitEleExpected);
  TGraph* GraphEleExpected1Sigma = new TGraph(2*nEle,massbandEle,limitEleExpected1S);
  TGraph* GraphEleExpected2Sigma = new TGraph(2*nEle,massbandEle,limitEleExpected2S);
  GraphEle->SetLineWidth(2);
  GraphEleExpected->SetLineWidth(2);
  
  // Theoretical cross section
  TGraphErrors* xsvsmass= new TGraphErrors(18,massneu,xsneu,0,sigmaxsneu);
  TGraph* xsvsmass2=new TGraph(18,massneu,xsneu);
  c1->cd();
  // Style issues
  xsvsmass->SetLineWidth(2);
  xsvsmass->SetTitle("");
  xsvsmass->GetXaxis()->SetTitle("W' mass (GeV)");
  xsvsmass2->SetMarkerStyle(1);
  xsvsmass->SetMarkerStyle(1);
  xsvsmass->GetYaxis()->SetTitleSize(0.05);
  xsvsmass->GetYaxis()->SetTitle("#sigma ^{.} BR(W' #rightarrow e + #nu) (pb)");
  
//xsvsmass->SetMarkerStyle(3);
  
  
  xsvsmass->SetFillColor(kGreen-8);
  xsvsmass->SetLineWidth(2);
  xsvsmass->GetXaxis()->SetRangeUser(1200,2500);
  
  xsvsmass->Draw("APE3L");
  xsvsmass->Draw("E3LSAME");
  // xsvsmass->Draw("AL");
  xsvsmass->GetXaxis()->SetRangeUser(1200.,2500.);
  xsvsmass->GetXaxis()->SetRangeUser(0.001,1.);
  xsvsmass2->SetLineColor(kBlack);
  xsvsmass2->SetLineWidth(2);
  xsvsmass2->Draw("LSAME");
  //xsvsmass2->Draw("AL");
  
  cmspre = new TText(1200,0.025,"CMS");
  cmspre->SetTextSize(0.033);
  cmspre -> Draw();
  ecm.SetTextAlign(12);
  ecm.SetTextSize(0.03);
  ecm.DrawLatex(1200,0.02,"#sqrt{s} = 7 TeV");
  intlumi.SetTextAlign(12);
  intlumi.SetTextSize(0.03);
  intlumi.DrawLatex(1200,0.013,"#int L dt = "+SlumiEle+" pb^{-1}");
  xsvsmass->GetXaxis()->SetRangeUser(1200.,2500.);
  xsvsmass->GetYaxis()->SetRangeUser(0.001,1.);
  
  TLegend *leg = new TLegend(0.5,0.65,0.95,0.95);
  leg->SetShadowColor(0);
  
  
  // ELECTRON LIMIT IN ONE CANVAS
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->SetLogy();
  c2->SetGrid(1);
  xsvsmass->Draw("APE3L");
  GraphEleExpected2Sigma->Draw("FLSAME");
  GraphEleExpected1Sigma->Draw("FLSAME");
  xsvsmass->Draw("E3LSAME");
  xsvsmass2->Draw("LSAME");
  GraphEle->Draw("LSAME");
  GraphEleExpected->Draw("LSAME");
  GraphEleExpected->SetLineColor(kRed);
  GraphEleExpected1Sigma->SetFillColor(kGreen);
  GraphEleExpected2Sigma->SetFillColor(kYellow);
  GraphEle->SetLineColor(kRed);
  GraphEleExpected->SetLineStyle(9);
  
  cmspre = new TText(1200,0.025,"CMS");
  cmspre->SetTextSize(0.033);
  cmspre -> Draw();
  ecm.SetTextAlign(12);
  ecm.SetTextSize(0.03);
  ecm.DrawLatex(1200,0.02,"#sqrt{s} = 7 TeV");
  intlumi.SetTextAlign(12);
  intlumi.SetTextSize(0.03);
  
  
  intlumi.DrawLatex(1200,0.013,"#int L dt = "+SlumiEle+" pb^{-1}");
  
  xsvsmass->GetXaxis()->SetRangeUser(1200.,2500.);
  xsvsmass->GetYaxis()->SetRangeUser(0.001,1.);
  TLegend *leg2 = new TLegend(0.5,0.65,0.95,0.95);
  leg2->SetShadowColor(0);
  
  leg2->AddEntry(GraphEle,"95% Observed Limit Electron","l");
  leg2->AddEntry(GraphEleExpected,"95% Expected Limit Electron","l");
  leg2->AddEntry(GraphEleExpected1Sigma,"Expected #pm1#sigma Electron","f");
  
  leg2->AddEntry(GraphEleExpected2Sigma,"Expected #pm2#sigma Electron","f");
  leg2->SetFillColor(kWhite);
  
  //   leg->AddEntry(l2,"Expected Limit at 95% CL","l");
  leg2->AddEntry(xsvsmass,"Theoretical Cross Section","fl");
  leg2->Draw();
  c2->Print("ExpectedEle.pdf");
  //~ 
  //~ 
  // EVERYTHING TOGETHER
  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetLogy();
  c3->SetGrid(1);
  
  xsvsmass->Draw("APE3L");
  xsvsmass2->Draw("LSAME");
  GraphEle->Draw("LSAME");
  GraphEleExpected->Draw("LSAME");
  GraphEleExpected->SetLineColor(kRed);
  GraphEle->SetLineColor(kRed);
  GraphEleExpected->SetLineStyle(9);
  cmspre = new TText(1200,0.135,"CMS");
  cmspre->SetTextSize(0.033);
  cmspre -> Draw();
  ecm.SetTextAlign(12);
  ecm.SetTextSize(0.03);
  ecm.DrawLatex(1200,0.108,"#sqrt{s} = 7 TeV");
  //~ intlumi.SetTextAlign(12);
  //~ intlumi.SetTextSize(0.03);
  //~ intlumi.DrawLatex(1200,0.013,"#int L dt = 153 pb^{-1}");
  
  xsvsmass->GetXaxis()->SetRangeUser(1200.,2500.);
  xsvsmass->GetYaxis()->SetRangeUser(0.001,1.);
  TLegend *leg3 = new TLegend(0.5,0.65,0.95,0.95);
  leg3->SetShadowColor(0);
  
  leg3->AddEntry(GraphEle,"95% Observed Limit Electron "+SlumiEle+" pb^{-1}","l");
  leg3->AddEntry(GraphEleExpected,"95% Expected Limit Electron "+SlumiEle+" pb^{-1}","l");
  leg3->SetFillColor(kWhite);
  
  //   leg->AddEntry(l2,"Expected Limit at 95% CL","l");
  leg3->AddEntry(xsvsmass,"Theoretical Cross Section","fl");
  leg3->Draw();
  c3->Print("Both.pdf");
  
}


void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("emr"); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.13);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

   //tdrStyle->SetBarOffset(Float_t baroff = 0.5);
   //tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
   //tdrStyle->SetPaintTextFormat(const char* format = "g");
   tdrStyle->SetPalette(1);
   //tdrStyle->SetTimeOffset(Double_t toffset);
   //tdrStyle->SetHistMinimumZero(kTRUE);






   const Int_t NRGBs = 5;
   const Int_t NCont = 255;

   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   tdrStyle->SetNumberContours(NCont);


   //TLatex *lab = new TLatex(0.70,0.85, "CMS 2008");
   //lab->SetNDC();
   //lab->SetTextFont(42);
   //lab->SetTextSize(0.05);
   //lab->Draw("same");
   



  gROOT -> ForceStyle();

  tdrStyle->cd();
}

Double_t Poisson(Double_t Mu, Int_t n)
// Calculate the Poission prob. of seeing n events given an expectation of Mu.
{
	if (Mu <= 0) return n > 0 ? 0. : 1.;
	
	Double_t logP;
	//
	// Tabulate values of -\sum log(i+2) up to n=1000 for faster calculation
	//
	
	//   Double_t sum = 0.;
	//   for (Int_t j = 0; j < 999; j++)
	//   {
	//		sum -= log(j+2.);
	//		printf("%10f, ",sum);
	//		if (((j+1)/10)*10 == j + 1) printf("\n");
	//   }
	
	Double_t logTable[999] = {-0.693147,  -1.791759,  -3.178054,  -4.787492,  -6.579251,  -8.525161, -10.604603, -12.801827, -15.104413, -17.502308, 
		-19.987214, -22.552164, -25.191221, -27.899271, -30.671860, -33.505073, -36.395445, -39.339884, -42.335616, -45.380139, 
		-48.471181, -51.606676, -54.784729, -58.003605, -61.261702, -64.557539, -67.889743, -71.257039, -74.658236, -78.092224, 
		-81.557959, -85.054467, -88.580828, -92.136176, -95.719695, -99.330612, -102.968199, -106.631760, -110.320640, -114.034212, 
		-117.771881, -121.533082, -125.317271, -129.123934, -132.952575, -136.802723, -140.673924, -144.565744, -148.477767, -152.409593, 
		-156.360836, -160.331128, -164.320112, -168.327445, -172.352797, -176.395848, -180.456291, -184.533829, -188.628173, -192.739047, 
		-196.866182, -201.009316, -205.168199, -209.342587, -213.532241, -217.736934, -221.956442, -226.190548, -230.439044, -234.701723, 
		-238.978390, -243.268849, -247.572914, -251.890402, -256.221136, -260.564941, -264.921650, -269.291098, -273.673124, -278.067573, 
		-282.474293, -286.893133, -291.323950, -295.766601, -300.220949, -304.686857, -309.164194, -313.652830, -318.152640, -322.663499, 
		-327.185288, -331.717887, -336.261182, -340.815059, -345.379407, -349.954118, -354.539086, -359.134205, -363.739376, -368.354496, 
		-372.979469, -377.614198, -382.258589, -386.912549, -391.575988, -396.248817, -400.930948, -405.622296, -410.322777, -415.032307, 
		-419.750806, -424.478193, -429.214392, -433.959324, -438.712914, -443.475088, -448.245773, -453.024896, -457.812388, -462.608179, 
		-467.412200, -472.224384, -477.044665, -481.872979, -486.709261, -491.553448, -496.405478, -501.265291, -506.132825, -511.008023, 
		-515.890825, -520.781174, -525.679014, -530.584288, -535.496943, -540.416924, -545.344178, -550.278652, -555.220294, -560.169054, 
		-565.124881, -570.087726, -575.057539, -580.034273, -585.017879, -590.008312, -595.005524, -600.009471, -605.020106, -610.037386, 
		-615.061266, -620.091704, -625.128657, -630.172082, -635.221938, -640.278184, -645.340779, -650.409683, -655.484857, -660.566261, 
		-665.653857, -670.747608, -675.847474, -680.953420, -686.065407, -691.183401, -696.307365, -701.437264, -706.573062, -711.714726, 
		-716.862220, -722.015512, -727.174567, -732.339353, -737.509837, -742.685987, -747.867770, -753.055156, -758.248113, -763.446610, 
		-768.650617, -773.860103, -779.075039, -784.295395, -789.521141, -794.752250, -799.988692, -805.230439, -810.477463, -815.729736, 
		-820.987232, -826.249922, -831.517780, -836.790780, -842.068894, -847.352098, -852.640365, -857.933670, -863.231987, -868.535292, 
		-873.843560, -879.156766, -884.474886, -889.797896, -895.125772, -900.458491, -905.796029, -911.138363, -916.485471, -921.837329, 
		-927.193915, -932.555207, -937.921183, -943.291821, -948.667100, -954.046997, -959.431492, -964.820564, -970.214191, -975.612354, 
		-981.015031, -986.422203, -991.833849, -997.249950, -1002.670485, -1008.095435, -1013.524780, -1018.958502, -1024.396582, -1029.838999, 
		-1035.285737, -1040.736775, -1046.192096, -1051.651682, -1057.115514, -1062.583574, -1068.055844, -1073.532308, -1079.012947, -1084.497744, 
		-1089.986681, -1095.479743, -1100.976911, -1106.478169, -1111.983501, -1117.492889, -1123.006318, -1128.523771, -1134.045232, -1139.570685, 
		-1145.100114, -1150.633503, -1156.170838, -1161.712101, -1167.257279, -1172.806355, -1178.359314, -1183.916142, -1189.476824, -1195.041344, 
		-1200.609689, -1206.181843, -1211.757792, -1217.337522, -1222.921018, -1228.508267, -1234.099254, -1239.693965, -1245.292387, -1250.894506, 
		-1256.500308, -1262.109780, -1267.722908, -1273.339679, -1278.960080, -1284.584097, -1290.211718, -1295.842930, -1301.477720, -1307.116075, 
		-1312.757982, -1318.403428, -1324.052403, -1329.704892, -1335.360884, -1341.020366, -1346.683326, -1352.349753, -1358.019634, -1363.692957, 
		-1369.369711, -1375.049884, -1380.733463, -1386.420439, -1392.110798, -1397.804530, -1403.501624, -1409.202067, -1414.905850, -1420.612960, 
		-1426.323387, -1432.037120, -1437.754148, -1443.474460, -1449.198045, -1454.924892, -1460.654992, -1466.388333, -1472.124906, -1477.864699, 
		-1483.607702, -1489.353905, -1495.103298, -1500.855871, -1506.611613, -1512.370515, -1518.132566, -1523.897757, -1529.666078, -1535.437519, 
		-1541.212071, -1546.989723, -1552.770467, -1558.554292, -1564.341189, -1570.131149, -1575.924163, -1581.720221, -1587.519313, -1593.321432, 
		-1599.126567, -1604.934709, -1610.745850, -1616.559981, -1622.377092, -1628.197175, -1634.020221, -1639.846221, -1645.675166, -1651.507049, 
		-1657.341860, -1663.179590, -1669.020232, -1674.863776, -1680.710215, -1686.559540, -1692.411742, -1698.266814, -1704.124747, -1709.985533, 
		-1715.849165, -1721.715633, -1727.584930, -1733.457047, -1739.331978, -1745.209714, -1751.090247, -1756.973569, -1762.859673, -1768.748551, 
		-1774.640196, -1780.534598, -1786.431752, -1792.331650, -1798.234283, -1804.139645, -1810.047728, -1815.958524, -1821.872027, -1827.788229, 
		-1833.707123, -1839.628702, -1845.552957, -1851.479884, -1857.409473, -1863.341718, -1869.276612, -1875.214148, -1881.154319, -1887.097119, 
		-1893.042539, -1898.990574, -1904.941217, -1910.894460, -1916.850298, -1922.808722, -1928.769728, -1934.733307, -1940.699454, -1946.668161, 
		-1952.639423, -1958.613233, -1964.589584, -1970.568470, -1976.549884, -1982.533820, -1988.520272, -1994.509233, -2000.500698, -2006.494659, 
		-2012.491111, -2018.490048, -2024.491463, -2030.495350, -2036.501703, -2042.510516, -2048.521784, -2054.535499, -2060.551656, -2066.570249, 
		-2072.591272, -2078.614720, -2084.640586, -2090.668864, -2096.699550, -2102.732636, -2108.768117, -2114.805988, -2120.846243, -2126.888876, 
		-2132.933881, -2138.981253, -2145.030987, -2151.083076, -2157.137515, -2163.194299, -2169.253423, -2175.314879, -2181.378665, -2187.444773, 
		-2193.513198, -2199.583936, -2205.656981, -2211.732327, -2217.809969, -2223.889902, -2229.972121, -2236.056620, -2242.143395, -2248.232440, 
		-2254.323750, -2260.417320, -2266.513144, -2272.611219, -2278.711537, -2284.814096, -2290.918889, -2297.025912, -2303.135160, -2309.246627, 
		-2315.360309, -2321.476201, -2327.594299, -2333.714596, -2339.837089, -2345.961772, -2352.088641, -2358.217692, -2364.348918, -2370.482316, 
		-2376.617881, -2382.755608, -2388.895493, -2395.037530, -2401.181716, -2407.328045, -2413.476513, -2419.627116, -2425.779849, -2431.934707, 
		-2438.091686, -2444.250781, -2450.411988, -2456.575303, -2462.740721, -2468.908238, -2475.077848, -2481.249549, -2487.423335, -2493.599202, 
		-2499.777146, -2505.957163, -2512.139248, -2518.323397, -2524.509606, -2530.697870, -2536.888185, -2543.080548, -2549.274953, -2555.471397, 
		-2561.669876, -2567.870385, -2574.072920, -2580.277478, -2586.484054, -2592.692644, -2598.903244, -2605.115850, -2611.330458, -2617.547065, 
		-2623.765665, -2629.986255, -2636.208831, -2642.433390, -2648.659926, -2654.888437, -2661.118919, -2667.351367, -2673.585777, -2679.822147, 
		-2686.060472, -2692.300747, -2698.542971, -2704.787138, -2711.033244, -2717.281287, -2723.531263, -2729.783166, -2736.036995, -2742.292745, 
		-2748.550413, -2754.809994, -2761.071486, -2767.334884, -2773.600185, -2779.867386, -2786.136482, -2792.407471, -2798.680348, -2804.955110, 
		-2811.231753, -2817.510275, -2823.790671, -2830.072937, -2836.357071, -2842.643070, -2848.930928, -2855.220644, -2861.512213, -2867.805632, 
		-2874.100898, -2880.398007, -2886.696957, -2892.997742, -2899.300361, -2905.604810, -2911.911085, -2918.219184, -2924.529102, -2930.840837, 
		-2937.154385, -2943.469743, -2949.786908, -2956.105876, -2962.426644, -2968.749209, -2975.073568, -2981.399718, -2987.727655, -2994.057376, 
		-3000.388877, -3006.722157, -3013.057211, -3019.394037, -3025.732631, -3032.072990, -3038.415112, -3044.758992, -3051.104629, -3057.452018, 
		-3063.801157, -3070.152043, -3076.504672, -3082.859042, -3089.215150, -3095.572992, -3101.932566, -3108.293868, -3114.656896, -3121.021647, 
		-3127.388118, -3133.756305, -3140.126206, -3146.497818, -3152.871137, -3159.246162, -3165.622889, -3172.001315, -3178.381438, -3184.763254, 
		-3191.146760, -3197.531955, -3203.918834, -3210.307396, -3216.697636, -3223.089553, -3229.483144, -3235.878406, -3242.275335, -3248.673930, 
		-3255.074188, -3261.476105, -3267.879679, -3274.284908, -3280.691788, -3287.100316, -3293.510491, -3299.922310, -3306.335768, -3312.750865, 
		-3319.167598, -3325.585963, -3332.005958, -3338.427580, -3344.850827, -3351.275696, -3357.702184, -3364.130290, -3370.560009, -3376.991340, 
		-3383.424280, -3389.858827, -3396.294977, -3402.732729, -3409.172079, -3415.613026, -3422.055566, -3428.499697, -3434.945417, -3441.392723, 
		-3447.841612, -3454.292083, -3460.744132, -3467.197757, -3473.652955, -3480.109725, -3486.568063, -3493.027968, -3499.489436, -3505.952465, 
		-3512.417053, -3518.883198, -3525.350897, -3531.820147, -3538.290947, -3544.763293, -3551.237184, -3557.712616, -3564.189589, -3570.668098, 
		-3577.148143, -3583.629720, -3590.112827, -3596.597463, -3603.083624, -3609.571308, -3616.060512, -3622.551236, -3629.043476, -3635.537230, 
		-3642.032495, -3648.529270, -3655.027552, -3661.527339, -3668.028629, -3674.531419, -3681.035707, -3687.541491, -3694.048769, -3700.557538, 
		-3707.067797, -3713.579542, -3720.092772, -3726.607485, -3733.123678, -3739.641349, -3746.160496, -3752.681117, -3759.203210, -3765.726773, 
		-3772.251802, -3778.778297, -3785.306255, -3791.835674, -3798.366551, -3804.898886, -3811.432675, -3817.967916, -3824.504607, -3831.042747, 
		-3837.582333, -3844.123363, -3850.665835, -3857.209747, -3863.755097, -3870.301882, -3876.850101, -3883.399752, -3889.950832, -3896.503340, 
		-3903.057274, -3909.612630, -3916.169409, -3922.727607, -3929.287222, -3935.848253, -3942.410697, -3948.974552, -3955.539817, -3962.106490, 
		-3968.674567, -3975.244049, -3981.814932, -3988.387214, -3994.960895, -4001.535970, -4008.112440, -4014.690301, -4021.269553, -4027.850192, 
		-4034.432217, -4041.015626, -4047.600417, -4054.186589, -4060.774139, -4067.363066, -4073.953367, -4080.545040, -4087.138085, -4093.732498, 
		-4100.328279, -4106.925425, -4113.523934, -4120.123804, -4126.725034, -4133.327622, -4139.931566, -4146.536864, -4153.143514, -4159.751515, 
		-4166.360864, -4172.971560, -4179.583601, -4186.196985, -4192.811711, -4199.427776, -4206.045179, -4212.663918, -4219.283991, -4225.905397, 
		-4232.528133, -4239.152198, -4245.777591, -4252.404308, -4259.032350, -4265.661713, -4272.292396, -4278.924398, -4285.557717, -4292.192350, 
		-4298.828297, -4305.465555, -4312.104122, -4318.743998, -4325.385180, -4332.027667, -4338.671457, -4345.316548, -4351.962938, -4358.610627, 
		-4365.259611, -4371.909890, -4378.561462, -4385.214325, -4391.868478, -4398.523918, -4405.180645, -4411.838656, -4418.497950, -4425.158525, 
		-4431.820380, -4438.483512, -4445.147921, -4451.813605, -4458.480562, -4465.148790, -4471.818288, -4478.489054, -4485.161087, -4491.834385, 
		-4498.508947, -4505.184770, -4511.861853, -4518.540196, -4525.219795, -4531.900649, -4538.582758, -4545.266119, -4551.950731, -4558.636592, 
		-4565.323700, -4572.012055, -4578.701654, -4585.392497, -4592.084580, -4598.777904, -4605.472466, -4612.168265, -4618.865299, -4625.563567, 
		-4632.263068, -4638.963799, -4645.665759, -4652.368947, -4659.073361, -4665.779001, -4672.485863, -4679.193947, -4685.903251, -4692.613774, 
		-4699.325515, -4706.038471, -4712.752642, -4719.468025, -4726.184620, -4732.902424, -4739.621438, -4746.341658, -4753.063083, -4759.785713, 
		-4766.509546, -4773.234579, -4779.960813, -4786.688244, -4793.416873, -4800.146697, -4806.877715, -4813.609926, -4820.343328, -4827.077919, 
		-4833.813700, -4840.550666, -4847.288819, -4854.028156, -4860.768675, -4867.510376, -4874.253256, -4880.997315, -4887.742552, -4894.488964, 
		-4901.236550, -4907.985310, -4914.735241, -4921.486343, -4928.238613, -4934.992051, -4941.746655, -4948.502424, -4955.259356, -4962.017451, 
		-4968.776706, -4975.537121, -4982.298694, -4989.061423, -4995.825308, -5002.590347, -5009.356539, -5016.123882, -5022.892375, -5029.662017, 
		-5036.432806, -5043.204742, -5049.977822, -5056.752046, -5063.527412, -5070.303919, -5077.081566, -5083.860351, -5090.640273, -5097.421330, 
		-5104.203522, -5110.986848, -5117.771305, -5124.556892, -5131.343609, -5138.131454, -5144.920426, -5151.710523, -5158.501745, -5165.294089, 
		-5172.087555, -5178.882142, -5185.677848, -5192.474671, -5199.272612, -5206.071668, -5212.871838, -5219.673121, -5226.475515, -5233.279021, 
		-5240.083635, -5246.889358, -5253.696187, -5260.504122, -5267.313161, -5274.123304, -5280.934548, -5287.746893, -5294.560338, -5301.374881, 
		-5308.190521, -5315.007257, -5321.825087, -5328.644011, -5335.464028, -5342.285135, -5349.107333, -5355.930619, -5362.754992, -5369.580452, 
		-5376.406998, -5383.234627, -5390.063339, -5396.893133, -5403.724007, -5410.555960, -5417.388992, -5424.223101, -5431.058286, -5437.894545, 
		-5444.731878, -5451.570283, -5458.409759, -5465.250306, -5472.091921, -5478.934605, -5485.778355, -5492.623170, -5499.469050, -5506.315993, 
		-5513.163998, -5520.013065, -5526.863191, -5533.714376, -5540.566618, -5547.419917, -5554.274272, -5561.129681, -5567.986143, -5574.843657, 
		-5581.702222, -5588.561837, -5595.422500, -5602.284212, -5609.146970, -5616.010773, -5622.875621, -5629.741512, -5636.608445, -5643.476419, 
		-5650.345434, -5657.215487, -5664.086579, -5670.958707, -5677.831871, -5684.706069, -5691.581301, -5698.457566, -5705.334862, -5712.213188, 
		-5719.092544, -5725.972928, -5732.854339, -5739.736777, -5746.620240, -5753.504726, -5760.390236, -5767.276768, -5774.164320, -5781.052893, 
		-5787.942484, -5794.833093, -5801.724719, -5808.617361, -5815.511017, -5822.405687, -5829.301370, -5836.198064, -5843.095769, -5849.994483, 
		-5856.894207, -5863.794937, -5870.696674, -5877.599417, -5884.503164, -5891.407915, -5898.313668, -5905.220423, -5912.128178 };
	
	logP = -Mu + n*log(Mu);
	if (n >= 2) logP += logTable[TMath::Min(n,1000)-2];
	
	for (Int_t i = 1001; i <= n; i++) logP -= log((Double_t) i);
	
	return exp(logP);
}


Double_t ExpectedCombined(Double_t ilum, Double_t flum, Double_t eff, Double_t feff, Double_t bck, Double_t fbck, Double_t ilum2, Double_t flum2, Double_t eff2, Double_t feff2, Double_t bck2, Double_t fbck2,Int_t type=0)
{
	
	Double_t CL95A = 0, precision = 1.e-4;
	Int_t i;
	Int_t j;
	for (j = bck2; j>=0; j--){
	for (i = bck; i >= 0; i--)
	{
		//
		//Double_t s95 = roostats_cl95_bc_single(ilum,flum,eff,feff,bck,fbck, i,type);
		Double_t s95 = roostats_cl95_bc_combined(ilum,flum,eff,feff,bck,fbck,i,ilum2,flum2,eff2,feff2,bck2,fbck2,j,type);
		Double_t s95w =s95*Poisson(bck,i)*Poisson(bck2,j);
		CL95A += s95w;
		cout << "n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;	
		cout <<endl << type << " " << j << " " << i << endl<< endl;
		//
		if (s95w < CL95A*precision) break;
	}
	cout << "Lower bound on n has been found at " << i+1 << endl;
	//
	for (i = bck+1; ; i++)
	{
		Double_t s95 = roostats_cl95_bc_combined(ilum,flum,eff,feff,bck,fbck,i,ilum2,flum2,eff2,feff2,bck2,fbck2,j,type);
		Double_t s95w =s95*Poisson(bck,i)*Poisson(bck2,j);
		CL95A += s95w;
		cout << "n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;
		cout <<endl << type << " " << j << " " << i << endl<< endl;
		//
		if (s95w < CL95A*precision) break;
	}
	cout << "Upper bound on n has been found at " << i << endl;
	//
	}

	for (j = bck2+1; ; j++){
	Double_t s95w2=0.;
	for (i = bck; i >= 0; i--)
	{
		//
		//Double_t s95 = roostats_cl95_bc_single(ilum,flum,eff,feff,bck,fbck, i,type);
		Double_t s95 = roostats_cl95_bc_combined(ilum,flum,eff,feff,bck,fbck,i,ilum2,flum2,eff2,feff2,bck2,fbck2,j,type);
		Double_t s95w =s95*Poisson(bck,i)*Poisson(bck2,j);
		CL95A += s95w;
		s95w2=s95*Poisson(bck2,j);
		cout << "n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;
		cout <<endl << type << " " << j << " " << i << endl<< endl;
		//
		if (s95w < CL95A*precision) break;
	}
	cout << "Lower bound on n has been found at " << i+1 << endl;
	//
	for (i = bck+1; ; i++)
	{
		Double_t s95 = roostats_cl95_bc_combined(ilum,flum,eff,feff,bck,fbck,i,ilum2,flum2,eff2,feff2,bck2,fbck2,j,type);
		Double_t s95w =s95*Poisson(bck,i)*Poisson(bck2,j);
		CL95A += s95w;
		s95w2=s95*Poisson(bck2,j);

		cout << "n = " << i << "; 95% C.L. = " << s95 << " pb; weighted 95% C.L. = " << s95w << " pb; running <s95> = " << CL95A << " pb" << endl;
		cout <<endl << type << " " << j << " " << i << endl<< endl;
		//
		if (s95w < CL95A*precision) break;
	}
	cout << "Upper bound on n has been found at " << i << endl;
	if (s95w2 < CL95A*precision) break;
	//
	}


	cout << "Average upper 95% C.L. limit = " << CL95A << " pb" << endl;
	return CL95A;
}

