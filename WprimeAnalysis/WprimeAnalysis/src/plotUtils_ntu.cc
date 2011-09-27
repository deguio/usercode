#include "plotUtils_ntu.h"
#include "setTDRStyle.h"
#include "WprimeTreeVariables.h"




double MyGetMinimum_ntu(const TH1F* histo, const double& minval, int binMin, int binMax)
{
  if(binMin == -1) binMin = 1;
  if(binMax == -1) binMax = histo -> GetNbinsX();
  
  double minimum = histo -> GetMaximum();
  double value = 0.;
  
  for(int bin = binMin; bin <= binMax; ++bin)
  {
    value = histo -> GetBinContent(bin);
    if( (value > minval) && (value < minimum) )
      minimum = value;
  }

  return minimum;
}






drawTStack_ntu::drawTStack_ntu(const std::string& inputDir,
                       const std::string& listFileName,
                       const std::string& baseRootFileName,
                       const std::string& outputDir):
 m_inputDir(inputDir),
 m_listFileName(listFileName),
 m_baseRootFileName(baseRootFileName),
 m_outputDir(outputDir),
 m_xAxisRange(false),
 m_xRangeMin(0.),
 m_xRangeMax(1.),
 m_xAxisTitle(false),
 m_xTitle(""),
 m_yAxisRange(false),
 m_yRangeMin(0.),
 m_yRangeMax(1.),
 m_yAxisTitle(false),
 m_yTitle(""),
 m_drawLegend(true),
 m_xLowLegend(0.68),
 m_yLowLegend(0.78),
 m_xHighLegend(0.99),
 m_yHighLegend(0.99)
{
  
  //----------------------------------------------------------
  // read the file with the list of samples and cross sections
  //----------------------------------------------------------
  
  //fede  std::string listFullFileName = inputDir+listFileName;
  std::string listFullFileName = listFileName;
  std::ifstream listFile( listFullFileName.c_str() );
  if(!listFile.is_open())
  {
    std::cout << "\n>>>plotUtils::Error opening file " << listFullFileName << std::endl;
    exit(-1);
  }
  else
    std::cout << "\n>>>plotUtils::Opening file " << listFullFileName << std::endl;
  
  
  
  while(!listFile.eof())
  {
    std::string sample;
    std::string sumName;
    int dataFlag;
    int isDD;
    double mH;
    double crossSection;
    int color;
    int recoil;
    int ptHatCut;
    
    listFile >> sample >> sumName >> dataFlag >> isDD >> mH >> crossSection >> color >> recoil >> ptHatCut;

    if(sample.size() == 0)
      continue;
    if(sample.at(0) == '#')
      continue;
    
    std::cout << std::fixed << std::setprecision(1) << std::setw(90)
	      << sample << " "
              << std::setw(20)
              << sumName << " "
	      << std::setw(5)
              << dataFlag << " "
	      << std::setw(5)
              << isDD << " "
              << std::setw(10)
              << mH << " "
              << std::setw(12)
              << crossSection << " " 
              << std::setw(8)
	      << color << " "
              << std::setw(8)
	      << recoil << " "
              << std::setw(8)
	      << ptHatCut << " "
              << std::endl;
    
    std::pair<std::string, std::string> dummyPair(sample, sumName);
    m_list.push_back(dummyPair);
    m_dataFlag[sample] = dataFlag;
    m_isDD[sample] = isDD;
    m_mH[sample] = mH;
    m_crossSection[sample] = crossSection;
    m_color[sample] = color;
    
  }
  
  listFile.close();

  std::cout << ">>>plotUtils::Read " << m_list.size() << " samples" << std::endl;
  std::cout << ">>>plotUtils::Closing file " << listFullFileName << "\n" << std::endl;
}






drawTStack_ntu::~drawTStack_ntu()
{
  delete c1;
}



void drawTStack_ntu::Draw(std::vector<std::string>& variableNames, const std::string& histoName,
			  const std::string& mode,
			  const float& lumi, const int& step,
			  const int& nBins, const bool& logy,
			  std::vector<std::string>* cut,
			  const bool& wantCumulative)
{ 
  std::cout << "\n>>>plotUtils::Draw::Drawing histogram " << histoName << std::endl;

  
  //---------------------------------------------
  // define the map with summed cross sections
  //---------------------------------------------
  
  std::map<std::string, double> crossSection_summed;
  std::map<std::string, int> color_summed;
  std::map<std::string, bool> isFirstSample_summed;
  std::map<std::string, int> dataFlag_summed;
  std::map<std::string, int> isDD_summed;
  std::map<std::string, TH1F*> histo_summed;
  std::map<std::string, double> m_mH_summed;
  
  //initialize summed vectors
  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    crossSection_summed[vecIt->second] = 0.;
    color_summed[vecIt->second] = 1;
    isFirstSample_summed[vecIt->second] = true;
    dataFlag_summed[vecIt->second] = false;
    isDD_summed[vecIt->second] = false;
    histo_summed[vecIt->second] = NULL;
    m_mH_summed[vecIt->second] = 0;
  }
  
  
  
  //---------------------------------------------
  // loop over all the samples and fill the stack
  //---------------------------------------------
  int binMin = -1;
  int binMax = -1;
  
  int i = 0;
  std::vector<TTree*> trees;
  std::vector<TFile*> rootFiles;
  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    // open root file
    std::string fullRootFileName = m_inputDir+vecIt->first+"/"+m_baseRootFileName+".root";
    //std::cout << "opening file: " << fullRootFileName << std::endl;
    rootFiles.push_back(new TFile(fullRootFileName.c_str(), "READ"));
    if(!(rootFiles.at(i))->IsOpen()) exit(-1);
    
    // Draw::get tree
    TTree* tree;
    char treeName[50];
    sprintf(treeName, "ntu_%d", step);
    rootFiles.at(i) -> GetObject(treeName, tree);
    trees.push_back(tree);
    
    // Draw:: dump tree into histogram
    TH1F* histo = new TH1F(histoName.c_str(), "", nBins, m_xRangeMin, m_xRangeMax);    
    for(unsigned int jj = 0; jj < variableNames.size(); ++jj)
    {
      //std::cout << "Draw::Dumping tree variable " << (variableNames.at(jj)+">>"+histoName).c_str() << std::endl;
      if(cut != NULL)
	tree -> Draw( (variableNames.at(jj)+" >>+ "+histoName).c_str(), (cut->at(jj)).c_str() );
      else
        tree -> Draw( (variableNames.at(jj)+" >>+ "+histoName).c_str());
    }
    
    //valuto la cumumativa
    if (wantCumulative == true)
      {
	TH1F* histo_cumulative = new TH1F("histo_cumulative", "", nBins, m_xRangeMin, m_xRangeMax); 
	for (int kk = 1; kk < histo->GetNbinsX()+1; ++kk)
	  histo_cumulative -> SetBinContent(kk, histo->Integral(kk,1000000));
	histo = (TH1F*)histo_cumulative->Clone("histo");
      }

    std::string eventsHistoName = "events";
    //std::cout << "getting histogram " << eventsHistoName << std::endl;
    
    TH1F* eventsHisto = NULL;
    rootFiles.at(i) -> GetObject(eventsHistoName.c_str(), eventsHisto);
    if(eventsHisto == NULL)
    {
      std::cout << ">>>plotUtils::Error in getting object " << eventsHistoName << " in file: " << fullRootFileName << std::endl;
      exit(-1);
    }
    
    
    
    // scale histograms normalizing to lumi (1. pb^-1)
    // if data do not apply any scale factor
    histo -> Sumw2();
    dataFlag_summed[vecIt->second] = m_dataFlag[vecIt->first];
    isDD_summed[vecIt->second] = m_isDD[vecIt->first];
    color_summed[vecIt->second] = m_color[vecIt->first];
    m_mH_summed[vecIt->second] = m_mH[vecIt->first];

    //se sono dati o se e' DD non scalo e non sommo le sezioni d'urto
    if( (histo->GetEntries() > 0.) && (m_dataFlag[vecIt->first] != 1) && (m_isDD[vecIt->first] != 1))
      {
	histo -> Scale(m_crossSection[vecIt->first]/eventsHisto->GetBinContent(1));
	crossSection_summed[vecIt->second] += m_crossSection[vecIt->first];
      }

    if( (histo->GetEntries() > 0.) && (m_dataFlag[vecIt->first] != 1) && (m_isDD[vecIt->first] == 1))
      {
	histo -> Scale(m_crossSection[vecIt->first] / lumi); //applico il fattore di scala del prescale (poi scalo *lumi)    
	crossSection_summed[vecIt->second] += m_crossSection[vecIt->first];
      }
    
    
    // sum histograms normalized to lumi (1. pb^-1)
    // if data, just sum histograms w/o normalization
    if( isFirstSample_summed[vecIt->second] == false )
    {
      histo_summed[vecIt->second] -> Add(histo);
    }
    if( isFirstSample_summed[vecIt->second] == true )
    {
      histo_summed[vecIt->second] = (TH1F*)(histo -> Clone());
      isFirstSample_summed[vecIt->second] = false;
      
      binMin = histo -> FindBin(m_xRangeMin);
      binMax = histo -> FindBin(m_xRangeMax);
    }
    
    
    
    ++i;
  }
  
  
  
  
  
  
  //---------------------------------------------
  // Define the final Stack/Histograms
  //---------------------------------------------
  
  THStack* hs = new THStack("hs", "hs");
  TH1F* histoData = NULL;
  std::vector<TH1F*> histoSig;

  TLegend legend(0.68, 0.78, 0.99, 0.99);
  legend.SetFillColor(kWhite);
  
  double globalMaximum = 0.;
  double globalMinimum = 999999999999.;
  
  
  char lumiBuffer[50];
  sprintf(lumiBuffer, "CMS Preliminary 2011");
  
  char lumiBuffer2[50];
  sprintf(lumiBuffer2, "#sqrt{s}=7 TeV   L=%.3f pb^{-1}", lumi);

  TLatex *latex = new TLatex(0.3, 0.91, lumiBuffer); 
  TLatex *latex2 = NULL;
  if( mode == "eventsScaled" )
    latex2 = new TLatex(0.3, 0.88, lumiBuffer2);
  if( mode == "sameAreaStack" )
    latex2 = new TLatex(0.3, 0.88, "norm. to same area");
  if( mode == "sameAreaNoStack" )
    latex2 = new TLatex(0.3, 0.88, "#splitline{each histo norm.}{to unit area}");
  //if( mode == "eventsScaled" )
  //  latex = new TLatex(0.76, 0.50, lumiBuffer);
  //if( mode == "sameAreaStack" )
  //  latex = new TLatex(0.76, 0.50, "norm. to same area");
  //if( mode == "sameAreaNoStack" )
  //  latex = new TLatex(0.76, 0.50, "#splitline{each histo norm.}{to unit area}");
  
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(0.03);
  latex2->SetNDC();
  latex2->SetTextFont(42);
  latex2->SetTextSize(0.03);
  
  
  // compute normalization factor
  double globalGlobalIntegral = 0.;
  double dataGlobalGlobalIntegral = 0.;
  if( mode == "sameAreaStack" )
  {
    TH1F* globalGlobalHisto = NULL;
    TH1F* dataGlobalGlobalHisto = NULL;
    
    bool isFirstSample = true;
    bool isFirstSample_data = true;
    for(std::map<std::string, double>::const_iterator mapIt = crossSection_summed.begin();
        mapIt != crossSection_summed.end(); ++mapIt)
    {    
      TH1F* globalHisto = histo_summed[mapIt->first];
      
      if( (isFirstSample == false) && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
      {
        globalGlobalHisto -> Add(globalHisto);
      }
      if( (isFirstSample == true) && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
      {
        globalGlobalHisto = (TH1F*)(globalHisto -> Clone());
        isFirstSample = false;
      }

      if( (isFirstSample_data == false) && (dataFlag_summed[mapIt->first] == 1) )
      {
        dataGlobalGlobalHisto -> Add(globalHisto);
      }
      if( (isFirstSample_data == true) && (dataFlag_summed[mapIt->first] == 1) )
      {
        dataGlobalGlobalHisto = (TH1F*)(globalHisto -> Clone());
        isFirstSample_data = false;
      }

    }
    
    globalGlobalIntegral = globalGlobalHisto -> Integral(1, globalGlobalHisto->GetNbinsX());
    dataGlobalGlobalIntegral = dataGlobalGlobalHisto -> Integral(1, dataGlobalGlobalHisto->GetNbinsX());
  }
  
  
  
  // loop over summed histograms
  i = 0;
  TH1F* globalGlobalHisto = NULL;
  bool isFirstSample = true;

  //  std::cout << "FEDE >>>> crossSection_summed.size() = " << crossSection_summed.size() << std::endl;

  for(std::map<std::string, double>::const_iterator mapIt = crossSection_summed.begin();
      mapIt != crossSection_summed.end(); ++mapIt)
  {    
    TH1F* globalHisto = histo_summed[mapIt->first];
    if(globalHisto -> GetEntries() == 0) continue;
    
    if(m_xAxisRange)
      globalHisto->GetXaxis()->SetRangeUser(m_xRangeMin, m_xRangeMax);
    
    
    
    
    if( (mode == "eventsScaled") && (dataFlag_summed[mapIt->first] == 1) )
    {
      globalHisto -> SetMarkerStyle(20);
      globalHisto -> SetMarkerSize(0.7);
      
      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax);
      
      histoData = globalHisto;
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "P");
    }
    
    //scalo anche il DD per la lumi che ho provveduto a scalare al contrario
    //setto i colori solo se non sono dati
    if( (mode == "eventsScaled") && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
    {
      globalHisto -> Scale(1. * lumi);
      //fede      globalHisto -> SetLineColor(getColor(i+1));
      //fede      globalHisto -> SetFillColor(getColor(i+1));
      //globalHisto -> SetLineColor(color_summed[mapIt->first]);
      globalHisto -> SetFillColor(color_summed[mapIt->first]);
      //fede globalHisto -> SetFillStyle(3003);
      //fede globalHisto -> SetLineWidth(2);
      
      hs -> Add(globalHisto);
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "F");
    }

    if( (mode == "eventsScaled") && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] > 0))
    {
      globalHisto -> Scale(1. * lumi);
      //fede      globalHisto -> SetLineColor(getColor(i+1));
      //fede      globalHisto -> SetFillColor(getColor(i+1));
      //globalHisto -> SetLineColor(color_summed[mapIt->first]);
      globalHisto -> SetLineColor(color_summed[mapIt->first]);
      //fede globalHisto -> SetFillStyle(3003);
      globalHisto -> SetLineWidth(2);
      
      histoSig.push_back(globalHisto);
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "L");
    }
    
    
    
    
    if( (mode == "sameAreaNoStack") )
    {
      globalHisto -> Scale(1./globalHisto->Integral(1, globalHisto->GetNbinsX()));
      //fede   globalHisto -> SetLineColor(getColor(i));
      //globalHisto -> SetLineColor(color_summed[mapIt->first]);
      globalHisto -> SetLineStyle(i+1);
      //fede globalHisto -> SetLineWidth(4);

      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax);
      
      hs -> Add(globalHisto);      
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "L");      
    }
    
    
    
    
    if( (mode == "sameAreaStack") && (dataFlag_summed[mapIt->first] == 1) )
    {

      //      std::cout << "FEDE >>> dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()) = " 
      //		<< dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()) << std::endl;

      globalHisto -> Scale(dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()));
      globalHisto -> SetMarkerStyle(20);
      globalHisto -> SetMarkerSize(0.7);
      
      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum_ntu(globalHisto, 1.E-15, binMin, binMax);
      
      histoData = globalHisto;
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "P");      
    }
    
    if( (mode == "sameAreaStack") && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
    {
      globalHisto -> Scale(dataGlobalGlobalIntegral/globalGlobalIntegral);
      //fede      globalHisto -> SetLineColor(getColor(i+1));
      //fede      globalHisto -> SetFillColor(getColor(i+1));
      //globalHisto -> SetLineColor(color_summed[mapIt->first]);
      globalHisto -> SetFillColor(color_summed[mapIt->first]);
      //fede globalHisto -> SetFillStyle(3003);
      //fede globalHisto -> SetLineWidth(2);
      
      hs -> Add(globalHisto);      
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "F");      
    }
    
    if( (mode == "sameAreaStack") && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] > 0))
    {
      globalHisto -> Scale(dataGlobalGlobalIntegral/globalGlobalIntegral);
      //fede      globalHisto -> SetLineColor(getColor(i+1));
      //fede      globalHisto -> SetFillColor(getColor(i+1));
      //globalHisto -> SetLineColor(color_summed[mapIt->first]);
      globalHisto -> SetLineColor(color_summed[mapIt->first]);
      //fede globalHisto -> SetFillStyle(3003);
      globalHisto -> SetLineWidth(2);
      
      histoSig.push_back(globalHisto);
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "L");
    }
    
    
    
    if( (isFirstSample == false) && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
    {
      globalGlobalHisto -> Add(globalHisto);
    }
    if( (isFirstSample == true) && (dataFlag_summed[mapIt->first] != 1) && (m_mH_summed[mapIt->first] < 0))
    {
      globalGlobalHisto = (TH1F*)(globalHisto -> Clone());
      isFirstSample = false;
    }
    
    
    ++i;
  }  // loop over summed histograms
  
  
  
  // if no histograms had entries, return
  if(i == 0)
  {
    std::cout << ">>>plotUtils::Error, histograms empty" << std::endl;

    // close root files
    for(unsigned int i = 0; i < rootFiles.size(); ++i)
    {
      rootFiles.at(i) -> Close();
      delete rootFiles.at(i);
    }
    return;
  }
  
  
  
  
  
  if( (mode != "sameAreaNoStack") )
    {
      if(globalGlobalHisto->GetMaximum() > globalMaximum)
	globalMaximum = globalGlobalHisto -> GetMaximum();
      if(MyGetMinimum_ntu(globalGlobalHisto, 1.E-15, binMin, binMax) < globalMinimum)
	globalMinimum = MyGetMinimum_ntu(globalGlobalHisto, 1.E-15, binMin, binMax);
    }
  
  
  
  // draw the stack and save file
  SaveTStack(hs, histoName, histoData); //salvo l'histo in ogni caso

  TCanvas* c1 = new TCanvas();
  c1 -> cd();
  TPad* p1 = new TPad("p1","p1",0., 0.25, 1., 1.);
  TPad* p2 = new TPad("p2","p2",0., 0., 1., 0.25);
  p1 -> Draw();
  p2 -> Draw();
  
  p1 -> cd();
  p1 -> SetGridx();
  p1 -> SetGridy();
  if(logy) p1 -> SetLogy();

  
  
  if( mode == "eventsScaled" )
  {
    hs -> Draw("HISTO");
    //DrawTStackError(hs);
    if(histoData != NULL) histoData -> Draw("P,same");
    for (unsigned int ii = 0; ii < histoSig.size(); ++ii)
      histoSig.at(ii)->Draw("HISTO,same");
     
    

    char buffer[50];
    sprintf(buffer, "events");
    hs->GetYaxis()->SetTitle(buffer);

    //draw second pad
    if(histoData != NULL)
    {
      TH1F* ratioHisto = (TH1F*)(histoData -> Clone());
      ratioHisto -> GetYaxis() -> SetTitleSize(0.09);
      ratioHisto -> GetYaxis() -> SetTitleOffset(1.0);
      ratioHisto -> GetYaxis() -> SetTitle("data/MC");
      
      TObjArray* histos = hs -> GetStack();
      int nHistos = histos -> GetEntries();
      TH1F* lastHisto = (TH1F*)(histos->At(nHistos-1))->Clone();
      
      int nPoints = lastHisto->GetNbinsX();
      TGraph* ratioGraph1s = new TGraph(2*nPoints);
      TGraph* ratioGraph2s = new TGraph(2*nPoints);

      for(int bin = 1; bin <= histoData->GetNbinsX(); ++bin)
      {
        if(globalGlobalHisto->GetBinContent(bin) == 0.) continue;
        ratioHisto -> SetBinContent(bin, 1.*histoData->GetBinContent(bin)/globalGlobalHisto->GetBinContent(bin));
        ratioHisto -> SetBinError(bin, 1.*histoData->GetBinError(bin)/globalGlobalHisto->GetBinContent(bin));
      }
      
      int point = 0;
      for(int bin = 1; bin <= lastHisto->GetNbinsX(); ++bin)
      {
        float binCenter = lastHisto->GetBinCenter(bin);
	//        float binError = sqrt(lastHisto->GetBinError(bin)*lastHisto->GetBinError(bin) + 
	//                         lastHisto->GetBinContent(bin));
        float binError = lastHisto->GetBinError(bin);

        float binErrorM;
        float binErrorP;
        float binError2M;
        float binError2P;
        if( lastHisto->GetBinContent(bin) != 0 )
        {
          binErrorM = 1. - binError/lastHisto->GetBinContent(bin);
          binErrorP = 1. + binError/lastHisto->GetBinContent(bin);
          binError2M = 1. - 2.*binError/lastHisto->GetBinContent(bin);
          binError2P = 1. + 2.*binError/lastHisto->GetBinContent(bin);
        }
        else
        {
          binErrorM = 1.;
          binErrorP = 1.;
          binError2M = 1.;
          binError2P = 1.;
        }
        ratioGraph1s -> SetPoint(point,binCenter,binErrorM);
        ratioGraph1s -> SetPoint(2*nPoints-point-1,binCenter,binErrorP);
        ratioGraph2s -> SetPoint(point,binCenter,binError2M);
        ratioGraph2s -> SetPoint(2*nPoints-point-1,binCenter,binError2P);
        
        ++point;
      }
      
      
      p2 -> cd();
      p2 -> SetGridx();
      p2 -> SetGridy();
      
      ratioHisto -> GetYaxis() -> SetRangeUser(0., 2.);
      ratioHisto -> Draw("P");
      
      ratioGraph2s -> SetFillColor(kYellow);
      ratioGraph2s -> SetLineColor(kYellow);
      //fede ratioGraph2s -> Draw("F,same");
      
      ratioGraph1s -> SetFillColor(kGreen);
      ratioGraph1s -> SetLineColor(kGreen);
      //fede ratioGraph1s -> Draw("F,same");
      
      ratioHisto -> Draw("P,same");
        
      TF1* line = new TF1("line", "1.", -1000000., 1000000.);
      line -> SetLineWidth(2);
      line -> SetLineColor(kRed);
      line -> Draw("same");
    }
    p1 -> cd();
  }

  if( mode == "sameAreaNoStack" )
  {
    hs -> Draw("HISTO,nostack");
    
    hs->GetYaxis()->SetTitle("event fraction");
  }
  
  if( mode == "sameAreaStack" )
  {
    hs -> Draw("HISTO");
    if(histoData != NULL) histoData -> Draw("P,same");
    for (unsigned int ii = 0; ii < histoSig.size(); ++ii)
      histoSig.at(ii)->Draw("HISTO,same");

    
    hs->GetYaxis()->SetTitle("event fraction");

    //draw second pad
    if(histoData != NULL)
    {
      TH1F* ratioHisto = (TH1F*)(histoData -> Clone());
      ratioHisto -> GetYaxis() -> SetTitleSize(0.09);
      ratioHisto -> GetYaxis() -> SetTitleOffset(1.0);
      ratioHisto -> GetYaxis() -> SetTitle("data/MC");
      
      TObjArray* histos = hs -> GetStack();
      int nHistos = histos -> GetEntries();
      TH1F* lastHisto = (TH1F*)(histos->At(nHistos-1))->Clone();
      
      int nPoints = lastHisto->GetNbinsX();
      TGraph* ratioGraph1s = new TGraph(2*nPoints);
      TGraph* ratioGraph2s = new TGraph(2*nPoints);

      for(int bin = 1; bin <= histoData->GetNbinsX(); ++bin)
      {
        if(globalGlobalHisto->GetBinContent(bin) == 0.) continue;
        ratioHisto -> SetBinContent(bin, 1.*histoData->GetBinContent(bin)/globalGlobalHisto->GetBinContent(bin));
        ratioHisto -> SetBinError(bin, 1.*histoData->GetBinError(bin)/globalGlobalHisto->GetBinContent(bin));
      }
      
      int point = 0;
      for(int bin = 1; bin <= lastHisto->GetNbinsX(); ++bin)
      {
        float binCenter = lastHisto->GetBinCenter(bin);
        // float binError = sqrt(lastHisto->GetBinError(bin)*lastHisto->GetBinError(bin) + 
        //                       lastHisto->GetBinContent(bin));
        float binError = lastHisto->GetBinError(bin);

        float binErrorM;
        float binErrorP;
        float binError2M;
        float binError2P;
        if( lastHisto->GetBinContent(bin) != 0 )
        {
          binErrorM = 1. - binError/lastHisto->GetBinContent(bin);
          binErrorP = 1. + binError/lastHisto->GetBinContent(bin);
          binError2M = 1. - 2.*binError/lastHisto->GetBinContent(bin);
          binError2P = 1. + 2.*binError/lastHisto->GetBinContent(bin);
        }
        else
        {
          binErrorM = 1.;
          binErrorP = 1.;
          binError2M = 1.;
          binError2P = 1.;
        }
        ratioGraph1s -> SetPoint(point,binCenter,binErrorM);
        ratioGraph1s -> SetPoint(2*nPoints-point-1,binCenter,binErrorP);
        ratioGraph2s -> SetPoint(point,binCenter,binError2M);
        ratioGraph2s -> SetPoint(2*nPoints-point-1,binCenter,binError2P);
        
        ++point;
      }
      
      
      p2 -> cd();
      p2 -> SetGridx();
      p2 -> SetGridy();
      
      ratioHisto -> GetYaxis() -> SetRangeUser(0., 2.);
      ratioHisto -> Draw("P");
      
      ratioGraph2s -> SetFillColor(kYellow);
      ratioGraph2s -> SetLineColor(kYellow);
      //fede ratioGraph2s -> Draw("F,same");
      
      ratioGraph1s -> SetFillColor(kGreen);
      ratioGraph1s -> SetLineColor(kGreen);
      //fede ratioGraph1s -> Draw("F,same");
      
      ratioHisto -> Draw("P,same");
        
      TF1* line = new TF1("line", "1.", -1000000., 1000000.);
      line -> SetLineWidth(2);
      line -> SetLineColor(kRed);
      line -> Draw("same");
    }
    p1 -> cd();
  }

  if(m_drawLegend == true)
    {
      legend.Draw("same");
      latex->Draw("same");
      latex2->Draw("same");
    }
  
  
  
  
  
  
  // set x-axis properties
  std::string fullTitle = histoName;
  char stepTitle[50];
  sprintf(stepTitle, "_%d", step);
  hs->GetXaxis()->SetTitle(fullTitle.c_str());
  if(m_xAxisTitle) hs->GetXaxis()->SetTitle(m_xTitle.c_str());
  hs -> GetXaxis() -> SetTitleSize(0.04);
  hs -> GetXaxis() -> SetLabelSize(0.03); //fede
  hs -> GetXaxis() -> SetTitleOffset(1.25);

  if(m_xAxisRange) hs->GetXaxis()->SetRangeUser(m_xRangeMin, m_xRangeMax);
  
  
  
  // set y-axis properties
  if(m_yAxisTitle) hs->GetYaxis()->SetTitle(m_yTitle.c_str());    
  hs -> GetYaxis() -> SetTitleSize(0.04);
  hs -> GetYaxis() -> SetLabelSize(0.03); //fede
  hs -> GetYaxis() -> SetTitleOffset(1.50);
  
  hs->SetMinimum(0.);
  hs->SetMaximum(globalMaximum+0.1*globalMaximum);

  if(logy)
  {
    hs->SetMinimum(pow(10., log10(globalMinimum) - 0.1));
    hs->SetMaximum(pow(10., log10(globalMaximum) + 0.1));
  } 
  
  if(m_yAxisRange)
  {
    hs->SetMinimum(m_yRangeMin);
    hs->SetMaximum(m_yRangeMax);
  }
  
  
  
  
  
  
  // write plots
  struct stat st;
  if(stat(m_outputDir.c_str(), &st) != 0)
  {
    std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
    exit(-1);
  }
  c1->Print((m_outputDir+fullTitle+".pdf").c_str(), "pdf");
  //c1->Print((m_outputDir+fullTitle+".root").c_str(), "root");
  
  
  
  // close root files
  for(unsigned int i = 0; i < rootFiles.size(); ++i)
  {
    rootFiles.at(i) -> Close();
    delete rootFiles.at(i);
  }
  
  m_xAxisRange = false;
  m_xAxisTitle = false;
  m_yAxisRange = false;
  m_yAxisTitle = false;
  
  delete c1;
  delete hs;
}




void drawTStack_ntu::DrawEvents(const std::string& mode, const float& lumi, const int& step, const bool& logy)
{ 
  std::cout << "\n>>>plotUtils::Drawing " << mode << std::endl;
  
  
  
  THStack* hs = new THStack("hs", "hs");
  TLegend legend(m_xLowLegend, m_yLowLegend, m_xHighLegend, m_yHighLegend);
  legend.SetFillColor(kWhite);
  
  double globalMaximum = 0.;
  double globalMinimum = 1.;
  
  
  //---------------------------------------------
  // define the map with summed cross sections
  //---------------------------------------------
 
  std::map<std::string, double> crossSection_summed;
  std::map<std::string, int> color_summed;
  std::map<std::string, bool> isFirstSample_summed;
  std::map<std::string, int> dataFlag_summed;
  std::map<std::string, bool> isSignal_summed;
  std::map<std::string, TH1F*> histo_summed;

  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    histo_summed[vecIt->second] = NULL;
    isFirstSample_summed[vecIt->second] = true;
    isSignal_summed[vecIt->second] = false;
  }
  
  
  
  //---------------------------------------------
  // loop over all the samples and fill the stack
  //---------------------------------------------

  std::vector<TFile*> rootFiles;
  TH1F* histoData = 0;
  TH1F* globalGlobalHisto = NULL;
  int i = 0;
  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    // open root file
    std::string fullRootFileName = m_inputDir+vecIt->first+"/"+m_baseRootFileName+".root";
    //std::cout << "opening file: " << fullRootFileName << std::endl;
    rootFiles.push_back(new TFile(fullRootFileName.c_str(), "READ"));
    
    
    // get histogram
    std::string fullHistoName = "events";
    //std::cout << "getting histogram " << fullHistoName << std::endl;
    TH1F* histo = NULL;
    rootFiles.at(i) -> GetObject(fullHistoName.c_str(), histo);
    if(histo == NULL)
    {
      std::cout << ">>>plotUtils::Error in getting object " << fullHistoName << " in file: " << fullRootFileName << std::endl;
      exit(-1);
    }
    
    crossSection_summed[vecIt->second] += m_crossSection[vecIt->first];
    color_summed[vecIt->second] = m_color[vecIt->first];
    if(m_mH[vecIt->first] > 0)
      isSignal_summed[vecIt->second] = true;
    
    
    if( (mode == "eventsScaled") ||
        (mode == "eventsScaledStack") ||
        (mode == "efficiencies") ||
        (mode == "efficienciesRelative") ||
        (mode == "significance") )
    {
      if( m_dataFlag[vecIt->first] != 1 )
	{
	  histo -> Scale(lumi*m_crossSection[vecIt->first]/histo->GetBinContent(1));
	}
     }
    
    
    
    // sum histograms
    if( isFirstSample_summed[vecIt->second] == false )
    {
      histo_summed[vecIt->second] -> Add(histo);
    }
    if( isFirstSample_summed[vecIt->second] == true )
    {
      histo_summed[vecIt->second] = (TH1F*)(histo -> Clone());
      isFirstSample_summed[vecIt->second] = false;
    }
    
    
    
    // save data histogram
    dataFlag_summed[vecIt->second] = m_dataFlag[vecIt->first];
    if(m_dataFlag[vecIt->first] == 1)
     histoData = histo;
    
    ++i;
  }
  
  
  
  std::ofstream* outFile = NULL;

  if(mode == "events")
    outFile = new std::ofstream((m_outputDir+"events.txt").c_str(), std::ios::out);  
  if(mode == "eventsScaled")
    outFile = new std::ofstream((m_outputDir+"eventsScaled.txt").c_str(), std::ios::out);
  if(mode == "eventsScaledStack")
    {
      outFile = new std::ofstream((m_outputDir+"eventsScaledStack.txt").c_str(), std::ios::out);
    }
  if(mode == "efficiencies")
    outFile = new std::ofstream((m_outputDir+"efficiencies.txt").c_str(), std::ios::out);
  if(mode == "efficienciesRelative")
    outFile = new std::ofstream((m_outputDir+"efficienciesRelative.txt").c_str(), std::ios::out);  
  if(mode == "significance")
    outFile = new std::ofstream((m_outputDir+"significance.txt").c_str(), std::ios::out);   
  
  
  
  
  
  int nSamples = 0;
  for(std::map<std::string, double>::const_iterator mapIt = crossSection_summed.begin();
      mapIt != crossSection_summed.end(); ++mapIt)
    ++nSamples;
  
  TH1F* stepHisto = new TH1F("stepHisto", "", nSamples, 0., 1.*nSamples);
  
  
  
  i = 0;
  for(std::map<std::string, double>::const_iterator mapIt = crossSection_summed.begin();
      mapIt != crossSection_summed.end(); ++mapIt)
  {    
    TH1F* globalHisto = histo_summed[mapIt->first];
    
    
    
    if(mode == "efficiencies")
    {
      float totalEvents = globalHisto -> GetBinContent(1);
      //fede
      //std::cout << "crossSection_summed = " << (*mapIt).first << std::endl;

      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
	{

	  //fede
	  // if (bin == 1)
	  //   {
	      // std::cout << std::endl;
	      // std::cout << "step = " << globalHisto -> GetXaxis() -> GetBinLabel(bin) << std::endl;
	      // std::cout << "bin = " << bin << std::endl;
	      // std::cout << "globalHisto->GetBinContent(bin) = " << globalHisto->GetBinContent(bin) << "; totalEvents = " << totalEvents << std::endl;
	      // std::cout << "ratio = " << globalHisto->GetBinContent(bin)/totalEvents << std::endl;
	      //	    }


	  globalHisto->SetBinContent(bin, globalHisto->GetBinContent(bin)/totalEvents);


	}
    }
    
    if(mode == "efficienciesRelative")
    {
      std::map<int, float > totalEvents;
      totalEvents[0] = globalHisto->GetBinContent(1);
      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
        totalEvents[bin] = globalHisto->GetBinContent(bin);
      
      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
      {        
        globalHisto->SetBinContent(bin, globalHisto->GetBinContent(bin)/totalEvents[bin-1]);  //-->> fede: err!
      }
    }
    
    if(mode == "significance")
    {
      if(isSignal_summed[mapIt->first] == false) continue;
    
      TH1F* totalBkgHisto = NULL;
      for(std::map<std::string, double>::const_iterator mapIt2 = crossSection_summed.begin();
          mapIt2 != crossSection_summed.end(); ++mapIt2)
      {
        if(m_mH[mapIt2->first] < 0) continue;
        TH1F* globalHisto2 = histo_summed[mapIt2->first];
        
        if(totalBkgHisto == NULL)
          totalBkgHisto = (TH1F*)(globalHisto2 -> Clone());
        else
          totalBkgHisto -> Add(globalHisto2);
      }
      
      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
        globalHisto->SetBinContent(bin, globalHisto->GetBinContent(bin) / 
                                        sqrt(totalBkgHisto->GetBinContent(bin)) );
    }
    
    
    //fede    globalHisto -> SetLineColor(getColor(i));
    globalHisto -> SetLineColor(color_summed[mapIt->first]);
    globalHisto -> SetLineStyle(i+1);
    if(m_xAxisRange)
      globalHisto->GetXaxis()->SetRangeUser(m_xRangeMin, m_xRangeMax);
    if(i == 0)
      globalHisto -> SetLineWidth(4);
    else
      globalHisto -> SetLineWidth(4);
        
    
    if(globalHisto->GetMaximum() > globalMaximum)
      globalMaximum = globalHisto -> GetMaximum();
    if( (globalHisto->GetMinimum() < globalMinimum) && (globalHisto->GetMinimum() > 0.) )
      globalMinimum = globalHisto -> GetMinimum();    
    
    
    
    
    
    
    if( ( (mode != "eventsScaled") && (mode != "eventsScaledStack") ) ||
        ( (mode == "eventsScaled")      && (dataFlag_summed[mapIt->first] != 1) ) || 
        ( (mode == "eventsScaledStack") && (dataFlag_summed[mapIt->first] != 1) ) )
    {
      hs -> Add(globalHisto);
      
      if(globalGlobalHisto != NULL)
        globalGlobalHisto -> Add(globalHisto);
      else
        globalGlobalHisto = (TH1F*)(globalHisto -> Clone());
    }
    
    legend.AddEntry(globalHisto, (mapIt->first).c_str(), "L");
    
    
    stepHisto -> SetBinContent(i+1, globalHisto->GetBinContent(step));
    stepHisto -> GetXaxis() -> SetBinLabel(i+1, (mapIt->first).c_str());
    
    
    
    
    
    
    if(outFile)
    {
      (*outFile) << "\n" << mapIt->first;
      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
      {
        const char* binLabel = globalHisto -> GetXaxis() -> GetBinLabel(bin);
        (*outFile) << "   " << binLabel << ":   " << std::fixed
                                                  << std::setprecision(4)
                                                  << std::scientific
                                                  << 1. * globalHisto -> GetBinContent(bin);
      }
      
      (*outFile) << "\n";
    }
  
    
    ++i;
  }
  
  
  
  
  
  
  // draw the stack and save file
  TCanvas* c1 = new TCanvas();
  c1 -> cd();
  c1 -> SetGridx();
  c1 -> SetGridy();
  
  if(mode != "eventsScaledStack")
  {
    hs -> Draw("nostack");
    //DrawTStackError(hs);
  }
  else
  {
    hs -> Draw("HISTO");
    //DrawTStackError(hs);
  }
  
  if( (mode != "significance") && 
      (mode != "efficiencies") &&
      (mode != "efficienciesRelative") )
    legend.Draw("same");
  
  
  if(logy)
    c1 -> SetLogy();
  if(logy)
    hs->GetYaxis()->SetRangeUser(globalMinimum - 0.1*globalMinimum,
                                 globalMaximum + 0.1*globalMaximum);
  
  
  hs -> GetXaxis() -> SetTitleSize(0.04);
  hs -> GetXaxis() -> SetLabelSize(0.03); //fede
  hs -> GetXaxis() -> SetTitleOffset(1.25);

  hs -> GetYaxis() -> SetTitleSize(0.04);
  hs -> GetYaxis() -> SetLabelSize(0.03); //fede
  hs -> GetYaxis() -> SetTitleOffset(1.50);
  
  
  if(m_xAxisRange)
    hs->GetXaxis()->SetRangeUser(m_xRangeMin, m_xRangeMax);
  
  hs->GetYaxis()->SetRangeUser(0., globalMaximum+0.1*globalMaximum);
  if(m_yAxisRange)
    hs->GetYaxis()->SetRangeUser(m_yRangeMin, m_yRangeMax);
  
  
  if(mode == "events")
  {
    hs->GetYaxis()->SetTitle("events left"); 
    stepHisto->GetYaxis()->SetTitle("events left");
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    c1->Print((m_outputDir+"events.pdf").c_str(), "pdf");
  }
  
  if(mode == "eventsScaled")
  {
    char buffer[50];
    sprintf(buffer, "events");
    hs->GetYaxis()->SetTitle(buffer);
    stepHisto->GetYaxis()->SetTitle(buffer);
    
    if(histoData != NULL)
    {
      histoData -> SetMarkerStyle(20);  
      histoData -> SetMarkerSize(0.7);
      histoData -> Draw("P,same");
    }
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    
    c1->Print((m_outputDir+"eventsScaled.pdf").c_str(), "pdf");
  }

  if(mode == "eventsScaledStack")
  {
    char buffer[50];
    sprintf(buffer, "events");
    hs->GetYaxis()->SetTitle(buffer);
    stepHisto->GetYaxis()->SetTitle(buffer);
    
    if(histoData != NULL)
    {
      histoData -> SetMarkerStyle(20);
      histoData -> SetMarkerSize(0.7);
      
      for(int bin = 1; bin <= histoData->GetNbinsX(); ++bin)
        histoData->SetBinError(bin, sqrt(histoData->GetBinContent(bin)));
      
      histoData -> Draw("P,same");
    }
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    
    c1->Print((m_outputDir+"eventsScaledStack.pdf").c_str(), "pdf");
    
    
    
    TCanvas* c2 = new TCanvas();
    for(int bin = 1; bin <= histoData->GetNbinsX(); ++bin)
    {
      histoData -> SetBinContent(bin, histoData->GetBinContent(bin)/globalGlobalHisto->GetBinContent(bin));
      histoData -> SetBinError(bin, histoData->GetBinError(bin)/globalGlobalHisto->GetBinContent(bin));
    }
    
    c2 -> cd();
    c2 -> SetGridx();
    c2 -> SetGridy();
    histoData -> Draw("P");
    c2->Print((m_outputDir+"eventsScaledStack_ratio.pdf").c_str(), "pdf");
    
    // save root file
    TFile* outRootFile = new TFile((m_outputDir+"eventsScaledStack_ratio.root").c_str(), "RECREATE");
    outRootFile -> cd();
    histoData -> Write();
    outRootFile -> Close();
    delete outRootFile;
    
    delete c2;
  }
    
  if(mode == "efficiencies")
  {
    hs->GetYaxis()->SetTitle("efficiency");
    stepHisto->GetYaxis()->SetTitle("efficiency");
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    c1->Print((m_outputDir+"efficiencies.pdf").c_str(), "pdf");
  }
  
  if(mode == "efficienciesRelative")
  {
    hs->GetYaxis()->SetTitle("relative efficiency");
    stepHisto->GetYaxis()->SetTitle("relative efficiency");
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    c1->Print((m_outputDir+"efficienciesRelative.pdf").c_str(), "pdf");  
  }

  if(mode == "significance")
  {
    hs->GetYaxis()->SetTitle("S / #sqrt{B}");
    
    struct stat st;
    if(stat(m_outputDir.c_str(), &st) != 0)
    {
      std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
      exit(-1);
    }
    c1->Print((m_outputDir+"significance.pdf").c_str(), "pdf");
  }  
  
  
  m_xAxisRange = false;
  m_xAxisTitle = false;
  m_yAxisRange = false;
  m_yAxisTitle = false;
  
  delete c1;
  delete hs;
  
  
  
  c1 = new TCanvas();
  c1 -> cd();
  c1 -> SetGridx();
  c1 -> SetGridy();
  c1 -> SetLogy();
  
  stepHisto -> GetXaxis() -> SetTitleSize(0.04);
  stepHisto -> GetXaxis() -> SetLabelSize(0.03); //fede
  stepHisto -> GetXaxis() -> SetTitleOffset(1.25);

  stepHisto -> GetYaxis() -> SetTitleSize(0.04);
  stepHisto -> GetYaxis() -> SetLabelSize(0.03); //fede
  stepHisto -> GetYaxis() -> SetTitleOffset(1.50);
  
  stepHisto -> SetLineWidth(3);
  
  stepHisto -> Draw();

  struct stat st;
  if(stat(m_outputDir.c_str(), &st) != 0)
  {
    std::cout << ">>>plotUtils::Error accessing directory " << m_outputDir << std::endl;
    exit(-1);
  }
  c1->Print( (m_outputDir+mode+"_step.pdf").c_str(), "pdf");
  
  delete c1;
  
  
  
  // save stepHisto root file
  TFile* outRootFile = new TFile((m_outputDir+mode+"_step.root").c_str(), "RECREATE");
  outRootFile -> cd();
  stepHisto -> Write();
  outRootFile -> Close();
  delete stepHisto;
  delete outRootFile;
  
  
  
  // close root files
  if(outFile)
    outFile -> close();
  
  for(unsigned int i = 0; i < rootFiles.size(); ++i)
  {
    rootFiles.at(i) -> Close();
    delete rootFiles.at(i);
  }
    
}

void drawTStack_ntu::SaveTStack(THStack* hs, const std::string& histoName, TH1F* histoData)
{
  TFile* fout = new TFile( (m_outputDir+histoName.c_str()+".root").c_str(),"RECREATE");

  TObjArray* histos = hs -> GetStack();
  int nHistos = histos -> GetEntries();
  TH1F* lastHisto = (TH1F*)(histos->At(nHistos-1))->Clone();
  lastHisto -> SetFillStyle(3005);
  lastHisto -> SetFillColor(kBlack);
  lastHisto -> SetMarkerSize(0);

  lastHisto -> Write( "histoMC" );
  histoData -> Write( "histoData" );
  
  fout -> Close();
}


void drawTStack_ntu::DrawTStackError(THStack* hs)
{
  TObjArray* histos = hs -> GetStack();
  int nHistos = histos -> GetEntries();
  TH1F* lastHisto = (TH1F*)(histos->At(nHistos-1))->Clone();

  lastHisto -> SetFillStyle(3005);
  lastHisto -> SetFillColor(kBlack);
  lastHisto -> SetMarkerSize(0);

  for(int bin = 0; bin <= lastHisto->GetNbinsX(); ++bin)
  {
    double oldErr = lastHisto -> GetBinError(bin);
    double poissonErr = sqrt(lastHisto -> GetBinContent(bin));
    lastHisto -> SetBinError(bin,sqrt(oldErr*oldErr + poissonErr*poissonErr));
    lastHisto -> SetBinError(bin,oldErr); //disegno solo la syst associata alla stat del MC

    //std::cout << "bin: " << bin << "   err: " << std::setprecision(4) << lastHisto -> GetBinError(bin) << std::endl;                                                             
  }

  lastHisto -> DrawClone("same,E2");
}


void drawTStack_ntu::SetXaxisRange(const double& xMin, const double& xMax)
{
  m_xAxisRange = true;
  m_xRangeMin = xMin;
  m_xRangeMax = xMax;
}

void drawTStack_ntu::SetXaxisTitle(const std::string& xTitle)
{
  m_xAxisTitle = true;
  m_xTitle = xTitle;
}



void drawTStack_ntu::SetYaxisRange(const double& yMin, const double& yMax)
{
  m_yAxisRange = true;
  m_yRangeMin = yMin;
  m_yRangeMax = yMax;
}

void drawTStack_ntu::SetYaxisTitle(const std::string& yTitle)
{
  m_yAxisTitle = true;
  m_yTitle = yTitle;
}



void drawTStack_ntu::SetDrawLegend(const bool& drawLegend)
{
  m_drawLegend = drawLegend;
}

void drawTStack_ntu::SetXLegend(const double& xLow, const double& xHigh)
{
  m_xLowLegend = xLow;
  m_xHighLegend = xHigh;
}

void drawTStack_ntu::SetYLegend(const double& yLow, const double& yHigh)
{
  m_yLowLegend = yLow;
  m_yHighLegend = yHigh;
}


//===========================
//==== limit calculation ====
//===========================
//esplicitamente per roostats_cl95:  https://indico.cern.ch:443/getFile.py/access?contribId=4&sessionId=3&resId=0&materialId=slides&confId=126178
void drawTStack_ntu::numbersForLimit(std::vector<std::string>& variableNames,
				     const float& lumi, const int& step,
				     std::vector<std::string>* cut,
				     const int& limitStep, 
				     const double& xLow, const double& xHigh)
{
  std::cout << "\n>>>plotUtils::Limit::evaluting numbers " << std::endl;

  
  //---------------------------------------------
  // loop over all the samples and fill the stack
  //---------------------------------------------
  int i = 0;
  std::vector<TTree*> trees;
  std::vector<TFile*> rootFiles;

  std::map<int,int> nEventsData;
  std::map<int,std::pair<float,float> > nEventsBkg;
  std::map<int, std::map<int,std::pair<float,float> > >nEventsSig; //per ogni segnale

  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    // open root file
    std::string fullRootFileName = m_inputDir+vecIt->first+"/"+m_baseRootFileName+".root";
    //std::cout << "opening file: " << fullRootFileName << std::endl;
    rootFiles.push_back(new TFile(fullRootFileName.c_str(), "READ"));
    if(!(rootFiles.at(i))->IsOpen()) exit(-1);
    
    // Draw::get tree
    TTree* tree;
    char treeName[50];
    sprintf(treeName, "ntu_%d", step);
    rootFiles.at(i) -> GetObject(treeName, tree);
    trees.push_back(tree);
    
    // Draw:: dump tree into histogram
    std::string histoName = "dummy";
    TH1F* histo = new TH1F(histoName.c_str(), "", 100000, 0., 100000);    //tutto il range
    for(unsigned int jj = 0; jj < variableNames.size(); ++jj)
    {
      //std::cout << "Draw::Dumping tree variable " << (variableNames.at(jj)+">>"+histoName).c_str() << std::endl;
      if(cut != NULL)
	tree -> Draw( (variableNames.at(jj)+" >>+ "+histoName).c_str(), (cut->at(jj)).c_str() );
      else
        tree -> Draw( (variableNames.at(jj)+" >>+ "+histoName).c_str());
    }
    
    std::string eventsHistoName = "events";
    //std::cout << "getting histogram " << eventsHistoName << std::endl;
    
    TH1F* eventsHisto = NULL;
    rootFiles.at(i) -> GetObject(eventsHistoName.c_str(), eventsHisto);
    if(eventsHisto == NULL)
    {
      std::cout << ">>>plotUtils::Error in getting object " << eventsHistoName << " in file: " << fullRootFileName << std::endl;
      exit(-1);
    }
    
    //std::cout << "sampleName = " << vecIt->first << std::endl;

    //calcolo integrali e numeri per limite ad ogni step
    for(float cut = xLow; cut <= xHigh; cut = cut+limitStep)
      {
	if (histo->GetEntries() <= 0.)
	  {
	    std::cout << "Limit Calculation: NO EVENTS! ERROR!" << std::endl;
	    continue;
	  }

	if(m_dataFlag[vecIt->first] == 1)
	    nEventsData[cut] = histo->Integral(histo->FindBin(cut), histo->GetNbinsX()+1);

	float scaleFactor = m_crossSection[vecIt->first]*lumi/eventsHisto->GetBinContent(1);
	if((m_dataFlag[vecIt->first] != 1) && (m_mH[vecIt->first] < 0))
	  {


	    float eff    =      histo->Integral(histo->FindBin(cut), histo->GetNbinsX()+1)  * scaleFactor;
	    float effErr = sqrt(histo->Integral(histo->FindBin(cut), histo->GetNbinsX()+1)) * scaleFactor;
	    
	    float effSummed    = nEventsBkg[cut].first  + eff;
	    float effErrSummed = nEventsBkg[cut].second + effErr*effErr;
	      
	    std::pair<float,float> dum (effSummed, effErrSummed);
	    nEventsBkg[cut] = dum;


	    //std::cout << "cutLow = " << cut << "; nEventsBkg[cut].first = " << nEventsBkg[cut].first << std::endl;

	  }

	if((m_dataFlag[vecIt->first] != 1) && (m_mH[vecIt->first] > 0))
	  {//per il segnale serve l'eff
	    float eff = histo->Integral(histo->FindBin(cut), histo->GetNbinsX()+1) / eventsHisto->GetBinContent(1);
	    float effErr = sqrt(histo->Integral(histo->FindBin(cut), histo->GetNbinsX()+1)) / eventsHisto->GetBinContent(1);

	    std::pair<float,float> dum (eff, effErr);
	    nEventsSig[cut][m_mH[vecIt->first]] = dum;
	  }
      }//loop over steps
    ++i;
  }

  //OUT LIMIT FILE
  std::ofstream* outFile = new std::ofstream((m_outputDir+"numbersForLimit.txt").c_str(), std::ios::out);  
  for (std::map<int, std::map<int,std::pair<float,float> > >::const_iterator cutItr = nEventsSig.begin(); cutItr != nEventsSig.end(); ++cutItr)
    {
      
      std::map<int,int>::const_iterator dataItr = nEventsData.find(cutItr->first);
      std::map<int,std::pair<float,float> >::const_iterator bkgItr = nEventsBkg.find(cutItr->first);
      (*outFile) << std::setprecision(4)
		 << std::scientific
		 << cutItr->first        << "       "
		 << dataItr->second      << "       "
		 << bkgItr->second.first << "       "
		 << sqrt(bkgItr->second.second)  << "       ";
      for(std::map<int,std::pair<float,float> >::const_iterator sigItr = cutItr->second.begin(); sigItr != cutItr->second.end(); ++sigItr)
	{
	  (*outFile) << std::setprecision(4)
		     << std::scientific
		     << sigItr->second.first        << "       "
		     << sqrt(sigItr->second.second)        << "       ";
	}//loop over signals
      (*outFile) << std::endl;
    }
  
  outFile->close();
}

//==================================
//==== print interesting events ====
//==================================
void drawTStack_ntu::printMtAboveThr(const int& step,
				     const double& thr)
{
  std::cout << "\n>>>plotUtils::Limit::print interesting events " << std::endl;

  
  //---------------------------------------------
  // loop over all the samples and fill the stack
  //---------------------------------------------
  int i = 0;
  std::vector<TTree*> trees;
  std::vector<TFile*> rootFiles;

  std::ofstream* outFile = NULL;
  outFile = new std::ofstream((m_outputDir+"interestingEvents.txt").c_str(), std::ios::out);  


  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
    {
      // open root file
      std::string fullRootFileName = m_inputDir+vecIt->first+"/"+m_baseRootFileName+".root";
      //std::cout << "opening file: " << fullRootFileName << std::endl;
      rootFiles.push_back(new TFile(fullRootFileName.c_str(), "READ"));
      if(!(rootFiles.at(i))->IsOpen()) exit(-1);
      
      // Draw::get tree
      TTree* tree;
      char treeName[50];
      sprintf(treeName, "ntu_%d", step);
      rootFiles.at(i) -> GetObject(treeName, tree);
      trees.push_back(tree);
      
      if(m_dataFlag[vecIt->first] == 1)
	{
	  WprimeVariables vars;
	  ClearWprimeVariables(vars);
	  SetBranchAddresses(vars, tree);

	  for (int entry = 0; entry < tree->GetEntries(); ++entry)
	    {
	      tree -> GetEntry(entry);
	      
	      if (vars.eleMet_mt < thr) continue;
	      if(!outFile) continue;
	      
	      (*outFile) << std::fixed << std::setprecision(5)
			 << "mt = " << std::setw(6) << vars.eleMet_mt << std::setw(18)
			 << "et/met = " << std::setw(6) << vars.ele.Et()/vars.met.Et() << std::setw(18)
			 << "eleMet dPhi = " << std::setw(6) << vars.eleMet_Dphi << std::setw(18)

			 << "metEt = " << std::setw(6) << vars.met.Et() << std::setw(18)
			 << "metPhi = " << std::setw(6) << vars.met.Phi() << std::setw(18)

			 << "eleEt = " << std::setw(6) << vars.ele.Et() << std::setw(18)
			 << "eleEta = " << std::setw(6) << vars.ele.eta() << std::setw(18)
			 << "elePhi = " << std::setw(6) << vars.ele.phi() << std::setw(18)
			 << "ele e/p = " << std::setw(6) << vars.ele_EOverP << std::setw(18)
			 << "ele h/e = " << std::setw(6) << vars.ele_HOverE << std::setw(18)
			 << "eleCharge = " << std::setw(6) << vars.ele_charge << std::setw(18)
			 << "ele dEta = " << std::setw(6) << vars.ele_DetaIn << std::setw(18)
			 << "ele dPhi = " << std::setw(6) << vars.ele_DphiIn << std::setw(18)
			 << "ele sIetaIeta = " << std::setw(6) << vars.ele_sigmaIetaIeta << std::setw(18)
			 << "ele EMiso = " << std::setw(6) << vars.ele_emIso << std::setw(18)
			 << "ele hadIso = " << std::setw(6) << vars.ele_hadIso << std::setw(18)
			 << "ele trkIso = " << std::setw(6) << vars.ele_tkIso << std::setw(18)

			 << "ele seedE = " << std::setw(6) << vars.ele_eSeed << std::setw(18)
			 << "ele seedTime = " << std::setw(6) << vars.ele_timeSeed << std::setw(18)
			 << "ele seedFlag = " << std::setw(6) << vars.ele_flagSeed << std::setw(18)
			 << "ele swissCross = " << std::setw(6) << vars.ele_swissCrossSeed << std::setw(18)
		
			 << "runId = " << std::setw(6) << vars.runId << std::setw(18)
			 << "lumiId = " << std::setw(6) << vars.lumiId << std::setw(18)
			 << "eventId = " << std::setw(6) << vars.eventId
			 << std::endl << std::endl;	    
	    }
	}
      ++i;
    }
  if(outFile)
    outFile -> close();
}
