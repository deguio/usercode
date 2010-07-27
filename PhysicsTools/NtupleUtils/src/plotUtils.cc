#include "PhysicsTools/NtupleUtils/interface/plotUtils.h"
#include "PhysicsTools/NtupleUtils/interface/setTDRStyle.h"






double MyGetMinimum(const TH1F* histo, const double& minval, int binMin, int binMax)
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






drawTStack::drawTStack(const std::string& inputDir,
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
  
  std::string listFullFileName = inputDir+listFileName;
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
    double mH;
    double crossSection;
    
    listFile >> sample >> sumName >> dataFlag >> mH >> crossSection;

    if(sample.size() == 0)
      continue;
    if(sample.at(0) == '#')
      continue;
    
    std::cout << sample << "\t"
              << sumName << "\t"
              << dataFlag << "\t"
              << mH << "\t"
              << crossSection << "\t" 
              << std::endl;
    
    std::pair<std::string, std::string> dummyPair(sample, sumName);
    m_list.push_back(dummyPair);
    m_dataFlag[sample] = dataFlag;
    m_mH[sample] = mH;
    m_crossSection[sample] = crossSection;
    
  }
  
  listFile.close();
  std::cout << ">>>plotUtils::Read " << m_list.size() << " samples" << std::endl;
  std::cout << ">>>plotUtils::Closing file " << listFullFileName << "\n" << std::endl;
}






drawTStack::~drawTStack()
{
  delete c1;
}






void drawTStack::Draw(const std::vector<std::string>& histoNames, const std::string& mode,
                      const float& lumi, const int& step,
                      const int& rebin, const bool& logy)
{ 
  std::cout << "\n>>>plotUtils::Drawing histogram " << histoNames.at(0);
  for(unsigned int j = 1; j < histoNames.size(); ++j)
    std::cout << " + " << histoNames.at(j);
  std::cout << std::endl;
  
  
  
  
  
  
  //---------------------------------------------
  // define the map with summed cross sections
  //---------------------------------------------
  
  std::map<std::string, double> crossSection_summed;
  std::map<std::string, bool> isFirstSample_summed;
  std::map<std::string, int> dataFlag_summed;
  std::map<std::string, TH1F*> histo_summed;
  
  //initialize summed vectors
  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    crossSection_summed[vecIt->second] = 0.;
    isFirstSample_summed[vecIt->second] = true;
    dataFlag_summed[vecIt->second] = false;
    histo_summed[vecIt->second] = NULL;
  }
  
  
  
  //---------------------------------------------
  // loop over all the samples and fill the stack
  //---------------------------------------------
  int binMin = -1;
  int binMax = -1;
  
  std::vector<TFile*> rootFiles;
  int i = 0;
  for(std::vector<std::pair<std::string, std::string> >::const_iterator vecIt = m_list.begin();
      vecIt != m_list.end(); ++vecIt)
  {
    // open root file
    std::string fullRootFileName = m_inputDir+vecIt->first+"/"+m_baseRootFileName+".root";
    //std::cout << "opening file: " << fullRootFileName << std::endl;
    rootFiles.push_back(new TFile(fullRootFileName.c_str(), "READ"));
    if(!(rootFiles.at(i))->IsOpen()) exit(-1);
    
    
    
    // get histograms
    TH1F* histo = NULL;
    for(unsigned int j = 0; j < histoNames.size(); ++j)
    {
      char buffer[10];
      sprintf(buffer, "%d", step);
      std::string fullHistoName = histoNames.at(j) + "/" + "h_" + buffer + "_" + histoNames.at(j);
      //std::cout << "getting histogram " << fullHistoName << std::endl;
      
      TH1F* tempHisto = NULL;
      rootFiles.at(i) -> GetObject(fullHistoName.c_str(), tempHisto);
      if(tempHisto == NULL)
      {
        std::cout << ">>>plotUtils::Error in getting object " << fullHistoName << std::endl;
        exit(-1);
      }
      
      if(histo != NULL) histo -> Add(tempHisto);
      else histo = tempHisto;
    }
    
    std::string eventsHistoName = "events";
    //std::cout << "getting histogram " << eventsHistoName << std::endl;
    
    TH1F* eventsHisto = NULL;
    rootFiles.at(i) -> GetObject(eventsHistoName.c_str(), eventsHisto);
    if(eventsHisto == NULL)
    {
      std::cout << ">>>plotUtils::Error in getting object " << eventsHistoName << std::endl;
      exit(-1);
    }
    
    
    
    // scale histograms normalizing to lumi (1. pb^-1)
    // if data do not apply any scale factor
    histo -> Sumw2();
    histo -> Rebin(rebin);
    dataFlag_summed[vecIt->second] = m_dataFlag[vecIt->first];
    
    if( (histo->GetEntries() > 0.) && (m_dataFlag[vecIt->first] != 1) )
    {
      histo -> Scale(m_crossSection[vecIt->first]/eventsHisto->GetBinContent(1));
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
  
  TLegend legend(0.68, 0.78, 0.99, 0.99);
  legend.SetFillColor(kWhite);
  
  double globalMaximum = 0.;
  double globalMinimum = 999999999999.;
  
  
  
  
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
      
      if( (isFirstSample == false) && (dataFlag_summed[mapIt->first] != 1) )
      {
        globalGlobalHisto -> Add(globalHisto);
      }
      if( (isFirstSample == true) && (dataFlag_summed[mapIt->first] != 1) )
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
      
      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum(globalHisto, 1.E-15, binMin, binMax);
      
      histoData = globalHisto;
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "P");
    }
    
    if( (mode == "eventsScaled") && (dataFlag_summed[mapIt->first] != 1) )
    {
      globalHisto -> Scale(1. * lumi);
      globalHisto -> SetLineColor(getColor(i+1));
      globalHisto -> SetFillColor(getColor(i+1));
      globalHisto -> SetFillStyle(3003);
      globalHisto -> SetLineWidth(2);
      
      hs -> Add(globalHisto);
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "F");
    }
    
    
    
    
    if( (mode == "sameAreaNoStack") )
    {
      globalHisto -> Scale(1./globalHisto->Integral(1, globalHisto->GetNbinsX()));
      globalHisto -> SetLineColor(getColor(i));
      globalHisto -> SetLineStyle(i+1);
      globalHisto -> SetLineWidth(4);

      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum(globalHisto, 1.E-15, binMin, binMax);
      
      hs -> Add(globalHisto);      
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "L");      
    }
    
    
    
    
    if( (mode == "sameAreaStack") && (dataFlag_summed[mapIt->first] == 1) )
    {

      //      std::cout << "FEDE >>> dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()) = " 
      //		<< dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()) << std::endl;

      globalHisto -> Scale(dataGlobalGlobalIntegral/globalHisto->Integral(1, globalHisto->GetNbinsX()));
      globalHisto -> SetMarkerStyle(20);
      
      if(globalHisto->GetMaximum() > globalMaximum)
        globalMaximum = globalHisto -> GetMaximum();
      if(MyGetMinimum(globalHisto, 1.E-15, binMin, binMax) < globalMinimum)
        globalMinimum = MyGetMinimum(globalHisto, 1.E-15, binMin, binMax);
      
      histoData = globalHisto;
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "P");      
    }
    
    if( (mode == "sameAreaStack") && (dataFlag_summed[mapIt->first] != 1) )
    {
      globalHisto -> Scale(dataGlobalGlobalIntegral/globalGlobalIntegral);
      globalHisto -> SetLineColor(getColor(i));
      globalHisto -> SetFillColor(getColor(i));
      globalHisto -> SetFillStyle(3003);
      globalHisto -> SetLineWidth(2);
      
      hs -> Add(globalHisto);      
      legend.AddEntry(globalHisto, (mapIt->first).c_str(), "F");      
    }
    
    
    
    
    if( (isFirstSample == false) && (dataFlag_summed[mapIt->first] != 1) )
    {
      globalGlobalHisto -> Add(globalHisto);
    }
    if( (isFirstSample == true) && (dataFlag_summed[mapIt->first] != 1) )
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
    if(MyGetMinimum(globalGlobalHisto, 1.E-15, binMin, binMax) < globalMinimum)
      globalMinimum = MyGetMinimum(globalGlobalHisto, 1.E-15, binMin, binMax);
  }
  
  
  
  
  // draw the stack and save file
  TCanvas* c1 = new TCanvas();
  c1 -> cd();
  c1 -> SetGridx();
  c1 -> SetGridy();
  if(logy) c1 -> SetLogy();
  
  
  
  if( mode == "eventsScaled" )
  {
    hs -> Draw("HISTO");
    if(histoData != NULL) histoData -> Draw("P,same");
    
    char buffer[50];
    sprintf(buffer, "events / %.3f pb^{-1}", lumi);
    hs->GetYaxis()->SetTitle(buffer);
  }

  if( mode == "sameAreaNoStack" )
  {
    hs -> Draw("nostack,HISTO");
    
    hs->GetYaxis()->SetTitle("event fraction");
  }
  
  if( mode == "sameAreaStack" )
  {
    hs -> Draw("HISTO");
    if(histoData != NULL) histoData -> Draw("P,same");
    
    hs->GetYaxis()->SetTitle("event fraction");
  }
  
  
  
  if(m_drawLegend == true)
    legend.Draw("same");
  
  
  
  
  
  
  // set x-axis properties
  std::string fullTitle = histoNames.at(0);
  for(unsigned int j = 1; j < histoNames.size(); ++j)
    fullTitle += "+" + histoNames.at(j); 
  hs->GetXaxis()->SetTitle(fullTitle.c_str());
  if(m_xAxisTitle) hs->GetXaxis()->SetTitle(m_xTitle.c_str());
  hs -> GetXaxis() -> SetTitleSize(0.04);
  hs -> GetXaxis() -> SetLabelSize(0.0); //fede
  hs -> GetXaxis() -> SetTitleOffset(1.25);

  if(m_xAxisRange) hs->GetXaxis()->SetRangeUser(m_xRangeMin, m_xRangeMax);
  
  
  
  // set y-axis properties
  if(m_yAxisTitle) hs->GetYaxis()->SetTitle(m_yTitle.c_str());    
  hs -> GetYaxis() -> SetTitleSize(0.04);
  hs -> GetYaxis() -> SetLabelSize(0.0); //fede
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


void drawTStack::DrawEvents(const std::string& mode, const float& lumi, const int& step, const bool& logy)
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
      std::cout << ">>>plotUtils::Error in getting object " << fullHistoName << std::endl;
      exit(-1);
    }
    
    crossSection_summed[vecIt->second] += m_crossSection[vecIt->first];
    if(m_mH[vecIt->first] > 0)
      isSignal_summed[vecIt->second] = true;
    
    
    if( (mode == "eventsScaled") ||
        (mode == "eventsScaledStack") ||
        (mode == "efficiencies") ||
        (mode == "efficienciesRelative") ||
        (mode == "significance") )
    {
      if( m_dataFlag[vecIt->first] != 1 )
        histo -> Scale(lumi*m_crossSection[vecIt->first]/histo->GetBinContent(1));
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
    outFile = new std::ofstream((m_outputDir+"eventsScaledStack.txt").c_str(), std::ios::out);
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
      int totalEvents = globalHisto -> GetBinContent(1);
      for(int bin = 1; bin <= globalHisto->GetNbinsX(); ++bin)
        globalHisto->SetBinContent(bin, globalHisto->GetBinContent(bin)/totalEvents);
    }
    
    if(mode == "efficienciesRelative")
    {
      std::map<int, int > totalEvents;
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
    
    
    globalHisto -> SetLineColor(getColor(i));
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
                                                  << std::setprecision(3)
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
  }
  else
  {
    hs -> Draw("HISTO");
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
    sprintf(buffer, "events / %.3f pb^{-1}", lumi);
    hs->GetYaxis()->SetTitle(buffer);
    stepHisto->GetYaxis()->SetTitle(buffer);
    
    if(histoData != NULL)
    {
      histoData -> SetMarkerStyle(20);  
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
    sprintf(buffer, "events / %.3f pb^{-1}", lumi);
    hs->GetYaxis()->SetTitle(buffer);
    stepHisto->GetYaxis()->SetTitle(buffer);
    
    if(histoData != NULL)
    {
      histoData -> SetMarkerStyle(20);
      
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






void drawTStack::SetXaxisRange(const double& xMin, const double& xMax)
{
  m_xAxisRange = true;
  m_xRangeMin = xMin;
  m_xRangeMax = xMax;
}

void drawTStack::SetXaxisTitle(const std::string& xTitle)
{
  m_xAxisTitle = true;
  m_xTitle = xTitle;
}



void drawTStack::SetYaxisRange(const double& yMin, const double& yMax)
{
  m_yAxisRange = true;
  m_yRangeMin = yMin;
  m_yRangeMax = yMax;
}

void drawTStack::SetYaxisTitle(const std::string& yTitle)
{
  m_yAxisTitle = true;
  m_yTitle = yTitle;
}



void drawTStack::SetDrawLegend(const bool& drawLegend)
{
  m_drawLegend = drawLegend;
}

void drawTStack::SetXLegend(const double& xLow, const double& xHigh)
{
  m_xLowLegend = xLow;
  m_xHighLegend = xHigh;
}

void drawTStack::SetYLegend(const double& yLow, const double& yHigh)
{
  m_yLowLegend = yLow;
  m_yHighLegend = yHigh;
}
