#ifndef plotUtils_ntu_h
#define plotUtils_ntu_h

#include <cstdlib>
#include <string>
#include <map>
#include <utility>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sys/stat.h> 

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"



double MyGetMinimum_ntu(const TH1F* histo, const double& minval, int binMin=-1, int binMax=-1);



class drawTStack_ntu
{
 public:
    
  //! ctor
  drawTStack_ntu(const std::string& inputDir,
             const std::string& listFileName,
             const std::string& baseRootFileName,
             const std::string& outputDir);
  
  //! dtor
  ~drawTStack_ntu();
  
  
  //! methods
  void Draw(std::vector<std::string>& variableNames, const std::string& histoName,
            const std::string& mode,
            const float& lumi, const int& step,
            const int& nBins, const bool& logy,
            std::vector<std::string>* cut = NULL,
	    const bool& wantCumulative = false);
  void DrawEvents(const std::string& mode, const float& lumi, const int& step, const bool& logy);
  void SaveTStack(THStack* hs, const std::string& histoName, TH1F* histoData);
  void DrawTStackError(THStack* hs);
  
  void SetXaxisRange(const double& xMin, const double& xMax);
  void SetXaxisTitle(const std::string& xTitle);

  void SetYaxisRange(const double& yMin, const double& yMax);
  void SetYaxisTitle(const std::string& yTitle);

  void SetDrawLegend(const bool& drawLegend);
  void SetXLegend(const double& xLow, const double& xHigh);
  void SetYLegend(const double& yLow, const double& yHigh);
  
  void numbersForLimit(std::vector<std::string>& variableNames,
		       const float& lumi, const int& step,
		       std::vector<std::string>* cut = NULL,
		       const int& limitStep = 25.,
		       const double& xLow = 200., const double& xHigh = 1500.);

  void printMtAboveThr(const int& step,
		     const double& thr = 600.);

 private:
  
  std::string m_inputDir;
  std::string m_listFileName;
  std::string m_baseRootFileName;
  std::string m_outputDir;
  
  std::vector<std::pair<std::string, std::string> > m_list;
  std::map<std::string, int> m_dataFlag;
  std::map<std::string, int> m_isDD;
  std::map<std::string, double> m_mH;
  std::map<std::string, double> m_crossSection;
  std::map<std::string, int> m_color;
  
  bool m_xAxisRange;
  double m_xRangeMin;
  double m_xRangeMax;
  bool m_xAxisTitle;
  std::string m_xTitle;
  
  bool m_yAxisRange;
  double m_yRangeMin;
  double m_yRangeMax;
  bool m_yAxisTitle;
  std::string m_yTitle;
  
  bool m_drawLegend;
  double m_xLowLegend;
  double m_yLowLegend;
  double m_xHighLegend;
  double m_yHighLegend;
  
  TCanvas* c1;
};

#endif
