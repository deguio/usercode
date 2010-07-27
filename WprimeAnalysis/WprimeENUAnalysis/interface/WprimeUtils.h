#ifndef WprimeUtils_h
#define WprimeUtils_h

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <utility>

#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TChain.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "TVectorF.h"


/** get the number of events from a list of files */
std::map<int, int> GetTotalEvents(const std::string& histoName, const std::string& inputFileList);


#endif
