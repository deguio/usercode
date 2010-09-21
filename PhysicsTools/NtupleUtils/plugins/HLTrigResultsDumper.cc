/** \class HLTrigResultsDumper
 *
 *  \author Leonardo Di Matteo
 *
 */

#include "PhysicsTools/NtupleUtils/plugins/HLTrigResultsDumper.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iomanip>


using namespace std;
using namespace edm;

//
// constructors and destructor
//
HLTrigResultsDumper::HLTrigResultsDumper(const edm::ParameterSet& iConfig) :
  hlTriggerResults_ (iConfig.getParameter<edm::InputTag> ("HLTriggerResults")),
  HLTPaths_ (iConfig.getParameter<std::vector<std::string> >("HLTPaths")),
  hltConfig_()
{
  Service<TFileService> fs ;
  HLTree = fs->make<TTree> ("HLTree","HLTree");
  HLTree -> Branch("HLTnpaths", &HLTnpaths, "HLTnpaths/i");
  HLTree -> Branch("HLTwasrun", HLTwasrun, "HLTwasrun[HLTnpaths]/i");
  HLTree -> Branch("HLTaccept", HLTaccept, "HLTaccept[HLTnpaths]/i");
  HLTree -> Branch("HLTerror", HLTerror, "HLTerror[HLTnpaths]/i");

  NameHLT = fs->make<TTree> ("NameHLT","NameHLT");
  NameHLT -> Branch("HLTTag_names",&HLTPaths_);
  NameHLT -> Fill();
}

HLTrigResultsDumper::~HLTrigResultsDumper()
{ }

//
// member functions
//

void
HLTrigResultsDumper::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;
  
  bool changed (true);
  if (hltConfig_.init(iRun,iSetup,hlTriggerResults_.process(),changed)) {
    if (changed) {
      // const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTR);
      hlNames_=hltConfig_.triggerNames();
    }
  } 
  else
    {
      // dump previous
      hlNames_.clear();
    }
  return;
}

    
// ------------ method called to produce the data  ------------
void
HLTrigResultsDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // accumulation of statistics event by event

  using namespace std;
  using namespace edm;

  //initialize branches
  HLTnpaths = 0;
  for (int kk = 0; kk < MAXPATH; ++kk)
    {
      HLTwasrun[kk] = 0;
      HLTaccept[kk] = 0;
      HLTerror[kk] = 0;
    }

  // get hold of TriggerResults
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResults_,HLTR);

  const unsigned int n(hlNames_.size());
  
  // decision for each HLTPaths
  HLTnpaths = (int) HLTPaths_.size();
  for ( int iPath = 0; iPath < (int) HLTPaths_.size(); iPath++ ) {
    //Find the index of the HLTPath under inspection
    int thePathIndex = -1;
    for (unsigned int i = 0; i < n; ++i) {
      if ( HLTPaths_[iPath] == hlNames_[i] ) thePathIndex = i;
    }
    
    //fill HLT info only for the HLT path which are in the HLT menu
    if ( thePathIndex > -1 ){
      if ( HLTR->wasrun((unsigned int)thePathIndex) ) HLTwasrun[iPath] = 1;
      if ( HLTR->accept((unsigned int)thePathIndex) ) HLTaccept[iPath] = 1;
      if ( HLTR->error((unsigned int)thePathIndex) )  HLTerror[iPath] = 1;
    }
  }
  
  //Fill HLT tree
  HLTree -> Fill();

  return;

}

void
HLTrigResultsDumper::endJob()
{ }

// ----------------------------------------------------------------

//define this as a plug-in
DEFINE_FWK_MODULE(HLTrigResultsDumper);

// ----------------------------------------------------------------
