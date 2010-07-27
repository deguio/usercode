#include "WprimeAnalysis/WprimeENUAnalysis/plugins/eventsCounter.h"






//! ctor
eventsCounter::eventsCounter(const edm::ParameterSet& iConfig)
{
  std::string histoName_ = iConfig.getParameter<std::string> ("histoName");

  edm::Service<TFileService> fs;
  
  m_totalEvents = fs -> make<TH1F>(histoName_.c_str(), histoName_.c_str(), 1,  0., 1.);

}

// ----------------------------------------------------------------






//! dtor
eventsCounter::~eventsCounter()
{}

// ----------------------------------------------------------------







bool eventsCounter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  
  m_totalEvents -> Fill(0.5);
  
  return true;
}
