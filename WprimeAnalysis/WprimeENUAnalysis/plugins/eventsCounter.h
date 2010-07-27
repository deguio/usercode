#ifndef EVENTSCOUNTER_h
#define EVENTSCOUNTER_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"






class eventsCounter : public edm::EDFilter
{
  
 public:
  
  //! ctor
  explicit eventsCounter (const edm::ParameterSet&);
  
  //! dtor 
  ~eventsCounter();
  
  
  
 private:
  
  //! the actual filter method 
  bool filter(edm::Event&, const edm::EventSetup&);
  
  
  
 private:
  
  TH1F* m_totalEvents;

  
};

#endif
