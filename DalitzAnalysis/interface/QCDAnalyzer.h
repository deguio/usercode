#ifndef QCDAnalyzer_h
#define QCDAnalyzer_h
// -*- C++ -*-
//
// Package:    QCDAnalyzer
// Class:      QCDAnalyzer
// 


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "TFile.h"
#include "TH1F.h"

class QCDAnalyzer : public edm::EDAnalyzer {
   public:
      explicit QCDAnalyzer(const edm::ParameterSet&);
      ~QCDAnalyzer();

   
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      
 private:

      std::string HepMC_;
      
      double nEle;

      //----------
      TH1F *ePt;
      TH1F *eP;
      TH1F *Egamma;
      TH1F *nPi;
      TH1F *PiDaughters;
      TH1F *InvMassee;
      TH1F *InvMassgee;

      //----------
      TH1F *eleEt;

};
#endif
