

#include "Dalitz/DalitzAnalysis/interface/QCDAnalyzer.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include <iostream>


using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

QCDAnalyzer::QCDAnalyzer(const edm::ParameterSet& iConfig)
{
  
  
  HepMC_ = iConfig.getUntrackedParameter<std::string>("HepMCLabel");
  
  
  
}

QCDAnalyzer::~QCDAnalyzer() {}

void QCDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  

  edm::Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  edm::Handle<reco::PixelMatchGsfElectron> electrons ;
  iEvent.getByLabel ("pixelMatchGsfElectrons",electrons) ;

  //------GEN PARTICLES-------

  int numPi = 0;
  for(GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) 
    {
      
      int nMo = p->numberOfMothers();
      int nDa = p->numberOfDaughters();
      
      int myMo=-99;
      int myDa=-99;
      
      //if( nMo == 1 ) myMo = p->mother(0)->pdgId();
      //if( nDa > 0 ) myDa = p->daughter(0)->pdgId();
      
      //DAUGHTERS-------------------------------------- 
      // search for e+e- couples with pi0 mother
      if (abs(p->pdgId()) == 11 && p->status()== 1 && p->mother(0)->pdgId() == 111) 
	{
	  ePt->Fill(p->pt());
	  
	}//end loop electrons from Pi
      
      //MOTHERS-----------------------------------------
      // search for pi0 (pdgId = 111)
      const reco::Candidate* ele = 0;
      const reco::Candidate* pos = 0;
      const reco::Candidate* gamma = 0;
      
      bool isEle = false;
      bool isPos = false;
      bool isGamma = false;
      
      if (p->pdgId() == 111 ) {
	for (int j=0; j<nDa; j++)
	  {
	    if(p->daughter(j)->pdgId() ==  11) 
	      {
		ele = p->daughter(j);
		isEle = true;
	      }
	    if(p->daughter(j)->pdgId() == -11) 
	      {
		pos = p->daughter(j);
		isPos = true;
	      }
	    if(p->daughter(j)->pdgId() == 22) 
	      {
		gamma = p->daughter(j);
		isGamma = true;
	      }
	    
	    PiDaughters->Fill(p->daughter(j)->pdgId());
	  }
	
	++numPi;
      }//end loop Pi
      
      //    massa inv e+ e-
      if(isEle == true && isPos == true && isGamma == true)
	{
	  double Minvee = (pos->p4() + ele->p4()).M();
	  double Minvgee = (gamma->p4() + pos->p4() + ele->p4()).M();
	  double Eg = gamma->energy();

	  InvMassee->Fill(Minvee);
	  InvMassgee->Fill(Minvgee);
	  Egamma->Fill(Eg);
	}
      
      
    } //end loop over particles
  nPi->Fill(numPi);
  
  
  //------RECO GSF ELECTRONS-------
  
  for(reco::PixelMatchGsfElectron::const_iterator ele = electrons->begin(); ele != electrons->end(); ++ele) 
    {
      eleEt->Fill(ele->et());
      
    }
  
  
  //   if(electrons->size() < 30 ){ nEle = electrons->size();}
  //   else {nEle = 30;}
  //   for(int i=0; i< nEle; i++)
  //     {
  //       eleEt->Fill(->et());
  //     }
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
QCDAnalyzer::beginJob(const edm::EventSetup&)
{ 
  edm::Service<TFileService> fs;

  // initialize histograms
  //------------------------
  ePt = fs->make<TH1F>("ePt","ePt",2000, 0., 200.);
  eP  = fs->make<TH1F>("eP","eP",2000, 0., 200.);
  Egamma =  fs->make<TH1F>("Egamma","Egamma",2000, 0., 200.);

  nPi = fs->make<TH1F>("nPi","nPi",200, 0., 200. );
  PiDaughters = fs->make<TH1F>("PiDaughters","PiDaughters",60, -30., 30. );
  
  InvMassee = fs->make<TH1F>("InvMassee","InvMassee",2000, 0., 20. );
  InvMassgee = fs->make<TH1F>("InvMassgee","InvMassgee",2000, 0., 20. );
  
  //------------------------
  eleEt = fs->make<TH1F>("eleEt","eleEt",2000, 0., 200. );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDAnalyzer::endJob() {

  
}


