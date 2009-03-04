#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "Calibration/Tools/interface/calibXMLwriter.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibTools.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/CaloMiscalibMapEcal.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalBarrel.h"
#include "CalibCalorimetry/CaloMiscalibTools/interface/MiscalibReaderFromXMLEcalEndcap.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

//#include "Calibration/EcalAlCaRecoProducers/interface/trivialParser.h"
//#include "Calibration/EcalAlCaRecoProducers/bin/trivialParser.h"

//---- new for IOV Iterator ----
#include "CondCore/Utilities/interface/CondIter.h"
#include "TSystem.h"

//---------------

#include "CaloOnlineTools/EcalTools/interface/EcalCosmicsTreeContent.h"
#include "CaloOnlineTools/EcalTools/interface/EcalCosmicsTreeUtils.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CRUZET/Calibration/interface/CRUtils.h"

#include <iostream>
#include <string>
#include <boost/foreach.hpp>

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TVector3.h"

#define PHI_MIN -3.1416
#define PHI_MAX +3.1416
#define PHI_BIN 360

#define ETA_MIN -1.5
#define ETA_MAX +1.5
#define ETA_BIN 170

#define IPHI_MIN 1.
#define IPHI_MAX +361.
#define IPHI_BIN 360

#define IETA_MIN -85.
#define IETA_MAX +86.
#define IETA_BIN 171

#define X_BIN 500
#define X_MIN -800
#define X_MAX +800

#define Y_BIN 500
#define Y_MIN -800
#define Y_MAX +800

#define P_BIN 500
#define P_MIN 0
#define P_MAX 1000

#define E_BIN 1000
#define E_MIN 0
#define E_MAX 250

#define PI 3.141592654



int main (int argc, char** argv)
{
  //------GET COEFF--------
  int EBetaStart = -85 ;
  int EBetaEnd = 86 ;
  int EBphiStart = 1 ;
  int EBphiEnd = 361 ;

  
  std::string NameDBOracle1 = "oracle://cms_orcoff_prod/CMS_COND_21X_ECAL";   //GAIN 50
  std::string TagDBOracle1 = "EcalIntercalibConstants_CosmGain50";  
  
  std::string NameDBOracle2 = "oracle://cms_orcoff_prod/CMS_COND_21X_ECAL";   //GAIN 200
  std::string TagDBOracle2 = "EcalIntercalibConstants_CosmGain200";  
  

  
  std::string Command2LineStr1 = "cmscond_export_iov -s " + NameDBOracle1 + " -d sqlite_file:Uno.db -D CondFormatsEcalObjects -t " + TagDBOracle1 + " -P /afs/cern.ch/cms/DB/conddb/";
  
  std::string Command2LineStr2 = "cmscond_export_iov -s " + NameDBOracle2 + " -d sqlite_file:Due.db -D CondFormatsEcalObjects -t " + TagDBOracle2 + " -P /afs/cern.ch/cms/DB/conddb/";
  
  std::cout << Command2LineStr1 << std::endl;
  std::cout << Command2LineStr2 << std::endl;
  
  
  gSystem->Exec(Command2LineStr1.c_str());
  gSystem->Exec(Command2LineStr2.c_str());
  
  
  std::string NameDB;
  std::string FileData;
  
  
  //----------------------------------
  //---- First Database Analyzed -----
  //----------------------------------
  NameDB = "sqlite_file:Uno.db";
  FileData = TagDBOracle1;
  CondIter<EcalIntercalibConstants> Iterator1;
  Iterator1.create(NameDB,FileData);
  
  //-----------------------------------
  //---- Second Database Analyzed -----
  //-----------------------------------
  NameDB = "sqlite_file:Due.db";  
  FileData = TagDBOracle2;
  CondIter<EcalIntercalibConstants> Iterator2;
  Iterator2.create(NameDB,FileData);

  //---------------------------------------------------------------------
  
  
  //-------------------------------------------------
  //---- Ottengo Mappe da entrambi gli Iterators ----
  //-------------------------------------------------

   
  const EcalIntercalibConstants* EBconstants1;
  EBconstants1 = Iterator1.next();
  EcalIntercalibConstantMap iEBcalibMap_1 = EBconstants1->getMap () ;
 
  
  
  const EcalIntercalibConstants* EBconstants2;
  EBconstants2 = Iterator2.next();
  EcalIntercalibConstantMap iEBcalibMap_2 = EBconstants2->getMap () ;

//   //LOOP OVER XTALS
//   //loop over eta
//   double BarrelAve_1 = 0.;
//   double BarrelAve_2 = 0.;
//   double BarrelNum = 0.;
//   for (int ieta = EBetaStart ; ieta < EBetaEnd ; ++ieta)
//   {
//       double phiSum = 0. ; 
//       double phiSumSq = 0. ; 
//       double N = 0. ;
//       //loop over phi
//       for (int iphi = EBphiStart ; iphi <= EBphiEnd ; ++iphi)
//       {
//           if (!EBDetId::validDetId (ieta,iphi)) continue ;
//           EBDetId det = EBDetId (ieta,iphi,EBDetId::ETAPHIMODE) ;

// 	  double coeff_1 = *(iEBcalibMap_1.find (det.rawId ()));
// 	  double coeff_2 = *(iEBcalibMap_2.find (det.rawId ()));


// 	  BarrelAve_1 += coeff_1;
// 	  BarrelAve_2 += coeff_2;
// 	  BarrelNum += 1.;
	  
//       } //loop over phi
//   } //loop over eta


  //-------END GET COEFF-------

  std::string fileName (argv[1]) ;
  boost::shared_ptr<edm::ProcessDesc> processDesc = edm::readConfigFile (fileName) ;
  boost::shared_ptr<edm::ParameterSet> parameterSet = processDesc->getProcessPSet () ;
  std::cout << parameterSet->dump () << std::endl ; // for testing
  
  edm::ParameterSet subPSetSelections =  parameterSet -> getParameter<edm::ParameterSet>("selections");
  int superClusterIPhiMIN = subPSetSelections.getUntrackedParameter<int>("superClusterIPhiMIN", 1);
  int superClusterIPhiMAX = subPSetSelections.getUntrackedParameter<int>("superClusterIPhiMAX", 361);
  int superClusterIEtaMIN = subPSetSelections.getUntrackedParameter<int>("superClusterIEtaMIN", -85);
  int superClusterIEtaMAX = subPSetSelections.getUntrackedParameter<int>("superClusterIEtaMAX", 86);
  std::string outputRootName = subPSetSelections.getUntrackedParameter<std::string> ("outputRootName","test.root") ;

  edm::ParameterSet subPSetInput = parameterSet->getParameter<edm::ParameterSet> ("inputNtuples") ;
  std::vector<std::string> inputFiles = subPSetInput.getParameter<std::vector<std::string> > ("inputFiles") ;
  
  TChain *chain = new TChain ("EcalCosmicsAnalysis") ;
  EcalCosmicsTreeContent treeVars ; 
  setBranchAddresses (chain, treeVars) ;
  
  
  // input files
  for (std::vector<std::string>::const_iterator listIt = inputFiles.begin () ;
       listIt != inputFiles.end () ; ++listIt)
    {
      std::cout << *listIt << " " << std::endl ;
      chain->Add (listIt->c_str ()) ;
    }
  
  int nEntries = chain->GetEntries () ;
  std::cout << "FOUND " << nEntries << " ENTRIES\n" ;    
  // input files
  
  
  // output file
  //  std::string outputRootName = outputRootName ;
  // output file  
  
  
  
  // output distributions
  TH2F SCoccupancy ("SCoccupancy", "SCoccupancy", PHI_BIN, PHI_MIN, PHI_MAX, ETA_BIN, ETA_MIN, ETA_MAX) ;
  TH1F SCenergy ("SCenergy", "SCenergy", E_BIN, E_MIN, E_MAX) ;
  TProfile2D SCenergy_ETAvsPHI ("SCenergy_ETAvsPHI", "SCenergy_ETAvsPHI", PHI_BIN, PHI_MIN, PHI_MAX, ETA_BIN, ETA_MIN, ETA_MAX) ;
  
  TProfile2D XtalEnergy_ETAvsPHI ("XtalEnergy_ETAvsPHI", "XtalEnergy_ETAvsPHI", IPHI_BIN, IPHI_MIN, IPHI_MAX, IETA_BIN, IETA_MIN, IETA_MAX);
  TH2F rapportoMap ("rapportoMap", "rapportoMap", IPHI_BIN, IPHI_MIN, IPHI_MAX, IETA_BIN, IETA_MIN, IETA_MAX);
  TH2F test ("test", "test", IPHI_BIN, IPHI_MIN, IPHI_MAX, IETA_BIN, IETA_MIN, IETA_MAX);
  TH1F rapporto_distr ("rapporto_distr", "rapporto_distr", 100000, 0., 100.);

  TH2F* EonESCvsCry = new TH2F ("EonESCvsCry", "EonESCvsCry", 61200, 0., 61200., 1000, 0., 10.);
  TProfile* EonESCvsCry_profile = new TProfile ("EonESCvsCry_profile", "EonESCvsCry_profile", 61200, 0., 61200.);
  TProfile* rapporto_profile = new TProfile ("rapporto_profile", "rapporto_profile", 61200, 0., 61200.);
  //--------------------------
  TH1F coeffDistr ("coeffDistr", "coeffDistr", 1000, 0., 2.);
  TH1F dEondxDistr ("dEondxDistr","dEondxDistr", 1000, 0., 0.1);

  TH2F dEondx_VS_ICmean ("dEondx_VS_ICmean","dEondx_VS_ICmean", 5000, 0., 5., 1000, 0., 0.1);
  TProfile* dEondx_VS_ICmean_pfx = new TProfile ("dEondx_VS_ICmean_pfx","dEondx_VS_ICmean_pfx", 5000, 0., 5.);
  TProfile* dEondx_VS_ICmean_pfy = new TProfile ("dEondx_VS_ICmean_pfy","dEondx_VS_ICmean_pfy", 1000, 0., 0.1);
  // output distributions
  

  
  // loop over entries
  for (int entry = 0; entry < nEntries; ++entry)
    //for (int entry = 0; entry < 20000; ++entry)
    {
      if ((entry % 100000) == 0)
        std::cout << "Reading entry " << entry << std::endl;
      
      chain->GetEntry (entry) ;
      
      // association MU-SC
      std::vector<ect::association> associations ;
      ect::fillAssocVector (associations, treeVars) ;
      ect::selectOnDR (associations, treeVars, 0.1) ;
      
      // numAssociations.Fill(associations.size()) ;
      
      
      
      //loop on associations vector
      for (unsigned int i = 0 ; i < associations.size () ; ++i)
	{
	  
	  int MUindex = associations.at (i).first  ;
	  int SCindex = associations.at (i).second ;

	  float muonP = treeVars.muonP[MUindex];
	  float muonPt = treeVars.muonPt[MUindex];
	  float muond0 = treeVars.muond0[MUindex];
	  float muondz = treeVars.muondz[MUindex];
	  float muonChi2 = treeVars.muonNChi2[MUindex];
	  float muonPhi = treeVars.muonPhi[MUindex];
	  float muonEta = treeVars.muonEta[MUindex];
	  float muonTkLengthInEcal = treeVars.muonTkLengthInEcalDetail[MUindex] ;
	  float muonTkLengthInEcalCurved = treeVars.muonTkLengthInEcalDetailCurved[MUindex] ;

	  float superClusterRawEnergy = treeVars.superClusterRawEnergy[SCindex] * 0.97 ;
	  float superClusterEta = treeVars.superClusterEta[SCindex] ;
	  float superClusterPhi = treeVars.superClusterPhi[SCindex] ;
	  int nClustersInSuperCluster = treeVars.nClustersInSuperCluster[SCindex] ;
	  int nXtalsInSuperCluster = treeVars.nXtalsInSuperCluster[SCindex] ;
	  int xtalIndexInSuperCluster = treeVars.xtalIndexInSuperCluster[SCindex] ;
	  

	  //SELECTIONS

	  //tutti gli evt nello stesso bin di momento
	  if (muonP <= 20. || muonP >= 40.) continue;
	  
	  //evt nel SM di interesse
	  double superClusterEtaMIN;
	  double superClusterEtaMAX;
	  double superClusterPhiMIN;
	  double superClusterPhiMAX;
	  
	  if (superClusterIEtaMIN >= 0) superClusterEtaMIN = (superClusterIEtaMIN - 0.5) * 0.0175;
	  else superClusterEtaMIN = (superClusterIEtaMIN + 0.5) * 0.0175;
	  superClusterPhiMIN = (superClusterIPhiMIN - 10) * 0.0175;
	  if(superClusterPhiMIN > PI) superClusterPhiMIN = - 2*PI + superClusterPhiMIN;        

	  if (superClusterIEtaMAX >= 0) superClusterEtaMAX = (superClusterIEtaMAX - 0.5) * 0.0175;
	  else superClusterEtaMAX = (superClusterIEtaMAX + 0.5) * 0.0175;
	  superClusterPhiMAX = (superClusterIPhiMAX - 10) * 0.0175;
	  if(superClusterPhiMAX > PI) superClusterPhiMAX = - 2*PI + superClusterPhiMAX;        

	  if ( (superClusterPhi < superClusterPhiMIN) ||
	       (superClusterPhi > superClusterPhiMAX) ) continue;
	  if ( (superClusterEta < superClusterEtaMIN) ||
	       (superClusterEta > superClusterEtaMAX) ) continue;

	  //END SELECTION

	  //test plots
	  SCoccupancy.Fill(superClusterPhi, superClusterEta);

	  //end test plot

	  if (muonTkLengthInEcal <= 0.) continue;
	  float dEdx = superClusterRawEnergy / muonTkLengthInEcal; // / 8.28 * 1000. ;
	  
	  //singoli cristalli
	  float ICmean = 0.;
	  int N = 0;
	  for (int XTLindex = treeVars.xtalIndexInSuperCluster[SCindex] ;
               XTLindex < treeVars.xtalIndexInSuperCluster[SCindex] + treeVars.nXtalsInSuperCluster[SCindex] ;
               ++XTLindex)
            {
	      EBDetId dummy = EBDetId::unhashIndex (treeVars.xtalHashedIndex[XTLindex]) ;
	      	    
	      double coeff_2 = *(iEBcalibMap_2.find (dummy.rawId ()));

	      coeffDistr.Fill(coeff_2);

	      ICmean += coeff_2;
	      ++N;
	      
	    }
	  //end singoli cristalli
	  
	  //if(N <= 0 || dEdx <= 0.) continue;
	  dEondxDistr.Fill(dEdx);
	  dEondx_VS_ICmean.Fill(ICmean/N, dEdx);
	  dEondx_VS_ICmean_pfx->Fill(ICmean/N, dEdx);
	  dEondx_VS_ICmean_pfy->Fill(dEdx, ICmean/N);

	 	  
	}
      
      
    } // loop over entries

//   //----BEGIN lavoro sulla mappa di energie-----
//   for(int ieta = (int)IETA_MIN; ieta <= (int)IETA_MAX; ++ieta)
//     for(int iphi = (int)IPHI_MIN; iphi <= (int)IPHI_MAX; ++iphi)
//       {
// 	if (!EBDetId::validDetId (ieta,iphi)) continue ;	
	
// 	double mean8 = 0.;
// 	int N = 0;

// 	//std::cerr << "test coord: ieta, iphi = " << ieta << " " << iphi << std::endl;

// 	for(int iieta = ieta -1; iieta <= ieta +1; ++iieta)
// 	  for(int iiphi = iphi -1; iiphi <= iphi +1; ++iiphi)

// 	    {
// 	      if (iieta == ieta && iiphi == iphi) continue;    //Xtal Centrale
// 	      if (!EBDetId::validDetId (iieta,iiphi)) continue ;

// 	      //std::cerr << "test coord: iieta, iiphi = " << iieta << " " << iiphi << std::endl;


// 	      mean8 += XtalEnergy_ETAvsPHI.GetBinContent(iiphi, iieta);
// 	      ++N;
// 	    }
	
	
// 	double rapporto = XtalEnergy_ETAvsPHI.GetBinContent(iphi, ieta) / mean8 * N;
// 	//rapportoMap.Fill(iphi, ieta, rapporto);
	
// 	//std::cerr << "test coord: ieta, iphi = " << ieta << " " << iphi << std::endl;
// 	//std::cerr << "rapporto, mean8  = " << rapporto << " " << mean8 << std::endl;
// 	//std::cerr << XtalEnergy_ETAvsPHI.GetBinContent(iphi, ieta) << std::endl;
//       }
  


  //----END lavoro sulla mappa di energie-----
  
  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;  
  
  SCoccupancy.Write () ;
  SCenergy.Write () ;
  SCenergy_ETAvsPHI.Write () ;

  XtalEnergy_ETAvsPHI.Write () ;
  rapportoMap.Write () ;
  test.Write () ;
  rapporto_distr.Write();

  EonESCvsCry->Write();
  EonESCvsCry_profile->Write();
  rapporto_profile->Write();
  //-----------------
  coeffDistr.Write();
  dEondxDistr.Write();
  dEondx_VS_ICmean.Write();
  dEondx_VS_ICmean_pfx->Write();
  dEondx_VS_ICmean_pfy->Write();
  
  saving.Close () ;
  //delete EonESCvsCry;
  //delete EonESCvsCry_profile;
  //delete rapporto_profile;
  
  return 0 ;
}


