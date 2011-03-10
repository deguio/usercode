//=====================+++AUTO-Reader+++============================
 std::vector<ROOT::Math::XYZTVector>* electrons = reader.Get4V("electrons");
 std::vector<ROOT::Math::XYZTVector>* electrons_SC = reader.Get4V("electrons_SC");
/*
 std::vector<float>* eidLoose = reader.GetFloat("eidLoose");
 std::vector<float>* eidRobustLoose = reader.GetFloat("eidRobustLoose");
 std::vector<float>* eidTight = reader.GetFloat("eidTight");
 std::vector<float>* eidRobustTight = reader.GetFloat("eidRobustTight");
*/

 std::vector<int>* electrons_mishits = reader.GetInt("electrons_mishits");
 std::vector<int>* electrons_nAmbiguousGsfTracks = reader.GetInt("electrons_nAmbiguousGsfTracks");
 std::vector<float>* electrons_eES = reader.GetFloat("electrons_eES");

 std::vector<float>* electrons_dxy_PV_noEle = reader.GetFloat("electrons_dxy_PV_noEle");
 std::vector<float>* electrons_dz_PV_noEle = reader.GetFloat("electrons_dz_PV_noEle");

 std::vector<ROOT::Math::XYZTVector>* muons = reader.Get4V("muons");
 std::vector<float>* muons_charge = reader.GetFloat("muons_charge");
 std::vector<float>* muons_dB = reader.GetFloat("muons_dB");
 std::vector<float>* muons_edB = reader.GetFloat("muons_edB");
 std::vector<float>* muons_dxy = reader.GetFloat("muons_dxy");
 std::vector<float>* muons_edxy = reader.GetFloat("muons_edxy");
 std::vector<float>* muons_dz = reader.GetFloat("muons_dz");
 std::vector<float>* muons_edz = reader.GetFloat("muons_edz");
 std::vector<float>* muons_dxy_PV = reader.GetFloat("muons_dxy_PV");
 std::vector<float>* muons_dz_PV = reader.GetFloat("muons_dz_PV");

 std::vector<float>* muonss_dxy_PV_noMuon = reader.GetFloat("muons_dxy_PV_noMuon");
 std::vector<float>* muons_dz_PV_noMuon = reader.GetFloat("muons_dz_PV_noMuon");

 std::vector<float>* muons_nTkIsoR03 = reader.GetFloat("muons_nTkIsoR03");
 std::vector<float>* muons_nTkIsoR05 = reader.GetFloat("muons_nTkIsoR05");
 std::vector<float>* muons_tkIsoR03 = reader.GetFloat("muons_tkIsoR03");
 std::vector<float>* muons_tkIsoR05 = reader.GetFloat("muons_tkIsoR05");
 std::vector<float>* muons_emIsoR03 = reader.GetFloat("muons_emIsoR03");
 std::vector<float>* muons_emIsoR05 = reader.GetFloat("muons_emIsoR05");
 std::vector<float>* muons_hadIsoR03 = reader.GetFloat("muons_hadIsoR03");
 std::vector<float>* muons_hadIsoR05 = reader.GetFloat("muons_hadIsoR05");
 std::vector<int>* muons_tracker = reader.GetInt("muons_tracker");
 std::vector<int>* muons_standalone = reader.GetInt("muons_standalone");
 std::vector<int>* muons_global = reader.GetInt("muons_global");
 std::vector<int>* muons_goodMuon = reader.GetInt("muons_goodMuon");
 std::vector<float>* muons_normalizedChi2 = reader.GetFloat("muons_normalizedChi2");
 std::vector<int>* muons_numberOfValidTrackerHits = reader.GetInt("muons_numberOfValidTrackerHits");
 std::vector<int>* muons_numberOfValidMuonHits = reader.GetInt("muons_numberOfValidMuonHits");

 std::vector<ROOT::Math::XYZTVector>* photons = reader.Get4V("photons");
// std::vector<float>* photons_hcalIso = reader.GetFloat("photons_hcalIso");
 std::vector<int>* photons_isGap = reader.GetInt("photons_isGap");
 std::vector<float>* photons_e1x5 = reader.GetFloat("photons_e1x5");
 std::vector<float>* photons_e2x5 = reader.GetFloat("photons_e2x5");
 std::vector<float>* photons_e3x3 = reader.GetFloat("photons_e3x3");
 std::vector<float>* photons_e5x5 = reader.GetFloat("photons_e5x5");
 std::vector<float>* photons_maxEnergyXtal = reader.GetFloat("photons_maxEnergyXtal");
 std::vector<float>* photons_sigmaEtaEta = reader.GetFloat("photons_sigmaEtaEta");
 std::vector<float>* photons_sigmaIetaIeta = reader.GetFloat("photons_sigmaIetaIeta");
 std::vector<float>* photons_r1x5 = reader.GetFloat("photons_r1x5");
 std::vector<float>* photons_r2x5 = reader.GetFloat("photons_r2x5");
 std::vector<float>* photons_r9 = reader.GetFloat("photons_r9");

std::vector<float>* photons_hadronicOverEm = reader.GetFloat("photons_hadronicOverEm");

std::vector<float>* photons_ecalIso = reader.GetFloat("photons_ecalIso");
std::vector<float>* photons_hcalIso = reader.GetFloat("photons_hcalIso");
std::vector<float>*  photons_trkSumPtHollowConeDR04 = reader.GetFloat("photons_trkSumPtHollowConeDR04");

std::vector<float>* photons_hasPixelSeed  = reader.GetFloat("photons_hasPixelSeed ");

std::vector<ROOT::Math::XYZTVector>* photons_SC = reader.Get4V("photons_SC");

/*
 std::vector<ROOT::Math::XYZTVector>* Met = reader.Get4V("Met");
 std::vector<ROOT::Math::XYZTVector>* TCMet = reader.Get4V("TCMet");
 std::vector<ROOT::Math::XYZTVector>* PFMet = reader.Get4V("PFMet");
 std::vector<ROOT::Math::XYZTVector>* jets = reader.Get4V("jets");
 std::vector<float>* jets_charge = reader.GetFloat("jets_charge");
 std::vector<float>* jets_corrFactor_raw = reader.GetFloat("jets_corrFactor_raw");
 std::vector<float>* jets_corrFactor_off = reader.GetFloat("jets_corrFactor_off");
 std::vector<float>* jets_corrFactor_rel = reader.GetFloat("jets_corrFactor_rel");
 std::vector<float>* jets_corrFactor_abs = reader.GetFloat("jets_corrFactor_abs");
 std::vector<float>* trackCountingHighEffBJetTags = reader.GetFloat("trackCountingHighEffBJetTags");
 std::vector<float>* trackCountingHighPurBJetTags = reader.GetFloat("trackCountingHighPurBJetTags");
 std::vector<float>* simpleSecondaryVertexHighPurBJetTags = reader.GetFloat("simpleSecondaryVertexHighPurBJetTags");
 std::vector<float>* simpleSecondaryVertexHighEffBJetTags = reader.GetFloat("simpleSecondaryVertexHighEffBJetTags");
 std::vector<float>* jets_etaetaMoment = reader.GetFloat("jets_etaetaMoment");
 std::vector<float>* jets_phiphiMoment = reader.GetFloat("jets_phiphiMoment");
 std::vector<float>* jets_etaphiMoment = reader.GetFloat("jets_etaphiMoment");
 std::vector<float>* jets_fHPD = reader.GetFloat("jets_fHPD");
 std::vector<float>* jets_fRBX = reader.GetFloat("jets_fRBX");
 std::vector<float>* jets_n90Hits = reader.GetFloat("jets_n90Hits");
 std::vector<float>* jets_nHCALTowers = reader.GetFloat("jets_nHCALTowers");
 std::vector<float>* jets_nECALTowers = reader.GetFloat("jets_nECALTowers");
 std::vector<float>* jets_towersArea = reader.GetFloat("jets_towersArea");
 std::vector<float>* jets_emEnergyFraction = reader.GetFloat("jets_emEnergyFraction");
 std::vector<float>* jets_chargedHadronEnergyFraction = reader.GetFloat("jets_chargedHadronEnergyFraction");
 std::vector<float>* jets_neutralHadronEnergyFraction = reader.GetFloat("jets_neutralHadronEnergyFraction");
 std::vector<float>* jets_chargedEmEnergyFraction = reader.GetFloat("jets_chargedEmEnergyFraction");
 std::vector<float>* jets_neutralEmEnergyFraction = reader.GetFloat("jets_neutralEmEnergyFraction");
 std::vector<float>* jets_photonEnergyFraction = reader.GetFloat("jets_photonEnergyFraction");
 std::vector<float>* jets_muonEnergyFraction = reader.GetFloat("jets_muonEnergyFraction");
 std::vector<int>* jets_chargedMultiplicity = reader.GetInt("jets_chargedMultiplicity");
 std::vector<int>* jets_neutralMultiplicity = reader.GetInt("jets_neutralMultiplicity");
 std::vector<int>* jets_muonMultiplicity = reader.GetInt("jets_muonMultiplicity");
*/
 std::vector<float>* PV_noEle_normalizedChi2 = reader.GetFloat("PV_noEle_normalizedChi2");
 std::vector<int>* PV_noEle_ndof = reader.GetInt("PV_noEle_ndof");
 std::vector<int>* PV_noEle_nTracks = reader.GetInt("PV_noEle_nTracks");
 std::vector<float>* PV_noEle_z = reader.GetFloat("PV_noEle_z");
 std::vector<float>* PV_noEle_d0 = reader.GetFloat("PV_noEle_d0");
 std::vector<float>* PV_noEle_SumPt2 = reader.GetFloat("PV_noEle_SumPt2");

 std::vector<float>* PV_noMuon_normalizedChi2 = reader.GetFloat("PV_noMuon_normalizedChi2");
 std::vector<int>* PV_noMuon_ndof = reader.GetInt("PV_noMuon_ndof");
 std::vector<int>* PV_noMuon_nTracks = reader.GetInt("PV_noMuon_nTracks");
 std::vector<float>* PV_noMuon_z = reader.GetFloat("PV_noMuon_z");
 std::vector<float>* PV_noMuon_d0 = reader.GetFloat("PV_noMuon_d0");
 std::vector<float>* PV_noMuon_SumPt2 = reader.GetFloat("PV_noMuon_SumPt2");

 std::vector<ROOT::Math::XYZVector>* tracks = reader.Get3V("tracks");
 std::vector<int>* tracks_PVindex = reader.GetInt("tracks_PVindex");
 std::vector<float>* tracks_PVtracks_dR = reader.GetFloat("tracks_PVtracks_dR");
 std::vector<float>* tracks_dz = reader.GetFloat("tracks_dz");
