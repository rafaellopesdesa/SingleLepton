#include <iostream>

#include "MetSelections.h"

#include "Math/VectorUtil.h"
#include "MetSelections.h"

#include "Tools/jetcorr/FactorizedJetCorrector.h"
#include "Tools/JetCorrector.h"

using namespace tas;

metStruct trackerMET(float deltaZCut, const std::vector<LorentzVector>* jets) {

  if ( vtxs_isFake().empty() ) return metStruct();
  double pX(0), pY(0);
  
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    if ( pfcands_charge().at(i)==0 ) continue;
    if ( jets ){
      bool matched = false;
      for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet )
    	if ( fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),*jet))<0.5 ) matched=true;
      if (matched) continue;
    }
    
    if ( fabs(pfcands_dz().at(i)) > deltaZCut) continue;
    
    pX -= pfcands_p4().at(i).px();
    pY -= pfcands_p4().at(i).py();
  }
  
  if (jets){
    for ( std::vector<LorentzVector>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ){
      pX -= jet->px();
      pY -= jet->py();
    }
  }
  metStruct met;
  met.met     = sqrt(pX * pX + pY * pY);
  met.metphi  = atan2(pY, pX);
  met.metx = pX;
  met.mety = pY;
  return met;
}

bool hbheNoiseFilter(int minZeros) {
    // http://cmslxr.fnal.gov/lxr/source/CommonTools/RecoAlgos/plugins/HBHENoiseFilterResultProducer.cc?v=CMSSW_7_4_1
    // by default (false --> reject), maxZeros() cut will never cause a lost event
    if(hcalnoise_maxHPDHits()>=17) return false;
    if(hcalnoise_maxHPDNoOtherHits()>=10) return false;
    if(hcalnoise_maxZeros()>=minZeros) return false;
    if(hcalnoise_HasBadRBXTS4TS5()) return false;
    return true;
}

bool hbheNoiseFilter_25ns(int minZeros) {
    // http://cmslxr.fnal.gov/lxr/source/CommonTools/RecoAlgos/plugins/HBHENoiseFilterResultProducer.cc?v=CMSSW_7_4_1
    // by default (false --> reject), maxZeros() cut will never cause a lost event
    if(hcalnoise_maxHPDHits()>=17) return false;
    if(hcalnoise_maxHPDNoOtherHits()>=10) return false;
    if(hcalnoise_maxZeros()>=minZeros) return false;
    if(hcalnoise_HasBadRBXRechitR45Loose()) return false;
    return true;
}

bool hbheIsoNoiseFilter() {
    // false = reject event
    if(hcalnoise_numIsolatedNoiseChannels()>=10) return false;
    if(hcalnoise_isolatedNoiseSumE()>=50.0) return false;
    if(hcalnoise_isolatedNoiseSumEt()>=25.0) return false;
    return true;
}

// takes in an already initialized FactorizedJetCorrector object
// and returns T1 Corrected MET using the CHS jet collection
pair <float, float> getT1CHSMET( FactorizedJetCorrector * jet_corrector ){
  float T1_met    = cms3.evt_METToolbox_pfmet_raw();
  float T1_metPhi = cms3.evt_METToolbox_pfmetPhi_raw();
  float T1_metx   = T1_met * cos(T1_metPhi);
  float T1_mety   = T1_met * sin(T1_metPhi);

  //Run over same jets that were produced with MET tools
  for(unsigned int iJet = 0; iJet < cms3.pfjets_METToolbox_p4().size(); iJet++){

	// // get uncorrected jet p4 to use as input for corrections
	LorentzVector jetp4_uncorr = cms3.pfjets_METToolbox_p4().at(iJet);		  
	float emfrac = (cms3.pfjets_METToolbox_chargedEmE().at(iJet) + cms3.pfjets_METToolbox_neutralEmE().at(iJet)) / jetp4_uncorr.E();

	if (emfrac > 0.9                  ) continue; // veto events with EM fraction > 0.9
	if( abs(jetp4_uncorr.eta()) > 9.9 ) continue; // veto jets with eta > 9.9

	// get L1FastL2L3 total correction
	jet_corrector->setRho   ( cms3.evt_fixgridfastjetMETTools_all_rho()      );
	jet_corrector->setJetA  ( cms3.pfjets_METToolbox_area().at(iJet) );
	jet_corrector->setJetPt ( jetp4_uncorr.pt()                      );
	jet_corrector->setJetEta( jetp4_uncorr.eta()                     );

	//Note the subcorrections are stored with corr_vals(N) = corr(N)*corr(N-1)*...*corr(1)
	vector<float> corr_vals = jet_corrector->getSubCorrections();

	double corr             = corr_vals.at(corr_vals.size()-1); // All corrections
	double corr_l1          = corr_vals.at(0);                  // offset correction
		  
	//	
	// remove SA or global muons from jets before correcting
	//
	// for (unsigned int imu = 0; imu < cms3.mus_p4().size(); imu++)
	//   {
	// 	int index = cms3.mus_pfidx().at(imu);
	// 	if (index < 0) continue;
	// 	bool is_global     = !(((cms3.mus_type().at(imu)) & (1<<1)) == 0);
	// 	bool is_standalone = !(((cms3.mus_type().at(imu)) & (1<<3)) == 0);
	// 	if (!(is_global || is_standalone)) continue;            
	// 	if (std::find(cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).begin(),
	// 				  cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).end(), index) == cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).end()) continue;
	// 	jetp4_uncorr -= cms3.pfcands_p4()   .at(index);
	//   }

	// Alternative way to do muon corrections, done by MET group
	for (unsigned int pfcind = 0; pfcind < cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).size(); pfcind++){
	  int index = cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).at(pfcind);
	  if( cms3.pfcands_isGlobalMuon()    .at(index) ||
		  cms3.pfcands_isStandAloneMuon().at(index)){
		jetp4_uncorr -= cms3.pfcands_p4()   .at(index);
	  }
	}
			  
	if (corr * jetp4_uncorr.pt() > 10.){		  
	  T1_metx += jetp4_uncorr.px() * ( corr_l1 - corr );
	  T1_mety += jetp4_uncorr.py() * ( corr_l1 - corr );
	}

  }
	  
  T1_met    = std::sqrt(pow(T1_metx, 2) + pow(T1_mety, 2));
  T1_metPhi = std::atan2(T1_mety, T1_metx);

  return make_pair(T1_met, T1_metPhi);
}

// takes in an already initialized FactorizedJetCorrector object
// and returns T1 Corrected MET using the CHS jet collection from miniAOD
pair <float, float> getT1CHSMET_fromMINIAOD( FactorizedJetCorrector * jet_corrector ){
  float T1_met    = cms3.evt_pfmet_raw();
  float T1_metPhi = cms3.evt_pfmetPhi_raw();
  float T1_metx   = T1_met * cos(T1_metPhi);
  float T1_mety   = T1_met * sin(T1_metPhi);

  //Run over same jets that were produced with MET tools
  for(unsigned int iJet = 0; iJet < cms3.pfjets_p4().size(); iJet++){

	// // get uncorrected jet p4 to use as input for corrections
	LorentzVector jetp4_uncorr = cms3.pfjets_p4().at(iJet)*cms3.pfjets_undoJEC().at(iJet);
	float emfrac = (cms3.pfjets_chargedEmE().at(iJet) + cms3.pfjets_neutralEmE().at(iJet)) / jetp4_uncorr.E();

	if (emfrac > 0.9                  ) continue; // veto events with EM fraction > 0.9
	if( abs(jetp4_uncorr.eta()) > 9.9 ) continue; // veto jets with eta > 9.9

	// get L1FastL2L3 total correction
	jet_corrector->setRho   ( cms3.evt_fixgridfastjet_all_rho()      );
	jet_corrector->setJetA  ( cms3.pfjets_area().at(iJet) );
	jet_corrector->setJetPt ( jetp4_uncorr.pt()                      );
	jet_corrector->setJetEta( jetp4_uncorr.eta()                     );

	//Note the subcorrections are stored with corr_vals(N) = corr(N)*corr(N-1)*...*corr(1)
	vector<float> corr_vals = jet_corrector->getSubCorrections();

	double corr             = corr_vals.at(corr_vals.size()-1); // All corrections
	double corr_l1          = corr_vals.at(0);                  // offset correction
		  
	//	
	// remove SA or global muons from jets before correcting
	//
	for (unsigned int pfcind = 0; pfcind < cms3.pfjets_pfcandIndicies().at(iJet).size(); pfcind++){
	  int index = cms3.pfjets_pfcandIndicies().at(iJet).at(pfcind);
	  if( cms3.pfcands_isGlobalMuon()    .at(index) ||
		  cms3.pfcands_isStandAloneMuon().at(index)){
		jetp4_uncorr -= cms3.pfcands_p4()   .at(index);
	  }
	}
			  
	if (corr * jetp4_uncorr.pt() > 15.){		  
	  T1_metx += jetp4_uncorr.px() * ( corr_l1 - corr );
	  T1_mety += jetp4_uncorr.py() * ( corr_l1 - corr );
	}

  }
	  
  T1_met    = std::sqrt(pow(T1_metx, 2) + pow(T1_mety, 2));
  T1_metPhi = std::atan2(T1_mety, T1_metx);

  return make_pair(T1_met, T1_metPhi);
}

pair<float,float> MET3p0() { 

  float met_pt  = 0.0;
  float met_phi = 0.0;
  
  LorentzVector pfcands3p0_p4(0,0,0,0);
  
  for( size_t pfind = 0; pfind < cms3.pfcands_p4().size(); pfind++ ){
	
	if( abs( cms3.pfcands_p4().at(pfind).eta() ) < 3.0 ){ 
	  pfcands3p0_p4 -= cms3.pfcands_p4().at(pfind);
	}
  }
  met_pt  = pfcands3p0_p4.pt();
  met_phi = pfcands3p0_p4.phi();

  return make_pair( met_pt, met_phi ); 
}

bool passesMETfilter(bool is25){

  //primary vertex filter (re-run by user)
  if (firstGoodVertex() == -1) return false;

  //CSC beam halo filter 
  if (!filt_cscBeamHalo()) return false;
  
  //HBHE noise filter (re-run by user)
  if (!is25 && !hbheNoiseFilter()) return false;
  if ( is25 && !hbheNoiseFilter_25ns()) return false;

  //ee badSC noise filter NEW
  if (!filt_eeBadSc()) return false;

  //Otherwise good
  return true;

}

// takes in an already initialized FactorizedJetCorrector object
// and returns T1 Corrected MET using the CHS jet collection
// THIS FUNCTION IS NOT VALIDATED CURRENTLY!
pair <float, float> getT1CHSMET3p0(   FactorizedJetCorrector * jet_corrector ){
  float T1_met    = cms3.evt_METToolboxNoHF_pfmet_raw();
  float T1_metPhi = cms3.evt_METToolboxNoHF_pfmetPhi_raw();
  float T1_metx   = T1_met * cos(T1_metPhi);
  float T1_mety   = T1_met * sin(T1_metPhi);

  //Run over same jets that were produced with MET tools
  for(unsigned int iJet = 0; iJet < cms3.pfjets_METToolbox_p4().size(); iJet++){

	if( cms3.pfjets_METToolbox_p4().at(iJet).eta() > 3.0 ) continue;// only use jets with eta<3.0 to correct MET 3.0
	
	// // get uncorrected jet p4 to use as input for corrections
	LorentzVector jetp4_uncorr = cms3.pfjets_METToolbox_p4().at(iJet);		  
	float emfrac = (cms3.pfjets_METToolbox_chargedEmE().at(iJet) + cms3.pfjets_METToolbox_neutralEmE().at(iJet)) / jetp4_uncorr.E();

	if (emfrac > 0.9                  ) continue; // veto events with EM fraction > 0.9
	if( abs(jetp4_uncorr.eta()) > 9.9 ) continue; // veto jets with eta > 9.9

	// get L1FastL2L3 total correction
	jet_corrector->setRho   ( cms3.evt_fixgridfastjet_all_rho()      );
	jet_corrector->setJetA  ( cms3.pfjets_METToolbox_area().at(iJet) );
	jet_corrector->setJetPt ( jetp4_uncorr.pt()                      );
	jet_corrector->setJetEta( jetp4_uncorr.eta()                     );

	//Note the subcorrections are stored with corr_vals(N) = corr(N)*corr(N-1)*...*corr(1)
	vector<float> corr_vals = jet_corrector->getSubCorrections();

	double corr             = corr_vals.at(corr_vals.size()-1); // All corrections
	double corr_l1          = corr_vals.at(0);                  // offset correction

	// Alternative way to do muon corrections, done by MET group
	for (unsigned int pfcind = 0; pfcind < cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).size(); pfcind++){
	  int index = cms3.pfjets_METToolbox_pfcandIndicies().at(iJet).at(pfcind);
	  if( cms3.pfcands_isGlobalMuon()    .at(index) ||
		  cms3.pfcands_isStandAloneMuon().at(index)){
		jetp4_uncorr -= cms3.pfcands_p4()   .at(index);
	  }
	}
			  
	if (corr * jetp4_uncorr.pt() > 10.){		  
	  T1_metx += jetp4_uncorr.px() * ( corr_l1 - corr );
	  T1_mety += jetp4_uncorr.py() * ( corr_l1 - corr );
	}

  }
	  
  T1_met    = std::sqrt(pow(T1_metx, 2) + pow(T1_mety, 2));
  T1_metPhi = std::atan2(T1_mety, T1_metx);

  return make_pair(T1_met, T1_metPhi);
}

