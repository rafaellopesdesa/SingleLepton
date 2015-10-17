#include "JetSelections.h"
#include "ElectronSelections.h"
#include "MuonSelections.h" 
#include "Math/VectorUtil.h"

using namespace tas;

bool isLoosePFJet(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = pfjets_chargedMultiplicity()[pfJetIdx];
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjets_pfcandIndicies()[pfJetIdx].size() < 2) return false;
  if (pfjet_nef_ >= 0.99) return false;
  if (pfjet_nhf_ >= 0.99) return false;

  if (pfjet_eta < 2.4){
    if (pfjet_cm_ < 1) return false;
    if (pfjet_chf_ < 1e-6) return false;
    if (pfjet_cef_ >= 0.99) return false;
  }

  return true;
}

bool isMediumPFJet(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());

  if (pfjet_nef_ >= 0.95) return false;
  if (pfjet_nhf_ >= 0.95) return false;

  if (!isLoosePFJet(pfJetIdx)) return false;

  return true;
}

bool isTightPFJet(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());

  if (pfjet_nef_ >= 0.90) return false;
  if (pfjet_nhf_ >= 0.90) return false;

  if (!isLoosePFJet(pfJetIdx)) return false;

  return true;
}

bool isLoosePFJetV2(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_muf_  = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = pfjets_chargedMultiplicity()[pfJetIdx];
  int   pfjet_nm_  = pfjets_neutralMultiplicity()[pfJetIdx];
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_cm_ + pfjet_nm_ < 2) return false;
  if (pfjet_nef_ >= 0.99) return false;
  if (pfjet_nhf_ >= 0.99) return false;
  if (pfjet_muf_ >= 0.8) return false;

  if (pfjet_eta < 2.4){
    if (pfjet_cm_ < 1) return false;
    if (!(pfjet_chf_ > 0.)) return false;
    if (pfjet_cef_ >= 0.99) return false;
  }

  return true;
}

bool isTightPFJetV2(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_nef_ >= 0.90) return false;
  if (pfjet_nhf_ >= 0.90) return false;
  if (pfjet_eta < 2.4){
    if (pfjet_cef_ >= 0.90) return false;
  }


  if (!isLoosePFJetV2(pfJetIdx)) return false;

  return true;
}


bool isLoosePFJet_50nsV1(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  int   pfjet_cm_  = pfjets_chargedMultiplicity()[pfJetIdx];
  int   pfjet_nm_  = pfjets_neutralMultiplicity()[pfJetIdx];
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nhf_ >= 0.99       ) return false;
	if (pfjet_nef_ >= 0.99       ) return false;
	if (pfjet_cm_ + pfjet_nm_ < 2) return false;
	// if (pfjet_muf_ >= 0.8        ) return false; removed again

	if (pfjet_eta < 2.4){
	  if (!(pfjet_cm_  >   0.  ) ) return false;
	  if (!(pfjet_chf_ >   0.  ) ) return false;
	  if (!(pfjet_cef_ <   0.99) ) return false;
	}
  }else if( pfjet_eta > 3.0 ){
	if (!(pfjet_nef_ < 0.9 ) ) return false;
	if (!(pfjet_nm_  > 10  ) ) return false;
  }

  return true;
}

bool isTightPFJet_50nsV1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
  }

  if (!isLoosePFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isTightPFJetLepVeto_50nsV1(unsigned int pfJetIdx){

  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_cef_  = pfjets_chargedEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_muf_  = pfjets_muonE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta <= 3.0){
	if (pfjet_nef_ >= 0.90) return false;
	if (pfjet_nhf_ >= 0.90) return false;
	if (pfjet_eta <= 2.4){
	  if (pfjet_cef_ >= 0.90) return false;
	}
	if (pfjet_muf_ >= 0.8        ) return false;
  }

  if (!isLoosePFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isMonoPFJet_MT2(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta > 3.0) return false;
  if (pfjet_chf_ <= 0.05) return false;
  if (pfjet_nhf_ >= 0.70) return false;
  if (pfjet_nef_ >= 0.80) return false;

  if (!isTightPFJet_50nsV1(pfJetIdx)) return false;

  return true;
}

bool isMonoPFJet_Monojet(unsigned int pfJetIdx){

  float pfjet_chf_  = pfjets_chargedHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nhf_  = pfjets_neutralHadronE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_nef_  = pfjets_neutralEmE()[pfJetIdx] / (pfjets_undoJEC().at(pfJetIdx)*pfjets_p4()[pfJetIdx].energy());
  float pfjet_eta  = fabs(pfjets_p4()[pfJetIdx].eta());

  if (pfjet_eta > 3.0) return false;
  if (pfjet_chf_ <= 0.20) return false;
  if (pfjet_nhf_ >= 0.70) return false;
  if (pfjet_nef_ >= 0.70) return false;

  if (!isMonoPFJet_MT2(pfJetIdx)) return false;

  return true;
}

bool loosePileupJetId(unsigned int pfJetIdx){

  float eta = fabs(pfjets_p4().at(pfJetIdx).eta());
  float value = pfjets_pileupJetId().at(pfJetIdx);

  if( (eta >= 0   ) && (eta <= 2.5 ) && (value > -0.63) ) return true;
  if( (eta > 2.5  ) && (eta <= 2.75) && (value > -0.60) ) return true;
  if( (eta > 2.75 ) && (eta <= 3.0 ) && (value > -0.55) ) return true;
  if( (eta > 3.0  ) && (eta <= 5.2 ) && (value > -0.45) ) return true;
    
  return false;
}

bool JetIsElectron(LorentzVector pfJet, id_level_t id_level, float ptcut, float deltaR){
  bool jetIsLep = false;
  for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++){
    LorentzVector electron = tas::els_p4().at(eidx);
    if (electron.pt() < ptcut) continue;
    if (!electronID(eidx,id_level)) continue;
    if (ROOT::Math::VectorUtil::DeltaR(pfJet, electron) > deltaR) continue;
    jetIsLep = true;
  }
  return jetIsLep;
}

bool JetIsMuon(LorentzVector pfJet, id_level_t id_level, float ptcut, float deltaR){
  bool jetIsLep = false;
  for (unsigned int muidx = 0; muidx < tas::mus_p4().size(); muidx++){
    LorentzVector muon = tas::mus_p4().at(muidx);
    if (muon.pt() < ptcut) continue;
    if (!muonID(muidx,id_level)) continue;
    if (ROOT::Math::VectorUtil::DeltaR(pfJet, muon) > deltaR) continue;
    jetIsLep = true;
  }
  return jetIsLep;
}
