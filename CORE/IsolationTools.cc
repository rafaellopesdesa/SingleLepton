#include "Math/VectorUtil.h"
#include "IsolationTools.h"
#include "Tools/JetCorrector.h"

using namespace std;
using namespace tas;

FactorizedJetCorrector *jetCorrAG3 = 0;
FactorizedJetCorrector *jetCorrAG2 = 0;

bool passMultiIsoCuts(float cutMiniIso, float cutPtRatio, float cutPtRel, float miniIsoValue, float ptRatioValue, float ptRelValue){
  return (miniIsoValue < cutMiniIso && (ptRatioValue>cutPtRatio || ptRelValue > cutPtRel));
}

bool passMultiIso(int id, int idx, float cutMiniIso, float cutPtRatio, float cutPtRel, int eaversion, int whichCorr){
  const LorentzVector& lep_p4 = abs(id)==11 ? els_p4().at(idx) : mus_p4().at(idx);
  const LorentzVector& jet_p4 = closestJet(lep_p4,0.4,2.4,whichCorr);
  float miniIso = abs(id)==11 ? elMiniRelIsoCMS3_EA(idx,eaversion) : muMiniRelIsoCMS3_EA(idx,eaversion);
  float closeJetPt = jet_p4.pt();
  float ptratio = ( closeJetPt>0. ? lep_p4.pt()/closeJetPt : 1.);
  float ptrel = ptRel(lep_p4, jet_p4, true);
  return passMultiIsoCuts(cutMiniIso, cutPtRatio, cutPtRel, miniIso, ptratio, ptrel);
}

bool passPtRel(int id, int idx, float cut, bool subtractLep, int whichCorr) {
  return getPtRel(id, idx, subtractLep, whichCorr) > cut;
}

float getPtRel(int id, int idx, bool subtractLep, int whichCorr) {
  LorentzVector lep_p4 = abs(id)==11 ? els_p4().at(idx) : mus_p4().at(idx);
  LorentzVector jet_p4 = closestJet(lep_p4,0.4,2.4,whichCorr);
  return ptRel(lep_p4, jet_p4, subtractLep);
}

float ptRel(const LorentzVector& lepp4, const LorentzVector& jetp4, bool subtractLep) {
  if (jetp4.pt()<=0.) return 0.;
  LorentzVector jp4 = jetp4;
  if (subtractLep) jp4-=lepp4;
  float dot = lepp4.Vect().Dot( jp4.Vect() );
  float ptrel = lepp4.P2() - dot*dot/jp4.P2();
  ptrel = ptrel>0 ? sqrt(ptrel) : 0.0;
  return ptrel;
}

int closestJetIdx(const LorentzVector& lep_p4, float dRmin, float maxAbsEta){
  int closestIdx = -1;
  for (unsigned int pfjidx=0; pfjidx<pfjets_p4().size(); ++pfjidx){
    LorentzVector jet_p4 = pfjets_p4().at(pfjidx);
    if (fabs(jet_p4.eta()) > maxAbsEta) continue;
    float tmp_dRmin = ROOT::Math::VectorUtil::DeltaR(jet_p4, lep_p4);
    if (tmp_dRmin < dRmin){
      closestIdx = pfjidx;
      dRmin = tmp_dRmin;
    }
  }
  return closestIdx;
}

LorentzVector closestJet(const LorentzVector& lep_p4, float dRmin, float maxAbsEta, int whichCorr){
  int closestIdx = closestJetIdx(lep_p4,dRmin,maxAbsEta);
  if (closestIdx < 0) return LorentzVector();
  LorentzVector jet = pfjets_p4().at(closestIdx);
  LorentzVector rawjet = pfjets_p4().at(closestIdx)*pfjets_undoJEC().at(closestIdx);
  
  if (whichCorr == 0) return jet; 

  //Calculate JEC
  if (jetCorrAG3 == 0){
    std::vector<std::string> filenames;
    if (tas::evt_isRealData() && whichCorr != 1) {
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV5_DATA_L1FastJet_AK4PFchs.txt");
    } else {
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
    }
    jetCorrAG3 = makeJetCorrector(filenames);
  }
  if (jetCorrAG2 == 0){
    std::vector<std::string> filenames;
    if (tas::evt_isRealData()) { 
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV5_DATA_L2Relative_AK4PFchs.txt");
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV5_DATA_L3Absolute_AK4PFchs.txt");
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV5_DATA_L2L3Residual_AK4PFchs.txt");
    } else {
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV2_MC_L2Relative_AK4PFchs.txt");
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
      filenames.push_back("CORE/Tools/jetcorr/data/run2_25ns/Summer15_25nsV2_MC_L2L3Residual_AK4PFchs.txt");
    }
    jetCorrAG2 = makeJetCorrector(filenames);
  }
  jetCorrAG3->setJetEta(rawjet.eta()); 
  jetCorrAG3->setJetPt(rawjet.pt()); 
  jetCorrAG3->setJetA(tas::pfjets_area().at(closestIdx)); 
  jetCorrAG3->setRho(tas::evt_fixgridfastjet_all_rho()); 
  jetCorrAG2->setJetEta(rawjet.eta()); 
  jetCorrAG2->setJetPt(rawjet.pt()); 
  jetCorrAG2->setJetA(tas::pfjets_area().at(closestIdx)); 
  jetCorrAG2->setRho(tas::evt_fixgridfastjet_all_rho()); 
  float JEC1 = jetCorrAG3->getCorrection(); 
  float JEC2 = jetCorrAG2->getCorrection(); 

  if (whichCorr == 1) return rawjet*JEC1;
  return (rawjet*JEC1 - lep_p4)*JEC2 + lep_p4;
}

float getMiniDR(float pt) {
  float dr = 0.2;
  if (pt>50) dr = 10./pt;
  if (pt>200) dr = 0.05;
  return dr;
}

float muRelIso03DB(unsigned int muIdx){
  float chiso     = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR03_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR03_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIso04DB(unsigned int muIdx){
  float chiso     = mus_isoR04_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR04_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR04_pf_PhotonEt().at(muIdx);
  float deltaBeta = mus_isoR04_pf_PUPt().at(muIdx);
  float absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIso03(unsigned int muIdx, analysis_t analysis){
  if (analysis == WW  ) return muRelIso03EA(muIdx);
  if (analysis == SS  ) return muRelIso03EA(muIdx,1);
  if (analysis == ZMET) return muRelIso03EA(muIdx);
  return muRelIso03DB(muIdx);
}

float muRelIso04(unsigned int muIdx, analysis_t analysis){
  return muRelIso04DB(muIdx);
}

float muEA03(unsigned int muIdx, int version){
  float ea = 0.;
  if (version==0) {
    //PHYS14 version
    if      (fabs(mus_p4().at(muIdx).eta())<=0.800) ea = 0.0913;
    else if (fabs(mus_p4().at(muIdx).eta())<=1.300) ea = 0.0765;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.000) ea = 0.0546;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.200) ea = 0.0728;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.500) ea = 0.1177;
  } else {
    //Spring15 version
    if      (fabs(mus_p4().at(muIdx).eta())<=0.800) ea = 0.0735;
    else if (fabs(mus_p4().at(muIdx).eta())<=1.300) ea = 0.0619;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.000) ea = 0.0465;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.200) ea = 0.0433;
    else if (fabs(mus_p4().at(muIdx).eta())<=2.500) ea = 0.0577;
  }
  return ea;
}

float muRelIso03EA(unsigned int muIdx, int eaversion){
  float chiso     = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso     = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso     = mus_isoR03_pf_PhotonEt().at(muIdx);
  float ea = muEA03(muIdx, eaversion);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea);
  return absiso/(mus_p4().at(muIdx).pt());
}
 
float muRelIso03_noCorr(unsigned int muIdx){
  float chiso  = mus_isoR03_pf_ChargedHadronPt().at(muIdx);
  float nhiso  = mus_isoR03_pf_NeutralHadronEt().at(muIdx);
  float emiso  = mus_isoR03_pf_PhotonEt().at(muIdx);
  float absiso = chiso + nhiso + emiso;
  return absiso/(mus_p4().at(muIdx).pt());
}

float muRelIsoCustomCone(unsigned int muIdx, float dr, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor, float mindr, int eaversion){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float correction = 0.;
  float deltaBeta = 0.;
  float deadcone_ch = 0.0001;
  float deadcone_pu = 0.01;
  float deadcone_ph = 0.01;
  float deadcone_nh = 0.01;

  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),mus_p4().at(muIdx)));
    if (thisDR<mindr) continue;
    if (thisDR>dr) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 ) {
      if (pfcands_fromPV().at(i) > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pfcands_p4().at(i).pt();
      else if (useDBcor && pfcands_fromPV().at(i) <= 1 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) deltaBeta+=pfcands_p4().at(i).pt();
    }
    if (fabs(pfcands_particleId().at(i))==130 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_nh)) nhiso+=pfcands_p4().at(i).pt();
    if (fabs(pfcands_particleId().at(i))==22 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_ph)) emiso+=pfcands_p4().at(i).pt();
  }

  if (useDBcor) correction = 0.5 * deltaBeta;
  else if (useEAcor) correction = evt_fixgridfastjet_all_rho() * muEA03(muIdx,eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - correction);
  return absiso/(mus_p4().at(muIdx).pt());
}

float muMiniRelIso(unsigned int idx, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor){
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return  muRelIsoCustomCone(idx,dr,useVetoCones,ptthresh,useDBcor,useEAcor);
}

float muMiniRelIsoCMS3_DB(unsigned int idx) {
  float correction = 0.5 * mus_miniIso_db().at(idx);
  float absiso = mus_miniIso_ch().at(idx) + std::max(float(0.0), mus_miniIso_nh().at(idx) + mus_miniIso_em().at(idx) - correction);
  return absiso/(mus_p4().at(idx).pt());
}

float muMiniRelIsoCMS3_EA(unsigned int idx, int eaversion) {
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  float rho = evt_fixgridfastjet_all_rho();
  if (eaversion!=0) rho = evt_fixgridfastjet_centralneutral_rho();
  float correction = rho * muEA03(idx,eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = mus_miniIso_ch().at(idx) + std::max(float(0.0), mus_miniIso_nh().at(idx) + mus_miniIso_em().at(idx) - correction);
  return absiso/(mus_p4().at(idx).pt());
}

float muRelIsoAn04(unsigned int idx, bool useDBcor) {
  float pt = mus_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return muRelIsoCustomCone(idx, 0.4, false, 0.0, useDBcor, false, dr);
}

float eleRelIso03(unsigned int elIdx, analysis_t analysis){
  if (analysis == HAD ) return eleRelIso03DB(elIdx);
  if (analysis == STOP) return eleRelIso03DB(elIdx);
  if (analysis == WW  ) return eleRelIso03_90ContEA(elIdx);
  if (analysis == SS  ) return eleRelIso03EA(elIdx,1);
  if (analysis == ZMET) return eleRelIso03EA(elIdx);

  else return eleRelIso03EA(elIdx,0);
}

float eleRelIso03DB(unsigned int elIdx){
  float chiso     = els_pfChargedHadronIso().at(elIdx);
  float nhiso     = els_pfNeutralHadronIso().at(elIdx);
  float emiso     = els_pfPhotonIso().at(elIdx);
  float deltaBeta = els_pfPUIso().at(elIdx);
  float absiso    = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return absiso/(els_p4().at(elIdx).pt());
}

float elEA03(unsigned int elIdx, int version) {
  float ea = 0.;
  if (version==0) {
    //PHYS14 version
    if      (fabs(els_p4().at(elIdx).eta())<=0.800) ea = 0.1013;
    else if (fabs(els_p4().at(elIdx).eta())<=1.300) ea = 0.0988;
    else if (fabs(els_p4().at(elIdx).eta())<=2.000) ea = 0.0572;
    else if (fabs(els_p4().at(elIdx).eta())<=2.200) ea = 0.0842;
    else if (fabs(els_p4().at(elIdx).eta())<=2.500) ea = 0.1530;
  } else {
    //Spring15 version
    if      (fabs(els_etaSC().at(elIdx))<=1.000) ea = 0.1752;
    else if (fabs(els_etaSC().at(elIdx))<=1.479) ea = 0.1862;
    else if (fabs(els_etaSC().at(elIdx))<=2.000) ea = 0.1411;
    else if (fabs(els_etaSC().at(elIdx))<=2.200) ea = 0.1534;
    else if (fabs(els_etaSC().at(elIdx))<=2.300) ea = 0.1903;
    else if (fabs(els_etaSC().at(elIdx))<=2.400) ea = 0.2243;
    else if (fabs(els_etaSC().at(elIdx))<=2.500) ea = 0.2687;
  }
  return ea;
}

float el90ContEA03(unsigned int elIdx) {
  float ea = 0.;
  if      (fabs(els_p4().at(elIdx).eta())<=1.000) ea = 0.1752;
  else if (fabs(els_p4().at(elIdx).eta())<=1.479) ea = 0.1862;
  else if (fabs(els_p4().at(elIdx).eta())<=2.000) ea = 0.1411;
  else if (fabs(els_p4().at(elIdx).eta())<=2.200) ea = 0.1534;
  else if (fabs(els_p4().at(elIdx).eta())<=2.300) ea = 0.1903;
  else if (fabs(els_p4().at(elIdx).eta())<=2.400) ea = 0.2243;
  else if (fabs(els_p4().at(elIdx).eta())<=2.500) ea = 0.2687;
  return ea;
}

float eleRelIso03EA(unsigned int elIdx, int eaversion){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float ea    = elEA03(elIdx, eaversion);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea);
  return absiso/(els_p4().at(elIdx).pt());
}

float eleRelIso03_90ContEA(unsigned int elIdx){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float ea    = el90ContEA03(elIdx);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - evt_fixgridfastjet_all_rho() * ea);
  return absiso/(els_p4().at(elIdx).pt());
}

float eleRelIso03_noCorr(unsigned int elIdx){
  float chiso = els_pfChargedHadronIso().at(elIdx);
  float nhiso = els_pfNeutralHadronIso().at(elIdx);
  float emiso = els_pfPhotonIso().at(elIdx);
  float absiso = chiso + nhiso + emiso;
  return absiso/(els_p4().at(elIdx).pt());
}

float elRelIsoCustomCone(unsigned int elIdx, float dr, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor, float mindr, int eaversion){
  float chiso     = 0.;
  float nhiso     = 0.;
  float emiso     = 0.;
  float correction = 0.;
  float deltaBeta = 0.;
  float deadcone_ch = 0.;
  float deadcone_pu = 0.;
  float deadcone_ph = 0.;
  // veto cones only in the endcap for electrons
  if (useVetoCones && fabs(els_etaSC().at(elIdx)) > 1.479) {
    deadcone_ch = 0.015;
    deadcone_pu = 0.015;
    deadcone_ph = 0.08;
  }
  for (unsigned int i=0; i<pfcands_particleId().size(); ++i){
    float thisDR = fabs(ROOT::Math::VectorUtil::DeltaR(pfcands_p4().at(i),els_p4().at(elIdx)));
    if ( thisDR<mindr ) continue;
    if ( thisDR>dr ) continue;  
    if ( fabs(pfcands_particleId().at(i))==211 ) {
      if (pfcands_fromPV().at(i) > 1 && (!useVetoCones || thisDR > deadcone_ch) ) chiso+=pfcands_p4().at(i).pt();
      else if (useDBcor && pfcands_fromPV().at(i) <= 1 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_pu)) deltaBeta+=pfcands_p4().at(i).pt();
    }
    if ( fabs(pfcands_particleId().at(i))==130 && (pfcands_p4().at(i).pt() > ptthresh) ) nhiso+=pfcands_p4().at(i).pt();
    if ( fabs(pfcands_particleId().at(i))==22 && (pfcands_p4().at(i).pt() > ptthresh) && (!useVetoCones || thisDR > deadcone_ph) ) emiso+=pfcands_p4().at(i).pt();
  }
  if (useDBcor) correction = 0.5 * deltaBeta;
  else if (useEAcor) correction = evt_fixgridfastjet_all_rho() * elEA03(elIdx, eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = chiso + std::max(float(0.0), nhiso + emiso - correction);
  return absiso/(els_p4().at(elIdx).pt());
}
float elMiniRelIso(unsigned int idx, bool useVetoCones, float ptthresh, bool useDBcor, bool useEAcor){

  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return elRelIsoCustomCone(idx,dr,useVetoCones,ptthresh,useDBcor,useEAcor);
}

float elMiniRelIsoCMS3_DB(unsigned int idx) {
  float correction = 0.5 * els_miniIso_db().at(idx);
  float absiso = els_miniIso_ch().at(idx) + std::max(float(0.0), els_miniIso_nh().at(idx) + els_miniIso_em().at(idx) - correction);
  return absiso/(els_p4().at(idx).pt());
}

float elMiniRelIsoCMS3_EA(unsigned int idx, int eaversion) {
  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  float rho = evt_fixgridfastjet_all_rho();
  if (eaversion!=0) rho = evt_fixgridfastjet_centralneutral_rho();
  float correction = rho * elEA03(idx, eaversion) * (dr/0.3) * (dr/0.3);
  float absiso = els_miniIso_ch().at(idx) + std::max(float(0.0), els_miniIso_nh().at(idx) + els_miniIso_em().at(idx) - correction);
  return absiso/(els_p4().at(idx).pt());
}

float elRelIsoAn04(unsigned int idx, bool useDBcor) {
  float pt = els_p4().at(idx).pt();
  float dr = getMiniDR(pt);
  return elRelIsoCustomCone(idx, 0.4, false, 0.0, useDBcor, false, dr);
}

//////////////////////
// Photon Isolation //
//////////////////////

// from https://indico.cern.ch/event/369239/contribution/2/attachments/1134693/1623149/spring15_pcb.pdf
float photon_CHEA03( int photonIdx )
{
  float eta = cms3.photons_p4().at(photonIdx).eta();
  float EA = -999;
  
  if(       abs(eta) < 1.0   ){ EA = 0.0157;
  }else if( abs(eta) < 1.479 ){ EA = 0.0143;
  }else if( abs(eta) < 2.0   ){ EA = 0.0115;
  }else if( abs(eta) < 2.2   ){ EA = 0.0094;
  }else if( abs(eta) < 2.3   ){ EA = 0.0095;
  }else if( abs(eta) < 2.4   ){ EA = 0.0068;
  }else if( abs(eta) > 2.4   ){ EA = 0.0053;
  }
  
  return EA;
}

float photon_NHEA03( int photonIdx )
{
  float eta = cms3.photons_p4().at(photonIdx).eta();
  float EA = -999;
  
  if(       abs(eta) < 1.0   ){ EA = 0.0143;
  }else if( abs(eta) < 1.479 ){ EA = 0.0210;
  }else if( abs(eta) < 2.0   ){ EA = 0.0148;
  }else if( abs(eta) < 2.2   ){ EA = 0.0082;
  }else if( abs(eta) < 2.3   ){ EA = 0.0124;
  }else if( abs(eta) < 2.4   ){ EA = 0.0186;
  }else if( abs(eta) > 2.4   ){ EA = 0.0320;
  }
  
  return EA;
}

float photon_EMEA03( int photonIdx )
{
  float eta = cms3.photons_p4().at(photonIdx).eta();
  float EA = -999;
  
  if(       abs(eta) < 1.0   ){ EA = 0.0725;
  }else if( abs(eta) < 1.479 ){ EA = 0.0604;
  }else if( abs(eta) < 2.0   ){ EA = 0.0320;
  }else if( abs(eta) < 2.2   ){ EA = 0.0512;
  }else if( abs(eta) < 2.3   ){ EA = 0.0766;
  }else if( abs(eta) < 2.4   ){ EA = 0.0949;
  }else if( abs(eta) > 2.4   ){ EA = 0.1160;
  }
  
  return EA;
}

float photonCHIso03EA( int photonIdx )
{
  float chiso = photons_recoChargedHadronIso().at(photonIdx);
  float ea    = photon_CHEA03(photonIdx);
  float CHIso = std::max(float(0.0), chiso - evt_fixgridfastjet_all_rho() * ea);
  return CHIso;
}

float photonNHIso03EA( int photonIdx )
{
  float nhiso = photons_recoNeutralHadronIso().at(photonIdx);
  float ea    = photon_NHEA03(photonIdx);
  float NHIso = std::max(float(0.0), nhiso - evt_fixgridfastjet_all_rho() * ea);
  return NHIso;
}

float photonEMIso03EA( int photonIdx )
{
  float emiso = photons_recoPhotonIso().at(photonIdx);
  float ea    = photon_EMEA03(photonIdx);
  float EMIso = std::max(float(0.0), emiso - evt_fixgridfastjet_all_rho() * ea);
  return EMIso;
}
