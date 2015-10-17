#include "WWSelections.h"

using namespace tas;
using namespace WWAnalysis;

int WWAnalysis::convertCMS3tag(TString tagName) 
{

  if (tagName.Length() == 0) return -1;
  if (!tagName.BeginsWith("CMS3")) return 9999999;
  TString t1(tagName[6]);
  TString t2(tagName[7]);
  TString t3(tagName[9]);
  TString t4(tagName[10]);
  TString t5(tagName[12]);
  TString t6(tagName[13]);
  return (t6.Atoi()+10*t5.Atoi()+100*t4.Atoi()+1000*t3.Atoi()+10000*t2.Atoi()+100000*t1.Atoi());
  
}

float WWAnalysis::getGenHT(bool is_b_a_jet)
{

  float ht_gen = 0.;
  if (tas::evt_isRealData()) return ht_gen;
  for (size_t gidx = 0; gidx < tas::genps_p4().size(); gidx++){
    if (tas::genps_status().at(gidx) != 23) continue;
    int id = abs(tas::genps_id().at(gidx));
    if ((id >= 1 && id<=4) || (id == 5 && is_b_a_jet) || (id == 21))
      ht_gen += tas::genps_p4().at(gidx).pt();
  }
  return ht_gen;

}

bool WWAnalysis::makesExtraGammaStar(int iHyp){

  std::vector<unsigned int> ele_idx;
  std::vector<unsigned int> mu_idx;

  int lt_id           = tas::hyp_lt_id().at(iHyp);
  int ll_id           = tas::hyp_ll_id().at(iHyp);
  unsigned int lt_idx = tas::hyp_lt_index().at(iHyp);
  unsigned int ll_idx = tas::hyp_ll_index().at(iHyp);

  (abs(lt_id) == 11) ? ele_idx.push_back(lt_idx) : mu_idx.push_back(lt_idx);
  (abs(ll_id) == 11) ? ele_idx.push_back(ll_idx) : mu_idx.push_back(ll_idx);

  if (ele_idx.size() + mu_idx.size() != 2) {
    std::cout << "ERROR: don't have 2 leptons in hypothesis!!!  Exiting..." << std::endl;
    return false;
  }
      
  if (ele_idx.size() > 0){
    for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++) {
      if (fabs(tas::els_p4().at(eidx).eta()) > 2.5) continue;
      if (tas::els_p4().at(eidx).pt() < 7.0) continue;
      if (!isGoodVetoElectron(eidx)) continue;

      for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {
          if (tas::els_charge().at(eidx) * tas::els_charge().at(ele_idx.at(vidx)) > 0) continue;
          LorentzVector gamma_p4 = tas::els_p4().at(eidx) + tas::els_p4().at(ele_idx.at(vidx));
          float gammacandmass = sqrt(fabs(gamma_p4.mass2()));
          if (gammacandmass < 12.0) return true;
      }
    } 
  }

  if (mu_idx.size() > 0){
    for (unsigned int midx = 0; midx < tas::mus_p4().size(); midx++) {
      if (fabs(tas::mus_p4().at(midx).eta()) > 2.4) continue;
      if (tas::mus_p4().at(midx).pt() < 5.0) continue;
      if (!isGoodVetoMuon(midx)) continue;

      for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {
          if (tas::mus_charge().at(midx) * tas::mus_charge().at(mu_idx.at(vidx)) > 0) continue;
          LorentzVector gamma_p4 = tas::mus_p4().at(midx) + tas::mus_p4().at(mu_idx.at(vidx));
          float gammacandmass = sqrt(fabs(gamma_p4.mass2()));
          if (gammacandmass < 12.0) return true;
      }
    }
  }

  return false;
}

Z_result_t WWAnalysis::makesExtraZ(int iHyp){

  Z_result_t result;
  result.result = false;
  result.id = 0;
  result.idx = 0; 

  std::vector<unsigned int> ele_idx;
  std::vector<unsigned int> mu_idx;

  int lt_id           = tas::hyp_lt_id().at(iHyp);
  int ll_id           = tas::hyp_ll_id().at(iHyp);
  unsigned int lt_idx = tas::hyp_lt_index().at(iHyp);
  unsigned int ll_idx = tas::hyp_ll_index().at(iHyp);

  (abs(lt_id) == 11) ? ele_idx.push_back(lt_idx) : mu_idx.push_back(lt_idx);
  (abs(ll_id) == 11) ? ele_idx.push_back(ll_idx) : mu_idx.push_back(ll_idx);

  if (ele_idx.size() + mu_idx.size() != 2) {
    std::cout << "ERROR: don't have 2 leptons in hypothesis!!!  Exiting" << std::endl;
    return result;
  }
      
  if (ele_idx.size() > 0) {
    for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++) {

      if (fabs(tas::els_p4().at(eidx).eta()) > 2.5) continue;
      if (tas::els_p4().at(eidx).pt() < 7) continue;

      if (!isGoodVetoElectron(eidx)) continue;

      for (unsigned int vidx = 0; vidx < ele_idx.size(); vidx++) {
        if (tas::els_charge().at(eidx) * tas::els_charge().at(ele_idx.at(vidx)) > 0) continue;
        LorentzVector zp4 = tas::els_p4().at(eidx) + tas::els_p4().at(ele_idx.at(vidx));
        float zcandmass = sqrt(fabs(zp4.mass2()));
        if (fabs(zcandmass-91.1875) < 15.){
          result.result = true; 
          result.id = -sgn(tas::els_charge().at(eidx))*11;
          result.idx = eidx;  
          return result;
        }
      }
    }
  }

  if (mu_idx.size() > 0) {
    for (unsigned int midx = 0; midx < tas::mus_p4().size(); midx++) {

      if (fabs(tas::mus_p4().at(midx).eta()) > 2.4) continue;
      if (tas::mus_p4().at(midx).pt() < 7.) continue;

      if (!isGoodVetoMuon(midx)) continue;

      for (unsigned int vidx = 0; vidx < mu_idx.size(); vidx++) {
        if (tas::mus_charge().at(midx) * tas::mus_charge().at(mu_idx.at(vidx)) > 0) continue;
        LorentzVector zp4 = tas::mus_p4().at(midx) + tas::mus_p4().at(mu_idx.at(vidx));
        float zcandmass = sqrt(fabs(zp4.mass2()));
        if (fabs(zcandmass-91.1875) < 15.){
          result.result = true; 
          result.id = -sgn(tas::mus_charge().at(midx))*13;
          result.idx = midx;  
          return result;
        }
      }
    }
  }

  return result;
}

std::pair <vector <Jet>, vector <Jet> > WWAnalysis::WWJetsCalculator(vector<LorentzVector> JetCollection){
  vector <Jet> result_jets;
  vector <Jet> result_btags;

  
  for (unsigned int i = 0; i < JetCollection.size(); i++){
    LorentzVector jet = JetCollection.at(i);
    
    //Kinematic jet cuts
    if (jet.pt() < 10.) continue;
    if (fabs(jet.eta()) > 4.7) continue;
    
    //Require loose jet ID
    if (!isLoosePFJet_50nsV1(i)) continue;
    
    //Jet cleaning -- electrons
    bool jetIsLep = false;
    for (unsigned int eidx = 0; eidx < tas::els_p4().size(); eidx++){
      LorentzVector electron = tas::els_p4().at(eidx);
      if (electron.pt() < 10) continue;
      if (!isGoodVetoElectron(eidx)) continue;
      if (ROOT::Math::VectorUtil::DeltaR(jet, electron) > 0.3) continue;
      jetIsLep = true;
    }
    if (jetIsLep == true) continue;
    
    //Jet cleaning -- muons
    for (unsigned int muidx = 0; muidx < tas::mus_p4().size(); muidx++){
      LorentzVector muon = tas::mus_p4().at(muidx);
      if (muon.pt() < 10) continue;
      if (!isGoodVetoMuon(muidx)) continue;
      if (ROOT::Math::VectorUtil::DeltaR(jet, muon) > 0.3) continue;
      jetIsLep = true;
    }
    if (jetIsLep == true) continue;      
        
    //Save jets that make it this far
    result_jets.push_back(Jet(i,
    			      tas::getbtagvalue("pfCombinedInclusiveSecondaryVertexV2BJetTags",i),
    			      tas::getbtagvalue("pfCombinedSecondaryVertexSoftLeptonBJetTags",i),
    			      tas::getbtagvalue("pfTrackCountingHighEffBJetTags",i)
    			      ));
    float disc = tas::getbtagvalue("pfCombinedInclusiveSecondaryVertexV2BJetTags",i);
    if (disc > 0. && disc < 0.605) continue;
    result_btags.push_back(Jet(i,
    			       tas::getbtagvalue("pfCombinedInclusiveSecondaryVertexV2BJetTags",i),
    			       tas::getbtagvalue("pfCombinedSecondaryVertexSoftLeptonBJetTags",i),
			       tas::getbtagvalue("pfTrackCountingHighEffBJetTags",i)
			       )); 
    
  }
  std::pair <vector <Jet>, vector <Jet> > result = std::make_pair(result_jets, result_btags);
  return result;
}

bool WWAnalysis::isLooseIsolatedLepton(int id, int idx){
  if (abs(id) == 11) return (elMiniRelIsoCMS3_EA(idx) < 0.4);
  else if (abs(id) == 13) return (muMiniRelIsoCMS3_EA(idx) < 0.4);
  return false;
}
bool WWAnalysis::isIsolatedLepton(int id, int idx){
  if (abs(id) == 11) return (elMiniRelIsoCMS3_EA(idx) < 0.1);
  else if (abs(id) == 13) return (muMiniRelIsoCMS3_EA(idx) < 0.1);
  return false;
}

bool WWAnalysis::isGoodLepton(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_medium_v2);
  else if (abs(id) == 13) return muonID(idx, WW_medium_v2);
  return false;
}

bool WWAnalysis::isGoodLeptonNoIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_medium_noiso_v2);
  else if (abs(id) == 13) return muonID(idx, WW_medium_noiso_v2);
  return false;
}

bool WWAnalysis::isGoodLeptonIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_medium_v2);
  else if (abs(id) == 13) return muonID(idx, WW_medium_v2);
  return false;
}

bool WWAnalysis::isInSituFRLepton(int id, int idx, bool expt){
  if (abs(id) == 11) return electronID(idx, WW_fo_v2);
  else if (abs(id) == 13) return muonID(idx, WW_fo_v2);
  return false;
}

bool WWAnalysis::isDenominatorLepton(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_fo_v2);
  else if (abs(id) == 13) return muonID(idx, WW_fo_v2);
  return false;
}

bool WWAnalysis::isDenominatorLeptonNoIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_fo_noiso_v2);
  else if (abs(id) == 13) return muonID(idx, WW_fo_noiso_v2);
  return false;
}

bool WWAnalysis::isDenominatorLeptonIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_fo_v2);
  else if (abs(id) == 13) return muonID(idx, WW_fo_v2);
  return false;
}

bool WWAnalysis::isVetoLepton(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_veto_v2);
  else if (abs(id) == 13) return muonID(idx, WW_veto_v2);
  return false;
}

bool WWAnalysis::isVetoLeptonNoIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_veto_noiso_v2);
  else if (abs(id) == 13) return muonID(idx, WW_veto_noiso_v2);
  return false;
}

bool WWAnalysis::isVetoLeptonIso(int id, int idx){
  if (abs(id) == 11) return electronID(idx, WW_veto_v2);
  else if (abs(id) == 13) return muonID(idx, WW_veto_v2);
  return false;
}

bool WWAnalysis::hypsFromFirstGoodVertex(size_t hypIdx, float dz_cut){

  int lt_idx = hyp_lt_index()[hypIdx];
  int ll_idx = hyp_ll_index()[hypIdx];

  int lt_id = hyp_lt_id()[hypIdx];
  int ll_id = hyp_ll_id()[hypIdx];

  float lt_dz = abs(lt_id) == 11 ? tas::els_dzPV().at(lt_idx) : tas::mus_dzPV().at(lt_idx);
  float ll_dz = abs(ll_id) == 11 ? tas::els_dzPV().at(ll_idx) : tas::mus_dzPV().at(ll_idx);

  if (fabs(lt_dz) < dz_cut && fabs(ll_dz) < dz_cut) return true;    
  return false;
}

anal_type_t WWAnalysis::analysisCategory(float lep1pt, float lep2pt){
  if      (lep1pt > ptCutHigh && lep2pt > ptCutHigh) return HighHigh;
  else if (lep1pt > ptCutHigh && lep2pt > ptCutLow)  return HighLow;
  else if (lep2pt > ptCutHigh && lep1pt > ptCutLow)  return HighLow;
  else if (lep1pt > ptCutLow  && lep2pt > ptCutLow)  return LowLow;
  return Undefined;
}

int WWAnalysis::baselineRegion(int njets, int nbtags, float met, float ht, float lep1_pt, float lep2_pt){
  //Kinematic cuts
  //  if (lep1_pt < 20) return -1;
  //  if (lep2_pt < 15) return -1;
  //  if (met < 25) return -1;
  return njets;
}

int WWAnalysis::signalRegion(int njets, int nbtags, float met, float ht, float mt_min, float lep1pt, float lep2pt){
  
  cout << "Still defined in private code" << endl;
  return -1;
}

bool WWAnalysis::isGoodVetoElectronNoIso(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 7.) return false;
  if (!electronID(elidx, WW_veto_noiso_v2)) return false;
  return true;
}

bool WWAnalysis::isGoodVetoElectron(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 7.) return false;
  if (!electronID(elidx, WW_veto_v2)) return false;
  return true;
}

bool WWAnalysis::isFakableElectronNoIso(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 10.) return false;
  if (!electronID(elidx, WW_fo_noiso_v2)) return false;
  return true;
}

bool WWAnalysis::isFakableElectron(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 10.) return false;
  if (!electronID(elidx, WW_fo_v2)) return false;
  return true;
}

bool WWAnalysis::isGoodElectronNoIso(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 10.) return false;
  if (!electronID(elidx, WW_medium_noiso_v2)) return false;
  return true;
}

bool WWAnalysis::isGoodElectron(unsigned int elidx){
  if (els_p4().at(elidx).pt() < 10.) return false;
  if (!electronID(elidx, WW_medium_v2)) return false;
  return true;
}

bool WWAnalysis::isGoodVetoMuonNoIso(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 5.)         return false;
  if (!muonID(muidx, WW_veto_noiso_v2))     return false;
  return true;
}

bool WWAnalysis::isGoodVetoMuon(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 5.)         return false;
  if (!muonID(muidx, WW_veto_v2))           return false;
  return true;
}

bool WWAnalysis::isFakableMuonNoIso(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 10.)        return false;
  if (!muonID(muidx, WW_fo_noiso_v2))       return false;
  return true;
}

bool WWAnalysis::isFakableMuon(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 10.)        return false;
  if (!muonID(muidx, WW_fo_v2))             return false;
  return true;
}

bool WWAnalysis::isGoodMuonNoIso(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 10.)        return false;
  if (!muonID(muidx, WW_medium_noiso_v2))    return false;
  return true;
}

bool WWAnalysis::isGoodMuon(unsigned int muidx){
  if (mus_p4().at(muidx).pt() < 10.)        return false;
  if (!muonID(muidx, WW_medium_v2))          return false;
  return true;
}

int WWAnalysis::lepMotherID(Lep lep){
  if (tas::evt_isRealData()) return 1;
  if (abs(lep.pdgId()) != abs(lep.mc_id())) return 0; 
  else if (isFromZ(lep.pdgId(),lep.idx()) || isFromW(lep.pdgId(),lep.idx())){
    if (sgn(lep.pdgId()) == sgn(lep.mc_id())) return 1;
    else return 2;
  }
  else if (isFromB(lep.pdgId(),lep.idx())) return -1;
  else if (isFromC(lep.pdgId(),lep.idx())) return -2;
  return 0;
}

int WWAnalysis::lepMotherID_inSituFR(Lep lep){
  if (tas::evt_isRealData()) return 1;
  else if (isFromW(lep.pdgId(),lep.idx())){
    if (sgn(lep.pdgId()) == sgn(lep.mc_id())) return 1;
    else return 2;
  }
  else if (isFromB(lep.pdgId(),lep.idx())) return -1;
  else if (isFromC(lep.pdgId(),lep.idx())) return -2;
  else if (isFromLight(lep.pdgId(),lep.idx())) return -3;
  return 0;
}

int WWAnalysis::isGoodHyp(int iHyp, bool expt, bool verbose){

  //Bunch o' variables

  float pt_ll = tas::hyp_ll_p4().at(iHyp).pt(); 
  float pt_lt = tas::hyp_lt_p4().at(iHyp).pt(); 
  float eta_ll = tas::hyp_ll_p4().at(iHyp).eta();
  float eta_lt = tas::hyp_lt_p4().at(iHyp).eta();
  int idx_ll = tas::hyp_ll_index().at(iHyp);
  int idx_lt = tas::hyp_lt_index().at(iHyp);
  int id_ll = tas::hyp_ll_id().at(iHyp);
  int id_lt = tas::hyp_lt_id().at(iHyp);
  bool isos = false;
  if (sgn(id_ll) != sgn(id_lt)) isos = true;  
  bool passed_id_numer_ll = isGoodLepton(id_ll, idx_ll);
  bool passed_id_numer_lt = isGoodLepton(id_lt, idx_lt);
  bool passed_id_denom_ll = isDenominatorLepton(id_ll, idx_ll);
  bool passed_id_denom_lt = isDenominatorLepton(id_lt, idx_lt);
  bool passed_id_inSituFR_ll = isInSituFRLepton(id_ll, idx_ll);
  bool passed_id_inSituFR_lt = isInSituFRLepton(id_lt, idx_lt);
  bool extraZ = makesExtraZ(iHyp).result;
  bool extraGammaStar = makesExtraGammaStar(iHyp);

  //pass in situ ID
  bool truth_inSituFR = passed_id_inSituFR_lt && passed_id_inSituFR_ll;

  //Verbose info:
  if (verbose && pt_ll > ptCutLow && pt_lt > ptCutLow){
    cout << "hyp " << iHyp << " leptons: " << id_ll << " " << pt_ll << " (idx: " << idx_ll << ") " << id_lt << " " << pt_lt << " (idx: " << idx_lt << ")" << endl;
    cout << "   isos: " << isos << endl;
    cout << "   extraZ: " << extraZ << endl;
    cout << "   extraG: " << extraGammaStar << endl;
    cout << "   invt mass: " << (tas::hyp_ll_p4().at(iHyp) + tas::hyp_lt_p4().at(iHyp)).M() << endl;
    cout << "   passes eta: " << ((abs(id_ll) == 11 ? fabs(eta_ll) < 2.5 : fabs(eta_ll) < 2.4) && (abs(id_lt) == 11 ? fabs(eta_lt) < 2.5 : fabs(eta_lt) < 2.4)) << " etas are " << eta_ll << " and " << eta_lt << endl;
    cout << "   passes hypsFromFirstGoodVertex: " << hypsFromFirstGoodVertex(iHyp) << endl;
    cout << "   lepton with pT " << pt_ll << " passes id: " << passed_id_numer_ll << endl;
    cout << "   lepton with pT " << pt_lt << " passes id: " << passed_id_numer_lt << endl;
    cout << "   lowMassVeto: " << ((tas::hyp_ll_p4().at(iHyp) + tas::hyp_lt_p4().at(iHyp)).M() < 8) << endl;
    if (abs(id_ll) == 11) cout << "   lepton with pT " << pt_ll << " passes 3chg: " << threeChargeAgree(idx_ll) << endl;
    if (abs(id_lt) == 11) cout << "   lepton with pT " << pt_lt << " passes 3chg: " << threeChargeAgree(idx_lt) << endl;
  }

  //Kinematic Cuts
  if (pt_ll < 10) return 0;
  if (pt_lt < 10) return 0;
  if (abs(id_ll) == 11 && fabs(eta_ll) > 2.5) return 0;
  if (abs(id_lt) == 11 && fabs(eta_lt) > 2.5) return 0;
  if (abs(id_ll) == 13 && fabs(eta_ll) > 2.4) return 0;
  if (abs(id_lt) == 13 && fabs(eta_lt) > 2.4) return 0;

  //Must pass at least one denominator selection
  if (!passed_id_inSituFR_ll && !passed_id_denom_ll && !passed_id_numer_ll) return 0;
  if (!passed_id_inSituFR_lt && !passed_id_denom_lt && !passed_id_numer_lt) return 0;

  //Other cuts
  if ((tas::hyp_ll_p4().at(iHyp) + tas::hyp_lt_p4().at(iHyp)).M() < 8) return 0; 
  if (!hypsFromFirstGoodVertex(iHyp)) return 0;

  //Finished for events that fail z veto
  if (extraZ) return 6;
  else if (extraGammaStar) return 7;

  //Results
  else if (passed_id_numer_lt && passed_id_numer_ll && isos) return 4;  // 4 if both numer pass, WW
  else if (!passed_id_numer_lt && !passed_id_numer_ll && passed_id_denom_lt && passed_id_denom_ll && isos == true) return 1; // 1 SS, if both denom and no numer pass
  else if (((passed_id_numer_ll && passed_id_denom_lt) || (passed_id_numer_lt && passed_id_denom_ll)) && isos == true) return 2; //2 SS, one numer and one denom not numer
  else if (isos && truth_inSituFR) return 5;  // 5 if both pass inSituFR
  else if (passed_id_numer_lt && passed_id_numer_ll && !isos) return 3;  // 3 if both numer pass, SS
  return 0; //non-highpass SS
}

hyp_result_t WWAnalysis::chooseBestHyp(bool expt, bool verbose){

  //List of good hyps
  vector <int> good_hyps_ss; //same sign, tight tight
  vector <int> good_hyps_fr; //opposite sign, inSituFR
  vector <int> good_hyps_sf; //opposite sign, single fail
  vector <int> good_hyps_df; //opposite sign, double fail
  vector <int> good_hyps_os; //opposite sign, tight tight
  vector <int> good_hyps_zv; //same sign, tight tight, fail Z veto
  vector <int> good_hyps_gs; //same sign, tight tight, fail gamma* veto
  for (unsigned int i = 0; i < tas::hyp_type().size(); i++){
    int good_hyp_result = isGoodHyp(i, verbose);
    if (good_hyp_result == 4) good_hyps_os.push_back(i); 
    else if (good_hyp_result == 2) good_hyps_sf.push_back(i); 
    else if (good_hyp_result == 1) good_hyps_df.push_back(i); 
    else if (good_hyp_result == 3) good_hyps_ss.push_back(i); 
    else if (good_hyp_result == 5) good_hyps_fr.push_back(i); 
    else if (good_hyp_result == 6) good_hyps_zv.push_back(i); 
    else if (good_hyp_result == 7) good_hyps_gs.push_back(i); 
  }

  //hyp_class_ to track SS(3), SF(2), DF(1), OS(4), or none(0)
  int hyp_class_;

  //Load good hyps in, OS then SF then DF then SS
  vector <int> good_hyps;
  if (good_hyps_os.size() != 0){
     good_hyps = good_hyps_os;
     hyp_class_ = 4;
  }
  else if (good_hyps_zv.size() != 0){
    good_hyps = good_hyps_zv;
    hyp_class_ = 6;
  }
  else if (good_hyps_gs.size() != 0){
    good_hyps = good_hyps_zv;
    hyp_class_ = 7;
  }
  else if (good_hyps_sf.size() != 0){
    good_hyps = good_hyps_sf;
    hyp_class_ = 2;
  }
  else if (good_hyps_df.size() != 0){
     good_hyps = good_hyps_df;
     hyp_class_ = 1;
  }
  else if (good_hyps_ss.size() != 0){
    good_hyps = good_hyps_ss;
    hyp_class_ = 3;
  }
  else if (good_hyps_fr.size() != 0){
    good_hyps = good_hyps_fr;
    hyp_class_ = 5;
  }
  else hyp_class_ = 0; 

  //If no hyps or one hyps, know what to do
  int best_hyp_ = -1;
  if (good_hyps.size() == 1) best_hyp_ = good_hyps.at(0);

  //Otherwise, pick ones with more muons, then highest pT  
  if (good_hyps.size() > 1){
    best_hyp_ = good_hyps.at(0); 
    for (unsigned int i = 1; i < good_hyps.size(); i++){
      int hyp = good_hyps.at(i);
      if (tas::hyp_type().at(hyp) < tas::hyp_type().at(best_hyp_)) best_hyp_ = hyp;
      else if (tas::hyp_type().at(hyp) == tas::hyp_type().at(best_hyp_) && (tas::hyp_ll_p4().at(hyp).pt()+tas::hyp_lt_p4().at(hyp).pt()) > (tas::hyp_ll_p4().at(best_hyp_).pt() + tas::hyp_lt_p4().at(best_hyp_).pt())) best_hyp_ = hyp;
    }
  }

  if (best_hyp_ < 0){
    hyp_result_t null = { -1, -1 };
    return null;
  }

  hyp_result_t temp;
  temp.best_hyp = best_hyp_;
  temp.hyp_class = hyp_class_; 
  return temp;
}

vector <particle_t> WWAnalysis::getGenPair(bool verbose){

  vector <particle_t> gen_particles;
  if (tas::evt_isRealData()) return gen_particles;
  
  //First get all gen leptons 
  for (unsigned int gidx = 0; gidx < tas::genps_p4().size(); gidx++){
    if (tas::genps_status().at(gidx) != 1) continue;
    int id = tas::genps_id().at(gidx);
    if (abs(id) != 11 && abs(id) != 13 && abs(id) != 15) continue;
    float eta = tas::genps_p4().at(gidx).eta();
    if (fabs(eta) > 2.4) continue;
    int did = -1;
    int didx = -1;
    if (abs(id) == 15){  
      float dpt = -1;
      for (unsigned int didx_ = 0; didx_ < tas::genps_lepdaughter_id().at(gidx).size(); didx_++){
        int did_ = tas::genps_lepdaughter_id().at(gidx).at(didx_);
        if (abs(did_) != 11 && abs(did_) != 13) continue;
        if (fabs(tas::genps_lepdaughter_p4().at(gidx).at(didx_).eta()) > 2.4) continue;
        float dpt_ = tas::genps_lepdaughter_p4().at(gidx).at(didx_).pt();
        if (dpt_ <= dpt) continue;
        dpt = dpt_;
        didx = didx_;
        did = did_;
      }
    }
    particle_t temp;
    temp.id = abs(id) == 15 ? did : id;
    if (abs(id) == 15 && didx < 0) continue;
    temp.idx = abs(id) == 15 ? tas::genps_lepdaughter_idx().at(gidx).at(didx) : gidx;
    temp.p4 = abs(id) == 15 ? tas::genps_lepdaughter_p4().at(gidx).at(didx) : tas::genps_p4().at(gidx);
    if (temp.id != -1) gen_particles.push_back(temp);
  }
 
  if (gen_particles.size() < 2) return gen_particles;

  //Now loop over gen hyps and pick the best
  int gen_hyp_class = 5; 
  int type = 5; 
  particle_t lep1; 
  particle_t lep2;
  lep1.id = -1;
  lep2.id = -1;
  for (unsigned int idx1 = 0; idx1 < gen_particles.size(); idx1++){
    if (verbose) cout << "gen lep " << idx1 << ": " << gen_particles.at(idx1).id << " " << gen_particles.at(idx1).p4.pt() << endl;
    for (unsigned int idx2 = idx1+1; idx2 < gen_particles.size(); idx2++){
      int id1 = gen_particles.at(idx1).id;
      int id2 = gen_particles.at(idx2).id;
      int gen_hyp_class_ = 5;
      int type_ = 5;
      LorentzVector p41 = gen_particles.at(idx1).p4;
      LorentzVector p42 = gen_particles.at(idx2).p4;
      if (id1*id2 < 0) gen_hyp_class_ = 4; 
      if (id1*id2 > 0) gen_hyp_class_ = 3; 
      if (min(p41.pt(), p42.pt()) < 5) continue;
      if (max(p41.pt(), p42.pt()) < 5) continue;
      if (abs(id1) == 11 && abs(id2) == 11) type_ = 3; 
      if (abs(id1) == 13 && abs(id2) == 13) type_ = 0; 
      if ((abs(id1) == 13 && abs(id2) == 11) || (abs(id1) == 11 && abs(id2) == 13)) type_ = 1; 
      if (gen_hyp_class_ < gen_hyp_class){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      }
      else if (gen_hyp_class_ == gen_hyp_class && type_ < type){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      } 
      else if (gen_hyp_class_ == gen_hyp_class && type_ == type && p41.pt()+p42.pt() > lep1.p4.pt()+lep2.p4.pt()){
        gen_hyp_class = gen_hyp_class_;
        type = type_;
        lep1 = gen_particles.at(idx1);
        lep2 = gen_particles.at(idx2);
      }
    }
  }
  
  //Return leptons you found
  vector <particle_t> selected_gen_leptons;
  if (lep1.p4.pt() > lep2.p4.pt()){
    selected_gen_leptons.push_back(lep1);
    selected_gen_leptons.push_back(lep2);
  }
  else{ 
    selected_gen_leptons.push_back(lep2);
    selected_gen_leptons.push_back(lep1);
  }

  return selected_gen_leptons;

}

pair<particle_t, int> WWAnalysis::getThirdLepton(int hyp){

  //If hyp_class == 6, save the lepton that triggered the Z veto (so long as it is veto or higher)
  Z_result_t Zresult = makesExtraZ(hyp);
  if (Zresult.result == true){
    particle_t result;
    result.id = Zresult.id;
    result.idx = Zresult.idx;
    result.p4 = abs(result.id) == 11 ? tas::els_p4().at(result.idx) : tas::mus_p4().at(result.idx);
    int quality = 0;
    if (abs(result.id) == 11 ? !isGoodVetoElectron(result.idx) : !isGoodVetoMuon(result.idx)) quality = -1;
    if (abs(result.id) == 11 ? isFakableElectron(result.idx) : isFakableMuon(result.idx)) quality = 1;
    if (abs(result.id) == 11 ? isGoodElectron(result.idx) : isGoodMuon(result.idx)) quality = 2;
    if (quality >= 0) return pair<particle_t, int>(result, quality);
  }
  
  //Otherwise, find the highest-quality lepton possible. 

  //Selected Lepton Information
  int ll_id = tas::hyp_ll_id().at(hyp);
  int lt_id = tas::hyp_lt_id().at(hyp);
  unsigned int ll_idx = tas::hyp_ll_index().at(hyp);
  unsigned int lt_idx = tas::hyp_lt_index().at(hyp);

  //Store best lepton
  int lep3_id_ = -1;
  int lep3_idx_ = -1;
  int quality = 0;
  LorentzVector lep3_p4_; 

  //Electron Loop 
  for (unsigned int i = 0; i < tas::els_p4().size(); i++){

    //Remove electrons already selected
    if (abs(ll_id) == 11 && ll_idx == i) continue; 
    if (abs(lt_id) == 11 && lt_idx == i) continue; 

    //Remove electrons that fail kinematically
    if (tas::els_p4().at(i).pt() < 7) continue;
    if (fabs(tas::els_p4().at(i).eta()) > 2.5) continue;

    //Remove electrons that fail loosest ID, determine tighter IDs
    int quality_ = 0;
    if (!isGoodVetoElectron(i)) continue;
    if (isFakableElectron(i)) quality_ = 1;
    if (isGoodElectron(i)) quality_ = 2;

    //Choose the highest-quality, highest-pT electron 
    if (quality_ > quality || (quality_ == quality && tas::els_p4().at(i).pt() > lep3_p4_.pt())){
       quality = quality_;
       lep3_p4_ = tas::els_p4().at(i); 
       lep3_id_ = -11*tas::els_charge().at(i);
    } 
  }
  
  //Muon Loop
  for (unsigned int i = 0; i < tas::mus_p4().size(); i++){

    //Remove electrons already selected
    if (abs(ll_id) == 13 && ll_idx == i) continue; 
    if (abs(lt_id) == 13 && lt_idx == i) continue; 
   
    //Remove electrons that fail kinematically
    if (tas::mus_p4().at(i).pt() < 7) continue;
    if (fabs(tas::mus_p4().at(i).eta()) > 2.4) continue;

    //Remove muons that fail ID
    int quality_ = 0; 
    if (!isGoodVetoMuon(i)) continue;
    if (isFakableMuon(i)) quality_ = 1;
    if (isGoodMuon(i)) quality_ = 2;

    //Choose the highest-quality, highest-pT electron 
    if (quality_ > quality || (quality_ == quality && tas::mus_p4().at(i).pt() > lep3_p4_.pt())){
       quality = quality_;
       lep3_p4_ = tas::mus_p4().at(i); 
       lep3_id_ = -11*tas::mus_charge().at(i);
    } 

  }//Muon loop

  particle_t result;
  result.id = lep3_id_;
  result.p4 = lep3_p4_;
  result.idx = lep3_idx_;

  return pair<particle_t, int>(result, quality);

}

bool WWAnalysis::ptsort (int i,int j) { return (genps_p4()[i].pt()>genps_p4()[j].pt()); }

bool WWAnalysis::lepsort (Lep i,Lep j) { 
  if ( abs(i.pdgId())==abs(j.pdgId()) ) return ( i.pt()>j.pt() ); //sort by pt if same flavor
  else return ( abs(i.pdgId())>abs(j.pdgId()) ); //prefer muons over electrons, but check that mu have pt>25//fixme, need to sync // && i.pt()>ptCutHigh
}

bool WWAnalysis::jetptsort (Jet i,Jet j) { return (i.pt()>j.pt()); }

float WWAnalysis::coneCorrPt(int id, int idx){
  float miniIso = abs(id)==11 ? elMiniRelIsoCMS3_EA(idx, 1) : muMiniRelIsoCMS3_EA(idx, 1);
  LorentzVector lep_p4 = abs(id)==11 ? els_p4().at(idx) : mus_p4().at(idx);
  LorentzVector jet_p4  = closestJet(lep_p4, 0.4, 2.4, 1);
  float ptrel = ptRel(lep_p4, jet_p4, true);
  float A = abs(id)==11 ? 0.12 : 0.16;
  float B = abs(id)==11 ? 0.80 : 0.76;
  float C = abs(id)==11 ? 7.20 : 7.20;
  return ((ptrel > C) ? lep_p4.pt()*(1 + std::max((float)0, miniIso - A)) : std::max(lep_p4.pt(), jet_p4.pt() * B));
}

int WWAnalysis::getHighPtTriggerPrescale(LorentzVector& p4, int& idx, int& id)  {

  //   int temp;
  //   int retval = 0;
  //   
  //   if (id == 13) {
  //     
  //     setHLTBranch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"             , p4, temp);
  //     if (temp >= 0) retval += temp; 
  //     setHLTBranch("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" , tas::mus_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_MuonLeg().at(idx), temp);
  //     if (temp >= 0) retval += 2*temp; 
  //     setHLTBranch("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"  , tas::mus_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_MuonLeg().at(idx) , temp);
  //     if (temp >= 0) retval += 4*temp; 
  //     setHLTBranch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"           , p4, temp);
  //     if (temp >= 0) retval += 8*temp; 
  // 
  //   } else if (id == 11) {
  //     
  //     setHLTBranch("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"       ,  p4, temp);
  //     if (temp >= 0) retval += temp; 
  //     setHLTBranch("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" , tas::els_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_ElectronLeg().at(idx), temp);
  //     if (temp >= 0) retval += 2*temp; 
  //     setHLTBranch("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"  , tas::els_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_ElectronLeg().at(idx) , temp);
  //     if (temp >= 0) retval += 4*temp; 
  //     
  //   }
  //   
  //   return retval;
  return 1;

}

int WWAnalysis::getLowPtTriggerPrescale(LorentzVector& p4, int& idx, int& id) {

  //   int temp;
  //   int retval = 0;
  //   
  //   if (id == 13) {
  //     
  //     setHLTBranch("HLT_DoubleMu8_Mass8_PFHT300_v"                     , tas::mus_HLT_DoubleMu8_Mass8_PFHT300_MuonLeg().at(idx)                      , temp);
  //     if (temp >= 0) retval += temp; 
  //     setHLTBranch("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v"     , tas::mus_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_MuonLeg().at(idx)      , temp);
  //     if (temp >= 0) retval += 2*temp; 
  // 
  //   } else if (id == 11) {
  //     
  //     setHLTBranch("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v"   , tas::els_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg().at(idx), temp);
  //     if (temp >= 0) retval += temp; 
  //     setHLTBranch("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v"     , tas::els_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_ElectronLeg().at(idx)  , temp);
  //     if (temp >= 0) retval += 2*temp; 
  //     
  //   }
  // 
  //   return retval;

  return 1;
}

void WWAnalysis::setHLTBranch(const char* pattern, const LorentzVector& p4, int& HLTbranch) {
  TString name_HLT = triggerName(pattern);
  if (name_HLT=="TRIGGER_NOT_FOUND") {
    HLTbranch=0;
    return;
  }
  if (tas::passHLTTrigger(name_HLT)) {
    HLTbranch = HLT_prescale(name_HLT);
    if (passHLTTrigger(name_HLT,p4)==0) HLTbranch*=-1;
  } else HLTbranch = 0;
}

void WWAnalysis::setHLTBranch(const char* pattern, bool legMatch, int& HLTbranch) {
  TString name_HLT = triggerName(pattern);
  if (name_HLT=="TRIGGER_NOT_FOUND"){
    HLTbranch=0;
    return;
  }
  if (tas::passHLTTrigger(name_HLT)) {
    HLTbranch = HLT_prescale(name_HLT);
    if (legMatch==0) HLTbranch*=-1;
  } else HLTbranch = 0;
}
