#ifndef WWSELECTIONS_H
#define WWSELECTIONS_H
#include "CMS3.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"
#include "MetSelections.h"
#include "JetSelections.h"
#include "VertexSelections.h"
#include "TriggerSelections.h"
#include "MCSelections.h"
#include "IsolationTools.h"
#include "Math/VectorUtil.h"

namespace WWAnalysis {
  
const static float ptCutHigh = 25.;
const static float ptCutLow = 10.;

//Enums
enum anal_type_t { HighHigh = 0, HighLow = 1, LowLow = 2, Undefined = -1 };
enum hyp_type_t { EE, MM, EM, UNASSIGNED }; 
enum est_type_t { MCWW, FLIPMC, SFAKEMC, DFAKEMC, FLIPDD, SFAKEDD, DFAKEDD, ALLMC }; 
//fixme: put WF and FSR in different categories
enum LeptonCategories { Prompt = 0, PromptWS = 1, PromptWF = 2, PromptFSR = 2, 
			FakeLightTrue = 3, FakeC = 4, FakeB = 5, FakeLightFake = 6, FakeHiPtGamma = 7, 
			FakeUnknown = 8, FakeLowPtGamma = 9, All9999 = 10,
			Other = 11, End = 12};

//Structs
struct hyp_result_t { int best_hyp; int hyp_class; };
struct particle_t { int id; LorentzVector p4; int idx; };
struct Lep;
struct DilepHyp;
struct Jet;
struct Z_result_t { bool result; int id; int idx; }; 

//helper function for sign
template <typename T> int sgn(T val){
  return (T(0) < val) - (val < T(0));
}

float coneCorrPt(int id, int idx);

//Main Object selections
 bool isGoodLepton(int id, int idx);
 bool isDenominatorLepton(int id, int idx);
 bool isVetoLepton(int id, int idx);

//Hyp selections
hyp_result_t chooseBestHyp(bool expt, bool verbose=false);
int isGoodHyp(int iHyp, bool expt, bool verbose=false);
bool makesExtraGammaStar(int iHyp);
Z_result_t makesExtraZ(int iHyp);
bool hypsFromFirstGoodVertex(size_t hypIdx, float dz_cut = 1.0);
std::pair<particle_t, int> getThirdLepton(int hyp);
std::vector<particle_t> getGenPair(bool verbose=false);

//Signal region selections
anal_type_t analysisCategory(float lep1pt, float lep2pt);
int baselineRegion(int njets, int nbtags, float met, float ht, float lep1_pt, float lep2_pt);
int signalRegion(int njets, int nbtags, float met, float ht, float mt_min, float lep1pt, float lep2pt);

//More Lepton selections
bool isGoodLeptonNoIso(int id, int idx);
bool isGoodLeptonIso(int id, int idx);
bool isDenominatorLeptonNoIso(int id, int idx);
bool isDenominatorLeptonIso(int id, int idx);
bool isVetoLeptonNoIso(int id, int idx);
bool isVetoLeptonIso(int id, int idx);
bool isGoodVetoElectronNoIso(unsigned int elidx);
bool isGoodVetoElectron(unsigned int elidx);
bool isFakableElectronNoIso(unsigned int elidx);
bool isFakableElectron(unsigned int elidx);
bool isGoodElectronNoIso(unsigned int elidx);
bool isGoodElectron(unsigned int elidx);
bool isGoodVetoMuonNoIso(unsigned int muidx);
bool isGoodVetoMuon(unsigned int muidx);
bool isFakableMuonNoIso(unsigned int muidx);
bool isFakableMuon(unsigned int muidx);
bool isGoodMuonNoIso(unsigned int muidx);
bool isGoodMuon(unsigned int muidx);
bool isLooseIsolatedLepton(int id, int idx);
bool isIsolatedLepton(int id, int idx);
bool isInSituFRLepton(int lep_id, int lep_idx, bool expt = false);

//MC truth functions
int lepMotherID(Lep lep);
int lepMotherID_inSituFR(Lep lep);

//Jet selection function
 std::pair <vector <Jet>, vector <Jet> > WWJetsCalculator(std::vector<LorentzVector> JetCollection);

// Calculate generator ht
float getGenHT(bool is_b_a_jet = true);

//Sorting functions
bool ptsort (int i,int j);
bool lepsort (Lep i,Lep j);
bool jetptsort (Jet i,Jet j);

//Definition of structs
struct Lep {
  Lep(int pdgid, int idxx):pdgid_(pdgid),idx_(idxx){}
  int charge() {return -1*pdgid_/abs(pdgid_);}
  int pdgId() {return pdgid_;}
  int idx() {return idx_;}
  float pt() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_).pt() : cms3.mus_p4().at(idx_).pt();}
  float eta() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_).eta() : cms3.mus_p4().at(idx_).eta();}
  LorentzVector p4() {return abs(pdgid_)==11 ? cms3.els_p4().at(idx_) : cms3.mus_p4().at(idx_);}
  float relIso03() { return abs(pdgid_)==11 ? eleRelIso03(idx_, WW) : muRelIso03(idx_, WW);}
  float miniRelIso() { return abs(pdgid_)==11 ? elMiniRelIsoCMS3_EA(idx_) : muMiniRelIsoCMS3_DB(idx_);}
  float dxyPV() { return abs(pdgid_)==11 ? cms3.els_dxyPV().at(idx_) : cms3.mus_dxyPV().at(idx_);}
  float dzPV() { return abs(pdgid_)==11 ? cms3.els_dzPV().at(idx_) : cms3.mus_dzPV().at(idx_);}
  float d0Err() { return abs(pdgid_)==11 ? cms3.els_d0Err().at(idx_) : cms3.mus_d0Err().at(idx_);}
  float ip3d() { return abs(pdgid_)==11 ? cms3.els_ip3d().at(idx_) : cms3.mus_ip3d().at(idx_);}
  float ip3dErr() { return abs(pdgid_)==11 ? cms3.els_ip3derr().at(idx_) : cms3.mus_ip3derr().at(idx_);}
  int mc3_id() {return abs(pdgid_)==11 ? cms3.els_mc3_id().at(idx_) : cms3.mus_mc3_id().at(idx_);}
  int mc3idx() {return abs(pdgid_)==11 ? cms3.els_mc3idx().at(idx_) : cms3.mus_mc3idx().at(idx_);}
  int mc3_motherid() {return abs(pdgid_)==11 ? cms3.els_mc3_motherid().at(idx_) : cms3.mus_mc3_motherid().at(idx_);}
  int mc3_motheridx() {return abs(pdgid_)==11 ? cms3.els_mc3_motheridx().at(idx_) : cms3.mus_mc3_motheridx().at(idx_);}
  int mc_id() { return abs(pdgid_)==11 ? cms3.els_mc_id().at(idx_) : cms3.mus_mc_id().at(idx_);}
  int mcidx() { return abs(pdgid_)==11 ? cms3.els_mcidx().at(idx_) : cms3.mus_mcidx().at(idx_);}
  int mc_motherid() {return abs(pdgid_)==11 ? cms3.els_mc_motherid().at(idx_) : cms3.mus_mc_motherid().at(idx_);}
  LorentzVector mc_p4() { return abs(pdgid_)==11 ? cms3.els_mc_p4().at(idx_) : cms3.mus_mc_p4().at(idx_);}
private:
  int pdgid_, idx_;
};

struct DilepHyp {
  DilepHyp(Lep lepone_, Lep leptwo_):
    leadlep(lepone_),trailep(leptwo_) {
    if (lepone_.pt()<leptwo_.pt()) {
      leadlep = leptwo_;
      trailep = lepone_;
    }
  }
  int charge() {return leadlep.charge()+trailep.charge();}
  LorentzVector p4() {return leadlep.p4()+trailep.p4();}
  Lep leadLep() {return leadlep;}
  Lep traiLep() {return trailep;}
private:
  Lep leadlep, trailep;
};

struct Jet {
Jet(int idxx,float __CSVv2,float __CSVsl,float __CSVtche):idx_(idxx),_CSVv2(__CSVv2),_CSVsl(__CSVsl),_CSVtche(__CSVtche) {}
  LorentzVector p4() {return cms3.pfjets_p4()[idx_];}
  float pt() {return p4().pt();}
  float eta() {return p4().eta();}
  float phi() {return p4().phi();}
  int   mc3_id() {return cms3.pfjets_mc3_id()[idx_];}
  LorentzVector genjet_p4() {return cms3.pfjets_mc_p4()[idx_];}
  LorentzVector genps_p4() {return cms3.pfjets_mc_gp_p4()[idx_];}
  float pileup_jet_id() {return cms3.pfjets_pileupJetId()[idx_];}
  int parton_flavor() {return cms3.pfjets_partonFlavour()[idx_];}
  float CSVv2() {return _CSVv2;}
  float CSVsl() {return _CSVsl;}
  float CSVtche() {return _CSVtche;}
  int idx() {return idx_;}
private:
  int idx_;
  float _CSVv2,_CSVsl,_CSVtche;
};

 int convertCMS3tag(TString tagName) ;
 int getHighPtTriggerPrescale(LorentzVector& p4, int& idx, int& id) ;
 int getLowPtTriggerPrescale(LorentzVector& p4, int& idx, int& id) ;
 void setHLTBranch(const char* pattern, const LorentzVector& p4, int& HLTbranch);
 void setHLTBranch(const char* pattern, bool legMatch, int& HLTbranch);

}


#endif
