// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "Math/GenVector/VectorUtil.h"
#include "TMath.h"

// CMS3
#include "CORE/CMS3.h"
#include "CORE/TriggerSelections.h"
#include "CORE/ElectronSelections.h"
#include "CORE/MuonSelections.h"
#include "CORE/JetSelections.h"
#include "CORE/IsolationTools.h"

using namespace std;
using namespace tas;

int ScanChain(TChain* chain, int nEvents, TString skimFilePrefix) {

  // Example Histograms
  TFile* outputFile = TFile::Open("baby.root", "UPDATE");
  TTree* baby = new TTree(skimFilePrefix, skimFilePrefix+"_tree");

  float scale1fb;
  baby->Branch("scale1fb", &scale1fb);

  LorentzVector lep_p4;
  baby->Branch("lep_p4", &lep_p4);
  int lep_id;
  baby->Branch("lep_id", &lep_id);
  LorentzVector lep_mc_p4;
  baby->Branch("lep_mc_p4", &lep_mc_p4);
  int lep_mc_id;
  baby->Branch("lep_mc_id", &lep_mc_id);

  LorentzVector met;
  baby->Branch("met", &met);

  int pass_trigger;
  baby->Branch("pass_trigger", &pass_trigger);

  vector<LorentzVector> jets_p4;
  baby->Branch("jets_p4", &jets_p4);
  vector<bool> jets_btag;
  baby->Branch("jets_btag", &jets_btag);
  vector<int> jets_mc_parton;
  baby->Branch("jets_mc_parton", &jets_mc_parton);
  vector<LorentzVector> jets_mc_p4;
  baby->Branch("jets_mc_p4", &jets_mc_p4);
  int njets;
  baby->Branch("njets", &njets);
  int nbtags;
  baby->Branch("nbtags", &nbtags);

  vector<LorentzVector> gen_jets_p4;
  baby->Branch("gen_jets_p4", &gen_jets_p4);

  int gen_els_n;
  baby->Branch("gen_els_n", &gen_els_n);
  vector<LorentzVector> gen_els_p4;
  baby->Branch("gen_els_p4", &gen_els_p4);
  vector<LorentzVector> gen_els_reco_p4;
  baby->Branch("gen_els_reco_p4", &gen_els_reco_p4);
  vector<int> gen_els_mother_id;
  baby->Branch("gen_els_mother_id", &gen_els_mother_id);
  vector<int> gen_els_grandma_id;
  baby->Branch("gen_els_grandma_id", &gen_els_grandma_id);
  vector<bool> gen_els_reco_id;
  baby->Branch("gen_els_reco_id", &gen_els_reco_id);
  vector<bool> gen_els_reco_iso;
  baby->Branch("gen_els_reco_iso", &gen_els_reco_iso);

  int gen_mus_n;
  baby->Branch("gen_mus_n", &gen_mus_n);
  vector<LorentzVector> gen_mus_p4;
  baby->Branch("gen_mus_p4", &gen_mus_p4);
  vector<LorentzVector> gen_mus_reco_p4;
  baby->Branch("gen_mus_reco_p4", &gen_mus_reco_p4);
  vector<int> gen_mus_mother_id;
  baby->Branch("gen_mus_mother_id", &gen_mus_mother_id);
  vector<int> gen_mus_grandma_id;
  baby->Branch("gen_mus_grandma_id", &gen_mus_grandma_id);
  vector<bool> gen_mus_reco_id;
  baby->Branch("gen_mus_reco_id", &gen_mus_reco_id);
  vector<bool> gen_mus_reco_iso;
  baby->Branch("gen_mus_reco_iso", &gen_mus_reco_iso);

  int gen_taus_n;
  baby->Branch("gen_taus_n", &gen_taus_n);
  vector<LorentzVector> gen_taus_p4;
  baby->Branch("gen_taus_p4", &gen_taus_p4);
  vector<int> gen_taus_decaymode;
  baby->Branch("gen_taus_decaymode", &gen_taus_decaymode);
  vector<int> gen_taus_mother_id;
  baby->Branch("gen_taus_mother_id", &gen_taus_mother_id);
  vector<int> gen_taus_grandma_id;
  baby->Branch("gen_taus_grandma_id", &gen_taus_grandma_id);

  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) { cout << "I am done!" << endl; continue; }
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) {cout << "I am done!" << endl; continue;}
      tree->LoadTree(event);
      cms3.GetEntry(event);
      ++nEventsTotal;

      // reset
      gen_els_n = -1.;
      gen_mus_n = -1.;
      gen_taus_n = -1.;

      scale1fb = -1;
      
      lep_p4 = LorentzVector(0,0,0,0);
      lep_id = -1;
      lep_mc_p4 = LorentzVector(0,0,0,0);
      lep_mc_id = -1;

      met = LorentzVector(0,0,0,0);

      njets = -1;
      nbtags = -1;
      pass_trigger = -1;

      jets_p4.clear();
      jets_btag.clear();
      jets_mc_parton.clear();
      jets_mc_p4.clear();

      gen_jets_p4.clear();

      gen_els_p4.clear();
      gen_els_reco_p4.clear();
      gen_els_mother_id.clear();
      gen_els_grandma_id.clear();
      gen_els_reco_id.clear();
      gen_els_reco_iso.clear();

      gen_mus_p4.clear();
      gen_mus_reco_p4.clear();
      gen_mus_mother_id.clear();
      gen_mus_grandma_id.clear();
      gen_mus_reco_id.clear();
      gen_mus_reco_iso.clear();

      gen_taus_p4.clear();
      gen_taus_decaymode.clear();
      gen_taus_mother_id.clear();
      gen_taus_grandma_id.clear();

      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );

      // Basic selection
      int theElectronIdx = -1;
      int nElectron = 0;
      for (uint i=0; i<cms3.els_p4().size(); i++) {
	if (!isTightElectronPOGspring15_v1(i)) continue;
	if (cms3.els_p4().at(i).pt() < 25.) continue;
	if (TMath::Abs(cms3.els_p4().at(i).eta()) > 2.4) continue;
	theElectronIdx = i;
	nElectron++;
      }
      if (nElectron > 1) continue;

      int theMuonIdx = -1;
      int nMuon = 0;
      for (uint i=0; i<cms3.mus_p4().size(); i++) {
	if (!isMediumMuonPOG(i)) continue;
	if (TMath::Abs(mus_dxyPV().at(i)) > 0.05 ) continue;
	if (TMath::Abs(mus_dzPV().at(i)) > 0.1 ) continue;
	if (cms3.mus_p4().at(i).pt() < 25.) continue;
	if (TMath::Abs(cms3.mus_p4().at(i).eta()) > 2.4) continue;
	if (muRelIso04DB(i) > 0.12) continue;
	theMuonIdx = i;
	nMuon++;
      }
      if (nMuon > 1) continue;
      if (!((nMuon == 0 && nElectron == 1) || (nMuon == 1 && nElectron == 0))) continue;
      vector<int> theJetsIdx;
      for (uint i=0; i<cms3.pfjets_p4().size(); i++) {
	if (!isLoosePFJet_50nsV1(i)) continue;
	bool jet_is_lep = false;
	for (uint j=0; j<cms3.els_p4().size(); j++) {
	  if (ROOT::Math::VectorUtil::DeltaR(cms3.els_p4().at(j), cms3.pfjets_p4().at(i)) < 0.3) jet_is_lep = true;
	} 
	for (uint j=0 ; j<cms3.mus_p4().size(); j++) {
	  if (ROOT::Math::VectorUtil::DeltaR(cms3.mus_p4().at(j), cms3.pfjets_p4().at(i)) < 0.3) jet_is_lep = true;
	} 
	if (jet_is_lep) continue;
	if (cms3.pfjets_p4().at(i).pt() < 20.) continue;
	if (TMath::Abs(cms3.pfjets_p4().at(i).eta()) > 2.5) continue; 
	theJetsIdx.push_back(i);
      }
      if (theJetsIdx.size() < 3) continue;
                
      // Analysis Code
      scale1fb = cms3.evt_scale1fb();

      if (nMuon == 1) {
	lep_p4 = cms3.mus_p4().at(theMuonIdx);
	lep_id = -cms3.mus_charge().at(theMuonIdx)*13;
	lep_mc_p4 = cms3.mus_mc_p4().at(theMuonIdx);
	lep_mc_id = cms3.mus_mc_id().at(theMuonIdx);
	setHLTBranch("HLT_IsoMu20_v", lep_p4, pass_trigger);
      } else if (nElectron == 1) {
	lep_p4 = cms3.els_p4().at(theElectronIdx);
	lep_id = -cms3.els_charge().at(theElectronIdx)*11;
	lep_mc_p4 = cms3.els_mc_p4().at(theElectronIdx);
	lep_mc_id = cms3.els_mc_id().at(theElectronIdx);
	setHLTBranch("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v", lep_p4, pass_trigger);
      }
      
      met = LorentzVector(cms3.evt_pfmet()*TMath::Cos(cms3.evt_pfmetPhi()), cms3.evt_pfmet()*TMath::Sin(cms3.evt_pfmetPhi()), 0., cms3.evt_pfmet());
      
      njets = theJetsIdx.size();
      nbtags = 0;
      for (auto idx : theJetsIdx) {
	jets_p4.push_back(cms3.pfjets_p4().at(idx));
	jets_btag.push_back(cms3.pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag().at(idx) > 0.890);
	if (cms3.pfjets_pfCombinedInclusiveSecondaryVertexV2BJetTag().at(idx) > 0.890) nbtags++;
	jets_mc_parton.push_back(cms3.pfjets_partonFlavour().at(idx));
	jets_mc_p4.push_back(cms3.pfjets_mc_p4().at(idx));
      }

      
      gen_jets_p4 = cms3.genjets_p4NoMuNoNu();

      for (int i=0; i<cms3.genps_id().size(); i++) {
	
		  
	if (cms3.genps_fromHardProcessBeforeFSR().at(i) && TMath::Abs(cms3.genps_id().at(i)) == 11) {
	  gen_els_p4.push_back(cms3.genps_p4().at(i));
	  gen_els_mother_id.push_back(cms3.genps_id_simplemother().at(i));
	  gen_els_grandma_id.push_back(cms3.genps_id_simplegrandma().at(i));

	  int reco_idx = -1;
	  float min_dr = 9999999.;
	  for (uint j=0; j<cms3.els_p4().size(); j++) {
	    float DR = ROOT::Math::VectorUtil::DeltaR(cms3.els_p4().at(j), cms3.genps_p4().at(i));
	    if (DR < min_dr) {
	      min_dr = DR;
	      reco_idx = j;
	    }
	  }
	  if (min_dr > 0.1) reco_idx = -1;

	  if (reco_idx > -1) {
	    gen_els_reco_p4.push_back(cms3.els_p4().at(reco_idx));
	    gen_els_reco_id.push_back(isTightElectronPOGspring15noIso_v1(reco_idx));
	    gen_els_reco_iso.push_back(isTightElectronPOGspring15_v1(reco_idx));
	  } else {
	    gen_els_reco_p4.push_back(LorentzVector(0,0,0,0));
	    gen_els_reco_id.push_back(false);
	    gen_els_reco_iso.push_back(false);
	  }	  
	  gen_els_n = gen_els_p4.size();
	}

	else if (cms3.genps_fromHardProcessBeforeFSR().at(i) && TMath::Abs(cms3.genps_id().at(i)) == 13) {
	  gen_mus_p4.push_back(cms3.genps_p4().at(i));
	  gen_mus_mother_id.push_back(cms3.genps_id_simplemother().at(i));
	  gen_mus_grandma_id.push_back(cms3.genps_id_simplegrandma().at(i));

	  int reco_idx = -1;
	  float min_dr = 9999999.;
	  for (uint j=0; j<cms3.mus_p4().size(); j++) {
	    float DR = ROOT::Math::VectorUtil::DeltaR(cms3.mus_p4().at(j), cms3.genps_p4().at(i));
	    if (DR < min_dr) {
	      min_dr = DR;
	      reco_idx = j;
	    }
	  }
	  if (min_dr > 0.1) reco_idx = -1;

	  if (reco_idx > -1) {
	    gen_mus_reco_p4.push_back(cms3.mus_p4().at(reco_idx));
	    gen_mus_reco_id.push_back(isTightMuonPOG(reco_idx) && TMath::Abs(mus_dxyPV().at(reco_idx)) < 0.05 && TMath::Abs(mus_dzPV().at(reco_idx)) < 0.1);
	    gen_mus_reco_iso.push_back(isTightMuonPOG(reco_idx) && TMath::Abs(mus_dxyPV().at(reco_idx)) < 0.05 && TMath::Abs(mus_dzPV().at(reco_idx)) < 0.1 && muRelIso04DB(reco_idx) > 0.12);
	  } else {
	    gen_mus_reco_p4.push_back(LorentzVector(0,0,0,0));
	    gen_mus_reco_id.push_back(false);
	    gen_mus_reco_iso.push_back(false);
	  }	  
	  gen_mus_n = gen_mus_p4.size();
	}

	else if (cms3.genps_fromHardProcessBeforeFSR().at(i) && TMath::Abs(cms3.genps_id().at(i)) == 15) {

	  gen_taus_p4.push_back(cms3.genps_p4().at(i));
          gen_taus_mother_id.push_back(cms3.genps_id_simplemother().at(i));
          gen_taus_grandma_id.push_back(cms3.genps_id_simplegrandma().at(i));

	  int lastcopy_idx = i;
	  bool isLastCopy = cms3.genps_isLastCopy().at(i);
	  while (!isLastCopy) {
	    for (int j=0; j<cms3.genps_id().size(); j++) {
	      if (cms3.genps_id().at(j) == cms3.genps_id().at(lastcopy_idx) && cms3.genps_idx_simplemother().at(j) == lastcopy_idx) {
		lastcopy_idx = j;
		isLastCopy = cms3.genps_isLastCopy().at(j);
		break;
	      }
	    }
	  }

	  int thisDecay = 0;
	  for (int j=0; j<cms3.genps_id().size(); j++) {
	    if (!cms3.genps_isDirectHardProcessTauDecayProductFinalState().at(j)) continue;
	    if (cms3.genps_idx_simplemother().at(j) != lastcopy_idx) continue;
	    if (cms3.genps_charge().at(j) == 0) continue;
	    if (TMath::Abs(cms3.genps_id().at(j)) == 11 || TMath::Abs(cms3.genps_id().at(j)) == 13) {
	      thisDecay = cms3.genps_id().at(j);
	      break;
	    } else {
	      thisDecay++;
	    }
	  }
	  gen_taus_decaymode.push_back(thisDecay);
	  gen_taus_n = gen_taus_p4.size();
	}
      }
      baby->Fill();
    }
  
  }
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }

  outputFile->cd();
  baby->Write();
  outputFile->Close();

  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << endl;
  
  // return
  return 0;
}
