#include <stdexcept>

#include "OSSelections.h"
#include "MuonSelections.h"
#include "ElectronSelections.h"
#include "PhotonSelections.h"

#include "Math/VectorUtil.h"

using namespace tas;

bool overlapMuon_ZMET_v1( int index , float ptcut ){
  for( unsigned int muind = 0; muind < cms3.mus_p4().size(); muind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.els_p4().at(index) , cms3.mus_p4().at(muind) );    
    if( dr > 0.05                                     ) continue;
    if( cms3.mus_p4().at(muind).pt() < ptcut          ) continue;
	if( !passMuonSelection_ZMET_v2(muind, true, true) ) continue;
    return true;
  }
  return false;
}

bool overlapElectron_ZMET_v1( int index , float ptcut ){
  for( unsigned int elind = 0; elind < cms3.els_p4().size(); elind++ ){
    float dr = ROOT::Math::VectorUtil::DeltaR( cms3.photons_p4().at(index) , cms3.els_p4().at(elind) );    
    if( dr > 0.2                                          ) continue;
    if( cms3.els_p4().at(elind).pt() < ptcut              ) continue;
	if( !passElectronSelection_ZMET_v2(elind, true, true) ) continue;
    return true;
  }
  return false;
}

//~-~-~-~-~-~-~-~-~-~//
//Electron selections//
//~-~-~-~-~-~-~-~-~-~//
bool passElectronSelection_ZMET(int index ){
  return passElectronSelection_ZMET_v4( index, true, true );
}

bool passElectronSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_tightMVA_v2 )       ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_tightMVA_v1 )       ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 10.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_v2 )          ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_NoIso_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapMuon_ZMET_v1( index, 15.0 )           ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_noiso_v2 )    ) return false; // Electron ID  
  return true;
}

bool passElectronSelection_ZMET_v1_NoIso(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMet2012_v1(index,10.0)          ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_noiso_v1 )    ) return false; // Electron ID  
  return true;
}


bool passElectronSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.els_p4().at(index).pt()) < 15.0    ) return false; // pT > 15 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.els_etaSC().at(index)) > 1.4442
	  && fabs(cms3.els_etaSC().at(index)) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.els_p4()[index].eta()) > 2.4    ) return false; // eta < 2.4
  // if( overlapMuon_ZMet2012_v1(index,10.0)          ) return false; // overlap removal
  if( !electronID( index, ZMET_loose_v1 )          ) return false; // Electron ID  
  return true;
}

//~-~-~-~-~-~-~-~//
//Muon selections//
//~-~-~-~-~-~-~-~//
bool passMuonSelection_ZMET(int index ){
  return passMuonSelection_ZMET_v4( index, true, true );
}

bool passMuonSelection_ZMET_v4(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v2 )              ) return false; // medium Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v3(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 15.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_mediumMu_v1 )              ) return false; // medium Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_v2 )                 ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_NoIso_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 15.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_noiso_v2 )           ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v1_NoIso(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_noiso_v1 )           ) return false; // tight Muon ID  
  return true;
}

bool passMuonSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.mus_p4().at(index).pt()) < 10.0       ) return false; // pT > 10 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.mus_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.mus_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.mus_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !muonID( index, ZMET_tight_v1 )                 ) return false; // tight Muon ID  
  return true;
}

bool passPhotonSelection_ZMET(int index ){
  return passPhotonSelection_ZMET_v2(index, true, true);
}

  bool passPhotonSelection_ZMET_v2(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 22 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( overlapElectron_ZMET_v1( index, 10.0 )              ) return false; // remove electrons from W
  if( !photonID(index, ZMET_photon_v1 )                   ) return false;
  return true;
}

bool passPhotonSelection_ZMET_v1(int index, bool vetoTransition, bool eta24 ){
  if( fabs(cms3.photons_p4().at(index).pt()) < 22.0       ) return false; // pT > 22 GeV - Minimum pT cut
  if( vetoTransition
	  && fabs(cms3.photons_p4().at(index).eta()) > 1.4442
	  && fabs(cms3.photons_p4().at(index).eta()) < 1.566  ) return false; // veto x-ition region
  if( eta24
	  && fabs(cms3.photons_p4().at(index).eta()) > 2.4    ) return false; // eta < 2.4
  if( !photonID(index, ZMET_photon_v1 ) ) return false;
  return true;
}
