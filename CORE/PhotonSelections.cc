#include <stdexcept>

#include "PhotonSelections.h"
#include "IsolationTools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

using namespace tas;

bool isLoosePhoton(unsigned int phIdx, analysis_t analysis, int version){
  // This is based on POG LOOSE, with an additional sieie sideband
  if (analysis == HAD){
    if (version == 1){
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.0559) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.049) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    }
    if (version == 2){
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.048) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.069) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    } 
    if (version == 3){ // Based on 50ns Spring15 IDs
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.015) return false; // looser than POG (0.01) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.05) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.035) return false; // larger than POG (0.0321) to get a sideband
        if(photons_full5x5_hOverEtowBC().at(phIdx) >= 0.05) return false; 
        if(photons_haspixelSeed().at(phIdx)) return false;
        return true;
      } 
      else return false;
    } 
  }
  else {
    std::cerr << "CORE::PhotonSelection -- Photon selection not implemented for analysis " << analysis << std::endl;
    std::cerr << "CORE::PhotonSelection -- I will return true." << std::endl;
  }

  return true;
}

bool isTightPhoton(unsigned int phIdx, analysis_t analysis, int version){
  // This should be the same as POG LOOSE: 
  //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  if (analysis == HAD){
    if (version == 1){
      if (!isLoosePhoton(phIdx, HAD)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.010) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0321) return false; 
        return true;
      } 
      else return false;
    }
    if (version == 2){
      if (!isLoosePhoton(phIdx, HAD, 2)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0106) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0266) return false; 
        return true;
      } 
      else return false;
    } 
    if (version == 3){
      if (!isLoosePhoton(phIdx, HAD, 3)) return false;
      if( fabs(photons_p4().at(phIdx).eta()) <= 1.479 ){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0103) return false; 
        return true;
      } 
      else if( fabs(photons_p4().at(phIdx).eta()) > 1.479 && fabs(photons_p4().at(phIdx).eta())  < 2.5){
        if(photons_full5x5_sigmaIEtaIEta().at(phIdx) >= 0.0277) return false; 
        return true;
      } 
      else return false;
    } 
  }
  else {
    std::cerr << "CORE::PhotonSelection -- Photon selection not implemented for analysis " << analysis << std::endl;
    std::cerr << "CORE::PhotonSelection -- I will return true." << std::endl;
  } 

  return true;
}

bool photonID(unsigned int phIdx, id_level_t id_level){
    
  // analysis_t analysis = whichAnalysis(id_level);

  switch (id_level){

	///////////////////
	// MET Templates //
	///////////////////
  
  case(ZMET_photon_v2):

	if( !isTemplatePhoton( phIdx ) ) return false;
	if( !isLoosePhoton_Spring15_50ns( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;
  
  case(ZMET_photon_v1):

	if( !isTemplatePhoton( phIdx ) ) return false;
	// Following cuts done in analysis code.
	// match to pfjet w/ neutral EM fraction > 70%
	// reject photons within electron with pT > 10 within cone of dR < 0.2
	// Reject photons aligned to MET within 0.14 radians in phi
	else return true;
	break;

	///////////////
	/// Default ///
	///////////////

  default:
	{
	  cout << "Warning! Photon ID not defined for this id_level!" << endl;
	  return false;
	}

  }//switch
  return true;
}

bool isTemplatePhoton( unsigned int phIdx ) {
  if( fabs(photons_p4()            .at(phIdx).pt())  <= 22.0 ) return false; // photon pT > 22 GeV
  if( fabs(photons_p4()            .at(phIdx).eta()) >= 2.4  ) return false; // photon eta < 2.4
  if(photons_full5x5_hOverEtowBC() .at(phIdx)        >= 0.1  ) return false; // hOverE < 0.1
  if(photons_haspixelSeed()        .at(phIdx)                ) return false; // veto pixel seed
  return true;
}

// from https://indico.cern.ch/event/369239/contribution/2/attachments/1134693/1623149/spring15_pcb.pdf
bool isLoosePhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0103                       ) return false;
	if( chiso  > 2.44                         ) return false;
	if( nhiso  > 2.57 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 1.92 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0277                       ) return false;
	if( chiso  > 1.84                         ) return false;
	if( nhiso  > 4.00 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 2.15 + 0.0041*pt             ) return false;
  }
  
  return true;
}

bool isMediumPhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0100                       ) return false;
	if( chiso  > 1.31                         ) return false;
	if( nhiso  > 0.60 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 1.33 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0267                       ) return false;
	if( chiso  > 1.25                         ) return false;
	if( nhiso  > 1.65 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 1.02 + 0.0041*pt             ) return false;
  }
  
  return true;
}

bool isTightPhoton_Spring15_50ns( int photonIdx )
{
  
  float eta = -999;  
  try{
	eta = cms3.photons_p4().at(photonIdx).eta();
  }
  catch( exception &e ){
	std::cout<<"Error! no photon with photonIdx: "<<photonIdx<<std::endl;
	return false;
  }

  float chiso = photonCHIso03EA(photonIdx);
  float nhiso = photonNHIso03EA(photonIdx);
  float emiso = photonEMIso03EA(photonIdx);

  float sieie  = cms3.photons_full5x5_sigmaIEtaIEta().at(photonIdx);
  float hovere = cms3.photons_full5x5_hOverE()       .at(photonIdx);
  float pt     = cms3.photons_p4()                   .at(photonIdx).pt();
  
  if(       abs(eta) < 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0010                       ) return false;
	if( chiso  > 0.91                         ) return false;
	if( nhiso  > 0.33 + exp(0.0044*pt+0.5809) ) return false;
	if( emiso  > 0.61 + 0.0043*pt             ) return false;
  }else if( abs(eta) > 1.479 ){
	if( hovere > 0.05                         ) return false;
	if( sieie  > 0.0267                       ) return false;
	if( chiso  > 0.65                         ) return false;
	if( nhiso  > 0.93 + exp(0.0040*pt+0.9402) ) return false;
	if( emiso  > 0.54 + 0.0041*pt             ) return false;
  }
  
  return true;
}
