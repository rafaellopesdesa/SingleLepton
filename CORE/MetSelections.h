#ifndef METSELECTIONS_H
#define METSELECTIONS_H

#include "CMS3.h"
#include "VertexSelections.h"

struct metStruct{
  metStruct() : met(-999.), metphi(-999.), metx(-999.), mety(-999.), sumet(-999.)  {}
  float met;
  float metphi;
  float metx;
  float mety;
  float sumet;
};

class FactorizedJetCorrector;

metStruct trackerMET(float deltaZCut = 0.2, const std::vector<LorentzVector>* jets = 0);
bool hbheNoiseFilter(int minZeros = 999);
bool hbheNoiseFilter_25ns(int minZeros = 999);
bool hbheIsoNoiseFilter();
pair <float, float> getT1CHSMET( FactorizedJetCorrector * jet_corrector );
pair <float, float> getT1CHSMET_fromMINIAOD( FactorizedJetCorrector * jet_corrector );
pair<float,float> MET3p0(); // tracker style MET with no particles with eta > 3.0
bool passesMETfilter(bool is25 = 1);
pair <float, float> getT1CHSMET3p0( FactorizedJetCorrector * jet_corrector );

std::pair<float,float> MET3p0();

#endif
