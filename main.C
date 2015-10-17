#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include <fstream>

using namespace std;

int ScanChain(TChain* chain, int nEvents, TString skimFilePrefix);

int main(int argc, char** argv) {

  gSystem->Exec("rm -rf baby.root");

  TString reader;
  ifstream fileList(argv[1]);
  
  TChain* ch = new TChain("Events");
  while (!fileList.eof()) {
    reader.ReadLine(fileList);
    ch->Add(reader);
  }
  ScanChain(ch, -1, argv[2]);
  return 0;
}
  
