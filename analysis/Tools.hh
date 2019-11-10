#include "Math/Vector4D.h"

using ROOT::Math::PtEtaPhiMVector;

#include "external/computePtEtaTable.h"
extern computePtEtaTable trigSFEl;
extern computePtEtaTable elecSF_;
extern computePtEtaTable elecSFRec_;

extern computePtEtaTable trigSFMu;
extern computePtEtaTable muonSF_;
extern computePtEtaTable muonSFIso_;

Bool_t LumiCheck(UInt_t run, UInt_t luminosityBlock);
float PileupGetWeight(int nTrueInt, int sys=0);

#if defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ function LumiCheck;
#pragma link C++ function PileupGetWeight;
#pragma link C++ global muonSF_;
#pragma link C++ global elecSF_;
#endif
