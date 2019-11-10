#include "ROOT/RVec.hxx"

#ifndef TOPSELECTIONS__
#define TOPSELECTIONS__

float DeltaPhi(float p1, float p2);
float DeltaR(float e1, float p1, float e2, float p2);
using namespace ROOT::VecOps;
RVec<int> fillOverlap(const RVec<float> &jeteta, const RVec<float> &jetphi, const RVec<float> &lepeta, const RVec<float> &lepphi);

#endif

#ifndef HADRONSELECTIONS__
#define HADRONSELECTIONS__

RVec<int> jetHadron(const RVec<float> &hadeta, const RVec<float> &hadphi, const RVec<float> &jeteta, const RVec<float> &jetphi);
RVec<float> GetX(const RVec<float> &hadpt, const RVec<int> &jeti, const RVec<float> &jetpt);
RVec<float> GetHadJetDR(const RVec<float> &hadeta, const RVec<float> &hadphi, const RVec<int> &jeti, const RVec<float> &jeteta, const RVec<float> &jetphi);

#endif

#if defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ function jetHadron;
#pragma link C++ function GetX;
#pragma link C++ function GetHadJetDR;
#pragma link C++ function fillOverlap;
#endif
