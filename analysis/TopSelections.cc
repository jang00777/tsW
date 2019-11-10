#include<cmath>
#include "ROOT/RVec.hxx"

#ifndef TOPSELECTIONS__
#define TOPSELECTIONS__

float DeltaPhi(float p1, float p2) {
  float dphi = p2-p1;
  if ( dphi > M_PI ) {
    dphi -= 2.0*M_PI;
  } else if ( dphi <= -M_PI ) {
    dphi += 2.0*M_PI;
  }
  return dphi;
}

float DeltaR(float e1, float p1, float e2, float p2) {
  float dphi = DeltaPhi(p1,p2);
  float deta = e2 - e1;
  return std::sqrt( dphi*dphi + deta*deta );
}

using namespace ROOT::VecOps;
RVec<int> fillOverlap(const RVec<float> &jeteta, const RVec<float> &jetphi, const RVec<float> &lepeta, const RVec<float> &lepphi)
{
  RVec<int> result;
  for (int i = 0; i < jeteta.size(); ++i) {
    int overlap = 0;
    for (int j = 0; j < lepeta.size(); ++j) {
      if (DeltaR(jeteta[i], jetphi[i], lepeta[j], lepphi[j]) < 0.4) overlap = 1;
    }
    result.push_back(overlap);
  }
  return result;
};

#endif
#ifndef HADRONSELECTIONS__
#define HADRONSELECTIONS__

RVec<int> jetHadron(const RVec<float> &hadeta, const RVec<float> &hadphi, const RVec<float> &jeteta, const RVec<float> &jetphi)
{
  RVec<int> result;
  for (int i = 0; i < hadeta.size(); ++i) {
    int overlap = -1;
    for (int j = 0; j < jeteta.size(); ++j) {
      if (DeltaR(jeteta[j], jetphi[j], hadeta[i], hadphi[i]) < 0.4) overlap = j;
    }
    result.push_back(overlap);
  }
  return result;
};

RVec<float> GetX(const RVec<float> &hadpt, const RVec<int> &jeti, const RVec<float> &jetpt)
{
  RVec<float> result;
  for (int i = 0; i < hadpt.size(); ++i) {
    result.push_back(hadpt[i] / jetpt[jeti[i]]);
  }
  return result;
};

RVec<float> GetHadJetDR(const RVec<float> &hadeta, const RVec<float> &hadphi, const RVec<int> &jeti, const RVec<float> &jeteta, const RVec<float> &jetphi)
{
  RVec<float> result;
  for (int i = 0; i < hadeta.size(); ++i) {
    result.push_back(DeltaR(hadeta[i], hadphi[i], jeteta[jeti[i]], jetphi[jeti[i]]));
  }
  return result;
};

#endif
