#include "WcsAnalysis.h"
#include "classes/DelphesClasses.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"

#include<string>
#include<iostream>
#include<vector>

using std::vector;
using std::string;

std::vector<GenParticle*> WDaughters(TClones<GenParticle>& Particles)
{
  std::vector<GenParticle*> Ws;
  for (auto& p : Particles) {
    if (p.M1 < 0) continue;
    if (abs(p.PID) == 24) continue;
    if (p.PID == 21) continue; // ignore gluons
    if (p.PID == 22) continue; // ignore photon
    GenParticle *mom = Particles[p.M1];
    if (abs(mom->PID) == 24)
      Ws.push_back(&p);
  }
  return std::move(Ws);
}

Jet* closestJet(GenParticle& g, TClones<Jet>& Jets, double min_closest=0.4)
{
  double closest_dr = min_closest;
  Jet *best = 0;
  TLorentzVector l = g.P4();
  for (auto& j : Jets) {
    double dr = l.DeltaR(j.P4());
    if (dr < closest_dr) {
      best = &j;
      closest_dr = dr;
    }
  }
  return best;
}

int findKS_Gen(Jet* j, TClones<GenParticle>& ps)
{
  int idx = -1, bestx = -999.;
  auto consts = TRefs<GenParticle>(&j->Constituents);
  for (auto c : consts) {
    if (abs(c->PID) == pi_pdg && c->M1 >= 0 && abs(ps[c->M1]->PID) == ks_pdg) {
      int ks_idx = c->M1;
      float x = ps[ks_idx]->PT / j->PT;
      if (x > bestx) {
	bestx = x;
	idx = ks_idx;
      }
    }
  }

  return idx;
}


/// find true ks
int findKS_Rec(Jet* j, TClones<GenParticle>& ps)
{
  int idx = -1, bestx = -999.;
  auto consts = TRefs<Track>(&j->Constituents);
  for (auto cr : consts) {
    if (!cr) continue; // skip over towers
    auto c = (GenParticle*) cr->Particle.GetObject();
    if (abs(c->PID) == pi_pdg && c->M1 >= 0 && abs(ps[c->M1]->PID) == ks_pdg) {
      int ks_idx = c->M1;
      float x = ps[ks_idx]->PT / j->PT;
      if (x > bestx) {
	bestx = x;
	idx = ks_idx;
      }
    }
  }

  return idx;
}

vector<TLorentzVector> ksCands(Jet& j)
{
  std::vector<TLorentzVector> cands;
  auto consts = TRefs<Track>(&j.Constituents);
  for (auto pip : consts) {
    if (!pip) continue;
    if (pip->Charge < 0) continue;
    if (pip->PT < 1.) continue;
    if ((pip->D0/pip->ErrorD0) < 5) continue;
    for (auto pim : consts) {
      if (!pim) continue;
      if (pim->Charge > 0) continue;
      if (pim->PT < 1.) continue;
      if ((pim->D0/pim->ErrorD0) < 5) continue;
      TLorentzVector cand = pip->P4()+pim->P4();
      if (fabs(cand.M() - 0.5) < 0.2)
	cands.push_back(cand);
    }
  }
  return std::move(cands);
}

int main()
{
  gROOT->SetBatch(true);
  TClones<GenParticle> Particles;
  TClones<Jet> GenJets, Jets;

  string name = "tsW";
  TChain chain("Delphes");
  chain.Add((name+".root").c_str());
  Particles.FromTree(&chain, "Particle");
  GenJets.FromTree(&chain, "GenJet");
  Jets.FromTree(&chain, "Jet");

  vector<TH1F*> vhists;
#define makeTH1F(name, title, n, b, e) name(#name, title, n, b, e); vhists.push_back(&name);
  TH1F makeTH1F(genjetPT, ";p_{T} (GeV);", 100, 0, 250);
  TH1F makeTH1F(jetPT, ";p_{T} (GeV);", 100, 0, 250);
  TH1F makeTH1F(sPT, "p_{T} (GeV)", 100, 0, 250);

  TH1F makeTH1F(ksPt, "p_{T} (GeV)", 100, 0, 250);
  TH1F makeTH1F(ksM, "M (GeV)", 50, 0.25, 0.75);

  TH1F makeTH1F(sgx, "x", 100, 0, 1.5);
  TH1F makeTH1F(srx, "x", 100, 0, 1.5);

  TH1F makeTH1F(sgx_withks, "x", 100, 0, 1.5);
  TH1F makeTH1F(ksgx, "x", 100, 0, 1.2);

  TH1F makeTH1F(srx_withks, "x", 100, 0, 1.5);
  TH1F makeTH1F(ksrx, "x", 100, 0, 1.2);

  Long64_t numberOfEntries = chain.GetEntries(), SemiLepEvents=0;
  std::cout << numberOfEntries << std::endl;
  
  for (size_t ii = 0; ii < numberOfEntries; ++ii) {
    chain.GetEntry(ii);

    auto WDs = WDaughters(Particles);
    if (WDs.size() != 4) {
      std::cout << "ERROR :: W daughters found " << WDs.size() << "  : ";
      for (auto &d : WDs)
	std::cout << d->PID << " ";
      std::cout << std::endl;
      continue;
    }

    // only semileptonic events:
    bool WudEvent = std::any_of(WDs.begin(), WDs.end(), [](GenParticle* g) {return abs(g->PID) == 1;});
    bool WcsEvent = std::any_of(WDs.begin(), WDs.end(), [](GenParticle* g) {return abs(g->PID) == 3;});
    bool LepEvent = std::any_of(WDs.begin(), WDs.end(), [](GenParticle* g) {return abs(g->PID) == 11 || abs(g->PID) == 13;});
    if (!(WcsEvent || WudEvent) || !LepEvent) continue;
    SemiLepEvents++;

    GenParticle* Ws = 0;
    Jet *WsGen=0, *WsRec=0;
    if (WcsEvent) {
      Ws = *std::find_if(WDs.begin(), WDs.end(), [](GenParticle* g) {return abs(g->PID) == 3;});
      sPT.Fill(Ws->PT);
      WsGen = closestJet(*Ws, GenJets);
      WsRec = closestJet(*Ws, Jets);

      int gidx=-1, ridx=-1;

      if (WsGen) {
	gidx = findKS_Gen(WsGen, Particles);
	if (gidx >= 0) {
	  ksgx.Fill(Particles[gidx]->PT / WsGen->PT);
	  sgx_withks.Fill(WsGen->PT/Ws->PT);
	}
	sgx.Fill(WsGen->PT/Ws->PT);
      }
      if (WsRec) {
	srx.Fill(WsRec->PT/Ws->PT);
	auto kss = ksCands(*WsRec);
	if (kss.size() > 0) {
	  std::sort(kss.begin(), kss.end(), [](const TLorentzVector &i, const TLorentzVector &j){return i.Pt() > j.Pt();});
	  
	  ksPt.Fill(kss[0].Pt());
	  ksM.Fill(kss[0].M());
	}

	ridx = findKS_Rec(WsRec, Particles);
	if (ridx >= 0) {
	  ksgx.Fill(Particles[ridx]->PT / WsRec->PT);
	  sgx_withks.Fill(WsRec->PT / Ws->PT);
	}
      }
    }
    
    for (auto& g : GenJets) {
      genjetPT.Fill(g.PT);
    }

    for (auto& g : Jets) {
      jetPT.Fill(g.PT);
    }
  }

  system("mkdir -p plots/");
  TCanvas cvs;
  for (auto& h : vhists) {
    h->Draw();
    cvs.Print(("plots/"+name+"_"+string(h->GetName())+".png").c_str());
  }

  std::cout << "SemiLepFrac: " << SemiLepEvents / (float) numberOfEntries << std::endl;
}
