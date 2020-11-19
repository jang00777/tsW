#include <stdio.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMVA/Reader.h"

#include "classes/DelphesClasses.h"
#include "vtsAnalyserDelphes.h"

using namespace std;

int main(int argc, char* argv[])
{
  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};
  m_decayChannel = std::string{argv[3]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  // read input file
  cout << "Input File : " << inf.c_str() << "\nOutput File : " << outf.c_str() << "\nDecay Channel : " << m_decayChannel.c_str() << endl;

  // set decay channel : di-leptonic (tt->qWqW->qqlvlv) or semi-leptonic (tt->qWqW->qqqqlv)
  int decay_ch = 0;
  if (m_decayChannel == "semilepton" || m_decayChannel == "semi") decay_ch = 1;
  if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   decay_ch = 2;

  // load inputs
  auto inFile = TFile::Open(inf.c_str(), "READ");
  auto inTree = (TTree*) inFile->Get("Delphes");
  inTree->SetBranchStatus("*", true);
  inTree->SetBranchAddress("Event",       &events);
  inTree->SetBranchAddress("GenJet",      &gen_jets);
  inTree->SetBranchAddress("Particle",    &particles);
  inTree->SetBranchAddress("Electron",    &electrons);
  inTree->SetBranchAddress("Muon",        &muons);
  inTree->SetBranchAddress("Jet",         &jets);
  inTree->SetBranchAddress("MissingET",   &missingET);
  inTree->SetBranchAddress("Track",       &tracks);

  // make out output file
  auto out   = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("event", "event");
  auto jettr = new TTree("MVA_jet", "MVA_jet");
  auto hadtr = new TTree("MVA_had", "MVA_had");

  DefBranch(outtr);
  SetJetBranch(jettr);
  SetHadBranch(hadtr);

  TString cutflow_title   = "cutflow" + inf;
  TH1F * cutflow          = new TH1F("cutflow", cutflow_title, 7,-1,6); // -1 : all events, 0 : events after lepton selection, 1~ : events after step

  SetMVAReader();
  //Event Loop Start!
  for (size_t iev = 0; iev < inTree->GetEntries(); ++iev){
    if (iev%1000 == 0 ) cout << "event check    iev    ----> " << iev << endl;
    inTree->GetEntry(iev);
    ResetBranch();
    EventSelection(cutflow, decay_ch);
    if (b_step < 4) continue; // Run events passing jet selection(step4) for di-leptonic case and b-jet selection for semi-leptonic case (see EventSelection() function)
    //cout << " ============================================ " << endl;
    MatchingGenJet();
    HadronReconstruction();
    HadronPreselection();//(m_recoHad);
    FindTruthMatchedHadron();
    FillJetTree(jets, jettr, "pT");
    FillHadTree(hadtr);
    //cout << " ============================================ " << endl;
    outtr->Fill();
  }
  cout << "total no. top : " << tot << endl;
  cout << "matched to tq1 : " << mat1 << endl;
  cout << "matched to tq2 : " << mat2 << endl;

  cout << "s + highest : " << sh << endl;
  cout << "b + highest : " << bh << endl;
  cout << "s + closest : " << sc << endl;
  cout << "b + closest : " << bc << endl;

  cout << "Ratio of (tq1+tq2)/tot : " << 100.*(mat1+mat2)*(1./tot) << endl;

  cout << " no. gen Had(Ks + Lambda) from s : " << setw(3) << testS << endl;
  cout << " no. gen Had(Ks + Lambda) from b : " << setw(3) << testB << endl;

  inFile->Close();
  cutflow->Write();
  h_genWeight->Write();

  out->Write();
  out->Close();
  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

std::map<int, Jet*> JetSelection(TClonesArray* jets, std::vector<struct Lepton> recoLep) {
  std::map<int, Jet*> selectedJets;
  for ( unsigned k = 0; k < jets->GetEntries(); ++k){
    auto jet = (Jet*) jets->At(k);
    if (jet->PT       < cut_JetPt) continue;
    if (abs(jet->Eta) > cut_JetEta) continue;
    bool hasOverlap = false;
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
    for (auto lep : recoLep) { if (jet_tlv.DeltaR(lep.tlv) < cut_JetConeSizeOverlap) hasOverlap = true; }
    if (hasOverlap) continue;
    selectedJets[k] = jet;
  }
  return selectedJets;
}

std::map<int, Jet*> BJetSelection(std::map<int, Jet*> selJets) {
  std::map<int, Jet*> selectedBJets;
  for (auto jet = selJets.begin(); jet != selJets.end(); ++jet){
    if (jet->second->BTag) selectedBJets[jet->first] = jet->second;
  }
  return selectedBJets;
}

void EventSelection(TH1F * cutflow, UInt_t decaychannel){
    for (unsigned i = 0; i < events->GetEntries(); ++i) {
        auto evt = (HepMCEvent*) events->At(i);
        b_genWeight = evt->Weight > 0 ? 1 : -1; // Correct manually
    }
    if (h_genWeight) h_genWeight->Fill(0.5, b_genWeight);

    cutflow->Fill(0);
    // object selection
    std::vector<struct Lepton> recoLep;
    std::vector<struct Lepton> recoEl;
    std::vector<struct Lepton> recoMu;
    SetCutValues(decaychannel); // Cut values are different by channels, decaychannel == 1 : Semileptonic decay channel / decaychannel == 2 : Dileptonic decay channel
    // get recoLep
    for (unsigned i = 0; i < muons->GetEntries(); ++i){
      auto mu = (Muon*) muons->At(i);
      if (abs(mu->Eta) > cut_MuonEta) continue;
      if (mu->PT       < cut_MuonPt)  continue;
      if (mu->IsolationVar > cut_MuonRelIso04All) continue;
      recoMu.push_back(toLepton(mu));
    }
    for (unsigned j = 0; j < electrons->GetEntries(); ++j){
      auto elec = (Electron*) electrons->At(j);
      if (abs(elec->Eta) > cut_ElectronEta) continue;
      if (elec->PT       < cut_ElectronPt)  continue;
      if (elec->IsolationVar > cut_ElectronRelIso03All) continue;
      recoEl.push_back(toLepton(elec));
    }
    // get MET
    auto met = ((MissingET*) missingET->At(0))->MET;
    auto met_eta = ((MissingET*) missingET->At(0))->Eta; // Just added, don't use
    auto met_phi = ((MissingET*) missingET->At(0))->Phi;
    b_MET     = met;
    b_MET_eta = met_eta;
    b_MET_phi = met_phi;

    // Semi-lepton case
    if (decaychannel == 1) {
      // get vetoLep
      std::vector<struct Lepton> vetoEl;
      std::vector<struct Lepton> vetoMu;
      for (unsigned j = 0; j < electrons->GetEntries(); ++j){
        auto elec = (Electron*) electrons->At(j);
        if (elec->PT < cut_VetoElectronPt) continue;
        if (abs(elec->Eta) > cut_VetoElectronEta) continue;
        if (elec->IsolationVar > cut_ElectronRelIso03All ) continue;
        vetoEl.push_back(toLepton(elec));
      }
      for (unsigned i = 0; i < muons->GetEntries(); ++i){
        auto mu = (Muon*) muons->At(i);
        if (mu->PT < cut_VetoMuonPt) continue;
        if (abs(mu->Eta) > cut_VetoMuonEta) continue;
        if (mu->IsolationVar > cut_MuonRelIso04All ) continue;
        vetoMu.push_back(toLepton(mu));
      }

      // step 1 : lepton reconstruction
      if (recoEl.size() + recoMu.size() != 1 ) {b_step = 0; return;}
      cutflow->Fill(1);
      if (recoEl.size() == 1 && recoMu.size() == 0) {b_recoLep1 = recoEl[0].tlv; b_recoLep1_pdgId = recoEl[0].pdgid; recoLep = recoEl;}
      if (recoEl.size() == 0 && recoMu.size() == 1) {b_recoLep1 = recoMu[0].tlv; b_recoLep1_pdgId = recoMu[0].pdgid; recoLep = recoMu;}

      // step 2 : veto additional lepton
      if ( (recoMu.size() == 0 && vetoMu.size() > 0) || (recoMu.size() == 1 && vetoMu.size() >1) ) {b_step = 1; return;}
      if ( (recoEl.size() == 0 && vetoEl.size() > 0) || (recoEl.size() == 1 && vetoEl.size() >1) ) {b_step = 1; return;}
      cutflow->Fill(2);

      // step 3 : jet selection 
      m_selectedJet  = JetSelection(jets, recoLep);
      m_selectedBJet = BJetSelection(m_selectedJet);
      b_nJet         = m_selectedJet.size();
      b_nBJet        = m_selectedBJet.size();
      if (m_selectedJet.size() < cut_nJet) {b_step = 2; return;}
      cutflow->Fill(3);

      // step 4 : b-jet selection
      if (m_selectedBJet.size() < cut_nBJet) {b_step = 3; return;}
      cutflow->Fill(4);
      b_step = 4;
    }

    // Di-lepton case
    if (decaychannel == 2) {
      // step 1 : lepton pair reconstruction
      if (recoEl.size() + recoMu.size() != 2 ) {b_step = 0; return;}
  
      if (recoEl.size() == 1 && recoMu.size() == 1) {
        recoLep.push_back(recoEl[0]);
        recoLep.push_back(recoMu[0]);
      } else if (recoEl.size() == 0 && recoMu.size() == 2) {
        recoLep.push_back(recoMu[0]);
        recoLep.push_back(recoMu[1]);
      } else if (recoEl.size() == 2 && recoMu.size() == 0) {
        recoLep.push_back(recoEl[0]);
        recoLep.push_back(recoEl[1]);
      }
      b_recoLep1       = recoLep[0].tlv;   b_recoLep2       = recoLep[1].tlv;
      b_recoLep1_pdgId = recoLep[0].pdgid; b_recoLep2_pdgId = recoLep[1].pdgid;

      b_recoLep1_MET_dphi = TVector2::Phi_mpi_pi(b_recoLep1.Phi() - b_MET_phi);
      b_recoLep2_MET_dphi = TVector2::Phi_mpi_pi(b_recoLep2.Phi() - b_MET_phi);

      auto dilepton = recoLep[0].tlv + recoLep[1].tlv;
      b_dilepton    = dilepton;
      b_dilepton_MET_dphi = TVector2::Phi_mpi_pi(b_dilepton.Phi() - b_MET_phi);

      auto mulpdg   = recoLep[0].charge * recoLep[1].charge;
      b_dilepton_ch   = abs(recoLep[0].pdgid) + abs(recoLep[1].pdgid); // 22 -> ee , 24 -> emu , 26 -> mumu
      b_dilepton_mass = dilepton.M();

      if (b_dilepton_mass < cut_Dilep_M || mulpdg > 0) {b_step = 0; return;}
      cutflow->Fill(1);
  
      // step 2 : veto Z boson
      if ( (b_dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) {b_step = 1; return;}
      cutflow->Fill(2);
  
      // step 3 : MET cut
      if ( (b_dilepton_ch != 24) && (met < cut_MET_Pt) ) {b_step = 2; return;}
      cutflow->Fill(3);
  
      // step 4 : jet selection
      m_selectedJet  = JetSelection(jets, recoLep);
      m_selectedBJet = BJetSelection(m_selectedJet);
      b_nJet         = m_selectedJet.size();
      b_nBJet        = m_selectedBJet.size();
      if (m_selectedJet.size() < cut_nJet) {b_step = 3; return;}
      cutflow->Fill(4);
  
      // step 5 : b-jet selection (not used)
      if (m_selectedBJet.size() < cut_nBJet) {b_step = 4; return;}
      cutflow->Fill(5);
      b_step = 5;
    }
}

void MatchingGenJet() {
  std::vector<const GenParticle*> genTops;
  genTops.clear();
  for (int i = 0; i < particles->GetEntries(); ++ i){ // Loop over genParticles and collect necessary information 
    auto p = (const GenParticle*) particles->At(i);
    if (p->Status   >  30) continue;
    if (abs(p->PID) != 6)  continue;
    auto top = getLast(particles,p);
    genTops.push_back(top);

    auto quark  = (const GenParticle*) particles->At(top->D1); // gen quark from top-quark
    auto Wboson = (const GenParticle*) particles->At(top->D2); // gen W-boson from top-quark

    if (abs(quark->PID) != 3 && abs(quark->PID) != 5) {
      quark     = (const GenParticle*) particles->At(top->D2);
      Wboson    = (const GenParticle*) particles->At(top->D1);
    }

    auto lastBoson = getLast(particles, Wboson); // Get a last copy of W-boson in the decay chain
    m_genQuark.push_back(quark);

    auto dau1_from_W = (const GenParticle*) particles->At(lastBoson->D1); // for dilep, always neutrino
    auto dau2_from_W = (const GenParticle*) particles->At(lastBoson->D2); // for dilep, always e/mu/tau

    struct Lepton lep1 = toLepton(dau1_from_W);
    struct Lepton lep2 = toLepton(dau2_from_W);

    //cout << " Gen.  TOP  (PID | STATUS) : ( " << setw(5) <<  top->PID         << " | " << setw(5) << top->Status         << endl;
    //cout << " Gen.  DAU1 (PID | STATUS) : ( " << setw(5) <<  quark->PID       << " | " << setw(5) << quark->Status       << endl;
    //cout << " Gen.  DAU2 (PID | STATUS) : ( " << setw(5) <<  Wboson->PID      << " | " << setw(5) << Wboson->Status      << endl;
    //cout << " Gen.  DDA1 (PID | STATUS) : ( " << setw(5) <<  dau1_from_W->PID << " | " << setw(5) << dau1_from_W->Status << endl;
    //cout << " Gen.  DDA2 (PID | STATUS) : ( " << setw(5) <<  dau2_from_W->PID << " | " << setw(5) << dau2_from_W->Status << endl;
    //cout << " Reco. DDA1 (PID | STATUS) : ( " << setw(5) <<  lep1.pdgid       << " | " << setw(5) << dau1_from_W->Status << endl;
    //cout << " Reco. DDA2 (PID | STATUS) : ( " << setw(5) <<  lep2.pdgid       << " | " << setw(5) << dau2_from_W->Status << endl;

    // Collect leptons from W-boson, for dilepton, the size of m_genLepton should be 2 and for semilepton, 1
    if (abs(lep1.pdgid) == 11 || abs(lep1.pdgid) == 13) m_genLepton.push_back(lep1);
    if (abs(lep2.pdgid) == 11 || abs(lep2.pdgid) == 13) m_genLepton.push_back(lep2);
  }

  if ( m_genQuark.size() < 2 ) {
    cout << "No two quarks from top quarks ... skip this event" << endl; 
    return;
  } else if ( (m_genLepton.size() != 2 && m_decayChannel == "di") || (m_genLepton.size() != 1 && m_decayChannel == "semi")) {
    //cout << "No two reco leptons in dileptonic channel or no reco lepton in semileptonic channel ===> Maybe because a case that W boson decay into tau / neutrino isn't included" << endl;
    return;
  }

  tot += m_genQuark.size();

  //std::cout << "TQ1 : " << m_genQuark[0]->PID << " TQ2 : " << m_genQuark[1]->PID << endl;

  TLorentzVector wl1_tlv;
  TLorentzVector wl2_tlv;
  if (m_decayChannel == "dilepton" || m_decayChannel == "di") {
    if (abs(m_genLepton[0].pdgid) == abs(m_genLepton[1].pdgid)) {
      if      (abs(m_genLepton[0].pdgid) == 11) b_channel = 1; // ee channel
      else if (abs(m_genLepton[0].pdgid) == 13) b_channel = 2; // mm channel
      else if (abs(m_genLepton[0].pdgid) == 15) b_channel = 3; // tau
    } else b_channel = 0;                                      // em channel
  } else if (m_decayChannel == "semilepton" || m_decayChannel == "semi") {
    if      (abs(m_genLepton[0].pdgid) == 11) b_channel = 1; // e+jet channel
    else if (abs(m_genLepton[0].pdgid) == 13) b_channel = 2; // mu+jet channel
    else if (abs(m_genLepton[0].pdgid) == 15) b_channel = 3; // tau+jet channel
  }

  int nMatchToQ1 = 0;  int nMatchToQ2 = 0;
  b_multimatchedS = 0; b_unmatchedS = 0;   b_jetOverlap = 0; 
  b_multimatchedB = 0; b_unmatchedB = 0;
  for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){ // Loop over jets that passed jet selection
    auto jIdx = jet->first; 
    auto selJet = jet->second;
    TLorentzVector jet_tlv, tq1_tlv, tq2_tlv;
    jet_tlv.SetPtEtaPhiM(selJet->PT, selJet->Eta, selJet->Phi, selJet->Mass);
    tq1_tlv.SetPtEtaPhiM(m_genQuark[0]->PT, m_genQuark[0]->Eta, m_genQuark[0]->Phi, m_genQuark[0]->Mass);
    tq2_tlv.SetPtEtaPhiM(m_genQuark[1]->PT, m_genQuark[1]->Eta, m_genQuark[1]->Phi, m_genQuark[1]->Mass);

    if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   { wl1_tlv = m_genLepton[0].tlv; wl2_tlv = m_genLepton[1].tlv;}
    if (m_decayChannel == "semilepton" || m_decayChannel == "semi") { wl1_tlv = m_genLepton[0].tlv;}

    //cout << "dR(jet, tq1) : " << setw(12) <<  jet_tlv.DeltaR(tq1_tlv) << " dR cut : " << cut_JetConeSizeOverlap << " jidx : " << jIdx << endl;
    //cout << "dR(jet, tq2) : " << setw(12) <<  jet_tlv.DeltaR(tq2_tlv) << " dR cut : " << cut_JetConeSizeOverlap << " jidx : " << jIdx << endl;
    //cout << "dR(jet, lep1): " << setw(12) <<  jet_tlv.DeltaR(wl1_tlv) << " dR cut : " << cut_JetConeSizeOverlap << " jidx : " << jIdx << endl;
    //cout << "dR(jet, lep2): " << setw(12) <<  jet_tlv.DeltaR(wl2_tlv) << " dR cut : " << cut_JetConeSizeOverlap << " jidx : " << jIdx << endl;

    // deltaR matching between reco. jet and quark from top-quark (for m_matchedJet, see struct RecoJet in header file)
    if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap && jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap ) { // Overlap case ==> we will not use this jet
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), -99,                          true,  -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), -99,                     -99,                          true,  -99, false, false}); 
      b_jetOverlap += 1;
    } else if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap && jet_tlv.DeltaR(tq2_tlv) > cut_JetConeSizeOverlap) { // jet matched to quark1 from top-quark
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[0]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[0]->PT, false, -99, false, false}); 
      mat1 += 1;
      nMatchToQ1 += 1;
    } else if (jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap && jet_tlv.DeltaR(tq1_tlv) > cut_JetConeSizeOverlap) { // jet matched to quark2 from to-quark
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
      mat2 += 1;
      nMatchToQ2 += 1;
    }
  }

  // Check ambiguity per event (mutlmatched/unmatched case will be shown from isSJet or isBJet and jetOverlap will be shown from isOverlap for the case of jet by jet)
  if (nMatchToQ1 != 0) {
    if (abs(m_genQuark[0]->PID) == 3) b_multimatchedS = nMatchToQ1;
    if (abs(m_genQuark[0]->PID) == 5) b_multimatchedB = nMatchToQ1;
  } else if (nMatchToQ1 == 0) {
    if (abs(m_genQuark[0]->PID) == 3) b_unmatchedS    = 1;
    if (abs(m_genQuark[0]->PID) == 5) b_unmatchedB    = 1;   
  }
  if (nMatchToQ2 != 0) {
    if (abs(m_genQuark[1]->PID) == 3) b_multimatchedS = nMatchToQ2;
    if (abs(m_genQuark[1]->PID) == 5) b_multimatchedB = nMatchToQ2;
  } else if (nMatchToQ2 == 0) {
    if (abs(m_genQuark[1]->PID) == 3) b_unmatchedS    = 1;
    if (abs(m_genQuark[1]->PID) == 5) b_unmatchedB    = 1;
  }

  if (m_matchedJet.size() != 0) {
    // Sort by pT and dR(lep, jet) and then set a flag for jets with highest top 2 pT or closest distance to lepton
    sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.j->PT > b.j->PT); } ); // Sort with pT ordering
    if (m_matchedJet.size() > 0) m_matchedJet[0].hasHighestPt  = true;
    if (m_matchedJet.size() > 1) m_matchedJet[1].hasHighestPt  = true;
    //for (auto j = 0; j < m_matchedJet.size(); ++j) cout << j << "th matched jet pT : " << m_matchedJet[j].j->PT << " / hasHighestPt : " << m_matchedJet[j].hasHighestPt << " selectedJet size : " << m_selectedJet.size() << " matchedJet size : " << m_matchedJet.size() << endl; 
    sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.drl1 < b.drl1); } ); // Sort with dR(lep, jet) ordering
    m_matchedJet[0].hasClosestLep = true;
    if (m_decayChannel == "dilepton" || m_decayChannel == "di") sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.drl2 < b.drl2); } ); // Sort with dR(lep, jet) ordering
    m_matchedJet[1].hasClosestLep = true;
  }
}

void HadronReconstruction() {
  //cout << " >>>>>>>>>>>>>>>>>>>>>>> HADRON RECO <<<<<<<<<<<<<<<<<<<<<<<< " << endl;
  // Pick two charged hadron pairs
  for (auto i = 0; i < tracks->GetEntries(); ++i) {
    auto hadCand1 = (Track*) tracks->At(i);
    if (hadCand1->Charge != 1) continue; // Only pick positive charge
    if (abs(hadCand1->PID) == 11 || abs(hadCand1->PID) == 13) continue; // Exclude leptons
    if (hadCand1->PT < tkPTCut_) continue;
    if (fabs(hadCand1->Eta) > tkEtaCut_) continue;
    if (fabs(hadCand1->D0/hadCand1->ErrorD0) < tkIPSigXYCut_) continue; // Cut of siginifcance of transverse impact parametr
    for (auto j = 0; j < tracks->GetEntries(); ++j) {
      auto hadCand2 = (Track*) tracks->At(j);
      if (hadCand2->Charge != -1) continue; // Only pick negative charge
      if (abs(hadCand2->PID) == 11 || abs(hadCand2->PID) == 13) continue; // Exclude leptons
      if (hadCand2->PT < tkPTCut_) continue;
      if (fabs(hadCand2->Eta) > tkEtaCut_) continue;
      if (fabs(hadCand2->D0/hadCand2->ErrorD0) < tkIPSigXYCut_) continue; // Cut of siginifcance of transverse impact parametr

      TLorentzVector dauCand1_pion_tlv;   
      TLorentzVector dauCand1_proton_tlv; 
      TLorentzVector dauCand2_pion_tlv;   
      TLorentzVector dauCand2_proton_tlv; 

      //dauCand1_pion_tlv.SetPtEtaPhiM(  hadCand1->PT, hadCand1->Eta, hadCand1->Phi, pionMass_);
      //dauCand1_proton_tlv.SetPtEtaPhiM(hadCand1->PT, hadCand1->Eta, hadCand1->Phi, protonMass_);
      //dauCand2_pion_tlv.SetPtEtaPhiM(  hadCand2->PT, hadCand2->Eta, hadCand2->Phi, pionMass_);
      //dauCand2_proton_tlv.SetPtEtaPhiM(hadCand2->PT, hadCand2->Eta, hadCand2->Phi, protonMass_);

      // Get not perigee momentum (see ParticlePropagator.cc) but original momentum (gen level pion) 
      auto genDau1     = (GenParticle*) hadCand1->Particle.GetObject();
      auto genDau2     = (GenParticle*) hadCand2->Particle.GetObject();
      dauCand1_pion_tlv.SetPtEtaPhiM(  hadCand1->PT, hadCand1->Eta, genDau1->Phi, pionMass_);
      dauCand1_proton_tlv.SetPtEtaPhiM(hadCand1->PT, hadCand1->Eta, genDau1->Phi, protonMass_);
      dauCand2_pion_tlv.SetPtEtaPhiM(  hadCand2->PT, hadCand2->Eta, genDau2->Phi, pionMass_);
      dauCand2_proton_tlv.SetPtEtaPhiM(hadCand2->PT, hadCand2->Eta, genDau2->Phi, protonMass_);

      TLorentzVector hadCand_pion_pion_tlv    = dauCand1_pion_tlv   + dauCand2_pion_tlv;
      TLorentzVector hadCand_pion_proton_tlv  = dauCand1_pion_tlv   + dauCand2_proton_tlv;
      TLorentzVector hadCand_proton_pion_tlv  = dauCand1_proton_tlv + dauCand2_pion_tlv;
      
      if ((fabs(hadCand_pion_pion_tlv.M() - ksMass_)/ksMass_) < 0.3) { // Invariant mass cut 
        m_recoHad.push_back({hadCand_pion_pion_tlv, dauCand1_pion_tlv, dauCand2_pion_tlv, -1, i, j, ksPdgId_, pionPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_pion_proton_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) { // Invariant mass cut
        m_recoHad.push_back({hadCand_pion_proton_tlv, dauCand1_pion_tlv, dauCand2_proton_tlv, -1, i, j, lambda0PdgId_, pionPdgId_, (-1)*protonPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_proton_pion_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) { // Invariant mass cut
        m_recoHad.push_back({hadCand_proton_pion_tlv, dauCand1_proton_tlv, dauCand2_pion_tlv, -1, i, j, lambda0PdgId_, protonPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
    }
  }
}

void HadronPreselection(){//(std::vector<RecoHad> recoHad) {
  // Preselection cuts for reconstructed hadrons
  // Not fully implemented yet
  for (auto i=0; i < m_recoHad.size(); ++i) {
    m_recoHad[i].isPreSelected = false;
    if (fabs(m_recoHad[i].tlv.Eta()) > 2.5) continue;
    if (m_recoHad[i].dau1_tlv.Pt() < 0.95 || m_recoHad[i].dau2_tlv.Pt() < 0.95) continue;
    //auto hadDau1 = (Track*) tracks->At(m_recoHad[i].dau1_idx);
    //auto hadDau2 = (Track*) tracks->At(m_recoHad[i].dau2_idx);
    //if (fabs(hadDau1->D0/hadDau1->ErrorD0) < 5. || fabs(hadDau2->D0/hadDau2->ErrorD0) < 5.) continue; // <- D0Sig cut( ==2sig) is applied already in HadronReconstruction !
    else m_recoHad[i].isPreSelected = true;
  }
}

void FindTruthMatchedHadron() {
  int n = 0;
  for (int i = 0; i < particles->GetEntries(); ++i){ // Loop over all genParticles
    bool  isFromTop = false;
    bool  isFromW   = false;
    int   isFrom    = -99;
    float x         = -99;
    auto p = (const GenParticle*) particles->At(i);
    if (abs(p->PID) != 310 && abs(p->PID) != 3122) continue; // We are only intrested in some hadrons (Ks(pdgId == 310) and Lambda0(pdgId == 3122))
    if (p->D1  == -1  || p->D2  == -1)   continue; // We need a gen. level link to their daughters
    // Get gen. daughter objects from gen hadrons
    auto d1 = (const GenParticle*) particles->At(p->D1);
    auto d2 = (const GenParticle*) particles->At(p->D2);
    if ( abs(p->PID == 310)  && (abs(d1->PID) != 211 || abs(d2->PID) != 211) ) continue; // Pick only charged pion pair, which is most dominant decay process (See http://pdg.lbl.gov/2019/listings/rpp2019-list-K-zero-S.pdf)
    if ( abs(p->PID == 3122) ) { // Pick only proton - pion pair, which is most dominant decay process (See http://pdg.lbl.gov/2019/listings/rpp2019-list-lambda.pdf) 
      if      ( (abs(d1->PID) != 211  && abs(d1->PID) != 2212) || (abs(d2->PID) != 211  && abs(d2->PID) != 2212) ) continue;
      else if ( abs(d1->PID) == abs(d2->PID) ) continue;
    }
    if (d1->PID*d2->PID > 0) continue; // Check opposite charge 
    //auto motherList = getMomList(particles, p); // Get All mothers of gen. hadron
    auto motherList = getMomIdxList(particles, p); // Get All mothers' idx of gen. hadron
    for (auto j=0; j < motherList.size()-1; ++j) {
      auto mom = (const GenParticle*)particles->At(motherList[j]);
      if ( abs(mom->PID) == 24 && !isFromW) { // Check whether gen. hadron comes from W-boson or not
        auto prev_mom = (const GenParticle*)particles->At(motherList[j-1]);
        isFromW = true;
        isFrom  = prev_mom->PID; // Check flavour of quark from which gen. hadron comes
        x       = p->PT/prev_mom->PT;
      }
      if ( abs(mom->PID) == 6 && mom->Status == 62 ) { // Check if gen. hadron comes from top-quark
        auto prev_mom = (const GenParticle*)particles->At(motherList[j-1]);
        isFromTop = true;
        if (isFromW == false) { isFrom = prev_mom->PID; x = p->PT/prev_mom->PT; }
        break;
      }
    }
    // Can't trust P4() function ...
    TLorentzVector p_tlv, d1_tlv, d2_tlv;
    if (isFromTop) {
    //  cout << "m_genHadron_CHK : " << "genHad pdgId      : " << setw(5) << p->PID   << " " 
    //                               << "genHad mother     : " << setw(5) << isFrom   << " "
    //                               << "genHad dau1 pdgId : " << setw(5) << d1->PID  << " "
    //                               << "genHad dau2 pdgId : " << setw(5) << d2->PID  << " "
    //                               << "genHad isFromTop  : " << setw(5) << isFromTop<< " " 
    //                               << "genHad isFromW    : " << setw(5) << isFromW  << " " 
    //       << endl;
      if (abs(isFrom) == 3) { testS += 1; }
      if (abs(isFrom) == 5) { testB += 1; }

    }
    p_tlv.SetPtEtaPhiM(p->PT, p->Eta, p->Phi, p->Mass);
    d1_tlv.SetPtEtaPhiM(d1->PT, d1->Eta, d1->Phi, d1->Mass);
    d2_tlv.SetPtEtaPhiM(d2->PT, d2->Eta, d2->Phi, d2->Mass);
    m_genHadron[i] = {p_tlv, d1_tlv, d2_tlv, x, i, p->PID, p->D1, d1->PID, p->D2, d2->PID, isFrom, isFromTop, isFromW};
  }
}

int FindMatchedHadron(TLorentzVector jet_tlv, TString method, TTree* jettr) {
  int   hIdx = -1;
  int   prevPdgId   = -99;
  float maxDeltaScore  = -99999.;
  float maxBDTScore = -1;
  float highestPt   = -1;
  float lowestDot   = 99999999.;
  float closestDr   = -1;
  for (auto i=0; i<m_recoHad.size(); ++i) {
    auto had_tlv = m_recoHad[i].tlv;
    auto dR      = had_tlv.DeltaR(jet_tlv);
    if (abs(m_recoHad[i].pdgid) != 310 && abs(m_recoHad[i].pdgid) != 3122) continue;
    if (m_hadPreCut && !m_recoHad[i].isPreSelected) continue;
    if (dR > cut_JetConeSizeOverlap) continue;

    if (method == "BDT") { // The method isn't constructed yet

      SetHadValues(jettr, i, jet_tlv); // Load hadron values temporarily for scoring

      float had_bdtScore = -99.;
      float deltaScore   = -99999.;
      if (abs(m_recoHad[i].pdgid) == 310)  { had_bdtScore = m_hadReader->EvaluateMVA("KS_BDT");     deltaScore = (had_bdtScore - KsBDTCut)/fabs(KsBDTCut); }
      if (abs(m_recoHad[i].pdgid) == 3122) { had_bdtScore = m_hadReader->EvaluateMVA("Lambda_BDT"); deltaScore  = (had_bdtScore - LambdaBDTCut)/fabs(LambdaBDTCut); }
      if (maxBDTScore < had_bdtScore) {
        if (abs(m_recoHad[i].pdgid) != abs(prevPdgId)) { // If pdgId is different, check how far the score is from the optimal cut of the particle type and choose more greater one that means more real-like
          if (maxDeltaScore < deltaScore) {
            prevPdgId     = m_recoHad[i].pdgid;
            maxDeltaScore = deltaScore;
            maxBDTScore   = had_bdtScore;
            hIdx = i;
          }
        } else {
          maxBDTScore = had_bdtScore;
          hIdx = i;
        }
      }
      //if (abs(m_recoHad[i].pdgid) == 310)  cout << setw(5) << i << " th had_bdtScore : " << setw(10) << had_bdtScore << " maxBDTScore : " << setw(10) << maxBDTScore << " pdgId : " << setw(6) << m_recoHad[i].pdgid << " prevPdgId : " << setw(6) << prevPdgId << " delta Optcut : " << setw(10) << (had_bdtScore - KsBDTCut)/fabs(KsBDTCut)         << " maxDelta : " << setw(10) << maxDeltaScore << endl;
      //if (abs(m_recoHad[i].pdgid) == 3122) cout << setw(5) << i << " th had_bdtScore : " << setw(10) << had_bdtScore << " maxBDTScore : " << setw(10) << maxBDTScore << " pdgId : " << setw(6) << m_recoHad[i].pdgid << " prevPdgId : " << setw(6) << prevPdgId << " delta Optcut : " << setw(10) << (had_bdtScore - LambdaBDTCut)/fabs(LambdaBDTCut) << " maxDelta : " << setw(10) << maxDeltaScore << endl;
    }
    if (method == "pT" || method == "pt") { // pT matching
      auto hadronPt = had_tlv.Pt();
      if (hadronPt > highestPt) {
        highestPt = hadronPt;
        hIdx = i;
      }
    }
    if (method == "dot") { // 4-vector matching
      auto dotProduct = had_tlv.Dot(jet_tlv);
      if (dotProduct < lowestDot) {
        lowestDot = dotProduct;
        hIdx = i;
      }
    }
    if (method == "dr" || method == "dR") { // deltaR matching
      if (dR < closestDr) {
        closestDr = dR;
        hIdx = i;
      }
    }

  }
  if (method == "BDT") ResetHadValues(); // clean hadron values loaded temporarily for evaluating a BDT score

  if (hIdx != -1) m_recoHad[hIdx].isJetMatched = true;
  return hIdx;
} 

void FillJetTree(TClonesArray* jets, TTree* jettr, TString matchingMethod) {
  int nHigh = 0;
  int nHighSel = 0;
  b_jet_start = jettr->GetEntries();

  std::vector<std::pair<int, float>> highestPtJet;
  for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){ // Loop over jets that passed jet selection
    auto jidx   = jet->first;
    auto selJet = jet->second;
    highestPtJet.push_back({jidx, selJet->PT});
    //cout << i << " th jet pT ===> : " << jet->PT << endl;
  }
  sort(highestPtJet.begin(), highestPtJet.end(), [](std::pair<int, float> a, std::pair<int, float> b) {return a.second > b.second;});
  //cout << "-----------------------------------------------------------------------------------------------" << endl;
  //cout << " Highest pT jet is : " << highestPtJet[0].first << " / " << highestPtJet[1].first << endl;

  for ( unsigned i = 0; i < jets->GetEntries(); ++i){
    ResetJetValues();
    auto jet     = (Jet*) jets->At(i);
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);

    if (m_selectedJet[i] != 0) { // check if the i-th jet is selected jet from event selection 
      b_isSelectedJet = true;
    } else continue;
    for ( auto j = 0; j < m_matchedJet.size(); ++j ) {
      if (m_matchedJet[j].idx == i) { // check if the i-th jet is jet which matched to gen s/b quark from top-quark
        b_isFrom        = m_matchedJet[j].gen_pdgid;
        b_hasHighestPt  = m_matchedJet[j].hasHighestPt;
        b_hasClosestLep = m_matchedJet[j].hasClosestLep;
        b_isOverlap     = m_matchedJet[j].isOverlap;
        //if (abs(b_isFrom) == 3 && b_hasHighestPt) cout << "CHK : " << j << " / " << m_matchedJet[j].idx <<endl;
      }
    }
    if (( i == highestPtJet[0].first || i == highestPtJet[1].first)) b_hasHighestPtFromSelected = true; // this is for checking jets with highest Pt but not matched to gen s/b quarks from top

    //cout << setw(2) << i << " / " << setw(2) << jets->GetEntries() 
    //     << " >>>> " 
    //     << " | FLAV : "              << setw(3)  << jet->Flavor 
    //     << " | isFrom : "            << setw(3)  << b_isFrom 
    //     << " | idx : "               << setw(2)  << i 
    //     << " | PT : "                << setw(10) << jet->PT 
    //     << " | hasHighest : "        << setw(2)  << b_hasHighestPt 
    //     << " | hasHighestFromSel : " << setw(2)  << b_hasHighestPtFromSelected 
    //     << " | bTag : "              << setw(2)  << jet->BTag 
    //     << " | dR(lep1, J) : "       << setw(10) << b_recoLep1.DeltaR(jet_tlv)
    //     << " | dR(lep2, J) : "       << setw(10) << b_recoLep2.DeltaR(jet_tlv)
    //     << endl;

    b_recoLep1_dr = b_recoLep1.DeltaR(jet_tlv);
    b_recoLep2_dr = b_recoLep2.DeltaR(jet_tlv);

    if (abs(b_isFrom) == 3 && b_hasHighestPt) sh +=1;
    if (abs(b_isFrom) == 5 && b_hasHighestPt) bh +=1;
    if (abs(b_isFrom) == 3 && b_hasClosestLep) sc +=1;
    if (abs(b_isFrom) == 5 && b_hasClosestLep) bc +=1;

    if (b_hasHighestPt == true) nHigh += 1;
    if (b_hasHighestPtFromSelected == true) nHighSel += 1;
    auto hIdx = FindMatchedHadron(jet_tlv, matchingMethod, jettr); // check if reco. hadron(Ks or lambda0) is inside the i-th jet according to a matching method, if you use "BDT" method, tree arg. should be imported into the function
    SetJetValues(jet);
    if (hIdx != -1) {
      SetHadValues(jettr, hIdx, jet_tlv); // if there is a reco. hadron within the i-th jet, then get information of reco. hadron
      if (abs(b_had_pdgId) == 310) {
        b_had_bdt_score       = m_hadReader->EvaluateMVA("KS_BDT");
        b_had_isSelectedByMVA = (b_had_bdt_score > KsBDTCut);
      } else if (abs(b_had_pdgId) == 3122) {
        b_had_bdt_score       = m_hadReader->EvaluateMVA("Lambda_BDT");
        b_had_isSelectedByMVA = (b_had_bdt_score > LambdaBDTCut);
      }
    }
    jettr->Fill();
  }
  b_jet_end = jettr->GetEntries();
  if (nHigh > 2)    cout << "no. hasHighestPt in the Event : "             << nHigh << endl;
  if (nHighSel > 2) cout << "no. hasHighestPtFromSelected in the Event : " << nHigh << endl;
  highestPtJet.clear();
}

void FillHadTree(TTree* hadtr) {
  b_had_start = hadtr->GetEntries();
  for ( unsigned i = 0; i < m_recoHad.size(); ++i){
    ResetHadValues();
    SetHadValues(hadtr, i);
    if (abs(b_had_pdgId) == 310) {
      b_had_bdt_score       = m_hadReader->EvaluateMVA("KS_BDT");
      b_had_isSelectedByMVA = (b_had_bdt_score > KsBDTCut);
    } else if (abs(b_had_pdgId) == 3122) {
      b_had_bdt_score       = m_hadReader->EvaluateMVA("Lambda_BDT");
      b_had_isSelectedByMVA = (b_had_bdt_score > LambdaBDTCut);
    }
    hadtr->Fill();
  }
  b_had_end = hadtr->GetEntries();
}

void SetJetValues(Jet* jet) {
  b_pt               = jet->PT;
  b_eta              = jet->Eta;
  b_phi              = jet->Phi;
  b_mass             = jet->Mass;
  b_radiusInEta      = jet->DeltaEta;
  b_radiusInPhi      = jet->DeltaPhi;
  //b_ptD              = jet->PTD;
  b_cMult            = jet->NCharged;
  b_nMult            = jet->NNeutrals;
  b_pdgId            = jet->Flavor;
  b_flavorAlgo       = jet->FlavorAlgo;
  b_flavorPhys       = jet->FlavorPhys;
  b_bTag             = jet->BTag;
  b_bTagAlgo         = jet->BTagAlgo;
  b_bTagPhys         = jet->BTagPhys;
  b_tauTag           = jet->TauTag;
  b_tauWeight        = jet->TauWeight;

  CollectJetConstituentInfo(jet);

  auto charged_partons = CollectJetConstituentPt(jet, "charged");
  auto neutral_partons = CollectJetConstituentPt(jet, "neutral");
  while (charged_partons.size() < 3) charged_partons.push_back(0.);
  while (neutral_partons.size() < 3) neutral_partons.push_back(0.);
  std::sort(charged_partons.begin(), charged_partons.end(), [] (float a, float b) { return (a > b); } );
  std::sort(neutral_partons.begin(), neutral_partons.end(), [] (float a, float b) { return (a > b); } );

  b_cpt1 = charged_partons[0]; 
  b_cpt2 = charged_partons[0] + charged_partons[1];
  b_cpt3 = charged_partons[0] + charged_partons[1] + charged_partons[2];
  b_cptA = std::accumulate(charged_partons.begin(), charged_partons.end(), 0.0);

  b_npt1 = neutral_partons[0]; 
  b_npt2 = neutral_partons[0] + neutral_partons[1];
  b_npt3 = neutral_partons[0] + neutral_partons[1] + neutral_partons[2];
  b_nptA = std::accumulate(neutral_partons.begin(), neutral_partons.end(), 0.0);

}

std::vector<float> CollectJetConstituentPt(Jet* jet, TString type) {
  std::vector<float> list;
  for (auto j=0; j < jet->Constituents.GetEntriesFast(); ++j) {
    auto object = jet->Constituents.At(j);
    if (object->IsA() == Track::Class() && type.Contains("charged")) {
      auto charged_parton = (Track*) object;
      list.push_back(charged_parton->PT); 
    }
    if (object->IsA() == Tower::Class() && type.Contains("neutral")) {
      auto neutral_parton = (Tower*) object;
      if (neutral_parton->ET < 1.0) continue;
      list.push_back(neutral_parton->ET);
    }
  }
  return list;
}

void CollectJetConstituentInfo(Jet* jet) {
  //cout << " nGenHadron : " << setw(4) << m_genHadron.size() << endl;
  //int n = 0;
  //for (auto it = m_genHadron.begin(); it != m_genHadron.end(); ++it) {
  //  auto i = it->second;
  //  n += 1;
  //  if (i.pdgid == 0) continue;
  //  cout << setw(4) << n << " th genHad : " << setw(6) << i.pdgid  << " dau1 : "      << setw(6) << i.dau1_pdgid << " dau2 : "    << setw(6) << i.dau2_pdgid 
  //                       << " isFrom : "    << setw(6) << i.isFrom << " isFromTop : " << setw(6) << i.isFromTop  << " isFromW : " << setw(6) << i.isFromW << endl;
  //}

  // Parameters for Quark gluon jet variables
  float sum_dpt    = 0, sum_dpt2    = 0, sum_deta   = 0, sum_deta2   = 0, sum_dphi   = 0, sum_dphi2   = 0, sum_detadphi   = 0;  
  float sum_c_dpt  = 0, sum_c_dpt2  = 0, sum_c_deta = 0, sum_c_deta2 = 0, sum_c_dphi = 0, sum_c_dphi2 = 0, sum_c_detadphi = 0; 
  float sum_n_dpt  = 0, sum_n_dpt2  = 0, sum_n_deta = 0, sum_n_deta2 = 0, sum_n_dphi = 0, sum_n_dphi2 = 0, sum_n_detadphi = 0; 
  float ave_deta   = 0, ave_deta2   = 0, ave_dphi   = 0, ave_dphi2   = 0;
  float ave_c_deta = 0, ave_c_deta2 = 0, ave_c_dphi = 0, ave_c_dphi2 = 0;
  float ave_n_deta = 0, ave_n_deta2 = 0, ave_n_dphi = 0, ave_n_dphi2 = 0;

  for (auto j=0; j < jet->Constituents.GetEntriesFast(); ++j) {
    auto object = jet->Constituents.At(j);
    if (object->IsA() == Track::Class()) {
      auto  cp   = (Track*) object;
      float pt2  = cp->PT * cp->PT;
      float deta = (cp->Eta - jet->Eta);
      float dphi = DeltaPhi(cp->Phi, jet->Phi);
      sum_dpt        += cp->PT;
      sum_dpt2       += pt2;
      sum_deta       += deta*pt2;
      sum_deta2      += deta*deta*pt2;
      sum_dphi       += dphi*pt2;
      sum_dphi2      += dphi*dphi*pt2;
      sum_detadphi   += deta*dphi*pt2;

      sum_c_dpt      += cp->PT;
      sum_c_dpt2     += pt2;
      sum_c_deta     += deta*pt2;
      sum_c_deta2    += deta*deta*pt2;
      sum_c_dphi     += dphi*pt2;
      sum_c_dphi2    += dphi*dphi*pt2;
      sum_c_detadphi += deta*dphi*pt2;

      b_constituent_pt.push_back(cp->PT); 
      b_constituent_eta.push_back(cp->Eta); 
      b_constituent_phi.push_back(cp->Phi); 
      b_constituent_charge.push_back(cp->Charge);
      b_constituent_D0.push_back(cp->D0);
      b_constituent_D0Err.push_back(cp->ErrorD0);
      b_constituent_DZ.push_back(cp->DZ);
      b_constituent_DZErr.push_back(cp->ErrorDZ);
      b_constituent_CtgTheta.push_back(cp->CtgTheta);
      if (abs(cp->PID) == 11 || abs(cp->PID) == 13) {
        b_constituent_type.push_back(abs(cp->PID));
      } else {
        b_constituent_type.push_back(1);
      }
      auto genp = (GenParticle*) cp->Particle.GetObject();
      if (abs(m_genHadron[genp->M1].pdgid) == 310) { 
        b_constituent_isFromKs.push_back(true);
        b_constituent_isKsFrom.push_back(m_genHadron[genp->M1].isFrom);
        if (m_genHadron[genp->M1].isFromTop)                                    b_constituent_isKsFromTop.push_back(true);
        else                                                                    b_constituent_isKsFromTop.push_back(false);
        if (m_genHadron[genp->M1].isFromW)                                      b_constituent_isKsFromW.push_back(true);
        else                                                                    b_constituent_isKsFromW.push_back(false);
        if (!m_genHadron[genp->M1].isFromTop && !m_genHadron[genp->M1].isFromW) b_constituent_isKsFromOther.push_back(true);
        else                                                                    b_constituent_isKsFromOther.push_back(false);
      } else { 
        b_constituent_isFromKs.push_back(false);
        b_constituent_isKsFrom.push_back(-99);
        b_constituent_isKsFromTop.push_back(false);
        b_constituent_isKsFromW.push_back(false);
        b_constituent_isKsFromOther.push_back(false);
      }
      if (abs(genp->PID) == 211) b_constituent_isChargedPion.push_back(true);
      else                       b_constituent_isChargedPion.push_back(false);
    }
    if (object->IsA() == Tower::Class()) {
      auto  np   = (Tower*) object;
      if (np->ET < 1.0) continue;
      float pt2  = np->ET * np->ET;
      float deta = (np->Eta - jet->Eta);
      float dphi = DeltaPhi(np->Phi, jet->Phi);
      sum_dpt        += np->ET;
      sum_dpt2       += pt2;
      sum_deta       += deta*pt2;
      sum_deta2      += deta*deta*pt2;
      sum_dphi       += dphi*pt2;
      sum_dphi2      += dphi*dphi*pt2;
      sum_detadphi   += deta*dphi*pt2;

      sum_n_dpt      += np->ET;
      sum_n_dpt2     += pt2;
      sum_n_deta     += deta*pt2;
      sum_n_deta2    += deta*deta*pt2;
      sum_n_dphi     += dphi*pt2;
      sum_n_dphi2    += dphi*dphi*pt2;
      sum_n_detadphi += deta*dphi*pt2;

      b_constituent_pt.push_back(np->ET); 
      b_constituent_eta.push_back(np->Eta);
      b_constituent_phi.push_back(np->Phi); 
      b_constituent_charge.push_back(0);
      b_constituent_type.push_back(0);
      b_constituent_D0.push_back(-99);
      b_constituent_D0Err.push_back(-99);
      b_constituent_DZ.push_back(-99);
      b_constituent_DZErr.push_back(-99);
      b_constituent_CtgTheta.push_back(-99);
      b_constituent_isFromKs.push_back(false);
      b_constituent_isChargedPion.push_back(false);
      b_constituent_isFromKs.push_back(false);
      b_constituent_isKsFrom.push_back(-99);
      b_constituent_isKsFromTop.push_back(false);
      b_constituent_isKsFromW.push_back(false);
      b_constituent_isKsFromOther.push_back(false);
    }
  }
  float a   = 0, b   = 0, c   = 0;
  float c_a = 0, c_b = 0, c_c = 0;
  float n_a = 0, n_b = 0, n_c = 0;

  if (sum_dpt2 > 0.) {
    b_ptD     = TMath::Sqrt(sum_dpt2)   / sum_dpt;
    b_c_ptD   = TMath::Sqrt(sum_c_dpt2) / sum_c_dpt;
    b_n_ptD   = TMath::Sqrt(sum_n_dpt2) / sum_n_dpt;

    // Calculation of major / minor axis of a jet, https://github.com/cms-sw/cmssw/blob/9834f5dc9ff342ddef08b73d6c294cad36575772/RecoJets/JetProducers/src/MVAJetPuId.cc#L437-L448
    ave_deta  = sum_deta  / sum_dpt2;
    ave_deta2 = sum_deta2 / sum_dpt2;
    ave_dphi  = sum_dphi  / sum_dpt2;
    ave_dphi2 = sum_dphi2 / sum_dpt2;
    a = ave_deta2 - ave_deta*ave_deta;
    b = ave_dphi2 - ave_dphi*ave_dphi;
    c = -(sum_detadphi/sum_dpt2 - ave_deta*ave_dphi);
    float delta = TMath::Sqrt(fabs( (a-b)*(a-b) + 4*c*c ));
    b_axis1 = (a+b+delta > 0) ? sqrt(0.5*(a+b+delta)) : 0;
    b_axis2 = (a+b-delta > 0) ? sqrt(0.5*(a+b-delta)) : 0;

    ave_c_deta  = sum_c_deta  / sum_c_dpt2;
    ave_c_deta2 = sum_c_deta2 / sum_c_dpt2;
    ave_c_dphi  = sum_c_dphi  / sum_c_dpt2;
    ave_c_dphi2 = sum_c_dphi2 / sum_c_dpt2;
    c_a = ave_c_deta2 - ave_c_deta*ave_c_deta;
    c_b = ave_c_dphi2 - ave_c_dphi*ave_c_dphi;
    c_c = -(sum_c_detadphi/sum_c_dpt2 - ave_c_deta*ave_c_dphi);
    float c_delta = TMath::Sqrt(fabs( (c_a-c_b)*(c_a-c_b) + 4*c_c*c_c ));
    b_c_axis1 = (c_a+c_b+c_delta > 0) ? sqrt(0.5*(c_a+c_b+c_delta)) : 0;
    b_c_axis2 = (c_a+c_b-c_delta > 0) ? sqrt(0.5*(c_a+c_b-c_delta)) : 0;

    ave_n_deta  = sum_n_deta  / sum_n_dpt2;
    ave_n_deta2 = sum_n_deta2 / sum_n_dpt2;
    ave_n_dphi  = sum_n_dphi  / sum_n_dpt2;
    ave_n_dphi2 = sum_n_dphi2 / sum_n_dpt2;
    n_a = ave_n_deta2 - ave_n_deta*ave_n_deta;
    n_b = ave_n_dphi2 - ave_n_dphi*ave_n_dphi;
    n_c = -(sum_n_detadphi/sum_n_dpt2 - ave_n_deta*ave_n_dphi);
    float n_delta = TMath::Sqrt(fabs( (n_a-n_b)*(n_a-n_b) + 4*n_c*n_c ));
    b_n_axis1 = (n_a+n_b+n_delta > 0) ? sqrt(0.5*(n_a+n_b+n_delta)) : 0;
    b_n_axis2 = (n_a+n_b-n_delta > 0) ? sqrt(0.5*(n_a+n_b-n_delta)) : 0;

  }
}

void SetHadValues(TTree* tr, int i, TLorentzVector jet_tlv) {
  if (string(tr->GetName()).find("jet") != string::npos) {
    b_had_x = m_recoHad[i].tlv.Pt()/jet_tlv.Pt();
    b_had_dr           = m_recoHad[i].tlv.DeltaR(jet_tlv);
  }
  b_had_pt            = m_recoHad[i].tlv.Pt();
  b_had_eta           = m_recoHad[i].tlv.Eta();
  b_had_phi           = m_recoHad[i].tlv.Phi();
  b_had_mass          = m_recoHad[i].tlv.M();
  b_had_pdgId         = m_recoHad[i].pdgid;
  b_had_isFrom        = m_recoHad[i].isFrom;
  b_had_isPreSelected = m_recoHad[i].isPreSelected;
  b_had_isJetMatched  = m_recoHad[i].isJetMatched;
  b_dau1_pdgId        = m_recoHad[i].dau1_pdgid;
  b_dau1_pt           = m_recoHad[i].dau1_tlv.Pt();
  b_dau1_eta          = m_recoHad[i].dau1_tlv.Eta();
  b_dau1_phi          = m_recoHad[i].dau1_tlv.Phi();
  b_dau1_mass         = m_recoHad[i].dau1_tlv.M();
  b_dau2_pdgId        = m_recoHad[i].dau2_pdgid;
  b_dau2_pt           = m_recoHad[i].dau2_tlv.Pt();
  b_dau2_eta          = m_recoHad[i].dau2_tlv.Eta();
  b_dau2_phi          = m_recoHad[i].dau2_tlv.Phi();
  b_dau2_mass         = m_recoHad[i].dau2_tlv.M();
  auto dau1           = (Track*) tracks->At(m_recoHad[i].dau1_idx);
  auto dau2           = (Track*) tracks->At(m_recoHad[i].dau2_idx);
  b_dau1_D0           = dau1->D0;
  b_dau1_DZ           = dau1->DZ;
  b_dau1_ctgTheta     = dau1->CtgTheta;
  b_dau1_ptErr        = dau1->ErrorPT;
  b_dau1_phiErr       = dau1->ErrorPhi;
  b_dau1_D0Err        = dau1->ErrorD0;
  b_dau1_DZErr        = dau1->ErrorDZ;
  b_dau1_ctgThetaErr  = dau1->ErrorCtgTheta;
  b_dau1_D0Sig        = dau1->D0/dau1->ErrorD0;
  b_dau1_DZSig        = dau1->DZ/dau1->ErrorDZ;
  b_dau2_D0           = dau2->D0;
  b_dau2_DZ           = dau2->DZ;
  b_dau2_ctgTheta     = dau2->CtgTheta;
  b_dau2_ptErr        = dau2->ErrorPT;
  b_dau2_phiErr       = dau2->ErrorPhi;
  b_dau2_D0Err        = dau2->ErrorD0;
  b_dau2_DZErr        = dau2->ErrorDZ;
  b_dau2_ctgThetaErr  = dau2->ErrorCtgTheta;
  b_dau2_D0Sig        = dau2->D0/dau2->ErrorD0;
  b_dau2_DZSig        = dau2->DZ/dau2->ErrorDZ;
  if (m_genHadron.size() != 0) { // Start filling truth matched information
    auto genDau1      = (GenParticle*) dau1->Particle.GetObject();
    auto genDau2      = (GenParticle*) dau2->Particle.GetObject();
    if ((genDau1->M1 == genDau2->M1) && m_genHadron[genDau1->M1].pdgid == m_recoHad[i].pdgid) { 
      // m_recoHad is a set of reconstructed hadron candidates wihch are reconstructed from 2 reco. daughters (Track obejcts, See HadronReconstruction() function) 
      // if reco. hadrons's 2 daughter's gen. information indicates same mother particle(genDau1->M1 == genDau2->M2), then we can say hadron paring is done well.
      // Also, if also gen. mother particle has same pdgId as reco. hadron, then we call this reco. hadron truth-matched hadron.
      auto genp = (const GenParticle*) particles->At(genDau1->M1);
      b_had_nMatched      = 2;
      b_genHad_x          = m_genHadron[genDau1->M1].x;
      b_genHad_pdgId      = m_genHadron[genDau1->M1].pdgid;
      b_genHad_pt         = m_genHadron[genDau1->M1].tlv.Pt();
      b_genHad_eta        = m_genHadron[genDau1->M1].tlv.Eta();
      b_genHad_phi        = m_genHadron[genDau1->M1].tlv.Phi();
      b_genHad_mass       = m_genHadron[genDau1->M1].tlv.M();
      b_genHad_isFrom     = m_genHadron[genDau1->M1].isFrom;
      b_genHad_isFromTop  = m_genHadron[genDau1->M1].isFromTop;
      b_genHad_isFromW    = m_genHadron[genDau1->M1].isFromW;
      b_genHad_isPU       = genp->IsPU;
      b_genHad_D0         = genp->D0;
      b_genHad_DZ         = genp->DZ;
      b_genHad_ctgTheta   = genp->CtgTheta;
      b_genDau1_pdgId     = m_genHadron[genDau1->M1].dau1_pdgid;
      b_genDau1_pt        = m_genHadron[genDau1->M1].dau1_tlv.Pt();
      b_genDau1_eta       = m_genHadron[genDau1->M1].dau1_tlv.Eta();
      b_genDau1_phi       = m_genHadron[genDau1->M1].dau1_tlv.Phi();
      b_genDau1_mass      = m_genHadron[genDau1->M1].dau1_tlv.M();
      b_genDau1_D0        = genDau1->D0;
      b_genDau1_DZ        = genDau1->DZ;
      b_genDau1_ctgTheta  = genDau1->CtgTheta;
      b_genDau2_pdgId     = m_genHadron[genDau1->M1].dau2_pdgid;
      b_genDau2_pt        = m_genHadron[genDau1->M1].dau2_tlv.Pt();
      b_genDau2_eta       = m_genHadron[genDau1->M1].dau2_tlv.Eta();
      b_genDau2_phi       = m_genHadron[genDau1->M1].dau2_tlv.Phi();
      b_genDau2_mass      = m_genHadron[genDau1->M1].dau2_tlv.M();
      b_genDau2_D0        = genDau2->D0;
      b_genDau2_DZ        = genDau2->DZ;
      b_genDau2_ctgTheta  = genDau2->CtgTheta;
    } else if ( (m_genHadron[genDau1->M1].idx == genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) || (m_genHadron[genDau1->M1].idx == genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) ) {
      b_had_nMatched = 1;
    } else if ( (m_genHadron[genDau1->M1].idx != genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) ) {
      b_had_nMatched = 0;
    }
  }
}


std::vector<float> collectHadron(std::vector<GenParticle> hadInJet, bool motherCheck){
}

void SetMVAReader() {
  #define TMVABranch_(reader, name) reader->AddVariable(#name, &(b_##name));
  #define hadTMVABranch(name) TMVABranch_(m_hadReader,name);
  #define jetTMVABranch(name) TMVABranch_(m_jetReader,name);
  #define jksTMVABranch(name) TMVABranch_(m_jksReader,name);

  m_hadReader = new TMVA::Reader();
  hadTMVABranch(had_pt);       
  hadTMVABranch(had_eta);      
  hadTMVABranch(had_phi);      
  hadTMVABranch(dau1_pt);      
  hadTMVABranch(dau1_eta);     
  hadTMVABranch(dau1_phi);     
  hadTMVABranch(dau1_ctgTheta);
  hadTMVABranch(dau1_D0);      
  hadTMVABranch(dau1_DZ);      
  hadTMVABranch(dau1_D0Sig);   
  hadTMVABranch(dau1_DZSig);   
  hadTMVABranch(dau2_pt);      
  hadTMVABranch(dau2_eta);     
  hadTMVABranch(dau2_phi);     
  hadTMVABranch(dau2_ctgTheta);
  hadTMVABranch(dau2_D0);      
  hadTMVABranch(dau2_DZ);      
  hadTMVABranch(dau2_D0Sig);   
  hadTMVABranch(dau2_DZSig);   
  m_hadReader->BookMVA("KS_BDT", tmvaWeight_Ks);
  m_hadReader->BookMVA("Lambda_BDT", tmvaWeight_Lambda);
  //m_jetReader = new TMVA::Reader();
  //jetTMVABranch(ptD); jetTMVABranch(area); jetTMVABranch(CSVV2);
  //m_jetReader->BookMVA("Jet_BDT_highest", tmvaWeight_Jet);

  //m_jksReader = new TMVA::Reader();
  //m_jksReader->BookMVA("JKS_BDT_highest", tmvaWeight_JKS);
}

