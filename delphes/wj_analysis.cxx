#include <stdio.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "classes/DelphesClasses.h"
#include "wj_analysis.h"

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
    HadronPreselection(m_recoHad);
    FindTruthMatchedHadron();
    FillJetTree(jets, jettr, "pT");
    FillHadTree(hadtr);
    //cout << " ============================================ " << endl;
    outtr->Fill();
    //cout << "chk 9 " << endl;
  }
  inFile->Close();
  cutflow->Write();
  //cout << "chk 11 " << endl;

  out->Write();
  out->Close();
  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

void ResetBranch(){
    m_genQuark.clear();    m_genLepton.clear();  m_genHadron.clear();
    m_selectedJet.clear(); m_matchedJet.clear();
    m_recoHad.clear();

    b_recoLep1.SetPtEtaPhiM(0,0,0,0); b_recoLep2.SetPtEtaPhiM(0,0,0,0);
    b_recoLep1_pdgId = -99; b_recoLep2_pdgId = -99;

    // MatchingGenJet()
    b_dilepton_mass = -99; b_dilepton_ch = 0; b_step = 0; 

    b_channel = -1;
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

void EventSelection(TH1F * cutflow, UInt_t decaychannel){
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

    // jet selection
    m_selectedJet = JetSelection(jets, recoLep);

    // b-jet selection
    std::vector<Jet*> selectedBJets;
    for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){
      if (jet->second->BTag) selectedBJets.push_back(jet->second);
    }

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
      if (recoEl.size() == 1 && recoMu.size() == 0) {b_recoLep1 = recoEl[0].tlv; b_recoLep1_pdgId = recoEl[0].pdgid ;}
      if (recoEl.size() == 0 && recoMu.size() == 1) {b_recoLep1 = recoMu[0].tlv; b_recoLep1_pdgId = recoMu[0].pdgid ;}

      // step 2 : veto additional lepton
      if ( (recoMu.size() == 0 && vetoMu.size() > 0) || (recoMu.size() == 1 && vetoMu.size() >1) ) {b_step = 1; return;}
      if ( (recoEl.size() == 0 && vetoEl.size() > 0) || (recoEl.size() == 1 && vetoEl.size() >1) ) {b_step = 1; return;}
      cutflow->Fill(2);

      // step 3 : jet selection 
      if (m_selectedJet.size() < cut_nJet) {b_step = 2; return;}
      cutflow->Fill(3);

      // step 4 : b-jet selection
      if (selectedBJets.size() < cut_nBJet) {b_step = 3; return;}
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

      auto dilepton = recoLep[0].tlv + recoLep[1].tlv;
      auto mulpdg   = recoLep[0].charge * recoLep[1].charge;
      b_dilepton_ch   = abs(recoLep[0].pdgid) + abs(recoLep[1].pdgid); // 22 -> ee , 24 -> emu , 26 -> mumu
      b_dilepton_mass = dilepton.M();

      if (b_dilepton_mass < cut_Dilep_M || mulpdg > 0) {b_step = 0; return;}
      cutflow->Fill(1);
  
      // step 2 : veto Z boson
      if ( (b_dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) {b_step = 1; return;}
      cutflow->Fill(2);
  
      // step 3 : MET cut
      auto met = ((MissingET*) missingET->At(0))->MET;
      if ( (b_dilepton_ch != 24) && (met < cut_MET_Pt) ) {b_step = 2; return;}
      cutflow->Fill(3);
  
      // step 4 : jet selection
      if (m_selectedJet.size() < cut_nJet) {b_step = 3; return;}
      cutflow->Fill(4);
  
      // step 5 : b-jet selection (not used)
      if (selectedBJets.size() < cut_nBJet) {b_step = 4; return;}
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

    auto quark     = (const GenParticle*) particles->At(top->D1); // gen quark from top-quark
    auto Wboson    = (const GenParticle*) particles->At(top->D2); // gen W-boson from top-quark
    auto lastBoson = getLast(particles, Wboson); // Get a last copy of W-boson in the decay chain
    m_genQuark.push_back(quark);

    auto dau1_from_W = (const GenParticle*) particles->At(lastBoson->D1); // for dilep, always neutrino
    auto dau2_from_W = (const GenParticle*) particles->At(lastBoson->D2); // for dilep, always e/mu/tau

    struct Lepton lep1 = toLepton(dau1_from_W);
    struct Lepton lep2 = toLepton(dau2_from_W);

    // Collect leptons from W-boson, for dilepton, the size of m_genLepton should be 2 and for semilepton, 1
    if (abs(lep1.pdgid) == 11 || abs(lep1.pdgid) == 13) m_genLepton.push_back(lep1);
    if (abs(lep2.pdgid) == 11 || abs(lep2.pdgid) == 13) m_genLepton.push_back(lep2);
  }

  if ( m_genQuark.size() < 2 ) {
    cout << "No two quarks from top quarks ... skip this event" << endl; 
    return;
  } else if ( (m_genLepton.size() != 2 && m_decayChannel == "di") || (m_genLepton.size() != 1 && m_decayChannel == "semi")) {
    cout << "No two reco leptons in dileptonic channel or no reco lepton in semileptonic channel ===> Maybe because a case that W boson decay into tau / neutrino isn't included" << endl;
    return;
  }

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

  for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){ // Loop over jets that passed jet selection
    auto jIdx = jet->first; 
    auto selJet = jet->second;
    TLorentzVector jet_tlv, tq1_tlv, tq2_tlv;
    jet_tlv.SetPtEtaPhiM(selJet->PT, selJet->Eta, selJet->Phi, selJet->Mass);
    tq1_tlv.SetPtEtaPhiM(m_genQuark[0]->PT, m_genQuark[0]->Eta, m_genQuark[0]->Phi, m_genQuark[0]->Mass);
    tq2_tlv.SetPtEtaPhiM(m_genQuark[1]->PT, m_genQuark[1]->Eta, m_genQuark[1]->Phi, m_genQuark[1]->Mass);

    if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   { wl1_tlv = m_genLepton[0].tlv; wl2_tlv = m_genLepton[1].tlv;}
    if (m_decayChannel == "semilepton" || m_decayChannel == "semi") { wl1_tlv = m_genLepton[0].tlv;}

    // deltaR matching between reco. jet and quark from top-quark (for m_matchedJet, see struct RecoJet in header file)
    if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap && jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap ) { // Overlap case ==> we will not use this jet
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), -99,                          true,  -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), -99,                     -99,                          true,  -99, false, false}); 
    } else if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap ) { // jet matched to quark1 from top-quark
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[0]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[0]->PT, false, -99, false, false}); 
    } else if (jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap ) { // jet matched to quark2 from to-quark
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
    }
  }

  if (m_matchedJet.size() != 0) {
    // Sort by pT and dR(lep, jet) and then set a flag for jets with highest top 2 pT or closest distance to lepton
    sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.j->PT > b.j->PT); } ); // Sort with pT ordering
    if (m_matchedJet.size() > 0) m_matchedJet[0].hasHighestPt  = true;
    if (m_matchedJet.size() > 1) m_matchedJet[1].hasHighestPt  = true;
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
    //if (hadCand1->PT < tkPTCut_) continue;
    //if (fabs(hadCand1->D0/hadCand1->ErrorD0) < tkIPSigXYCut_) continue; // Cut of siginifcance of transverse impact parametr
    for (auto j = 0; j < tracks->GetEntries(); ++j) {
      auto hadCand2 = (Track*) tracks->At(j);
      if (hadCand2->Charge != -1) continue; // Only pick negative charge
      if (abs(hadCand2->PID) == 11 || abs(hadCand2->PID) == 13) continue; // Exclude leptons
      //if (hadCand2->PT < tkPTCut_) continue;
      //if (fabs(hadCand2->D0/hadCand2->ErrorD0) < tkIPSigXYCut_) continue; // Cut of siginifcance of transverse impact parametr

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

void HadronPreselection(std::vector<RecoHad> recoHad) {
  // Preselection cuts for reconstructed hadrons
  // Not fully implemented yet
  for (auto i=0; i < recoHad.size(); ++i) {
    if (recoHad[i].tlv.Eta() > 2.5) continue;
    if (recoHad[i].dau1_tlv.Pt() < 0.95 || recoHad[i].dau2_tlv.Pt() < 0.95) continue;
    //if (recoHad[i].dau1_D0Sig    < 5.   || recoHad[i].dau2_D0Sig    < 5.)   continue;
    else recoHad[i].isSelected = true;
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
    if (p->PID != 310 && p->PID != 3122) continue; // We are only intrested in some hadrons (Ks(pdgId == 310) and Lambda0(pdgId == 3122))
    if (p->D1  == -1  || p->D2  == -1)   continue; // We need a gen. level link to their daughters
    // Get gen. daughter objects from gen hadrons
    auto d1 = (const GenParticle*) particles->At(p->D1);
    auto d2 = (const GenParticle*) particles->At(p->D2);
    if ( p->PID == 310  && (abs(d1->PID) != 211 || abs(d2->PID) != 211) ) continue; // Pick only charged pion pair, which is most dominant decay process (See http://pdg.lbl.gov/2019/listings/rpp2019-list-K-zero-S.pdf)
    if ( p->PID == 3122 ) { // Pick only proton - pion pair, which is most dominant decay process (See http://pdg.lbl.gov/2019/listings/rpp2019-list-lambda.pdf) 
      if      ( (abs(d1->PID) != 211  && abs(d1->PID) != 2212) || (abs(d2->PID) != 211  && abs(d2->PID) != 2212) ) continue;
      else if ( abs(d1->PID) == abs(d2->PID) ) continue;
    }
    if (d1->PID*d2->PID > 0) continue; // Check opposite charge 

    auto motherList = getMlist(particles, p); // Get All mothers of gen. hadron
    for (auto j=0; j < motherList.size()-1; ++j) {
      if ( abs(motherList[j]->PID) == 24 && !isFromW) { // Check whether gen. hadron comes from W-boson or not
        isFromW = true;
        isFrom  = motherList[j-1]->PID; // Check flavour of quark from which gen. hadron comes
        x       = p->PT/motherList[j-1]->PT;
      }
      if ( abs(motherList[j]->PID) == 6 && motherList[j]->Status == 62 ) { // Check if gen. hadron comes from top-quark
        isFromTop = true; 
        if (isFromW == false) { isFrom = motherList[j-1]->PID; x = p->PT/motherList[j-1]->PT; }
        break;
      }
    }
    // Can't trust P4() function ...
    TLorentzVector p_tlv, d1_tlv, d2_tlv;
    p_tlv.SetPtEtaPhiM(p->PT, p->Eta, p->Phi, p->Mass);
    d1_tlv.SetPtEtaPhiM(d1->PT, d1->Eta, d1->Phi, d1->Mass);
    d2_tlv.SetPtEtaPhiM(d2->PT, d2->Eta, d2->Phi, d2->Mass);
    
    m_genHadron[i] = {p_tlv, d1_tlv, d2_tlv, x, i, p->PID, p->D1, d1->PID, p->D2, d2->PID, isFrom, isFromTop, isFromW};
  }
}

int FindMatchedHadron(TLorentzVector jet_tlv, TString method) {
  int hIdx = -1;
  float maxBDTScore = -1;
  float highestPt   = -1;
  float lowestDot   = 99999999.;
  float closestDr   = -1;
  for (auto i=0; i<m_recoHad.size(); ++i) {
    auto had_tlv = m_recoHad[i].tlv;
    auto dR      = had_tlv.DeltaR(jet_tlv);
    if (abs(m_recoHad[i].pdgid) != 310 && abs(m_recoHad[i].pdgid) != 3122) continue;
    if (m_hadPreCut && !m_recoHad[i].isSelected) continue;
    if (dR > cut_JetConeSizeOverlap) continue;

    if (method == "BDT") { // The method isn't constructed yet
      cout << "BDT ==> Work in progress" << endl;
      // auto hadBDTScore = m_hadReader->EvaluateMVA("KS_BDT");
      // if (hadBDTScore > maxBDTScore) {
      //   maxBDTScore = hadBDTScore;
      //   hIdx = i;
      // }
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
  if (hIdx != -1) m_recoHad[hIdx].isJetMatched = true;
  return hIdx;
} 

void FillJetTree(TClonesArray* jets, TTree* jettr, TString matchingMethod) {
  for ( unsigned i = 0; i < jets->GetEntries(); ++i){
    ResetJetValues();
    auto jet     = (Jet*) jets->At(i);
    TLorentzVector jet_tlv;
    jet_tlv.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);

    if (m_selectedJet[i] != 0) { // check if the i-th jet is selected jet from event selection 
      b_isSelectedJet = true;
    }
    for ( auto j = 0; j < m_matchedJet.size(); ++j ) {
      if (m_matchedJet[j].idx == i) { // check if the i-th jet is jet which matched to gen s/b quark from top-quark
        b_isFrom        = m_matchedJet[j].gen_pdgid;
        b_hasHighestPt  = m_matchedJet[j].hasHighestPt;
        b_hasClosestLep = m_matchedJet[j].hasClosestLep;
      }
    }
    auto hIdx = FindMatchedHadron(jet_tlv, matchingMethod); // check if reco. hadron(Ks or lambda0) is inside the i-th jet according to a matching method
    SetJetValues(jet);
    if (hIdx != -1) {
      SetHadValues(jettr, hIdx, jet_tlv); // if there is a reco. hadron within the i-th jet, then get information of reco. hadron
    }
    jettr->Fill();
  }
}

void FillHadTree(TTree* hadtr) {
  for ( unsigned i = 0; i < m_recoHad.size(); ++i){
    ResetHadValues();
    SetHadValues(hadtr, i);
    hadtr->Fill();
  }
}

void SetJetValues(Jet* jet) {
  b_pt               = jet->PT;
  b_eta              = jet->Eta;
  b_phi              = jet->Phi;
  b_mass             = jet->Mass;
  b_radiusInEta      = jet->DeltaEta;
  b_radiusInPhi      = jet->DeltaPhi;
  b_ptD              = jet->PTD;
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
}

void SetHadValues(TTree* tr, int i, TLorentzVector jet_tlv) {
  if (string(tr->GetName()).find("jet") != string::npos) {
    b_had_x = m_recoHad[i].tlv.Pt()/jet_tlv.Pt();
    b_had_dr           = m_recoHad[i].tlv.DeltaR(jet_tlv);
  }
  b_had_pt           = m_recoHad[i].tlv.Pt();
  b_had_eta          = m_recoHad[i].tlv.Eta();
  b_had_phi          = m_recoHad[i].tlv.Phi();
  b_had_mass         = m_recoHad[i].tlv.M();
  b_had_pdgId        = m_recoHad[i].pdgid;
  b_had_isFrom       = m_recoHad[i].isFrom;
  b_had_isSelected   = m_recoHad[i].isSelected;
  b_had_isJetMatched = m_recoHad[i].isJetMatched;
  b_dau1_pdgId       = m_recoHad[i].dau1_pdgid;
  b_dau1_pt          = m_recoHad[i].dau1_tlv.Pt();
  b_dau1_eta         = m_recoHad[i].dau1_tlv.Eta();
  b_dau1_phi         = m_recoHad[i].dau1_tlv.Phi();
  b_dau1_mass        = m_recoHad[i].dau1_tlv.M();
  b_dau2_pdgId       = m_recoHad[i].dau2_pdgid;
  b_dau2_pt          = m_recoHad[i].dau2_tlv.Pt();
  b_dau2_eta         = m_recoHad[i].dau2_tlv.Eta();
  b_dau2_phi         = m_recoHad[i].dau2_tlv.Phi();
  b_dau2_mass        = m_recoHad[i].dau2_tlv.M();
  auto dau1          = (Track*) tracks->At(m_recoHad[i].dau1_idx);
  auto dau2          = (Track*) tracks->At(m_recoHad[i].dau2_idx);
  b_dau1_D0          = dau1->D0;
  b_dau1_DZ          = dau1->DZ;
  b_dau1_ctgTheta    = dau1->CtgTheta;
  b_dau1_ptErr       = dau1->ErrorPT;
  b_dau1_phiErr      = dau1->ErrorPhi;
  b_dau1_D0Err       = dau1->ErrorD0;
  b_dau1_DZErr       = dau1->ErrorDZ;
  b_dau1_ctgThetaErr = dau1->ErrorCtgTheta;
  b_dau1_D0Sig       = dau1->D0/dau1->ErrorD0;
  b_dau1_DZSig       = dau1->DZ/dau1->ErrorDZ;
  b_dau2_D0          = dau2->D0;
  b_dau2_DZ          = dau2->DZ;
  b_dau2_ctgTheta    = dau2->CtgTheta;
  b_dau2_ptErr       = dau2->ErrorPT;
  b_dau2_phiErr      = dau2->ErrorPhi;
  b_dau2_D0Err       = dau2->ErrorD0;
  b_dau2_DZErr       = dau2->ErrorDZ;
  b_dau2_ctgThetaErr = dau2->ErrorCtgTheta;
  b_dau2_D0Sig       = dau2->D0/dau2->ErrorD0;
  b_dau2_DZSig       = dau2->DZ/dau2->ErrorDZ;
  if (m_genHadron.size() != 0) { // Start filling truth matched information
    auto genDau1     = (GenParticle*) dau1->Particle.GetObject();
    auto genDau2     = (GenParticle*) dau2->Particle.GetObject();
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
