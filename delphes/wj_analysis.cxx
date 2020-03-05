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

  //read input file
  cout << "Input File : " << inf.c_str() << "\nOutput File : " << outf.c_str() << "\nDecay Channel : " << m_decayChannel.c_str() << endl;
  int decay_ch = 0;
  if (m_decayChannel == "semilepton" || m_decayChannel == "semi") decay_ch = 1;
  if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   decay_ch = 2;

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

  //make out file
  auto out   = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("event", "event");
  auto jettr = new TTree("MVA_jet", "MVA_jet");
  auto hadtr = new TTree("MVA_had", "MVA_had");


  DefBranch(outtr);
  SetJetBranch(jettr);
  SetHadBranch(hadtr);

  TString cutflow_title   = "cutflow" + inf;
  TH1F * cutflow          = new TH1F("cutflow", cutflow_title, 7,-1,6); // -1 : all events, 0 : events after lepton selection, 1~ : events after step
  TH1F * histo_diHadron_S = new TH1F("diHadron_mass_S", "diHadron_mass_S", 300, 0, 3000);
  TH1F * histo_diHadron_B = new TH1F("diHadron_mass_B", "diHadron_mass_B", 300, 0, 3000);
  TH1F * histo_x_Ks       = new TH1F("hist_x_Ks", "hist_x_Ks", 100, 0, 1);
  TH1F * histo_rho_Ks     = new TH1F("hist_rho_Ks", "hist_rho_Ks", 1000, 0, 300);
  TH1F * histo_d_Ks       = new TH1F("hist_d_Ks", "hist_d_Ks", 1000, 0, 20);
  TH1F * histo_x_lamb     = new TH1F("hist_x_lamb", "hist_x_lamb", 100, 0, 1);
  TH1F * histo_rho_lamb   = new TH1F("hist_rho_lamb", "hist_rho_lamb", 1800, 0, 18000);
  TH1F * histo_d_lamb     = new TH1F("hist_d_lamb", "hist_d_lamb", 1000, 0, 100);

  //Event Loop Start!
  for (size_t iev = 0; iev < inTree->GetEntries(); ++iev){
    if (iev%1000 == 0 ) cout << "event check    iev    ----> " << iev << endl;
    inTree->GetEntry(iev);
    ResetBranch();

    EventSelection(cutflow, decay_ch);
    if (b_step < 4) continue;
    //cout << " ============================================ " << endl;
    MatchingGenJet(histo_diHadron_S, histo_diHadron_B);
    HadronReconstruction();
    HadronPreselection(m_recoHad);
    FindTruthMatchedHadron();
    FillJetTree(jets, jettr, "pT");
    FillHadTree(hadtr);
    //cout << " ============================================ " << endl;
    /*
    // pion track test
    for (unsigned i = 0; i < tracks->GetEntries(); ++ i){
      auto track = (Track*) tracks->At(i);
      if (abs(track->PID) == 321){
        auto genTrack = (GenParticle*) track->Particle.GetObject();
        auto genTrack_mom = (const GenParticle*) particles->At(genTrack->M1);
        cout << genTrack_mom->PID << endl;
        float pionR = sqrt(pow(track->X,2)+pow(track->Y,2)+pow(track->Z,2));
        float pionRd = sqrt(pow(track->Xd,2)+pow(track->Yd,2)+pow(track->Zd,2));
        float pionROuter = sqrt(pow(track->XOuter,2)+pow(track->YOuter,2)+pow(track->ZOuter,2));
      }
    }
    */
    outtr->Fill();
  }

  inFile->Close();

  //outtr->Write();

  //outtr->SetBranchAddress("x_Ks", &b_Ks_x);
  //outtr->SetBranchAddress("x_lamb", &b_lamb_x);
  //outtr->SetBranchAddress("rho_Ks", &b_Ks_rho);
  //outtr->SetBranchAddress("rho_lamb", &b_lamb_rho);
  //outtr->SetBranchAddress("d_Ks", &b_Ks_d);
  //outtr->SetBranchAddress("d_lamb", &b_lamb_d);
  //for (size_t ent = 0; ent < outtr->GetEntries(); ++ent){
  //  outtr->GetEntry(ent);
  //  histo_x_Ks->Fill(b_Ks_x);
  //  histo_rho_Ks->Fill(b_Ks_rho);
  //  histo_d_Ks->Fill(b_Ks_d);
  //  histo_x_lamb->Fill(b_lamb_x);
  //  histo_rho_lamb->Fill(b_lamb_rho);
  //  histo_d_lamb->Fill(b_lamb_d);
  //}
  //histo_x_Ks->Write();
  //histo_rho_Ks->Write();
  //histo_d_Ks->Write();
  //histo_x_lamb->Write();
  //histo_rho_lamb->Write();
  //histo_d_lamb->Write();

  cutflow->Write();

  //histo_diHadron_S->Write();
  //histo_diHadron_B->Write();

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

    b_Ks_x   = -1; b_Ks_x_S   = -1; b_Ks_x_B   = -1; b_Ks_rho   = -99;  b_Ks_rho_S   = -99;  b_Ks_rho_B   = -99;  b_Ks_d   = -99; b_Ks_d_S   = -99; b_Ks_d_B   = -99;
    b_lamb_x = -1; b_lamb_x_S = -1; b_lamb_x_B = -1; b_lamb_rho = -999; b_lamb_rho_S = -999; b_lamb_rho_B = -999; b_lamb_d = -99; b_lamb_d_S = -99; b_lamb_d_B = -99;

    b_significance_S = -99; b_significance_B = -99;

    b_channel = -1;

    b_gen_step0 = false; b_gen_step1 = false;

    b_jet_pt.clear(); b_jet_eta.clear(); b_jet_phi.clear(); b_jet_energy.clear();
    b_jet_pdgId.clear();

    b_KsInJet_pt.clear();   b_KsInJet_eta.clear();   b_KsInJet_phi.clear();   b_KsInJet_energy.clear();   b_KsInJet_R.clear();   b_KsInJet_outR.clear();   b_KsInJet_rho.clear();   b_KsInJet_d.clear();
    b_lambInJet_pt.clear(); b_lambInJet_eta.clear(); b_lambInJet_phi.clear(); b_lambInJet_energy.clear(); b_lambInJet_R.clear(); b_lambInJet_outR.clear(); b_lambInJet_rho.clear(); b_lambInJet_d.clear();
    b_lepInJet_pt.clear();  b_lepInJet_eta.clear();  b_lepInJet_phi.clear();  b_lepInJet_energy.clear();  b_lepInJet_R.clear();

    b_nKsInJet.clear(); b_nLambInJet.clear(); b_nLepInJet.clear();

    b_jet1_diHadron_mass.clear(); b_jet2_diHadron_mass.clear();
}

//std::vector<Jet*> 
std::map<int, Jet*> JetSelection(TClonesArray* jets, std::vector<struct Lepton> recoLep) {
  //std::vector<Jet*>
  std::map<int, Jet*> selectedJets;
  for ( unsigned k = 0; k < jets->GetEntries(); ++k){
    auto jet = (Jet*) jets->At(k);
    if (jet->PT < cut_JetPt) continue;
    if (abs(jet->Eta) > cut_JetEta) continue;
    bool hasOverlap = false;
    for (auto lep : recoLep) { if ((jet->P4()).DeltaR(lep.tlv) < 0.4) hasOverlap = true; }
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

    SetCutValues(decaychannel); // decaychannel == 1 : Semileptonic decay channel / decaychannel == 2 : Dileptonic decay channel

    // get recoLep vector
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
    //for (auto jet : m_selectedJet){
    for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){
      if (jet->second->BTag) selectedBJets.push_back(jet->second);
    }

    // Semi-lepton case
    if (decaychannel == 1) {
      // get vetoLep vector
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

void MatchingGenJet(TH1F * histo_diHadron_S, TH1F * histo_diHadron_B) {
  std::vector<const GenParticle*> genTops;
  std::vector<std::vector<float> > jet_diHadron_mass;

  genTops.clear();
  jet_diHadron_mass.clear();

  for (int i = 0; i < particles->GetEntries(); ++ i){
    auto p = (const GenParticle*) particles->At(i);
    if (p->Status   >  30) continue;
    if (abs(p->PID) != 6)  continue;
    auto top = getLast(particles,p);
    genTops.push_back(top);

    auto quark     = (const GenParticle*) particles->At(top->D1);
    auto Wboson    = (const GenParticle*) particles->At(top->D2);
    auto lastBoson = getLast(particles, Wboson);
    m_genQuark.push_back(quark);

    auto dau1_from_W = (const GenParticle*) particles->At(lastBoson->D1); // for dilep, always neutrino
    auto dau2_from_W = (const GenParticle*) particles->At(lastBoson->D2); // for dilep, always e/mu/tau

    struct Lepton lep1 = toLepton(dau1_from_W);
    struct Lepton lep2 = toLepton(dau2_from_W);

    if (abs(lep1.pdgid) == 11 || abs(lep1.pdgid) == 13 || abs(lep1.pdgid) == 15) m_genLepton.push_back(lep1);
    if (abs(lep2.pdgid) == 11 || abs(lep2.pdgid) == 13 || abs(lep2.pdgid) == 15) m_genLepton.push_back(lep2);

    //cout << m_genLepton.size() << " : quark from top : " << setw(10) << quark->PID << " | W boson : " << setw(10) << lastBoson->PID << " | dau1 : " << setw(10) << dau1_from_W->PID << " | dau2 : " << setw(10) << dau2_from_W->PID << endl;
  }

  if ( m_genQuark.size() < 2 ) {
    cout << "No two quarks from top quarks ... skip this event" << endl; 
    return;
  }

  TLorentzVector wl1_tlv;
  TLorentzVector wl2_tlv;
  if (m_decayChannel == "dilepton" || m_decayChannel == "di") {//m_genLepton.size() == 2){
    if (abs(m_genLepton[0].pdgid) == abs(m_genLepton[1].pdgid)) {
      if      (abs(m_genLepton[0].pdgid) == 11) b_channel = 1; // ee channel
      else if (abs(m_genLepton[0].pdgid) == 13) b_channel = 2; // mm channel
      else if (abs(m_genLepton[0].pdgid) == 15) b_channel = 3; // tau
    } else b_channel = 0; // em channel
  } else if (m_decayChannel == "semilepton" || m_decayChannel == "semi") {//else if (m_genLepton.size() == 1) {
    if      (abs(m_genLepton[0].pdgid) == 11) b_channel = 1; // e+jet channel
    else if (abs(m_genLepton[0].pdgid) == 13) b_channel = 2; // mu+jet channel
    else if (abs(m_genLepton[0].pdgid) == 15) b_channel = 3; // tau+jet channel
  }

  //for (unsigned int i = 0; i < m_selectedJet.size(); ++i) {
  for (auto jet = m_selectedJet.begin(); jet != m_selectedJet.end(); ++jet){
    auto jIdx = jet->first; 
    auto selJet = jet->second;
    TLorentzVector jet_tlv = selJet->P4();
    auto tq1_tlv = m_genQuark[0]->P4();
    auto tq2_tlv = m_genQuark[1]->P4();

    if (m_decayChannel == "dilepton"   || m_decayChannel == "di") { wl1_tlv = m_genLepton[0].tlv; wl2_tlv = m_genLepton[1].tlv;}
    if (m_decayChannel == "semilepton" || m_decayChannel == "semi") { wl1_tlv = m_genLepton[0].tlv;}

    if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap && jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap ) {
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), -99,                          true,  -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, -99,                -99,                     jet_tlv.DeltaR(wl1_tlv), -99,                     -99,                          true,  -99, false, false}); 
    } else if (jet_tlv.DeltaR(tq1_tlv) <= cut_JetConeSizeOverlap ) {
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[0]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[0]->PID, jet_tlv.DeltaR(tq1_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[0]->PT, false, -99, false, false}); 
    } else if (jet_tlv.DeltaR(tq2_tlv) <= cut_JetConeSizeOverlap ) {
      if (m_decayChannel == "dilepton"   || m_decayChannel == "di")   m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), jet_tlv.DeltaR(wl2_tlv), selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
      if (m_decayChannel == "semilepton" || m_decayChannel == "semi") m_matchedJet.push_back({jIdx, selJet, selJet->Flavor, m_genQuark[1]->PID, jet_tlv.DeltaR(tq2_tlv), jet_tlv.DeltaR(wl1_tlv), -99,                     selJet->PT/m_genQuark[1]->PT, false, -99, false, false});
    }

    //cout << "selected Jet idx : "                                     << setw(2)  << jIdx 
    //     << " dR (tq1 = " << setw(2) << m_genQuark[0]->PID << " ) : " << setw(10) << jet_tlv.DeltaR(tq1_tlv) 
    //     << " dR (tq2 = " << setw(2) << m_genQuark[1]->PID << " ) : " << setw(10) << jet_tlv.DeltaR(tq2_tlv) 
    //     << " jet PID : "                                             << setw(2)  << selJet->Flavor 
    //     << " Gen PT 1 : "                                            << setw(10) << m_genQuark[0]->PT 
    //     << " Gen PT 2 : "                                            << setw(10) << m_genQuark[1]->PT 
    //     << " Jet PT : "                                              << setw(10) << selJet->PT 
    //     << " Jet PT / Gen1 PT : "                                    << setw(10) << selJet->PT/m_genQuark[0]->PT 
    //     << " Jet PT / Gen2 PT : "                                    << setw(10) << selJet->PT/m_genQuark[1]->PT  
    //     << " no. matched jet : "                                     << setw(2)  << m_matchedJet.size() 
    //     << endl;
  }

  sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.j->PT > b.j->PT); } ); // sort with pT ordering
  if (m_matchedJet.size() > 0) m_matchedJet[0].hasHighestPt  = true;
  if (m_matchedJet.size() > 1) m_matchedJet[1].hasHighestPt  = true;
  sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.drl1 < b.drl1); } ); // sort with dR(lep, jet) ordering
  m_matchedJet[0].hasClosestLep = true;
  if (m_decayChannel == "dilepton" || m_decayChannel == "di") sort(m_matchedJet.begin(), m_matchedJet.end(), [] (RecoJet a, RecoJet b) { return (a.drl2 < b.drl2); } ); // sort with dR(lep, jet) ordering
  m_matchedJet[1].hasClosestLep = true;

  //for (auto i =0; i < m_matchedJet.size(); ++i) {
  //  cout << " >>>>> matched Jet idx : " << setw(2)  << m_matchedJet[i].idx 
  //       << " dR : "                    << setw(10) << m_matchedJet[i].dr 
  //       << " pT : "                    << setw(10) << m_matchedJet[i].j->PT 
  //       << " x : "                     << setw(10) << m_matchedJet[i].x 
  //       << " jet PID : "               << setw(2)  << m_matchedJet[i].pdgid 
  //       << " no. matched jet : "       << setw(2)  << m_matchedJet.size() << endl;
  //}
}

void HadronReconstruction() {
  //cout << " >>>>>>>>>>>>>>>>>>>>>>> HADRON RECO <<<<<<<<<<<<<<<<<<<<<<<< " << endl;
  // Pick two charged hadron pairs
  for (auto i = 0; i < tracks->GetEntries(); ++i) {
    auto hadCand1 = (Track*) tracks->At(i);
    if (hadCand1->Charge != 1 || ( abs(hadCand1->PID) == 11 || abs(hadCand1->PID) == 13) ) continue;
    for (auto j = 0; j < tracks->GetEntries(); ++j) {
      auto hadCand2 = (Track*) tracks->At(j);
      if (hadCand2->Charge != -1 || ( abs(hadCand2->PID) == 11 || abs(hadCand2->PID) == 13) ) continue; 

      TLorentzVector dauCand1_pion_tlv;   
      TLorentzVector dauCand1_proton_tlv; 
      TLorentzVector dauCand2_pion_tlv;   
      TLorentzVector dauCand2_proton_tlv; 

      dauCand1_pion_tlv.SetPtEtaPhiM(  hadCand1->PT, hadCand1->Eta, hadCand1->Phi, pionMass_);
      dauCand1_proton_tlv.SetPtEtaPhiM(hadCand1->PT, hadCand1->Eta, hadCand1->Phi, protonMass_);
      dauCand2_pion_tlv.SetPtEtaPhiM(  hadCand2->PT, hadCand2->Eta, hadCand2->Phi, pionMass_);
      dauCand2_proton_tlv.SetPtEtaPhiM(hadCand2->PT, hadCand2->Eta, hadCand2->Phi, protonMass_);

      TLorentzVector hadCand_pion_pion_tlv    = dauCand1_pion_tlv   + dauCand2_pion_tlv;
      TLorentzVector hadCand_pion_proton_tlv  = dauCand1_pion_tlv   + dauCand2_proton_tlv;
      TLorentzVector hadCand_proton_pion_tlv  = dauCand1_proton_tlv + dauCand2_pion_tlv;
      
      if ((fabs(hadCand_pion_pion_tlv.M() - ksMass_)/ksMass_) < 0.3) {// && (abs(hadCand1->PID * hadCand2->PID) == 211*211))  { 
        /*
        cout << "Pairing : "
             << " Mass diff "       << setw(14) << fabs(hadCand_pion_pion_tlv.M() - ksMass_)
             << " Reco Mass : "     << setw(8)  << hadCand_pion_pion_tlv.M()
             << " GEN Mass : "      << setw(8)  << ksMass_ 
             << " GEN PDG : "       << setw(5)  << ksPdgId_ 
             << " dau1 PID : "      << setw(5)  << hadCand1->PID
             << " dau2 PID : "      << setw(5)  << hadCand2->PID
             << " GEN idx : "       << setw(5)  << "-1"
             << " dau1 idx : "      << setw(5)  << i
             << " dau2 idx : "      << setw(5)  << j
             << endl;
        */
        m_recoHad.push_back({hadCand_pion_pion_tlv, dauCand1_pion_tlv, dauCand2_pion_tlv, -1, i, j, ksPdgId_, pionPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_pion_proton_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) {// && (abs(hadCand1->PID * hadCand2->PID) == 211*2212)) {
        /*
        cout << "Pairing : "
             << " Mass diff "      << setw(14) << fabs(hadCand_pion_proton_tlv.M() - lambda0Mass_)
             << " Reco Mass : "    << setw(8)  << hadCand_pion_proton_tlv.M()
             << " GEN Mass : "     << setw(8)  << lambda0Mass_ 
             << " GEN PDG : "      << setw(5)  << lambda0PdgId_ 
             << " dau1 PID : "     << setw(5)  << hadCand1->PID
             << " dau2 PID : "     << setw(5)  << hadCand2->PID
             << " GEN idx : "      << setw(5)  << "-1"
             << " dau1 idx : "     << setw(5)  << i
             << " dau2 idx : "     << setw(5)  << j
             << endl;
        */
        m_recoHad.push_back({hadCand_pion_proton_tlv, dauCand1_pion_tlv, dauCand2_proton_tlv, -1, i, j, lambda0PdgId_, pionPdgId_, (-1)*protonPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_proton_pion_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) {// && (abs(hadCand1->PID * hadCand2->PID) == 211*2212)) {
        /*
        cout << "Pairing : "
             << " Mass diff "      << setw(14) << fabs(hadCand_proton_pion_tlv.M() - lambda0Mass_)
             << " Reco Mass : "    << setw(8)  << hadCand_proton_pion_tlv.M()
             << " GEN Mass : "     << setw(8)  << lambda0Mass_ 
             << " GEN PDG : "      << setw(5)  << lambda0PdgId_ 
             << " dau1 PID : "     << setw(5)  << hadCand1->PID
             << " dau2 PID : "     << setw(5)  << hadCand2->PID
             << " GEN idx : "      << setw(5)  << "-1"
             << " dau1 idx : "     << setw(5)  << i
             << " dau2 idx : "     << setw(5)  << j
             << endl;
        */
        m_recoHad.push_back({hadCand_proton_pion_tlv, dauCand1_proton_tlv, dauCand2_pion_tlv, -1, i, j, lambda0PdgId_, protonPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
    }
  }
}

void HadronPreselection(std::vector<RecoHad> recoHad) {
  for (auto i=0; i < recoHad.size(); ++i) {
    if (recoHad[i].tlv.Eta() > 2.5) continue;
    else recoHad[i].isSelected = true;
  }
}

void FindTruthMatchedHadron() {
  int n = 0;
  for (int i = 0; i < particles->GetEntries(); ++i){
    bool  isFromTop = false;
    bool  isFromW   = false;
    int   isFrom    = -99;
    float x         = -99;
    auto p = (const GenParticle*) particles->At(i);
    if (p->PID != 310 && p->PID != 3122) continue;
    if (p->D1  == -1  || p->D2  == -1) continue;
    auto d1 = (const GenParticle*) particles->At(p->D1);
    auto d2 = (const GenParticle*) particles->At(p->D2);
    if ( p->PID == 310  && (abs(d1->PID) != 211 || abs(d2->PID) != 211) ) continue;
    if ( p->PID == 3122 ) {
      if      ( (abs(d1->PID) != 211  && abs(d2->PID) == 2212) || (abs(d1->PID) == 211  && abs(d2->PID) != 2212) ) continue;
      else if ( (abs(d1->PID) != 2212 && abs(d2->PID) == 211 ) || (abs(d1->PID) == 2212 && abs(d2->PID) != 211 ) ) continue;
    }

    auto motherList = getMlist(particles, p);
    for (auto j=0; j < motherList.size()-1; ++j) {
      //if ( j == 0) {
      //cout << "[ " << setw(3) << j   << " / " << setw(3) << motherList.size() 
      //     << " ] th starting point (Index)  ===> "   << setw(5) << p->PID                << " ( " << setw(4) << i                   << " ) " 
      //     << " dau1 (Index): "                       << setw(5) << d1->PID               << " ( " << setw(4) << p->D1               << " ) " 
      //     << " dau2 (Index): "                       << setw(5) << d2->PID               << " ( " << setw(4) << p->D2               << " ) " 
      //     << endl;
      //}
      //cout << "[ " << setw(3) << j+1 << " / " << setw(3) << motherList.size() 
      //     << " ] th mother particle (Index) ===> "   << setw(5) << motherList[j]->PID    << " ( " << setw(4) << motherList[j+1]->D1 << " ) "
      //     << " status : "                            << setw(5) << motherList[j]->Status
      //     << endl;
      if ( abs(motherList[j]->PID) == 24) { 
        isFromW = true;
        isFrom  = motherList[j-1]->PID;
        x       = p->PT/motherList[j-1]->PT;
        if ( abs(motherList[j-1]->PID) == 24) { isFrom = motherList[j-2]->PID; x = p->PT/motherList[j-2]->PT; } 
        if ( abs(motherList[j-2]->PID) == 24) { isFrom = motherList[j-3]->PID; x = p->PT/motherList[j-3]->PT; } 
      }
      if ( abs(motherList[j]->PID) == 6 && motherList[j]->Status == 62 ) {
        isFromTop = true; 
        if (isFromW == false) { isFrom = motherList[j-1]->PID; x = p->PT/motherList[j-1]->PT; }
        break;
      }
    }
    m_genHadron[i] = {p->P4(), d1->P4(), d2->P4(), x, i, p->PID, p->D1, d1->PID, p->D2, d2->PID, isFrom, isFromTop, isFromW};
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

    if (method == "BDT") {
      cout << "BDT ==> Work in progress" << endl;
      // auto hadBDTScore = m_hadReader->EvaluateMVA("KS_BDT");
      // if (hadBDTScore > maxBDTScore) {
      //   maxBDTScore = hadBDTScore;
      //   hIdx = i;
      // }
    }
    if (method == "pT" || method == "pt") {
      auto hadronPt = had_tlv.Pt();
      if (hadronPt > highestPt) {
        highestPt = hadronPt;
        hIdx = i;
      }
    }
    if (method == "dot") {
      auto dotProduct = had_tlv.Dot(jet_tlv);
      if (dotProduct < lowestDot) {
        lowestDot = dotProduct;
        hIdx = i;
      }
    }
    if (method == "dr" || method == "dR") {
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
    auto jet_tlv = jet->P4();

    if (m_selectedJet[i] != 0) { 
      b_isSelectedJet = true;
    }
    for ( auto j = 0; j < m_matchedJet.size(); ++j ) {
      if (m_matchedJet[j].idx == i) {
        b_isFrom        = m_matchedJet[j].gen_pdgid;
        b_hasHighestPt  = m_matchedJet[j].hasHighestPt;
        b_hasClosestLep = m_matchedJet[j].hasClosestLep;
      }
      //if (m_genQuark.size() > 1) cout << " Gen Quark : " << setw(4) << m_genQuark[0]->PID << " | " << setw(4) << m_genQuark[1]->PID << " isFrom ====> : " << setw(4) << m_matchedJet[j].gen_pdgid << " /// " << setw(4) << m_matchedJet[j].pdgid << endl;
      //else cout << "Less than 2 gen quark" << endl;
    }
    auto hIdx = FindMatchedHadron(jet_tlv, matchingMethod);
    SetJetValues(jet);
    if (hIdx != -1) {
      SetHadValues(jettr, hIdx, jet_tlv);
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
  if (m_genHadron.size() != 0) {
    auto dau1        = (Track*) tracks->At(m_recoHad[i].dau1_idx);
    auto dau2        = (Track*) tracks->At(m_recoHad[i].dau2_idx);
    auto genDau1     = (GenParticle*) dau1->Particle.GetObject();
    auto genDau2     = (GenParticle*) dau2->Particle.GetObject();
    //cout << " genHad[genDau1->M1].idx : "   << setw(5) << m_genHadron[genDau1->M1].idx 
    //     << " genDau1->M1 : "               << setw(5) << genDau1->M1
    //     << " genHad[genDau2->M1].idx : "   << setw(5) << m_genHadron[genDau2->M1].idx 
    //     << " genDau2->M1 : "               << setw(5) << genDau2->M1
    //     << " genHad[genDau1->M1].pdgid : " << setw(5) << m_genHadron[genDau1->M1].pdgid
    //     << " genDau1->pdgid : "            << setw(5) << genDau1->PID
    //     << " genHad[genDau2->M1].pdgid : " << setw(5) << m_genHadron[genDau2->M1].pdgid
    //     << " genDau2->pdgid : "            << setw(5) << genDau2->PID
    //     << endl;
    if ((genDau1->M1 == genDau2->M1) && m_genHadron[genDau1->M1].pdgid != 0) {
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
      b_genDau1_pdgId     = m_genHadron[genDau1->M1].dau1_pdgid;
      b_genDau1_pt        = m_genHadron[genDau1->M1].dau1_tlv.Pt();
      b_genDau1_eta       = m_genHadron[genDau1->M1].dau1_tlv.Eta();
      b_genDau1_phi       = m_genHadron[genDau1->M1].dau1_tlv.Phi();
      b_genDau1_mass      = m_genHadron[genDau1->M1].dau1_tlv.M();
      b_genDau2_pdgId     = m_genHadron[genDau1->M1].dau2_pdgid;
      b_genDau2_pt        = m_genHadron[genDau1->M1].dau2_tlv.Pt();
      b_genDau2_eta       = m_genHadron[genDau1->M1].dau2_tlv.Eta();
      b_genDau2_phi       = m_genHadron[genDau1->M1].dau2_tlv.Phi();
      b_genDau2_mass      = m_genHadron[genDau1->M1].dau2_tlv.M();
    } else if ( (m_genHadron[genDau1->M1].idx == genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) || (m_genHadron[genDau1->M1].idx == genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) ) {
      b_had_nMatched = 1;
    } else if ( (m_genHadron[genDau1->M1].idx != genDau1->M1 && m_genHadron[genDau2->M1].idx != genDau2->M1) ) {
      b_had_nMatched = 0;
    }
  }
}


std::vector<float> collectHadron(std::vector<GenParticle> hadInJet, bool motherCheck){
}
