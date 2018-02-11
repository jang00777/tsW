#include <stdio.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
#include "wj_analysis.h"

int main(int argc, char* argv[])
{
  auto inf = std::string{argv[1]};
  auto outf = std::string{argv[2]};

  // check cpu time (start)
  std::clock_t c_start = std::clock();

  //read input file
  auto tfiles = TFile::Open(inf.c_str(), "READ");
  auto trees = (TTree*) tfiles->Get("Delphes");
  trees->SetBranchStatus("*", true);
  trees->SetBranchAddress("GenJet",      &gen_jets);
  trees->SetBranchAddress("Particle",    &particles);
  trees->SetBranchAddress("Electron",    &electrons);
  trees->SetBranchAddress("Muon",        &muons);
  trees->SetBranchAddress("Jet",         &jets);
  trees->SetBranchAddress("MissingET",   &missingET);
  trees->SetBranchAddress("Track",       &tracks);

  //make out file
  auto out = TFile::Open(outf.c_str(), "RECREATE");
  auto outtr = new TTree("tsW", "tsW");

  defBranchFucns(outtr);

  TString cutflow_title = "cutflow" + inf;
  TH1F * cutflow = new TH1F("cutflow", cutflow_title, 7,-1,6); // -1 : all events, 0 : events after lepton selection, 1~ : events after step

  //Event Loop Start!
  for (size_t iev = 0; iev < trees->GetEntries(); ++iev){
    if (iev%1000 == 0 ) std::cout << "event check    iev    ----> " << iev << std::endl;
    trees->GetEntry(iev);
    initValues();

    recoParticle(cutflow);
    genParticle();

    /*
    // pion track test
    for (unsigned i = 0; i < tracks->GetEntries(); ++ i){
      auto track = (Track*) tracks->At(i);
      if (abs(track->PID) == 321){
        auto genTrack = (GenParticle*) track->Particle.GetObject();
        auto genTrack_mom = (const GenParticle*) particles->At(genTrack->M1);
        std::cout << genTrack_mom->PID << std::endl;
        float pionR = sqrt(pow(track->X,2)+pow(track->Y,2)+pow(track->Z,2));
        float pionRd = sqrt(pow(track->Xd,2)+pow(track->Yd,2)+pow(track->Zd,2));
        float pionROuter = sqrt(pow(track->XOuter,2)+pow(track->YOuter,2)+pow(track->ZOuter,2));
      }
    }
    */

    outtr->Fill();  
  }
  tfiles->Close();
  outtr->Write();
  cutflow->Write();

  out->Close();

  //check cpu time (end)
  std::clock_t c_end = std::clock();
  long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  std::cout << "CPU time used(sec): " << time_elapsed_ms/1000 << " sec\n";
  return 0;
}

void initValues(){
    dilepton_mass = -99;
    dilepton_ch = 0; 

    gen_step0 = false; gen_step1 = false;

    jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_energy.clear();
    jet_pid.clear();

    kshortsInjet_pt.clear(); kshortsInjet_eta.clear(); kshortsInjet_phi.clear(); kshortsInjet_energy.clear(); kshortsInjet_R.clear(); kshortsInjet_outR.clear();
    lambdasInjet_pt.clear(); lambdasInjet_eta.clear(); lambdasInjet_phi.clear(); lambdasInjet_energy.clear(); lambdasInjet_R.clear(); lambdasInjet_outR.clear();
    leptonsInjet_pt.clear(); leptonsInjet_eta.clear(); leptonsInjet_phi.clear(); leptonsInjet_energy.clear(); leptonsInjet_R.clear();
    nkshortsInjet.clear(); nlambdasInjet.clear(); nleptonsInjet.clear();

    jet1_diHadron_mass.clear(); jet2_diHadron_mass.clear();
}

void recoParticle(TH1F * cutflow){
    cutflow->Fill(-1);
    // object selection
    std::vector<struct Lepton> recolep;
    for (unsigned i = 0; i < muons->GetEntries(); ++i){
      auto mu = (Muon*) muons->At(i);
      if (abs(mu->Eta) > 2.4) continue;
      if (mu->PT < 20.) continue;
      recolep.push_back(toLepton(mu));
    }
    for (unsigned j = 0; j < electrons->GetEntries(); ++j){
      auto elec = (Electron*) electrons->At(j);
      if (abs(elec->Eta) > 2.4) continue;
      if (elec->PT < 20.) continue;
      recolep.push_back(toLepton(elec));
    }

    // pick highest pt lep pair
    if (recolep.size() < 2 ) return;

    sort(recolep.begin(), recolep.end(), [](struct Lepton a, struct Lepton b){return a.tlv.Pt() > b.tlv.Pt();});
    recolep.erase(recolep.begin()+2,recolep.end());

    cutflow->Fill(0);

    auto dilepton = recolep[0].tlv + recolep[1].tlv;
    dilepton_ch = abs(recolep[0].pdgid) + abs(recolep[1].pdgid); // 22 -> ee , 24 -> emu , 26 -> mumu
    dilepton_mass = dilepton.M();

    // step 1
    if(dilepton.M() < 20. || recolep[0].charge * recolep[1].charge > 0 ) return;
    cutflow->Fill(1);

    // step2
    if ( (dilepton_ch != 24) && ( dilepton.M() > 76. && dilepton.M() < 106.) ) return;
    cutflow->Fill(2);

    // step3
    auto met = ((MissingET*) missingET->At(0))->MET;
    if ( (dilepton_ch != 24) && (met < 40.) ) return;
    cutflow->Fill(3);

    // step4  
    std::vector<Jet*> selectedJets;
    for ( unsigned k = 0; k < jets->GetEntries(); ++k){
      auto jet = (Jet*) jets->At(k);
      if (jet->PT < 30.) continue;
      if (abs(jet->Eta) > 2.4) continue;
      bool hasOverlap = false;
      for (auto lep : recolep) { if ((jet->P4()).DeltaR(lep.tlv) < 0.4) hasOverlap = true; }
      if (hasOverlap) continue;
      selectedJets.push_back(jet);
    }
    if (selectedJets.size() < 2) return;
    cutflow->Fill(4);

    // step5
    std::vector<Jet*> selectedBJets;
    for (auto jet : selectedJets){
      if (jet->BTag) selectedBJets.push_back(jet);
    }
    if (selectedBJets.size() < 1) return;
    cutflow->Fill(5);

}

void genParticle(){
    std::vector<const GenParticle*> genTops;
    std::vector<const GenParticle*> genJets;
    std::vector<struct Lepton> genLeps;
    std::vector<std::vector<float> > jet_diHadron_mass;
    for (unsigned i = 0; i < particles->GetEntries(); ++ i){
      auto p = (const GenParticle*) particles->At(i);
      if (p->Status > 30) continue;
      if (abs(p->PID) != 6) continue;
      auto top = getLast(particles,p);
      genTops.push_back(top);

      auto jet = (const GenParticle*) particles->At(top->D2);
      auto Wboson = (const GenParticle*) particles->At(top->D1);
      auto lastBoson = getLast(particles, Wboson);
      genJets.push_back(jet);

      auto lep_tmp = (const GenParticle*) particles->At(lastBoson->D1);
      auto neu = (const GenParticle*) particles->At(lastBoson->D2);
      struct Lepton lep = toLepton(lep_tmp);
      genLeps.push_back(lep);
      jet_pid.push_back(jet->PID);
      auto jet_tlv = jet->P4();
      jet_pt.push_back(jet_tlv.Pt());
      jet_eta.push_back(jet_tlv.Eta());
      jet_phi.push_back(jet_tlv.Phi());
      jet_energy.push_back(jet_tlv.E());

      // collect jet constitues
      std::vector<GenParticle> jetConstitues;
      std::vector<int> pdgList = {211, 321, 2212, 310, 3122, 11, 13};
      for (unsigned i2 = 0; i2 < particles->GetEntries(); ++ i2){
        auto p2 = (const GenParticle*) particles->At(i2);

        if ( !(std::find(pdgList.begin(), pdgList.end(), abs(p2->PID)) != pdgList.end()) ) continue;
	    if ( jet_tlv.DeltaR(p2->P4()) > 0.5) continue;
        bool fromJet = false;
        std::vector<const GenParticle*> mlist = getMlist(particles, p2);
        for (auto m : mlist){ if (m == jet) fromJet = true; }
        if (!fromJet) continue;
        
        jetConstitues.push_back(*p2);
      }
      /*
      // collect jet constitues
      std::vector<GenParticle> jetConstitues;
      std::vector<int> pdgList = {211, 321, 2212, 310, 3122, 11, 13};
      for ( unsigned k = 0; k < gen_jets->GetEntries(); ++k){
        auto genjet = (Jet*) gen_jets->At(k);
        auto test = (genjet->Constituents).At(0);
        std::cout << test->PID << std::endl;
        for (c : jet.Constituents){
          if ( !(std::find(pdgList.begin(), pdgList.end(), abs(p2->PID)) != pdgList.end()) ) continue;
          jetConstitues.push_back(*p2);
        }
      }
      */

      std::vector<GenParticle> kshortsInjet;
      std::vector<GenParticle> lambdasInjet;
      std::vector<GenParticle> hadronsInjet;
      std::vector<GenParticle> leptonsInjet;
      for (auto c : jetConstitues){
        int absPid = abs(c.PID);
        if (absPid == 310)  kshortsInjet.push_back(c);
        if (absPid == 3122) lambdasInjet.push_back(c);
        if (absPid == 2212 || absPid == 311 || absPid == 211)  hadronsInjet.push_back(c);
        if (absPid == 11 || absPid == 13)  leptonsInjet.push_back(c);
      }

      // save highest pT Ks only
      nkshortsInjet.push_back(kshortsInjet.size());
      if (kshortsInjet.size() >0){
        sort(kshortsInjet.begin(), kshortsInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});
        kshortsInjet_pt.push_back(kshortsInjet[0].PT);
        kshortsInjet_eta.push_back(kshortsInjet[0].Eta);
        kshortsInjet_phi.push_back(kshortsInjet[0].Phi);
        kshortsInjet_energy.push_back(kshortsInjet[0].E);
        auto R = sqrt(pow(kshortsInjet[0].X,2)+pow(kshortsInjet[0].Y,2)+pow(kshortsInjet[0].Z,2));
        kshortsInjet_R.push_back(R);

        auto kshortDau = (const GenParticle*) particles->At(kshortsInjet[0].D1);
        auto outR = sqrt(pow(kshortDau->X,2)+pow(kshortDau->Y,2)+pow(kshortDau->Z,2));
        kshortsInjet_outR.push_back(outR);
      }
      else {
        kshortsInjet_pt.push_back(-99);
        kshortsInjet_eta.push_back(-99);
        kshortsInjet_phi.push_back(-99);
        kshortsInjet_energy.push_back(-99);
        kshortsInjet_R.push_back(-99);
        kshortsInjet_outR.push_back(-99);
      }

      // save highest pT lambda only
      nlambdasInjet.push_back(lambdasInjet.size());
      if (lambdasInjet.size() >0){
        sort(lambdasInjet.begin(), lambdasInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});
        lambdasInjet_pt.push_back(lambdasInjet[0].PT);
        lambdasInjet_eta.push_back(lambdasInjet[0].Eta);
        lambdasInjet_phi.push_back(lambdasInjet[0].Phi);
        lambdasInjet_energy.push_back(lambdasInjet[0].E);
        auto R = sqrt(pow(lambdasInjet[0].X,2)+pow(lambdasInjet[0].Y,2)+pow(lambdasInjet[0].Z,2));
        lambdasInjet_R.push_back(R);

        auto lambdaDau = (const GenParticle*) particles->At(lambdasInjet[0].D1);
        auto outR = sqrt(pow(lambdaDau->X,2)+pow(lambdaDau->Y,2)+pow(lambdaDau->Z,2));
        lambdasInjet_outR.push_back(outR);
      }
      else {
        lambdasInjet_pt.push_back(-99);
        lambdasInjet_eta.push_back(-99);
        lambdasInjet_phi.push_back(-99);
        lambdasInjet_energy.push_back(-99);
        lambdasInjet_R.push_back(-99);
        lambdasInjet_outR.push_back(-99);
      }

      // save highest pT lepton only
      nleptonsInjet.push_back(leptonsInjet.size());
      if (leptonsInjet.size() >0){
        sort(leptonsInjet.begin(), leptonsInjet.end(), [](GenParticle a, GenParticle b){return a.PT > b.PT;});
        leptonsInjet_pt.push_back(leptonsInjet[0].PT);
        leptonsInjet_eta.push_back(leptonsInjet[0].Eta);
        leptonsInjet_phi.push_back(leptonsInjet[0].Phi);
        leptonsInjet_energy.push_back(leptonsInjet[0].E);
        auto R = sqrt(pow(leptonsInjet[0].X,2)+pow(leptonsInjet[0].Y,2)+pow(leptonsInjet[0].Z,2));
        leptonsInjet_R.push_back(R);
      }
      else {
        leptonsInjet_pt.push_back(-99);
        leptonsInjet_eta.push_back(-99);
        leptonsInjet_phi.push_back(-99);
        leptonsInjet_energy.push_back(-99);
        leptonsInjet_R.push_back(-99);
      }

      std::vector<float> diHadron_mass;
      if (hadronsInjet.size() >= 2){
        //std::vector<float> diHadron_energy;
        //std::vector<bool> diHadron_sameMother;
        for (unsigned ih1 = 0; ih1 < hadronsInjet.size(); ++ih1){
          auto hadron1 = hadronsInjet[ih1].P4();
          for (unsigned ih2 = 0; ih2 < hadronsInjet.size(); ++ih2){
            if (ih1 >= ih2) continue;
            auto hadron2 = hadronsInjet[ih2].P4();

            auto diHadron = hadron1+hadron2;
            diHadron_mass.push_back(diHadron.M());
            //diHadron_energy.push_back(diHadron.E());
            //diHadron_sameMother.push_back(hadronsInjet[ih1].M1 == hadronsInjet[ih2].M1 );
          }
        }
      }
      jet_diHadron_mass.push_back(diHadron_mass);

    }

    jet1_diHadron_mass = jet_diHadron_mass[0];
    jet2_diHadron_mass = jet_diHadron_mass[1];

    // accepted leptons
    int selectedGenLep = 0;
    if (abs(genLeps[0].tlv.Eta()) < 2.4 && genLeps[0].tlv.Pt() > 20.) selectedGenLep++;
    if (abs(genLeps[1].tlv.Eta()) < 2.4 && genLeps[1].tlv.Pt() > 20.) selectedGenLep++;

    if (selectedGenLep == 2){
      gen_step0 = true;
      auto genDilep = genLeps[0].tlv+genLeps[1].tlv;
      if (genDilep.M() > 20. && genLeps[0].charge * genLeps[1].charge < 0 ) gen_step1 = true;
    }
}

std::vector<float> collectHadron(std::vector<GenParticle> hadronsInjet, bool motherCheck){
}
