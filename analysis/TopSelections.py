from NanoDataFrame import NanoDataFrame
import ROOT as r
import PyRDF


class TopDataFrame(NanoDataFrame):

    muonTriggers = "HLT_IsoTkMu24 || HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"

    def pvFilter(self):
        return self.Filter("PV_npvs!=0 && PV_ndof >=4 && fabs(PV_z) < 24.", "PV")

    def electronSelection(self):
        return self.Define("el_scEta", "Electron_deltaEtaSC + Electron_eta")\
                   .Define("ElectronSelection", "Electron_cutBased>=3 && Electron_pt>20. && abs(Electron_eta)<2.4 && (abs(el_scEta) < 1.4442 || abs(el_scEta) > 1.566)")\
                   .Define("ElectronPt", "Electron_pt[ElectronSelection]")\
                   .Define("ElectronEta", "Electron_eta[ElectronSelection]")\
                   .Define("ElectronPhi", "Electron_phi[ElectronSelection]")\
                   .Define("ElectronM", "Electron_mass[ElectronSelection]")\
                   .Define("ElectronCharge", "Electron_charge[ElectronSelection]")\
                   .Define("nElectrons", "ElectronPt.size()")

    def muonSelection(self):
        return self.Define("MuonSelection", "Muon_tightId && Muon_pt>20. && abs(Muon_eta)<2.4 && Muon_pfRelIso04_all<0.15")\
                   .Define("MuonPt", "Muon_pt[MuonSelection]")\
                   .Define("MuonEta", "Muon_eta[MuonSelection]")\
                   .Define("MuonPhi", "Muon_phi[MuonSelection]")\
                   .Define("MuonM", "Muon_mass[MuonSelection]")\
                   .Define("MuonCharge", "Muon_charge[MuonSelection]")\
                   .Define("nMuons", "MuonPt.size()")

    overlap_code ='''
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
    '''
    r.gInterpreter.Declare(overlap_code)

    def jetSelection(self):
        return self.Define("JetMuOverlap", "fillOverlap(Jet_eta, Jet_phi, MuonEta, MuonPhi)")\
                   .Define("JetElOverlap", "fillOverlap(Jet_eta, Jet_phi, ElectronEta, ElectronPhi)")\
                   .Define("JetSelection", "Jet_jetId>=1 && Jet_pt>30. && abs(Jet_eta)<2.4 && !JetMuOverlap && !JetElOverlap")\
                   .Define("BTagged", "JetSelection && Jet_btagCSVV2 > 0.8484")\
                   .Define("JetPt", "Jet_pt[JetSelection]")\
                   .Define("JetEta", "Jet_eta[JetSelection]")\
                   .Define("JetPhi", "Jet_phi[JetSelection]")\
                   .Define("JetM", "Jet_mass[JetSelection]")\
                   .Define("BJetPt", "Jet_pt[BTagged]")\
                   .Define("BTag", "BTagged[JetSelection]")\
                   .Define("nJets", "JetPt.size()")\
                   .Define("nBJets", "BJetPt.size()")

    hadron_code ='''
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
    '''

    r.gInterpreter.Declare(hadron_code)

    def hadronSelection(self):
        return self.Define("JetHadron", "jetHadron(had_eta, had_phi, JetEta, JetPhi)")\
            .Define("JetHadronOverlap", "had_pdgId == 310 && JetHadron != -1")\
            .Define("HadJet", "JetHadron[JetHadronOverlap]")\
            .Define("HadEta", "had_eta[JetHadronOverlap]")\
            .Define("HadPhi", "had_phi[JetHadronOverlap]")\
            .Define("HadPt", "had_pt[JetHadronOverlap]")\
            .Define("HadMass", "had_mass[JetHadronOverlap]")\
            .Define("HadX", "GetX(HadPt, HadJet, JetPt)")\
            .Define("HadJetDR", "GetHadJetDR(HadEta, HadPhi, HadJet, JetEta, JetPhi)")

    def __init__(self, tree=None, files=None, spark=False, npartitions=172, initialize=None, weight=None):
        if initialize is None: initialize = []
        super(TopDataFrame, self).__init__(tree, files, spark, npartitions, initialize+[self.hadron_code, self.overlap_code], weight)
