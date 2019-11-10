from NanoDataFrame import NanoDataFrame
import ROOT as r
import PyRDF


class TopDataFrame(NanoDataFrame):

    """Filters and selections for running a basic Ttbar analysis.

Currently, only supports the mumu dilepton channel.
    """

    mumuTriggers = "HLT_IsoTkMu24 || HLT_IsoMu24 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"

    def mumuTriggerFilter(self):
        return self.Filter(self.mumuTriggers, "MuMu Trigger")

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

    #define TOPSELECTIONS__

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

    #define HADRONSELECTIONS__

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
    def __init__(self, tree=None, files=None, spark=False, npartitions=172, weight=None):
        super(TopDataFrame, self).__init__(tree, files, spark, npartitions, weight)
