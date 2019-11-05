# import PyRDF
# PyRDF.use('spark', {'npartitions':64})
#PyRDF.use('local')

import os
import Tools
import ROOT as r
from glob import glob
import TopSelections
r.ROOT.EnableImplicitMT()

fns = glob("/hdfs/store/group/nanoAOD/run2_2016v5/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_115926/**/*root")
fs = Tools.toVector('string', fns[:1])
# print(len(fns))
powheg = r.RDataFrame("Events", fs)

# Based on topEventSelectionDL.cc. Currently mumu only
powheg = powheg.Define("GenWeight", "copysign(1.,genWeight)")
wcount = powheg.Sum("GenWeight") # Need the weighted sum for normalizing
dimu = powheg.pvFilter()\
             .Filter(TopSelections.muonTriggers, "Trigger")\
             .Define("PUWeight", "PileupGetWeight(Pileup_nTrueInt)")\
             .electronSelection()\
             .muonSelection()\
             .Filter("nElectrons==0", "NoElectron")\
             .Filter("nMuons==2 && (MuonCharge[0]!=MuonCharge[1])", "DiMuon")\
             .Define("Lep0", "ROOT::Math::PtEtaPhiMVector(MuonPt[0],MuonEta[0],MuonPhi[0],MuonM[0])")\
             .Define("Lep1", "ROOT::Math::PtEtaPhiMVector(MuonPt[1],MuonEta[1],MuonPhi[1],MuonM[1])")\
             .Define("DiLep", "Lep0+Lep1")\
             .Define("DiLep_Mass", "DiLep.M()")
cutflow =  dimu.Filter("DiLep.M()>20.","DiLep M")\
               .Filter("MET_pt>40.","MET")\
               .jetSelection()\
               .Filter("nJets>=2", "2JetCut")\
               .hadronSelection()

hdim_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "DiLep_Mass")
hmupt_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "MuonPt")

os.system('mkdir -p skim2016')
cutflow.Snapshot('Events', 'skim2016/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root',
                 Tools.toVector('string', [
                     'Lep0', 'Lep1', 'MET_pt', 'MET_phi',
                     'JetPt', 'JetEta', 'JetPhi', 'JetM', 'BTag',
                     'HadJet', 'HadPt', 'HadEta', 'HadPhi', 'HadMass', 'HadX', 'HadJetDR',
                     'PUWeight', 'genWeight'
                 ]))

# print([c for c in cutflow.GetColumnNames() if not c.startswith('HLT')])

# hdim_precut.Draw()
print "Weighted entries:", wcount.GetValue()
report = cutflow.Report()
report.Print()

# add cuts histo to the output
f = r.TFile('skim2016/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', 'update')
cuts = [cut for cut in report.GetValue()]
hcuts = r.TH1F("cuts", "cuts", 1+len(cuts), -0.5, len(cuts)+0.5)
hcuts.GetXaxis().SetBinLabel(1,'weighted')
hcuts.SetBinContent(1,wcount.GetValue())
[hcuts.GetXaxis().SetBinLabel(i+2,cut.GetName()) for i, cut in enumerate(cuts)]
[hcuts.SetBinContent(i+2,cut.GetPass()) for i, cut in enumerate(cuts)]
hcuts.Write()
f.Close()
