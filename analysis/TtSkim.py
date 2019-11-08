# files -> time (local)
# nfil = [1,   2,   3,    10,   20,   ]
# time = [27, 28.1, 30.9, 83.1, 170,  ]

# files(npart) -> time (spark)
# nfil = [  64(32), 128(64), 982(172) ]
# time = [ 305,     373    , 5754]

import os
import Tools
import ROOT as r
from glob import glob
import TopSelections
r.ROOT.EnableImplicitMT()


def run_skim(fglob, output, spark=False, npartitions=172):
    fns = glob(fglob)
#    print(len(fns))
    print fglob, "---", output
    fs = Tools.toVector('string', fns[:5])
    # print(len(fns))
    powheg = None
    powheg = TopSelections.TopDataFrame("Events", fs, spark, npartitions, Tools.ini)

    # Based on topEventSelectionDL.cc. Currently mumu only
    powheg = powheg.Define("GenWeight", "copysign(1.,genWeight)").SetWeight('GenWeight')
    wcount = powheg.Sum("GenWeight") # Need the weighted sum for normalizing
    dimu = powheg.pvFilter()\
                 .Filter(powheg.muonTriggers, "Trigger")\
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
#                   .Define('LepWeight', 'muonSF_.getScaleFactor({Lep0.Eta(), Lep0.Pt()}, 0)*muonSF_.getScaleFactor({Lep1.Eta(), Lep1.Pt()}, 0)')
    hdim_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "DiLep_Mass")
    hmupt_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "MuonPt")

    os.system('mkdir -p skim2016')
    cutflow.Snapshot('Events', output,
                     Tools.toVector('string', [
                         'Lep0', 'Lep1', 'MET_pt', 'MET_phi',
                         'JetPt', 'JetEta', 'JetPhi', 'JetM', 'BTag',
                         'HadJet', 'HadPt', 'HadEta', 'HadPhi', 'HadMass', 'HadX', 'HadJetDR',
#                         'LepWeight', 'PUWeight', 'genWeight',
                     ]))

    cutflow.Report()
    cutflow.AppendCounts(output)

# print([c for c in cutflow.GetColumnNames() if not c.startswith('HLT')])

# # hdim_precut.Draw()
# print "Weighted entries:", wcount.GetValue()
# report = cutflow.Report()
# report.Print()

# # add cuts histo to the output
# f = r.TFile('skim2016/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', 'update')
# cuts = [cut for cut in report.GetValue()]
# hcuts = r.TH1F("cuts", "cuts", 1+len(cuts), -0.5, len(cuts)+0.5)
# hcuts.GetXaxis().SetBinLabel(1,'weighted')
# hcuts.SetBinContent(1,wcount.GetValue())
# [hcuts.GetXaxis().SetBinLabel(i+2,cut.GetName()) for i, cut in enumerate(cuts)]
# [hcuts.SetBinContent(i+2,cut.GetPass()) for i, cut in enumerate(cuts)]
# hcuts.Write()
# f.Close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run the (for now mumu only) TTbar skim')
    parser.add_argument('--fglob', '-f', default="/hdfs/store/group/nanoAOD/run2_2016v5/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_115926/**/*root", help='File glob for input files to run')
    parser.add_argument('--output', '-o', default='skim2016/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', help='filename for output root file')
    parser.add_argument('--spark', '-s', default=False, action='store_true', help='Run on a spark cluster instead of locally')
    parser.add_argument('--npartitions', '-n', type=int, default=32, help='Number of data partitions (for spark)')
    args = parser.parse_args()
    run_skim(args.fglob, args.output, args.spark, args.npartitions)
