import os
import Tools
import ROOT as r
from glob import glob
import TopSelections
import atexit
r.ROOT.EnableImplicitMT()

# files -> time (local)
# nfil = [10, 20 ]
# time = [57, 119]

# files(npart) -> time (spark)
# nfil = [ 20(32) ]
# time = [ 93 ]

def run_skim(fglob, output, spark=False, npartitions=172, merge=False, hadron=False, isMc=None):
    fns = glob(fglob)
    fs = Tools.toVector('string', fns)
    powheg = TopSelections.TopDataFrame("Events", fs, spark, npartitions)
    if isMc is None:
        # Spark doesn't seem to have GetColumnNames...
        isMc = True if spark else ("GenPart_pt" in powheg.GetColumnNames())
    print fglob, " [ {} ] ---".format(len(fns)), output, isMc

    # Based on topEventSelectionDL.cc. Currently mumu only
    if isMc:
        powheg = powheg.Define("GenWeight", "copysign(1.,genWeight)")\
                       .SetWeight('GenWeight')\
                       .Define("PUWeight", "PileupGetWeight(Pileup_nTrueInt)")
    dimu = powheg.mumuTriggerFilter()\
                 .pvFilter()\
                 .electronSelection()\
                 .muonSelection()\
                 .Filter("nElectrons==0", "NoElectron")\
                 .Filter("nMuons==2 && (MuonCharge[0]!=MuonCharge[1])", "DiMuon")\
                 .Define("Lep0", "PtEtaPhiMVector(MuonPt[0],MuonEta[0],MuonPhi[0],MuonM[0])")\
                 .Define("Lep1", "PtEtaPhiMVector(MuonPt[1],MuonEta[1],MuonPhi[1],MuonM[1])")\
                 .Define("DiLep", "Lep0+Lep1")\
                 .Define("DiLep_Mass", "DiLep.M()")
    cutflow =  dimu.Filter("DiLep.M()>20.","DiLep M")\
                   .Filter("MET_pt>40.","MET")\
                   .jetSelection()\
                   .Filter("nJets>=2", "2JetCut")\
                   .Define('LepWeight', 'muonSF_.getScaleFactor({Lep0.Eta(), Lep0.Pt()}, 0)*muonSF_.getScaleFactor({Lep1.Eta(), Lep1.Pt()}, 0)')

    # hdim_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "DiLep_Mass")
    # hmupt_precut = dimu.Histo1D(r.RDF.TH1DModel("", "", 100, 0, 300), "MuonPt")

    if hadron:
        cutflow = cutflow.hadronSelection()

    # Note that currently, the PyRDF backend doesn't support selecting branches to save :-(
    columns = ['Lep0', 'Lep1', 'MET_pt', 'MET_phi', 'JetPt', 'JetEta',
               'JetPhi', 'JetM', 'BTag',]
    if isMc:
        columns += ['PUWeight', 'GenWeight','LepWeight']
    if hadron:
        columns += ['HadJet', 'HadPt', 'HadEta', 'HadPhi', 'HadMass',
                    'HadX', 'HadJetDR']
    ch = cutflow.Snapshot('Events', output, Tools.toVector('string', columns)) 
    
    if spark and isinstance(ch.args[0], r.TChain) and merge:
        print "Merging TChain into output file"
        fs = ch.get_inputfiles()
        ch.args[0].Merge(output)
        print "Merging finished"
        def clean(filelist):
            print "Cleaning up"
            for f in filelist:
                os.system("rm {}".format(f))
            print "Finished clean"
        print dir(ch)
        ch.__del__()
        atexit.register(clean, [f.GetTitle() for f in ch.args[0].GetListOfFiles()])
    cutflow.Report()
    cutflow.AppendCounts(output)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run the (for now mumu only) TTbar skim')
    parser.add_argument('--fglob', '-f', default="/hdfs/store/group/nanoAOD/run2_2016v5/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_115926/**/*root", help='File glob for input files to run')
    parser.add_argument('--output', '-o', default='skim2016/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root', help='filename for output root file')
    parser.add_argument('--hadron', '-a', default=False, action='store_true', help='Add hadron information (requires the Hadron version of the NanoAOD skim)')
    parser.add_argument('--spark', '-s', default=False, action='store_true', help='Run on a spark cluster instead of locally')
    parser.add_argument('--npartitions', '-n', type=int, default=32, help='Number of data partitions (for spark)')
    parser.add_argument('--merge', '-m', default=False, action='store_true', help="Merge root file at the end of spark run (this will give an error but seems to work). If this option isn't used, the spark output will be one file per partition")
    args = parser.parse_args()
    outdir = os.path.split(os.path.realpath(args.output))[0]
    os.system('mkdir -p {}'.format(outdir))
    run_skim(args.fglob, args.output, args.spark, args.npartitions, args.merge, args.hadron)
