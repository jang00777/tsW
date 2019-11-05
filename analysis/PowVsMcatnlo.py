# import PyRDF
# PyRDF.use('spark', {'npartitions':64})
#PyRDF.use('local')

import ROOT as r
from glob import glob
import TopSelections
r.ROOT.EnableImplicitMT()

def toVector(typ, l):
    v = r.vector(typ)()
    for i in l: v.push_back(i)
    return v

fns = glob("/hdfs/store/group/nanoAOD/run2_2016v5/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180607_115926/**/*root")
fs = toVector('string', fns[:2])
# print(len(fns))
powheg = r.RDataFrame("Events", fs)

# Based on topEventSelectionDL.cc. Currently mumu only
dimu = powheg.pvFilter()\
             .Filter(TopSelections.muonTriggers, "Trigger")\
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

cutflow.Snapshot('Events', 'test.root', toVector('string', ['Lep0', 'Lep1', 'MET_pt', 'MET_phi',
                                                            'JetPt', 'JetEta', 'JetPhi', 'JetM', 'BTag',
                                                            'HadJet', 'HadPt', 'HadEta', 'HadPhi', 'HadMass', 'HadX', 'HadJetDR',
]))

# print([c for c in cutflow.GetColumnNames() if not c.startswith('HLT')])

# hdim_precut.Draw()
report = cutflow.Report()
report.Print()

r.gInterpreter.Declare('using namespace ROOT::VecOps; float Maximum(const RVec<float> & a) {float ret=-9.9e99; for(auto i : a){if(i>ret) ret=i;} return ret;}')

# disp=r.vector('string')()
# disp.push_back('Muon_pt')
# # dimu.Display(disp).Print()
