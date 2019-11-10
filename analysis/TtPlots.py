from glob import glob
import Tools
import ROOT as r
r.ROOT.EnableImplicitMT()
r.gStyle.SetOptStat(0)

lumi = 35.9
mc = [
    ("DYJets", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8", 6225420., r.kMagenta),
    ("Single antitop", "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1", 38060., r.kBlue),
    ("Single top", "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1", 38090., r.kCyan),
    ("TTbar", "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8", 831760., r.kRed),
    ("WW", "WW_TuneCUETP8M1_13TeV-pythia8", 75800., r.kGreen),
    ("WZ", "WZ_TuneCUETP8M1_13TeV-pythia8", 23430., r.kYellow+3),
    ("ZZ", "ZZ_TuneCUETP8M1_13TeV-pythia8", 10160., r.kOrange),
]

cvs = r.TCanvas()

mcCounts = [r.TFile('skim2016/{0}/{0}.root'.format(d)) for _, d, _, _ in mc]
mcDF = [r.RDataFrame('Events', Tools.toVector('string', ['/hdfs/store/user/iawatson/TT/skim2016/{0}/{0}_*.root'.format(d)])) for _, d, _, _ in mc]
dt = r.RDataFrame('Events', Tools.toVector('string', glob('/hdfs/store/user/iawatson/TT/skim2016/DoubleMuon**/*_*.root')))

mcDF = [h.Define('Weight', '{}*GenWeight*PUWeight*LepWeight'.format(float(xsec)*lumi/float(c.Get("weights").GetBinContent(1)))) for c,h,(_,_,xsec,_) in zip(mcCounts, mcDF, mc)]
dt, mcDF = dt.Define('DiLep_Pt', 'DiLep.Pt()'), [h.Define('DiLep_Pt', 'DiLep.Pt()') for h in mcDF]

histos = [
    ("DiLep_Mass",
     [df.Histo1D(r.RDF.TH1DModel('',';M(#ell#ell) [GeV];Events / 5 GeV',100,0,500), 'DiLep_Mass', 'Weight') for df in mcDF],
     dt.Histo1D(r.RDF.TH1DModel('',';M(#ell#ell) [GeV];Events / 5 GeV',100,0,500), 'DiLep_Mass')),
    ("DiLep_Pt",
     [df.Histo1D(r.RDF.TH1DModel('',';p_{T}(#ell#ell) [GeV];Events / 5 GeV',100,0,500), 'DiLep_Pt', 'Weight') for df in mcDF],
     dt.Histo1D(r.RDF.TH1DModel('',';p_{T}(#ell#ell) [GeV];Events / 5 GeV',100,0,500), 'DiLep_Pt')),
    ("MuonPt",
     [df.Histo1D(r.RDF.TH1DModel('',';p_{T}(#ell) [GeV];Events / 5 GeV',100,0,500), 'MuonPt', 'Weight') for df in mcDF],
     dt.Histo1D(r.RDF.TH1DModel('',';p_{T}(#ell) [GeV];Events / 5 GeV',100,0,500), 'MuonPt')),
]   

for name, hmc, hdt in histos:
    leg = r.TLegend(0.65,0.45,0.83,0.87, "", "brNDC")
    leg.SetBorderSize(0)
    leg.SetTextSize(0.046)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)

    hdt.SetMarkerStyle(8)
    hdt.SetMarkerSize(1)
    hdt.SetMarkerColor(r.kBlack)
    #hdt.Draw('p')
    leg.AddEntry(hdt.GetPtr(), "#mu#mu", "")
    leg.AddEntry(hdt.GetPtr(), "Data", "p")

    stk = r.THStack()
    # hmc = [h.GetValue() for h in hmc]
    [h.SetFillColor(col) for (_,_,_,col), h in zip(mc, hmc)]
    [h.SetLineColor(col) for (_,_,_,col), h in zip(mc, hmc)]
    [h.SetMarkerColor(col) for (_,_,_,col), h in zip(mc, hmc)]
    [h.SetMarkerStyle(21) for (_,_,_,col), h in zip(mc, hmc)]
    [h.SetMarkerSize(3) for (_,_,_,col), h in zip(mc, hmc)]
    [leg.AddEntry(h.GetPtr(), n, "p") for (n,_,_,_), h in zip(mc, hmc)]
    [stk.Add(h.GetPtr()) for h in hmc]

    mchist = stk.GetHists().At(0).Clone()
    mchist.Reset()
    for i in range(stk.GetHists().GetSize()):
        mchist.Add(stk.GetHists().At(i))

    rat = r.TRatioPlot(hdt.GetPtr(), mchist)
    rat.SetH1DrawOpt("p")
    rat.SetH2DrawOpt("hist")
    rat.Draw()

    scale = 1.1;
    stk.SetMaximum(scale*max(stk.GetMaximum(), d.GetMaximum()));
    d.SetMaximum(scale*max(stk.GetMaximum(), d.GetMaximum()));
    cvs.Update()

    up = rat.GetUpperPad()
    up.cd()
    stk.Draw("same hist")
    hdt.Draw("same p")
    leg.Draw()

    rat_range=0.6
    rat.GetLowerRefYaxis().SetTitleSize(0.045)
    rat.GetLowerRefYaxis().SetTitleOffset(1.3)
    rat.GetLowerRefXaxis().SetTitleSize(0.045)
    rat.GetLowerRefYaxis().SetLabelSize(0.035)
    rat.GetLowerRefXaxis().SetLabelSize(0.035)
    rat.GetLowerRefXaxis().CenterTitle(True)
    rat.GetLowerRefYaxis().SetRangeUser(1-rat_range, 1+rat_range)
    rat.GetLowerRefYaxis().SetNdivisions(103)
    rat.GetLowerRefXaxis().SetTitleOffset(0.)
    rat.GetLowerPad().Update()
    rat.GetLowerRefYaxis().SetRangeUser(1-rat_range, 1+rat_range)
    rat.GetUpperRefYaxis().SetTitleFont(42)
    rat.GetUpperRefYaxis().SetLabelSize(0.035)
    rat.GetUpperRefYaxis().SetTitleSize(0.045)
    rat.GetUpperRefYaxis().SetTitleOffset(1.3)
    rat.SetRightMargin(.03)
    rat.SetLeftMargin(.12)

    up.SetLogy()
    cvs.Update()
    cvs.Print("{}.png".format(name))
