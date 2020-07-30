#!/usr/bin/env python3
# https://twiki.cern.ch/twiki/bin/view/RooStats/RooStatsExercisesMarch2015

import ROOT, time, os, sys, copy
from ROOT import RooStats, RooFit, RooWorkspace

import json

#datasets = json.load(open("%s/src/nano/nanoAOD/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))
datasets = json.load(open("/home/wjang/CMSSW_9_4_4/src/nano/nanoAOD/data/dataset/dataset.json"))
bkglist = ['WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']

# !! Background here is non-sjet

# Configuration
IntegLumi_2016    = 35.92 # fb-1 
IntegLumi_2017    = 41.86 # fb-1 
IntegLumi_2018    = 58.83 # fb-1 
IntegLumi         = IntegLumi_2016 + IntegLumi_2017 + IntegLumi_2018 # fb-1 (Run2 Total) 

nPoint = 40
poiMin = 0.
poiMax = 4.

for ntimes in [1, 2, 3, 4, 5, 6, 7, 10, 20]:
#for ntimes in [1]:
    xsec_bkgs = {}
    for d in datasets:
        name = str(d["name"])
        if name == 'TTJets_aMC':
            xsec_TTJets_aMC = d["xsec"]
        if name in bkglist:
            xsec_bkgs[name] = d["xsec"]

    TotalLumi = ntimes * IntegLumi
     
    nttbar            = ntimes * IntegLumi * xsec_TTJets_aMC * 1000
    nbkgprocess       = ntimes * IntegLumi * sum(xsec_bkgs.values()) * 1000
    
    BR_w_to_e         = 0.1071
    BR_w_to_mu        = 0.1063
    BR_w_to_tau       = 0.1138
    BR_tau_to_e       = 0.1782
    BR_tau_to_mu      = 0.1739
    BR_ttbar_to_dilep = (BR_w_to_e + BR_w_to_mu + BR_w_to_tau)**2 # ratio of dilepton (e, mu, tau)
    BR_ttbar_to_emu   = (BR_w_to_e + BR_w_to_mu + BR_w_to_tau*(BR_tau_to_e + BR_tau_to_mu))**2 # ratio of dilepton channel (e, mu)
    
    
    # from 'getConfigForDelphes.py' 
  
    ep_eventSel_ttbar = 0.198968947927
    ep_sjetReco = 0.724568906197
    ep_bjetReco = 0.706045762461
    ep_KsReco_BR_Ks_inSJet = 0.0626274375162
    ep_KsReco_BR_Ks_inBJet = 0.0261087962463
    ep_KsReco_BR_Ks_inOtherJet = 0.0449888641425

    ep_eventSel_bkgs = {'SingleTop_tW': 0.007754, # from seulgi's bkgs
                        'TT_powheg': 0.01178,
                        'SingleTbar_tW': 0.007598,
                        'WW': 0.0029722602423566153,
                        'WJets': 0.0, # delphes
                        #'WJets': 2.0494579371971637e-06, # cmssw
                        'ZZ': 0.00033473373635007936,
                        'SingleTop_tch': 6e-06,
                        'SingleTbar_tch': 8e-06,
                        'DYJets': 0.0, # delphes
                        #'DYJets': 0.0004163625516465179, # cmssw
                        #'DYJets_10to50': 3.589531159199871e-05, # cmssw
                        'DYJets_10to50' : 0.0, # no sample in seulgi's bkg
                        'WZ': 0.0007759261865037088}
     
     
    # http://pdg.lbl.gov/2018/reviews/rpp2018-rev-ckm-matrix.pdf
    Vts               = 0.04133
    Vtb               = 0.999105
    Vts2              = Vts**2
    Vtb2              = Vtb**2
    
    N_selected_ttbar  = nttbar * BR_ttbar_to_emu * ep_eventSel_ttbar
    N_selected_bkg    = IntegLumi * sum([xsec_bkgs[name] * ep_eventSel_bkgs[name] for name in bkglist]) * 1000
    N_bWbW            = N_selected_ttbar * (1 - Vts2) * (1 - Vts2) * ep_bjetReco
    N_sWbW            = N_selected_ttbar * Vts2 * (1 - Vts2) * ep_sjetReco * 2 # sWbW + bWsW
    N_sWsW            = N_selected_ttbar * Vts2 * Vts2 * ep_sjetReco * ep_sjetReco
    
    expNsig           = N_sWbW + N_sWsW
    expNbkg           = N_bWbW + N_selected_bkg

    expNsig_wKs       = expNsig * ep_KsReco_BR_Ks_inSJet
    expNbkg_wKs       = expNbkg * (ep_KsReco_BR_Ks_inBJet+ep_KsReco_BR_Ks_inOtherJet)
    expNsig_woKs      = expNsig * (1-ep_KsReco_BR_Ks_inSJet)
    expNbkg_woKs      = expNbkg * (1-ep_KsReco_BR_Ks_inBJet-ep_KsReco_BR_Ks_inOtherJet)

    if "all" in sys.argv[1] :
        expNsig_wKs       = expNsig * ep_KsReco_BR_Ks_inSJet
        expNbkg_wKs       = expNbkg * (ep_KsReco_BR_Ks_inBJet+ep_KsReco_BR_Ks_inOtherJet)

    if "tt" in sys.argv[1] : 
        expNsig_wKs       = ntimes*(209)         # expNsig_wKs from https://docs.google.com/spreadsheets/d/1WTB4gv29uJpDs_ft6oG_WONj8n9G-u_KJrzV2PdVN0o/edit#gid=810703199 
        expNbkg_wKs       = ntimes*(85 + 26812) # expNsig_wKs from https://docs.google.com/spreadsheets/d/1WTB4gv29uJpDs_ft6oG_WONj8n9G-u_KJrzV2PdVN0o/edit#gid=810703199 
        #expNsig_wKs       = expNsig * ep_KsReco_BR_Ks_inSJet
        #expNbkg_wKs       = N_bWbW * (ep_KsReco_BR_Ks_inBJet+ep_KsReco_BR_Ks_inOtherJet)

    #For sigma line
    sigmaLevels = {2:0.0455, 3:0.0027, 4:0.000063, 5:0.00000057}
     
    def setConst(fitFunc):
        argList = ROOT.RooArgList(fitFunc.getVariables())
        for i in range(argList.getSize()):
            vname = argList.at(i).GetName()
            if all(s not in vname for s in ["bdt", "nsig", "nbkg"]) :
                argList.at(i).setConstant()
    
    def getExpValue(inverter, cl, useCLs):
        inverter.SetConfidenceLevel(cl)
        inverter.UseCLs(useCLs)
        r = inverter.GetInterval()
        return r.GetExpectedUpperLimit()
    
    def drawSigmaLines(color, poiMin = 0, poiMax = 0):
        line = ROOT.TLine()
        line.SetLineColor(color)
        sigma = ROOT.TLatex()
        sigma.SetTextColor(color)
        sigma.SetTextAlign(12)
        for n in sigmaLevels:
            y = sigmaLevels[n]
            if poiMin == 0 and poiMax == 0 : 
                line.DrawLine(0, y, 2, y)
                sigma.DrawLatex(poiMax + 0.05, y, "%d #sigma"%n)
            else : 
                line.DrawLine(poiMin, y, poiMax, y)
                sigma.DrawLatex(poiMax + 0.05, y, "%d #sigma"%n)
    
    def correctErr(h, i):
        x = h.get(i)
        y = h.weight(x)
        h.set(x, y, ROOT.TMath.Sqrt(y))
    
    def getFitFunc(w, trainType):
        sigName = "sig_" + trainType
        bkgName = "bkg_" + trainType
    
        w.factory("Gaussian::{0}1(bdt, {0}_mean1[-0.03,-1, 1.], {0}_sigma1[0.06,0, 1])".format(sigName))
        w.factory("Gaussian::{0}2(bdt, {0}_mean2[-0.1, -1, 1.], {0}_sigma2[0.06, 0, 1])".format(sigName))
        w.factory("Gaussian::{0}3(bdt, {0}_mean3[-0.1, -1, 1.], {0}_sigma3[0.1, 0, 1])".format(sigName))
        fit_sig = w.factory("SUM::{0}({0}1, {0}_frac1[.33, 0, 1]*{0}2, {0}_frac2[.33, 0, 1]*{0}3)".format(sigName))
    
        #Two Gaussian
        w.factory("Gaussian::{0}1(bdt, {0}_mean1[-0.15, -1, 1.], {0}_sigma1[0.033, 0, 1.] )".format(bkgName))
        w.factory("Gaussian::{0}2(bdt, {0}_mean2[-0.05, -1, 1.], {0}_sigma2[0.033, 0, 1.] )".format(bkgName))
        w.factory("Gaussian::{0}3(bdt, {0}_mean3[ 0.03, -1, 1.], {0}_sigma2[0.033, 0, 1.] )".format(bkgName))
        fit_bkg = w.factory("SUM::{0}({0}1, {0}_frac1[0.33, 0, 1]*{0}2, {0}_frac2[.33, 0, 1]*{0}3)".format(bkgName))
    
        return fit_sig, fit_bkg
    
    
    
    #category = "s_vs_b_JKS_BDT_highest"
    category = "PTMAT_s_vs_n_JKS_BDT_highest" 
    mu_val = 0.;
    
    if len(sys.argv) > 2 :
        category = str(sys.argv[2])
    if len(sys.argv) > 3 :
        mu_val = float(sys.argv[3])
    
    #fdir = os.environ['CMSSW_BASE']+"/src/tsW/delphes/TMVA/output/PTMAT_TMVA_JetDiscrimination.root"
    fdir = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/TMVA/output/PTMAT_TMVA_JetDiscrimination.root"
    if not os.path.exists(fdir):
        print "WARNING: Run tmva/jetDiscriminator.C first to get BDT distribution."
        print "WARNING: Stopping running."
        sys.exit()
    
    f = ROOT.TFile(fdir)
    t = f.Get(category+"/Method_BDT/BDT")
    h_sig = t.Get("MVA_BDT_S")
    h_bkg = t.Get("MVA_BDT_B")
    
    w = RooWorkspace("combined_"+str(ntimes))
    
    w.factory("bdt[-1,1]")
    bdt = w.var("bdt")
    rooh_sig = ROOT.RooDataHist("sigHist", "sigHist", ROOT.RooArgList(bdt), h_sig)
    rooh_bkg = ROOT.RooDataHist("bkgHist", "bkgHist", ROOT.RooArgList(bdt), h_bkg)
    getattr(w, 'import')(rooh_sig)
    getattr(w, 'import')(rooh_bkg)
    
    trainType = "JKS"#category.split("_")[0]
    fit_sig, fit_bkg = getFitFunc(w, "JKS")
    fit_sig.fitTo(rooh_sig)
    fit_bkg.fitTo(rooh_bkg)
     
    # ------------------ pad1: BDT distribution -------------------
    frBDT = bdt.frame()
    rooh_sig.plotOn(frBDT, RooFit.MarkerColor(2))
    rooh_bkg.plotOn(frBDT, RooFit.MarkerColor(4))
    fit_sig.plotOn(frBDT, RooFit.LineColor(2))
    fit_bkg.plotOn(frBDT, RooFit.LineColor(4))
    frBDT.SetTitle("BDT distribution")
     
    # ------------------- pad2: Generated Dataset ----------------
    setConst(fit_sig)
    setConst(fit_bkg)
     
    if "JKS" in category:
        print(">>>>>>>>>>>>>>>> Category : JKS")
        expected_sig = expNsig_wKs
        expected_bkg = expNbkg_wKs
    elif "Jet" in category:
        print(">>>>>>>>>>>>>>>> Category : Jet")
        expected_sig = expNsig_woKs
        expected_bkg = expNbkg_woKs
    else:
        print(">>>>>>>>>>>>>>>> Category : " + category)
        expected_sig = expNsig
        expected_bkg = expNbkg
     
    w.factory("mu[1, 0, 4]")
    w.factory("expected_sig[10000]")
    w.var("expected_sig").setVal(expected_sig)
    w.factory("expr::nsig('mu*expected_sig', {expected_sig, mu})")
     
    w.factory("nu[1, 0, 10]")
    w.factory("expected_bkg[1000000]")
    w.var("expected_bkg").setVal(expected_bkg)
    w.factory("expr::nbkg('nu*expected_bkg', {expected_bkg, nu})")
     
    w.factory("SUM::model(nsig*sig_%s, nbkg*bkg_%s)"%(trainType, trainType))
    model = w.pdf("model")
    modelData = model.generate(ROOT.RooArgSet(bdt))
    #modelData = copy.deepcopy(data)

    frDATA = bdt.frame(RooFit.Bins(nPoint))
    #frDATA = bdt.frame(RooFit.Bins(400))
    modelData.plotOn(frDATA, RooFit.MarkerColor(1))
    model.plotOn(frDATA, RooFit.LineColor(1))
    frDATA.SetTitle("Generated Dataset")
    
    # ------------------- pad3: Asimov -------------------
    mu = w.var("mu")
    nu = w.var("nu")
    expected_sig = w.var("expected_sig")
    expected_bkg = w.var("expected_bkg")

    mcfg = RooStats.ModelConfig("mcfg", w)
    mcfg.SetPdf(model)
    mcfg.SetParametersOfInterest("mu")
    mcfg.SetObservables("bdt")
    mcfg.SetNuisanceParameters("nu")
    
    data = ROOT.RooStats.AsymptoticCalculator.MakeAsimovData(modelData, mcfg, ROOT.RooArgSet(mu), ROOT.RooArgSet())
    dataHist = ROOT.RooDataHist("","",w.argSet("bdt"), data)
    for i in range(0, dataHist.numEntries()): correctErr(dataHist, i) #error correction
    
    frAsimov = bdt.frame()
    dataHist.plotOn(frAsimov)
    frAsimov.SetTitle("Asimov Dataset")
    
    # ------------------ pad4: p-value (HypoInv) ----------------
    valueCL = 0.95
    if len(sys.argv)>4:
        valueCL = float(sys.argv[4])
    print("Set Confidence Level : " + str(valueCL))
    useCLs = False
    
    pod = w.var("mu").getVal()
    print(">>> NS    ", w.var("expected_sig").getVal(),     " >>> NB    ", w.var("expected_bkg").getVal())
    print(">>> MU    ", pod,                                " >>> NU    ", w.var("nu").getVal())
    print(">>> Mu*Ns ", pod*w.var("expected_sig").getVal(), " >>> Nu*Nb ", w.var("nu").getVal()*w.var("expected_bkg").getVal())
    
    w.var("mu").setVal(pod)
    mcfg.SetSnapshot(w.argSet("mu"))
    
    bMod = mcfg.Clone()
    bMod.SetName("bkg")
    w.var("mu").setVal(mu_val)
    #bMod.SetParametersOfInterest("mu")
    #bMod.SetObservables("bdt")
    #bMod.SetNuisanceParameters("nu")
    bMod.SetSnapshot(w.argSet("mu")) # Configure bModel to encode current poi=0 scenario as its hypothesis
    
    getattr(w, 'import')(mcfg)
    getattr(w, 'import')(bMod)
    getattr(w, 'import')(data)
    getattr(w, 'import')(modelData)

    #fileName = "s_vs_b_JKS_with_modelData.root"
    fileName = "JKS_with_modelData.root"#"results/JKS_with_modelDta.root"
    if "all" in sys.argv[1] : fileName = "All_JKS_with_modelData.root"
    if "tt"  in sys.argv[1] : fileName = "TT_JKS_with_modelData.root"

    w.writeToFile(fileName,False)#True)

    dirPath = time.strftime("%Y%m%d")
    try: 
        os.mkdir(dirPath)
        print(dirPath+" folder created")
    except OSError:
        print(dirPath+" folder already exists")

    ####################### ProfileLikelihood plotting #####################
    w.Print()
   
    """
    mcfg_profLLC = RooStats.ProfileLikelihoodCalculator(data, mcfg)
    mcfg_profLLC.SetConfidenceLevel(valueCL) # 95% interval
    mcfg_interval = mcfg_profLLC.GetInterval()
    # print out the interval on the first Parameter of Interest
    mcfg_firstPOI = mcfg.GetParametersOfInterest().first()
    print("\n>>>> RESULT : ", valueCL * 100, " % interval on ", mcfg_firstPOI.GetName(), " is : [ Lower : ", mcfg_interval.LowerLimit(mcfg_firstPOI), ", Upper : ", mcfg_interval.UpperLimit(mcfg_firstPOI), "]\n ")
    # make a plot
    mcfg_plot = RooStats.LikelihoodIntervalPlot(mcfg_interval)
    mcfg_plot.SetNPoints(nPoint)  #do not use too many points, it could become very slow for some models
    mcfg_plot.SetRange(0, 4)

    c0 = ROOT.TCanvas()

    mcfg_plot.Draw("TF1")  # use option TF1 if too slow (plot.Draw("tf1")
    mcfg_profLLC_outName = "mcfg_"+str(ntimes)+"_ProfileLikelihood.png"
    c0.SaveAs("%s/%s.png"%(dirPath, mcfg_profLLC_outName))
    c0.SaveAs("%s/%s.pdf"%(dirPath, mcfg_profLLC_outName))

    bMod_profLLC = RooStats.ProfileLikelihoodCalculator(data, bMod)
    bMod_profLLC.SetConfidenceLevel(valueCL) # 95% interval
    bMod_interval = bMod_profLLC.GetInterval()
    # print out the interval on the first Parameter of Interest
    bMod_firstPOI = bMod.GetParametersOfInterest().first()
    print("\n>>>> RESULT : ", valueCL * 100, " % interval on ", bMod_firstPOI.GetName(), " is : [ Lower : ", bMod_interval.LowerLimit(bMod_firstPOI), ", Upper : ", bMod_interval.UpperLimit(bMod_firstPOI), "]\n ")
    # make a plot
    bMod_plot = RooStats.LikelihoodIntervalPlot(bMod_interval)
    bMod_plot.SetNPoints(nPoint)  #do not use too many points, it could become very slow for some models
    bMod_plot.SetRange(0, 4)
    bMod_plot.Draw("TF1")  # use option TF1 if too slow (plot.Draw("tf1")
    #bMod_plot.Print("bMod_"+str(ntimes)+"_ProfileLikelihood.png")
    bMod_profLLC_outName = "bMod_"+str(ntimes)+"_ProfileLikelihood.png"
    c0.SaveAs("%s/%s.png"%(dirPath, bMod_profLLC_outName))
    c0.SaveAs("%s/%s.pdf"%(dirPath, bMod_profLLC_outName))

    #mcfg_bMod_ratio_plot = mcfg_plot.GetPlottedObject().Divide(bMod_plot.GetPlottedObject())
    #mcfg_bMod_ratio_profLLC_outName = "ratio_"+str(ntimes)+"_ProfileLikelihood.png"
    #c0.SaveAs("%s/%s.png"%(dirPath, mcfg_bMod_ratio_profLLC_outName))
    #c0.SaveAs("%s/%s.pdf"%(dirPath, mcfg_bMod_ratio_profLLC_outName))

    #################
    """


    #fcalc = RooStats.FrequentistCalculator(dataHist, bMod, mcfg)
    #fcalc.SetToys(300,300)

    # Asymp. calc
    fcalc = RooStats.AsymptoticCalculator(dataHist, bMod, mcfg)
    #fcalc.SetOneSided(True)
    fcalc.SetOneSidedDiscovery(True)
    #fcalc.SetPrintLevel(-1)
    #fcalc.SetToys(300,300)

    result = fcalc.GetHypoTest()

#    # profile likelihood test statistics 
#    profLL = RooStats.ProfileLikelihoodTestStat(mcfg.GetPdf()) # here we have to tell it that we want the profile likelihood ratio test statistic
#    profLL.SetOneSided(True) # Configure calculator for a limit (=one-sided interval)
#    # for CLs (bounded intervals) use one-sided profile likelihood
#    if (useCLs): profLL.SetOneSided(True) # Configure calculator for a limit (=one-sided interval)
#    
#    # ratio of profile likelihood - need to pass snapshot for the alt 
#    # ropl = RooStats.RatioOfProfiledLikelihoodsTestStat(mcfg.GetPdf(), bMod.GetPdf(), bMod.GetSnapshot())
#    # set the test statistic to use
#
#    toymcs = fcalc.GetTestStatSampler()
#    #toymcs = inverter.GetHypoTestCalculator().GetTestStatSampler()
#    toymcs.SetTestStatistic(profLL)
#
    # a tool that can calculate the POI value takes a certain value. This inversion requires a scan over possible values of mu
    inverter = RooStats.HypoTestInverter(fcalc)
    inverter.SetConfidenceLevel(valueCL) # Statistical configuration of hypothesis test inverter
    inverter.UseCLs(useCLs)              # Statistical configuration of hypothesis test inverter
    inverter.SetFixedScan(nPoint, poiMin, poiMax) # Technical configuration of hypothesis test inverter

    r = inverter.GetInterval() # Perform calculation of limit

    print(str(100.*inverter.ConfidenceLevel())+ " %  upper limit : " + str(r.UpperLimit()))
    print("Expected upper limits, using the B (alternate) model : " )
    print(" expected limit (median) " +str(r.GetExpectedUpperLimit(0)) ) 
    print(" expected limit (-1 sig) " +str(r.GetExpectedUpperLimit(-1)))
    print(" expected limit (+1 sig) " +str(r.GetExpectedUpperLimit(1)) ) 
    print(" expected limit (-2 sig) " +str(r.GetExpectedUpperLimit(-2)))  
    print(" expected limit (+2 sig) " +str(r.GetExpectedUpperLimit(2)) )


    plot = RooStats.HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result with Asymptotic Interval at L = "+str(TotalLumi)+" fb^-1",r) 

    obs = plot.MakePlot()
    exp = plot.MakeExpectedPlot()
    # ------------------ Saving Canvas ------------------
    ww = int(698*1.6)
    wh = int(474*1.6)
    c = ROOT.TCanvas("plots","plots",ww, wh)
    c.Divide(2,2)
    c.cd(1)
    frBDT.Draw()
    c.cd(2)
    frDATA.Draw()
    c.cd(3)
    frAsimov.Draw()
    c.cd(4)
    ROOT.gPad.SetLogy()

    #obs.Draw()
    #exp.Draw("same")
    #exp.GetYaxis().SetRangeUser(10**-8, 1.)
    #exp.GetXaxis().SetRangeUser(poiMin, poiMax)
    #obs.Draw("samepl")
    #drawSigmaLines(ROOT.kRed, poiMin, poiMax)

    plot.Draw("CLb 2CL")
    plot4 = c.GetListOfPrimitives().At(3)
    plot4XMin = plot4.GetUxmin()
    plot4XMax = plot4.GetUxmax()
    drawSigmaLines(ROOT.kRed, plot4XMin, plot4XMax)

    ROOT.gPad.RedrawAxis()
    c.Draw()
    
    outName = "Lumi_"+str(TotalLumi).replace(".", "_")+"fb_"+"pvalue_%s_mu%d_CL%.2f"%(category,mu_val,valueCL)
    if "all" in sys.argv[1] : outName = "All_Lumi_"+str(TotalLumi).replace(".", "_")+"fb_"+"pvalue_%s_mu%d_CL%.2f"%(category,mu_val,valueCL)
    if "tt"  in sys.argv[1] : outName = "TT_Lumi_"+str(TotalLumi).replace(".", "_")+"fb_"+"pvalue_%s_mu%d_CL%.2f"%(category,mu_val,valueCL)

    c.SaveAs("%s/%s.png"%(dirPath, outName))
    c.SaveAs("%s/%s.pdf"%(dirPath, outName))
    
    #pad4 = c.GetPad(4)
    #pad4.SaveAs("%s/%s_pad4.png"%(dirPath, outName))
    #
    #if os.path.exists("vts.wiki"):
    #    c.SaveAs("vts.wiki/plots/%s.png"%(outName))
    #    c.SaveAs("vts.wiki/plots/%s.pdf"%(outName))
    #    pad4.SaveAs("vts.wiki/plots/%s_pad4.png"%(outName))
    
    # ----------------- Save numbers in text file -------------------
    f = open(dirPath+"/model_output.txt", "a")
    if "all" in sys.argv[1] : f = open(dirPath+"/All_model_output.txt", "a")
    if "tt"  in sys.argv[1] : f = open(dirPath+"/TT_model_output.txt", "a")
    
    expValueCL95 = getExpValue(inverter, 0.9545, False)
    expValueCL99 = getExpValue(inverter, 0.9973, False)
    pValueAtMuZero = inverter.GetInterval().CLsplusb(0)
    
    f.write(str(TotalLumi) + " fb^-1 >>> " + outName+", p-value at mu = 0, " + str(pValueAtMuZero) + ", mean 95% upper limit, " + str(expValueCL95) + ", mean 99% upper limit, " + str(expValueCL99) + "\n")
    f.write("\n")
    f.close()
    continue 
