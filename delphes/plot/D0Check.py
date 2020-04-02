from math import sqrt, log10
from ROOT import *
import os, sys, array
from datetime import*
from pytz import timezone

now_time     = datetime.now().strftime('%Y%m%d')

gSystem.Load("libDelphes")
gInterpreter.Declare('#include "classes/DelphesClasses.h"')
gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

outf_path = "./outroot/"

inFile = TFile("input/sum_tt012j_bbars_2l_FxFx_elIso03_muIso04_newResForm_noVtxSmearing_20200319.root")
inTree = inFile.Get("MVA_had")
inFile = TFile("input/delphes_tt012j_bbars_2l_FxFx_elIso03_muIso04_newResForm_noVtxSmearing_20200316.root")
inTree = inFile.Get("Delphes")

def getMotherList(tree, gen) :
    motherList = []
    idx = gen.M1
    if idx == -1 :
        return motherList
    while True :
        mother = tree.Particle.At(idx)
        motherList.append(mother)
        idx = mother.M1
        if idx == -1 : break
    return motherList

def checkFromTop(moms) :
    for imom in moms :
        if imom.PID == 6 and imom.Status == 62 :
            return True
    return False

opt_draw_first = "norm hist"
opt_draw_other = "norm hist same"

binning = [100, 0, 300]
binEta  = [100, 0, 3.5]
if inFile.GetName().find("delphes_") != -1 :
    if sys.argv[1].find("pion") != -1:
        outf = TFile(outf_path+now_time+"_D0Check.root", "RECREATE")
        c = TCanvas("c", "c", 1920, 1080)
        c.Divide(2,2)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        nentries = inTree.GetEntries()
        h1  = TH1D("RxyGenPion",          "Rxy of Gen Pions;R_{xy};Entries",                                    *binning)
        h2  = TH1D("RxyGenPionFromKs",    "Rxy of Gen Pions From Ks;R_{xy};Entries",                            *binning)
        h3  = TH1D("D0RecoPion",          "D0 of Reco Pions;|D0|;Entries",                                      *binning)
        h4  = TH1D("D0RecoPionFromKs",    "D0 of Reco Pions From Ks;|D0|;Entries",                              *binning)
        h5  = TH1D("RxyRecoPionFromKs",   "Rxy of Reco Pions;R_{xy};Entries",                                   *binning)
        h6  = TH1D("RxyRecoPionFromKs",   "Rxy of Reco Pions From Ks;R_{xy};Entries",                           *binning)
        h7  = TH1D("D0SigRecoPionFromKs", "Significance of D0 of Reco Pions;|D0 Significance|;Entries",         *binning)
        h8  = TH1D("D0SigRecoPionFromKs", "Significance of D0 of Reco Pions From Ks;|D0 Significance|;Entries", *binning)

        for iev in range(1000):#range(nentries):
            print("[ {0:>6} / {1:>6} ] th Entry".format(str(iev), str(nentries)))
            inTree.GetEntry(iev)
            # Gen
            for i in range(len(inTree.Particle)):
                if inTree.Particle[i] == None : break
                Rxy = sqrt(inTree.Particle[i].X*inTree.Particle[i].X +inTree.Particle[i].Y*inTree.Particle[i].Y)
                if abs(inTree.Particle[i].PID) != 211 : continue
                isFromKs  = False
                moms = getMotherList(inTree, inTree.Particle[i])
                for im in range(len(moms)):
                    if moms[0].PID == 310 : 
                        isFromKs = True
                        break
                    else: break
                if isFromKs == True : 
                    h2.Fill(Rxy)
                else :
                    h1.Fill(Rxy)
            # Reco
            for i in range(len(inTree.Track)):
                if inTree.Track[i] == None : break
                isFromKs  = False
                genp = inTree.Track[i].Particle.GetObject()
                if abs(genp.PID) != 211 : continue 
                Rxy = sqrt(inTree.Track[i].X*inTree.Track[i].X +inTree.Track[i].Y*inTree.Track[i].Y) 
                moms = getMotherList(inTree, genp)
                for im in range(len(moms)):
                    if moms[0].PID == 310 :
                        isFromKs = True
                        break
                    else: break
                if isFromKs == True :
                    h4.Fill(abs(inTree.Track[i].D0))
                    h6.Fill(Rxy)
                    if inTree.Track[i].ErrorD0 != 0 : h8.Fill(abs(inTree.Track[i].D0/float(inTree.Track[i].ErrorD0)))
                else :
                    h3.Fill(abs(inTree.Track[i].D0))
                    h5.Fill(Rxy)
                    if inTree.Track[i].ErrorD0 != 0 : h7.Fill(abs(inTree.Track[i].D0/float(inTree.Track[i].ErrorD0)))

        c.cd(1)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        h1.SetLineColor(4)
        h1.Draw()
        h2.SetLineColor(2)
        h2.Draw("same")
        tl1 = TLegend(0.6, 0.6, 0.9, 0.9)
        tl1.AddEntry(h1, "Gen Pions not from Ks")
        tl1.AddEntry(h2, "Gen Pions from Ks")
        tl1.Draw("same")

        c.cd(3)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        h3.SetLineColor(4)
        h3.Draw()
        h4.SetLineColor(2)
        h4.Draw("same")
        tl2 = TLegend(0.6, 0.6, 0.9, 0.9)
        tl2.AddEntry(h3, "Reco Pions not from Ks")
        tl2.AddEntry(h4, "Reco Pions from Ks")
        tl2.Draw("same")

        c.cd(2)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        h5.SetLineColor(4)
        h5.Draw()
        h6.SetLineColor(2)
        h6.Draw("same")
        tl2.Draw("same")

        c.cd(4)
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)
        h7.SetLineColor(4)
        h7.Draw()
        h8.SetLineColor(2)
        h8.Draw("same")
        tl2.Draw("same")
        c.SaveAs("CheckD0.png")
        outf.Write()
        outf.Close()
    elif sys.argv[1].find("ks") != -1 :
        outf = TFile(outf_path+now_time+"_MassCheck.root", "RECREATE")
        #c = TCanvas("c", "c", 1920, 1080)
        #c.Divide(2,2)
        #gPad.SetLogy(1)
        #gStyle.SetOptStat(0)
        nentries = inTree.GetEntries()

        ksMass_   = 0.49761
        pionMass_ = 0.13967

        binningPt   = [100,    0,  100]
        binningEta  = [35,     0,  3.5]
        binningPhi  = [35,     0,  3.5]
        binningMass = [100,  0.3,  0.7]
        binningD0   = [100,    0,  200]
        binningDS   = [100,    0,  200]
        binningDE   = [100,    0,    1]

        binningPtP  = [100, -0.5,  0.5]
        binningEtaP = [100, -0.5,  0.5]
        binningPhiP = [100,   -5,    5]

        binningPtR  = [100, -0.5,  0.5]
        binningEtaR = [100, -0.5,  0.5]
        binningPhiR = [100,   -5,    5]


        binningPt2    = binningPt   + binningPt
        binningEta2   = binningEta  + binningEta
        binningPhi2   = binningPhi  + binningPhi
        binningMass2  = binningMass + binningMass
        binningMassD0 = binningMass + binningD0
        binningMassDS = binningMass + binningDS
        binningMassDE = binningMass + binningDE
        binningD0DS   = binningD0   + binningDS
        binningD0DE   = binningD0   + binningDE

        binningPtDS   = binningPt   + binningDS
        binningEtaDS  = binningEta  + binningDS
        binningPhiDS  = binningPhi  + binningDS

        binningPtPDS  = binningPtP   + binningDS
        binningEtaPDS = binningEtaP  + binningDS
        binningPhiPDS = binningPhiP  + binningDS

        binningPtRDS  = binningPtR   + binningDS
        binningEtaRDS = binningEtaR  + binningDS
        binningPhiRDS = binningPhiR  + binningDS


        hGenMAD0  = TH2D("hGenMAD0",  "Gen Ks D0 vs mass     w/o D0 sig > 5 cut;Gen mass (GeV);D0",            *binningMassD0)
        hGenMADS  = TH2D("hGenMADS",  "Gen Ks D0 sig vs mass w/o D0 sig > 5 cut;Gen mass (GeV);D0_{sig}",      *binningMassDS)
        hGenMADE  = TH2D("hGenMADE",  "Gen Ks D0 Err vs mass w/o D0 sig > 5 cut;Gen mass (GeV);D0_{Err}",      *binningMassDE)

        hPT    = TH2D("hPT",    "Gen vs Reco Ks pt   w/o D0 sig > 5 cut;Reco p_{T} (GeV);Gen p_{T} (GeV)", *binningPt2) 
        hETA   = TH2D("hETA",   "Gen vs Reco Ks eta  w/o D0 sig > 5 cut;Reco |#eta|;Gen |#eta|",           *binningEta2)  
        hPHI   = TH2D("hPHI",   "Gen vs Reco Ks phi  w/o D0 sig > 5 cut;Reco |#phi|;Gen |#phi|",           *binningPhi2) 
        hMASS  = TH2D("hMASS",  "Gen vs Reco Ks mass w/o D0 sig > 5 cut;Reco mass (GeV);Gen mass (GeV)",   *binningMass2) 
        hMAD0  = TH2D("hMAD0",  "Reco Ks D0 vs mass     w/o D0 sig > 5 cut;Reco mass (GeV);D0",            *binningMassD0)
        hMADS  = TH2D("hMADS",  "Reco Ks D0 sig vs mass w/o D0 sig > 5 cut;Reco mass (GeV);D0_{sig}",      *binningMassDS)
        hMADE  = TH2D("hMADE",  "Reco Ks D0 Err vs mass w/o D0 sig > 5 cut;Reco mass (GeV);D0_{Err}",      *binningMassDE)
        hD0DS  = TH2D("hD0DS",  "Reco Ks D0 sig vs D0   w/o D0 sig > 5 cut;D0;D0_{sig}",                   *binningD0DS)
        hD0DE  = TH2D("hD0DE",  "Reco Ks D0 Err vs D0   w/o D0 sig > 5 cut;D0;D0_{Err}",                   *binningD0DE)

        hPTDS    = TH2D("hPTDS",    "Reco Ks D0 sig vs pT w/o D0 sig > 5 cut;p_{T} (GeV);D0_{sig}",     *binningPtDS)
        hETADS   = TH2D("hETADS",   "Reco Ks D0 sig vs eta w/o D0 sig > 5 cut;|#eta|;D0_{sig}",         *binningEtaDS)
        hPHIDS   = TH2D("hPHIDS",   "Reco Ks D0 sig vs phi w/o D0 sig > 5 cut;|#phi|;D0_{sig}",         *binningPhiDS)
        hPTPDS   = TH2D("hPTPDS",   "Reco Ks D0 sig vs Pull of pT w/o D0 sig > 5 cut;Pull_{p_{T}};D0_{sig}",         *binningPtPDS)
        hPTRDS   = TH2D("hPTRDS",   "Reco Ks D0 sig vs Residual of pT w/o D0 sig > 5 cut;Residual_{p_{T}};D0_{sig}", *binningPtRDS)
        hETAPDS  = TH2D("hETAPDS",  "Reco Ks D0 sig vs Pull of eta w/o D0 sig > 5 cut;Pull_{#eta};D0_{sig}",         *binningEtaPDS)
        hETARDS  = TH2D("hETARDS",  "Reco Ks D0 sig vs Residual of eta w/o D0 sig > 5 cut;Residual_{#eta};D0_{sig}", *binningEtaRDS)
        hPHIPDS  = TH2D("hPHIPDS",  "Reco Ks D0 sig vs Pull of phi w/o D0 sig > 5 cut;Pull_{#phi};D0_{sig}",         *binningPhiPDS)
        hPHIRDS  = TH2D("hPHIRDS",  "Reco Ks D0 sig vs Residual of phi w/o D0 sig > 5 cut;Residual_{#phi};D0_{sig}", *binningPhiRDS)

        hPTc   = TH2D("hPTc",   "Gen vs Reco Ks pt   w D0 sig > 5 cut;Reco p_{T} (GeV);Gen p_{T} (GeV)",    *binningPt2)    
        hETAc  = TH2D("hETAc",  "Gen vs Reco Ks eta  w D0 sig > 5 cut;Reco |#eta|;Gen |#eta|",              *binningEta2) 
        hPHIc  = TH2D("hPHIc",  "Gen vs Reco Ks phi  w D0 sig > 5 cut;Reco |#phi|;Gen |#phi|",              *binningPhi2) 
        hMASSc = TH2D("hMASSc", "Gen vs Reco Ks mass w D0 sig > 5 cut;Reco mass (GeV);Gen mass (GeV)",      *binningMass2)  
        hMAD0c = TH2D("hMAD0c", "Reco Ks D0 vs mass     w D0 sig > 5 cut;Reco mass (GeV);D0",               *binningMassD0)
        hMADSc = TH2D("hMADSc", "Reco Ks D0 sig vs mass w D0 sig > 5 cut;Reco mass (GeV);D0_{sig}",         *binningMassDS)
        hMADEc = TH2D("hMADEc", "Reco Ks D0 Err vs mass w D0 sig > 5 cut;Reco mass (GeV);D0_{Err}",         *binningMassDE)
        hD0DSc = TH2D("hD0DSc", "Reco Ks D0 sig vs D0   w D0 sig > 5 cut;D0;D0_{sig}",                      *binningD0DS)
        hD0DEc = TH2D("hD0DEc", "Reco Ks D0 Err vs D0   w D0 sig > 5 cut;D0;D0_{Err}",                      *binningD0DE)

        hPTDSc   = TH2D("hPTDSc",    "Reco Ks D0 sig vs pT w D0 sig > 5 cut;p_{T} (GeV);D0_{sig}",     *binningPtDS)
        hETADSc  = TH2D("hETADSc",   "Reco Ks D0 sig vs eta w D0 sig > 5 cut;|#eta|;D0_{sig}",         *binningEtaDS)
        hPHIDSc  = TH2D("hPHIDSc",   "Reco Ks D0 sig vs phi w D0 sig > 5 cut;|#phi|;D0_{sig}",         *binningPhiDS)
        hPTPDSc  = TH2D("hPTPDSc",   "Reco Ks D0 sig vs Pull of pT w D0 sig > 5 cut;Pull_{p_{T}};D0_{sig}",         *binningPtPDS)
        hPTRDSc  = TH2D("hPTRDSc",   "Reco Ks D0 sig vs Residual of pT w D0 sig > 5 cut;Residual_{p_{T}};D0_{sig}", *binningPtRDS)
        hETAPDSc = TH2D("hETAPDSc",  "Reco Ks D0 sig vs Pull of eta w D0 sig > 5 cut;Pull_{#eta};D0_{sig}",         *binningEtaPDS)
        hETARDSc = TH2D("hETARDSc",  "Reco Ks D0 sig vs Residual of eta w D0 sig > 5 cut;Residual_{#eta};D0_{sig}", *binningEtaRDS)
        hPHIPDSc = TH2D("hPHIPDSc",  "Reco Ks D0 sig vs Pull of phi w D0 sig > 5 cut;Pull_{#phi};D0_{sig}",         *binningPhiPDS)
        hPHIRDSc = TH2D("hPHIRDSc",  "Reco Ks D0 sig vs Residual of phi w D0 sig > 5 cut;Residual_{#phi};D0_{sig}", *binningPhiRDS)


        h1  = TH1D("h1",  "TruthMatched Ks pt     w/o D0 sig cut > 5;Ks p_{T} (GeV);Entries", *binningPt)
        h2  = TH1D("h2",  "TruthMatched Ks eta    w/o D0 sig cut > 5;Ks |#eta|;Entries",      *binningEta)
        h3  = TH1D("h3",  "TruthMatched Ks phi    w/o D0 sig cut > 5;Ks |#phi|;Entries",      *binningPhi)
        h4  = TH1D("h4",  "TruthMatched Ks mass   w/o D0 sig cut > 5;Ks mass (GeV);Entries",  *binningMass)
        h5  = TH1D("h5",  "TruthMatched Ks D0;D0  w/o D0 sig cut > 5;Entries",                *binningD0)
        h6  = TH1D("h6",  "TruthMatched Ks D0 sig w/o D0 sig cut > 5;D0_{sig};Entries",       *binningDS)

        h1c = TH1D("h1c", "TruthMatched Ks pt     w D0 sig > 5 cut;Ks p_{T} (GeV);Entries",   *binningPt)
        h2c = TH1D("h2c", "TruthMatched Ks eta    w D0 sig > 5 cut;Ks |#eta|;Entries",        *binningEta)
        h3c = TH1D("h3c", "TruthMatched Ks phi    w D0 sig > 5 cut;Ks |#phi|;Entries",        *binningPhi)
        h4c = TH1D("h4c", "TruthMatched Ks mass   w D0 sig > 5 cut;Ks mass (GeV);Entries",    *binningMass)
        h5c = TH1D("h5c", "TruthMatched Ks D0     w D0 sig > 5 cut;D0;Entries",               *binningD0)
        h6c = TH1D("h6c", "TruthMatched Ks D0 sig w D0 sig > 5 cut;D0_{sig};Entries",         *binningDS)

        for iev in range(10000):#range(nentries):
            print("[ {0:>6} / {1:>6} ] th Entry".format(str(iev), str(nentries)))
            inTree.GetEntry(iev)
            # Reco
            for i in range(len(inTree.Track)):
                if inTree.Track[i] == None : break
                dau1 = inTree.Track[i]
                if dau1.Charge != 1 : continue
                if abs(dau1.PID) == 11 or abs(dau1.PID) == 13 : continue
                for j in range(len(inTree.Track)):
                    if inTree.Track[j] == None : break
                    dau2 = inTree.Track[j]
                    if dau2.Charge != -1 : continue
                    if abs(dau2.PID) == 11 or abs(dau2.PID) == 13 : continue
                    dau1_tlv = TLorentzVector()
                    dau2_tlv = TLorentzVector()
                    dau1_tlv.SetPtEtaPhiM(dau1.PT, dau1.Eta, dau1.Phi, pionMass_)
                    dau2_tlv.SetPtEtaPhiM(dau2.PT, dau2.Eta, dau2.Phi, pionMass_)
                    had_tlv  = dau1_tlv + dau2_tlv
                    if abs(had_tlv.M() - ksMass_)/float(ksMass_) < 0.3 : 
                        #print("Invariant mass cut (0.3) : " + str(abs(had_tlv.M() - ksMass_)/float(ksMass_)))
                        genDau1 = dau1.Particle.GetObject()
                        genDau2 = dau2.Particle.GetObject()
                        if genDau1.M1 == genDau2.M1 and inTree.Particle[genDau1.M1].PID == 310:
                            genDau1_tlv = genDau1.P4()
                            genDau2_tlv = genDau2.P4()
                            genHad_tlv = genDau1_tlv + genDau2_tlv
                            hGenMAD0.Fill(genHad_tlv.M(), abs(dau1.D0)) 
                            hGenMAD0.Fill(genHad_tlv.M(), abs(dau2.D0))
                            hGenMADE.Fill(genHad_tlv.M(), abs(dau1.ErrorD0))
                            hGenMADE.Fill(genHad_tlv.M(), abs(dau2.ErrorD0)) 
                            #print("Truth matching")
                            genKs = inTree.Particle[genDau1.M1]
                            hPT.Fill(had_tlv.Pt(),    genKs.PT)
                            hETA.Fill(had_tlv.Eta(),  genKs.Eta)
                            hPHI.Fill(had_tlv.Phi(),  genKs.Phi)
                            hMASS.Fill(had_tlv.M(),   genKs.Mass)
                            hMAD0.Fill(had_tlv.M(),   abs(dau1.D0))
                            hMAD0.Fill(had_tlv.M(),   abs(dau2.D0))
                            hMADE.Fill(had_tlv.M(),   abs(dau1.ErrorD0))
                            hMADE.Fill(had_tlv.M(),   abs(dau2.ErrorD0))
                            hD0DE.Fill(abs(dau1.D0),  abs(dau1.ErrorD0))
                            hD0DE.Fill(abs(dau1.D0),  abs(dau2.ErrorD0))
                            h1.Fill(had_tlv.Pt())
                            h2.Fill(had_tlv.Eta())
                            h3.Fill(had_tlv.Phi())
                            h4.Fill(had_tlv.M())
                            h5.Fill(dau1.D0)
                            h5.Fill(dau2.D0)
                            if dau1.ErrorD0 != 0 : 
                                h6.Fill(abs(dau1.D0/dau1.ErrorD0))
                                hGenMADS.Fill(genHad_tlv.M(),  abs(dau1.D0/dau1.ErrorD0))
                                hMADS.Fill(had_tlv.M(),  abs(dau1.D0/dau1.ErrorD0))
                                hD0DS.Fill(abs(dau1.D0), abs(dau1.D0/dau1.ErrorD0))

                                hPTDS.Fill(had_tlv.Pt(),   abs(dau1.D0/dau1.ErrorD0))  
                                hETADS.Fill(abs(had_tlv.Eta()), abs(dau1.D0/dau1.ErrorD0))  
                                hPHIDS.Fill(abs(had_tlv.Phi()), abs(dau1.D0/dau1.ErrorD0))  
                                hPTPDS.Fill((had_tlv.Pt()-genHad_tlv.Pt()),                     abs(dau1.D0/dau1.ErrorD0))  
                                hPTRDS.Fill((had_tlv.Pt()-genHad_tlv.Pt())/genHad_tlv.Pt(),     abs(dau1.D0/dau1.ErrorD0))   
                                hETAPDS.Fill((had_tlv.Eta()-genHad_tlv.Eta()),                  abs(dau1.D0/dau1.ErrorD0))  
                                hETARDS.Fill((had_tlv.Eta()-genHad_tlv.Eta())/genHad_tlv.Eta(), abs(dau1.D0/dau1.ErrorD0))  
                                hPHIPDS.Fill((had_tlv.Phi()-genHad_tlv.Phi()),                  abs(dau1.D0/dau1.ErrorD0))  
                                hPHIRDS.Fill((had_tlv.Phi()-genHad_tlv.Phi())/genHad_tlv.Phi(), abs(dau1.D0/dau1.ErrorD0))   
                            if dau2.ErrorD0 != 0 : 
                                h6.Fill(abs(dau2.D0/dau2.ErrorD0))
                                hGenMADS.Fill(genHad_tlv.M(),  abs(dau2.D0/dau2.ErrorD0))
                                hMADS.Fill(had_tlv.M(),  abs(dau2.D0/dau2.ErrorD0))
                                hD0DS.Fill(abs(dau2.D0), abs(dau2.D0/dau2.ErrorD0))

                                hPTDS.Fill(had_tlv.Pt(),   abs(dau2.D0/dau2.ErrorD0))
                                hETADS.Fill(abs(had_tlv.Eta()), abs(dau2.D0/dau2.ErrorD0))
                                hPHIDS.Fill(abs(had_tlv.Phi()), abs(dau2.D0/dau2.ErrorD0))
                                hPTPDS.Fill((had_tlv.Pt()-genHad_tlv.Pt()),                     abs(dau2.D0/dau2.ErrorD0))
                                hPTRDS.Fill((had_tlv.Pt()-genHad_tlv.Pt())/genHad_tlv.Pt(),     abs(dau2.D0/dau2.ErrorD0))
                                hETAPDS.Fill((had_tlv.Eta()-genHad_tlv.Eta()),                  abs(dau2.D0/dau2.ErrorD0))
                                hETARDS.Fill((had_tlv.Eta()-genHad_tlv.Eta())/genHad_tlv.Eta(), abs(dau2.D0/dau2.ErrorD0))
                                hPHIPDS.Fill((had_tlv.Phi()-genHad_tlv.Phi()),                  abs(dau2.D0/dau2.ErrorD0))
                                hPHIRDS.Fill((had_tlv.Phi()-genHad_tlv.Phi())/genHad_tlv.Phi(), abs(dau2.D0/dau2.ErrorD0))
                            if dau1.ErrorD0 != 0 and dau2.ErrorD0 != 0 :
                                if abs(dau1.D0/dau1.ErrorD0) > 5 and abs(dau2.D0/dau2.ErrorD0) > 5 :
                                    hPTc.Fill(had_tlv.Pt(),   genKs.PT)
                                    hETAc.Fill(had_tlv.Eta(), genKs.Eta)
                                    hPHIc.Fill(had_tlv.Phi(), genKs.Phi)
                                    hMASSc.Fill(had_tlv.M(),  genKs.Mass)
                                    hMAD0c.Fill(had_tlv.M(),  abs(dau1.D0))
                                    hMAD0c.Fill(had_tlv.M(),  abs(dau2.D0))
                                    hMADSc.Fill(had_tlv.M(),  abs(dau1.D0/dau1.ErrorD0))
                                    hMADSc.Fill(had_tlv.M(),  abs(dau2.D0/dau2.ErrorD0))
                                    hMADEc.Fill(had_tlv.M(),  abs(dau1.ErrorD0))
                                    hMADEc.Fill(had_tlv.M(),  abs(dau2.ErrorD0))
                                    hD0DSc.Fill(abs(dau1.D0), abs(dau1.D0/dau1.ErrorD0))
                                    hD0DSc.Fill(abs(dau2.D0), abs(dau2.D0/dau2.ErrorD0))
                                    hD0DEc.Fill(abs(dau1.D0), abs(dau1.ErrorD0))
                                    hD0DEc.Fill(abs(dau2.D0), abs(dau2.ErrorD0))

                                    hPTDSc.Fill(had_tlv.Pt(),   abs(dau1.D0/dau1.ErrorD0))
                                    hETADSc.Fill(abs(had_tlv.Eta()), abs(dau1.D0/dau1.ErrorD0))
                                    hPHIDSc.Fill(abs(had_tlv.Phi()), abs(dau1.D0/dau1.ErrorD0))
                                    hPTPDSc.Fill((had_tlv.Pt()-genHad_tlv.Pt()),                     abs(dau1.D0/dau1.ErrorD0))
                                    hPTRDSc.Fill((had_tlv.Pt()-genHad_tlv.Pt())/genHad_tlv.Pt(),     abs(dau1.D0/dau1.ErrorD0))
                                    hETAPDSc.Fill((had_tlv.Eta()-genHad_tlv.Eta()),                  abs(dau1.D0/dau1.ErrorD0))
                                    hETARDSc.Fill((had_tlv.Eta()-genHad_tlv.Eta())/genHad_tlv.Eta(), abs(dau1.D0/dau1.ErrorD0))
                                    hPHIPDSc.Fill((had_tlv.Phi()-genHad_tlv.Phi()),                  abs(dau1.D0/dau1.ErrorD0))
                                    hPHIRDSc.Fill((had_tlv.Phi()-genHad_tlv.Phi())/genHad_tlv.Phi(), abs(dau1.D0/dau1.ErrorD0))
                                    hPTDSc.Fill(had_tlv.Pt(),   abs(dau2.D0/dau2.ErrorD0))
                                    hETADSc.Fill(abs(had_tlv.Eta()), abs(dau2.D0/dau2.ErrorD0))
                                    hPHIDSc.Fill(abs(had_tlv.Phi()), abs(dau2.D0/dau2.ErrorD0))
                                    hPTPDSc.Fill((had_tlv.Pt()-genHad_tlv.Pt()),                     abs(dau2.D0/dau2.ErrorD0))
                                    hPTRDSc.Fill((had_tlv.Pt()-genHad_tlv.Pt())/genHad_tlv.Pt(),     abs(dau2.D0/dau2.ErrorD0))
                                    hETAPDSc.Fill((had_tlv.Eta()-genHad_tlv.Eta()),                  abs(dau2.D0/dau2.ErrorD0))
                                    hETARDSc.Fill((had_tlv.Eta()-genHad_tlv.Eta())/genHad_tlv.Eta(), abs(dau2.D0/dau2.ErrorD0))
                                    hPHIPDSc.Fill((had_tlv.Phi()-genHad_tlv.Phi()),                  abs(dau2.D0/dau2.ErrorD0))
                                    hPHIRDSc.Fill((had_tlv.Phi()-genHad_tlv.Phi())/genHad_tlv.Phi(), abs(dau2.D0/dau2.ErrorD0))

                                    h1c.Fill(had_tlv.Pt())
                                    h2c.Fill(had_tlv.Eta())
                                    h3c.Fill(had_tlv.Phi())
                                    h4c.Fill(had_tlv.M())
                                    h5c.Fill(dau1.D0)
                                    h5c.Fill(dau2.D0)
                                    h6c.Fill(abs(dau1.D0/dau1.ErrorD0))
                                    h6c.Fill(abs(dau2.D0/dau2.ErrorD0))


        c = TCanvas()
        gPad.SetLogy(1)
        gStyle.SetOptStat(0)

        hGenMAD0.Draw("colz") 
        hGenMADS.Draw("colz") 
        hGenMADE.Draw("colz") 

        hPT.Draw("colz")
        hETA.Draw("colz")
        hPHI.Draw("colz")
        hMASS.Draw("colz")
        hMAD0.Draw("colz")

        hPTc.Draw("colz")
        hETAc.Draw("colz")
        hPHIc.Draw("colz") 
        hMASSc.Draw("colz") 
        hMAD0c.Draw("colz") 

        h1.Draw() 
        h2.Draw() 
        h3.Draw() 
        h4.Draw() 
        h5.Draw() 
        h6.Draw() 

        h1c.Draw()
        h2c.Draw() 
        h3c.Draw() 
        h4c.Draw() 
        h5c.Draw() 
        h6c.Draw() 

        outf.Write()
        outf.Close()

    #        # Gen
    #        for i in range(len(inTree.Particle)):
    #            if inTree.Particle[i] == None : break
    #            if inTree.Particle[i].PID != 310 : continue
    #            dau1_idx = inTree.Particle[i].D1
    #            dau2_idx = inTree.Particle[i].D2
    #            if inTree.Particle[dau1_idx] == -1 or inTree.Particle[dau2_idx] == -1 : continue
    #            if abs(inTree.Particle[dau1_idx].PID) != 211 or abs(inTree.Particle[dau2_idx].PID) != 211 : continue
    #            if inTree.Particle[dau1_idx].PID * inTree.Particle[dau2_idx].PID > 0 : continue

    #            isFromW   = False
    #            isFromTop = False
    #            isFrom    = -99
    #            moms = getMotherList(inTree, inTree.Particle[i])
    #            for j in range(len(moms)):
    #                if abs(moms[j].PID) == 24 and isFromW == False : 
    #                    isFromW = True
    #                    isFrom  = moms[j-1].PID
    #                if abs(moms[j].PID) == 6 and moms[j].Status == 62 : 
    #                   isFromTop == True
    #                   if isFromW == False : isFrom = moms[j-1].PID
    #                   break

    else :
        nentries = inTree.GetEntries()
        for iev in range(nentries):
            inTree.GetEntry(iev)
            for i in range(len(inTree.Particle)):
                if inTree.Particle[i] == None : break
                if inTree.Particle[i].PID != 310 : continue
                if inTree.Particle[i].D1 == -1 and inTree.Particle[i].D2 == -1 : continue
                np = 0
                charge = 1
                for idau in range(inTree.Particle[i].D1,inTree.Particle[i].D2+1) :
                  if abs(inTree.Particle[idau].PID) == 211 : 
                    np += 1
                    charge = charge*inTree.Particle[idau].Charge/(211.)
                    print(str(np) + " th charged pion index : " + str(idau))
                print("no. daughter : " + str(len(inTree.Particle[inTree.Particle[i].D1:inTree.Particle[i].D2+1])) + " and charged pion : " + str(np) + " charge mult. : " + str(charge))    
else :
    outf = TFile(outf_path+now_time+"_RecoKsMassCheck.root", "RECREATE")
    binningMass = [40, 0.3, 0.7]

    hAllNoCut    = TH1D("hAllNoCut",    "All reco Ks mass without d0sig > 5;Ks Mass (GeV);Entries",          *binningMass)
    hAllSigCut   = TH1D("hAllSigCut",   "All reco Ks mass with d0sig > 5;Ks Mass (Gev);Entries",             *binningMass)
    hTruthNoCut  = TH1D("hTruthNoCut",  "TruthMatched reco Ks mass without d0sig > 5;Ks Mass (GeV);Entries", *binningMass)
    hTrtuhSigCut = TH1D("hTruthSigCut", "TruthMatched reco Ks mass with d0sig > 5;Ks Mass (Gev);Entries",    *binningMass)


    cut_match     = TCut("had_nMatched == 2")
    cut_not_match = TCut("had_nMatched != 2")
    cut_d0sig5    = TCut("fabs(dau1_D0/dau1_D0Err) > 5 && fabs(dau2_D0/dau2_D0Err) > 5")
    cut_ks        = TCut("had_pdgId == 310")
    cut_hadgen    = TCut("had_pdgId == genHad_pdgId")

    inTree.Draw("had_mass >> hAllNoCut",    cut_ks, "")
    inTree.Draw("had_mass >> hAllSigCut",   cut_ks + cut_d0sig5, "") 
    inTree.Draw("had_mass >> hTruthNoCut",  cut_match + cut_ks, "")
    inTree.Draw("had_mass >> hTruthSigCut", cut_match + cut_ks + cut_d0sig5, "")

    c = TCanvas("c", "c", 1920, 1080)
    c.Divide(2,1)

    c.cd(1)
    gPad.SetLogy(1)
    gStyle.SetOptStat(0)
    hAllNoCut.SetLineColor(4)   
    hAllSigCut.SetLineColor(2)
    hAllNoCut.GetYaxis().SetRangeUser(1, 200000)

    hAllNoCut.Draw()   
    hAllSigCut.Draw("same")  
    tl1 = TLegend(0.6, 0.6, 0.9, 0.9)
    tl1.AddEntry(hAllNoCut,  "All reco Ks w/o D0_{sig} > 5")
    tl1.AddEntry(hAllSigCut, "All reco Ks w D0_{sig} > 5")
    tl1.Draw("same")

    c.cd(2)
    #gPad.SetLogy(1)
    gStyle.SetOptStat(0)
    hTruthNoCut.SetLineColor(4) 
    hTrtuhSigCut.SetLineColor(2)
    hTruthNoCut.Draw() 
    hTrtuhSigCut.Draw("same")
    tl2 = TLegend(0.6, 0.6, 0.9, 0.9)
    tl2.AddEntry(hTruthNoCut,  "TruthMatched reco Ks w/o D0_{sig} > 5")
    tl2.AddEntry(hTruthSigCut, "TruthMatched reco Ks w D0_{sig} > 5")
    tl2.Draw("same")

    c.SaveAs("RecoKsMass.png")
    outf.Write()
    outf.Close()


    
