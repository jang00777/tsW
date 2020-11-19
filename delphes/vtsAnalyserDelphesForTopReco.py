import os, sys, copy
import xml.etree.ElementTree as ET
from array import array
from math import sqrt
from ROOT import *


def listdir_full(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

addCut    = ""
if len(sys.argv) > 1 :  
    addCut    = sys.argv[1]

base_path_bkg       = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result_bkg/"
base_path_sig_bbars = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result/tt012j_bbars_2l_FxFx/"
base_path_sig_bsbar = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result/tt012j_bsbar_2l_FxFx/"
base_path_bkg_bbbar = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result/tt012j_bbbar_2l_FxFx/"

out_path            = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result_topReco/" 

sig_path  = listdir_full(base_path_sig_bbars) + listdir_full(base_path_sig_bsbar)
bkg_path  = listdir_full(base_path_bkg_bbbar) + listdir_full(base_path_bkg+"TT_s/") + listdir_full(base_path_bkg+"TT_ss/") + listdir_full(base_path_bkg+"Z_2j/") + listdir_full(base_path_bkg+"Z_bb/")

channel   = "DL"

lumi2016  = 35.922 # fb-1 
lumi2017  = 41.86 # fb-1 
lumi2018  = 58.83 # fb-1 
totalLumi = lumi2016 + lumi2017 + lumi2018 # fb-1

xsec_ttjj = 831.36 # pb
xsec_sig  = 0.2815844982 #0.1689786308 # pb
xsec_bkg  = 88.76650802  #53.26871002 # pb
xsec_tts  = 9.407 # pb
xsec_ttss = 1.567 # pb
xsec_zjj  = 5653  # pb
xsec_zbb  = 820   # pb
xsec_zcc  = 60.78 # pb


for ibkg in sig_path + bkg_path:
#for ibkg in bkg_path[:1]:
#for ibkg in sig_path[:1]:
    if not "tt012j_bbars_2l" in ibkg and not "tt012j_bsbar_2l" in ibkg and not "tt012j_bbbar_2l" in ibkg and not "dilep_" in ibkg : continue
    if   channel == "DL" and ( "semi_"  in ibkg or "_1l" in ibkg) : continue
    elif channel == "SL" and ( "dilep_" in ibkg or "_2l" in ibkg) : continue
    tree_file  = TFile(ibkg)
    tree_file2 = TFile(ibkg)
    tree_evt  = tree_file.Get("event")
    tree_jet1 = tree_file.Get("MVA_jet")
    tree_jet2 = tree_file2.Get("MVA_jet")

    if tree_file.IsZombie() : continue
    try :
        nentries = tree_jet1.GetEntries()
    except AttributeError:
        print(">>>> "+ ibkg + " has an error ... SKIP this file ...")
        continue

    nentries = tree_jet1.GetEntries()

    outName = out_path + ibkg.split("/")[-2] + "/TopReco"+ibkg.split("Out")[-1]
    if os.path.isdir(outName.split("TopReco")[0]) == False : 
        print(outName.split("TopReco")[0] + " IS NOT EXSIST ===> LET'S MAKE THE FOLDER")
        os.system("mkdir -p " + outName.split("TopReco")[0])
    #outf    = TFile(outName, "RECREATE")
    #outtr   = TTree("MVA_jet_pair", "MVA_jet_pair")

    evt_step                      = array('i', [-99])
    evt_nJet                      = array('i', [-99])
    evt_nBJet                     = array('i', [-99])

    isFromTop                     = array('I', [False])

    jet1_jetFlavour               = array('i', [-99]) 
    jet1_isFrom                   = array('i', [-99]) 
    jet1_hasHighestPt             = array('I', [False]) 
    jet1_hasHighestPtFromSelected = array('I', [False]) 
    jet1_bTag                     = array('i', [-99])

    jet2_jetFlavour               = array('i', [-99]) 
    jet2_isFrom                   = array('i', [-99]) 
    jet2_hasHighestPt             = array('I', [False])
    jet2_hasHighestPtFromSelected = array('I', [False])
    jet2_bTag                     = array('i', [-99])

    jet1_jet2_dr                  = array('d', [-99.])

    jet1_pt                       = array('d', [-99.])
    jet1_eta                      = array('d', [-99.])
    jet1_phi                      = array('d', [-99.])
    jet1_mass                     = array('d', [-99.])
    jet1_lep1_deta                = array('d', [-99.])
    jet1_lep1_dphi                = array('d', [-99.])
    jet1_lep1_dr                  = array('d', [-99.])
    jet1_lep2_deta                = array('d', [-99.])
    jet1_lep2_dphi                = array('d', [-99.])
    jet1_lep2_dr                  = array('d', [-99.])
    jet1_dlep_deta                = array('d', [-99.])
    jet1_dlep_dphi                = array('d', [-99.])
    jet1_dlep_dr                  = array('d', [-99.])
    jet1_MET_deta                 = array('d', [-99.])
    jet1_MET_dphi                 = array('d', [-99.])
    jet1_lep_closest_deta         = array('d', [-99.]) 
    jet1_lep_closest_dphi         = array('d', [-99.]) 
    jet1_lep_closest_dr           = array('d', [-99.]) 
    jet1_lep_matched_idx          = array('i', [-99])
    jet1_lep_matched_pt           = array('d', [-99.]) 
    jet1_lep_matched_eta          = array('d', [-99.]) 
    jet1_lep_matched_phi          = array('d', [-99.]) 
    jet1_lep_matched_mass         = array('d', [-99.]) 
    jet1_lep_matched_deta         = array('d', [-99.]) 
    jet1_lep_matched_dphi         = array('d', [-99.]) 
    jet1_lep_matched_dr           = array('d', [-99.]) 

    jet2_pt                       = array('d', [-99.])
    jet2_eta                      = array('d', [-99.])
    jet2_phi                      = array('d', [-99.])
    jet2_mass                     = array('d', [-99.])
    jet2_lep1_deta                = array('d', [-99.])
    jet2_lep1_dphi                = array('d', [-99.])
    jet2_lep1_dr                  = array('d', [-99.])
    jet2_lep2_deta                = array('d', [-99.])
    jet2_lep2_dphi                = array('d', [-99.])
    jet2_lep2_dr                  = array('d', [-99.])
    jet2_dlep_deta                = array('d', [-99.])
    jet2_dlep_dphi                = array('d', [-99.])
    jet2_dlep_dr                  = array('d', [-99.])
    jet2_MET_deta                 = array('d', [-99.])
    jet2_MET_dphi                 = array('d', [-99.])
    jet2_lep_closest_deta         = array('d', [-99.])
    jet2_lep_closest_dphi         = array('d', [-99.])
    jet2_lep_closest_dr           = array('d', [-99.])
    jet2_lep_matched_idx          = array('i', [-99])
    jet2_lep_matched_pt           = array('d', [-99.])
    jet2_lep_matched_eta          = array('d', [-99.])
    jet2_lep_matched_phi          = array('d', [-99.])
    jet2_lep_matched_mass         = array('d', [-99.])
    jet2_lep_matched_deta         = array('d', [-99.])
    jet2_lep_matched_dphi         = array('d', [-99.])
    jet2_lep_matched_dr           = array('d', [-99.])

    lep1_pt                       = array('d', [-99.])
    lep1_eta                      = array('d', [-99.])
    lep1_phi                      = array('d', [-99.])
    lep1_mass                     = array('d', [-99.])
    lep1_MET_deta                 = array('d', [-99.])
    lep1_MET_dphi                 = array('d', [-99.])
    lep1_pdgId                    = array('i', [-99])

    lep2_pt                       = array('d', [-99.])
    lep2_eta                      = array('d', [-99.])
    lep2_phi                      = array('d', [-99.])
    lep2_mass                     = array('d', [-99.])
    lep2_MET_deta                 = array('d', [-99.])
    lep2_MET_dphi                 = array('d', [-99.])
    lep2_pdgId                    = array('i', [-99])

    dlep_pt                       = array('d', [-99.])
    dlep_eta                      = array('d', [-99.])
    dlep_phi                      = array('d', [-99.])
    dlep_mass                     = array('d', [-99.])
    dlep_MET_deta                 = array('d', [-99.])
    dlep_MET_dphi                 = array('d', [-99.])

    MET_pt                        = array('d', [-99.])
    MET_eta                       = array('d', [-99.])
    MET_phi                       = array('d', [-99.])


    outf    = TFile(outName, "RECREATE")
    outtr   = TTree("MVA_jet_pair", "MVA_jet_pair")

    outtr.Branch("evt_step",                      evt_step                     , "evt_step/I")
    outtr.Branch("evt_nJet",                      evt_nJet                     , "evt_nJet/I")     
    outtr.Branch("evt_nBJet",                     evt_nBJet                    , "evt_nBJet/I")    

    outtr.Branch("isFromTop",                     isFromTop                    , "isFromTop/O")
                                                                                                               
    outtr.Branch("jet1_jetFlavour",               jet1_jetFlavour              , "jet1_jetFlavour/I")    
    outtr.Branch("jet1_isFrom",                   jet1_isFrom                  , "jet1_isFrom/I")    
    outtr.Branch("jet1_hasHighestPt",             jet1_hasHighestPt            , "jet1_hasHighestPt/O")    
    outtr.Branch("jet1_hasHighestPtFromSelected", jet1_hasHighestPtFromSelected, "jet1_hasHighestPtFromSelected/O")    
    outtr.Branch("jet1_bTag",                     jet1_bTag                    , "jet1_bTag/I");
                                                                                                              
    outtr.Branch("jet2_jetFlavour",               jet2_jetFlavour              , "jet2_jetFlavour/I")  
    outtr.Branch("jet2_isFrom",                   jet2_isFrom                  , "jet2_isFrom/I")  
    outtr.Branch("jet2_hasHighestPt",             jet2_hasHighestPt            , "jet2_hasHighestPt/O")  
    outtr.Branch("jet2_hasHighestPtFromSelected", jet2_hasHighestPtFromSelected, "jet2_hasHighestPtFromSelected/O")  
    outtr.Branch("jet2_bTag",                     jet2_bTag                    , "jet2_bTag/I");

    outtr.Branch("jet1_jet2_dr",                  jet1_jet2_dr                 , "jet1_jet2_dr/D")  
                                                                                                
    outtr.Branch("jet1_pt",                       jet1_pt                      , "jet1_pt/D")  
    outtr.Branch("jet1_eta",                      jet1_eta                     , "jet1_eta/D")  
    outtr.Branch("jet1_phi",                      jet1_phi                     , "jet1_phi/D")  
    outtr.Branch("jet1_mass",                     jet1_mass                    , "jet1_mass/D")  
    outtr.Branch("jet1_lep1_deta",                jet1_lep1_deta               , "jet1_lep1_deta/D")  
    outtr.Branch("jet1_lep1_dphi",                jet1_lep1_dphi               , "jet1_lep1_dphi/D")  
    outtr.Branch("jet1_lep1_dr",                  jet1_lep1_dr                 , "jet1_lep1_dr/D")  
    outtr.Branch("jet1_lep2_deta",                jet1_lep2_deta               , "jet1_lep2_deta/D")  
    outtr.Branch("jet1_lep2_dphi",                jet1_lep2_dphi               , "jet1_lep2_dphi/D")  
    outtr.Branch("jet1_lep2_dr",                  jet1_lep2_dr                 , "jet1_lep2_dr/D")  
    outtr.Branch("jet1_dlep_deta",                jet1_dlep_deta               , "jet1_dlep_deta/D")  
    outtr.Branch("jet1_dlep_dphi",                jet1_dlep_dphi               , "jet1_dlep_dphi/D")  
    outtr.Branch("jet1_dlep_dr",                  jet1_dlep_dr                 , "jet1_dlep_dr/D")  
    outtr.Branch("jet1_MET_deta",                 jet1_MET_deta                , "jet1_MET_deta/D")  
    outtr.Branch("jet1_MET_dphi",                 jet1_MET_dphi                , "jet1_MET_dphi/D")  
    outtr.Branch("jet1_lep_closest_deta",         jet1_lep_closest_deta        , "jet1_lep_closest_deta/D")
    outtr.Branch("jet1_lep_closest_dphi",         jet1_lep_closest_dphi        , "jet1_lep_closest_dphi/D")
    outtr.Branch("jet1_lep_closest_dr",           jet1_lep_closest_dr          , "jet1_lep_closest_dr/D")
    outtr.Branch("jet1_lep_matched_idx",          jet1_lep_matched_idx         , "jet1_lep_matched_idx/I")
    outtr.Branch("jet1_lep_matched_pt",           jet1_lep_matched_pt          , "jet1_lep_matched_pt/D")
    outtr.Branch("jet1_lep_matched_eta",          jet1_lep_matched_eta         , "jet1_lep_matched_eta/D")
    outtr.Branch("jet1_lep_matched_phi",          jet1_lep_matched_phi         , "jet1_lep_matched_phi/D")
    outtr.Branch("jet1_lep_matched_mass",         jet1_lep_matched_mass        , "jet1_lep_matched_mass/D")
    outtr.Branch("jet1_lep_matched_deta",         jet1_lep_matched_deta        , "jet1_lep_matched_deta/D")
    outtr.Branch("jet1_lep_matched_dphi",         jet1_lep_matched_dphi        , "jet1_lep_matched_dphi/D")
    outtr.Branch("jet1_lep_matched_dr",           jet1_lep_matched_dr          , "jet1_lep_matched_dr/D")

    outtr.Branch("jet2_pt",                       jet2_pt                      , "jet2_pt/D")  
    outtr.Branch("jet2_eta",                      jet2_eta                     , "jet2_eta/D")  
    outtr.Branch("jet2_phi",                      jet2_phi                     , "jet2_phi/D")  
    outtr.Branch("jet2_mass",                     jet2_mass                    , "jet2_mass/D")  
    outtr.Branch("jet2_lep1_deta",                jet2_lep1_deta               , "jet2_lep1_deta/D")  
    outtr.Branch("jet2_lep1_dphi",                jet2_lep1_dphi               , "jet2_lep1_dphi/D")  
    outtr.Branch("jet2_lep1_dr",                  jet2_lep1_dr                 , "jet2_lep1_dr/D")  
    outtr.Branch("jet2_lep2_deta",                jet2_lep2_deta               , "jet2_lep2_deta/D")  
    outtr.Branch("jet2_lep2_dphi",                jet2_lep2_dphi               , "jet2_lep2_dphi/D")  
    outtr.Branch("jet2_lep2_dr",                  jet2_lep2_dr                 , "jet2_lep2_dr/D")  
    outtr.Branch("jet2_dlep_deta",                jet2_dlep_deta               , "jet2_dlep_deta/D")  
    outtr.Branch("jet2_dlep_dphi",                jet2_dlep_dphi               , "jet2_dlep_dphi/D")  
    outtr.Branch("jet2_dlep_dr",                  jet2_dlep_dr                 , "jet2_dlep_dr/D")  
    outtr.Branch("jet2_MET_deta",                 jet2_MET_deta                , "jet2_MET_deta/D")  
    outtr.Branch("jet2_MET_dphi",                 jet2_MET_dphi                , "jet2_MET_dphi/D")  
    outtr.Branch("jet2_lep_closest_deta",         jet2_lep_closest_deta        , "jet2_lep_closest_deta/D")
    outtr.Branch("jet2_lep_closest_dphi",         jet2_lep_closest_dphi        , "jet2_lep_closest_dphi/D")
    outtr.Branch("jet2_lep_closest_dr",           jet2_lep_closest_dr          , "jet2_lep_closest_dr/D")
    outtr.Branch("jet2_lep_matched_idx",          jet2_lep_matched_idx         , "jet2_lep_matched_idx/I")
    outtr.Branch("jet2_lep_matched_pt",           jet2_lep_matched_pt          , "jet2_lep_matched_pt/D")
    outtr.Branch("jet2_lep_matched_eta",          jet2_lep_matched_eta         , "jet2_lep_matched_eta/D")
    outtr.Branch("jet2_lep_matched_phi",          jet2_lep_matched_phi         , "jet2_lep_matched_phi/D")
    outtr.Branch("jet2_lep_matched_mass",         jet2_lep_matched_mass        , "jet2_lep_matched_mass/D")
    outtr.Branch("jet2_lep_matched_deta",         jet2_lep_matched_deta        , "jet2_lep_matched_deta/D")
    outtr.Branch("jet2_lep_matched_dphi",         jet2_lep_matched_dphi        , "jet2_lep_matched_dphi/D")
    outtr.Branch("jet2_lep_matched_dr",           jet2_lep_matched_dr          , "jet2_lep_matched_dr/D")
                                                                                               
    outtr.Branch("lep1_pt",                       lep1_pt                      , "lep1_pt/D")  
    outtr.Branch("lep1_eta",                      lep1_eta                     , "lep1_eta/D")  
    outtr.Branch("lep1_phi",                      lep1_phi                     , "lep1_phi/D")  
    outtr.Branch("lep1_mass",                     lep1_mass                    , "lep1_mass/D")  
    outtr.Branch("lep1_MET_deta",                 lep1_MET_deta                , "lep1_MET_deta/D")  
    outtr.Branch("lep1_MET_dphi",                 lep1_MET_dphi                , "lep1_MET_dphi/D")  
    outtr.Branch("lep1_pdgId",                    lep1_pdgId                   , "lep1_pdgId/I")                                                                                                
    outtr.Branch("lep2_pt",                       lep2_pt                      , "lep2_pt/D")  
    outtr.Branch("lep2_eta",                      lep2_eta                     , "lep2_eta/D")  
    outtr.Branch("lep2_phi",                      lep2_phi                     , "lep2_phi/D")  
    outtr.Branch("lep2_mass",                     lep2_mass                    , "lep2_mass/D")  
    outtr.Branch("lep2_MET_deta",                 lep2_MET_deta                , "lep2_MET_deta/D")  
    outtr.Branch("lep2_MET_dphi",                 lep2_MET_dphi                , "lep2_MET_dphi/D")  
    outtr.Branch("lep2_pdgId",                    lep2_pdgId                   , "lep2_pdgId/I")                                                                                                
    outtr.Branch("dlep_pt",                       dlep_pt                      , "dlep_pt/D")  
    outtr.Branch("dlep_eta",                      dlep_eta                     , "dlep_eta/D")  
    outtr.Branch("dlep_phi",                      dlep_phi                     , "dlep_phi/D")  
    outtr.Branch("dlep_mass",                     dlep_mass                    , "dlep_mass/D")  
    outtr.Branch("dlep_MET_deta",                 dlep_MET_deta                , "dlep_MET_deta/D")  
    outtr.Branch("dlep_MET_dphi",                 dlep_MET_dphi                , "dlep_MET_dphi/D")  
                                                                                                
    outtr.Branch("MET_pt",                        MET_pt                       , "MET_pt/D")  
    outtr.Branch("MET_eta",                       MET_eta                      , "MET_eta/D")  
    outtr.Branch("MET_phi",                       MET_phi                      , "MET_phi/D")  


    print(">>>>>>>>>>>>>>>>>>>>>>>",ibkg)
    nevent = tree_evt.GetEntries()

    for i in range(nevent) :
        tree_evt.GetEntry(i)
        if tree_evt.step < 4: continue
        if tree_evt.jet_start == tree_evt.jet_end : continue
        if not "noCut" in str(addCut) and "nJet" in str(addCut) and tree_evt.nJet != int(addCut.split("nJet")[-1]) : continue # nJet Cut

        lep1_tlv = tree_evt.recoLep1
        lep2_tlv = tree_evt.recoLep2
        dlep_tlv = lep1_tlv + lep2_tlv

        evt_step[0]       = tree_evt.step
        evt_nJet[0]       = tree_evt.nJet
        evt_nBJet[0]      = tree_evt.nBJet

        lep1_pt[0]        = lep1_tlv.Pt()
        lep1_eta[0]       = lep1_tlv.Eta()
        lep1_phi[0]       = lep1_tlv.Phi()
        lep1_mass[0]      = lep1_tlv.M()
        lep1_MET_deta[0]  = lep1_tlv.Eta() - tree_evt.MET_eta
        lep1_MET_dphi[0]  = TVector2.Phi_mpi_pi(lep1_tlv.Phi() - tree_evt.MET_phi)
        lep1_pdgId[0]     = tree_evt.recoLep1_pdgId

        lep2_pt[0]        = lep2_tlv.Pt()
        lep2_eta[0]       = lep2_tlv.Eta()
        lep2_phi[0]       = lep2_tlv.Phi()
        lep2_mass[0]      = lep2_tlv.M()
        lep2_MET_deta[0]  = lep2_tlv.Eta() - tree_evt.MET_eta
        lep2_MET_dphi[0]  = TVector2.Phi_mpi_pi(lep2_tlv.Phi() - tree_evt.MET_phi)
        lep2_pdgId[0]     = tree_evt.recoLep2_pdgId

        dlep_pt[0]        = dlep_tlv.Pt()
        dlep_eta[0]       = dlep_tlv.Eta()
        dlep_phi[0]       = dlep_tlv.Phi()
        dlep_mass[0]      = dlep_tlv.M()
        dlep_MET_deta[0]  = dlep_tlv.Eta() - tree_evt.MET_eta
        dlep_MET_dphi[0]  = TVector2.Phi_mpi_pi(dlep_tlv.Phi() - tree_evt.MET_phi)

        MET_pt[0]         = tree_evt.MET
        MET_eta[0]        = tree_evt.MET_eta
        MET_phi[0]        = tree_evt.MET_phi

        for iev in range(tree_evt.jet_start, tree_evt.jet_end) :
            tree_jet1.GetEntry(iev)
            if "btag" in str(addCut) :
                if tree_jet1.bTag != 1 : continue
                jet2_start = bkg_evt.jet_start
            else :
                jet2_start = iev+1

            jet1_tlv = TLorentzVector()
            jet1_tlv.SetPtEtaPhiM(tree_jet1.pt, tree_jet1.eta, tree_jet1.phi, tree_jet1.mass)

            jet1_jetFlavour[0]               = tree_jet1.pdgId
            jet1_isFrom[0]                   = tree_jet1.isFrom
            jet1_hasHighestPt[0]             = tree_jet1.hasHighestPt
            jet1_hasHighestPtFromSelected[0] = tree_jet1.hasHighestPtFromSelected
            jet1_bTag[0]                     = tree_jet1.bTag

            jet1_pt[0]                       = tree_jet1.pt
            jet1_eta[0]                      = tree_jet1.eta
            jet1_phi[0]                      = tree_jet1.phi
            jet1_mass[0]                     = tree_jet1.mass
            jet1_lep1_deta[0]                = jet1_tlv.Eta() - lep1_tlv.Eta()
            jet1_lep1_dphi[0]                = jet1_tlv.DeltaPhi(lep1_tlv)
            jet1_lep1_dr[0]                  = jet1_tlv.DeltaR(lep1_tlv)
            jet1_lep2_deta[0]                = jet1_tlv.Eta() - lep2_tlv.Eta()
            jet1_lep2_dphi[0]                = jet1_tlv.DeltaPhi(lep2_tlv)
            jet1_lep2_dr[0]                  = jet1_tlv.DeltaR(lep2_tlv)
            jet1_dlep_deta[0]                = jet1_tlv.Eta() - dlep_tlv.Eta()
            jet1_dlep_dphi[0]                = jet1_tlv.DeltaPhi(dlep_tlv)
            jet1_dlep_dr[0]                  = jet1_tlv.DeltaR(dlep_tlv)
            jet1_MET_deta[0]                 = jet1_tlv.Eta() - tree_evt.MET_eta
            jet1_MET_dphi[0]                 = TVector2.Phi_mpi_pi(jet1_tlv.Phi() - tree_evt.MET_phi)

            if jet1_lep1_dr[0] < jet1_lep2_dr[0] :
                jet1_lep_closest_deta[0] = jet1_lep1_deta[0] 
                jet1_lep_closest_dphi[0] = jet1_lep1_dphi[0]
                jet1_lep_closest_dr[0]   = jet1_lep1_dr[0]
            elif jet1_lep1_dr[0] > jet1_lep2_dr[0] :
                jet1_lep_closest_deta[0] = jet1_lep2_deta[0]
                jet1_lep_closest_dphi[0] = jet1_lep2_dphi[0]
                jet1_lep_closest_dr[0]   = jet1_lep2_dr[0]
            else :
                jet1_lep_closest_deta[0] = -99. 
                jet1_lep_closest_dphi[0] = -99. 
                jet1_lep_closest_dr[0]   = -99. 

            if abs(jet1_isFrom[0]) == 5 :
                if   jet1_isFrom[0] * lep1_pdgId[0] < 0 :
                    jet1_lep_matched_idx[0]  = 1
                    jet1_lep_matched_tlv     = jet1_tlv + lep1_tlv
                    jet1_lep_matched_pt[0]   = jet1_lep_matched_tlv.Pt()   
                    jet1_lep_matched_eta[0]  = jet1_lep_matched_tlv.Eta() 
                    jet1_lep_matched_phi[0]  = jet1_lep_matched_tlv.Phi()  
                    jet1_lep_matched_mass[0] = jet1_lep_matched_tlv.M() 
                    jet1_lep_matched_deta[0] = jet1_lep1_deta[0] 
                    jet1_lep_matched_dphi[0] = jet1_lep1_dphi[0] 
                    jet1_lep_matched_dr[0]   = jet1_lep1_dr[0]  
                elif jet1_isFrom[0] * lep2_pdgId[0] < 0 : 
                    jet1_lep_matched_idx[0]  = 2
                    jet1_lep_matched_tlv     = jet1_tlv + lep2_tlv
                    jet1_lep_matched_pt[0]   = jet1_lep_matched_tlv.Pt()
                    jet1_lep_matched_eta[0]  = jet1_lep_matched_tlv.Eta()
                    jet1_lep_matched_phi[0]  = jet1_lep_matched_tlv.Phi()
                    jet1_lep_matched_mass[0] = jet1_lep_matched_tlv.M() 
                    jet1_lep_matched_deta[0] = jet1_lep2_deta[0]
                    jet1_lep_matched_dphi[0] = jet1_lep2_dphi[0]
                    jet1_lep_matched_dr[0]   = jet1_lep2_dr[0]
                else :
                    jet1_lep_matched_idx[0]  = -99
                    jet1_lep_matched_pt[0]   = -99. 
                    jet1_lep_matched_eta[0]  = -99. 
                    jet1_lep_matched_phi[0]  = -99. 
                    jet1_lep_matched_mass[0] = -99. 
                    jet1_lep_matched_deta[0] = -99. 
                    jet1_lep_matched_dphi[0] = -99. 
                    jet1_lep_matched_dr[0]   = -99.
            else :
                jet1_lep_matched_idx[0]  = -99
                jet1_lep_matched_pt[0]   = -99.
                jet1_lep_matched_eta[0]  = -99.
                jet1_lep_matched_phi[0]  = -99.
                jet1_lep_matched_mass[0] = -99.
                jet1_lep_matched_deta[0] = -99.
                jet1_lep_matched_dphi[0] = -99.
                jet1_lep_matched_dr[0]   = -99.
 
            for jev in range(jet2_start, tree_evt.jet_end) :
                tree_jet2.GetEntry(jev)  
                if "btag" in str(addCut) :
                    if jev == iev : continue # Avoid self-pairing
                    if tree_jet2.bTag and jev < iev : continue # Avoid double-counting

                jet2_tlv = TLorentzVector()
                jet2_tlv.SetPtEtaPhiM(tree_jet2.pt, tree_jet2.eta, tree_jet2.phi, tree_jet2.mass)
                
                jet2_lep_matched_pt  
                jet2_lep_matched_eta 
                jet2_lep_matched_phi 
                jet2_lep_matched_mass
                jet2_lep_matched_deta
                jet2_lep_matched_dphi
                jet2_lep_matched_dr  

                jet1_jet2_dr[0]                   = jet1_tlv.DeltaR(jet2_tlv) 

                jet2_jetFlavour[0]               = tree_jet2.pdgId
                jet2_isFrom[0]                   = tree_jet2.isFrom
                jet2_hasHighestPt[0]             = tree_jet2.hasHighestPt
                jet2_hasHighestPtFromSelected[0] = tree_jet2.hasHighestPtFromSelected
                jet2_bTag[0]                     = tree_jet2.bTag
 
                jet2_pt[0]                       = tree_jet2.pt 
                jet2_eta[0]                      = tree_jet2.eta 
                jet2_phi[0]                      = tree_jet2.phi
                jet2_mass[0]                     = tree_jet2.mass 
                jet2_lep1_deta[0]                = jet2_tlv.Eta() - lep1_tlv.Eta() 
                jet2_lep1_dphi[0]                = jet2_tlv.DeltaPhi(lep1_tlv) 
                jet2_lep1_dr[0]                  = jet2_tlv.DeltaR(lep1_tlv) 
                jet2_lep2_deta[0]                = jet2_tlv.Eta() - lep2_tlv.Eta() 
                jet2_lep2_dphi[0]                = jet2_tlv.DeltaPhi(lep2_tlv) 
                jet2_lep2_dr[0]                  = jet2_tlv.DeltaR(lep2_tlv) 
                jet2_dlep_deta[0]                = jet2_tlv.Eta() - dlep_tlv.Eta() 
                jet2_dlep_dphi[0]                = jet2_tlv.DeltaPhi(dlep_tlv) 
                jet2_dlep_dr[0]                  = jet2_tlv.DeltaR(dlep_tlv)
                jet2_MET_deta[0]                 = jet2_tlv.Eta() - tree_evt.MET_eta
                jet2_MET_dphi[0]                 = TVector2.Phi_mpi_pi(jet2_tlv.Phi() - tree_evt.MET_phi)

                if jet2_lep1_dr[0] < jet2_lep2_dr[0] :
                    jet2_lep_closest_deta[0] = jet2_lep1_deta[0]
                    jet2_lep_closest_dphi[0] = jet2_lep1_dphi[0]
                    jet2_lep_closest_dr[0]   = jet2_lep1_dr[0]
                elif jet2_lep1_dr[0] > jet2_lep2_dr[0] :
                    jet2_lep_closest_deta[0] = jet2_lep2_deta[0]
                    jet2_lep_closest_dphi[0] = jet2_lep2_dphi[0]
                    jet2_lep_closest_dr[0]   = jet2_lep2_dr[0]
                else :
                    jet2_lep_closest_deta[0] = -99.
                    jet2_lep_closest_dphi[0] = -99.
                    jet2_lep_closest_dr[0]   = -99.
   
                if abs(jet1_isFrom[0]) != 99 and abs(jet2_isFrom[0]) != 99 : isFromTop[0] = True
                else : isFromTop[0] = False

                bjet_lep_tlv1 = TLorentzVector()
                if abs(jet2_isFrom[0]) == 5 :
                    if   jet2_isFrom[0] * lep1_pdgId[0] < 0 :
                        jet2_lep_matched_idx[0]  = 1
                        jet2_lep_matched_tlv     = jet2_tlv + lep1_tlv
                        jet2_lep_matched_pt[0]   = jet2_lep_matched_tlv.Pt()   
                        jet2_lep_matched_eta[0]  = jet2_lep_matched_tlv.Eta() 
                        jet2_lep_matched_phi[0]  = jet2_lep_matched_tlv.Phi()  
                        jet2_lep_matched_mass[0] = jet2_lep_matched_tlv.M() 
                        jet2_lep_matched_deta[0] = jet2_lep1_deta[0] 
                        jet2_lep_matched_dphi[0] = jet2_lep1_dphi[0] 
                        jet2_lep_matched_dr[0]   = jet2_lep1_dr[0]  
                    elif jet2_isFrom[0] * lep2_pdgId[0] < 0 : 
                        jet2_lep_matched_idx[0]  = 2
                        jet2_lep_matched_tlv     = jet2_tlv + lep2_tlv
                        jet2_lep_matched_pt[0]   = jet2_lep_matched_tlv.Pt()
                        jet2_lep_matched_eta[0]  = jet2_lep_matched_tlv.Eta()
                        jet2_lep_matched_phi[0]  = jet2_lep_matched_tlv.Phi()
                        jet2_lep_matched_mass[0] = jet2_lep_matched_tlv.M() 
                        jet2_lep_matched_deta[0] = jet2_lep2_deta[0]
                        jet2_lep_matched_dphi[0] = jet2_lep2_dphi[0]
                        jet2_lep_matched_dr[0]   = jet2_lep2_dr[0]
                    else :
                        jet2_lep_matched_idx[0]  = -99
                        jet2_lep_matched_pt[0]   = -99. 
                        jet2_lep_matched_eta[0]  = -99. 
                        jet2_lep_matched_phi[0]  = -99. 
                        jet2_lep_matched_mass[0] = -99. 
                        jet2_lep_matched_deta[0] = -99. 
                        jet2_lep_matched_dphi[0] = -99. 
                        jet2_lep_matched_dr[0]   = -99.
                else : 
                    jet2_lep_matched_idx[0]  = -99
                    jet2_lep_matched_pt[0]   = -99.
                    jet2_lep_matched_eta[0]  = -99.
                    jet2_lep_matched_phi[0]  = -99.
                    jet2_lep_matched_mass[0] = -99.
                    jet2_lep_matched_deta[0] = -99.
                    jet2_lep_matched_dphi[0] = -99.
                    jet2_lep_matched_dr[0]   = -99.
 
                #print("IEV : ", i, " | jet1 idx : ", iev, " | jet2 idx : ", jev, " | jet1 isFrom : ", jet1_isFrom[0], " | jet2 isFrom : ", jet2_isFrom[0], " | jet1 bTag : ", jet1_bTag[0], " | jet2 bTag : ", jet2_bTag[0], " | isFromTop : ", isFromTop[0])

                outtr.Fill()

    outf.Write()
    outf.Close()
