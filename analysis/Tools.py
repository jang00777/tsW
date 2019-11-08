import ROOT as r

def toVector(typ, l):
    v = r.vector(typ)()
    for i in l: v.push_back(i)
    return v

# GOOD LUMIBLOCK
r.gInterpreter.LoadFile("external/lumiTool.cc")
lumi_code ='''
#ifndef LUMITOOL__
#define LUMITOOL__
static lumiTool lumitool("data/lumi/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt");

Bool_t LumiCheck(UInt_t run, UInt_t luminosityBlock) {
  return lumitool.LumiCheck(run, luminosityBlock);
}
#endif
'''
r.gInterpreter.Declare(lumi_code)

# PILEUP REWEIGHTING
r.gInterpreter.LoadFile("external/pileUpTool.cc")
pileup_code ='''
#ifndef PILEUPTOOL__
#define PILEUPTOOL__
static pileUpTool puTool;

float PileupGetWeight(int nTrueInt, int sys=0) {
  return puTool.getWeight(nTrueInt, sys);
}
#endif
'''
r.gInterpreter.Declare(pileup_code)

# LEPTON SCALE FACTORS
r.gInterpreter.LoadFile("external/computePtEtaTable.cc")
scalefactor_code = '''
#ifndef LEPSFTOOL__
#define LEPSFTOOL__
  string m_strTrigSFEl   = "data/SF/SF_El_Trig_2016_Ele32_WPTight.csv";
  string m_strLeptonSFEl = "data/SF/SF_El_ID_2016_Tight.csv";
  string m_strRecSFEl    = "data/SF/SF_El_Rec_2016.csv";

  string m_strTrigSFMu   = "data/SF/SF_Mu_Trig_2016_IsoTrk24_OR_Iso24.csv";
  string m_strLeptonSFMu = "data/SF/SF_Mu_ID_2016_Tight.csv";
  string m_strIsoSFMu    = "data/SF/SF_Mu_Iso_2016.csv";

  computePtEtaTable trigSFEl(m_strTrigSFEl);
  computePtEtaTable elecSF_(m_strLeptonSFEl);
  computePtEtaTable elecSFRec_(m_strRecSFEl);

  computePtEtaTable trigSFMu(m_strTrigSFMu);
  computePtEtaTable muonSF_(m_strLeptonSFMu);
  computePtEtaTable muonSFIso_(m_strIsoSFMu);
#endif
'''
r.gInterpreter.Declare(scalefactor_code)

ini = [lumi_code, pileup_code, scalefactor_code]
