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
