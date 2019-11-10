import ROOT as r

def toVector(typ, l):
    v = r.vector(typ)()
    for i in l: v.push_back(i)
    return v

# GOOD LUMIBLOCK
r.gInterpreter.LoadFile("external/lumiTool.cc")
#define LUMITOOL__

# PILEUP REWEIGHTING
r.gInterpreter.LoadFile("external/pileUpTool.cc")
#define PILEUPTOOL__

# LEPTON SCALE FACTORS
r.gInterpreter.LoadFile("external/computePtEtaTable.cc")
#define LEPSFTOOL__
