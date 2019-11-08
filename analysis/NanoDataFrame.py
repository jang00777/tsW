import ROOT as r
import Tools

initialize = Tools.ini

def ini():
    from glob import glob
    # for f in glob('external/*cc'):
    #     r.gInterpreter.LoadFile(f)
    r.gInterpreter.LoadFile('external/lumiTool.cc')
#    r.gInterpreter.LoadFile('external/computePtEtaTable.cc')
    r.gInterpreter.LoadFile('external/pileUpTool.cc')
    for code in initialize:
        r.gInterpreter.Declare(code)

class NanoDataFrame(object):

    """Switch between PyRDF and ROOTs RDataFrame, plus a bunch of helpers
for e.g. keeping track of weights.

    """
    
    def __init__(self, tree=None, files=None, spark=False, npartitions=172, initialize=None, weight=None):
        if tree is None:
            return
        if spark:
            import PyRDF
            # According to https://jaceklaskowski.gitbooks.io/mastering-apache-spark/spark-rdd-partitions.html
            # we want 2-3x the number of partitions as executors, and according to
            # PyRDF.current_backend.sparkContext.defaultParallelism
            # we have 88
            PyRDF.use('spark', {'npartitions':npartitions})
            PyRDF.initialize(ini)
            self.df = PyRDF.RDataFrame("Events", files)
        else:
            self.df = r.RDataFrame("Events", files)
        self.weight = weight
        self.filters = [["Overall", self.df.Count(),  self.df.Sum(self.weight) if self.weight else self.df.Count()]]

    def copy(self, new_df):
        new = type(self)()
        new.df = new_df
        new.filters = self.filters
        new.weight = self.weight
        return new

    def Define(self, *args):
        return self.copy(self.df.Define(*args))

    def Filter(self, cut, name=None, weight=None):
        filt = self.df.Filter(cut)
        if name:
            count = filt.Count()
            weight = filt.Sum(weight) if weight else (filt.Sum(self.weight) if self.weight else filt.Count())
            new = self.copy(filt)
            new.filters = self.filters+[[name, count, weight]]
            return new
        return self.copy(filt)

    def Sum(self, *args):
        return self.copy(self.df.Sum(*args))

    def Count(self, *args):
        return self.copy(self.df.Count(*args))

    def Histo1D(self, *args):
        return self.copy(self.df.Histo1D(*args))

    def Snapshot(self, *args):
        return self.copy(self.df.Snapshot(*args))

    def SetWeight(self, weight):
        """Set the weight for the event count of current and subsequent
filters (doesn't affect previous filters)
        """
        self.weight = weight
        self.filters[-1][2] = self.df.Sum(weight)
        return self

    def AppendCounts(self, filename):
        "Add a counts and weights histogram to the file of the given filename"
        f = r.TFile(filename, "update")
        counts = [c.GetValue() for _,c,_ in self.filters]
        weights = [w.GetValue() for _,_,w in self.filters]
        names = [n for n,_,_ in self.filters]
        h = r.TH1D("counts", "", len(counts), 0.5, len(counts)+0.5)
        h.SetDirectory(f)
        [h.SetBinContent(i+1, c) for i,c in enumerate(counts)]
        [h.GetXaxis().SetBinLabel(i+1, n) for i,n in enumerate(names)]
        h.Write()
        h = r.TH1D("weights", "", len(weights), 0.5, len(weights)+0.5)
        h.SetDirectory(f)
        [h.SetBinContent(i+1, c) for i,c in enumerate(weights)]
        [h.GetXaxis().SetBinLabel(i+1, n) for i,n in enumerate(names)]
        h.Write()
        f.Close()

    def Report(self):
        filters = [(n, "{}".format(c.GetValue()), "{}".format(w.GetValue())) for n, c, w in self.filters]
        nmaxi = max([len(n) for n, _, _ in filters]+[4])
        cmaxi = max([len(c) for _, c, _ in filters]+[5])
        wmaxi = max([len(w) for _, _, w in filters]+[6])
        print "| {{:{}s}} | {{:{}s}} | {{:{}s}} |".format(nmaxi, cmaxi, wmaxi).format("Name", "Count", "Weight")
        print "| {{:{}s}} | {{:{}s}} | {{:{}s}} |".format(nmaxi, cmaxi, wmaxi).format("-"*nmaxi, "-"*cmaxi, "-"*wmaxi)
        for name, count, weight in filters:
            print "| {{:{}s}} | {{:{}s}} | {{:{}s}} |".format(nmaxi, cmaxi, wmaxi).format(name, count, weight)
