# Compare every numeric branch of every tree in two ROOT files, value by value (exact).
# Usage: compare_trees.py a.root b.root
import ROOT, sys, numpy as np
ROOT.gROOT.SetBatch(True)
fa, fb = ROOT.TFile.Open(sys.argv[1]), ROOT.TFile.Open(sys.argv[2])
names = sorted(k.GetName() for k in fa.GetListOfKeys() if k.GetClassName() == "TTree")
assert names == sorted(k.GetName() for k in fb.GetListOfKeys() if k.GetClassName() == "TTree"), "tree lists differ"
bad = 0
def column(tree, leaf):
    n = tree.Draw(leaf, "", "goff")
    if n < 0: return None
    if n == 0: return np.zeros(0)
    return np.array(tree.GetV1().reshape((n,)), dtype=np.float64, copy=True)
for name in names:
    ta, tb = fa.Get(name), fb.Get(name)
    if ta.GetEntries() != tb.GetEntries():
        print(f"{name}: entries differ {ta.GetEntries()} vs {tb.GetEntries()}"); bad += 1; continue
    leaves = [l.GetName() for l in ta.GetListOfLeaves()]
    nchecked = nskipped = ndiff = 0
    for leaf in leaves:
        a, b = column(ta, leaf), column(tb, leaf)
        if a is None or b is None: nskipped += 1; continue
        nchecked += 1
        if a.shape != b.shape or not np.array_equal(a, b, equal_nan=True):
            ndiff += 1; bad += 1
            extra = "" if a.shape != b.shape else f", {np.sum(~((a == b) | (np.isnan(a) & np.isnan(b))))} of {a.size} values"
            print(f"  DIFF {name}.{leaf}: sizes {a.size} vs {b.size}{extra}")
    print(f"{name}: {ta.GetEntries()} entries, {nchecked} branches compared, {nskipped} skipped (non-numeric), {ndiff} differ")
print("IDENTICAL" if bad == 0 else f"DIFFERENCES: {bad}")
