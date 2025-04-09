import sys

import ROOT
ROOT.gROOT.SetBatch(True)


if len(sys.argv) < 2:
    print("Usage: python print_event_rates.py filename.py")
    exit(0)

filename = sys.argv[1]
tf = ROOT.TFile(filename)
assert not tf.IsZombie(), f"Had issue opening file {filename}"
hists = ["event_rates__tms__tms_event_rates_with_shell_cut", "event_rates__tms__tms_event_rates",
         "event_rates__tms__tms_event_rates_lar_active", "event_rates__tms__tms_event_rates_by_vertex",
         "event_rates__tms__tms_muon_reco_rates_with_shell_cut", "event_rates__tms__tms_muon_reco_rates",
         "event_rates__tms__tms_muon_reco_rates_lar_active_region", "event_rates__reco_pdg__pdg_with_shell_cut",
         "event_rates__reco_pdg__pdg", "event_rates__reco_pdg__pdg_lar_active_region"]
for name in hists:
    h = tf.Get(name)
    assert h != None, f"Couldn't find hist named {name} in {filename}"
    #print(h.GetXaxis().GetModifiedLabels())
    #print(h.GetXaxis().GetModifiedLabels()[0].GetText())
    #print(name, h.GetXaxis().GetNlabels())
    title = h.GetTitle()
    second_title = name.rsplit("__", 1)[-1]
    print(f"\n{title} ({second_title}):")
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        I = i-1 # get index starting from zero
        label = f"bin {i}"
        if h.GetXaxis().GetNlabels() > I: label = h.GetXaxis().GetModifiedLabels()[I].GetText()
        print(f"{label}:\t{c:g}")
