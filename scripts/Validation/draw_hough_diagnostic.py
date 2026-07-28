#!/usr/bin/env python3
"""Draw one rejected 2D Hough attempt from the optional diagnostic trees.

Example:
  python scripts/Validation/draw_hough_diagnostic.py reco.root \
      --event 0 --slice 0 --view X --attempt 0

The plot is in the plane/bar coordinates used by the Hough post-processing
DBSCAN.  Grey points are all cleaned hits in the selected view, orange points
are the Hough seed, and blue points are the retained largest DBSCAN cluster.
"""

import argparse
import os

import ROOT


VIEW_TO_BAR_TYPE = {"X": 0, "Y": 1, "U": 2, "V": 3}
STAGE_NAMES = {
    0: "accepted",
    1: "input empty",
    2: "empty after cleaning",
    3: "below per-view minimum hits",
    4: "no Hough seed hits",
    5: "no DBSCAN cluster",
    6: "A* emptied candidate",
    7: "too few final hits",
    8: "too short",
    9: "extrapolation emptied candidate",
    10: "post-Hough A* emptied candidate",
    11: "event below minimum hits",
}


def matching_entry(tree, event, slice_no, view, attempt):
    for entry_number in range(tree.GetEntries()):
        tree.GetEntry(entry_number)
        if (int(tree.EventNo) == event and int(tree.SliceNo) == slice_no and
                int(tree.View) == ord(view) and int(tree.Attempt) == attempt):
            return entry_number
    return None


def matching_line_entry(tree, event, slice_no, spill_no, run_no):
    for entry_number in range(tree.GetEntries()):
        tree.GetEntry(entry_number)
        if (int(tree.EventNo) == event and int(tree.SliceNo) == slice_no and
                int(tree.SpillNo) == spill_no and int(tree.RunNo) == run_no):
            return entry_number
    return None


def graph(points, color, marker_style, marker_size):
    result = ROOT.TGraph(len(points))
    for index, (plane, bar) in enumerate(points):
        result.SetPoint(index, plane, bar)
    result.SetMarkerColor(color)
    result.SetMarkerStyle(marker_style)
    result.SetMarkerSize(marker_size)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="dune-tms reconstruction ROOT output")
    parser.add_argument("--event", type=int, required=True, help="EventNo")
    parser.add_argument("--slice", dest="slice_no", type=int, required=True, help="SliceNo")
    parser.add_argument("--view", choices=sorted(VIEW_TO_BAR_TYPE), required=True, help="2D view")
    parser.add_argument("--attempt", type=int, default=0, help="Hough attempt number (default: 0)")
    parser.add_argument("--output", help="PNG output path (default: beside input)")
    args = parser.parse_args()

    input_file = ROOT.TFile.Open(args.input)
    if not input_file or input_file.IsZombie():
        raise RuntimeError("Could not open {}".format(args.input))

    diagnostics = input_file.Get("Hough_Diagnostics")
    snapshots = input_file.Get("Hough_Diagnostic_Hits")
    lines = input_file.Get("Line_Candidates")
    if not diagnostics or not snapshots or not lines:
        raise RuntimeError("Input needs Hough_Diagnostics, Hough_Diagnostic_Hits, and Line_Candidates")

    diagnostic_entry = matching_entry(diagnostics, args.event, args.slice_no, args.view, args.attempt)
    if diagnostic_entry is None:
        raise RuntimeError("No matching diagnostic row")
    diagnostics.GetEntry(diagnostic_entry)
    if int(diagnostics.RejectStage) == 0:
        raise RuntimeError("Selected attempt was accepted; it has no rejected-candidate snapshot")

    snapshot_entry = matching_entry(snapshots, args.event, args.slice_no, args.view, args.attempt)
    if snapshot_entry is None:
        raise RuntimeError("No matching hit snapshot; rerun with WriteDiagnosticHitSnapshots = true")
    snapshots.GetEntry(snapshot_entry)

    line_entry = matching_line_entry(lines, args.event, args.slice_no,
                                     int(diagnostics.SpillNo), int(diagnostics.RunNo))
    if line_entry is None:
        raise RuntimeError("Could not find the matching Line_Candidates entry")
    lines.GetEntry(line_entry)

    selected_bar_type = VIEW_TO_BAR_TYPE[args.view]
    all_hits = []
    for index in range(int(lines.nHits)):
        if int(lines.RecoHitBarType[index]) == selected_bar_type:
            all_hits.append((int(lines.RecoHitPlane[index]), int(lines.RecoHitBar[index])))

    def hit_points(indices):
        return [(int(lines.RecoHitPlane[index]), int(lines.RecoHitBar[index])) for index in indices]

    seed_hits = hit_points(list(snapshots.SeedHitIndices))
    post_dbscan_hits = hit_points(list(snapshots.PostDBSCANHitIndices))
    if not all_hits:
        raise RuntimeError("No cleaned hits found in {} view".format(args.view))

    planes = [point[0] for point in all_hits]
    bars = [point[1] for point in all_hits]
    x_padding = max(1.0, 0.05 * (max(planes) - min(planes) + 1))
    y_padding = max(1.0, 0.05 * (max(bars) - min(bars) + 1))
    frame = ROOT.TH2D("frame", "", 10, min(planes) - x_padding, max(planes) + x_padding,
                      10, min(bars) - y_padding, max(bars) + y_padding)
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("TMS plane number")
    frame.GetYaxis().SetTitle("TMS bar number")

    output = args.output
    if not output:
        stem = os.path.splitext(os.path.basename(args.input))[0]
        output = os.path.join(os.path.dirname(os.path.abspath(args.input)),
                              "{}_event{}_slice{}_{}{}_hough_diagnostic.png".format(
                                  stem, args.event, args.slice_no, args.view, args.attempt))

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    canvas = ROOT.TCanvas("hough_diagnostic", "Hough diagnostic", 1100, 850)
    frame.Draw()
    all_graph = graph(all_hits, ROOT.kGray + 1, 20, 0.7)
    seed_graph = graph(seed_hits, ROOT.kOrange + 7, 24, 1.25)
    post_dbscan_graph = graph(post_dbscan_hits, ROOT.kAzure + 2, 21, 1.4)
    all_graph.Draw("P SAME")
    seed_graph.Draw("P SAME")
    post_dbscan_graph.Draw("P SAME")

    stage = int(diagnostics.RejectStage)
    title = "Event {} slice {} {} attempt {}: {}".format(
        args.event, args.slice_no, args.view, args.attempt, STAGE_NAMES.get(stage, "unknown"))
    label = ROOT.TPaveText(0.12, 0.83, 0.62, 0.92, "NDC")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.AddText(title)
    label.AddText("seed={}  largest DBSCAN={}  final={}  endpoint distance={:.2f}".format(
        int(diagnostics.nSeed), int(diagnostics.nLargestDBSCAN),
        int(diagnostics.nFinal), float(diagnostics.EndpointDistance)))
    label.Draw()

    legend = ROOT.TLegend(0.68, 0.77, 0.90, 0.91)
    legend.AddEntry(all_graph, "all cleaned view hits", "p")
    legend.AddEntry(seed_graph, "Hough seed", "p")
    legend.AddEntry(post_dbscan_graph, "post-DBSCAN", "p")
    legend.Draw()
    canvas.SaveAs(output)
    print("Wrote {}".format(output))


if __name__ == "__main__":
    main()
