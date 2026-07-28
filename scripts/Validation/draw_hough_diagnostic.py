#!/usr/bin/env python3
"""Render rejected 2D Hough attempts from the optional diagnostic trees.

Example:
  python scripts/Validation/draw_hough_diagnostic.py reco.root \
      --output-dir hough_diagnostic_plots

Optional --event, --slice, --view, and --attempt filters select a subset.

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


def draw_attempt(lines, diagnostic, seed_indices, post_dbscan_indices, output):
    line_entry = matching_line_entry(lines, diagnostic["event"], diagnostic["slice"],
                                     diagnostic["spill"], diagnostic["run"])
    if line_entry is None:
        print("Skipping event {} slice {}: no Line_Candidates entry".format(
            diagnostic["event"], diagnostic["slice"]))
        return False
    lines.GetEntry(line_entry)

    selected_bar_type = VIEW_TO_BAR_TYPE[diagnostic["view"]]
    all_hits = []
    for index in range(int(lines.nHits)):
        if int(lines.RecoHitBarType[index]) == selected_bar_type:
            all_hits.append((int(lines.RecoHitPlane[index]), int(lines.RecoHitBar[index])))

    def hit_points(indices):
        return [(int(lines.RecoHitPlane[index]), int(lines.RecoHitBar[index])) for index in indices]

    seed_hits = hit_points(seed_indices)
    post_dbscan_hits = hit_points(post_dbscan_indices)
    if not all_hits:
        print("Skipping event {} slice {} {}{}: no cleaned view hits".format(
            diagnostic["event"], diagnostic["slice"], diagnostic["view"], diagnostic["attempt"]))
        return False

    planes = [point[0] for point in all_hits]
    bars = [point[1] for point in all_hits]
    x_padding = max(1.0, 0.05 * (max(planes) - min(planes) + 1))
    y_padding = max(1.0, 0.05 * (max(bars) - min(bars) + 1))
    frame = ROOT.TH2D("frame", "", 10, min(planes) - x_padding, max(planes) + x_padding,
                      10, min(bars) - y_padding, max(bars) + y_padding)
    frame.SetDirectory(0)
    frame.GetXaxis().SetTitle("TMS plane number")
    frame.GetYaxis().SetTitle("TMS bar number")

    canvas = ROOT.TCanvas("hough_diagnostic", "Hough diagnostic", 1100, 850)
    frame.Draw()
    all_graph = graph(all_hits, ROOT.kGray + 1, 20, 0.7)
    seed_graph = graph(seed_hits, ROOT.kOrange + 7, 24, 1.25)
    post_dbscan_graph = graph(post_dbscan_hits, ROOT.kAzure + 2, 21, 1.4)
    all_graph.Draw("P SAME")
    seed_graph.Draw("P SAME")
    post_dbscan_graph.Draw("P SAME")

    stage = diagnostic["stage"]
    title = "Event {} slice {} {} attempt {}: {}".format(
        diagnostic["event"], diagnostic["slice"], diagnostic["view"], diagnostic["attempt"],
        STAGE_NAMES.get(stage, "unknown"))
    label = ROOT.TPaveText(0.12, 0.83, 0.62, 0.92, "NDC")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.AddText(title)
    label.AddText("seed={}  largest DBSCAN={}  final={}  endpoint distance={:.2f}".format(
        diagnostic["n_seed"], diagnostic["largest_dbscan"],
        diagnostic["n_final"], diagnostic["endpoint_distance"]))
    label.Draw()

    legend = ROOT.TLegend(0.68, 0.77, 0.90, 0.91)
    legend.AddEntry(all_graph, "all cleaned view hits", "p")
    legend.AddEntry(seed_graph, "Hough seed", "p")
    legend.AddEntry(post_dbscan_graph, "post-DBSCAN", "p")
    legend.Draw()
    canvas.SaveAs(output)
    print("Wrote {}".format(output))
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="dune-tms reconstruction ROOT output")
    parser.add_argument("--event", type=int, help="only this EventNo")
    parser.add_argument("--slice", dest="slice_no", type=int, help="only this SliceNo")
    parser.add_argument("--view", choices=sorted(VIEW_TO_BAR_TYPE), help="only this 2D view")
    parser.add_argument("--attempt", type=int, help="only this Hough attempt number")
    parser.add_argument("--output-dir", help="directory for PNGs (default: beside input)")
    args = parser.parse_args()

    input_file = ROOT.TFile.Open(args.input)
    if not input_file or input_file.IsZombie():
        raise RuntimeError("Could not open {}".format(args.input))

    diagnostics = input_file.Get("Hough_Diagnostics")
    snapshots = input_file.Get("Hough_Diagnostic_Hits")
    lines = input_file.Get("Line_Candidates")
    if not diagnostics or not snapshots or not lines:
        raise RuntimeError("Input needs Hough_Diagnostics, Hough_Diagnostic_Hits, and Line_Candidates")

    diagnostic_rows = {}
    for entry_number in range(diagnostics.GetEntries()):
        diagnostics.GetEntry(entry_number)
        key = (int(diagnostics.EventNo), int(diagnostics.SliceNo), int(diagnostics.View), int(diagnostics.Attempt))
        diagnostic_rows[key] = {
            "event": int(diagnostics.EventNo),
            "slice": int(diagnostics.SliceNo),
            "spill": int(diagnostics.SpillNo),
            "run": int(diagnostics.RunNo),
            "view": chr(int(diagnostics.View)),
            "attempt": int(diagnostics.Attempt),
            "stage": int(diagnostics.RejectStage),
            "n_seed": int(diagnostics.nSeed),
            "largest_dbscan": int(diagnostics.nLargestDBSCAN),
            "n_final": int(diagnostics.nFinal),
            "endpoint_distance": float(diagnostics.EndpointDistance),
        }

    stem = os.path.splitext(os.path.basename(args.input))[0]
    output_dir = args.output_dir or os.path.join(os.path.dirname(os.path.abspath(args.input)),
                                                  stem + "_hough_diagnostic_plots")
    os.makedirs(output_dir, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    rendered = 0
    selected = 0
    for entry_number in range(snapshots.GetEntries()):
        snapshots.GetEntry(entry_number)
        event = int(snapshots.EventNo)
        slice_no = int(snapshots.SliceNo)
        view = chr(int(snapshots.View))
        attempt = int(snapshots.Attempt)
        if ((args.event is not None and event != args.event) or
                (args.slice_no is not None and slice_no != args.slice_no) or
                (args.view is not None and view != args.view) or
                (args.attempt is not None and attempt != args.attempt)):
            continue
        selected += 1
        if view not in VIEW_TO_BAR_TYPE:
            print("Skipping event {} slice {} view {}: no plane/bar projection".format(
                event, slice_no, view))
            continue
        key = (event, slice_no, int(snapshots.View), attempt)
        diagnostic = diagnostic_rows.get(key)
        if not diagnostic:
            print("Skipping event {} slice {} {}{}: no diagnostic row".format(
                event, slice_no, view, attempt))
            continue
        output = os.path.join(output_dir, "event{:06d}_slice{:04d}_{}{:02d}_stage{:02d}.png".format(
            event, slice_no, view, attempt, diagnostic["stage"]))
        if draw_attempt(lines, diagnostic, list(snapshots.SeedHitIndices),
                        list(snapshots.PostDBSCANHitIndices), output):
            rendered += 1
    print("Rendered {} of {} selected rejected attempts into {}".format(rendered, selected, output_dir))


if __name__ == "__main__":
    main()
