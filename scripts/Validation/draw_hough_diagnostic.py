#!/usr/bin/env python3
"""Render 2D Hough attempts from the optional diagnostic trees.

Example:
  python scripts/Validation/draw_hough_diagnostic.py reco.root \
      --output-dir hough_diagnostic_plots

Optional --event, --slice, --view, and --attempt filters select a subset.
Use --failed-only to omit accepted candidates.

Each image has a physical Z/not-Z Hough panel with the fitted line and a
plane/bar panel in the coordinates used by post-Hough DBSCAN. Grey points are
all cleaned hits in the selected view; orange points are the Hough seed, green
points are hits added by walking, and blue points are the retained largest
DBSCAN cluster.
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


def snapshot_points(plane_bar, z_not_z):
    if len(plane_bar) != len(z_not_z) or len(plane_bar) % 2:
        raise RuntimeError("Malformed diagnostic hit snapshot")
    plane_bar_points = []
    physical_points = []
    for index in range(0, len(plane_bar), 2):
        plane_bar_points.append((int(plane_bar[index]), int(plane_bar[index + 1])))
        physical_points.append((float(z_not_z[index]) / 1000.0, float(z_not_z[index + 1]) / 1000.0))
    return plane_bar_points, physical_points


def draw_attempt(lines, diagnostic, seed_snapshot, walked_snapshot, post_dbscan_snapshot, output):
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
            not_z_component = 1 if selected_bar_type == 0 else 0
            all_hits.append({
                "plane": int(lines.RecoHitPlane[index]),
                "bar": int(lines.RecoHitBar[index]),
                "z": float(lines.RecoHitPos[index * 4 + 2]),
                "not_z": float(lines.RecoHitPos[index * 4 + not_z_component]),
            })

    seed_hits, seed_physical_hits = snapshot_points(*seed_snapshot)
    walked_all_hits, walked_all_physical_hits = snapshot_points(*walked_snapshot)
    post_dbscan_hits, post_dbscan_physical_hits = snapshot_points(*post_dbscan_snapshot)
    seed_keys = set(zip(seed_hits, seed_physical_hits))
    walked_pairs = list(zip(walked_all_hits, walked_all_physical_hits))
    walked_added_pairs = [pair for pair in walked_pairs if pair not in seed_keys]
    walked_hits = [pair[0] for pair in walked_added_pairs]
    walked_physical_hits = [pair[1] for pair in walked_added_pairs]
    all_plane_bar_hits = [(hit["plane"], hit["bar"]) for hit in all_hits]
    all_physical_hits = [(hit["z"] / 1000.0, hit["not_z"] / 1000.0) for hit in all_hits]
    if not all_hits:
        print("Skipping event {} slice {} {}{}: no cleaned view hits".format(
            diagnostic["event"], diagnostic["slice"], diagnostic["view"], diagnostic["attempt"]))
        return False

    z_values = [point[0] for point in all_physical_hits]
    not_z_values = [point[1] for point in all_physical_hits]
    z_padding = max(0.1, 0.05 * (max(z_values) - min(z_values) + 0.1))
    not_z_padding = max(0.1, 0.05 * (max(not_z_values) - min(not_z_values) + 0.1))

    object_name = "hough_event{}_slice{}_{}{}".format(
        diagnostic["event"], diagnostic["slice"], diagnostic["view"], diagnostic["attempt"])
    canvas = ROOT.TCanvas(object_name, "Hough diagnostic", 1800, 1000)
    plot_pad = ROOT.TPad(object_name + "_plot", "", 0.0, 0.0, 0.68, 1.0)
    info_pad = ROOT.TPad(object_name + "_info", "", 0.68, 0.0, 1.0, 1.0)
    plot_pad.SetLeftMargin(0.11)
    plot_pad.SetRightMargin(0.03)
    plot_pad.SetBottomMargin(0.10)
    info_pad.SetMargin(0.04, 0.04, 0.04, 0.04)
    plot_pad.Draw()
    info_pad.Draw()

    plot_pad.cd()
    physical_frame = ROOT.TH2D(object_name + "_physical", "", 10, min(z_values) - z_padding, max(z_values) + z_padding,
                               10, min(not_z_values) - not_z_padding, max(not_z_values) + not_z_padding)
    physical_frame.SetDirectory(0)
    physical_frame.GetXaxis().SetTitle("Z (m)")
    physical_frame.GetYaxis().SetTitle("not-Z coordinate (m)")
    physical_frame.Draw()
    all_physical_graph = graph(all_physical_hits, ROOT.kGray + 1, 20, 0.65)
    seed_physical_graph = graph(seed_physical_hits, ROOT.kOrange + 7, 24, 1.1)
    walked_physical_graph = graph(walked_physical_hits, ROOT.kGreen + 2, 25, 1.0)
    post_dbscan_physical_graph = graph(post_dbscan_physical_hits, ROOT.kAzure + 2, 21, 1.3)
    line_graph = ROOT.TGraph(2)
    line_graph.SetPoint(0, min(z_values) - z_padding,
                        (diagnostic["slope"] * (min(z_values) - z_padding) * 1000.0 + diagnostic["intercept"]) / 1000.0)
    line_graph.SetPoint(1, max(z_values) + z_padding,
                        (diagnostic["slope"] * (max(z_values) + z_padding) * 1000.0 + diagnostic["intercept"]) / 1000.0)
    line_graph.SetLineColor(ROOT.kRed + 1)
    line_graph.SetLineStyle(2)
    line_graph.SetLineWidth(2)
    all_physical_graph.Draw("P SAME")
    line_graph.Draw("L SAME")
    seed_physical_graph.Draw("P SAME")
    walked_physical_graph.Draw("P SAME")
    post_dbscan_physical_graph.Draw("P SAME")

    info_pad.cd()
    stage = diagnostic["stage"]
    title = "Event {} slice {} {} attempt {}: {}".format(
        diagnostic["event"], diagnostic["slice"], diagnostic["view"], diagnostic["attempt"],
        STAGE_NAMES.get(stage, "unknown"))
    label = ROOT.TPaveText(0.06, 0.48, 0.94, 0.94, "NDC")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.SetTextAlign(12)
    label.SetTextSize(0.037)
    label.AddText(title)
    label.AddText("input: {}    after cleaning: {}    projected: {}".format(
        diagnostic["n_input"], diagnostic["n_after_clean"], diagnostic["n_projected"]))
    label.AddText("Hough seed: {}    after walking: {}".format(
        diagnostic["n_seed"], diagnostic["n_after_walk"]))
    label.AddText("DBSCAN: {} clusters, largest {}, retained {}".format(
        diagnostic["n_dbscan_clusters"], diagnostic["largest_dbscan"], diagnostic["n_after_dbscan"]))
    label.AddText("after A*: {}    after extrapolation: {}    final: {}".format(
        diagnostic["n_after_astar"], diagnostic["n_after_extrapolation"], diagnostic["n_final"]))
    label.AddText("endpoint distance={:.2f}".format(diagnostic["endpoint_distance"]))
    label.Draw()

    legend = ROOT.TLegend(0.08, 0.12, 0.92, 0.42)
    legend.AddEntry(all_physical_graph, "all cleaned view hits", "p")
    legend.AddEntry(line_graph, "Hough line", "l")
    legend.AddEntry(seed_physical_graph, "Hough seed", "p")
    legend.AddEntry(walked_physical_graph, "added by walking", "p")
    legend.AddEntry(post_dbscan_physical_graph, "post-DBSCAN", "p")
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
    parser.add_argument("--failed-only", action="store_true", help="render rejected attempts only")
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
    if not snapshots.GetBranch("SeedPlaneBar"):
        raise RuntimeError("Input has an older diagnostic snapshot schema; rerun reconstruction")

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
            "n_input": int(diagnostics.nInput),
            "n_after_clean": int(diagnostics.nAfterClean),
            "n_projected": int(diagnostics.nProjected),
            "n_seed": int(diagnostics.nSeed),
            "n_after_walk": int(diagnostics.nAfterWalk),
            "n_dbscan_clusters": int(diagnostics.nDBSCANClusters),
            "largest_dbscan": int(diagnostics.nLargestDBSCAN),
            "n_after_dbscan": int(diagnostics.nAfterDBSCAN),
            "n_after_astar": int(diagnostics.nAfterAStar),
            "n_after_extrapolation": int(diagnostics.nAfterExtrapolation),
            "n_final": int(diagnostics.nFinal),
            "slope": float(diagnostics.Slope),
            "intercept": float(diagnostics.Intercept),
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
        if args.failed_only and diagnostic["stage"] == 0:
            continue
        selected += 1
        output = os.path.join(output_dir, "event{:06d}_slice{:04d}_{}{:02d}_stage{:02d}.png".format(
            event, slice_no, view, attempt, diagnostic["stage"]))
        seed_snapshot = (list(snapshots.SeedPlaneBar), list(snapshots.SeedZNotZ))
        walked_snapshot = (list(snapshots.WalkedPlaneBar), list(snapshots.WalkedZNotZ))
        post_dbscan_snapshot = (list(snapshots.PostDBSCANPlaneBar), list(snapshots.PostDBSCANZNotZ))
        if draw_attempt(lines, diagnostic, seed_snapshot, walked_snapshot, post_dbscan_snapshot, output):
            rendered += 1
    print("Rendered {} of {} selected attempts into {}".format(rendered, selected, output_dir))


if __name__ == "__main__":
    main()
