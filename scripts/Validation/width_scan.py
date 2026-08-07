#!/usr/bin/env python3
"""Run a reproducible transverse-width reconstruction scan.

Each configuration is run over the same input ROOT files.  The script keeps
the produced ROOT files, generated TOMLs, machine-readable summary, and a
compact CSV together under one output directory.

The reported stage efficiencies are reconstruction-defined: their denominator
is the union of hit-truth primary IDs seen in any reconstructed object.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


DEFAULT_POINT = ("corrected_nominal", 1.5, 1.0, False)
PR289_SUITE = [
    ("legacy_ab", 1.5, 1.0, True),
    ("corrected_nominal", 1.5, 1.0, False),
    ("hough_4x", 6.0, 1.0, False),
    ("extrapolation_4x", 1.5, 4.0, False),
    ("retuned_4x", 6.0, 4.0, False),
]
STAGE_NAMES = {
    0: "accepted", 1: "input_empty", 2: "post_clean_empty",
    3: "below_minimum_hits", 4: "seed_empty", 5: "no_dbscan_cluster",
    6: "astar_empty", 7: "final_too_few_hits", 8: "final_too_short",
    9: "post_extrapolation_empty", 10: "post_hough_astar_empty",
    11: "event_below_minimum_hits",
}


@dataclass
class Point:
    name: str
    hough_half_width: float
    extrapolation_multiplier: float
    legacy_x: bool


def parse_point(value):
    fields = value.split(",")
    if len(fields) not in (3, 4):
        raise argparse.ArgumentTypeError(
            "point must be NAME,HOUGH_HALF_WIDTH,EXTRAPOLATION_MULTIPLIER[,LEGACY_X]"
        )
    try:
        point = Point(fields[0], float(fields[1]), float(fields[2]),
                      len(fields) == 4 and fields[3].lower() in ("1", "true", "yes"))
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error))
    if not point.name or min(point.hough_half_width, point.extrapolation_multiplier) <= 0:
        raise argparse.ArgumentTypeError("point values must all be positive")
    return point


def replace_toml_value(text, key, value):
    pattern = re.compile(r"^(\s*" + re.escape(key) + r"\s*=\s*)[^#\n]*(.*)$", re.MULTILINE)
    text, replacements = pattern.subn(r"\g<1>" + str(value) + r"\2", text)
    if replacements != 1:
        raise RuntimeError("expected exactly one TOML setting named " + key)
    return text


def write_config(template, destination, point, max_events, diagnostics):
    text = template.read_text()
    values = {
        "ContainmentHalfWidth": point.hough_half_width,
        "ContainmentWidthMultiplier": point.extrapolation_multiplier,
        "UseLegacyXBarContainment": str(point.legacy_x).lower(),
        "WriteDiagnostics": str(diagnostics).lower(),
        "WriteDiagnosticHitSnapshots": "false",
    }
    if max_events is not None:
        values["MaximumNEvents"] = max_events
    for key, value in values.items():
        text = replace_toml_value(text, key, value)
    destination.write_text(text)


def get_required_tree(root_file, name):
    tree = root_file.Get(name)
    if not tree:
        raise RuntimeError("missing " + name + " tree")
    return tree


def get_branch(tree, name):
    if not tree.GetBranch(name):
        raise RuntimeError("missing " + name + " branch in " + tree.GetName())


def fixed_array_value(entry, branch_name, row, column, row_width):
    """Read a fixed-size 2-D ROOT leaf across PyROOT array representations."""
    values = getattr(entry, branch_name)
    try:
        return values[row][column]
    except TypeError:
        # Some PyROOT versions expose TrackHitBarType[nTracks][200] as one
        # flat buffer rather than an array of 200-element rows.
        return values[row * row_width + column]


def entry_key(entry):
    """Identity shared by the Line_Candidates and Reco_Tree slice entries."""
    return tuple(int(getattr(entry, name))
                 for name in ("RunNo", "SpillNo", "EventNo", "SliceNo"))


def global_particle_key(key):
    """Identity of a hit-truth primary, independent of reconstruction slice."""
    return key[0], key[1], key[2], key[4], key[5]


def candidate_primary(line_data, view, index):
    """Return the leading truth-particle identity for one 2-D candidate."""
    candidates = line_data[view]
    if index < 0 or index >= len(candidates):
        return None
    return candidates[index]


def add_view_truth_metrics(lines, summary, view, source_id):
    """Summarise completeness, cleanliness, and primary-particle duplication."""
    suffix = view
    count_branch = "nLines" + suffix
    branches = (
        count_branch, "PrimaryVertexId" + suffix, "PrimaryTrackId" + suffix,
        "PrimaryTrueVisibleEnergy" + suffix,
        "PrimaryTrueVisibleEnergyInView" + suffix,
        "TotalTrueVisibleEnergy" + suffix,
    )
    for branch in branches:
        get_branch(lines, branch)

    prefix = "candidate_" + view.lower() + "_"
    completeness = []
    cleanliness = []
    primary_counts = {}
    for entry in lines:
        for index in range(int(getattr(entry, count_branch))):
            vertex = int(getattr(entry, "PrimaryVertexId" + suffix)[index])
            track = int(getattr(entry, "PrimaryTrackId" + suffix)[index])
            primary_energy = float(getattr(entry, "PrimaryTrueVisibleEnergy" + suffix)[index])
            view_energy = float(getattr(entry, "PrimaryTrueVisibleEnergyInView" + suffix)[index])
            total_energy = float(getattr(entry, "TotalTrueVisibleEnergy" + suffix)[index])
            if vertex < 0 or track < 0 or primary_energy < 0:
                continue
            if view_energy > 0:
                completeness.append(primary_energy / view_energy)
            if total_energy > 0:
                cleanliness.append(primary_energy / total_energy)
            primary_counts[(source_id, int(entry.RunNo), int(entry.SpillNo), int(entry.EventNo),
                            int(entry.SliceNo), vertex, track)] = (
                primary_counts.get((source_id, int(entry.RunNo), int(entry.SpillNo), int(entry.EventNo),
                                    int(entry.SliceNo), vertex, track), 0) + 1
            )

    multiplicities = list(primary_counts.values())
    summary[prefix + "count"] = len(completeness)
    summary[prefix + "completeness_sum"] = sum(completeness)
    summary[prefix + "cleanliness_sum"] = sum(cleanliness)
    summary[prefix + "unique_primaries"] = len(multiplicities)
    summary[prefix + "multi_reco_primaries"] = sum(count > 1 for count in multiplicities)
    summary[prefix + "multiplicity_mean"] = (
        sum(multiplicities) / len(multiplicities) if multiplicities else None
    )
    summary[prefix + "from_multi_reco_fraction"] = (
        sum(count for count in multiplicities if count > 1) / sum(multiplicities)
        if multiplicities else None
    )


def add_stage_metrics(summary, stage, view, candidates, truth_energy):
    """Record object truth summaries for one reconstruction stage/view.

    The scan's denominator is intentionally not a Truth_Spill selection.  It
    is built later from the union of primary IDs observed in *any* recorded
    reconstruction object.  Keep the per-object primary-energy information
    here, and make that union only after all stages have been seen.
    """
    prefix = "stage_" + stage + "_" + view.lower() + "_"
    counts = summary.setdefault("_stage_primary_counts", {}).setdefault(
        stage + ":" + view.upper(), Counter())
    completeness = []
    cleanliness = []
    for key, primary_energy, total_energy in candidates:
        if primary_energy <= 0:
            continue
        counts[global_particle_key(key)] += 1
        denominator = truth_energy.get(key, 0.0)
        if denominator > 0:
            completeness.append(primary_energy / denominator)
        if total_energy > 0:
            cleanliness.append(primary_energy / total_energy)
    summary[prefix + "candidates"] = sum(counts.values())
    summary[prefix + "completeness_sum"] = sum(completeness)
    summary[prefix + "cleanliness_sum"] = sum(cleanliness)
    summary[prefix + "quality_count"] = len(completeness)


def derive_stage_metrics(root_file, lines, summary):
    truth = root_file.Get("Hough_View_Truth")
    diagnostics = root_file.Get("Hough_Diagnostics")
    if not truth or not diagnostics:
        return
    for tree, names in ((truth, ("EventNo", "SliceNo", "SpillNo", "RunNo", "View", "VertexId", "TrackId", "VisibleEnergy")),
                        (diagnostics, ("EventNo", "SliceNo", "SpillNo", "RunNo", "View", "PrimaryVertexId", "PrimaryTrackId", "PrimaryEnergy", "TotalEnergy"))):
        for name in names:
            get_branch(tree, name)
    energies = {"X": {}, "Y": {}, "XY": {}}
    for entry in truth:
        view = "XY" if int(entry.View) == 3 else chr(int(entry.View))
        if view in energies:
            key = (int(entry.RunNo), int(entry.SpillNo), int(entry.EventNo), int(entry.SliceNo),
                   int(entry.VertexId), int(entry.TrackId))
            energies[view][key] = float(entry.VisibleEnergy)
    stage_indices = (("seed", 0), ("dbscan", 1), ("hough_final", 2))
    stage_candidates = {(stage, view): [] for stage, _ in stage_indices for view in ("X", "Y")}
    for entry in diagnostics:
        view = chr(int(entry.View))
        if view not in energies:
            continue
        base = (int(entry.RunNo), int(entry.SpillNo), int(entry.EventNo), int(entry.SliceNo))
        for stage, index in stage_indices:
            vertex = int(entry.PrimaryVertexId[index])
            track = int(entry.PrimaryTrackId[index])
            if vertex >= 0 and track >= 0:
                stage_candidates[(stage, view)].append(
                    (base + (vertex, track), float(entry.PrimaryEnergy[index]), float(entry.TotalEnergy[index])))
    for stage, _ in stage_indices:
        for view in ("X", "Y"):
            add_stage_metrics(summary, stage, view, stage_candidates[(stage, view)], energies[view])
    final_candidates = {"X": [], "Y": []}
    for entry in lines:
        base = entry_key(entry)
        for view in final_candidates:
            for index in range(int(getattr(entry, "nLines" + view))):
                vertex = int(getattr(entry, "PrimaryVertexId" + view)[index])
                track = int(getattr(entry, "PrimaryTrackId" + view)[index])
                if vertex >= 0 and track >= 0:
                    final_candidates[view].append((base + (vertex, track),
                        float(getattr(entry, "PrimaryTrueVisibleEnergy" + view)[index]),
                        float(getattr(entry, "TotalTrueVisibleEnergy" + view)[index])))
    for view in energies:
        if view in final_candidates:
            add_stage_metrics(summary, "final_2d", view, final_candidates[view], energies[view])
    reco = get_required_tree(root_file, "Reco_Tree")
    for branch in ("nTracks", "PrimaryVertexId", "PrimaryTrackId", "PrimaryVisibleEnergy", "TotalVisibleEnergy"):
        get_branch(reco, branch)
    final_tracks = []
    for entry in reco:
        base = entry_key(entry)
        for index in range(int(entry.nTracks)):
            vertex = int(entry.PrimaryVertexId[index])
            track = int(entry.PrimaryTrackId[index])
            if vertex >= 0 and track >= 0:
                final_tracks.append((base + (vertex, track), float(entry.PrimaryVisibleEnergy[index]), float(entry.TotalVisibleEnergy[index])))
    add_stage_metrics(summary, "final_3d", "XY", final_tracks, energies["XY"])


def summarise_output(path):
    try:
        import ROOT
    except ImportError as error:
        raise RuntimeError("PyROOT is required to summarise outputs; set up ROOT first") from error

    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("could not open output ROOT file " + str(path))

    lines = get_required_tree(root_file, "Line_Candidates")
    reco = get_required_tree(root_file, "Reco_Tree")
    for branch in ("nLinesX", "nLinesY", "nHitsInTrackX", "nHitsInTrackY"):
        get_branch(lines, branch)
    for branch in ("RunNo", "SpillNo", "EventNo", "SliceNo",
                   "PrimaryVertexIdX", "PrimaryTrackIdX",
                   "PrimaryVertexIdY", "PrimaryTrackIdY"):
        get_branch(lines, branch)
    for branch in ("nTracks", "nHits", "TrackHitBarType", "MatchedCandidateU",
                   "MatchedCandidateV", "MatchedCandidateX", "MatchedCandidateY"):
        get_branch(reco, branch)

    summary = {
        "line_entries": int(lines.GetEntries()),
        "lines_x": 0,
        "lines_y": 0,
        "line_hits_x": 0,
        "line_hits_y": 0,
        "slices_with_x_line": 0,
        "slices_with_y_line": 0,
        "reco_slices": int(reco.GetEntries()),
        "slices_with_tracks": 0,
        "reco_tracks": 0,
        "reco_track_hits": 0,
        "reco_track_hits_x": 0,
        "reco_track_hits_y": 0,
        "tracks_with_x_hits": 0,
        "tracks_with_y_hits": 0,
        "tracks_with_xy_hits": 0,
        "reco_tracks_from_uvx_match": 0,
        "reco_tracks_from_xy_match": 0,
        "xy_source_pairs_with_truth": 0,
        "xy_source_pairs_same_primary": 0,
        "xy_source_pairs_different_primary": 0,
        "xy_source_pairs_missing_truth": 0,
        "x_candidates_used_in_xy_match": 0,
        "y_candidates_used_in_xy_match": 0,
        "reco_tracks_without_x_candidate": 0,
        "reco_tracks_without_y_candidate": 0,
    }

    line_entries = {}
    for entry in lines:
        n_x = int(entry.nLinesX)
        n_y = int(entry.nLinesY)
        line_entries[entry_key(entry)] = {
            view: [
                (int(getattr(entry, "PrimaryVertexId" + view)[index]),
                 int(getattr(entry, "PrimaryTrackId" + view)[index]))
                if int(getattr(entry, "PrimaryVertexId" + view)[index]) >= 0
                and int(getattr(entry, "PrimaryTrackId" + view)[index]) >= 0 else None
                for index in range(int(getattr(entry, "nLines" + view)))
            ]
            for view in ("X", "Y")
        }
        summary["lines_x"] += n_x
        summary["lines_y"] += n_y
        summary["line_hits_x"] += sum(int(entry.nHitsInTrackX[index]) for index in range(n_x))
        summary["line_hits_y"] += sum(int(entry.nHitsInTrackY[index]) for index in range(n_y))
        summary["slices_with_x_line"] += n_x > 0
        summary["slices_with_y_line"] += n_y > 0

    add_view_truth_metrics(lines, summary, "X", str(path))
    add_view_truth_metrics(lines, summary, "Y", str(path))
    derive_stage_metrics(root_file, lines, summary)

    used_x_candidates = set()
    used_y_candidates = set()
    for entry in reco:
        key = entry_key(entry)
        line_entry = line_entries.get(key)
        if line_entry is None:
            raise RuntimeError("Reco_Tree entry has no matching Line_Candidates entry: "
                               + repr(key))
        n_tracks = int(entry.nTracks)
        summary["reco_tracks"] += n_tracks
        summary["slices_with_tracks"] += n_tracks > 0
        for track in range(n_tracks):
            n_hits = int(entry.nHits[track])
            summary["reco_track_hits"] += n_hits
            has_x = False
            has_y = False
            for hit in range(n_hits):
                bar_type = int(fixed_array_value(entry, "TrackHitBarType", track, hit, 200))
                if bar_type == 0:  # TMS_Bar::kXBar
                    summary["reco_track_hits_x"] += 1
                    has_x = True
                elif bar_type == 1:  # TMS_Bar::kYBar
                    summary["reco_track_hits_y"] += 1
                    has_y = True
            summary["tracks_with_x_hits"] += has_x
            summary["tracks_with_y_hits"] += has_y
            summary["tracks_with_xy_hits"] += has_x and has_y
            source_u = int(entry.MatchedCandidateU[track]) >= 0
            source_v = int(entry.MatchedCandidateV[track]) >= 0
            source_x = int(entry.MatchedCandidateX[track]) >= 0
            source_y = int(entry.MatchedCandidateY[track]) >= 0
            summary["reco_tracks_from_uvx_match"] += source_u and source_v and source_x
            summary["reco_tracks_from_xy_match"] += source_x and source_y
            if source_x and source_y:
                x_index = int(entry.MatchedCandidateX[track])
                y_index = int(entry.MatchedCandidateY[track])
                used_x_candidates.add((key, x_index))
                used_y_candidates.add((key, y_index))
                x_primary = candidate_primary(line_entry, "X", x_index)
                y_primary = candidate_primary(line_entry, "Y", y_index)
                if x_primary is None or y_primary is None:
                    summary["xy_source_pairs_missing_truth"] += 1
                else:
                    summary["xy_source_pairs_with_truth"] += 1
                    if x_primary == y_primary:
                        summary["xy_source_pairs_same_primary"] += 1
                    else:
                        summary["xy_source_pairs_different_primary"] += 1
            summary["reco_tracks_without_x_candidate"] += not source_x
            summary["reco_tracks_without_y_candidate"] += not source_y

    summary["x_candidates_used_in_xy_match"] = len(used_x_candidates)
    summary["y_candidates_used_in_xy_match"] = len(used_y_candidates)

    diagnostics = root_file.Get("Hough_Diagnostics")
    if diagnostics:
        for branch in ("View", "RejectStage", "nSeed", "nAfterDBSCAN", "nFinal"):
            get_branch(diagnostics, branch)
        summary["hough_attempts"] = int(diagnostics.GetEntries())
        summary["hough_attempts_x"] = 0
        summary["hough_accepted_x"] = 0
        summary["hough_seed_hits_x"] = 0
        summary["hough_post_dbscan_hits_x"] = 0
        summary["hough_final_hits_x"] = 0
        for stage_name in STAGE_NAMES.values():
            summary["hough_x_" + stage_name] = 0
        for entry in diagnostics:
            if chr(int(entry.View)) != "X":
                continue
            summary["hough_attempts_x"] += 1
            stage = int(entry.RejectStage)
            stage_name = STAGE_NAMES.get(stage, "unknown_" + str(stage))
            summary.setdefault("hough_x_" + stage_name, 0)
            summary["hough_x_" + stage_name] += 1
            summary["hough_accepted_x"] += stage == 0
            summary["hough_seed_hits_x"] += int(entry.nSeed)
            summary["hough_post_dbscan_hits_x"] += int(entry.nAfterDBSCAN)
            summary["hough_final_hits_x"] += int(entry.nFinal)

    root_file.Close()
    return summary


def merge_summaries(summaries):
    merged = {}
    stage_counts = {}
    for summary in summaries:
        for key, value in summary.items():
            if key == "_stage_primary_counts":
                for stage, counts in value.items():
                    stage_counts.setdefault(stage, Counter()).update(counts)
                continue
            merged[key] = merged.get(key, 0) + value

    # This is the deliberately reconstruction-defined denominator: an ID is
    # eligible if it was the energy-leading primary of at least one object at
    # any recorded stage.  There is no Truth_Spill population or truth-match.
    denominator_ids = set()
    for counts in stage_counts.values():
        denominator_ids.update(counts)
    merged["reconstruction_primary_denominator"] = len(denominator_ids)
    for view in ("x", "y"):
        prefix = "candidate_" + view + "_"
        count = merged[prefix + "count"]
        merged[prefix + "completeness_mean"] = (
            merged[prefix + "completeness_sum"] / count if count else None
        )
        merged[prefix + "cleanliness_mean"] = (
            merged[prefix + "cleanliness_sum"] / count if count else None
        )
    for stage, views in (("seed", ("x", "y")), ("dbscan", ("x", "y")),
                         ("hough_final", ("x", "y")), ("final_2d", ("x", "y")),
                         ("final_3d", ("xy",))):
        for view in views:
            prefix = "stage_" + stage + "_" + view + "_"
            counts = stage_counts.get(stage + ":" + view.upper(), Counter())
            found = len(counts)
            quality_count = merged.get(prefix + "quality_count", 0)
            candidates = sum(counts.values())
            merged[prefix + "truth"] = len(denominator_ids)
            merged[prefix + "found"] = found
            merged[prefix + "candidates"] = candidates
            merged[prefix + "reco_efficiency"] = found / len(denominator_ids) if denominator_ids else None
            # Fraction of represented primary IDs that appear in more than one
            # object at this stage.  The raw object count remains printed too.
            merged[prefix + "multi_reco_rate"] = (
                sum(count > 1 for count in counts.values()) / found if found else None
            )
            merged[prefix + "completeness"] = merged.get(prefix + "completeness_sum", 0) / quality_count if quality_count else None
            merged[prefix + "cleanliness"] = merged.get(prefix + "cleanliness_sum", 0) / quality_count if quality_count else None
    paired_truth = merged.get("xy_source_pairs_with_truth", 0)
    merged["xy_pair_same_primary_rate"] = merged.get("xy_source_pairs_same_primary", 0) / paired_truth if paired_truth else None
    return merged


def run_reconstruction(executable, repo, config, input_file, output_file):
    environment = os.environ.copy()
    environment["TMS_DIR"] = str(repo)
    environment["TMS_TOML"] = str(config)
    environment["ND_PRODUCTION_TMSRECO_OUTFILE"] = str(output_file)
    command = [str(executable), str(input_file), str(output_file)]
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=environment)


def default_reco_executable(repo):
    candidates = [repo / "bin" / "ConvertToTMSTree.exe", repo / "bin" / "ConvertToTMSTree"]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def print_report(results):
    columns = [
        ("point", "point"),
        ("reco_tracks", "tracks"),
        ("lines_x", "X lines"),
        ("lines_y", "Y lines"),
        ("reconstruction_primary_denominator", "primary-ID union"),
    ]
    widths = [max(len(title), *(len(format_value(row.get(key))) for row in results))
              for key, title in columns]
    print("\n" + "  ".join(title.ljust(width) for (_, title), width in zip(columns, widths)))
    print("  ".join("-" * width for width in widths))
    for row in results:
        print("  ".join(format_value(row.get(key)).ljust(width)
                        for (key, _), width in zip(columns, widths)))
    print("\nAll stage rows use the same reconstruction-defined denominator: "
          "the union of global primary IDs leading any object in any recorded stage.")
    for stage, title in (("seed", "Hough seed"), ("dbscan", "Post-DBSCAN"),
                         ("hough_final", "Post-Hough/A*/extrapolation"),
                         ("final_2d", "Final 2-D candidates"),
                         ("final_3d", "Final XY 3-D tracks")):
        print_stage_table(results, stage, title)
    print_xy_pairing_table(results)


def print_stage_table(results, stage, title):
    print("\n" + title)
    columns = [("point", "point"), ("view", "view"), ("truth", "truth"),
               ("candidates", "candidates"), ("reco_efficiency", "reco eff"),
               ("multi_reco_rate", "multi-reco"), ("completeness", "completeness"),
               ("cleanliness", "cleanliness")]
    rows = []
    for result in results:
        for view in (("xy",) if stage == "final_3d" else ("x", "y")):
            prefix = "stage_" + stage + "_" + view + "_"
            if prefix + "truth" not in result:
                continue
            rows.append({"point": result["point"], "view": view.upper(),
                         **{name: result.get(prefix + name) for name in ("truth", "candidates", "reco_efficiency", "multi_reco_rate", "completeness", "cleanliness")}})
    if not rows:
        print("(requires fresh outputs with Hough_View_Truth)")
        return
    widths = [max(len(title), *(len(format_value(row.get(key))) for row in rows)) for key, title in columns]
    print("  ".join(title.ljust(width) for (_, title), width in zip(columns, widths)))
    print("  ".join("-" * width for width in widths))
    for row in rows:
        print("  ".join(format_value(row.get(key)).ljust(width) for (key, _), width in zip(columns, widths)))


def print_xy_pairing_table(results):
    print("\nXY pairing")
    columns = [("point", "point"), ("reco_tracks_from_xy_match", "pairs"),
               ("xy_source_pairs_same_primary", "same primary"),
               ("xy_source_pairs_different_primary", "wrong primary"),
               ("xy_source_pairs_missing_truth", "no truth"),
               ("xy_pair_same_primary_rate", "same-primary rate"),
               ("x_candidates_used_in_xy_match", "X paired"),
               ("y_candidates_used_in_xy_match", "Y paired")]
    widths = [max(len(title), *(len(format_value(row.get(key))) for row in results)) for key, title in columns]
    print("  ".join(title.ljust(width) for (_, title), width in zip(columns, widths)))
    print("  ".join("-" * width for width in widths))
    for row in results:
        print("  ".join(format_value(row.get(key)).ljust(width) for (key, _), width in zip(columns, widths)))


def format_value(value):
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, nargs="+", type=Path,
                        help="one or more identical edep-sim ROOT inputs for every scan point")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="directory for generated configs, ROOT outputs, and summary files")
    parser.add_argument("--reco-exe", type=Path,
                        help="ConvertToTMSTree executable (default: detect .exe or bare binary in <repo>/bin)")
    parser.add_argument("--config", type=Path,
                        help="base TOML (default: <repo>/config/TMS_Default_Config.toml)")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2],
                        help="dune-tms checkout containing the executable and config")
    parser.add_argument("--point", type=parse_point, action="append", default=[], metavar="NAME,H,E[,LEGACY]",
                        help="additional point: common Hough half-width, extrapolation multiplier, optional legacy-X true/false")
    parser.add_argument("--suite", choices=("pr289",),
                        help="run the diagnostic PR-289 A/B suite (legacy, nominal, Hough-only, extrapolation-only, retuned)")
    parser.add_argument("--no-diagnostics", action="store_true",
                        help="do not enable compact Hough_Diagnostics output in generated TOMLs")
    parser.add_argument("--skip-default", action="store_true",
                        help="do not add corrected_nominal,1.5,1 automatically")
    parser.add_argument("--max-events", type=int,
                        help="override Applications.MaximumNEvents in generated TOMLs")
    parser.add_argument("--force", action="store_true", help="replace existing output ROOT files")
    parser.add_argument("--analyze-only", action="store_true",
                        help="do not run reconstruction; summarise existing ROOT outputs")
    arguments = parser.parse_args()

    repo = arguments.repo.resolve()
    template = (arguments.config or repo / "config" / "TMS_Default_Config.toml").resolve()
    executable = (arguments.reco_exe or default_reco_executable(repo)).resolve()
    inputs = [path.resolve() for path in arguments.input]
    if not template.is_file():
        parser.error("base config does not exist: " + str(template))
    if not arguments.analyze_only and not executable.is_file():
        parser.error("reconstruction executable does not exist: " + str(executable))
    for input_file in inputs:
        if not input_file.is_file():
            parser.error("input does not exist: " + str(input_file))
    if arguments.max_events is not None and arguments.max_events < 1:
        parser.error("--max-events must be positive")
    points = [] if arguments.skip_default else [Point(*DEFAULT_POINT)]
    if arguments.suite == "pr289":
        points = [Point(*values) for values in PR289_SUITE]
    points.extend(arguments.point)
    if not points:
        parser.error("provide at least one --point or omit --skip-default")
    names = [point.name for point in points]
    if len(names) != len(set(names)):
        parser.error("every scan point needs a unique name")

    output_dir = arguments.output_dir.resolve()
    config_dir = output_dir / "configs"
    root_dir = output_dir / "root"
    config_dir.mkdir(parents=True, exist_ok=True)
    root_dir.mkdir(parents=True, exist_ok=True)

    report_rows = []
    for point in points:
        config_path = config_dir / (point.name + ".toml")
        write_config(template, config_path, point, arguments.max_events,
                     diagnostics=not arguments.no_diagnostics)
        summaries = []
        output_paths = []
        for input_index, input_file in enumerate(inputs):
            output_path = root_dir / point.name / f"{input_index:02d}_{input_file.stem}.root"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            if not arguments.analyze_only:
                if output_path.exists() and not arguments.force:
                    raise RuntimeError("output exists (use --force): " + str(output_path))
                if output_path.exists():
                    output_path.unlink()
                run_reconstruction(executable, repo, config_path, input_file, output_path)
            if not output_path.is_file():
                raise RuntimeError("expected reconstruction output was not created: " + str(output_path))
            summaries.append(summarise_output(output_path))
            output_paths.append(str(output_path))

        row = {
            "point": point.name,
            "hough_half_width": point.hough_half_width,
            "extrapolation_width_multiplier": point.extrapolation_multiplier,
            "legacy_x_bar_containment": point.legacy_x,
            "config": str(config_path),
            "outputs": output_paths,
        }
        row.update(merge_summaries(summaries))
        report_rows.append(row)

    (output_dir / "summary.json").write_text(json.dumps(report_rows, indent=2) + "\n")
    csv_fields = [key for key in report_rows[0] if key != "outputs"]
    with (output_dir / "summary.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=csv_fields)
        writer.writeheader()
        for row in report_rows:
            writer.writerow({key: row[key] for key in csv_fields})
    print_report(report_rows)
    print("\nWrote", output_dir / "summary.csv")
    print("Wrote", output_dir / "summary.json")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as error:
        print("error:", error, file=sys.stderr)
        sys.exit(1)
