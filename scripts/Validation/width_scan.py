#!/usr/bin/env python3
"""Run a reproducible transverse-width reconstruction scan.

Each configuration is run over the same input ROOT files.  The script keeps
the produced ROOT files, generated TOMLs, machine-readable summary, and a
compact CSV together under one output directory.

The reported track-per-truth-muon number is a scan diagnostic, not a matched
physics efficiency: one truth muon can produce zero, one, or several tracks.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


DEFAULT_POINT = ("default", 1.5, 1.5, 1.0, 1.0)


@dataclass
class Point:
    name: str
    hough_x: float
    hough_y: float
    extrapolation_x: float
    extrapolation_y: float


def parse_point(value):
    fields = value.split(",")
    if len(fields) != 5:
        raise argparse.ArgumentTypeError(
            "point must be NAME,HOUGH_X,HOUGH_Y,EXTRAPOLATION_X,EXTRAPOLATION_Y"
        )
    try:
        point = Point(fields[0], *[float(field) for field in fields[1:]])
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error))
    if not point.name or min(point.hough_x, point.hough_y,
                             point.extrapolation_x, point.extrapolation_y) <= 0:
        raise argparse.ArgumentTypeError("point values must all be positive")
    return point


def replace_toml_value(text, key, value):
    pattern = re.compile(r"^(\s*" + re.escape(key) + r"\s*=\s*)[^#\n]*(.*)$", re.MULTILINE)
    text, replacements = pattern.subn(r"\g<1>" + str(value) + r"\2", text)
    if replacements != 1:
        raise RuntimeError("expected exactly one TOML setting named " + key)
    return text


def write_config(template, destination, point, max_events):
    text = template.read_text()
    values = {
        "ContainmentHalfWidthX": point.hough_x,
        "ContainmentHalfWidthY": point.hough_y,
        "ContainmentWidthScaleX": point.extrapolation_x,
        "ContainmentWidthScaleY": point.extrapolation_y,
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
    truth_spill = root_file.Get("Truth_Spill")
    for branch in ("nLinesX", "nLinesY", "nLines3D", "nHitsInTrackX", "nHitsInTrackY"):
        get_branch(lines, branch)
    for branch in ("nTracks", "nHits", "TrackHitBarType"):
        get_branch(reco, branch)

    summary = {
        "line_entries": int(lines.GetEntries()),
        "lines_x": 0,
        "lines_y": 0,
        "lines_3d": 0,
        "line_hits_x": 0,
        "line_hits_y": 0,
        "slices_with_x_line": 0,
        "slices_with_y_line": 0,
        "slices_with_3d_line": 0,
        "reco_slices": int(reco.GetEntries()),
        "slices_with_tracks": 0,
        "reco_tracks": 0,
        "reco_track_hits": 0,
        "reco_track_hits_x": 0,
        "reco_track_hits_y": 0,
        "tracks_with_x_hits": 0,
        "tracks_with_y_hits": 0,
        "tracks_with_xy_hits": 0,
        "truth_spills": 0,
        "truth_primary_muons_touching_tms": 0,
        "truth_spills_with_primary_muon_touching_tms": 0,
    }

    for entry in lines:
        n_x = int(entry.nLinesX)
        n_y = int(entry.nLinesY)
        n_3d = int(entry.nLines3D)
        summary["lines_x"] += n_x
        summary["lines_y"] += n_y
        summary["lines_3d"] += n_3d
        summary["line_hits_x"] += sum(int(entry.nHitsInTrackX[index]) for index in range(n_x))
        summary["line_hits_y"] += sum(int(entry.nHitsInTrackY[index]) for index in range(n_y))
        summary["slices_with_x_line"] += n_x > 0
        summary["slices_with_y_line"] += n_y > 0
        summary["slices_with_3d_line"] += n_3d > 0

    for entry in reco:
        n_tracks = int(entry.nTracks)
        summary["reco_tracks"] += n_tracks
        summary["slices_with_tracks"] += n_tracks > 0
        for track in range(n_tracks):
            n_hits = int(entry.nHits[track])
            summary["reco_track_hits"] += n_hits
            has_x = False
            has_y = False
            for hit in range(n_hits):
                bar_type = int(entry.TrackHitBarType[track][hit])
                if bar_type == 0:  # TMS_Bar::kXBar
                    summary["reco_track_hits_x"] += 1
                    has_x = True
                elif bar_type == 1:  # TMS_Bar::kYBar
                    summary["reco_track_hits_y"] += 1
                    has_y = True
            summary["tracks_with_x_hits"] += has_x
            summary["tracks_with_y_hits"] += has_y
            summary["tracks_with_xy_hits"] += has_x and has_y

    if truth_spill:
        for branch in ("nTrueParticles", "PDG", "IsPrimary", "TMSFiducialTouch"):
            get_branch(truth_spill, branch)
        summary["truth_spills"] = int(truth_spill.GetEntries())
        for entry in truth_spill:
            n_target_muons = 0
            for particle in range(int(entry.nTrueParticles)):
                if (abs(int(entry.PDG[particle])) == 13 and bool(entry.IsPrimary[particle])
                        and bool(entry.TMSFiducialTouch[particle])):
                    n_target_muons += 1
            summary["truth_primary_muons_touching_tms"] += n_target_muons
            summary["truth_spills_with_primary_muon_touching_tms"] += n_target_muons > 0

    root_file.Close()
    return summary


def merge_summaries(summaries):
    merged = {}
    for summary in summaries:
        for key, value in summary.items():
            merged[key] = merged.get(key, 0) + value
    denominator = merged["truth_primary_muons_touching_tms"]
    merged["tracks_per_truth_primary_muon_touching_tms"] = (
        merged["reco_tracks"] / denominator if denominator else None
    )
    return merged


def run_reconstruction(executable, repo, config, input_file, output_file):
    environment = os.environ.copy()
    environment["TMS_DIR"] = str(repo)
    environment["TMS_TOML"] = str(config)
    environment["ND_PRODUCTION_TMSRECO_OUTFILE"] = str(output_file)
    command = [str(executable), str(input_file), str(output_file)]
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=environment)


def print_report(results):
    columns = [
        ("point", "point"),
        ("reco_tracks", "tracks"),
        ("slices_with_tracks", "track slices"),
        ("lines_x", "X lines"),
        ("lines_y", "Y lines"),
        ("lines_3d", "3D lines"),
        ("reco_track_hits_x", "X track hits"),
        ("reco_track_hits_y", "Y track hits"),
        ("truth_primary_muons_touching_tms", "truth #mu touch"),
        ("tracks_per_truth_primary_muon_touching_tms", "tracks / truth #mu"),
    ]
    widths = [max(len(title), *(len(format_value(row.get(key))) for row in results))
              for key, title in columns]
    print("\n" + "  ".join(title.ljust(width) for (_, title), width in zip(columns, widths)))
    print("  ".join("-" * width for width in widths))
    for row in results:
        print("  ".join(format_value(row.get(key)).ljust(width)
                        for (key, _), width in zip(columns, widths)))
    print("\n'tracks / truth #mu' is a relative scan diagnostic, not a matched efficiency.")


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
                        help="ConvertToTMSTree executable (default: <repo>/bin/ConvertToTMSTree)")
    parser.add_argument("--config", type=Path,
                        help="base TOML (default: <repo>/config/TMS_Default_Config.toml)")
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2],
                        help="dune-tms checkout containing the executable and config")
    parser.add_argument("--point", type=parse_point, action="append", default=[], metavar="NAME,HX,HY,EX,EY",
                        help="additional point: Hough half-widths and extrapolation width scales")
    parser.add_argument("--skip-default", action="store_true",
                        help="do not add default,1.5,1.5,1,1 automatically")
    parser.add_argument("--max-events", type=int,
                        help="override Applications.MaximumNEvents in generated TOMLs")
    parser.add_argument("--force", action="store_true", help="replace existing output ROOT files")
    parser.add_argument("--analyze-only", action="store_true",
                        help="do not run reconstruction; summarise existing ROOT outputs")
    arguments = parser.parse_args()

    repo = arguments.repo.resolve()
    template = (arguments.config or repo / "config" / "TMS_Default_Config.toml").resolve()
    executable = (arguments.reco_exe or repo / "bin" / "ConvertToTMSTree").resolve()
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
        write_config(template, config_path, point, arguments.max_events)
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
            "hough_half_width_x": point.hough_x,
            "hough_half_width_y": point.hough_y,
            "extrapolation_width_scale_x": point.extrapolation_x,
            "extrapolation_width_scale_y": point.extrapolation_y,
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
