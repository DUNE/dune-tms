import os
import collections
import logging
import argparse
from dataclasses import dataclass
from typing import Optional, Dict

import ROOT

ROOT.gROOT.SetBatch(True)

import dunestyle.root as dunestyle

ROOT.gStyle.SetPaintTextFormat(".2f")

import sys


# Basic logging setup (keeps prints for key user messages)
logging.basicConfig(level=logging.INFO)

# Visual constants
LEGEND_TEXT_SIZE = 0.04
# Progress logging cadence (base plot saves)
PROGRESS_EVERY = 50
# Stack drawing constants
STACK_COLORS = [
    ROOT.kBlue,
    ROOT.kRed,
    ROOT.kGreen,
    ROOT.kMagenta,
    ROOT.kBlack,
    ROOT.kOrange,
]
STACK_LINE_STYLES = [1, 2, 7]


@dataclass
class HistInfo:
    """Lightweight description of a histogram inferred from its name.

    Captures subdirectory, image filename, and flags used to control
    drawing (stacking, log-y, normalization, etc.).
    """
    raw_name: str
    subdir: str
    image_name: str
    # Stack-related
    is_stack: bool = False
    stack_key: Optional[str] = None
    stack_item_name: Optional[str] = None
    stack_mode: Optional[str] = None  # "stack" or "nostack"
    # Rendering flags inferred from name
    log_y: bool = False
    area_norm: bool = False
    optsmooth: bool = False
    normalized: bool = False
    scale_width: bool = False
    # Reconstruction efficiency role
    recoeff_role: Optional[str] = None  # "numerator" | "denominator" | None
    recoeff_key: Optional[str] = None


def parse_hist_name(hist_name: str) -> HistInfo:
    """Parse a ROOT object name into a HistInfo structure.

    Name conventions supported:
    - Subdirectories separated by "__" (e.g., a__b__hist -> a/b with image "hist").
    - Stack grouping via "_stack_" or "_nostack_".
    - Flags via substrings: width, normalized, log, area_norm, optsmooth,
      numerator, denominator.
    """
    name = hist_name.strip()
    split = name.split("__")
    subdir = os.path.join(*split[:-1]) if len(split) > 1 else ""
    image_name = split[-1]

    # Flags derived from the full name
    scale_width = ("width" in name)
    normalized = ("normalized" in image_name)
    optsmooth = ("optsmooth" in name)
    area_norm = ("area_norm" in name)
    log_y = ("log" in name)

    recoeff_role = None
    recoeff_key = None
    if "numerator" in name:
        recoeff_role = "numerator"
        recoeff_key = name.replace("_numerator", "")
    elif "denominator" in name:
        recoeff_role = "denominator"
        recoeff_key = name.replace("_denominator", "")

    # Stack parsing
    is_stack = False
    stack_key = None
    stack_item_name = None
    stack_mode = None
    if "_stack_" in name:
        parts = name.split("_stack_")
        stack_key, stack_item_name = parts[0], parts[1]
        is_stack = True
        stack_mode = "stack"
    elif "_nostack_" in name:
        parts = name.split("_nostack_")
        stack_key, stack_item_name = parts[0], parts[1]
        is_stack = True
        stack_mode = "nostack"

    info = HistInfo(
        raw_name=name,
        subdir=subdir,
        image_name=image_name,
        is_stack=is_stack,
        stack_key=stack_key,
        stack_item_name=stack_item_name,
        stack_mode=stack_mode,
        log_y=log_y,
        area_norm=area_norm,
        optsmooth=optsmooth,
        normalized=normalized,
        scale_width=scale_width,
        recoeff_role=recoeff_role,
        recoeff_key=recoeff_key,
    )
    return info


def get_subdir_and_name(hist_name):
    """Backward-compatible helper returning (subdir, image_name).

    Prefer using parse_hist_name directly; kept to minimize diff.
    """
    # Backward compatible helper; prefer parse_hist_name
    info = parse_hist_name(hist_name)
    logging.debug("Name parse: %s -> (%s, %s)", hist_name, info.subdir, info.image_name)
    return info.subdir, info.image_name


def ensure_subdir(base_dir: str, subdir: str, extra: Optional[str] = None) -> str:
    """Ensure output subdirectory exists and return its path.

    If extra is provided, it is appended (e.g., "additional_plots").
    """
    path = os.path.join(base_dir, subdir) if subdir else base_dir
    if extra:
        path = os.path.join(path, extra)
    os.makedirs(path, exist_ok=True)
    return path


def setup_th2_style(h: ROOT.TH2, use_text_overlay: bool) -> None:
    """Apply consistent 2D histogram styling and optional text overlay tweaks."""
    h.GetYaxis().SetTitleOffset(1.4)
    h.GetZaxis().SetTitleOffset(0.5)
    if use_text_overlay:
        h.SetMarkerColor(ROOT.kRed + 1)
        h.GetZaxis().SetRangeUser(0.001, h.GetMaximum())


def setup_th1_style(h: ROOT.TH1) -> None:
    """Apply consistent 1D histogram styling and headroom."""
    top = h.GetMaximum() * 1.2
    h.GetYaxis().SetRangeUser(0, top)
    h.SetLineColor(ROOT.kBlack)


def save_canvas(canvas: ROOT.TCanvas, outdir: str, image_name: str) -> None:
    """Save the current canvas to outdir/image_name.png.

    Temporarily raises ROOT's gErrorIgnoreLevel to suppress noisy
    "Info in <TCanvas::Print>" messages emitted on save.
    """
    path = os.path.join(outdir, image_name + ".png")
    old_level = int(ROOT.gErrorIgnoreLevel)
    try:
        ROOT.gErrorIgnoreLevel = ROOT.kWarning  # suppress Info, keep warnings/errors
        canvas.Print(path)
    finally:
        ROOT.gErrorIgnoreLevel = old_level


# TODO erase the output dir subdirs (except event displays) because sometimes the hists don't exist if they weren't filled at least once
def draw_histograms(input_file):
    # Open the input ROOT file
    root_file = ROOT.TFile.Open(input_file)

    # Validate file open
    if not root_file or root_file.IsZombie():
        logging.error("Failed to open ROOT file: %s", input_file)
        return

    # Create a directory to save images
    output_dir = os.path.splitext(input_file)[0] + "_images"
    os.makedirs(output_dir, exist_ok=True)
    logging.info("Output directory: %s", output_dir)

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    # Save special plots here
    recoeff_plots_numerators = dict()
    recoeff_plots_denominators = dict()
    stack_plots = collections.defaultdict(dict)

    # Loop over all keys in the ROOT file
    base_plots_rendered = 0
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        name = obj.GetName()
        info = parse_hist_name(name)

        output_subdir = ensure_subdir(
            output_dir,
            info.subdir,
            "additional_plots" if (info.recoeff_role is not None or info.is_stack) else None,
        )

        if info.scale_width:
            obj.Scale(1, "width")

        if isinstance(obj, ROOT.TH2):
            # For 2D histograms, draw with "colz" option and save as png
            use_text_overlay = info.normalized and obj.GetXaxis().GetNbins() < 20
            setup_th2_style(obj, use_text_overlay)
            obj.Draw("colz text" if use_text_overlay else "colz")
            dunestyle.Simulation()
            logging.debug("%s integral: %s", obj.GetName(), obj.Integral())
            save_canvas(canvas, output_subdir, info.image_name)
            base_plots_rendered += 1
            if base_plots_rendered % PROGRESS_EVERY == 0:
                logging.info("Rendered %d base plots...", base_plots_rendered)
        elif isinstance(obj, ROOT.TH1):
            # For 1D histograms, draw and save as png
            setup_th1_style(obj)
            obj.Draw()
            dunestyle.Simulation()
            # print(f"{obj.GetName()} integral: {obj.Integral()}")
            save_canvas(canvas, output_subdir, info.image_name)
            base_plots_rendered += 1
            if base_plots_rendered % PROGRESS_EVERY == 0:
                logging.info("Rendered %d base plots...", base_plots_rendered)
        # Collect for reco efficiency if relevant
        if info.recoeff_role is not None and info.recoeff_key is not None:
            key_base = info.recoeff_key
            if info.recoeff_role == "numerator":
                recoeff_plots_numerators[key_base] = obj
            elif info.recoeff_role == "denominator":
                recoeff_plots_denominators[key_base] = obj
        # Collect for stacks if relevant
        if info.is_stack and info.stack_key is not None and info.stack_item_name is not None:
            stack_plots[info.stack_key][info.stack_item_name] = obj
            # stack_plots[stack_key + "_log"][split_stack[1]] = obj

    logging.info("Base plots rendered: %d", base_plots_rendered)

    # Draw the stacks
    def draw_stack_group(canvas: ROOT.TCanvas, name: str, hist_and_name: Dict[str, ROOT.TH1], output_dir: str) -> None:
        """Draw a group of stacked or overlaid histograms and save the image.

        Behavior is inferred from the group name: supports log y, area normalization,
        and smoothing options consistent with the original script.
        """
        if not hist_and_name:
            logging.debug("Empty stack group: %s (skipping)", name)
            return

        index = 0
        first_hist = None
        ymax = 0
        hist_stack = ROOT.THStack()
        leg = ROOT.TLegend(0.2, 0.67, 0.8, 0.82)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetNColumns(2)
        leg.SetTextSize(LEGEND_TEXT_SIZE)

        log = ("log" in name)
        canvas.SetLogy(True if log else False)

        area_norm = ("area_norm" in name)

        headroom = 1.3
        if log:
            headroom = 5

        items = list(hist_and_name.items())
        items.sort()
        for item_name, hist in items:
            if area_norm:
                integral = hist.Integral()
                if integral != 0:
                    hist.Scale(1 / integral)
            hist.SetLineColor(STACK_COLORS[index % len(STACK_COLORS)])
            hist.SetLineStyle(STACK_LINE_STYLES[index % len(STACK_LINE_STYLES)])
            if index == 0:
                first_hist = hist
            ymax = max(ymax, hist.GetMaximum())
            hist_stack.Add(hist)
            name_for_legend = item_name
            if ":" in hist.GetTitle():
                name_for_legend = hist.GetTitle().split(":")[1].strip()
            leg.AddEntry(hist, name_for_legend, "lep")
            index += 1

        opts = "hist"
        if first_hist is not None:
            if "nostack" in first_hist.GetName():
                opts += " nostack"
            else:
                opts += " stack"
            if "optsmooth" in first_hist.GetName():
                opts += " C"
        hist_stack.Draw(opts)

        if first_hist is not None:
            ytitle = first_hist.GetYaxis().GetTitle()
            hist_stack.SetTitle(first_hist.GetTitle().split(":")[0])
            hist_stack.GetXaxis().SetTitle(first_hist.GetXaxis().GetTitle())
            hist_stack.GetYaxis().SetTitle(ytitle)

        ymin = 0
        if log:
            ymin = 10
        for hist in hist_and_name.values():
            hist.GetYaxis().SetRangeUser(ymin, ymax * headroom)

        hist_stack.Draw(opts)
        leg.Draw()
        dunestyle.Simulation()
        subdir, image_name = get_subdir_and_name(name)
        output_subdir = ensure_subdir(output_dir, subdir)
        logging.debug("Saving in %s", os.path.join(output_subdir, image_name + ".png"))
        save_canvas(canvas, output_subdir, image_name)

    # Loop over each stack
    logging.info("Drawing stacks for %d groups...", len(stack_plots))
    stack_groups_rendered = 0
    for name, hist_and_name in stack_plots.items():
        draw_stack_group(canvas, name, hist_and_name, output_dir)
        stack_groups_rendered += 1
    logging.info("Stacks rendered: %d", stack_groups_rendered)

    # Turn off log y if it's already on
    canvas.SetLogy(False)

    # Draw reco eff
    def draw_recoefficiencies(canvas: ROOT.TCanvas, num_map: Dict[str, ROOT.TH1], den_map: Dict[str, ROOT.TH1], output_dir: str) -> None:
        """Compute and draw reconstruction efficiencies from numerator/denominator pairs."""
        denominator_only = set(den_map.keys()) - set(num_map.keys())
        numerator_only = set(num_map.keys()) - set(den_map.keys())

        for name in sorted(denominator_only):
            logging.warning(
                "Missing reco-eff numerator for %s; drawing zero-efficiency plot from denominator",
                name,
            )

        for name in sorted(numerator_only):
            logging.error(
                "!!! RECO-EFF NUMERATOR WITHOUT DENOMINATOR: %s. "
                "This should not happen; skipping this efficiency plot.",
                name,
            )

        all_names = set(num_map.keys()) & set(den_map.keys())
        all_names |= denominator_only
        for name in sorted(all_names):
            if name not in den_map:
                continue

            denominator = den_map[name]
            if name in num_map:
                numerator = num_map[name]
            else:
                numerator = denominator.Clone(name + "_numerator")
                numerator.Reset()
                numerator.SetTitle(
                    denominator.GetTitle().replace(": Denominator", "").strip()
                )

            # Clone before divide to avoid altering originals
            eff = numerator.Clone()
            newtitle = numerator.GetTitle()
            newtitle = newtitle.replace(": Numerator", "").strip()
            eff.SetTitle(newtitle)
            eff.Divide(denominator)

            eff.GetYaxis().SetRangeUser(0, 1.2)
            eff.GetYaxis().SetTitle("Reconstruction Efficiency")
            xmin = eff.GetXaxis().GetXmin()
            xmax = eff.GetXaxis().GetXmax()
            eff.Draw()
            dunestyle.Simulation()
            line = ROOT.TLine(xmin, 1, xmax, 1)
            line.SetLineStyle(2)
            line.Draw()

            out_name = numerator.GetName().replace("_numerator", "")
            subdir, image_name = get_subdir_and_name(out_name)
            output_subdir = ensure_subdir(output_dir, subdir)
            save_canvas(canvas, output_subdir, image_name)

    # Reco-eff: compute intersection once and log counts
    recoeff_pairs = set(recoeff_plots_numerators.keys()) & set(recoeff_plots_denominators.keys())
    recoeff_denominator_only = set(recoeff_plots_denominators.keys()) - set(recoeff_plots_numerators.keys())
    logging.info(
        "Rendering reconstruction efficiencies for %d pairs and %d denominator-only plots...",
        len(recoeff_pairs),
        len(recoeff_denominator_only),
    )
    draw_recoefficiencies(canvas, recoeff_plots_numerators, recoeff_plots_denominators, output_dir)
    logging.info(
        "Reconstruction efficiencies rendered: %d",
        len(recoeff_pairs) + len(recoeff_denominator_only),
    )

    # Close the input ROOT file
    root_file.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw all histograms from a ROOT file and save images.")
    parser.add_argument("input_root_file", help="Path to the input ROOT file")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose debug logging")
    args = parser.parse_args()

    # Configure logging
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    input_file = args.input_root_file
    if not os.path.isfile(input_file):
        print("Error: Input file does not exist!")
        sys.exit(1)

    # Initialize ROOT
    ROOT.gROOT.SetBatch(True)  # Prevent ROOT from trying to open X11 windows
    ROOT.gStyle.SetOptStat(0)  # Hide statistics box in histograms

    # Call function to draw histograms
    draw_histograms(input_file)
