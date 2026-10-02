#!/usr/bin/env python3
"""Detector-simulation validation: hit-level metrics and plots comparing ConvertToTMSTree outputs.

Compares two or more runs over the same input files (the first is the reference) and writes a Markdown report, a
JSON file with every number in it, and a few plots. Only simulation-level quantities: no tracking or slicing.

Sections (muon-dominated hits = TrueLeptonicEnergy / TrueHitE > 0.95, as in the 2026-09 segmentation benchmark):
  1. Totals: readout hits and PE above threshold per spill; light yield (PE after fibers per MeV) for muon- and
     hadron-dominated hits.
  2. Hit survival (not pedestal-suppressed) vs true path length in the bar (TrueHitDx, if the run has it) and vs true
     deposited energy.
  3. Photostatistics: reco PE mean and relative spread per path-length (or energy) band; PE distribution at partial fill.
  4. Timing: mean (reco - true) hit time vs the number of Geant4 contributions (TrueNTrueParticles) at full crossings.
  5. Coincident-hit merging: groups of more than one final true hit for one particle in one bar and slice, by time gap.
     Within the readout window these are merge failures; beyond it, separate readouts by design.
  6. Plane coverage: per (particle, plane) with >= 0.5 MeV deposited by a muon-dominated particle, the fraction with a
     surviving hit, by deposited energy.
  7. Light provenance (runs with Sim.DetSim.UseResponseElements = true): photon count vs detected PE, most-light vs
     most-energy particle, first photon vs most-light particle.

Inputs per run: a directory with <NNNNNNN>_RecoCandidates.root (Truth_Info, Reco_Tree) and <NNNNNNN>_Readout.root
(TMS tree) for each file number.

Usage:
  detsim_validation.py --out <dir> [--files 1-15] <label>=<run dir> [<label>=<run dir> ...]
Needs uproot, awkward, numpy, matplotlib.
"""
import argparse
import json
import os
import sys

import awkward as ak
import numpy as np
import uproot

LEPTONIC_CUT = 0.95
THRESHOLD_PE = 3.0  # default pedestal threshold, drawn on the PE plot only
DX_EDGES = [0, 2, 4, 6, 8, 12, 20, 50]
E_EDGES = [0, 0.25, 0.5, 1, 2, 4, 1e9]
PLANE_E_EDGES = [0.5, 1, 2, 4, 1e9]
NSEG_BINS = [1, 2, 3, 4, 5]  # last bin is 5+
DT_EDGES = [0, 10, 120, 1000, 1e5, np.inf]  # ns; 120 ns = default readout window
FULL_DX = (8, 12)  # mm, full crossing of a 1 cm bar
FULL_E = (1.5, 2.5)  # MeV, used instead when the run has no TrueHitDx
PARTIAL_DX = (4, 6)
PARTIAL_E = (0.5, 1.0)

TRUTH = ["NTrueHits", "TrueHitE", "TrueLeptonicEnergy", "TrueHitT", "TrueRecoHitT", "TrueRecoHitPE",
         "TrueRecoHitIsPedSupped", "TrueNTrueParticles", "TrueHitPlane", "TrueHitView", "TrueHitBar",
         "TrueHitPrimaryId", "TrueHitVertexId", "TrueHitPEAfterFibersShortPath", "TrueHitPEAfterFibersLongPath"]
OPTIONAL = ["TrueHitDx", "TrueHitNPhotons", "TrueHitPrimaryIdByLight", "TrueHitVertexIdByLight", "TrueHitLightShare",
            "TrueHitFirstPhotonPrimaryId", "TrueHitFirstPhotonVertexId"]

# Plot style: validated reference categorical slots in fixed order; marker and dash give a second encoding
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
DASHES = ["-", "--", "-.", ":"]
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"


def parse_files(spec):
    files = []
    for part in spec.split(","):
        lo, _, hi = part.partition("-")
        files += list(range(int(lo), int(hi or lo) + 1))
    return files


def load(run_dir, files):
    """Per-hit columns (flattened over files and slices) plus per-file totals for one run."""
    cols, slice_entry, file_idx = {}, [], []
    tot = dict(spills=0, readout_hits=0, hits_above_thr=0, pe_above_thr=0.0)
    present = None
    for fi, n in enumerate(files):
        with uproot.open(f"{run_dir}/{n:07d}_RecoCandidates.root") as f:
            tree = f["Truth_Info"]
            if present is None:
                present = [c for c in OPTIONAL if c in tree.keys()]
            a = tree.arrays(TRUTH + present, library="ak")
        nh = a["NTrueHits"]
        mask = ak.local_index(a["TrueHitE"]) < nh
        for c in TRUTH[1:] + present:
            cols.setdefault(c, []).append(ak.to_numpy(ak.flatten(a[c][mask])))
        counts = ak.to_numpy(nh)
        slice_entry.append(np.repeat(np.arange(len(counts)), counts))
        file_idx.append(np.full(int(counts.sum()), fi))
        with uproot.open(f"{run_dir}/{n:07d}_Readout.root") as f:
            r = f["TMS"].arrays(["NRecoHits", "RecoHitPE", "RecoHitIsPedSupped"], library="ak")
        sup = ak.values_astype(r["RecoHitIsPedSupped"], bool)
        tot["spills"] += len(r["NRecoHits"])
        tot["readout_hits"] += int(ak.sum(r["NRecoHits"]))
        tot["hits_above_thr"] += int(ak.sum(~sup))
        tot["pe_above_thr"] += float(ak.sum(r["RecoHitPE"][~sup]))
    d = {c: np.concatenate(v) for c, v in cols.items()}
    d["entry"] = np.concatenate(slice_entry)
    d["file"] = np.concatenate(file_idx)
    d["survived"] = ~d["TrueRecoHitIsPedSupped"].astype(bool)
    e = d["TrueHitE"]
    d["lep_frac"] = np.divide(d["TrueLeptonicEnergy"], e, out=np.zeros_like(e), where=e > 0)
    d["pe_fibers"] = d["TrueHitPEAfterFibersShortPath"] + d["TrueHitPEAfterFibersLongPath"]
    return d, tot


def band_mask(x, lo, hi):
    return (x >= lo) & (x < hi)


def groups(keys):
    """Sort order, group start indices and group sizes for rows with identical keys (list of int arrays)."""
    order = np.lexsort(keys[::-1])
    k = np.stack([a[order] for a in keys])
    new = np.ones(len(order), bool)
    new[1:] = np.any(k[:, 1:] != k[:, :-1], axis=0)
    starts = np.nonzero(new)[0]
    return order, starts, np.diff(np.append(starts, len(order)))


def analyze(d, tot, has_dx):
    """All metrics for one run; has_dx = use TrueHitDx bands (only if every compared run has it)."""
    m = {}
    mu = d["lep_frac"] > LEPTONIC_CUT
    had = (d["TrueHitE"] > 0) & (d["lep_frac"] < 1 - LEPTONIC_CUT)
    # 1. totals
    m["totals"] = dict(
        spills=tot["spills"],
        readout_hits_per_spill=tot["readout_hits"] / tot["spills"],
        hits_above_thr_per_spill=tot["hits_above_thr"] / tot["spills"],
        pe_above_thr_per_spill=tot["pe_above_thr"] / tot["spills"],
        pe_per_mev_muon=float(d["pe_fibers"][mu].sum() / d["TrueHitE"][mu].sum()),
        pe_per_mev_hadronic=float(d["pe_fibers"][had].sum() / d["TrueHitE"][had].sum()),
        true_hits=int(len(d["TrueHitE"])),
    )
    # 2. survival
    def survival(x, edges):
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = mu & band_mask(x, lo, hi)
            n = int(s.sum())
            p = float(d["survived"][s].mean()) if n else float("nan")
            rows.append(dict(lo=lo, hi=hi, n=n, eff=p, err=float(np.sqrt(p * (1 - p) / n)) if n else float("nan")))
        return rows
    if has_dx:
        m["survival_vs_dx"] = survival(d["TrueHitDx"], DX_EDGES)
    m["survival_vs_e"] = survival(d["TrueHitE"], E_EDGES)
    # 3. photostatistics per band
    var, edges = ("TrueHitDx", DX_EDGES) if has_dx else ("TrueHitE", E_EDGES)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = mu & band_mask(d[var], lo, hi)
        pe = d["TrueRecoHitPE"][s]
        rows.append(dict(lo=lo, hi=hi, n=int(s.sum()), mean_pe=float(pe.mean()) if len(pe) else float("nan"),
                         rel_rms=float(pe.std() / pe.mean()) if len(pe) and pe.mean() > 0 else float("nan"),
                         below_thr=float(np.mean(pe < THRESHOLD_PE)) if len(pe) else float("nan")))
    m["pe_by_band"] = dict(variable=var, rows=rows)
    m["integer_pe_after_fibers"] = float(np.mean(np.isclose(d["pe_fibers"], np.round(d["pe_fibers"]))))
    # 4. timing at full crossings
    full = mu & (band_mask(d["TrueHitDx"], *FULL_DX) if has_dx else band_mask(d["TrueHitE"], *FULL_E))
    dt = d["TrueRecoHitT"] - d["TrueHitT"]
    trows = []
    for k in NSEG_BINS:
        s = full & ((d["TrueNTrueParticles"] >= k) if k == NSEG_BINS[-1] else (d["TrueNTrueParticles"] == k))
        x = dt[s]
        trows.append(dict(nseg=f"{k}+" if k == NSEG_BINS[-1] else str(k), n=int(s.sum()),
                          mean=float(x.mean()) if len(x) else float("nan"),
                          sem=float(x.std() / np.sqrt(len(x))) if len(x) > 1 else float("nan"),
                          rms=float(x.std()) if len(x) else float("nan")))
    m["time_vs_nseg"] = dict(selection=f"TrueHitDx {FULL_DX} mm" if has_dx else f"TrueHitE {FULL_E} MeV", rows=trows)
    # 5. coincident-hit merging
    idx = np.nonzero(mu)[0]
    keys = [d[c][idx].astype(np.int64) for c in ("file", "entry", "TrueHitPlane", "TrueHitView", "TrueHitBar",
                                                  "TrueHitPrimaryId", "TrueHitVertexId")]
    order, starts, sizes = groups(keys)
    t = d["TrueHitT"][idx][order]
    span = np.maximum.reduceat(t, starts) - np.minimum.reduceat(t, starts)
    multi = sizes > 1
    gaps = span[multi]
    m["merging"] = dict(groups=int(len(sizes)), multi_groups=int(multi.sum()),
                        rate=float(multi.mean()),
                        by_gap=[dict(lo=lo, hi=hi, n=int(band_mask(gaps, lo, hi).sum()))
                                for lo, hi in zip(DT_EDGES[:-1], DT_EDGES[1:])],
                        gaps=gaps.tolist())
    # 6. plane coverage
    keys = [d[c][idx].astype(np.int64) for c in ("file", "entry", "TrueHitPrimaryId", "TrueHitVertexId",
                                                  "TrueHitView", "TrueHitPlane")]
    order, starts, sizes = groups(keys)
    esum = np.add.reduceat(d["TrueHitE"][idx][order], starts)
    ok = np.maximum.reduceat(d["survived"][idx][order].astype(np.int8), starts).astype(bool)
    prow = []
    for lo, hi in zip(PLANE_E_EDGES[:-1], PLANE_E_EDGES[1:]):
        s = band_mask(esum, lo, hi)
        prow.append(dict(lo=lo, hi=hi, n=int(s.sum()), eff=float(ok[s].mean()) if s.any() else float("nan")))
    s = esum >= PLANE_E_EDGES[0]
    m["plane_coverage"] = dict(rows=prow, all=dict(n=int(s.sum()), eff=float(ok[s].mean())))
    # 7. light provenance
    if "TrueHitNPhotons" in d and np.any(d["TrueHitNPhotons"] >= 0):
        f = d["TrueHitNPhotons"] >= 0
        same = (d["TrueHitPrimaryIdByLight"] == d["TrueHitPrimaryId"]) & (d["TrueHitVertexIdByLight"] == d["TrueHitVertexId"])
        first = ((d["TrueHitFirstPhotonPrimaryId"] == d["TrueHitPrimaryIdByLight"])
                 & (d["TrueHitFirstPhotonVertexId"] == d["TrueHitVertexIdByLight"]))
        multi_c = d["TrueNTrueParticles"] > 1
        unfilled = ~f
        m["provenance"] = dict(
            filled=float(f.mean()),
            unfilled_with_zero_pe=float(np.mean(d["pe_fibers"][unfilled] == 0)) if unfilled.any() else 1.0,
            nphotons_equals_pe=float(np.mean(np.isclose(d["TrueHitNPhotons"][f], d["pe_fibers"][f]))),
            light_eq_energy_multi=float(same[f & multi_c].mean()),
            first_eq_light_multi=float(first[f & multi_c].mean()),
            light_share_below_half_multi=float(np.mean(d["TrueHitLightShare"][f & multi_c] < 0.5)),
            light_share_multi=d["TrueHitLightShare"][f & multi_c].tolist(),
        )
    # kept for plots only
    pb = mu & (band_mask(d["TrueHitDx"], *PARTIAL_DX) if has_dx else band_mask(d["TrueHitE"], *PARTIAL_E))
    m["_partial_pe"] = d["TrueRecoHitPE"][pb]
    m["_partial_label"] = f"TrueHitDx {PARTIAL_DX[0]}-{PARTIAL_DX[1]} mm" if has_dx else f"TrueHitE {PARTIAL_E[0]}-{PARTIAL_E[1]} MeV"
    return m


# ---------------------------------------------------------------- report
def fmt(x, nd=4):
    return "nan" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def rel(a, b):
    return f"{100 * (b - a) / a:+.2f}%" if a else "n/a"


def band_label(r):
    return f">= {r['lo']:g}" if r["hi"] >= 1e8 else f"{r['lo']:g}-{r['hi']:g}"


def report(res, labels, files):
    ref = labels[0]
    L = [f"# Detector-simulation validation\n", f"Runs: " + ", ".join(f"`{l}`" for l in labels) +
         f" (reference `{ref}`); files {files[0]}-{files[-1]} ({len(files)}).\n",
         "Muon-dominated hits: TrueLeptonicEnergy / TrueHitE > 0.95.\n"]

    def table(title, header, rows):
        L.append(f"\n## {title}\n")
        L.append("| " + " | ".join(header) + " |")
        L.append("|" + "---|" * len(header))
        L.extend("| " + " | ".join(r) + " |" for r in rows)

    tk = [("spills", "spills", 0), ("readout_hits_per_spill", "readout hits / spill", 1),
          ("hits_above_thr_per_spill", "hits above threshold / spill", 1), ("pe_above_thr_per_spill", "PE above threshold / spill", 0),
          ("pe_per_mev_muon", "PE after fibers / MeV, muon hits", 4), ("pe_per_mev_hadronic", "PE after fibers / MeV, hadronic hits", 4)]
    rows = []
    for k, name, nd in tk:
        vals = [res[l]["totals"][k] for l in labels]
        rows.append([name] + [f"{v:.{nd}f}" if isinstance(v, float) else str(v) for v in vals] +
                    [rel(vals[0], v) for v in vals[1:]])
    table("1. Totals", ["quantity"] + labels + [f"{l} vs {ref}" for l in labels[1:]], rows)

    for key, var in (("survival_vs_dx", "TrueHitDx (mm)"), ("survival_vs_e", "TrueHitE (MeV)")):
        if all(key in res[l] for l in labels):
            rows = []
            for i, r0 in enumerate(res[ref][key]):
                rows.append([band_label(r0)] + [f"{res[l][key][i]['eff']:.4f} ({res[l][key][i]['n']})" for l in labels] +
                            [f"{100 * (res[l][key][i]['eff'] - r0['eff']):+.2f} pp" for l in labels[1:]])
            table(f"2. Hit survival vs {var}", [var] + [f"{l} eff (N)" for l in labels] + [f"{l} - {ref}" for l in labels[1:]], rows)

    var = res[ref]["pe_by_band"]["variable"]
    rows = []
    for i, r0 in enumerate(res[ref]["pe_by_band"]["rows"]):
        cells = [band_label(r0)]
        for l in labels:
            r = res[l]["pe_by_band"]["rows"][i]
            cells.append(f"{fmt(r['mean_pe'], 2)} / {fmt(r['rel_rms'], 3)} / {fmt(r['below_thr'], 4)}")
        rows.append(cells)
    table(f"3. Reco PE by {var} band: mean / relative rms / fraction below {THRESHOLD_PE:g} PE",
          [var] + labels, rows)
    L.append("\nFraction of true hits with integer PE after fibers (Poisson after attenuation gives 1): " +
             ", ".join(f"`{l}` {res[l]['integer_pe_after_fibers']:.4f}" for l in labels) + "\n")

    sel = res[ref]["time_vs_nseg"]["selection"]
    rows = []
    for i, r0 in enumerate(res[ref]["time_vs_nseg"]["rows"]):
        rows.append([r0["nseg"]] + [f"{fmt(res[l]['time_vs_nseg']['rows'][i]['mean'], 2)} +- "
                                    f"{fmt(res[l]['time_vs_nseg']['rows'][i]['sem'], 2)} ({res[l]['time_vs_nseg']['rows'][i]['n']})"
                                    for l in labels])
    table(f"4. Mean reco - true hit time (ns) vs TrueNTrueParticles, {sel}", ["contributions"] + labels, rows)
    L.append("\nIn real spills the contribution count is not the only thing that varies: hits with more contributions "
             "carry more energy (more photons, so an earlier first photon) and sit at different distances from the "
             "readout. A segmentation artifact shows as a large trend (several ns); for a clean test at fixed physics "
             "use app/ArtificialResegmentationTest.\n")

    rows = [["groups (particle x bar x slice)"] + [str(res[l]["merging"]["groups"]) for l in labels],
            ["with > 1 final hit"] + [f"{res[l]['merging']['multi_groups']} ({100 * res[l]['merging']['rate']:.4f}%)" for l in labels]]
    for i, g in enumerate(res[ref]["merging"]["by_gap"]):
        hi = "inf" if np.isinf(g["hi"]) else f"{g['hi']:g}"
        rows.append([f"  time gap {g['lo']:g}-{hi} ns"] + [str(res[l]["merging"]["by_gap"][i]["n"]) for l in labels])
    table("5. Coincident-hit merging (muon-dominated)", ["quantity"] + labels, rows)
    L.append("\nGaps under the readout window (120 ns by default) are merge failures; longer gaps are separate readouts.\n")

    rows = []
    for i, r0 in enumerate(res[ref]["plane_coverage"]["rows"]):
        rows.append([band_label(r0)] + [f"{res[l]['plane_coverage']['rows'][i]['eff']:.4f} ({res[l]['plane_coverage']['rows'][i]['n']})" for l in labels])
    rows.append(["all >= 0.5"] + [f"{res[l]['plane_coverage']['all']['eff']:.4f} ({res[l]['plane_coverage']['all']['n']})" for l in labels])
    table("6. Plane coverage by deposited energy (MeV) per particle and plane", ["deposit (MeV)"] + labels, rows)

    prov = [l for l in labels if "provenance" in res[l]]
    if prov:
        keys = [("filled", "hits with light provenance"), ("unfilled_with_zero_pe", "unfilled hits with 0 PE"),
                ("nphotons_equals_pe", "photon count == PE after fibers"),
                ("light_eq_energy_multi", "most-light == most-energy particle (multi-contributor)"),
                ("first_eq_light_multi", "first photon from most-light particle (multi-contributor)"),
                ("light_share_below_half_multi", "most-light share < 0.5 (multi-contributor)")]
        table("7. Light provenance", ["quantity"] + prov, [[n] + [f"{res[l]['provenance'][k]:.4f}" for l in prov] for k, n in keys])
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------- plots
def style(ax, xlabel, ylabel, title):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlabel(xlabel, color=INK2, fontsize=10)
    ax.set_ylabel(ylabel, color=INK2, fontsize=10)
    ax.set_title(title, color=INK, fontsize=11, loc="left")


def plots(res, labels, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []

    def save(fig, name):
        fig.patch.set_facecolor(SURFACE)
        fig.tight_layout()
        p = os.path.join(out, name)
        fig.savefig(p, dpi=150)
        plt.close(fig)
        paths.append(p)

    def series(i):
        return dict(color=COLORS[i % 4], marker=MARKERS[i % 4], linestyle=DASHES[i % 4], linewidth=1.5, markersize=5)

    # survival
    key, xl = ("survival_vs_dx", "True path length in bar, TrueHitDx (mm)") if all("survival_vs_dx" in res[l] for l in labels) \
        else ("survival_vs_e", "True deposited energy (MeV)")
    # upper panel: survival; lower panel: difference from the reference in percentage points
    xmax = 20 if key == "survival_vs_dx" else 4
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(6.4, 5.4), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    ref = [r for r in res[labels[0]][key] if r["n"] and r["hi"] <= xmax]
    for i, l in enumerate(labels):
        rows = [r for r in res[l][key] if r["n"] and r["hi"] <= xmax]
        x = [0.5 * (r["lo"] + r["hi"]) for r in rows]
        ax.errorbar(x, [r["eff"] for r in rows], yerr=[r["err"] for r in rows], label=l, capsize=2, **series(i))
        if i:
            diff = [100 * (r["eff"] - r0["eff"]) for r, r0 in zip(rows, ref)]
            err = [100 * np.hypot(r["err"], r0["err"]) for r, r0 in zip(rows, ref)]
            axd.errorbar(x, diff, yerr=err, capsize=2, **series(i))
    ax.set_ylim(0, 1.04)
    axd.axhline(0, color="#c3c2b7", linewidth=0.8)
    style(ax, "", "Fraction above pedestal threshold", "Hit survival, muon-dominated hits")
    style(axd, xl, f"minus {labels[0]} (pp)", "")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "survival.png")

    # PE distribution at partial fill
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    hi = np.percentile(np.concatenate([res[l]["_partial_pe"] for l in labels]), 99.5)
    bins = np.arange(0, max(hi, 10) + 1, 1.0)
    for i, l in enumerate(labels):
        pe = res[l]["_partial_pe"]
        ax.hist(pe, bins=bins, histtype="step", color=COLORS[i % 4], linestyle=DASHES[i % 4], linewidth=1.5,
                density=True, label=f"{l} (mean {pe.mean():.1f}, rms {pe.std():.1f})")
    ax.axvline(THRESHOLD_PE, color=MUTED, linewidth=1, linestyle=":")
    ax.text(THRESHOLD_PE, ax.get_ylim()[1] * 0.97, f" {THRESHOLD_PE:g} PE threshold", color=MUTED, fontsize=8, va="top")
    style(ax, "Reco PE", "Fraction of hits per PE", f"Reco PE at partial fill ({res[labels[0]]['_partial_label']})")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "pe_partial_fill.png")

    # timing vs contributions
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, l in enumerate(labels):
        rows = res[l]["time_vs_nseg"]["rows"]
        x = np.arange(len(rows)) + (i - (len(labels) - 1) / 2) * 0.08
        ax.errorbar(x, [r["mean"] for r in rows], yerr=[r["sem"] for r in rows], label=l, capsize=2, **series(i))
    ax.set_xticks(range(len(NSEG_BINS)))
    ax.set_xticklabels([r["nseg"] for r in res[labels[0]]["time_vs_nseg"]["rows"]])
    ax.axhline(0, color="#c3c2b7", linewidth=0.8)
    style(ax, "Geant4 contributions to the hit (TrueNTrueParticles)", "Mean reco - true time (ns)",
          f"Hit time bias at full crossings ({res[labels[0]]['time_vs_nseg']['selection']})")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "time_bias_vs_contributions.png")

    # merge gaps
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    bins = np.linspace(-1, 6, 57)
    for i, l in enumerate(labels):
        g = np.asarray(res[l]["merging"]["gaps"])
        ax.hist(np.log10(np.maximum(g, 0.1)), bins=bins, histtype="step", color=COLORS[i % 4], linestyle=DASHES[i % 4],
                linewidth=1.5, label=f"{l} ({len(g)} groups)")
    ax.axvline(np.log10(120), color=MUTED, linewidth=1, linestyle=":")
    ax.text(np.log10(120), ax.get_ylim()[1] * 0.97, " 120 ns window", color=MUTED, fontsize=8, va="top")
    style(ax, "log10(time gap between the hits, ns)", "Groups", "Same particle, bar and slice: more than one final hit")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "merge_gaps.png")

    # plane coverage
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for i, l in enumerate(labels):
        rows = res[l]["plane_coverage"]["rows"]
        ax.plot(range(len(rows)), [r["eff"] for r in rows], label=l, **series(i))
    ax.set_xticks(range(len(PLANE_E_EDGES) - 1))
    ax.set_xticklabels([band_label(r) for r in res[labels[0]]["plane_coverage"]["rows"]])
    style(ax, "Energy deposited by the particle in the plane (MeV)", "Fraction with a surviving hit", "Plane coverage, muon-dominated particles")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "plane_coverage.png")

    # light share
    prov = [l for l in labels if "provenance" in res[l]]
    if prov:
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        for i, l in enumerate(labels):
            if l not in prov:
                continue
            ax.hist(res[l]["provenance"]["light_share_multi"], bins=np.linspace(0, 1, 41), histtype="step",
                    color=COLORS[i % 4], linestyle=DASHES[i % 4], linewidth=1.5, label=l)
        ax.set_yscale("log")
        style(ax, "Light share of the most-light particle", "Hits", "Light provenance, multi-contributor hits")
        ax.legend(frameon=False, fontsize=9)
        save(fig, "light_share.png")
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=run directory; the first is the reference")
    ap.add_argument("--out", required=True, help="output directory for report.md, metrics.json and plots")
    ap.add_argument("--files", default="1-15", help="file numbers, e.g. 1-15 or 1,4,7-9")
    args = ap.parse_args()
    files = parse_files(args.files)
    runs = dict(r.split("=", 1) for r in args.runs)
    labels = list(runs)
    os.makedirs(args.out, exist_ok=True)
    data = {}
    for l in labels:
        print(f"loading {l}: {runs[l]}", file=sys.stderr)
        data[l] = load(runs[l], files)
    has_dx = all("TrueHitDx" in data[l][0] for l in labels)
    res = {l: analyze(*data[l], has_dx) for l in labels}
    del data
    with open(os.path.join(args.out, "report.md"), "w") as f:
        f.write(report(res, labels, files))
    plots(res, labels, args.out)
    clean = {l: {k: v for k, v in res[l].items() if not k.startswith("_")} for l in labels}
    for l in labels:
        clean[l]["merging"] = {k: v for k, v in clean[l]["merging"].items() if k != "gaps"}
        if "provenance" in clean[l]:
            clean[l]["provenance"] = {k: v for k, v in clean[l]["provenance"].items() if k != "light_share_multi"}
    with open(os.path.join(args.out, "metrics.json"), "w") as f:
        json.dump(dict(runs=runs, files=files, metrics=clean), f, indent=1, default=float)
    print(f"wrote {args.out}/report.md, metrics.json and plots", file=sys.stderr)


if __name__ == "__main__":
    main()
