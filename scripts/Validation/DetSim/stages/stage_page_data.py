#!/usr/bin/env python3
"""Reduce the stage scans (run_stage_scans.sh) and, optionally, a real-spill run to the data of the stage-by-stage
detector-simulation page (page_data.json, read by build_stage_page.py).

    python3 stage_page_data.py --scans <scan dir> [--spills <flag-on run dir> --files 1-20]
                               [--metrics <detsim_validation.py metrics.json>] --out page_data.json

--spills: a ConvertToTMSTree run with the new pipeline (light-provenance branches in Truth_Info), for stage 6.
--metrics: metrics.json of detsim_validation.py over a legacy/new pair of real-spill runs, for the real-spill
cross-checks drawn next to the controlled scans.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

LEGACY, NEW = "legacy", "new"
THRESHOLDS = ["0.5", "1.0", "1.5", "2.0", "2.5", "3.0"]
FULL_PATH = 16.0  # mm, the full thickness of the bar


def sem(x):
    return float(np.std(x, ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0


def r(x, nd=4):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


def path_scan(df):
    """Per-path rows at minimum-ionizing dE/dx (the brighter full crossings are left out)."""
    mip = df[np.isclose(df.edep_MeV / df.path_mm, df.edep_MeV.iloc[0] / df.path_mm.iloc[0], rtol=1e-3)]
    rows = []
    for p, g in mip.groupby("path_mm"):
        surv = g.n_hits_surviving > 0
        s = g[surv]
        dt = s.hit_time - s.true_time
        rows.append(dict(path=r(p), n=len(g), surv=r(surv.mean()), surv_err=r(np.sqrt(surv.mean() * (1 - surv.mean()) / len(g))),
                         pe=r(g.pe_all.mean()), pe_err=r(sem(g.pe_all)), var=r(g.pe_all.var()),
                         fano=r(g.pe_all.var() / g.pe_all.mean()) if g.pe_all.mean() > 0 else None,
                         nph=r(g.n_photons.clip(lower=0).mean()) if (g.n_photons >= 0).any() else None,
                         dt=r(dt.mean()) if len(s) > 20 else None, dt_rms=r(dt.std()) if len(s) > 20 else None,
                         dt_err=r(sem(dt)) if len(s) > 20 else None))
    return rows


def position_scan(df):
    """Full MIP crossing vs distance from the readout end (DetSimStageScan position)."""
    rows = []
    for d, g in df.groupby("distance_mm"):
        surv = g.n_hits_surviving > 0
        s = g[surv]
        dt = s.hit_time - s.true_time
        rows.append(dict(d=r(d, 1), n=len(g), bar=f"{g.bar_type.iloc[0]}{int(g.bar_number.iloc[0])}", length=r(g.bar_length_mm.iloc[0], 0),
                         pe=r(g.pe_all.mean()), pe_err=r(sem(g.pe_all)), surv=r(surv.mean()),
                         surv_err=r(np.sqrt(surv.mean() * (1 - surv.mean()) / len(g))),
                         dt=r(dt.mean()), dt_err=r(sem(dt)), dt_rms=r(dt.std())))
    return rows


def time_vs_pe(df, edges):
    """Hit time - true time of surviving throws, binned in the light of the surviving hit (all scan cells)."""
    s = df[df.n_hits_surviving > 0]
    dt = s.hit_time - s.true_time
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s.pe_surviving >= lo) & (s.pe_surviving < hi)
        if m.sum() < 30:
            continue
        rows.append(dict(lo=lo, hi=hi, n=int(m.sum()), x=r(s.pe_surviving[m].mean(), 2), dt=r(dt[m].mean()),
                         dt_err=r(sem(dt[m])), dt_rms=r(dt[m].std())))
    return rows


def reseg(df):
    out = {}
    for sc, gs in df.groupby("scenario"):
        rows = []
        for nseg, g in gs.groupby("nseg"):
            surv = g.n_hits_surviving > 0
            dt = (g.min_hit_time - g.true_time)[surv]
            rows.append(dict(nseg=int(nseg), n=len(g), surv=r(surv.mean()), pe=r(g.total_pe.mean()), pe_err=r(sem(g.total_pe)),
                             hits=r(g.n_hits_surviving[surv].mean()), hits_total=r(g.n_hits_total.mean()),
                             efrac=r((g.true_energy_in_hits / g.injected_energy).mean(), 5),
                             dt=r(dt.mean()), dt_err=r(sem(dt)), dt_rms=r(dt.std())))
        out[sc] = rows
    return out


def pairs(scans):
    out = {}
    for name in ["legacy", "new", "new_dead500", "new_dead500_zombie100"]:
        f = f"{scans}/pair_{name}.csv"
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f)
        same = df[df.case == "same_bar"].groupby("dt").n_readouts
        out[name] = dict(
            scan=[dict(dt=float(dt), mean=r(g.mean()), two=r((g == 2).mean())) for dt, g in same],
            xbar=dict(readouts=r(df[df.case == "xbar_halves"].n_readouts.mean()),
                      channels=r(df[df.case == "xbar_halves"].n_channels.mean())),
            single=r(df[df.case == "single"].n_readouts.mean()))
    return out


def readout_params(conf):
    """ReadoutTime, Deadtime, ZombieTime of a pinned readout config, for the expected-readouts curves."""
    vals = {}
    for line in open(conf):
        k, _, v = line.partition("=")
        k = k.strip()
        if k in ("ReadoutTime", "Deadtime", "ZombieTime"):
            vals[k] = float(v.split("#")[0])
    return vals


def timing(scans):
    out = dict(thresholds=[], survival={}, tot={})
    for t in THRESHOLDS:
        f = f"{scans}/path_timing_thr{t}.csv"
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f)
        full = df[(df.path_mm == FULL_PATH) & (df.edep_MeV == df[df.path_mm == FULL_PATH].edep_MeV.min())]
        s = full[full.n_hits_surviving > 0]
        dt = s.hit_time - s.true_time
        q16, q84 = np.percentile(dt, [16, 84])
        out["thresholds"].append(dict(thr=float(t), n=len(full), surv=r(full.n_hits_surviving.gt(0).mean()),
                                      dt=r(dt.mean()), dt_rms=r(dt.std()), dt_hw68=r((q84 - q16) / 2),
                                      on_grid=r(np.mean(np.isclose(np.mod(s.hit_time * 2, 1), 0) | np.isclose(np.mod(s.hit_time * 2, 1), 1)))))
        out["survival"][t] = [dict(path=p["path"], surv=p["surv"]) for p in path_scan(df)]
        g = df[(df.n_hits_surviving > 0) & (df.n_photons >= 0) & (df.tot > 0)]
        rows = []
        for lo in range(0, 100, 4):
            m = (g.n_photons >= lo) & (g.n_photons < lo + 4)
            if m.sum() < 30:
                continue
            rows.append(dict(lo=lo, hi=lo + 4, n=int(m.sum()), x=r(g.n_photons[m].mean(), 2), tot=r(g.tot[m].mean()),
                             q16=r(np.percentile(g.tot[m], 16)), q84=r(np.percentile(g.tot[m], 84))))
        out["tot"][t] = rows
    return out


def provenance(run_dir, files):
    import awkward as ak
    import uproot
    cols = ["NTrueHits", "TrueNTrueParticles", "TrueHitPrimaryId", "TrueHitVertexId", "TrueRecoHitIsPedSupped",
            "TrueHitNPhotons", "TrueHitPrimaryIdByLight", "TrueHitVertexIdByLight", "TrueHitLightShare",
            "TrueHitFirstPhotonPrimaryId", "TrueHitFirstPhotonVertexId"]
    parts = {c: [] for c in cols[1:]}
    for n in files:
        with uproot.open(f"{run_dir}/{n:07d}_RecoCandidates.root") as f:
            a = f["Truth_Info"].arrays(cols, library="ak")
        mask = ak.local_index(a["TrueHitPrimaryId"]) < a["NTrueHits"]
        for c in cols[1:]:
            parts[c].append(ak.to_numpy(ak.flatten(a[c][mask])))
    d = {c: np.concatenate(v) for c, v in parts.items()}
    sel = (d["TrueHitNPhotons"] > 0) & (d["TrueNTrueParticles"] > 1) & ~d["TrueRecoHitIsPedSupped"].astype(bool)
    share = d["TrueHitLightShare"][sel]
    same = ((d["TrueHitPrimaryIdByLight"] == d["TrueHitPrimaryId"]) & (d["TrueHitVertexIdByLight"] == d["TrueHitVertexId"]))[sel]
    first = ((d["TrueHitFirstPhotonPrimaryId"] == d["TrueHitPrimaryIdByLight"])
             & (d["TrueHitFirstPhotonVertexId"] == d["TrueHitVertexIdByLight"]))[sel]
    edges = np.round(np.linspace(0.3, 1.0, 15), 3)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (share >= lo) & ((share < hi) | (hi == 1.0))
        if m.sum() < 30:
            continue
        rows.append(dict(lo=float(lo), hi=float(hi), n=int(m.sum()), frac=r(m.mean()),
                         same=r(same[m].mean()), first=r(first[m].mean())))
    return dict(files=f"{files[0]}-{files[-1]}", n_multi=int(sel.sum()),
                n_surviving=int((~d["TrueRecoHitIsPedSupped"].astype(bool) & (d["TrueHitNPhotons"] > 0)).sum()),
                same=r(same.mean()), first=r(first.mean()), below_half=r((share < 0.5).mean()), bins=rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scans", required=True)
    ap.add_argument("--spills")
    ap.add_argument("--files", default="1-20")
    ap.add_argument("--metrics")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S = a.scans
    page = dict(provenance=open(f"{S}/provenance.txt").read().strip())
    paths = {k: pd.read_csv(f"{S}/path_{k}.csv") for k in (LEGACY, NEW)}
    page["reseg"] = {k: reseg(pd.read_csv(f"{S}/reseg_{k}.csv")) for k in (LEGACY, NEW)}
    page["path"] = {k: path_scan(v) for k, v in paths.items()}
    pe_edges = [3, 5, 7, 9, 11, 13, 15, 18, 21, 25, 30, 36, 42, 50, 60, 70, 80, 90, 105]
    page["time_vs_pe"] = {k: time_vs_pe(v, pe_edges) for k, v in paths.items()}
    if all(os.path.exists(f"{S}/position_{k}.csv") for k in (LEGACY, NEW)):
        page["position"] = {k: position_scan(pd.read_csv(f"{S}/position_{k}.csv")) for k in (LEGACY, NEW)}
    page["pairs"] = pairs(S)
    page["readout"] = {k: readout_params(f"{S}/configs/readout_{k}.toml") for k in page["pairs"]}
    page["timing"] = timing(S)
    if a.spills:
        lo, _, hi = a.files.partition("-")
        page["provenance_spills"] = provenance(a.spills, list(range(int(lo), int(hi or lo) + 1)))
    if a.metrics:
        m = json.load(open(a.metrics))["metrics"]
        page["spills"] = {k: dict(survival_vs_dx=v["survival_vs_dx"], pe_by_band=v["pe_by_band"], time_vs_nseg=v["time_vs_nseg"],
                                  spills=v["totals"]["spills"]) for k, v in m.items()}
    json.dump(page, open(a.out, "w"), separators=(",", ":"))
    print("wrote", a.out, os.path.getsize(a.out), "bytes")


if __name__ == "__main__":
    main()
