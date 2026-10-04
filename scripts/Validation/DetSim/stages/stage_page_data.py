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


def bin_convergence(scans):
    """Mean PE and survival vs path length for deposit bin lengths 0.5 / 1 / 2 mm (new pipeline), as the difference to
    the default 1 mm, plus the largest differences of the resegmentation scenarios."""
    cells = {"0.5": "path_new_bin0.5", "1.0": "path_new", "2.0": "path_new_bin2.0"}
    rows = {b: {p["path"]: p for p in path_scan(pd.read_csv(f"{scans}/{f}.csv"))} for b, f in cells.items()}
    out = dict(paths=[], reseg={})
    for p, ref in rows["1.0"].items():
        row = dict(path=p, n=ref["n"])
        for b in ("0.5", "2.0"):
            o = rows[b][p]
            row[b] = dict(dpe=r(100 * (o["pe"] / ref["pe"] - 1), 3),
                          dpe_err=r(100 * o["pe"] / ref["pe"] * np.hypot(o["pe_err"] / o["pe"], ref["pe_err"] / ref["pe"]), 3),
                          dsurv=r(100 * (o["surv"] - ref["surv"]), 3), dsurv_err=r(100 * np.hypot(o["surv_err"], ref["surv_err"]), 3))
        out["paths"].append(row)
    rs = {b: reseg(pd.read_csv(f"{scans}/{f}.csv")) for b, f in (("0.5", "reseg_new_bin0.5"), ("1.0", "reseg_new"), ("2.0", "reseg_new_bin2.0"))}
    for b in ("0.5", "2.0"):
        dpe, ddt, dsv, zpe, zdt = [], [], [], [], []
        for sc in rs["1.0"]:
            for a, o in zip(rs["1.0"][sc], rs[b][sc]):
                zpe.append((o["pe"] - a["pe"]) / np.hypot(o["pe_err"], a["pe_err"]))
                zdt.append((o["dt"] - a["dt"]) / np.hypot(o["dt_err"], a["dt_err"]))
                dpe.append(100 * (o["pe"] / a["pe"] - 1))
                ddt.append(o["dt"] - a["dt"])
                dsv.append(100 * (o["surv"] - a["surv"]))
        out["reseg"][b] = dict(max_dpe=r(max(abs(x) for x in dpe), 2), rms_dpe=r(float(np.sqrt(np.mean(np.square(dpe)))), 2),
                               max_ddt=r(max(abs(x) for x in ddt), 2), max_dsurv=r(max(abs(x) for x in dsv), 2), n=len(dpe),
                               max_zpe=r(max(abs(x) for x in zpe), 2), rms_zpe=r(float(np.sqrt(np.mean(np.square(zpe)))), 2),
                               max_zdt=r(max(abs(x) for x in zdt), 2), rms_zdt=r(float(np.sqrt(np.mean(np.square(zdt)))), 2))
    return out


def timing_grid(scans):
    """Timing-mode resolution vs threshold for the 9 and 16 mm crossings at the reference position and at the far end."""
    rows = []
    for t in THRESHOLDS:
        for pos, suf in (("ref", ""), ("far", "_far")):
            f = f"{scans}/path_timing_thr{t}{suf}.csv"
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f)
            for p in (9.0, 16.0):
                g = df[(df.path_mm == p) & (df.edep_MeV == df[df.path_mm == p].edep_MeV.min())]
                s = g[g.n_hits_surviving > 0]
                dt = s.hit_time - s.true_time
                q16, q84 = np.percentile(dt, [16, 84])
                rows.append(dict(thr=float(t), pos=pos, path=p, n=len(g), surv=r(g.n_hits_surviving.gt(0).mean()), pe=r(g.pe_all.mean(), 2),
                                 dt=r(dt.mean()), dt_rms=r(dt.std()), dt_hw68=r((q84 - q16) / 2)))
    return rows


def expected_readouts(n, dt, window):
    """Readouts of one channel for n crossings dt apart: a window opens at the first and takes everything within it."""
    count, start = 0, None
    for i in range(n):
        t = i * dt
        if start is None or t >= start + window:
            count, start = count + 1, t
    return count


def pileup(scans, window):
    out = {}
    scale = {"equal": [1, 1, 1, 1], "bright_first": [1, .5, .25, .125], "bright_last": [.125, .25, .5, 1]}
    for k in (LEGACY, NEW):
        f = f"{scans}/pileup_{k}.csv"
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f)
        single = df[df.n_particles == 1]
        pe1, e1 = single.pe_surviving.mean(), single.reco_energy.mean()
        rows = []
        for (pat, n, dt), g in df[df.n_particles > 1].groupby(["pattern", "n_particles", "dt"]):
            tot = sum(scale[pat][:int(n)])
            alive = g[g.n_readouts > 0]
            bright = alive[alive.top_light_share > 0]  # light provenance exists (new pipeline only)
            row = dict(pattern=pat, n=int(n), dt=float(dt), throws=len(g), readouts=r(g.n_readouts.mean(), 3),
                       expected=expected_readouts(int(n), float(dt), window),
                       pe_ratio=r(g.pe_surviving.mean() / (pe1 * tot), 4), e_ratio=r(g.reco_energy.mean() / (e1 * tot), 4),
                       lost=r((g.n_readouts == 0).mean(), 4))
            if len(bright):
                ex = max(scale[pat][:int(n)]) / tot
                row.update(top_is_max_energy=r((bright.top_light_id == bright.max_energy_id).mean(), 3),
                           top_share=r(bright.top_light_share.mean(), 3), top_share_expected=r(ex, 3),
                           first_is_earliest=r((bright.first_photon_id == bright.min_contrib_id).mean(), 3),
                           n_contrib=r(bright.n_contrib.mean(), 2))
            rows.append(row)
        out[k] = dict(single_pe=r(pe1, 2), single_e=r(e1, 3), rows=rows)
    return out


def contribution_stats(readout_root, tree="TMS"):
    """Consistency of the per-particle light breakdown in a ConvertToTMSTree readout file."""
    import uproot
    cols = ["NTrueHits", "TrueHitNPhotons", "TrueHitPrimaryIdByLight", "TrueHitLightShare", "TrueHitLightContribOffset",
            "TrueHitNLightContrib", "TrueHitLightContribPrimaryId", "TrueHitLightContribPhotons", "TrueHitLightContribShare"]
    a = uproot.open(readout_root)[tree].arrays(cols, library="np")
    n_hits = n_light = bad_sum = multi = bad_top = 0
    counts = {}
    for ev in range(len(a["NTrueHits"])):
        for i in range(int(a["NTrueHits"][ev])):
            n_hits += 1
            nc = int(a["TrueHitNLightContrib"][ev][i])
            if nc == 0:
                continue
            n_light += 1
            off = int(a["TrueHitLightContribOffset"][ev][i])
            ph = a["TrueHitLightContribPhotons"][ev][off:off + nc]
            sh = a["TrueHitLightContribShare"][ev][off:off + nc]
            pid = a["TrueHitLightContribPrimaryId"][ev][off:off + nc]
            key = str(nc) if nc < 4 else "4+"
            counts[key] = counts.get(key, 0) + 1
            if ph.sum() != a["TrueHitNPhotons"][ev][i] or abs(sh.sum() - 1) > 1e-4:
                bad_sum += 1
            if nc > 1:
                multi += 1
                j = int(np.argmax(ph))
                if abs(sh[j] - a["TrueHitLightShare"][ev][i]) > 1e-4 or pid[j] != a["TrueHitPrimaryIdByLight"][ev][i]:
                    bad_top += 1
    return dict(hits=n_hits, with_light=n_light, by_n=counts, multi=multi, bad_sum=bad_sum, bad_top=bad_top)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scans", required=True)
    ap.add_argument("--spills")
    ap.add_argument("--files", default="1-20")
    ap.add_argument("--metrics")
    ap.add_argument("--contrib", help="a ConvertToTMSTree *_Readout.root file, for the light-breakdown consistency check")
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
    if os.path.exists(f"{S}/path_new_bin0.5.csv"):
        page["bins"] = bin_convergence(S)
    page["timing_grid"] = timing_grid(S)
    page["pileup"] = pileup(S, page["readout"]["new"]["ReadoutTime"])
    if a.contrib:
        page["contrib"] = contribution_stats(a.contrib)
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
