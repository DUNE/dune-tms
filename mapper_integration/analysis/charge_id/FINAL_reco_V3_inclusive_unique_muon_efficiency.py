import ROOT
import numpy as np
import math
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

FN = (
    "/exp/dune/app/users/sushils/dune-tms/bin/"
    "neutrino.1001_1005.merged.edep_TMS_RecoCandidates_Hough_Cluster1.root"
)

MAP = (
    "/exp/dune/app/users/sushils/TMS_MAPPER_RUN/"
    "Mapper_global_DIAGNOSTIC.txt"
)

OUT = (
    "/exp/dune/app/users/sushils/AUTO_SANDBOX/"
    "RECO_FIELDAWARE_SIGNED_DISTANCE"
)

os.makedirs(OUT, exist_ok=True)

FIELD_STEP_MM = 10.0

# V2:
# estimate the entrance x-z direction from the first upstream
# reconstructed hits instead of using StartDirection directly.
LOCAL_FIT_NHITS = 6


# ============================================================
# LOAD REALISTIC GLOBAL FIELD MAP
# ============================================================

print("Loading mapper:", MAP)

with open(MAP, "r") as fh:

    metadata = None

    for line in fh:
        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        metadata = np.fromstring(line, sep=" ")
        break

    if metadata is None or len(metadata) != 6:
        raise RuntimeError(
            "Could not find six-column mapper metadata row"
        )

    x0, y0, z0, dx, dy, dz = metadata

    data = np.loadtxt(
        fh,
        comments="#",
        usecols=(0,1,2,3,4,5),
        dtype=np.float32
    )


xmin = float(np.min(data[:,0]))
xmax = float(np.max(data[:,0]))

ymin = float(np.min(data[:,1]))
ymax = float(np.max(data[:,1]))

zmin = float(np.min(data[:,2]))
zmax = float(np.max(data[:,2]))

nx = int(round((xmax-x0)/dx)) + 1
ny = int(round((ymax-y0)/dy)) + 1
nz = int(round((zmax-z0)/dz)) + 1

expected = nx * ny * nz

if len(data) != expected:
    raise RuntimeError(
        f"Mapper rows={len(data)} "
        f"but expected={expected}"
    )

Bgrid = data[:,3:6].reshape(
    nx,
    ny,
    nz,
    3
)

print("Mapper origin [mm] :", x0, y0, z0)
print("Mapper spacing [mm]:", dx, dy, dz)
print("Mapper dimensions   :", nx, ny, nz)


# ============================================================
# FIELD INTERPOLATION
# ============================================================

def field_at(x, y, z):

    if (
        x < xmin or x > xmax or
        y < ymin or y > ymax or
        z < zmin or z > zmax
    ):
        return np.zeros(3, dtype=np.float64)

    fx = (x-x0)/dx
    fy = (y-y0)/dy
    fz = (z-z0)/dz

    ix = int(math.floor(fx))
    iy = int(math.floor(fy))
    iz = int(math.floor(fz))

    tx = fx - ix
    ty = fy - iy
    tz = fz - iz

    if ix >= nx-1:
        ix = nx-2
        tx = 1.0

    if iy >= ny-1:
        iy = ny-2
        ty = 1.0

    if iz >= nz-1:
        iz = nz-2
        tz = 1.0

    c000 = Bgrid[ix,   iy,   iz]
    c001 = Bgrid[ix,   iy,   iz+1]
    c010 = Bgrid[ix,   iy+1, iz]
    c011 = Bgrid[ix,   iy+1, iz+1]

    c100 = Bgrid[ix+1, iy,   iz]
    c101 = Bgrid[ix+1, iy,   iz+1]
    c110 = Bgrid[ix+1, iy+1, iz]
    c111 = Bgrid[ix+1, iy+1, iz+1]

    c00 = c000*(1-tx) + c100*tx
    c01 = c001*(1-tx) + c101*tx

    c10 = c010*(1-tx) + c110*tx
    c11 = c011*(1-tx) + c111*tx

    c0 = c00*(1-ty) + c10*ty
    c1 = c01*(1-ty) + c11*ty

    return np.asarray(
        c0*(1-tz) + c1*tz,
        dtype=np.float64
    )


# ============================================================
# ROOT FILE
# ============================================================

f = ROOT.TFile.Open(FN)

if not f or f.IsZombie():
    raise RuntimeError("Could not open ROOT file")

r = f.Get("Reco_Tree")
t = f.Get("Truth_Info")

if not r:
    raise RuntimeError("Reco_Tree not found")

if not t:
    raise RuntimeError("Truth_Info not found")

print("Reco entries :", r.GetEntries())
print("Truth entries:", t.GetEntries())

if r.GetEntries() != t.GetEntries():
    raise RuntimeError(
        "Reco_Tree and Truth_Info entry counts differ"
    )



# ============================================================
# V3 — INCLUSIVE UNIQUE-MUON CHARGE-ID ANALYSIS
# ============================================================
#
# Truth definition:
#   1. Match MuonP4 to BirthMomentum[nTrueParticles]
#   2. Use matched particle PDG (NOT LeptonPDG)
#   3. Deduplicate repeated Truth_Info slices
#   4. Require true muon BirthPosition in current ND-LAr fiducial
#   5. Require valid truth MomentumTMSStart
#
# Reconstruction:
#   - search all slices belonging to the unique true muon
#   - truth-match reco track via RecoTrackPrimaryParticleIndex
#   - choose highest-quality VALID V2 candidate
#       quality = Length_3D, then nHits
#   - use exact frozen V2 six-hit local-fit + Kx algorithm
#
# Report:
#   SYSTEM efficiency =
#       correct charge / all unique true muons entering TMS
#
#   CONDITIONAL charge-ID efficiency =
#       correct charge / unique entering muons with usable V2 reco
#
# ============================================================

from collections import Counter, defaultdict

MU_MASS_MEV = 105.6583755

# ============================================================
# CURRENT 1.1.0 LAr FIDUCIAL VOLUME [mm]
#
# Active LAr:
# x = -3478.48 ... +3478.48
# y = -2166.71 ... +829.282
# z = +4179.24 ... +9135.88
#
# Current config:
# XYCut = 500 mm
# DownstreamZCut = 1500 mm
# ============================================================

LAR_XMIN = -3478.48 + 500.0
LAR_XMAX =  3478.48 - 500.0

LAR_YMIN = -2166.71 + 500.0
LAR_YMAX =   829.282 - 500.0

LAR_ZMIN =  4179.24 + 500.0
LAR_ZMAX =  9135.88 - 1500.0


OUTDIR_V3 = (
    "/exp/dune/app/users/sushils/AUTO_SANDBOX/"
    "RECO_V3_INCLUSIVE_UNIQUE_MUON"
)

os.makedirs(OUTDIR_V3, exist_ok=True)


# ============================================================
# GENERIC ROOT HELPERS
# ============================================================

def get4(arr, i=None):

    if i is None:
        return [
            float(arr[k])
            for k in range(4)
        ]

    try:
        return [
            float(arr[i][k])
            for k in range(4)
        ]

    except Exception:
        return [
            float(arr[i*4+k])
            for k in range(4)
        ]


def valid4(v):

    return (
        all(
            math.isfinite(x)
            for x in v
        )
        and
        all(
            abs(x) < 1.0e7
            for x in v[:3]
        )
        and
        math.isfinite(v[3])
        and
        abs(v[3]) < 1.0e8
    )


def dist4(a, b):

    return math.sqrt(
        sum(
            (a[k]-b[k])**2
            for k in range(4)
        )
    )


def p3mag(v):

    return math.sqrt(
        v[0]**2 +
        v[1]**2 +
        v[2]**2
    )


def in_lar_fiducial(p):

    return (
        LAR_XMIN <= p[0] <= LAR_XMAX
        and
        LAR_YMIN <= p[1] <= LAR_YMAX
        and
        LAR_ZMIN <= p[2] <= LAR_ZMAX
    )


# ============================================================
# TRUE-MUON MATCHING
# ============================================================

def find_true_muon_index(truth):

    mu4 = get4(
        truth.MuonP4
    )

    if not valid4(mu4):
        return None

    best_j = -1
    best_d = float("inf")

    for j in range(
        int(truth.nTrueParticles)
    ):

        birth4 = get4(
            truth.BirthMomentum,
            j
        )

        if not valid4(birth4):
            continue

        d = dist4(
            birth4,
            mu4
        )

        if d < best_d:
            best_d = d
            best_j = j

    if best_j < 0:
        return None

    # Audited: matched values are essentially exact.
    if best_d > 1.0e-2:
        return None

    pdg = int(
        truth.PDG[best_j]
    )

    if abs(pdg) != 13:
        return None

    return best_j


def make_unique_muon_key(
    truth,
    mu_index
):

    pdg = int(
        truth.PDG[mu_index]
    )

    mu4 = get4(
        truth.MuonP4
    )

    birth = get4(
        truth.BirthPosition,
        mu_index
    )

    try:
        vertex_id = int(
            truth.VertexID[
                mu_index
            ]
        )
    except Exception:
        vertex_id = -999999

    track_id = int(
        truth.TrackId[
            mu_index
        ]
    )

    # EventNo deliberately excluded:
    # repeated slices can have different EventNo.
    return (
        int(truth.RunNo),
        int(truth.SpillNo),
        vertex_id,
        track_id,
        pdg,

        round(birth[0], 3),
        round(birth[1], 3),
        round(birth[2], 3),

        round(mu4[0], 3),
        round(mu4[1], 3),
        round(mu4[2], 3),
        round(mu4[3], 3),
    )


# ============================================================
# EXACT FROZEN V2 RECO SIGNED-DISTANCE ESTIMATOR
# ============================================================

def reco_track_length3d(reco, j):

    try:
        return float(
            reco.Length_3D[j]
        )
    except Exception:
        return float(
            reco.Length[j]
        )


def compute_v2_for_track(
    reco,
    j
):
    """
    Exact physics logic copied from frozen V2:
      StartPos
      EndPos
      StartDirection sanity/orientation
      TrackHitPos
      first-six-hit x(z) fit
      same reconstructed path
      same Kx field integral
      same field-aware SD definition
    """

    nreco = int(
        reco.nTracks
    )

    if j < 0 or j >= nreco:
        return None

    start = np.frombuffer(
        reco.StartPos,
        dtype=np.float32,
        count=nreco*4
    ).reshape(
        nreco,
        4
    )

    end = np.frombuffer(
        reco.EndPos,
        dtype=np.float32,
        count=nreco*4
    ).reshape(
        nreco,
        4
    )

    direction = np.frombuffer(
        reco.StartDirection,
        dtype=np.float32,
        count=nreco*3
    ).reshape(
        nreco,
        3
    )

    nhits = np.frombuffer(
        reco.nHits,
        dtype=np.int32,
        count=nreco
    )

    hits_all = np.frombuffer(
        reco.TrackHitPos,
        dtype=np.float32,
        count=nreco*200*4
    ).reshape(
        nreco,
        200,
        4
    )

    p0 = np.asarray(
        start[j,:3],
        dtype=np.float64
    )

    p3 = np.asarray(
        end[j,:3],
        dtype=np.float64
    )

    d0 = np.asarray(
        direction[j,:3],
        dtype=np.float64
    )

    if not (
        np.all(np.isfinite(p0))
        and
        np.all(np.isfinite(p3))
        and
        np.all(np.isfinite(d0))
    ):
        return None

    x1 = p0[0]
    z1 = p0[2]

    x3 = p3[0]
    z3 = p3[2]

    dxdir = d0[0]
    dzdir = d0[2]

    # Same orientation requirement as frozen V2.
    if z3 <= z1:
        return None

    if dzdir <= 1.0e-6:
        return None

    # --------------------------------------------------------
    # V1 quantity preserved only as diagnostic
    # --------------------------------------------------------

    slope_v1 = (
        dxdir /
        dzdir
    )

    x_extrap_v1 = (
        x1 +
        slope_v1 * (z3-z1)
    )

    sd_raw_v1 = (
        x3 -
        x_extrap_v1
    )

    # --------------------------------------------------------
    # RECONSTRUCTED HIT PATH
    # --------------------------------------------------------

    nh = int(
        nhits[j]
    )

    if nh < 2:
        return None

    nh = min(
        nh,
        200
    )

    hits = np.asarray(
        hits_all[
            j,
            :nh,
            :3
        ],
        dtype=np.float64
    )

    good = (
        np.all(
            np.isfinite(hits),
            axis=1
        )
        &
        np.all(
            np.abs(hits) < 1.0e7,
            axis=1
        )
    )

    hits = hits[good]

    if len(hits) < 2:
        return None

    # TrackHitPos observed downstream -> upstream.
    # Orient robustly relative to StartPos.
    d_first = np.linalg.norm(
        hits[0] - p0
    )

    d_last = np.linalg.norm(
        hits[-1] - p0
    )

    if d_last < d_first:
        hits = hits[::-1]

    # --------------------------------------------------------
    # V2 LOCAL ENTRANCE x(z) FIT
    # --------------------------------------------------------

    nfit = min(
        LOCAL_FIT_NHITS,
        len(hits)
    )

    if nfit < 3:
        return None

    fit_hits = np.asarray(
        hits[:nfit],
        dtype=np.float64
    )

    zfit = fit_hits[:,2]
    xfit = fit_hits[:,0]

    if np.ptp(zfit) < 100.0:
        return None

    slope_v2, intercept_v2 = (
        np.polyfit(
            zfit,
            xfit,
            1
        )
    )

    if not (
        math.isfinite(slope_v2)
        and
        math.isfinite(intercept_v2)
    ):
        return None

    x_extrap_v2 = (
        x1 +
        slope_v2 * (z3-z1)
    )

    sd_raw = (
        x3 -
        x_extrap_v2
    )

    # --------------------------------------------------------
    # Full reconstructed path
    # --------------------------------------------------------

    path = np.vstack([
        p0,
        hits,
        p3
    ])

    clean = [
        path[0]
    ]

    for point in path[1:]:

        if np.linalg.norm(
            point-clean[-1]
        ) > 1.0e-3:

            clean.append(
                point
            )

    path = np.asarray(
        clean
    )

    if len(path) < 2:
        return None

    # --------------------------------------------------------
    # SAME FIELD-AWARE Kx INTEGRAL AS FROZEN V2
    # --------------------------------------------------------

    Kx = 0.0

    for iseg in range(
        len(path)-1
    ):

        a = path[iseg]
        b = path[iseg+1]

        delta = b-a

        seglen = np.linalg.norm(
            delta
        )

        if seglen < 1.0e-6:
            continue

        u = (
            delta /
            seglen
        )

        uy = u[1]
        uz = u[2]

        nsub = max(
            1,
            int(
                math.ceil(
                    seglen /
                    FIELD_STEP_MM
                )
            )
        )

        dsub = (
            seglen /
            nsub
        )

        for k in range(
            nsub
        ):

            frac = (
                k + 0.5
            ) / nsub

            p = (
                a +
                frac*delta
            )

            lever = (
                z3 -
                p[2]
            )

            if lever <= 0:
                continue

            B = field_at(
                p[0],
                p[1],
                p[2]
            )

            By = B[1]
            Bz = B[2]

            orientation = (
                uz*By -
                uy*Bz
            )

            Kx += (
                lever *
                orientation *
                dsub
            )

    if not math.isfinite(Kx):
        return None

    if abs(Kx) < 1.0e-9:
        return None

    sd_field = (
        np.sign(Kx) *
        sd_raw
    )

    sd_field_v1 = (
        np.sign(Kx) *
        sd_raw_v1
    )

    return {
        "sd_field": float(
            sd_field
        ),
        "sd_raw": float(
            sd_raw
        ),
        "sd_field_v1": float(
            sd_field_v1
        ),
        "Kx": float(
            Kx
        ),
        "nHits": int(
            nhits[j]
        ),
        "length3d": float(
            reco_track_length3d(
                reco,
                j
            )
        ),
        "track_index": int(j),
    }


# ============================================================
# PASS 1:
# BUILD UNIQUE TRUTH-MUON DENOMINATOR
# ============================================================

unique_muons = {}

truth_entries_total = int(
    t.GetEntries()
)

for ievt in range(
    truth_entries_total
):

    t.GetEntry(
        ievt
    )

    mu_index = (
        find_true_muon_index(
            t
        )
    )

    if mu_index is None:
        continue

    pdg = int(
        t.PDG[
            mu_index
        ]
    )

    birth = get4(
        t.BirthPosition,
        mu_index
    )

    if not valid4(birth):
        continue

    if not in_lar_fiducial(
        birth
    ):
        continue

    ptms4 = get4(
        t.MomentumTMSStart,
        mu_index
    )

    if not valid4(ptms4):
        continue

    ptms = p3mag(
        ptms4
    )

    if not (
        math.isfinite(ptms)
        and
        ptms > 0.0
    ):
        continue

    mu4 = get4(
        t.MuonP4
    )

    ke = (
        float(mu4[3]) -
        MU_MASS_MEV
    )

    if not (
        math.isfinite(ke)
        and
        ke > 0.0
    ):
        continue

    key = (
        make_unique_muon_key(
            t,
            mu_index
        )
    )

    if key not in unique_muons:

        unique_muons[key] = {
            "pdg": pdg,
            "ke": ke,
            "ptms": ptms,
            "birth": tuple(
                birth[:3]
            ),
            "slice_entries": [],
            "candidate_count": 0,
            "valid_candidate_count": 0,
            "best": None,
        }

    unique_muons[
        key
    ][
        "slice_entries"
    ].append(
        int(ievt)
    )


# ============================================================
# PASS 2:
# SEARCH ALL SLICES FOR TRUTH-MATCHED RECO CANDIDATES
#
# Candidate choice is independent of charge correctness:
#   highest Length_3D,
#   then highest nHits,
#   then lowest entry number / track index.
#
# We do NOT pick the candidate whose sign agrees with truth.
# ============================================================

for key, rec in (
    unique_muons.items()
):

    candidates = []

    for ievt in rec[
        "slice_entries"
    ]:

        t.GetEntry(
            ievt
        )

        r.GetEntry(
            ievt
        )

        mu_index = (
            find_true_muon_index(
                t
            )
        )

        if mu_index is None:
            continue

        nreco = int(
            r.nTracks
        )

        ntruth_reco = int(
            t.RecoTrackN
        )

        # Same structural sanity condition used by frozen V2.
        if nreco <= 0:
            continue

        if nreco != ntruth_reco:
            continue

        for j in range(
            nreco
        ):

            try:
                primary_index = int(
                    t.
                    RecoTrackPrimaryParticleIndex[
                        j
                    ]
                )
            except Exception:
                continue

            if primary_index != mu_index:
                continue

            rec[
                "candidate_count"
            ] += 1

            result = (
                compute_v2_for_track(
                    r,
                    j
                )
            )

            if result is None:
                continue

            rec[
                "valid_candidate_count"
            ] += 1

            result[
                "entry"
            ] = int(
                ievt
            )

            candidates.append(
                result
            )

    if candidates:

        # Deterministic and charge-independent.
        candidates.sort(
            key=lambda x: (
                -x[
                    "length3d"
                ],
                -x[
                    "nHits"
                ],
                x[
                    "entry"
                ],
                x[
                    "track_index"
                ],
            )
        )

        rec[
            "best"
        ] = candidates[0]


# ============================================================
# CLASSIFICATION
# ============================================================

for rec in (
    unique_muons.values()
):

    best = rec[
        "best"
    ]

    if best is None:

        rec[
            "reco_available"
        ] = False

        rec[
            "correct"
        ] = False

        rec[
            "sd"
        ] = float(
            "nan"
        )

        continue

    rec[
        "reco_available"
    ] = True

    sd = float(
        best[
            "sd_field"
        ]
    )

    rec[
        "sd"
    ] = sd

    pdg = int(
        rec[
            "pdg"
        ]
    )

    # Same convention validated in V2:
    # mu- (PDG +13) should have positive field-aware SD.
    # mu+ (PDG -13) should have negative field-aware SD.
    rec[
        "correct"
    ] = bool(
        (
            pdg == 13
            and
            sd > 0.0
        )
        or
        (
            pdg == -13
            and
            sd < 0.0
        )
    )


# ============================================================
# SUMMARY HELPERS
# ============================================================

def species_rows(
    pdg
):

    return [
        rec
        for rec in
        unique_muons.values()
        if int(
            rec[
                "pdg"
            ]
        ) == pdg
    ]


def print_species_summary(
    pdg,
    label
):

    rows = species_rows(
        pdg
    )

    total = len(
        rows
    )

    reco_rows = [
        x
        for x in rows
        if x[
            "reco_available"
        ]
    ]

    correct_rows = [
        x
        for x in rows
        if x[
            "correct"
        ]
    ]

    nreco = len(
        reco_rows
    )

    ncorrect = len(
        correct_rows
    )

    print("\n" + "="*78)
    print(label)
    print("="*78)

    print(
        "unique true entering      :",
        total
    )

    print(
        "usable V2 reco candidate  :",
        nreco
    )

    print(
        "no usable V2 candidate    :",
        total-nreco
    )

    print(
        "correct charge            :",
        ncorrect
    )

    print(
        "wrong charge among reco   :",
        nreco-ncorrect
    )

    system_eff = (
        ncorrect/total
        if total
        else float("nan")
    )

    conditional_eff = (
        ncorrect/nreco
        if nreco
        else float("nan")
    )

    reco_fraction = (
        nreco/total
        if total
        else float("nan")
    )

    print(
        "reco availability         :",
        reco_fraction
    )

    print(
        "SYSTEM efficiency         :",
        system_eff
    )

    print(
        "CHARGE|RECO efficiency    :",
        conditional_eff
    )

    return {
        "total": total,
        "reco": nreco,
        "correct": ncorrect,
        "system": system_eff,
        "conditional": conditional_eff,
        "reco_fraction": reco_fraction,
    }


# ============================================================
# BINNED EFFICIENCY
# ============================================================

def binned_metrics(
    rows,
    edges,
    variable
):

    out = []

    for ibin in range(
        len(edges)-1
    ):

        lo = float(
            edges[ibin]
        )

        hi = float(
            edges[ibin+1]
        )

        selected = [
            x
            for x in rows
            if (
                lo <= float(
                    x[
                        variable
                    ]
                ) < hi
            )
        ]

        ntotal = len(
            selected
        )

        if ntotal == 0:
            continue

        reco_rows = [
            x
            for x in selected
            if x[
                "reco_available"
            ]
        ]

        nreco = len(
            reco_rows
        )

        ncorrect = sum(
            1
            for x in selected
            if x[
                "correct"
            ]
        )

        system = (
            ncorrect /
            ntotal
        )

        conditional = (
            ncorrect /
            nreco
            if nreco
            else float("nan")
        )

        reco_fraction = (
            nreco /
            ntotal
        )

        out.append({
            "lo": lo,
            "hi": hi,
            "x": 0.5*(
                lo+hi
            ),
            "total": ntotal,
            "reco": nreco,
            "correct": ncorrect,
            "system": system,
            "conditional": conditional,
            "reco_fraction": reco_fraction,
        })

    return out


def print_table(
    title,
    rows,
    unit
):

    print("\n" + "="*100)
    print(title)
    print("="*100)

    print(
        "range                  "
        "correct / total   "
        "correct / reco    "
        "reco / total"
    )

    for x in rows:

        cond = (
            f"{x['conditional']:.4f}"
            if math.isfinite(
                x[
                    "conditional"
                ]
            )
            else "nan"
        )

        print(
            f"{x['lo']:7.0f} - "
            f"{x['hi']:7.0f} {unit:<6s}  "
            f"{x['correct']:3d}/{x['total']:3d} "
            f"= {x['system']:.4f}     "
            f"{x['correct']:3d}/{x['reco']:3d} "
            f"= {cond:>6s}     "
            f"{x['reco']:3d}/{x['total']:3d} "
            f"= {x['reco_fraction']:.4f}"
        )


# ============================================================
# PRINT CORE V3 RESULTS
# ============================================================

print("\n")
print("="*78)
print("V3 INCLUSIVE UNIQUE-MUON ANALYSIS")
print("="*78)

print(
    "LAr fiducial x [mm]:",
    LAR_XMIN,
    LAR_XMAX
)

print(
    "LAr fiducial y [mm]:",
    LAR_YMIN,
    LAR_YMAX
)

print(
    "LAr fiducial z [mm]:",
    LAR_ZMIN,
    LAR_ZMAX
)

print(
    "\nUnique LAr-fiducial true muons entering TMS:",
    len(
        unique_muons
    )
)

mult = Counter(
    len(
        x[
            "slice_entries"
        ]
    )
    for x in
    unique_muons.values()
)

print(
    "\nSlice multiplicity:"
)

for n in sorted(
    mult
):
    print(
        f"  {n} slices :",
        mult[n]
    )


mu_minus_summary = (
    print_species_summary(
        13,
        "mu-  (matched true PDG +13)"
    )
)

mu_plus_summary = (
    print_species_summary(
        -13,
        "mu+  (matched true PDG -13)"
    )
)


# ============================================================
# 0.5–5 GeV PHYSICS-RANGE TABLES
# ============================================================

KE_EDGES_PHYS = np.arange(
    500.0,
    5000.0 + 500.0,
    500.0
)

P_EDGES_PHYS = np.arange(
    500.0,
    5000.0 + 500.0,
    500.0
)


for pdg, label in [
    (13, "mu-"),
    (-13, "mu+")
]:

    rows = species_rows(
        pdg
    )

    ke_table = binned_metrics(
        rows,
        KE_EDGES_PHYS,
        "ke"
    )

    p_table = binned_metrics(
        rows,
        P_EDGES_PHYS,
        "ptms"
    )

    print_table(
        f"{label} : TRUE KE 0.5–5 GeV",
        ke_table,
        "MeV"
    )

    print_table(
        f"{label} : TRUE pTMS 0.5–5 GeV/c",
        p_table,
        "MeV/c"
    )


# ============================================================
# WRITE EVENT-LEVEL CSV FOR COMPLETE AUDITABILITY
# ============================================================

csv_path = os.path.join(
    OUTDIR_V3,
    "v3_unique_muon_results.csv"
)

with open(
    csv_path,
    "w"
) as fout:

    fout.write(
        "pdg,ke_mev,ptms_mevc,"
        "n_slices,n_candidates,"
        "n_valid_candidates,"
        "reco_available,correct,"
        "sd_field_mm,"
        "chosen_entry,"
        "chosen_track,"
        "chosen_length3d,"
        "chosen_nhits\n"
    )

    for rec in (
        unique_muons.values()
    ):

        best = rec[
            "best"
        ]

        if best is None:

            values = [
                rec["pdg"],
                rec["ke"],
                rec["ptms"],
                len(
                    rec[
                        "slice_entries"
                    ]
                ),
                rec[
                    "candidate_count"
                ],
                rec[
                    "valid_candidate_count"
                ],
                0,
                0,
                "nan",
                -1,
                -1,
                "nan",
                -1,
            ]

        else:

            values = [
                rec["pdg"],
                rec["ke"],
                rec["ptms"],
                len(
                    rec[
                        "slice_entries"
                    ]
                ),
                rec[
                    "candidate_count"
                ],
                rec[
                    "valid_candidate_count"
                ],
                1,
                int(
                    rec[
                        "correct"
                    ]
                ),
                rec[
                    "sd"
                ],
                best[
                    "entry"
                ],
                best[
                    "track_index"
                ],
                best[
                    "length3d"
                ],
                best[
                    "nHits"
                ],
            ]

        fout.write(
            ",".join(
                str(v)
                for v in values
            )
            + "\n"
        )


print("\n")
print(
    "V3 CSV written to:"
)

print(
    csv_path
)

print("\nIMPORTANT:")
print(
    "Do NOT freeze V3 yet."
)

print(
    "First validate its denominators, "
    "reco availability and overlap behaviour."
)
print()
