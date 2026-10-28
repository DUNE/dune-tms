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
# STORAGE
# ============================================================

sd_minus = []
sd_plus = []

raw_minus = []
raw_plus = []

# V1 values for the exact same tracks that survive V2.
sd_minus_v1_same_tracks = []
sd_plus_v1_same_tracks = []

n_mu = 0
n_lar = 0
n_touch = 0
n_contained = 0
n_track_mismatch = 0
n_valid_reco = 0
n_nonzero_field = 0


# ============================================================
# EVENT LOOP
# ============================================================

for ievt in range(r.GetEntries()):

    r.GetEntry(ievt)
    t.GetEntry(ievt)

    nreco = int(r.nTracks)
    ntruth = int(t.RecoTrackN)

    if nreco <= 0:
        continue

    # Truth_Info is organized per reconstructed track.
    # Require exact track multiplicity correspondence.
    if nreco != ntruth:
        n_track_mismatch += 1
        continue


    start = np.frombuffer(
        r.StartPos,
        dtype=np.float32,
        count=nreco*4
    ).reshape(nreco,4)

    end = np.frombuffer(
        r.EndPos,
        dtype=np.float32,
        count=nreco*4
    ).reshape(nreco,4)

    direction = np.frombuffer(
        r.StartDirection,
        dtype=np.float32,
        count=nreco*3
    ).reshape(nreco,3)

    nhits = np.frombuffer(
        r.nHits,
        dtype=np.int32,
        count=nreco
    )

    hits_all = np.frombuffer(
        r.TrackHitPos,
        dtype=np.float32,
        count=nreco*200*4
    ).reshape(nreco,200,4)


    for j in range(nreco):

        # ====================================================
        # TRUTH USED ONLY FOR SELECTION / LABEL
        # ====================================================

        pdg = int(
            t.RecoTrackPrimaryParticlePDG[j]
        )

        if abs(pdg) != 13:
            continue

        n_mu += 1


        if not bool(
            t.RecoTrackPrimaryParticleLArFiducialStart[j]
        ):
            continue

        n_lar += 1


        if not bool(
            t.RecoTrackPrimaryParticleTMSFiducialTouch[j]
        ):
            continue

        n_touch += 1


        if not bool(
            t.RecoTrackPrimaryParticleTMSFiducialEnd[j]
        ):
            continue

        n_contained += 1


        # ====================================================
        # RECONSTRUCTED TRACK QUANTITIES
        # ====================================================

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
            np.all(np.isfinite(p0)) and
            np.all(np.isfinite(p3)) and
            np.all(np.isfinite(d0))
        ):
            continue


        x1 = p0[0]
        z1 = p0[2]

        x3 = p3[0]
        z3 = p3[2]

        dxdir = d0[0]
        dzdir = d0[2]


        # Reco tracks observed above are oriented
        # upstream StartPos -> downstream EndPos.
        if z3 <= z1:
            continue

        if dzdir <= 1e-6:
            continue


        # ====================================================
        # V1 RAW SD — StartDirection
        #
        # Kept only so V1 and V2 can be compared on exactly
        # the same tracks.
        # ====================================================

        slope_v1 = dxdir / dzdir

        x_extrap_v1 = (
            x1 +
            slope_v1 * (z3-z1)
        )

        sd_raw_v1 = (
            x3 -
            x_extrap_v1
        )


        # ====================================================
        # RECONSTRUCTED HIT PATH
        # ====================================================

        nh = int(nhits[j])

        if nh < 2:
            continue

        nh = min(nh,200)

        hits = np.asarray(
            hits_all[j,:nh,:3],
            dtype=np.float64
        )


        good = (
            np.all(
                np.isfinite(hits),
                axis=1
            )
            &
            np.all(
                np.abs(hits) < 1e7,
                axis=1
            )
        )

        hits = hits[good]

        if len(hits) < 2:
            continue


        # We established experimentally that TrackHitPos
        # is stored downstream -> upstream.
        #
        # Still orient robustly by choosing the end
        # closest to reconstructed StartPos.

        d_first = np.linalg.norm(
            hits[0] - p0
        )

        d_last = np.linalg.norm(
            hits[-1] - p0
        )

        if d_last < d_first:
            hits = hits[::-1]


        # ====================================================
        # V2 LOCAL ENTRANCE-DIRECTION FIT
        #
        # Fit x(z) using only the first upstream reconstructed
        # hits. No truth position, momentum, or PDG enters
        # this fit.
        # ====================================================

        nfit = min(
            LOCAL_FIT_NHITS,
            len(hits)
        )

        if nfit < 3:
            continue

        fit_hits = np.asarray(
            hits[:nfit],
            dtype=np.float64
        )

        zfit = fit_hits[:,2]
        xfit = fit_hits[:,0]

        # Protect against pathological / duplicate-z fits.
        if np.ptp(zfit) < 100.0:
            continue

        slope_v2, intercept_v2 = np.polyfit(
            zfit,
            xfit,
            1
        )

        if not (
            math.isfinite(slope_v2) and
            math.isfinite(intercept_v2)
        ):
            continue

        # Keep the same reconstructed StartPos used in V1.
        # Only the entrance slope is changed.
        x_extrap_v2 = (
            x1 +
            slope_v2 * (z3-z1)
        )

        sd_raw = (
            x3 -
            x_extrap_v2
        )


        # Complete reconstructed path:
        # StartPos -> TrackHitPos -> EndPos

        path = np.vstack([
            p0,
            hits,
            p3
        ])


        # Remove duplicate adjacent points

        clean = [path[0]]

        for point in path[1:]:

            if np.linalg.norm(
                point-clean[-1]
            ) > 1e-3:

                clean.append(point)

        path = np.asarray(clean)

        if len(path) < 2:
            continue

        n_valid_reco += 1


        # ====================================================
        # FIELD-AWARE RECO BENDING INTEGRAL
        #
        # NO truth PDG is used in Kx.
        #
        # Kx =
        # integral
        # (zend-z) (uz By - uy Bz) ds
        # ====================================================

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

            if seglen < 1e-6:
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


            for k in range(nsub):

                frac = (
                    k+0.5
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
            continue

        if abs(Kx) < 1e-9:
            continue

        n_nonzero_field += 1


        # ====================================================
        # FIELD-AWARE RECONSTRUCTED SD
        # ====================================================

        # V2 field-aware SD
        sd_field = (
            np.sign(Kx) *
            sd_raw
        )

        # V1 field-aware SD for the SAME accepted track
        sd_field_v1 = (
            np.sign(Kx) *
            sd_raw_v1
        )


        if pdg == 13:

            raw_minus.append(
                sd_raw
            )

            sd_minus.append(
                sd_field
            )

            sd_minus_v1_same_tracks.append(
                sd_field_v1
            )


        elif pdg == -13:

            raw_plus.append(
                sd_raw
            )

            sd_plus.append(
                sd_field
            )

            sd_plus_v1_same_tracks.append(
                sd_field_v1
            )


# ============================================================
# PRINT RESULTS
# ============================================================

def summary(
    name,
    raw,
    corrected,
    expected_positive
):

    raw = np.asarray(
        raw,
        dtype=float
    )

    corrected = np.asarray(
        corrected,
        dtype=float
    )

    print("\n" + "="*65)
    print(name)
    print("="*65)

    print("N =", len(corrected))

    if len(corrected) == 0:
        return


    print("\nRAW RECO SD")

    print(
        "mean   =",
        np.mean(raw),
        "mm"
    )

    print(
        "median =",
        np.median(raw),
        "mm"
    )


    print("\nFIELD-AWARE RECO SD")

    print(
        "mean   =",
        np.mean(corrected),
        "mm"
    )

    print(
        "median =",
        np.median(corrected),
        "mm"
    )

    print(
        "std    =",
        np.std(corrected),
        "mm"
    )

    print(
        "q16/q50/q84 =",
        np.percentile(
            corrected,
            [16,50,84]
        )
    )


    npos = int(
        np.sum(corrected > 0)
    )

    nneg = int(
        np.sum(corrected < 0)
    )


    print(
        "SD > 0:",
        npos,
        "/",
        len(corrected),
        "=",
        npos/len(corrected)
    )

    print(
        "SD < 0:",
        nneg,
        "/",
        len(corrected),
        "=",
        nneg/len(corrected)
    )


    correct = (
        npos
        if expected_positive
        else nneg
    )


    print(
        "expected charge side:",
        correct,
        "/",
        len(corrected),
        "=",
        correct/len(corrected)
    )


print("\n")
print("="*70)
print("FIELD-AWARE RECONSTRUCTED SIGNED DISTANCE — V2 LOCAL ENTRANCE FIT")
print("MERGED MAPPED RUNS 1001-1005")
print("="*70)

print(
    "truth-matched muons          :",
    n_mu
)

print(
    "truth starts in LAr          :",
    n_lar
)

print(
    "truth enters TMS             :",
    n_touch
)

print(
    "truth endpoint inside TMS    :",
    n_contained
)

print(
    "track-count mismatched events:",
    n_track_mismatch
)

print(
    "valid reconstructed paths    :",
    n_valid_reco
)

print(
    "non-zero field integral      :",
    n_nonzero_field
)


summary(
    "mu- RECO SD",
    raw_minus,
    sd_minus,
    True
)

summary(
    "mu+ RECO SD",
    raw_plus,
    sd_plus,
    False
)



# ============================================================
# EXACT SAME-TRACK V1 vs V2 COMPARISON
# ============================================================

def expected_side_fraction(values, positive_side):

    a = np.asarray(values, dtype=float)

    if len(a) == 0:
        return 0, 0, float("nan")

    if positive_side:
        good = int(np.sum(a > 0))
    else:
        good = int(np.sum(a < 0))

    return good, len(a), good/len(a)


print("\n")
print("="*70)
print("EXACT SAME-TRACK COMPARISON: V1 vs V2")
print("="*70)

g1, n1, f1 = expected_side_fraction(
    sd_minus_v1_same_tracks,
    True
)

g2, n2, f2 = expected_side_fraction(
    sd_minus,
    True
)

print("\nmu-:")
print(
    "V1 StartDirection :",
    g1, "/", n1, "=",
    f1
)

print(
    "V2 local x(z) fit :",
    g2, "/", n2, "=",
    f2
)


g1p, n1p, f1p = expected_side_fraction(
    sd_plus_v1_same_tracks,
    False
)

g2p, n2p, f2p = expected_side_fraction(
    sd_plus,
    False
)

print("\nmu+:")
print(
    "V1 StartDirection :",
    g1p, "/", n1p, "=",
    f1p
)

print(
    "V2 local x(z) fit :",
    g2p, "/", n2p, "=",
    f2p
)

print(
    "\nLocal-fit hits used:",
    LOCAL_FIT_NHITS
)


# ============================================================
# PLOT
# ============================================================

NBINS = 80
XMIN = -2000.0
XMAX = 2000.0

hminus = ROOT.TH1D(
    "hminus",
    "",
    NBINS,
    XMIN,
    XMAX
)

hplus = ROOT.TH1D(
    "hplus",
    "",
    NBINS,
    XMIN,
    XMAX
)


for value in sd_minus:
    hminus.Fill(value)

for value in sd_plus:
    hplus.Fill(value)


hminus.SetLineColor(
    ROOT.kBlue+1
)

hminus.SetLineWidth(2)

hplus.SetLineColor(
    ROOT.kRed+1
)

hplus.SetLineWidth(2)


c = ROOT.TCanvas(
    "c",
    "",
    1050,
    750
)

c.SetTopMargin(0.10)
c.SetBottomMargin(0.13)
c.SetLeftMargin(0.11)
c.SetRightMargin(0.05)


hminus.GetXaxis().SetTitle(
    "Field-aware reconstructed signed distance [mm]"
)

hminus.GetYaxis().SetTitle(
    "Muon tracks"
)


maximum = max(
    hminus.GetMaximum(),
    hplus.GetMaximum()
)

if maximum <= 0:
    maximum = 1

hminus.SetMaximum(
    1.25*maximum
)

hminus.SetMinimum(0)

hminus.Draw("HIST")
hplus.Draw("HIST SAME")


zero = ROOT.TLine(
    0,
    0,
    0,
    hminus.GetMaximum()
)

zero.SetLineStyle(2)
zero.Draw()


title = ROOT.TLatex()

title.SetNDC(True)
title.SetTextAlign(22)
title.SetTextSize(0.045)

title.DrawLatex(
    0.50,
    0.955,
    "Field-aware reconstructed signed distance - V2 local fit"
)


leg = ROOT.TLegend(
    0.72,
    0.80,
    0.92,
    0.90
)

leg.SetBorderSize(0)
leg.SetFillStyle(0)

leg.AddEntry(
    hminus,
    "#mu^{-}",
    "l"
)

leg.AddEntry(
    hplus,
    "#mu^{+}",
    "l"
)

leg.Draw()


png = os.path.join(
    OUT,
    "reco_fieldaware_signed_distance_V2_localfit.png"
)

c.SaveAs(png)

print("\nPlot:")
print(png)
