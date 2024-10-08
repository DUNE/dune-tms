import glob
import argparse
import os
import math
import ROOT
import array
import logging

# truth_signed_distance_mu_momentum_slices.py
# This script creates histograms of the signed distance of muons and anti-muons
# using the Truth_Info TTree. The histograms sliced based in muon kinetic energy


# TODO: add the B Field info...from the inlist.

# Set ROOT to batch mode and configure styles
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

# setup the logger
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Define muon mass and fiducial cuts
MUON_MASS = 105.7  # MeV/c^2
FUDICIAL_CUT = 50
LAR_START = (-3478.48, -2166.71, 4179.24)
LAR_END = (3478.48, 829.282, 9135.88)

class Momentum:
    def __init__(self, kinetic_energy, classification="muon"):
        self.ke = kinetic_energy
        self.classification = classification
        self.ranges = [(0, 250), (250, 500), (500, 750), (750, 1000),
                       (1000, 2000), (2000, 3000), (3000, 4000),
                       (4000, 4250), (4250, 4500), (4500, 4750), (4750, 5000)]  # MeV

    def momentum_from_kinetic_energy(self, ke):
        return math.sqrt((ke + MUON_MASS) ** 2 - MUON_MASS ** 2)

    def get_muon_ke_index(self):
        for i, r in enumerate(self.ranges):
            if r[0] <= self.ke <= r[1]:
                return i
        return None


class Neutrino:
    def __init__(self, energy, classification='nu'):
        self.energy = energy
        self.classification = classification
        self.ranges = [(i, i + 500) for i in range(0, 5000, 500)]  # in 500 MeV slices

    def get_enu_index(self):
        for i, r in enumerate(self.ranges):
            if r[0] <= self.energy <= r[1]:
                return i
        return None


def get_muon_ke_entering_tms(momentum_tms_start):
    return math.sqrt(momentum_tms_start.Mag2() + MUON_MASS ** 2) - MUON_MASS

def calc_signed_distance(p_x, p_z, x_end_death, z_end_death, x_start_tms, z_start_tms):
    return x_end_death - (p_x / p_z * z_end_death + (x_start_tms - p_x / p_z * z_start_tms))

# NOTE: all values here are in mm.
def inside_tms(x, y, z):
    return -3500 < x < 3500 and -3700 < y < 1000 and 11000 < z < 18200

def inside_lar(x, y, z):
    return -4500 < x < 3700 and -3200 < y < 1000 and 4100 < z < 9200

def region1(x):
    return -4000 < x < -2500

def region2(x):
    return -1500 < x < 1500

def region3(x):
    return 2500 < x < 4000

def run(c, truth, outfilename, nmax=-1):
    bin_edges = array.array('d', [0, 250, 500, 750, 1000, 2000, 3000, 4000, 4250, 4500, 4750, 5000])
    bin_edges_gev = [edge // 1000 for edge in bin_edges]  # integer division to get GeV

    # create histogram, and bin edges are from -2m -- 2m. This is a list of TH1s
    # These two lists of TH1s will use the KE from the True_MuonKE (from its birth)
    hist_signed_distance_muon_lar_ke = [ROOT.TH1D(f"muon_lar_ke_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]
    hist_signed_distance_amuon_lar_ke = [ROOT.TH1D(f"amuon_lar_ke_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]

    # create histograms that will be sliced based on their entering TMS KE
    hist_signed_distance_muon_tms_ke = [ROOT.TH1D(f"muon_tms_ke_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]
    hist_signed_distance_amuon_tms_ke = [ROOT.TH1D(f"amuon_tms_ke_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]

    # create histograms that will be sliced based on the Neutrino Energy
    hist_signed_distance_enu_lar = [ROOT.TH1D(f"nu_lar_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges_gev)-1)]
    hist_signed_distance_enubar_lar = [ROOT.TH1D(f"nubar_lar_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges_gev)-1)]

    # Set axis labels for each histogram
    edge_counter = 0
    for hist in (hist_signed_distance_muon_lar_ke + hist_signed_distance_amuon_lar_ke +
                 hist_signed_distance_muon_tms_ke + hist_signed_distance_amuon_tms_ke +
                 hist_signed_distance_enu_lar + hist_signed_distance_enubar_lar):
        hist.SetXTitle("Signed Distance (mm)")
        hist.SetYTitle("Events")
        hist.GetYaxis().SetTitleOffset(0.95)
        if hist.GetName().startswith("muon") and "tms_ke" in hist.GetName():  # only want to set the title for the first histogram drawn
            hist.SetTitle(rf"{bin_edges[edge_counter]} < {'KE_{#mu}'} < {bin_edges[edge_counter + 1]} MeV Entering TMS")
        elif hist.GetName().startswith("amuon") and "tms_ke" in hist.GetName():
            hist.SetTitle("")  # Remove the histogram title
        elif hist.GetName().startswith("muon") and "lar_ke" in hist.GetName():
            hist.SetTitle(rf"{bin_edges[edge_counter]} < {'KE_{#mu}'} < {bin_edges[edge_counter + 1]} MeV Inside LAr")
        elif hist.GetName().startswith("amuon") and "lar_ke" in hist.GetName():
            hist.SetTitle("")  # Remove the histogram title
        elif hist.GetName().startswith("nu") and "lar" in hist.GetName():
            hist.SetTitle(rf"{bin_edges[edge_counter]} < E_{{#nu}} < {bin_edges[edge_counter + 1]} GeV Inside LAr")
        elif hist.GetName().startswith("nubar") and "lar" in hist.GetName():
            hist.SetTitle("")
        else:
            raise ValueError("Histogram naming error")
        edge_counter += 1
        if edge_counter == len(bin_edges) - 1:
            edge_counter = 0  # rest the counter

    nevents = min(c.GetEntries(), nmax if nmax >= 0 else float('inf'))

    for i in range(nevents):
        if i % 10000 == 0:
            logging.info(f"Processing event {i} / {nevents} ({i/nevents*100:.1f}%)")
        c.GetEntry(i)
        truth.GetEntry(i)

        for index, pdg in enumerate(truth.PDG):
            signed_dist = None
            if pdg in [13, -13]:
                x_start = truth.BirthPosition[4*index]
                y_start = truth.BirthPosition[4*index+1]
                z_start = truth.BirthPosition[4*index+2]
                x_end = truth.DeathPosition[4*index]
                y_end = truth.DeathPosition[4*index+1]
                z_end = truth.DeathPosition[4*index+2]
                x_start_tms = truth.PositionTMSStart[4*index]
                z_start_tms = truth.PositionTMSStart[4*index+2]

                enu = truth.NeutrinoP4[3]  # fourth component is the energy
                print('index -- pdg -- enu:', index, pdg, enu)

                if isinstance(truth.Muon_TrueKE, (list, tuple, ROOT.vector('float'))):
                    KE_muon = truth.Muon_TrueKE[index]
                else:
                    KE_muon = truth.Muon_TrueKE

                if inside_lar(x_start, y_start, z_start) and inside_tms(x_end, y_end, z_end):
                    p_tms_start = ROOT.TVector3(truth.MomentumTMSStart[4 * index], truth.MomentumTMSStart[4 * index + 1], truth.MomentumTMSStart[4 * index + 2])
                    ke_muon_tms_start = get_muon_ke_entering_tms(p_tms_start)
                    pz_tms_start, px_tms_start = truth.MomentumTMSStart[4 * index + 2], truth.MomentumTMSStart[4 * index]

                    if region1(x_start_tms) and region1(x_end):
                        if pz_tms_start != 0:
                            signed_dist = calc_signed_distance(px_tms_start, pz_tms_start, x_end, z_end, x_start_tms, z_start_tms)
                    elif region2(x_start_tms) and region2(x_end):
                        if pz_tms_start != 0:
                            signed_dist = - calc_signed_distance(px_tms_start, pz_tms_start, x_end, z_end, x_start_tms, z_start_tms)
                    elif region3(x_start_tms) and region3(x_end):
                        if pz_tms_start != 0:
                            signed_dist = calc_signed_distance(px_tms_start, pz_tms_start, x_end, z_end, x_start_tms, z_start_tms)

                    # this is the Muon_TrueKE from the birth of the muon (inside LAr)
                    p = Momentum(KE_muon, classification="muon" if pdg == 13 else "amuon")
                    range_index = p.get_muon_ke_index()  # find which set of KE ranges the muon KE falls into
                    if range_index is not None and signed_dist is not None:
                        if pdg == 13:
                            hist_signed_distance_muon_lar_ke[range_index].Fill(signed_dist)
                        else:
                            hist_signed_distance_amuon_lar_ke[range_index].Fill(signed_dist)

                    # this is the MuonKE entering the TMS
                    p_tms = Momentum(ke_muon_tms_start, classification="muon" if pdg == 13 else "amuon")
                    range_index_tms = p_tms.get_muon_ke_index()
                    if range_index_tms is not None and signed_dist is not None:
                        if pdg == 13:
                            hist_signed_distance_muon_tms_ke[range_index_tms].Fill(signed_dist)
                        else:
                            hist_signed_distance_amuon_tms_ke[range_index_tms].Fill(signed_dist)

                    # fill the neutrino histograms
                    e_nu = Neutrino(enu, classification="nu" if pdg == 14 else "nubar")
                    range_index_nu = e_nu.get_enu_index()
                    if range_index_nu is not None and signed_dist is not None:
                        if pdg == 14:
                            hist_signed_distance_enu_lar[range_index_nu].Fill(signed_dist)
                        else:
                            hist_signed_distance_enubar_lar[range_index_nu].Fill(signed_dist)

    logging.info("Done filling histograms.")


    tf = ROOT.TFile(outfilename, "recreate")
    canvas = ROOT.TCanvas("canvas", "", 800, 600)  # Set an empty string for the title to avoid unwanted titles

    # make the legend only once, since they will be the same
    legend = ROOT.TLegend(0.65, 0.75, 0.95, 0.9)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.AddEntry(hist_signed_distance_muon_lar_ke[0], "#mu^{-}", "l")  # just pick some index
    legend.AddEntry(hist_signed_distance_amuon_lar_ke[0], "#mu^{+}", "l")

    edge_counter = 0
    # loop through the mu and amu hists together for LAr KE and TMS KE.
    for index, (hists_mu_list, hists_amu_list) in enumerate([[hist_signed_distance_muon_lar_ke, hist_signed_distance_amuon_lar_ke], [hist_signed_distance_muon_tms_ke, hist_signed_distance_amuon_tms_ke], [hist_signed_distance_enu_lar, hist_signed_distance_enubar_lar]]):
        for hist_mu, hist_amu in zip(hists_mu_list, hists_amu_list):  # pair these together
            hist_mu.SetLineColor(ROOT.kRed)
            hist_amu.SetLineColor(ROOT.kBlue)

            hist_mu.SetLineWidth(3)
            hist_amu.SetLineWidth(3)

            hist_mu.GetXaxis().CenterTitle()
            hist_amu.GetXaxis().CenterTitle()
            hist_mu.GetYaxis().CenterTitle()
            hist_amu.GetYaxis().CenterTitle()

            hist_mu.Draw("hist")
            hist_amu.Draw("hist same")

            max_content = max(hist_mu.GetMaximum(), hist_amu.GetMaximum())
            hist_mu.SetMaximum(max_content * 1.25)
            hist_amu.SetMaximum(max_content * 1.25)

            legend.Draw("same")

            # vertical line for easier reading
            line0 = ROOT.TLine(0, 0, 0, max_content)
            line0.SetLineColor(ROOT.kBlack)
            line0.SetLineStyle(2)
            line0.Draw("same")

            # get the integrals, be sure to not double count the 0 mm bin
            muon_integral = hist_mu.Integral()
            events_mu_gt_0 = hist_mu.Integral(hist_mu.FindBin(0), hist_mu.GetNbinsX())
            events_mu_lt_0 = hist_mu.Integral(1, (hist_mu.FindBin(0) - 1))
            assert muon_integral == events_mu_gt_0 + events_mu_lt_0, f"Integral mu calculation failed, {hist_mu.GetTitle()}: {muon_integral} != {events_mu_gt_0} + {events_mu_lt_0}"

            amuon_integral = hist_amu.Integral()
            events_am_gt_0 = hist_amu.Integral(hist_amu.FindBin(0), hist_amu.GetNbinsX())
            events_am_lt_0 = hist_amu.Integral(1, (hist_amu.FindBin(0) - 1))
            assert amuon_integral == events_am_gt_0 + events_am_lt_0, f"Integral amuon calculation failed, {hist_amu.GetTitle()}: {amuon_integral} != {events_am_gt_0} + {events_am_lt_0}"

            # efficiency for Mu and AMu
            # purity for events of signed distance > or < 0 mm.
            # we may have zero division error if the integral is zero -- when testing over only a few events.
            try:
                efficiency_mu_gt_0 = events_mu_gt_0 / muon_integral
                efficiency_amuon_lt_0 = events_am_lt_0 / amuon_integral
                purity_mu = events_mu_gt_0 / (events_mu_gt_0 + events_am_gt_0)
                purity_amuon = events_am_lt_0 / (events_mu_lt_0 + events_am_lt_0)
            except ZeroDivisionError:
                print("Zero division error. Probably due to empty histograms. Setting to -5.0.")
                efficiency_mu_gt_0 = -5
                efficiency_amuon_lt_0 = -5
                purity_mu = -5
                purity_amuon = -5


            pt = ROOT.TPaveText(0.18, 0.7, 0.48, 0.85, "NDC")
            pt.AddText(f"Efficiency S.D. > 0 mm: {efficiency_mu_gt_0*100:.2f} %")
            pt.AddText(f"Efficiency S.D. < 0 mm: {efficiency_amuon_lt_0*100:.2f} %")
            pt.AddText(f"Purity S.D. > 0 mm: {purity_mu*100:.2f} %")
            pt.AddText(f"Purity S.D. < 0 mm: {purity_amuon*100:.2f} %")
            pt.SetTextSize(0.03)
            pt.SetTextFont(102)
            pt.SetBorderSize(0)
            pt.SetFillStyle(0)
            pt.Draw("same")

            # use the index to determine if inside LAr or TMS
            label = ''
            if index == 0:
                label = "lar" # inside LAr
            elif index == 1:
                label = "tms"
            else:
                label = "enu_lar"
            print(f'Saving plots of {label}')

            canvas.Write()
            for ext in ['png', 'pdf']:
                canvas.Print(f"{outfilename.replace('.root', '')}_{label}_ke_{bin_edges[edge_counter]}MeV_{bin_edges[edge_counter+1]}MeV." + ext)
            edge_counter += 1
            if edge_counter == len(bin_edges) - 1:
                edge_counter = 0  # reset the counter
            # end of loop through mu and amu hists

    for hist in (hist_signed_distance_muon_lar_ke + hist_signed_distance_amuon_lar_ke +
                 hist_signed_distance_muon_tms_ke + hist_signed_distance_amuon_tms_ke +
                 hist_signed_distance_enu_lar + hist_signed_distance_enubar_lar):
        hist.Write()
    tf.Close()
    logging.info("Done saving.")

def validate_then_run(args):
    indir = args.indir
    inlist = args.inlist
    infile = args.filename
    files_to_use = [infile] if infile else glob.glob(f"{indir}/**/*tmsreco.root", recursive=True) if indir else open(inlist).read().splitlines()
    if not files_to_use: 
        raise ValueError("No input files found")

    outdir = args.outdir or f"/exp/dune/data/users/{os.environ['USER']}/dune-tms_hists/truth_signed_distance_mu_momentum_slices"
    os.makedirs(outdir, exist_ok=True)
    outfilename = os.path.join(outdir, args.out_rootfile_name)
    if os.path.exists(outfilename) and not args.allow_overwrite:
        raise ValueError(f"Output file {outfilename} already exists")

    # get the TTrees: Line_Candidates (Reco) and Truth_Info
    c, truth = ROOT.TChain("Line_Candidates"), ROOT.TChain("Truth_Info")
    for f in files_to_use:
        c.Add(f)
        truth.Add(f)
    assert c.GetEntries() > 0 and truth.GetEntries() > 0, "No entries in input files"

    run(c, truth, outfilename, args.nevents)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', type=str, default="")
    parser.add_argument('--out_rootfile_name', type=str, default="dune-tms_sd_mu_momentum_slices.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', type=str, default="")
    group.add_argument('--inlist', type=str, default="")
    group.add_argument('--filename', '-f', type=str, default="")
    parser.add_argument('--nevents', '-n', type=int, default=-1)
    parser.add_argument('--allow_overwrite', action='store_true')
    parser.add_argument('--preview', action='store_true')

    args = parser.parse_args()
    validate_then_run(args)

