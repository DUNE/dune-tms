import glob
import argparse
import os
import math
import ROOT
import array
import logging

# run this script:
#   python signed_distance_efficiency_vs_muon_true_ke.py
# This script creates histograms of the efficiency of signed distance of muons and anti-muons
# using the Truth_Info TTree, and saves them to a ROOT file.
# S.D. = signed distance


# TODO: add the B Field info...from the inlist to outfile name.

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
MUON_KE_BINNING = [(i, i+100) for i in range(0, 5001, 100)]  # MeV

class Momentum:
    def __init__(self, kinetic_energy, classification="muon"):
        self.ke = kinetic_energy
        self.classification = classification
        self.ranges = MUON_KE_BINNING  # MeV

    def momentum_from_kinetic_energy(self, ke):
        return math.sqrt((ke + MUON_MASS) ** 2 - MUON_MASS ** 2)

    def get_muon_ke_bin(self):
        for i, r in enumerate(self.ranges):
            if r[0] <= self.ke < r[1]:  # [lower, upper) bounds
                return i + 1  # bin number starts at 1 (enumerate idx starts at 0)
            elif self.ke >= self.ranges[-1][1]:  # if the KE is greater than the last bin, return the last bin, which is an overflow.
                return len(self.ranges)
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

def run(truth, outfilename, nmax=-1):
    bin_edges = array.array('d', [i for i in range(0, 5001, 100)])  # 1 bin is 100 MeV
    # bin_edges_gev = [edge // 1000 for edge in bin_edges]  # integer division to get GeV

    # create histogram, for 'True_MuonKE'.
    # These two lists of TH1s will use the KE from the True_MuonKE (from its birth)
    hist_sd_eff_muon_lar_ke = ROOT.TH1D("muon_sd_eff_lar_ke", "", len(bin_edges)-1, bin_edges[0], bin_edges[-1])  # 0-5 GeV.
    hist_sd_eff_amuon_lar_ke = ROOT.TH1D("amuon_sd_eff_lar_ke", "", len(bin_edges)-1, bin_edges[0], bin_edges[-1])

    # create histogram, for muon TMS KE.
    hist_sd_eff_muon_tms_ke = ROOT.TH1D("muon_sd_eff_tms_ke", "", len(bin_edges)-1, bin_edges[0], bin_edges[-1])
    hist_sd_eff_amuon_tms_ke = ROOT.TH1D("amuon_sd_eff_tms_ke", "", len(bin_edges)-1, bin_edges[0], bin_edges[-1])


    # Set axis labels for each histogram
    edge_counter = 0
    for hist in (hist_sd_eff_muon_lar_ke , hist_sd_eff_amuon_lar_ke , hist_sd_eff_muon_tms_ke , hist_sd_eff_amuon_tms_ke):
        hist.SetXTitle("Muon Kinetic Energy (MeV)")
        hist.SetYTitle("Efficiency")
        hist.GetYaxis().SetTitleOffset(0.95)
    hist_sd_eff_muon_lar_ke.SetTitle("Signed Distance Efficiency, KE Inside LAr")
    hist_sd_eff_amuon_lar_ke.SetTitle("Signed Distance Efficiency, KE Inside LAr")
    hist_sd_eff_muon_tms_ke.SetTitle("Signed Distance Efficiency, KE Entering TMS")
    hist_sd_eff_amuon_tms_ke.SetTitle("Signed Distance Efficiency, KE Entering Inside TMS")

    # dictionary of bin number (for KE) and events: {key=bin_number, value=[sd > 0, sd < 0, sd = 0]}
    # Need separate containers for muons and anti-muons.
    data_lar_mu, data_tms_mu, data_lar_amu, data_tms_amu = {}, {}, {}, {}
    for binIdx in range(len(MUON_KE_BINNING)):
        data_lar_mu[binIdx + 1] = [0, 0, 0]  # first index is S.D. > 0, second is S.D. < 0, third is S.D. = 0.
        data_tms_mu[binIdx + 1] = [0, 0, 0]
        data_lar_amu[binIdx + 1] = [0, 0, 0]
        data_tms_amu[binIdx + 1] = [0, 0, 0]

    # add up the garbage KE values from 'tmsreco.root' Truth_Info TTree
    garbage_ke_count = 0

    n_events = min(truth.GetEntries(), nmax if nmax >= 0 else float('inf'))
    for i in range(n_events):
        if i % 10000 == 0:
            logging.info(f"Processing event {i} / {n_events} ({i/n_events*100:.1f}%)")
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

                # check if muon starts in LAr and ends in TMS, before creating any vars.
                if inside_lar(x_start, y_start, z_start) and inside_tms(x_end, y_end, z_end):
                    # muon KE at birth (in LAr)
                    if isinstance(truth.Muon_TrueKE, (list, tuple, ROOT.vector('float'))):
                        ke_muon_lar = truth.Muon_TrueKE[index]
                    else:
                        ke_muon_lar = truth.Muon_TrueKE
                    muon_ke_lar_bin = Momentum(ke_muon_lar, classification="muon" if pdg == 13 else "amuon").get_muon_ke_bin()

                    # muon KE at TMS start
                    p_tms_start = ROOT.TVector3(truth.MomentumTMSStart[4 * index], truth.MomentumTMSStart[4 * index + 1], truth.MomentumTMSStart[4 * index + 2])
                    ke_muon_tms_start = get_muon_ke_entering_tms(p_tms_start)
                    muon_ke_tms_bin = Momentum(ke_muon_tms_start, classification="muon" if pdg == 13 else "amuon").get_muon_ke_bin()

                    # NOTE: bin number may be None if the KE is outside the binning range.
                    # this is bc some of the Truth_Info have garbage values (173204975.05692) for KE.
                    if muon_ke_lar_bin is None or muon_ke_tms_bin is None:
                        garbage_ke_count += 1
                        logger.warning(f"Garbage KE value found: {ke_muon_lar} or {ke_muon_tms_start}")
                        continue
                    # starting and ending momenta in TMS
                    pz_tms_start, px_tms_start = truth.MomentumTMSStart[4 * index + 2], truth.MomentumTMSStart[4 * index]

                    # must be in the same region to be considered. A good # todo for future...multiple regions.
                    if (region1(x_start) and region1(x_end)) or (region2(x_start) and region2(x_end)) or (region3(x_start) and region3(x_end)):
                        if pz_tms_start != 0:
                            signed_dist = calc_signed_distance(px_tms_start, pz_tms_start, x_end, z_end, x_start_tms, z_start_tms)
                            # Recall: [sd > 0, sd < 0, sd = 0]
                            if signed_dist > 0:
                                if pdg == 13:
                                    data_lar_mu[muon_ke_lar_bin][0] += 1
                                    data_tms_mu[muon_ke_tms_bin][0] += 1
                                else:
                                    data_lar_amu[muon_ke_lar_bin][0] += 1
                                    data_tms_amu[muon_ke_tms_bin][0] += 1
                            elif signed_dist < 0:
                                if pdg == 13:
                                    data_lar_mu[muon_ke_lar_bin][1] += 1
                                    data_tms_mu[muon_ke_tms_bin][1] += 1
                                else:
                                    data_lar_amu[muon_ke_lar_bin][1] += 1
                                    data_tms_amu[muon_ke_tms_bin][1] += 1
                            elif signed_dist == 0:
                                if pdg == 13:
                                    data_lar_mu[muon_ke_lar_bin][2] += 1
                                    data_tms_mu[muon_ke_tms_bin][2] += 1
                                else:
                                    data_lar_amu[muon_ke_lar_bin][2] += 1
                                    data_tms_amu[muon_ke_tms_bin][2] += 1
                            else:
                                print('Unknown sign distance', signed_dist)
                    else:
                        continue
    logging.info("Done looping through events.")
    logging.warning(f"Garbage KE count: {garbage_ke_count}")

    # fill the histograms via SetBinContent()  # key is the bin number, value is the sign distance counts [SD>0, SD<0, SD=0].
    # NOTE: if --n is small, you may get a divide by zero error: just set denominator to 'inf' in that case.
    # NOTE: number of muon events is not necessarily the same as anti-muon events.
    for bin_num, sign_dists in data_lar_mu.items():
        total_countable_events_in_bin = (sign_dists[0] + sign_dists[1] + sign_dists[2])
        hist_sd_eff_muon_lar_ke.SetBinContent(bin_num, sign_dists[0] / total_countable_events_in_bin if total_countable_events_in_bin != 0 else float('inf'))
    for bin_num, sign_dists in data_lar_amu.items():
        total_countable_events_in_bin = (sign_dists[0] + sign_dists[1] + sign_dists[2])
        hist_sd_eff_amuon_lar_ke.SetBinContent(bin_num, sign_dists[1] / total_countable_events_in_bin if total_countable_events_in_bin != 0 else float('inf'))
    for bin_num, sign_dists in data_tms_mu.items():
        total_countable_events_in_bin = (sign_dists[0] + sign_dists[1] + sign_dists[2])
        hist_sd_eff_muon_tms_ke.SetBinContent(bin_num, sign_dists[0] / total_countable_events_in_bin if total_countable_events_in_bin != 0 else float('inf'))
    for bin_num, sign_dists in data_tms_amu.items():
        total_countable_events_in_bin = (sign_dists[0] + sign_dists[1] + sign_dists[2])
        hist_sd_eff_amuon_tms_ke.SetBinContent(bin_num, sign_dists[1] / total_countable_events_in_bin if total_countable_events_in_bin != 0 else float('inf'))


    tf = ROOT.TFile(outfilename, "recreate")
    canvas = ROOT.TCanvas("canvas", "", 800, 600)  # Set an empty string for the title to avoid unwanted titles

    # make the legend only once, since they will be the same
    legend = ROOT.TLegend(0.65, 0.75, 0.95, 0.9)
    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetFillColor(0)
    legend.AddEntry(hist_sd_eff_muon_lar_ke, "#mu^{-}", "l")  # just pick one mu-amu pair
    legend.AddEntry(hist_sd_eff_amuon_lar_ke, "#mu^{+}", "l")

    edge_counter = 0
    # loop through the mu and amu hists together for LAr KE and TMS KE.
    for index, (hist_mu, hist_amu) in enumerate([[hist_sd_eff_muon_lar_ke, hist_sd_eff_amuon_lar_ke], [hist_sd_eff_muon_tms_ke, hist_sd_eff_amuon_tms_ke]]):
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
        hist_mu.SetMaximum(max_content * 1.4)
        hist_amu.SetMaximum(max_content * 1.4)

        legend.Draw("same")

        # vertical line for easier reading
        line_half = ROOT.TLine(0, 0.5, hist_mu.GetBinLowEdge(hist_mu.GetNbinsX()) + hist_mu.GetBinWidth(hist_mu.GetNbinsX()), 0.5)
        line_half.SetLineColor(ROOT.kBlack)
        line_half.SetLineStyle(2)
        line_half.Draw("same")

        # use the index to determine if inside LAr or TMS
        label = ''
        if index == 0:
            label = "lar"  # inside LAr
        elif index == 1:
            label = "tms"
        else:
            print("unknown")
        print(f'Saving plots of {label}')

        canvas.Write()
        for ext in ['png', 'pdf']:
            canvas.Print(f"{outfilename.replace('.root', '')}_{label}_ke." + ext)

    for hist in (hist_sd_eff_muon_lar_ke, hist_sd_eff_amuon_lar_ke, hist_sd_eff_muon_tms_ke, hist_sd_eff_amuon_tms_ke):
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

    outdir = args.outdir or f"/exp/dune/data/users/{os.environ['USER']}/dune-tms_hists/sign_distance_efficiency_vs_muon_true_ke"
    os.makedirs(outdir, exist_ok=True)
    outfilename = os.path.join(outdir, args.out_rootfile_name)
    if os.path.exists(outfilename) and not args.allow_overwrite:
        raise ValueError(f"Output file {outfilename} already exists")

    # get the TTrees: Line_Candidates (Reco) and Truth_Info, though we don't use Line_Candidates here.
    c, truth = ROOT.TChain("Line_Candidates"), ROOT.TChain("Truth_Info")
    for f in files_to_use:
        c.Add(f)
        truth.Add(f)
    assert c.GetEntries() > 0 and truth.GetEntries() > 0, "No entries in input files"

    run(truth, outfilename, args.nevents)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', type=str, default="")
    parser.add_argument('--out_rootfile_name', type=str, default="dune-tms_sd_eff_vs_mu_true_ke.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', type=str, default="")
    group.add_argument('--inlist', type=str, default="")
    group.add_argument('--filename', '-f', type=str, default="")
    parser.add_argument('--nevents', '-n', type=int, default=-1)
    parser.add_argument('--allow_overwrite', action='store_true')
    parser.add_argument('--preview', action='store_true')

    args = parser.parse_args()
    validate_then_run(args)
