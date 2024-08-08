import glob
import argparse
import os
import math
import ROOT
import array
import logging

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

    def get_range_index(self):
        for i, r in enumerate(self.ranges):
            if r[0] <= self.ke <= r[1]:
                return i
        return None

def calc_signed_distance(p_x, p_z, x_end_death, z_end_death, x_start_tms, z_start_tms):
    return x_end_death - (p_x / p_z * z_end_death + (x_start_tms - p_x / p_z * z_start_tms))

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
    hist_signed_distance_muon = [ROOT.TH1D(f"muon_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]
    hist_signed_distance_amuon = [ROOT.TH1D(f"amuon_{i}", "", 100, -2000, 2000) for i in range(len(bin_edges)-1)]

    # Set axis labels for each histogram
    edge_counter = 0
    for hist in hist_signed_distance_muon + hist_signed_distance_amuon:
        hist.SetXTitle("Signed Distance (mm)")
        hist.SetYTitle("Events")
        hist.GetYaxis().SetTitleOffset(0.95)
        if hist.GetName().startswith("muon"):  # only want to set the title for the first histogram drawn
            hist.SetTitle(rf"{bin_edges[edge_counter]} < {'KE_{#mu}'} < {bin_edges[edge_counter + 1]} MeV")
        elif hist.GetName().startswith("amuon"):
            hist.SetTitle("")  # Remove the histogram title
        edge_counter += 1

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

                if isinstance(truth.Muon_TrueKE, (list, tuple, ROOT.vector('float'))):
                    KE_muon = truth.Muon_TrueKE[index]
                else:
                    KE_muon = truth.Muon_TrueKE

                if inside_lar(x_start, y_start, z_start) and inside_tms(x_end, y_end, z_end):
                    p_z, p_x = truth.MomentumTMSStart[4 * index + 2], truth.MomentumTMSStart[4 * index]
                    if region1(x_start_tms) and region1(x_end):
                        if p_z != 0:
                            signed_dist = calc_signed_distance(p_x, p_z, x_end, z_end, x_start_tms, z_start_tms)
                    elif region2(x_start_tms) and region2(x_end):
                        if p_z != 0: 
                            signed_dist = - calc_signed_distance(p_x, p_z, x_end, z_end, x_start_tms, z_start_tms)
                    elif region3(x_start_tms) and region3(x_end):
                        if p_z != 0: 
                            signed_dist = calc_signed_distance(p_x, p_z, x_end, z_end, x_start_tms, z_start_tms)

                    p = Momentum(KE_muon, classification="muon" if pdg == 13 else "amuon")
                    range_index = p.get_range_index()  # find which set of KE ranges the muon KE falls into
                    if range_index is not None and signed_dist is not None:
                        if pdg == 13:
                            hist_signed_distance_muon[range_index].Fill(signed_dist)
                        else:
                            hist_signed_distance_amuon[range_index].Fill(signed_dist)

    tf = ROOT.TFile(outfilename, "recreate")
    canvas = ROOT.TCanvas("canvas", "", 800, 600)  # Set an empty string for the title to avoid unwanted titles

    for i in range(len(bin_edges) - 1):
        hist_signed_distance_muon[i].SetLineColor(ROOT.kRed)
        hist_signed_distance_amuon[i].SetLineColor(ROOT.kBlue)

        hist_signed_distance_muon[i].SetLineWidth(3)
        hist_signed_distance_amuon[i].SetLineWidth(3)

        hist_signed_distance_muon[i].GetXaxis().CenterTitle()
        hist_signed_distance_amuon[i].GetXaxis().CenterTitle()
        hist_signed_distance_muon[i].GetYaxis().CenterTitle()
        hist_signed_distance_amuon[i].GetYaxis().CenterTitle()

        hist_signed_distance_muon[i].Draw("hist")
        hist_signed_distance_amuon[i].Draw("hist same")

        legend = ROOT.TLegend(0.65, 0.75, 0.95, 0.9)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetFillColor(0)
        legend.AddEntry(hist_signed_distance_muon[i], "#mu^{-}", "l")
        legend.AddEntry(hist_signed_distance_amuon[i], "#mu^{+}", "l")
        legend.Draw("same")

        # vertical line for easier reading
        max_content = max(hist_signed_distance_muon[i].GetMaximum(), hist_signed_distance_amuon[i].GetMaximum())
        line0 = ROOT.TLine(0, 0, 0, max_content)
        line0.SetLineColor(ROOT.kBlack)
        line0.SetLineStyle(2)
        line0.Draw("same")

        canvas.Write()
        for ext in ['png', 'pdf']:
            canvas.Print(f"{outfilename.replace('.root', '')}_ke_{bin_edges[i]}MeV_{bin_edges[i+1]}MeV." + ext)

    for hist in hist_signed_distance_muon + hist_signed_distance_amuon:
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

    outdir = args.outdir or f"/exp/dune/data/users/{os.environ['USER']}/dune-tms_hists/signed_distance_momentum_slices"
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
    parser.add_argument('--out_rootfile_name', type=str, default="dune-tms_hists.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', type=str, default="")
    group.add_argument('--inlist', type=str, default="")
    group.add_argument('--filename', '-f', type=str, default="")
    parser.add_argument('--nevents', '-n', type=int, default=-1)
    parser.add_argument('--allow_overwrite', action='store_true')
    parser.add_argument('--preview', action='store_true')

    args = parser.parse_args()
    validate_then_run(args)

