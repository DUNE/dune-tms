import sys
import os
import array
import random
random.seed(1234)
import math
import ROOT
ROOT.gROOT.SetBatch(True)
import dunestyle.root as dunestyle
ROOT.gStyle.SetOptStat(1111)

canvas = ROOT.TCanvas()
canvas.UseCurrentStyle()

def make_box(x1, y1, x2, y2, text_size=0.04, text_align=12):
    text_box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    text_box.SetFillColor(0)  # Transparent background
    text_box.SetFillStyle(0)  # Transparent background
    text_box.SetBorderSize(0)
    text_box.SetTextSize(text_size)
    text_box.SetTextAlign(text_align)
    return text_box
    

def add_text(x1, y1, x2, y2, lines, text_size=0.04):
    """
    Add multi-line text to a ROOT canvas.

    Parameters:
        x1, y1 (float): Bottom-left corner of the text box (NDC coordinates, 0.0 to 1.0).
        x2, y2 (float): Top-right corner of the text box (NDC coordinates, 0.0 to 1.0).
        lines (list of str): List of lines of text to display.
        text_size (float): Size of the text.
        text_align (int): Text alignment (12 = left-center, 22 = center-center, etc.).
    """
    text_box_a = make_box(x1, y1, x2, y2, text_align = 11)
    text_box_b = make_box(x1, y1, x2, y2, text_align = 21)
    text_box_c = make_box(x1, y1, x2, y2, text_align = 31)
    
    for line in lines.split("\n"):
        foo = line.split("\t")
        a = foo[0]
        b = None
        c = None
        if len(foo) > 1: b = foo[1]
        if len(foo) > 2: c = foo[2]
        if c == None and b != None: 
            c = b
            b = None
        text_box_a.AddText(a)
        if b != None: text_box_b.AddText(b)
        if c != None: text_box_c.AddText(c)
    
    text_box_a.Draw("same")
    text_box_b.Draw("same")
    text_box_c.Draw("same")
    return [text_box_a, text_box_b, text_box_c] 

filename = sys.argv[1]
assert os.path.exists(filename)
tf = ROOT.TFile(filename)
if not os.path.exists("output"):
    os.makedirs("output")
n_start = 1
n_end = 19
bin_edges = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0]
centers = [0.5*(bin_edges[i] + bin_edges[i+1]) for i in range(len(bin_edges) - 1)]
means = []
errors = []
gauss = ROOT.TF1("gauss", "gaus", -0.4, 0.4)
gauss.SetLineColor(ROOT.kRed)
crystalball = ROOT.TF1("crystalball", "crystalball", -0.4, 0.4)


def do_multi_fit(h, seeds):
    print(f"Got {len(seeds)} seeds")
    out_chi2 = 1e12
    out_result = -1
    out_ndf = 0
    out_params = [0] * 5
    out_err = [1e12] * 5
    for seed in seeds:
        #print(f"Testing {seed}")
        crystalball.SetParameters(*seed)
        fit_result = h.Fit(crystalball, "RIESQ", "")
        chi2 = crystalball.GetChisquare()
        #print(f"Got chi2 {chi2}. Need lower than {out_chi2}", int(fit_result), fit_result.IsValid(), fit_result.Prob(), fit_result.Status(), fit_result.Edm())
        if fit_result.IsValid() and chi2 < out_chi2:
            out_chi2 = chi2
            out_result = fit_result
            out_ndf = crystalball.GetNDF()
            for i in range(len(out_params)):
                out_params[i] = crystalball.GetParameter(i)
                out_err[i] = crystalball.GetParError(i)
    for i in range(len(out_params)):
        crystalball.SetParameter(i, out_params[i])
        crystalball.SetChisquare(out_chi2)
        crystalball.SetNDF(out_ndf)
    return out_result
    
def create_n_seeds(n, *seed_arrays):
    out = set()
    n_iterations = 0
    while len(out) < n and n_iterations < n * 100:
        foo = []
        for seed_array in seed_arrays:
            seed = random.choice(seed_array)
            foo.append(seed)
        foo = tuple(foo)
        out.add(foo)
        n_iterations += 1
    return list(out)

    
def do_fit(h, outfilename, i):
    #h.SetStats(True)
    mean = h.GetMean()
    stddev = h.GetStdDev()
    simple_mean = mean
    simple_stddev = stddev
    n_entries = h.Integral()
    norm = h.Integral()
    #h.Sumw2()
    #h.Scale(1/h.Integral())
    max_with_headroom = 1.2 * h.GetMaximum()
    h.GetYaxis().SetRangeUser(0, max_with_headroom)
    print(f"{i}\t{mean:0.4f}\t{stddev:0.4f}\t{norm:0.4f}")
    
    #crystalball.SetParLimits(0, 0, 30) # N
    crystalball.SetParLimits(1, -0.4, 0.4) # mean
    crystalball.SetParLimits(2, 0.0, 2) # stddev
    crystalball.SetParLimits(3, 0.0, 50) # alpha
    crystalball.SetParLimits(4, 0, 100) # n
    
    if True:
        n = 3
        alpha = 0.1
        sigma = 0.1
        seeds = [(n, mean, sigma, alpha, norm)]
        aN = [10, 100, 500, 1000]
        aMean = [0, 0.1, -0.1, -0.2, 0.2]
        aSigma = [0.05, 0.1, 0.2]
        aAlpha = [0.05, 0.1, 0.2, 0.6, 1, 2, 5, 10, 20]
        an = [1, 2, 3, 4, 10, 20, 40, 80]
        additional_seeds = create_n_seeds(100, aN, aMean, aSigma, aAlpha, an)
        seeds.extend(additional_seeds)
        fit_result = do_multi_fit(h, seeds)
    if False:
        # Parameters: alpha, n, mean, sigma, norm
        #crystalball.SetParameters(1.0, 3.0, mean, stddev, norm)
        # params are n, mean, sigma, alpha, N
        # power
        n = 3
        alpha = 0.1
        sigma = 0.1
        C = n/abs(alpha) * 1/(n-1) * math.exp(-0.5 * alpha**2)
        D = math.sqrt(math.pi*0.5)*(1 + ROOT.TMath.Erf(alpha)/math.sqrt(2.0))
        #norm = (sigma*(C+D))
        norm = 1 # h.Integral() * h.Integral()
        crystalball.SetParameters(n, mean, sigma, alpha, norm)
        print(n, mean, sigma, alpha, norm)
        #crystalball.FixParameter(0, 2.0)
        #crystalball.FixParameter(1, 0.0) 
        #crystalball.FixParameter(2, stddev)
        #crystalball.FixParameter(3, 0.1) 
        #crystalball.FixParameter(4, 1) # norm
        crystalball.SetParLimits(0, 0, 30) # n
        crystalball.SetParLimits(1, mean - 0.2, mean + 0.2) # mean
        crystalball.SetParLimits(2, 0.05, 1) # stddev
        crystalball.SetParLimits(3, 0.05, 3) # alpha
        crystalball.SetParLimits(4, 0, 40) # norm
        fit_result = h.Fit(crystalball, "RIES", "")
        #crystalball.FixParameter(4, 1)
        #crystalball.SetParLimits(4, 0, 10) # norm
        #fit_result = h.Fit(crystalball, "RIES", "")
        fit_result.Print("V")
    chi2 = crystalball.GetChisquare()
    ndf = crystalball.GetNDF()
    if ndf == 0: ndf = 1
    mean = crystalball.GetParameter(1)
    mean_error = crystalball.GetParError(1)
    stddev = crystalball.GetParameter(2)
    stddev_error = crystalball.GetParError(2)
    out_mean = mean * centers[i] + centers[i]
    out_error = stddev * centers[i]
    h.Draw("hist fit")
    stats = f"N Muons:\t{n_entries:0.0f}" # \nSimple mean:\t{simple_mean: 0.2f}\nSimple stddev:\t{simple_stddev: 0.2f}"
    stats += f"\n#chi^{{2}}/dof:\t{chi2:0.2f}/{ndf:0.0f}\n#bar{{x}}:\t{mean:0.4f}\n#sigma:\t{stddev: 0.4f}"
    stats += f"\nN:\t{crystalball.GetParameter(0):0.0f}\n#alpha:\t{crystalball.GetParameter(3):0.2f}\nn:\t{crystalball.GetParameter(4):0.2f}"
    textbox = add_text(0.575, 0.55, 0.83, 0.85, stats)
    crystalball.Draw("same")
    dunestyle.Simulation()
    print(f"#{{Chi}}^{2}/dof: {chi2:0.3f}\tndf: {ndf}\tchi2/ndf: {chi2/ndf:0.4f}\t" \
          f"mean: {mean:0.2f} +/- {mean_error:0.2f}\tsigma: {stddev:0.2f} +/- {stddev_error:0.2f}\t{int(fit_result)}")
    canvas.Print(outfilename)
    return out_mean, out_error
    


selected_sample_hist = tf.Get("energy_resolution__resolution__muon_starting_ke_fractional_resolution")
if selected_sample_hist != None:
    selected_sample_hist.SetLineColor(ROOT.kBlack)
    do_fit(selected_sample_hist, "selected_sample_resolution.png", 1)


selected_sample_hist = tf.Get("energy_resolution__lar_resolution__lar_muon_fractional_resolution")
if selected_sample_hist != None:
    selected_sample_hist.SetLineColor(ROOT.kBlack)
    do_fit(selected_sample_hist, "selected_sample_full_muon_resolution.png", 1)
    
for i in range(n_start, n_end + 1):
    histname = f"energy_resolution__resolution__slices__muon_{i}"
    print(f"\n\n## {histname} ##")
    h = tf.Get(histname)
    assert h != None, f"Couldn't get {histname}"
    title = h.GetTitle().replace("true KE", "True KE / GeV")
    title = title.replace(".000000", "").replace(".500000", ".5").replace(".750000", ".75").replace(".250000", ".25")
    h.SetTitle(title)
    h.SetLineColor(ROOT.kBlack)
    
    mean, error = do_fit(h, f"output/slice_{i:02d}.png", i) 
    means.append(mean)
    errors.append(error)
    
centers = centers[n_start:n_end+1]
assert len(centers) == len(means)
n = len(centers)
zero_error = n * [0.0]
as_array = lambda x: array.array("d", x)
col = tf.Get("energy_resolution__resolution__muon_starting_ke_resolution_column_normalized")
g = ROOT.TGraphErrors(n, as_array(centers), as_array(means), as_array(zero_error), as_array(errors))
fit_min = 0.5 # GeV
fit_max = 4 # GeV
f1 = ROOT.TF1("fit", "[0]*x+[1]", fit_min, fit_max)
g.Fit(f1, "", "", fit_min, fit_max);

p0 = f1.GetParameter(0)
p1 = f1.GetParameter(1)
m = 1/p0
b = -1000.0 * p1 / p0
print(f"Fit function:\n#define length_to_energy(l) (1.75*l*{m:0.4f} + {b:0.4f})*1e-3")

f1.SetLineColor(ROOT.kRed)
col.GetZaxis().SetRangeUser(0.001, col.GetMaximum()*1.0001)
col.Draw("colz")
dunestyle.Simulation()
g.Draw("P same")
canvas.Print("output.png")















