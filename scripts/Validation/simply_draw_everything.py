import os

import ROOT
ROOT.gROOT.SetBatch(True)

import sys

def draw_histograms(input_file):
    # Open the input ROOT file
    root_file = ROOT.TFile.Open(input_file)

    # Create a directory to save images
    output_dir = os.path.splitext(input_file)[0] + "_images"
    os.makedirs(output_dir, exist_ok=True)
    
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    
    recoeff_plots_numerators = dict()
    recoeff_plots_denominators = dict()

    # Loop over all keys in the ROOT file
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        output_subdir = output_dir
        name = obj.GetName()
        split = name.split("_")
        if len(split) > 2:
            output_subdir = os.path.join(output_dir, split[0])
        reco_eff = False
        if "numerator" in name or "denominator" in name: reco_eff = True
        if reco_eff:
            output_subdir = os.path.join(output_subdir, "additional_plots")
        os.makedirs(output_subdir, exist_ok=True)  
        if isinstance(obj, ROOT.TH2):
            # For 2D histograms, draw with "colz" option and save as png
            obj.GetYaxis().SetTitleOffset(1.4)
            obj.GetZaxis().SetTitleOffset(0.5)
            obj.Draw("colz")
            print(f"{obj.GetName()} integral: {obj.Integral()}")
            canvas.Print(os.path.join(output_subdir, obj.GetName() + ".png"))
        elif isinstance(obj, ROOT.TH1):
            # For 1D histograms, draw and save as png
            top = obj.GetMaximum()*1.2
            obj.GetYaxis().SetRangeUser(0, top)
            obj.Draw()
            #print(f"{obj.GetName()} integral: {obj.Integral()}")
            canvas.Print(os.path.join(output_subdir, obj.GetName() + ".png"))
        if reco_eff:
            key = name.replace("_numerator", "").replace("_denominator", "")
            if "numerator" in name: recoeff_plots_numerators[key] = obj 
            if "denominator" in name: recoeff_plots_denominators[key] = obj 
            
    all_names = set(recoeff_plots_numerators.keys()) & set(recoeff_plots_denominators.keys())
    for name in all_names:
        error = False
        if name not in recoeff_plots_numerators:
            print(f"Didn't find {name} in {recoeff_plots_numerators}")
            error = True
        if name not in recoeff_plots_denominators:
            print(f"Didn't find {name} in {recoeff_plots_denominators}")
            error = True
        if error:
            print(f"Had one or more errors, skipping {name}")
            continue
        numerator = recoeff_plots_numerators[name]
        denominator = recoeff_plots_denominators[name]
        newtitle = numerator.GetTitle()
        newtitle = newtitle.replace("Numerator", "").strip()
        numerator.SetTitle(newtitle)
        numerator.Divide(denominator)
        numerator.GetYaxis().SetRangeUser(0, 1.2)
        numerator.GetYaxis().SetTitle("Reconstruction Efficiency")
        numerator.Draw()
        
        
        output_subdir = output_dir
        name = numerator.GetName()
        name = name.replace("_numerator", "")
        split = name.split("_")
        if len(split) > 1:
            output_subdir = os.path.join(output_dir, split[0])
        canvas.Print(os.path.join(output_subdir, name + ".png"))

    # Close the input ROOT file
    root_file.Close()



if __name__ == "__main__":
    # Check if input file is provided as argument
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_root_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.isfile(input_file):
        print("Error: Input file does not exist!")
        sys.exit(1)

    # Initialize ROOT
    ROOT.gROOT.SetBatch(True)  # Prevent ROOT from trying to open X11 windows
    ROOT.gStyle.SetOptStat(0)  # Hide statistics box in histograms

    # Call function to draw histograms
    draw_histograms(input_file)
