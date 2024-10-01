import os
import collections

import ROOT
ROOT.gROOT.SetBatch(True)

import sys

def get_subdir_and_name(hist_name):
    subdir = ""
    name = hist_name.strip()
    split = name.split("__")
    if len(split) > 2:
        # Subdir
        subdir = split[0].replace("_", "/")
        name = split[1]
    else:
        # Simple subdir
        split = name.split("_")
        subdir = split[0]
        # And the rest is the name
        name = "_".join(split[1:])
    print(hist_name, subdir, name)
    return subdir, name
    

def draw_histograms(input_file):
    # Open the input ROOT file
    root_file = ROOT.TFile.Open(input_file)

    # Create a directory to save images
    output_dir = os.path.splitext(input_file)[0] + "_images"
    os.makedirs(output_dir, exist_ok=True)
    
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    
    recoeff_plots_numerators = dict()
    recoeff_plots_denominators = dict()
    
    stack_plots = collections.defaultdict(dict)

    # Loop over all keys in the ROOT file
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        output_subdir = output_dir
        name = obj.GetName()
            
        subdir, image_name = get_subdir_and_name(name)
        output_subdir = os.path.join(output_dir, subdir)
            
        # Can add reco eff
        reco_eff = False
        if "numerator" in name or "denominator" in name: reco_eff = True
        stack = False
        if "stack" in name: stack = True
        if reco_eff or stack:
            output_subdir = os.path.join(output_subdir, "additional_plots")
        os.makedirs(output_subdir, exist_ok=True)  
        if isinstance(obj, ROOT.TH2):
            # For 2D histograms, draw with "colz" option and save as png
            obj.GetYaxis().SetTitleOffset(1.4)
            obj.GetZaxis().SetTitleOffset(0.5)
            obj.Draw("colz")
            print(f"{obj.GetName()} integral: {obj.Integral()}")
            canvas.Print(os.path.join(output_subdir, image_name + ".png"))
        elif isinstance(obj, ROOT.TH1):
            # For 1D histograms, draw and save as png
            top = obj.GetMaximum()*1.2
            obj.GetYaxis().SetRangeUser(0, top)
            obj.Draw()
            #print(f"{obj.GetName()} integral: {obj.Integral()}")
            canvas.Print(os.path.join(output_subdir, image_name + ".png"))
        if reco_eff:
            key = name.replace("_numerator", "").replace("_denominator", "")
            if "numerator" in name: recoeff_plots_numerators[key] = obj 
            if "denominator" in name: recoeff_plots_denominators[key] = obj 
        if stack:
            split_stack = name.split("_stack_")
            stack_key = split_stack[0]
            stack_plots[stack_key][split_stack[1]] = obj
            stack_plots[stack_key + "_log"][split_stack[1]] = obj
    
    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen, ROOT.kMagenta, ROOT.kBlack]
    line_styles = [1, 2, 7, 9]
    for name, hist_and_name in stack_plots.items():
        print(f"Doing stack plot {name}")
        index = 0
        first_hist = None
        ymax = 0
        hist_stack = ROOT.THStack()
        
        log = False
        if "log" in name: log = True
        if log: canvas.SetLogy(True)
        else: canvas.SetLogy(False)
        
        headroom = 1.1
        if log: headroom = 6
        for item_name, hist in hist_and_name.items():
            print(f"Doing stack plot subhist {item_name}. Integral: {hist.Integral()}")
            hist.SetLineColor(colors[index % len(colors)])
            hist.SetLineStyle(line_styles[index % len(line_styles)])
            #hist.DrawCopy("" if index == 0 else "same")
            if index == 0: first_hist = hist
            ymax = max(ymax, hist.GetMaximum())
            hist_stack.Add(hist)
            if log: hist.GetYaxis().SetRangeUser(10, ymax*headroom)
            index += 1
        #if first_hist != None:
        #    if not log: first_hist.GetYaxis().SetRangeUser(0, ymax*headroom)
        print(f"ymax: {ymax}")
        hist_stack.Draw("nostack")
        subdir, image_name = get_subdir_and_name(name)
        output_subdir = os.path.join(output_dir, subdir)
        
        outfilename = os.path.join(output_subdir, image_name + ".png")
        print(f"Saving in {outfilename}")
        canvas.Print(outfilename)
        
    # Turn off log y if it's already on
    canvas.SetLogy(False)
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
        
        
        name = numerator.GetName()
        name = name.replace("_numerator", "")
        subdir, image_name = get_subdir_and_name(name)
        output_subdir = os.path.join(output_dir, subdir)
        canvas.Print(os.path.join(output_subdir, image_name + ".png"))

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
