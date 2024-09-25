import os

import ROOT
ROOT.gROOT.SetBatch(True)

import sys

def draw_histograms(input_file1, input_file2):
    # The filename can include the key separated with comma, so check for that
    foo = input_file1.split(",")
    print(input_file1)
    if len(foo) > 1: 
        input_file1 = foo[0]
        key1 = foo[1]
    foo = input_file2.split(",")
    if len(foo) > 1: 
        input_file2 = foo[0]
        key2 = foo[1]
        
    # Open the input ROOT files
    root_file = ROOT.TFile.Open(input_file1)
    root_file2 = ROOT.TFile.Open(input_file2)

    # Create a directory to save images
    output_dir = os.path.join(os.path.split(input_file1)[0], f"comparison_{key1}_{key2}_images")
    os.makedirs(output_dir, exist_ok=True)
    

    # Loop over all keys in the ROOT file
    for key in root_file.GetListOfKeys():
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        obj = key.ReadObj()
        output_subdir = output_dir
        name = obj.GetName()
        obj2 = root_file2.Get(name)
        if obj2 == None:
            print(f"Could not find {name} in {input_file2}. Skipping")
            continue
        
        # Create a uniq subdirectory if it's a complex name
        split = name.split("_")
        if len(split) > 2:
            output_subdir = os.path.join(output_dir, split[0])
        os.makedirs(output_subdir, exist_ok=True)  
        
        if isinstance(obj, ROOT.TH2):
            canvas.cd(0)
            canvas.Divide(2, 1)
            canvas.cd(1)
            # For 2D histograms, draw with "colz" option and save as png
            obj.GetYaxis().SetTitleOffset(1.4)
            obj.GetZaxis().SetTitleOffset(0.5)
            obj.Draw("colz")
            canvas.cd(2)
            # For 2D histograms, draw with "colz" option and save as png
            obj2.GetYaxis().SetTitleOffset(1.4)
            obj2.GetZaxis().SetTitleOffset(0.5)
            obj2.Draw("colz")
            canvas.Print(os.path.join(output_subdir, obj.GetName() + ".png"))
        elif isinstance(obj, ROOT.TH1):
            canvas.cd(0)
            # For 1D histograms, draw and save as png
            leg = ROOT.TLegend(0.2,0.8,0.8,0.9)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetNColumns(2)
            leg.AddEntry(obj, key1, "lep")
            leg.AddEntry(obj2, key2, "lep")
            
            top = max(obj.GetMaximum(), obj2.GetMaximum())*1.3
            obj.GetYaxis().SetRangeUser(0, top)
            color1 = ROOT.kBlack
            obj.SetLineColor(color1)
            obj.Draw("e")
            color2 = 46
            obj2.SetLineColor(color2)
            obj2.SetMarkerColor(color2)
            obj2.SetMarkerStyle(22)
            obj2.SetLineStyle(2)
            obj2.Draw("e same")
            leg.Draw()
            #print(f"{obj.GetName()} integral: {obj.Integral()}")
            canvas.Print(os.path.join(output_subdir, obj.GetName() + ".png"))

    # Close the input ROOT file
    root_file.Close()



if __name__ == "__main__":
    # Check if input file is provided as argument
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_root_file>,<name (optional)> <input_root_file>,<name (optional)>")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]
    error = False
    if not os.path.isfile(input_file1.split(",")[0]):
        print(f"Error: Input file {input_file1} does not exist!")
        error = True
    if not os.path.isfile(input_file2.split(",")[0]):
        print(f"Error: Input file {input_file2} does not exist!")
        error = True
    if error:
        sys.exit(1)

    # Initialize ROOT
    ROOT.gROOT.SetBatch(True)  # Prevent ROOT from trying to open X11 windows
    ROOT.gStyle.SetOptStat(0)  # Hide statistics box in histograms

    # Call function to draw histograms
    draw_histograms(input_file1, input_file2)
