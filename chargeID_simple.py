import glob
import argparse
import os
import math
import ROOT # type: ignore
import array
# Tells root to be in batch mode so it doesn't try to create a canvas on your screen, which is slow
ROOT.gROOT.SetBatch(True)
# Don't draw the stats box on histograms
ROOT.gStyle.SetOptStat(0)
# Python's wrapper around ROOT sometimes crashes when it automatically adds hists to its garbage collector. This fixed that
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

def inside_tms(x, y, z, only_thin_section = False):
  
    is_inside = True
    if not -3730 < x < 3730: is_inside = False
    if not -3702 < y < 998: is_inside = False
    if only_thin_section:
        if not 11333 < z < 13500: is_inside = False
    else:
        if not 11333 < z < 18535: is_inside = False
    return is_inside 

def inside_lar(x, y, z):
   
    is_inside = True
    if not -4500 < x < 3700: is_inside = False
    if not -3200 < y < 1000: is_inside = False
    if not 4100 < z < 9200: is_inside = False
    return is_inside 


def region1(x):
    is_region1= True
    if not -3520 < x < -1750: is_region1= False
    return is_region1

# def region1_2(x):
#     is_region1_2= True
#     if not -2500 < x < -1500: is_region1_2= False
#     return is_region1_2

def region2(x):
    is_region2= True
    if not -1750 < x < 1750: is_region2= False
    return is_region2

# def region2_3(x):
#     is_region2_3= True
#     if not 1500 < x < 2500: is_region2_3= False
#     return is_region2_3

def region3(x):
    is_region3= True
    if not 1750 < x < 3520: is_region3= False
    return is_region3


def run(c, truth, reco, outfilename, nmax=-1):
  
    hist_correct_charge_id = ROOT.TH1D("hist_correct_charge_id", "Muons correct charge number (using reco);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_incorrect_charge_id = ROOT.TH1D("hist_incorrect_charge_id", "Muons incorrect charge number plot(using reco);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_total_charge_id = ROOT.TH1D("hist_total_charge_id", "Muons total charge number (using reco);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_correct_charge_id_percentage = ROOT.TH1D("hist_correct_charge_id_percentage", "Muons correct charge ID percentage (using reco);True muon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage.SetMaximum(1.2)  # Upper limit of Y-axis

    hist_correct_charge_id_antimuon = ROOT.TH1D("hist_correct_charge_id_antimuon", "Antimuons correct charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_incorrect_charge_id_antimuon = ROOT.TH1D("hist_incorrect_charge_id_antimuon", "Antimuons incorrect charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_total_charge_id_antimuon = ROOT.TH1D("hist_total_charge_id_antimuon", "Antimuons total charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon", "Antimuons correct charge ID percentage (using reco);True antimuon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_antimuon.SetMaximum(1.2)  # Upper limit of Y-axis

    nevents = c.GetEntries()
    if nevents < 100: print_every = 1 # Print every event if < 100 events
    elif 100 <= nevents < 1000: print_every = 20
    elif 1000 <= nevents < 10000: print_every = 100
    elif 10000 <= nevents < 100000: print_every = 1000
    else: print_every = 10000
    if nmax >= 0 and nevents > nmax: nevents = nmax
    
    n_true_muons=0
    n_true_antimuons=0

    n_correct_charge_muon_geo = 0
    n_incorrect_charge_muon_geo = 0
    
    n_correct_charge_antimuon_geo = 0
    n_incorrect_charge_antimuon_geo = 0
    
    n_correct_charge_muon_geo_per = 0
    n_correct_charge_antimuon_geo_per = 0
    
    
    # Now loop over all events
    lastspill = -1
    for i, event in enumerate(c):
        if i > nevents: break
        if i % print_every == 0 : print(f"On {i} / {nevents}")
        if truth != None: 
            truth.GetEntry(i)
            
        if reco != None: 
            reco.GetEntry(i)

        if not truth.IsCC:
            continue

        # on_new_spill = False
        
        # if lastspill != reco.SpillNo:
        #     on_new_spill = True
        # lastspill = reco.SpillNo

        # if on_new_spill:
        #    continue    
        index = 0
        #for index, value in enumerate(truth.RecoTrackPrimaryParticlePDG):
        # Pick out the muons.
        if len(truth.RecoTrackPrimaryParticlePDG) <= 0: continue
        if truth.RecoTrackPrimaryParticlePDG[index] == 13:
            n_true_muons+=1
            x_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+0]
            y_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+1]
            z_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+2]
            x_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+0]
            y_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+1]
            z_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+2]
            p_z = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+2] 
            p_y = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+1]
            p_x = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+0]
            
            
            KE = math.sqrt(p_x**2 + p_y**2 +p_z**2 +105.7**2) - 105.7
            #KE = math.sqrt((truth.Muon_TrueKE + 105.7) ** 2 - 105.7 ** 2)
            
            #if inside_lar(x_start,y_start,z_start) and inside_tms(truth.Muon_Death[0], truth.Muon_Death[1], truth.Muon_Death[2]):
            if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                if KE <10:
                    continue
                if len(reco.TrackHitPos) < 3:               
                    continue
                nTracks = reco.nTracks
                if nTracks >0:
                    charge_value = reco.Charge[index]
                    # Selected Muon Energy Spectrum
                    hist_total_charge_id.Fill(KE)
                    
                    if charge_value == 13: 
                        hist_correct_charge_id.Fill(KE)
                        n_correct_charge_muon_geo+=1
                    else: 
                        n_incorrect_charge_muon_geo += 1 
                        hist_incorrect_charge_id.Fill(KE)   

        # Pick out the antimuons.            
        if truth.RecoTrackPrimaryParticlePDG[index] == -13:
            n_true_antimuons+=1
            x_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+0]
            y_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+1]
            z_start = truth.RecoTrackPrimaryParticleTruePositionStart[4*index+2]
            x_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+0]
            y_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+1]
            z_end = truth.RecoTrackPrimaryParticleTruePositionEnd[4*index+2]
            p_z = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+2] 
            p_y = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+1]
            p_x = truth.RecoTrackPrimaryParticleTrueMomentumEnteringTMS[4*index+0]
            
            KE = math.sqrt(p_x**2 + p_y**2 +p_z**2 +105.7**2) - 105.7
            #KE = math.sqrt((truth.Muon_TrueKE + 105.7) ** 2 - 105.7 ** 2)
            
            #if inside_lar(x_start,y_start,z_start) and inside_tms(truth.Muon_Death[0], truth.Muon_Death[1], truth.Muon_Death[2]):
            if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                if KE <10:
                    continue
                if len(reco.TrackHitPos) < 3:               
                    continue
                
                nTracks = reco.nTracks
                if nTracks >0:
                    charge_value = reco.Charge[index]
                    hist_total_charge_id_antimuon.Fill(KE)
                    
                    if charge_value == -13: 
                        hist_correct_charge_id_antimuon.Fill(KE)
                        n_correct_charge_antimuon_geo+=1

                    else:
                        n_incorrect_charge_antimuon_geo +=1
                        hist_incorrect_charge_id_antimuon.Fill(KE)
   
    if n_correct_charge_muon_geo >0:  n_correct_charge_muon_geo_per = n_correct_charge_muon_geo/(n_incorrect_charge_muon_geo +n_correct_charge_muon_geo)
    if n_correct_charge_antimuon_geo >0:  n_correct_charge_antimuon_geo_per = n_correct_charge_antimuon_geo/(n_incorrect_charge_antimuon_geo+ n_correct_charge_antimuon_geo)

    
    # Divide histograms to compute percentage histograms
    hist_correct_charge_id_percentage.Divide(hist_correct_charge_id, hist_total_charge_id, 1.0, 1.0, "B")
    hist_correct_charge_id_percentage_antimuon.Divide(hist_correct_charge_id_antimuon, hist_total_charge_id_antimuon, 1.0, 1.0, "B")
   

    print(f"Completed loop of {nevents} events. Saving to {outfilename}")
    hists_to_save = []
   
    for name, item in locals().items():
        if "hist_" in name:
            hists_to_save.append(item)
    
    tf = ROOT.TFile(outfilename, "recreate")
    for hist in hists_to_save:
        hist.Write()
    tf.Close()
    print("Done saving.")
    
    

    outdir = args.outdir
    if outdir == "":
        # No output directory was specified, so use the default
        username = os.environ["USER"]
        outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
    else:
        # Check if it follows the correct conventions
        good_locations = ["/dune/data/users", "/dune/data2/users", "/pnfs/dune/persistent", "/pnfs/dune/scratch"]
        if not any(location in outdir for location in good_locations):
            print(f"Warning: outdir is not in list of good locations. Don't want to write root files to app area. {outdir}")

    # Make sure the output directory exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Build output file path
    outname = args.name
    base_name = outname[:-5]  # Remove the last 5 characters (.root)
    output_filename = base_name + '.txt'
    file_path = os.path.join(outdir, output_filename)

    print(f"Output text file will be: {file_path}")

    # Open the file for writing
    with open(file_path, 'w') as file:
    # Write the output to the file
        file.write(f"N true muons: {n_true_muons}\n")
        file.write(f"N true antimuons: {n_true_antimuons}\n")
        
        file.write(f"correct chargeID in reco (for muon): {n_correct_charge_muon_geo}\n")
        file.write(f"incorrect chargeID in reco (for muon): {n_incorrect_charge_muon_geo}\n")
        file.write(f"correct percentage geoPDG in reco (for muon): {n_correct_charge_muon_geo_per}\n")

        file.write(f"correct chargeID in reco (for antimuon): {n_correct_charge_antimuon_geo}\n")
        file.write(f"incorrect chargeID in reco (for antimuon): {n_incorrect_charge_antimuon_geo}\n")
        file.write(f"correct percentage geoPDG in reco (for antimuon): {n_correct_charge_antimuon_geo_per}\n")

    
    # Return this hists if the user requested previews
    return hists_to_save 



def validate_then_run(args):
    ''' This function has all the code to validate the arguments. Usually users don't need to edit this code. '''
    # First we validate the args
    indir = args.indir
    inlist = args.inlist
    infile = args.filename
    files_to_use = []
    if indir != "":
        # User specified an input directory. So loop through all the files to find tmsreco files.
        print(f"Finding all tmsreco.root files in {indir}")
        for root, dirs, files in os.walk(indir):
            filenum = 0
            for name in files:
                if "root" in name and filenum < 100:
                    # This is a tmsreco file. So add it to the list
                    fullfilename = os.path.join(root, name)
                    files_to_use.append(fullfilename)
                    filenum += 1
        nfiles = len(files_to_use)
        if nfiles == 0:
            raise ValueError(f"Did not find any files in {indir}")
        print(f"Found {nfiles} files in {indir}")
    if inlist != "":
        # In this case the user specified a text file with the full paths
        with open(inlist) as f:
            file_data = f.read()
            files_to_use = file_data.splitlines()
        nfiles = len(files_to_use)
        if nfiles == 0:
            raise ValueError(f"Did not find any files in {inlist}")
        print(f"Found {nfiles} files in {inlist}")
    if infile != "":
        # In this case, the user specified exactly one file. Usually they'd hadd many files together.
        files_to_use = [infile]
        
    outdir = args.outdir
    if outdir == "":
        # No output directory was specified so use the default
        # First we need the username
        username = os.environ["USER"]
        outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
    else:
        # Check if it follows the correct conventions
        good_locations = ["/dune/data/users", "/dune/data2/users", "/pnfs/dune/persistent", "/pnfs/dune/scratch"]
        if not any(location in outdir for location in good_locations):
            print(f"Warning: outdir is not in list of good locations. Don't want to write root files to app area. {outdir}")
    # Make sure the output directory exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outname = args.name
    if ".root" not in outname:
        raise ValueError(f"The output file should be a root file. {outname}")
    outfilename = os.path.join(outdir, outname)
    if os.path.exists(outfilename):
        if args.allow_overwrite:
            print(f"Warning: A file already exists at {outfilename}, but you specified allow_overwrite. It will be overwritten")
        else:
            raise ValueError(f"The output file already exists at {outfilename}. Please specify another location, allow_overwrite or delete the file yourself.")
    print(f"Output will be in {outfilename}")
    
    has_truth = True
    has_reco = True
    
    # Make the TChain objects. One for truth information and one for reconstructed information.
    c = ROOT.TChain("Line_Candidates")
    truth = None
    if has_truth: truth = ROOT.TChain("Truth_Info")
    for f in files_to_use:
        c.Add(f)
        if has_truth: truth.Add(f)
    assert c.GetEntries() > 0, "Didn't get any entries in Line_Candidates TChain." 
    if has_truth: assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info TChain."

    reco = None
    if has_reco: reco = ROOT.TChain("Reco_Tree")
    for f in files_to_use:
        c.Add(f)
        if has_reco: reco.Add(f)
    assert c.GetEntries() > 0, "Didn't get any entries in Line_Candidates TChain." 
    if has_reco: assert reco.GetEntries() > 0, "Didn't get any entries in Reco_Tree TChain."
    
    nevents = args.nevents
    assert nevents >= -1, f"nevents <= -1, why? {nevents}"
    
    # Now finally run
    hists = run(c, truth, reco, outfilename, nevents)

    preview = args.preview
    if preview:
        print("Making preview histograms.")
        print("WARNING: Eventually you'll want to make plots using a separate plotting script")
        print("That way you can format the plots without having to rerun the make_hists.py script")
        # Make the output directory if needed
        outname_without_root = outname.replace(".root", "")
        preview_dir = os.path.join(outdir, "previews", outname_without_root)
        if not os.path.exists(preview_dir): 
            os.makedirs(preview_dir)
        print(f"Saving previews in {preview_dir}/<hist_name>.png")
        canvas = ROOT.TCanvas()
        for hist in hists:
            name = hist.GetName().replace("hist_", "")
            opts = "colz"
            hist.Draw(opts)
            # hist.Draw("HIST E1")
            canvas.Print(os.path.join(preview_dir, f"{name}.png"))

# This is how python knows whether it's a script that's being imported into another script or running as the main script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', type=str, help="The output dir. Will be made if it doesn't exist.", default="")
    parser.add_argument('--name', type=str, help="The name of the output files.", default="dune-tms_hists.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', '-dir', type=str, help="The location of the input tmsreco files", default="")
    group.add_argument('--inlist', type=str, help="The input filelist", default="")
    group.add_argument('--filename', '-f', type=str, help="The input file, if you have a single file", default="")
    parser.add_argument('--nevents', '-n', type=int, help="The maximum number of events to loop over", default=-1)
    parser.add_argument('--allow_overwrite', help="Allow the output file to overwrite", action=argparse.BooleanOptionalAction)
    parser.add_argument('--preview', help="Save preview images of the histograms", action=argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    validate_then_run(args)
    
