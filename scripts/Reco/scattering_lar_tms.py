 # use MomentumLArEnd PositionLArEnd MomentumTMSStart PositionTMSStart PositionLArStart
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
    """ Returns true if x,y,z are inside the TMS
    Currently rough estimates of the positions """
    is_inside = True
    if not -3500 < x < 3500: is_inside = False
    if not -3700 < y < 1000: is_inside = False
    if only_thin_section:
        if not 11000 < z < 13500: is_inside = False
    else:
        if not 11000 < z < 18200: is_inside = False
    return is_inside 

def inside_lar(x, y, z):
    """ Returns true if x,y,z are inside the TMS
    Currently rough estimates of the positions """
    is_inside = True
    if not -4500 < x < 3700: is_inside = False
    if not -3200 < y < 1000: is_inside = False
    if not 4100 < z < 9200: is_inside = False
    return is_inside 


def region1(x):
    is_region1= True
    if not -4000 < x < -2500: is_region1= False
    return is_region1

def region1_2(x):
    is_region1_2= True
    if not -2500 < x < -1500: is_region1_2= False
    return is_region1_2

def region2(x):
    is_region2= True
    if not -1500 < x < 1500: is_region2= False
    return is_region2

def region2_3(x):
    is_region2_3= True
    if not 1500 < x < 2500: is_region2_3= False
    return is_region2_3

def region3(x):
    is_region3= True
    if not 2500 < x < 4000: is_region3= False
    return is_region3


def run(c, truth, reco, outfilename, nmax=-1):
    
    # PositionTMSStart, MomentumTMSStart, and PositionTMSEnd
   
    bin_edges = array.array('d', [0, 200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,3000,4000,5000,6000,7000,8000,9000])
    
    if truth != None:
        hist_delta_theta = ROOT.TH1D("hist_delta_theta", "True absolute delta theta: abs(theta_tms - theta_lar) (using truth_info);True absolute angle difference(radian) ;Number of muons", 100, 0, 1)
        hist_delta_x = ROOT.TH1D("hist_delta_x", "True absolute delta x: abs(x_tms_start - x_extrapolate) (using truth_info);True absolute distance(mm) ;Number of muons", 100, 0, 2000)
        hist_delta_theta_signed = ROOT.TH1D("hist_delta_theta_signed", "Delta theta: (theta_tms - theta_lar) (using truth_info);True angle difference(radian) ;Number of muons", 100, -1, 1)
        hist_delta_x_signed = ROOT.TH1D("hist_delta_x_signed", "True delta x: (x_tms_start - x_extrapolate) (using truth_info);True distance(mm) ;Number of muons", 100, -2000, 2000)
        hist_theta_tms_signed = ROOT.TH1D("hist_theta_tms_signed", "True tms angle;True tms angle (radian) ;Number of muons", 100, -1, 1)
        hist_theta_lar_signed = ROOT.TH1D("hist_theta_lar_signed", "True lar angle;True angle (radian) ;Number of muons", 100, -1, 1)
        hist_delta_theta_antimuons = ROOT.TH1D("hist_delta_theta_antimuons", "True absolute delta theta: abs(theta_tms - theta_lar) (using truth_info);True absolute angle difference(radian) ;Number of antimuons", 100, 0, 1)
        hist_delta_x_antimuons = ROOT.TH1D("hist_delta_x_antimuons", "True absolute delta x: abs(x_tms_start - x_extrapolate) (using truth_info);True absolute distance(mm) ;Number of antimuons", 100, 0, 2000)
        hist_delta_theta_signed_antimuons = ROOT.TH1D("hist_delta_theta_signed_antimuons", "Delta theta: (theta_tms - theta_lar) (using truth_info);True angle difference(radian) ;Number of antimuons", 100, -1, 1)
        hist_delta_x_signed_antimuons = ROOT.TH1D("hist_delta_x_signed_antimuons", "True delta x: (x_tms_start - x_extrapolate) (using truth_info);True distance(mm) ;Number of antimuons", 100, -2000, 2000)
        hist_theta_tms_signed_antimuons = ROOT.TH1D("hist_theta_tms_signed_antimuons", "True tms angle;True tms angle (radian) ;Number of antimuons", 100, -1, 1)
        hist_theta_lar_signed_antimuons = ROOT.TH1D("hist_theta_lar_signed_antimuons", "True lar angle;True angle (radian) ;Number of antimuons", 100, -1, 1)
        
        hist_matched_muons_in_lar_xz = ROOT.TH2D("hist_matched_muons_in_lar_xz", 
           "Matched muons in lar xz plane (using truth_info);True muon start position Z (mm);True muon start position X (mm)", 50, 3600, 10000, 50, -5000, 4000)
        hist_matched_muons_in_lar_yz = ROOT.TH2D("hist_matched_muons_in_lar_yz", 
           "Matched muons in lar yz plane (using truth_info);True muon start position Z (mm);True muon start position Y (mm)", 50, 3600, 10000, 50, -4000, 1300)
        
        hist_all_muons_in_lar_xz = ROOT.TH2D("hist_all_muons_in_lar_xz", 
           "All muons in lar xz plane (using truth_info);True muon start position Z (mm);True muon start position X (mm)", 50, 3600, 10000, 50, -5000, 4000)
        hist_all_muons_in_lar_yz = ROOT.TH2D("hist_all_muons_in_lar_yz", 
           "All muons in lar yz plane (using truth_info);True muon start position Z (mm);True muon start position Y (mm)", 50, 3600, 10000, 50, -4000, 1300)
        

        hist_ratio_of_matched_muons_vs_all_muons_in_lar_xz = ROOT.TH2D("hist_ratio_of_matched_muons_vs_all_muons_in_lar_xz", 
           "Ratio of matched muons vs all muons in lar xz plane (using truth_info);True muon start position Z (mm);True muon start position X (mm)", 50, 3600, 10000, 50, -5000, 4000)
        hist_ratio_of_matched_muons_vs_all_muons_in_lar_yz = ROOT.TH2D("hist_ratio_of_matched_muons_vs_all_muons_in_lar_yz", 
           "Ratio of matched muons vs all muons in lar yz plane (using truth_info);True muon start position Z (mm);True muon start position Y (mm)", 50, 3600, 10000, 50, -4000, 1300)
        

        hist_all_muons_in_lar_xz_using_positionlarstart = ROOT.TH2D("hist_all_muons_in_lar_xz_using_positionlarstart", 
           "All muons in lar xz plane (using truth_info)(using_positionlarstart);True muon start position Z (mm);True muon start position X (mm)", 50, 3600, 10000, 50, -5000, 4000)
        hist_all_muons_in_lar_yz_using_positionlarstart = ROOT.TH2D("hist_all_muons_in_lar_yz_using_positionlarstart", 
           "All muons in lar yz plane (using truth_info)(using_positionlarstart);True muon start position Z (mm);True muon start position Y (mm)", 50, 3600, 10000, 50, -4000, 1300)
        
    
    # User can request fewer events, so check how many we're looping over.
    nevents = c.GetEntries()
    if nmax >= 0 and nevents > nmax: nevents = nmax
    # Figure out how often to print progress information.
    # Setting carriage = True will use a carriage return which keeps the progress on a single line
    # But if you add print statements, it will be ugly so it's not default
    carriage = False
    if carriage:
        if nevents <= 100: print_every = 1 
        elif nevents <= 1000: print_every = 10
        elif nevents <= 10000: print_every = 100
        else: print_every = 1000
    else:
        if nevents < 100: print_every = 1 # Print every event if < 100 events
        elif 100 <= nevents < 1000: print_every = 20
        elif 1000 <= nevents < 10000: print_every = 100
        elif 10000 <= nevents < 100000: print_every = 1000
        else: print_every = 10000
    
    n_true_muons = 0
    n_true_antimuons=0
    n_matched_muons =0
    n_matched_antimuons =0
    matched_percentage_muons=0
    matched_percentage_antimuons=0
    n_muons_lar_start=0
    


    # Now loop over all events
    for i, event in enumerate(c):
        if i > nevents: break
        if truth != None: 
            truth.GetEntry(i)
            
        if reco != None: 
            reco.GetEntry(i)
        # Print current progress, with carriage return \r to prevent long list of progress and have everything on a singe line.
        if i % print_every == 0 and carriage: print(f"\rOn {i} / {nevents}  {i/nevents*100}%", end='')
        if i % print_every == 0 and not carriage: print(f"On {i} / {nevents}   {i/nevents*100:.1f}%")
        
       
        
        
        # use MomentumLArStart PositionLArEnd MomentumTMSStart PositionTMSStart PositionLArStart
        for index, particle in enumerate(truth.PDG):
            if truth.PDG[index] == 13:

                x_start_lar = truth.PositionLArEnd[4*index+0]
                z_start_lar = truth.PositionLArEnd[4*index+2]
                
                x_start_tms = truth.PositionTMSStart[4*index+0]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                
               
                x_birth = truth.BirthPosition[4*index+0]
                y_birth = truth.BirthPosition[4*index+1]
                z_birth = truth.BirthPosition[4*index+2]
                
                x_death = truth.DeathPosition[4*index+0]
                y_death = truth.DeathPosition[4*index+1]
                z_death = truth.DeathPosition[4*index+2]
                
                
                
                if inside_lar(x_birth,y_birth,z_birth):
                    hist_all_muons_in_lar_xz.Fill(z_birth,x_birth)
                    hist_all_muons_in_lar_yz.Fill(z_birth,y_birth)

                    hist_all_muons_in_lar_xz_using_positionlarstart.Fill(truth.PositionLArStart[4*index+2],truth.PositionLArStart[4*index+0])
                    hist_all_muons_in_lar_yz_using_positionlarstart.Fill(truth.PositionLArStart[4*index+2],truth.PositionLArStart[4*index+1])
                    
                    n_muons_lar_start+=1


                if inside_lar(x_birth,y_birth,z_birth) and inside_tms(x_death,y_death,z_death):
                    p_z = truth.MomentumLArEnd[4*index+2]
                    p_x = truth.MomentumLArEnd[4*index+0]
                    p_z_tms = truth.MomentumTMSStart[4*index+2]
                    p_x_tms = truth.MomentumTMSStart[4*index+0]
                    if p_z ==0 or p_z_tms ==0 : continue
                    if p_z <0 or p_z_tms <0 : continue
                    m= p_x / p_z 
                    m2= p_x_tms/p_z_tms
                    theta_tms = math.atan(m2)
                    theta_lar = math.atan(m)
                    #if abs(theta_tms)> 0.7 or abs(theta_lar)>0.7 : continue

                    n_true_muons+=1
                    x_extrapolate =m*(z_start_tms-z_start_lar) + x_start_lar
                    delta_x= abs(x_start_tms-x_extrapolate)
                    delta_x_signed= x_start_tms-x_extrapolate
                    hist_delta_x_signed.Fill(delta_x_signed)
                    hist_delta_x.Fill(delta_x)
                    

                    hist_theta_tms_signed.Fill(theta_tms)
                    hist_theta_lar_signed.Fill(theta_lar)
                    delta_theta =abs(theta_tms-theta_lar)
                    delta_theta_signed =theta_tms-theta_lar
                    hist_delta_theta.Fill(delta_theta)
                    hist_delta_theta_signed.Fill(delta_theta_signed)
                    if delta_x < 400 and delta_theta < 0.31 :  
                        n_matched_muons += 1
                        hist_matched_muons_in_lar_xz.Fill(z_birth,x_birth)
                        hist_matched_muons_in_lar_yz.Fill(z_birth,y_birth)

            
            
            
            if truth.PDG[index] == -13:

                x_start_lar = truth.PositionLArEnd[4*index+0]
                z_start_lar = truth.PositionLArEnd[4*index+2]
                
                x_start_tms = truth.PositionTMSStart[4*index+0]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                
               
                x_birth = truth.BirthPosition[4*index+0]
                y_birth = truth.BirthPosition[4*index+1]
                z_birth = truth.BirthPosition[4*index+2]
                
                x_death = truth.DeathPosition[4*index+0]
                y_death = truth.DeathPosition[4*index+1]
                z_death = truth.DeathPosition[4*index+2]


                if inside_lar(x_birth,y_birth,z_birth) and inside_tms(x_death,y_death,z_death):
                    p_z = truth.MomentumLArEnd[4*index+2]
                    p_x = truth.MomentumLArEnd[4*index+0]
                    p_z_tms = truth.MomentumTMSStart[4*index+2]
                    p_x_tms = truth.MomentumTMSStart[4*index+0]
                    if p_z ==0 or p_z_tms ==0 : continue
                    if p_z <0 or p_z_tms <0 : continue
                    m= p_x / p_z 
                    m2= p_x_tms/p_z_tms
                    theta_tms = math.atan(m2)
                    theta_lar = math.atan(m)
                    #if abs(theta_tms)> 0.7 or abs(theta_lar)>0.7 : continue

                    n_true_antimuons+=1
                    x_extrapolate =m*(z_start_tms-z_start_lar) + x_start_lar
                    delta_x= abs(x_start_tms-x_extrapolate)
                    delta_x_signed= x_start_tms-x_extrapolate
                    hist_delta_x_signed_antimuons.Fill(delta_x_signed)
                    hist_delta_x_antimuons.Fill(delta_x)
                    

                    hist_theta_tms_signed_antimuons.Fill(theta_tms)
                    hist_theta_lar_signed_antimuons.Fill(theta_lar)
                    delta_theta =abs(theta_tms-theta_lar)
                    delta_theta_signed =theta_tms-theta_lar
                    hist_delta_theta_antimuons.Fill(delta_theta)
                    hist_delta_theta_signed_antimuons.Fill(delta_theta_signed)
                    if delta_x < 400 and delta_theta < 0.31 :  n_matched_antimuons += 1           


                     

   
    matched_percentage_muons = n_matched_muons/n_true_muons
    matched_percentage_antimuons = n_matched_antimuons/n_true_antimuons
    hist_ratio_of_matched_muons_vs_all_muons_in_lar_xz.Divide(hist_matched_muons_in_lar_xz,hist_all_muons_in_lar_xz)
    hist_ratio_of_matched_muons_vs_all_muons_in_lar_yz.Divide(hist_matched_muons_in_lar_yz,hist_all_muons_in_lar_yz)
    
    ## Now save all the histograms
    if carriage: print("\r", end="") # Erase the current line
    print(f"Completed loop of {nevents} events. Saving to {outfilename}")
    hists_to_save = []
    # This is a pythonistic cheat to quickly find the histograms
    # and add them to the list of files to save
    # Otherwise we'd add them one by one, like: hists_to_save.append(hist1), etc.
    for name, item in locals().items():
        if "hist_" in name:
            hists_to_save.append(item)
    
    tf = ROOT.TFile(outfilename, "recreate")
    for hist in hists_to_save:
        hist.Write()
    tf.Close()
    print("Done saving.")
    print(f"N true muons: {n_true_muons}")
    print(f"N true antimuons: {n_true_antimuons}")
    
    


    
    # Return this hists if the user requested previews
    return hists_to_save 

def new_func(n_true_muons, n_true_antimuons, correct_percentage_total):
    print(f"N true muons: {n_true_muons}")
    print(f"N true antimuons: {n_true_antimuons}")
    print(f"total correct percentage for muons: {correct_percentage_total}")
    

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
            for name in files:
                if "tmsreco.root" in name:
                    # This is a tmsreco file. So add it to the list
                    fullfilename = os.path.join(root, name)
                    files_to_use.append(fullfilename)
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
            canvas.Print(os.path.join(preview_dir, f"{name}.png"))
    
    
        
        
        
  
    

# This is how python knows whether it's a script that's being imported into another script or running as the main script
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', type=str, help="The output dir. Will be made if it doesn't exist.", default="")
    parser.add_argument('--name', type=str, help="The name of the output files.", default="dune-tms_hists.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', type=str, help="The location of the input tmsreco files", default="")
    group.add_argument('--inlist', type=str, help="The input filelist", default="")
    group.add_argument('--filename', '-f', type=str, help="The input file, if you have a single file", default="")
    parser.add_argument('--nevents', '-n', type=int, help="The maximum number of events to loop over", default=-1)
    parser.add_argument('--allow_overwrite', help="Allow the output file to overwrite", action=argparse.BooleanOptionalAction)
    parser.add_argument('--preview', help="Save preview images of the histograms", action=argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    validate_then_run(args)
    

