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
    """ This code does 3 things:
    Makes histograms
    Loops over all events and fills histograms
    Saves histograms in outfilename
    """
    # use PositionTMSStart, MomentumTMSStart, and PositionTMSEnd
    
    
    bin_edges = array.array('d', [0, 200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,3000,4000,5000,6000,7000,8000,9000])
    
    if truth != None and reco != None:
       
        hist_correct_charge_id = ROOT.TH1D("hist_correct_charge_id", "Muons correct charge number (using reco);True muon KE (MeV);Number of muons", 100, 0, 5000)
        hist_incorrect_charge_id = ROOT.TH1D("hist_incorrect_charge_id", "Muons incorrect charge number plot(using reco);True muon KE (MeV);Number of muons", 100, 0, 5000)
        hist_total_charge_id = ROOT.TH1D("hist_total_charge_id", "Muons total charge number (using reco);True muon KE (MeV);Number of muons", 100, 0, 5000)
        hist_correct_charge_id_percentage = ROOT.TH1D("hist_correct_charge_id_percentage", "Muons correct charge percentage (using reco);True muon KE (MeV);Fraction", 100, 0, 5000)
        hist_signed_distance = ROOT.TH1D("hist_signed_distance", "Muons signed distance: (x_extropolate - x_reco) (using reco);True signed distance(mm) ;Number of muons", 100, -2000, 2000)
        
        hist_correct_charge_id_antimuon = ROOT.TH1D("hist_correct_charge_id_antimuon", "Antimuons correct charge percentage (using reco);True antimuon KE (MeV);Number of antimuons", 100, 0, 5000)
        hist_incorrect_charge_id_antimuon = ROOT.TH1D("hist_incorrect_charge_id_antimuon", "Antimuons incorrect charge percentage (using reco);True antimuon KE (MeV);Number of antimuons", 100, 0, 5000)
        hist_total_charge_id_antimuon = ROOT.TH1D("hist_total_charge_id_antimuon", "Antimuons total charge number (using reco);True antimuon KE (MeV);Number of antimuons", 100, 0, 5000)
        hist_correct_charge_id_percentage_antimuon = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon", "Antimuons correct charge percentage (using reco);True antimuon KE (MeV);Fraction", 100, 0, 5000)
        hist_signed_distance_antimuon = ROOT.TH1D("hist_signed_distance_antimuon", "Antimuons signed distance: (x_extropolate - x_reco) (using reco) ;True signed distance(mm); Number of antimuons", 100, -2000, 2000)
        
        
        hist_tracklength_vs_starting_KE = ROOT.TH2D("hist_tracklength_vs_starting_KE", "Muons tracklength vs starting KE (using truth_info);True starting KE (MeV);True tracklength(mm) ", 100, 0, 6000,100,0,10000)
        


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
    
    

    n_true_muons=0
    n_true_antimuons=0
    n_muon_total_lar_start_tms_end =0
    n_antimuon_total_lar_start_tms_end =0
    n_correct_muon=0
    n_incorrect_muon =0 
    n_correct_antimuon=0
    n_incorrect_antimuon=0
    correct_percentage_total =0
    correct_percentage_total_antimuon =0
    no_reco = 0
    reco_ttree =0


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
        
       
        
        
        # use PositionTMSStart, MomentumTMSStart, and PositionTMSEnd, Muon_TrueKE ,truth level study
        for index, particle in enumerate(truth.PDG):
            if truth.PDG[index] == 13:
                n_true_muons+=1
                x_start = truth.BirthPosition[4*index+0]
                y_start = truth.BirthPosition[4*index+1]
                z_start = truth.BirthPosition[4*index+2]
                x_end = truth.DeathPosition[4*index+0]
                y_end = truth.DeathPosition[4*index+1]
                z_end = truth.DeathPosition[4*index+2]
                x_start_tms = truth.PositionTMSStart[4*index+0]
                y_start_tms = truth.PositionTMSStart[4*index+1]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                tracklength = math.sqrt((z_end-z_start_tms)**2 + (y_end-y_start_tms)**2 +(x_end-x_start_tms)**2)

                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    n_muon_total_lar_start_tms_end+=1
                    if len(reco.Direction) <= 0: 
                        no_reco+=1
                        break
                    if len(reco.Direction) > 0: reco_ttree +=1
                    p_z = reco.Direction[2]  
                    p_y = reco.Direction[1]  
                    p_x = reco.Direction[0]  
                    
                    if p_z <=0 : break
                    m= p_x / p_z 
                 
                    x_extrapolate =m*(reco.EndPos[2] - reco.StartPos[2]) + reco.StartPos[0]
                    
                    signed_dist = reco.EndPos[0] - x_extrapolate

                    KE = math.sqrt(p_x**2 + p_y**2 +p_z**2 +105.7**2) - 105.7
                    
                    
                    
                    if region1(x_start_tms) and region1(x_end):
                        
                        hist_signed_distance.Fill(signed_dist)
                        hist_total_charge_id.Fill(KE)
                        hist_tracklength_vs_starting_KE.Fill(KE,tracklength)
                        
                        if signed_dist > 0 :
                            hist_correct_charge_id.Fill(KE)
                            n_correct_muon+=1
                        else:  
                            hist_incorrect_charge_id.Fill(KE)
                            n_incorrect_muon +=1
                
                
                
               
                    if region2(x_start_tms) and region2(x_end):
                        
                        hist_signed_distance.Fill(-signed_dist)
                        hist_total_charge_id.Fill(KE)
                        hist_tracklength_vs_starting_KE.Fill(KE,tracklength)
                        if signed_dist < 0 :
                            hist_correct_charge_id.Fill(KE)
                            n_correct_muon+=1
                        else:  
                            hist_incorrect_charge_id.Fill(KE) 
                            n_incorrect_muon +=1
                
                
                
                
                    if region3(x_start_tms) and region3(x_end):
                        
                        hist_signed_distance.Fill(signed_dist)
                        hist_total_charge_id.Fill(KE)
                        hist_tracklength_vs_starting_KE.Fill(KE,tracklength)
                        if signed_dist > 0 :
                            hist_correct_charge_id.Fill(KE)
                            n_correct_muon+=1
                        else:  
                            hist_incorrect_charge_id.Fill(KE) 
                            n_incorrect_muon +=1
 
                           
            
            
            
            if truth.PDG[index] == -13:
                n_true_antimuons+=1
                x_start = truth.BirthPosition[4*index+0]
                y_start = truth.BirthPosition[4*index+1]
                z_start = truth.BirthPosition[4*index+2]
                x_end = truth.DeathPosition[4*index+0]
                y_end = truth.DeathPosition[4*index+1]
                z_end = truth.DeathPosition[4*index+2]
                x_start_tms = truth.PositionTMSStart[4*index+0]
                y_start_tms = truth.PositionTMSStart[4*index+1]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                tracklength = math.sqrt((z_end-z_start_tms)**2 + (y_end-y_start_tms)**2 +(x_end-x_start_tms)**2)

                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    n_antimuon_total_lar_start_tms_end+=1
                    if len(reco.Direction) <= 0:  
                        no_reco+=1
                        break
                    if len(reco.Direction) > 0: reco_ttree +=1
                    p_z = reco.Direction[2]  
                    p_y = reco.Direction[1]  
                    p_x = reco.Direction[0]  
                    
                    if p_z <=0 : break
                    m= p_x / p_z 
                 
                    x_extrapolate =m*(reco.EndPos[2] - reco.StartPos[2]) + reco.StartPos[0]
                    
                    signed_dist = reco.EndPos[0] - x_extrapolate

                    KE = math.sqrt(p_x**2 + p_y**2 +p_z**2 +105.7**2) - 105.7
                    if region1(x_start_tms) and region1(x_end):
                       
                        hist_signed_distance_antimuon.Fill(signed_dist)
                        hist_total_charge_id_antimuon.Fill(KE)
                        
                        if signed_dist < 0 :
                            hist_correct_charge_id_antimuon.Fill(KE)
                            n_correct_antimuon+=1
                        else:  
                            hist_incorrect_charge_id_antimuon.Fill(KE)
                            n_incorrect_antimuon +=1
                
                
                
               
                    if region2(x_start_tms) and region2(x_end):
                       
                        hist_signed_distance_antimuon.Fill(-signed_dist)
                        hist_total_charge_id_antimuon.Fill(KE)
                        
                        if signed_dist > 0 :
                            hist_correct_charge_id_antimuon.Fill(KE)
                            n_correct_antimuon+=1
                        else:  
                            hist_incorrect_charge_id_antimuon.Fill(KE) 
                            n_incorrect_antimuon +=1
                
                
                
                
                    if region3(x_start_tms) and region3(x_end):
                        
                        hist_signed_distance_antimuon.Fill(signed_dist)
                        hist_total_charge_id_antimuon.Fill(KE)
                       
                        if signed_dist < 0 :
                            hist_correct_charge_id_antimuon.Fill(KE)
                            n_correct_antimuon+=1
                        else:  
                            hist_incorrect_charge_id_antimuon.Fill(KE) 
                            n_incorrect_antimuon +=1
 
            
            
   
    hist_correct_charge_id_percentage.Divide(hist_correct_charge_id,hist_total_charge_id)
    hist_correct_charge_id_percentage_antimuon.Divide(hist_correct_charge_id_antimuon,hist_total_charge_id_antimuon) 
    correct_percentage_total = n_correct_muon/(n_correct_muon+n_incorrect_muon)
    correct_percentage_total_antimuon= n_correct_antimuon/(n_correct_antimuon+n_incorrect_antimuon)
   
    #n_region1_contained_percentage=   (n_region1_total-n_region1_not_contained)/n_region1_total   
    #n_region2_contained_percentage=   (n_region2_total-n_region2_not_contained)/n_region2_total        
    #n_region3_contained_percentage=   (n_region3_total-n_region3_not_contained)/n_region3_total
    #correct_percentage=n_correct_muon/(n_region1_total+n_region2_total+n_region3_total)    
       


   
    
    
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
    print(f"total correct percentage for muons: {correct_percentage_total}")
    print(f"total correct percentage for antimuons: {correct_percentage_total_antimuon}")
    print(f"correct signed distance muons: {n_correct_muon}")
    print(f"incorrect signed distance muons: {n_incorrect_muon}")
    print(f"correct signed distance antimuons: {n_correct_antimuon}")
    print(f"incorrect signed distance antimuons: {n_incorrect_antimuon}")
    print(f"N lar start tms end muons: {n_muon_total_lar_start_tms_end}")
    print(f"N lar start tms end antimuons: {n_antimuon_total_lar_start_tms_end}")
    print(f"empty reco ttree : {no_reco}")
    print(f"reco ttree : {reco_ttree}")
    


    
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
    


