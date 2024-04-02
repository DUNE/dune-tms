import glob
import argparse
import os
import math
import ROOT
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

def region2(x):
    is_region2= True
    if not -1500 < x < 1500: is_region2= False
    return is_region2

def region3(x):
    is_region3= True
    if not 2500 < x < 4000: is_region3= False
    return is_region3


def run(c, truth, outfilename, nmax=-1):
    """ This code does 3 things:
    Makes histograms
    Loops over all events and fills histograms
    Saves histograms in outfilename
    """
    
    # First make some histograms
    # You can add any histograms you want. Just make sure that histogram name hasn't been used before.
    
    # Hists related to track length
    
    bin_edges = array.array('d', [0, 200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,3000,4000,5000,6000,7000,8000,9000])
    
    if truth != None:
        # KE estimators
        hist_correct_charge_id = ROOT.TH1D("hist_correct_charge_id", "hist_correct_charge_id;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_incorrect_charge_id = ROOT.TH1D("hist_incorrect_charge_id", "hist_incorrect_charge_id;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_total_charge_id = ROOT.TH1D("hist_total_charge_id", "hist_total_charge_id;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_correct_charge_id_percentage = ROOT.TH1D("hist_correct_charge_id_percentage", "hist_correct_charge_id_percentage;True muon KE (MeV);percentage", 100, 0, 5000)
        hist_signed_distance = ROOT.TH1D("hist_signed_distance", "hist_signed_distance;signed distance(mm) ;Number of events", 100, -200000, 1000000)
        hist_correct_charge_id_antimuon = ROOT.TH1D("hist_correct_charge_id_antimuon", "hist_correct_charge_id_antimuon;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_incorrect_charge_id_antimuon = ROOT.TH1D("hist_incorrect_charge_id_antimuon", "hist_incorrect_charge_id_antimuon;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_total_charge_id_antimuon = ROOT.TH1D("hist_total_charge_id_antimuon", "hist_total_charge_id_antimuon;True muon KE (MeV);Number of events", 100, 0, 5000)
        hist_correct_charge_id_percentage_antimuon = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon", "hist_correct_charge_id_percentage_antimuon;True muon KE (MeV);percentage", 100, 0, 5000)
        hist_signed_distance_antimuon = ROOT.TH1D("hist_signed_distance_antimuon", "hist_signed_distance_antimuon ;signed distance(mm); Number of events", 100, -2000, 2000)
        # Can also calculate efficiency
       
      
        
        
        
        # We can also do some simple counts
        n_start_inside_tms = 0 
        n_end_inside_tms = 0
        n_start_and_end_inside_tms = 0
        n_true_muons = 0
        
    
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
    
    n_region1_total=0
    n_region1_not_contained=0
    n_region2_total=0
    n_region2_not_contained=0
    n_region3_total=0
    n_region3_not_contained=0
    n_region1_contained_percentage=0
    n_region2_contained_percentage=0
    n_region3_contained_percentage=0
    n_total=0
    n_correct=0
    correct_percentage=0



    # Now loop over all events
    for i, event in enumerate(c):
        if i > nevents: break
        if truth != None: truth.GetEntry(i)
        # Print current progress, with carriage return \r to prevent long list of progress and have everything on a singe line.
        if i % print_every == 0 and carriage: print(f"\rOn {i} / {nevents}", end='')
        if i % print_every == 0 and not carriage: print(f"On {i} / {nevents}")
        
        # Get some basic info
        nLines = event.nLines
        ntracks = nLines
        nhits = event.nHits
        hitpos = event.TrackHitPos
        
        # Loop over all tracks
        # Find the track with highest areal density (event.TrackLength)
        # This will be our primary muon candidate, longtrack
        # Also find the track length by first and last hough hit positions
        # This is not as good of an estimator because it's not taking into
        # account the various amounts of material in the path.
        # It's also not taking into account the curvature due to magnetic field
        # This code also plots all track lengths in a histogram
        # And saves them in an array for later use
        longtrack_length = 0
        longtrack = -1
        track_lengths = [0] * ntracks
        areal_densities = [0] * ntracks
        for track in range(ntracks):
           
            
            TrackLength = event.TrackLength[track]
            # Check if this is the longest track by areal density
            if TrackLength > longtrack_length:
                longtrack = track
                longtrack_length = TrackLength
            
            # Plot them in histograms
            
            
            # Save them in these arrays
        
        if longtrack >= 0:
            
        
        
        # If there's a muon candidate, plot its occupancy, and start and stop positions
        
           
            
            track_start_z = event.FirstHoughHit[2*longtrack+0]
            track_end_z = event.LastHoughHit[2*longtrack+0]
            track_start_x = event.FirstHoughHit[2*longtrack+1]
            track_end_x = event.LastHoughHit[2*longtrack+1]
            if truth != None:
             has_true_muon = truth.Muon_TrueKE 
             Muon_TrueKE = truth.Muon_TrueKE
            #if has_true_muon: print(Muon_TrueKE, has_true_muon)
            
            # Only fill these hists if there is a true muon
             if has_true_muon: 
                n_true_muons+=1
                
                # Now check if the muon started and ended in the TMS
                mx = truth.Muon_Vertex[0]
                my = truth.Muon_Vertex[1]
                mz = truth.Muon_Vertex[2]
                mdx = truth.Muon_Death[0]
                mdy = truth.Muon_Death[1]
                mdz = truth.Muon_Death[2]
                
                

               
                end_inside_tms = inside_tms(mdx, mdy, mdz)
                start_inside_lar = inside_lar(mx, my, mz)
                if start_inside_lar and end_inside_tms and truth.LeptonPDG==13 :
                    if region1(track_start_x) and region1(mdx):
                        n_region1_total+=1
                        n_total+=1
                        
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z
                        
                        signed_dist=-(k*track_end_z+b-track_end_x )*abs(k*track_end_z+b-track_end_x )
                        for i in range(0, len(event.TrackHitPos), 2):
                            if event.TrackHitPos[i]==-999: continue
                            signed_dist=signed_dist-(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])*abs(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])
                        hist_signed_distance.Fill(signed_dist)
                        hist_total_charge_id.Fill(Muon_TrueKE)
                        if signed_dist > 0 :
                            hist_correct_charge_id.Fill(Muon_TrueKE)
                            n_correct+=1
                        else:  
                            hist_incorrect_charge_id.Fill(Muon_TrueKE) 
                        for i in range(1, len(hitpos), 2):
                          if not region1(hitpos[i]): 
                              n_region1_not_contained+=1
                              break


                    if region2(track_start_x) and region2(mdx):
                        n_total+=1
                        n_region2_total+=1
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z                  
                           
                        signed_dist=(k*track_end_z+b-track_end_x )*abs(k*track_end_z+b-track_end_x )
                        for i in range(0, len(event.TrackHitPos), 2):
                            if event.TrackHitPos[i]==-999: continue
                            signed_dist=signed_dist+(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])*abs(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])
                        hist_signed_distance.Fill(signed_dist)
                        hist_total_charge_id.Fill(Muon_TrueKE)
                        if signed_dist > 0 :
                            hist_correct_charge_id.Fill(Muon_TrueKE)
                            n_correct+=1
                        else:  
                            hist_incorrect_charge_id.Fill(Muon_TrueKE)  
                        for i in range(1, len(hitpos), 2):
                          if  not region2(hitpos[i]): 
                              n_region2_not_contained+=1
                              break   
                
                    if region3(track_start_x) and region3(mdx):
                        n_total+=1
                        n_region3_total+=1
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z
                        signed_dist=-(k*track_end_z+b-track_end_x )*abs(k*track_end_z+b-track_end_x )
                        for i in range(0, len(event.TrackHitPos), 2):
                            if event.TrackHitPos[i]==-999: continue
                            signed_dist=signed_dist-(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])*abs(k*event.TrackHitPos[i]+b-event.TrackHitPos[i+1])
                        hist_signed_distance.Fill(signed_dist)
                        hist_total_charge_id.Fill(Muon_TrueKE)
                        if signed_dist > 0 :
                            hist_correct_charge_id.Fill(Muon_TrueKE)
                            n_correct+=1
                        else:  
                            hist_incorrect_charge_id.Fill(Muon_TrueKE) 
                        for i in range(1, len(hitpos), 2):
                          if  not region3(hitpos[i]): 
                              n_region3_not_contained+=1
                              break
                
                
                if start_inside_lar and end_inside_tms and truth.LeptonPDG==-13 :
                    if region1(track_start_x) and region1(mdx):
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z
                        signed_dist=-(k*track_end_z+b-mdx )
                        hist_signed_distance_antimuon.Fill(signed_dist)
                        hist_total_charge_id_antimuon.Fill(Muon_TrueKE)
                        if signed_dist < 0 :
                            hist_correct_charge_id_antimuon.Fill(Muon_TrueKE)
                        else:
                            hist_incorrect_charge_id_antimuon.Fill(Muon_TrueKE)    

                    if region2(track_start_x) and region2(mdx):
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z
                        signed_dist=(k*track_end_z+b-mdx )
                        hist_signed_distance_antimuon.Fill(signed_dist)
                        hist_total_charge_id_antimuon.Fill(Muon_TrueKE)
                        if signed_dist < 0 :
                            hist_correct_charge_id_antimuon.Fill(Muon_TrueKE)
                        else:
                            hist_incorrect_charge_id_antimuon.Fill(Muon_TrueKE)    
    
                
                    if region3(track_start_x) and region3(mdx):
                        k=(mx-track_start_x)/(mz-track_start_z)
                        b=track_start_x-k*track_start_z
                        signed_dist=-(k*track_end_z+b-mdx )
                        hist_signed_distance_antimuon.Fill(signed_dist)
                        hist_total_charge_id_antimuon.Fill(Muon_TrueKE)
                        if signed_dist < 0 :
                            hist_correct_charge_id_antimuon.Fill(Muon_TrueKE)  
                        else:
                            hist_incorrect_charge_id_antimuon.Fill(Muon_TrueKE)    

                
                
                          


   
    hist_correct_charge_id_percentage.Divide(hist_correct_charge_id,hist_total_charge_id)
    hist_correct_charge_id_percentage_antimuon.Divide(hist_correct_charge_id_antimuon,hist_total_charge_id_antimuon)
   
    # Before the division operation
    if n_region1_total != 0:
        n_region1_contained_percentage = (n_region1_total - n_region1_not_contained) / n_region1_total
    else:
    # Handle the case where n_region1_total is zero
        n_region1_contained_percentage = 0  
  
        n_region2_contained_percentage=   (n_region2_total-n_region2_not_contained)/n_region2_total        
        n_region3_contained_percentage=   (n_region3_total-n_region3_not_contained)/n_region3_total  
        correct_percentage=n_correct/n_total
        # Now check if the muon is a good candidate
        # Does this by checking a list of conditions. If any of them are false, don't fill the reco vs true ke plots.
       


   
    
    
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
    
    if truth != None:
        print(f"N true muons: {n_true_muons}")
        print(f"n_region1_total: {n_region1_total}")
        print(f"n_region2_total: {n_region2_total}")
        print(f"n_region3_total: {n_region3_total}")
        print(f"n_region1_not_contained: {n_region1_not_contained}")
        print(f"n_region2_not_contained: {n_region2_not_contained}")
        print(f"n_region3_not_contained: {n_region3_not_contained}")
        print(f"n_region1_contained_percentage: {n_region1_contained_percentage}")
        print(f"n_region2_contained_percentage: {n_region2_contained_percentage}")
        print(f"n_region3_contained_percentage: {n_region3_contained_percentage}")
        print(f"correct_percentage: {correct_percentage}")

    
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
        # In this case, the user specified exactly one file. Usually they'd hadd many files together. # Return this hists if the user requested previews        
        files_to_use = [infile]
        
    outdir = args.outdir
    if outdir == "":
        # No output directory was specified so use the default
        # First we need the username
        username = os.environ["USER"]
        outdir = f"/exp/dune/app/users/{username}/dune-tms_hists2"
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
    
    # Make the TChain objects. One for truth information and one for reconstructed information.
    c = ROOT.TChain("Line_Candidates")
    truth = None
    if has_truth: truth = ROOT.TChain("Truth_Info")
    for f in files_to_use:
        c.Add(f)
        if has_truth: truth.Add(f)
    assert c.GetEntries() > 0, "Didn't get any entries in Line_Candidates TChain." 
    if has_truth: assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info TChain."
    
    nevents = args.nevents
    assert nevents >= -1, f"nevents <= -1, why? {nevents}"
    
    # Now finally run
    hists = run(c, truth, outfilename, nevents)
    
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
