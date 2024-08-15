import glob
import argparse
import os
import math
import ROOT # type: ignore
import array
import pickle
# from mat_helper import Hists_Graph
# import numpy as np

# Tells root to be in batch mode so it doesn't try to create a canvas on your screen, which is slow
ROOT.gROOT.SetBatch(True)
# Don't draw the stats box on histograms
ROOT.gStyle.SetOptStat(0)
# Python's wrapper around ROOT sometimes crashes when it automatically adds hists to its garbage collector. This fixed that
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

FUDICIAL_CUT = 50
LAR_START = (-3478.48,-2166.71,4179.24)
LAR_END =(3478.48,829.282,9135.88)

def in_between(x_lar_start,x_lar_end,y_lar_start,y_lar_end,z_lar_start,z_lar_end): 
    if LAR_START[0]< x_lar_start and LAR_START[1]<y_lar_start and LAR_START[2]<z_lar_start and LAR_END[0]>x_lar_end and LAR_END[1]>y_lar_end and LAR_END[2]>z_lar_end: 
        return True
    else:
        return False


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


def run(c, truth, f,outfilename, nmax=-1):
    """ This code does 3 things:
    Makes histograms
    Loops over all events and fills histograms
    Saves histograms in outfilename
    """
    # use PositionTMSStart, MomentumTMSStart, and PositionTMSEnd


    bin_edges = array.array('d', [0, 200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,3000,4000,5000,6000,7000,8000,9000])

    if truth != None:
        # Define histograms for muons
        # hist_signed_distance = Hists_Graph("hist_signed_distance", "Muons signed distance: (x_extropolate - x_truth) (using truth_info);True signed distance(mm) ;Number of muons", 100, -2000, 2000)
        pass








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
    n_muon_total_lar_start_tms_end =0
    n_correct=0
    n_true_amuons = 0
    n_acorrect= 0
    correct_percentage_total =0
    muon_signed_distances = []
    amuon_signed_distances = []
    correct_muon_ke = []
    total_muon_ke = []
    correct_amuon_ke = []
    total_amuon_ke = []


    # Now loop over all events
    for i, event in enumerate(c):
        if i > nevents: break
        #if i>100000: break
        if truth != None: truth.GetEntry(i)
        # Print current progress, with carriage return \r to prevent long list of progress and have everything on a singe line.
        if i % print_every == 0 and carriage: print(f"\rOn {i} / {nevents}  {i/nevents*100}%", end='')
        if i % print_every == 0 and not carriage: print(f"On {i} / {nevents}   {i/nevents*100:.1f}%")

        
        # use PositionTMSStart, MomentumTMSStart, and PositionTMSEnd, truth level study
        for index, particle in enumerate(truth.PDG):
            print(i)
            signed_dist=None
            if truth.PDG[index] == 13:
                n_true_muons+=1
                x_start = truth.BirthPosition[4*index+0]
                y_start = truth.BirthPosition[4*index+1]
                z_start = truth.BirthPosition[4*index+2]
                x_end = truth.DeathPosition[4*index+0]
                y_end = truth.DeathPosition[4*index+1]
                z_end = truth.DeathPosition[4*index+2]
                x_start_tms = truth.PositionTMSStart[4*index+0]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                KE_muon = truth.Muon_TrueKE
                x_lar_start = truth.PositionLArStart[4*index+0]+FUDICIAL_CUT
                y_lar_start = truth.PositionLArStart[4*index+1]+FUDICIAL_CUT
                z_lar_start = truth.PositionLArStart[4*index+2]+FUDICIAL_CUT

                x_lar_end = truth.PositionLArEnd[4*index+0]-FUDICIAL_CUT
                y_lar_end = truth.PositionLArEnd[4*index+1]-FUDICIAL_CUT
                z_lar_end = truth.PositionLArEnd[4*index+2]-FUDICIAL_CUT

                # p = Momentum(KE_muon,classification="muon")
                
            
                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    n_muon_total_lar_start_tms_end+=1
                    if region1(x_start_tms) and region1(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]
                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = x_end - x_extrapolate
                        total_muon_ke.append(truth.Muon_TrueKE)
                        if signed_dist > 0 :
                            correct_muon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1
                            # n_correct+=1


                    if region2(x_start_tms) and region2(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]
                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = -(x_end - x_extrapolate)
                        total_muon_ke.append(truth.Muon_TrueKE)
                        if signed_dist > 0 :
                            correct_muon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1

                        
                    if region3(x_start_tms) and region3(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]
                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = x_end - x_extrapolate
                        total_muon_ke.append(truth.Muon_TrueKE)
                        if signed_dist > 0 :
                            correct_muon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1
                        
                    if signed_dist!=None:
                        muon_signed_distances.append(signed_dist)
                        # print(f"{n_correct} / {n_true_muons}")
                        
                    
                    
            if truth.PDG[index] == -13:
                n_true_amuons+=1
                x_start = truth.BirthPosition[4*index+0]
                y_start = truth.BirthPosition[4*index+1]
                z_start = truth.BirthPosition[4*index+2]
                x_end = truth.DeathPosition[4*index+0]
                y_end = truth.DeathPosition[4*index+1]
                z_end = truth.DeathPosition[4*index+2]
                x_start_tms = truth.PositionTMSStart[4*index+0]
                z_start_tms = truth.PositionTMSStart[4*index+2]
                KE_muon = truth.Muon_TrueKE
                # p = Momentum(KE_muon,classification="amuon")
                x_lar_start = truth.PositionLArStart[4*index+0]+FUDICIAL_CUT
                y_lar_start = truth.PositionLArStart[4*index+1]+FUDICIAL_CUT
                z_lar_start = truth.PositionLArStart[4*index+2]+FUDICIAL_CUT

                x_lar_end = truth.PositionLArEnd[4*index+0]-FUDICIAL_CUT
                y_lar_end = truth.PositionLArEnd[4*index+1]-FUDICIAL_CUT
                z_lar_end = truth.PositionLArEnd[4*index+2]-FUDICIAL_CUT

                

                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    if region1(x_start_tms) and region1(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]
                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = x_end - x_extrapolate
                        total_amuon_ke.append(truth.Muon_TrueKE)
                        if signed_dist > 0 :
                            correct_amuon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1
                        
                    if region2(x_start_tms) and region2(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]

                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = -(x_end - x_extrapolate)
                        total_amuon_ke.append(truth.Muon_TrueKE)
                        if signed_dist < 0 :
                            correct_amuon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1
                        



                    if region3(x_start_tms) and region3(x_end):
                        p_z = truth.MomentumTMSStart[4*index+2]
                        p_x = truth.MomentumTMSStart[4*index+0]
                        if p_z ==0: break
                        m= p_x / p_z
                        b= x_start_tms-m*z_start_tms
                        x_extrapolate =m*z_end +b
                        signed_dist = x_end - x_extrapolate
                        total_amuon_ke.append(truth.Muon_TrueKE)
                        if signed_dist > 0 :
                            correct_amuon_ke.append(truth.Muon_TrueKE)
                            n_correct+=1
                        
                    if signed_dist!=None:
                        amuon_signed_distances.append(signed_dist)
                            # print("Amuon Identified")
                            # print(n_acorrect/n_true_amuons)
                            # amuon_correct_charge[1].append(n_acorrect/n_true_amuons)
                            # amuon_correct_charge[0].append(KE_muon)
                            
                            


                    
            
    correct_percentage_total = n_correct/n_muon_total_lar_start_tms_end
    #n_region1_contained_percentage=   (n_region1_total-n_region1_not_contained)/n_region1_total
    #n_region2_contained_percentage=   (n_region2_total-n_region2_not_contained)/n_region2_total
    #n_region3_contained_percentage=   (n_region3_total-n_region3_not_contained)/n_region3_total
    #correct_percentage=n_correct/(n_region1_total+n_region2_total+n_region3_total)






    ## Now save all the histograms
    if carriage: print("\r", end="") # Erase the current line
    print(f"Completed loop of {nevents} events.From Filename : {f}.Now, Saving to {outfilename}")
    


    # Return this hists if the user requested previews
    return muon_signed_distances,amuon_signed_distances,[correct_muon_ke,total_muon_ke],[correct_amuon_ke,total_amuon_ke]


def validate_then_run(args):
    ''' This function has all the code to validate the arguments. Usually users don't need to edit this code. '''
    # First we validate the args
    infile = args.filename
    if infile != "":
        # In this case, the user specified exactly one file. Usually they'd hadd many files together.
        files_to_use = infile

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

    # muon_graph = Hists_Graph("hist_signed_distance", "Muons signed distance: (x_extropolate - x_truth) (using truth_info);True signed distance(mm) ;Number of muons")
    # amuon_graph = Hists_Graph("hist_signed_distance_antimuon", "Antimuons signed distance: (x_extropolate - x_truth) (using truth_info) ;True signed distance(mm); Number of antimuons")
    entire_array = {}
    entire_array["MUONS"] = {}
    entire_array["AMUONS"] = {}
    entire_array["MUONS_CORRECT"]={}
    entire_array["AMUONS_CORRECT"]={}
    for f in files_to_use:
    # Make the TChain objects. One for truth information and one for reconstructed information.
        c = ROOT.TChain("Line_Candidates")
        truth = None
        if has_truth: truth = ROOT.TChain("Truth_Info")
    # for f in files_to_use:
        c.Add(f)
        if has_truth: truth.Add(f)
        assert c.GetEntries() > 0, "Didn't get any entries in Line_Candidates TChain."
        if has_truth: assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info TChain."

        nevents = args.nevents
        assert nevents >= -1, f"nevents <= -1, why? {nevents}"

        # Now finally run
        muon_arr,amuon_arr,muon_correct_charge,amuon_correct_charge = run(c, truth,f, outfilename, nevents)
        entire_array["MUONS"][f] = muon_arr
        entire_array["AMUONS"][f] = amuon_arr
        entire_array["MUONS_CORRECT"][f] = muon_correct_charge
        entire_array["AMUONS_CORRECT"][f] = amuon_correct_charge


        # muon_graph.add(muon_arr)
        # amuon_graph.add(amuon_arr)
        
    # muon_graph.finish(f)
    # amuon_graph.finish(f)
    with open("python_object.sushil_dai", 'wb') as file:
        pickle.dump(entire_array, file)

    # preview = args.preview
    # if preview:
    #     print("Making preview histograms.")
    #     print("WARNING: Eventually you'll want to make plots using a separate plotting script")
    #     print("That way you can format the plots without having to rerun the make_hists.py script")
    #     # Make the output directory if needed
    #     outname_without_root = outname.replace(".root", "")
    #     preview_dir = os.path.join(outdir, "previews", outname_without_root)
    #     if not os.path.exists(preview_dir):
    #         os.makedirs(preview_dir)
    #     print(f"Saving previews in {preview_dir}/<hist_name>.png")
    #     # canvas = ROOT.TCanvas()
    #     # for hist in hists:
    #         # name = hist.GetName().replace("hist_", "")
    #         # opts = "colz"
    #         # hist.Draw(opts)
    #         canvas.Print(os.path.join(preview_dir, f"{name}.png"))








# This is how python knows whether it's a script that's being imported into another script or running as the main script
if __name__ == "__main__":
    def list_of_strings(arg):
        return arg.split(',')

    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', type=str, help="The output dir. Will be made if it doesn't exist.", default="")
    parser.add_argument('--name', type=str, help="The name of the output files.", default="dune-tms_hists.root")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--indir', type=str, help="The location of the input tmsreco files", default="")
    group.add_argument('--inlist', type=str, help="The input filelist", default="")
    group.add_argument('--filename', '-f', type=list_of_strings, help="The input file, if you have a single file", default="")
    parser.add_argument('--nevents', '-n', type=int, help="The maximum number of events to loop over", default=-1)
    parser.add_argument('--allow_overwrite', help="Allow the output file to overwrite", action=argparse.BooleanOptionalAction)
    parser.add_argument('--preview', help="Save preview images of the histograms", action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    validate_then_run(args)
