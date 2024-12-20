import glob
import argparse
import os
import math
import ROOT # type: ignore
import array
import dunestyle.root as dunestyle
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
    if not -3300 < x < 3300: is_inside = False
    if not -2850 < y < 160: is_inside = False
    if only_thin_section:
        if not 11362 < z < 13500: is_inside = False
    else:
        if not 11362 < z < 18314: is_inside = False
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
    if not -3300 < x < -1750: is_region1= False
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
    if not 1750 < x < 3300: is_region3= False
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

        hist_muon_signed_distance = ROOT.TH1D("hist_muon_signed_distance", "Muons' signed distance (using reco, positive number means muon-like);signed distance (mm);Number of muons", 100, -15000, 15000)
        hist_antimuon_signed_distance = ROOT.TH1D("hist_antimuon_signed_distance", "Antimuons' signed distance (using reco, negative number means antimuon-like);signed distance (mm);Number of antimuons", 100, -15000, 15000)
        hist_muon_signed_distance_0_to_300 = ROOT.TH1D("hist_muon_signed_distance_0_to_300", "Muons' signed distance (0 to 300MeV)(using reco, positive number means muon-like);signed distance (mm);Number of muons", 100,  -15000, 15000)
        hist_antimuon_signed_distance_0_to_300 = ROOT.TH1D("hist_antimuon_signed_distance_0_to_300", "Antimuons' signed distance(0 to 300MeV) (using reco, negative number means antimuon-like);signed distance (mm);Number of antimuons", 100,  -15000, 15000)
        hist_muon_signed_distance_300_to_1000 = ROOT.TH1D("hist_muon_signed_distance_300_to_1000", "Muons' signed distance (300 to 1000MeV)(using reco, positive number means muon-like);signed distance (mm);Number of muons", 100, -15000, 15000)
        hist_antimuon_signed_distance_300_to_1000 = ROOT.TH1D("hist_antimuon_signed_distance_300_to_1000", "Antimuons' signed distance(300 to 1000MeV) (using reco, negative number means antimuon-like);signed distance (mm);Number of antimuons", 100,  -15000, 15000)
        hist_muon_signed_distance_1000_to_4000 = ROOT.TH1D("hist_muon_signed_distance_1000_to_4000", "Muons' signed distance (1000 to 4000MeV)(using reco, positive number means muon-like);signed distance (mm);Number of muons", 100,  -15000, 15000)
        hist_antimuon_signed_distance_1000_to_4000 = ROOT.TH1D("hist_antimuon_signed_distance_1000_to_4000", "Antimuons' signed distance(1000 to 4000MeV) (using reco, negative number means antimuon-like);signed distance (mm);Number of antimuons", 100,  -15000, 15000)
        hist_muon_signed_distance_4000_to_5000 = ROOT.TH1D("hist_muon_signed_distance_4000_to_5000", "Muons' signed distance (4000 to 5000MeV)(using reco, positive number means muon-like);signed distance (mm);Number of muons", 100,  -45000, 45000)
        hist_antimuon_signed_distance_4000_to_5000 = ROOT.TH1D("hist_antimuon_signed_distance_4000_to_5000", "Antimuons' signed distance(4000 to 5000MeV) (using reco, negative number means antimuon-like);signed distance (mm);Number of antimuons", 100,  -45000, 45000)
        

        hist_kemu_vs_thetamu_correct_chargeID_fraction_muon = ROOT.TH2D("hist_kemu_vs_thetamu_correct_chargeID_fraction_muon", "Muon correct chargeID fraction (RHC); TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_correct_chargeID_muon = ROOT.TH2D("hist_kemu_vs_thetamu_correct_chargeID_muon", "Muon correct chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_incorrect_chargeID_muon = ROOT.TH2D("hist_kemu_vs_thetamu_incorrect_chargeID_muon", "Muon incorrect chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_total_chargeID_muon = ROOT.TH2D("hist_kemu_vs_thetamu_total_chargeID_muon", "Muon total chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        
        
        
        hist_kemu_vs_thetamu_correct_chargeID_fraction_antimuon = ROOT.TH2D("hist_kemu_vs_thetamu_correct_chargeID_fraction_antimuon", "Antimuon correct chargeID fraction (RHC); TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_correct_chargeID_antimuon = ROOT.TH2D("hist_kemu_vs_thetamu_correct_chargeID_antimuon", "Antimuon correct chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_incorrect_chargeID_antimuon = ROOT.TH2D("hist_kemu_vs_thetamu_incorrect_chargeID_antimuon", "Antimuon incorrect chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        hist_kemu_vs_thetamu_total_chargeID_antimuon = ROOT.TH2D("hist_kemu_vs_thetamu_total_chargeID_antimuon", "Antimuon total chargeID, TMS Entry Energy vs entrance angle; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 65)
        

        hist_kemu_vs_thetamu_correct_chargeID_comparison_muon_antimuon = ROOT.TH2D("hist_kemu_vs_thetamu_correct_chargeID_comparison_muon_antimuon", "Muon correct chargeID over antimuon correct chargeID; TMS Entry KE (MeV);TMS Entry Angle (in x-z plane)", 50, 0, 5000, 50, -50, 50)
        hist_correct_charge_id_percentage_muon_over_antimuon = ROOT.TH1D("hist_correct_charge_id_percentage_muon_over_antimuon", "Muons correct charge percentage over antimuon (using reco);True muon TMS Entry KE (MeV);Fraction", 25, 0, 5000)



        hist_hits_vs_energy = ROOT.TH2D("hist_hits_vs_energy", "Muon hit number vs Energy;TMS Entry Energy (MeV);Hits leave in TMS", 50, 0, 5000, 50, 0, 100)
        hist_hits_vs_energy_antimuon = ROOT.TH2D("hist_hits_vs_energy_antimuon", "Antimuon hit number vs TMS Entry Energy;Energy (MeV);Hits leave in TMS", 50, 0, 5000, 50, 0, 100)
        
        
        
        
        
        
        hist_unidentifiable_muon = ROOT.TH1D("hist_unidentifiable_muon", "unidentifiable muons hits number (using reco);Hits leave in TMS;Muon number", 100, 0, 100)
        hist_unidentifiable_muon_less_than_20_hits = ROOT.TH1D("hist_unidentifiable_muon_less_than_20_hits", "unidentifiable muons hits number(less than 20 hits) (using reco);Hits leave in TMS;Muon number", 21, 0, 20)
        
        hist_correct_charge_id_with_hits = ROOT.TH1D("hist_correct_charge_id_with_hits", "correct chargeID muons hits number (using reco);Hits leave in TMS;Muon number", 100, 0, 100)
        hist_correct_charge_id_with_hits_less_than_20 = ROOT.TH1D("hist_correct_charge_id_with_hits_less_than_20", "correct chargeID muons hits number(less than 20) (using reco);Hits leave in TMS;Muon number", 20, 0, 20)
        
        hist_incorrect_charge_id_with_hits = ROOT.TH1D("hist_incorrect_charge_id_with_hits", "incorrect chargeID muons hits number (using reco);Hits leave in TMS;Muon number", 100, 0, 100)
        hist_incorrect_charge_id_with_hits_less_than_20 = ROOT.TH1D("hist_incorrect_charge_id_with_hits_less_than_20", "incorrect chargeID muons hits number(less than 20) (using reco);Hits leave in TMS;Muon number", 20, 0, 20)
        
        hist_muon_hits_in_TMS = ROOT.TH1D("hist_muon_hits_in_TMS", " muons hits number (using reco);Hits leave in TMS;Muon number", 100, 0, 100)
        hist_reconstructed_muons = ROOT.TH1D("hist_reconstructed_muons", "reconstructed muons distribution with energy;True muon TMS Entry KE (MeV);Number of renconstructed muons", 25, 0, 5000)
        hist_unreconstructed_muons = ROOT.TH1D("hist_unreconstructed_muons", "unreconstructed muons distribution with energy;True muon TMS Entry KE (MeV);Number of unrenconstructed muons", 25, 0, 5000)
        hist_total_muons = ROOT.TH1D("hist_total_muons", "total muons distribution with energy(lar to tms);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
      
        
        hist_reconstructed_muons_percentage = ROOT.TH1D("hist_reconstructed_muons_percentage", "reconstructed muons percentage;True muon TMS Entry KE (MeV);Percentage of renconstructed muons", 25, 0, 5000)
   
        hist_endpoint_diff_truth_reco_in_z = ROOT.TH1D("hist_endpoint_diff_truth_reco_in_z", "muon endpoint difference between truth info and reco in z; Difference (mm); Number of muons", 100, -700, 700)
        hist_endpoint_diff_truth_reco_in_z_antimuon = ROOT.TH1D("hist_endpoint_diff_truth_reco_in_z_antimuon", "antimuon endpoint difference between truth info and reco in z; Difference (mm); Number of antimuons", 100, -700, 700)
        dunestyle.Simulation()
   
   
   
   
   
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
    region1_muon_number =0
    region2_muon_number =0
    region3_muon_number =0
    num_unidentifiable_muon = 0
    num_unidentifiable_antimuon = 0
    gross_correct_percentage_total =0
    gross_correct_percentage_total_antimuon =0
    n_default_number =0
    reco_muons = 0
    unreco_muons = 0
    event_display_num_unidentifiable =0
    event_display_num_unidentifiable_misidentified =0
    event_display_num = 0
    n_correct_muon_100_cut=0
    n_incorrect_muon_100_cut =0 
    n_correct_muon_300_cut=0
    n_incorrect_muon_300_cut =0 
    n_correct_muon_500_cut=0
    n_incorrect_muon_500_cut =0 
    n_correct_muon_700_cut=0
    n_incorrect_muon_700_cut =0 
    n_correct_antimuon_100_cut=0
    n_incorrect_antimuon_100_cut =0 
    n_correct_antimuon_300_cut=0
    n_incorrect_antimuon_300_cut =0 
    n_correct_antimuon_500_cut=0
    n_incorrect_antimuon_500_cut =0 
    n_correct_antimuon_700_cut=0
    n_incorrect_antimuon_700_cut =0 
    correct_percentage_total_100_cut = 0
    correct_percentage_total_300_cut = 0
    correct_percentage_total_500_cut = 0
    correct_percentage_total_700_cut = 0
    correct_percentage_total_100_cut_antimuon = 0
    correct_percentage_total_300_cut_antimuon = 0
    correct_percentage_total_500_cut_antimuon = 0
    correct_percentage_total_700_cut_antimuon = 0
    
    
    

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
            if truth.PDG[index] == 13 or truth.PDG[index] == -13:
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
                p_z = truth.MomentumTMSStart[4*index+2] 
                p_y = truth.MomentumTMSStart[4*index+1]
                p_x = truth.MomentumTMSStart[4*index+0]
                tracklength = math.sqrt((z_end-z_start_tms)**2 + (y_end-y_start_tms)**2 +(x_end-x_start_tms)**2)
                KE = math.sqrt(p_x**2 + p_y**2 +p_z**2 +105.7**2) - 105.7
                if p_z > 0:
                    thetamu = math.degrees(math.atan(p_x / p_z))
                
                
                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end) and p_z > 0:
                    chargeID = -999
                    muon_hits = []
                    region1_hits =[]
                    region2_hits =[]
                    region3_hits =[]
                    n_changeregion = 0
                    dist_plus_overweighted =0
                    dist_minus_overweighted =0
                    dist_plus_underweighted =0
                    dist_minus_underweighted =0
                    dist_plus =0
                    dist_minus =0
                    n_hit_total=0
                    n_hits=0
                    markers = []
                    z_reco =0
                    n_last_valid_hit = 0
                    z_diff = -999

                    if len(reco.TrackHitPos) < 3:               
                        break
                    
                    # This line can be ignored
                    n_muon_total_lar_start_tms_end+=1
                    
                    #This loop takes in all the muon hits and put them in a list that I named
                    #Check if the first hit is correctly filled up

                    if  reco.TrackHitPos[0] == 0 or reco.TrackHitPos[0] == -999:
                        for i in range(3,len(reco.TrackHitPos)):
                            muon_hits.append(reco.TrackHitPos[i])
                    else:
                        for i in range(len(reco.TrackHitPos)):
                            muon_hits.append(reco.TrackHitPos[i])
                    for i in range(len(muon_hits)):
                       if not ( muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0):  
                           n_last_valid_hit = (i//3)*3
                    z_reco = muon_hits[n_last_valid_hit + 2] 
                    z_diff = z_end -  z_reco
                    if truth.PDG[index] == 13 :
                        hist_endpoint_diff_truth_reco_in_z.Fill(z_diff)
                    if truth.PDG[index] == -13:
                        hist_endpoint_diff_truth_reco_in_z_antimuon.Fill(z_diff)

                    #If muon starts in region 1 
                    if region1(muon_hits[0]):
                        
                        for i in range(len(muon_hits)): 
                            # Check if the muon is still in the same region as the starting region. 
                            # Mark down the point where the muon crosses into a different region, 
                            # or if there is an invalid entry. 
                            if not region1(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0: 
                                n_changeregion = i
                                break
                            # Now we have all the hits in region 1 
                            region1_hits.append(muon_hits[i])
                            
                        if region2(muon_hits[n_changeregion]):
                            # Fill region 2 hits if there are any
                            for i in range(n_changeregion, len(muon_hits)): 
                                if not region2(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0 : break
                                region2_hits.append(muon_hits[i])
                    
                    #If muon starts in region 3 ,Similar to the above algorithm
                    if region3(muon_hits[0]):
                        for i in range(len(muon_hits)): 
                            if not region3(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0: 
                                n_changeregion = i
                                break
                            region3_hits.append(muon_hits[i])
                        if region2(muon_hits[n_changeregion]):
                            for i in range(n_changeregion, len(muon_hits)): 
                                if not region2(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0 : break
                                region2_hits.append(muon_hits[i])

                    
                    #If muon starts in region 2, Similar to the above algorithm
                    if region2(muon_hits[0]):
                       
                        for i in range(len(muon_hits)): 
                            if not region2(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0: 
                                n_changeregion = i
                                break
                            region2_hits.append(muon_hits[i])
                        if region1(muon_hits[n_changeregion]):
                            
                            for i in range(n_changeregion, len(muon_hits)): 
                                if not region1(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0 : break
                                region1_hits.append(muon_hits[i])
                                
                        if region3(muon_hits[n_changeregion]):
                            
                            for i in range(n_changeregion, len(muon_hits)): 
                                if not region3(muon_hits[(i//3)*3]) or muon_hits[(i//3)*3]==-999 or muon_hits[(i//3)*3]==0 : break
                                region3_hits.append(muon_hits[i])
                    
                    total_hit_region1 = len(region1_hits)//3
                    total_hit_region2 = len(region2_hits)//3
                    total_hit_region3 = len(region3_hits)//3

                    n_hit_total = total_hit_region1 + total_hit_region2 + total_hit_region3

                
                    
                    
                    #No need to do any calculation if fewer than four hits are available 
                    if (total_hit_region1 + total_hit_region2 + total_hit_region3) <4: break
                   
                    #Now the muon hits are collected in three different regions, do the calculation
                    
                    #Region 1 calculation 
                    if total_hit_region1>2 and region1_hits[2]- region1_hits[3*total_hit_region1-1] <0:
                        
                        m = (region1_hits[0] - region1_hits[3*total_hit_region1-3])/(region1_hits[2]-region1_hits[3*total_hit_region1-1] )  
                        for i in range(1,total_hit_region1-2):
                            x_interpolation =m*(region1_hits[3*i+2] - region1_hits[2]) + region1_hits[0]
                            signed_dist = region1_hits[3*i] - x_interpolation
                            if signed_dist > 0: 
                                dist_plus +=abs(signed_dist)
                                dist_plus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_plus_underweighted += math.sqrt(abs(signed_dist))

                            if signed_dist < 0: 
                                dist_minus +=abs(signed_dist)
                                dist_minus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_minus_underweighted += math.sqrt(abs(signed_dist))
                                

                    #Region 2 calculation 
                    if total_hit_region2>2 and region2_hits[2]-region2_hits[3*total_hit_region2-1] <0:
                       
                        m = (region2_hits[0] - region2_hits[3*total_hit_region2-3])/(region2_hits[2]-region2_hits[3*total_hit_region2-1] )  
                        for i in range(1,total_hit_region2-2):
                            x_interpolation =m*(region2_hits[3*i+2] - region2_hits[2]) + region2_hits[0]
                            signed_dist = region2_hits[3*i] - x_interpolation
                            if signed_dist < 0: 
                                dist_plus +=abs(signed_dist)
                                dist_plus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_plus_underweighted += math.sqrt(abs(signed_dist))
                                
                            if signed_dist > 0: 
                                dist_minus +=abs(signed_dist)
                                dist_minus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_minus_underweighted += math.sqrt(abs(signed_dist))
    
                    #Region 3 calculation 
                    if total_hit_region3>2 and region3_hits[2]-region3_hits[3*total_hit_region3-1] <0:
                        m = (region3_hits[0] - region3_hits[3*total_hit_region3-3])/(region3_hits[2]-region3_hits[3*total_hit_region3-1] )  
                        for i in range(1,total_hit_region3-2):
                            x_interpolation =m*(region3_hits[3*i+2] - region3_hits[2]) + region3_hits[0]
                            signed_dist = region3_hits[3*i] - x_interpolation
                            if signed_dist > 0: 
                                dist_plus +=abs(signed_dist)
                                dist_plus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_plus_underweighted += math.sqrt(abs(signed_dist))

                            if signed_dist < 0: 
                                dist_minus +=abs(signed_dist)
                                dist_minus_overweighted += abs(signed_dist)*abs(signed_dist)
                                dist_minus_underweighted += math.sqrt(abs(signed_dist))
                    
                    
                    #Now judge whether the particle is a muon or an antimuon 
                   
                    
                    #normal
                    if truth.PDG[index] == 13 :
                        hist_hits_vs_energy.Fill(KE, n_hit_total)
                        hist_muon_signed_distance.Fill(dist_minus-dist_plus)
                        if KE < 300: hist_muon_signed_distance_0_to_300.Fill(dist_minus-dist_plus)
                        if KE > 300 and KE < 1000: hist_muon_signed_distance_300_to_1000.Fill(dist_minus-dist_plus)
                        if KE > 1000 and KE < 4000: hist_muon_signed_distance_1000_to_4000.Fill(dist_minus-dist_plus)
                        if KE > 4000 and KE < 5000: hist_muon_signed_distance_4000_to_5000.Fill(dist_minus-dist_plus)
                        if dist_plus < dist_minus : 
                            
                            chargeID = 13
                            n_correct_muon +=1
                            hist_correct_charge_id.Fill(KE)
                            hist_total_charge_id.Fill(KE)
                            hist_kemu_vs_thetamu_correct_chargeID_muon.Fill(KE,thetamu)
                            hist_kemu_vs_thetamu_total_chargeID_muon.Fill(KE,thetamu)
                            hist_correct_charge_id_with_hits.Fill(n_hits)
                            hist_correct_charge_id_with_hits_less_than_20.Fill(n_hits)
                            if KE >100:
                                n_correct_muon_100_cut+=1
                            if KE >300:
                                n_correct_muon_300_cut+=1
                            if KE >500:
                                n_correct_muon_500_cut+=1    
                            if KE >700:
                                n_correct_muon_700_cut+=1
                             #draw some spills pictures
                            if n_hit_total > 50 and event_display_num<100:
                                for i in range(n_hit_total):
                                    marker = ROOT.TMarker(muon_hits[3*i+2], muon_hits[3*i], 21)
                                    marker.SetMarkerSize(0.5)
                                    markers.append(marker)
                                
                                event_display_num+=1
                                # Create a canvas
                                canvas = ROOT.TCanvas(f"canvas{event_display_num}", "Markers", 800, 600)

                                # Optionally set up the axes if needed (e.g., using a TH2F histogram as a backdrop)
                                hist = ROOT.TH2F("hist", "Hits Positions", 100, 10000, 20000, 100, -5000, 5000)
                                hist.SetXTitle("Z Axis")  # Label the horizontal axis as Z
                                hist.SetYTitle("X Axis")  # Label the vertical axis as X
                                hist.Draw()
                                dunestyle.Simulation()

                                # Draw each marker on the canvas
                                for marker in markers:
                                    marker.Draw()
           
                                # Draw the rectangular box
                                box = ROOT.TBox(11362, -3300, 18314, 3300)
                                box.SetLineColor(ROOT.kBlack)  # Set box edge color to black
                                box.SetLineWidth(2)  # Set box edge thickness
                                box.SetFillStyle(0)  # Transparent fill
                                box.Draw("same")

                                # Draw the red dashed lines
                                line1 = ROOT.TLine(11362, -1750, 18314, -1750)  # Line at x = -1750
                                line2 = ROOT.TLine(11362, 1750, 18314, 1750)    # Line at x = 1750
                                line1.SetLineColor(ROOT.kRed)  # Set line color to red
                                line2.SetLineColor(ROOT.kRed)
                                line1.SetLineStyle(2)  # Dashed line
                                line2.SetLineStyle(2)
                                line1.Draw("same")
                                line2.Draw("same")
                                
                                # Update the canvas
                                canvas.Update()

                                # Define the path where the file should be saved
                                username = os.getenv('USER')  # Get the username from environment variables
                                outname = args.name
                                outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
                                outname_without_root = outname.replace(".root", "")
                                output_path = os.path.join(outdir, "eventdisplay/example", outname_without_root)
                                

                                # Ensure the directory exists
                                os.makedirs(output_path, exist_ok=True)

                                # Save the canvas as a PNG file in the specified directory
                                output_file = os.path.join(output_path, f"eventdisplay{event_display_num}.png")
                                canvas.SaveAs(output_file)

                        if dist_plus > dist_minus :
                            
                            chargeID = -13
                            n_incorrect_muon +=1
                            hist_incorrect_charge_id.Fill(KE)
                            hist_total_charge_id.Fill(KE)
                            hist_kemu_vs_thetamu_incorrect_chargeID_muon.Fill(KE,thetamu)
                            hist_kemu_vs_thetamu_total_chargeID_muon.Fill(KE,thetamu)
                            hist_incorrect_charge_id_with_hits.Fill(n_hits)
                            hist_incorrect_charge_id_with_hits_less_than_20.Fill(n_hits)
                            if KE >100:
                                n_incorrect_muon_100_cut+=1
                            if KE >300:
                                n_incorrect_muon_300_cut+=1
                            if KE >500:
                                n_incorrect_muon_500_cut+=1    
                            if KE >700:
                                n_incorrect_muon_700_cut+=1
                            
                            
                            #draw some spills pictures
                            if n_hit_total > 10 and event_display_num_unidentifiable_misidentified<100:
                                for i in range(n_hit_total):
                                    marker = ROOT.TMarker(muon_hits[3*i+2], muon_hits[3*i], 21)
                                    marker.SetMarkerSize(0.5)
                                    markers.append(marker)
                                
                                event_display_num_unidentifiable_misidentified+=1
                                # Create a canvas
                                canvas = ROOT.TCanvas(f"canvas{event_display_num_unidentifiable_misidentified}", "Markers", 800, 600)

                                # Optionally set up the axes if needed (e.g., using a TH2F histogram as a backdrop)
                                hist = ROOT.TH2F("hist", "Hits Positions", 100, 10000, 20000, 100, -5000, 5000)
                                hist.SetXTitle("Z Axis")  # Label the horizontal axis as Z
                                hist.SetYTitle("X Axis")  # Label the vertical axis as X
                                hist.Draw()
                                dunestyle.Simulation()

                                # Draw each marker on the canvas
                                for marker in markers:
                                    marker.Draw()
           
                                # Draw the rectangular box
                                box = ROOT.TBox(11362, -3300, 18314, 3300)
                                box.SetLineColor(ROOT.kBlack)  # Set box edge color to black
                                box.SetLineWidth(2)  # Set box edge thickness
                                box.SetFillStyle(0)  # Transparent fill
                                box.Draw("same")

                                # Draw the red dashed lines
                                line1 = ROOT.TLine(11362, -1750, 18314, -1750)  # Line at x = -1750
                                line2 = ROOT.TLine(11362, 1750, 18314, 1750)    # Line at x = 1750
                                line1.SetLineColor(ROOT.kRed)  # Set line color to red
                                line2.SetLineColor(ROOT.kRed)
                                line1.SetLineStyle(2)  # Dashed line
                                line2.SetLineStyle(2)
                                line1.Draw("same")
                                line2.Draw("same")
                                
                                # Update the canvas
                                canvas.Update()

                                # Define the path where the file should be saved
                                username = os.getenv('USER')  # Get the username from environment variables
                                outname = args.name
                                outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
                                outname_without_root = outname.replace(".root", "")
                                output_path = os.path.join(outdir, "eventdisplay/misidentified", outname_without_root)
                                

                                # Ensure the directory exists
                                os.makedirs(output_path, exist_ok=True)

                                # Save the canvas as a PNG file in the specified directory
                                output_file = os.path.join(output_path, f"eventdisplay{event_display_num_unidentifiable_misidentified}.png")
                                canvas.SaveAs(output_file)

                

                        if chargeID == -999: 
                            num_unidentifiable_muon += 1
                            hist_unidentifiable_muon.Fill(n_hits)
                            hist_unidentifiable_muon_less_than_20_hits.Fill(n_hits)
                            if n_hits > 20 and event_display_num_unidentifiable<100:
                                for i in range(n_hits):
                                    marker = ROOT.TMarker(muon_hits[3*i+2], muon_hits[3*i], 21)
                                    marker.SetMarkerSize(0.5)
                                    markers.append(marker)
                                
                                event_display_num_unidentifiable+=1
                                # Create a canvas
                                canvas = ROOT.TCanvas(f"canvas{event_display_num_unidentifiable}", "Markers", 800, 600)

                                # Optionally set up the axes if needed (e.g., using a TH2F histogram as a backdrop)
                                hist = ROOT.TH2F("hist", "Hits Positions", 100, 10000, 20000, 100, -4000, 4000)
                                hist.SetXTitle("Z Axis")  # Label the horizontal axis as Z
                                hist.SetYTitle("X Axis")  # Label the vertical axis as X
                                hist.Draw()
                                dunestyle.Simulation()

                                # Draw each marker on the canvas
                                for marker in markers:
                                    marker.Draw()

                                # Update the canvas
                                canvas.Update()

                                # Define the path where the file should be saved
                                username = os.getenv('USER')  # Get the username from environment variables
                                outname = args.name
                                outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
                                outname_without_root = outname.replace(".root", "")
                                output_path = os.path.join(outdir, "eventdisplay", outname_without_root)
                                

                                # Ensure the directory exists
                                os.makedirs(output_path, exist_ok=True)

                                # Save the canvas as a PNG file in the specified directory
                                output_file = os.path.join(output_path, f"eventdisplay{event_display_num_unidentifiable}.png")
                                canvas.SaveAs(output_file)
                             #Now judge whether the particle is a muon or an antimuon 
                    
                    
                    if truth.PDG[index] == -13:
                        hist_hits_vs_energy_antimuon.Fill(KE, n_hit_total)
                        hist_antimuon_signed_distance.Fill(dist_minus-dist_plus)
                        if KE < 300: hist_antimuon_signed_distance_0_to_300.Fill(dist_minus-dist_plus)
                        if KE > 300 and KE < 1000: hist_antimuon_signed_distance_300_to_1000.Fill(dist_minus-dist_plus)
                        if KE > 1000 and KE < 4000: hist_antimuon_signed_distance_1000_to_4000.Fill(dist_minus-dist_plus)
                        if KE > 4000 and KE < 5000: hist_antimuon_signed_distance_4000_to_5000.Fill(dist_minus-dist_plus)
                        
                        if dist_plus < dist_minus : 
                            
                            chargeID = 13
                            n_incorrect_antimuon +=1
                            hist_incorrect_charge_id_antimuon.Fill(KE)
                            hist_total_charge_id_antimuon.Fill(KE)
                            hist_kemu_vs_thetamu_incorrect_chargeID_antimuon.Fill(KE,thetamu)
                            hist_kemu_vs_thetamu_total_chargeID_antimuon.Fill(KE,thetamu)
                            if KE >100:
                                n_incorrect_antimuon_100_cut+=1
                            if KE >300:
                                n_incorrect_antimuon_300_cut+=1
                            if KE >500:
                                n_incorrect_antimuon_500_cut+=1    
                            if KE >700:
                                n_incorrect_antimuon_700_cut+=1
                            
                        if dist_plus > dist_minus : 
                            
                            chargeID = -13
                            n_correct_antimuon +=1
                            hist_correct_charge_id_antimuon.Fill(KE)
                            hist_total_charge_id_antimuon.Fill(KE)
                            hist_kemu_vs_thetamu_correct_chargeID_antimuon.Fill(KE,thetamu)
                            hist_kemu_vs_thetamu_total_chargeID_antimuon.Fill(KE,thetamu)
                            if KE >100:
                                n_correct_antimuon_100_cut+=1
                            if KE >300:
                                n_correct_antimuon_300_cut+=1
                            if KE >500:
                                n_correct_antimuon_500_cut+=1    
                            if KE >700:
                                n_correct_antimuon_700_cut+=1
                            

                        if chargeID == -999: num_unidentifiable_antimuon += 1
            
                  


            
            
   
    hist_correct_charge_id_percentage.Divide(hist_correct_charge_id,hist_total_charge_id)
    hist_correct_charge_id_percentage_antimuon.Divide(hist_correct_charge_id_antimuon,hist_total_charge_id_antimuon)

    hist_kemu_vs_thetamu_correct_chargeID_fraction_muon.Divide(hist_kemu_vs_thetamu_correct_chargeID_muon,hist_kemu_vs_thetamu_total_chargeID_muon)
    hist_kemu_vs_thetamu_correct_chargeID_fraction_antimuon.Divide(hist_kemu_vs_thetamu_correct_chargeID_antimuon,hist_kemu_vs_thetamu_total_chargeID_antimuon)

    hist_total_muons.Add(hist_reconstructed_muons)
    hist_total_muons.Add(hist_unreconstructed_muons)
    hist_reconstructed_muons_percentage.Divide(hist_reconstructed_muons,hist_total_muons)

    hist_kemu_vs_thetamu_correct_chargeID_comparison_muon_antimuon.Divide(hist_kemu_vs_thetamu_correct_chargeID_fraction_muon,hist_kemu_vs_thetamu_correct_chargeID_fraction_antimuon)
    hist_correct_charge_id_percentage_muon_over_antimuon.Divide(hist_correct_charge_id_percentage,hist_correct_charge_id_percentage_antimuon)

    if (n_correct_muon+n_incorrect_muon) > 0:
        correct_percentage_total = n_correct_muon/(n_correct_muon+n_incorrect_muon)
    if (n_correct_antimuon+n_incorrect_antimuon) >0:
        correct_percentage_total_antimuon= n_correct_antimuon/(n_correct_antimuon+n_incorrect_antimuon)
    if (n_correct_muon + n_incorrect_muon + num_unidentifiable_muon)>0:
        gross_correct_percentage_total = n_correct_muon/(n_correct_muon + n_incorrect_muon + num_unidentifiable_muon)
    if (n_correct_antimuon + n_incorrect_antimuon + num_unidentifiable_antimuon)>0:
        gross_correct_percentage_total_antimuon = n_correct_antimuon/(n_correct_antimuon + n_incorrect_antimuon + num_unidentifiable_antimuon)

    # correct_percentage_total_100_cut = n_correct_muon_100_cut /(n_correct_muon_100_cut + n_incorrect_muon_100_cut)
    # correct_percentage_total_300_cut = n_correct_muon_300_cut /(n_correct_muon_300_cut + n_incorrect_muon_300_cut)
    # correct_percentage_total_500_cut = n_correct_muon_500_cut /(n_correct_muon_500_cut + n_incorrect_muon_500_cut)
    # correct_percentage_total_700_cut = n_correct_muon_700_cut /(n_correct_muon_700_cut + n_incorrect_muon_700_cut)

    # correct_percentage_total_100_cut_antimuon = n_correct_antimuon_100_cut /(n_correct_antimuon_100_cut + n_incorrect_antimuon_100_cut)
    # correct_percentage_total_300_cut_antimuon = n_correct_antimuon_300_cut /(n_correct_antimuon_300_cut + n_incorrect_antimuon_300_cut)
    # correct_percentage_total_500_cut_antimuon = n_correct_antimuon_500_cut /(n_correct_antimuon_500_cut + n_incorrect_antimuon_500_cut)
    # correct_percentage_total_700_cut_antimuon = n_correct_antimuon_700_cut /(n_correct_antimuon_700_cut + n_incorrect_antimuon_700_cut)



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
        dunestyle.Simulation()
    tf.Close()
    print("Done saving.")
    
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
    base_name = outname[:-5]  # Remove the last 5 characters (.root)
    output_filename = base_name + '.txt'
    file_path = f"/exp/dune/app/users/{username}/dune-tms_hists/{output_filename}"

    # Open the file for writing
    with open(file_path, 'w') as file:
    # Write the output to the file
        file.write(f"N true muons: {n_true_muons}\n")
        file.write(f"N true antimuons: {n_true_antimuons}\n")
        file.write(f"Total correct percentage for muons: {correct_percentage_total}\n")
        file.write(f"Total correct percentage for antimuons: {correct_percentage_total_antimuon}\n")
        file.write(f"Correct chargeID muons: {n_correct_muon}\n")
        file.write(f"correct percentage for muons(100MeV cut): {correct_percentage_total_100_cut}\n")
        file.write(f"correct percentage for muons(300MeV cut): {correct_percentage_total_300_cut}\n")
        file.write(f"correct percentage for muons(500MeV cut): {correct_percentage_total_500_cut}\n")
        file.write(f"correct percentage for muons(700MeV cut): {correct_percentage_total_700_cut}\n")
        file.write(f"correct percentage for antimuons(100MeV cut): {correct_percentage_total_100_cut_antimuon}\n")
        file.write(f"correct percentage for antimuons(300MeV cut): {correct_percentage_total_300_cut_antimuon}\n")
        file.write(f"correct percentage for antimuons(500MeV cut): {correct_percentage_total_500_cut_antimuon}\n")
        file.write(f"correct percentage for antimuons(700MeV cut): {correct_percentage_total_700_cut_antimuon}\n")
        file.write(f"Incorrect chargeID muons: {n_incorrect_muon}\n")
        file.write(f"Correct chargeID antimuons: {n_correct_antimuon}\n")
        file.write(f"Incorrect chargeID antimuons: {n_incorrect_antimuon}\n")
        file.write(f"N lar start tms end muons: {n_muon_total_lar_start_tms_end}\n")
        file.write(f"N lar start tms end antimuons: {n_antimuon_total_lar_start_tms_end}\n")
        # file.write(f"region1 contained muons: {region1_muon_number}\n")
        # file.write(f"region2 contained muons: {region2_muon_number}\n")
        # file.write(f"region3 contained muons: {region3_muon_number}\n")
        file.write(f"unidentifiable muons: {num_unidentifiable_muon}\n")
        file.write(f"unidentifiable antimuons: {num_unidentifiable_antimuon}\n")
        file.write(f"gross correct percentage for muons: {gross_correct_percentage_total}\n")
        file.write(f"gross correct percentage for antimuons: {gross_correct_percentage_total_antimuon}\n")
        file.write(f"default to -999 times: {n_default_number}\n")
        file.write(f"reco muons: {reco_muons}\n")
        file.write(f"unreco muons: {unreco_muons}\n")
        if n_muon_total_lar_start_tms_end:
            file.write(f"reco percentage: {reco_muons/n_muon_total_lar_start_tms_end}\n")

        

    
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
    dunestyle.Simulation()


    
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
            dunestyle.Simulation()
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
    









