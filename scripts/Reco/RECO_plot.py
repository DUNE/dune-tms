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
    if not -3520 < x < 3520: is_inside = False
    if not -3864 < y < 1159: is_inside = False
    if only_thin_section:
        if not 11362 < z < 13500: is_inside = False
    else:
        if not 11362 < z < 18314: is_inside = False
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

    #kalman chargeID
    hist_correct_charge_id_kalman = ROOT.TH1D("hist_correct_charge_id_kalman", "Muons correct charge number (using reco, kalman);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_incorrect_charge_id_kalman = ROOT.TH1D("hist_incorrect_charge_id_kalman", "Muons incorrect charge number plot(using reco, kalman);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_total_charge_id_kalman = ROOT.TH1D("hist_total_charge_id_kalman", "Muons total charge number (using reco, kalman);True muon TMS Entry KE (MeV);Number of muons", 25, 0, 5000)
    hist_correct_charge_id_percentage_kalman = ROOT.TH1D("hist_correct_charge_id_percentage_kalman", "Muons correct charge ID percentage (using reco, kalman);True muon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_kalman.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_kalman.SetMaximum(1.2)  # Upper limit of Y-axis
    
    hist_correct_charge_id_antimuon_kalman = ROOT.TH1D("hist_correct_charge_id_antimuon_kalman", "Antimuons correct charge number (using reco, kalman);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_incorrect_charge_id_antimuon_kalman = ROOT.TH1D("hist_incorrect_charge_id_antimuon_kalman", "Antimuons incorrect charge number (using reco, kalman;True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_total_charge_id_antimuon_kalman = ROOT.TH1D("hist_total_charge_id_antimuon_kalman", "Antimuons total charge number (using reco, kalman);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_kalman = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon_kalman", "Antimuons correct charge ID percentage (using reco, kalman);True antimuon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_kalman.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_antimuon_kalman.SetMaximum(1.2)  # Upper limit of Y-axis
    
    hist_correct_charge_id_antimuon_kalman_central = ROOT.TH1D("hist_correct_charge_id_antimuon_kalman_central", "Antimuons correct charge number, central region (using reco, kalman);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_total_charge_id_antimuon_kalman_central = ROOT.TH1D("hist_total_charge_id_antimuon_kalman_central", "Antimuons total charge number, central region (using reco, kalman);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_kalman_central = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon_kalman_central", "Antimuons correct charge ID percentage, central region (using reco, kalman);True antimuon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_kalman_central.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_antimuon_kalman_central.SetMaximum(1.2)  # Upper limit of Y-axis

    hist_correct_charge_id_antimuon_central = ROOT.TH1D("hist_correct_charge_id_antimuon_central", "Antimuons correct charge number, central region (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_total_charge_id_antimuon_central = ROOT.TH1D("hist_total_charge_id_antimuon_central", "Antimuons total charge number, central region (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_central = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon_central", "Antimuons correct charge ID percentage, central region (using reco);True antimuon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon_central.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_antimuon_central.SetMaximum(1.2)  # Upper limit of Y-axis

    hist_incorrect_charge_id_antimuon_endpoint_reco = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_reco", "Antimuons incorrect charge, endpoint distribution in z(using reco);endpoint pos (mm);Number of antimuons", 25, 11000, 19000)
    hist_incorrect_charge_id_antimuon_endpoint_truth = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_truth", "Antimuons correct charge, endpoint distribution in z(using truth);endpoint pos (mm);Number of antimuons", 25, 11000, 19000)
    hist_incorrect_charge_id_antimuon_endpoint_reco_kalman = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_reco_kalman", "Antimuons incorrect charge, endpoint distribution in z(using reco,kalman);endpoint pos (mm);Number of antimuons", 25, 11000, 19000)
    hist_incorrect_charge_id_antimuon_endpoint_truth_kalman = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_truth_kalman", "Antimuons correct charge, endpoint distribution in z(using truth,kalman);endpoint pos (mm);Number of antimuons", 25, 11000, 19000)
    hist_incorrect_charge_id_antimuon_endpoint_truth_x = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_truth_x", "Antimuons correct charge, endpoint distribution in x(using truth);endpoint pos (mm);Number of antimuons", 25, -3000, 3000)
    hist_incorrect_charge_id_antimuon_endpoint_truth_kalman_x = ROOT.TH1D("hist_incorrect_charge_id_antimuon_endpoint_truth_kalman_x", "Antimuons correct charge, endpoint distribution in x(using truth,kalman);endpoint pos (mm);Number of antimuons", 25, -3000, 3000)

    hist_correct_charge_id_antimuon = ROOT.TH1D("hist_correct_charge_id_antimuon", "Antimuons correct charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_incorrect_charge_id_antimuon = ROOT.TH1D("hist_incorrect_charge_id_antimuon", "Antimuons incorrect charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_total_charge_id_antimuon = ROOT.TH1D("hist_total_charge_id_antimuon", "Antimuons total charge number (using reco);True antimuon TMS Entry KE (MeV);Number of antimuons", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon = ROOT.TH1D("hist_correct_charge_id_percentage_antimuon", "Antimuons correct charge ID percentage (using reco);True antimuon TMS Entry KE (MeV);Fraction", 25, 0, 5000)
    hist_correct_charge_id_percentage_antimuon.SetMinimum(0)  # Lower limit of Y-axis
    hist_correct_charge_id_percentage_antimuon.SetMaximum(1.2)  # Upper limit of Y-axis

    hist_endpoint_resolution_with_kalman = ROOT.TH1D("hist_endpoint_resolution_with_kalman", "muons endpoint resolution (with kalman);reco_z - true_z (mm);Number of muons", 20, -1000, 1000)
    hist_endpoint_resolution_with_kalman_antimuon = ROOT.TH1D("hist_endpoint_resolution_with_kalman_antimuon", "antimuons endpoint resolution (with kalman);reco_z - true_z (mm);Number of antimuons", 20, -1000, 1000)
    hist_endpoint_resolution_no_kalman = ROOT.TH1D("hist_endpoint_resolution_no_kalman", "muons endpoint resolution (without kalman);reco_z - true_z (mm);Number of muons", 20, -1000, 1000)
    hist_endpoint_resolution_no_kalman_antimuon = ROOT.TH1D("hist_endpoint_resolution_no_kalman_antimuon", "antimuons endpoint resolution (without kalman);reco_z - true_z (mm);Number of antimuons", 20, -1000, 1000)
    hist_endpoint_resolution_with_kalman_antimuon1 = ROOT.TH1D("hist_endpoint_resolution_with_kalman_antimuon1", "antimuons endpoint resolution (with kalman);reco_z - true_z (mm);Number of antimuons", 20, -1000, 1000)
    hist_endpoint_resolution_with_kalman_antimuon2 = ROOT.TH1D("hist_endpoint_resolution_with_kalman_antimuon2", "antimuons endpoint resolution (with kalman);reco_z - true_z (mm);Number of antimuons", 20, -1000, 1000)
    
    hist_kalman_chi2_plus_diff = ROOT.TH1D("hist_kalman_chi2_plus_diff", "muons chi2_plus diff with chi2;chi2_plus - chi2;Number of muons", 50, -20, 20)
    hist_kalman_chi2_minus_diff_muon = ROOT.TH1D("hist_kalman_chi2_minus_diff_muon", "muons chi2_minus diff with chi2;chi2_minus - chi2;Number of muons", 50, -20, 20)
    hist_kalman_chi2_minus_diff = ROOT.TH1D("hist_kalman_chi2_minus_diff", "antimuons chi2_minus diff with chi2;chi2_minus - chi2;Number of antimuons", 50, -20, 20)
    hist_kalman_chi2_plus_diff_antimuon = ROOT.TH1D("hist_kalman_chi2_plus_diff_antimuon", "antimuons chi2_plus diff with chi2;chi2_plus - chi2;Number of antimuons", 50, -20, 20)

    hist_kalman_chi2_plus_diff_per_dof = ROOT.TH1D("hist_kalman_chi2_plus_diff_per_dof", "muons chi2_plus diff with chi2 per node;(chi2_plus - chi2)/nkalmannode;Number of muons", 50, -5, 5)
    hist_kalman_chi2_minus_diff_per_dof = ROOT.TH1D("hist_kalman_chi2_minus_diff_per_dof", "muons chi2_minus diff with chi2 per node;(chi2_minus - chi2)/nkalmannode;Number of muons", 50, -5, 5)
    hist_kalman_chi2_minus_diff_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_minus_diff_per_dof_antimuon", "antimuons chi2_minus diff with chi2 per node;(chi2_minus - chi2)/nkalmannode;Number of antimuons", 50, -5, 5)
    hist_kalman_chi2_plus_diff_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_plus_diff_per_dof_antimuon", "antimuons chi2_plus diff with chi2 per node;(chi2_plus - chi2)/nkalmannode;Number of antimuons", 50, -5, 5)

    hist_kalman_chi2_plus_minus_diff_per_dof = ROOT.TH1D("hist_kalman_chi2_plus_minus_diff_per_dof", "muons chi2_plus diff with chi2_minus per node;(chi2_plus - chi2_minus)/nkalmannode;Number of muons", 50, -1, 1)
    hist_kalman_chi2_plus_minus_diff_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_plus_minus_diff_per_dof_antimuon", "antimuons chi2_plus diff with chi2_minus per node;(chi2_plus - chi2_minus)/nkalmannode;Number of antimuons", 50, -1, 1)
    
    
    hist_kalmannode_plus_diff_muon = ROOT.TH1D("hist_kalmannode_plus_diff_muon", "muons kalmannode_plus diff with kalmannode;kalmannode_plus - kalmannode;Number of muons", 50, -10, 10)
    hist_kalmannode_minus_diff_muon = ROOT.TH1D("hist_kalmannode_minus_diff_muon", "muons kalmannode_minus diff with kalmannode;kalmannode_minus - kalmannode;Number of muons", 50, -10, 10)
    hist_kalmannode_plus_diff_antimuon = ROOT.TH1D("hist_kalmannode_plus_diff_antimuon", "antimuons kalmannode_plus diff with kalmannode;kalmannode_plus - kalmannode;Number of antimuons", 50, -10, 10)
    hist_kalmannode_minus_diff_antimuon = ROOT.TH1D("hist_kalmannode_minus_diff_antimuon", "antimuons kalmannode_minus diff with kalmannode;kalmannode_minus - kalmannode;Number of antimuons", 50, -10, 10)

    hist_kalmannode_muon = ROOT.TH1D("hist_kalmannode_muon", "muons kalmannode; Number of kalmannode;Number of muons", 50, -10, 100)
    hist_kalmannode_antimuon = ROOT.TH1D("hist_kalmannode_antimuon", "antimuons kalmannode; Number of kalmannode;Number of antimuons", 50, -10, 100)
    
    hist_kalman_chi2 = ROOT.TH1D("hist_kalman_chi2", "muons chi2;chi2;Number of muons", 50, -50, 50)
    hist_kalman_chi2_plus = ROOT.TH1D("hist_kalman_chi2_plus", "muons chi2_plus;chi2_plus;Number of muons", 50, -50, 50)
    hist_kalman_chi2_minus = ROOT.TH1D("hist_kalman_chi2_minus", "muons chi2_minus;chi2_minus;Number of muons", 50, -50, 50)

    hist_kalman_chi2_per_dof = ROOT.TH1D("hist_kalman_chi2_per_dof", "muons chi2 per node;chi2_per_node;Number of muons", 50, -5, 5)
    hist_kalman_chi2_plus_per_dof = ROOT.TH1D("hist_kalman_chi2_plus_per_dof", "muons chi2_plus per node;chi2_plus_per_node;Number of muons", 50, -5, 5)
    hist_kalman_chi2_minus_per_dof = ROOT.TH1D("hist_kalman_chi2_minus_per_dof", "muons chi2_minus per node;chi2_minus_per_node;Number of muons", 50, -5, 5)
 
    hist_kalman_chi2_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_per_dof_antimuon", "antimuons chi2 per node;chi2_per_node;Number of antimuons", 50, -5, 5)
    hist_kalman_chi2_plus_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_plus_per_dof_antimuon", "antimuons chi2_plus per node;chi2_plus_per_node;Number of antimuons", 50, -5, 5)
    hist_kalman_chi2_minus_per_dof_antimuon = ROOT.TH1D("hist_kalman_chi2_minus_per_dof_antimuon", "antimuons chi2_minus per node;chi2_minus_per_node;Number of antimuons", 50, -5, 5)

    hist_kalman_chi2_antimuon = ROOT.TH1D("hist_kalman_chi2_antimuon", "antimuons chi2;chi2;Number of antimuons", 50, -50, 50)
    hist_kalman_chi2_plus_antimuon = ROOT.TH1D("hist_kalman_chi2_plus_antimuon", "antimuons chi2_plus;chi2_plus;Number of antimuons", 50, -50, 50)
    hist_kalman_chi2_minus_antimuon = ROOT.TH1D("hist_kalman_chi2_minus_antimuon", "antimuons chi2_minus;chi2_minus;Number of antimuons", 50, -50, 50)
    
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
    n_unidentifiable_muon = 0
    n_correct_charge_antimuon_geo = 0
    n_incorrect_charge_antimuon_geo = 0
    n_correct_charge_muon_geo_per = 0
    n_correct_charge_antimuon_geo_per = 0
    
    n_correct_charge_muon = 0
    n_incorrect_charge_muon = 0
    n_unidentifiable_muon = 0
    n_correct_charge_antimuon = 0
    n_incorrect_charge_antimuon = 0
    n_unidentifiable_antimuon = 0
    n=0
    n_kalmanpdg_muon = 0
    n_kalmanpdg_correct_muon = 0
    n_kalmanpdg_incorrect_muon = 0
    n_unidentifiable_kalmanpdg_muon = 0
    n_kalmanpdg_muon_correct_percentage = 0
    n_kalmanpdg_antimuon = 0
    n_kalmanpdg_correct_antimuon = 0
    n_kalmanpdg_incorrect_antimuon = 0
    n_unidentifiable_kalmanpdg_antimuon = 0
    n_kalmanpdg_antimuon_correct_percentage = 0
    n_chi_neg = 0
    n_chi_pos = 0
    chi_diff_muon = -9999
    chi_diff_antimuon = -9999
    event_display_num = 0
    
    # Now loop over all events
    for i, event in enumerate(c):
        if i > nevents: break
        if i % print_every == 0 : print(f"On {i} / {nevents}")
        if truth != None: 
            truth.GetEntry(i)
            
        if reco != None: 
            reco.GetEntry(i)
      
        for index, value in enumerate(truth.RecoTrackPrimaryParticlePDG):
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
                

                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    # if KE <10:
                    #     break
                    # if len(reco.TrackHitPos) < 3:               
                    #     break
                    nTracks = reco.nTracks
                    if nTracks >0:
                        z_resolution = reco.EndPos[3*index + 2] - z_end
                        # hist_endpoint_resolution_with_kalman.Fill(z_resolution)
                        # hist_kalmannode_plus_diff_muon.Fill(reco.nKalmanNodes_plus[index] - reco.nKalmanNodes[index])
                        # hist_kalmannode_minus_diff_muon.Fill(reco.nKalmanNodes_minus[index] - reco.nKalmanNodes[index])
                        # hist_kalmannode_muon.Fill(reco.nKalmanNodes[index])
                        # hist_kalman_chi2.Fill(reco.Chi2[index])
                        # hist_kalman_chi2_plus.Fill(reco.Chi2_plus[index])
                        # hist_kalman_chi2_minus.Fill(reco.Chi2_minus[index])
                        
                        # hist_endpoint_resolution_no_kalman.Fill(z_resolution)
                        #print(f"true PDG = {truth.RecoTrackPrimaryParticlePDG[index]}, reco PDG_kalman = {reco.Charge_Kalman[index]}, reco PDG = {reco.Charge[index]}, chi2={reco.Chi2[index]:.1f}, chi2_plus={reco.Chi2_plus[index]:.1f}, chi2_minus={reco.Chi2_minus[index]:.1f}")
                        # chi_diff_muon = reco.Chi2_plus[index] - reco.Chi2[index]
                        # hist_kalman_chi2_plus_diff.Fill(chi_diff_muon)
                        # hist_kalman_chi2_minus_diff_muon.Fill(reco.Chi2_minus[index] - reco.Chi2[index])
                        # charge_value_kalman = reco.Charge_Kalman[index]
                        charge_value = reco.Charge[index]
                        hist_total_charge_id.Fill(KE)
                        hist_total_charge_id_kalman.Fill(KE)
                        # if reco.nKalmanNodes[index] > 0:
                        #    hist_kalman_chi2_per_dof.Fill((reco.Chi2[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_plus_per_dof.Fill((reco.Chi2_plus[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_minus_per_dof.Fill((reco.Chi2_minus[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_plus_diff_per_dof.Fill((reco.Chi2_plus[index] - reco.Chi2[index])/reco.nKalmanNodes[index]) 
                        #    hist_kalman_chi2_minus_diff_per_dof.Fill((reco.Chi2_minus[index] - reco.Chi2[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_plus_minus_diff_per_dof.Fill((reco.Chi2_plus[index] - reco.Chi2_minus[index])/reco.nKalmanNodes[index]) 
                           

 

                        if charge_value == 13: 
                            hist_correct_charge_id.Fill(KE)
                            n_correct_charge_muon_geo+=1
                        else: 
                            n_incorrect_charge_muon_geo += 1    
                        
                        # if charge_value_kalman == 13: 
                        #     hist_correct_charge_id_kalman.Fill(KE)
                        #     n_correct_charge_muon += 1

                        # elif charge_value_kalman == -13: n_incorrect_charge_muon += 1
                        # elif charge_value_kalman == -999999999 : n_unidentifiable_muon += 1
                               #draw some spills pictures
                            num_hits = len(reco.TrackHitTruePos)//3   
                            if num_hits > 30 and event_display_num<200  and nTracks == 1:
                                for i in range(num_hits):
                                    marker = ROOT.TMarker(reco.TrackHitPos[3*i+2], reco.TrackHitPos[3*i], 21)
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
                

                if inside_lar(x_start,y_start,z_start) and inside_tms(x_end,y_end,z_end):
                    # if KE <10:
                    #     break
                    # if len(reco.TrackHitPos) < 3:               
                    #     break
                    markers = []
                    nTracks = reco.nTracks
                    if nTracks >0:
                        z_resolution = reco.EndPos[3*index + 2] - z_end
                        z_resolution1 = reco.TrackHitPos[3*index +2] - z_end
                        z_resolution2 = reco.KalmanPos[3*index +2] - z_end
                    
                        # hist_kalmannode_plus_diff_antimuon.Fill(reco.nKalmanNodes_plus[index] - reco.nKalmanNodes[index])
                        # hist_kalmannode_minus_diff_antimuon.Fill(reco.nKalmanNodes_minus[index] - reco.nKalmanNodes[index])
                        # hist_kalmannode_antimuon.Fill(reco.nKalmanNodes[index])
                        
                        # hist_kalman_chi2_antimuon.Fill(reco.Chi2[index])
                        # hist_kalman_chi2_plus_antimuon.Fill(reco.Chi2_plus[index])
                        # hist_kalman_chi2_minus_antimuon.Fill(reco.Chi2_minus[index])

                        # hist_endpoint_resolution_with_kalman_antimuon1.Fill(z_resolution1)
                        # hist_endpoint_resolution_with_kalman_antimuon2.Fill(z_resolution2)

                        # hist_endpoint_resolution_with_kalman_antimuon.Fill(z_resolution)
                        # hist_endpoint_resolution_no_kalman_antimuon.Fill(z_resolution)
                        #print(f"true PDG = {truth.RecoTrackPrimaryParticlePDG[index]}, reco PDG_kalman = {reco.Charge_Kalman[index]}, reco PDG = {reco.Charge[index]}, chi2={reco.Chi2[index]:.1f}, chi2_plus={reco.Chi2_plus[index]:.1f}, chi2_minus={reco.Chi2_minus[index]:.1f}")

                        if_contained_in_central_region = True
                        num_hits = len(reco.TrackHitTruePos)//3
                        for i in range(num_hits):
                            if not region2(reco.TrackHitTruePos[3 * i]): if_contained_in_central_region = False
                        if if_contained_in_central_region == True:
                            # hist_total_charge_id_antimuon_kalman_central.Fill(KE)
                            hist_total_charge_id_antimuon_central.Fill(KE)
                            if reco.Charge[index]==13:
                                hist_correct_charge_id_antimuon_central.Fill(KE)
                            # if reco.Charge_Kalman[index]==-13:
                            #     hist_correct_charge_id_antimuon_kalman_central.Fill(KE)


                        # chi_diff_antimuon = reco.Chi2_minus[index] - reco.Chi2[index]
                        # hist_kalman_chi2_minus_diff.Fill(chi_diff_antimuon)
                        # hist_kalman_chi2_plus_diff_antimuon.Fill(reco.Chi2_plus[index] - reco.Chi2[index])
                        # charge_value_kalman = reco.Charge_Kalman[index]
                        charge_value = reco.Charge[index]
                        # hist_total_charge_id_antimuon.Fill(KE)
                        # hist_total_charge_id_antimuon_kalman.Fill(KE)
                        # if reco.nKalmanNodes[index] > 0:
                        #    hist_kalman_chi2_per_dof_antimuon.Fill((reco.Chi2[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_plus_per_dof_antimuon.Fill((reco.Chi2_plus[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_minus_per_dof_antimuon.Fill((reco.Chi2_minus[index])/reco.nKalmanNodes[index])
                        #    hist_kalman_chi2_plus_diff_per_dof_antimuon.Fill((reco.Chi2_plus[index] - reco.Chi2[index])/reco.nKalmanNodes[index]) 
                        #    hist_kalman_chi2_minus_diff_per_dof_antimuon.Fill((reco.Chi2_minus[index] - reco.Chi2[index])/reco.nKalmanNodes[index]) 
                        #    hist_kalman_chi2_plus_minus_diff_per_dof_antimuon.Fill((reco.Chi2_plus[index] - reco.Chi2_minus[index])/reco.nKalmanNodes[index]) 
                        
                        if charge_value == -13: 
                            hist_correct_charge_id_antimuon.Fill(KE)
                            n_correct_charge_antimuon_geo+=1

                            #draw some spills pictures
                            # if num_hits > 30 and event_display_num<200 and charge_value_kalman == 13 and nTracks == 1:
                            #     for i in range(num_hits):
                            #         marker = ROOT.TMarker(reco.TrackHitPos[3*i+2], reco.TrackHitPos[3*i], 21)
                            #         marker.SetMarkerSize(0.5)
                            #         markers.append(marker)
                                
                            #     event_display_num+=1
                            #     # Create a canvas
                            #     canvas = ROOT.TCanvas(f"canvas{event_display_num}", "Markers", 800, 600)

                            #     # Optionally set up the axes if needed (e.g., using a TH2F histogram as a backdrop)
                            #     hist = ROOT.TH2F("hist", "Hits Positions", 100, 10000, 20000, 100, -5000, 5000)
                            #     hist.SetXTitle("Z Axis")  # Label the horizontal axis as Z
                            #     hist.SetYTitle("X Axis")  # Label the vertical axis as X
                            #     hist.Draw()
                              

                            #     # Draw each marker on the canvas
                            #     for marker in markers:
                            #         marker.Draw()
           
                            #     # Draw the rectangular box
                            #     box = ROOT.TBox(11362, -3300, 18314, 3300)
                            #     box.SetLineColor(ROOT.kBlack)  # Set box edge color to black
                            #     box.SetLineWidth(2)  # Set box edge thickness
                            #     box.SetFillStyle(0)  # Transparent fill
                            #     box.Draw("same")

                            #     # Draw the red dashed lines
                            #     line1 = ROOT.TLine(11362, -1750, 18314, -1750)  # Line at x = -1750
                            #     line2 = ROOT.TLine(11362, 1750, 18314, 1750)    # Line at x = 1750
                            #     line1.SetLineColor(ROOT.kRed)  # Set line color to red
                            #     line2.SetLineColor(ROOT.kRed)
                            #     line1.SetLineStyle(2)  # Dashed line
                            #     line2.SetLineStyle(2)
                            #     line1.Draw("same")
                            #     line2.Draw("same")
                                
                            #     # Update the canvas
                            #     canvas.Update()

                            #     # Define the path where the file should be saved
                            #     username = os.getenv('USER')  # Get the username from environment variables
                            #     outname = args.name
                            #     outdir = f"/exp/dune/app/users/{username}/dune-tms_hists"
                            #     outname_without_root = outname.replace(".root", "")
                            #     output_path = os.path.join(outdir, "eventdisplay/example", outname_without_root)
                                

                            #     # Ensure the directory exists
                            #     os.makedirs(output_path, exist_ok=True)

                            #     # Save the canvas as a PNG file in the specified directory
                            #     output_file = os.path.join(output_path, f"eventdisplay{event_display_num}.png")
                            #     canvas.SaveAs(output_file)
                        else:
                            n_incorrect_charge_antimuon_geo +=1 
                            hist_incorrect_charge_id_antimuon_endpoint_reco.Fill(reco.EndPos[2])
                            hist_incorrect_charge_id_antimuon_endpoint_truth.Fill(reco.TrackHitTruePos[2])
                            hist_incorrect_charge_id_antimuon_endpoint_truth_x.Fill(reco.TrackHitTruePos[0])
                        # if charge_value_kalman == -13: 
                        #     hist_correct_charge_id_antimuon_kalman.Fill(KE)
                        #     n_correct_charge_antimuon += 1

                        # elif charge_value_kalman == 13: 
                        #     n_incorrect_charge_antimuon += 1
                        #     hist_incorrect_charge_id_antimuon_endpoint_reco_kalman.Fill(reco.EndPos[2])
                        #     hist_incorrect_charge_id_antimuon_endpoint_truth_kalman.Fill(reco.TrackHitTruePos[2])
                        #     hist_incorrect_charge_id_antimuon_endpoint_truth_kalman_x.Fill(reco.TrackHitTruePos[0])
                            
                        # elif charge_value_kalman == -999999999 : n_unidentifiable_antimuon += 1
    
    n_kalmanpdg_muon = n_correct_charge_muon + n_incorrect_charge_muon
    n_kalmanpdg_antimuon = n_incorrect_charge_antimuon + n_correct_charge_antimuon
    if n_kalmanpdg_muon != 0: n_kalmanpdg_muon_correct_percentage = n_correct_charge_muon / n_kalmanpdg_muon
    if n_kalmanpdg_antimuon != 0: n_kalmanpdg_antimuon_correct_percentage = n_correct_charge_antimuon / n_kalmanpdg_antimuon

    if n_correct_charge_muon_geo >0:  n_correct_charge_muon_geo_per = n_correct_charge_muon_geo/(n_incorrect_charge_muon_geo +n_correct_charge_muon_geo)
    if n_correct_charge_antimuon_geo >0:  n_correct_charge_antimuon_geo_per = n_correct_charge_antimuon_geo/(n_incorrect_charge_antimuon_geo+ n_correct_charge_antimuon_geo)

    
    # Divide histograms to compute percentage histograms
    hist_correct_charge_id_percentage.Divide(hist_correct_charge_id, hist_total_charge_id)
    hist_correct_charge_id_percentage_antimuon.Divide(hist_correct_charge_id_antimuon, hist_total_charge_id_antimuon)
    hist_correct_charge_id_percentage_kalman.Divide(hist_correct_charge_id_kalman, hist_total_charge_id_kalman)
    hist_correct_charge_id_percentage_antimuon_kalman.Divide(hist_correct_charge_id_antimuon_kalman, hist_total_charge_id_antimuon_kalman)
    
    hist_correct_charge_id_percentage_antimuon_kalman_central.Divide(hist_correct_charge_id_antimuon_kalman_central,hist_total_charge_id_antimuon_kalman_central)
    hist_correct_charge_id_percentage_antimuon_central.Divide(hist_correct_charge_id_antimuon_central,hist_total_charge_id_antimuon_central)
    # Add error bars using binomial statistics for each percentage histogram
    # for hist_percentage, hist_correct, hist_total in [
    #     (hist_correct_charge_id_percentage, hist_correct_charge_id, hist_total_charge_id),
    #     (hist_correct_charge_id_percentage_antimuon, hist_correct_charge_id_antimuon, hist_total_charge_id_antimuon),
    #     (hist_correct_charge_id_percentage_kalman, hist_correct_charge_id_kalman, hist_total_charge_id_kalman),
    #     (hist_correct_charge_id_percentage_antimuon_kalman, hist_correct_charge_id_antimuon_kalman, hist_total_charge_id_antimuon_kalman)
    # ]:
    #     for bin in range(1, hist_percentage.GetNbinsX() + 1):
    #         n_correct = hist_correct.GetBinContent(bin)
    #         n_total = hist_total.GetBinContent(bin)
            
    #         if n_total > 0:  # Avoid division by zero
    #             error = math.sqrt((n_correct / n_total) * (1 - (n_correct / n_total)) / n_total)
    #             hist_percentage.SetBinError(bin, error)
    #         else:
    #             hist_percentage.SetBinError(bin, 0)  # No data, no error
           

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
        
        file.write(f"correct chargeID in reco (for muon): {n_correct_charge_muon_geo}\n")
        file.write(f"incorrect chargeID in reco (for muon): {n_incorrect_charge_muon_geo}\n")
        file.write(f"correct percentage geoPDG in reco (for muon): {n_correct_charge_muon_geo_per}\n")

        file.write(f"correct chargeID in reco (for antimuon): {n_correct_charge_antimuon_geo}\n")
        file.write(f"incorrect chargeID in reco (for antimuon): {n_incorrect_charge_antimuon_geo}\n")
        file.write(f"correct percentage geoPDG in reco (for antimuon): {n_correct_charge_antimuon_geo_per}\n")

        file.write(f"correct KalmanPDG in reco (for muon): {n_correct_charge_muon}\n")
        file.write(f"incorrect KalmanPDG in reco (for muon): {n_incorrect_charge_muon}\n")
        file.write(f"unidentifiable KalmanPDG in reco (for muon): {n_unidentifiable_muon}\n")
        file.write(f"correct percentage KalmanPDG in reco (for muon): {n_kalmanpdg_muon_correct_percentage}\n")

        file.write(f"correct KalmanPDG in reco (for antimuon): {n_correct_charge_antimuon}\n")
        file.write(f"incorrect KalmanPDG in reco (for antimuon): {n_incorrect_charge_antimuon}\n")
        file.write(f"unidentifiable KalmanPDG in reco (for antimuon): {n_unidentifiable_antimuon}\n")
        file.write(f"correct percentage KalmanPDG in reco (for antimuon): {n_kalmanpdg_antimuon_correct_percentage}\n")

        file.write(f"number of chisquareneg: {n_chi_neg}\n")
        file.write(f"number of chisquarepos: {n_chi_pos}\n")

    
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
    

