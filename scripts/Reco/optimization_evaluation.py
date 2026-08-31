import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll

import math # filter out nans

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-v0_8-poster')

mp.rc('axes', labelsize = 15)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 15) # fontsize of the tick labels
mp.rc('ytick', labelsize = 15) # fontsize of the tick labels

# dyslexia friendly background
mp.rcParams['axes.facecolor'] = '#fafafa'
mp.rcParams["figure.facecolor"] = '#fafafa'
mp.rcParams["savefig.facecolor"] = '#fafafa'

# off-black text
mp.rcParams['text.color'] = '#424242'
mp.rcParams['axes.labelcolor'] = '#424242'
mp.rcParams['xtick.color'] = '#424242'
mp.rcParams['ytick.color'] = '#424242'

### Actual function that loops through the spills
def draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained):
    if not os.path.exists(input_filename): raise ValueError(f"Cannor find input_filename {input_filename}")
    
    # Make sure we read in the correct file and have the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir):
        raise ValueError(f"Could not make out_dir {out_dir}")
            
    # Read in the Reco_Tree that contains the TMS_Tracks
    r = ROOT.TChain("Reco_Tree")
    r.Add(input_filename)
    print("N entries:", r.GetEntries())
    if not r.GetEntries() > 0:
        print("Didn't get any entries, are you sure the input_filename is right?\n", input_filename)
    
    truth = ROOT.TChain("Truth_Info")
    truth.Add(input_filename)
    if not truth.GetEntries() > 0:
        print("Didn't get any entries in Truth_Info, are you sure the input_filename is right?\n", input_filename)
            
    max_n_spills = 129000 # TODO (old) add some meta info to output file with n spill info for file
    
    spill_number_cache = dict()
    n_events = r.GetEntries()
    
    # prepare arrays to be filled
    if plot_Start:
        Reco_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Primary_True_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
    if plot_End:
        Reco_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Primary_True_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
    if plot_Charge:
        Reco_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_Charge_stop = np.ones((n_events, 5), dtype = float) * -9999.
    if plot_Charge or plot_Angle:
        True_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_KE = np.ones((n_events, 5), dtype = float) * -9999.
        True_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.
        Reco_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.        
    True_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true
    Reco_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true

    # some simple counters for efficiency evaluation
    correct_tracks_reco = 0
    correct_reco_hits = 0
    correct_true_hits = 0
    count_muons = 0
    counter_leaving = 0
    
    # now fill the arrays
    for current_spill_number in range(max_n_spills):
        for i in range(n_events):
            try:
                spill_number = spill_number_cache[i]
                event = None
                true_event = None
            except KeyError:
                r.GetEntry(i)
                event = r
                truth.GetEntry(i)
                true_event = truth
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event == None:
                r.GetEntry(i)
                event = r
            if true_event == None:
                truth.GetEntry(i)
                true_event = truth

            # track start and end position
            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)
            True_Position_TMS_Start = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionEnteringTMS, dtype = np.float32)
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            True_Position_TMS_End = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionLeavingTMS, dtype = np.float32)
            # track direction/momentum at start
            if plot_Angle or plot_Charge:
                Reco_Track_StartDirection = np.frombuffer(event.StartDirection, dtype = np.float32)
                MomentumTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumTrackStart, dtype = np.float32)
            # charge related stuff
            if plot_Charge:
                Reco_Track_Charge = event.Charge
            True_PDG = true_event.PDG
            # classify if true muon actually enough distance in TMS and start in LAr
            LArFiducialTouch = true_event.RecoTrackPrimaryParticleLArFiducialStart
            Particle_PDG = true_event.LeptonPDG
            Muon_Start = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionStart, dtype = np.float32)
            Muon_End = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionEnd, dtype = np.float32)
            True_Momentum_Leaving = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumLeavingTMS, dtype = np.float32)
            # actual hits
            Reco_Hits = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype = np.float32)
            sum_reco_hits = event.nHits
            sum_true_hits = true_event.RecoTrackNHits            
            
            if (abs(Particle_PDG) == 13):
                if 4179.24 < Muon_Start[2] < 9135.88 and 11134 < Muon_End[2] < 18535:
                    if 371.77 > Muon_End[1] > -3076.23 and abs(Muon_End[0]) < 3491:
                        count_muons += 1
            
            nTracks = event.nTracks
            if nTracks <= 0: continue
            if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break

                # check if (anti-)muon as true primary particle and if origin in LAr
                if True: #np.abs(True_PDG[j]) == 13 and LArFiducialTouch[j]:
                    # check if (anti-)muon travesers at least 4 planes in TMS
                    if (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) >= 440.:
                        # if so, then this is a true muon: 0 -> 1
                        True_Muon_Track[i, j] = 1.
                        counter_correct = 0
                        for true_hits in range(int(sum_true_hits[j])):
                            for reco_hits in range(int(sum_reco_hits[j])):
                                if True_Hits[j*600 + true_hits*4 + 2] == Reco_Hits[j*600 + reco_hits*3 + 2]:
                                    if np.abs(True_Hits[j*600 + true_hits*4 + 0] - Reco_Hits[j*600 + reco_hits*3 + 0]) <= 2 * 36:
                                        counter_correct += 1

                        if counter_correct / sum_true_hits[j] >= 0.1:
                            correct_tracks_reco += 1
                            correct_reco_hits += counter_correct
                            correct_true_hits += sum_true_hits[j]
                        
                        # check if reconstructed tracks exist for this event
                        if StartPos.size != 0:
                            # if so, then add all identified as muons: 0 -> number tracks
                            Reco_Muon_Track[i, j] = 1. #StartPos.size / 3
                
                if True_Position_TMS_Start[j*4 + 0] > -8000. and not StartPos.size == 0:
                    # checking for muon tracks (minimal length for this are 20 planes traversed -> 890 mm in thin area
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 4 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 4 * 65.:
                            # normalization of true direction
                            if plot_Angle or plot_Charge:
                                magnitude = MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 + 1]**2 + MomentumTrackStart[j*4 + 2]**2
                                True_TrackDirection[i, j, 0] = dir_to_angle(MomentumTrackStart[j*4 + 0] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                True_TrackDirection[i, j, 1] = dir_to_angle(MomentumTrackStart[j*4 + 1] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                Reco_TrackDirection[i, j, 0] = dir_to_angle(Reco_Track_StartDirection[j*3 + 0], Reco_Track_StartDirection[j*3 + 2])
                                Reco_TrackDirection[i, j, 1] = dir_to_angle(Reco_Track_StartDirection[j*3 + 1], Reco_Track_StartDirection[j*3 + 2])
                                True_KE[i, j] = np.sqrt((MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 +1]**2 + MomentumTrackStart[j*4 + 2]**2 + 105.7**2) - 105.7)
                                True_Charge[i, j] = 13  #True_PDG[j]                            
                            if plot_Charge:
                                # Make sure that true muon stopped in TMS and is not leaving
                                if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                    Reco_Charge[i, j] = Reco_Track_Charge[j]
                                    True_Charge_stop[i, j] = 13 #True_PDG[j]
                            if plot_Start:
                                if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) >= 1:
                                    Reco_Start[i, j, 0] = StartPos[j*3 + 0]
                                    Reco_Start[i, j, 1] = StartPos[j*3 + 1]
                                    Reco_Start[i, j, 2] = StartPos[j*3 + 2]
                                    Primary_True_Start[i, j, 0] = True_Position_TMS_Start[j*4 + 0]
                                    Primary_True_Start[i, j, 1] = True_Position_TMS_Start[j*4 + 1]
                                    Primary_True_Start[i, j, 2] = True_Position_TMS_Start[j*4 + 2]
                if plot_End and Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 4 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 4 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
                            if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                Reco_End[i, j, 0] = EndPos[j*3 + 0]
                                Reco_End[i, j, 1] = EndPos[j*3 + 1]
                                Reco_End[i, j, 2] = EndPos[j*3 + 2]
                                Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                                Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                                Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                            else:
                                counter_leaving += 1
                if plot_End and not Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 4 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 4 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
                            if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) > 1:
                                counter_leaving += 1
                            Reco_End[i, j, 0] = EndPos[j*3 + 0]
                            Reco_End[i, j, 1] = EndPos[j*3 + 1]
                            Reco_End[i, j, 2] = EndPos[j*3 + 2]
                            Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                            Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                            Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]

    
    # filter out not filled indice
    if plot_Start:
        boolean_Reco_Start = (Reco_Start[:, :, 0] != -9999.)
        Reco_Start = Reco_Start[boolean_Reco_Start]
        boolean_Primary_Start = (Primary_True_Start[:, :, 0] != -9999.)
        Primary_True_Start = Primary_True_Start[boolean_Primary_Start]
    if plot_End:
        boolean_Reco_End = (Reco_End[:, :, 0] != -9999.)
        Reco_End = Reco_End[boolean_Reco_End]
        boolean_Primary_End = (Primary_True_End[:, :, 0] != -9999.)
        Primary_True_End = Primary_True_End[boolean_Primary_End]
    if plot_Charge:
        boolean_Reco_Charge = (Reco_Charge != -9999.)
        Reco_Charge = Reco_Charge[boolean_Reco_Charge]
        boolean_True_Charge_stop = (True_Charge_stop != -9999.)
        True_Charge_stop = True_Charge_stop[boolean_True_Charge_stop]
        True_KE_stop = True_KE[boolean_True_Charge_stop & (True_KE != -9999.)]
    if plot_Charge or plot_Angle:
        boolean_True_Charge = (True_Charge != -9999.)
        True_Charge = True_Charge[boolean_True_Charge]    
        boolean_True_TrackDirection = (True_TrackDirection[:, :, 0] != -9999.)
        True_TrackDirection = True_TrackDirection[boolean_True_TrackDirection]
        boolean_Reco_TrackDirection = (Reco_TrackDirection[:, :, 0] != -9999.)
        Reco_TrackDirection = Reco_TrackDirection[boolean_Reco_TrackDirection]
        boolean_True_KE = (True_KE != -9999.)
        True_KE = True_KE[boolean_True_KE]
    
    # flatten arrays
    if plot_Start:
        Reco_Start_x = Reco_Start[:, 0]
        Reco_Start_y = Reco_Start[:, 1]
        Reco_Start_z = Reco_Start[:, 2]
        Primary_True_Start_x = Primary_True_Start[:, 0]
        Primary_True_Start_y = Primary_True_Start[:, 1]
        Primary_True_Start_z = Primary_True_Start[:, 2]
    if plot_End:
        Reco_End_x = Reco_End[:, 0]
        Reco_End_y = Reco_End[:, 1]
        Reco_End_z = Reco_End[:, 2]
        Primary_True_End_x = Primary_True_End[:, 0]
        Primary_True_End_y = Primary_True_End[:, 1]
        Primary_True_End_z = Primary_True_End[:, 2]
    if plot_Angle or plot_Charge:
        True_TrackDirection_xz = True_TrackDirection[:, 0]
        True_TrackDirection_yz = True_TrackDirection[:, 1]
        Reco_TrackDirection_xz = Reco_TrackDirection[:, 0]
        Reco_TrackDirection_yz = Reco_TrackDirection[:, 1]

    
    # total number of events after filtering
    if plot_End:
        print("#events reconstruction: ", len(Reco_End), " + not stopped in TMS:", counter_leaving)
    print("true (anti-)muons: ", sum(True_Muon_Track), "  vs. reconstructed particles: ", sum(Reco_Muon_Track))
    print("correctly identified tracks: ", correct_tracks_reco)
    #print("  correct hits reco: ", correct_reco_hits, " vs. true hits: ", correct_true_hits, " -> ", correct_reco_hits / correct_true_hits * 100)
    print("total muons (outside of reco): ", count_muons)
    
    # subtract reconstruction from truth for all directions
    if plot_End:
        Diff_End_x = Primary_True_End_x - Reco_End_x
        Diff_End_y = Primary_True_End_y - Reco_End_y
        Diff_End_z = Primary_True_End_z - Reco_End_z
    
    if plot_Start:
        Diff_Start_x = Primary_True_Start_x - Reco_Start_x
        Diff_Start_y = Primary_True_Start_y - Reco_Start_y
        Diff_Start_z = Primary_True_Start_z - Reco_Start_z

    # create histograms for the differences
    if plot_End:
        Diff_End_x_hist, Diff_End_x_bins = np.histogram(Diff_End_x, bins = 50)
        Diff_End_y_hist, Diff_End_y_bins = np.histogram(Diff_End_y, bins = 50)
        Diff_End_z_hist, Diff_End_z_bins = np.histogram(Diff_End_z, bins = 50)
    
        Diff_End_x_histX, Diff_End_x_histY = histogram_arr_handle(Diff_End_x_hist, Diff_End_x_bins)
        Diff_End_y_histX, Diff_End_y_histY = histogram_arr_handle(Diff_End_y_hist, Diff_End_y_bins)
        Diff_End_z_histX, Diff_End_z_histY = histogram_arr_handle(Diff_End_z_hist, Diff_End_z_bins)

        mean_End_x = np.mean(Diff_End_x)
        mean_End_y = np.mean(Diff_End_y)
        mean_End_z = np.mean(Diff_End_z)
        std_End_x = np.std(Diff_End_x)
        std_End_y = np.std(Diff_End_y)
        std_End_z = np.std(Diff_End_z)
        
        print("End   x: ", mean_End_x, "±", std_End_x, "    y: ", mean_End_y, "±", std_End_y, "    z: ", mean_End_z, "±", std_End_z)

    if plot_Start:    
        Diff_Start_x_hist, Diff_Start_x_bins = np.histogram(Diff_Start_x, bins = 50)
        Diff_Start_y_hist, Diff_Start_y_bins = np.histogram(Diff_Start_y, bins = 50)
        Diff_Start_z_hist, Diff_Start_z_bins = np.histogram(Diff_Start_z, bins = 50)
    
        Diff_Start_x_histX, Diff_Start_x_histY = histogram_arr_handle(Diff_Start_x_hist, Diff_Start_x_bins)
        Diff_Start_y_histX, Diff_Start_y_histY = histogram_arr_handle(Diff_Start_y_hist, Diff_Start_y_bins)
        Diff_Start_z_histX, Diff_Start_z_histY = histogram_arr_handle(Diff_Start_z_hist, Diff_Start_z_bins)

        mean_Start_x = np.mean(Diff_Start_x)
        mean_Start_y = np.mean(Diff_Start_y)
        mean_Start_z = np.mean(Diff_Start_z)
        std_Start_x = np.std(Diff_Start_x)
        std_Start_y = np.std(Diff_Start_y)
        std_Start_z = np.std(Diff_Start_z)
    
    # plot
    if plot_Start:
        mp.plot(Diff_Start_x_histX, Diff_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_x_histX, 0, Diff_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_x, 0, max(Diff_Start_x_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_Start_x)
        mp.axvspan(mean_Start_x - std_Start_x, mean_Start_x + std_Start_x, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_Start_x)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('Start x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/Difference_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_y, 0, max(Diff_Start_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_Start_y)
        mp.axvspan(mean_Start_y - std_Start_y, mean_Start_y + std_Start_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_Start_y)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('Start y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/Difference_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_z, 0, max(Diff_Start_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_Start_z)
        mp.axvspan(mean_Start_z - std_Start_z, mean_Start_z + std_Start_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_Start_z)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('Start z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/Difference_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()    
    
    if plot_End:    
        mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_x, 0, max(Diff_End_x_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_End_x)
        mp.axvspan(mean_End_x - std_End_x, mean_End_x + std_End_x, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_End_x)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('End x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_X_hist.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_X_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_y, 0, max(Diff_End_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_End_y)
        mp.axvspan(mean_End_y - std_End_y, mean_End_y + std_End_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_End_y)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('End y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_z, 0, max(Diff_End_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f mm)' % mean_End_z)
        mp.axvspan(mean_End_z - std_End_z, mean_End_z + std_End_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f mm)' % std_End_z)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [mm]')
        mp.ylabel('#')
        mp.title('End z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
       
    ### Charge ID
    if plot_Charge:
        # Filter out all events without truth +/-13 as PDG
        boolean_true_muon_stop = (np.abs(True_Charge_stop) == 13)
        True_Charge_stop = True_Charge_stop[boolean_true_muon_stop]
        True_KE_stop = True_KE_stop[boolean_true_muon_stop]
        Reco_Charge = Reco_Charge[boolean_true_muon_stop]
        boolean_issues = (Reco_Charge == -9.99999999e+08)
        not_identified = Reco_Charge[boolean_issues]
        
        # Prepare the KE arrays
        Muons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        AMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
    
        # Counter for both positive, both negative, different reco/truth positive and negative
        true_muons = len(True_Charge[True_Charge_ctop == 13])
        true_antimuons = len(True_Charge[True_Charge_stop == -13])
        reco_muons = len(Reco_Charge[Reco_Charge == 13])
        reco_antimuons = len(Reco_Charge[Reco_Charge == -13])
        counter_true_positive = 0   #truth and reco agree on muon
        counter_true_negative = 0   #truth and reco agree on antimuon
        counter_false_positive = 0  #truth and reco disagree, reco -> muon
        counter_false_negative = 0  #truth and reco disagree, reco -> antimuon
    
        True_Muons_KE = True_KE_stop[True_Charge_stop == 13]
        True_Antimuons_KE = True_KE_stop[True_Charge_stop == -13]
    
        # Compare sign of truth and reco    
        for i in range(len(True_Charge_stop)):
            if True_Charge[i] == Reco_Charge[i]:
                if True_Charge[i] > 0:
                    counter_true_positive += 1
                    Muons_KE[i] = True_KE_stop[i]
                else:
                    counter_true_negative += 1
                    AMuons_KE[i] = True_KE_stop[i]
            else:
                if Reco_Charge[i] > 0:
                    counter_false_positive += 1
                else:
                    counter_false_negative += 1
        
        if (counter_true_positive + counter_true_negative + counter_false_positive + counter_false_negative) != len(True_Charge_stop):
            print("Counters don't add up")
            print("Length: ", len(True_Charge_stop))
        
        print("Charge ID numbers")
        print("  Not identified:  ", len(not_identified))
        print("  True muons:      ", counter_true_positive)
        print("  True antimuons:  ", counter_true_negative)
        print("  False muons:     ", counter_false_positive)
        print("  False antimuons: ", counter_false_negative, " here the not identified go into!")
        
        # Performance evaluation
        if counter_true_positive > 0:
            efficiency_muons = counter_true_positive / (counter_true_positive + (counter_false_negative - len(not_identified)))
            purity_muons = counter_true_positive / (counter_true_positive + counter_false_positive)
            print("  Muons (efficiency | purity):     ", efficiency_muons, " | ", purity_muons)
        if counter_true_negative > 0:
            efficiency_antimuons = counter_true_negative / (counter_true_negative + counter_false_positive)
            purity_antimuons = counter_true_negative / (counter_true_negative + (counter_false_negative - len(not_identified)))
            print("  Antimuons (efficiency | purity): ", efficiency_antimuons, " | ", purity_antimuons)

        accuracy_anti_muons = (counter_true_positive + counter_true_negative) / (counter_true_positive + counter_true_negative + counter_false_positive + (counter_false_negative - len(not_identified)))    
        print("  Accuracy (both): ", accuracy_anti_muons)

        # Energy dependent      
        Muons_KE = Muons_KE[Muons_KE != -9999.]
        AMuons_KE = AMuons_KE[AMuons_KE != -9999.]
        
        muons_ke_hist, muons_ke_bins = np.histogram(Muons_KE, bins = 50, range = (0, 5000))
        amuons_ke_hist, amuons_ke_bins = np.histogram(AMuons_KE, bins = muons_ke_bins)
    
        true_muons_ke_hist, muons_ke_bins = np.histogram(True_Muons_KE, bins = muons_ke_bins)
        true_antimuons_ke_hist, amuons_ke_bins = np.histogram(True_Antimuons_KE, bins = amuons_ke_bins)
    
        muons_ke_hist_x, muons_ke_hist_y = histogram_arr_handle(muons_ke_hist, muons_ke_bins)
        amuons_ke_hist_x, amuons_ke_hist_y = histogram_arr_handle(amuons_ke_hist, amuons_ke_bins)
        true_muons_ke_hist_x, true_muons_ke_hist_y = histogram_arr_handle(true_muons_ke_hist, muons_ke_bins)
        true_antimuons_ke_hist_x, true_antimuons_ke_hist_y = histogram_arr_handle(true_antimuons_ke_hist, amuons_ke_bins)
    
        if counter_true_positive > 0:
            mp.plot(muons_ke_hist_x, muons_ke_hist_y / true_muons_ke_hist_y, color = blue_cbf, linestyle = '-', linewidth = 2)
            mp.fill_between(muons_ke_hist_x, 0, muons_ke_hist_y / true_muons_ke_hist_y, color = blue_cbf, alpha = 0.3, hatch = '//', label = '$\\mu$ correct')
        if counter_true_negative > 0:
            mp.plot(amuons_ke_hist_x, amuons_ke_hist_y / true_antimuons_ke_hist_y, color = orange_cbf, linestyle = '--', linewidth = 2)
            mp.fill_between(amuons_ke_hist_x, 0, amuons_ke_hist_y / true_antimuons_ke_hist_y, color = orange_cbf, alpha = 0.3, hatch = '\\\\', label = '$\\bar{\\mu}$ correct')
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Fraction')
        mp.legend(loc = 'lower center', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Charge_ID_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
    ### Angular resolution
    if plot_Angle:
        # Filter out all events without truth +/-13 as PDG    
        boolean_true_muon = (np.abs(True_Charge) == 13)
        True_KE = True_KE[boolean_true_muon]
        True_TrackDirection_xz_KE = True_TrackDirection_xz[boolean_true_muon]
        True_TrackDirection_yz_KE = True_TrackDirection_yz[boolean_true_muon]
        Reco_TrackDirection_xz_KE = Reco_TrackDirection_xz[boolean_true_muon]
        Reco_TrackDirection_yz_KE = Reco_TrackDirection_yz[boolean_true_muon]
    
        # calculate differences in angles
        xz_difference = True_TrackDirection_xz - Reco_TrackDirection_xz
        yz_difference = True_TrackDirection_yz - Reco_TrackDirection_yz
        
        xz_difference = xz_difference[~np.isnan(xz_difference)] #filter out potential nan's
        yz_difference = yz_difference[~np.isnan(yz_difference)] #filter out potential nan's
        
        # calculate mean and std of differences
        mean_xz_difference = np.mean(xz_difference)
        mean_yz_difference = np.mean(yz_difference)
        std_xz_difference = np.std(xz_difference)
        std_yz_difference = np.std(yz_difference)
        
        print("XZ angle difference: ", mean_xz_difference, " ± ", std_xz_difference)
        print("YZ angle difference: ", mean_yz_difference, " ± ", std_yz_difference)
        
        # calculate histograms of differences
        Diff_xz_hist, Diff_xz_bins = np.histogram(xz_difference, bins = 50)
        Diff_yz_hist, Diff_yz_bins = np.histogram(yz_difference, bins = 50)
        
        Diff_xz_histX, Diff_xz_histY = histogram_arr_handle(Diff_xz_hist, Diff_xz_bins)
        Diff_yz_histX, Diff_yz_histY = histogram_arr_handle(Diff_yz_hist, Diff_yz_bins)
        
        # plot
        mp.plot(Diff_xz_histX, Diff_xz_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(Diff_xz_histX, 0, Diff_xz_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_xz_difference, 0, max(Diff_xz_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f°)' % mean_xz_difference)
        mp.axvspan(mean_xz_difference - std_xz_difference, mean_xz_difference + std_xz_difference, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f°)' % std_xz_difference)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('#')
        mp.title('xz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/Difference_XZ_angle_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        mp.plot(Diff_yz_histX, Diff_yz_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(Diff_yz_histX, 0, Diff_yz_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_yz_difference, 0, max(Diff_yz_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f°)' % mean_yz_difference)
        mp.axvspan(mean_yz_difference - std_yz_difference, mean_yz_difference + std_yz_difference, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f°)' % std_yz_difference)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('#')
        mp.title('yz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/Difference_YZ_angle_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        # energy dependent
        bins = np.histogram_bin_edges(True_KE, bins = 50, range = (0, 5000))
        
        xz_mean_KE = np.zeros(len(bins), dtype = float)
        yz_mean_KE = np.zeros(len(bins), dtype = float)
        xz_std_KE = np.zeros(len(bins), dtype = float)
        yz_std_KE = np.zeros(len(bins), dtype = float)
        
        for i in range(len(bins) - 1):
            boolean_KE = (True_KE >= bins[i]) & (True_KE < bins[i+1])
            True_xz = True_TrackDirection_xz_KE[boolean_KE]
            True_yz = True_TrackDirection_yz_KE[boolean_KE]
            Reco_xz = Reco_TrackDirection_xz_KE[boolean_KE]
            Reco_yz = Reco_TrackDirection_yz_KE[boolean_KE]
            
            xz_diff = True_xz - Reco_xz
            yz_diff = True_yz - Reco_yz
            
            xz_mean_KE[i] = np.mean(xz_diff)
            xz_std_KE[i] = np.std(xz_diff)
            yz_mean_KE[i] = np.mean(yz_diff)
            yz_std_KE[i] = np.std(yz_diff)
    
        # plot
        mp.errorbar(bins, xz_mean_KE, yerr = xz_std_KE, linestyle = '--', linewidth = 1.5, color = red_cbf)
        mp.scatter(bins, xz_mean_KE, marker = '.', color = red_cbf, label = 'mean')
        mp.fill_between(bins, xz_mean_KE - xz_std_KE, xz_mean_KE + xz_std_KE, color = red_cbf, alpha = 0.1, label = '1 $\\sigma$')
        mp.hlines(0, 0, 5000, color = black_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('xz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_XZ_angle_True_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        mp.errorbar(bins, yz_mean_KE, yerr = yz_std_KE, linestyle = '--', linewidth = 1.5, color = blue_cbf)
        mp.scatter(bins, yz_mean_KE, marker = '.', color = blue_cbf, label = 'mean')
        mp.fill_between(bins, yz_mean_KE - yz_std_KE, yz_mean_KE + yz_std_KE, color = blue_cbf, alpha = 0.1, label = '1 $\\sigma$')
        mp.hlines(0, 0, 5000, color = black_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('yz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_YZ_angle_True_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
    return

def dir_to_angle(xy, z):
    return np.degrees(np.arctan(xy / z))

def histogram_arr_handle(bins, edges):#(u bins,v edges)
    """Function to calculated x- and y-values from np.histogram

	@param[in] bins: bin-values from np.histogram to be transformed into y-values
	@param[in] edges: edge-values from np.histogram to be transformed into x-values

	@return: 1d-array containing the x-values and 1d-array containing the y-values
	"""
    x_output = np.zeros(2 * len(bins))
    y_output = np.zeros(2 * len(bins))
    #v edges to u*2 x-values
    x_output[0] = edges[0]
    counter = 1
    for i in range(1, len(x_output)-1, 2):
        x_output[i] = edges[counter]
        x_output[i + 1] = edges[counter]
        counter = counter + 1
    x_output[-1] = edges[-1]
	
    #u bins to u*2 y-values
    counter = 0
    for i in range(0, len(y_output), 2):
        y_output[i] = bins[counter]
        y_output[i + 1] = bins[counter]
        counter = counter + 1

    return x_output, y_output

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--Start', "-s", help = "Do you want to plot the Start point resolution? Yes -> --Start, No -> --no-Start", action = argparse.BooleanOptionalAction)
    parser.add_argument('--End', "-e", help = "Do you want to plot the End point resolution? Yes -> --End, No -> --no-End", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Charge', "-c", help = "Do you want to plot the Charge identification efficiency plots? Yes -> --Charge, No -> --no-Charge", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Angle', "-a", help = "Do you want to plot the Angular resolution? Yes -> --Angle, No -> --no-Angle", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Contained', "-con", help = "Do you want to plot the end resolution only for contained muons or all? Yes -> --Contained, No -> --no-Contained", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    #print(layer_dict)
    out_dir = args.outdir
    input_filename = args.input_filename
    plot_Start = args.Start
    plot_End = args.End
    plot_Charge = args.Charge
    plot_Angle = args.Angle
    Contained = args.Contained
    draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained)
