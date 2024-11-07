import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll
import matplotlib.cm as cm
from matplotlib_venn import venn2, venn2_circles

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-poster')

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
def draw_performance(out_dir, input_filename, Xlayers):
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
    
    Reco_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
    Reco_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
    Primary_True_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
    Primary_True_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
    True_TrackLength = np.ones((n_events, 5), dtype = float) * -9999.
    Reco_TrackLength = np.ones((n_events, 5), dtype = float) * -9999.
    Reco_Charge = np.ones((n_events, 5), dtype = float) * -9999.
    True_Charge = np.ones((n_events, 5), dtype = float) * -9999.
    True_KE = np.ones((n_events, 5), dtype = float) * -9999.
    True_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true
    Reco_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true
    
    correct_tracks_reco = 0
    correct_reco_hits = 0
    correct_true_hits = 0
    
    count_muons = 0
    
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

             
            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
            Muon_TrueTrackLength = true_event.Muon_TrueTrackLength
            Reco_Track_Length = np.frombuffer(event.Length, dtype = np.float32)
            #Reco_Track_Charge = event.Charge
            MomentumTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumTrackStart, dtype = np.float32)
            True_PDG = true_event.PDG
            LArFiducialTouch = true_event.RecoTrackPrimaryParticleLArFiducialStart
            True_Position_TMS_Start = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionEnteringTMS, dtype = np.float32)
            True_Position_TMS_End = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionLeavingTMS, dtype = np.float32)
            
            Reco_Hits = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype = np.float32)
            sum_reco_hits = event.nHits
            sum_true_hits = true_event.RecoTrackNHits
            
            Particle_PDG = true_event.LeptonPDG
            Muon_Start = np.frombuffer(true_event.Muon_Vertex, dtype = np.float32)
            Muon_End = np.frombuffer(true_event.Muon_Death, dtype = np.float32)
            
            if (abs(Particle_PDG) == 13):
                if 4179.24 < Muon_Start[2] < 9135.88 and 11185 < Muon_End[2] < 18535:    #12462
                    if 1159 > Muon_End[1] > -3864 and abs(Muon_End[0]) < 3520:
                        count_muons += 1
            
            nTracks = event.nTracks
            if nTracks <= 0: continue
            if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break
                
                # check if (anti-)muon as true primary particle and if origin in LAr
                if np.abs(True_PDG[j]) == 13 and LArFiducialTouch[j]:
                    # check if (anti-)muon travesers at least 4 planes in TMS
                    if (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) >= 440.:
                        # if so, then this is a true muon: 0 -> 1
                        True_Muon_Track[i, j] = 1.
                        counter_correct = 0
                        for true_hits in range(int(sum_true_hits[j])):
                            for reco_hits in range(int(sum_reco_hits[j])):
                                if True_Hits[j*600 + true_hits*4 + 2] == Reco_Hits[j*600 + reco_hits*3 + 2]:
                                    if np.abs(True_Hits[j*600 + true_hits*4 + 0] - Reco_Hits[j*600 + reco_hits*3 + 0]) <= 2 * 35.42:
                                        counter_correct += 1

                        if counter_correct / sum_true_hits[j] >= 0.1:
                            correct_tracks_reco += 1
                            correct_reco_hits += counter_correct
                            correct_true_hits += sum_true_hits[j]
                        
                        # check if reconstructed tracks exist for this event
                        if StartPos.size != 0:
                            # if so, then add all identified as muons: 0 -> number tracks
                            Reco_Muon_Track[i, j] = 1. #StartPos.size / 3
                
                if RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] > -8000. and not StartPos.size == 0:
                    # checking for muon tracks (minimal length for this are 20 planes traversed -> 890 mm in thin area
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:
                        Reco_Start[i, j, 0] = StartPos[j*3 + 0]
                        Reco_Start[i, j, 1] = StartPos[j*3 + 1]
                        Reco_Start[i, j, 2] = StartPos[j*3 + 2]
                        Primary_True_Start[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0]
                        Primary_True_Start[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1]
                        Primary_True_Start[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]
                        Primary_True_px = MomentumTrackStart[j*4 + 0]
                        Primary_True_py = MomentumTrackStart[j*4 + 1]
                        Primary_True_pz = MomentumTrackStart[j*4 + 2]
                        True_KE[i, j] = np.sqrt((Primary_True_px**2 + Primary_True_py**2 + Primary_True_pz**2 + 105.7**2) - 105.7)
                        True_TrackLength[i, j] = Muon_TrueTrackLength
                        Reco_TrackLength[i, j] = Reco_Track_Length[j]
                        True_Charge[i, j] = True_PDG[j]
                        Reco_Charge[i, j] = Reco_Track_Charge[j]
                        
                if RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] > -8000. and not EndPos.size == 0:
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:            
                        Reco_End[i, j, 0] = EndPos[j*3 + 0]
                        Reco_End[i, j, 1] = EndPos[j*3 + 1]
                        Reco_End[i, j, 2] = EndPos[j*3 + 2]
                        Primary_True_End[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0]
                        Primary_True_End[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1]
                        Primary_True_End[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2]
    
    # filter out not filled indice
    boolean_Reco_Start = (Reco_Start[:, :, 0] != -9999.)
    Reco_Start = Reco_Start[boolean_Reco_Start]
    boolean_Reco_End = (Reco_End[:, :, 0] != -9999.)
    Reco_End = Reco_End[boolean_Reco_End]
    boolean_Primary_Start = (Primary_True_Start[:, :, 0] != -9999.)
    Primary_True_Start = Primary_True_Start[boolean_Primary_Start]
    boolean_Primary_End = (Primary_True_End[:, :, 0] != -9999.)
    Primary_True_End = Primary_True_End[boolean_Primary_End]
    boolean_True_TrackLength = (True_TrackLength != -9999.)
    True_TrackLength = True_TrackLength[boolean_True_TrackLength]
    boolean_Reco_TrackLength = (Reco_TrackLength != -9999.)
    Reco_TrackLength = Reco_TrackLength[boolean_Reco_TrackLength]
    boolean_Reco_Charge = (Reco_Charge != -9999.)
    Reco_Charge = Reco_Charge[boolean_Reco_Charge]
    boolean_True_Charge = (True_Charge != -9999.)
    True_Charge = True_Charge[boolean_True_Charge]
    boolean_True_KE = (True_KE != -9999.)
    True_KE = True_KE[boolean_True_KE]
    
    # flatten arrays
    Reco_Start_x = Reco_Start[:, 0]
    Reco_Start_y = Reco_Start[:, 1]
    Reco_Start_z = Reco_Start[:, 2]
    Reco_End_x = Reco_End[:, 0]
    Reco_End_y = Reco_End[:, 1]
    Reco_End_z = Reco_End[:, 2]
    Primary_True_Start_x = Primary_True_Start[:, 0]
    Primary_True_Start_y = Primary_True_Start[:, 1]
    Primary_True_Start_z = Primary_True_Start[:, 2]
    Primary_True_End_x = Primary_True_End[:, 0]
    Primary_True_End_y = Primary_True_End[:, 1]
    Primary_True_End_z = Primary_True_End[:, 2]
    
    # total number of events after filtering
    print("#events reconstruction: ", len(Reco_Start), "# events truth: ", len(Primary_True_Start))
    print("tracklength truth: ", len(True_TrackLength), "reco: ", len(Reco_TrackLength))
    print("true (anti-)muons: ", sum(True_Muon_Track), "  vs. reconstructed particles: ", sum(Reco_Muon_Track))
    print("correctly identified tracks: ", correct_tracks_reco)
    print("  correct hits reco: ", correct_reco_hits, " vs. true hits: ", correct_true_hits, " -> ", correct_reco_hits / correct_true_hits * 100)
    print("total muons (outside of reco): ", count_muons)
    
    # subtract reconstruction from truth for all directions
    Diff_Start_x = Primary_True_Start_x - Reco_Start_x
    Diff_Start_y = Primary_True_Start_y - Reco_Start_y
    Diff_Start_z = Primary_True_Start_z - Reco_Start_z
    Diff_End_x = Primary_True_End_x - Reco_End_x
    Diff_End_y = Primary_True_End_y - Reco_End_y
    Diff_End_z = Primary_True_End_z - Reco_End_z

    # create histograms for the differences
    Diff_Start_x_hist, Diff_Start_x_bins = np.histogram(Diff_Start_x, bins = 50)
    Diff_Start_y_hist, Diff_Start_y_bins = np.histogram(Diff_Start_y, bins = 50)
    Diff_Start_z_hist, Diff_Start_z_bins = np.histogram(Diff_Start_z, bins = 50)
    Diff_End_x_hist, Diff_End_x_bins = np.histogram(Diff_End_x, bins = 50)
    Diff_End_y_hist, Diff_End_y_bins = np.histogram(Diff_End_y, bins = 50)
    Diff_End_z_hist, Diff_End_z_bins = np.histogram(Diff_End_z, bins = 50)
    
    Diff_Start_x_histX, Diff_Start_x_histY = histogram_arr_handle(Diff_Start_x_hist, Diff_Start_x_bins)
    Diff_Start_y_histX, Diff_Start_y_histY = histogram_arr_handle(Diff_Start_y_hist, Diff_Start_y_bins)
    Diff_Start_z_histX, Diff_Start_z_histY = histogram_arr_handle(Diff_Start_z_hist, Diff_Start_z_bins)
    Diff_End_x_histX, Diff_End_x_histY = histogram_arr_handle(Diff_End_x_hist, Diff_End_x_bins)
    Diff_End_y_histX, Diff_End_y_histY = histogram_arr_handle(Diff_End_y_hist, Diff_End_y_bins)
    Diff_End_z_histX, Diff_End_z_histY = histogram_arr_handle(Diff_End_z_hist, Diff_End_z_bins)
    
    mean_Start_x = np.mean(Diff_Start_x)
    mean_Start_y = np.mean(Diff_Start_y)
    mean_Start_z = np.mean(Diff_Start_z)
    mean_End_x = np.mean(Diff_End_x)
    mean_End_y = np.mean(Diff_End_y)
    mean_End_z = np.mean(Diff_End_z)
    std_Start_x = np.std(Diff_Start_x)
    std_Start_y = np.std(Diff_Start_y)
    std_Start_z = np.std(Diff_Start_z)
    std_End_x = np.std(Diff_End_x)
    std_End_y = np.std(Diff_End_y)
    std_End_z = np.std(Diff_End_z)
    
    print("Start x: ", mean_Start_x, "±", std_Start_x, "    y: ", mean_Start_y, "±", std_Start_y, "    z: ", mean_Start_z, "±", std_Start_z)
    print("End   x: ", mean_End_x, "±", std_End_x, "    y: ", mean_End_y, "±", std_End_y, "    z: ", mean_End_z, "±", std_End_z)
    
    # plot
    mp.plot(Diff_Start_x_histX, Diff_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_Start_x_histX, 0, Diff_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_Start_x, 0, max(Diff_Start_x_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_Start_x - std_Start_x, mean_Start_x + std_Start_x, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('Start x', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_Start_y, 0, max(Diff_Start_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_Start_y - std_Start_y, mean_Start_y + std_Start_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('Start y', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_Start_z, 0, max(Diff_Start_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_Start_z - std_Start_z, mean_Start_z + std_Start_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('Start z', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_x, 0, max(Diff_End_x_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_x - std_End_x, mean_End_x + std_End_x, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End x', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_y, 0, max(Diff_End_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_y - std_End_y, mean_End_y + std_End_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End y', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_z, 0, max(Diff_End_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_z - std_End_z, mean_End_z + std_End_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End z', fontsize = 'xx-large')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
        
    # create histograms for all directions for truth and reconstruction 
    z_bins = np.zeros(93, dtype = float)
    z_bins[0] = 11185
    for i in range(1, 93):
        if i < 50:
            z_bins[i] = z_bins[0] + i * 65.
        if i >= 50:
            z_bins[i] = 14435. + (i - 50) * 90.
        if i >= 84:
            z_bins[i] = 17495. + (i - 84) * 130.
    
    Reco_Start_x_hist, Reco_Start_x_bins = np.histogram(Reco_Start_x, bins = 48, range = (-3520., 3520.))
    Reco_Start_y_hist, Reco_Start_y_bins = np.histogram(Reco_Start_y, bins = 20, range = (-2949., 244.)) 
    Reco_Start_z_hist, Reco_Start_z_bins = np.histogram(Reco_Start_z, bins = z_bins)
    Reco_End_x_hist, Reco_End_x_bins = np.histogram(Reco_End_x, bins = Reco_Start_x_bins)
    Reco_End_y_hist, Reco_End_y_bins = np.histogram(Reco_End_y, bins = Reco_Start_y_bins)
    Reco_End_z_hist, Reco_End_z_bins = np.histogram(Reco_End_z, bins = Reco_Start_z_bins)
    Primary_True_Start_x_hist, Primary_True_Start_x_bins = np.histogram(Primary_True_Start_x, bins = Reco_Start_x_bins)
    Primary_True_Start_y_hist, Primary_True_Start_y_bins = np.histogram(Primary_True_Start_y, bins = Reco_Start_y_bins)
    Primary_True_Start_z_hist, Primary_True_Start_z_bins = np.histogram(Primary_True_Start_z, bins = Reco_Start_z_bins)
    Primary_True_End_x_hist, Primary_True_End_x_bins = np.histogram(Primary_True_End_x, bins = Reco_Start_x_bins)
    Primary_True_End_y_hist, Primary_True_End_y_bins = np.histogram(Primary_True_End_y, bins = Reco_Start_y_bins)
    Primary_True_End_z_hist, Primary_True_End_z_bins = np.histogram(Primary_True_End_z, bins = Reco_Start_z_bins)
    
    # filter histograms for tracks with and without X hits for Xlayers == true
    if (Xlayers):
        start_helper = np.array([check_orientation(int(Reco_Start_z[i])) == 'kXBar' for i in range(len(Reco_Start_z))], dtype = bool)
        Reco_Start_x_WX = Reco_Start_x[start_helper]
        Reco_Start_y_WX = Reco_Start_y[start_helper]
        Reco_Start_z_WX = Reco_Start_z[start_helper]
        end_helper = np.array([check_orientation(int(Reco_End_z[i])) == 'kXBar' for i in range(len(Reco_End_z))], dtype = bool)
        Reco_End_x_WX = Reco_End_x[end_helper]
        Reco_End_y_WX = Reco_End_y[end_helper]
        Reco_End_z_WX = Reco_End_z[end_helper]
        
        Reco_Start_x_WX_hist, Reco_Start_x_bins = np.histogram(Reco_Start_x_WX, bins = Reco_Start_x_bins)
        Reco_Start_y_WX_hist, Reco_Start_y_bins = np.histogram(Reco_Start_y_WX, bins = Reco_Start_y_bins)
        Reco_Start_z_WX_hist, Reco_Start_z_bins = np.histogram(Reco_Start_z_WX, bins = Reco_Start_z_bins)
        
        Reco_End_x_WX_hist, Reco_End_x_bins = np.histogram(Reco_End_x_WX, bins = Reco_End_x_bins)
        Reco_End_y_WX_hist, Reco_End_y_bins = np.histogram(Reco_End_y_WX, bins = Reco_End_y_bins)
        Reco_End_z_WX_hist, Reco_End_z_bins = np.histogram(Reco_End_z_WX, bins = Reco_End_z_bins)
        
        Reco_Start_x_WX_histX, Reco_Start_x_WX_histY = histogram_arr_handle(Reco_Start_x_WX_hist, Reco_Start_x_bins)
        Reco_Start_y_WX_histX, Reco_Start_y_WX_histY = histogram_arr_handle(Reco_Start_y_WX_hist, Reco_Start_y_bins)
        Reco_Start_z_WX_histX, Reco_Start_z_WX_histY = histogram_arr_handle(Reco_Start_z_WX_hist, Reco_Start_z_bins)
        
        Reco_End_x_WX_histX, Reco_End_x_WX_histY = histogram_arr_handle(Reco_End_x_WX_hist, Reco_End_x_bins)
        Reco_End_y_WX_histX, Reco_End_y_WX_histY = histogram_arr_handle(Reco_End_y_WX_hist, Reco_End_y_bins)
        Reco_End_z_WX_histX, Reco_End_z_WX_histY = histogram_arr_handle(Reco_End_z_WX_hist, Reco_End_z_bins)
    
    # make the histograms usable
    Reco_Start_x_histX, Reco_Start_x_histY = histogram_arr_handle(Reco_Start_x_hist, Reco_Start_x_bins)
    Reco_Start_y_histX, Reco_Start_y_histY = histogram_arr_handle(Reco_Start_y_hist, Reco_Start_y_bins)
    Reco_Start_z_histX, Reco_Start_z_histY = histogram_arr_handle(Reco_Start_z_hist, Reco_Start_z_bins)
    Reco_End_x_histX, Reco_End_x_histY = histogram_arr_handle(Reco_End_x_hist, Reco_End_x_bins)
    Reco_End_y_histX, Reco_End_y_histY = histogram_arr_handle(Reco_End_y_hist, Reco_End_y_bins)
    Reco_End_z_histX, Reco_End_z_histY = histogram_arr_handle(Reco_End_z_hist, Reco_End_z_bins)
    Primary_True_Start_x_histX, Primary_True_Start_x_histY = histogram_arr_handle(Primary_True_Start_x_hist, Primary_True_Start_x_bins)
    Primary_True_Start_y_histX, Primary_True_Start_y_histY = histogram_arr_handle(Primary_True_Start_y_hist, Primary_True_Start_y_bins)
    Primary_True_Start_z_histX, Primary_True_Start_z_histY = histogram_arr_handle(Primary_True_Start_z_hist, Primary_True_Start_z_bins)
    Primary_True_End_x_histX, Primary_True_End_x_histY = histogram_arr_handle(Primary_True_End_x_hist, Primary_True_End_x_bins)
    Primary_True_End_y_histX, Primary_True_End_y_histY = histogram_arr_handle(Primary_True_End_y_hist, Primary_True_End_y_bins)
    Primary_True_End_z_histX, Primary_True_End_z_histY = histogram_arr_handle(Primary_True_End_z_hist, Primary_True_End_z_bins)        

    # now plot this
    mp.plot(Primary_True_Start_x_histX, Primary_True_Start_x_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_x_histX, Reco_Start_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_Start_x_WX_histX, Reco_Start_x_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')
    mp.fill_between(Primary_True_Start_x_histX, 0, Primary_True_Start_x_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_x_histX, 0, Reco_Start_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start x', fontsize = 'xx-large')
    mp.savefig('%s/Track_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_y_histX, Primary_True_Start_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_y_histX, Reco_Start_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_Start_y_WX_histX, Reco_Start_y_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')    
    mp.fill_between(Primary_True_Start_y_histX, 0, Primary_True_Start_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_y_histX, 0, Reco_Start_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start y', fontsize = 'xx-large')   
    mp.savefig('%s/Track_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_z_histX, Primary_True_Start_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_z_histX, Reco_Start_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_Start_z_WX_histX, Reco_Start_z_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')
    mp.fill_between(Primary_True_Start_z_histX, 0, Primary_True_Start_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_z_histX, 0, Reco_Start_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start z', fontsize = 'xx-large')
    mp.savefig('%s/Track_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_x_histX, Primary_True_End_x_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_x_histX, Reco_End_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_End_x_WX_histX, Reco_End_x_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')
    mp.fill_between(Primary_True_End_x_histX, 0, Primary_True_End_x_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_x_histX, 0, Reco_End_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    if (Xlayers): mp.fill_between(Reco_End_x_WX_histX, 0, Reco_End_x_WX_histY, color = black_cbf, alpha = 0.2)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End x', fontsize = 'xx-large')   
    mp.savefig('%s/Track_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_y_histX, Primary_True_End_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_y_histX, Reco_End_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_End_y_WX_histX, Reco_End_y_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')
    mp.fill_between(Primary_True_End_y_histX, 0, Primary_True_End_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_y_histX, 0, Reco_End_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    if (Xlayers): mp.fill_between(Reco_End_y_WX_histX, 0, Reco_End_y_WX_histY, color = black_cbf, alpha = 0.2)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End y', fontsize = 'xx-large')
    mp.savefig('%s/Track_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_z_histX, Primary_True_End_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_z_histX, Reco_End_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    if (Xlayers): mp.plot(Reco_End_z_WX_histX, Reco_End_z_WX_histY, color = black_cbf, linestyle = '-', linewidth = 2, label = 'just X (reco)')
    mp.fill_between(Primary_True_End_z_histX, 0, Primary_True_End_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_z_histX, 0, Reco_End_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    if (Xlayers): mp.fill_between(Reco_End_z_WX_histX, 0, Reco_End_z_WX_histY, color = black_cbf, alpha = 0.2)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End z', fontsize = 'xx-large')
    mp.savefig('%s/Track_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    ### plot difference in dependence of hit position
    # create 2d histograms
    
    dependence_Start_x_hist, dependence_Start_x_binsX, dependence_Start_x_binsY = np.histogram2d(Primary_True_Start_x, Diff_Start_x, bins = [Primary_True_Start_x_bins, Diff_Start_x_bins])
    dependence_Start_y_hist, dependence_Start_y_binsX, dependence_Start_y_binsY = np.histogram2d(Primary_True_Start_y, Diff_Start_y, bins = [Primary_True_Start_y_bins, Diff_Start_y_bins])
    dependence_Start_z_hist, dependence_Start_z_binsX, dependence_Start_z_binsY = np.histogram2d(Primary_True_Start_z, Diff_Start_z, bins = [Primary_True_Start_z_bins, Diff_Start_z_bins])
    dependence_End_x_hist, dependence_End_x_binsX, dependence_End_x_binsY = np.histogram2d(Primary_True_End_x, Diff_End_x, bins = [Primary_True_End_x_bins, Diff_End_x_bins])
    dependence_End_y_hist, dependence_End_y_binsX, dependence_End_y_binsY = np.histogram2d(Primary_True_End_y, Diff_End_y, bins = [Primary_True_End_y_bins, Diff_End_y_bins])
    dependence_End_z_hist, dependence_End_z_binsX, dependence_End_z_binsY = np.histogram2d(Primary_True_End_z, Diff_End_z, bins = [Primary_True_End_z_bins, Diff_End_z_bins])
    
    cmap = cm.get_cmap('cividis');
    
    im = mp.pcolormesh(dependence_Start_x_binsX, dependence_Start_x_binsY, np.transpose(dependence_Start_x_hist), cmap = cmap);
    mp.xlabel('Start in TMS (X) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start x', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_X_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_y_binsX, dependence_Start_y_binsY, np.transpose(dependence_Start_y_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start y', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_Y_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_z_binsX, dependence_Start_z_binsY, np.transpose(dependence_Start_z_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start z', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_Z_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_x_binsX, dependence_End_x_binsY, np.transpose(dependence_End_x_hist), cmap = cmap);
    mp.xlabel('End in TMS (X) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End x', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/End_X_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();
    
    im = mp.pcolormesh(dependence_End_y_binsX, dependence_End_y_binsY, np.transpose(dependence_End_y_hist), cmap = cmap);
    mp.xlabel('End in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End y', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/End_Y_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_z_binsX, dependence_End_z_binsY, np.transpose(dependence_End_z_hist), cmap = cmap);
    mp.xlabel('End in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End z', fontsize = 'xx-large')
    mp.colorbar(im);
    mp.savefig('%s/End_Z_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    # Track length
    # histogram
    
    # single histogram as well
    reco_tracklength_hist, reco_tracklength_bins = np.histogram(True_TrackLength - Reco_TrackLength / 10, bins = 40)#, range = (min(min(True_TrackLength - Reco_TrackLength / 10), min(True_TrackLength)), max(max(True_TrackLength - Reco_TrackLength / 10), max(True_TrackLength))))#(50, 5000))#
    true_tracklength_hist, true_tracklength_bins = np.histogram(True_TrackLength, bins = 40)#reco_tracklength_bins)
    
    #true_tracklength_0 = True_TrackLength[True_TrackLength == 0]
    #print(len(true_tracklength_0))
    
    reco_tracklength_histX, reco_tracklength_histY = histogram_arr_handle(reco_tracklength_hist, reco_tracklength_bins)
    true_tracklength_histX, true_tracklength_histY = histogram_arr_handle(true_tracklength_hist, true_tracklength_bins)
    
    mp.plot(true_tracklength_histX, true_tracklength_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(reco_tracklength_histX, reco_tracklength_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(true_tracklength_histX, 0, true_tracklength_histY, alpha = 0.6, color = blue_cbf, hatch = '\\\\')
    mp.fill_between(reco_tracklength_histX, 0, reco_tracklength_histY, alpha = 0.6, color = orange_cbf, hatch = '//')
    mp.xlabel('True track length - Reco track length [$\\frac{g}{cm^2}$]')
    mp.ylabel('#')
    mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.savefig('%s/Reco_tracklength.png' % out_dir, bbox_inches = 'tight')
    mp.close()    
    
    tracklength_hist, tracklength_binsX, tracklength_binsY = np.histogram2d(True_TrackLength, True_TrackLength - Reco_TrackLength / 10, bins = [true_tracklength_bins, reco_tracklength_bins])
    
    cmap = cm.get_cmap('cividis');
    
    im = mp.pcolormesh(tracklength_binsX, tracklength_binsY, np.transpose(tracklength_hist), cmap = cmap)
    mp.xlabel('True track length [$\\frac{g}{cm^2}$]')
    mp.ylabel('True track length - Reco track length [$\\frac{g}{cm^2}$]')
    mp.colorbar(im)
    mp.savefig('%s/Tracklength.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    ### single slices for comparison of geometries
    # low energy slice
    # find all events with energy between lower and upper bound
    combined_TrackLength = True_TrackLength - Reco_TrackLength / 10
    low_helper = (True_TrackLength > 150.) & (True_TrackLength < 250.)
    high_helper = (True_TrackLength > 400.)
    low_energy_slice = combined_TrackLength[low_helper]
    high_energy_slice = combined_TrackLength[high_helper]    
    
    # calculate mean and std of this slice
    low_energy_mean = np.mean(low_energy_slice)
    low_energy_std = np.std(low_energy_slice)
    print('Low energy: ', low_energy_mean, ' ± ', low_energy_std)
    
    high_energy_mean = np.mean(high_energy_slice)
    high_energy_std = np.std(high_energy_slice)
    print('High energy: ', high_energy_mean, ' ± ', high_energy_std)
    
    # plot slice distribution
    low_energy_slice_hist, low_energy_slice_bins = np.histogram(low_energy_slice, bins = reco_tracklength_bins)
    high_energy_slice_hist, high_energy_slice_bins = np.histogram(high_energy_slice, bins = reco_tracklength_bins)
    
    low_energy_slice_histX, low_energy_slice_histY = histogram_arr_handle(low_energy_slice_hist, low_energy_slice_bins)
    high_energy_slice_histX, high_energy_slice_histY = histogram_arr_handle(high_energy_slice_hist, high_energy_slice_bins)

    
    mp.plot(low_energy_slice_histX, low_energy_slice_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.plot(high_energy_slice_histX, high_energy_slice_histY, color = orange_cbf, linestyle = ':', linewidth = 2)
    mp.fill_between(low_energy_slice_histX, 0, low_energy_slice_histY, alpha = 0.6, color = blue_cbf, hatch = '//', label = 'low energy [150 - 250]')
    mp.fill_between(high_energy_slice_histX, 0, high_energy_slice_histY, alpha = 0.6, color = orange_cbf, hatch = '\\\\', label = 'high energy [> 400]')
    mp.xlabel('True track length - Reco track length [$\\frac{g}{cm^2}$]')
    mp.ylabel('#')
    mp.legend(loc = 'upper left', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.savefig('%s/Low_energy_slice.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    
    #reco_pure_length = np.sqrt((Reco_End[:, 0] - Reco_Start[:, 0])**2 + (Reco_End[:, 1] - Reco_Start[:, 1])**2 + (Reco_End[:, 2] - Reco_Start[:, 2])**2)
    #true_pure_length = np.sqrt((Primary_True_End[:, 0] - Primary_True_Start[:, 0])**2 + (Primary_True_End[:, 1] - Primary_True_Start[:, 1])**2 + (Primary_True_End[:, 2] - Primary_True_Start[:, 2])**2)
    
    #reco_pure_length_hist, reco_pure_length_bins = np.histogram(reco_pure_length, bins = 40, range = (min(min(reco_pure_length), min(true_pure_length)), max(max(reco_pure_length), max(true_pure_length))))
    #true_pure_length_hist, true_pure_length_bins = np.histogram(true_pure_length, bins = reco_pure_length_bins)
    
    #reco_pure_length_histX, reco_pure_length_histY = histogram_arr_handle(reco_pure_length_hist, reco_pure_length_bins)
    #true_pure_length_histX, true_pure_length_histY = histogram_arr_handle(true_pure_length_hist, true_pure_length_bins)
    
    #mp.plot(true_pure_length_histX, true_pure_length_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    #mp.plot(reco_pure_length_histX, reco_pure_length_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    #mp.fill_between(true_pure_length_histX, 0, true_pure_length_histY, alpha = 0.6, color = blue_cbf, hatch = '\\\\')
    #mp.fill_between(reco_pure_length_histX, 0, reco_pure_length_histY, alpha = 0.6, color = orange_cbf, hatch = '//')
    #mp.xlabel('Track length [mm]')
    #mp.ylabel('#')
    #mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    #mp.grid(True, linestyle = '--', alpha = 0.2)
    #mp.savefig('%s/Reco_pure_length.png' % out_dir, bbox_inches = 'tight')
    #mp.close();
    
    length_density_hist, length_density_binsX, length_density_binsY = np.histogram2d(reco_pure_length, Reco_TrackLength / 10, bins = [30, 30], range = [[min(reco_pure_length), 6500], [min(Reco_TrackLength / 10), max(Reco_TrackLength / 10)]])
    
    im = mp.pcolormesh(length_density_binsX, length_density_binsY, np.transpose(length_density_hist), cmap = cmap)
    mp.xlabel('Reco track length [mm]')
    mp.ylabel('Reco track density length [$\\frac{g}{cm^2}$]')
    mp.colorbar(im)
    mp.savefig('%s/TrackLength_Density.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    length_density_true_hist, length_density_true_binsX, length_density_true_binsY = np.histogram2d(true_pure_length, True_TrackLength, bins = [30, 30], range = [[min(reco_pure_length), 6500], [min(Reco_TrackLength / 10), max(Reco_TrackLength / 10)]])
    
    im = mp.pcolormesh(length_density_true_binsX, length_density_true_binsY, np.transpose(length_density_true_hist), cmap = cmap)
    mp.xlabel('True track length [mm]')
    mp.ylabel('True track density length [$\\frac{g}{cm^2}$]')
    mp.colorbar(im)
    mp.savefig('%s/TrackLength_Density_truth.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    ### Stopping vs exiting evaluation
    opposite_direction_counter = 0
    stopping_and_exiting_counter = 0
    stopping_counter = 0
    exiting_counter = 0
    exiting_truth = 0
    exiting_reco = 0
    stopping_truth = 0
    stopping_reco = 0
    
    for i in range(len(Reco_Start)):
        true_slope, true_intercept = line_slope_intercept(Primary_True_Start[i, 1], Primary_True_End[i, 1], Primary_True_Start[i, 2], Primary_True_End[i, 2])
        reco_slope, reco_intercept = line_slope_intercept(Reco_Start[i, 1], Reco_End[i, 1], Reco_Start[i, 2], Reco_End[i, 2])
        
        truth_exiting = False
        truth_stopping = False
        reco_exiting = False
        reco_stopping = False
        
        # check if line crosses detector volume in yz
        if (true_slope * 18535 + true_intercept) < -2949 and Primary_True_End[i, 1] <= (-2939):# + 338.3):
            truth_exiting = True
        elif Primary_True_End[i, 2] < 17750:
            truth_stopping = True
        if (reco_slope * 18535 + reco_intercept) < -2949 and Reco_End[i, 1] < (-2949 + 338.3):
            reco_exiting = True
        elif Reco_End[i, 2] < 17750:
            reco_stopping = True
        
        if truth_exiting:
            exiting_truth += 1
        elif truth_stopping:
            stopping_truth += 1
        
        if reco_exiting:
            exiting_reco += 1
        elif reco_stopping:
            stopping_reco += 1
        
        if truth_exiting and reco_exiting:
            stopping_and_exiting_counter += 1
            exiting_counter += 1
        elif truth_stopping and reco_stopping:
            stopping_and_exiting_counter += 1
            stopping_counter += 1
        
        if (true_slope < 0 and reco_slope > 0) or (true_slope > 0 and reco_slope < 0):
            opposite_direction_counter += 1
    
    efficiency_stop, purity_stop, accuracy_stop = eff_pur_acc(stopping_counter, stopping_truth, stopping_reco)
    efficiency_exit, purity_exit, accuracy_exit = eff_pur_acc(exiting_counter, exiting_truth, exiting_reco)
    
    
    print("Number of opposite directions: ", opposite_direction_counter, " (", opposite_direction_counter / len(Reco_Start) * 100, "%)")
    print("Stopping and exiting total: ", stopping_and_exiting_counter, " (", stopping_and_exiting_counter / len(Reco_Start) * 100, "%)")
    print("Exiting: ", exiting_counter, " (", exiting_counter / len(Reco_Start) * 100, "%) | eff.: ", efficiency_exit * 100, " , pur.: ", purity_exit * 100, " , acc.: ", accuracy_exit * 100)
    print("Stopping: ", stopping_counter, " (", stopping_counter / len(Reco_Start) * 100, "%) | eff.: ", efficiency_stop * 100, " , pur.: ", purity_stop * 100, " , acc.: ", accuracy_stop * 100)
    print("Truth: exiting: ", exiting_truth, " stopping: ", stopping_truth)
    print("Reco: exiting: ", exiting_reco, " stopping: ", stopping_reco)
    
    v = venn2(subsets = (stopping_truth - stopping_counter, stopping_reco - stopping_counter, stopping_counter), set_labels = ('Truth', 'Reco'), set_colors = ('#03558e', '#e95125', '#f39e92'), alpha = 1)
    c = venn2_circles(subsets = (stopping_truth - stopping_counter, stopping_reco - stopping_counter, stopping_counter), color = 'grey', linestyle = 'dashed', linewidth = 1.5)
    v.get_label_by_id('10').set_text('%i' % stopping_truth)
    v.get_label_by_id('10').set_color('lightgrey')
    v.get_label_by_id('10').set_fontsize('x-large')
    v.get_label_by_id('A').set_color('#03558e')
    v.get_label_by_id('A').set_fontsize('xx-large')
    v.get_patch_by_id('10').set_color('#03558e')
    v.get_label_by_id('01').set_text('%i' % stopping_reco)
    v.get_label_by_id('01').set_color('black')
    v.get_label_by_id('01').set_fontsize('x-large')
    v.get_label_by_id('B').set_color('#e95125')
    v.get_label_by_id('B').set_fontsize('xx-large')
    v.get_patch_by_id('01').set_color('#e95125')
    if stopping_counter > 0:
        v.get_label_by_id('11').set_text('%i' % stopping_counter)
        v.get_label_by_id('11').set_color('black')
        v.get_label_by_id('11').set_fontsize('x-large')
        v.get_patch_by_id('11').set_color('#f39e92')
    mp.title('Stopping')
    mp.savefig('%s/Stopping_venn.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    v = venn2(subsets = (exiting_truth - exiting_counter, exiting_reco - exiting_counter, exiting_counter), set_labels = ('Truth', 'Reco'), set_colors = ('#03558e', '#e95125', '#f39e92'), alpha = 1)
    c = venn2_circles(subsets = (exiting_truth - exiting_counter, exiting_reco - exiting_counter, exiting_counter), color = 'grey', linestyle = 'dashed', linewidth = 1.5)
    v.get_label_by_id('10').set_text('%i' % exiting_truth)
    v.get_label_by_id('10').set_color('lightgrey')
    v.get_label_by_id('10').set_fontsize('x-large')
    v.get_label_by_id('A').set_color('#03558e')
    v.get_label_by_id('A').set_fontsize('xx-large')
    v.get_patch_by_id('10').set_color('#03558e')
    v.get_label_by_id('01').set_text('%i' % exiting_reco)
    v.get_label_by_id('01').set_color('black')
    v.get_label_by_id('01').set_fontsize('x-large')
    v.get_label_by_id('B').set_color('#e95125')
    v.get_label_by_id('B').set_fontsize('xx-large')
    v.get_patch_by_id('01').set_color('#e95125')
    if exiting_counter > 0:
        v.get_label_by_id('11').set_text('%i' % exiting_counter)
        v.get_label_by_id('11').set_color('black')
        v.get_label_by_id('11').set_fontsize('x-large')
        v.get_patch_by_id('11').set_color('#f39e92')
    mp.title('Exiting')
    mp.savefig('%s/Exiting_venn.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    
    ### Charge ID
    # Filter out all events without truth +/-13 as PDG
    boolean_true_muon = (np.abs(True_Charge) == 13)
    True_Charge = True_Charge[boolean_true_muon]
    Reco_Charge = Reco_Charge[boolean_true_muon]
    True_KE = True_KE[boolean_true_muon]
    boolean_issues = (Reco_Charge == -9.99999999e+08)
    not_identified = Reco_Charge[boolean_issues]
    #True_KE = True_KE[Reco_Charge != -9.99999999e+08]

    # Counter for both positive, both negative, different reco/truth positive and negative
    true_muons = len(True_Charge[True_Charge == 13])
    true_antimuons = len(True_Charge[True_Charge == -13])
    reco_muons = len(Reco_Charge[Reco_Charge == 13])
    reco_antimuons = len(Reco_Charge[Reco_Charge == -13])
    counter_true_positive = 0   #truth and reco agree on muon
    counter_true_negative = 0   #truth and reco agree on antimuon
    counter_false_positive = 0  #truth and reco disagree, reco -> muon
    counter_false_negative = 0  #truth and reco disagree, reco -> antimuon
    
    # Prepare the KE arrays
    Muons_KE = np.ones(len(True_KE), dtype = float) * -9999.
    AMuons_KE = np.ones(len(True_KE), dtype = float) * -9999.
    True_Muons_KE = True_KE[True_Charge == 13]
    True_Antimuons_KE = True_KE[True_Charge == -13]
    
    # Compare sign of truth and reco    
    for i in range(len(True_Charge)):
        if True_Charge[i] == Reco_Charge[i]:
            if True_Charge[i] > 0:
                counter_true_positive += 1
                Muons_KE[i] = True_KE[i]
            else:
                counter_true_negative += 1
                AMuons_KE[i] = True_KE[i]
        else:
            if Reco_Charge[i] > 0:
                counter_false_positive += 1
            else:
                counter_false_negative += 1
    
    if (counter_true_positive + counter_true_negative + counter_false_positive + counter_false_negative) != len(True_Charge):
        print("Counters don't add up")
        print("Length: ", len(True_Charge))
    
    print("Charge ID numbers")
    print("  Not identified:  ", len(not_identified))
    print("  True muons:      ", counter_true_positive)
    print("  True antimuons:  ", counter_true_negative)
    print("  False muons:     ", counter_false_positive)
    print("  False antimuons: ", counter_false_negative, " here the not identified go into!")

    #TODO look into 'edge' cases with region jumps
    
    # Performance evaluation
    efficiency_muons = counter_true_positive / (counter_true_positive + (counter_false_negative - len(not_identified)))
    efficiency_antimuons = counter_true_negative / (counter_true_negative + counter_false_positive)
    purity_muons = counter_true_positive / (counter_true_positive + counter_false_positive)
    purity_antimuons = counter_true_negative / (counter_true_negative + (counter_false_negative - len(not_identified)))
    accuracy_anti_muons = (counter_true_positive + counter_true_negative) / (counter_true_positive + counter_true_negative + counter_false_positive + (counter_false_negative - len(not_identified)))
    
    print("  Muons (efficiency | purity):     ", efficiency_muons, " | ", purity_muons)
    print("  Antimuons (efficiency | purity): ", efficiency_antimuons, " | ", purity_antimuons)
    print("  Accuracy (both): ", accuracy_anti_muons)
    
    # Plotting
    # Total
    v = venn2(subsets = (true_muons - counter_true_positive, reco_muons- counter_true_positive, counter_true_positive), set_labels = ('Truth', 'Reco'), set_colors = ('#03558e', '#e95125', '#f39e92'), alpha = 1)
    c = venn2_circles(subsets = (true_muons - counter_true_positive, reco_muons - counter_true_positive, counter_true_positive), color = 'grey', linestyle = 'dashed', linewidth = 1.5)
    v.get_label_by_id('10').set_text('%i' % true_muons)
    v.get_label_by_id('10').set_color('lightgrey')
    v.get_label_by_id('10').set_fontsize('x-large')
    v.get_label_by_id('A').set_color('#03558e')
    v.get_label_by_id('A').set_fontsize('xx-large')
    v.get_patch_by_id('10').set_color('#03558e')
    v.get_label_by_id('01').set_text('%i' % reco_muons)
    v.get_label_by_id('01').set_color('black')
    v.get_label_by_id('01').set_fontsize('x-large')
    v.get_label_by_id('B').set_color('#e95125')
    v.get_label_by_id('B').set_fontsize('xx-large')
    v.get_patch_by_id('01').set_color('#e95125')
    if counter_true_positive > 0:
        v.get_label_by_id('11').set_text('%i' % counter_true_positive)
        v.get_label_by_id('11').set_color('black')
        v.get_label_by_id('11').set_fontsize('x-large')
        v.get_patch_by_id('11').set_color('#f39e92')
    mp.title('Charge ID $\\mu$')
    mp.savefig('%s/Muons_venn.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    v = venn2(subsets = (true_antimuons - counter_true_negative, reco_antimuons - counter_true_negative, counter_true_negative), set_labels = ('Truth', 'Reco'), set_colors = ('#03558e', '#e95125', '#f39e92'), alpha = 1)
    c = venn2_circles(subsets = (true_antimuons - counter_true_negative, reco_antimuons - counter_true_negative, counter_true_negative), color = 'grey', linestyle = 'dashed', linewidth = 1.5)
    v.get_label_by_id('10').set_text('%i' % true_antimuons)
    v.get_label_by_id('10').set_color('lightgrey')
    v.get_label_by_id('10').set_fontsize('x-large')
    v.get_label_by_id('A').set_color('#03558e')
    v.get_label_by_id('A').set_fontsize('xx-large')
    v.get_patch_by_id('10').set_color('#03558e')
    v.get_label_by_id('01').set_text('%i' % reco_antimuons)
    v.get_label_by_id('01').set_color('black')
    v.get_label_by_id('01').set_fontsize('x-large')
    v.get_label_by_id('B').set_color('#e95125')
    v.get_label_by_id('B').set_fontsize('xx-large')
    v.get_patch_by_id('01').set_color('#e95125')
    if counter_true_negative > 0:
        v.get_label_by_id('11').set_text('%i' % counter_true_negative)
        v.get_label_by_id('11').set_color('black')
        v.get_label_by_id('11').set_fontsize('x-large')
        v.get_patch_by_id('11').set_color('#f39e92')
    mp.title('Charge ID $\\bar{\\mu}$')
    mp.savefig('%s/Antimuons_venn.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
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
    
    mp.plot(muons_ke_hist_x, muons_ke_hist_y / true_muons_ke_hist_y, color = blue_cbf, linestyle = '-', linewidth = 2)
    mp.plot(amuons_ke_hist_x, amuons_ke_hist_y / true_antimuons_ke_hist_y, color = orange_cbf, linestyle = '-', linewidth = 2)
    mp.fill_between(muons_ke_hist_x, 0, muons_ke_hist_y / true_muons_ke_hist_y, color = blue_cbf, alpha = 0.3, hatch = '//', label = '$\\mu$ correct')
    mp.fill_between(amuons_ke_hist_x, 0, amuons_ke_hist_y / true_antimuons_ke_hist_y, color = orange_cbf, alpha = 0.3, hatch = '\\\\', label = '$\\bar{\\mu}$ correct')
    mp.xlabel('True Muon KE [MeV]')
    mp.ylabel('Fraction')
    mp.legend(loc = 'lower center', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.savefig('%s/Charge_ID_KE.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    return

def eff_pur_acc(overlap, truth, reco):
    if truth == 0 or reco == 0:
        efficiency = 0
        purity = 0
        accuracy = 0
    else:
        efficiency = overlap / truth
        purity = overlap / reco
        accuracy = overlap / (truth + reco - overlap)
    return efficiency, purity, accuracy


def line_slope_intercept(y_start, y_end, z_start, z_end):
    slope = (y_end - y_start) / (z_end - z_start)
    intercept = y_start - slope * z_start
    return slope, intercept

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

def check_orientation(hit_z):
    return layer_dict["%s" % hit_z]

first_z = 11185
layer_dict = { "%s" % first_z : "kUBar" }
        
def calculate_layers(Xlayers):
    increment = 2
    if Xlayers:
        increment = 3
    thin_layers = 49
    thick_layers = 34
    double_layers = 9
    # Calculate the z position for each layer for the thin section
    for i in range(thin_layers):
        hit_z = first_z + i * 65
        # even layers
        if not Xlayers:
            if (((hit_z - first_z) / 65) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
            # odd layers
            elif (((hit_z - first_z) / 65) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "kVBar" })
        # x layers
        if Xlayers:
            # even layers
            if (((hit_z - first_z) / 65) % 2) == 0:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
            # odd layers
            elif (((hit_z - first_z) / 65) % 2) == 1:
                layer_dict.update({ "%s" % hit_z : "kVBar" })
            if (((hit_z - first_z) / 65) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "kXBar" })
    # Calculate the z position for each layers the thick section
    start_thick = first_z + thin_layers * 65
    for i in range(thick_layers):
        hit_z = start_thick + i * 90
        if not Xlayers:
            # even layers
            if (((hit_z - start_thick) / 90) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "kVBar" })
            # odd layers
            elif (((hit_z - start_thick) / 90) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
        # x layers
        if Xlayers:
            # even layers
            if (((hit_z - start_thick) / 90) % 2) == 0:
                layer_dict.update({ "%s" % hit_z : "kVBar" })
            # odd layers
            elif (((hit_z - start_thick) / 90) % 2) == 1:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
            if (((hit_z - start_thick) / 90) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "kXBar" })
    # Calculate the z position for each layers the double section
    start_double = first_z + thin_layers * 65 + thick_layers * 90
    for i in range(double_layers):
        hit_z = start_double + i * 130
        if not Xlayers:
            # even layers
            if (((hit_z - start_double) / 130) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
            # odd layers
            elif (((hit_z - start_thick) / 90) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "kVBar" })
        # x layers
        if Xlayers:
            # even layers
            if (((hit_z - start_double) / 130) % 2) == 0:
                layer_dict.update({ "%s" % hit_z : "kUBar" })
            # odd layers
        elif (((hit_z - start_double) / 130) % 2) == 1:
            layer_dict.update({ "%s" % hit_z : "kVBar" })
        if ((( hit_z - start_double) / 130) % increment) == 0:
            layer_dict.update({ "%s" % hit_z : "kXBar" })

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--Xlayers', "-X", help = "Does the geometry use X (90 degree orientated) scintillator layers? Yes -> --Xlayers, No -> --no-Xlayers", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    Xlayers = args.Xlayers
    
    calculate_layers(Xlayers)
    #print(layer_dict)
    out_dir = args.outdir
    input_filename = args.input_filename
    draw_performance(out_dir, input_filename, Xlayers)

