import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll
import matplotlib.cm as cm

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-poster')

mp.rc('axes', labelsize = 12)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 12) # fontsize of the tick labels
mp.rc('ytick', labelsize = 12) # fontsize of the tick labels

### Actual function that loops through the spills
def draw_comparison(out_dir, input_filename1, input_filename2):
    if not os.path.exists(input_filename1): raise ValueError(f"Cannor find input_filename1 {input_filename1}")
    if not os.path.exists(input_filename2): raise ValueError(f"Cannor find input_filename2 {input_filename2}")
    
    # Make sure we read in the correct file and have the output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir):
        raise ValueError(f"Could not make out_dir {out_dir}")
            
    # Read in the Reco_Tree that contains the TMS_Tracks
    r1 = ROOT.TChain("Reco_Tree")
    r1.Add(input_filename1)
    print("N entries:", r1.GetEntries())
    if not r1.GetEntries() > 0:
        print("Didn't get any entries, are you sure the input_filename1 is right?\n", input_filename1)
    
    truth1 = ROOT.TChain("Truth_Info")
    truth1.Add(input_filename1)
    if not truth1.GetEntries() > 0:
        print("Didn't get any entries in Truth_Info, are you sure the input_filename1 is right?\n", input_filename1)
        
        # Read in the Reco_Tree that contains the TMS_Tracks
    r2 = ROOT.TChain("Reco_Tree")
    r2.Add(input_filename2)
    print("N2 entries:", r2.GetEntries())
    if not r2.GetEntries() > 0:
        print("Didn't get any entries, are you sure the input_filename2 is right?\n", input_filename2)
    
    truth2 = ROOT.TChain("Truth_Info")
    truth2.Add(input_filename2)
    if not truth2.GetEntries() > 0:
        print("Didn't get any entries in Truth_Info, are you sure the input_filename2 is right?\n", input_filename2)
            
    max_n_spills = 10000 # TODO (old) add some meta info to output file with n spill info for file
    
    spill_number_cache = dict()
    n_events_1 = r1.GetEntries()
    
    n_events_2 = r2.GetEntries()
    
    Reco_Start_1 = np.ones((n_events_1, 5, 3), dtype = float) * -9999.
    Reco_End_1 = np.ones((n_events_1, 5, 3), dtype = float) * -9999.
    Primary_True_Start_1 = np.ones((n_events_1, 5, 3), dtype = float) * -9999.
    Primary_True_End_1 = np.ones((n_events_1, 5, 3), dtype = float) * -9999.  
    
    Reco_Start_2 = np.ones((n_events_2, 5, 3), dtype = float) * -9999.
    Reco_End_2 = np.ones((n_events_2, 5, 3), dtype = float) * -9999.
    Primary_True_Start_2 = np.ones((n_events_2, 5, 3), dtype = float) * -9999.
    Primary_True_End_2 = np.ones((n_events_2, 5, 3), dtype = float) * -9999.  
    
    for current_spill_number in range(max_n_spills):
        for i in range(n_events_1):
            try:
                spill_number = spill_number_cache[i]
                event = None
                true_event = None
            except KeyError:
                r1.GetEntry(i)
                event = r1
                truth1.GetEntry(i)
                true_event = truth1
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event == None:
                r1.GetEntry(i)
                event = r1
            if true_event == None:
                truth1.GetEntry(i)
                true_event = truth1

            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
            
            nTracks = event.nTracks
            if nTracks <= 0: continue
            if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break
            
                if RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] > -8000. and not StartPos.size == 0:
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:
                        Reco_Start_1[i, j, 0] = StartPos[j*3 + 0]
                        Reco_Start_1[i, j, 1] = StartPos[j*3 + 1]
                        Reco_Start_1[i, j, 2] = StartPos[j*3 + 2]
                        Primary_True_Start_1[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0]
                        Primary_True_Start_1[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1]
                        Primary_True_Start_1[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]
                if RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] > -8000. and not EndPos.size == 0:
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:
                        Reco_End_1[i, j, 0] = EndPos[j*3 + 0]
                        Reco_End_1[i, j, 1] = EndPos[j*3 + 1]
                        Reco_End_1[i, j, 2] = EndPos[j*3 + 2]
                        Primary_True_End_1[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0]
                        Primary_True_End_1[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1]
                        Primary_True_End_1[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2]
        
        for i in range(n_events_2):
            try:
                spill_number = spill_number_cache[i]
                event2 = None
                true_event2 = None
            except KeyError:
                r2.GetEntry(i)
                event2 = r2
                truth2.GetEntry(i)
                true_event2 = truth2
                spill_number = event2.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event2 == None:
                r2.GetEntry(i)
                event2 = r2
            if true_event2 == None:
                truth2.GetEntry(i)
                true_event2 = truth2

            StartPos = np.frombuffer(event2.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event2.EndPos, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event2.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event2.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
            
            nTracks = event2.nTracks
            if nTracks <= 0: continue
            if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break
            
                if RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] > -8000. and not StartPos.size == 0:
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:
                        Reco_Start_2[i, j, 0] = StartPos[j*3 + 0]
                        Reco_Start_2[i, j, 1] = StartPos[j*3 + 1]
                        Reco_Start_2[i, j, 2] = StartPos[j*3 + 2]
                        Primary_True_Start_2[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0]
                        Primary_True_Start_2[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1]
                        Primary_True_Start_2[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]
                if RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] > -8000. and not EndPos.size == 0:
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] - RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2]) > 890.:
                        Reco_End_2[i, j, 0] = EndPos[j*3 + 0]
                        Reco_End_2[i, j, 1] = EndPos[j*3 + 1]
                        Reco_End_2[i, j, 2] = EndPos[j*3 + 2]
                        Primary_True_End_2[i, j, 0] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0]
                        Primary_True_End_2[i, j, 1] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1]
                        Primary_True_End_2[i, j, 2] = RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2]
    
    # filter out not filled indice
    boolean_Reco_Start_1 = (Reco_Start_1[:, :, 0] != -9999.)
    Reco_Start_1 = Reco_Start_1[boolean_Reco_Start_1]
    boolean_Reco_End_1 = (Reco_End_1[:, :, 0] != -9999.)
    Reco_End_1 = Reco_End_1[boolean_Reco_End_1]
    boolean_Primary_Start_1 = (Primary_True_Start_1[:, :, 0] != -9999.)
    Primary_True_Start_1 = Primary_True_Start_1[boolean_Primary_Start_1]
    boolean_Primary_End_1 = (Primary_True_End_1[:, :, 0] != -9999.)
    Primary_True_End_1 = Primary_True_End_1[boolean_Primary_End_1]
    
    boolean_Reco_Start_2 = (Reco_Start_2[:, :, 0] != -9999.)
    Reco_Start_2 = Reco_Start_2[boolean_Reco_Start_2]
    boolean_Reco_End_2 = (Reco_End_2[:, :, 0] != -9999.)
    Reco_End_2 = Reco_End_2[boolean_Reco_End_2]
    boolean_Primary_Start_2 = (Primary_True_Start_2[:, :, 0] != -9999.)
    Primary_True_Start_2 = Primary_True_Start_2[boolean_Primary_Start_2]
    boolean_Primary_End_2 = (Primary_True_End_2[:, :, 0] != -9999.)
    Primary_True_End_2 = Primary_True_End_2[boolean_Primary_End_2]
    
    # subtract reconstruction from truth for all directions
    Diff_Start_x_1 = Primary_True_Start_1[:, 0] - Reco_Start_1[:, 0]
    Diff_Start_y_1 = Primary_True_Start_1[:, 1] - Reco_Start_1[:, 1]
    Diff_Start_z_1 = Primary_True_Start_1[:, 2] - Reco_Start_1[:, 2]
    Diff_End_x_1 = Primary_True_End_1[:, 0] - Reco_End_1[:, 0]
    Diff_End_y_1 = Primary_True_End_1[:, 1] - Reco_End_1[:, 1]
    Diff_End_z_1 = Primary_True_End_1[:, 2] - Reco_End_1[:, 2]
    
    Diff_Start_x_2 = Primary_True_Start_2[:, 0] - Reco_Start_2[:, 0]
    Diff_Start_y_2 = Primary_True_Start_2[:, 1] - Reco_Start_2[:, 1]
    Diff_Start_z_2 = Primary_True_Start_2[:, 2] - Reco_Start_2[:, 2]
    Diff_End_x_2 = Primary_True_End_2[:, 0] - Reco_End_2[:, 0]
    Diff_End_y_2 = Primary_True_End_2[:, 1] - Reco_End_2[:, 1]
    Diff_End_z_2 = Primary_True_End_2[:, 2] - Reco_End_2[:, 2]

    # create histograms for the differences
    Diff_Start_x_1_hist, Diff_Start_x_bins = np.histogram(Diff_Start_x_1, bins = 50, range = (min(min(Diff_Start_x_1), min(Diff_Start_x_2)), max(max(Diff_Start_x_1), max(Diff_Start_x_2))))
    Diff_Start_y_1_hist, Diff_Start_y_bins = np.histogram(Diff_Start_y_1, bins = 50, range = (min(min(Diff_Start_y_1), min(Diff_Start_y_2)), max(max(Diff_Start_y_1), max(Diff_Start_y_2))))
    Diff_Start_z_1_hist, Diff_Start_z_bins = np.histogram(Diff_Start_z_1, bins = 50, range = (min(min(Diff_Start_z_1), min(Diff_Start_z_2)), max(max(Diff_Start_z_1), max(Diff_Start_z_2))))
    Diff_End_x_1_hist, Diff_End_x_bins = np.histogram(Diff_End_x_1, bins = 50, range = (min(min(Diff_End_x_1), min(Diff_End_x_2)), max(max(Diff_End_x_1), max(Diff_End_x_2))))
    Diff_End_y_1_hist, Diff_End_y_bins = np.histogram(Diff_End_y_1, bins = 50, range = (min(min(Diff_End_y_1), min(Diff_End_y_2)), max(max(Diff_End_y_1), max(Diff_End_y_2))))
    Diff_End_z_1_hist, Diff_End_z_bins = np.histogram(Diff_End_z_1, bins = 50, range = (min(min(Diff_End_z_1), min(Diff_End_z_2)), max(max(Diff_End_z_1), max(Diff_End_z_2))))
    
    Diff_Start_x_1_histX, Diff_Start_x_1_histY = histogram_arr_handle(Diff_Start_x_1_hist, Diff_Start_x_bins)
    Diff_Start_y_1_histX, Diff_Start_y_1_histY = histogram_arr_handle(Diff_Start_y_1_hist, Diff_Start_y_bins)
    Diff_Start_z_1_histX, Diff_Start_z_1_histY = histogram_arr_handle(Diff_Start_z_1_hist, Diff_Start_z_bins)
    Diff_End_x_1_histX, Diff_End_x_1_histY = histogram_arr_handle(Diff_End_x_1_hist, Diff_End_x_bins)
    Diff_End_y_1_histX, Diff_End_y_1_histY = histogram_arr_handle(Diff_End_y_1_hist, Diff_End_y_bins)
    Diff_End_z_1_histX, Diff_End_z_1_histY = histogram_arr_handle(Diff_End_z_1_hist, Diff_End_z_bins)
    
    Diff_Start_x_2_hist, Diff_Start_x_2_bins = np.histogram(Diff_Start_x_2, bins = Diff_Start_x_bins)
    Diff_Start_y_2_hist, Diff_Start_y_2_bins = np.histogram(Diff_Start_y_2, bins = Diff_Start_y_bins)
    Diff_Start_z_2_hist, Diff_Start_z_2_bins = np.histogram(Diff_Start_z_2, bins = Diff_Start_z_bins)
    Diff_End_x_2_hist, Diff_End_x_2_bins = np.histogram(Diff_End_x_2, bins = Diff_End_x_bins)
    Diff_End_y_2_hist, Diff_End_y_2_bins = np.histogram(Diff_End_y_2, bins = Diff_End_y_bins)
    Diff_End_z_2_hist, Diff_End_z_2_bins = np.histogram(Diff_End_z_2, bins = Diff_End_z_bins)
    
    Diff_Start_x_2_histX, Diff_Start_x_2_histY = histogram_arr_handle(Diff_Start_x_2_hist, Diff_Start_x_bins)
    Diff_Start_y_2_histX, Diff_Start_y_2_histY = histogram_arr_handle(Diff_Start_y_2_hist, Diff_Start_y_bins)
    Diff_Start_z_2_histX, Diff_Start_z_2_histY = histogram_arr_handle(Diff_Start_z_2_hist, Diff_Start_z_bins)
    Diff_End_x_2_histX, Diff_End_x_2_histY = histogram_arr_handle(Diff_End_x_2_hist, Diff_End_x_bins)
    Diff_End_y_2_histX, Diff_End_y_2_histY = histogram_arr_handle(Diff_End_y_2_hist, Diff_End_y_bins)
    Diff_End_z_2_histX, Diff_End_z_2_histY = histogram_arr_handle(Diff_End_z_2_hist, Diff_End_z_bins)
    
    # plot
    mp.plot(Diff_Start_x_1_histX, Diff_Start_x_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_Start_x_1_histX, 0, Diff_Start_x_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_Start_x_2_histX, Diff_Start_x_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_Start_x_2_histX, 0, Diff_Start_x_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('x')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_1_histX, Diff_Start_y_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_Start_y_1_histX, 0, Diff_Start_y_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_Start_y_2_histX, Diff_Start_y_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_Start_y_2_histX, 0, Diff_Start_y_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_1_histX, Diff_Start_z_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_Start_z_1_histX, 0, Diff_Start_z_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_Start_z_2_histX, Diff_Start_z_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_Start_z_2_histX, 0, Diff_Start_z_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8) 
    mp.savefig('%s/Difference_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_1_histX, Diff_End_x_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_End_x_1_histX, 0, Diff_End_x_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_End_x_2_histX, Diff_End_x_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_End_x_2_histX, 0, Diff_End_x_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('x')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_1_histX, Diff_End_y_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_End_y_1_histX, 0, Diff_End_y_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_End_y_2_histX, Diff_End_y_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_End_y_2_histX, 0, Diff_End_y_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_1_histX, Diff_End_z_1_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'stereo')
    mp.fill_between(Diff_End_z_1_histX, 0, Diff_End_z_1_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.plot(Diff_End_z_2_histX, Diff_End_z_2_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'hybrid')
    mp.fill_between(Diff_End_z_2_histX, 0, Diff_End_z_2_histY, alpha = 0.6, hatch = '\\\\', color = orange_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8) 
    mp.savefig('%s/Difference_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()  
        
    # create histograms for all directions for truth and reconstruction 


    Reco_Start_x_1_hist, Reco_Start_x_1_bins = np.histogram(Reco_Start_1[:, 0], bins = int((max(Reco_Start_1[:, 0]) - min(Reco_Start_1[:, 0])) / 35.42))#3, range = (-35.42, 35.42)) # #48
    Reco_Start_y_1_hist, Reco_Start_y_1_bins = np.histogram(Reco_Start_1[:, 1], bins = 40)#, range = (-4110., -910.))
    Reco_Start_z_1_hist, Reco_Start_z_1_bins = np.histogram(Reco_Start_1[:, 2], bins = 100)#, range = (11360., 18320.))
    Reco_End_x_1_hist, Reco_End_x_1_bins = np.histogram(Reco_End_1[:, 0], bins = int((max(Reco_End_1[:, 0]) - min(Reco_End_1[:, 0])) / 35.42))#48)#, range = (-3400., 3400.))
    Reco_End_y_1_hist, Reco_End_y_1_bins = np.histogram(Reco_End_1[:, 1], bins = 40)#, range = (-4110., -910.))
    Reco_End_z_1_hist, Reco_End_z_1_bins = np.histogram(Reco_End_1[:, 2], bins = 100)#, range = (11360., 18320.))
    Primary_True_Start_x_1_hist, Primary_True_Start_x_1_bins = np.histogram(Primary_True_Start_1[:, 0], bins = Reco_Start_x_1_bins) #48
    Primary_True_Start_y_1_hist, Primary_True_Start_y_1_bins = np.histogram(Primary_True_Start_1[:, 1], bins = Reco_Start_y_1_bins)#20)
    Primary_True_Start_z_1_hist, Primary_True_Start_z_1_bins = np.histogram(Primary_True_Start_1[:, 2], bins = 100)
    Primary_True_End_x_1_hist, Primary_True_End_x_1_bins = np.histogram(Primary_True_End_1[:, 0], bins = Reco_End_x_1_bins) #48)
    Primary_True_End_y_1_hist, Primary_True_End_y_1_bins = np.histogram(Primary_True_End_1[:, 1], bins = Reco_End_y_1_bins)#20)
    Primary_True_End_z_1_hist, Primary_True_End_z_1_bins = np.histogram(Primary_True_End_1[:, 2], bins = 100)
    
    Reco_Start_x_2_hist, Reco_Start_x_2_bins = np.histogram(Reco_Start_2[:, 0], bins = int((max(Reco_Start_2[:, 0]) - min(Reco_Start_2[:, 0])) / 35.42))#3, range = (-35.42, 35.42))#10)# #48)#
    Reco_Start_y_2_hist, Reco_Start_y_2_bins = np.histogram(Reco_Start_2[:, 1], bins = 40)#, range = (-4110., -910.))
    Reco_Start_z_2_hist, Reco_Start_z_2_bins = np.histogram(Reco_Start_2[:, 2], bins = 100)#, range = (11360., 18320.))
    Reco_End_x_2_hist, Reco_End_x_2_bins = np.histogram(Reco_End_2[:, 0], bins = int((max(Reco_End_2[:, 0]) - min(Reco_End_2[:, 0])) / 35.42)) #48)#, range = (-3400., 3400.))
    Reco_End_y_2_hist, Reco_End_y_2_bins = np.histogram(Reco_End_2[:, 1], bins = 40)#, range = (-4110., -910.))
    Reco_End_z_2_hist, Reco_End_z_2_bins = np.histogram(Reco_End_2[:, 2], bins = 100)#, range = (11360., 18320.))
    Primary_True_Start_x_2_hist, Primary_True_Start_x_2_bins = np.histogram(Primary_True_Start_2[:, 0], bins = Reco_Start_x_2_bins) #48)
    Primary_True_Start_y_2_hist, Primary_True_Start_y_2_bins = np.histogram(Primary_True_Start_2[:, 1], bins = Reco_Start_y_2_bins) #20)
    Primary_True_Start_z_2_hist, Primary_True_Start_z_2_bins = np.histogram(Primary_True_Start_2[:, 2], bins = 100)
    Primary_True_End_x_2_hist, Primary_True_End_x_2_bins = np.histogram(Primary_True_End_2[:, 0], bins = Reco_End_x_2_bins) #48)
    Primary_True_End_y_2_hist, Primary_True_End_y_2_bins = np.histogram(Primary_True_End_2[:, 1], bins = Reco_End_y_2_bins) #20)
    Primary_True_End_z_2_hist, Primary_True_End_z_2_bins = np.histogram(Primary_True_End_2[:, 2], bins = 100)
    
    # make the histograms usable
    Reco_Start_x_1_histX, Reco_Start_x_1_histY = histogram_arr_handle(Reco_Start_x_1_hist, Reco_Start_x_1_bins)
    Reco_Start_y_1_histX, Reco_Start_y_1_histY = histogram_arr_handle(Reco_Start_y_1_hist, Reco_Start_y_1_bins)
    Reco_Start_z_1_histX, Reco_Start_z_1_histY = histogram_arr_handle(Reco_Start_z_1_hist, Reco_Start_z_1_bins)
    Reco_End_x_1_histX, Reco_End_x_1_histY = histogram_arr_handle(Reco_End_x_1_hist, Reco_End_x_1_bins)
    Reco_End_y_1_histX, Reco_End_y_1_histY = histogram_arr_handle(Reco_End_y_1_hist, Reco_End_y_1_bins)
    Reco_End_z_1_histX, Reco_End_z_1_histY = histogram_arr_handle(Reco_End_z_1_hist, Reco_End_z_1_bins)
    Primary_True_Start_x_1_histX, Primary_True_Start_x_1_histY = histogram_arr_handle(Primary_True_Start_x_1_hist, Primary_True_Start_x_1_bins)
    Primary_True_Start_y_1_histX, Primary_True_Start_y_1_histY = histogram_arr_handle(Primary_True_Start_y_1_hist, Primary_True_Start_y_1_bins)
    Primary_True_Start_z_1_histX, Primary_True_Start_z_1_histY = histogram_arr_handle(Primary_True_Start_z_1_hist, Primary_True_Start_z_1_bins)
    Primary_True_End_x_1_histX, Primary_True_End_x_1_histY = histogram_arr_handle(Primary_True_End_x_1_hist, Primary_True_End_x_1_bins)
    Primary_True_End_y_1_histX, Primary_True_End_y_1_histY = histogram_arr_handle(Primary_True_End_y_1_hist, Primary_True_End_y_1_bins)
    Primary_True_End_z_1_histX, Primary_True_End_z_1_histY = histogram_arr_handle(Primary_True_End_z_1_hist, Primary_True_End_z_1_bins)        

    Reco_Start_x_2_histX, Reco_Start_x_2_histY = histogram_arr_handle(Reco_Start_x_2_hist, Reco_Start_x_2_bins)
    Reco_Start_y_2_histX, Reco_Start_y_2_histY = histogram_arr_handle(Reco_Start_y_2_hist, Reco_Start_y_2_bins)
    Reco_Start_z_2_histX, Reco_Start_z_2_histY = histogram_arr_handle(Reco_Start_z_2_hist, Reco_Start_z_2_bins)
    Reco_End_x_2_histX, Reco_End_x_2_histY = histogram_arr_handle(Reco_End_x_2_hist, Reco_End_x_2_bins)
    Reco_End_y_2_histX, Reco_End_y_2_histY = histogram_arr_handle(Reco_End_y_2_hist, Reco_End_y_2_bins)
    Reco_End_z_2_histX, Reco_End_z_2_histY = histogram_arr_handle(Reco_End_z_2_hist, Reco_End_z_2_bins)
    Primary_True_Start_x_2_histX, Primary_True_Start_x_2_histY = histogram_arr_handle(Primary_True_Start_x_2_hist, Primary_True_Start_x_2_bins)
    Primary_True_Start_y_2_histX, Primary_True_Start_y_2_histY = histogram_arr_handle(Primary_True_Start_y_2_hist, Primary_True_Start_y_2_bins)
    Primary_True_Start_z_2_histX, Primary_True_Start_z_2_histY = histogram_arr_handle(Primary_True_Start_z_2_hist, Primary_True_Start_z_2_bins)
    Primary_True_End_x_2_histX, Primary_True_End_x_2_histY = histogram_arr_handle(Primary_True_End_x_2_hist, Primary_True_End_x_2_bins)
    Primary_True_End_y_2_histX, Primary_True_End_y_2_histY = histogram_arr_handle(Primary_True_End_y_2_hist, Primary_True_End_y_2_bins)
    Primary_True_End_z_2_histX, Primary_True_End_z_2_histY = histogram_arr_handle(Primary_True_End_z_2_hist, Primary_True_End_z_2_bins)        

    # now plot this
    mp.fill_between(Reco_Start_x_2_histX, 0, Reco_Start_x_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)
    mp.fill_between(Reco_Start_x_1_histX, 0, Reco_Start_x_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_Start_x_1_histX, Primary_True_Start_x_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_Start_x_1_histX, Reco_Start_x_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_Start_x_2_histX, Primary_True_Start_x_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_Start_x_2_histX, Reco_Start_x_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('x')   
    mp.savefig('%s/Track_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.fill_between(Reco_Start_y_2_histX, 0, Reco_Start_y_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)
    mp.fill_between(Reco_Start_y_1_histX, 0, Reco_Start_y_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_Start_y_1_histX, Primary_True_Start_y_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_Start_y_1_histX, Reco_Start_y_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_Start_y_2_histX, Primary_True_Start_y_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_Start_y_2_histX, Reco_Start_y_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')    
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('y')   
    mp.savefig('%s/Track_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.fill_between(Reco_Start_z_2_histX, 0, Reco_Start_z_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)
    mp.fill_between(Reco_Start_z_1_histX, 0, Reco_Start_z_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_Start_z_1_histX, Primary_True_Start_z_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_Start_z_1_histX, Reco_Start_z_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_Start_z_2_histX, Primary_True_Start_z_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_Start_z_2_histX, Reco_Start_z_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')    
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('%s/Track_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.fill_between(Reco_End_x_2_histX, 0, Reco_End_x_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)
    mp.fill_between(Reco_End_x_1_histX, 0, Reco_End_x_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_End_x_1_histX, Primary_True_End_x_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_End_x_1_histX, Reco_End_x_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_End_x_2_histX, Primary_True_End_x_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_End_x_2_histX, Reco_End_x_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')    
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('x')   
    mp.savefig('%s/Track_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.fill_between(Reco_End_y_2_histX, 0, Reco_End_y_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)
    mp.fill_between(Reco_End_y_1_histX, 0, Reco_End_y_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_End_y_1_histX, Primary_True_End_y_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_End_y_1_histX, Reco_End_y_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_End_y_2_histX, Primary_True_End_y_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_End_y_2_histX, Reco_End_y_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('%s/Track_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    print('Truth stereo: ', np.mean(Primary_True_End_1[:, 1]), ' +/- ', np.std(Primary_True_End_1[:, 1]))
    print('Reco stereo:  ', np.mean(Reco_End_1[:, 1]), ' +/- ', np.std(Reco_End_1[:, 1]))
    print('Truth hybrid: ', np.mean(Primary_True_End_2[:, 1]), ' +/- ', np.std(Primary_True_End_2[:, 1]))
    print('Reco hybrid:  ', np.mean(Reco_End_2[:, 1]), ' +/- ', np.std(Reco_End_2[:, 1]))
    
    mp.fill_between(Reco_End_z_2_histX, 0, Reco_End_z_2_histY, color = blue_cbf, hatch = '\\\\', alpha = 0.3)    
    mp.fill_between(Reco_End_z_1_histX, 0, Reco_End_z_1_histY, color = orange_cbf, hatch = '//', alpha = 0.3)
    mp.plot(Primary_True_End_z_1_histX, Primary_True_End_z_1_histY, color = green_cbf, linestyle = '-', linewidth = 3, label = 'truth stereo')
    mp.plot(Reco_End_z_1_histX, Reco_End_z_1_histY, color = orange_cbf, linestyle = '--', linewidth = 2, label = 'reco stereo')
    mp.plot(Primary_True_End_z_2_histX, Primary_True_End_z_2_histY, color = black_cbf, linestyle = '-', linewidth = 3, label = 'truth hybrid')
    mp.plot(Reco_End_z_2_histX, Reco_End_z_2_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'reco hybrid')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('%s/Track_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
           
    return

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
    parser.add_argument('--input_filename1', "-f1", type = str, help = "The stereo file with the events to draw.")
    parser.add_argument('--input_filename2', "-f2", type = str, help = "The hybrid file with the events to draw.")
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    input_filename1 = args.input_filename1
    input_filename2 = args.input_filename2
    draw_comparison(out_dir, input_filename1, input_filename2)

