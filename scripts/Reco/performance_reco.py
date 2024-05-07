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

mp.rc('axes', labelsize = 12)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 12) # fontsize of the tick labels
mp.rc('ytick', labelsize = 12) # fontsize of the tick labels

### Actual function that loops through the spills
def draw_performance(out_dir, input_filename):
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
            
    max_n_spills = 121000 # TODO (old) add some meta info to output file with n spill info for file
    
    spill_number_cache = dict()
    n_events = r.GetEntries()
    
    Reco_Start = np.ones((n_events, 3), dtype = float) * -9999.
    Reco_End = np.ones((n_events, 3), dtype = float) * -9999.
    Primary_True_Start = np.ones((n_events, 3), dtype = float) * -9999.
    Primary_True_End = np.ones((n_events, 3), dtype = float) * -9999.
    True_TrackLength = np.ones(n_events, dtype = float) * -9999.
    Reco_TrackLength = np.ones(n_events, dtype = float) * -9999.
    
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
            
            nTracks = event.nTracks
            if nTracks <= 0: continue

            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
            #Muon_TrueTrackLength = true_event.Muon_TrueTrackLength
            #Reco_Track_Length = np.frombuffer(event.Length, dtype = np.float32)
            
            if RecoTrackPrimaryParticleTruePositionTrackStart[0] > -8000. and not StartPos.size == 0:
                # checking for muon tracks (minimal length for this are 20 planes traversed -> 890 mm in thin area
                if (EndPos[2] - StartPos[2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[2] - RecoTrackPrimaryParticleTruePositionTrackStart[2]) > 890.:
                    Reco_Start[i, 0] = StartPos[0]
                    Reco_Start[i, 1] = StartPos[1]
                    Reco_Start[i, 2] = StartPos[2]
                    #if RecoTrackPrimaryParticleTruePositionTrackStart[2] < 11362:
                    #    print('PDG: ', true_event.PDG[np.frombuffer(true_event.RecoTrackPrimaryParticleIndex, dtype = np.uint8)[0]], ' End: ', RecoTrackPrimaryParticleTruePositionTrackEnd[2])
                    Primary_True_Start[i, 0] = RecoTrackPrimaryParticleTruePositionTrackStart[0]
                    Primary_True_Start[i, 1] = RecoTrackPrimaryParticleTruePositionTrackStart[1]
                    Primary_True_Start[i, 2] = RecoTrackPrimaryParticleTruePositionTrackStart[2]
                    #True_TrackLength[i] = Muon_TrueTrackLength
                    #Reco_TrackLength[i] = Reco_Track_Length
            if RecoTrackPrimaryParticleTruePositionTrackEnd[0] > -8000. and not EndPos.size == 0:
                if (EndPos[2] - StartPos[2]) > 890. and (RecoTrackPrimaryParticleTruePositionTrackEnd[2] - RecoTrackPrimaryParticleTruePositionTrackStart[2]) > 890.:            
                    Reco_End[i, 0] = EndPos[0]
                    Reco_End[i, 1] = EndPos[1]
                    Reco_End[i, 2] = EndPos[2]
                    Primary_True_End[i, 0] = RecoTrackPrimaryParticleTruePositionTrackEnd[0]
                    Primary_True_End[i, 1] = RecoTrackPrimaryParticleTruePositionTrackEnd[1]
                    Primary_True_End[i, 2] = RecoTrackPrimaryParticleTruePositionTrackEnd[2]
    
    # filter out not filled indice
    boolean_Reco_Start = (Reco_Start[:, 0] != -9999.)
    Reco_Start = Reco_Start[boolean_Reco_Start]
    boolean_Reco_End = (Reco_End[:, 0] != -9999.)
    Reco_End = Reco_End[boolean_Reco_End]
    boolean_Primary_Start = (Primary_True_Start[:, 0] != -9999.)
    Primary_True_Start = Primary_True_Start[boolean_Primary_Start]
    boolean_Primary_End = (Primary_True_End[:, 0] != -9999.)
    Primary_True_End = Primary_True_End[boolean_Primary_End]
    boolean_True_TrackLength = (True_TrackLength != -9999.)
    True_TrackLength = True_TrackLength[boolean_True_TrackLength]
    boolean_Reco_TrackLength = (Reco_TrackLength != -9999.)
    Reco_TrackLength = Reco_TrackLength[boolean_Reco_TrackLength]
    
    # total number of events after filtering
    print("#events reconstruction: ", len(Reco_Start), "# events truth: ", len(Primary_True_Start))
    
    # subtract reconstruction from truth for all directions
    Diff_Start_x = Primary_True_Start[:, 0] - Reco_Start[:, 0]
    Diff_Start_y = Primary_True_Start[:, 1] - Reco_Start[:, 1]# - 915
    Diff_Start_z = Primary_True_Start[:, 2] - Reco_Start[:, 2]
    Diff_End_x = Primary_True_End[:, 0] - Reco_End[:, 0]
    Diff_End_y = Primary_True_End[:, 1] - Reco_End[:, 1]# - 915
    Diff_End_z = Primary_True_End[:, 2] - Reco_End[:, 2]

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
    mp.title('Start x', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_Start_y, 0, max(Diff_Start_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_Start_y - std_Start_y, mean_Start_y + std_Start_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('Start y', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_Start_z, 0, max(Diff_Start_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_Start_z - std_Start_z, mean_Start_z + std_Start_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('Start z', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_x, 0, max(Diff_End_x_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_x - std_End_x, mean_End_x + std_End_x, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End x', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_y, 0, max(Diff_End_y_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_y - std_End_y, mean_End_y + std_End_y, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End y', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
    mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.vlines(mean_End_z, 0, max(Diff_End_z_histY), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean')
    mp.axvspan(mean_End_z - std_End_z, mean_End_z + std_End_z, alpha = 0.2, color = red_cbf, label = '1$\\sigma$')
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('End z', fontsize = 'x-large')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.savefig('%s/Difference_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
        
    # create histograms for all directions for truth and reconstruction 
    z_bins = np.zeros(102, dtype = float)
    z_bins[0] = 11362
    for i in range(1, 102):
        if i < 40:
            z_bins[i] = z_bins[0] + i * 55.
        if i >= 40:
            z_bins[i] = 13507. + (i - 40) * 80.
    
    Reco_Start_x_hist, Reco_Start_x_bins = np.histogram(Reco_Start[:, 0], bins = 48, range = (-3520., 3520.))
    Reco_Start_y_hist, Reco_Start_y_bins = np.histogram(Reco_Start[:, 1], bins = 20, range = (-2949., 244.))    #+915
    Reco_Start_z_hist, Reco_Start_z_bins = np.histogram(Reco_Start[:, 2], bins = z_bins)
    Reco_End_x_hist, Reco_End_x_bins = np.histogram(Reco_End[:, 0], bins = Reco_Start_x_bins)
    Reco_End_y_hist, Reco_End_y_bins = np.histogram(Reco_End[:, 1], bins = Reco_Start_y_bins)   #+915
    Reco_End_z_hist, Reco_End_z_bins = np.histogram(Reco_End[:, 2], bins = Reco_Start_z_bins)
    Primary_True_Start_x_hist, Primary_True_Start_x_bins = np.histogram(Primary_True_Start[:, 0], bins = Reco_Start_x_bins)
    Primary_True_Start_y_hist, Primary_True_Start_y_bins = np.histogram(Primary_True_Start[:, 1], bins = Reco_Start_y_bins)
    Primary_True_Start_z_hist, Primary_True_Start_z_bins = np.histogram(Primary_True_Start[:, 2], bins = Reco_Start_z_bins)
    Primary_True_End_x_hist, Primary_True_End_x_bins = np.histogram(Primary_True_End[:, 0], bins = Reco_Start_x_bins)
    Primary_True_End_y_hist, Primary_True_End_y_bins = np.histogram(Primary_True_End[:, 1], bins = Reco_Start_y_bins)
    Primary_True_End_z_hist, Primary_True_End_z_bins = np.histogram(Primary_True_End[:, 2], bins = Reco_Start_z_bins)
    
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
    mp.fill_between(Primary_True_Start_x_histX, 0, Primary_True_Start_x_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_x_histX, 0, Reco_Start_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start x', fontsize = 'x-large')
    mp.savefig('%s/Track_Start_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_y_histX, Primary_True_Start_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_y_histX, Reco_Start_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_Start_y_histX, 0, Primary_True_Start_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_y_histX, 0, Reco_Start_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start y', fontsize = 'x-large')   
    mp.savefig('%s/Track_Start_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_z_histX, Primary_True_Start_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_z_histX, Reco_Start_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_Start_z_histX, 0, Primary_True_Start_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_z_histX, 0, Reco_Start_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('Start z', fontsize = 'x-large')
    mp.savefig('%s/Track_Start_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_x_histX, Primary_True_End_x_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_x_histX, Reco_End_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_x_histX, 0, Primary_True_End_x_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_x_histX, 0, Reco_End_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End x', fontsize = 'x-large')   
    mp.savefig('%s/Track_End_X_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_y_histX, Primary_True_End_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_y_histX, Reco_End_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_y_histX, 0, Primary_True_End_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_y_histX, 0, Reco_End_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End y', fontsize = 'x-large')
    mp.savefig('%s/Track_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_z_histX, Primary_True_End_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_z_histX, Reco_End_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_z_histX, 0, Primary_True_End_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_z_histX, 0, Reco_End_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('End z', fontsize = 'x-large')
    mp.savefig('%s/Track_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    ### plot difference in dependence of hit position
    # create 2d histograms
    #print(min(min(Primary_True_Start[:, 1]), min(Primary_True_End[:, 1])), max(max(Primary_True_Start[:, 1]), max(Primary_True_End[:, 1])), min(min(Reco_Start[:, 1]), min(Reco_End[:, 1])), max(max(Reco_Start[:, 1]), max(Reco_End[:, 1])))
    
    dependence_Start_x_hist, dependence_Start_x_binsX, dependence_Start_x_binsY = np.histogram2d(Primary_True_Start[:, 0], Diff_Start_x, bins = [Primary_True_Start_x_bins, Diff_Start_x_bins])
    dependence_Start_y_hist, dependence_Start_y_binsX, dependence_Start_y_binsY = np.histogram2d(Primary_True_Start[:, 1], Diff_Start_y, bins = [Primary_True_Start_y_bins, Diff_Start_y_bins])
    dependence_Start_z_hist, dependence_Start_z_binsX, dependence_Start_z_binsY = np.histogram2d(Primary_True_Start[:, 2], Diff_Start_z, bins = [Primary_True_Start_z_bins, Diff_Start_z_bins])
    dependence_End_x_hist, dependence_End_x_binsX, dependence_End_x_binsY = np.histogram2d(Primary_True_End[:, 0], Diff_End_x, bins = [Primary_True_End_x_bins, Diff_End_x_bins])
    dependence_End_y_hist, dependence_End_y_binsX, dependence_End_y_binsY = np.histogram2d(Primary_True_End[:, 1], Diff_End_y, bins = [Primary_True_End_y_bins, Diff_End_y_bins])
    dependence_End_z_hist, dependence_End_z_binsX, dependence_End_z_binsY = np.histogram2d(Primary_True_End[:, 2], Diff_End_z, bins = [Primary_True_End_z_bins, Diff_End_z_bins])
    
    cmap = cm.get_cmap('cividis');
    
    im = mp.pcolormesh(dependence_Start_x_binsX, dependence_Start_x_binsY, np.transpose(dependence_Start_x_hist), cmap = cmap);
    mp.xlabel('Start in TMS (X) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start x', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_X_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_y_binsX, dependence_Start_y_binsY, np.transpose(dependence_Start_y_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start y', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_Y_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_z_binsX, dependence_Start_z_binsY, np.transpose(dependence_Start_z_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('Start z', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/Start_Z_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_x_binsX, dependence_End_x_binsY, np.transpose(dependence_End_x_hist), cmap = cmap);
    mp.xlabel('End in TMS (X) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End x', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/End_X_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();
    
    im = mp.pcolormesh(dependence_End_y_binsX, dependence_End_y_binsY, np.transpose(dependence_End_y_hist), cmap = cmap);
    mp.xlabel('End in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End y', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/End_Y_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_z_binsX, dependence_End_z_binsY, np.transpose(dependence_End_z_hist), cmap = cmap);
    mp.xlabel('End in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('End z', fontsize = 'x-large')
    mp.colorbar(im);
    mp.savefig('%s/End_Z_real_and_difference.png' % out_dir, bbox_inches = 'tight');
    mp.close();

    # Track length
    # histogram
    
    # single histogram as well
    reco_tracklength_hist, reco_tracklength_bins = np.histogram(Reco_TrackLength / 10, bins = 40, range = (min(min(Reco_TrackLength / 10), min(True_TrackLength)), max(max(Reco_TrackLength / 10), max(True_TrackLength))))
    true_tracklength_hist, true_tracklength_bins = np.histogram(True_TrackLength, bins = reco_tracklength_bins)
    
    reco_tracklength_histX, reco_tracklength_histY = histogram_arr_handle(reco_tracklength_hist, reco_tracklength_bins)
    true_tracklength_histX, true_tracklength_histY = histogram_arr_handle(true_tracklength_hist, true_tracklength_bins)
    
    mp.plot(true_tracklength_histX, true_tracklength_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(reco_tracklength_histX, reco_tracklength_histY, color = orange_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(true_tracklength_histX, 0, true_tracklength_histY, alpha = 0.6, color = blue_cbf, hatch = '\\\\')
    mp.fill_between(reco_tracklength_histX, 0, reco_tracklength_histY, alpha = 0.6, color = orange_cbf, hatch = '//')
    mp.xlabel('Reco track length [$\\frac{g}{cm^2}$]')
    mp.ylabel('#')
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.savefig('%s/Reco_tracklength.png' % out_dir, bbox_inches = 'tight')
    mp.close()    
    
    tracklength_hist, tracklength_binsX, tracklength_binsY = np.histogram2d(True_TrackLength, Reco_TrackLength / 10, bins = [true_tracklength_bins, reco_tracklength_bins])
    
    im = mp.pcolormesh(tracklength_binsX, tracklength_binsY, np.transpose(tracklength_hist), cmap = cmap)
    mp.xlabel('True track length [$\\frac{g}{cm^2}$]')
    mp.ylabel('Reco track length [$\\frac{g}{cm^2}$]')
    mp.colorbar(im)
    mp.savefig('%s/Tracklength.png' % out_dir, bbox_inches = 'tight')
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
        if (true_slope * 18294 + true_intercept) < -2949 and Primary_True_End[i, 1] < (-2949 + 338.3):
            truth_exiting = True
        elif Primary_True_End[i, 2] < 17750:
            truth_stopping = True
        if (reco_slope * 18294 + reco_intercept) < -2949 and Reco_End[i, 1] < (-2949 + 338.3):
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
    print("Exiting: ", exiting_counter, " (", exiting_counter / len(Reco_Start) * 100, "%) | eff.: ", efficiency_stop * 100, " , pur.: ", purity_stop * 100, " , acc.: ", accuracy_stop * 100)
    print("Stopping: ", stopping_counter, " (", stopping_counter / len(Reco_Start) * 100, "%) | eff.: ", efficiency_exit * 100, " , pur.: ", purity_exit * 100, " , acc.: ", accuracy_exit * 100)
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
    v.get_label_by_id('11').set_text('%i' % exiting_counter)
    v.get_label_by_id('11').set_color('black')
    v.get_label_by_id('11').set_fontsize('x-large')
    v.get_patch_by_id('11').set_color('#f39e92')
    mp.title('Exiting')
    mp.savefig('%s/Exiting_venn.png' % out_dir, bbox_inches = 'tight')
    mp.close()
    
    return

def eff_pur_acc(overlap, truth, reco):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    input_filename = args.input_filename
    draw_performance(out_dir, input_filename)

