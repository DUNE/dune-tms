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
def draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons = False):
    if not os.path.exists(input_filename): raise ValueError(f"Cannor find input_filename {input_filename}")
    if readout_filename != "" and not os.path.exists(readout_filename): raise ValueError(f"Cannot find readout_filename {readout_filename}")
    if spill_number < -1: raise ValueError(f"Got spill_number = {spill_number}")
    if time_slice < -1: raise ValueError(f"Got time_slice = {time_slice}")
    
    # Make sure we read in the correct file and have the output directory
    use_readout = True
    if readout_filename == "": use_readout = False
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
    
    # Not used yet
    readout = None
    if use_readout:
        readout = ROOT.TChain("TMS")
        readout.Add(readout_filename)
        if not readout.GetEntries() > 0:
            print("Didnt't get any entries in TMS, are you sure the readout_filename is right?\n", readout_filename)
            
    max_n_spills = 10000 # TODO (old) add some meta info to output file with n spill info for file
    
    simplify_tracks = False
    
    spill_number_cache = dict()
    n_events = r.GetEntries()
    
    Reco_Start = np.ones((n_events, 3), dtype = float) * -9999.
    Reco_End = np.ones((n_events, 3), dtype = float) * -9999.
    Primary_True_Start = np.ones((n_events, 3), dtype = float) * -9999.
    Primary_True_End = np.ones((n_events, 3), dtype = float) * -9999.    
    
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
            
            if RecoTrackPrimaryParticleTruePositionTrackStart[0] > -8000. and not StartPos.size == 0:
                Reco_Start[i, 0] = StartPos[0]
                Reco_Start[i, 1] = StartPos[1]
                Reco_Start[i, 2] = StartPos[2]
                Primary_True_Start[i, 0] = RecoTrackPrimaryParticleTruePositionTrackStart[0]
                Primary_True_Start[i, 1] = RecoTrackPrimaryParticleTruePositionTrackStart[1]
                Primary_True_Start[i, 2] = RecoTrackPrimaryParticleTruePositionTrackStart[2]
            if RecoTrackPrimaryParticleTruePositionTrackEnd[0] > -8000. and not EndPos.size == 0:
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
    
    # subtract reconstruction from truth for all directions
    Diff_Start_x = Primary_True_Start[:, 0] - Reco_Start[:, 0]
    Diff_Start_y = Primary_True_Start[:, 1] - Reco_Start[:, 1]
    Diff_Start_z = Primary_True_Start[:, 2] - Reco_Start[:, 2]
    Diff_End_x = Primary_True_End[:, 0] - Reco_End[:, 0]
    Diff_End_y = Primary_True_End[:, 1] - Reco_End[:, 1]
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
    
    # plot
    mp.plot(Diff_Start_x_histX, Diff_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_x_histX, 0, Diff_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('x')
    mp.savefig('performance_plots_NERSC/Difference_Start_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Difference_Start_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('z')   
    mp.savefig('performance_plots_NERSC/Difference_Start_Z_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('x')
    mp.savefig('performance_plots_NERSC/Difference_End_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Difference_End_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction [mm]')
    mp.ylabel('#')
    mp.title('z')   
    mp.savefig('performance_plots_NERSC/Difference_End_Z_hist.png', bbox_inches = 'tight')
    mp.close()  
        
    # create histograms for all directions for truth and reconstruction 
    Reco_Start_x_hist, Reco_Start_x_bins = np.histogram(Reco_Start[:, 0], bins = 48)#, range = (-3400., 3400.))
    Reco_Start_y_hist, Reco_Start_y_bins = np.histogram(Reco_Start[:, 1], bins = 20)#, range = (-4110., -910.))
    Reco_Start_z_hist, Reco_Start_z_bins = np.histogram(Reco_Start[:, 2], bins = 100)#, range = (11360., 18320.))
    Reco_End_x_hist, Reco_End_x_bins = np.histogram(Reco_End[:, 0], bins = 48)#, range = (-3400., 3400.))
    Reco_End_y_hist, Reco_End_y_bins = np.histogram(Reco_End[:, 1], bins = 20)#, range = (-4110., -910.))
    Reco_End_z_hist, Reco_End_z_bins = np.histogram(Reco_End[:, 2], bins = 100)#, range = (11360., 18320.))
    Primary_True_Start_x_hist, Primary_True_Start_x_bins = np.histogram(Primary_True_Start[:, 0], bins = 48)
    Primary_True_Start_y_hist, Primary_True_Start_y_bins = np.histogram(Primary_True_Start[:, 1], bins = 20)
    Primary_True_Start_z_hist, Primary_True_Start_z_bins = np.histogram(Primary_True_Start[:, 2], bins = 100)
    Primary_True_End_x_hist, Primary_True_End_x_bins = np.histogram(Primary_True_End[:, 0], bins = 48)
    Primary_True_End_y_hist, Primary_True_End_y_bins = np.histogram(Primary_True_End[:, 1], bins = 20)
    Primary_True_End_z_hist, Primary_True_End_z_bins = np.histogram(Primary_True_End[:, 2], bins = 100)
    
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
    mp.title('x')   
    mp.savefig('performance_plots_NERSC/Track_Start_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_y_histX, Primary_True_Start_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_y_histX, Reco_Start_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_Start_y_histX, 0, Primary_True_Start_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_y_histX, 0, Reco_Start_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('y')   
    mp.savefig('performance_plots_NERSC/Track_Start_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_Start_z_histX, Primary_True_Start_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_z_histX, Reco_Start_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_Start_z_histX, 0, Primary_True_Start_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_Start_z_histX, 0, Reco_Start_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('performance_plots_NERSC/Track_Start_Z_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_x_histX, Primary_True_End_x_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_x_histX, Reco_End_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_x_histX, 0, Primary_True_End_x_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_x_histX, 0, Reco_End_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('x')   
    mp.savefig('performance_plots_NERSC/Track_End_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_y_histX, Primary_True_End_y_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_y_histX, Reco_End_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_y_histX, 0, Primary_True_End_y_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_y_histX, 0, Reco_End_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Track_End_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Primary_True_End_z_histX, Primary_True_End_z_histY, color = blue_cbf, linestyle = '-.', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_z_histX, Reco_End_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Primary_True_End_z_histX, 0, Primary_True_End_z_histY, color = blue_cbf, alpha = 0.6, hatch = '//')
    mp.fill_between(Reco_End_z_histX, 0, Reco_End_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS [mm]')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('performance_plots_NERSC/Track_End_Z_hist.png', bbox_inches = 'tight')
    mp.close()
    
    ### plot difference in dependence of hit position
    # create 2d histograms
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
    mp.title('x')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/Start_X_real_and_difference.png', bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_y_binsX, dependence_Start_y_binsY, np.transpose(dependence_Start_y_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('y')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/Start_Y_real_and_difference.png', bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_Start_z_binsX, dependence_Start_z_binsY, np.transpose(dependence_Start_z_hist), cmap = cmap);
    mp.xlabel('Start in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('z')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/Start_Z_real_and_difference.png', bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_x_binsX, dependence_End_x_binsY, np.transpose(dependence_End_x_hist), cmap = cmap);
    mp.xlabel('End in TMS (X) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('x')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/End_X_real_and_difference.png', bbox_inches = 'tight');
    mp.close();
    
    im = mp.pcolormesh(dependence_End_y_binsX, dependence_End_y_binsY, np.transpose(dependence_End_y_hist), cmap = cmap);
    mp.xlabel('End in TMS (Y) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('y')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/End_Y_real_and_difference.png', bbox_inches = 'tight');
    mp.close();

    im = mp.pcolormesh(dependence_End_z_binsX, dependence_End_z_binsY, np.transpose(dependence_End_z_hist), cmap = cmap);
    mp.xlabel('End in TMS (Z) [mm]')
    mp.ylabel('Difference truth - reco [mm]')
    mp.title('z')
    mp.colorbar(im);
    mp.savefig('performance_plots_NERSC/End_Z_real_and_difference.png', bbox_inches = 'tight');
    mp.close();
    
            
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
    parser.add_argument('--name', "-n", type = str, help = "The name of the output files. Will be <name>_<spill>_<slice>.png. Default = spill", default = "spill")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--spillnum', "-s", type = int, help = "The spill to draw. -1 for all", default = -1)
    parser.add_argument('--timeslice', "-t", type = int, help = "The time slice to draw. -1 for all", default = -1)
    parser.add_argument('--readout_filename', "-r", type = str, help = "(optional) A file with the raw readout.", default = "")
    parser.add_argument('--only_true_tms_muons', help = "Only draw events with true muons inside the TMS", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    readout_filename =  args.readout_filename
    only_true_tms_muons = args.only_true_tms_muons
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons)

