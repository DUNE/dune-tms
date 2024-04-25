import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll

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
    
    Position_TMS_Start = np.empty((n_events, 3), dtype = float) * -9999.
    Position_TMS_End = np.empty((n_events, 3), dtype = float) * -9999.
    Reco_Start = np.empty((n_events, 3), dtype = float) * -9999.
    Reco_End = np.empty((n_events, 3), dtype = float) * -9999.
    
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
            
            PositionTMSStart = np.frombuffer(true_event.PositionTMSStart, dtype = np.float32)
            PositionTMSEnd = np.frombuffer(true_event.PositionTMSEnd, dtype = np.float32)
            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            
            if PositionTMSStart[0] > -8000. and not StartPos.size == 0:
                Position_TMS_Start[i, 0] = PositionTMSStart[0]
                Position_TMS_Start[i, 1] = PositionTMSStart[1]
                Position_TMS_Start[i, 2] = PositionTMSStart[2]
                Reco_Start[i, 0] = StartPos[0]
                Reco_Start[i, 1] = StartPos[1]
                Reco_Start[i, 2] = StartPos[2]                
            if PositionTMSStart[0] > -8000. and not EndPos.size == 0:
                Position_TMS_End[i, 0] = PositionTMSEnd[0]
                Position_TMS_End[i, 1] = PositionTMSEnd[1]
                Position_TMS_End[i, 2] = PositionTMSEnd[2]
                Reco_End[i, 0] = EndPos[0]
                Reco_End[i, 1] = EndPos[1]
                Reco_End[i, 2] = EndPos[2]
    
    # subtract reconstruction from truth for all directions
    Diff_Start_x = Position_TMS_Start[:, 0] - Reco_Start[:, 0]
    Diff_Start_y = Position_TMS_Start[:, 1] - Reco_Start[:, 1]
    Diff_Start_z = Position_TMS_Start[:, 2] - Reco_Start[:, 2]
    Diff_End_x = Position_TMS_End[:, 0] - Reco_End[:, 0]
    Diff_End_y = Position_TMS_End[:, 1] - Reco_End[:, 1]
    Diff_End_z = Position_TMS_End[:, 2] - Reco_End[:, 2]

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
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('x')
    mp.savefig('performance_plots_NERSC/Difference_Start_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Difference_Start_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('z')   
    mp.savefig('performance_plots_NERSC/Difference_Start_Z_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('x')
    mp.savefig('performance_plots_NERSC/Difference_End_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Difference_End_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('z')   
    mp.savefig('performance_plots_NERSC/Difference_End_Z_hist.png', bbox_inches = 'tight')
    mp.close()    
    
    mp.plot(Diff_Start_x_histX, Diff_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_x_histX, 0, Diff_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('x')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Difference_Start_X_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_y_histX, Diff_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_y_histX, 0, Diff_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('y')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Difference_Start_Y_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_Start_z_histX, Diff_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_Start_z_histX, 0, Diff_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('z')
    mp.yscale('log')    
    mp.savefig('performance_plots_NERSC/Difference_Start_Z_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_x_histX, Diff_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_x_histX, 0, Diff_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('x')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Difference_End_X_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_y_histX, Diff_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_y_histX, 0, Diff_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('y')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Difference_End_Y_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Diff_End_z_histX, Diff_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2)
    mp.fill_between(Diff_End_z_histX, 0, Diff_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.xlabel('Difference truth – reconstruction')
    mp.ylabel('#')
    mp.title('z')
    mp.yscale('log')    
    mp.savefig('performance_plots_NERSC/Difference_End_Z_hist_log.png', bbox_inches = 'tight')
    mp.close()
        
    # create histograms for all directions for truth and reconstruction 
    Position_TMS_Start_x_hist, Position_TMS_Start_x_bins = np.histogram(Position_TMS_Start[:, 0], bins = 192)#, range = (-3400., 3400.))
    Position_TMS_Start_y_hist, Position_TMS_Start_y_bins = np.histogram(Position_TMS_Start[:, 1], bins = 100)#, range = (-4110., -910.))
    Position_TMS_Start_z_hist, Position_TMS_Start_z_bins = np.histogram(Position_TMS_Start[:, 2], bins = 100)#, range = (11360., 18320.))
    Position_TMS_End_x_hist, Position_TMS_End_x_bins = np.histogram(Position_TMS_End[:, 0], bins = 192)#, range = (-3400., 3400.))
    Position_TMS_End_y_hist, Position_TMS_End_y_bins = np.histogram(Position_TMS_End[:, 1], bins = 100)#, range = (-4110., -910.))
    Position_TMS_End_z_hist, Position_TMS_End_z_bins = np.histogram(Position_TMS_End[:, 2], bins = 100)#, range = (11360., 18320.))
    Reco_Start_x_hist, Reco_Start_x_bins = np.histogram(Reco_Start[:, 0], bins = 192)#, range = (-3400., 3400.))
    Reco_Start_y_hist, Reco_Start_y_bins = np.histogram(Reco_Start[:, 1], bins = 20)#, range = (-4110., -910.))
    Reco_Start_z_hist, Reco_Start_z_bins = np.histogram(Reco_Start[:, 2], bins = 100)#, range = (11360., 18320.))
    Reco_End_x_hist, Reco_End_x_bins = np.histogram(Reco_End[:, 0], bins = 192)#, range = (-3400., 3400.))
    Reco_End_y_hist, Reco_End_y_bins = np.histogram(Reco_End[:, 1], bins = 20)#, range = (-4110., -910.))
    Reco_End_z_hist, Reco_End_z_bins = np.histogram(Reco_End[:, 2], bins = 100)#, range = (11360., 18320.))
    
    # make the histograms usable
    Position_TMS_Start_x_histX, Position_TMS_Start_x_histY = histogram_arr_handle(Position_TMS_Start_x_hist, Position_TMS_Start_x_bins)
    Position_TMS_Start_y_histX, Position_TMS_Start_y_histY = histogram_arr_handle(Position_TMS_Start_y_hist, Position_TMS_Start_y_bins)
    Position_TMS_Start_z_histX, Position_TMS_Start_z_histY = histogram_arr_handle(Position_TMS_Start_z_hist, Position_TMS_Start_z_bins)
    Position_TMS_End_x_histX, Position_TMS_End_x_histY = histogram_arr_handle(Position_TMS_End_x_hist, Position_TMS_End_x_bins)
    Position_TMS_End_y_histX, Position_TMS_End_y_histY = histogram_arr_handle(Position_TMS_End_y_hist, Position_TMS_End_y_bins)
    Position_TMS_End_z_histX, Position_TMS_End_z_histY = histogram_arr_handle(Position_TMS_End_z_hist, Position_TMS_End_z_bins)
    Reco_Start_x_histX, Reco_Start_x_histY = histogram_arr_handle(Reco_Start_x_hist, Reco_Start_x_bins)
    Reco_Start_y_histX, Reco_Start_y_histY = histogram_arr_handle(Reco_Start_y_hist, Reco_Start_y_bins)
    Reco_Start_z_histX, Reco_Start_z_histY = histogram_arr_handle(Reco_Start_z_hist, Reco_Start_z_bins)
    Reco_End_x_histX, Reco_End_x_histY = histogram_arr_handle(Reco_End_x_hist, Reco_End_x_bins)
    Reco_End_y_histX, Reco_End_y_histY = histogram_arr_handle(Reco_End_y_hist, Reco_End_y_bins)
    Reco_End_z_histX, Reco_End_z_histY = histogram_arr_handle(Reco_End_z_hist, Reco_End_z_bins)
    
    # now plot this
    mp.plot(Position_TMS_Start_x_histX, Position_TMS_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_x_histX, Reco_Start_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_x_histX, 0, Position_TMS_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_x_histX, 0, Reco_Start_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('x')   
    mp.savefig('performance_plots_NERSC/Track_Start_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_Start_y_histX, Position_TMS_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_y_histX, Reco_Start_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_y_histX, 0, Position_TMS_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_y_histX, 0, Reco_Start_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('y')   
    mp.savefig('performance_plots_NERSC/Track_Start_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_Start_z_histX, Position_TMS_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_z_histX, Reco_Start_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_z_histX, 0, Position_TMS_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_z_histX, 0, Reco_Start_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('performance_plots_NERSC/Track_Start_Z_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_x_histX, Position_TMS_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_x_histX, Reco_End_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_x_histX, 0, Position_TMS_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_x_histX, 0, Reco_End_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('x')   
    mp.savefig('performance_plots_NERSC/Track_End_X_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_y_histX, Position_TMS_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_y_histX, Reco_End_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_y_histX, 0, Position_TMS_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_y_histX, 0, Reco_End_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('y')
    mp.savefig('performance_plots_NERSC/Track_End_Y_hist.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_z_histX, Position_TMS_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_z_histX, Reco_End_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_z_histX, 0, Position_TMS_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_z_histX, 0, Reco_End_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('z')
    mp.savefig('performance_plots_NERSC/Track_End_Z_hist.png', bbox_inches = 'tight')
    mp.close()

    mp.plot(Position_TMS_Start_x_histX, Position_TMS_Start_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_x_histX, Reco_Start_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_x_histX, 0, Position_TMS_Start_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_x_histX, 0, Reco_Start_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('x')
    mp.yscale('log')    
    mp.savefig('performance_plots_NERSC/Track_Start_X_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_Start_y_histX, Position_TMS_Start_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_y_histX, Reco_Start_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_y_histX, 0, Position_TMS_Start_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_y_histX, 0, Reco_Start_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('y')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Track_Start_Y_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_Start_z_histX, Position_TMS_Start_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_Start_z_histX, Reco_Start_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_Start_z_histX, 0, Position_TMS_Start_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_Start_z_histX, 0, Reco_Start_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('Start in TMS')
    mp.ylabel('#')
    mp.title('z')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Track_Start_Z_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_x_histX, Position_TMS_End_x_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_x_histX, Reco_End_x_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_x_histX, 0, Position_TMS_End_x_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_x_histX, 0, Reco_End_x_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('x')
    mp.yscale('log')    
    mp.savefig('performance_plots_NERSC/Track_End_X_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_y_histX, Position_TMS_End_y_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_y_histX, Reco_End_y_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_y_histX, 0, Position_TMS_End_y_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_y_histX, 0, Reco_End_y_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('y')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Track_End_Y_hist_log.png', bbox_inches = 'tight')
    mp.close()
    
    mp.plot(Position_TMS_End_z_histX, Position_TMS_End_z_histY, color = blue_cbf, linestyle = '--', linewidth = 2, label = 'truth')
    mp.plot(Reco_End_z_histX, Reco_End_z_histY, color = red_cbf, linestyle = ':', linewidth = 2, label = 'reco')
    mp.fill_between(Position_TMS_End_z_histX, 0, Position_TMS_End_z_histY, alpha = 0.6, hatch = '//', color = blue_cbf)
    mp.fill_between(Reco_End_z_histX, 0, Reco_End_z_histY, alpha = 0.6, hatch = '\\\\', color = red_cbf)
    mp.grid(True, linestyle = '--', alpha = 0.2)
    mp.legend(loc = 'best', fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
    mp.xlabel('End in TMS')
    mp.ylabel('#')
    mp.title('z')
    mp.yscale('log')
    mp.savefig('performance_plots_NERSC/Track_End_Z_hist_log.png', bbox_inches = 'tight')
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

