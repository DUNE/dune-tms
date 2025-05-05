import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll
from math import sqrt

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-v0_8-poster')

mp.rc('axes', labelsize = 12)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 12) # fontsize of the tick labels
mp.rc('ytick', labelsize = 12) # fontsize of the tick labels

cos_3 = 0.99863
sin_3 = 0.05234
tan_87 = 19.08114
tan_3 = 0.05241

### Widths of hits
delta_x = 0.018    # half a bar width
delta_y = 0.3383    # uncertainty from +/-3 degree tilted bars
delta_z = 0.025      # space of scintilator with air gap

MUON_MASS = 105.7  # MeV/c^2

tms_top_stereo = 0.29777
tms_bottom_stereo = -3.00223
tms_top_hybrid = 0.37177
tms_bottom_hybrid = -3.07623

### Function for upper limit of tilted bar 'hit'
def upper_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'VBar':  # assumption VBar is tilted in positive way and UBar then in negative
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y            
        if x < r:
            return_value = tan_3 * x - tan_3 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
        elif x >= r:
            return_value = -tan_87 * x + tan_87 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
    elif orientation_bar == 'UBar':
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y
        if x < r:
            return_value = tan_87 * x - tan_87 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
        elif x >= r:
            return_value = -tan_3 * x + tan_3 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value

### Function for lower limit of tilted bar 'hit'
def lower_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'VBar':
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return_value = -tan_87 * x + tan_87 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
        elif x >= r:
            return_value = tan_3 * x - tan_3 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
    elif orientation_bar == 'UBar':
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return_value = -tan_3 * x + tan_3 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value
        elif x >= r:
            return_value = tan_87 * x - tan_87 * r + s
            if return_value >= tms_top_stereo: return tms_top_stereo
            elif return_value <= tms_bottom_stereo: return tms_bottom_stereo
            else: return return_value

### Function for hits to appear in size
def hit_size(hit_x, hit_y, orientation, orientation_bar):
    if orientation == 'xy':
        if orientation_bar == 'XBar':
            size_array = np.zeros((2,2))
            size_array[0, 0] = hit_x / 1000.0 + delta_x
            size_array[0, 1] = hit_x / 1000.0 - delta_x
            if (hit_y / 1000.0 + delta_x) >= tms_top_hybrid:
                size_array[1, 0] = tms_top_hybrid
            else:
                size_array[1, 0] = hit_y / 1000.0 + delta_x
            if (hit_y / 1000.0 - delta_x) <= tms_bottom_hybrid:
                size_array[1, 1] = tms_bottom_hybrid
            else:
                size_array[1, 1] = hit_y / 1000.0 - delta_x
            return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]
        else:
            left_top = hit_x / 1000.0 - cos_3 * delta_x - sin_3 * delta_y
            right_bottom = hit_x / 1000.0 + cos_3 * delta_x + sin_3 * delta_y
            x_array = np.linspace(left_top, right_bottom, num = 50)
            return x_array, np.array([lower_limit(hit_x / 1000.0, hit_y / 1000.0, i, orientation_bar) for i in x_array]), np.array([upper_limit(hit_x / 1000.0, hit_y / 1000.0, i, orientation_bar) for i in x_array])
                       
    elif orientation == 'zy':
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x / 1000.0 + delta_z
        size_array[0, 1] = hit_x / 1000.0 - delta_z
        if orientation_bar == 'XBar':
            if (hit_y / 1000.0 + delta_x) >= tms_top_hybrid:
                size_array[1, 0] = tms_top_hybrid
            else:
                size_array[1, 0] = hit_y / 1000.0 + delta_x
            if (hit_y / 1000.0 - delta_x) <= tms_bottom_hybrid:
                size_array[1, 1] = tms_bottom_hybrid
            else:
                size_array[1, 1] = hit_y / 1000.0 - delta_x
        else:
            if (hit_y / 1000.0 + delta_y) >= tms_top_stereo:
                size_array[1, 0] = tms_top_stereo
            else:
                size_array[1, 0] = hit_y / 1000.0 + delta_y
            if (hit_y / 1000.0 - delta_y) <= tms_bottom_stereo:
                size_array[1, 1] = tms_bottom_stereo
            else:
                size_array[1, 1] = hit_y / 1000.0 - delta_y
        return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]   
             
    elif orientation == 'xz':
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x / 1000.0 + delta_z
        size_array[0, 1] = hit_x / 1000.0 - delta_z
        size_array[1, 0] = hit_y / 1000.0 + delta_x
        size_array[1, 1] = hit_y / 1000.0 - delta_x
        return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]        

### Actual function that loops through the spills
def draw_spill(out_dir, name, input_filename, spill_number, time_slice, histograms = False, lines2D = False, report_true_ke = False):
    if not os.path.exists(input_filename): raise ValueError(f"Cannor find input_filename {input_filename}")
    #if readout_filename != "" and not os.path.exists(readout_filename): raise ValueError(f"Cannot find readout_filename {readout_filename}")
    if spill_number < -1: raise ValueError(f"Got spill_number = {spill_number}")
    if time_slice < -1: raise ValueError(f"Got time_slice = {time_slice}")
 
    # Make sure we read in the correct file and have the output directory
    #use_readout = True
    #if readout_filename == "": use_readout = False
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

    if lines2D:
        b = ROOT.TChain("Line_Candidates")
        b.Add(input_filename)
        print("Line_Candidates: ", b.GetEntries())
        if not b.GetEntries() > 0:
            print("Didn't get any entries in Branch_Lines, are you sure the input_filename is right?\n", input_filename)

    DrawKalmanTrack = False
    if hasattr(r,"KalmanPos"):
        print("Kalman Filter info present in input file, will draw Kalman tracks.\n")
        DrawKalmanTrack = True
    
    truth = ROOT.TChain("Truth_Info")
    truth.Add(input_filename)
    if not truth.GetEntries() > 0:
        print("Didn't get any entries in Truth_Info, are you sure the input_filename is right?\n", input_filename)
    
    ## Not used yet
    #readout = None
    #if use_readout:
    #    readout = ROOT.TChain("TMS")
    #    readout.Add(readout_filename)
    #    if not readout.GetEntries() > 0:
    #        print("Didnt't get any entries in TMS, are you sure the readout_filename is right?\n", readout_filename)
            
    max_n_spills = 10000 # TODO (old) add some meta info to output file with n spill info for file
    
    start_spill = 0
    if spill_number != 0:
        start_spill = spill_number
        max_n_spills = spill_number + 1

    simplify_tracks = False
    
    spill_number_cache = dict()
    n_events = r.GetEntries()
    
    # First loop through all the slices and draw one overall spill
    for current_spill_number in range(start_spill, max_n_spills):
        for i in range(n_events):
            try:
                spill_number = spill_number_cache[i]
                event = None
                true_event = None
                if lines2D:
                    branch = None
            except KeyError:
                r.GetEntry(i)
                event = r
                truth.GetEntry(i)
                true_event = truth
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
                if lines2D:
                    b.GetEntry(i)
                    branch = b
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event == None:
                r.GetEntry(i)
                event = r
            if lines2D and branch == None:
                b.GetEntry(i)
                branch = b
            # Sync up the readout info if it's there. Note that it has one entry per spill, not timeslice
            #if readout != None: readout.GetEntry(current_spill_number)
            if true_event == None:
                truth.GetEntry(i)
                true_event = truth
            
            ### Check if a track exists in the event/spill, otherwise skip it
            nTracks = event.nTracks
            if nTracks <= 0: continue
                
            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])
            TrackHitPos = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            TrackHitBarType = np.frombuffer(event.TrackHitBarType, dtype = np.uint8)

            if histograms:
                TrackHitEnergies = np.frombuffer(event.TrackHitEnergies, dtype = np.float32)
            
            # I want the number of true hits.
            # nTrueHits = true_event.RecoTrackNHits
            # True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype=np.float32)

            if lines2D:
                nLinesU = branch.nLinesU
                nLinesV = branch.nLinesV
                nLinesX = branch.nLinesX
                nLinesY = branch.nLinesY
                nHitsU = np.frombuffer(branch.nHitsInTrackU, dtype = np.uint8)
                nHitsU = np.array([nHitsU[i] for i in range(0, nLinesU * 4, 4)])
                nHitsV = np.frombuffer(branch.nHitsInTrackV, dtype = np.uint8)
                nHitsV = np.array([nHitsV[i] for i in range(0, nLinesV * 4, 4)])
                nHitsX = np.frombuffer(branch.nHitsInTrackX, dtype = np.uint8)
                nHitsX = np.array([nHitsX[i] for i in range(0, nLinesX * 4, 4)])
                nHitsY = np.frombuffer(branch.nHitsInTrackY, dtype = np.uint8)
                nHitsY = np.array([nHitsY[i] for i in range(0, nLinesY * 4, 4)])
                TrackHitPosU = np.frombuffer(branch.TrackHitPosU, dtype = np.float32)
                TrackHitPosV = np.frombuffer(branch.TrackHitPosV, dtype = np.float32)
                TrackHitPosX = np.frombuffer(branch.TrackHitPosX, dtype = np.float32)
                TrackHitPosY = np.frombuffer(branch.TrackHitPosY, dtype = np.float32)

            # Kalman stuff
            nKalmanNodes = np.frombuffer(event.nKalmanNodes, dtype = np.uint8)
            KalmanPos = np.frombuffer(event.KalmanPos, dtype = np.float32)
            KalmanTruePos = np.frombuffer(event.KalmanTruePos, dtype = np.float32)

            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)            
            #Direction = np.frombuffer(event.Direction, dtype = np.float32)
            
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
            
            for j in range(nTracks):
                ### Create subplots
                fig = mp.figure(constrained_layout = False)
                if histograms:
                    gs = fig.add_gridspec(ncols=2, nrows=3, hspace = 0.3, wspace = 0.0)
                    x_y = fig.add_subplot(gs[0, 0])
                    z_y = fig.add_subplot(gs[1, 0])
                    x_z = fig.add_subplot(gs[0:2, 1:])
                    time = fig.add_subplot(gs[2, 0])
                    energy = fig.add_subplot(gs[2, 1:])
                else:
                    gs = fig.add_gridspec(2, 2, hspace = 0.25, wspace = 0.15)
                    x_y = fig.add_subplot(gs[0, 0])
                    z_y = fig.add_subplot(gs[1, 0])
                    x_z = fig.add_subplot(gs[0:, 1:])
    
                ### Set labels and ticks
                x_y.set(xlabel = 'x [m]', ylabel = 'y [m]', xticks = [3, 2, 1, 0, -1, -2, -3, -4], yticks = [-3, -2, -1, 0])
                z_y.set(xlabel = 'z [m]', ylabel = 'y [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0])
                x_z.set(xlabel = 'z [m]', ylabel = 'x [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0, 1, 2, 3])
                x_y.text(3.79, -2, 'front view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                z_y.text(18.6, -2, 'side view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                x_z.text(18.6, -0.5, 'top view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                if histograms:
                    time.set(xlabel = 'Time [ns]', ylabel = 'Total Hit E [MeV]')
                    energy.set(xlabel = 'Hit E [MeV]', ylabel = 'N Hits')
            
                ### Set TMS name
                x_y.text(-3.73, tms_top_hybrid, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8) #0.1
                z_y.text(11.2, tms_top_hybrid, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8) #0.1
                x_z.text(11.2, 3.79, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)

                ### Position plots efficient/nice in subplots
                x_z.axis('equal')
                x_z.axes.set_box_aspect(1)
                x_z.axes.set_anchor('W')
                z_y.axis('equal')
                z_y.axes.set_box_aspect(0.5)
                z_y.axes.set_anchor('NW')
                x_y.axis('equal')
                x_y.axes.set_box_aspect(0.5)
                x_y.axes.set_anchor('SW')
                if histograms:
                    time.axes.set_box_aspect(0.5)
                    time.axes.set_anchor('W')
                    energy.axes.set_box_aspect(0.5)
                    energy.axes.set_anchor('C')

                ### Put in outlines of scintillator parts
                x_z.hlines(-3.73, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')   # outer steel plate
                x_z.hlines(3.73, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')    # outer steel plate
                x_z.hlines(-1.87, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')   # inner steel plate
                x_z.hlines(0, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')       # middle of steel
                x_z.hlines(1.87, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')    # inner steel plate
                x_z.vlines(11.124, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.vlines(18.544, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.vlines(14.460, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5))) # to thick steel
                x_z.vlines(17.520, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5))) # to double thick steel

                z_y.hlines(tms_bottom_hybrid, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.hlines(tms_top_hybrid, 11.124, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.vlines(11.124, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.vlines(18.544, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.vlines(14.460, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5)))   # to thick steel
                z_y.vlines(17.520, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5)))   # to double thick steel

                x_y.hlines(tms_bottom_hybrid, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.hlines(tms_top_hybrid, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(-3.73, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(3.73, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(-1.87, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(0, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(1.87, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')


                # This isn't working with current file, Aug 2024. File doesn't have true_event.RecoTrackNHits
                # # loop through those true hits in the TMS
                # for hit in range(nTrueHits[i]):
                #     thit_x = True_Hits[i * 600 + hit * 4 + 0]
                #     thit_y = True_Hits[i * 600 + hit * 4 + 1]
                #     thit_z = True_Hits[i * 600 + hit * 4 + 2]
                #
                #     # print("THit", thit_x, thit_y, thit_z)
                #
                #     if thit_z < 11000.: continue
                #
                #     x_z.scatter(thit_z / 1000.0, thit_x / 1000.0, color=black_cbf, marker='.', alpha=0.5)
                #     z_y.scatter(thit_z / 1000.0, thit_y / 1000.0, color=black_cbf, marker='.', alpha=0.5)
                #     x_y.scatter(thit_x / 1000.0, thit_y / 1000.0, color=black_cbf, marker='.', alpha=0.5)
                #
                #     output_filename_thits = os.path.join(out_dir, f"{name}_truth_{current_spill_number:03d}")
                #     mp.savefig(output_filename_thits + ".png", bbox_inches='tight')
                #     mp.close()
                
                if histograms:
                    times = np.zeros(nHits[j], dtype = float)
                    energies = np.zeros(nHits[j], dtype = float)
                ### Hit positions of the hits in the track
                for hit in range(nHits[j]):
                    hit_x = TrackHitPos[j*800 + hit*4 + 0]
                    hit_y = TrackHitPos[j*800 + hit*4 + 1]
                    hit_z = TrackHitPos[j*800 + hit*4 + 2]
                    if histograms:
                        times[hit] = TrackHitPos[j*800 + hit*4 + 3]
                        energies[hit] = TrackHitEnergies[j*200 + hit]
                
                    #print('(%s)' %check_orientation(int(hit_z)), hit_x, hit_y, hit_z)
                
                    #temporary fix
                    if hit_z < 11000.: continue 
                    if np.abs(hit_x) > 10000. or np.abs(hit_y) > 10000. or np.abs(hit_z) > 20000.: continue 
                    
                    orientation_bar = check_orientation(TrackHitBarType[j*800 + hit*4])
                    if not lines2D:
                        color_cbf = red_cbf
                        if orientation_bar == 'VBar':
                            color_cbf = blue_cbf
                        elif orientation_bar == 'XBar':
                            color_cbf = black_cbf
                        elif orientation_bar == 'YBar':
                            color_cbf = blue_cbf
                            orientation_bar = 'XBar'
                    else:
                        color_cbf = black_cbf
                    
                    # check if gap with two hits successively have the same BarType
                    helper = 0
                    if TrackHitPos[j*800 + (hit + 1)*4 + 2] == hit_z: helper = 1    # if this is the case set helper to 1, to check from the next hit
                    # check for close by X hit for V hits
                    if orientation_bar == 'VBar':
                        if hit + helper + 2 < nHits[j]:
                            if check_orientation(TrackHitBarType[j*800 + (hit + helper + 2)*4]) == 'XBar':
                                # if close by X hit, set BarType for plotting (not color) to X type
                                orientation_bar = 'XBar'
                        elif hit + helper + 1 < nHits[j]:
                            if check_orientation(TrackHitBarType[j*800 + (hit + helper + 1)*4]) == 'XBar':
                                orientation_bar = 'XBar'
                        else:
                            if check_orientation(TrackHitBarType[j*800 + (hit - helper - 1)*4]) == 'XBar' or check_orientation(TrackHitBarType[j*800 + (hit - helper - 2)*4]) == 'XBar':
                                orientation_bar = 'XBar'
                    # check fro close by X hit for U hits
                    if orientation_bar == 'UBar':
                        if hit + helper + 1 < nHits[j]:
                            if check_orientation(TrackHitBarType[j*800 + (hit + helper + 1)*4]) == 'XBar':
                                orientation_bar = 'XBar'
                        elif hit + helper + 2 < nHits[j]:
                            if check_orientation(TrackHitBarType[j*800 + (hit + helper + 2)*4]) == 'XBar':
                                orientation_bar = 'XBar'
                        else:
                            if check_orientation(TrackHitBarType[j*800 + (hit - helper - 2)*4]) == 'XBar' or check_orientation(TrackHitBarType[j*800 + (hit - helper - 1)*4]) == 'XBar':
                                orientation_bar = 'XBar'

                    if hit + 3 >= nHits[j]:
                        x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', orientation_bar), color = color_cbf, label = 'hit area %s' % check_orientation(TrackHitBarType[j*800 + hit*4]))
                    else:
                        x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', orientation_bar), color = color_cbf)
                    z_y.fill_between(*hit_size(hit_z, hit_y, 'zy', orientation_bar), color = color_cbf)
                    x_y.fill_between(*hit_size(hit_x, hit_y, 'xy', orientation_bar), color = color_cbf, alpha = 0.5, linewidth = 0.5)
        
                if DrawKalmanTrack:
                    print("Track: ", j, "\t Hits: ", nHits[j], "\t Nodes: ", nKalmanNodes[j])

                    prev_kal_x = -1E100
                    prev_kal_y = -1E100
                    prev_kal_z = -1E100
                    kal_x = np.zeros(nKalmanNodes[j])
                    kal_y = np.zeros(nKalmanNodes[j])
                    kal_z = np.zeros(nKalmanNodes[j])
                    kal_true_x = np.zeros(nKalmanNodes[j])
                    kal_true_y = np.zeros(nKalmanNodes[j])

                    for node in range(nKalmanNodes[j]):
                        kal_x[node] = KalmanPos[j*600 + node*3 + 0]/1000.0 # from mm to m
                        kal_y[node] = KalmanPos[j*600 + node*3 + 1]/1000.0
                        kal_z[node] = KalmanPos[j*600 + node*3 + 2]/1000.0
                        kal_true_x[node] = KalmanTruePos[j*600 + node*3 + 0]/1000.0 # from mm to m
                        kal_true_y[node] = KalmanTruePos[j*600 + node*3 + 1]/1000.0

                    x_z.plot(kal_z[0:], kal_x[0:], ls=':', lw = 1.3, color = green_cbf, label = 'Kalman reco')
                    z_y.plot(kal_z[0:], kal_y[0:], ls=':', lw = 1.3, color = green_cbf)
                    x_y.plot(kal_x[0:], kal_y[0:], ls=':', lw = 1.3, color = green_cbf)

                    x_z.plot(kal_z[0:], kal_true_x[0:], ls='--', lw = 1.3, color = magenta_cbf, label = 'Kalman true')
                    z_y.plot(kal_z[0:], kal_true_y[0:], ls='--', lw = 1.3, color = magenta_cbf)
                    x_y.plot(kal_true_x[0:], kal_true_y[0:], ls='--', lw = 1.3, color = magenta_cbf)

                if histograms:
                    # create hit times histogram and plot
                    # weighted with hit energy as done by Jeffrey at https://github.com/DUNE/dune-tms/blob/kleykamp_validation/scripts/Reco/draw_spill.py#L188
                    time.hist(times % 1.2e9, bins = int(max(times % 1.2e9) - min(times % 1.2e9)),color = black_cbf, align = 'mid', weights = energies)  
                    # create hit energies histogram and plot
                    energy.hist(energies, bins = int(max(energies) * 3),color = black_cbf, align = 'mid')
                
                if lines2D:
                    if nLinesU != 0:
                        Uhits_x = np.zeros(nHitsU[j])
                        Uhits_z = np.zeros(nHitsU[j])
                        for hit in range(nHitsU[j]):
                            Uhits_z[hit] = TrackHitPosU[j*400 + hit*2 + 0] / 1000.0
                            Uhits_x[hit] = TrackHitPosU[j*400 + hit*2 + 1] / 1000.0
                            #print(Uhits_x[hit], Uhits_z[hit])
                        x_z.plot(Uhits_z, Uhits_x, ls = '-', lw = 1.3, color = red_cbf, label = '2D U')

                    if nLinesV != 0:
                        Vhits_x = np.zeros(nHitsV[j])
                        Vhits_z = np.zeros(nHitsV[j])
                        for hit in range(nHitsV[j]):
                            Vhits_z[hit] = TrackHitPosV[j*400 + hit*2 + 0] / 1000.0
                            Vhits_x[hit] = TrackHitPosV[j*400 + hit*2 + 1] / 1000.0
                        x_z.plot(Vhits_z, Vhits_x, ls = '-', lw = 1.3, color = green_cbf, label = '2D V')

                    if nLinesX != 0:
                        Xhits_y = np.zeros(nHitsX[j])   # TODO how are X simple tracks saved? Is the same iteration valid (j) or not for hybrid cases???
                        Xhits_z = np.zeros(nHitsX[j])
                        for hit in range(nHitsX[j]):
                            Xhits_z[hit] = TrackHitPosX[j*400 + hit*2 + 0] / 1000.0
                            Xhits_y[hit] = TrackHitPosX[j*400 + hit*2 + 1] / 1000.0
                        z_y.plot(Xhits_z, Xhits_y, ls = '-', lw = 1.3, color = blue_cbf, label = '2D X')

                    if nLinesY != 0:
                        Yhits_x = np.zeros(nHitsY[j])
                        Yhits_z = np.zeros(nHitsY[j])
                        for hit in range(nHitsY[j]):
                            Yhits_z[hit] = TrackHitPosY[j*400 + hit*2 + 0] / 1000.0
                            Yhits_x[hit] = TrackHitPosY[j*400 + hit*2 + 1] / 1000.0
                        x_z.plot(Yhits_z, Yhits_x, ls = '-', lw = 1.3, color = red_cbf, label = '2D Y')

                ### Track start
                #temporary fix
                if not (StartPos[j*4 + 2] < 11000.): 
                
                    if not StartPos[j*4 + 1] == 0.0:
                        orientation_bar = check_orientation(TrackHitBarType[j*800 + (nHits[j] - 1)*4])
                        if not orientation_bar == 'XBar':
                            if orientation_bar == 'VBar':
                                if check_orientation(TrackHitBarType[j*800 + (nHits[j] - 1 - 1)*4]) == 'XBar':
                                    orientation_bar = 'XBar'
                            if orientation_bar == 'UBar':
                                if check_orientation(TrackHitBarType[j*800 + (nHits[j] - 1 - 2)*4]) == 'XBar':
                                    orientation_bar = 'XBar'
                            if orientation_bar == 'YBar': orientation_bar = 'XBar'
                        x_z.fill_between(*hit_size(StartPos[j*4 + 2], StartPos[j*4 + 0], 'xz', orientation_bar), color = green_cbf, label = 'Start/End reco')
                        z_y.fill_between(*hit_size(StartPos[j*4 + 2], StartPos[j*4 + 1], 'zy', orientation_bar), color = green_cbf)
                        x_y.fill_between(*hit_size(StartPos[j*4 + 0], StartPos[j*4 + 1], 'xy', orientation_bar), color = green_cbf, alpha = 0.5, linewidth = 0.5)

            
                ### Track end               
                #temporary fix
                if not (EndPos[j*4 + 2] < 11000.):  
    
                    if not EndPos[j*4 + 1] == 0.0:
                        orientation_bar = check_orientation(TrackHitBarType[j*800 + 0])
                        if not orientation_bar == 'XBar':
                            if orientation_bar == 'VBar':
                                if check_orientation(TrackHitBarType[j*800 + 2*4]) == 'XBar':
                                    orientation_bar = 'XBar'
                            if orientation_bar == 'UBar':
                                if check_orientation(TrackHitBarType[j*800 + 1*4]) == 'XBar':
                                    orientation_bar = 'XBar'
                            if orientation_bar == 'YBar': orientation_bar = 'XBar'
                        x_z.fill_between(*hit_size(EndPos[j*4 + 2], EndPos[j*4 + 0], 'xz', orientation_bar), color = green_cbf)
                        z_y.fill_between(*hit_size(EndPos[j*4 + 2], EndPos[j*4 + 1], 'zy', orientation_bar), color = green_cbf)
                        x_y.fill_between(*hit_size(EndPos[j*4 + 0], EndPos[j*4 + 1], 'xy', orientation_bar), color = green_cbf, alpha = 0.5, linewidth = 0.5)
                
                ### Track direction
              
                #temporary fix
                # Add check on DrawKalmanTrack so we draw the true kalman info instead of a line
                if not DrawKalmanTrack and not (StartPos[j*4 + 2] < 11000. or EndPos[j*4 + 2] < 11000.): 

                    if not StartPos[j*4 + 1] == 0.0 or EndPos[j*4 + 1] == 0.0:
                        x_z.plot([StartPos[j*4 + 2] / 1000.0, EndPos[j*4 + 2] / 1000.0], [StartPos[j*4 + 0] / 1000.0, EndPos[j*4 + 0] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--', label = 'Direction')
                        z_y.plot([StartPos[j*4 + 2] / 1000.0, EndPos[j*4 + 2] / 1000.0], [StartPos[j*4 + 1] / 1000.0, EndPos[j*4 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
                        x_y.plot([StartPos[j*4 + 0] / 1000.0, EndPos[j*4 + 0] / 1000.0], [StartPos[j*4 + 1] / 1000.0, EndPos[j*4 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
                
                if RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] > 11000.: 
                    x_z.scatter(RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] / 1000.0, c = magenta_cbf, marker = '2', alpha = 0.5, label = 'Start true')
                    x_z.scatter(RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] / 1000.0, c = magenta_cbf, marker = '1', alpha = 0.5, label = 'End true')
                    
                    z_y.scatter(RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '2', alpha = 0.5)
                    z_y.scatter(RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '1', alpha = 0.5)
                    
                    x_y.scatter(RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '2', alpha = 0.5)
                    x_y.scatter(RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '1', alpha = 0.5)


                # Write the True Muon KE to each spill plot.
                if report_true_ke:
                    for idx, pdg in enumerate(true_event.PDG):
                        if pdg != abs(13): continue
    
                        muon_ke_lar = true_event.Muon_TrueKE / 1000.0
                        p_tms_start = ROOT.TVector3(truth.MomentumTMSStart[4 * idx], truth.MomentumTMSStart[4 * idx + 1], truth.MomentumTMSStart[4 * idx + 2])
                        muon_ke_tms_start = sqrt(p_tms_start.Mag2() + MUON_MASS ** 2) - MUON_MASS
                        muon_ke_tms_start /= 1000.0
                        x_z.text(11, 4, f'Muon KE at birth (LAr): {muon_ke_lar:.2f} GeV', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                        x_z.text(11, 5, f'Muon KE entering TMS: {muon_ke_tms_start:.2f} GeV', fontsize = 12, fontweight = 'bold', color = orange_cbf)
    
                        if muon_ke_tms_start > 5.0 or muon_ke_lar > 5.0:  # GeV
                            print(f'Event: {i}, Spill {spill_number}, Muon KE at birth (LAr): {muon_ke_lar}, Muon KE entering TMS: {muon_ke_tms_start}, GeV.')
    
                # add a legend
                fig.legend(loc = 7, fontsize = 'x-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
                fig.tight_layout()
                fig.subplots_adjust(right = 0.84)
                # save plot
                output_filename = os.path.join(out_dir, f"{name}_{current_spill_number:03d}_{i:03d}_{j:02d}")
                print("plotted ", output_filename)
                mp.savefig(output_filename + ".png", bbox_inches = 'tight')
                mp.close()
        
    return

### This is for plotting the hits according to their different orientations
def check_orientation(BarType):#hit_z):
    if BarType == 0:
        return 'XBar'
    elif BarType == 2:
        return 'UBar'
    elif BarType == 3:
        return 'VBar'
    elif BarType == 1:
        return 'YBar'
#    return layer_dict["%s" % hit_z]

### Dictionary that after calculate_layers contains for each z-coordinate the orientation str
first_z = 11185
layer_dict = { "%s" % 11133 : "VBar" }
layer_dict.update({ "%s" % first_z : "XBar" })
        
def calculate_layers(Xlayers):
    increment = 2
    if Xlayers:
        increment = 3
    thin_layers = 50
    thick_layers = 34
    double_layers = 9
    # Calculate the z position for each layer for the thin section
    for i in range(thin_layers):
        hit_z = first_z + i * 65
        # even layers
        if (((hit_z - first_z) / 65) % increment) == 0:
            layer_dict.update({ "%s" % hit_z : "UBar" })
        # odd layers
        elif (((hit_z - first_z) / 65) % increment) == 1:
            layer_dict.update({ "%s" % hit_z : "VBar" })
        # x layers
        if Xlayers:
            if (((hit_z - first_z) / 65) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "XBar" })
            # even layers
            if (((hit_z - first_z) / 65) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "UBar" })
            # odd layers
            elif (((hit_z - first_z) / 65) % increment) == 2:
                layer_dict.update({ "%s" % hit_z : "VBar" })
    # Calculate the z position for each layer for the thick section
    start_thick = first_z + thin_layers * 65
    for i in range(thick_layers):
        hit_z = start_thick + i * 90
        # even layers
        if (((hit_z - start_thick) / 90) % increment) == 0:
            layer_dict.update({ "%s" % hit_z : "UBar" })
        # odd layers
        elif (((hit_z - start_thick) / 90) % increment) == 1:
            layer_dict.update({ "%s" % hit_z : "VBar" })
        # x layers
        if Xlayers:
            if (((hit_z - start_thick) / 90) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "VBar" })
            # even layers
            elif (((hit_z - start_thick) / 90) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "XBar" })
            # odd layers
            elif (((hit_z - start_thick) / 90) % increment) == 2:
                layer_dict.update({ "%s" % hit_z : "UBar" })
    # Calculate the z position for each layer for the double section
    start_double = first_z + thin_layers * 65 + thick_layers * 90
    for i in range(double_layers):
        hit_z = start_double + i * 130
        # even layers
        if (((hit_z - start_double) / 130) % increment) == 0:
            layer_dict.update({ "%s" % hit_z : "UBar" })
        # odd layers
        elif (((hit_z - start_double) / 130) % increment) == 1:
            layer_dict.update({ "%s" % hit_z : "VBar" })
        # x layers
        if Xlayers:
            if (((hit_z - start_double) / 130) % increment) == 0:
                layer_dict.update({ "%s" % hit_z : "XBar" })
            # even layers
            elif (((hit_z - start_double) / 130) % increment) == 1:
                layer_dict.update({ "%s" % hit_z : "UBar" })
            # odd layers
            elif (((hit_z - start_double) / 130) % increment) == 2:
                layer_dict.update({ "%s" % hit_z : "VBar" })
                
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--name', "-n", type = str, help = "The name of the output files. Will be <name>_<spill>_<slice>.png. Default = spill", default = "spill")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--spillnum', "-s", type = int, help = "The spill to draw. -1 for all", default = -1)
    parser.add_argument('--timeslice', "-t", type = int, help = "The time slice to draw. -1 for all", default = -1)
    #parser.add_argument('--readout_filename', "-r", type = str, help = "(optional) A file with the raw readout.", default = "")
    parser.add_argument('--report_true_ke', help = "Add the true KE of muon to plot.", action = argparse.BooleanOptionalAction)
    #parser.add_argument('--Xlayers', "-X", help = "Does the geometry use X (90 degree orientated) scintillator layers? Yes -> --Xlayers, No -> --no-Xlayers", action = argparse.BooleanOptionalAction)
    parser.add_argument('--hists', "-H", help = "Plot hit times and energies histogram. Yes -> --hists", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Lines2D', "-l", help = "Plot low level 2D lines for single orientations (debugging). Yes -> --Lines2D", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    #readout_filename =  args.readout_filename
    report_true_ke = args.report_true_ke
    #Xlayers = args.Xlayers
    #print(Xlayers)
    #calculate_layers(Xlayers)
    #print(layer_dict)
    histograms = args.hists
    lines2D = args.Lines2D
    
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, histograms, lines2D, report_true_ke)

