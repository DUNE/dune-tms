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
delta_z = 0.025      # space of scintilattor with air gap

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
def hit_size(hit_x, hit_y, orientation, hit_z):
    if orientation == 'xy':
        orientation_bar = check_orientation(int(hit_z))
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
        orientation_bar = check_orientation(int(hit_z))
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
def draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, report_true_ke = False):
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

    DrawKalmanTrack = False
    #if hasattr(r,"KalmanPos"):
    #    print("Kalman Filter info present in input file, will draw Kalman tracks.\n")
    #    DrawKalmanTrack = True
    
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
    
    # First loop through all the slices and draw one overall spill
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
            # Sync up the readout info if it's there. Note that it has one entry per spill, not timeslice
            if readout != None: readout.GetEntry(current_spill_number)
            if true_event == None:
                truth.GetEntry(i)
                true_event = truth
            
            ### Check if a track exists in the event/spill, otherwise skip it
            nTracks = event.nTracks
            if nTracks <= 0: continue
    
#            ### Create subplots
#            fig = mp.figure(constrained_layout = False)
#            gs = fig.add_gridspec(2, 2, hspace = 0.25, wspace = 0.15)
#            x_y = fig.add_subplot(gs[0, 0])
#            z_y = fig.add_subplot(gs[1, 0])
#            x_z = fig.add_subplot(gs[0:, 1:])
    
#            ### Set labels and ticks
#            x_y.set(xlabel = 'x [m]', ylabel = 'y [m]', xticks = [4, 3, 2, 1, 0, -1, -2, -3, -4], yticks = [-3, -2, -1, 0])
#            z_y.set(xlabel = 'z [m]', ylabel = 'y [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0])
#            x_z.set(xlabel = 'z [m]', ylabel = 'x [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0, 1, 2, 3])
#            x_y.text(3.55, -2, 'front view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
#            z_y.text(18.6, -2, 'side view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
#            x_z.text(18.6, -0.5, 'top view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            
#            ### Set TMS name
#            x_y.text(-3.45, 0.1, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
#            z_y.text(11.22, 0.1, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
#            x_z.text(11.2, 3.55, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
    
#            ### Position plots efficient/nice in subplots
#            x_z.axis('equal')
#            x_z.axes.set_box_aspect(1)
#            z_y.axis('equal')
#            z_y.axes.set_box_aspect(0.5)
#            z_y.axes.set_anchor('NW')
#            x_y.axis('equal')
#            x_y.axes.set_box_aspect(0.5)
#            x_y.axes.set_anchor('SW')
    
#            ### Put in outlines of scintillator parts
#            x_z.hlines(-3.49, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.hlines(3.49, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.hlines(-1.75, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.hlines(0, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.hlines(1.75, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.vlines(11.176, -3.49, 3.52, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_z.vlines(18.544, -3.49, 3.52, color = orange_cbf, linewidth = 1, linestyle = ':')
    
#            z_y.hlines(tms_bottom_hybrid, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            z_y.hlines(tms_top_hybrid, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
#            z_y.vlines(11.176, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
#            z_y.vlines(18.544, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
    
#            x_y.hlines(tms_bottom_hybrid, -3.49, 3.49, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.hlines(tms_top_hybrid, -3.49, 3.49, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.vlines(-3.49, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.vlines(3.49, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.vlines(-1.75, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.vlines(0, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
#            x_y.vlines(1.75, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
            
            #print("number of tracks: ", nTracks)
            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])
            #print("number of hits: ", nHits)
            TrackHitPos = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            
            # I want the number of true hits.
            # nTrueHits = true_event.RecoTrackNHits
            # True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype=np.float32)

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
                gs = fig.add_gridspec(2, 2, hspace = 0.25, wspace = 0.15)
                x_y = fig.add_subplot(gs[0, 0])
                z_y = fig.add_subplot(gs[1, 0])
                x_z = fig.add_subplot(gs[0:, 1:])
    
                ### Set labels and ticks
                x_y.set(xlabel = 'x [m]', ylabel = 'y [m]', xticks = [4, 3, 2, 1, 0, -1, -2, -3, -4], yticks = [-3, -2, -1, 0])
                z_y.set(xlabel = 'z [m]', ylabel = 'y [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0])
                x_z.set(xlabel = 'z [m]', ylabel = 'x [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0, 1, 2, 3])
                x_y.text(3.79, -2, 'front view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                z_y.text(18.6, -2, 'side view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                x_z.text(18.6, -0.5, 'top view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            
                ### Set TMS name
                x_y.text(-3.69, 0.1, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
                z_y.text(11.22, 0.1, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
                x_z.text(11.2, 3.79, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)
    
                ### Position plots efficient/nice in subplots
                x_z.axis('equal')
                x_z.axes.set_box_aspect(1)
                z_y.axis('equal')
                z_y.axes.set_box_aspect(0.5)
                z_y.axes.set_anchor('NW')
                x_y.axis('equal')
                x_y.axes.set_box_aspect(0.5)
                x_y.axes.set_anchor('SW')
    
                ### Put in outlines of scintillator parts
                x_z.hlines(-3.73, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.hlines(3.73, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.hlines(-1.75, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.hlines(0, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.hlines(1.75, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.vlines(11.176, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_z.vlines(18.544, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
    
                z_y.hlines(tms_bottom_hybrid, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.hlines(tms_top_hybrid, 11.176, 18.544, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.vlines(11.176, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                z_y.vlines(18.544, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
    
                x_y.hlines(tms_bottom_hybrid, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.hlines(tms_top_hybrid, -3.73, 3.73, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(-3.73, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(3.73, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(-1.75, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(0, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                x_y.vlines(1.75, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')


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

                ### Hit positions of the hits in the track
                for hit in range(nHits[j]):

                    hit_x = TrackHitPos[j*600 + hit*3 + 0]
                    hit_y = TrackHitPos[j*600 + hit*3 + 1]
                    hit_z = TrackHitPos[j*600 + hit*3 + 2]
                
                    #print('(%s)' %check_orientation(int(hit_z)), hit_x, hit_y, hit_z)
                
                    #temporary fix
                    if hit_z < 11000.: continue #hit_y > -2000.0 or 
                    if np.abs(hit_x) > 10000. or np.abs(hit_y) > 10000. or np.abs(hit_z) > 20000.: continue
                
                    #print("hit", hit_y, hit_z)
                    color_cbf = red_cbf
                    if check_orientation(int(hit_z)) == 'VBar':
                        color_cbf = blue_cbf
                    elif check_orientation(int(hit_z)) == 'XBar':
                        color_cbf = black_cbf
                    
                    if hit < 3:
                        x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', hit_z), color = color_cbf, label = 'hit area %s' % check_orientation(int(hit_z)))
                    else:
                        x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', hit_z), color = color_cbf)
                    z_y.fill_between(*hit_size(hit_z, hit_y, 'zy', hit_z), color = color_cbf)
                    x_y.fill_between(*hit_size(hit_x, hit_y, 'xy', hit_z), color = color_cbf, alpha = 0.5, linewidth = 0.5)
                
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

                    x_z.plot(kal_z[1:], kal_x[1:], ls='-', lw=2, color=black_cbf, label = 'Kalman reco')
                    z_y.plot(kal_z[1:], kal_y[1:], ls='-', lw=2, color=black_cbf)
                    x_y.plot(kal_x[1:], kal_y[1:], ls='-', lw=2, color=black_cbf)

                    x_z.plot(kal_z[1:], kal_true_x[1:], ls='--', lw=2, color=green_cbf, label = 'Kalman true')
                    z_y.plot(kal_z[1:], kal_true_y[1:], ls='--', lw=2, color=green_cbf)
                    x_y.plot(kal_true_x[1:], kal_true_y[1:], ls='--', lw=2, color=green_cbf)

                ### Track start
                #print(StartPos)
                #temporary fix
                if not (StartPos[j*3 + 2] < 11000.):    #StartPos[i*3 + 1] > -2000.0 or 
                    
                    #print("Start", StartPos[i*3 + 0], StartPos[i*3 + 1], StartPos[i*3 + 2])
                
                    if not StartPos[j*3 + 1] == 0.0:
                        x_z.fill_between(*hit_size(StartPos[j*3 + 2], StartPos[j*3 + 0], 'xz', StartPos[j*3 + 2]), color = green_cbf, label = 'Start/End reco')
                        z_y.fill_between(*hit_size(StartPos[j*3 + 2], StartPos[j*3 + 1], 'zy', StartPos[j*3 + 2]), color = green_cbf)
                        x_y.fill_between(*hit_size(StartPos[j*3 + 0], StartPos[j*3 + 1], 'xy', StartPos[j*3 + 2]), color = green_cbf, alpha = 0.5, linewidth = 0.5)

            
                ### Track end
                #print(EndPos)               
                #temporary fix
                if not (EndPos[j*3 + 2] < 11000.):  #EndPos[i*3 + 1] > -2000.0 or 
                    #print(check_orientation(int(EndPos[i*3 + 2])))
                    #print("End", EndPos[i*3 + 0], EndPos[i*3 + 1], EndPos[i*3 + 2])
    
                    if not EndPos[j*3 + 1] == 0.0:
                        x_z.fill_between(*hit_size(EndPos[j*3 + 2], EndPos[j*3 + 0], 'xz', EndPos[j*3 + 2]), color = green_cbf)
                        z_y.fill_between(*hit_size(EndPos[j*3 + 2], EndPos[j*3 + 1], 'zy', EndPos[j*3 + 2]), color = green_cbf)
                        x_y.fill_between(*hit_size(EndPos[j*3 + 0], EndPos[j*3 + 1], 'xy', EndPos[j*3 + 2]), color = green_cbf, alpha = 0.5, linewidth = 0.5)
                
                ### Track direction
                #print(Direction)              
                #temporary fix
                # Add check on DrawKalmanTrack so we draw the true kalman info instead of a line
                if not DrawKalmanTrack and not (StartPos[j*3 + 2] < 11000. or EndPos[j*3 + 2] < 11000.):   #StartPos[i*3 + 1] > -2000.0 or EndPos[i*3 + 1] > -2000.0 or 

                    #print("Direction", StartPos[i*3 + 1], StartPos[i*3 + 2])

                    if not StartPos[j*3 + 1] == 0.0 or EndPos[j*3 + 1] == 0.0:
                        x_z.plot([StartPos[j*3 + 2] / 1000.0, EndPos[j*3 + 2] / 1000.0], [StartPos[j*3 + 0] / 1000.0, EndPos[j*3 + 0] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--', label = 'Direction')
                        z_y.plot([StartPos[j*3 + 2] / 1000.0, EndPos[j*3 + 2] / 1000.0], [StartPos[j*3 + 1] / 1000.0, EndPos[j*3 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
                        x_y.plot([StartPos[j*3 + 0] / 1000.0, EndPos[j*3 + 0] / 1000.0], [StartPos[j*3 + 1] / 1000.0, EndPos[j*3 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
        
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
                mp.savefig(output_filename + ".png", bbox_inches = 'tight')
                mp.close()
        
    return

### This is for cplotting the hits according to their different orientations
def check_orientation(hit_z):
    return layer_dict["%s" % hit_z]

### Dictionary that after calculate_layers contains for each z-coordinate the orientation str
first_z = 11185
layer_dict = { "%s" % first_z : "UBar" }
        
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
    parser.add_argument('--readout_filename', "-r", type = str, help = "(optional) A file with the raw readout.", default = "")
    parser.add_argument('--report_true_ke', help = "Add the true KE of muon to plot.", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Xlayers', "-X", help = "Does the geometry use X (90 degree orientated) scintillator layers? Yes -> --Xlayers, No -> --no-Xlayers", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    readout_filename =  args.readout_filename
    report_true_ke = args.report_true_ke
    Xlayers = args.Xlayers
    print(Xlayers)
    calculate_layers(Xlayers)
    print(layer_dict)
    
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, report_true_ke)

