import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
#import cppyy.ll
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

# dyslexia friendly background
mp.rcParams['axes.facecolor'] = '#fafafa'
mp.rcParams["figure.facecolor"] = '#fafafa'
mp.rcParams["savefig.facecolor"] = '#fafafa'

# off-black text
mp.rcParams['text.color'] = '#424242'
mp.rcParams['axes.labelcolor'] = '#424242'
mp.rcParams['xtick.color'] = '#424242'
mp.rcParams['ytick.color'] = '#424242'

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

tms_start = 11.124
tms_end = 18.544
tms_start_thick = 14.460
tms_start_double = 17.520
tms_inner_steel = 1.87
tms_outer_steel = 3.73

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

### Function which compares truth to reco - in order to only plot extremes (toggle True/Fase in draw_spill)
def draw_performance(out_dir, input_filename,i,j):

    # Make sure we read in the correct file and have the output directory (we try to make an output directory if one doesn't exist)
    if not os.path.exists(input_filename): raise ValueError(f"Cannor find input_filename {input_filename}")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir):
        raise ValueError(f"Could not make out_dir {out_dir}")
            
    # Read in the Reco_Tree that contains the TMS_Tracks, checking if we have entries
    r = ROOT.TChain("Reco_Tree")
    r.Add(input_filename)
    if not r.GetEntries() > 0:
        print("Didn't get any entries, are you sure the input_filename is right?\n", input_filename)
    
    # Then read in the Truth_Info with the truth tracks, checking if we have entries
    truth = ROOT.TChain("Truth_Info")
    truth.Add(input_filename)
    if not truth.GetEntries() > 0:
        print("Didn't get any entries in Truth_Info, are you sure the input_filename is right?\n", input_filename)

    # Additionally, we read in the Line_Candidates, which is closer to the extrapolation than Reco_Tree which was steps after the extrapolation like 3d track matching and Kalman filter
    b=ROOT.TChain("Line_Candidates")
    b.Add(input_filename)
    if not b.GetEntries() > 0:
        print("Didn't get any entries, are you sure the input_filename is right?\n", input_filename)
    
    # make a dictionary for the spill number and number of events
    spill_number_cache = dict()

    #take spill number and event number as an input to look at a specific event
    try:
        event = None
        true_event = None
        branch = None
    except KeyError:
        r.GetEntry(i)
        event = r
        b.GetEntry(i)
        branch=b
        truth.GetEntry(i)
        true_event = truth
    #here, we fill events, true_events and "branches" with our truth, reco, and Line_Candidates
    if event == None:
        r.GetEntry(i)
        event = r
    if branch == None:
        b.GetEntry(i)
        branch=b
    if true_event == None:
        truth.GetEntry(i)
        true_event = truth

    spill_number = event.SpillNo
    spill_number_cache[i] = spill_number
    #end position for true_event (truth) using StartPos and EndPos which return start/end coordinates for each event
    EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
    True_Position_TMS_End = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionLeavingTMS, dtype = np.float32)

    TrueFiducialEnd = true_event.RecoTrackPrimaryParticleTMSFiducialEnd

    #we need nHits to loop over when we look at reco events
    nTracks=event.nTracks
    nHits = np.frombuffer(event.nHits, dtype = np.uint8)
    nHits = np.array([nHits[j] for j in range(0, nTracks * 4, 4)])

    #from our line candidates, we take a look at the U branch - the no. lines and hits in each line, and we get the TrackHits
    nLinesU = branch.nLinesU
    nHitsU = np.frombuffer(branch.nHitsInTrackU, dtype = np.uint8)
    nHitsU = np.array([nHitsU[i] for i in range(0, nLinesU * 4, 4)])
    #Track hits are in the format [z,x,z,x...0,0,...] so there are nHits lots of [x,z] and then a bunch of zeros
    TrackHitPosU = np.frombuffer(branch.TrackHitPosU, dtype = np.float32)
    #V:
    nLinesV = branch.nLinesV
    nHitsV = np.frombuffer(branch.nHitsInTrackV, dtype = np.uint8)
    nHitsV = np.array([nHitsV[i] for i in range(0, nLinesV * 4, 4)])
    TrackHitPosV = np.frombuffer(branch.TrackHitPosV, dtype = np.float32)
    #X:
    nLinesX = branch.nLinesX
    nHitsX = np.frombuffer(branch.nHitsInTrackX, dtype = np.uint8)
    nHitsX = np.array([nHitsX[i] for i in range(0, nLinesX * 4, 4)])
    TrackHitPosX = np.frombuffer(branch.TrackHitPosX, dtype = np.float32)

    #we take j as an input to look at a specific track
    #we want our reco track to have a start position and for the true start position to have the a non-invalid (where it would be saved as -9999) value
    #the "+0" could be "+1" or "+2", too, it doesn't matter: if one is -9999, they all will be
    #as above, we want our reco track to have an end position and valid end point (contained)
    UX_separation=0
    UZ_separation=0
    VX_separation=0
    VZ_separation=0
    XY_separation=0
    XZ_separation=0

    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
        #when it comes to the LineCandidates, let's still only look at contained tracks with a valid end point, but not worry about the length of track as above
        if (TrueFiducialEnd[j]):
            Primary_True_End_x = True_Position_TMS_End[j*4 + 0]
            Primary_True_End_y = True_Position_TMS_End[j*4 + 1]
            Primary_True_End_z = True_Position_TMS_End[j*4 + 2]
            #make sure there are actual U lines
            if nLinesU != 0:
                Uhits_x = np.zeros(nHitsU[j])
                Uhits_z = np.zeros(nHitsU[j])
                #we prepare variables which will mark start and end
                max_z=-9999
                max_z_index=0
                #these will find the end point (most downstream)
                #we loop over the number of hits, k for a particular track (labelled j)
                for k in range(nHitsU[j]):
                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                    # finally, /1000 to put in mm. 
                    Uhits_x[k] = TrackHitPosU[j*400 + k*2 + 1] / 1000.0 
                    Uhits_z[k] = TrackHitPosU[j*400 + k*2 + 0] / 1000.0
                    #we iterate over our data to save the most downstream and most upstream hits per event in U_end
                    if Uhits_z[k]>max_z:
                        max_z=Uhits_z[k]
                        max_z_index=k
                #U_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                U_End_x=Uhits_x[max_z_index]
                U_End_z=Uhits_z[max_z_index]
                if U_End_x!=-9999.:
                    UX_separation=U_End_x-Primary_True_End_x
                    UZ_separation=U_End_z-Primary_True_End_z
            
            #make sure there are actual V lines
            if nLinesV != 0:
                Vhits_x = np.zeros(nHitsV[j])
                Vhits_z = np.zeros(nHitsV[j])
                #we prepare variables which will mark start and end
                max_z=-9999
                max_z_index=0
                #we loop over the number of hits, k for a particular track (labelled j)
                for k in range(nHitsV[j]):
                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                    # finally, /1000 to put in mm. 
                    Vhits_x[k] = TrackHitPosV[j*400 + k*2 + 1] / 1000.0 
                    Vhits_z[k] = TrackHitPosV[j*400 + k*2 + 0] / 1000.0
                    #we iterate over our data to save the most downstream hit per event in V_end
                    if Vhits_z[k]>max_z:
                        max_z=Vhits_z[k]
                        max_z_index=k
                #V_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                V_End_x=Vhits_x[max_z_index]
                V_End_z=Vhits_z[max_z_index]
                if V_End_x!=-9999.:
                    VX_separation=V_End_x-Primary_True_End_x
                    VZ_separation=V_End_z-Primary_True_End_z
            #make sure there are actual X lines
            if nLinesX != 0:
                Xhits_y = np.zeros(nHitsX[j])
                Xhits_z = np.zeros(nHitsX[j])
                #we prepare variables which will mark start and end
                max_z=-9999
                max_z_index=0
                #we loop over the number of hits, k for a particular track (labelled j)
                for k in range(nHitsX[j]):
                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                    # finally, /1000 to put in mm. 
                    Xhits_y[k] = TrackHitPosX[j*400 + k*2 + 1] / 1000.0 
                    Xhits_z[k] = TrackHitPosX[j*400 + k*2 + 0] / 1000.0
                    #we iterate over our data to save the most downstream hit per event in V_end
                    if Xhits_z[k]>max_z:
                        max_z=Xhits_z[k]
                        max_z_index=k
                    #V_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                X_End_y=Xhits_y[max_z_index]
                X_End_z=Xhits_z[max_z_index]
                if X_End_y!=-9999.:
                    XY_separation=X_End_y-Primary_True_End_y
                    XZ_separation=X_End_z-Primary_True_End_z
    
    output = False
    if  UZ_separation<-12000:
        output = True
    return output

### Actual function that loops through the spills
def draw_spill(out_dir, name, input_filename, spill_number, time_slice, report_true_ke = False):
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
            
    max_n_spills = 10000 # TODO (old) add some meta info to output file with n spill info for file

    start_spill = 0
    if spill_number != -1:
        start_spill = spill_number
        max_n_spills = spill_number + 1
    
    spill_number_cache = dict()
    n_events = r.GetEntries()

    # First loop through all the slices and draw one overall spill
    for current_spill_number in range(start_spill, max_n_spills):
        counter_tracks = 0
        for i in range(n_events):
            try:
                spill_number = spill_number_cache[i]
                event = None
                true_event = None
                branch = None
            except KeyError:
                r.GetEntry(i)
                event = r
                truth.GetEntry(i)
                true_event = truth
                b.GetEntry(i)
                branch = b
                
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            
            if event == None:
                r.GetEntry(i)
                event = r
            if branch == None:
                b.GetEntry(i)
                branch = b
            if true_event == None:
                truth.GetEntry(i)
                true_event = truth
            
            ### Check if a track exists in the event/spill, otherwise skip it
            nTracks = event.nTracks

            if nTracks <= 0: continue
            else: counter_tracks += nTracks
                
            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])

            
            # 2d lines
            nLinesU = branch.nLinesU
            nLinesV = branch.nLinesV
            nLinesX = branch.nLinesX
            nLinesY = branch.nLinesY
            print("Lines: ", nLinesU, nLinesV, nLinesX, nLinesY)
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
            # 2d clusters
            nClustersU = branch.nClustersU
            nClustersV = branch.nClustersV
            nClustersX = branch.nClustersX
            nClustersY = branch.nClustersY
            print("Clusters: ", nClustersU, nClustersV, nClustersX, nClustersY)
            nHitsInClusterU = np.frombuffer(branch.nHitsInClusterU, dtype = np.uint8)
            nHitsInClusterU = np.array([nHitsInClusterU[i] for i in range(0, nClustersU * 4, 4)])
            nHitsInClusterV = np.frombuffer(branch.nHitsInClusterV, dtype = np.uint8)
            nHitsInClusterV = np.array([nHitsInClusterV[i] for i in range(0, nClustersV * 4, 4)])
            nHitsInClusterX = np.frombuffer(branch.nHitsInClusterX, dtype = np.uint8)
            nHitsInClusterX = np.array([nHitsInClusterX[i] for i in range(0, nClustersX * 4, 4)])
            nHitsInClusterY = np.frombuffer(branch.nHitsInClusterY, dtype = np.uint8)
            nHitsInClusterY = np.array([nHitsInClusterY[i] for i in range(0, nClustersY * 4, 4)])
            ClusterHitPosU = np.frombuffer(branch.ClusterHitPosU, dtype = np.float32)
            ClusterHitPosV = np.frombuffer(branch.ClusterHitPosV, dtype = np.float32)
            ClusterHitPosX = np.frombuffer(branch.ClusterHitPosX, dtype = np.float32)
            ClusterHitPosY = np.frombuffer(branch.ClusterHitPosY, dtype = np.float32)
            # true hit information
            nRecoTrackN = true_event.RecoTrackN
            nTrueHits = np.frombuffer(true_event.RecoTrackNHits, dtype = np.uint8)
            nTrueHits = np.array([nTrueHits[i] for i in range(0, nRecoTrackN * 4, 4)])
            True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype = np.float32)

            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)            
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)            
            #Direction = np.frombuffer(event.Direction, dtype = np.float32)
            
            RecoTrackPrimaryParticleTruePositionTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackStart, dtype = np.float32)
            RecoTrackPrimaryParticleTruePositionTrackEnd = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionTrackEnd, dtype = np.float32)
  
            for j in range(nTracks):
                plot_true=draw_performance(out_dir,input_filename,i,j)
                if plot_true==True:
                    ### Create subplots
                    fig = mp.figure(constrained_layout = False)
                    gs = fig.add_gridspec(2, 2, hspace = 0.25, wspace = 0.15)
                    x_y = fig.add_subplot(gs[0, 0])
                    x_z = fig.add_subplot(gs[0:, 1:])
                    z_y = fig.add_subplot(gs[1, 0])
    
                    ### Set labels and ticks
                    x_y.set(xlabel = 'x [m]', ylabel = 'y [m]', xticks = [3, 2, 1, 0, -1, -2, -3, -4], yticks = [-3, -2, -1, 0])
                    z_y.set(xlabel = 'z [m]', ylabel = 'y [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0])
                    x_z.set(xlabel = 'z [m]', ylabel = 'x [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0, 1, 2, 3])
                    x_y.text(tms_outer_steel + 0.06, -2, 'front view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                    z_y.text(tms_end + 0.056, -2, 'side view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
                    x_z.text(tms_end + 0.056, -0.5, 'top view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            
                    ### Set TMS name
                    x_y.text(-tms_outer_steel, tms_top_hybrid, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8) #0.1
                    z_y.text(tms_start - 0.044, tms_top_hybrid, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8) #0.1
                    x_z.text(tms_start - 0.044, tms_outer_steel + 0.06, 'TMS', fontsize = 14, fontweight = 'bold', color = orange_cbf, alpha = 0.8)

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

                    ### Put in outlines of scintillator parts
                    x_z.hlines(-tms_outer_steel, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')   # outer steel plate
                    x_z.hlines(tms_outer_steel, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')    # outer steel plate
                    x_z.hlines(-tms_inner_steel, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')   # inner steel plate
                    x_z.hlines(0, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')       # middle of steel
                    x_z.hlines(tms_inner_steel, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')    # inner steel plate
                    x_z.vlines(tms_start, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_z.vlines(tms_end, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_z.vlines(tms_start_thick, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5))) # to thick steel
                    x_z.vlines(tms_start_double, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5))) # to double thick steel

                    z_y.hlines(tms_bottom_hybrid, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')
                    z_y.hlines(tms_top_hybrid, tms_start, tms_end, color = orange_cbf, linewidth = 1, linestyle = ':')
                    z_y.vlines(tms_start, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    z_y.vlines(tms_end, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    z_y.vlines(tms_start_thick, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5)))   # to thick steel
                    z_y.vlines(tms_start_double, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = (0, (1, 5)))   # to double thick steel

                    x_y.hlines(tms_bottom_hybrid, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.hlines(tms_top_hybrid, -tms_outer_steel, tms_outer_steel, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.vlines(-tms_outer_steel, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.vlines(tms_outer_steel, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.vlines(-tms_inner_steel, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.vlines(0, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')
                    x_y.vlines(tms_inner_steel, tms_top_hybrid, tms_bottom_hybrid, color = orange_cbf, linewidth = 1, linestyle = ':')

                        # plot true hits
                    for hit in range(nTrueHits[j]):
                        hit_x = True_Hits[j*800 + hit*4 + 0]
                        hit_y = True_Hits[j*800 + hit*4 + 1]
                        hit_z = True_Hits[j*800 + hit*4 + 2]

                        if hit_z < 11000.: continue
                        
                        if hit + 1 >= nTrueHits[j]:
                            x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', 'XBar'), color = black_cbf, alpha = 0.3, label = 'hit area')
                        else:
                            x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', 'XBar'), color = black_cbf, alpha = 0.3)
                        z_y.fill_between(*hit_size(hit_z, hit_y, 'zy', 'XBar'), color = black_cbf, alpha = 0.3)
                        x_y.fill_between(*hit_size(hit_x, hit_y, 'xy', 'XBar'), color = black_cbf, alpha = 0.2)
                
                    if nLinesU != 0:
                        Uhits_x = np.zeros(nHitsU[j])
                        Uhits_z = np.zeros(nHitsU[j])
                        U_End = np.zeros(2)
                        max_z=-9999
                        max_z_index=0
                        for hit in range(nHitsU[j]):
                            Uhits_z[hit] = TrackHitPosU[j*400 + hit*2 + 0] / 1000.0
                            Uhits_x[hit] = TrackHitPosU[j*400 + hit*2 + 1] / 1000.0
                            if Uhits_z[hit]>max_z:
                                max_z=Uhits_z[hit]
                                max_z_index=hit
                        U_End[0]=Uhits_x[max_z_index]
                        U_End[1]=Uhits_z[max_z_index]
                        x_z.plot(Uhits_z, Uhits_x, ls = '-', lw = 1.3, color = red_cbf, label = '2D U')
                        x_z.scatter(U_End[1],U_End[0],c = red_cbf, marker = '1', alpha = 1,label='U End')

                    if nLinesV != 0:
                        Vhits_x = np.zeros(nHitsV[j])
                        Vhits_z = np.zeros(nHitsV[j])
                        V_End = np.zeros(2)
                        max_z=-9999
                        max_z_index=0
                        for hit in range(nHitsV[j]):
                            Vhits_z[hit] = TrackHitPosV[j*400 + hit*2 + 0] / 1000.0
                            Vhits_x[hit] = TrackHitPosV[j*400 + hit*2 + 1] / 1000.0
                            if Vhits_z[hit]>max_z:
                                max_z=Vhits_z[hit]
                                max_z_index=hit
                        V_End[0]=Vhits_x[max_z_index]
                        V_End[1]=Vhits_z[max_z_index]
                        x_z.plot(Vhits_z, Vhits_x, ls = '-', lw = 1.3, color = green_cbf, label = '2D V')
                        x_z.scatter(V_End[1],V_End[0],c = green_cbf, marker = '1', alpha = 1,label='V End')

                    if nLinesX != 0:
                        helper = 0
                        X_End = np.zeros(2)
                        max_z=-9999
                        max_z_index=0
                        if nLinesX != nLinesU: 
                            if j >= nLinesX:
                                helper = j - (nLinesX - 1)
                        Xhits_y = np.zeros(nHitsX[j - helper])   
                        Xhits_z = np.zeros(nHitsX[j - helper])
                        for hit in range(nHitsX[j - helper]):
                            Xhits_z[hit] = TrackHitPosX[(j - helper)*400 + hit*2 + 0] / 1000.0
                            Xhits_y[hit] = TrackHitPosX[(j - helper)*400 + hit*2 + 1] / 1000.0
                            if Xhits_z[hit]>max_z:
                                max_z=Xhits_z[hit]
                                max_z_index=hit
                        X_End[0]=Xhits_y[max_z_index]
                        X_End[1]=Xhits_z[max_z_index]
                        z_y.plot(Xhits_z, Xhits_y, ls = '-', lw = 1.3, color = blue_cbf, label = '2D X')
                        z_y.scatter(X_End[1],X_End[0],c = blue_cbf, marker = '1', alpha = 1,label='X End')


                    if nLinesY != 0:
                        Yhits_x = np.zeros(nHitsY[j])
                        Yhits_z = np.zeros(nHitsY[j])
                        for hit in range(nHitsY[j]):
                            Yhits_z[hit] = TrackHitPosY[j*400 + hit*2 + 0] / 1000.0
                            Yhits_x[hit] = TrackHitPosY[j*400 + hit*2 + 1] / 1000.0
                        x_z.plot(Yhits_z, Yhits_x, ls = '-', lw = 1.3, color = red_cbf, label = '2D Y')

                    if nClustersU != 0:
                        helper = 0
                        if j >= nClustersU:
                            helper = j - (nClustersU - 1)
                        for hit in range(nHitsInClusterU[j - helper]):
                            Cluster_z = ClusterHitPosU[(j - helper)*400 + hit*2 + 0] / 1000.0
                            Cluster_x = ClusterHitPosU[(j - helper)*400 + hit*2 + 1] / 1000.0
                            #if hit == 0:
                               # x_z.scatter(Cluster_z, Cluster_x, c = red_cbf, marker = '1', alpha = 0.5, label = 'U Cluster Hits')
                           # else:
                                #x_z.scatter(Cluster_z, Cluster_x, c = red_cbf, marker = '1', alpha = 0.5)

                    if nClustersV != 0:
                        helper = 0
                        if j >= nClustersV:
                            helper = j - (nClustersV - 1)
                        for hit in range(nHitsInClusterV[j - helper]):
                            Cluster_z = ClusterHitPosV[(j - helper)*400 + hit*2 + 0] / 1000.0
                            Cluster_x = ClusterHitPosV[(j - helper)*400 + hit*2 + 1] / 1000.0
                            #if hit == 0:
                                #x_z.scatter(Cluster_z, Cluster_x, c = green_cbf, marker = '1', alpha = 0.5, label = 'V Cluster Hits')
                            #else:
                               # x_z.scatter(Cluster_z, Cluster_x, c = green_cbf, marker = '1', alpha = 0.5)

                    if nClustersX != 0:
                        helper = 0
                        if j >= nClustersX:
                            helper = j - (nClustersX - 1)
                        for hit in range(nHitsInClusterX[j - helper]):
                            Cluster_z = ClusterHitPosX[(j - helper)*400 + hit*2 + 0] / 1000.0
                            Cluster_y = ClusterHitPosX[(j - helper)*400 + hit*2 + 1] / 10000
                            #if hit == 0:
                                #z_y.scatter(Cluster_z, Cluster_y, c = blue_cbf, marker = '1', alpha = 0.5, label = 'X Cluster Hits')
                            #else:
                                #z_y.scatter(Cluster_z, Cluster_y, c = blue_cbf, marker = '1', alpha = 0.5)

                    if nClustersY != 0:
                        helper = 0
                        if j >= nClustersY:
                            helper = j - (nClustersY - 1)
                        for hit in range(nHitsInClusterY[j - helper]):
                            Cluster_z = ClusterHitPosY[(j - helper)*400 + hit*2 + 0] / 1000.0
                            Cluster_x = ClusterHitPosY[(j - helper)*400 + hit*2 + 1] / 1000.0
                            #if hit == 0:
                            #    x_z.scatter(Cluster_z, Cluster_x, c = red_cbf, marker = '1', alpha = 0.5,  label = 'Y Cluster Hits')
                           # else:
                                #x_z.scatter(Cluster_z, Cluster_x, c = red_cbf, marker = '1', alpha = 0.5)

                    ### Track start
                    #temporary fix
                    if not (StartPos[j*4 + 2] < 11000.): 
                        
                        if not StartPos[j*4 + 1] == 0.0:
                            x_z.scatter(StartPos[j*4 + 2] / 1000.0, StartPos[j*4 + 0] / 1000.0, c = black_cbf, alpha = 0.5, marker = '2')
                            z_y.scatter(StartPos[j*4 + 2] / 1000.0, StartPos[j*4 + 1] / 1000.0, c = black_cbf, alpha = 0.5, marker = '2', label = 'Start reco')
                            x_y.scatter(StartPos[j*4 + 0] / 1000.0, StartPos[j*4 + 1] / 1000.0, c = black_cbf, alpha = 0.5, marker = '2')

                
                    ### Track end               
                    #temporary fix
                    if not (EndPos[j*4 + 2] < 11000.):  
        
                        if not EndPos[j*4 + 1] == 0.0:
                            x_z.scatter(EndPos[j*4 + 2] / 1000.0, EndPos[j*4 + 0] / 1000.0, c = black_cbf, alpha = 0.5, marker = '1')
                            z_y.scatter(EndPos[j*4 + 2] / 1000.0, EndPos[j*4 + 1] / 1000.0, c = black_cbf, alpha = 0.5, marker = '1', label = 'End reco')
                            x_y.scatter(EndPos[j*4 + 0] / 1000.0, EndPos[j*4 + 1] / 1000.0, c = black_cbf, alpha = 0.5, marker = '1')
                    
                    ### Track direction
                    #temporary fix
                    # Add check on DrawKalmanTrack so we draw the true kalman info instead of a line
                    if not DrawKalmanTrack and not (StartPos[j*4 + 2] < 11000. or EndPos[j*4 + 2] < 11000.): 

                        if not StartPos[j*4 + 1] == 0.0 or EndPos[j*4 + 1] == 0.0:
                            x_z.plot([StartPos[j*4 + 2] / 1000.0, EndPos[j*4 + 2] / 1000.0], [StartPos[j*4 + 0] / 1000.0, EndPos[j*4 + 0] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
                            z_y.plot([StartPos[j*4 + 2] / 1000.0, EndPos[j*4 + 2] / 1000.0], [StartPos[j*4 + 1] / 1000.0, EndPos[j*4 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--', label = 'Direction')
                            x_y.plot([StartPos[j*4 + 0] / 1000.0, EndPos[j*4 + 0] / 1000.0], [StartPos[j*4 + 1] / 1000.0, EndPos[j*4 + 1] / 1000.0], color = black_cbf, linewidth = 1.5, linestyle = '--')
                    
                    if RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] > 11000.: 
                        x_z.scatter(RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 0] / 1000.0, c = magenta_cbf, marker = '2', alpha = 0.5)
                        x_z.scatter(RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 0] / 1000.0, c = magenta_cbf, marker = '1', alpha = 0.5)
                    
                        z_y.scatter(RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackStart[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '2', alpha = 0.5, label = 'Start true')
                        z_y.scatter(RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 2] / 1000.0, RecoTrackPrimaryParticleTruePositionTrackEnd[j*4 + 1] / 1000.0, c = magenta_cbf, marker = '1', alpha = 0.5, label = 'End true')
                
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
    parser.add_argument('--report_true_ke', help = "Add the true KE of muon to plot.", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    report_true_ke = args.report_true_ke
    
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, report_true_ke)
