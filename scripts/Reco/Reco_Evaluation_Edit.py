import ROOT
import numpy as np
import os
import argparse

# function which compares truth to reco - input the file, what we want to do, and the name and directory for output file
def draw_performance(out_dir, input_filename, output_name):

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

    #take a look at this number below!        
    max_n_spills = 129000
    
    # make a dictionary for the spill number and number of events
    spill_number_cache = dict()
    n_events = r.GetEntries()

    # prepare arrays to be filled (if we tell it to do the Start and End)
    #start:
    Reco_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
    Primary_True_Start = np.ones((n_events, 5, 3), dtype = float) * -9999.
    StartX = np.ones((n_events, 5), dtype = float) * -9999.
    #end:
    Reco_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
    Primary_True_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
    EndX = np.ones((n_events, 5), dtype = float) * -9999.
    #energy:
    True_Muon_Energy = np.ones((n_events, 5), dtype = float) * -9999.
    Reco_Energy = np.ones((n_events, 5), dtype = float) * -9999.
    EnergyX=np.ones((n_events, 5), dtype = float) * -9999.
    HitsPerTrack = np.ones((n_events, 5), dtype = float) * -9999.
    #U:
    U_Start=np.ones((n_events,5,2),dtype=float)*-9999.
    U_End=np.ones((n_events,5,2),dtype=float)*-9999.
    #V:
    V_Start=np.ones((n_events,5,2),dtype=float)*-9999.
    V_End=np.ones((n_events,5,2),dtype=float)*-9999.
    #X:
    X_Start=np.ones((n_events,5,2),dtype=float)*-9999.
    X_End=np.ones((n_events,5,2),dtype=float)*-9999.

    # loop through spills, then events, obtaining events to compare against each other
    for current_spill_number in range(max_n_spills):
        for i in range(n_events):
            try:
                spill_number = spill_number_cache[i]
                event = None
                true_event = None
            except KeyError:
                r.GetEntry(i)
                event = r
                b.GetEntry(i)
                branch=b
                truth.GetEntry(i)
                true_event = truth
            #just dw about this stuff with the spill number
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
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
            # here is the number of tracks which we loop through later to look at each track
            nTracks = event.nTracks

            # start and end position for event (reco) and true_event (truth) using StartPos and EndPos which return start/end coordinates for each event
            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)
            True_Position_TMS_Start = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionEnteringTMS, dtype = np.float32)
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            True_Position_TMS_End = np.frombuffer(true_event.RecoTrackPrimaryParticleTruePositionLeavingTMS, dtype = np.float32)

            #define particle identification
            True_PDG=true_event.RecoTrackPrimaryParticlePDG
            #we get initial momentum of muon
            MomentumEnteringTMS = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumEnteringTMS, dtype = np.float32)
            #the true length of track
            length=event.Length
            TrueFiducialEnd = true_event.RecoTrackPrimaryParticleTMSFiducialEnd

            #we need nHits to loop over when we look at reco events
            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[j] for j in range(0, nTracks * 4, 4)])

            # X hits in reconstructed track - a test we perform later on
            BarTypes = np.frombuffer(event.TrackHitBarType, dtype = np.uint8)

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

            #already looping over spills and events, we loop over tracks within the events here
            if nTracks <= 0: continue
            for j in range(nTracks):
                #we don't want to operate on events with too many tracks going on
                if j > 4: break
                #we only concern ourselves with muons
                if (np.abs(True_PDG[j])!=13): break
                # check if X hits in reconstructed track
                ContainsX = False
                for hit in range(nHits[j]):
                    if BarTypes[j*800 + hit*4] == 0:
                        ContainsX = True
                        break
                #we want our reco track to have a start position and for the true start position to have the a non-invalid (where it would be saved as -9999) value
                #the "+0" could be "+1" or "+2", too, it doesn't matter: if one is -9999, they all will be
                if True_Position_TMS_Start[j*4 + 0] > -8000. and not StartPos.size == 0 and True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                    #we want our tracks to have length (end z - start z co-ord for both truth and reco) greater than 14x65, corresponding to more than 14 potential hits (over 65mm scintillator separation)
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                        # we want an end in fiducial volume
                        if (TrueFiducialEnd[j]):
                            EnergyX[i, j] = ContainsX
                            true_muon_starting_ke = MomentumEnteringTMS[j*4 + 3] * 1e-3
                            True_Muon_Energy[i, j] = true_muon_starting_ke
                            starting_ke = (length[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                            Reco_Energy[i, j] = starting_ke
                            #getting the start points for each track and saving them (x,y,z)
                            Reco_Start[i, j, 0] = StartPos[j*3 + 0]
                            Reco_Start[i, j, 1] = StartPos[j*3 + 1]
                            Reco_Start[i, j, 2] = StartPos[j*3 + 2]
                            Primary_True_Start[i, j, 0] = True_Position_TMS_Start[j*4 + 0]
                            Primary_True_Start[i, j, 1] = True_Position_TMS_Start[j*4 + 1]
                            Primary_True_Start[i, j, 2] = True_Position_TMS_Start[j*4 + 2]
                            StartX[i, j] = ContainsX
                            #doing the same for the end position
                            Reco_End[i, j, 0] = EndPos[j*3 + 0]
                            Reco_End[i, j, 1] = EndPos[j*3 + 1]
                            Reco_End[i, j, 2] = EndPos[j*3 + 2]
                            Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                            Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                            Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                            HitsPerTrack[i, j] = event.nHits[j]
                            EndX[i, j] = ContainsX
                
                            #when it comes to the LineCandidates, let's still only look at contained tracks with a valid end point, but not worry about the length of track as above
                            #make sure there are actual U lines
                            if nLinesU != 0:
                                Uhits_x = np.zeros(nHitsU[j])
                                Uhits_z = np.zeros(nHitsU[j])
                                #we prepare variables which will mark start and end
                                max_z=-9999
                                min_z=999999
                                max_z_index=0
                                min_z_index=0
                                #these will find the end point (most downstream)
                                #we loop over the number of hits, k for a particular track (labelled j)
                                for k in range(nHitsU[j]):
                                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                                    Uhits_x[k] = TrackHitPosU[j*400 + k*2 + 1]
                                    Uhits_z[k] = TrackHitPosU[j*400 + k*2 + 0]
                                    #we iterate over our data to save the most downstream and most upstream hits per event in U_end
                                    if Uhits_z[k]>max_z:
                                        max_z=Uhits_z[k]
                                        max_z_index=k
                                    if Uhits_z[k]<min_z:
                                        min_z=Uhits_z[k]
                                        min_z_index=k
                                #U_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                                U_End[i,j,0]=Uhits_x[max_z_index]
                                U_End[i,j,1]=Uhits_z[max_z_index]
                                U_Start[i,j,0]=Uhits_x[min_z_index]
                                U_Start[i,j,1]=Uhits_z[min_z_index]
                        
                            #make sure there are actual V lines
                            if nLinesV != 0:
                                Vhits_x = np.zeros(nHitsV[j])
                                Vhits_z = np.zeros(nHitsV[j])
                                #we prepare variables which will mark start and end
                                max_z=-9999
                                min_z=999999
                                max_z_index=0
                                min_z_index=0
                                #we loop over the number of hits, k for a particular track (labelled j)
                                for k in range(nHitsV[j]):
                                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                                    # finally, /1000 to put in mm. 
                                    Vhits_x[k] = TrackHitPosV[j*400 + k*2 + 1]
                                    Vhits_z[k] = TrackHitPosV[j*400 + k*2 + 0]
                                    #we iterate over our data to save the most downstream hit per event in V_end
                                    if Vhits_z[k]>max_z:
                                        max_z=Vhits_z[k]
                                        max_z_index=k
                                    if Vhits_z[k]<min_z:
                                        min_z=Vhits_z[k]
                                        min_z_index=k
                                #V_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                                V_End[i,j,0]=Vhits_x[max_z_index]
                                V_End[i,j,1]=Vhits_z[max_z_index]
                                V_Start[i,j,0]=Vhits_x[min_z_index]
                                V_Start[i,j,1]=Vhits_z[min_z_index]

                            #make sure there are actual X lines
                            if nLinesX != 0:
                                Xhits_y = np.zeros(nHitsX[j])
                                Xhits_z = np.zeros(nHitsX[j])
                                #we prepare variables which will mark start and end
                                max_z=-9999
                                min_z=999999
                                max_z_index=0
                                min_z_index=0
                                #we loop over the number of hits, k for a particular track (labelled j)
                                for k in range(nHitsX[j]):
                                    #j*400 because for each event we store a max of 200 events (400 data points), so j*400 moves to look at a new event when we increase j
                                    # k*2 iterates over the hits, each hit having an x and z data point, where the +0 or +1 selects z and x respectively
                                    # finally, /1000 to put in mm. 
                                    Xhits_y[k] = TrackHitPosX[j*400 + k*2 + 1]
                                    Xhits_z[k] = TrackHitPosX[j*400 + k*2 + 0]
                                    #we iterate over our data to save the most downstream hit per event in V_end
                                    if Xhits_z[k]>max_z:
                                        max_z=Xhits_z[k]
                                        max_z_index=k
                                    if Xhits_z[k]<min_z:
                                        min_z=Xhits_z[k]
                                        min_z_index=k
                                #X_end will save an x and z data point for the most downstream hit in j: each track in i: each event   
                                X_End[i,j,0]=Xhits_y[max_z_index]
                                X_End[i,j,1]=Xhits_z[max_z_index]
                                X_Start[i,j,0]=Xhits_y[min_z_index]
                                X_Start[i,j,1]=Xhits_z[min_z_index]
                   
                    
    #start point for reco
    boolean_Reco_Start = (Reco_Start[:, :, 0] != -9999.)
    Reco_Start = Reco_Start[boolean_Reco_Start]
    Primary_True_Start = Primary_True_Start[boolean_Reco_Start]
    StartX = StartX[boolean_Reco_Start]
    Reco_Start_x = Reco_Start[:, 0]
    Reco_Start_y = Reco_Start[:, 1]
    Reco_Start_z = Reco_Start[:, 2]
    Primary_True_Start_x = Primary_True_Start[:, 0]
    Primary_True_Start_y = Primary_True_Start[:, 1]
    Primary_True_Start_z = Primary_True_Start[:, 2]

    #end point for reco
    boolean_Reco_End = (Reco_End[:, :, 0] != -9999.)
    Reco_End = Reco_End[boolean_Reco_End]
    Primary_True_End = Primary_True_End[boolean_Reco_End]
    EndX = EndX[boolean_Reco_End]
    Reco_End_x = Reco_End[:, 0]
    Reco_End_y = Reco_End[:, 1]
    Reco_End_z = Reco_End[:, 2]
    Primary_True_End_x = Primary_True_End[:, 0]
    Primary_True_End_y = Primary_True_End[:, 1]
    Primary_True_End_z = Primary_True_End[:, 2]
    
    #start point for U
    boolean_U_Start = (U_Start[:, :, 0] != -9999.)
    U_Start = U_Start[boolean_U_Start]
    U_Start_x=U_Start[:,0]
    U_Start_z=U_Start[:,1]
    #end point for U
    boolean_U_End = (U_End[:, :, 0] != -9999.)
    U_End = U_End[boolean_U_End]
    U_End_x=U_End[:,0]
    U_End_z=U_End[:,1]

    #start point for V
    boolean_V_Start = (V_Start[:, :, 0] != -9999.)
    V_Start = V_Start[boolean_V_Start]
    V_Start_x=V_Start[:,0]
    V_Start_z=V_Start[:,1]
    #end point for V
    boolean_V_End = (V_End[:, :, 0] != -9999.)
    V_End = V_End[boolean_V_End]
    V_End_x=V_End[:,0]
    V_End_z=V_End[:,1]
    
    #start point for X
    boolean_X_Start = (X_Start[:, :, 0] != -9999.)
    X_Start = X_Start[boolean_X_Start]
    X_Start_y=X_Start[:,0]
    X_Start_z=X_Start[:,1]
    #end point for X
    boolean_X_End = (X_End[:, :, 0] != -9999.)
    X_End = X_End[boolean_X_End]
    X_End_y=X_End[:,0]
    X_End_z=X_End[:,1]

    # energy
    boolean_energy = (True_Muon_Energy != -9999.)
    True_Muon_Energy = True_Muon_Energy[boolean_energy]
    Reco_Energy = Reco_Energy[boolean_energy]
    EnergyX = EnergyX[boolean_energy]
    #double check hits per track:
    HitsPerTrack = HitsPerTrack[boolean_energy]

    # we name our output file and give it a directory 
    output = open('%s/%s.txt' % (out_dir, output_name), 'w')
    # we tell our output file what to include
    #reco info - start
    output.write(str([Reco_Start_x[i] for i in range(len(Reco_Start_x))])) #0
    output.write('\n')
    output.write(str([Reco_Start_y[i] for i in range(len(Reco_Start_y))]))
    output.write('\n')
    output.write(str([Reco_Start_z[i] for i in range(len(Reco_Start_z))]))
    output.write('\n')
    #truth info - start
    output.write(str([Primary_True_Start_x[i] for i in range(len(Primary_True_Start_x))])) #3
    output.write('\n')
    output.write(str([Primary_True_Start_y[i] for i in range(len(Primary_True_Start_y))]))
    output.write('\n')
    output.write(str([Primary_True_Start_z[i] for i in range(len(Primary_True_Start_z))]))
    output.write('\n')
    output.write(str([StartX[i] for i in range(len(StartX))]))
    output.write('\n')
    #end - reco vs truth
    output.write(str([Reco_End_x[i] for i in range(len(Reco_End_x))])) #7
    output.write('\n')
    output.write(str([Reco_End_y[i] for i in range(len(Reco_End_y))]))
    output.write('\n')
    output.write(str([Reco_End_z[i] for i in range(len(Reco_End_z))]))
    output.write('\n')
    output.write(str([Primary_True_End_x[i] for i in range(len(Primary_True_End_x))])) #10
    output.write('\n')
    output.write(str([Primary_True_End_y[i] for i in range(len(Primary_True_End_y))]))
    output.write('\n')
    output.write(str([Primary_True_End_z[i] for i in range(len(Primary_True_End_z))]))
    output.write('\n')
    output.write(str([EndX[i] for i in range(len(EndX))]))
    output.write('\n')
    #line candidates info - U
    output.write(str([U_Start_x[i] for i in range(len(U_Start_x))])) #14
    output.write('\n')
    output.write(str([U_Start_z[i] for i in range(len(U_Start_z))]))
    output.write('\n')
    output.write(str([U_End_x[i] for i in range(len(U_End_x))])) #16
    output.write('\n')
    output.write(str([U_End_z[i] for i in range(len(U_End_z))]))
    output.write('\n')
    #line candidates info - V
    output.write(str([V_Start_x[i] for i in range(len(V_Start_x))])) #18
    output.write('\n')
    output.write(str([V_Start_z[i] for i in range(len(V_Start_z))]))
    output.write('\n')
    output.write(str([V_End_x[i] for i in range(len(V_End_x))])) #20
    output.write('\n')
    output.write(str([V_End_z[i] for i in range(len(V_End_z))]))
    output.write('\n')
    #line candidates info - X
    output.write(str([X_Start_y[i] for i in range(len(X_Start_y))])) #22
    output.write('\n')
    output.write(str([X_Start_z[i] for i in range(len(X_Start_z))]))
    output.write('\n')
    output.write(str([X_End_y[i] for i in range(len(X_End_y))])) #24
    output.write('\n')
    output.write(str([X_End_z[i] for i in range(len(X_End_z))]))
    output.write('\n')
    # energy - reco vs. truth
    output.write(str([HitsPerTrack[i] for i in range(len(HitsPerTrack))])) #26
    output.write('\n')
    output.write(str([EndX[i] for i in range(len(EndX))]))
    output.write('\n')
    output.write(str([True_Muon_Energy[i] for i in range(len(True_Muon_Energy))]))
    output.write('\n')
    output.write(str([Reco_Energy[i] for i in range(len(Reco_Energy))]))
    output.write('\n')
    output.write(str([EnergyX[i] for i in range(len(EnergyX))]))
    output.close()

    #we close our function, returning the ourput file
    return

#below, we actually run the function above
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    #give it a place to save
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    #give it a data file
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    #give it a filename 
    parser.add_argument('--output', type = str, help = "The output file name.")
    #parse the arguments above as variables for the function
    args = parser.parse_args() 
    out_dir = args.outdir
    input_filename = args.input_filename
    output = args.output
    #run the function
    draw_performance(out_dir, input_filename, output)





