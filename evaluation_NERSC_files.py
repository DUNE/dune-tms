import ROOT
import numpy as np
import os
import argparse
import cppyy.ll

import math # filter out nans

import csv
import matplotlib.pyplot as mp

### Actual function that loops through the spills
def draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained, output_name):
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
        Reco_End_EscX = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Reco_End_EscY = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Primary_True_End_EscX = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Primary_True_End_EscY = np.ones((n_events, 5, 3), dtype = float) * -9999.
    if plot_Charge:
        Reco_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_Charge_stop = np.ones((n_events, 5), dtype = float) * -9999.
    if plot_Charge or plot_Angle:
        True_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_KE = np.ones((n_events, 5), dtype = float) * -9999.
        True_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.
        Reco_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.
        True_ZDirection = np.ones((n_events, 5), dtype = float) * -9999.
        Reco_ZDirection = np.ones((n_events, 5), dtype = float) * -9999.
    True_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true
    Reco_Muon_Track = np.zeros((n_events, 5), dtype = float)    # treat as boolean array: 0 -> false, 1 -> true
    HitsPerTrack = np.ones((n_events, 5), dtype = float) * -9999.

    # some simple counters for efficiency evaluation
    correct_tracks_reco = 0
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

            nTracks = event.nTracks
            
            # track direction/momentum at start
            if plot_Angle or plot_Charge:
                Reco_Track_StartDirection = np.frombuffer(event.StartDirection, dtype = np.float32)
                MomentumTrackStart = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumTrackStart, dtype = np.float32)
            # charge related stuff
            if plot_Charge:
                Reco_Track_Charge = event.Charge
            True_PDG = true_event.RecoTrackPrimaryParticlePDG
            # classify if true muon actually enough distance in TMS and start in LAr
            LArFiducialTouch = true_event.RecoTrackPrimaryParticleLArFiducialStart
            Particle_PDG = true_event.LeptonPDG
            Muon_Start = np.frombuffer(true_event.Muon_Vertex, dtype = np.float32)
            Muon_End = np.frombuffer(true_event.Muon_Death, dtype = np.float32)
            True_Momentum_Leaving = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumLeavingTMS, dtype = np.float32)
            # actual hits
            Reco_Hits = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            True_Hits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype = np.float32)
            sum_reco_hits = event.nHits
            sum_true_hits = true_event.RecoTrackNHits            
            
            if (abs(Particle_PDG) == 13):
                if 4179.24 < Muon_Start[2] < 9135.88 and 11185 < Muon_End[2] < 18535:
                    if 371.77 > Muon_End[1] > -3076.23 and abs(Muon_End[0]) < 3491:
                        count_muons += 1

            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])
            nTrueHits = true_event.RecoTrackNHits
            #nTrueHits = np.array([nTrueHits[i] for i in range(0, nTracks * 4, 4)])

            if nTracks <= 0: continue
            #if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break

                # check if (anti-)muon as true primary particle and if origin in LAr
                #if True: #np.abs(True_PDG[j]) == 13 and LArFiducialTouch[j]:
                #    # check if (anti-)muon travesers at least 4 planes in TMS
                #    if (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) >= 440.:
                #        # if so, then this is a true muon: 0 -> 1
                #        True_Muon_Track[i, j] = 1.
                #        counter_correct = 0
                #        for true_hits in range(int(sum_true_hits[j])):
                #            for reco_hits in range(int(sum_reco_hits[j])):
                #                if True_Hits[j*600 + true_hits*4 + 2] == Reco_Hits[j*600 + reco_hits*3 + 2]:
                #                    if np.abs(True_Hits[j*600 + true_hits*4 + 0] - Reco_Hits[j*600 + reco_hits*3 + 0]) <= 2 * 36:
                #                        counter_correct += 1

                #        if counter_correct / sum_true_hits[j] >= 0.1:
                #            correct_tracks_reco += 1
                        
                #        # check if reconstructed tracks exist for this event
                #        if StartPos.size != 0:
                #            # if so, then add all identified as muons: 0 -> number tracks
                #            Reco_Muon_Track[i, j] = 1. #StartPos.size / 3
                
                if True_Position_TMS_Start[j*4 + 0] > -8000. and not StartPos.size == 0:
                    # checking for muon tracks (minimal length for this are 20 planes traversed -> 890 mm in thin area
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # normalization of true direction
                            if plot_Angle or plot_Charge:
                                magnitude = np.sqrt(MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 + 1]**2 + MomentumTrackStart[j*4 + 2]**2)
                                magnitude_r = np.sqrt(Reco_Track_StartDirection[j*3 + 0]**2 + Reco_Track_StartDirection[j*3 + 1]**2 + Reco_Track_StartDirection[j*3 + 2]**2)
                                True_TrackDirection[i, j, 0] = dir_to_angle(MomentumTrackStart[j*4 + 0] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                True_TrackDirection[i, j, 1] = dir_to_angle(MomentumTrackStart[j*4 + 1] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                Reco_TrackDirection[i, j, 0] = dir_to_angle(Reco_Track_StartDirection[j*3 + 0], Reco_Track_StartDirection[j*3 + 2])
                                Reco_TrackDirection[i, j, 1] = dir_to_angle(Reco_Track_StartDirection[j*3 + 1], Reco_Track_StartDirection[j*3 + 2])
                                True_ZDirection[i, j] = MomentumTrackStart[j*4 + 2] / magnitude
                                Reco_ZDirection[i, j] = Reco_Track_StartDirection[j*3 + 2] / magnitude_r
                                True_KE[i, j] = np.sqrt((MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 +1]**2 + MomentumTrackStart[j*4 + 2]**2 + 105.7**2) - 105.7)
                                True_Charge[i, j] = True_PDG[j]  #13                          
                            if plot_Charge:
                                # Make sure that true muon stopped in TMS and is not leaving
                                if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                    Reco_Charge[i, j] = Reco_Track_Charge[j]
                                    True_Charge_stop[i, j] = True_PDG[j]    #13
                                    if Reco_Track_Charge[j] != True_PDG[j]:
                                        TrackHitPos = np.frombuffer(event.TrackHitPos, dtype = np.float32)
                                        TrueHits = np.frombuffer(true_event.RecoTrackTrueHitPosition, dtype = np.float32)
                                        mp.hlines(-3.49, 11.176, 18.544, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.hlines(3.49, 11.176, 18.544, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.hlines(-1.75, 11.176, 18.544, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.hlines(0, 11.176, 18.544, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.hlines(1.75, 11.176, 18.544, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.vlines(11.176, -3.49, 3.52, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.vlines(18.544, -3.49, 3.52, color = 'blue', linewidth = 1, linestyle = ':')
                                        mp.axis('equal')
                                        for hit in range(nHits[j]):
                                            hit_x = TrackHitPos[j*600 + hit*3 + 0]
                                            hit_z = TrackHitPos[j*600 + hit*3 + 2]
                                            mp.scatter(hit_z / 1000, hit_x / 1000, marker = 'd', c = 'black', alpha = 0.5)
                                        for hit in range(nTrueHits[j]):
                                            thit_x = TrueHits[j*600 + hit*4 + 0]
                                            thit_z = TrueHits[j*600 + hit*4 + 2]
                                            if thit_z < 11000.: continue
                                            mp.scatter(thit_z / 1000, thit_x / 1000, marker = '+', c = 'yellow')
                                        mp.xlabel('z [m]')
                                        mp.ylabel('x [m]')
                                        mp.scatter(True_Position_TMS_End[j*4 + 2] / 1000, True_Position_TMS_End[j*4 + 0] / 1000, c = 'red', marker = '2')
                                        mp.scatter(True_Position_TMS_Start[j*4 + 2] / 1000, True_Position_TMS_Start[j*4 + 0] / 1000, c = 'red', marker = '1')
                                        mp.scatter(EndPos[j*3 + 2] / 1000, EndPos[j*3 + 0] / 1000, c = 'green', marker = '3')
                                        mp.scatter(StartPos[j*3 + 2] / 1000, StartPos[j*3 + 0] / 1000, c = 'green', marker = '4')
                                        mp.savefig('events/event_%03d_%03d_%02d.png' % (current_spill_number, i, j), bbox_inches = 'tight')
                                        mp.close()

                            if plot_Start:
                                Reco_Start[i, j, 0] = StartPos[j*3 + 0]
                                Reco_Start[i, j, 1] = StartPos[j*3 + 1]
                                Reco_Start[i, j, 2] = StartPos[j*3 + 2]
                                Primary_True_Start[i, j, 0] = True_Position_TMS_Start[j*4 + 0]
                                Primary_True_Start[i, j, 1] = True_Position_TMS_Start[j*4 + 1]
                                Primary_True_Start[i, j, 2] = True_Position_TMS_Start[j*4 + 2]
                if plot_End and Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
                            if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                Reco_End[i, j, 0] = EndPos[j*3 + 0]
                                Reco_End[i, j, 1] = EndPos[j*3 + 1]
                                Reco_End[i, j, 2] = EndPos[j*3 + 2]
                                Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                                Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                                Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                                HitsPerTrack[i, j] = event.nHits[j]
                            else:
                                counter_leaving += 1
                if plot_End and not Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
                            if np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) >= 1:
                                counter_leaving += 1
                            Reco_End[i, j, 0] = EndPos[j*3 + 0]
                            Reco_End[i, j, 1] = EndPos[j*3 + 1]
                            Reco_End[i, j, 2] = EndPos[j*3 + 2]
                            if np.abs(EndPos[j*3 + 0]) > 3419.:
                                if EndPos[j*3 + 0] > 0:
                                    Reco_End[i, j, 0] = 3419.
                                else:
                                    Reco_End[i, j, 0] = -3419.
                            if EndPos[j*3 + 1] > 291.77:
                                Reco_End[i, j, 1] = 291.77
                            if EndPos[j*3 + 1] < -2976.23:
                                Reco_End[i, j, 1] = -2976.23
                            Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                            Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                            Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                            if np.abs(True_Position_TMS_End[j*4 + 0]) >= 3219:
                                Reco_End_EscX[i, j, 0] = EndPos[j*3 + 0]
                                Reco_End_EscX[i, j, 1] = EndPos[j*3 + 1]
                                Reco_End_EscX[i, j, 2] = EndPos[j*3 + 2]
                                Primary_True_End_EscX[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                                Primary_True_End_EscX[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                                Primary_True_End_EscX[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                            elif True_Position_TMS_End[j*4 + 1] >= 271.77 or True_Position_TMS_End[j*4 + 1] <= -2776.23:
                                Reco_End_EscY[i, j, 0] = EndPos[j*3 + 0]
                                Reco_End_EscY[i, j, 1] = EndPos[j*3 + 1]
                                Reco_End_EscY[i, j, 2] = EndPos[j*3 + 2]
                                Primary_True_End_EscY[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                                Primary_True_End_EscY[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                                Primary_True_End_EscY[i, j, 2] = True_Position_TMS_End[j*4 + 2]
    
    boolean_True_KE = (True_KE != -9999.)
    True_KE = True_KE[boolean_True_KE]
    boolean_Reco_Start = (Reco_Start[:, :, 0] != -9999.)
    Reco_Start = Reco_Start[boolean_Reco_Start]
    Primary_True_Start = Primary_True_Start[boolean_Reco_Start]
    Reco_End = Reco_End[boolean_Reco_Start]
    Primary_True_End = Primary_True_End[boolean_Reco_Start]
    Reco_Charge = Reco_Charge[boolean_True_KE]
    True_Charge = True_Charge[boolean_True_KE]
    True_Charge_stop = True_Charge_stop[boolean_True_KE]
    True_TrackDirection = True_TrackDirection[boolean_Reco_Start]
    Reco_TrackDirection = Reco_TrackDirection[boolean_Reco_Start]
    True_ZDirection = True_ZDirection[boolean_True_KE]
    Reco_ZDirection = Reco_ZDirection[boolean_True_KE]
    HitsPerTrack = HitsPerTrack[boolean_True_KE]
    
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
        HitsPerTrack
    if plot_Angle or plot_Charge:
        True_TrackDirection_xz = True_TrackDirection[:, 0]
        True_TrackDirection_yz = True_TrackDirection[:, 1]
        Reco_TrackDirection_xz = Reco_TrackDirection[:, 0]
        Reco_TrackDirection_yz = Reco_TrackDirection[:, 1]
    
    output = open('%s/%s.txt' % (out_dir, output_name), 'w')
    output.write(str([Reco_Start_x[i] for i in range(len(Reco_Start_x))]))
    output.write('\n')
    output.write(str([Reco_Start_y[i] for i in range(len(Reco_Start_y))]))
    output.write('\n')
    output.write(str([Reco_Start_z[i] for i in range(len(Reco_Start_z))]))
    output.write('\n')
    output.write(str([Primary_True_Start_x[i] for i in range(len(Primary_True_Start_x))]))
    output.write('\n')
    output.write(str([Primary_True_Start_y[i] for i in range(len(Primary_True_Start_y))]))
    output.write('\n')
    output.write(str([Primary_True_Start_z[i] for i in range(len(Primary_True_Start_z))]))
    output.write('\n')
    output.write(str([Reco_End_x[i] for i in range(len(Reco_End_x))]))
    output.write('\n')
    output.write(str([Reco_End_y[i] for i in range(len(Reco_End_y))]))
    output.write('\n')
    output.write(str([Reco_End_z[i] for i in range(len(Reco_End_z))]))
    output.write('\n')
    output.write(str([Primary_True_End_x[i] for i in range(len(Primary_True_End_x))]))
    output.write('\n')
    output.write(str([Primary_True_End_y[i] for i in range(len(Primary_True_End_y))]))
    output.write('\n')
    output.write(str([Primary_True_End_z[i] for i in range(len(Primary_True_End_z))]))
    output.write('\n')
    output.write(str([Reco_TrackDirection_xz[i] for i in range(len(Reco_TrackDirection_xz))]))
    output.write('\n')
    output.write(str([Reco_TrackDirection_yz[i] for i in range(len(Reco_TrackDirection_yz))]))
    output.write('\n')
    output.write(str([True_TrackDirection_xz[i] for i in range(len(True_TrackDirection_xz))]))
    output.write('\n')
    output.write(str([True_TrackDirection_yz[i] for i in range(len(True_TrackDirection_yz))]))
    output.write('\n')
    output.write(str([Reco_Charge[i] for i in range(len(Reco_Charge))]))
    output.write('\n')
    output.write(str([True_Charge[i] for i in range(len(True_Charge))]))
    output.write('\n')
    output.write(str([True_Charge_stop[i] for i in range(len(True_Charge_stop))]))
    output.write('\n')
    output.write(str([True_KE[i] for i in range(len(True_KE))]))
    output.write('\n')
    output.write(str([True_ZDirection[i] for i in range(len(True_ZDirection))]))
    output.write('\n')
    output.write(str([Reco_ZDirection[i] for i in range(len(Reco_ZDirection))]))
    output.write('\n')
    output.write(str([HitsPerTrack[i] for i in range(len(HitsPerTrack))]))
    output.close()
        
    return

def dir_to_angle(xy, z):
    return np.degrees(np.arctan(xy / z))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--Start', "-s", help = "Do you want to plot the Start point resolution? Yes -> --Start, No -> --no-Start", action = argparse.BooleanOptionalAction)
    parser.add_argument('--End', "-e", help = "Do you want to plot the End point resolution? Yes -> --End, No -> --no-End", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Charge', "-c", help = "Do you want to plot the Charge identification efficiency plots? Yes -> --Charge, No -> --no-Charge", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Angle', "-a", help = "Do you want to plot the Angular resolution? Yes -> --Angle, No -> --no-Angle", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Contained', "-con", help = "Do you want to plot the end resolution only for contained muons or all? Yes -> --Contained, No -> --no-Contained", action = argparse.BooleanOptionalAction)
    parser.add_argument('--output', type = str, help = "The output file name.")
    
    args = parser.parse_args()
    
    #print(layer_dict)
    out_dir = args.outdir
    input_filename = args.input_filename
    plot_Start = args.Start
    plot_End = args.End
    plot_Charge = args.Charge
    plot_Angle = args.Angle
    Contained = args.Contained
    output = args.output
    draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained, output)
