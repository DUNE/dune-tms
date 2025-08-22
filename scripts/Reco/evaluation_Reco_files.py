"""
- This script takes reconstructed files as input, loops through the spills, events and tracks and extracts the
Start and End position of the reconstructed tracks, as well as the kinetic energy, start direction and charge.
- This is done for the reconstruction and the truth allowing a comparison between the two.
- The extracted values are output into txt files for further analysis with plotting_RecoEvaluation_files.py.

- Usage: 
python3 evaluation_Reco_files.py --outdir /output/directory/ --input_filename /input/file.root --Start --End --Charge --Angle --Contained --output /file/name/for/the/output
- As only one file at a time is read in, either combine all reco files into one, or run it in a bash script looping through the files in a directory.

- Requirements: ROOT, numpy
-> at the gpvms:
$ source setup_FNAL.sh (in dune-tms)
$ spack load py-packaging
$ spack load py-numpy

- The output file(s) include [ and ] at the start and end of each line. To delete them: 
sed -i 's/[][]//g'
afterwards or in the bash script as well.
"""
import ROOT
import numpy as np
import os
import argparse
import cppyy.ll

import math # filter out nans

import csv
#import matplotlib.pyplot as mp

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
        StartX = np.ones((n_events, 5), dtype = float) * -9999.
    if plot_End:
        Reco_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
        Primary_True_End = np.ones((n_events, 5, 3), dtype = float) * -9999.
        #Reco_End_EscX = np.ones((n_events, 5, 3), dtype = float) * -9999.
        #Reco_End_EscY = np.ones((n_events, 5, 3), dtype = float) * -9999.
        #Primary_True_End_EscX = np.ones((n_events, 5, 3), dtype = float) * -9999.
        #Primary_True_End_EscY = np.ones((n_events, 5, 3), dtype = float) * -9999.
        EndX = np.ones((n_events, 5), dtype = float) * -9999.
    if plot_Charge:
        Reco_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_Charge_stop = np.ones((n_events, 5), dtype = float) * -9999.
        ChargeX = np.ones((n_events, 5), dtype = float) * -9999.
    if plot_Charge or plot_Angle:
        True_Charge = np.ones((n_events, 5), dtype = float) * -9999.
        True_KE = np.ones((n_events, 5), dtype = float) * -9999.
        True_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.
        Reco_TrackDirection = np.ones((n_events, 5, 2), dtype = float) * -9999.
        AngleX = np.ones((n_events, 5), dtype = float) * -9999.
    HitsPerTrack = np.ones((n_events, 5), dtype = float) * -9999.

    ## reco efficiency
    #TRUETrack = np.ones((max_n_spills, 200), dtype = float) * -9999.
    #TrueTrack = np.ones((n_events, 5), dtype = float) * -9999.

    # energy resolution
    True_Muon_Energy = np.ones((n_events, 5), dtype = float) * -9999.
    Reco_Muon_Energy = np.ones((n_events, 5), dtype = float) * -9999.
    AltReco_Energy = np.ones((n_events, 5), dtype = float) * -9999.
    EnergyX = np.ones((n_events, 5), dtype = float) * -9999.
    True_Muon_Energy_sameZ = np.ones((n_events, 5), dtype = float) * -9999.
    Reco_Muon_Energy_sameZ = np.ones((n_events, 5), dtype = float) * -9999.
    AltReco_Energy_sameZ = np.ones((n_events, 5), dtype = float) * -9999.
    EnergyStrictX = np.ones((n_events, 5), dtype = float) * -9999.
    
    # now fill the arrays
    for current_spill_number in range(max_n_spills):
    #    # true tracks for reco efficiency
    #    true_track_it = 0
    #    for j in range(truth.nTrueParticles):
    #        # true muon?
    #        ismuon = False
    #        if (np.abs(truth.PDG[j]) == 13):
    #            ismuon = True
    #        # muon touching TMS?
    #        tms_touch = False
    #        if (truth.TrueVisibleEnergy[j] >= 5):
    #            tms_touch = True # MeV

    #        if (ismuon and tms_touch):
    #            TRUETrack[current_spill_number, true_track_it] = truth.MomentumTMSStart[j*4 + 3] * 1e-3
    #            true_track_it += 1
    #            if true_track_it > 200: break 


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
            True_Momentum_Leaving = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumLeavingTMS, dtype = np.float32)

            # reco efficiency and energy resolution
            MomentumEnteringTMS = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueMomentumEnteringTMS, dtype = np.float32)
            length_to_use = event.Length
            alt_length = event.Length_3D
            TrueFiducialEnd = true_event.RecoTrackPrimaryParticleTMSFiducialEnd
            #RecoVisibleEnergy = np.frombuffer(true_event.RecoTrackPrimaryParticleTrueVisibleEnergy, dtype = np.float32)

            nHits = np.frombuffer(event.nHits, dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])
            
            # X hits in reconstructed track
            BarTypes = np.frombuffer(event.TrackHitBarType, dtype = np.uint8)

            if nTracks <= 0: continue
            #if nTracks > 4: print("Too many tracks in event. Limit to first 5")
            for j in range(nTracks):
                if j > 4: break

                # check if X hits in reconstructed track
                ContainsX = False
                for hit in range(nHits[j]):
                    if BarTypes[j*800 + hit*4] == 0: # 0 == XBar
                        ContainsX = True
                        break
                
                if True_Position_TMS_Start[j*4 + 0] > -8000. and not StartPos.size == 0:
                    # checking for muon tracks (minimal length for this are 20 planes traversed -> 890 mm in thin area
                    if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # reco efficiency
                            #ismuon = False
                            #if (np.abs(True_PDG[j]) == 13):
                            #    ismuon = True
                            #tms_touch = False
                            #if (RecoVisibleEnergy[j] >= 5):
                            #    tms_touch = True # MeV
                            #if (ismuon and tms_touch):
                            #    TrueTrack[i, j] = MomentumEnteringTMS[j*4 + 3] * 1e-3

                            # energy resolution
                            if (np.abs(True_PDG[j]) == 13):
                                if (TrueFiducialEnd[j]):
                                    EnergyX[i, j] = ContainsX

                                    true_muon_starting_ke = MomentumEnteringTMS[j*4 + 3] * 1e-3
                                    True_Muon_Energy[i, j] = true_muon_starting_ke
    
                                    reco_muon_starting_ke = (length_to_use[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                    Reco_Muon_Energy[i, j] = reco_muon_starting_ke
                                    alt_starting_ke = (alt_length[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                    AltReco_Energy[i, j] = alt_starting_ke
                                    
                                    # roughly same end z for reco and truth
                                    if (True_Position_TMS_End[j*4 + 2] <= 14435):
                                        if np.abs(EndPos[j*3 + 2] - True_Position_TMS_End[j*4 + 2]) <= 65.:
                                            True_Muon_Energy_sameZ[i, j] = MomentumEnteringTMS[j*4 + 3] * 1e-3
                                            Reco_Muon_Energy_sameZ[i, j] = (length_to_use[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            AltReco_Energy_sameZ[i, j] = (alt_length[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            EnergyStrictX[i, j] = ContainsX
                                    elif (True_Position_TMS_End[j*4 + 2] >= 14525 and True_Position_TMS_End[j*4 + 2] <= 17495):
                                        if np.abs(EndPos[j*3 + 2] - True_Position_TMS_End[j*4 + 2]) <= 90.:
                                            True_Muon_Energy_sameZ[i, j] = MomentumEnteringTMS[j*4 + 3] * 1e-3
                                            Reco_Muon_Energy_sameZ[i, j] = (length_to_use[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            AltReco_Energy_sameZ[i, j] = (alt_length[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            EnergyStrictX[i, j] = ContainsX
                                    elif (True_Position_TMS_End[j*4 + 2] > 17495):
                                        if np.abs(EndPos[j*3 + 2] - True_Position_TMS_End[j*4 + 2]) <= 130.:
                                            True_Muon_Energy_sameZ[i, j] = MomentumEnteringTMS[j*4 + 3] * 1e-3
                                            Reco_Muon_Energy_sameZ[i, j] = (length_to_use[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            AltReco_Energy_sameZ[i, j] = (alt_length[j] * 1.75 * 0.9428 + 18.73) * 1e-3
                                            EnergyStrictX[i, j] = ContainsX

                            # normalization of true direction
                            if plot_Angle or plot_Charge:
                                magnitude = np.sqrt(MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 + 1]**2 + MomentumTrackStart[j*4 + 2]**2)
                                magnitude_r = np.sqrt(Reco_Track_StartDirection[j*3 + 0]**2 + Reco_Track_StartDirection[j*3 + 1]**2 + Reco_Track_StartDirection[j*3 + 2]**2)
                                True_TrackDirection[i, j, 0] = dir_to_angle(MomentumTrackStart[j*4 + 0] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                True_TrackDirection[i, j, 1] = dir_to_angle(MomentumTrackStart[j*4 + 1] / magnitude, MomentumTrackStart[j*4 + 2] / magnitude)
                                Reco_TrackDirection[i, j, 0] = dir_to_angle(Reco_Track_StartDirection[j*3 + 0], Reco_Track_StartDirection[j*3 + 2])
                                Reco_TrackDirection[i, j, 1] = dir_to_angle(Reco_Track_StartDirection[j*3 + 1], Reco_Track_StartDirection[j*3 + 2])
                                True_KE[i, j] = np.sqrt((MomentumTrackStart[j*4 + 0]**2 + MomentumTrackStart[j*4 +1]**2 + MomentumTrackStart[j*4 + 2]**2 + 105.7**2) - 105.7)
                                True_Charge[i, j] = True_PDG[j]
                                AngleX[i, j] = ContainsX
                            if plot_Charge:
                                # Make sure that true muon stopped in TMS and is not leaving
                                if TrueFiducialEnd[j]:  #np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                    Reco_Charge[i, j] = Reco_Track_Charge[j]
                                    True_Charge_stop[i, j] = True_PDG[j]
                                    ChargeX[i, j] = ContainsX

                            if plot_Start:
                                Reco_Start[i, j, 0] = StartPos[j*3 + 0]
                                Reco_Start[i, j, 1] = StartPos[j*3 + 1]
                                Reco_Start[i, j, 2] = StartPos[j*3 + 2]
                                Primary_True_Start[i, j, 0] = True_Position_TMS_Start[j*4 + 0]
                                Primary_True_Start[i, j, 1] = True_Position_TMS_Start[j*4 + 1]
                                Primary_True_Start[i, j, 2] = True_Position_TMS_Start[j*4 + 2]
                                StartX[i, j] = ContainsX
                if plot_End and Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
                            if TrueFiducialEnd[j]:  #np.sqrt(True_Momentum_Leaving[j*4 + 0]**2 + True_Momentum_Leaving[j*4 + 1]**2 + True_Momentum_Leaving[j*4 + 2]**2) < 1:
                                Reco_End[i, j, 0] = EndPos[j*3 + 0]
                                Reco_End[i, j, 1] = EndPos[j*3 + 1]
                                Reco_End[i, j, 2] = EndPos[j*3 + 2]
                                Primary_True_End[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                                Primary_True_End[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                                Primary_True_End[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                                HitsPerTrack[i, j] = event.nHits[j]
                                EndX[i, j] = ContainsX

                if plot_End and not Contained:
                    if True_Position_TMS_End[j*4 + 0] > -8000. and not EndPos.size == 0:
                        # Make sure that at least 4 potential hits in TMS
                        if (EndPos[j*3 + 2] - StartPos[j*3 + 2]) > 14 * 65. and (True_Position_TMS_End[j*4 + 2] - True_Position_TMS_Start[j*4 + 2]) > 14 * 65.:
                            # Make sure that true muon stopped in TMS and is not leaving
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
                            #if np.abs(True_Position_TMS_End[j*4 + 0]) >= 3219:
                            #    Reco_End_EscX[i, j, 0] = EndPos[j*3 + 0]
                            #    Reco_End_EscX[i, j, 1] = EndPos[j*3 + 1]
                            #    Reco_End_EscX[i, j, 2] = EndPos[j*3 + 2]
                            #    Primary_True_End_EscX[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                            #    Primary_True_End_EscX[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                            #    Primary_True_End_EscX[i, j, 2] = True_Position_TMS_End[j*4 + 2]
                            #elif True_Position_TMS_End[j*4 + 1] >= 271.77 or True_Position_TMS_End[j*4 + 1] <= -2776.23:
                            #    Reco_End_EscY[i, j, 0] = EndPos[j*3 + 0]
                            #    Reco_End_EscY[i, j, 1] = EndPos[j*3 + 1]
                            #    Reco_End_EscY[i, j, 2] = EndPos[j*3 + 2]
                            #    Primary_True_End_EscY[i, j, 0] = True_Position_TMS_End[j*4 + 0]
                            #    Primary_True_End_EscY[i, j, 1] = True_Position_TMS_End[j*4 + 1]
                            #    Primary_True_End_EscY[i, j, 2] = True_Position_TMS_End[j*4 + 2]
    
    boolean_True_KE = (True_KE != -9999.)
    True_KE = True_KE[boolean_True_KE]
    boolean_Reco_Start = (Reco_Start[:, :, 0] != -9999.)
    Reco_Start = Reco_Start[boolean_Reco_Start]
    Primary_True_Start = Primary_True_Start[boolean_Reco_Start]
    StartX = StartX[boolean_Reco_Start]
    Reco_End = Reco_End[boolean_Reco_Start]
    Primary_True_End = Primary_True_End[boolean_Reco_Start]
    EndX = EndX[boolean_Reco_Start]
    Reco_Charge = Reco_Charge[boolean_True_KE]
    True_Charge = True_Charge[boolean_True_KE]
    ChargeX = ChargeX[boolean_True_KE]
    True_Charge_stop = True_Charge_stop[boolean_True_KE]
    True_TrackDirection = True_TrackDirection[boolean_Reco_Start]
    Reco_TrackDirection = Reco_TrackDirection[boolean_Reco_Start]
    AngleX = AngleX[boolean_Reco_Start]
    HitsPerTrack = HitsPerTrack[boolean_True_KE]
    #boolean_TRUETrack = (TRUETrack != -9999.)
    #TRUETrack = TRUETrack[boolean_TRUETrack]
    #boolean_TrueTrack = (TrueTrack != -9999.)
    #TrueTrack = TrueTrack[boolean_TrueTrack]
    boolean_energy = (True_Muon_Energy != -9999.)
    True_Muon_Energy = True_Muon_Energy[boolean_energy]
    Reco_Muon_Energy = Reco_Muon_Energy[boolean_energy]
    AltReco_Energy = AltReco_Energy[boolean_energy]
    EnergyX = EnergyX[boolean_energy]
    boolean_energy_sameZ = (True_Muon_Energy_sameZ != -9999.)
    True_Muon_Energy_sameZ = True_Muon_Energy_sameZ[boolean_energy_sameZ]
    Reco_Muon_Energy_sameZ = Reco_Muon_Energy_sameZ[boolean_energy_sameZ]
    AltReco_Energy_sameZ = AltReco_Energy_sameZ[boolean_energy_sameZ]
    EnergyStrictX = EnergyStrictX[boolean_energy_sameZ]
    
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
    output.write(str([StartX[i] for i in range(len(StartX))]))
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
    output.write(str([HitsPerTrack[i] for i in range(len(HitsPerTrack))]))
    output.write('\n')
    output.write(str([EndX[i] for i in range(len(EndX))]))
    output.write('\n')
    output.write(str([Reco_TrackDirection_xz[i] for i in range(len(Reco_TrackDirection_xz))]))
    output.write('\n')
    output.write(str([Reco_TrackDirection_yz[i] for i in range(len(Reco_TrackDirection_yz))]))
    output.write('\n')
    output.write(str([True_TrackDirection_xz[i] for i in range(len(True_TrackDirection_xz))]))
    output.write('\n')
    output.write(str([True_TrackDirection_yz[i] for i in range(len(True_TrackDirection_yz))]))
    output.write('\n')
    output.write(str([True_Charge[i] for i in range(len(True_Charge))]))
    output.write('\n')
    output.write(str([True_KE[i] for i in range(len(True_KE))]))
    output.write('\n')
    output.write(str([AngleX[i] for i in range(len(AngleX))]))
    output.write('\n')
    output.write(str([Reco_Charge[i] for i in range(len(Reco_Charge))]))
    output.write('\n')
    output.write(str([True_Charge_stop[i] for i in range(len(True_Charge_stop))]))
    output.write('\n')
    output.write(str([ChargeX[i] for i in range(len(ChargeX))]))
    output.close()

    #output.write(str([TRUETrack[i] for i in range(len(TRUETrack))]))
    #output.write('\n')
    #output.write(str([TrueTrack[i] for i in range(len(TrueTrack))]))
    #output.write('\n')
    output1 = open('%s/energy_%s.txt' % (out_dir, output_name), 'w')
    output1.write(str([True_Muon_Energy[i] for i in range(len(True_Muon_Energy))]))
    output1.write('\n')
    output1.write(str([Reco_Muon_Energy[i] for i in range(len(Reco_Muon_Energy))]))
    output1.write('\n')
    output1.write(str([AltReco_Energy[i] for i in range(len(AltReco_Energy))]))
    output1.write('\n')
    output1.write(str([EnergyX[i] for i in range(len(EnergyX))]))
    output1.close()
    output2 = open('%s/strict_%s.txt' % (out_dir, output_name), 'w')
    output2.write(str([True_Muon_Energy_sameZ[i] for i in range(len(True_Muon_Energy_sameZ))]))
    output2.write('\n')
    output2.write(str([Reco_Muon_Energy_sameZ[i] for i in range(len(Reco_Muon_Energy_sameZ))]))
    output2.write('\n')
    output2.write(str([AltReco_Energy_sameZ[i] for i in range(len(AltReco_Energy_sameZ))]))
    output2.write('\n')
    output2.write(str([EnergyStrictX[i] for i in range(len(EnergyStrictX))]))
    output2.close()
        
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

    out_dir = args.outdir
    input_filename = args.input_filename
    plot_Start = args.Start
    plot_End = args.End
    plot_Charge = args.Charge
    plot_Angle = args.Angle
    Contained = args.Contained
    output = args.output
    draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained, output)
