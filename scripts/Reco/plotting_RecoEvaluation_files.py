import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll

import math # filter out nans
#from scipy.optimize import curve_fit
#from scipy.stats import crystalball
#import scipy.special as sc
#from scipy._lib._util import _lazywhere

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-v0_8-poster')

mp.rc('axes', labelsize = 15)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 15) # fontsize of the tick labels
mp.rc('ytick', labelsize = 15) # fontsize of the tick labels

# dyslexia friendly background
mp.rcParams['axes.facecolor'] = '#fafafa'
mp.rcParams["figure.facecolor"] = '#fafafa'
mp.rcParams["savefig.facecolor"] = '#fafafa'

# off-black text
mp.rcParams['text.color'] = '#424242'
mp.rcParams['axes.labelcolor'] = '#424242'
mp.rcParams['xtick.color'] = '#424242'
mp.rcParams['ytick.color'] = '#424242'

def fill(input_array):
    """Fill arrays in usable manner from output of Auswertung
    
    @param[in] input_array: array read in from output file of Auswertung
    
    @return: 2 1d-arrays containing each one line from output file of Auswertung
    """
    output_1 = input_array[:, 0]
    output_2 = input_array[:, 1]
    output_3 = input_array[:, 2]
    output_4 = input_array[:, 3]
    output_5 = input_array[:, 4]
    output_6 = input_array[:, 5]
    output_7 = input_array[:, 6]
    output_8 = input_array[:, 7]
    output_9 = input_array[:, 8]
    output_10 = input_array[:, 9]
    output_11 = input_array[:, 10]
    output_12 = input_array[:, 11]
    output_13 = input_array[:, 12]
    output_14 = input_array[:, 13]
    output_15 = input_array[:, 14]
    output_16 = input_array[:, 15]
    output_17 = input_array[:, 16]
    output_18 = input_array[:, 17]
    output_19 = input_array[:, 18]
    output_20 = input_array[:, 19]
    output_21 = input_array[:, 20]
    output_22 = input_array[:, 21]
    output_23 = input_array[:, 22]
    output_24 = input_array[:, 23]
    output_25 = input_array[:, 24]
    return output_1, output_2, output_3, output_4, output_5, output_6, output_7, output_8, output_9, output_10, output_11, output_12, output_13, output_14, output_15, output_16, output_17, output_18, output_19, output_20, output_21, output_22, output_23, output_24, output_25

def fill_short(input_array):
    output_1 = input_array[:, 0]
    output_2 = input_array[:, 1]
    output_3 = input_array[:, 2]
    output_4 = input_array[:, 3]
    return output_1, output_2, output_3, output_4

### Actual function that loops through the spills
def draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained):
    # Define empty arrays to get filled in the read-in loop
    Reco_Start_x = np.array([], dtype = float)
    Reco_Start_y = np.array([], dtype = float)
    Reco_Start_z = np.array([], dtype = float)
    Primary_True_Start_x = np.array([], dtype = float)
    Primary_True_Start_y = np.array([], dtype = float)
    Primary_True_Start_z = np.array([], dtype = float)
    StartX = np.array([], dtype = int)
    Reco_End_x = np.array([], dtype = float)
    Reco_End_y = np.array([], dtype = float)
    Reco_End_z = np.array([], dtype = float)
    Primary_True_End_x = np.array([], dtype = float)
    Primary_True_End_y = np.array([], dtype = float)
    Primary_True_End_z = np.array([], dtype = float)
    HitsPerTrack = np.array([], dtype = float)
    EndX = np.array([], dtype = int)
    Reco_TrackDirection_xz = np.array([], dtype = float)
    Reco_TrackDirection_yz = np.array([], dtype = float)
    True_TrackDirection_xz = np.array([], dtype = float)
    True_TrackDirection_yz = np.array([], dtype = float)
    True_Charge = np.array([], dtype = float)
    True_KE = np.array([], dtype = float)
    AngleX = np.array([], dtype = int)
    Reco_Charge = np.array([], dtype = float)
    True_Charge_stop = np.array([], dtype = float)
    ChargeX = np.array([], dtype = int)
    #TRUETracks = np.array([], dtype = float)
    #TrueTracks = np.array([], dtype = float)
    TrueEnergy = np.array([], dtype = float)
    RecoEnergy = np.array([], dtype = float)
    RecoEnergy3D = np.array([], dtype = float)
    EnergyX = np.array([], dtype = int)
    TrueEnergySame = np.array([], dtype = float)
    RecoEnergySame = np.array([], dtype = float)
    RecoEnergy3DSame = np.array([], dtype = float)
    EnergyStrictX = np.array([], dtype = int)
    
    # read in files
    os.chdir(input_filename)
    
    for files in os.listdir():
        # Check whether file is in text format
        if files.endswith('.txt'):
            file_path = '%s' % files
            #print("reading in: ", file_path)
            # Read in data
            #open(file_path, 'r')
            if file_path.startswith('output', 0, 6):
                data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
                
                # Fill arrays in correctly
                RSx, RSy, RSz, TSx, TSy, TSz, SX, REx, REy, REz, TEx, TEy, TEz, HPT, EX, RTxz, RTyz, TTxz, TTyz, TC, TK, AX, RC, TCs, CX = fill(data)
                # Concatenate arrays
                Reco_Start_x = np.concatenate([Reco_Start_x, RSx])
                Reco_Start_y = np.concatenate([Reco_Start_y, RSy])
                Reco_Start_z = np.concatenate([Reco_Start_z, RSz])
                Primary_True_Start_x = np.concatenate([Primary_True_Start_x, TSx])
                Primary_True_Start_y = np.concatenate([Primary_True_Start_y, TSy])
                Primary_True_Start_z = np.concatenate([Primary_True_Start_z, TSz])
                StartX = np.concatenate([StartX, SX])
                Reco_End_x = np.concatenate([Reco_End_x, REx])
                Reco_End_y = np.concatenate([Reco_End_y, REy])
                Reco_End_z = np.concatenate([Reco_End_z, REz])
                Primary_True_End_x = np.concatenate([Primary_True_End_x, TEx])
                Primary_True_End_y = np.concatenate([Primary_True_End_y, TEy])
                Primary_True_End_z = np.concatenate([Primary_True_End_z, TEz])
                HitsPerTrack = np.concatenate([HitsPerTrack, HPT])
                EndX = np.concatenate([EndX, EX])
                Reco_TrackDirection_xz = np.concatenate([Reco_TrackDirection_xz, RTxz])
                Reco_TrackDirection_yz = np.concatenate([Reco_TrackDirection_yz, RTyz])
                True_TrackDirection_xz = np.concatenate([True_TrackDirection_xz, TTxz])
                True_TrackDirection_yz = np.concatenate([True_TrackDirection_yz, TTyz])
                True_Charge = np.concatenate([True_Charge, TC])
                True_KE = np.concatenate([True_KE, TK])
                AngleX = np.concatenate([AngleX, AX])
                Reco_Charge = np.concatenate([Reco_Charge, RC])
                True_Charge_stop = np.concatenate([True_Charge_stop, TCs])
                ChargeX = np.concatenate([ChargeX, CX])

            #elif file_path.startswith('res_t', 0, 5):
            #    data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
            #    #TT = data

            #    TRUETracks = np.concatenate([TRUETracks, data])
            #elif file_path.startswith('res_r', 0, 5):
            #    data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
            #    #Tt = data

            #    TrueTracks = np.concatenate([TrueTracks, data])
            elif file_path.startswith('energy', 0, 6):
                data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
                TE, RE, RE3D, EX = fill_short(data)

                TrueEnergy = np.concatenate([TrueEnergy, TE])
                RecoEnergy = np.concatenate([RecoEnergy, RE])
                RecoEnergy3D = np.concatenate([RecoEnergy3D, RE3D])
                EnergyX = np.concatenate([EnergyX, EX])
            elif file_path.startswith('strict', 0, 6):
                data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
                TE, RE, RE3D, EX = fill_short(data)

                TrueEnergySame = np.concatenate([TrueEnergySame, TE])
                RecoEnergySame = np.concatenate([RecoEnergySame, RE])
                RecoEnergy3DSame = np.concatenate([RecoEnergy3DSame, RE3D])
                EnergyStrictX = np.concatenate([EnergyStrictX, EX])
            
    # filter out not filled indice
    if plot_Start:
        boolean_Reco_Start = (Reco_Start_x != -9999.)
        Reco_Start_x = Reco_Start_x[boolean_Reco_Start]
        Reco_Start_y = Reco_Start_y[boolean_Reco_Start]
        Reco_Start_z = Reco_Start_z[boolean_Reco_Start]
        boolean_Primary_Start = (Primary_True_Start_x != -9999.)
        Primary_True_Start_x = Primary_True_Start_x[boolean_Primary_Start]
        Primary_True_Start_y = Primary_True_Start_y[boolean_Primary_Start]
        Primary_True_Start_z = Primary_True_Start_z[boolean_Primary_Start]
        StartX = StartX[boolean_Reco_Start].astype(bool)
    if plot_End:
        boolean_Reco_End = (Reco_End_x != -9999.)
        Reco_End_x = Reco_End_x[boolean_Reco_End]
        Reco_End_y = Reco_End_y[boolean_Reco_End]
        Reco_End_z = Reco_End_z[boolean_Reco_End]
        boolean_Primary_End = (Primary_True_End_x != -9999.)
        Primary_True_End_x = Primary_True_End_x[boolean_Primary_End]
        Primary_True_End_y = Primary_True_End_y[boolean_Primary_End]
        Primary_True_End_z = Primary_True_End_z[boolean_Primary_End]
        HitsPerTrack = HitsPerTrack[boolean_Reco_End]
        EndX = EndX[boolean_Reco_End].astype(bool)
    if plot_Charge:
        boolean_Reco_Charge = (Reco_Charge != -9999.)
        Reco_Charge = Reco_Charge[boolean_Reco_Charge]
        boolean_True_Charge_stop = (True_Charge_stop != -9999.)
        True_Charge_stop = True_Charge_stop[boolean_True_Charge_stop]
        True_KE_stop = True_KE[boolean_True_Charge_stop & (True_KE != -9999.)]
        ChargeX = ChargeX[boolean_Reco_Charge].astype(bool)
    if plot_Charge or plot_Angle:
        boolean_True_Charge = (True_Charge != -9999.)
        True_Charge = True_Charge[boolean_True_Charge]    
        boolean_True_TrackDirection = (True_TrackDirection_xz != -9999.)
        True_TrackDirection_xz = True_TrackDirection_xz[boolean_True_TrackDirection]
        True_TrackDirection_yz = True_TrackDirection_yz[boolean_True_TrackDirection]
        boolean_Reco_TrackDirection = (Reco_TrackDirection_xz != -9999.)
        Reco_TrackDirection_xz = Reco_TrackDirection_xz[boolean_Reco_TrackDirection]
        Reco_TrackDirection_yz = Reco_TrackDirection_yz[boolean_Reco_TrackDirection]
        boolean_True_KE = (True_KE != -9999.)
        True_KE = True_KE[boolean_True_KE]
        AngleX = AngleX[boolean_Reco_TrackDirection].astype(bool)
    #boolean_TRUETracks = (TRUETracks != -9999.)
    #TRUETracks = TRUETracks[boolean_TRUETracks]
    #TrueTracks = TrueTracks[boolean_TRUETracks]
    boolean_energy = (RecoEnergy != -9999.)
    TrueEnergy = TrueEnergy[boolean_energy]
    RecoEnergy = RecoEnergy[boolean_energy]
    RecoEnergy3D = RecoEnergy3D[boolean_energy]
    EnergyX = EnergyX[boolean_energy].astype(bool)
    boolean_energy_strict = (RecoEnergySame != -9999.)
    TrueEnergySame = TrueEnergySame[boolean_energy_strict]
    RecoEnergySame = RecoEnergySame[boolean_energy_strict]
    RecoEnergy3DSame = RecoEnergy3DSame[boolean_energy_strict]
    EnergyStrictX = EnergyStrictX[boolean_energy_strict].astype(bool)
    
    # total number of events after filtering
    if plot_End:
        print("------------------------------------------")
        print("     #events reconstruction: ", len(Reco_End_x))
        print("------------------------------------------")
    
    # subtract reconstruction from truth for all directions
    if plot_End:
        ### x
        # getting rid of nan's in arrays
        boolean_EndX_nan = ~np.isnan(Reco_End_x) & ~np.isnan(Primary_True_End_x)
        Reco_End_x = Reco_End_x[boolean_EndX_nan]
        Primary_True_End_x = Primary_True_End_x[boolean_EndX_nan]
        EndX_x = EndX[boolean_EndX_nan]
        # actual difference
        Diff_End_x = Reco_End_x - Primary_True_End_x
        # now taking out very unreasonable values
        sanity_Endx = np.abs(Reco_End_x) < 4000.
        EndX_x = EndX[sanity_Endx]
        Diff_End_x = Diff_End_x[sanity_Endx]

        ### y
        # getting rid of nan's in arrays
        boolean_EndY_nan = ~np.isnan(Reco_End_y) & ~np.isnan(Primary_True_End_y)
        Reco_End_y = Reco_End_y[boolean_EndY_nan]
        Primary_True_End_y = Primary_True_End_y[boolean_EndY_nan]
        EndX_y = EndX[boolean_EndY_nan]
        # actual difference
        Diff_End_y = Reco_End_y - Primary_True_End_y
        # now taking out very unreasonable values
        sanity_Endy = (Reco_End_y < 500.) & (Reco_End_y > -3500.)
        EndX_y = EndX[sanity_Endy]
        Diff_End_y = Diff_End_y[sanity_Endy]

        ### z
        Diff_End_z = Reco_End_z - Primary_True_End_z
        EndX_z = EndX[~np.isnan(Diff_End_z)]
        Diff_End_z = Diff_End_z[~np.isnan(Diff_End_z)]

    if plot_Start:
        ### x
        # getting rid of nan's in arrays
        boolean_StartX_nan = ~np.isnan(Reco_Start_x) & ~np.isnan(Primary_True_Start_x)
        Reco_Start_x = Reco_Start_x[boolean_StartX_nan]
        Primary_True_Start_x = Primary_True_Start_x[boolean_StartX_nan]
        StartX_x = StartX[boolean_StartX_nan]
        # actual difference
        Diff_Start_x = Reco_Start_x - Primary_True_Start_x
        # now taking out very unreasonable values
        sanity_Startx = np.abs(Reco_Start_x) < 4000.
        StartX_x = StartX[sanity_Startx]
        Diff_Start_x = Diff_Start_x[sanity_Startx]
    
        ### y
        # getting rid of nan's in arrays
        boolean_StartY_nan = ~np.isnan(Reco_Start_y) & ~np.isnan(Primary_True_Start_y)
        Reco_Start_y = Reco_Start_y[boolean_StartY_nan]
        Primary_True_Start_y = Primary_True_Start_y[boolean_StartY_nan]
        StartX_y = StartX[boolean_StartY_nan]
        # actual difference
        Diff_Start_y = Reco_Start_y - Primary_True_Start_y
        # now taking out very unreasonable values
        sanity_Starty = (Reco_Start_y < 500.) & (Reco_Start_y > -3500.)
        StartX_y = StartX[sanity_Starty]
        Diff_Start_y = Diff_Start_y[sanity_Starty]
    
        ### z
        Diff_Start_z = Reco_Start_z - Primary_True_Start_z
        StartX_z = StartX[~np.isnan(Diff_Start_z)]
        Diff_Start_z = Diff_Start_z[~np.isnan(Diff_Start_z)]

    # create histograms for the differences
    if plot_Start:    
        Diff_Start_x_hist, Diff_Start_x_bins = np.histogram(Diff_Start_x, bins = 100, range = (-500, 500))
        Diff_Start_y_hist, Diff_Start_y_bins = np.histogram(Diff_Start_y, bins = 100, range = (-2000, 2000))
        Diff_Start_z_hist, Diff_Start_z_bins = np.histogram(Diff_Start_z, bins = 100, range = (-3000, 2000))
    
        Diff_Start_x_histX, Diff_Start_x_histY = histogram_arr_handle(Diff_Start_x_hist, Diff_Start_x_bins)
        Diff_Start_y_histX, Diff_Start_y_histY = histogram_arr_handle(Diff_Start_y_hist, Diff_Start_y_bins)
        Diff_Start_z_histX, Diff_Start_z_histY = histogram_arr_handle(Diff_Start_z_hist, Diff_Start_z_bins)

        mean_Start_x = np.nanmean(Diff_Start_x)
        mean_Start_y = np.nanmean(Diff_Start_y)
        mean_Start_z = np.nanmean(Diff_Start_z)
        std_Start_x = np.nanstd(Diff_Start_x)
        std_Start_y = np.nanstd(Diff_Start_y)
        std_Start_z = np.nanstd(Diff_Start_z)
        Start_x_q_25 = np.nanquantile(Diff_Start_x, 0.25)
        Start_y_q_25 = np.nanquantile(Diff_Start_y, 0.25)
        Start_z_q_25 = np.nanquantile(Diff_Start_z, 0.25)
        Start_x_q_75 = np.nanquantile(Diff_Start_x, 0.75)
        Start_y_q_75 = np.nanquantile(Diff_Start_y, 0.75)
        Start_z_q_75 = np.nanquantile(Diff_Start_z, 0.75)

        Diff_StartX_x = Diff_Start_x[StartX_x]
        Diff_StartX_y = Diff_Start_y[StartX_y]
        Diff_StartX_z = Diff_Start_z[StartX_z]

        Diff_StartX_x_hist, Diff_StartX_x_bins = np.histogram(Diff_StartX_x, bins = 100, range = (-500, 500))
        Diff_StartX_y_hist, Diff_StartX_y_bins = np.histogram(Diff_StartX_y, bins = 100, range = (-2000, 2000))
        Diff_StartX_z_hist, Diff_StartX_z_bins = np.histogram(Diff_StartX_z, bins = 100, range = (-3000, 2000))
    
        Diff_StartX_x_histX, Diff_StartX_x_histY = histogram_arr_handle(Diff_StartX_x_hist, Diff_StartX_x_bins)
        Diff_StartX_y_histX, Diff_StartX_y_histY = histogram_arr_handle(Diff_StartX_y_hist, Diff_StartX_y_bins)
        Diff_StartX_z_histX, Diff_StartX_z_histY = histogram_arr_handle(Diff_StartX_z_hist, Diff_StartX_z_bins)

        mean_StartX_x = np.nanmean(Diff_StartX_x)
        mean_StartX_y = np.nanmean(Diff_StartX_y)
        mean_StartX_z = np.nanmean(Diff_StartX_z)
        StartX_x_q_25 = np.nanquantile(Diff_StartX_x, 0.25)
        StartX_y_q_25 = np.nanquantile(Diff_StartX_y, 0.25)
        StartX_z_q_25 = np.nanquantile(Diff_StartX_z, 0.25)
        StartX_x_q_75 = np.nanquantile(Diff_StartX_x, 0.75)
        StartX_y_q_75 = np.nanquantile(Diff_StartX_y, 0.75)
        StartX_z_q_75 = np.nanquantile(Diff_StartX_z, 0.75)

        # fit distributions
        #popt_startX, pcov_startX = curve_fit(double_gauss, Diff_Start_x_histX / 10, Diff_Start_x_histY / len(Diff_Start_x), p0 = [Start_x_q_25 / 10, Start_x_q_75 / 10, mean_Start_x / 10, 0.5, 0.5], sigma = np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x), absolute_sigma = False)
        #popt_startY, pcov_startY = curve_fit(double_gauss, Diff_Start_y_histX / 10, Diff_Start_y_histY / len(Diff_Start_y), p0 = [Start_y_q_25 / 10, Start_y_q_75 / 10, mean_Start_y / 10, 0.5, 0.5], sigma = np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y), absolute_sigma = False)
    
        print("------------------------------------------")
        print("     Start   x: ", mean_Start_x, "±", std_Start_x, "    y: ", mean_Start_y, "±", std_Start_y, "    z: ", mean_Start_z, "±", std_Start_z)
        print("     Quantile (25%): x: ", Start_x_q_25, "     y: ", Start_y_q_25, "     z: ", Start_z_q_25)
        print("     Quantile (75%): x: ", Start_x_q_75, "     y: ", Start_y_q_75, "     z: ", Start_z_q_75)
        print("     Only tracks with X hits")
        print("             x: ", mean_StartX_x, "  y: ", mean_StartX_y, "  z: ", mean_StartX_z)
        print("     Quantile (25%): x: ", StartX_x_q_25, "     y: ", StartX_y_q_25, "     z: ", StartX_z_q_25)
        print("     Quantile (75%): x: ", StartX_x_q_75, "     y: ", StartX_y_q_75, "     z: ", StartX_z_q_75)
        #print("     Fit results:")
        #print("         x: ", popt_startX, " ± ", np.sqrt(np.diagonal(pcov_startX)))
        #print("           chi_square: ", chi_square(Diff_Start_x_histY / len(Diff_Start_x), double_gauss(Diff_Start_x_histX / 10, *popt_startX), np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x)), " / ", len(Diff_Start_x_histX) / 2, " : ", chi_square(Diff_Start_x_histY / len(Diff_Start_x), double_gauss(Diff_Start_x_histX / 10, *popt_startX), np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x)) / (len(Diff_Start_x_histX) / 2 - len(popt_startX)))
        #print("         y: ", popt_startY, " ± ", np.sqrt(np.diagonal(pcov_startY)))
        #print("           chi_square: ", chi_square(Diff_Start_y_histY / len(Diff_Start_y), double_gauss(Diff_Start_y_histX / 10, *popt_startY), np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y)), " / ", len(Diff_Start_y_histX) / 2, " : ", chi_square(Diff_Start_y_histY / len(Diff_Start_y), double_gauss(Diff_Start_y_histX / 10, *popt_startY), np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y)) / (len(Diff_Start_y_histX) / 2 - len(popt_startY)))
        output = open('%s/shapes.txt' % out_dir, 'w')
        output.write(str([Diff_Start_x_histX[i] / 10 for i in range(len(Diff_Start_x_histX))]))
        output.write('\n')
        output.write(str([Diff_Start_x_histY[i] / len(Diff_Start_x) for i in range(len(Diff_Start_y_histY))]))
        output.write('\n')
        output.write(str([Diff_Start_y_histX[i] / 10 for i in range(len(Diff_Start_y_histX))]))
        output.write('\n')
        output.write(str([Diff_Start_y_histY[i] / len(Diff_Start_y) for i in range(len(Diff_Start_y_histY))]))
        output.write('\n')
        print("------------------------------------------")

    if plot_End:
        Diff_End_x_hist, Diff_End_x_bins = np.histogram(Diff_End_x, bins = 100, range = (-2000, 2000))
        Diff_End_y_hist, Diff_End_y_bins = np.histogram(Diff_End_y, bins = 100, range = (-3000, 3000))
        Diff_End_z_hist, Diff_End_z_bins = np.histogram(Diff_End_z, bins = 100, range = (-3000, 2000))
   
        Diff_End_x_histX, Diff_End_x_histY = histogram_arr_handle(Diff_End_x_hist, Diff_End_x_bins)
        Diff_End_y_histX, Diff_End_y_histY = histogram_arr_handle(Diff_End_y_hist, Diff_End_y_bins)
        Diff_End_z_histX, Diff_End_z_histY = histogram_arr_handle(Diff_End_z_hist, Diff_End_z_bins)               

        mean_End_x = np.nanmean(Diff_End_x)
        mean_End_y = np.nanmean(Diff_End_y)
        mean_End_z = np.nanmean(Diff_End_z)
        std_End_x = np.nanstd(Diff_End_x)
        std_End_y = np.nanstd(Diff_End_y)
        std_End_z = np.nanstd(Diff_End_z)
        End_x_q_25 = np.nanquantile(Diff_End_x, 0.25)
        End_y_q_25 = np.nanquantile(Diff_End_y, 0.25)
        End_z_q_25 = np.nanquantile(Diff_End_z, 0.25)
        End_x_q_75 = np.nanquantile(Diff_End_x, 0.75)
        End_y_q_75 = np.nanquantile(Diff_End_y, 0.75)
        End_z_q_75 = np.nanquantile(Diff_End_z, 0.75)

        Diff_EndX_x = Diff_End_x[EndX_x]
        Diff_EndX_y = Diff_End_y[EndX_y]
        Diff_EndX_z = Diff_End_z[EndX_z]

        Diff_EndX_x_hist, Diff_EndX_x_bins = np.histogram(Diff_EndX_x, bins = 100, range = (-2000, 2000))
        Diff_EndX_y_hist, Diff_EndX_y_bins = np.histogram(Diff_EndX_y, bins = 100, range = (-3000, 3000))
        Diff_EndX_z_hist, Diff_EndX_z_bins = np.histogram(Diff_EndX_z, bins = 100, range = (-3000, 2000))
   
        Diff_EndX_x_histX, Diff_EndX_x_histY = histogram_arr_handle(Diff_EndX_x_hist, Diff_EndX_x_bins)
        Diff_EndX_y_histX, Diff_EndX_y_histY = histogram_arr_handle(Diff_EndX_y_hist, Diff_EndX_y_bins)
        Diff_EndX_z_histX, Diff_EndX_z_histY = histogram_arr_handle(Diff_EndX_z_hist, Diff_EndX_z_bins)               

        mean_EndX_x = np.nanmean(Diff_EndX_x)
        mean_EndX_y = np.nanmean(Diff_EndX_y)
        mean_EndX_z = np.nanmean(Diff_EndX_z)
        EndX_x_q_25 = np.nanquantile(Diff_EndX_x, 0.25)
        EndX_y_q_25 = np.nanquantile(Diff_EndX_y, 0.25)
        EndX_z_q_25 = np.nanquantile(Diff_EndX_z, 0.25)
        EndX_x_q_75 = np.nanquantile(Diff_EndX_x, 0.75)
        EndX_y_q_75 = np.nanquantile(Diff_EndX_y, 0.75)
        EndX_z_q_75 = np.nanquantile(Diff_EndX_z, 0.75)

        # fit distributions
        #popt_endY, pcov_endY = curve_fit(end_y_gauss, Diff_End_y_histX / 10, Diff_End_y_histY / len(Diff_End_y), p0 = [End_y_q_25 / 10, End_y_q_75 / 10, -130., 20., mean_End_y / 10, 20., 80., 0.8, 1.5, 0.1, 1.], sigma = np.sqrt(Diff_End_y_histY) / len(Diff_End_y), absolute_sigma = False)
        
        print("------------------------------------------")
        print("     End   x: ", mean_End_x, "±", std_End_x, "    y: ", mean_End_y, "±", std_End_y, "    z: ", mean_End_z, "±", std_End_z)
        print("         Quantile (25%): x: ", End_x_q_25, "     y: ", End_y_q_25, "     z: ", End_z_q_25)
        print("         Quantile (75%): x: ", End_x_q_75, "     y: ", End_y_q_75, "     z: ", End_z_q_75)
        print("     Only tracks with X hits")
        print("           x: ", mean_EndX_x, "  y: ", mean_EndX_y, "  z: ", mean_EndX_z)
        print("         Quantile (25%): x: ", EndX_x_q_25, "    y: ", EndX_y_q_25, "    z: ", EndX_z_q_25)
        print("         Quantile (75%): x: ", EndX_x_q_75, "    y: ", EndX_y_q_75, "    z: ", EndX_z_q_75)
        #print("     Fit results:")
        #print("         y: ", popt_endY, " ± ", np.sqrt(np.diagonal(pcov_endY)))
        #print("           chi_square: ", chi_square(Diff_End_y_histY / len(Diff_End_y), end_y_gauss(Diff_End_y_histX / 10, *popt_endY), np.sqrt(Diff_End_y_histY) / len(Diff_End_y)), " / ", len(Diff_End_y_histY) / 2, " : ", chi_square(Diff_End_y_histY / len(Diff_End_y), end_y_gauss(Diff_End_y_histX / 10, *popt_endY), np.sqrt(Diff_End_y_histY) / len(Diff_End_y)) / (len(Diff_End_y_histX) / 2 - len(popt_endY)))
        
        output.write(str([Diff_End_x_histX[i] / 10 for i in range(len(Diff_End_x_histX))]))
        output.write('\n')
        output.write(str([Diff_End_x_histY[i] / len(Diff_End_x) for i in range(len(Diff_End_x_histY))]))
        output.write('\n')
        output.write(str([Diff_End_y_histX[i] / 10 for i in range(len(Diff_End_y_histX))]))
        output.write('\n')
        output.write(str([Diff_End_y_histY[i] / len(Diff_End_y) for i in range(len(Diff_End_y_histY))]))
        output.close()
        print("------------------------------------------")

    # plot
    if plot_Start:
        start_x_array = np.linspace(Diff_Start_x_histX[0] / 10, Diff_Start_x_histX[-1] / 10, len(Diff_Start_x_histX) * 2)
        Diff_StartX_x = Diff_Start_x[StartX_x]
        mp.plot(Diff_Start_x_histX / 10, Diff_Start_x_histY / len(Diff_Start_x), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_x_histX / 10, 0, Diff_Start_x_histY / len(Diff_Start_x), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_x / 10, 0, max(Diff_Start_x_histY / len(Diff_Start_x)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_x / 10))
        mp.axvspan(Start_x_q_25 / 10, Start_x_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_x_q_25 / 10, Start_x_q_75 / 10, (Start_x_q_75 - Start_x_q_25) / 10))
        mp.plot(Diff_StartX_x_histX / 10, Diff_StartX_x_histY / len(Diff_StartX_x), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks w X hits')
        mp.vlines(mean_StartX_x / 10, 0, max(Diff_StartX_x_histY / len(Diff_StartX_x)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm)' % (mean_StartX_x / 10))
        mp.axvspan(StartX_x_q_25 / 10, StartX_x_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((StartX_x_q_75 - StartX_x_q_25) / 10))
        #mp.plot(start_x_array, double_gauss(start_x_array, *popt_startX), color = black_cbf, linestyle = '-', linewidth = 2, label = 'fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('Start x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_X.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        start_y_array = np.linspace(Diff_Start_y_histX[0] / 10, Diff_Start_y_histX[-1] / 10, len(Diff_Start_y_histX) * 2)
        Diff_StartX_y = Diff_Start_y[StartX_y]
        mp.plot(Diff_Start_y_histX / 10, Diff_Start_y_histY / len(Diff_Start_y), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_y_histX / 10, 0, Diff_Start_y_histY / len(Diff_Start_y), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_y / 10, 0, max(Diff_Start_y_histY / len(Diff_Start_y)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_y / 10))
        mp.axvspan(Start_y_q_25 / 10, Start_y_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_y_q_25 / 10, Start_y_q_75 / 10, (Start_y_q_75 - Start_y_q_25) / 10))
        mp.plot(Diff_StartX_y_histX / 10, Diff_StartX_y_histY / len(Diff_StartX_y), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks w X hits')
        mp.vlines(mean_StartX_y / 10, 0, max(Diff_StartX_y_histY / len(Diff_StartX_y)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm)' % (mean_StartX_y / 10))
        mp.axvspan(StartX_y_q_25 / 10, StartX_y_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((StartX_y_q_75 - StartX_y_q_25) / 10))
        #mp.plot(start_y_array, double_gauss(start_y_array, *popt_startY), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('Start y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_Y.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        Diff_StartX_z = Diff_Start_z[StartX_z]
        mp.plot(Diff_Start_z_histX / 10, Diff_Start_z_histY / len(Diff_Start_z), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_z_histX / 10, 0, Diff_Start_z_histY / len(Diff_Start_z), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_z / 10, 0, max(Diff_Start_z_histY / len(Diff_Start_z)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_z / 10))
        mp.axvspan(Start_z_q_25 / 10, Start_z_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_z_q_25 / 10, Start_z_q_75 / 10, (Start_z_q_75 - Start_z_q_25) / 10))
        mp.plot(Diff_StartX_z_histX / 10, Diff_StartX_z_histY / len(Diff_StartX_z), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks w X hits')
        mp.vlines(mean_StartX_z / 10, 0, max(Diff_StartX_z_histY / len(Diff_StartX_z)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm)' % (mean_StartX_z / 10))
        mp.axvspan(StartX_z_q_25 / 10, StartX_z_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((StartX_z_q_75 - StartX_z_q_25) / 10))
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('Start z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_Z.png' % out_dir, bbox_inches = 'tight')
        mp.close()    
    
    if plot_End:
        Diff_EndX_x = Diff_End_x[EndX_x]
        mp.plot(Diff_End_x_histX / 10, Diff_End_x_histY / len(Diff_End_x), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_x_histX / 10, 0, Diff_End_x_histY / len(Diff_End_x), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_x / 10, 0, max(Diff_End_x_histY / len(Diff_End_x)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_x / 10))
        mp.axvspan(End_x_q_25 / 10, End_x_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_x_q_25 / 10, End_x_q_75 / 10, (End_x_q_75 - End_x_q_25) / 10))
        mp.plot(Diff_EndX_x_histX / 10, Diff_EndX_x_histY / len(Diff_EndX_x), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks w X hits')
        mp.vlines(mean_EndX_x / 10, 0, max(Diff_EndX_x_histY / len(Diff_EndX_x)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm)' % (mean_EndX_x / 10))
        mp.axvspan(EndX_x_q_25 / 10, EndX_x_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((EndX_x_q_75 - EndX_x_q_25) / 10)) 
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('End x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_X.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_X.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        end_y_array = np.linspace(Diff_End_y_histX[0] / 10, Diff_End_y_histX[-1] / 10, len(Diff_End_x_histX) * 2)
        Diff_EndX_y = Diff_End_y[EndX_y]
        mp.plot(Diff_End_y_histX / 10, Diff_End_y_histY / len(Diff_End_y), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')   # / 10
        mp.fill_between(Diff_End_y_histX / 10, 0, Diff_End_y_histY / len(Diff_End_y), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_y / 10, 0, max(Diff_End_y_histY / len(Diff_End_y)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_y / 10))
        mp.axvspan(End_y_q_25 / 10, End_y_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_y_q_25 / 10, End_y_q_75 / 10, (End_y_q_75 - End_y_q_25) / 10))
        mp.plot(Diff_EndX_y_histX / 10, Diff_EndX_y_histY / len(Diff_EndX_y), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks w X hits')
        mp.vlines(mean_EndX_y / 10, 0, max(Diff_EndX_y_histY / len(Diff_EndX_y)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm' % (mean_EndX_y / 10))
        mp.axvspan(EndX_y_q_25 / 10, EndX_y_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((EndX_y_q_75 - EndX_y_q_25) / 10))
        #mp.plot(end_y_array, end_y_gauss(end_y_array, *popt_endY), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('End y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Y.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        Diff_EndX_z = Diff_End_z[EndX_z]
        mp.plot(Diff_End_z_histX / 10, Diff_End_z_histY / len(Diff_End_z), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_z_histX / 10, 0, Diff_End_z_histY / len(Diff_End_z), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_z / 10, 0, max(Diff_End_z_histY / len(Diff_End_z)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_z / 10))
        mp.axvspan(End_z_q_25 / 10, End_z_q_75 / 10, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_z_q_25 / 10, End_z_q_75 / 10, (End_z_q_75 - End_z_q_25) / 10))
        mp.plot(Diff_EndX_z_histX / 10, Diff_EndX_z_histY / len(Diff_EndX_z), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only track w X hits')
        mp.vlines(mean_EndX_z / 10, 0, max(Diff_EndX_z_histY / len(Diff_EndX_z)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f cm)' % (mean_EndX_z / 10))
        mp.axvspan(EndX_z_q_25 / 10, EndX_z_q_75 / 10, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f cm)' % ((EndX_z_q_75 - EndX_z_q_25) / 10))
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('Normalized #')
        mp.title('End z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Z.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Z.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        start_end_hist, start_end_bins = np.histogram(HitsPerTrack[~np.isnan(HitsPerTrack)], bins = 50)
        start_end_hist_x, start_end_hist_y = histogram_arr_handle(start_end_hist, start_end_bins)

        mp.plot(start_end_hist_x, start_end_hist_y, color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(start_end_hist_x, 0, start_end_hist_y, alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('# Hits per track')
        mp.ylabel('#')
        mp.savefig('%s/HitsPerTrack.png' % out_dir, bbox_inches = 'tight')
        mp.close()
       
    ### Charge ID
    if plot_Charge:
        # Filter out all events without truth +/-13 as PDG
        boolean_true_muon_stop = (np.abs(True_Charge_stop) == 13)
        True_Charge_stop = True_Charge_stop[boolean_true_muon_stop]
        True_KE_stop = True_KE_stop[boolean_true_muon_stop]
        Reco_Charge = Reco_Charge[boolean_true_muon_stop]
        boolean_issues = (Reco_Charge == -9.99999999e+08)
        not_identified = Reco_Charge[boolean_issues]
    
        ChargeX = ChargeX[boolean_true_muon_stop]
        True_Charge_stopX = True_Charge_stop[ChargeX]
        True_KE_stopX = True_KE_stop[ChargeX]
        Reco_ChargeX = Reco_Charge[ChargeX]
        boolean_issuesX = (Reco_ChargeX == -9.99999999e+08)
        not_identifiedX = Reco_ChargeX[boolean_issuesX]
        
        # Prepare the KE arrays
        Muons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        AMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        FMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        FAMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        Muons_KEX = np.ones(len(True_KE_stopX), dtype = float) * -9999.
        AMuons_KEX = np.ones(len(True_KE_stopX), dtype = float) * -9999.
        FMuons_KEX = np.ones(len(True_KE_stopX), dtype = float) * -9999.
        FAMuons_KEX = np.ones(len(True_KE_stopX), dtype = float) * -9999.
    
        # Counter for both positive, both negative, different reco/truth positive and negative
        true_muons = len(True_Charge_stop[True_Charge_stop == 13])
        true_muonsX = len(True_Charge_stopX[True_Charge_stopX == 13])
        true_antimuons = len(True_Charge_stop[True_Charge_stop == -13])
        true_antimuonsX = len(True_Charge_stopX[True_Charge_stopX == -13])
        reco_muons = len(Reco_Charge[Reco_Charge == 13])
        reco_muonsX = len(Reco_ChargeX[Reco_ChargeX == 13])
        reco_antimuons = len(Reco_Charge[Reco_Charge == -13])
        reco_antimuonsX = len(Reco_ChargeX[Reco_ChargeX == -13])
        counter_true_positive = 0   #truth and reco agree on muon
        counter_true_negative = 0   #truth and reco agree on antimuon
        counter_false_positive = 0  #truth and reco disagree, reco -> muon
        counter_false_negative = 0  #truth and reco disagree, reco -> antimuon
        counter_TPX = 0
        counter_TNX = 0
        counter_FPX = 0
        counter_FNX = 0
    
        True_Muons_KE = True_KE_stop[True_Charge_stop == 13]
        True_Muons_KEX = True_KE_stopX[True_Charge_stopX == 13]
        True_Antimuons_KE = True_KE_stop[True_Charge_stop == -13]
        True_Antimuons_KEX = True_KE_stopX[True_Charge_stopX == -13]
    
        # Compare sign of truth and reco   
        for i in range(len(True_Charge_stop)):
            if True_Charge_stop[i] == Reco_Charge[i] and True_KE_stop[i] >= 500:
                if True_Charge_stop[i] == 13:
                    counter_true_positive += 1
                    Muons_KE[i] = True_KE_stop[i]
                else:
                    counter_true_negative += 1
                    AMuons_KE[i] = True_KE_stop[i]
            elif True_Charge_stop[i] != Reco_Charge[i] and True_KE_stop[i] >= 500:
                if Reco_Charge[i] == 13:
                    counter_false_positive += 1
                    FMuons_KE[i] = True_KE_stop[i]
                else:
                    counter_false_negative += 1
                    FAMuons_KE[i] = True_KE_stop[i]
        
        for i in range(len(True_Charge_stopX)):
            if True_Charge_stopX[i] == Reco_ChargeX[i] and True_KE_stopX[i] >= 500:
                if True_Charge_stopX[i] == 13:
                    counter_TPX += 1
                    Muons_KEX[i] = True_KE_stopX[i]
                else:
                    counter_TNX += 1
                    AMuons_KEX[i] = True_KE_stopX[i]
            elif True_Charge_stopX[i] != Reco_ChargeX[i] and True_KE_stopX[i] >= 500:
                if Reco_ChargeX[i] == 13:
                    counter_FPX += 1
                    FMuons_KEX[i] = True_KE_stopX[i]
                else:
                    counter_FNX += 1
                    FAMuons_KEX[i] = True_KE_stopX[i]

        
        if (counter_true_positive + counter_true_negative + counter_false_positive + counter_false_negative) != len(True_Charge_stop):
            print("     Counters don't add up")
            print("     Length: ", len(True_Charge_stop))

        print("------------------------------------------")        
        print("     Charge ID numbers")
        print("         Not identified:  ", len(not_identified), " (", len(not_identifiedX), ")")
        print("         True muons:      ", counter_true_positive, " (", counter_TPX, ")")
        print("         True antimuons:  ", counter_true_negative, " (", counter_TNX, ")")
        print("         False muons:     ", counter_false_positive, " (", counter_FPX, ")")
        print("         False antimuons: ", counter_false_negative, " (", counter_FNX, ")  here the not identified go into!")
        print("         true_muons: ", true_muons, " true_antimuons: ", true_antimuons, " reco_muons: ", reco_muons, " reco_antimuons: ", reco_antimuons)
        print("         true m (X): ", true_muonsX, "    true am (X): ", true_antimuonsX, " reco m (X): ", reco_muonsX, "    reco am (X): ", reco_antimuonsX)
        print("------------------------------------------")
        
        # Performance evaluation
        if counter_true_positive > 0:
            efficiency_muons = counter_true_positive / (counter_true_positive + (counter_false_negative - len(not_identified)))
            efficiency_muonsX = counter_TPX / (counter_TPX + (counter_FNX - sum(boolean_issuesX)))
            purity_muons = counter_true_positive / (counter_true_positive + counter_false_positive)
            purity_muonsX = counter_TPX / (counter_TPX + counter_FPX)
            print("  Muons (efficiency | purity):     ", efficiency_muons, " | ", purity_muons, " (", efficiency_muonsX, " | ", purity_muonsX, ")")
        if counter_true_negative > 0:
            efficiency_antimuons = counter_true_negative / (counter_true_negative + counter_false_positive)
            efficiency_antimuonsX = counter_TNX / (counter_TNX + counter_FPX)
            purity_antimuons = counter_true_negative / (counter_true_negative + (counter_false_negative - len(not_identified)))
            purity_antimuonsX = counter_TNX / (counter_TNX + (counter_FNX - sum(boolean_issuesX)))
            print("  Antimuons (efficiency | purity): ", efficiency_antimuons, " | ", purity_antimuons, " (", efficiency_antimuonsX, " | ", purity_antimuonsX, ")")

        accuracy_anti_muons = (counter_true_positive + counter_true_negative) / (counter_true_positive + counter_true_negative + counter_false_positive + (counter_false_negative - len(not_identified)))
        accuracy_anti_muonsX = (counter_TPX + counter_TNX) / (counter_TPX + counter_TNX + counter_FPX + (counter_FNX - len(not_identifiedX)))
        print("  Accuracy (both): ", accuracy_anti_muons, " (", accuracy_anti_muonsX, " ± ", accuracy_uncert(counter_TPX, counter_TNX, counter_FPX, counter_FNX - len(not_identifiedX)), ")")

        # Energy dependent      
        Muons_KE = Muons_KE[Muons_KE != -9999.]
        AMuons_KE = AMuons_KE[AMuons_KE != -9999.]
        FMuons_KE = FMuons_KE[FMuons_KE != -9999.]
        FAMuons_KE = FAMuons_KE[FAMuons_KE != -9999.]
        Muons_KEX = Muons_KEX[Muons_KEX != -9999.]
        AMuons_KEX = AMuons_KEX[AMuons_KEX != -9999.]
        FMuons_KEX = FMuons_KEX[FMuons_KEX != -9999.]
        FAMuons_KEX = FAMuons_KEX[FAMuons_KEX != -9999.]
        
        muons_ke_hist, muons_ke_bins = np.histogram(Muons_KE, bins = 60, range = (0, 6000))
        amuons_ke_hist, amuons_ke_bins = np.histogram(AMuons_KE, bins = muons_ke_bins)
        fmuons_ke_hist, fmuons_ke_bins = np.histogram(FMuons_KE, bins = muons_ke_bins)
        famuons_ke_hist, famuons_ke_bins = np.histogram(FAMuons_KE, bins = muons_ke_bins)
        muons_keX_hist, muons_keX_bins = np.histogram(Muons_KEX, bins = 45, range = (0, 4500))
        amuons_keX_hist, amuons_keX_bins = np.histogram(AMuons_KEX, bins = muons_keX_bins)
        fmuons_keX_hist, fmuons_keX_bins = np.histogram(FMuons_KEX, bins = muons_keX_bins)
        famuons_keX_hist, famuons_keX_bins = np.histogram(FAMuons_KEX, bins = muons_keX_bins)
    
        true_muons_ke_hist, muons_ke_bins = np.histogram(True_Muons_KE, bins = muons_ke_bins)
        true_antimuons_ke_hist, amuons_ke_bins = np.histogram(True_Antimuons_KE, bins = amuons_ke_bins)
        true_muons_keX_hist, muons_keX_bins = np.histogram(True_Muons_KEX, bins = muons_keX_bins)
        true_antimuons_keX_hist, amuons_keX_bins = np.histogram(True_Antimuons_KEX, bins = amuons_keX_bins) 
        
        muons_2 = muons_ke_hist[19] + muons_ke_hist[20]
        amuons_2 = amuons_ke_hist[19] + amuons_ke_hist[20]
        fmuons_2 = fmuons_ke_hist[19] + fmuons_ke_hist[20]
        famuons_2 = famuons_ke_hist[19] + famuons_ke_hist[20]
        muonsX_2 = muons_keX_hist[19] + muons_keX_hist[20]
        amuonsX_2 = amuons_keX_hist[19] + amuons_keX_hist[20]
        fmuonsX_2 = fmuons_keX_hist[19] + fmuons_keX_hist[20]
        famuonsX_2 = famuons_keX_hist[19] + famuons_keX_hist[20]

        print("--------------------------------")
        print("     Charge at 2GeV")
        print("         true muons:     ", muons_2, " (", muonsX_2, ")")
        print("         true amuons:    ", amuons_2, " (", amuonsX_2, ")")
        print("         false muons:    ", fmuons_2, " (", fmuonsX_2, ")")
        print("         false amuons:   ", famuons_2, " (", famuonsX_2, ")")
        print("         Muons (efficiency | purity):    ", muons_2 / (muons_2 + famuons_2), " | ", muons_2 / (muons_2 + fmuons_2), " (", muonsX_2 / (muonsX_2 + famuonsX_2), " | ", muonsX_2 / (muonsX_2 + fmuonsX_2), ")")
        print("         Antimuons (efficiency | purity: ", amuons_2 / (amuons_2 + fmuons_2), " | ", amuons_2 / (amuons_2 + famuons_2), " (", amuonsX_2 / (amuonsX_2 + famuonsX_2), " | ", amuonsX_2 / (amuonsX_2 + famuonsX_2), ")")
        print("         Accuracy (both): ", (muons_2 + amuons_2) / (muons_2 + amuons_2 + fmuons_2 + famuons_2), " (", (muonsX_2 + amuonsX_2) / (muonsX_2 + amuonsX_2 + fmuonsX_2 + famuonsX_2), " ± ", accuracy_uncert(muonsX_2, amuonsX_2, fmuonsX_2, famuonsX_2), ")")

        muons_ke_hist_x, muons_ke_hist_y = histogram_arr_handle(muons_ke_hist, muons_ke_bins)
        amuons_ke_hist_x, amuons_ke_hist_y = histogram_arr_handle(amuons_ke_hist, amuons_ke_bins)
        fmuons_ke_hist_x, fmuons_ke_hist_y = histogram_arr_handle(fmuons_ke_hist, fmuons_ke_bins)
        famuons_ke_hist_x, famuons_ke_hist_y = histogram_arr_handle(famuons_ke_hist, famuons_ke_bins)
        true_muons_ke_hist_x, true_muons_ke_hist_y = histogram_arr_handle(true_muons_ke_hist, muons_ke_bins)
        true_antimuons_ke_hist_x, true_antimuons_ke_hist_y = histogram_arr_handle(true_antimuons_ke_hist, amuons_ke_bins)
        
        muons_keX_hist_x, muons_keX_hist_y = histogram_arr_handle(muons_keX_hist, muons_keX_bins)
        amuons_keX_hist_x, amuons_keX_hist_y = histogram_arr_handle(amuons_keX_hist, amuons_keX_bins)
        fmuons_keX_hist_x, fmuons_keX_hist_y = histogram_arr_handle(fmuons_keX_hist, fmuons_keX_bins)
        famuons_keX_hist_x, famuons_keX_hist_y = histogram_arr_handle(famuons_keX_hist, famuons_keX_bins)
        true_muons_keX_hist_x, true_muons_keX_hist_y = histogram_arr_handle(true_muons_keX_hist, muons_keX_bins)
        true_antimuons_keX_hist_x, true_antimuons_keX_hist_y = histogram_arr_handle(true_antimuons_keX_hist, amuons_keX_bins)
        
        output1 = open('%s/charge_shapes.txt' % out_dir, 'w')
        output1.write(str([muons_ke_hist_x[i] for i in range(len(muons_ke_hist_x))]))
        output1.write('\n')
        output1.write(str([muons_ke_hist_y[i] / (muons_ke_hist_y[i] + famuons_ke_hist_y[i]) for i in range(len(muons_ke_hist_y))]))
        output1.write('\n')
        output1.write(str([amuons_ke_hist_y[i] / (amuons_ke_hist_y[i] + fmuons_ke_hist_y[i]) for i in range(len(amuons_ke_hist_y))]))
        output1.write('\n')
        output1.write(str([np.sqrt( famuons_ke_hist_y[i]**2 * muons_ke_hist_y[i] + muons_ke_hist_y[i]**2 * famuons_ke_hist_y[i]) / (muons_ke_hist_y[i] + famuons_ke_hist_y[i])**2 for i in range(len(muons_ke_hist_y))]))
        output1.write('\n')
        output1.write(str([np.sqrt( fmuons_ke_hist_y[i]**2 * amuons_ke_hist_y[i] + amuons_ke_hist_y[i]**2 * fmuons_ke_hist_y[i]) / (muons_ke_hist_y[i] + famuons_ke_hist_y[i])**2 for i in range(len(muons_ke_hist_y))]))
        output1.close()
        print("--------------------------------")

        if counter_true_positive > 0:
            #uncert = 1 / (muons_ke_hist_y + famuons_ke_hist_y)**2 * np.sqrt( (muons_ke_hist_y**3 + famuons_ke_hist_y**3) / (muons_ke_hist_y * famuons_ke_hist_y))
            uncert = 1 / (muons_ke_hist_y + famuons_ke_hist_y)**2 * np.sqrt( famuons_ke_hist_y**2 * muons_ke_hist_y + muons_ke_hist_y**2 * famuons_ke_hist_y)
            mp.plot(muons_ke_hist_x, muons_ke_hist_y / (muons_ke_hist_y + famuons_ke_hist_y), color = blue_cbf, linestyle = '-', linewidth = 2) 
            mp.fill_between(muons_ke_hist_x, muons_ke_hist_y / (muons_ke_hist_y + famuons_ke_hist_y) - uncert, muons_ke_hist_y / (muons_ke_hist_y + famuons_ke_hist_y) + uncert, color = blue_cbf, alpha = 0.3, label = '$\\mu$ correct')
            mp.plot(muons_ke_hist_x, muons_ke_hist_y / counter_true_positive, color = blue_cbf, linestyle = ':', linewidth = 1, label = 'Fraction $\\mu$ (TP)')
        if counter_true_negative > 0:
            #uncert = 1 / (amuons_ke_hist_y + fmuons_ke_hist_y)**2 * np.sqrt( (amuons_ke_hist_y**3 + fmuons_ke_hist_y**3) / (amuons_ke_hist_y * fmuons_ke_hist_y))
            uncert = 1 / (amuons_ke_hist_y + fmuons_ke_hist_y)**2 * np.sqrt( fmuons_ke_hist_y**2 * amuons_ke_hist_y + amuons_ke_hist_y**2 * fmuons_ke_hist_y)
            mp.plot(amuons_ke_hist_x, amuons_ke_hist_y / (amuons_ke_hist_y + fmuons_ke_hist_y), color = orange_cbf, linestyle = '--', linewidth = 2) 
            mp.fill_between(amuons_ke_hist_x, amuons_ke_hist_y / (amuons_ke_hist_y + fmuons_ke_hist_y) - uncert, amuons_ke_hist_y / (amuons_ke_hist_y + fmuons_ke_hist_y) + uncert, color = orange_cbf, alpha = 0.3, label = '$\\bar{\\mu}$ correct')
            mp.plot(amuons_ke_hist_x, amuons_ke_hist_y / counter_true_negative, color = orange_cbf, linestyle = ':', linewidth = 1, label = 'Fraction $\\bar{\\mu}$ (TN)')
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Fraction')
        mp.title('Charge ID efficiency', fontsize = 'xx-large')
        mp.legend(loc = 'center', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Charge_ID_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        if counter_TPX > 0:
            uncertX = 1 / (muons_keX_hist_y + famuons_keX_hist_y)**2 * np.sqrt( famuons_keX_hist_y**2 + muons_keX_hist_y + muons_keX_hist_y**2 + famuons_keX_hist_y)
            mp.plot(muons_keX_hist_x, muons_keX_hist_y / (muons_keX_hist_y + famuons_keX_hist_y), color = blue_cbf, linestyle = '-', linewidth = 2)
            mp.fill_between(muons_keX_hist_x, muons_keX_hist_y / (muons_keX_hist_y + famuons_keX_hist_y) - uncertX, muons_keX_hist_y / (muons_keX_hist_y + famuons_keX_hist_y) + uncertX, color = blue_cbf, alpha = 0.3, label = '$\\mu$ correct')
            mp.plot(muons_keX_hist_x, muons_keX_hist_y / counter_TPX, color = blue_cbf, linestyle = ':', linewidth = 1, label = 'Fraction $\\mu$ (TP)')
        if counter_TNX > 0:
            uncertX = 1 / (amuons_keX_hist_y + fmuons_keX_hist_y)**2 * np.sqrt( fmuons_keX_hist_y**2 * amuons_keX_hist_y + amuons_keX_hist_y**2 * fmuons_keX_hist_y)
            mp.plot(amuons_keX_hist_x, amuons_keX_hist_y / (amuons_keX_hist_y + fmuons_keX_hist_y), color = orange_cbf, linestyle = '--', linewidth = 2)
            mp.fill_between(amuons_keX_hist_x, amuons_keX_hist_y / (amuons_keX_hist_y + fmuons_keX_hist_y) - uncertX, amuons_keX_hist_y / ( amuons_keX_hist_y + fmuons_keX_hist_y) + uncertX, color = orange_cbf, alpha = 0.3, label = '$\\bar{\\mu}$ correct')
            mp.plot(amuons_keX_hist_x, amuons_keX_hist_y / counter_TNX, color = orange_cbf, linestyle = ':', linewidth = 1, label = 'Fraction $\\bar{\\mu}$ (TN)')
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Fraction')
        mp.title('Charge ID efficiency (only tracks with X hits)', fontsize = 'xx-large')
        mp.legend(loc = 'center', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Charge_ID_KEX.png' % out_dir, bbox_inches = 'tight')
        mp.close()

    
    ### Angular resolution
    if plot_Angle:
        # Filter out all events without truth +/-13 as PDG    
        boolean_true_muon = (np.abs(True_Charge) == 13)
        True_KE = True_KE[boolean_true_muon]
        True_TrackDirection_xz_KE = True_TrackDirection_xz[boolean_true_muon]
        True_TrackDirection_yz_KE = True_TrackDirection_yz[boolean_true_muon]
        Reco_TrackDirection_xz_KE = Reco_TrackDirection_xz[boolean_true_muon]
        Reco_TrackDirection_yz_KE = Reco_TrackDirection_yz[boolean_true_muon]

        AngleX_KE = AngleX[boolean_true_muon]
        True_TrackDirection_xz_KEX = True_TrackDirection_xz_KE[AngleX_KE]
        True_TrackDirection_yz_KEX = True_TrackDirection_yz_KE[AngleX_KE]
        Reco_TrackDirection_xz_KEX = Reco_TrackDirection_xz_KE[AngleX_KE]
        Reco_TrackDirection_yz_KEX = Reco_TrackDirection_yz_KE[AngleX_KE]
    
        # calculate differences in angles
        xz_difference = True_TrackDirection_xz - Reco_TrackDirection_xz
        yz_difference = True_TrackDirection_yz - Reco_TrackDirection_yz
        xz_differenceX = True_TrackDirection_xz[AngleX] - Reco_TrackDirection_xz[AngleX]
        yz_differenceX = True_TrackDirection_yz[AngleX] - Reco_TrackDirection_yz[AngleX]
        
        xz_difference = xz_difference[~np.isnan(xz_difference)] #filter out potential nan's
        yz_difference = yz_difference[~np.isnan(yz_difference)] #filter out potential nan's
        xz_differenceX = xz_differenceX[~np.isnan(xz_differenceX)]
        yz_differenceX = yz_differenceX[~np.isnan(yz_differenceX)]
        
        # calculate mean and std of differences
        mean_xz_difference = np.mean(xz_difference)
        mean_yz_difference = np.mean(yz_difference)
        std_xz_difference = np.std(xz_difference)
        std_yz_difference = np.std(yz_difference)
        xz_q_25 = np.quantile(xz_difference, 0.25)
        yz_q_25 = np.quantile(yz_difference, 0.25)
        xz_q_75 = np.quantile(xz_difference, 0.75)
        yz_q_75 = np.quantile(yz_difference, 0.75)

        mean_xz_differenceX = np.mean(xz_differenceX)
        mean_yz_differenceX = np.mean(yz_differenceX)
        xzX_q_25 = np.quantile(xz_differenceX, 0.25)
        yzX_q_25 = np.quantile(yz_differenceX, 0.25)
        xzX_q_75 = np.quantile(xz_differenceX, 0.75)
        yzX_q_75 = np.quantile(yz_differenceX, 0.75)

        print("------------------------------------------")        
        print("     XZ angle difference: ", mean_xz_difference, " ± ", std_xz_difference)
        print("         25% quantile: ", xz_q_25, "  75% quantile: ", xz_q_75)
        print("     YZ angle difference: ", mean_yz_difference, " ± ", std_yz_difference)
        print("         25% quantile: ", yz_q_25, "  75% quantile: ", yz_q_75)
        print("     Only tracks with X hits")
        print("     XZ mean: ", mean_xz_differenceX, " 25%: ", xzX_q_25, " 75%: ", xzX_q_75)
        print("     YZ mean: ", mean_yz_differenceX, " 25%: ", yzX_q_25, " 75%: ", yzX_q_75)
        #print("------------------------------------------")
        
        # calculate histograms of differences
        Diff_xz_hist, Diff_xz_bins = np.histogram(xz_difference, bins = 200, range = (-20, 20))
        Diff_yz_hist, Diff_yz_bins = np.histogram(yz_difference, bins = 200, range = (-75, 75))
        Diff_xzX_hist, Diff_xzX_bins = np.histogram(xz_differenceX, bins = 200, range = (-20, 20))
        Diff_yzX_hist, Diff_yzX_bins = np.histogram(yz_differenceX, bins = 200, range = (-75, 75))
       
        Diff_xz_histX, Diff_xz_histY = histogram_arr_handle(Diff_xz_hist, Diff_xz_bins)
        Diff_yz_histX, Diff_yz_histY = histogram_arr_handle(Diff_yz_hist, Diff_yz_bins)
        Diff_xzX_histX, Diff_xzX_histY = histogram_arr_handle(Diff_xzX_hist, Diff_xzX_bins)
        Diff_yzX_histX, Diff_yzX_histY = histogram_arr_handle(Diff_yzX_hist, Diff_yzX_bins)

        output2 = open("%s/angle_shapes.txt" % out_dir, 'w')
        output2.write(str([Diff_xz_histX[i] for i in range(len(Diff_xz_histX))]))
        output2.write('\n')
        output2.write(str([Diff_xz_histY[i] / len(xz_difference) for i in range(len(Diff_xz_histY))]))
        output2.write('\n')
        output2.write(str([Diff_yz_histX[i] for i in range(len(Diff_yz_histX))]))
        output2.write('\n')
        output2.write(str([Diff_yz_histY[i] / len(yz_difference) for i in range(len(Diff_yz_histY))]))
        output2.close()
        
        # fit distributions
        #popt_xz, pcov_xz = curve_fit(double_gauss_sep, Diff_xz_histX, Diff_xz_histY / len(xz_difference), p0 = [xz_q_25, xz_q_75, mean_xz_difference - 4, mean_xz_difference + 2, 0.5, 0.5], sigma = np.sqrt(Diff_xz_histY) / len(xz_difference), absolute_sigma = False)
        #popt_yz, pcov_yz = curve_fit(double_gauss_sep, Diff_yz_histX, Diff_yz_histY / len(yz_difference), p0 = [yz_q_25, yz_q_75, mean_yz_difference - 8, mean_yz_difference + 1, 0.5, 0.5], sigma = np.sqrt(Diff_yz_histY) / len(yz_difference), absolute_sigma = False)

        #print("     xz fit: ", popt_xz, " ± ", np.sqrt(np.diagonal(pcov_xz)))
        #print("       chi_square: ", chi_square(Diff_xz_histY / len(xz_difference), double_gauss_sep(Diff_xz_histX, *popt_xz), np.sqrt(Diff_xz_histY) / len(xz_difference)), " / ", len(Diff_xz_histX) / 2 - len(popt_xz), " : ", chi_square(Diff_xz_histY / len(xz_difference), double_gauss_sep(Diff_xz_histX, *popt_xz), np.sqrt(Diff_xz_histY) / len(xz_difference)) / (len(Diff_xz_histX) / 2 - len(popt_xz)))
        #print("     yz fit: ", popt_yz, " ± ", np.sqrt(np.diagonal(pcov_yz)))
        #print("       chi_square: ", chi_square(Diff_yz_histY / len(yz_difference), double_gauss_sep(Diff_yz_histX, *popt_yz), np.sqrt(Diff_yz_histY) / len(yz_difference)), " / ", len(Diff_yz_histX) / 2 - len(popt_yz), " : ", chi_square(Diff_yz_histY / len(yz_difference), double_gauss_sep(Diff_yz_histX, *popt_yz), np.sqrt(Diff_yz_histY) / len(yz_difference)) / (len(Diff_yz_histX) / 2 - len(popt_yz)))
        print("------------------------------------------")
        
        # plot
        mp.plot(Diff_xz_histX, Diff_xz_histY / len(xz_difference), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(Diff_xz_histX, 0, Diff_xz_histY / len(xz_difference), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_xz_difference, 0, max(Diff_xz_histY / len(xz_difference)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f°)' % mean_xz_difference)
        mp.axvspan(xz_q_25, xz_q_75, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f°)$\\to\\ Q_3$ (%2.2f°)\n(= %2.2f°)' % (xz_q_25, xz_q_75, xz_q_75 - xz_q_25))
        mp.plot(Diff_xzX_histX, Diff_xzX_histY / len(xz_differenceX), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks with X hits')
        mp.vlines(mean_xz_differenceX, 0, max(Diff_xzX_histY / len(xz_differenceX)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f°)' % mean_xz_differenceX)
        mp.axvspan(xzX_q_25, xzX_q_75, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f°)' % (xzX_q_75 - xzX_q_25))
        #mp.plot(Diff_xz_histX, double_gauss_sep(Diff_xz_histX, *popt_xz), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('Normalized #')
        mp.title('xz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_XZ_angle.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        mp.plot(Diff_yz_histX, Diff_yz_histY / len(yz_difference), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(Diff_yz_histX, 0, Diff_yz_histY / len(yz_difference), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_yz_difference, 0, max(Diff_yz_histY / len(yz_difference)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f°)' % mean_yz_difference)
        mp.axvspan(yz_q_25, yz_q_75, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f°)$\\to\\ Q_3$ (%2.2f°)\n(= %2.2f°)' % (yz_q_25, yz_q_75, yz_q_75 - yz_q_25))
        mp.plot(Diff_yzX_histX, Diff_yzX_histY / len(yz_difference), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks with X hits')
        mp.vlines(mean_yz_differenceX, 0, max(Diff_yzX_histY / len(yz_differenceX)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.2f°)' % mean_yz_differenceX)
        mp.axvspan(yzX_q_25, yzX_q_75, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f°)' % (yzX_q_75 - yzX_q_25))
        #mp.plot(Diff_yz_histX, double_gauss_sep(Diff_yz_histX, *popt_yz), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('Normalized #')
        mp.title('yz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_YZ_angle.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        # energy dependent
        bins = np.histogram_bin_edges(True_KE, bins = 50, range = (500, 5000))
        
        xz_mean_KE = np.zeros(len(bins), dtype = float)
        yz_mean_KE = np.zeros(len(bins), dtype = float)
        xz_std_KE = np.zeros(len(bins), dtype = float)
        yz_std_KE = np.zeros(len(bins), dtype = float)
        xz_q_25_KE = np.zeros(len(bins), dtype = float)
        xz_q_75_KE = np.zeros(len(bins), dtype = float)
        yz_q_25_KE = np.zeros(len(bins), dtype = float)
        yz_q_75_KE = np.zeros(len(bins), dtype = float)

        xz_mean_KEX = np.zeros(len(bins), dtype = float)
        yz_mean_KEX = np.zeros(len(bins), dtype = float)
        xz_q_25_KEX = np.zeros(len(bins), dtype = float)
        yz_q_25_KEX = np.zeros(len(bins), dtype = float)
        xz_q_75_KEX = np.zeros(len(bins), dtype = float)
        yz_q_75_KEX = np.zeros(len(bins), dtype = float)
        
        for i in range(len(bins) - 1):
            boolean_KE = (True_KE >= bins[i]) & (True_KE < bins[i+1])
            True_xz = True_TrackDirection_xz_KE[boolean_KE]
            True_yz = True_TrackDirection_yz_KE[boolean_KE]
            Reco_xz = Reco_TrackDirection_xz_KE[boolean_KE]
            Reco_yz = Reco_TrackDirection_yz_KE[boolean_KE]
            True_KEX = True_KE[AngleX_KE]
            boolean_KEX = (True_KEX >= bins[i]) & (True_KEX < bins[i+1])
            True_xzX = True_TrackDirection_xz_KEX[boolean_KEX]
            True_yzX = True_TrackDirection_yz_KEX[boolean_KEX]
            Reco_xzX = Reco_TrackDirection_xz_KEX[boolean_KEX]
            Reco_yzX = Reco_TrackDirection_yz_KEX[boolean_KEX]
            
            xz_diff = True_xz - Reco_xz
            yz_diff = True_yz - Reco_yz
            xz_diffX = True_xzX - Reco_xzX
            yz_diffX = True_yzX - Reco_yzX
            
            xz_mean_KE[i] = np.nanmean(xz_diff)
            xz_std_KE[i] = np.nanstd(xz_diff)
            yz_mean_KE[i] = np.nanmean(yz_diff)
            yz_std_KE[i] = np.nanstd(yz_diff)
            xz_q_25_KE[i] = np.nanquantile(xz_diff, 0.25)
            xz_q_75_KE[i] = np.nanquantile(xz_diff, 0.75)
            yz_q_25_KE[i] = np.nanquantile(yz_diff, 0.25)
            yz_q_75_KE[i] = np.nanquantile(yz_diff, 0.75)

            xz_mean_KEX[i] = np.nanmean(xz_diffX)
            yz_mean_KEX[i] = np.nanmean(yz_diffX)
            xz_q_25_KEX[i] = np.nanquantile(xz_diffX, 0.25)
            xz_q_75_KEX[i] = np.nanquantile(xz_diffX, 0.75)
            yz_q_25_KEX[i] = np.nanquantile(yz_diffX, 0.25)
            yz_q_75_KEX[i] = np.nanquantile(yz_diffX, 0.75)
    
        # plot
        #mp.errorbar(bins, xz_mean_KE, yerr = xz_std_KE, linestyle = '--', linewidth = 1.5, color = red_cbf)
        mp.plot(bins, xz_mean_KE, linestyle = '--', linewidth = 1.5, color = red_cbf)
        mp.scatter(bins, xz_mean_KE, marker = '.', color = red_cbf, label = 'mean')
        mp.fill_between(bins, xz_q_25_KE, xz_q_75_KE, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1\\ \\to\\ Q_3$')
        mp.plot(bins - 30, xz_mean_KEX, linestyle = ':', linewidth = 1.5, color = black_cbf)
        mp.scatter(bins - 30, xz_mean_KEX, marker = '.', color = black_cbf, label = 'Only tracks with X hits\n(-30 MeV)')
        mp.fill_between(bins - 30, xz_q_25_KEX, xz_q_75_KEX, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$')
        mp.hlines(0, 0, 5000, color = red_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('xz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_XZ_angle_True_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        #mp.errorbar(bins, yz_mean_KE, yerr = yz_std_KE, linestyle = '--', linewidth = 1.5, color = blue_cbf)
        mp.plot(bins, yz_mean_KE, linestyle = '--', linewidth = 1.5, color = blue_cbf)
        mp.scatter(bins, yz_mean_KE, marker = '.', color = blue_cbf, label = 'mean')
        mp.fill_between(bins, yz_q_25_KE, yz_q_75_KE, color = blue_cbf, alpha = 0.2, hatch = '//', label = '$Q_1\\ \\to\\ Q_3$')
        mp.plot(bins - 30, yz_mean_KEX, linestyle = ':', linewidth = 1.5, color = black_cbf)
        mp.scatter(bins - 30, yz_mean_KEX, marker = '.', color = black_cbf, label = 'Only tracks with X hits\n(-30 MeV)')
        mp.fill_between(bins - 30, yz_q_25_KEX, yz_q_75_KEX, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$')
        mp.hlines(0, 0, 5000, color = blue_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('yz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_YZ_angle_True_KE.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        ### reconstruction efficiency
        #boolean_reco_tracks = (TrueTracks != -9999.)
        #reco_tracks = TrueTracks[boolean_reco_tracks]
        #print("RecoTracks Min:", np.min(reco_tracks), " Max: ", np.max(reco_tracks))
        #print("------------------------------------------")
        #print("     total reconstruction efficiency: ", len(reco_tracks), " / ", len(TrueTracks), " : ", len(reco_tracks) / len(TrueTracks))
        #print("------------------------------------------")
        #print('TRUE Min: ', np.min(TRUETracks), ' Max: ', np.max(TRUETracks))

        ## calculate histogram
        #TRUE_tracks_hist, efficiency_bins = np.histogram(TRUETracks, bins = 50)
        #True_tracks_hist, efficiency_bins = np.histogram(TrueTracks, bins = efficiency_bins)

        #TRUE_tracks_x, TRUE_tracks_y = histogram_arr_handle(TRUE_tracks_hist, efficiency_bins)
        #True_tracks_x, True_tracks_y = histogram_arr_handle(True_tracks_hist, efficiency_bins)

        ## plot
        #mp.plot(TRUE_tracks_x, True_tracks_y / TRUE_tracks_y, color = blue_cbf, linestyle = '--', linewidth = 2)
        #mp.fill_between(TRUE_tracks_x, 0, True_tracks_y / TRUE_tracks_y, alpha = 0.6, hatch = '//', color = blue_cbf)
        #mp.xlabel('True Muon KE Entering TMS [GeV]')
        #mp.ylabel('Reconstruction Efficiency')
        #mp.grid(True, linestyle = '--', alpha = 0.2)
        #mp.title('Reco Efficiency vs True TMS-Entering KE')
        #mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        #mp.savefig('%s/Reco_Efficiency.png' % out_dir, bbox_inches = 'tight')
        #mp.close()

        ### Energy resolution
        #print("RecoEnergy Min: ", np.nanmin(RecoEnergy), " Max: ", np.nanmax(RecoEnergy))
        #print("RecoEnergy strict Min: ", np.nanmin(RecoEnergySame), " Max: ", np.nanmax(RecoEnergySame))
        #print("TrueEnergy Min: ", np.nanmin(TrueEnergy), " Max: ", np.nanmax(TrueEnergy))
        #print("TrueEnergy strict Min: ", np.nanmin(TrueEnergySame), " Max: ", np.nanmax(TrueEnergySame))
        difference_energy = RecoEnergy3D - TrueEnergy
        resolution_energy = difference_energy / TrueEnergy        
        #print("Resolution Min: ", np.nanmin(resolution_energy), " Max: ", np.nanmax(resolution_energy))
        difference_energy_strict = RecoEnergy3DSame - TrueEnergySame
        resolution_energy_strict = difference_energy_strict / TrueEnergySame
        #print("Resolution strict Min: ", np.nanmin(resolution_energy_strict), " Max: ", np.nanmax(resolution_energy_strict))

        difference_energyX = RecoEnergy3D[EnergyX] - TrueEnergy[EnergyX]
        resolution_energyX = difference_energyX / TrueEnergy[EnergyX]
        difference_energy_strictX = RecoEnergy3DSame[EnergyStrictX] - TrueEnergySame[EnergyStrictX]
        resolution_energy_strictX = difference_energy_strictX / TrueEnergySame[EnergyStrictX]

        energy_res_mean = np.nanmean(resolution_energy)
        energy_res_std = np.nanmean(resolution_energy)
        energy_res_strict_mean = np.nanmean(resolution_energy_strict)
        energy_res_strict_std = np.nanstd(resolution_energy_strict)
        energy_res_q_25 = np.nanquantile(resolution_energy, 0.25)
        energy_res_q_50 = np.nanquantile(resolution_energy, 0.5)
        energy_res_q_75 = np.nanquantile(resolution_energy, 0.75)
        energy_res_strict_q25 = np.nanquantile(resolution_energy_strict, 0.25)
        energy_res_strict_q50 = np.nanquantile(resolution_energy_strict, 0.5)
        energy_res_strict_q75 = np.nanquantile(resolution_energy_strict, 0.75)

        energy_q_25 = np.nanquantile(resolution_energyX, 0.25)
        energy_q_50 = np.nanquantile(resolution_energyX, 0.5)
        energy_q_75 = np.nanquantile(resolution_energyX, 0.75)
        energy_strict_q_25 = np.nanquantile(resolution_energy_strictX, 0.25)
        energy_strict_q_50 = np.nanquantile(resolution_energy_strictX, 0.5)
        energy_strict_q_75 = np.nanquantile(resolution_energy_strictX, 0.75)

        print("------------------------------------------")
        print("     Energy resolution: ")
        print("         median: ", energy_res_q_50, " (", energy_q_50, ")")
        print("         25% quantile: ", energy_res_q_25, " (", energy_q_25, ")\t75% quantile: ", energy_res_q_75, " (", energy_q_75, ")")
        print("     Strict resolution: ")
        print("         median: ", energy_res_strict_q50, " (", energy_strict_q_50, ")")
        print("         25% quantile: ", energy_res_strict_q25, " (", energy_strict_q_25, ")\t75% quantile: ", energy_res_strict_q75, " (", energy_strict_q_75, ")")

        boolean_resolution = (-0.75 < resolution_energy) & (resolution_energy < 0.5)
        #print(" res (cut) min: ", np.nanmin(resolution_energy[boolean_resolution]), " max: ", np.nanmax(resolution_energy[boolean_resolution]))
        boolean_resolution_strict = (-0.75 < resolution_energy_strict) & (resolution_energy_strict < 0.5)
        #print(" res strict (cut) min: ", np.nanmin(resolution_energy_strict[boolean_resolution_strict]), " max: ", np.nanmax(resolution_energy_strict[boolean_resolution_strict]))
        boolean_resolutionX = (-0.75 < resolution_energyX) & (resolution_energyX < 0.5)
        boolean_resolution_strictX = (-0.75 < resolution_energy_strictX) & (resolution_energy_strictX < 0.5)

        # calculcate histogram
        resolution_energy_hist, resolution_energy_bins = np.histogram(resolution_energy[boolean_resolution], bins = 50)
        resolution_energy_strict_hist, resolution_energy_strict_bins = np.histogram(resolution_energy_strict[boolean_resolution_strict], bins = 50)
        resolution_energyX_hist, resolution_energyX_bins = np.histogram(resolution_energyX[boolean_resolutionX], bins = 50)
        resolution_energy_strictX_hist, resolution_energy_strictX_bins = np.histogram(resolution_energy_strictX[boolean_resolution_strictX], bins = 50)

        res_energy_hist_x, res_energy_hist_y = histogram_arr_handle(resolution_energy_hist, resolution_energy_bins)
        res_energy_strict_hist_x, res_energy_strict_hist_y = histogram_arr_handle(resolution_energy_strict_hist, resolution_energy_strict_bins)
        res_energyX_hist_x, res_energyX_hist_y = histogram_arr_handle(resolution_energyX_hist, resolution_energyX_bins)
        res_energy_strictX_hist_x, res_energy_strictX_hist_y = histogram_arr_handle(resolution_energy_strictX_hist, resolution_energy_strictX_bins)
        
        output3 = open('%s/energy3D_shape.txt' % out_dir, 'w')
        output3.write(str([res_energy_hist_x[i] for i in range(len(res_energy_hist_x))]))
        output3.write('\n')
        output3.write(str([res_energy_hist_y[i] / len(resolution_energy) for i in range(len(res_energy_hist_y))]))
        output3.write('\n')
        output3.write(str([res_energy_strict_hist_x[i] for i in range(len(res_energy_strict_hist_x))]))
        output3.write('\n')
        output3.write(str([res_energy_strict_hist_y[i] / len(resolution_energy_strict) for i in range(len(res_energy_strict_hist_y))]))
        output3.close()

        # plot
        mp.plot(res_energy_hist_x, res_energy_hist_y / len(resolution_energy), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(res_energy_hist_x, 0, res_energy_hist_y / len(resolution_energy), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(energy_res_q_50, 0, np.nanmax(res_energy_hist_y) / len(resolution_energy), color = red_cbf, linestyle = ':', linewidth = 3, label = 'median (%2.2f)' % energy_res_q_50)
        mp.axvspan(energy_res_q_25, energy_res_q_75, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f)$\\to\\ Q_3$ (%2.2f)\n(= %2.2f)' % (energy_res_q_25, energy_res_q_75, energy_res_q_75 - energy_res_q_25))
        mp.plot(res_energyX_hist_x, res_energyX_hist_y / len(resolution_energyX), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks with X hits')
        mp.vlines(energy_q_50, 0, np.nanmax(res_energyX_hist_y) / len(resolution_energyX), color = black_cbf, linestyle = ':', linewidth = 2, label = 'median (%2.2f)' % energy_q_50)
        mp.axvspan(energy_q_25, energy_q_75, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f)' % (energy_q_75 - energy_q_25))
        mp.xlabel('Energy Resolution (Reco - True) / True')
        mp.ylabel('Normalized #')
        mp.title('Energy Resolution')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/EnergyResolution.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        mp.plot(res_energy_strict_hist_x, res_energy_strict_hist_y / len(resolution_energy_strict), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(res_energy_strict_hist_x, 0, res_energy_strict_hist_y / len(resolution_energy_strict), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(energy_res_strict_q50, 0, np.nanmax(res_energy_strict_hist_y) / len(resolution_energy_strict), color = red_cbf, linestyle = ':', linewidth = 3, label = 'median (%2.2f)' % energy_res_strict_q50)
        mp.axvspan(energy_res_strict_q25, energy_res_strict_q75, color = red_cbf, alpha = 0.2, hatch = '//', label = '$Q_1$ (%2.2f)$\\to\\ Q_3$ (%2.2f)\n(= %2.2f)' % (energy_res_strict_q25, energy_res_strict_q75, energy_res_strict_q75 - energy_res_strict_q25))
        mp.plot(res_energy_strictX_hist_x, res_energy_strictX_hist_y / len(resolution_energy_strictX), color = black_cbf, linestyle = '--', linewidth = 2, label = 'Only tracks with X hits')
        mp.vlines(energy_strict_q_50, 0, np.nanmax(res_energy_strictX_hist_y) / len(resolution_energy_strictX), color = black_cbf, linestyle = ':', linewidth = 2, label = 'median (%2.2f)' % energy_strict_q_50)
        mp.axvspan(energy_strict_q_25, energy_strict_q_75, color = black_cbf, alpha = 0.1, hatch = '\\\\', label = '$Q_1\\ \\to\\ Q_3$ (%2.2f)' % (energy_strict_q_75 - energy_strict_q_25))
        mp.xlabel('Energy Resolution (Reco - True) / True')
        mp.ylabel('Normalized #')
        mp.title('Strict Energy Resolution (roughly same z end position)')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.savefig('%s/EnergyResolutionStrict.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        
def double_gauss_sep(x, sigma1, sigma2, mean1, mean2, norm1, norm2):
    return norm1 / np.sqrt(2 * np.pi * sigma1**2) * np.exp(-(x - mean1)**2 / (2 * sigma1**2)) + norm2 / np.sqrt(2 * np.pi * sigma2**2) * np.exp(-(x - mean2)**2 / (2 * sigma2**2))

def double_gauss(x, sigma1, sigma2, mean, norm1, norm2):
    return norm1 / np.sqrt(2 * np.pi * sigma1**2) * np.exp(-(x - mean)**2 / (2 * sigma1**2)) + norm2 / np.sqrt(2 * np.pi * sigma2**2) * np.exp(-(x - mean)**2 / (2 * sigma2**2))

def end_y_gauss(x, sigma1, sigma2, sigma3, sigma4, mean, mean3, mean4, norm1, norm2, norm3, norm4):
    return norm1 / np.sqrt(2 * np.pi * sigma1**2) * np.exp(-(x - mean)**2 / (2 * sigma1**2)) +  norm2 / np.sqrt(2 * np.pi + sigma2**2) * np.exp(-(x - mean)**2 / (2 * sigma2**2)) + norm3 / np.sqrt(2 * np.pi * sigma3**2) * np.exp(-(x - mean3)**2 / (2 * sigma3**2)) + norm4 / np.sqrt(2 * np.pi * sigma4**2) * np.exp(-(x - mean4)**2 / (2 * sigma4**2))

def chi_square(data, fit_y, yerr):
    diff = data - fit_y
    chi_square = (diff / yerr)**2
    return sum(chi_square)

def dir_to_angle(xy, z):
    return np.degrees(np.arctan(xy / z))

def accuracy_uncert(TP, TN, FP, FN):
    return np.sqrt( (TP + TN) * (FP + FN) / (TP + TN + FP + FN)**3 )

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
    parser.add_argument('--input_folder', "-f", type = str, help = "The folder with the files with the events to draw.")
    parser.add_argument('--Start', "-s", help = "Do you want to plot the Start point resolution? Yes -> --Start, No -> --no-Start", action = argparse.BooleanOptionalAction)
    parser.add_argument('--End', "-e", help = "Do you want to plot the End point resolution? Yes -> --End, No -> --no-End", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Charge', "-c", help = "Do you want to plot the Charge identification efficiency plots? Yes -> --Charge, No -> --no-Charge", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Angle', "-a", help = "Do you want to plot the Angular resolution? Yes -> --Angle, No -> --no-Angle", action = argparse.BooleanOptionalAction)
    parser.add_argument('--Contained', "-con", help = "Do you want to plot the end resolution only for contained muons or all? Yes -> --Contained, No -> --no-Contained", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()
    
    #print(layer_dict)
    out_dir = args.outdir
    input_filename = args.input_folder
    plot_Start = args.Start
    plot_End = args.End
    plot_Charge = args.Charge
    plot_Angle = args.Angle
    Contained = args.Contained
    draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained)
