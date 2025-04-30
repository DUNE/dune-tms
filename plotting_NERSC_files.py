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
    output_1 = input_array[:, 0]#[input_array[i, 0] for i in range(len(input_array))]
    output_2 = input_array[:, 1]#[input_array[i, 1] for i in range(len(input_array))]
    output_3 = input_array[:, 2]#[input_array[i, 2] for i in range(len(input_array))]
    output_4 = input_array[:, 3]#[input_array[i, 3] for i in range(len(input_array))]
    output_5 = input_array[:, 4]#[input_array[i, 4] for i in range(len(input_array))]
    output_6 = input_array[:, 5]#[input_array[i, 5] for i in range(len(input_array))]
    output_7 = input_array[:, 6]#[input_array[i, 6] for i in range(len(input_array))]
    output_8 = input_array[:, 7]#[input_array[i, 7] for i in range(len(input_array))]
    output_9 = input_array[:, 8]#[input_array[i, 8] for i in range(len(input_array))]
    output_10 = input_array[:, 9] #[input_array[i, 9] for i in range(len(input_array))]
    output_11 = input_array[:, 10] #[input_array[i, 10] for i in range(len(input_array))]
    output_12 = input_array[:, 11] #[input_array[i, 11] for i in range(len(input_array))]
    output_13 = input_array[:, 12] #[input_array[i, 12] for i in range(len(input_array))]
    output_14 = input_array[:, 13] #[input_array[i, 13] for i in range(len(input_array))]
    output_15 = input_array[:, 14] #[input_array[i, 14] for i in range(len(input_array))]
    output_16 = input_array[:, 15]#[input_array[i, 15] for i in range(len(input_array))]
    output_17 = input_array[:, 16]#[input_array[i, 16] for i in range(len(input_array))]
    output_18 = input_array[:, 17]#[input_array[i, 17] for i in range(len(input_array))]
    output_19 = input_array[:, 18]#[input_array[i, 18] for i in range(len(input_array))]
    output_20 = input_array[:, 19]#[input_array[i, 19] for i in range(len(input_array))]
    output_21 = input_array[:, 20]
    output_22 = input_array[:, 21]
    output_23 = input_array[:, 22]

    return output_1, output_2, output_3, output_4, output_5, output_6, output_7, output_8, output_9, output_10, output_11, output_12, output_13, output_14, output_15, output_16, output_17, output_18, output_19, output_20, output_21, output_22, output_23

### Actual function that loops through the spills
def draw_performance(out_dir, input_filename, plot_Start, plot_End, plot_Charge, plot_Angle, Contained):
    # Define empty arrays to get filled in the read-in loop
    Reco_Start_x = np.array([], dtype = float)
    Reco_Start_y = np.array([], dtype = float)
    Reco_Start_z = np.array([], dtype = float)
    Primary_True_Start_x = np.array([], dtype = float)
    Primary_True_Start_y = np.array([], dtype = float)
    Primary_True_Start_z = np.array([], dtype = float)
    Reco_End_x = np.array([], dtype = float)
    Reco_End_y = np.array([], dtype = float)
    Reco_End_z = np.array([], dtype = float)
    Primary_True_End_x = np.array([], dtype = float)
    Primary_True_End_y = np.array([], dtype = float)
    Primary_True_End_z = np.array([], dtype = float)
    Reco_TrackDirection_xz = np.array([], dtype = float)
    Reco_TrackDirection_yz = np.array([], dtype = float)
    True_TrackDirection_xz = np.array([], dtype = float)
    True_TrackDirection_yz = np.array([], dtype = float)
    Reco_Charge = np.array([], dtype = float)
    True_Charge = np.array([], dtype = float)
    True_Charge_stop = np.array([], dtype = float)
    True_KE = np.array([], dtype = float)
    True_ZDirection = np.array([], dtype = float)
    Reco_ZDirection = np.array([], dtype = float)
    HitsPerTrack = np.array([], dtype = float)
    
    # read in files
    os.chdir(input_filename)
    
    for files in os.listdir():
        # Check whether file is in text format
        if files.endswith('.txt'):
            file_path = '%s' % files
            print("reading in: ", file_path)
            # Read in data
            #open(file_path, 'r')
            data = np.genfromtxt(file_path, dtype = float, delimiter = ', ', unpack = True)
            
            # Fill arrays in correctly
            RSx, RSy, RSz, TSx, TSy, TSz, REx, REy, REz, TEx, TEy, TEz, RTxz, RTyz, TTxz, TTyz, RC, TC, TCs, TK, TZ, RZ, HPT = fill(data)
            
            # Concatenate arrays
            Reco_Start_x = np.concatenate([Reco_Start_x, RSx])
            Reco_Start_y = np.concatenate([Reco_Start_y, RSy])
            Reco_Start_z = np.concatenate([Reco_Start_z, RSz])
            Primary_True_Start_x = np.concatenate([Primary_True_Start_x, TSx])
            Primary_True_Start_y = np.concatenate([Primary_True_Start_y, TSy])
            Primary_True_Start_z = np.concatenate([Primary_True_Start_z, TSz])
            Reco_End_x = np.concatenate([Reco_End_x, REx])
            Reco_End_y = np.concatenate([Reco_End_y, REy])
            Reco_End_z = np.concatenate([Reco_End_z, REz])
            Primary_True_End_x = np.concatenate([Primary_True_End_x, TEx])
            Primary_True_End_y = np.concatenate([Primary_True_End_y, TEy])
            Primary_True_End_z = np.concatenate([Primary_True_End_z, TEz])
            Reco_TrackDirection_xz = np.concatenate([Reco_TrackDirection_xz, RTxz])
            Reco_TrackDirection_yz = np.concatenate([Reco_TrackDirection_yz, RTyz])
            True_TrackDirection_xz = np.concatenate([True_TrackDirection_xz, TTxz])
            True_TrackDirection_yz = np.concatenate([True_TrackDirection_yz, TTyz])
            Reco_Charge = np.concatenate([Reco_Charge, RC])
            True_Charge = np.concatenate([True_Charge, TC])
            True_Charge_stop = np.concatenate([True_Charge_stop, TCs])
            True_KE = np.concatenate([True_KE, TK])
            True_ZDirection = np.concatenate([True_ZDirection, TZ])
            Reco_ZDirection = np.concatenate([Reco_ZDirection, RZ])
            HitsPerTrack = np.concatenate([HitsPerTrack, HPT])
            
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
    if plot_Charge:
        boolean_Reco_Charge = (Reco_Charge != -9999.)
        Reco_Charge = Reco_Charge[boolean_Reco_Charge]
        boolean_True_Charge_stop = (True_Charge_stop != -9999.)
        True_Charge_stop = True_Charge_stop[boolean_True_Charge_stop]
        True_KE_stop = True_KE[boolean_True_Charge_stop & (True_KE != -9999.)]
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
        boolean_ZDirection = (True_ZDirection != -9999.)
        True_ZDirection = True_ZDirection[boolean_ZDirection]
        Reco_ZDirection = Reco_ZDirection[boolean_ZDirection]
    
    # total number of events after filtering
    if plot_End:
        print("------------------------------------------")
        print("     #events reconstruction: ", len(Reco_End_x))#, " + not stopped in TMS:", counter_leaving)
        print("------------------------------------------")
    #print("true (anti-)muons: ", sum(True_Muon_Track), "  vs. reconstructed particles: ", sum(Reco_Muon_Track))
    #print("correctly identified tracks: ", correct_tracks_reco)
    #print("total muons (outside of reco): ", count_muons)
    
    # subtract reconstruction from truth for all directions
    if plot_End:
        Diff_End_x = Reco_End_x - Primary_True_End_x
        Diff_End_y = Reco_End_y - Primary_True_End_y
        Diff_End_z = Reco_End_z - Primary_True_End_z
    
    if plot_Start:
        Diff_Start_x = Reco_Start_x - Primary_True_Start_x
        Diff_Start_y = Reco_Start_y - Primary_True_Start_y
        Diff_Start_z = Reco_Start_z - Primary_True_Start_z

    # create histograms for the differences
    if plot_End:
        Diff_End_x_hist, Diff_End_x_bins = np.histogram(Diff_End_x[~np.isnan(Diff_End_x)], bins = 100, range = (-2000, 2000))#range = (min([min(Diff_End_x), min(Diff_End_y), min(Diff_End_z)]), max([max(Diff_End_x), max(Diff_End_y), max(Diff_End_z)])))
        Diff_End_y_hist, Diff_End_y_bins = np.histogram(Diff_End_y[~np.isnan(Diff_End_y)], bins = 100, range = (-3000, 3000))#Diff_End_x_bins)
        Diff_End_z_hist, Diff_End_z_bins = np.histogram(Diff_End_z[~np.isnan(Diff_End_z)], bins = 100, range = (-3000, 2000))#Diff_End_x_bins)
   
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

        # fit distributions
        #popt_endY, pcov_endY = curve_fit(end_y_gauss, Diff_End_y_histX / 10, Diff_End_y_histY / len(Diff_End_y), p0 = [End_y_q_25 / 10, End_y_q_75 / 10, -130., 20., mean_End_y / 10, 20., 80., 0.8, 1.5, 0.1, 1.], sigma = np.sqrt(Diff_End_y_histY) / len(Diff_End_y), absolute_sigma = False)
        
        print("------------------------------------------")
        print("     End   x: ", mean_End_x, "±", std_End_x, "    y: ", mean_End_y, "±", std_End_y, "    z: ", mean_End_z, "±", std_End_z)
        print("     Quantile (25%): x: ", End_x_q_25, "     y: ", End_y_q_25, "     z: ", End_z_q_25)
        print("     Quantile (75%): x: ", End_x_q_75, "     y: ", End_y_q_75, "     z: ", End_z_q_75)
        #print("     Fit results:")
        #print("         y: ", popt_endY, " ± ", np.sqrt(np.diagonal(pcov_endY)))
        #print("           chi_square: ", chi_square(Diff_End_y_histY / len(Diff_End_y), end_y_gauss(Diff_End_y_histX / 10, *popt_endY), np.sqrt(Diff_End_y_histY) / len(Diff_End_y)), " / ", len(Diff_End_y_histY) / 2, " : ", chi_square(Diff_End_y_histY / len(Diff_End_y), end_y_gauss(Diff_End_y_histX / 10, *popt_endY), np.sqrt(Diff_End_y_histY) / len(Diff_End_y)) / (len(Diff_End_y_histX) / 2 - len(popt_endY)))
        print("------------------------------------------")

    if plot_Start:    
        Diff_Start_x_hist, Diff_Start_x_bins = np.histogram(Diff_Start_x[~np.isnan(Diff_Start_x)], bins = 100, range = (-500, 500))#range = (min([min(Diff_Start_x), min(Diff_Start_y), min(Diff_Start_z)]), max([max(Diff_Start_x), max(Diff_Start_y), max(Diff_Start_z)])))
        Diff_Start_y_hist, Diff_Start_y_bins = np.histogram(Diff_Start_y[~np.isnan(Diff_Start_y)], bins = 100, range = (-2000, 2000))#Diff_Start_x_bins)
        Diff_Start_z_hist, Diff_Start_z_bins = np.histogram(Diff_Start_z[~np.isnan(Diff_Start_z)], bins = 100, range = (-3000, 2000))#Diff_Start_x_bins)
    
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

        # fit distributions
        #popt_startX, pcov_startX = curve_fit(double_gauss, Diff_Start_x_histX / 10, Diff_Start_x_histY / len(Diff_Start_x), p0 = [Start_x_q_25 / 10, Start_x_q_75 / 10, mean_Start_x / 10, 0.5, 0.5], sigma = np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x), absolute_sigma = False)
        #popt_startY, pcov_startY = curve_fit(double_gauss, Diff_Start_y_histX / 10, Diff_Start_y_histY / len(Diff_Start_y), p0 = [Start_y_q_25 / 10, Start_y_q_75 / 10, mean_Start_y / 10, 0.5, 0.5], sigma = np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y), absolute_sigma = False)
    
        print("------------------------------------------")
        print("     Start   x: ", mean_Start_x, "±", std_Start_x, "    y: ", mean_Start_y, "±", std_Start_y, "    z: ", mean_Start_z, "±", std_Start_z)
        print("     Quantile (25%): x: ", Start_x_q_25, "     y: ", Start_y_q_25, "     z: ", Start_z_q_25)
        print("     Quantile (75%): x: ", Start_x_q_75, "     y: ", Start_y_q_75, "     z: ", Start_z_q_75)
        #print("     Fit results:")
        #print("         x: ", popt_startX, " ± ", np.sqrt(np.diagonal(pcov_startX)))
        #print("           chi_square: ", chi_square(Diff_Start_x_histY / len(Diff_Start_x), double_gauss(Diff_Start_x_histX / 10, *popt_startX), np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x)), " / ", len(Diff_Start_x_histX) / 2, " : ", chi_square(Diff_Start_x_histY / len(Diff_Start_x), double_gauss(Diff_Start_x_histX / 10, *popt_startX), np.sqrt(Diff_Start_x_histY) / len(Diff_Start_x)) / (len(Diff_Start_x_histX) / 2 - len(popt_startX)))
        #print("         y: ", popt_startY, " ± ", np.sqrt(np.diagonal(pcov_startY)))
        #print("           chi_square: ", chi_square(Diff_Start_y_histY / len(Diff_Start_y), double_gauss(Diff_Start_y_histX / 10, *popt_startY), np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y)), " / ", len(Diff_Start_y_histX) / 2, " : ", chi_square(Diff_Start_y_histY / len(Diff_Start_y), double_gauss(Diff_Start_y_histX / 10, *popt_startY), np.sqrt(Diff_Start_y_histY) / len(Diff_Start_y)) / (len(Diff_Start_y_histX) / 2 - len(popt_startY)))
        print("------------------------------------------")

    # plot
    if plot_Start:
        start_x_array = np.linspace(Diff_Start_x_histX[0] / 10, Diff_Start_x_histX[-1] / 10, len(Diff_Start_x_histX) * 2)
        mp.plot(Diff_Start_x_histX / 10, Diff_Start_x_histY / len(Diff_Start_x), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_x_histX / 10, 0, Diff_Start_x_histY / len(Diff_Start_x), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_x / 10, 0, max(Diff_Start_x_histY / len(Diff_Start_x)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_x / 10))
        #mp.axvspan(mean_Start_x / 10 - std_Start_x / 10, mean_Start_x / 10 + std_Start_x / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_Start_x / 10))
        #mp.vlines(Start_x_q_25 / 10, 0, max(Diff_Start_x_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (Start_x_q_25 / 10))
        #mp.vlines(Start_x_q_75 / 10, 0, max(Diff_Start_x_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (Start_x_q_75 / 10))
        mp.axvspan(Start_x_q_25 / 10, Start_x_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_x_q_25 / 10, Start_x_q_75 / 10, (Start_x_q_75 - Start_x_q_25) / 10))
        #mp.plot(start_x_array, double_gauss(start_x_array, *popt_startX), color = black_cbf, linestyle = '-', linewidth = 2, label = 'fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('Start x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_X_hist_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        start_y_array = np.linspace(Diff_Start_y_histX[0] / 10, Diff_Start_y_histX[-1] / 10, len(Diff_Start_y_histX) * 2)
        mp.plot(Diff_Start_y_histX / 10, Diff_Start_y_histY / len(Diff_Start_y), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_y_histX / 10, 0, Diff_Start_y_histY / len(Diff_Start_y), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_y / 10, 0, max(Diff_Start_y_histY / len(Diff_Start_y)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_y / 10))
        #mp.axvspan(mean_Start_y / 10 - std_Start_y / 10, mean_Start_y / 10 + std_Start_y / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_Start_y / 10))
        #mp.vlines(Start_y_q_25 / 10, 0, max(Diff_Start_y_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (Start_y_q_25 / 10))
        #mp.vlines(Start_y_q_75 / 10, 0, max(Diff_Start_y_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (Start_y_q_75 / 10))
        mp.axvspan(Start_y_q_25 / 10, Start_y_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_y_q_25 / 10, Start_y_q_75 / 10, (Start_y_q_75 - Start_y_q_25) / 10))
        #mp.plot(start_y_array, double_gauss(start_y_array, *popt_startY), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('Start y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_Y_hist_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_Start_z_histX / 10, Diff_Start_z_histY / len(Diff_Start_z), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_Start_z_histX / 10, 0, Diff_Start_z_histY / len(Diff_Start_z), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_Start_z / 10, 0, max(Diff_Start_z_histY / len(Diff_Start_z)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_Start_z / 10))
        #mp.axvspan(mean_Start_z / 10 - std_Start_z / 10, mean_Start_z / 10 + std_Start_z / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_Start_z / 10))
        #mp.vlines(Start_z_q_25 / 10, 0, max(Diff_Start_z_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (Start_z_q_25 / 10))
        #mp.vlines(Start_z_q_75 / 10, 0, max(Diff_Start_z_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (Start_z_q_75 / 10))
        mp.axvspan(Start_z_q_25 / 10, Start_z_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (Start_z_q_25 / 10, Start_z_q_75 / 10, (Start_z_q_75 - Start_z_q_25) / 10))
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('Start z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_Start_Z_hist_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()    
    
    if plot_End:
        mp.plot(Diff_End_x_histX / 10, Diff_End_x_histY / len(Diff_End_x), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_x_histX / 10, 0, Diff_End_x_histY / len(Diff_End_x), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_x / 10, 0, max(Diff_End_x_histY / len(Diff_End_x)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_x / 10))
        #mp.axvspan(mean_End_x / 10 - std_End_x / 10, mean_End_x / 10 + std_End_x / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_End_x / 10))
        #mp.vlines(End_x_q_25 / 10, 0, max(Diff_End_x_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (End_x_q_25 / 10))
        #mp.vlines(End_x_q_75 / 10, 0, max(Diff_End_x_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (End_x_q_75 / 10))
        mp.axvspan(End_x_q_25 / 10, End_x_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_x_q_25 / 10, End_x_q_75 / 10, (End_x_q_75 - End_x_q_25) / 10)) 
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('End x', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_X_hist_6_quant.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_X_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        end_y_array = np.linspace(Diff_End_y_histX[0] / 10, Diff_End_y_histX[-1] / 10, len(Diff_End_x_histX) * 2)
        mp.plot(Diff_End_y_histX / 10, Diff_End_y_histY / len(Diff_End_y), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')   # / 10
        mp.fill_between(Diff_End_y_histX / 10, 0, Diff_End_y_histY / len(Diff_End_y), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_y / 10, 0, max(Diff_End_y_histY / len(Diff_End_y)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_y / 10))
        #mp.axvspan(mean_End_y / 10 - std_End_y / 10, mean_End_y / 10 + std_End_y / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_End_y / 10))
        #mp.vlines(End_y_q_25 / 10, 0, max(Diff_End_y_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (End_y_q_25 / 10))
        #mp.vlines(End_y_q_75 / 10, 0, max(Diff_End_y_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (End_y_q_75 / 10))
        mp.axvspan(End_y_q_25 / 10, End_y_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_y_q_25 / 10, End_y_q_75 / 10, (End_y_q_75 - End_y_q_25) / 10))
        #mp.plot(end_y_array, end_y_gauss(end_y_array, *popt_endY), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('End y', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Y_hist_6_quant.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Y_hist.png' % out_dir, bbox_inches = 'tight')
        mp.close()
    
        mp.plot(Diff_End_z_histX / 10, Diff_End_z_histY / len(Diff_End_z), color = blue_cbf, linestyle = '--', linewidth = 2, label = 'Difference')
        mp.fill_between(Diff_End_z_histX / 10, 0, Diff_End_z_histY / len(Diff_End_z), alpha = 0.5, hatch = '//', color = blue_cbf)
        mp.vlines(mean_End_z / 10, 0, max(Diff_End_z_histY / len(Diff_End_z)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f cm)' % (mean_End_z / 10))
        #mp.axvspan(mean_End_z / 10 - std_End_z / 10, mean_End_z / 10 + std_End_z / 10, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f cm)' % (std_End_z / 10))
        #mp.vlines(End_z_q_25 / 10, 0, max(Diff_End_z_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f cm)' % (End_z_q_25 / 10))
        #mp.vlines(End_z_q_75 / 10, 0, max(Diff_End_z_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f cm)' % (End_z_q_75 / 10))
        mp.axvspan(End_z_q_25 / 10, End_z_q_75 / 10, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f cm)$\\to\\ Q_3$ (%2.2f cm)\n(= %2.2f cm)' % (End_z_q_25 / 10, End_z_q_75 / 10, (End_z_q_75 - End_z_q_25) / 10))
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth [cm]')
        mp.ylabel('#')
        mp.title('End z', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        if Contained:
            mp.savefig('%s/Difference_Contained_End_Z_hist_6_quant.png' % out_dir, bbox_inches = 'tight')
        if not Contained:
            mp.savefig('%s/Difference_End_Z_hist.png' % out_dir, bbox_inches = 'tight')
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
        
        # Prepare the KE arrays
        Muons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        AMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        FMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
        FAMuons_KE = np.ones(len(True_KE_stop), dtype = float) * -9999.
    
        # Counter for both positive, both negative, different reco/truth positive and negative
        true_muons = len(True_Charge_stop[True_Charge_stop == 13])
        true_antimuons = len(True_Charge_stop[True_Charge_stop == -13])
        reco_muons = len(Reco_Charge[Reco_Charge == 13])
        reco_antimuons = len(Reco_Charge[Reco_Charge == -13])
        counter_true_positive = 0   #truth and reco agree on muon
        counter_true_negative = 0   #truth and reco agree on antimuon
        counter_false_positive = 0  #truth and reco disagree, reco -> muon
        counter_false_negative = 0  #truth and reco disagree, reco -> antimuon
    
        True_Muons_KE = True_KE_stop[True_Charge_stop == 13]
        True_Antimuons_KE = True_KE_stop[True_Charge_stop == -13]
    
        # Compare sign of truth and reco    
        for i in range(len(True_Charge_stop)):
            if True_Charge[i] == Reco_Charge[i]:
                if True_Charge[i] > 0:
                    counter_true_positive += 1
                    Muons_KE[i] = True_KE_stop[i]
                else:
                    counter_true_negative += 1
                    AMuons_KE[i] = True_KE_stop[i]
            else:
                if Reco_Charge[i] > 0:
                    counter_false_positive += 1
                    FMuons_KE[i] = True_KE_stop[i]
                else:
                    counter_false_negative += 1
                    FAMuons_KE[i] = True_KE_stop[i]
        
        if (counter_true_positive + counter_true_negative + counter_false_positive + counter_false_negative) != len(True_Charge_stop):
            print("     Counters don't add up")
            print("     Length: ", len(True_Charge_stop))

        print("------------------------------------------")        
        print("     Charge ID numbers")
        print("         Not identified:  ", len(not_identified))
        print("         True muons:      ", counter_true_positive)
        print("         True antimuons:  ", counter_true_negative)
        print("         False muons:     ", counter_false_positive)
        print("         False antimuons: ", counter_false_negative, " here the not identified go into!")
        print("------------------------------------------")
        
        # Performance evaluation
        if counter_true_positive > 0:
            efficiency_muons = counter_true_positive / (counter_true_positive + (counter_false_negative - len(not_identified)))
            purity_muons = counter_true_positive / (counter_true_positive + counter_false_positive)
            print("  Muons (efficiency | purity):     ", efficiency_muons, " | ", purity_muons)
        if counter_true_negative > 0:
            efficiency_antimuons = counter_true_negative / (counter_true_negative + counter_false_positive)
            purity_antimuons = counter_true_negative / (counter_true_negative + (counter_false_negative - len(not_identified)))
            print("  Antimuons (efficiency | purity): ", efficiency_antimuons, " | ", purity_antimuons)

        accuracy_anti_muons = (counter_true_positive + counter_true_negative) / (counter_true_positive + counter_true_negative + counter_false_positive + (counter_false_negative - len(not_identified)))    
        print("  Accuracy (both): ", accuracy_anti_muons)

        # Energy dependent      
        Muons_KE = Muons_KE[Muons_KE != -9999.]
        AMuons_KE = AMuons_KE[AMuons_KE != -9999.]
        FMuons_KE = FMuons_KE[FMuons_KE != -9999.]
        FAMuons_KE = FAMuons_KE[FAMuons_KE != -9999.]
        
        muons_ke_hist, muons_ke_bins = np.histogram(Muons_KE, bins = 50, range = (0, 5000))
        amuons_ke_hist, amuons_ke_bins = np.histogram(AMuons_KE, bins = muons_ke_bins)
        fmuons_ke_hist, fmuons_ke_bins = np.histogram(FMuons_KE, bins = muons_ke_bins)
        famuons_ke_hist, famuons_ke_bins = np.histogram(FAMuons_KE, bins = muons_ke_bins)
    
        true_muons_ke_hist, muons_ke_bins = np.histogram(True_Muons_KE, bins = muons_ke_bins)
        true_antimuons_ke_hist, amuons_ke_bins = np.histogram(True_Antimuons_KE, bins = amuons_ke_bins)
        
        muons_2 = muons_ke_hist[19] + muons_ke_hist[20]
        amuons_2 = amuons_ke_hist[19] + amuons_ke_hist[20]
        fmuons_2 = fmuons_ke_hist[19] + fmuons_ke_hist[20]
        famuons_2 = famuons_ke_hist[19] + famuons_ke_hist[20]

        print("--------------------------------")
        print("     Charge at 2GeV")
        print("         true muons:     ", muons_2)
        print("         true amuons:    ", amuons_2)
        print("         false muons:    ", fmuons_2)
        print("         false amuons:   ", famuons_2)
        print("         eff. muon:      ", muons_2 / (muons_2 + famuons_2))
        print("         pur. muon:      ", muons_2 / (muons_2 + fmuons_2))
        print("         eff. amuon:     ", amuons_2 / (amuons_2 + fmuons_2))
        print("         pur. amuon:     ", amuons_2 / (amuons_2 + famuons_2))
        print("         accuracy:       ", (muons_2 + amuons_2) / (muons_2 + amuons_2 + fmuons_2 + famuons_2))

        muons_ke_hist_x, muons_ke_hist_y = histogram_arr_handle(muons_ke_hist, muons_ke_bins)
        amuons_ke_hist_x, amuons_ke_hist_y = histogram_arr_handle(amuons_ke_hist, amuons_ke_bins)
        fmuons_ke_hist_x, fmuons_ke_hist_y = histogram_arr_handle(fmuons_ke_hist, fmuons_ke_bins)
        famuons_ke_hist_x, famuons_ke_hist_y = histogram_arr_handle(famuons_ke_hist, famuons_ke_bins)
        true_muons_ke_hist_x, true_muons_ke_hist_y = histogram_arr_handle(true_muons_ke_hist, muons_ke_bins)
        true_antimuons_ke_hist_x, true_antimuons_ke_hist_y = histogram_arr_handle(true_antimuons_ke_hist, amuons_ke_bins)
    
        if counter_true_positive > 0:
            mp.plot(muons_ke_hist_x, muons_ke_hist_y / (true_muons_ke_hist_y + famuons_ke_hist_y), color = blue_cbf, linestyle = '-', linewidth = 2)
            mp.fill_between(muons_ke_hist_x, 0, muons_ke_hist_y / (true_muons_ke_hist_y + famuons_ke_hist_y), color = blue_cbf, alpha = 0.3, hatch = '//', label = '$\\mu$ correct')
        if counter_true_negative > 0:
            mp.plot(amuons_ke_hist_x, amuons_ke_hist_y / (true_antimuons_ke_hist_y + fmuons_ke_hist_y), color = orange_cbf, linestyle = '--', linewidth = 2)
            mp.fill_between(amuons_ke_hist_x, 0, amuons_ke_hist_y / (true_antimuons_ke_hist_y + fmuons_ke_hist_y), color = orange_cbf, alpha = 0.3, hatch = '\\\\', label = '$\\bar{\\mu}$ correct')
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Fraction')
        mp.legend(loc = 'lower center', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Charge_ID_KE.png' % out_dir, bbox_inches = 'tight')
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
    
        # calculate differences in angles
        xz_difference = True_TrackDirection_xz - Reco_TrackDirection_xz
        yz_difference = True_TrackDirection_yz - Reco_TrackDirection_yz
        
        xz_difference = xz_difference[~np.isnan(xz_difference)] #filter out potential nan's
        yz_difference = yz_difference[~np.isnan(yz_difference)] #filter out potential nan's
        
        # calculate mean and std of differences
        mean_xz_difference = np.mean(xz_difference)
        mean_yz_difference = np.mean(yz_difference)
        std_xz_difference = np.std(xz_difference)
        std_yz_difference = np.std(yz_difference)
        xz_q_25 = np.quantile(xz_difference, 0.25)
        yz_q_25 = np.quantile(yz_difference, 0.25)
        xz_q_75 = np.quantile(xz_difference, 0.75)
        yz_q_75 = np.quantile(yz_difference, 0.75)

        print("------------------------------------------")        
        print("     XZ angle difference: ", mean_xz_difference, " ± ", std_xz_difference)
        print("         25% quantile: ", xz_q_25, "  75% quantile: ", xz_q_75)
        print("     YZ angle difference: ", mean_yz_difference, " ± ", std_yz_difference)
        print("         25% quantile: ", yz_q_25, "  75% quantile: ", yz_q_75)
        #print("------------------------------------------")
        
        # calculate histograms of differences
        Diff_xz_hist, Diff_xz_bins = np.histogram(xz_difference, bins = 200, range = (-20, 20))
        Diff_yz_hist, Diff_yz_bins = np.histogram(yz_difference, bins = 400, range = (-75, 75))
       
        Diff_xz_histX, Diff_xz_histY = histogram_arr_handle(Diff_xz_hist, Diff_xz_bins)
        Diff_yz_histX, Diff_yz_histY = histogram_arr_handle(Diff_yz_hist, Diff_yz_bins)
        
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
        #mp.axvspan(mean_xz_difference - std_xz_difference, mean_xz_difference + std_xz_difference, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f°)' % std_xz_difference)
        #mp.vlines(xz_q_25, 0, max(Diff_xz_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f°)' % xz_q_25)
        #mp.vlines(xz_q_75, 0, max(Diff_xz_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f°)' % xz_q_75)
        mp.axvspan(xz_q_25, xz_q_75, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f°)$\\to\\ Q_3$ (%2.2f°)\n(= %2.2f°)' % (xz_q_25, xz_q_75, xz_q_75 - xz_q_25))
        #mp.plot(Diff_xz_histX, double_gauss_sep(Diff_xz_histX, *popt_xz), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('#')
        mp.title('xz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_XZ_angle_hist_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        mp.plot(Diff_yz_histX, Diff_yz_histY / len(yz_difference), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(Diff_yz_histX, 0, Diff_yz_histY / len(yz_difference), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(mean_yz_difference, 0, max(Diff_yz_histY / len(yz_difference)), color = red_cbf, linestyle = ':', linewidth = 3, label = 'mean (%2.2f°)' % mean_yz_difference)
        #mp.axvspan(mean_yz_difference - std_yz_difference, mean_yz_difference + std_yz_difference, alpha = 0.2, color = red_cbf, label = '1$\\sigma$ (%2.2f°)' % std_yz_difference)
        #mp.vlines(yz_q_25, 0, max(Diff_yz_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.2f°)' % yz_q_25)
        #mp.vlines(yz_q_75, 0, max(Diff_yz_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.2f°)' % yz_q_75)
        mp.axvspan(yz_q_25, yz_q_75, color = red_cbf, alpha = 0.2, label = '$Q_1$ (%2.2f°)$\\to\\ Q_3$ (%2.2f°)\n(= %2.2f°)' % (yz_q_25, yz_q_75, yz_q_75 - yz_q_25))
        #mp.plot(Diff_yz_histX, double_gauss_sep(Diff_yz_histX, *popt_yz), color = black_cbf, linestyle = '-', linewidth = 2, label = 'Fit')
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference truth – reconstruction [°]')
        mp.ylabel('#')
        mp.title('yz angle difference', fontsize = 'xx-large')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_YZ_angle_hist_quant.png' % out_dir, bbox_inches = 'tight')
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
        
        for i in range(len(bins) - 1):
            boolean_KE = (True_KE >= bins[i]) & (True_KE < bins[i+1])
            True_xz = True_TrackDirection_xz_KE[boolean_KE]
            True_yz = True_TrackDirection_yz_KE[boolean_KE]
            Reco_xz = Reco_TrackDirection_xz_KE[boolean_KE]
            Reco_yz = Reco_TrackDirection_yz_KE[boolean_KE]
            
            xz_diff = True_xz - Reco_xz
            yz_diff = True_yz - Reco_yz
            
            xz_mean_KE[i] = np.nanmean(xz_diff)
            xz_std_KE[i] = np.nanstd(xz_diff)
            yz_mean_KE[i] = np.nanmean(yz_diff)
            yz_std_KE[i] = np.nanstd(yz_diff)
            xz_q_25_KE[i] = np.nanquantile(xz_diff, 0.25)
            xz_q_75_KE[i] = np.nanquantile(xz_diff, 0.75)
            yz_q_25_KE[i] = np.nanquantile(yz_diff, 0.25)
            yz_q_75_KE[i] = np.nanquantile(yz_diff, 0.75)
    
        # plot
        mp.errorbar(bins, xz_mean_KE, yerr = xz_std_KE, linestyle = '--', linewidth = 1.5, color = red_cbf)
        mp.scatter(bins, xz_mean_KE, marker = '.', color = red_cbf, label = 'mean')
        #mp.fill_between(bins, xz_mean_KE - xz_std_KE, xz_mean_KE + xz_std_KE, color = red_cbf, alpha = 0.1, label = '1 $\\sigma$')
        #mp.plot(bins, xz_q_25_KE, color = black_cbf, linestyle = ':', linewidth = 2, label = '0.25 quantile')
        #mp.plot(bins, xz_q_75_KE, color = black_cbf, linestyle = '-', linewidth = 2, label = '0.75 quantile')
        mp.fill_between(bins, xz_q_25_KE, xz_q_75_KE, color = red_cbf, alpha = 0.1, label = '$Q_1\\ \\to\\ Q_3$')
        mp.hlines(0, 0, 5000, color = black_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('xz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_XZ_angle_True_KE_500MeV_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()
        
        mp.errorbar(bins, yz_mean_KE, yerr = yz_std_KE, linestyle = '--', linewidth = 1.5, color = blue_cbf)
        mp.scatter(bins, yz_mean_KE, marker = '.', color = blue_cbf, label = 'mean')
        #mp.fill_between(bins, yz_mean_KE - yz_std_KE, yz_mean_KE + yz_std_KE, color = blue_cbf, alpha = 0.1, label = '1 $\\sigma$')
        #mp.plot(bins, yz_q_25_KE, color = black_cbf, linestyle = ':', linewidth = 2, label = '0.25 quantile')
        #mp.plot(bins, yz_q_75_KE, color = black_cbf, linestyle = '-', linewidth = 2, label = '0.75 quantile')
        mp.fill_between(bins, yz_q_25_KE, yz_q_75_KE, color = blue_cbf, alpha = 0.1, label = '$Q_1\\ \\to\\ Q_3$')
        mp.hlines(0, 0, 5000, color = black_cbf, linewidth = 1)
        mp.xlabel('True Muon KE [MeV]')
        mp.ylabel('Difference truth - reco [°]')
        mp.title('yz')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.savefig('%s/Difference_YZ_angle_True_KE_500MeV_quant.png' % out_dir, bbox_inches = 'tight')
        mp.close()

        ### energy resolution estimation by cos(theta) = d_z / |d| (z-direction of angle)
        True_ZDirection = True_ZDirection[boolean_true_muon]
        Reco_ZDirection = Reco_ZDirection[boolean_true_muon]

        # calculate differences
        z_difference = -Reco_ZDirection - True_ZDirection

        z_difference = z_difference[~np.isnan(z_difference)]

        z_difference_mean = np.nanmean(z_difference)
        z_difference_std = np.nanstd(z_difference)
        z_difference_q_25 = np.nanquantile(z_difference, 0.25)
        z_difference_q_50 = np.nanquantile(z_difference, 0.5)
        z_difference_q_75 = np.nanquantile(z_difference, 0.75)

        print("------------------------------------------")        
        print("     Z difference: ", z_difference_mean, " ± ", z_difference_std)
        print("         25% quantile: ", z_difference_q_25, "  75% quantile: ", z_difference_q_75)
        print("         50% quantile: ", z_difference_q_50, "  50% quantile: ", z_difference_q_50)
        print("------------------------------------------")


        # calculate histogram
        z_difference_hist, z_difference_bins = np.histogram(z_difference, bins = 100, range = (-0.5, 0.3))

        z_difference_histX, z_difference_histY = histogram_arr_handle(z_difference_hist, z_difference_bins)

        # plot
        mp.plot(z_difference_histX, z_difference_histY / len(z_difference), color = blue_cbf, linestyle = '--', linewidth = 2)
        mp.fill_between(z_difference_histX, 0, z_difference_histY / len(z_difference), alpha = 0.6, hatch = '//', color = blue_cbf)
        mp.vlines(z_difference_mean, 0, max(z_difference_histY / len(z_difference)), color = black_cbf, linestyle = ':', linewidth = 2, label = 'mean (%2.4f)' % z_difference_mean)
        mp.vlines(z_difference_q_50, 0, max(z_difference_histY / len(z_difference)), color = red_cbf, linestyle = ':', linewidth = 3, label = '$Q_2 (%2.4f)$' % z_difference_q_50)
        mp.axvspan(z_difference_q_25, z_difference_q_75, alpha = 0.2, color = red_cbf, label = '$Q_1$ (%1.4f)$\\to\\ Q_3$ (%1.4f)\n(= %1.4f)' % (z_difference_q_25, z_difference_q_75, z_difference_q_75 - z_difference_q_25))
        #mp.vlines(z_difference_q_25, 0, max(z_difference_histY), color = black_cbf, linestyle = ':', linewidth = 3, label = '0.25 quant. (%2.4f)' % z_difference_q_25)
        #mp.vlines(z_difference_q_75, 0, max(z_difference_histY), color = black_cbf, linestyle = '-', linewidth = 3, label = '0.75 quant. (%2.4f)' % z_difference_q_75)
        mp.grid(True, linestyle = '--', alpha = 0.2)
        mp.xlabel('Difference reconstruction – truth')
        mp.ylabel('#')
        mp.title('$\\cos(\\theta)$')
        mp.legend(loc = 'best', fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
        mp.text(0.05, 0.90, r"$\mathdefault{\bf{DUNE}}$" + " Simulation", horizontalalignment = 'left', color = "black", fontsize = 18, transform = mp.gca().transAxes)
        mp.savefig('%s/Difference_ZDirection.png' % out_dir, bbox_inches = 'tight')
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
