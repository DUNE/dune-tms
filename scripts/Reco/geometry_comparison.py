import numpy as np
import matplotlib.pyplot as mp

red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
green_cbf = '#009e73'
mp.style.use('seaborn-poster')

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


### mean and std arrays for directions for all geometries
same3_start_mean = np.array([0.47, -1.61, 1.50])
same3_end_mean = np.array([38.25, -1.72, -3.69])
same3_start_std = np.array([9.69, 53.45, 55.91])
same3_end_std = np.array([66.50, 243.23, 19.39])

same4_start_mean = np.array([1.65, 2.90, 29.68])
same4_start_std = np.array([10.13, 63.62, 59.62])
same4_end_mean = np.array([30.80, 480.49, 8.70])
same4_end_std = np.array([90.22, 765.52, 33.63])

nsame4_start_mean = np.array([0.26, -1.33, -0.28])
nsame4_start_std = np.array([15.86, 88.58, 58.52])
nsame4_end_mean = np.array([84.75, -11.22, -3.46])
nsame4_end_std = np.array([182.73, 267.52, 20.62])

same5_start_mean = np.array([-0.29, 0.16, 45.49])
same5_start_std = np.array([19.24, 92.14, 78.75])
same5_end_mean = np.array([75.33, -68.41, -4.33])
same5_end_std = np.array([156.13, 673.49, 22.74])

same6_start_mean = np.array([-0.51, -32.21, 24.05])
same6_start_std = np.array([5.87, 106.94, 58.38])
same6_end_mean = np.array([55.07, 297.38, -6.63])
same6_end_std = np.array([65-33, 519.39, 18.96])

nsame6_start_mean = np.array([-0.33, -11.58, -1.15])
nsame6_start_std = np.array([16.21, 108.03, 56.53])
nsame6_end_mean = np.array([36.27, 75.38, -3.85])
nsame6_end_std = np.array([71.94, 442.97, 18.41])

same7_start_mean = np.array([-0.42, 6.86, 44.62])
same7_start_std = np.array([17.76, 108.29, 78.98])
same7_end_mean = np.array([65.07, -24.76, -4.10])
same7_end_std = np.array([147.83, 869.13, 20.74])

noX_start_mean = np.array([-0.51, -3.23, 28.70])
noX_start_std = np.array([7.81, 152.34, 63.89])
noX_end_mean = np.array([16.00, 34.42, -4.34])
noX_end_std = np.array([16.57, 413.97, 18.17])

### re-arranging by direction
x_start_mean = np.array([same3_start_mean[0], nsame4_start_mean[0], same5_start_mean[0], nsame6_start_mean[0], same7_start_mean[0], noX_start_mean[0]])
y_start_mean = np.array([same3_start_mean[1], nsame4_start_mean[1], same5_start_mean[1], nsame6_start_mean[1], same7_start_mean[1], noX_start_mean[1]])
z_start_mean = np.array([same3_start_mean[2], nsame4_start_mean[2], same5_start_mean[2], nsame6_start_mean[2], same7_start_mean[2], noX_start_mean[2]])
x_start_std = np.array([same3_start_std[0], nsame4_start_std[0], same5_start_std[0], nsame6_start_std[0], same7_start_std[0], noX_start_std[0]])
y_start_std = np.array([same3_start_std[1], nsame4_start_std[1], same5_start_std[1], nsame6_start_std[1], same7_start_std[1], noX_start_std[1]])
z_start_std = np.array([same3_start_std[2], nsame4_start_std[2], same5_start_std[2], nsame6_start_std[2], same7_start_std[2], noX_start_std[2]])
x_end_mean = np.array([same3_end_mean[0], nsame4_end_mean[0], same5_end_mean[0], nsame6_end_mean[0], same7_end_mean[0], noX_end_mean[0]])
y_end_mean = np.array([same3_end_mean[1], nsame4_end_mean[1], same5_end_mean[1], nsame6_end_mean[1], same7_end_mean[1], noX_end_mean[1]])
z_end_mean = np.array([same3_end_mean[2], nsame4_end_mean[2], same5_end_mean[2], nsame6_end_mean[2], same7_end_mean[2], noX_end_mean[2]])
x_end_std = np.array([same3_end_std[0], nsame4_end_std[0], same5_end_std[0], nsame6_end_std[0], same7_end_std[0], noX_end_std[0]])
y_end_std = np.array([same3_end_std[1], nsame4_end_std[1], same5_end_std[1], nsame6_end_std[1], same7_end_std[1], noX_end_std[1]])
z_end_std = np.array([same3_end_std[2], nsame4_end_std[2], same5_end_std[2], nsame6_end_std[2], same7_end_std[2], noX_end_std[2]])

x_array = np.array([3, 4, 5, 6, 7, 8])

extra_x_start_mean = np.array([nsame4_start_mean[0], same6_start_mean[0]])
extra_y_start_mean = np.array([nsame4_start_mean[1], same6_start_mean[1]])
extra_z_start_mean = np.array([nsame4_start_mean[2], same6_start_mean[2]])
extra_x_start_std = np.array([nsame4_start_std[0], same6_start_std[0]])
extra_y_start_std = np.array([nsame4_start_std[1], same6_start_std[1]])
extra_z_start_std = np.array([nsame4_start_std[2], same6_start_std[2]])
extra_x_end_mean = np.array([nsame4_end_mean[0], same6_end_mean[0]])
extra_y_end_mean = np.array([nsame4_end_mean[1], same6_end_mean[1]])
extra_z_end_mean = np.array([nsame4_end_mean[2], same6_end_mean[2]])
extra_x_end_std = np.array([nsame4_end_std[0], same6_end_std[0]])
extra_y_end_std = np.array([nsame4_end_std[1], same6_end_std[1]])
extra_z_end_std = np.array([nsame4_end_std[2], same6_end_std[2]])

extra_x = np.array([4, 6])

### plot for directions
mp.errorbar(x_array - 0.1, z_start_mean, yerr = z_start_std, linestyle = '', linewidth = 1.5, color = black_cbf)
mp.errorbar(x_array, x_start_mean, yerr = x_start_std, linestyle = '', linewidth = 1.5, color = red_cbf)
mp.errorbar(x_array + 0.1, y_start_mean, yerr = y_start_std, linestyle = '', linewidth = 1.5, color = blue_cbf)
#mp.errorbar(extra_x, extra_x_start_mean, yerr = extra_x_start_std, linestyle = '', linewidth = 1.5, color = red_cbf)
#mp.errorbar(extra_x + 0.1, extra_y_start_mean, yerr = extra_y_start_std, linestyle = '', linewidth = 1.5, color = blue_cbf)
#mp.errorbar(extra_x - 0.1, extra_z_start_mean, yerr = extra_z_start_std, linestyle = '', linewidth = 1.5, color = black_cbf)
mp.scatter(x_array - 0.1, z_start_mean, marker = '1', color = black_cbf, label = 'z (-0.1)')
mp.scatter(x_array, x_start_mean, marker = '3', color = red_cbf, label = 'x')
mp.scatter(x_array + 0.1, y_start_mean, marker = '2', color = blue_cbf, label = 'y (+0.1)')
#mp.scatter(extra_x - 0.1, extra_z_start_mean, marker = '1', color = black_cbf)
#mp.scatter(extra_x, extra_x_start_mean, marker = '3', color = red_cbf)
#mp.scatter(extra_x + 0.1, extra_y_start_mean, marker = '2', color = blue_cbf)
mp.fill_between(x_array - 0.1, z_start_mean - z_start_std, z_start_mean + z_start_std, alpha = 0.1, color = black_cbf, label = 'z $1\\sigma$')
mp.fill_between(x_array, x_start_mean - x_start_std, x_start_mean + x_start_std, alpha = 0.2, color = red_cbf, label = 'x $1\\sigma$')
mp.fill_between(x_array + 0.1, y_start_mean - y_start_std, y_start_mean + y_start_std, alpha = 0.2, color = blue_cbf, label = 'y $1\\sigma$')
mp.xticks([3, 4, 5, 6, 7, 8], ['3layer', '4layer', '5layer', '6layer', '7layer', 'no x'])
mp.xlabel('Every Xth layer is x layer')
mp.ylabel('Start point resolution')
mp.title('Start')
mp.legend(loc = 'lower left', ncol = 2, fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
mp.grid(True, linestyle = '--', alpha = 0.2)
mp.savefig('geometry_comparison/start_point_resolution.png', bbox_inches = 'tight')
mp.close()

mp.errorbar(x_array - 0.1, z_end_mean, yerr = z_end_std, linestyle = '', linewidth = 1.5, color = black_cbf)
mp.errorbar(x_array, x_end_mean, yerr = x_end_std, linestyle = '', linewidth = 1.5, color = red_cbf)
mp.errorbar(x_array + 0.1, y_end_mean, yerr = y_end_std, linestyle = '', linewidth = 1.5, color = blue_cbf)
#mp.errorbar(extra_x - 0.1, extra_z_end_mean, yerr = extra_z_end_std, linestyle = '', linewidth = 1.5, color = black_cbf)
#mp.errorbar(extra_x, extra_x_end_mean, yerr = extra_x_end_std, linestyle = '', linewidth = 1.5, color = red_cbf)
#mp.errorbar(extra_x + 0.1, extra_y_end_mean, yerr = extra_y_end_std, linestyle = '', linewidth = 1.5, color = blue_cbf)
mp.scatter(x_array - 0.1, z_end_mean, marker = '1', color = black_cbf, label = 'z (-0.1)')
mp.scatter(x_array, x_end_mean, marker = '3', color = red_cbf, label = 'x')
mp.scatter(x_array + 0.1, y_end_mean, marker = '2', color = blue_cbf, label = 'y (+0.1)')
#mp.scatter(extra_x - 0.1, extra_z_end_mean, marker = '1', color = black_cbf)
#mp.scatter(extra_x, extra_x_end_mean, marker = '3', color = red_cbf)
#mp.scatter(extra_x + 0.1, extra_y_end_mean, marker = '2', color = blue_cbf)
mp.fill_between(x_array - 0.1, z_end_mean - z_end_std, z_end_mean + z_end_std, alpha = 0.1, color = black_cbf, label = 'z $1\\sigma$')
mp.fill_between(x_array, x_end_mean - x_end_std, x_end_mean + x_end_std, alpha = 0.2, color = red_cbf, label = 'x $1\\sigma$')
mp.fill_between(x_array + 0.1, y_end_mean - y_end_std, y_end_mean + y_end_std, alpha = 0.2, color = blue_cbf, label = 'y $1\\sigma$')
mp.xticks([3, 4, 5, 6, 7, 8], ['3layer', '4layer', '5layer', '6layer', '7layer', 'no x'])
mp.xlabel('Every Xth layer is x layer')
mp.ylabel('End point resolution')
mp.title('End')
mp.legend(loc = 'lower left', ncol = 2, fontsize = 'xx-large', markerscale = 1.0, columnspacing = 0.5, handlelength = 0.8)
mp.grid(True, linestyle = '--', alpha = 0.2)
mp.savefig('geometry_comparison/end_point_resolution.png', bbox_inches = 'tight')
mp.close()
