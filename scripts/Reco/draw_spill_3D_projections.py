import ROOT
import numpy as np
import matplotlib.pyplot as mp

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
mp.style.use('seaborn-poster')

mp.rc('axes', labelsize = 12)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 12) # fontsize of the tick labels
mp.rc('ytick', labelsize = 12) # fontsize of the tick labels

#def draw_spill(

#TODO path to file

#TODO create loop over all spills/events
root_file = ROOT.TFile.Open("path/to/file.root");
for spill in root_file.Reco_Tree: #or TChain???


#TODO access ROOT Tree with hit information


cos_3 = 0.99863
sin_3 = 0.05234
tan_87 = 19.08114
tan_3 = 0.05241

### Widths of hits
delta_x = 0.0177    # half a bar width
delta_y = 0.3383    # uncertainty from +/-3 degree tilted bars
delta_z = 0.02      # space of scintilattor with air gap

### Function for lower limit of tilted bar 'hit'
def lower_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'plus':   # dummy one here TODO exchange for correct one
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return -tan_87 * x + tan_87 * r + s
        elif x >= r:
            return tan_3 * x - tan_3 * r + s
    elif orientation_bar == 'minus':
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return -tan_3 * x + tan_3 * r + s
        elif x >= r:
            return tan_87 * x - tan_87 * r + s

### Function for upper limit of tilted bar 'hit'
def upper_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'plus':
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y
        if x < r:
            return tan_3 * x - tan_3 * r + s
        elif x >= r:
            return -tan_87 * x + tan_87 * r + s
    elif orientation_bar == 'minus':
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y
        if x < r:
            return tan_87 * x - tan_87 * r + s
        elif x >= r:
            return -tan_3 * x + tan_3 * r + s

### Function for hits to appear in size
def hit_size(hit_x, hit_y, orientation, orientation_bar):
    if orientation == 'xy':
        left_top = hit_x - cos_3 * delta_x - sin_3 * delta_y
        right_bottom = hit_x + cos_3 * delta_x + sin_3 * delta_y
        x_array = np.linspace(left_top, right_bottom, num = 50)
        return x_array, np.array([lower_limit(hit_x, hit_y, i, orientation_bar) for i in x_array]), np.array([upper_limit(hit_x, hit_y, i, orientation_bar) for i in x_array])
                       
    elif orientation == 'zy':
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x + delta_z
        size_array[0, 1] = hit_x - delta_z
        size_array[1, 0] = hit_y + delta_y
        size_array[1, 1] = hit_y - delta_y
        return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]        
    elif orientation == 'xz':
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x + delta_z
        size_array[0, 1] = hit_x - delta_z
        size_array[1, 0] = hit_y + delta_x
        size_array[1, 1] = hit_y - delta_x
        return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]      



### Create subplots
fig = mp.figure(constrained_layout = False)
gs = fig.add_gridspec(2, 2, hspace = 0.25, wspace = 0.15)
x_y = fig.add_subplot(gs[0, 0])
z_y = fig.add_subplot(gs[1, 0])
x_z = fig.add_subplot(gs[0:, 1:])

### Set labels and ticks
x_y.set(xlabel = 'x [m]', ylabel = 'y [m]', xticks = [-4, -3, -2, -1, 0, 1, 2, 3, 4], yticks = [-5, -4, -3])
z_y.set(xlabel = 'z [m]', ylabel = 'y [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-5, -4, -3])
x_z.set(xlabel = 'z [m]', ylabel = 'x [m]', xticks = [11, 12, 13, 14, 15, 16, 17, 18], yticks = [-3, -2, -1, 0, 1, 2, 3])
x_y.text(2.2, -6, 'front view', fontsize = 12, fontweight = 'bold', color = orange_cbf)
z_y.text(16.8, -6, 'side view', fontsize = 12, fontweight = 'bold', color = orange_cbf)
x_z.text(16.9, -3.8, 'top view', fontsize = 12, fontweight = 'bold', color = orange_cbf)

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
x_z.hlines(-3.5, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
x_z.hlines(3.5, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
x_z.hlines(-1.75, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
x_z.hlines(1.75, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
x_z.vlines(11, -3.5, 3.5, color = orange_cbf, linewidth = 1, linestyle = ':')
x_z.vlines(18, -3.5, 3.5, color = orange_cbf, linewidth = 1, linestyle = ':')

z_y.hlines(-5.71, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
z_y.hlines(-2.51, 11, 18, color = orange_cbf, linewidth = 1, linestyle = ':')
z_y.vlines(11, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':')
z_y.vlines(18, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':')

x_y.hlines(-5.71, -3.5, 3.5, color = orange_cbf, linewidth = 1, linestyle = ':')
x_y.hlines(-2.51, -3.5, 3.5, color = orange_cbf, linewidth = 1, linestyle = ':')
x_y.vlines(-3.5, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':')  #TODO this is simplified without tilt of modules
x_y.vlines(3.5, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':')   #TODO this is simplified without tilt of modules
x_y.vlines(-1.75, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':') #TODO this is simplified without tilt of modules
x_y.vlines(1.75, -2.51, -5.71, color = orange_cbf, linewidth = 1, linestyle = ':')  #TODO this is simplified without tilt of modules

### Test hit
x_z.fill_between(*hit_size(15, -2, 'xz', 'plus'), color = blue_cbf)
z_y.fill_between(*hit_size(15, -3, 'zy', 'plus'), color = blue_cbf)
x_y.fill_between(*hit_size(-2, -3, 'xy', 'plus'), color = blue_cbf, alpha = 0.5, linewidth = 0.5)
x_y.fill_between(*hit_size(-2, -3, 'xy', 'minus'), color = red_cbf, alpha = 0.5, linewidth = 0.5)

mp.savefig('event_layout_test.png', bbox_inches = 'tight')
mp.close()

#TODO create subplots with all 2D projections and time?
#TODO fill subplots with hit information necessary
  #TODO make sure size of bars are reflected
#TODO fill subplots with track information necessary
#TODO save subplot in numeric variable name file




#TODO make sure all arguments from old draw spill file also work here?
