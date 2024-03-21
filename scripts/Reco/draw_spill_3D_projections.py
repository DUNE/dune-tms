import ROOT
import numpy as np
import matplotlib.pyplot as mp
import os
import argparse
import cppyy.ll

# plotstyle
red_cbf = '#d55e00'
blue_cbf = '#0072b2'
orange_cbf = '#e69f00'
magenta_cbf = '#cc79a7'
black_cbf = '#000000'
mp.style.use('seaborn-poster')

mp.rc('axes', labelsize = 12)  # fontsize of the x and y labels
mp.rc('xtick', labelsize = 12) # fontsize of the tick labels
mp.rc('ytick', labelsize = 12) # fontsize of the tick labels

cos_3 = 0.99863
sin_3 = 0.05234
tan_87 = 19.08114
tan_3 = 0.05241

### Widths of hits
delta_x = 0.0177    # half a bar width
delta_y = 0.3383    # uncertainty from +/-3 degree tilted bars
delta_z = 0.02      # space of scintilattor with air gap

### Function for upper limit of tilted bar 'hit'
def upper_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'kVBar':  # assumption VBar is tilted in positive way and UBar then in negative
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y
        if x < r:
            return_value = tan_3 * x - tan_3 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
        elif x >= r:
            return_value = -tan_87 * x + tan_87 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
    elif orientation_bar == 'kUBar':
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y + sin_3 * delta_x + cos_3 * delta_y
        if x < r:
            return_value = tan_87 * x - tan_87 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
        elif x >= r:
            return_value = -tan_3 * x + tan_3 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value

### Function for lower limit of tilted bar 'hit'
def lower_limit(hit_x, hit_y, x, orientation_bar):
    if orientation_bar == 'kVBar':
        r = hit_x - cos_3 * delta_x + sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return_value = -tan_87 * x + tan_87 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
        elif x >= r:
            return_value = tan_3 * x - tan_3 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
    elif orientation_bar == 'kUBar':
        r = hit_x + cos_3 * delta_x - sin_3 * delta_y
        s = hit_y - sin_3 * delta_x - cos_3 * delta_y
        if x < r:
            return_value = -tan_3 * x + tan_3 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value
        elif x >= r:
            return_value = tan_87 * x - tan_87 * r + s
            if return_value >= -2.51: return -2.51
            elif return_value <= -5.71: return -5.71
            else: return return_value

### Function for hits to appear in correct size (according to bar size, or reconstructed hit area size)
def hit_size(hit_x, hit_y, orientation, hit_z):
    if orientation == 'xy':   # here it is the reconstructed hit area
        orientation_bar = check_orientation(int(hit_z))
        if orientation_bar == 'kXBar':
            size_array = np.zeros((2,2))
            size_array[0, 0] = hit_x + delta_x
            size_array[0, 1] = hit_x - delta_x
            if (hit_y + delta_x) >= -2.51:
                size_array[1, 0] = -2.51
            else: 
                size_array[1, 0] = hit_y + delta_x
            if (hit_y - delta_x) <= -5.71:
                size_array[1, 1] = -5.71
            else:
                size_array[1, 1] = hit_y - delta_x
            return = np.array(size_array[0]), size_array[1, 0], size_array[1, 1]
        else: 
            left_top = hit_x / 1000.0 - cos_3 * delta_x - sin_3 * delta_y
            right_bottom = hit_x / 1000.0 + cos_3 * delta_x + sin_3 * delta_y
            x_array = np.linspace(left_top, right_bottom, num = 50)
            return x_array, np.array([lower_limit(hit_x / 1000.0, hit_y / 1000.0, i, orientation_bar) for i in x_array]), np.array([upper_limit(hit_x / 1000.0, hit_y / 1000.0, i, orientation_bar) for i in x_array])
                            
    elif orientation == 'zy': # here it is the reconstructed hit area
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x / 1000.0 + delta_z
        size_array[0, 1] = hit_x / 1000.0 - delta_z
        orientation_bar = check_orientation(int(hit_z))
        if orientation_bar == 'kXBar':
            if (hit_y + delta_x) >= -2.51:
                size_array[1, 0] = -2.51
            else:
                size_array[1, 0] = hit_y / 1000.0 + delta_x
            if (hit_y - delta_x) <= -5.71:
                size_array[1, 1] = -5.71
            else:
                size_array[1, 1] = hit_y / 1000.0 - delta_x
        else:
            if (hit_y / 1000.0 + delta_y) >= -2.51:
                size_array[1, 0] = -2.51
            else:
                size_array[1, 0] = hit_y / 1000.0 + delta_y
            if (hit_y / 1000.0 - delta_y) <= -5.71:
                size_array[1, 1] = -5.71
            else:
                size_array[1, 1] = hit_y / 1000.0 - delta_y
       return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]        

    elif orientation == 'xz': # here it is only the bar size
        size_array = np.zeros((2,2))
        size_array[0, 0] = hit_x / 1000.0 + delta_z
        size_array[0, 1] = hit_x / 1000.0 - delta_z
        size_array[1, 0] = hit_y / 1000.0 + delta_x
        size_array[1, 1] = hit_y / 1000.0 - delta_x
        return np.array(size_array[0]), size_array[1, 0], size_array[1, 1]        

### Actual function that loops through the spills
def draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons = False):
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
            except KeyError:
                r.GetEntry(i)
                event = r
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event == None:
                r.GetEntry(i)
                event = r
            # Sync up the readout info if it's there. Note that it has one entry per spill, not timeslice
            if readout != None: readout.GetEntry(current_spill_number)
           
            ### Check if a track exists in the event/spill, otherwise skip it
            nTracks = event.nTracks
            if nTracks <= 0: continue
            
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
            x_y.text(3.6, -5, 'front view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            z_y.text(18.1, -4.8, 'side view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            x_z.text(18.1, -1, 'top view', rotation = 'vertical', fontsize = 12, fontweight = 'bold', color = orange_cbf)
            
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
            
            print("number of tracks: ", nTracks)
            nHits = np.frombuffer(event.nHits), dtype = np.uint8)
            nHits = np.array([nHits[i] for i in range(0, nTracks * 4, 4)])
            print("number of hits: ", nHits)
            
            TrackHitPos = np.frombuffer(event.TrackHitPos, dtype = np.float32)
            StartPos = np.frombuffer(event.StartPos, dtype = np.float32)
            EndPos = np.frombuffer(event.EndPos, dtype = np.float32)
            Direction = np.frombuffer(event.Direction, dtype = np.float32)

            for i in range(nTracks):

                ###  Now fill with the hits
                for hit in range(nHits[i]):
                    #print(TrackHitPos[hit])
                    hit_x = TrackHitPos[i*600 + hit*3 + 0]
                    hit_y = TrackHitPos[i*600 + hit*3 + 1]
                    hit_z = TrackHitPos[i*600 + hit*3 + 2]
                
                    #print(hit_x, hit_y, hit_z)

                    #temporary fix
                    if hit_y > -2000.0 or hit_z < 11000.: continue
                    if np.abs(hit_x) > 10000. or np.aby(hit_y) > 10000. or np.abs(hit_z) > 20000.: continue

                    color_cbf = red_cbf
                    if check_orientation(int(hit_z)) == 'kVBar':
                        color_cbf = blue_cbf
                    elif check_orientation(int(hit_z)) == 'kXBar':
                        color_cbf = black_cbf
                
                    if hit_x == -999.0 and hit_y == -999.0 and hit_z == -999.0: continue
                    x_z.fill_between(*hit_size(hit_z, hit_x, 'xz', hit_z), color = color_cbf)
                    z_y.fill_between(*hit_size(hit_z, hit_y, 'zy', hit_z), color = color_cbf)
                    x_y.fill_between(*hit_size(hit_x, hit_y, 'xy', hit_z), color = color_cbf, alpha = 0.5, linewidth = 0.5)
                 
                ### Track start
                #print(StartPos)
                
                #temporary fix
                if not (StartPos[i*3 + 1] > -2000.0 or StartPos[i*3 + 2] < 11000.0):
                    if not StartPos[i*3 + 1] == 0.0:
                        x_z.fill_between(*hit_size(StartPos[i*3 + 2], StartPos[i*3 + 0], 'xz', StartPos[i*3 + 2]), color = green_cbf)
                        z_y.fill_between(*hit_size(StartPos[i*3 + 2], StartPos[i*3 + 1], 'zy', StartPos[i*3 + 2]), color = green_cbf)
                        x_y.fill_between(*hit_size(StartPos[i*3 + 0], StartPos[i*3 + 1], 'xy', StartPos[i*3 + 2]), color = green_cbf, alpha = 0.5, linewidth = 0.5)
                 
                ### Track end
                #print(EndPos)

                #temporary fix
                if not (EndPos[i*3 + 1] > -2000.0 or EndPos[i*3 + 2] < 11000.0):
                    if not EndPos[i*3 + 1] == 0.0: 
                        x_z.fill_between(*hit_size(EndPos[i*3 + 2], EndPos[i*3 + 0], 'xz', EndPos[i*3 + 2]), color = green_cbf)
                        z_y.fill_between(*hit_size(EndPos[i*3 + 2], EndPos[i*3 + 1], 'zy', EndPos[i*3 + 2]), color = green_cbf)
                        x_y.fill_between(*hit_size(EndPos[i*3 + 0], EndPos[i*3 + 1], 'xy', EndPos[i*3 + 2]), color = green_cbf, alpha = 0.5, linewidth = 0.5)
                 
                ### Track direction
                #print(Direction)

                #temporary fix
                if not (StartPos[i*3 + 1] > -2000.0 or StartPos[i*3 + 2] < 11000.0 or EndPos[i*3 + 1] > -2000.0 or EndPos[i*3 + 2] < 11000.0):
                    if not StartPos[i*3 + 1] == 0.0:
                        x_z.plot([StartPos[i*3 + 2] / 1000.0, (StartPos[i*3 + 2] + Direction[i*3 + 2]) / 1000.0], [StartPos[i*3 + 0] / 1000.0, (StartPos[i*3 + 0] + Direction[i*3 + 0]) / 1000.0], color = green_cbf, linewidth = 1.5, linestyle = '--')
                        z_y.plot([StartPos[i*3 + 2] / 1000.0, (StartPos[i*3 + 2] + Direction[i*3 + 2]) / 1000.0], [StartPos[i*3 + 1] / 1000.0, (StartPos[i*3 + 1] + Direction[i*3 + 1]) / 1000.0], color = green_cbf, linewidth = 1.5, linestyle = '--')
                        x_y.plot([StartPos[i*3 + 0] / 1000.0, (StartPos[i*3 + 0] + Direction[i*3 + 0]) / 1000.0], [StartPos[i*3 + 1] / 1000.0, (StartPos[i*3 + 1] + Direction[i*3 + 1]) / 1000.0], color = green_cbf, linewidth = 1.5, linestyle = '--')

            output_filename = os.path.join(out_dir, f"{name}_{current_spill_number:03d}")
            mp.savefig(output_filename + ".png", bbox_inches = 'tight')
            mp.close()
             
    return

### This is for plotting the hits according to their different orientations
def check_orientation(hit_z):
    return layer_dict["%s" % hit_z]

### Dictionary that after calculate_layers contains for each z-coordinate the orientation str
first_z = 11368
layer_dict = { "%s" % first_z : "kUBar" }

def calculate_layers():
    thin_layers = 39
    thick_layers = 61
    # Calculate the z position for each layer for the thin section
    for i in range(thin_layers):
        hit_z = first_z + i * 55
        if ((hit_z - first_z) / 55) % 2 == 0: # even layers
            layer_dict.update({ "%s" % hit_z : "kUBar" })
        elif ((hit_z - first_z) / 55) % 2 == 1: # odd layers
            layer_dict.update({ "%s" % hit_z : "kVBar" })
    
    # Calculate the z position for each layer for the thick section
    start_thick = first_z + thin_layers * 55
    for i in range(thick_layers):
        hit_z = start_thick + i * 80
        if ((hit_z - start_thick) / 80) % 2 == 0: # even layers
            layer_dict.update({ "%s" % hit_z : "kVBar" })
        elif ((hit_z - start_thick) / 80) % 2 == 1: # odd layers
            layer_dict.update({ "%s" % hit_z : "kUBar" })

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Draws spills.")
    parser.add_argument('--outdir', "-o", type = str, help = "The output dir. Will be made if it doesn't exist. Default = spills/", default = "spills")
    parser.add_argument('--name', "-n", type = str, help = "The name of the output files. Will be <name>_<spill>_<slice>.png. Default = spill", default = "spill")
    parser.add_argument('--input_filename', "-f", type = str, help = "The file with the events to draw.")
    parser.add_argument('--spillnum', "-s", type = int, help = "The spill to draw. -1 for all", default = -1)
    parser.add_argument('--timeslice', "-t", type = int, help = "The time slice to draw. -1 for all", default = -1)
    parser.add_argument('--readout_filename', "-r", type = str, help = "(optional) A file with the raw readout.", default = "")
    parser.add_argument('--only_true_tms_muons', help = "Only draw events with true muons inside the TMS", action = argparse.BooleanOptionalAction)
    
    args = parser.parse_args()

    calculate_layers()
    
    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    readout_filename =  args.readout_filename
    only_true_tms_muons = args.only_true_tms_muons
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons)

