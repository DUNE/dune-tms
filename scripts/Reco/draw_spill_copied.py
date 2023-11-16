import glob
import argparse
import math
import os
import collections
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
#ROOT.TImage.AddDirectory(False)
ROOT.TH2.AddDirectory(False)


ROOT.gStyle.SetPalette(52) # greyscale
ROOT.TColor.InvertPalette();

canvas = ROOT.TCanvas("c", "c", 1200, 800)
#canvas.Divide(1, 2)
pad1 = ROOT.TPad("p1", "p1", 0.0,0.2,1,1.)
pad1.SetRightMargin(0.13) # Needed to see z axis label
pad1.SetRightMargin(0.025)
pad2 = ROOT.TPad("p2", "p2", 0.0,0.,0.6,0.21)
pad2.SetBottomMargin(0.15)
pad2.SetRightMargin(0)
pad2.SetTopMargin(0.05)
pad3 = ROOT.TPad("p2", "p2", 0.6,0.,1,0.21)
pad3.SetBottomMargin(0.15)
pad3.SetLeftMargin(0)
pad3.SetTopMargin(0.05)
pad4 = ROOT.TPad("p4", "p4", 0.0,0.,1,0.21)
pad4.SetBottomMargin(0.35)
pad4.SetRightMargin(0.025)
pad4.SetTopMargin(0.05)

def draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons = False):
    if not os.path.exists(input_filename): raise ValueError(f"Cannot find input_filename {input_filename}")
    if readout_filename != "" and not os.path.exists(readout_filename): raise ValueError(f"Cannot find readout_filename {readout_filename}")
    if spill_number < -1: raise ValueError(f"Got spill_number = {spill_number}")
    if time_slice < -1: raise ValueError(f"Got time_slice = {time_slice}")
    
    # Now configure some stuff
    use_readout = True
    if readout_filename == "": use_readout = False
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        if not os.path.exists(out_dir):
            raise ValueError(f"Could not make out_dir {out_dir}")
            
    bound_z_min = 11.0
    bound_z_max = 18.5000
    bound_x_min = -4.0
    bound_x_max = 4.0
    
    font_size = 0.15
    
    markers = []
    markers_not_in_slice = []
    lines = []
    text_to_draw = []
    
    c = ROOT.TChain("Line_Candidates")
    c.Add(input_filename)
    print("N entries:", c.GetEntries())
    assert c.GetEntries() > 0, "Didn't get any entries, are you sure the input_filename is right?\n" + input_filename
    
    truth = ROOT.TChain("Truth_Info")
    truth.Add(input_filename)
    assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info, are you sure the input_filename is right?\n" + input_filename
    
    readout = None
    if use_readout:
        readout = ROOT.TChain("TMS")
        readout.Add(readout_filename)
        assert readout.GetEntries() > 0, "Didn't get any entries in TMS, are you sure the readout_filename is right?\n" + readout_filename
        
        
    max_n_spills = 10000 # TODO add some meta info to output file with n spill info for file
    
    haxis = ROOT.TH2D("hits_in_event", ";Position Z (m);Position X (m);Total Energy (MeV)", 100, bound_z_min, bound_z_max, 100, bound_x_min, bound_x_max)
    
    simplify_tracks = False
    
    spill_number_cache = dict()
    n_events = c.GetEntries()
        
    # First loop through all the slices and draw one overall spill 
    for current_spill_number in range(max_n_spills):
        htime = ROOT.TH1D(f"htime_{current_spill_number}", ";Time (ns);Total Hit E (MeV)", 500, -200, 12000)
        htime.GetYaxis().SetTitleOffset(0.2)
        htime.GetXaxis().SetTitleSize(font_size)
        htime.GetXaxis().SetLabelSize(font_size)
        htime.GetYaxis().SetTitleSize(font_size*0.8)
        htime.GetYaxis().SetLabelSize(font_size*0.8)
        htime.SetFillColor(ROOT.kBlack)
        htime.SetLineColor(ROOT.kBlack)
        
        total_energy = 0
        time_high = float("-inf")
        time_low = float("inf")
        
        markers = []
        text_to_draw = []
        for i in range(n_events):
            try:
                spill_number = spill_number_cache[i]
                event = None
            except KeyError:
                c.GetEntry(i)
                event = c
                spill_number = event.SpillNo
                spill_number_cache[i] = spill_number
            if spill_number < current_spill_number: continue
            if spill_number > current_spill_number: break
            if event == None:
                c.GetEntry(i)
                event = c
            # Sync up the truth info if it's there
            if truth != None: truth.GetEntry(i)
            # Sync up the readout info if it's there. Note that it has one entry per spill, not timeslice
            if readout != None: readout.GetEntry(current_spill_number)
            
            if only_true_tms_muons:
                assert truth != None, "Truth shouldn't be None inside only_true_tms_muons"
                def inside_tms(x, y, z):
                    is_inside = True
                    if not -4000 < x < 4000: is_inside = False
                    if not -1548-4000 < y < -1548+4000: is_inside = False
                    if not 11000 < z < 18000: is_inside = False
                    return is_inside
                mx = truth.Muon_Vertex[0]
                my = truth.Muon_Vertex[1]
                mz = truth.Muon_Vertex[2]
                mdx = truth.Muon_Death[0]
                mdy = truth.Muon_Death[1]
                mdz = truth.Muon_Death[2]
                start_inside_tms = inside_tms(mx, my, mz)
                end_inside_tms = inside_tms(mdx, mdy, mdz)
                # Skip any events that don't start and stop inside the TMS
                if not start_inside_tms: continue
                if not end_inside_tms: continue
                print(f"Muon Start XYZ: ({mx:0.2f}, {my:0.2f}, {mz:0.2f})\tMuon End XYZ: ({mdx:0.2f}, {mdy:0.2f}, {mdz:0.2f})\tstart_inside_tms: {start_inside_tms}\tend_inside_tms: {end_inside_tms},\tPDG: {truth.LeptonPDG}")
            
            print(f"Event {i} has {event.nHits} hits, {event.nClustersOne} clusters(one) | {event.nClustersOther} clusters(other), and {event.nLinesOne} lines(one) | {event.nLinesOther} lines(other).")
            
            
            for hit in range(event.nHits):
                if event.RecoHitEnergy[hit] < 0: continue
                # Seems like pyROOT maps RecoHitPos[nHits][4] onto a RecoHitPos[nHits*4] array
                reco_hit_x = event.RecoHitPos[hit*4 + 0] / 1000.0
                reco_hit_y = event.RecoHitPos[hit*4 + 1] / 1000.0
                reco_hit_z = event.RecoHitPos[hit*4 + 2] / 1000.0
                reco_hit_time = event.RecoHitPos[hit*4 + 3]
                reco_hit_energy = event.RecoHitEnergy[hit]
                
                time_high = max(reco_hit_time, time_high)
                time_low = min(reco_hit_time, time_low)
                
                
                e = int(min(255, 255 * reco_hit_energy / 10.0))
                t = int(min(255, 255 * reco_hit_time / 10000.0))
                
                htime.Fill(reco_hit_time, reco_hit_energy)
                x = reco_hit_z
                y = reco_hit_x
                
                marker = ROOT.TMarker(x, y, 21)
                color = ROOT.TColor.GetColor(e, 32, 32)
                marker.SetMarkerColor(color)
                markers.append(marker)
                marker = 0
                

                if reco_hit_energy > 0:
                    total_energy += reco_hit_energy
                elif reco_hit_energy < 0:
                    print(f"Found unexpected reco_hit_energy < 0: {reco_hit_energy}")
                
                for cluster in range(min(25, event.nClustersOne)):
                    cluster_energy = event.ClusterEnergyOne[cluster]
                    cluster_time = event.ClusterTimeOne[cluster]
                    cluster_pos_z = event.ClusterPosMeanOne[cluster*2 + 0] / 1000.0
                    cluster_pos_x = event.ClusterPosMeanOne[cluster*2 + 1] / 1000.0
                    cluster_pos_z_std_dev = event.ClusterPosStdDevOne[cluster*2 + 0] / 1000.0
                    cluster_pos_x_std_dev = event.ClusterPosStdDevOne[cluster*2 + 1] / 1000.0
                    #cluster_total_std_dev = math.sqrt(cluster_pos_z_std_dev**2 + cluster_pos_x_std_dev**2)
                    cluster_total_std_dev = max(cluster_pos_z_std_dev, cluster_pos_x_std_dev)
                    e = int(min(255, 255 * cluster_energy / 10.0))
                    t = int(min(255, 255 * cluster_time / 10000.0))
                    color = ROOT.kAzure # ROOT.TColor.GetColor(e, t, 128)
                    x = cluster_pos_z
                    y = cluster_pos_x
                    #print("cluster:", cluster, e, t, 255 * cluster_energy / 10.0, 255 * cluster_time / 10000.0, color, x, y, cluster_total_std_dev)
                    
                    marker = ROOT.TMarker(x, y, 53)
                    marker.SetMarkerColor(color)
                    marker.SetMarkerSize(20 * cluster_total_std_dev)
                    markers.append(marker)
                    marker = 0
                    
                    for ClusterHit in range(event.nHitsInClusterOne[cluster]):
                        hit_z = event.ClusterHitPosOne[ClusterHit*2 + 0] / 1000.0
                        hit_x = event.ClusterHitPosOne[ClusterHit*2 + 1] / 1000.0

                        marker = ROOT.TMarker(hit_z, hit_x, 21)
                        marker.SetMarkerColor(ROOT.kAzure-8)
                        markers.append(marker)
                        marker = 0

                for cluster in range(min(25, event.nClustersOther)):
                    cluster_energy = event.ClusterEnergyOther[cluster]
                    cluster_time = event.ClusterTimeOther[cluster]
                    cluster_pos_z = event.ClusterPosMeanOther[cluster*2 + 0] / 1000.0
                    cluster_pos_x = event.ClusterPosMeanOther[cluster*2 + 1] / 1000.0
                    cluster_pos_z_std_dev = event.ClusterPosStdDevOther[cluster*2 + 0] / 1000.0
                    cluster_pos_x_std_dev = event.ClusterPosStdDevOther[cluster*2 + 1] / 1000.0

                    cluster_total_std_dev = max(cluster_pos_z_std_dev, cluster_pos_x_std_dev)
                    e = int(min(255, 255*cluster_energy / 10.0))
                    t = int(min(255, 255*cluster_time / 10000.0))
                    color = ROOT.kMagenta
                    x = cluster_pos_z
                    y = cluster_pos_x

                    marker = ROOT.TMarker(x, y, 53)
                    marker.SetMarkerColor(color)
                    marker.SetMarkerSize(20 * cluster_total_std_dev)
                    markers.append(marker)
                    marker = 0

                    for ClusterHit in range(event.nHitsInClusterOther[cluster]):
                        hit_z = event.ClusterHitPosOther[ClusterHit*2 + 0] / 1000.0
                        hit_x = event.ClusterHitPosOther[ClusterHit*2 + 1] / 1000.0

                        marker = ROOT.TMarker(hit_z, hit_x, 21)
                        marker.SetMarkerColor(ROOT.kMagenta-8)
                        markers.append(marker)
                        marker = 0

                for line in range(event.nLinesOne):
                    #track_z = event.TrackHitPos[line*2 + 0] / 1000.0
                    #track_x = event.TrackHitPos[line*2 + 1] / 1000.0
                    track_z = event.FirstHoughHitOne[line*2 + 0] / 1000.0
                    track_x = event.FirstHoughHitOne[line*2 + 1] / 1000.0
                    track_z2 = event.LastHoughHitOne[line*2 + 0] / 1000.0
                    track_x2 = event.LastHoughHitOne[line*2 + 1] / 1000.0
                    track_length = event.TrackLengthOne[line] / 1000.0
                    direction_z = event.DirectionZOne[line]
                    direction_x = event.DirectionXOne[line]
                    
                    
                    if track_length == 0: continue
                    
                    x = track_z
                    y = track_x
                    x2 = track_z2
                    y2 = track_x2
                    
                    diffx = x2 - x
                    diffy = y2 - y
                    track_length_from_first_last = math.sqrt(diffx**2 + diffy**2)
                    
                    if not simplify_tracks: 
                        line = ROOT.TLine(x, y, x2, y2)
                        line.SetLineWidth(4)
                        line.SetLineStyle(ROOT.kDotted)
                        line.SetLineColorAlpha(ROOT.kBlue-6, 0.15)
                        markers.append(line)
                        line = 0
                    
                    # Add markers for front and end of track
                    offset = 0.075
                    marker_start = ROOT.TMarker(x, y - offset, 22)
                    marker_start.SetMarkerColor(ROOT.kBlue-10)
                    marker_start.SetMarkerSize(2)
                    markers.append(marker_start)
                    marker_start = 0
                    marker_end = ROOT.TMarker(x2, y2 + offset, 23)
                    marker_end.SetMarkerColor(ROOT.kBlue+3)
                    marker_end.SetMarkerSize(2)
                    markers.append(marker_end)
                    marker_end = 0

                    for TrackHit in range(event.nHitsInTrackOne[line]):
                        hit_z = event.TrackHitPosOne[TrackHit*2 + 0] / 1000.0
                        hit_x = event.TrackHitPosOne[TrackHit*2 + 1] / 1000.0

                        marker = ROOT.TMarker(hit_z, hit_x, 21)
                        marker.SetMarkerColor(ROOT.kAzure-3)
                        markers.append(marker)
                        marker = 0

                for line in range(event.nLinesOther):
                    
                    track_z = event.FirstHoughHitOther[line*2 + 0] / 1000.0
                    track_x = event.FirstHoughHitOther[line*2 + 1] / 1000.0
                    track_z2 = event.LastHoughHitOther[line*2 + 0] / 1000.0
                    track_x2 = event.LastHoughHitOther[line*2 + 1] / 1000.0
                    track_length = event.TrackLengthOther[line] / 1000.0
                    direction_z = event.DirectionZOther[line]
                    direction_x = event.DirectionXOther[line]

                    if track_length == 0: continue

                    x = track_z
                    y = track_x
                    x2 = track_z2
                    y2 = track_x2

                    diffx = x2 - x
                    diffy = y2 - y
                    track_length_from_first_last = math.sqrt(diffx**2 + diffy**2)

                    if not simplify_tracks:
                        line = ROOT.TLine(x, y, x2, y2)
                        line.SetLineWidth(4)
                        line.SetLineStyle(ROOT.kDashed)
                        line.SetLineColorAlpha(ROOT.kMagenta-6, 0.15)
                        markers.append(line)
                        line = 0

                    # Add markers for front and end of track
                    offset = 0.075
                    marker_start = ROOT.TMarker(x, y - offset, 22)
                    marker_start.SetMarkerColor(ROOT.kMagenta-10)
                    marker_start.SetMarkerSize(2)
                    markers.append(marker_start)
                    marker_start = 0
                    marker_end = ROOT.TMarker(x2, y2 + offset, 23)
                    marker_end.SetMarkerColor(ROOT.kMagenta+3)
                    marker_end.SetMarkerSize(2)
                    markers.append(marker_end)
                    marker_end = 0
                 
                    for TrackHit in range(event.nHitsInTrackOther[line]):
                        hit_z = event.TrackHitPosOther[TrackHit*2 + 0] / 1000.0
                        hit_x = event.TrackHitPosOther[TrackHit*2 + 1] / 1000.0
                        
                        marker = ROOT.TMarker(hit_z, hit_x, 21)
                        marker.SetMarkerColor(ROOT.kMagenta-3)
                        markers.append(marker)
                        marker = 0
            
        minimum_energy_to_print = 0
        if total_energy > minimum_energy_to_print:
            print(f"Time range: {time_low:0.0f}-{time_high:0.0f} ns,\tEnergy: {total_energy:0.2f} MeV")
            pad1.cd()
            haxis.Draw("colz axis")
            for marker in markers:
                marker.Draw()
            for text in text_to_draw:
                text.Draw()
            # Draw the time slice
            pad4.cd()
            htime.Draw("hist")
            
            canvas.cd()
            pad1.Draw()
            pad4.Draw()
            output_filename = os.path.join(out_dir, f"{name}_{current_spill_number:03d}")
            canvas.Print(output_filename + ".png")
            canvas.Clear("d")
            markers=[]
        
        
        
    return
        
    


    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', "-o", type=str, help="The output dir. Will be made if it doesn't exist. Default = spills/", default="spills")
    parser.add_argument('--name', "-n", type=str, help="The name of the output files. Will be <name>_<spill>_<slice>.png. Default = spill", default="spill")
    parser.add_argument('--input_filename', "-f", type=str, help="The file with the events to draw")
    parser.add_argument('--spillnum', "-s", type=int, help="The spill to draw. -1 for all", default=-1)
    parser.add_argument('--timeslice', "-t", type=int, help="The time slice to draw. -1 for all", default=-1)
    parser.add_argument('--readout_filename', "-r", type=str, help="(optional) A file with the raw readout.", default="")
    parser.add_argument('--only_true_tms_muons', help="Only draw events with true muons inside the TMS", action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    readout_filename = args.readout_filename
    only_true_tms_muons = args.only_true_tms_muons
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, only_true_tms_muons)

    
