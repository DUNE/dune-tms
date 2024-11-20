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
import numpy


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
pad4 = ROOT.TPad("p4", "p4", 0.0,0.,0.7,0.21)
pad4.SetBottomMargin(0.35)
pad4.SetRightMargin(0.05)
pad4.SetTopMargin(0.05)
pad5 = ROOT.TPad("p5", "p5", 0.7,0.,1,0.21)
pad5.SetBottomMargin(0.35)
pad5.SetRightMargin(0.05)
pad5.SetLeftMargin(0.15)
pad5.SetTopMargin(0.05)


lar_start = 4179.24
lar_end = 9135.88
tms_start = 11362.0
tms_end = 18314.0
lar_x_start = -3478.48
lar_x_end = 3478.48
tms_x_start = -3300
tms_x_end = 3300
def inside_tms(x, y, z):
    is_inside = True
    if not -4000 < x < 4000: is_inside = False
    if not -1548-4000 < y < -1548+4000: is_inside = False
    if not 11000 < z < 18000: is_inside = False
    return is_inside


def inside_lar(x, y, z):
    is_inside = True
    if not -3478.48 < x < 3478.48: is_inside = False
    if not -2166 < y < -829.282: is_inside = False
    if not 4179.24 < z < 9135.88: is_inside = False
    return is_inside

def draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, interactions):
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
            
    only_true_tms_muons = False
    only_true_lar_muons = False
    if interactions.lower() == "tms":
        only_true_tms_muons = True
    if interactions.lower() == "lar":
        only_true_lar_muons = True
            
    bound_z_min = 4.0
    bound_z_max = 18.5000
    bound_x_min = -5.0
    bound_x_max = 5.0
    
    slice_colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kPink, ROOT.kCyan, ROOT.kMagenta, ROOT.kAzure]
    
    # TODO save this info in the file
    spill_time = 1.2e9
    
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
    
    truth_spill = ROOT.TChain("Truth_Spill")
    truth_spill.Add(input_filename)
    assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Spill, are you sure the input_filename is right?\n" + input_filename
    
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
        htime.GetYaxis().SetTitleOffset(0.35)
        htime.GetXaxis().SetTitleSize(font_size)
        htime.GetXaxis().SetLabelSize(font_size)
        htime.GetYaxis().SetTitleSize(font_size*0.8)
        htime.GetYaxis().SetLabelSize(font_size*0.8)
        htime.SetFillColor(ROOT.kBlack)
        htime.SetLineColor(ROOT.kBlack)
        
        max_hit_energy = 50
        henergy = ROOT.TH1D(f"henergy_{current_spill_number}", ";Hit E (MeV);N Hits", 50, 0, max_hit_energy)
        henergy.GetYaxis().SetTitleOffset(0.5)
        henergy.GetXaxis().SetTitleSize(font_size)
        henergy.GetXaxis().SetLabelSize(font_size)
        henergy.GetYaxis().SetTitleSize(font_size*0.8)
        henergy.GetYaxis().SetLabelSize(font_size*0.8)
        henergy.SetFillColor(ROOT.kBlack)
        henergy.SetLineColor(ROOT.kBlack)
        
        vertex_id_mapping = dict()
        slice_mapping = dict()
        
        total_energy = 0
        time_high = float("-inf")
        time_low = float("inf")
        n_hits = 0
        
        markers = []
        markers_in_slice = collections.defaultdict(list)
        markers_from_vertex_id = collections.defaultdict(list)
        text_to_draw = []
        
        if truth_spill != None: truth_spill.GetEntry(current_spill_number)
        # Sync up the readout info if it's there. Note that it has one entry per spill, not timeslice
        if readout != None: readout.GetEntry(current_spill_number)
        
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
            
            if only_true_tms_muons:
                assert truth != None, "Truth shouldn't be None inside only_true_tms_muons"
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
            
            if only_true_lar_muons:
                assert truth != None, "Truth shouldn't be None inside only_true_lar_muons"
                mx = truth.Muon_Vertex[0]
                my = truth.Muon_Vertex[1]
                mz = truth.Muon_Vertex[2]
                mdx = truth.Muon_Death[0]
                mdy = truth.Muon_Death[1]
                mdz = truth.Muon_Death[2]
                start_inside_lar = inside_lar(mx, my, mz)
                end_inside_tms = inside_tms(mdx, mdy, mdz)
                # Skip any events that don't start and stop inside the TMS
                if not start_inside_lar: continue
                if not end_inside_tms: continue
            
            print(f"Event {i} has {event.nHits} hits, and {event.nLines3D} lines.")
            
            for hit in range(event.nHits):
                if event.RecoHitEnergy[hit] < 0: continue
                # Seems like pyROOT maps RecoHitPos[nHits][4] onto a RecoHitPos[nHits*4] array
                reco_hit_x = event.RecoHitPos[hit*4 + 0] / 1000.0
                reco_hit_y = event.RecoHitPos[hit*4 + 1] / 1000.0
                reco_hit_z = event.RecoHitPos[hit*4 + 2] / 1000.0
                reco_hit_time = event.RecoHitPos[hit*4 + 3]
                reco_hit_energy = event.RecoHitEnergy[hit]
                reco_hit_slice = event.RecoHitSlice[hit]
                
            
                color_for_slice = ROOT.kGray if reco_hit_slice == 0 else slice_colors[reco_hit_slice % len(slice_colors)] + reco_hit_slice // len(slice_colors) - 3
                vertex_id_mapping[truth.VertexIdOfMostEnergyInEvent] = color_for_slice
                slice_mapping[truth.VertexIdOfMostEnergyInEvent] = reco_hit_slice
                
                n_hits += 1
                
                time_high = max(reco_hit_time, time_high)
                time_low = min(reco_hit_time, time_low)
                
                
                e = int(min(255, 255 * reco_hit_energy / 10.0))
                t = int(min(255, 255 * reco_hit_time / 10000.0))
                
                htime.Fill(reco_hit_time % spill_time, reco_hit_energy)
                henergy.Fill(reco_hit_energy if reco_hit_energy < max_hit_energy else 0.95 * max_hit_energy)
                x = reco_hit_z
                y = reco_hit_x
                
                marker = ROOT.TMarker(x, y, 21)
                #color = ROOT.TColor.GetColor(e, 32, 32)
                color = color_for_slice
                marker.SetMarkerColor(color)
                markers.append(marker)
                markers_in_slice[reco_hit_slice].append(marker)
                markers_from_vertex_id[truth.VertexIdOfMostEnergyInEvent].append(marker)
                
                if reco_hit_energy > 0:
                    total_energy += reco_hit_energy
                elif reco_hit_energy < 0:
                    print(f"Found unexpected reco_hit_energy < 0: {reco_hit_energy}")
                     
                    
        for hit in range(truth_spill.TrueNonTMSNHits):
            x = truth_spill.TrueNonTMSHitPos[hit*4 + 2] / 1000.0
            y = truth_spill.TrueNonTMSHitPos[hit*4 + 0] / 1000.0
            e = truth_spill.TrueNonTMSHitEnergy[hit]
            h = truth_spill.TrueNonTMSHitHadronicEnergy[hit]
            dedx = truth_spill.TrueNonTMSHitdEdx[hit]
            vid = truth_spill.TrueNonTMSHitVertexID[hit]
            a = e / 5
            if a > 1: a = 1
            if e < 3: a = 0.8
            if e < 1: a = 0.6
            if e < 0.5: a = 0.4
            if e < 0.1: a = 0.2
            if vid in vertex_id_mapping:
                marker = ROOT.TMarker(x, y, 21)
                color = ROOT.kGray if vid not in vertex_id_mapping else vertex_id_mapping[vid]
                marker.SetMarkerColorAlpha(color, a)
                #marker.SetMarkerColor(color)
                if h > 0.1: marker.SetMarkerStyle(20)
                markers.append(marker)
                vertex_slice = slice_mapping[vid]
                markers_in_slice[vertex_slice].append(marker)
                markers_from_vertex_id[vid].append(marker)
                #print(x, z)
                
        def make_box(bounds, color = ROOT.kMagenta+2, style = 1):
            out = ROOT.TBox(*bounds)
            out.SetLineColor(color)
            out.SetFillStyle(0)
            out.SetLineStyle(style)
            return out
        
        def draw(output_filename, markers, text_to_draw):
            minimum_energy_to_print = 0
            if total_energy > minimum_energy_to_print:
                print(f"Spill {current_spill_number},\tTime range: {time_low:0.0f}-{time_high:0.0f} ns,\tEnergy: {total_energy:0.2f} MeV,\tN hits: {n_hits}")
                pad1.cd()
                haxis.Draw("colz axis")
                for marker in markers:
                    marker.Draw()
                for text in text_to_draw:
                    text.Draw()
                    
                M = 1e-3
                box_lar = make_box((lar_start * M, lar_x_start * M, lar_end * M, lar_x_end * M))
                box_lar_had = make_box((lar_start * M + 0.3, lar_x_start * M + 0.3, lar_end * M - 0.3, lar_x_end * M - 0.3), style = 2)
                box_lar_fiducial = make_box((lar_start * M + 1.50, lar_x_start * M + 0.5, lar_end * M - 1.50, lar_x_end * M - 0.5), ROOT.kRed)
                box_tms = make_box((tms_start * M, tms_x_start * M, tms_end * M, tms_x_end * M))
                box_lar.Draw()
                box_lar_had.Draw()
                box_lar_fiducial.Draw()
                box_tms.Draw()
                    
                # Draw the time slice
                pad4.cd()
                htime.Draw("hist")
                # Draw the time slice
                pad5.cd()
                henergy.Draw("hist")
                
                
                canvas.cd()
                pad1.Draw()
                pad4.Draw()
                pad5.Draw()
                canvas.Print(output_filename)
                
        output_filename = os.path.join(out_dir, f"{name}_{current_spill_number:03d}.png")
        draw(output_filename, markers, text_to_draw)
        #for s in markers_in_slice:
        #    output_filename_slice = os.path.join(out_dir, f"{name}_{current_spill_number:03d}_s{s:02d}.png")
        #    draw(output_filename_slice, markers_in_slice[s], text_to_draw)
        #for v in markers_from_vertex_id:
        #    output_filename_vertex = os.path.join(out_dir, f"{name}_{current_spill_number:03d}_v{v:02d}.png")
        #    draw(output_filename_vertex, markers_from_vertex_id[v], text_to_draw)
        
        
        
        
    return
        
    


    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draws spills.')
    parser.add_argument('--outdir', "-o", type=str, help="The output dir. Will be made if it doesn't exist. Default = spills/", default="spills")
    parser.add_argument('--name', "-n", type=str, help="The name of the output files. Will be <name>_<spill>_<slice>.png. Default = spill", default="spill")
    parser.add_argument('--input_filename', "-f", type=str, help="The file with the events to draw")
    parser.add_argument('--spillnum', "-s", type=int, help="The spill to draw. -1 for all", default=-1)
    parser.add_argument('--timeslice', "-t", type=int, help="The time slice to draw. -1 for all", default=-1)
    parser.add_argument('--readout_filename', "-r", type=str, help="(optional) A file with the raw readout.", default="")
    parser.add_argument('--interactions', "-i", type=str, help="(optional) The subdet with interactions", default="")

    args = parser.parse_args()

    out_dir = args.outdir
    name = args.name
    input_filename = args.input_filename
    spill_number = args.spillnum
    time_slice = args.timeslice
    readout_filename = args.readout_filename
    interactions = args.interactions
    draw_spill(out_dir, name, input_filename, spill_number, time_slice, readout_filename, interactions)

    
