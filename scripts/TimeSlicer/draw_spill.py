import glob
import math
import os
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

def draw_spill(output_filename, filename, target_spill_number, target_slice_number = -1):
    #image = ROOT.TImage(1000, 1000)
    
    markers = []
    markers_not_in_slice = []
    lines = []
    text_to_draw = []
    
    c = ROOT.TChain("Line_Candidates")
    c.Add(filename)
    print("N entries:", c.GetEntries())
    assert c.GetEntries() > 0, "Didn't get any entries, are you sure the filename is right?\n" + filename
    
    truth = ROOT.TChain("Truth_Info")
    truth.Add(filename)
    assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info, are you sure the filename is right?\n" + filename
    
    bound_z_min = 11.0
    bound_z_max = 18.5000
    bound_x_min = -4.0
    bound_x_max = 4.0
    h2 = ROOT.TH2D("hits_in_event", ";Position Z (m);Position X (m);Total Energy (MeV)", 100, bound_z_min, bound_z_max, 100, bound_x_min, bound_x_max)
    
    htime = ROOT.TH1D("htime", ";Time (ns);Energy (MeV)", 500, -200, 10000)
    htime_highres = ROOT.TH1D("htime_highres", ";Time (ns);Energy (MeV)", 5100, -200, 10000)
    htime_spill = ROOT.TH1D("htime_spill", ";Time (ns);Energy (MeV)", 500, -200, 10000)
    htime_highres_spill = ROOT.TH1D("htime_highres_spill", ";Time (ns);Energy (MeV)", 5100, -200, 10000)
    htime.SetLineColor(ROOT.kGray)
    htime_highres.SetLineColor(ROOT.kGray)
    htime.SetFillColor(ROOT.kGray)
    htime_highres.SetFillColor(ROOT.kGray)
    htime_spill.SetLineColor(ROOT.kBlack)
    htime_highres_spill.SetLineColor(ROOT.kBlack)
    htime_spill.SetFillColor(ROOT.kBlack)
    htime_highres_spill.SetFillColor(ROOT.kBlack)
    htime_truth = ROOT.TH1D("htime_truth", ";Time (ns);Energy (MeV)", 500*2, -200, 10000)
    htime_highres_truth = ROOT.TH1D("htime_highres_truth", ";Time (ns);Energy (MeV)", 5100*2, -200, 10000)
    truth_color = ROOT.kMagenta + 1
    htime_truth.SetFillColor(truth_color)
    htime_truth.SetLineColor(truth_color)
    htime_highres_truth.SetFillColor(truth_color)
    htime_highres_truth.SetLineColor(truth_color)
    
    
    font_size = 0.15
    hists = (htime, htime_highres, htime_spill, htime_highres_spill)
    for hist in hists:
        hist.GetXaxis().SetTitleSize(font_size)
        hist.GetXaxis().SetLabelSize(font_size)
        hist.GetYaxis().SetTitleSize(font_size*0.8)
        hist.GetYaxis().SetLabelSize(font_size*0.8)
    htime_highres.GetXaxis().SetNdivisions(504)
    htime_highres_spill.GetXaxis().SetNdivisions(504)
    htime_highres.GetYaxis().SetNdivisions(504)
    htime_highres_spill.GetYaxis().SetNdivisions(504)
    htime.GetYaxis().SetTitleOffset(0.4)
    htime_spill.GetYaxis().SetTitleOffset(0.4)
    
    event_time = 1e10
    
    for i, event in enumerate(c):
        spill_number = event.SpillNo
        if spill_number < target_spill_number: continue
        if spill_number > target_spill_number: break
        if truth != None: truth.GetEntry(i)
        slice_number = event.SliceNo
            
        in_slice = False
        if target_slice_number == -1 or target_slice_number == slice_number: in_slice = True
        
        if truth != None and in_slice and slice_number != 0 and target_slice_number != -1:
            text = "Nu PDG = %s" % truth.NeutrinoPDG
            write_this = ROOT.TText(0.2,0.92, text);
            write_this.SetNDC()
            #write_this.SetTextColor(ROOT.kBlack)
            #write_this.SetTextFont(43);
            #write_this.SetTextSize(40);
            write_this.SetTextAlign(21);
            print("Writing:", text)
            text_to_draw.append(write_this)
            
            # Draw a star at the start of the true interaction, but stay within bounds
            true_x = truth.LeptonX4[0] / 1000.0
            true_z = truth.LeptonX4[2] / 1000.0
            true_t = truth.LeptonX4[3]
            bound_true_x = true_x
            bound_true_z = true_z
            buffer_amount = 0.03
            buffer_x = (bound_x_max - bound_x_min) * buffer_amount
            buffer_z = (bound_z_max - bound_z_min) * buffer_amount
            if true_x < bound_x_min - buffer_x: bound_true_x = bound_x_min - buffer_x
            if true_x > bound_x_max + buffer_x: bound_true_x = bound_x_min + buffer_x
            if true_z < bound_z_min - buffer_z: bound_true_z = bound_z_min - buffer_z
            if true_z > bound_z_max + buffer_z: bound_true_z = bound_z_min + buffer_z
            true_interaction_vtx_marker = ROOT.TMarker(bound_true_z, bound_true_x, 31)
            true_interaction_vtx_marker.SetMarkerColor(truth_color)
            true_interaction_vtx_marker.SetMarkerSize(3)
            markers.append(true_interaction_vtx_marker)
            
            htime_truth.Fill(true_t, 100)
            htime_highres_truth.Fill(true_t, 100)
        
        for hit in range(event.nHits):
            if event.RecoHitEnergy[hit] < 0: continue
            # Seems like pyROOT maps RecoHitPos[nHits][4] onto a RecoHitPos[nHits*4] array
            reco_hit_x = event.RecoHitPos[hit*4 + 0] / 1000.0
            reco_hit_y = event.RecoHitPos[hit*4 + 1] / 1000.0
            reco_hit_z = event.RecoHitPos[hit*4 + 2] / 1000.0
            reco_hit_time = event.RecoHitPos[hit*4 + 3]
            reco_hit_energy = event.RecoHitEnergy[hit]
            
            
            e = int(min(255, 255 * reco_hit_energy / 10.0))
            t = int(min(255, 255 * reco_hit_time / 10000.0))
            
            if in_slice:
                event_time = min(event_time, reco_hit_time)
            
            htime.Fill(reco_hit_time, reco_hit_energy)
            htime_highres.Fill(reco_hit_time, reco_hit_energy)
            if in_slice:
                htime_spill.Fill(reco_hit_time, reco_hit_energy)
                htime_highres_spill.Fill(reco_hit_time, reco_hit_energy)
            
            if in_slice: color = ROOT.TColor.GetColor(e, 32, 32)
            else: color = ROOT.kGray
            #print("hit:", e, t, 255 * reco_hit_energy / 10.0, 255 * reco_hit_time / 10000.0, color)
            x = reco_hit_z
            y = reco_hit_x
            
            marker = ROOT.TMarker(x, y, 21)
            marker.SetMarkerColor(color)
            if in_slice: markers.append(marker)
            else: markers_not_in_slice.append(marker)
            
        if target_slice_number == -1 or target_slice_number == slice_number:
            
            for cluster in range(min(25, event.nClusters)):
                cluster_energy = event.ClusterEnergy[cluster]
                cluster_time = event.ClusterTime[cluster]
                cluster_pos_z = event.ClusterPosMean[cluster*2 + 0] / 1000.0
                cluster_pos_x = event.ClusterPosMean[cluster*2 + 1] / 1000.0
                cluster_pos_z_std_dev = event.ClusterPosStdDev[cluster*2 + 0] / 1000.0
                cluster_pos_x_std_dev = event.ClusterPosStdDev[cluster*2 + 1] / 1000.0
                #cluster_total_std_dev = math.sqrt(cluster_pos_z_std_dev**2 + cluster_pos_x_std_dev**2)
                cluster_total_std_dev = max(cluster_pos_z_std_dev, cluster_pos_x_std_dev)
                e = int(min(255, 255 * cluster_energy / 10.0))
                t = int(min(255, 255 * cluster_time / 10000.0))
                color = ROOT.kBlack # ROOT.TColor.GetColor(e, t, 128)
                x = cluster_pos_z
                y = cluster_pos_x
                #print("cluster:", cluster, e, t, 255 * cluster_energy / 10.0, 255 * cluster_time / 10000.0, color, x, y, cluster_total_std_dev)
                
                marker = ROOT.TMarker(x, y, 53)
                marker.SetMarkerColor(color)
                marker.SetMarkerSize(20 * cluster_total_std_dev)
                markers.append(marker)
                
            for line in range(event.nLines):
                #track_z = event.TrackHitPos[line*2 + 0] / 1000.0
                #track_x = event.TrackHitPos[line*2 + 1] / 1000.0
                track_z = event.FirstHoughHit[line*2 + 0] / 1000.0
                track_x = event.FirstHoughHit[line*2 + 1] / 1000.0
                track_z2 = event.LastHoughHit[line*2 + 0] / 1000.0
                track_x2 = event.LastHoughHit[line*2 + 1] / 1000.0
                track_length = event.TrackLength[line] / 1000.0
                direction_z = event.DirectionZ[line]
                direction_x = event.DirectionX[line]
                
                
                if track_length == 0: continue
                
                #x = track_z
                #y = track_x
                #x2 = track_z + direction_z * track_length
                #y2 = track_x + direction_x * track_length
                x = track_z
                y = track_x
                x2 = track_z2
                y2 = track_x2
                
                diffx = x2 - x
                diffy = y2 - y
                track_length_from_first_last = math.sqrt(diffx**2 + diffy**2)
                
                #print("Line:", track_z, track_x, track_length, track_length_from_first_last, track_length_from_first_last / track_length, x, y, x2, y2)
                #print("Line lengths:", track_length, track_length_from_first_last, track_length_from_first_last / track_length)
                
                line = ROOT.TLine(x, y, x2, y2)
                line.SetLineWidth(4)
                line.SetLineColorAlpha(ROOT.kBlue, 0.15)
                markers.append(line)
                
                # Add markers for front and end of track
                offset = 0.075
                marker_start = ROOT.TMarker(x, y - offset, 22)
                marker_start.SetMarkerColor(ROOT.kGreen)
                marker_start.SetMarkerSize(2)
                markers.append(marker_start)
                marker_end = ROOT.TMarker(x2, y2 + offset, 23)
                marker_end.SetMarkerColor(ROOT.kRed)
                marker_end.SetMarkerSize(2)
                markers.append(marker_end)
            
            
            
    pad1.cd()
    h2.Draw("colz axis")
    for marker in markers_not_in_slice:
        marker.Draw()
    for marker in markers:
        marker.Draw()
    for text in text_to_draw:
        #print("Drawing", text, "at x =", text.GetX(), " y =", text.GetY())
        text.Draw()
    pad2.cd()
    htime.Draw("hist")
    htime_spill.Draw("hist same")
    htime_truth.Draw("hist same")
    htime_spill.Draw("hist same axis")
    pad3.cd()
    htime_highres.GetXaxis().SetRangeUser(event_time - 150, event_time + 150)
    htime_highres_spill.GetXaxis().SetRangeUser(event_time - 150, event_time + 150)
    htime_highres.Draw("hist Y+")
    htime_highres_spill.Draw("hist same")
    htime_highres_truth.Draw("hist same")
    htime_highres_spill.Draw("hist same axis")
    
    canvas.cd()
    pad1.Draw()
    pad2.Draw()
    pad3.Draw()
    output_dir = os.path.dirname(output_filename)
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    if target_slice_number < 0:
        canvas.Print(output_filename + ".png")
    else:
        canvas.Print(output_filename + "_slice{:03d}.png".format(target_slice_number))
        
    



def draw_spill_old(filename, target_spill_number, output_filename, target_slice_number = -1):
    assert False, "This code is no longer supported"
    c = ROOT.TChain("Line_Candidates")
    c.Add(filename)
    print("N entries:", c.GetEntries())
    assert c.GetEntries() > 0, "Didn't get any entries, are you sure the filename is right?\n" + filename
    
    t = ROOT.TChain("Truth_Info")
    t.Add(filename)
    assert t.GetEntries() > 0, "Didn't get any entries in Truth_Info, are you sure the filename is right?\n" + filename
    
    hits_in_event = ROOT.TH2D("hits_in_event", "Hits X vs Z;Position Z (m);Position X (m);Total Energy (MeV)", 100, 11.0, 18.5000, 100, -4.000, 4.000)
    hits_in_slice = ROOT.TH2D("hits_in_slice", "Hits X vs Z;Position Z (m);Position X (m);Total Energy (MeV)", 100, 11.0, 18.5000, 100, -4.000, 4.000)
    hit_time = ROOT.TH1D("hit_time", "Hit Time (ns);Hit Time (ns);N hits", 1000, 0, 10000)
    font_size = 0.15
    hit_time.GetXaxis().SetTitleSize(font_size)
    hit_time.GetXaxis().SetLabelSize(font_size)
    hit_time.GetYaxis().SetTitleSize(font_size*0.8)
    hit_time.GetYaxis().SetLabelSize(font_size*0.8)
    hit_time_in_slice = ROOT.TH1D("hit_time_in_slice", "Hit Time;Hit Time (ns);N hits", 1000, 0, 10000)
    hit_time_in_slice.SetFillColor(ROOT.kBlack)
    true_lepton_times = ROOT.TH1D("true_lepton_times", "True Lepton Time;Time (ns);Lepton Code", 1000*5, 0, 10000)
    true_lepton_times.SetFillColor(ROOT.kRed)
    true_lepton_times.SetLineColor(ROOT.kRed)
    true_lepton_times_nc = ROOT.TH1D("true_lepton_times_nc", "True Lepton Time;Time (ns);Lepton Code", 1000*5, 0, 10000)
    true_lepton_times_nc.SetFillColor(ROOT.kBlue)
    true_lepton_times_nc.SetLineColor(ROOT.kBlue)
    
    
    max_time_found = -1e9
    min_time_found = 1e9
    skip_slice_zero = True
    if target_slice_number == 0: skip_slice_zero = False
    found_target_slice = False
    if target_slice_number == -1: found_target_slice = True
    for i, event in enumerate(c):
        if t != None: t.GetEntry(i)
        spill_number = event.SpillNo
        if spill_number < target_spill_number: continue
        if spill_number > target_spill_number: break
        slice_number = event.SliceNo
        if slice_number == target_slice_number or target_slice_number == -1: found_target_slice = True
        # Let's just draw nearby slices for now
        should_skip = abs(slice_number - target_slice_number) > 3
        if skip_slice_zero and slice_number == 0: should_skip = True
        if not skip_slice_zero and slice_number == 0: should_skip = False
        if target_slice_number == 0: should_skip = False
        if target_slice_number == -1: should_skip = False
        if should_skip: continue
        for hit in range(event.nHits):
            if event.RecoHitEnergy[hit] < 0: continue
            # Seems like pyROOT maps RecoHitPos[nHits][4] onto a RecoHitPos[nHits*4] array
            reco_hit_x = event.RecoHitPos[hit*4 + 0] / 1000.0
            reco_hit_y = event.RecoHitPos[hit*4 + 1] / 1000.0
            reco_hit_z = event.RecoHitPos[hit*4 + 2] / 1000.0
            reco_hit_time = event.RecoHitPos[hit*4 + 3]
            reco_hit_energy = event.RecoHitEnergy[hit]
            
            # Slice 0 can have any time so don't track that
            if slice_number != 0:
                max_time_found = max(reco_hit_time, max_time_found)
                min_time_found = min(reco_hit_time, min_time_found)
            
            hits_in_event.Fill(reco_hit_z, reco_hit_x, reco_hit_energy)
            if slice_number == target_slice_number: hits_in_slice.Fill(reco_hit_z, reco_hit_x, reco_hit_energy)
            if slice_number == target_slice_number: hit_time_in_slice.Fill(reco_hit_time, reco_hit_energy)
            hit_time.Fill(reco_hit_time, reco_hit_energy)
        if t != None:
            true_lepton_time = t.LeptonX4[3];
            true_lepton_pdg = t.LeptonPDG;
            if abs(true_lepton_pdg) == 13: # Only muons
                true_lepton_times.Fill(true_lepton_time, 100)
            if abs(true_lepton_pdg) == 14: # Only muon neutrinos
                true_lepton_times_nc.Fill(true_lepton_time, 100)
            if abs(true_lepton_pdg) != 13 and abs(true_lepton_pdg) != 14:
                print("Found non-muon lepton pdg: ", true_lepton_pdg)
    
    #canvas.cd(1)
    pad1.cd()
    if hits_in_slice.Integral() > 0: 
        nbinsx = hits_in_event.GetNbinsX()
        nbinsy = hits_in_event.GetNbinsY()
        #m = ROOT.TMatrix(nbinsx+2, nbinsy+2, hits_in_event.GetArray(), "F")
        m = hits_in_event
        for x in range(0, nbinsx+2):
            for y in range(0, nbinsy+2):
                b = m.GetBin(x, y)
                value = m.GetBinContent(b)
                if value > 0:
                    m.SetBinContent(b, 0.1)
                value2 = hits_in_slice.GetBinContent(b)
                if value2 > 0:
                    hits_in_slice.SetBinContent(b, 1)
        m.Draw("col")
        m.GetZaxis().SetRangeUser(0, 1)
        hits_in_slice.Draw("col same")
    else: hits_in_event.Draw("colz")
    #canvas.cd(2)
    pad2.cd()
    max_to_use = max(true_lepton_times.GetMaximum(), true_lepton_times.GetMaximum(), hit_time.GetMaximum()) * 1.1
    hit_time.GetXaxis().SetRangeUser(min_time_found, max_time_found)
    hit_time.GetYaxis().SetRangeUser(0, max_to_use)
    hit_time.Draw("hist")
    if hit_time_in_slice.Integral() > 0: hit_time_in_slice.Draw("hist same")
    if t != None:
        true_lepton_times.Draw("hist same")
        true_lepton_times_nc.Draw("hist same")
        hit_time.Draw("axis same")
    
    canvas.cd()
    pad1.Draw()
    pad2.Draw()
    
    # Require some amount of minimum energy
    if target_slice_number == -1 and hits_in_event.Integral() < 10: return False
    
    if found_target_slice:
        if target_slice_number >= 0: print_filename = output_filename + "_{}_slice_{}.png".format(target_spill_number, target_slice_number)
        else: print_filename = output_filename + "_{}.png".format(target_spill_number)
        canvas.Print(print_filename)
        return True
    else: return False
    
    
if __name__ == "__main__":
    output_filename = sys.argv[1]
    filename = sys.argv[2]
    for i in range(100):
        draw_spill(output_filename, filename, 1, i)
    if False:
        if len(sys.argv) > 3:
            spill_number = int(sys.argv[3])
            saw_anything = draw_spill(filename, spill_number, output_filename)
            if saw_anything and False:
                for i in range(1000):
                    if not draw_spill(filename, spill_number, output_filename, i):
                        break
        else:
            for spill_number in range(1000):
                saw_anything = draw_spill(filename, spill_number, output_filename)
                if saw_anything:
                    for i in range(1000):
                        if not draw_spill(filename, spill_number, output_filename, i):
                            break

    
