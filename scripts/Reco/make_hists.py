import glob
import argparse
import os
import math
import ROOT
import array
# Tells root to be in batch mode so it doesn't try to create a canvas on your screen, which is slow
ROOT.gROOT.SetBatch(True)
# Don't draw the stats box on histograms
ROOT.gStyle.SetOptStat(0)
# Python's wrapper around ROOT sometimes crashes when it automatically adds hists to its garbage collector. This fixed that
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)

def inside_tms(x, y, z, only_thin_section = False):
  """ Returns true if x,y,z are inside the TMS
  Currently rough estimates of the positions """
  is_inside = True
  if not -3000 < x < 3000: is_inside = False
  if not -2800 < y < 0: is_inside = False
  if only_thin_section:
    if not 11400 < z < 13500: is_inside = False
  else:
    if not 11400 < z < 18200: is_inside = False
  return is_inside 

def run(c, truth, outfilename, nmax=-1):
  """ This code does 3 things:
  Makes histograms
  Loops over all events and fills histograms
  Saves histograms in outfilename
  """
  
  # First make some histograms
  # You can add any histograms you want. Just make sure that histogram name hasn't been used before.
    
  # Hists related to track length
  hist_TrackLengthOne = ROOT.TH1D("hist_TrackLength_one", "TrackLength;Track Length (g/cm^{2});N Tracks", 50, 0, 2500)
  hist_track_lengthOne = ROOT.TH1D("hist_track_length_one", "Track Length;Track Length (cm);N Tracks", 50, 0, 2500)
  hist_track_length_longestOne = ROOT.TH1D("hist_track_length_longest_one", "Track Length of Muon Candidate;Track Length (cm);N Tracks", 50, 10, 10000)
  hist_track_length_max_dz = ROOT.TH1D("hist_track_length_max_dz", "Longest Track Length by Max dz;Track Length (cm);N Tracks", 50, 10, 10000)

  hist_TrackLengthOther = ROOT.TH1D("hist_TrackLength_other", "TrackLength;TrackLength (g/cm^{2};N Tracks", 50, 0, 2500)
  hist_track_lengthOther = ROOT.TH1D("hist_track_length_other", "Track Length;Track Length (cm);N Tracks", 50, 0, 2500);
  hist_track_length_longestOther = ROOT.TH1D("hist_track_length_longest_other", "Track Length of Muon Candidate;Track Length (cm);N Tracks", 50, 10, 10000)
#  hist_track_length_max_dzOther = ROOT.TH1D("hist_track_length_max_dz_other", "Longest Track Length by Max dz;Track Length (cm);N Tracks", 50, 10, 10000)
  # Other hists
  hist_n_tracksOne = ROOT.TH1D("hist_n_tracks_one", "N Tracks;N Tracks;N Events", 6, -0.5, 5.5)
  hist_n_hits = ROOT.TH1D("hist_n_hits", "N Hits;N Hits;N Events", 101, 0, 100)
  hist_occupancyOne = ROOT.TH1D("hist_occupancy_one", "Occupancy;Occupancy of longest track;N Events", 110, 0, 1.1)

  hist_n_tracksOther = ROOT.TH1D("hist_n_tracks_other", "N Tracks;N Tracks; N Events", 6, -0.5, 5.5)
  hist_occupancyOther = ROOT.TH1D("hist_occupancy_other", "Occupancy;Occupancy of longest track;N Events", 110, 0, 1.1)
  # Track start and stop positions
  hist_track_start_xOne = ROOT.TH1D("hist_track_start_x_one", "Track Start Position X;X (mm);N Events", 100, -4000, 4000)
  hist_track_start_zOne = ROOT.TH1D("hist_track_start_z_one", "Track Start Position Z;X (mm);N Events", 100, 11000, 18000)
  hist_track_end_xOne = ROOT.TH1D("hist_track_end_x_one", "Track End Position X;X (mm);N Events", 100, -4000, 4000)
  hist_track_end_zOne = ROOT.TH1D("hist_track_end_z_one", "Track End Position Z;X (mm);N Events", 100, 11000, 18000)
  hist_track_startOne = ROOT.TH2D("hist_track_start_one", 
    "Track Start Position;Track Start Position X (mm);Track Start Position Z (mm)", 100, 11000, 18000, 100, -4000, 4000)
  hist_track_endOne = ROOT.TH2D("hist_track_end_one", 
    "Track End Position;Track End Position X (mm);Track End Position Z (mm)", 100, 11000, 18000, 100, -4000, 4000)
  hist_track_angleOne = ROOT.TH1D("hist_track_angle_one", "Track Angle;Angle (deg);N Events", 50, -180, 180)

  hist_track_start_xOther = ROOT.TH1D("hist_track_start_x_other", "Track Start Position X;X (mm);N Events", 100, -4000, 4000)
  hist_track_start_zOther = ROOT.TH1D("hist_track_start_z_other", "Track Start Position Z;X (mm);N Events", 100, 11000, 18000)
  hist_track_end_xOther = ROOT.TH1D("hist_track_end_x_other", "Track End Position X;X (mm);N Events", 100, -4000, 4000)
  hist_track_end_zOther = ROOT.TH1D("hist_track_end_z_other", "Track End Position Z;X (mm);N Events", 100, 11000, 18000)
  hist_track_startOther = ROOT.TH2D("hist_track_start_other",
    "Track Start Position;Track Start Position X (mm);Track Start Position Z (mm)", 100, 11000, 18000, 100, -4000, 4000)
  hist_track_endOther = ROOT.TH2D("hist_track_end_other",
    "Track End Position;Track End Position X (mm); Track End Position Z (mm)", 100, 11000, 18000, 100, -4000, 4000)
  hist_track_angleOther = ROOT.TH1D("hist_track_angle_other", "Track Angle;Angle (deg);N Events", 50, -180, 180)

  # Now make histograms that rely on truth information as well
  if truth != None:
    # KE estimators
    hist_KEOne = ROOT.TH2D("hist_KE_one", "KE;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 2500)
    hist_KE_estimatedOne = ROOT.TH2D("hist_KE_estimated_one", "KE estimator;True muon KE (MeV);Reco muon KE (MeV)", 100, 0, 5000, 50, 0, 5000)
    hist_KE_vs_track_length_max_dz = ROOT.TH2D("hist_KE_vs_track_length_max_dz", "KE vs Track Length;True muon KE (MeV);Track Length (mm)", 100, 0, 5000, 50, 0, 10000)
    hist_track_length_vs_max_dz_distOne = ROOT.TH2D("hist_track_length_vs_max_dz_dist_one", "Track Dist;Track Length (mm);Max dz Dist (mm)", 50, 0, 10000, 50, 0, 10000)
    hist_KE_inside_TMS = ROOT.TH1D("hist_KE_inside_TMS", "KE of Muons Starting in TMS;True muon KE (MeV);N events", 100, 0, 5000)

    hist_KEOther = ROOT.TH2D("hist_KE_other", "KE;True muon KE (MeV);Track length of best track (g/cm^{2}=)", 100, 0, 5000, 50, 0, 2500)
    hist_KE_estimatedOther = ROOT.TH2D("hist_KE_estimated_other", "KE estimator;True muon KE (MeV);Reco muon KE (MeV)", 100, 0, 5000, 50, 0, 5000)
#    hist_KE_vs_track_length_max_dzOther = ROOT.TH2D("hist_KE_vs_track_length_max_dz_other", "KE vs Track Length;True muon KE (MeV);Track Length (mm)", 100, 0, 5000, 50, 0, 10000)
    hist_track_length_vs_max_dz_distOther = ROOT.TH2D("hist_track_length_vs_max_dz_dist_other", "Track Dist;Track Length (mm);Max dz Dist (mm)", 50, 0, 10000, 50, 0, 10000)
#    hist_KE_inside_TMSOther = ROOT.TH1D("hist_KE_inside_TMS_other", "KE of Muons Starting in TMS;True muon KE (MeV);N events", 100, 0, 5000)
   
    # Additional estimators that only require a true muon start and stop inside the TMS (it's "interior")
    # These are not ideal because they're using truth information, but it's a nice way to test for issues with reco
    hist_true_interior_muon_ke_vs_track_length_dz = ROOT.TH2D("hist_true_interior_muon_ke_vs_track_length_dz", 
      "True KE vs Dist Between Widest dz Hits;True muon KE (MeV);Dist between two widest dz hits (mm)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_max_distOne = ROOT.TH2D("hist_true_interior_muon_ke_vs_max_dist_one", 
      "True KE vs Longest Track;True muon KE (MeV);Longest Track (mm)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_max_areal_densityOne = ROOT.TH2D("hist_true_interior_muon_ke_vs_max_areal_density_one", 
      "True KE vs Largest Areal Density Reco Track;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_estimated_keOne = ROOT.TH2D("hist_true_interior_muon_ke_vs_estimated_ke_one", 
      "True KE vs 3.5*(areal density);True muon KE (MeV);Estimated KE (MeV)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_estimated_ke_originalOne = ROOT.TH2D("hist_true_interior_muon_ke_vs_estimated_ke_original_one", 
      "True KE vs 82+1.75*(areal density);True muon KE (MeV);Estimated KE (MeV)", 100, 0, 5000, 50, 0, 5000)
   
#    hist_true_interior_muon_ks_vs_track_length_dzOther = ROOT.TH2D("hist_true_interior_muon_ke_vs_track_length_dz_other",
#      "True KE vs Dist Between Widest dz Hits;True muon KE (MeV);Dist between two widest dz hits (mm)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_max_distOther = ROOT.TH2D("hist_true_interior_muon_ke_vs_max_dist_other",
      "True KE vs Longest Track;True muon KE (MeV);Longest Track (mm)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_max_areal_densityOther = ROOT.TH2D("hist_true_interior_muon_ke_vs_max_areal_density_other",
      "True KE vs Largest Areal Density Reco Track;True muon KE (MeV);Track length of best track (g/cm^{2})", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_estimated_keOther = ROOT.TH2D("hist_true_interior_muon_ke_vs_estimated_ke_other",
      "True KE vs 3.5*(areal density);True muon KE (MeV);Estimated KE (MeV)", 100, 0, 5000, 50, 0, 5000)
    hist_true_interior_muon_ke_vs_estimated_ke_originalOther = ROOT.TH2D("hist_true_interior_muon_ke_vs_estimated_ke_original_other",
      "True KE vs 82+1.75*(areal density);True muon KE (MeV);Estimated KE (MeV)", 100, 0, 5000, 50, 0, 5000)                                                                   

    # Vertex Resolution
    hist_track_start_vtx_z_resolutionOne = ROOT.TH1D("hist_track_start_vtx_z_resolution_one", 
      "Track Start Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_start_vtx_z_resolution_using_spanOne = ROOT.TH1D("hist_track_start_vtx_z_resolution_using_span_one", 
      "Track Start Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_end_vtx_z_resolutionOne = ROOT.TH1D("hist_track_end_vtx_z_resolution_one", "End Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_end_vtx_z_resolution_using_spanOne = ROOT.TH1D("hist_track_end_vtx_z_resolution_using_span_one", 
      "Track End Vtx Resolution Z;Reco - True Vtx Z (mm);N events", 51, -1000, 1000)
    hist_track_start_vtx_x_resolutionOne = ROOT.TH1D("hist_track_start_vtx_x_resolution_one", 
      "Track Start Vtx Resolution X;Reco - True Vtx X (mm);N events", 51, -1000, 1000)
    hist_track_end_vtx_x_resolutionOne = ROOT.TH1D("hist_track_end_vtx_x_resolution_one", 
      "Track End Vtx Resolution X;Reco - True Vtx X (mm); N events", 51, -1000, 1000)
     
    hist_track_start_vtx_z_resolutionOther = ROOT.TH1D("hist_track_start_vtx_z_resolution_other",
      "Track Start Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_start_vtx_z_resolution_using_spanOther = ROOT.TH1D("hist_track_start_vtx_z_resolution_using_span_other",
      "Track Start Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_end_vtx_z_resolutionOther = ROOT.TH1D("hist_track_end_vtx_z_resolution_other", "End Vtx Resolution Z;Reco - True Vtx Z (mm); N events", 51, -1000, 1000)
    hist_track_end_vtx_z_resolution_using_spanOther = ROOT.TH1D("hist_track_end_vtx_z_resolution_using_span_other",
      "Track End Vtx Resolution Z;Reco - True Vtx Z (mm);N events", 51, -1000, 1000)
    hist_track_start_vtx_x_resolutionOther = ROOT.TH1D("hist_track_start_vtx_x_resolution_other",
      "Track Start Vtx Resolution X;Reco - True Vtx X (mm);N events", 51, -1000, 1000)
    hist_track_end_vtx_x_resolutionOther = ROOT.TH1D("hist_track_end_vtx_x_resolution_other",
      "Track End Vtx Resolution X; Reco - True Vtx X (mm); N events", 51, -1000, 1000)                                               
    # Can also calculate efficiency
    bin_edges = array.array('d', [0, 200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,3000,4000,5000])
    hist_eff_track_finding_numeratorOne = ROOT.TH1D("hist_eff_track_finding_numerator_one", "N Tracks Found;True KE (MeV);N Muons", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_after_cuts_numeratorOne = ROOT.TH1D("hist_eff_track_finding_after_cuts_numerator_one", "N Tracks After Cuts;True KE (MeV);N Muons", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_denominatorOne = ROOT.TH1D("hist_eff_track_finding_denominator_one", "N Tracks;True KE (MeV);N Muons",len(bin_edges) - 1, bin_edges)
    hist_eff_track_findingOne = ROOT.TH1D("hist_eff_track_finding_one", "Eff. of Reco'ing TMS-Starting Muons;True KE (MeV);Eff", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_after_cutsOne = ROOT.TH1D("hist_eff_track_finding_after_cuts_one", 
      "Eff. of Reco'ing TMS-Starting Muons After Cuts;True KE (MeV);Eff",len(bin_edges) - 1, bin_edges)
    
    hist_eff_track_finding_numeratorOther = ROOT.TH1D("hist_eff_track_finding_numerator_other", "N Tracks Found;True KE (MeV);N Muons", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_after_cuts_numeratorOther = ROOT.TH1D("hist_eff_track_finding_after_cuts_numerator_other", "N Tracks After Cuts; True KE (MeV);N Muons", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_denominatorOther = ROOT.TH1D("hist_eff_track_finding_denominator_other", "N Tracks;True KE (MeV);N Muons", len(bin_edges) - 1, bin_edges)
    hist_eff_track_findingOther = ROOT.TH1D("hist_eff_track_finding_other", "Eff. of Reco'ing TMS-Starting Muons;True KE (MeV);Eff", len(bin_edges) - 1, bin_edges)
    hist_eff_track_finding_after_cutsOther = ROOT.TH1D("hist_eff_track_finding_after_cuts_other",
      "Eff. of Reco'ing TMS-Starting Muons After Cuts;True KE (MeV);Eff", len(bin_edges) - 1, bin_edges)                                                 

    # We can also do some simple counts
    n_start_inside_tms = 0 
    n_end_inside_tms = 0
    n_start_and_end_inside_tms = 0
    n_true_muons = 0

    # User can request fewer events, so check how many we're looping over.
    nevents = c.GetEntries()
    if nmax >= 0 and nevents > nmax: nevents = nmax
    # Figure out how often to print progress information.
    # Setting carriage = True will use a carriage return which keeps the progress on a single line
    # But if you add print statements, it will be ugly so it's not default
    carriage = False
    if carriage:
      if nevents <= 100: print_every = 1 
      elif nevents <= 1000: print_every = 10
      elif nevents <= 10000: print_every = 100
      else: print_every = 1000
    else:
      if nevents < 100: print_every = 1 # Print every event if < 100 events
      elif 100 <= nevents < 1000: print_every = 20
      elif 1000 <= nevents < 10000: print_every = 100
      elif 10000 <= nevents < 100000: print_every = 1000
      else: print_every = 10000
     
    # Now loop over all events
    for i, event in enumerate(c):
      if i > nevents: break
      if truth != None: truth.GetEntry(i)
      # Print current progress, with carriage return \r to prevent long list of progress and have everything on a singe line.
      if i % print_every == 0 and carriage: print(f"\rOn {i} / {nevents}", end='')
      if i % print_every == 0 and not carriage: print(f"On {i} / {nevents}")
        
      # Get some basic info
      nLinesOne = event.nLinesOne
      ntracksOne = nLinesOne
      nhits = event.nHits
       
      nLinesOther = event.nLinesOther
      ntracksOther = nLinesOther

      # And fill that basic info
      hist_n_tracksOne.Fill(ntracksOne)
      hist_n_hits.Fill(nhits)
      
      hist_n_tracksOther.Fill(ntracksOther)

      # This finds the two hits with the largest dz distance between them.
      # Then it calculates the distance between the two using dz and dx
      # If there was a nice muon, this would be roughly equivalent to the track length
      # So this is used as a reconstruction-independent track length estimator
      # in case there's an issue with the reconstruction.
      max_dz = -999
      track_length_max_dz = -999
      min_z_hit = -999
      max_z_hit = -999
      max_z = -1e9
      min_z = 1e9
  
      for hit in range(nhits):
        x = event.RecoHitPos[hit * 4 + 0]
        z = event.RecoHitPos[hit * 4 + 2]
        if z > max_z:
          max_z_hit = hit
          max_z = z
        if z < min_z:
          min_z_hit = hit
          min_z = z
      # Now if we found two hits, calculate the distance
      if max_z_hit >= 0 and min_z_hit < 1e8:
        x1 = event.RecoHitPos[min_z_hit * 4 + 0]
        z1 = event.RecoHitPos[min_z_hit * 4 + 2]
        x2 = event.RecoHitPos[max_z_hit * 4 + 0]
        z2 = event.RecoHitPos[max_z_hit * 4 + 2]
        track_length_max_dz = math.sqrt((x2 - x1)**2 + (z2 - z1)**2)
     
        # Loop over all tracks
        # Find the track with highest areal density (event.TrackLength)
        # This will be our primary muon candidate, longtrack
        # Also find the track length by first and last hough hit positions
        # This is not as good of an estimator because it's not taking into
        # account the various amounts of material in the path.
        # It's also not taking into account the curvature due to magnetic field
        # This code also plots all track lengths in a histogram
        # And saves them in an array for later use
        longtrack_lengthOne = 0
        longtrackOne = -1
        track_lengthsOne = [0] * ntracksOne
        areal_densitiesOne = [0] * ntracksOne

        longtrack_lengthOther = 0
        longtrackOther = -1
        track_lengthsOther = [0] * ntracksOther
        areal_densitiesOther = [0] * ntracksOther

        for track in range(ntracksOne):
          track_start_x = event.FirstHoughHitOne[2*track+0]
          track_start_y = event.FirstHoughHitOne[2*track+1]
          track_end_x = event.LastHoughHitOne[2*track+0]
          track_end_y = event.LastHoughHitOne[2*track+1]
          xdist = track_end_x - track_start_x
          ydist = track_end_y - track_start_y
          dist = math.sqrt(xdist*xdist+ydist*ydist)
          TrackLength = event.TrackLengthOne[track]
          
          # Check if this is the longest track by areal density
          if TrackLength > longtrack_lengthOne:
            longtrackOne = track
            longtrack_lengthOne = TrackLength
            
          # Plot them in histograms
          hist_track_lengthOne.Fill(dist)
          hist_TrackLengthOne.Fill(TrackLength)
            
          # Save them in these arrays
          track_lengthsOne[track] = dist
          areal_densitiesOne[track] = event.TrackLengthOne[track]

        for track in range(ntracksOther):
          track_start_x = event.FirstHoughHitOther[2*track+0]
          track_start_y = event.FirstHoughHitOther[2*track+1]
          track_end_x = event.LastHoughHitOther[2*track+0]
          track_end_y = event.LastHoughHitOther[2*track+1]
          xdist = track_end_x - track_start_x
          ydist = track_end_y - track_start_y
          dist = math.sqrt(xdist*xdist+ydist*ydist)
          TrackLength = event.TrackLengthOther[track]

          # Check if this is the longest track by areal density
          if TrackLength > longtrack_lengthOther:
            longtrackOther = track
            longtrack_lengthOther = TrackLength

          # Plot them in histograms
          hist_track_lengthOther.Fill(dist)
          hist_TrackLengthOther.Fill(TrackLength)

          # Save them in these arrays
          track_lengthsOther[track] = dist
          areal_densitiesOther[track] = event.TrackLengthOther[track]

        best_tracklengthOne = -999.0
        reconstructed_muon_keOne = -999.0
        if longtrackOne >= 0:
          hist_track_length_longestOne.Fill(track_lengthsOne[longtrackOne])
        hist_track_length_max_dz.Fill(track_length_max_dz)
         
        best_tracklengthOther = -999.9
        reconstructed_muon_keOther = -999.0
        if longtrackOther >= 0:
          hist_track_length_longestOther.Fill(track_lengthsOther[longtrackOther])
#        hist_track_length_max_dzOther.Fill(track_length_max_dzOther)
          
        # If there's a muon candidate, plot its occupancy, and start and stop positions
        if longtrackOne >= 0:
          occupancy = event.OccupancyOne[longtrackOne]
          hist_occupancyOne.Fill(occupancy)
          track_start_z = event.FirstHoughHitOne[2*longtrackOne+0]
          track_end_z = event.LastHoughHitOne[2*longtrackOne+0]
          track_start_x = event.FirstHoughHitOne[2*longtrackOne+1]
          track_end_x = event.LastHoughHitOne[2*longtrackOne+1]
          hist_track_start_xOne.Fill(track_start_x)
          hist_track_start_zOne.Fill(track_start_z)
          hist_track_startOne.Fill(track_start_z, track_start_x)
          hist_track_end_xOne.Fill(track_end_x)
          hist_track_end_zOne.Fill(track_end_z)
          hist_track_endOne.Fill(track_end_z, track_end_x)
            
          slope = event.SlopeOne[longtrackOne]
          angle = math.atan(slope) * 180 / math.pi
          hist_track_angleOne.Fill(angle)
        
        if longtrackOther >= 0:
          occupancy = event.OccupancyOther[longtrackOther]
          hist_occupancyOther.Fill(occupancy)
          track_start_z = event.FirstHoughHitOther[2*longtrackOther+0]
          track_end_z = event.LastHoughHitOther[2*longtrackOther+0]
          track_start_x = event.FirstHoughHitOther[2*longtrackOther+1]
          track_end_x = event.LastHoughHitOther[2*longtrackOther+1]
          hist_track_start_xOther.Fill(track_start_x)
          hist_track_start_zOther.Fill(track_start_z)
          hist_track_startOther.Fill(track_start_z, track_start_x)
          hist_track_end_xOther.Fill(track_end_x)
          hist_track_end_zOther.Fill(track_end_z)
          hist_track_endOther.Fill(track_end_z, track_end_x)

          slope = event.SlopeOther[longtrackOther]
          angle = math.atan(slope) * 180 / math.pi
          hist_track_angleOther.Fill(angle)

        # Now check if the muon is a good candidate
        # Does this by checking a list of conditions. If any of them are false, don't fill the reco vs true ke plots.
        best_muon_candidateOne = longtrackOne
        best_muon_candidateOther = longtrackOther
        if longtrackOne >= 0:
          OccupancyCut = 0.0  #0.5
          alldet = True
          track_start_z = event.FirstHoughHitOne[2*longtrackOne+0]
          track_end_z = event.LastHoughHitOne[2*longtrackOne+0]
          track_start_x = event.FirstHoughHitOne[2*longtrackOne+1]
          track_end_x = event.LastHoughHitOne[2*longtrackOne+1]
          # Check if the muon starts within the correct z range
          if alldet and track_start_z < 11362+55*2: best_muon_candidateOne = -1
          if not alldet and (track_start_z < 11362+55*2 or track_start_z > 13600): best_muon_candidateOne = -1 
          if track_end_z > 18294-80*2: best_muon_candidateOne = -1
          # Check if the muon is in the right x position
          if abs(track_start_x) > 3520-200: best_muon_candidateOne = -1
          if abs(track_end_x) > 3520-200: best_muon_candidateOne = -1
          # Look only at events with true muons that die inside the detector
          if truth != None and (truth.Muon_Death[1] > 1159 or truth.Muon_Death[1] > -3864): best_muon_candidateOne = -1 
          # Make sure that at least OccupancyCut of the hits are from the candidate
          if event.OccupancyOne[longtrackOne] < OccupancyCut: best_muon_candidateOne = -1
          
        if longtrackOther >= 0:
          OccupancyCut = 0.0  #0.5
          alldet = True
          track_start_z = event.FirstHoughHitOther[2*longtrackOther+0]
          track_end_z = event.LastHoughHitOther[2*longtrackOther+0]
          track_start_x = event.FirstHoughHitOther[2*longtrackOther+1]
          track_end_x = event.LastHoughHitOther[2*longtrackOther+1]
          # Check if the muon starts within the correct z range
          if alldet and track_start_z < 11362+55*2: best_muon_candidateOther = -1
          if not alldet and (track_start_z < 11362+55*2 or track_start_z > 13600): best_muon_candidateOther = -1
          if track_end_z > 18294-80*2: best_muon_candidateOther = -1
          # Check if the muon is in the right x position
          if abs(track_start_x) > 3520-200: best_muon_candidateOther = -1
          if abs(track_end_x) > 3520-200: best_muon_candidateOther = -1
          # Look only at events with true muons that die inside the detector
          if truth != None and (truth.Muon_Death[1] > 1159 or truth.Muon_Death[1] > -3864): best_muon_candidateOther = -1
          # Make sure that at least OccupancyCut of the hits are from the candidate
          if event.OccupancyOther[longtrackOther] < OccupancyCut: best_muon_candidateOther = -1
            
        # Only fill if we found a good candidate
        if best_muon_candidateOne >= 0:
          best_tracklengthOne = event.TrackLengthOne[best_muon_candidateOne]
          reconstructed_muon_keOne = 82+1.75*best_tracklengthOne
        
        if best_muon_candidateOther >= 0:
          best_tracklengthOther = event.TrackLengthOther[best_muon_candidateOther]
          reconstructed_muon_keOther = 82+1.75*best_tracklengthOther

        track_lengthOne = -999
        track_distOne = -999
        track_lengthOther = -999
        track_distOther = -999

        if best_muon_candidateOne >= 0: 
          track_lengthOne = event.TrackLengthOne[best_muon_candidateOne]
          track_distOne = track_lengthsOne[best_muon_candidateOne]
        elif event.nLinesOne > 0: 
          track_lengthOne = event.TrackLengthOne[0]
          track_distOne = track_lengthsOne[0]

        if best_muon_candidateOther >= 0:
          track_lengthOther = event.TrackLengthOther[best_muon_candidateOther]
          track_distOther = track_lengthsOther[best_muon_candidateOther]
        elif event.nLinesOther > 0:
          track_lengthOther = event.TrackLengthOther[0]
          track_distOther = track_lengthsOther[0]
        #if track_length >= 0 or track_length_max_dz >= 0:
        #    print(f"True muon KE: {truth.Muon_TrueKE:0.2f}\tPath length: {track_length:0.2f}\tTrack Dist: {track_dist:0.2f}\tMax dz: {track_length_max_dz:0.2f}")
        hist_track_length_vs_max_dz_distOne.Fill(track_distOne, track_length_max_dz)
        hist_track_length_vs_max_dz_distOther.Fill(track_distOther, track_length_max_dz)
      
        # Get information that relies on truth information as well
        if truth != None:
          has_true_muon = truth.Muon_TrueKE >= 0
          Muon_TrueKE = truth.Muon_TrueKE
          #if has_true_muon: print(Muon_TrueKE, has_true_muon)
            
          # Only fill these hists if there is a true muon
          if has_true_muon: 
            hist_KEOne.Fill(Muon_TrueKE, best_tracklengthOne)
            hist_KEOther.Fill(Muon_TrueKE, best_tracklengthOther)
            hist_KE_estimatedOne.Fill(Muon_TrueKE, reconstructed_muon_keOne)
            hist_KE_estimatedOther.Fill(Muon_TrueKE, reconstructed_muon_keOther)
            hist_KE_vs_track_length_max_dz.Fill(Muon_TrueKE, track_length_max_dz)
#            hist_KE_vs_track_length_max_dzOther.Fill(Muon_TrueKE, track_length_max_dzOther)
          
            # Now check if the muon started and ended in the TMS
            mx = truth.Muon_Vertex[0]
            my = truth.Muon_Vertex[1]
            mz = truth.Muon_Vertex[2]
            mdx = truth.Muon_Death[0]
            mdy = truth.Muon_Death[1]
            mdz = truth.Muon_Death[2]
            start_inside_tms = inside_tms(mx, my, mz, True)
            end_inside_tms = inside_tms(mdx, mdy, mdz)
            #print(f"Muon Start XYZ: ({mx:0.2f}, {my:0.2f}, {mz:0.2f})\tMuon End XYZ: ({mdx:0.2f}, {mdy:0.2f}, {mdz:0.2f})\tstart_inside_tms: {start_inside_tms}\tend_inside_tms: {end_inside_tms}")
          
            # Does some simple counts which it prints later
            n_true_muons += 1
            if start_inside_tms: n_start_inside_tms += 1
            if end_inside_tms: n_end_inside_tms += 1
            if start_inside_tms and end_inside_tms: n_start_and_end_inside_tms += 1
         
            if start_inside_tms:
              hist_KE_inside_TMS.Fill(Muon_TrueKE)
            # Plot KE estimators for muons which start and end in the detector
            if start_inside_tms and end_inside_tms:
              hist_true_interior_muon_ke_vs_track_length_dz.Fill(Muon_TrueKE, track_length_max_dz)
              if len(track_lengthsOne) > 0:
                hist_true_interior_muon_ke_vs_max_distOne.Fill(Muon_TrueKE, max(track_lengthsOne))
              if len(track_lengthsOther) > 0:
                hist_true_interior_muon_ke_vs_max_distOther.Fill(Muon_TrueKE, max(track_lengthsOther))
              if len(areal_densitiesOne) > 0:
                hist_true_interior_muon_ke_vs_max_areal_densityOne.Fill(Muon_TrueKE, max(areal_densitiesOne))
                hist_true_interior_muon_ke_vs_estimated_ke_originalOne.Fill(Muon_TrueKE, 82+1.75*max(areal_densitiesOne))
                hist_true_interior_muon_ke_vs_estimated_keOne.Fill(Muon_TrueKE, 2*1.75*max(areal_densitiesOne))
              if len(areal_densitiesOther) > 0:
                hist_true_interior_muon_ke_vs_max_areal_densityOther.Fill(Muon_TrueKE, max(areal_densitiesOther))
                hist_true_interior_muon_ke_vs_estimated_ke_originalOther.Fill(Muon_TrueKE, 82+1.75*max(areal_densitiesOther))
                hist_true_interior_muon_ke_vs_estimated_keOther.Fill(Muon_TrueKE, 2*1.75*max(areal_densitiesOther))
             
            # Plot vertex resolution of muon
            # But only if it starts or ends in detector respectively
            if longtrackOne >= 0:
              dz_start = event.FirstHoughHitOne[2*longtrackOne+0] - mz 
              dz_start_span = z1 - mz
              dz_end = event.LastHoughHitOne[2*longtrackOne+0] - mdz
              dz_end_span = z2 - mdz
              dx_start = event.FirstHoughHitOne[2*longtrackOne+1] - mx
              dx_end = event.LastHoughHitOne[2*longtrackOne+1] - mdx
              if start_inside_tms:
                hist_track_start_vtx_z_resolutionOne.Fill(dz_start)
                hist_track_start_vtx_z_resolution_using_spanOne.Fill(dz_start_span)
                hist_track_start_vtx_x_resolutionOne.Fill(dx_start)
              if end_inside_tms: 
                hist_track_end_vtx_z_resolutionOne.Fill(dz_end)
                hist_track_end_vtx_z_resolution_using_spanOne.Fill(dz_end_span)
                hist_track_end_vtx_x_resolutionOne.Fill(dx_end)
            
            if longtrackOther >= 0:
              dz_start = event.FirstHoughHitOther[2*longtrackOther+0] - mz
              dz_start_span = z1 - mz
              dz_end = event.LastHoughHitOther[2*longtrackOther+0] - mdz
              dz_end_span = z2 - mdz
              dx_start = event.FirstHoughHitOther[2*longtrackOther+1] - mx
              dx_end = event.LastHoughHitOther[2*longtrackOther+1] - mdx
              if start_inside_tms:
                hist_track_start_vtx_z_resolutionOther.Fill(dz_start)
                hist_track_start_vtx_z_resolution_using_spanOther.Fill(dz_start_span)
                hist_track_start_vtx_x_resolutionOther.Fill(dx_start)
              if end_inside_tms:
                hist_track_end_vtx_z_resolutionOther.Fill(dz_end)
                hist_track_end_vtx_z_resolution_using_spanOther.Fill(dz_end_span)
                hist_track_end_vtx_x_resolutionOther.Fill(dx_end)
              
            # For finding efficiency
            # Plot true KE in all cases
            # But for numerator, look if we found a candidate
            # For now only consider muons starting in the TMS
            if start_inside_tms:
              if longtrackOne >= 0:
                hist_eff_track_finding_numeratorOne.Fill(Muon_TrueKE)
              #if longtrackOther >= 0:
              #  hist_eff_track_finding_numeratorOther.Fill(Muon_TrueKE)
              if best_muon_candidateOne >= 0:
                hist_eff_track_finding_after_cuts_numeratorOne.Fill(Muon_TrueKE)
              hist_eff_track_finding_denominatorOne.Fill(Muon_TrueKE)
              if longtrackOther >= 0:
                hist_eff_track_finding_numeratorOther.Fill(Muon_TrueKE)
              #if longtrackOther >= 0:
              #  hist_eff_track_finding_numeratorOther.Fill(Muon_TrueKE)
              if best_muon_candidateOther >= 0:
                hist_eff_track_finding_after_cuts_numeratorOther.Fill(Muon_TrueKE)
              hist_eff_track_finding_denominatorOther.Fill(Muon_TrueKE)
      
      
    ## Calculate the efficiency
    hist_eff_track_findingOne.Divide(hist_eff_track_finding_numeratorOne, hist_eff_track_finding_denominatorOne)
    hist_eff_track_finding_after_cutsOne.Divide(hist_eff_track_finding_after_cuts_numeratorOne, hist_eff_track_finding_denominatorOne)

    hist_eff_track_findingOther.Divide(hist_eff_track_finding_numeratorOther, hist_eff_track_finding_denominatorOther)
    hist_eff_track_finding_after_cutsOther.Divide(hist_eff_track_finding_after_cuts_numeratorOther, hist_eff_track_finding_denominatorOther)
  
    ## Now save all the histograms
    if carriage: print("\r", end="") # Erase the current line
    print(f"Completed loop of {nevents} events. Saving to {outfilename}")
    hists_to_save = []
    # This is a pythonistic cheat to quickly find the histograms
    # and add them to the list of files to save
    # Otherwise we'd add them one by one, like: hists_to_save.append(hist1), etc.
    for name, item in locals().items():
      if "hist_" in name:
        hists_to_save.append(item)
       
    tf = ROOT.TFile(outfilename, "recreate")
    for hist in hists_to_save:
      hist.Write()
    tf.Close()
    print("Done saving.")
  
    if truth != None:
      print(f"N true muons: {n_true_muons}")
      print(f"N true muons that start inside TMS: {n_start_inside_tms}")
      print(f"N true muons that end inside TMS: {n_end_inside_tms}")
      print(f"N true muons that start and end inside TMS: {n_start_and_end_inside_tms}")
   
    # Return this hists if the user requested previews
    return hists_to_save 
  

def validate_then_run(args):
  ''' This function has all the code to validate the arguments. Usually users don't need to edit this code. '''
  # First we validate the args
  indir = args.indir
  inlist = args.inlist
  infile = args.filename
  files_to_use = []
  if indir != "":
    # User specified an input directory. So loop through all the files to find tmsreco files.
    print(f"Finding all tmsreco.root files in {indir}")
    for root, dirs, files in os.walk(indir):
      for name in files:
        if "tmsreco.root" in name:
          # This is a tmsreco file. So add it to the list
          fullfilename = os.path.join(root, name)
          files_to_use.append(fullfilename)
    nfiles = len(files_to_use)
    if nfiles == 0:
      raise ValueError(f"Did not find any files in {indir}")
    print(f"Found {nfiles} files in {indir}")
  if inlist != "":
    # In this case the user specified a text file with the full paths
    with open(inlist) as f:
      file_data = f.read()
      files_to_use = file_data.splitlines()
    nfiles = len(files_to_use)
    if nfiles == 0:
      raise ValueError(f"Did not find any files in {inlist}")
    print(f"Found {nfiles} files in {inlist}")
  if infile != "":
    # In this case, the user specified exactly one file. Usually they'd hadd many files together.
    files_to_use = [infile]
  
  outdir = args.outdir
  if outdir == "":
    # No output directory was specified so use the default
    # First we need the username
    username = os.environ["USER"]
    outdir = f"/dune/data/users/{username}/dune-tms_hists"
  else:
    # Check if it follows the correct conventions
    good_locations = ["/dune/data/users", "/dune/data2/users", "/pnfs/dune/persistent", "/pnfs/dune/scratch"]
    if not any(location in outdir for location in good_locations):
      print(f"Warning: outdir is not in list of good locations. Don't want to write root files to app area. {outdir}")
  # Make sure the output directory exists
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  outname = args.name
  if ".root" not in outname:
    raise ValueError(f"The output file should be a root file. {outname}")
  outfilename = os.path.join(outdir, outname)
  if os.path.exists(outfilename):
    if args.allow_overwrite:
      print(f"Warning: A file already exists at {outfilename}, but you specified allow_overwrite. It will be overwritten")
    else:
      raise ValueError(f"The output file already exists at {outfilename}. Please specify another location, allow_overwrite or delete the file yourself.")
  print(f"Output will be in {outfilename}")
    
  has_truth = True
   
  # Make the TChain objects. One for truth information and one for reconstructed information.
  c = ROOT.TChain("Line_Candidates")
  truth = None
  if has_truth: truth = ROOT.TChain("Truth_Info")
  for f in files_to_use:
    c.Add(f)
    if has_truth: truth.Add(f)
  assert c.GetEntries() > 0, "Didn't get any entries in Line_Candidates TChain." 
  if has_truth: assert truth.GetEntries() > 0, "Didn't get any entries in Truth_Info TChain."
    
  nevents = args.nevents
  assert nevents >= -1, f"nevents <= -1, why? {nevents}"
    
  # Now finally run
  hists = run(c, truth, outfilename, nevents)
  
  preview = args.preview
  if preview:
    print("Making preview histograms.")
    print("WARNING: Eventually you'll want to make plots using a separate plotting script")
    print("That way you can format the plots without having to rerun the make_hists.py script")
    # Make the output directory if needed
    outname_without_root = outname.replace(".root", "")
    preview_dir = os.path.join(outdir, "previews", outname_without_root)
    if not os.path.exists(preview_dir): 
      os.makedirs(preview_dir)
    print(f"Saving previews in {preview_dir}/<hist_name>.png")
    canvas = ROOT.TCanvas()
    for hist in hists:
      name = hist.GetName().replace("hist_", "")
      opts = "colz"
      hist.Draw(opts)
      canvas.Print(os.path.join(preview_dir, f"{name}.png"))








# This is how python knows whether it's a script that's being imported into another script or running as the main script
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Draws spills.')
  parser.add_argument('--outdir', type=str, help="The output dir. Will be made if it doesn't exist.", default="")
  parser.add_argument('--name', type=str, help="The name of the output files.", default="dune-tms_hists.root")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--indir', type=str, help="The location of the input tmsreco files", default="")
  group.add_argument('--inlist', type=str, help="The input filelist", default="")
  group.add_argument('--filename', '-f', type=str, help="The input file, if you have a single file", default="")
  parser.add_argument('--nevents', '-n', type=int, help="The maximum number of events to loop over", default=-1)
  parser.add_argument('--allow_overwrite', help="Allow the output file to overwrite", action=argparse.BooleanOptionalAction)
  parser.add_argument('--preview', help="Save preview images of the histograms", action=argparse.BooleanOptionalAction)
   
  args = parser.parse_args()
   
  validate_then_run(args)
