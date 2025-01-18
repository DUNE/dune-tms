
namespace DrawSliceN {
  enum max_prints {
      disabled, 
      handfull,
      few,
      many,
      all
  };


  int GetMaxDrawSlicePrints(max_prints n) {
    switch (n) {
      case disabled: return 0;
      case handfull: return 5;
      case few: return 20;
      case many: return 50;
      case all: return -1;
      default: return 0;
    }
  }
}

void ConfigureHist(TH1& hist) {
    double textSize = 0.07;
    double titleOffset = 0.8;
    double titleOffsetX = 1.1;
    // Set axis title sizes
    hist.GetXaxis()->SetTitleSize(textSize);
    hist.GetYaxis()->SetTitleSize(textSize);

    // Set axis label sizes
    hist.GetXaxis()->SetLabelSize(textSize);
    hist.GetYaxis()->SetLabelSize(textSize);

    // Set axis title offsets
    hist.GetXaxis()->SetTitleOffset(titleOffsetX);
    hist.GetYaxis()->SetTitleOffset(titleOffset);
}

void DrawSlice(std::string outfilename, std::string reason, std::string message, Reco_Tree& reco, 
                Line_Candidates& lc, Truth_Info& truth, DrawSliceN::max_prints max_n_prints = DrawSliceN::all) { 
  // Quit early if we already drew n copies of slices that have this reason
  static std::map<std::string, int> nSlicesDrawn;
  int nmax = GetMaxDrawSlicePrints(max_n_prints);
  if (nmax >= 0 && nSlicesDrawn[reason] >= nmax) return;
  nSlicesDrawn[reason] += 1;

  if (save_location == "") {
    throw std::runtime_error("Do not have a save directory for the event displays in DrawSlice");
  }
  std::string directoryPath = save_location + "EventDisplays/" + reason + "/";
  createDirectory(directoryPath);
  static TCanvas* canvas;
  if (canvas == NULL) {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE); 
    canvas = new TCanvas();
    canvas->SetLeftMargin(0.12);
    canvas->Divide(1, 2);
  }
  
  float buffer = 0;
  TH2D hist("", "X vs Z View;Z (mm);X (mm);N Hits", 100, 11000 - buffer, 19000 + buffer, 100, -4000 - buffer, 4000 + buffer);
  TH2D histy("", "Y vs Z View;Z (mm);Y (mm);N Hits", 100, 11000 - buffer, 19000 + buffer, 100, -5000 - buffer, 1000 + buffer);
  ConfigureHist(hist);
  ConfigureHist(histy);
  const double titleh = 0.09;
  gStyle->SetTitleH(titleh);
  
  std::vector<TMarker> markers;
  std::vector<TMarker> markersy;
  // Get approx time slice range
  float hit_time_earliest = 1e9;
  float hit_time_latest = -1e9;
  for (int ih = 0; ih < lc.nHits; ih++) {
    if (lc.RecoHitPos[ih][3] > hit_time_latest) hit_time_latest = lc.RecoHitPos[ih][3];
    if (lc.RecoHitPos[ih][3] < hit_time_earliest) hit_time_earliest = lc.RecoHitPos[ih][3];
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][0];
      TMarker marker(mx, my, 21);
      marker.SetMarkerStyle(21);
      markers.push_back(marker);
    }
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][1];
      TMarker marker(mx, my, 21);
      markersy.push_back(marker);
    }
  }
  
  bool draw_true_hits = true;
  if (draw_true_hits) {
    for (int ih = 0; ih < truth.NTrueHits; ih++) {
      if (truth.TrueRecoHitIsPedSupped[ih]) continue;
      {
        float mx = truth.TrueHitZ[ih];
        float my = truth.TrueHitX[ih];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen+2);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markers.push_back(marker);
      }
      {
        float mx = truth.TrueHitZ[ih];
        float my = truth.TrueHitY[ih];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen+2);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markersy.push_back(marker);
      }
    }
  }
  
  // This section draws reco tracks
  int track_colors[] = {kRed, kRed + 2, kRed - 2, kRed + 4, kRed - 5};
  int n_track_colors = sizeof(track_colors) / sizeof(track_colors[0]);
  for (int it = 0; it < reco.nTracks; it++) {
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      {
        float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
        float my = truth.RecoTrackTrueHitPosition[it][ih][0];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markers.push_back(marker);
      }
      {
        float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
        float my = truth.RecoTrackTrueHitPosition[it][ih][1];
        TMarker marker(mx, my, 21);
        marker.SetMarkerColor(kGreen);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markersy.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][0];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(track_colors[it % n_track_colors]);
        marker.SetMarkerSize(0.5);
        markers.push_back(marker);
      }
      {
        float mx = reco.TrackHitPos[it][ih][2];
        float my = reco.TrackHitPos[it][ih][1];
        TMarker marker(mx, my, 20);
        marker.SetMarkerColor(track_colors[it % n_track_colors]);
        marker.SetMarkerSize(0.5);
        markersy.push_back(marker);
      }
    }
  }
  // This section draws true particles
  bool should_draw_true_particles = true;
  if (should_draw_true_particles) {
    int colors[] = {kGreen, kBlue, kOrange, kYellow, kMagenta, kGreen + 2};
    int ncolors = sizeof(colors)/sizeof(colors[0]);
    int n_offset = 0;
    for (int ip = 0; ip < reco.nTracks; ip++) {
      int it = truth.RecoTrackPrimaryParticleIndex[ip];
      if (it < 0 || it > truth.nTrueParticles) continue; // Outside valid range
      // Only draw charged particles
      bool should_draw_true_particle = false;
      if (std::abs(truth.PDG[it]) == 13) should_draw_true_particle = true;
      // Pion
      if (std::abs(truth.PDG[it]) == 211) should_draw_true_particle = true;
      // Kaon
      if (std::abs(truth.PDG[it]) == 321) should_draw_true_particle = true;
      // Proton
      if (std::abs(truth.PDG[it]) == 2212) should_draw_true_particle = true;
      if (should_draw_true_particle) {
        {
          float mx = truth.BirthPosition[it][2];
          float my = truth.BirthPosition[it][0];
          if (mx < 11000) {
            mx = truth.PositionZIsTMSStart[it][2];
            my = truth.PositionZIsTMSStart[it][0];
          }
          TMarker marker(mx, my, 22);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markers.push_back(marker);
        }
        {
          float mx = truth.BirthPosition[it][2];
          float my = truth.BirthPosition[it][1];
          if (mx < 11000) {
            mx = truth.PositionZIsTMSStart[it][2];
            my = truth.PositionZIsTMSStart[it][1];
          }
          TMarker marker(mx, my, 22);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markersy.push_back(marker);
        }
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][0];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][0];
          }
          TMarker marker(mx, my, 23);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markers.push_back(marker);
        }
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][1];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][1];
          }
          TMarker marker(mx, my, 23);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          markersy.push_back(marker);
        }
        n_offset += 1;
      }
    }
  }
  
  // Now draw
  canvas->cd(1);
  hist.Draw("Axis");
  for (auto& marker : markers) {
    marker.Draw();
  }
  canvas->cd(2);
  histy.Draw("Axis");
  for (auto& marker : markersy) {
    marker.Draw();
  }
  
  canvas->Print((directoryPath + outfilename + ".png").c_str());
}
