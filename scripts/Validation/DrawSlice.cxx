
namespace DrawSliceN {
enum max_prints { disabled, handfull, few, many, tons, all };

static int max_slices = -1;

int GetMaxDrawSlicePrints(max_prints n) {
  int out = 0;
  switch (n) {
  case disabled:
    out = 0;
    break;
  case handfull:
    out = 5;
    break;
  case few:
    out = 20;
    break;
  case many:
    out = 50;
    break;
  case tons:
    out = 1000;
    break;
  case all:
    out = -1;
    break;
  default:
    out = 0;
    break;
  }
  // User can optionally limit the number of slices
  if (max_slices != -1 && (out > max_slices || out < 0))
    out = max_slices;
  return out;
}
} // namespace DrawSliceN

void ConfigureHist(TH1 &hist) {
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

TPaveText MakeTextBox(double x1, double y1, double x2, double y2,
                      double text_size = 0.08, int text_align = 12) {
  /*
  Create a transparent text box

  Parameters:
      x1, y1, x2, y2: Coordinates for the box (normalized to canvas size)
      text_size: Size of the text
      text_align: Text alignment (12 = left-center, etc.)
  */
  TPaveText text_box(x1, y1, x2, y2, "NDC");
  text_box.SetFillColor(0);
  text_box.SetFillStyle(0);
  text_box.SetBorderSize(0);
  text_box.SetTextSize(text_size);
  text_box.SetTextAlign(text_align);
  return text_box;
}

// Function to split a string by a delimiter
std::vector<std::string> split(const std::string &str, char delimiter) {
  std::vector<std::string> tokens;
  std::stringstream ss(str);
  std::string token;
  while (std::getline(ss, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

void AddMarkerToLegend(TLegend &leg, int marker_style, int marker_size,
                       int marker_color, TString text) {
  auto reco_entry = leg.AddEntry("", text, "p");
  reco_entry->SetMarkerStyle(marker_style);
  reco_entry->SetMarkerSize(marker_size);
  reco_entry->SetMarkerColor(marker_color);
}

void DrawSlice(std::string outfilename, std::string reason, std::string message,
               Reco_Tree &reco, Line_Candidates &lc, Truth_Info &truth,
               DrawSliceN::max_prints max_n_prints = DrawSliceN::all) {
  // Quit early if we already drew n copies of slices that have this reason
  static std::map<std::string, int> nSlicesDrawn;
  int nmax = GetMaxDrawSlicePrints(max_n_prints);
  if (nmax >= 0 && nSlicesDrawn[reason] >= nmax)
    return;
  nSlicesDrawn[reason] += 1;

  if (save_location == "") {
    throw std::runtime_error(
        "Do not have a save directory for the event displays in DrawSlice");
  }
  std::string directoryPath = save_location + "EventDisplays/" + reason + "/";
  createDirectory(directoryPath);
  static TCanvas *canvas;
  if (canvas == NULL) {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);
    canvas = new TCanvas();
    // canvas->SetLeftMargin(0.12);
    canvas->Divide(1, 2);
    for (size_t i = 1; i <= 2; i++) {
      auto pad = canvas->cd(i);
      pad->SetLeftMargin(0.42);
      pad->SetRightMargin(0.02);
    }
  }

  float buffer = 0;
  TH2D hist("", "X vs Z View;Z (cm);X (cm);N Hits", 100, 1100 - buffer,
            1900 + buffer, 100, -400 - buffer, 400 + buffer);
  TH2D histy("", "Y vs Z View;Z (cm);Y (cm);N Hits", 100, 1100 - buffer,
             1900 + buffer, 100, -500 - buffer, 100 + buffer);
  ConfigureHist(hist);
  ConfigureHist(histy);
  const double titleh = 0.09;
  gStyle->SetTitleH(titleh);

  std::vector<TMarker> markers;
  std::vector<TMarker> markersy;
  std::vector<TLine> lines;
  std::vector<TLine> linesy;

  bool draw_true_hits = true;
  if (draw_true_hits) {
    for (int ih = 0; ih < truth.NTrueHits; ih++) {
      if (truth.TrueRecoHitIsPedSupped[ih])
        continue;
      {
        float mx = truth.TrueHitZ[ih];
        float my = truth.TrueHitX[ih];
        TMarker marker(mx * CM, my * CM, 21);
        marker.SetMarkerColor(kGreen + 2);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markers.push_back(marker);
      }
      {
        float mx = truth.TrueHitZ[ih];
        float my = truth.TrueHitY[ih];
        TMarker marker(mx * CM, my * CM, 21);
        marker.SetMarkerColor(kGreen + 2);
        marker.SetMarkerStyle(33);
        marker.SetMarkerSize(0.75);
        markersy.push_back(marker);
      }
    }
  }

  // Get approx time slice range
  float hit_time_earliest = 1e9;
  float hit_time_latest = -1e9;
  for (int ih = 0; ih < lc.nHits; ih++) {
    if (lc.RecoHitPos[ih][3] > hit_time_latest)
      hit_time_latest = lc.RecoHitPos[ih][3];
    if (lc.RecoHitPos[ih][3] < hit_time_earliest)
      hit_time_earliest = lc.RecoHitPos[ih][3];
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][0];
      TMarker marker(mx * CM, my * CM, 21);
      marker.SetMarkerStyle(21);
      markers.push_back(marker);
    }
    {
      float mx = lc.RecoHitPos[ih][2];
      float my = lc.RecoHitPos[ih][1];
      TMarker marker(mx * CM, my * CM, 21);
      markersy.push_back(marker);
    }
  }

  int track_colors[] = {kRed, kBlue, kOrange - 3, kMagenta - 6, kGray};
  int n_track_colors = sizeof(track_colors) / sizeof(track_colors[0]);
  // This section draws reco tracks
  // Not as clean as line-base system in the next chunk
  bool draw_reco_track_nodes = false;
  if (draw_reco_track_nodes) {
    for (int it = 0; it < reco.nTracks; it++) {
      for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
        {
          float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
          float my = truth.RecoTrackTrueHitPosition[it][ih][0];
          TMarker marker(mx * CM, my * CM, 21);
          marker.SetMarkerColor(kGreen);
          marker.SetMarkerStyle(33);
          marker.SetMarkerSize(0.75);
          markers.push_back(marker);
        }
        {
          float mx = truth.RecoTrackTrueHitPosition[it][ih][2];
          float my = truth.RecoTrackTrueHitPosition[it][ih][1];
          TMarker marker(mx * CM, my * CM, 21);
          marker.SetMarkerColor(kGreen);
          marker.SetMarkerStyle(33);
          marker.SetMarkerSize(0.75);
          markersy.push_back(marker);
        }
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][0];
          TMarker marker(mx * CM, my * CM, 20);
          marker.SetMarkerColor(track_colors[it % n_track_colors]);
          marker.SetMarkerSize(0.5);
          markers.push_back(marker);
        }
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][1];
          TMarker marker(mx * CM, my * CM, 20);
          marker.SetMarkerColor(track_colors[it % n_track_colors]);
          marker.SetMarkerSize(0.5);
          markersy.push_back(marker);
        }
      }
    }
  }

  // This section draws reco tracks as lines between two points
  if (true) {
    for (int it = 0; it < reco.nTracks; it++) {
      auto color_to_use = track_colors[it % n_track_colors];
      for (int ih = 0; ih < reco.nKalmanNodes[it] - 1; ih++) {
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][0];
          float mx2 = reco.TrackHitPos[it][ih + 1][2];
          float my2 = reco.TrackHitPos[it][ih + 1][0];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(2);
          lines.push_back(line);
        }
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][1];
          float mx2 = reco.TrackHitPos[it][ih + 1][2];
          float my2 = reco.TrackHitPos[it][ih + 1][1];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(2);
          linesy.push_back(line);
        }
      }
    }
  }
  int empty_up_triangle = 55;
  int empty_down_triangle = 59;
  double marker_size = 1.5;
  // This section draws true particles
  bool should_draw_true_particles = true;
  if (should_draw_true_particles) {
    int colors[] = {kRed, kBlue, kOrange - 3, kMagenta - 6, kGray};
    int ncolors = sizeof(colors) / sizeof(colors[0]);
    int n_offset = 0;
    for (int ip = 0; ip < reco.nTracks; ip++) {
      int it = truth.RecoTrackPrimaryParticleIndex[ip];
      if (it < 0 || it > truth.nTrueParticles)
        continue; // Outside valid range
      // Only draw charged particles
      bool should_draw_true_particle = false;
      if (std::abs(truth.PDG[it]) == 13)
        should_draw_true_particle = true;
      // Pion
      if (std::abs(truth.PDG[it]) == 211)
        should_draw_true_particle = true;
      // Kaon
      if (std::abs(truth.PDG[it]) == 321)
        should_draw_true_particle = true;
      // Proton
      if (std::abs(truth.PDG[it]) == 2212)
        should_draw_true_particle = true;
      if (should_draw_true_particle) {
        int marker_to_use_start = 22;
        if (!truth.TMSFiducialStart[it])
          marker_to_use_start = empty_up_triangle;
        {
          float mx = truth.PositionZIsTMSStart[it][2];
          float my = truth.PositionZIsTMSStart[it][0];
          TMarker marker(mx * CM, my * CM, marker_to_use_start);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.PositionZIsTMSStart[it][2];
          float my = truth.PositionZIsTMSStart[it][1];
          TMarker marker(mx * CM, my * CM, marker_to_use_start);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markersy.push_back(marker);
        }
        int marker_to_use_end = empty_down_triangle;
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][0];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][0];
          }
          TMarker marker(mx * CM, my * CM, marker_to_use_end);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.PositionTMSEnd[it][2];
          float my = truth.PositionTMSEnd[it][1];
          if (mx > 19000) {
            mx = truth.PositionZIsTMSEnd[it][2];
            my = truth.PositionZIsTMSEnd[it][1];
          }
          TMarker marker(mx * CM, my * CM, marker_to_use_end);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markersy.push_back(marker);
        }
        int marker_to_use_death = 23;
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][0];
          TMarker marker(mx * CM, my * CM, marker_to_use_death);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.DeathPosition[it][2];
          float my = truth.DeathPosition[it][1];
          TMarker marker(mx * CM, my * CM, marker_to_use_death);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markersy.push_back(marker);
        }
        n_offset += 1;
      }
    }
  }

  // Create text boxes for each alignment
  double text_size = 0.08;
  double x1 = 0.01;
  double x2 = 0.28;
  double y1 = 0.02;
  double y2 = 0.9;
  auto textBoxA = MakeTextBox(x1, y1, x2, y2, text_size, 12); // Left-aligned
  auto textBoxC = MakeTextBox(x1, y1, x2, y2, text_size, 32); // Right-aligned

  textBoxA.AddText("Event Info:");
  textBoxC.AddText("");
  textBoxA.AddText(
      TString::Format("Entry: %d, Run: %d", reco.EventNo, truth.RunNo));
  textBoxC.AddText("");
  textBoxA.AddText(
      TString::Format("Spill: %d, Slice: %d", reco.SpillNo, reco.SliceNo));
  textBoxC.AddText("");
  for (auto foo : split(truth.Interaction->c_str(), ';')) {
    textBoxA.AddText(foo.c_str());
    textBoxC.AddText("");
  }

  y1 = 0.6;
  y2 = 1;
  auto textBoxD = MakeTextBox(x1, y1, x2, y2, text_size, 12); // Left-aligned
  auto textBoxE = MakeTextBox(x1, y1, x2, y2, text_size, 32); // Right-aligned
  textBoxD.AddText("N reco tracks:");
  textBoxE.AddText(TString::Format("%d", reco.nTracks));
  int time_slice_start_time = reco.TimeSliceStartTime;
  int time_slice_end_time = reco.TimeSliceEndTime;
  int slice_width = time_slice_end_time - time_slice_start_time;
  textBoxD.AddText("Slice start:");
  textBoxE.AddText(TString::Format("%dns", time_slice_start_time));
  textBoxD.AddText("Slice width:");
  textBoxE.AddText(TString::Format("%dns", slice_width));
  // Add message but only if it's not the default message
  if (message.find("n tracks =") == std::string::npos) {
    textBoxD.AddText(message.c_str());
    textBoxE.AddText("");
  }

  TLegend leg(x1, 0, x2 + 0.1, y1);
  leg.SetFillStyle(0);
  AddMarkerToLegend(leg, 21, 1.5, kBlack, "Reco Hits");
  AddMarkerToLegend(leg, 33, 1, kGreen + 2, "True Hits");
  if (draw_reco_track_nodes) {
    AddMarkerToLegend(leg, 20, 1.5, kRed, "Reco Track Nodes");
    AddMarkerToLegend(leg, 33, 1, kGreen + 2, "True Hits on Track");
  }
  auto reco_entry = leg.AddEntry("", "Reco Track", "l");
  reco_entry->SetLineColor(kRed);
  AddMarkerToLegend(leg, 22, 1, kBlack, "True Particle Start");
  AddMarkerToLegend(leg, 32, 1, kBlack, "True Part. Last TMS Pos");
  AddMarkerToLegend(leg, 23, 1, kBlack, "True Particle End");

  // Now draw
  canvas->cd(1);
  hist.Draw("Axis");
  for (auto &marker : markers) {
    marker.Draw();
  }
  for (auto &line : lines) {
    line.Draw();
  }
  textBoxA.Draw("same");
  textBoxC.Draw("same");
  canvas->cd(2);
  histy.Draw("Axis");
  for (auto &marker : markersy) {
    marker.Draw();
  }
  for (auto &line : linesy) {
    line.Draw();
  }
  textBoxD.Draw("same");
  textBoxE.Draw("same");
  leg.Draw();

  canvas->Print((directoryPath + outfilename + ".png").c_str());
}
