
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
    out = 500;
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
  auto legend_marker = new TMarker(0, 0, marker_style);
  auto reco_entry = leg.AddEntry(legend_marker, text, "p");
  reco_entry->SetMarkerStyle(marker_style);
  reco_entry->SetMarkerSize(marker_size);
  reco_entry->SetMarkerColor(marker_color);
}

void AddLineToLegend(TLegend &leg, int line_width, int line_color,
                     TString text) {
  auto legend_line = new TLine(0, 0, 1, 0);
  auto entry = leg.AddEntry(legend_line, text, "l");
  entry->SetLineWidth(line_width);
  entry->SetLineColor(line_color);
}

const double BOUNDS_X_START = -400;
const double BOUNDS_X_END = 400;
const double BOUNDS_Y_START = -500;
const double BOUNDS_Y_END = 100;
const double BOUNDS_Z_START = 1100;
const double BOUNDS_Z_END = 1900;

void ConstrainMarkerX(float &mx, float &my) {
  double buffer = 5;
  if (mx > BOUNDS_Z_END - buffer)
    mx = BOUNDS_Z_END - buffer;
  if (mx < BOUNDS_Z_START + buffer)
    mx = BOUNDS_Z_START + buffer;
  if (my > BOUNDS_X_END - buffer)
    my = BOUNDS_X_END - buffer;
  if (my < BOUNDS_X_START + buffer)
    my = BOUNDS_X_START + buffer;
}

void ConstrainMarkerY(float &mx, float &my) {
  double buffer = 5;
  if (mx > BOUNDS_Z_END - buffer)
    mx = BOUNDS_Z_END - buffer;
  if (mx < BOUNDS_Z_START + buffer)
    mx = BOUNDS_Z_START + buffer;
  if (my > BOUNDS_Y_END - buffer)
    my = BOUNDS_Y_END - buffer;
  if (my < BOUNDS_Y_START + buffer)
    my = BOUNDS_Y_START + buffer;
}

void AddStageGuide(std::vector<TLine> &lines, double y) {
  TLine line(BOUNDS_Z_START, y, BOUNDS_Z_END, y);
  line.SetLineColor(kGray + 1);
  line.SetLineStyle(3);
  line.SetLineWidth(1);
  lines.push_back(line);
}

void DrawSliceReco3DStages(std::string outfilename, std::string reason,
                           std::string message, Reco_Tree &reco,
                           Line_Candidates &lc, Truth_Info &truth,
                           DrawSliceN::max_prints max_n_prints =
                               DrawSliceN::all) {
  static std::map<std::string, int> nSlicesDrawn;
  int nmax = GetMaxDrawSlicePrints(max_n_prints);
  if (nmax >= 0 && nSlicesDrawn[reason] >= nmax)
    return;
  nSlicesDrawn[reason] += 1;

  if (save_location == "") {
    throw std::runtime_error("Do not have a save directory for the event "
                             "displays in DrawSliceReco3DStages");
  }

  std::string directoryPath =
      save_location + "EventDisplays/" + reason + "_reco3d/";
  createDirectory(directoryPath);

  static TCanvas *canvas;
  if (canvas == NULL) {
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);
    canvas = new TCanvas("draw_slice_reco3d_stages",
                         "draw_slice_reco3d_stages", 1200, 900);
    canvas->Divide(1, 2);
    for (size_t i = 1; i <= 2; i++) {
      auto pad = canvas->cd(i);
      pad->SetLeftMargin(0.30);
      pad->SetRightMargin(0.02);
    }
  }

  const double raw_offset = 320;
  const double track2d_offset = 0;
  const double track3d_offset = -320;
  const double zmin = BOUNDS_Z_START;
  const double zmax = BOUNDS_Z_END;

  TH2D hist_x("", "X vs Z Reconstruction Stages;Z (cm);X + stage offset (cm)",
              100, zmin, zmax, 100, -760, 760);
  TH2D hist_y("", "Y vs Z Reconstruction Stages;Z (cm);Y + stage offset (cm)",
              100, zmin, zmax, 100, -860, 520);
  ConfigureHist(hist_x);
  ConfigureHist(hist_y);
  hist_x.GetYaxis()->SetTitleOffset(1.1);
  hist_y.GetYaxis()->SetTitleOffset(1.1);
  gStyle->SetTitleH(0.09);

  std::vector<TMarker> markers_x;
  std::vector<TMarker> markers_y;
  std::vector<TLine> lines_x;
  std::vector<TLine> lines_y;

  AddStageGuide(lines_x, raw_offset);
  AddStageGuide(lines_x, track2d_offset);
  AddStageGuide(lines_x, track3d_offset);
  AddStageGuide(lines_y, raw_offset);
  AddStageGuide(lines_y, track2d_offset);
  AddStageGuide(lines_y, track3d_offset);

  int color_2d = kBlue + 1;
  int color_3d = kRed;
  int color_kalman = kMagenta + 2;
  int color_2d_pale = TColor::GetColor("#b9d6ff");
  int color_3d_pale = TColor::GetColor("#ffc6c6");
  int color_kalman_pale = TColor::GetColor("#d8c8ff");

  for (int it = 0; it < lc.nLinesU; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackU[it]; ih++) {
      TMarker marker(lc.TrackHitPosU[it][ih][0] * CM,
                     lc.TrackHitPosU[it][ih][1] * CM + raw_offset, 21);
      marker.SetMarkerColor(color_2d_pale);
      marker.SetMarkerSize(1.35);
      markers_x.push_back(marker);
    }
  }
  for (int it = 0; it < lc.nLinesV; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackV[it]; ih++) {
      TMarker marker(lc.TrackHitPosV[it][ih][0] * CM,
                     lc.TrackHitPosV[it][ih][1] * CM + raw_offset, 25);
      marker.SetMarkerColor(color_2d_pale);
      marker.SetMarkerSize(1.35);
      markers_x.push_back(marker);
    }
  }
  for (int it = 0; it < lc.nLinesX; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackX[it]; ih++) {
      TMarker marker(lc.TrackHitPosX[it][ih][0] * CM,
                     lc.TrackHitPosX[it][ih][1] * CM + raw_offset, 21);
      marker.SetMarkerColor(color_2d_pale);
      marker.SetMarkerSize(1.35);
      markers_y.push_back(marker);
    }
  }

  for (int ih = 0; ih < lc.nHits; ih++) {
    {
      TMarker marker(lc.RecoHitPos[ih][2] * CM,
                     lc.RecoHitPos[ih][0] * CM + raw_offset, 20);
      marker.SetMarkerColor(kBlack);
      marker.SetMarkerSize(0.45);
      markers_x.push_back(marker);
    }
    {
      TMarker marker(lc.RecoHitPos[ih][2] * CM,
                     lc.RecoHitPos[ih][1] * CM + raw_offset, 20);
      marker.SetMarkerColor(kBlack);
      marker.SetMarkerSize(0.45);
      markers_y.push_back(marker);
    }
  }

  for (int it = 0; it < reco.nTracks; it++) {
    for (int ih = 0; ih < reco.nHits[it]; ih++) {
      {
        TMarker marker(reco.TrackHitPos[it][ih][2] * CM,
                       reco.TrackHitPos[it][ih][0] * CM + track2d_offset, 20);
        marker.SetMarkerColor(color_3d_pale);
        marker.SetMarkerSize(1.35);
        markers_x.push_back(marker);
      }
      {
        TMarker marker(reco.TrackHitPos[it][ih][2] * CM,
                       reco.TrackHitPos[it][ih][1] * CM + track2d_offset, 20);
        marker.SetMarkerColor(color_3d_pale);
        marker.SetMarkerSize(1.35);
        markers_y.push_back(marker);
      }
    }
  }

  for (int it = 0; it < lc.nLinesU; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackU[it]; ih++) {
      TMarker marker(lc.TrackHitPosU[it][ih][0] * CM,
                     lc.TrackHitPosU[it][ih][1] * CM + track2d_offset, 21);
      marker.SetMarkerColor(color_2d);
      marker.SetMarkerSize(0.75);
      markers_x.push_back(marker);
    }
  }
  for (int it = 0; it < lc.nLinesV; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackV[it]; ih++) {
      TMarker marker(lc.TrackHitPosV[it][ih][0] * CM,
                     lc.TrackHitPosV[it][ih][1] * CM + track2d_offset, 25);
      marker.SetMarkerColor(color_2d);
      marker.SetMarkerSize(0.75);
      markers_x.push_back(marker);
    }
  }
  for (int it = 0; it < lc.nLinesX; it++) {
    for (int ih = 0; ih < lc.nHitsInTrackX[it]; ih++) {
      TMarker marker(lc.TrackHitPosX[it][ih][0] * CM,
                     lc.TrackHitPosX[it][ih][1] * CM + track2d_offset, 21);
      marker.SetMarkerColor(color_2d);
      marker.SetMarkerSize(0.75);
      markers_y.push_back(marker);
    }
  }

  int n_3d_hits = 0;
  for (int it = 0; it < reco.nTracks; it++) {
    for (int ih = 0; ih < reco.nKalmanNodes[it]; ih++) {
      {
        TMarker marker(reco.KalmanPos[it][ih][2] * CM,
                       reco.KalmanPos[it][ih][0] * CM + track3d_offset, 20);
        marker.SetMarkerColor(color_kalman_pale);
        marker.SetMarkerSize(1.35);
        markers_x.push_back(marker);
      }
      {
        TMarker marker(reco.KalmanPos[it][ih][2] * CM,
                       reco.KalmanPos[it][ih][1] * CM + track3d_offset, 20);
        marker.SetMarkerColor(color_kalman_pale);
        marker.SetMarkerSize(1.35);
        markers_y.push_back(marker);
      }
    }

    n_3d_hits += reco.nHits[it];
    for (int ih = 0; ih < reco.nHits[it]; ih++) {
      {
        TMarker marker(reco.TrackHitPos[it][ih][2] * CM,
                       reco.TrackHitPos[it][ih][0] * CM + track3d_offset, 20);
        marker.SetMarkerColor(color_3d);
        marker.SetMarkerSize(0.75);
        markers_x.push_back(marker);
      }
      {
        TMarker marker(reco.TrackHitPos[it][ih][2] * CM,
                       reco.TrackHitPos[it][ih][1] * CM + track3d_offset, 20);
        marker.SetMarkerColor(color_3d);
        marker.SetMarkerSize(0.75);
        markers_y.push_back(marker);
      }
    }
    for (int ih = 0; ih < reco.nHits[it] - 1; ih++) {
      {
        TLine line(reco.TrackHitPos[it][ih][2] * CM,
                   reco.TrackHitPos[it][ih][0] * CM + track3d_offset,
                   reco.TrackHitPos[it][ih + 1][2] * CM,
                   reco.TrackHitPos[it][ih + 1][0] * CM + track3d_offset);
        line.SetLineColor(color_3d);
        line.SetLineWidth(1);
        lines_x.push_back(line);
      }
      {
        TLine line(reco.TrackHitPos[it][ih][2] * CM,
                   reco.TrackHitPos[it][ih][1] * CM + track3d_offset,
                   reco.TrackHitPos[it][ih + 1][2] * CM,
                   reco.TrackHitPos[it][ih + 1][1] * CM + track3d_offset);
        line.SetLineColor(color_3d);
        line.SetLineWidth(1);
        lines_y.push_back(line);
      }
    }
    for (int ih = 0; ih < reco.nKalmanNodes[it] - 1; ih++) {
      {
        TLine line(reco.KalmanPos[it][ih][2] * CM,
                   reco.KalmanPos[it][ih][0] * CM + track3d_offset,
                   reco.KalmanPos[it][ih + 1][2] * CM,
                   reco.KalmanPos[it][ih + 1][0] * CM + track3d_offset);
        line.SetLineColor(color_kalman);
        line.SetLineWidth(2);
        lines_x.push_back(line);
      }
      {
        TLine line(reco.KalmanPos[it][ih][2] * CM,
                   reco.KalmanPos[it][ih][1] * CM + track3d_offset,
                   reco.KalmanPos[it][ih + 1][2] * CM,
                   reco.KalmanPos[it][ih + 1][1] * CM + track3d_offset);
        line.SetLineColor(color_kalman);
        line.SetLineWidth(2);
        lines_y.push_back(line);
      }
    }
  }

  int n_2d_hits = 0;
  for (int it = 0; it < lc.nLinesU; it++)
    n_2d_hits += lc.nHitsInTrackU[it];
  for (int it = 0; it < lc.nLinesV; it++)
    n_2d_hits += lc.nHitsInTrackV[it];
  for (int it = 0; it < lc.nLinesX; it++)
    n_2d_hits += lc.nHitsInTrackX[it];
  int n_2d_lines = lc.nLinesU + lc.nLinesV + lc.nLinesX;

  auto textBox = MakeTextBox(0.01, 0.28, 0.28, 0.90, 0.045, 12);
  textBox.AddText("Reco stages");
  textBox.AddText(
      TString::Format("Entry: %d, Run: %d", reco.EventNo, truth.RunNo));
  textBox.AddText(TString::Format("Spill: %d", reco.SpillNo));
  textBox.AddText(TString::Format("Slice: %d", reco.SliceNo));
  textBox.AddText(TString::Format("All hits: %d", lc.nHits));
  textBox.AddText(TString::Format("2D hits: %d", n_2d_hits));
  textBox.AddText(TString::Format("3D hits: %d", n_3d_hits));
  textBox.AddText(TString::Format("2D lines: %d", n_2d_lines));
  textBox.AddText(TString::Format("3D tracks: %d", reco.nTracks));
  int n_chi2_lines = reco.nTracks;
  if (n_chi2_lines > 4)
    n_chi2_lines = 4;
  for (int it = 0; it < n_chi2_lines; it++) {
    int ndf = reco.nKalmanNodes[it];
    double chi2_ndf = ndf > 0 ? reco.Chi2[it] / ndf : 0;
    textBox.AddText(TString::Format("T%d #chi^{2}/ndf: %.1f/%d = %.2f", it,
                                    reco.Chi2[it], ndf, chi2_ndf));
  }
  if (reco.nTracks > n_chi2_lines)
    textBox.AddText("...");
  if (message != "") {
    for (auto foo : split(message.c_str(), ';')) {
      textBox.AddText(foo.c_str());
    }
  }

  auto labels_x = MakeTextBox(0.01, 0.06, 0.28, 0.24, 0.060, 12);
  labels_x.AddText(TString::Format("All hits +%.0f cm", raw_offset));
  labels_x.AddText("2D reco hits +0 cm");
  labels_x.AddText(TString::Format("3D reco hits %.0f cm", track3d_offset));

  TLegend leg(0.01, 0.00, 0.28, 0.27);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.052);
  AddMarkerToLegend(leg, 20, 1.5, kBlack, "All hits");
  AddMarkerToLegend(leg, 21, 1.5, color_2d, "2D reco hits");
  AddMarkerToLegend(leg, 20, 1.5, color_3d, "3D reco hits");
  AddLineToLegend(leg, 2, color_kalman, "Kalman fit");

  canvas->cd(1);
  hist_x.Draw("Axis");
  for (auto &line : lines_x)
    line.Draw();
  for (auto &marker : markers_x)
    marker.Draw();
  textBox.Draw("same");
  labels_x.Draw("same");

  canvas->cd(2);
  hist_y.Draw("Axis");
  for (auto &line : lines_y)
    line.Draw();
  for (auto &marker : markers_y)
    marker.Draw();
  leg.Draw();

  canvas->Print((directoryPath + outfilename + ".png").c_str());
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
             1900 + buffer, 100, -350 - buffer, 100 + buffer);
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
  // Draw reco hits
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

  bool draw_2d_tracks_u = true;
  bool draw_2d_tracks_v = true;
  bool draw_2d_tracks_x = true;
  int color_adjustment = -2;
  std::cout << "lc.nLinesU: " << lc.nLinesU << ", lc.nLinesV: " << lc.nLinesV
            << ", lc.nLinesX: " << lc.nLinesX
            << ", reco.nTracks: " << reco.nTracks << std::endl;
  if (draw_2d_tracks_u) {
    for (int it = 0; it < lc.nLinesU; it++) {
      auto color_to_use = track_colors[it % n_track_colors] + color_adjustment;
      for (int ih = 0; ih < lc.nHitsInTrackU[it]; ih++) {
        float mx = lc.TrackHitPosU[it][ih][0] * CM;
        float my = lc.TrackHitPosU[it][ih][1] * CM;
        ConstrainMarkerX(mx, my);
        TMarker marker(mx, my, 21);
        marker.SetMarkerStyle(21);
        marker.SetMarkerColor(color_to_use);
        markers.push_back(marker);
      }
      /*for (int ih = 0; ih < lc.nHitsInTrackU[it] - 1; ih++) {
        {
          float mx = lc.TrackHitPosU[it][ih][0];
          float my = lc.TrackHitPosU[it][ih][1];
          float mx2 = lc.TrackHitPosU[it][ih + 1][0];
          float my2 = lc.TrackHitPosU[it][ih + 1][1];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(1);
          lines.push_back(line);
        }
      }*/
    }
  }
  if (draw_2d_tracks_v) {
    for (int it = 0; it < lc.nLinesV; it++) {
      std::cout << it << " lc.nHitsInTrackV[it]: " << lc.nHitsInTrackV[it]
                << std::endl;
      auto color_to_use = track_colors[it % n_track_colors] - color_adjustment;
      for (int ih = 0; ih < lc.nHitsInTrackV[it]; ih++) {
        float mx = lc.TrackHitPosV[it][ih][0] * CM;
        float my = lc.TrackHitPosV[it][ih][1] * CM;
        ConstrainMarkerX(mx, my);
        TMarker marker(mx, my, 21);
        marker.SetMarkerStyle(21);
        marker.SetMarkerColor(color_to_use);
        markers.push_back(marker);
      }
      /*for (int ih = 0; ih < lc.nHitsInTrackV[it] - 1; ih++) {
        {
          float mx = lc.TrackHitPosV[it][ih][0];
          float my = lc.TrackHitPosV[it][ih][1];
          float mx2 = lc.TrackHitPosV[it][ih + 1][0];
          float my2 = lc.TrackHitPosV[it][ih + 1][1];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(1);
          lines.push_back(line);
        }
      }*/
    }
  }
  if (draw_2d_tracks_x) {
    for (int it = 0; it < lc.nLinesX; it++) {
      // std::cout<<"lc.nHitsInTrackX[it]: "<<lc.nHitsInTrackX[it]<<std::endl;
      auto color_to_use = track_colors[it % n_track_colors] + color_adjustment;
      for (int ih = 0; ih < lc.nHitsInTrackX[it]; ih++) {
        float mx = lc.TrackHitPosX[it][ih][0] * CM;
        float my = lc.TrackHitPosX[it][ih][1] * CM;
        ConstrainMarkerY(mx, my);
        TMarker marker(mx, my, 21);
        marker.SetMarkerStyle(21);
        marker.SetMarkerColor(color_to_use);
        markersy.push_back(marker);
      }
      /*for (int ih = 0; ih < lc.nHitsInTrackX[it] - 1; ih++) {
        {
          float mx = lc.TrackHitPosX[it][ih][0];
          float my = lc.TrackHitPosX[it][ih][1];
          float mx2 = lc.TrackHitPosX[it][ih + 1][0];
          float my2 = lc.TrackHitPosX[it][ih + 1][1];
          //std::cout<<mx<<","<<my<<","<<mx2<<","<<my2<<std::endl;
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(1);
          linesy.push_back(line);
        }
      }*/
    }
  }

  // This section draws reco tracks as lines between two points
  if (true) {
    for (int it = 0; it < reco.nTracks; it++) {
      auto color_to_use = track_colors[(it + 1) % n_track_colors] - 2;
      for (int ih = 0; ih < reco.nHits[it] - 1; ih++) {
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][0];
          float mx2 = reco.TrackHitPos[it][ih + 1][2];
          float my2 = reco.TrackHitPos[it][ih + 1][0];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(1);
          lines.push_back(line);
        }
        {
          float mx = reco.TrackHitPos[it][ih][2];
          float my = reco.TrackHitPos[it][ih][1];
          float mx2 = reco.TrackHitPos[it][ih + 1][2];
          float my2 = reco.TrackHitPos[it][ih + 1][1];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(1);
          linesy.push_back(line);
        }
      }
    }
  }

  // Same but for kalman
  if (true) {
    for (int it = 0; it < reco.nTracks; it++) {
      auto color_to_use = track_colors[it % n_track_colors];
      for (int ih = 0; ih < reco.nKalmanNodes[it] - 1; ih++) {
        {
          float mx = reco.KalmanPos[it][ih][2];
          float my = reco.KalmanPos[it][ih][0];
          float mx2 = reco.KalmanPos[it][ih + 1][2];
          float my2 = reco.KalmanPos[it][ih + 1][0];
          TLine line(mx * CM, my * CM, mx2 * CM, my2 * CM);
          line.SetLineColor(color_to_use);
          line.SetLineWidth(2);
          lines.push_back(line);
        }
        {
          float mx = reco.KalmanPos[it][ih][2];
          float my = reco.KalmanPos[it][ih][1];
          float mx2 = reco.KalmanPos[it][ih + 1][2];
          float my2 = reco.KalmanPos[it][ih + 1][1];
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
  
  
  bool should_draw_reco_endpoint = true;
  if (should_draw_reco_endpoint) {
    for (int it = 0; it < reco.nTracks; it++) {
      double reco_startpoint_x = reco.StartPos[it][0];
      double reco_startpoint_y = reco.StartPos[it][1];
      double reco_startpoint_z = reco.StartPos[it][2];
      double kalman_startpoint_x = reco.KalmanPos[it][0][0];
      double kalman_startpoint_y = reco.KalmanPos[it][0][1];
      double kalman_startpoint_z = reco.KalmanPos[it][0][2];
      double reco_endpoint_x = reco.EndPos[it][0];
      double reco_endpoint_y = reco.EndPos[it][1];
      double reco_endpoint_z = reco.EndPos[it][2];
      int last_kalman_node = reco.nKalmanNodes[it] - 1;
      double kalman_endpoint_x = reco.KalmanPos[it][last_kalman_node][0];
      double kalman_endpoint_y = reco.KalmanPos[it][last_kalman_node][1];
      double kalman_endpoint_z = reco.KalmanPos[it][last_kalman_node][2];
      
      int color_to_use = kOrange - 3;
      int marker_to_use = 22;
      {
        float mx = reco_startpoint_z;
        float my = reco_startpoint_x;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markers.push_back(marker);
      }
      {
        float mx = reco_startpoint_z;
        float my = reco_startpoint_y;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markersy.push_back(marker);
      }
      {
        float mx = reco_endpoint_z;
        float my = reco_endpoint_x;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markers.push_back(marker);
      }
      {
        float mx = reco_endpoint_z;
        float my = reco_endpoint_y;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markersy.push_back(marker);
      }
      color_to_use = kBlue;
      {
        float mx = kalman_startpoint_z;
        float my = kalman_startpoint_x;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markers.push_back(marker);
      }
      {
        float mx = kalman_startpoint_z;
        float my = kalman_startpoint_y;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markersy.push_back(marker);
      }
      {
        float mx = kalman_endpoint_z;
        float my = kalman_endpoint_x;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markers.push_back(marker);
      }
      {
        float mx = kalman_endpoint_z;
        float my = kalman_endpoint_y;
        TMarker marker(mx * CM, my * CM, marker_to_use);
        marker.SetMarkerColor(color_to_use);
        marker.SetMarkerSize(marker_size);
        markersy.push_back(marker);
      }
    }
  }
  
  
  // This section draws true particles
  bool should_draw_true_particles = true;
  if (should_draw_true_particles) {
    int colors[] = {kRed, kBlue, kOrange - 3, kMagenta - 6, kGray};
    int ncolors = sizeof(colors) / sizeof(colors[0]);
    int n_offset = 0;
    for (int ip = 0; ip < reco.nTracks; ip++) {
      // Only draw charged particles
      bool should_draw_true_particle = false;
      if (std::abs(truth.RecoTrackPrimaryParticlePDG[ip]) == 13)
        should_draw_true_particle = true;
      // Pion
      if (std::abs(truth.RecoTrackPrimaryParticlePDG[ip]) == 211)
        should_draw_true_particle = true;
      // Kaon
      if (std::abs(truth.RecoTrackPrimaryParticlePDG[ip]) == 321)
        should_draw_true_particle = true;
      // Proton
      if (std::abs(truth.RecoTrackPrimaryParticlePDG[ip]) == 2212)
        should_draw_true_particle = true;
      
      if (should_draw_true_particle) {
        int marker_to_use_start = 22;
        if (!truth.RecoTrackPrimaryParticleTMSFiducialStart[ip])
          marker_to_use_start = empty_up_triangle;
        {
          float mx = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[ip][2];
          float my = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[ip][0];
          TMarker marker(mx * CM, my * CM, marker_to_use_start);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[ip][2];
          float my = truth.RecoTrackPrimaryParticleTruePositionEnteringTMS[ip][1];
          TMarker marker(mx * CM, my * CM, marker_to_use_start);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markersy.push_back(marker);
        }
        int marker_to_use_end = empty_down_triangle;
        {
          float mx = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][2];
          float my = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][0];
          if (mx > 19000) {
            mx = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][2];
            my = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][0];
          }
          TMarker marker(mx * CM, my * CM, marker_to_use_end);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][2];
          float my = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][1];
          if (mx > 19000) {
            mx = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][2];
            my = truth.RecoTrackPrimaryParticleTrueMomentumLeavingTMS[ip][1];
          }
          TMarker marker(mx * CM, my * CM, marker_to_use_end);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markersy.push_back(marker);
        }
        int marker_to_use_death = 23;
        {
          float mx = truth.RecoTrackPrimaryParticleTruePositionEnd[ip][2];
          float my = truth.RecoTrackPrimaryParticleTruePositionEnd[ip][0];
          TMarker marker(mx * CM, my * CM, marker_to_use_death);
          marker.SetMarkerColor(colors[n_offset % ncolors]);
          marker.SetMarkerSize(marker_size);
          markers.push_back(marker);
        }
        {
          float mx = truth.RecoTrackPrimaryParticleTruePositionEnd[ip][2];
          float my = truth.RecoTrackPrimaryParticleTruePositionEnd[ip][1];
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
  y2 = 1.0;
  auto textBoxD =
      MakeTextBox(x1, y1 - 0.1, x2, y2, text_size, 12); // Left-aligned
  auto textBoxE =
      MakeTextBox(x1, y1 - 0.1, x2, y2, text_size, 32); // Right-aligned
  // Add message but only if it's not the default message
  if (message.find("n tracks =") == std::string::npos) {
    // textBoxD.AddText("Display info:");
    // textBoxE.AddText("");
    int n_lines = 0;
    for (auto foo : split(message.c_str(), ';')) {
      textBoxD.AddText(foo.c_str());
      textBoxE.AddText("");
      n_lines += 1;
    }
    // Fill out to make it a bit higher
    while (n_lines < 3) {
      textBoxD.AddText("");
      textBoxE.AddText("");
      n_lines += 1;
    }
  } else {
    textBoxD.AddText("");
    textBoxE.AddText("");
  }
  textBoxA.AddText("N reco tracks:");
  textBoxC.AddText(TString::Format("%d", reco.nTracks));
  int time_slice_start_time = reco.TimeSliceStartTime;
  int time_slice_end_time = reco.TimeSliceEndTime;
  int slice_width = time_slice_end_time - time_slice_start_time;
  textBoxA.AddText("Slice start:");
  textBoxC.AddText(TString::Format("%dns", time_slice_start_time));
  textBoxA.AddText("Slice width:");
  textBoxC.AddText(TString::Format("%dns", slice_width));

  TLegend leg(x1, 0, x2 + 0.1, y1 - 0.1);
  leg.SetFillStyle(0);
  AddMarkerToLegend(leg, 21, 1.5, kBlack, "Reco Hits");
  AddMarkerToLegend(leg, 33, 1, kGreen + 2, "True Hits");
  if (draw_reco_track_nodes) {
    AddMarkerToLegend(leg, 20, 1.5, kRed, "Reco Track Nodes");
    AddMarkerToLegend(leg, 33, 1, kGreen + 2, "True Hits on Track");
  }
  auto kalman_entry = leg.AddEntry("", "Kalman Track", "l");
  kalman_entry->SetLineColor(kRed);
  auto reco_entry = leg.AddEntry("", "3D Reco Track", "l");
  reco_entry->SetLineColor(kRed - 1);
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
