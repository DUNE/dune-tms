int PDGtoIndex(int pdgCode) {
  // const char *pdg[] = {"e^{+/-}, #gamma", "#mu^{-}", "#mu^{+}", "#pi^{+}",
  // "#pi^{-}", "K", "n", "p", "other", "unknown"}; Unknown is -999999999
  if (pdgCode < -999999990)
    return 9;
  switch (pdgCode) {
  case 11:
    return 0; // e-
  case -11:
    return 0; // e+
  case 22:
    return 0; // gamma
  case 13:
    return 1; // mu-
  case -13:
    return 2; // mu+
  case 211:
    return 3; // pi+
  case -211:
    return 4; // pi-
  case 321:
    return 5; // K+
  case -321:
    return 5; // K-
  case 310:
    return 5; // K0
  case 130:
    return 5; // K0_L
  case 311:
    return 5; // K0_S
  case 2112:
    return 6; // Neutron
  case 2212:
    return 7; // Proton
  case -2212:
    return 7; // anti-Proton
  default:
    return 8; // other
  }
}

int PDGtoIndexReduced(int pdgCode) {
  // const char *particle_types[] = {"electron", "muon", "pion", "kaon",
  // "neutron", "proton", "other", "unknown"};
  if (pdgCode < -999999990)
    return 7;
  switch (pdgCode) {
  case 11:
    return 0; // e-
  case -11:
    return 0; // e+
  case 22:
    return 0; // gamma
  case 13:
    return 1; // mu-
  case -13:
    return 1; // mu+
  case 211:
    return 2; // pi+
  case -211:
    return 2; // pi-
  case 321:
    return 3; // K+
  case -321:
    return 3; // K-
  case 310:
    return 3; // K0
  case 130:
    return 3; // K0_L
  case 311:
    return 3; // K0_S
  case 2112:
    return 4; // Neutron
  case 2212:
    return 5; // Proton
  case -2212:
    return 5; // anti-Proton
  default:
    return 6; // other
  }
}

int NuPDGtoIndex(int pdg) {
  switch (pdg) {
  case -12:
    return 1; // nue bar
  case -14:
    return 2; // numu bar
  case -16:
    return 3; // nutau bar
  case 12:
    return 4; // nue bar
  case 14:
    return 5; // numu bar
  case 16:
    return 6; // nutau bar
  default:
    return 0; // other
  }
}

std::unordered_map<std::string, TH1 *> mapForGetHist;
std::unordered_map<std::string, std::string> specialHists;

std::unordered_map<std::string, std::tuple<std::string, int, double, double>>
    registeredAxes;
inline void
RegisterAxis(std::string axis_name,
             std::tuple<std::string, int, double, double> axis_tuple) {
  // Only register the first time
  if (registeredAxes.find(axis_name) == registeredAxes.end()) {
    registeredAxes[axis_name] = axis_tuple;
  }
}

class AutoregisterAxis {
public:
  AutoregisterAxis(const std::string &name,
                   std::tuple<std::string, int, double, double> axis_tuple) {
    std::cout << "Registering axis " << name << std::endl;
    RegisterAxis(name, axis_tuple);
  }
};

#define REGISTER_AXIS(name, axis_tuple)                                        \
  static AutoregisterAxis reg_axis_##name(#name, axis_tuple)

std::tuple<std::string, int, double, double> GetBinning(std::string axis_name) {
  if (axis_name == "ntracks")
    return std::make_tuple("N Tracks", 10, -0.5, 9.5);
  if (axis_name == "n0-120")
    return std::make_tuple("N", 24, 0, 120);
  if (axis_name == "n0-500")
    return std::make_tuple("N", 25, 0, 500);
  if (axis_name == "EventNo")
    return std::make_tuple("Event Number", 100, 0, 5000);
  if (axis_name == "SliceNo")
    return std::make_tuple("Slice Number", 61, -0.5, 60.5);
  if (axis_name == "SpillNo")
    return std::make_tuple("Spill Number", 100, 0, 300);
  if (axis_name == "X")
    return std::make_tuple("X (cm)", 100, -400, 400);
  if (axis_name == "Y")
    return std::make_tuple("Y (cm)", 100, -500, 100);
  if (axis_name == "Y_full")
    return std::make_tuple("Y (cm)", 100, -500, 500);
  if (axis_name == "Z")
    return std::make_tuple("Z (cm)", 100, 1100, 1900);
  if (axis_name == "Z_full")
    return std::make_tuple("Z (cm)", 100, 0, 2300);
  if (axis_name == "T")
    return std::make_tuple("T (ns)", 100, 0, 10000);
  if (axis_name == "direction_xz")
    return std::make_tuple("XZ Direction (dx/dz)", 31, -2, 2);
  if (axis_name == "direction_yz")
    return std::make_tuple("YZ Direction (dy/dz)", 31, -2, 2);
  if (axis_name == "dx")
    return std::make_tuple("dX", 51, -1, 1);
  if (axis_name == "dy")
    return std::make_tuple("dY", 51, -1, 1);
  if (axis_name == "dz")
    return std::make_tuple("dZ", 51, -1, 1);
  if (axis_name == "pdg")
    return std::make_tuple("Particle", 10, -0.5, 9.5);
  if (axis_name == "nu_pdg")
    return std::make_tuple("Neutrino", 7, -0.5, 6.5);
  if (axis_name == "angle_tms_enter")
    return std::make_tuple("Angle (deg)", 30, -60, 60);
  if (axis_name == "unusual_hit_locations")
    return std::make_tuple("Hit Location", 4, -0.5, 3.5);
  if (axis_name == "charge")
    return std::make_tuple("Charge", 50, -100, 100);
  if (axis_name == "areal_density")
    return std::make_tuple("Areal Density (g/cm^2)", 50, 0, 3000);
  if (axis_name == "momentum")
    return std::make_tuple("Momentum (GeV)", 50, 0, 5);
  if (axis_name == "energy_range")
    return std::make_tuple("Energy Range (GeV)", 50, 0, 5);
  if (axis_name == "energy_deposit")
    return std::make_tuple("Energy Deposit (MeV)", 50, 0, 3000);
  if (axis_name == "yesno")
    return std::make_tuple("Yes or No", 2, -0.5, 1.5);
  if (axis_name == "falsetrue")
    return std::make_tuple("", 2, -0.5, 1.5);
  if (axis_name == "rock_muon_event_rate")
    return std::make_tuple("Rock Muon Event Rates", 4, -0.5, 3.5);
  if (axis_name == "tms_event_rate")
    return std::make_tuple("TMS Event Rates", 4, -0.5, 3.5);
  if (axis_name == "tms_event_rate_by_vertex")
    return std::make_tuple("TMS Event Rates", 3, -0.5, 2.5);
  if (axis_name == "tms_muon_reco_rates")
    return std::make_tuple("TMS Muon Rates", 4, -0.5, 3.5);

  // Allow for registered axis so we don't need to add them away from where
  // they're used
  if (registeredAxes.find(axis_name) != registeredAxes.end())
    return registeredAxes[axis_name];

  std::cerr << "Fatal: Add axis to GetBinning. Did not understand axis name "
            << axis_name << std::endl;
  throw std::runtime_error("Unable to understand axis name");
}

double muon_ke_bins[] = {0.0,  0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
                         2.25, 2.5,  2.75, 3.0,  3.5, 4.0,  4.5, 5.0};
int n_muon_ke_bins = sizeof(muon_ke_bins) / sizeof(double) - 1;

std::tuple<bool, std::string, int, double *>
GetComplexBinning(std::string axis_name) {
  if (axis_name == "ke_tms_enter")
    return std::make_tuple(true,
                           "True Muon KE Entering TMS (GeV);N Muons / GeV",
                           n_muon_ke_bins, muon_ke_bins);
  if (axis_name == "ke_tms_enter_true")
    return std::make_tuple(true, "True Muon KE Entering TMS (GeV)",
                           n_muon_ke_bins, muon_ke_bins);
  if (axis_name == "ke_tms_enter_reco")
    return std::make_tuple(true, "Reco Muon KE Entering TMS (GeV)",
                           n_muon_ke_bins, muon_ke_bins);
  return std::make_tuple(false, "", 0, (double *)NULL);
}

void AdjustAxis(TH1 *hist, std::string xaxis, std::string yaxis = "",
                std::string zaxis = "") {
  (void)zaxis;
  (void)yaxis;
  if (xaxis == "pdg") {
    const char *pdg[] = {"e^{+/-}, #gamma", "#mu^{-}", "#mu^{+}", "#pi^{+}",
                         "#pi^{-}",         "K",       "n",       "p",
                         "other",           "unknown"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    hist->SetNdivisions(npdg);
    for (int ib = 0; ib < npdg; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, pdg[ib]);
    }
  }
  if (xaxis == "nu_pdg") {
    const char *pdg[] = {
        "other",   "#bar{#nu_{e}}", "#bar{#nu_{#mu}}", "#bar{#nu_{#tau}}",
        "#nu_{e}", "#nu_{#mu}",     "#nu_{#tau}"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    hist->SetNdivisions(npdg);
    for (int ib = 0; ib < npdg; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, pdg[ib]);
    }
  }

  if (xaxis == "yesno") {
    const char *pdg[] = {"yes", "no"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    hist->SetNdivisions(npdg);
    for (int ib = 0; ib < npdg; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, pdg[ib]);
    }
  }

  if (xaxis == "falsetrue") {
    const char *pdg[] = {"false", "true"};
    const int npdg = sizeof(pdg) / sizeof(pdg[0]);
    hist->SetNdivisions(npdg);
    for (int ib = 0; ib < npdg; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, pdg[ib]);
    }
  }

  if (xaxis == "unusual_hit_locations") {
    const char *labels[] = {"Hit at Zero", "Hit at -999", "Hit Below Fiducial",
                            "Hit Above Fiducial"};
    const int nlabels = sizeof(labels) / sizeof(labels[0]);
    hist->SetNdivisions(nlabels);
    for (int ib = 0; ib < nlabels; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, labels[ib]);
    }
  }

  if (xaxis == "rock_muon_event_rate") {
    const char *labels[] = {"Rock Muons", "Ending in ND-LAr Active",
                            "Ending in TMS Active", "Ending Other"};
    const int nlabels = sizeof(labels) / sizeof(labels[0]);
    hist->SetNdivisions(nlabels);
    for (int ib = 0; ib < nlabels; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, labels[ib]);
    }
  }

  if (xaxis == "tms_event_rate") {
    const char *labels[] = {"N Total", "Touch", "End", "Exit"};
    const int nlabels = sizeof(labels) / sizeof(labels[0]);
    hist->SetNdivisions(nlabels);
    for (int ib = 0; ib < nlabels; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, labels[ib]);
    }
  }

  if (xaxis == "tms_event_rate_by_vertex") {
    const char *labels[] = {"N Total", "Any Energy", "5+ MeV Particle(s)"};
    const int nlabels = sizeof(labels) / sizeof(labels[0]);
    hist->SetNdivisions(nlabels);
    for (int ib = 0; ib < nlabels; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, labels[ib]);
    }
  }

  if (xaxis == "tms_muon_reco_rates") {
    const char *labels[] = {"N Total", "Touch", "End", "Exit"};
    const int nlabels = sizeof(labels) / sizeof(labels[0]);
    hist->SetNdivisions(nlabels);
    for (int ib = 0; ib < nlabels; ib++) {
      hist->GetXaxis()->ChangeLabel(ib + 1, -1, -1, -1, -1, -1, labels[ib]);
    }
  }
}

TH1 *MakeHist(std::string directory_and_name, std::string title,
              std::string xaxis, std::string yaxis = "",
              std::string zaxis = "") {
  if (zaxis != "") {
    // 3d hist case
    throw std::runtime_error("3d hists are not implemented yet");
  } else if (yaxis != "" && yaxis != "Probability" &&
             yaxis != "Normalized Frequency" && yaxis.find("#") != 0) {
    // 2d hist case
    TH2D *out;

    std::string xaxis_title;
    int xaxis_nbins;
    double *xaxis_bins;
    bool has_complex_binning_x;
    std::tie(has_complex_binning_x, xaxis_title, xaxis_nbins, xaxis_bins) =
        GetComplexBinning(xaxis);

    std::string yaxis_title;
    int yaxis_nbins;
    double *yaxis_bins;
    bool has_complex_binning_y;
    std::tie(has_complex_binning_y, yaxis_title, yaxis_nbins, yaxis_bins) =
        GetComplexBinning(yaxis);
    // Can't do one and not the other, but could manually create the binning if
    // we wanted
    if (has_complex_binning_y != has_complex_binning_x)
      throw std::runtime_error(
          "2d hists have to use complex y and x bins simultaneously");
    if (has_complex_binning_x) {
      auto complete_title = title + ";" + xaxis_title + ";" + yaxis_title;
      out = new TH2D(directory_and_name.c_str(), complete_title.c_str(),
                     xaxis_nbins, xaxis_bins, yaxis_nbins, yaxis_bins);
    } else {
      double xaxis_start;
      double xaxis_end;
      std::tie(xaxis_title, xaxis_nbins, xaxis_start, xaxis_end) =
          GetBinning(xaxis);

      double yaxis_start;
      double yaxis_end;
      std::tie(yaxis_title, yaxis_nbins, yaxis_start, yaxis_end) =
          GetBinning(yaxis);
      auto complete_title = title + ";" + xaxis_title + ";" + yaxis_title;
      out = new TH2D(directory_and_name.c_str(), complete_title.c_str(),
                     xaxis_nbins, xaxis_start, xaxis_end, yaxis_nbins,
                     yaxis_start, yaxis_end);
    }
    // Add special naming here
    AdjustAxis(out, xaxis, yaxis);
    return out;
  } else {
    // 1d hist case
    TH1D *out;
    std::string xaxis_title;
    int xaxis_nbins;
    double *xaxis_bins;
    bool has_complex_binning;
    std::tie(has_complex_binning, xaxis_title, xaxis_nbins, xaxis_bins) =
        GetComplexBinning(xaxis);
    if (has_complex_binning) {

      auto complete_title = title + ";" + xaxis_title;
      out = new TH1D(directory_and_name.c_str(), complete_title.c_str(),
                     xaxis_nbins, xaxis_bins);
    } else {
      double xaxis_start;
      double xaxis_end;
      std::tie(xaxis_title, xaxis_nbins, xaxis_start, xaxis_end) =
          GetBinning(xaxis);
      auto complete_title = title + ";" + xaxis_title;
      out = new TH1D(directory_and_name.c_str(), complete_title.c_str(),
                     xaxis_nbins, xaxis_start, xaxis_end);
    }
    // Add special naming here
    AdjustAxis(out, xaxis);
    if (yaxis.find("#") == 0) yaxis = yaxis.substr(1, yaxis.size()-1);
    out->GetYaxis()->SetTitle(yaxis.c_str());
    return out;
  }
}

void NormalizeColumns(TH2 *hist, bool colMax = false) {
  if (!hist) {
    std::cerr << "Error: Null histogram passed to NormalizeColumns."
              << std::endl;
    return;
  }

  int nBinsX = hist->GetNbinsX();
  int nBinsY = hist->GetNbinsY();

  // Iterate over each column (x-bin)
  for (int binX = 1; binX <= nBinsX; ++binX) {
    double columnSum = 0.0;
    double columnMax = 0.0;

    // Calculate the sum of the column
    for (int binY = 1; binY <= nBinsY; ++binY) {
      columnSum += hist->GetBinContent(binX, binY);
      if (columnMax < hist->GetBinContent(binX, binY))
        columnMax = hist->GetBinContent(binX, binY);
    }

    double weight = 1.0;
    if (colMax && columnMax > 0)
      weight = 1 / columnMax;
    if (!colMax && columnSum > 0)
      weight = 1 / columnSum;
    for (int binY = 1; binY <= nBinsY; ++binY) {
      double value = hist->GetBinContent(binX, binY);
      // Only set the value if nonzero
      if (value > 0.001)
        hist->SetBinContent(binX, binY, value * weight);
    }
  }
}

void NormalizeRows(TH2 *hist) {
  if (!hist) {
    std::cerr << "Error: Null histogram passed to NormalizeColumns."
              << std::endl;
    return;
  }

  int nBinsX = hist->GetNbinsX();
  int nBinsY = hist->GetNbinsY();

  // Iterate over each column (x-bin)
  for (int binY = 1; binY <= nBinsY; ++binY) {
    double columnSum = 0.0;

    // Calculate the sum of the column
    for (int binX = 1; binX <= nBinsX; ++binX) {
      columnSum += hist->GetBinContent(binX, binY);
    }

    // Normalize the column if the sum is not zero
    if (columnSum > 0) {
      for (int binX = 1; binX <= nBinsX; ++binX) {
        double value = hist->GetBinContent(binX, binY);
        // Only set the value if nonzero
        if (value > 0.001)
          hist->SetBinContent(binX, binY, value / columnSum);
      }
    }
  }
}

void FinalizeHists() {
  // Could check here for certain hists and automatically make copies with
  // column norm for example
  for (auto &hist : mapForGetHist) {
    if (hist.first.find("column_normalize") != std::string::npos) {
      // This hist wants column normalization
      NormalizeColumns((TH2 *)hist.second);
    }
    if (hist.first.find("column_maximize") != std::string::npos) {
      // This hist wants column normalization such that the max is one
      NormalizeColumns((TH2 *)hist.second, true);
    }
    if (hist.first.find("row_normalize") != std::string::npos) {
      // This hist wants column normalization
      NormalizeRows((TH2 *)hist.second);
    }
    if (hist.first.find("normalize_to_one") != std::string::npos) {
      // Normalize max to 1
      hist.second->Scale(1 / hist.second->GetMaximum());
    }
    if (hist.first.find("normalize_to_probability") != std::string::npos) {
      // Normalize to 1, as in a probability for example
      hist.second->Scale(1 / hist.second->Integral());
    }
    if (hist.first.find("normalize_to_one") != std::string::npos) {
      // Normalize max to 1
      hist.second->Scale(1 / hist.second->GetMaximum());
    }
  }
  
  // Now save copies of all the special hists
  for (auto& is : specialHists) {
    auto special_name = is.first;
    auto regular_name = is.second;
    // Save a copy under the new special name
    mapForGetHist[special_name] = dynamic_cast<TH1*>(mapForGetHist[regular_name]->Clone(special_name.c_str()));
  }
}

double total_lookup_time = 0;
double total_no_make_time = 0;
double total_make_time = 0;
int n_lookups = 0;
int n_makes = 0;
int n_no_makes = 0;

TH1 *GetHist(std::string directory_and_name, std::string title,
             std::string xaxis, std::string yaxis = "",
             std::string zaxis = "") {
  auto time_start = std::chrono::high_resolution_clock::now();
  bool did_make = false;
  if (mapForGetHist.find(directory_and_name) == mapForGetHist.end()) {
    // Object doesn't exist, create it and store it
    mapForGetHist[directory_and_name] =
        MakeHist(directory_and_name, title, xaxis, yaxis, zaxis);
    auto time_stop_make_only = std::chrono::high_resolution_clock::now();
    auto duration_make_only =
        std::chrono::duration_cast<std::chrono::microseconds>(
            time_stop_make_only - time_start)
            .count();
    total_make_time += duration_make_only;
    n_makes += 1;
    did_make = true;
  }
  auto out = mapForGetHist[directory_and_name];
  auto time_stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                      time_stop - time_start)
                      .count();
  total_lookup_time += duration;
  n_lookups += 1;
  if (!did_make) {
    total_no_make_time += duration;
    n_no_makes += 1;
  }
#ifdef PRINT_STATS
  if (n_lookups % 1000000 == 0) {
    double avg_lookup_time = total_lookup_time / n_lookups;
    double avg_make_time = total_make_time / n_makes;
    double avg_no_make_time = total_no_make_time / n_no_makes;

    std::cout << "Avg lookup time: " << avg_lookup_time << " ("
              << total_lookup_time << "/" << n_lookups << ")" << std::endl;
    std::cout << "Avg make time: " << avg_make_time << " (" << total_make_time
              << "/" << n_makes << ")" << std::endl;
    std::cout << "Avg no make time: " << avg_no_make_time << " ("
              << total_no_make_time << "/" << n_no_makes << ")" << std::endl;
  }
#endif // end ifdef PRINT_STATS
  return out;
}

TH1 *GetSpecialHist(std::string special_dir_and_name,
             std::string directory_and_name, std::string title,
             std::string xaxis, std::string yaxis = "",
             std::string zaxis = "") {
  TH1* out = GetHist(directory_and_name, title, xaxis, yaxis, zaxis);
  // Save a copy
  specialHists[special_dir_and_name] = directory_and_name;
  return out;
}
