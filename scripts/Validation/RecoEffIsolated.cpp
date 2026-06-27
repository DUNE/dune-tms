#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLine.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "Reco_Tree.h"
#include "Truth_Info.h"

namespace {

constexpr double kGeV = 1e-3;
constexpr double kMinVisibleEnergy = 5.0; // MeV
constexpr double kDefaultLike = -1e8;

double muon_ke_bins[] = {0.0,  0.25, 0.5,  0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
                         2.25, 2.5,  2.75, 3.0,  3.5, 4.0,  4.5, 5.0};
constexpr int n_muon_ke_bins = sizeof(muon_ke_bins) / sizeof(double) - 1;

struct Counters {
  long long truth_spill_entries = 0;
  long long truth_info_entries = 0;
  long long reco_entries = 0;
  long long truth_spill_particles = 0;
  long long denom_muons = 0;
  long long denom_tms_touch = 0;
  long long denom_tms_touch_lar_start = 0;
  long long denom_tms_touch_tms_end = 0;
  long long denom_tms_touch_lar_only = 0;
  long long denom_tms_touch_tms_only = 0;
  long long denom_tms_touch_neither_strict_flag = 0;
  long long denom_strict = 0;
  long long denom_negative_ke = 0;
  long long denom_default_ke = 0;
  long long denom_out_of_bounds_ke = 0;
  long long truth_info_reco_tracks = 0;
  long long reco_tree_tracks = 0;
  long long reco_track_count_mismatch = 0;
  long long numerator_muons = 0;
  long long numerator_tms_touch = 0;
  long long numerator_all_unique = 0;
  long long numerator_unique_lar_start = 0;
  long long numerator_unique_tms_end = 0;
  long long numerator_unique_lar_only = 0;
  long long numerator_unique_tms_only = 0;
  long long numerator_unique_neither_strict_flag = 0;
  long long numerator_strict_unique = 0;
  long long numerator_duplicates = 0;
  long long numerator_negative_ke = 0;
  long long numerator_default_ke = 0;
  long long numerator_out_of_bounds_ke = 0;
  long long numerator_reco_trackn_over_max = 0;
};

Long64_t FindMaxTruthSpillParticles(TChain *chain) {
  TTreeReader reader(chain);
  TTreeReaderValue<Int_t> n_true_particles(reader, "nTrueParticles");
  Long64_t max_particles = 0;
  while (reader.Next()) {
    max_particles = std::max(max_particles, static_cast<Long64_t>(*n_true_particles));
  }
  return std::max<Long64_t>(max_particles, 1);
}

struct TruthSpill {
  TChain *chain = nullptr;
  Int_t SpillNo = 0;
  Int_t nTrueParticles = 0;
  Long64_t max_particles = 0;
  std::vector<Int_t> PDG;
  std::vector<Float_t> TrueVisibleEnergy;
  std::unique_ptr<Bool_t[]> LArFiducialStart;
  std::unique_ptr<Bool_t[]> TMSFiducialEnd;
  std::vector<Float_t> MomentumTMSStart;

  TruthSpill(TChain *input, Long64_t max_particles_in)
      : chain(input), max_particles(max_particles_in), PDG(max_particles),
        TrueVisibleEnergy(max_particles),
        LArFiducialStart(std::make_unique<Bool_t[]>(max_particles)),
        TMSFiducialEnd(std::make_unique<Bool_t[]>(max_particles)),
        MomentumTMSStart(max_particles * 4) {
    chain->SetMakeClass(1);
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("SpillNo", 1);
    chain->SetBranchStatus("nTrueParticles", 1);
    chain->SetBranchStatus("PDG", 1);
    chain->SetBranchStatus("TrueVisibleEnergy", 1);
    chain->SetBranchStatus("LArFiducialStart", 1);
    chain->SetBranchStatus("TMSFiducialEnd", 1);
    chain->SetBranchStatus("MomentumTMSStart", 1);
    chain->SetBranchAddress("SpillNo", &SpillNo);
    chain->SetBranchAddress("nTrueParticles", &nTrueParticles);
    chain->SetBranchAddress("PDG", PDG.data());
    chain->SetBranchAddress("TrueVisibleEnergy", TrueVisibleEnergy.data());
    chain->SetBranchAddress("LArFiducialStart", LArFiducialStart.get());
    chain->SetBranchAddress("TMSFiducialEnd", TMSFiducialEnd.get());
    chain->SetBranchAddress("MomentumTMSStart", MomentumTMSStart.data());
  }

  Long64_t GetEntries() const { return chain->GetEntries(); }
  void GetEntry(Long64_t entry) { chain->GetEntry(entry); }
  Float_t MomentumTMSStartE(int particle) const {
    return MomentumTMSStart[particle * 4 + 3];
  }
};

bool HasTree(const std::string &filename, const std::string &tree_name) {
  TFile file(filename.c_str(), "READ");
  if (file.IsZombie()) return false;
  return file.Get(tree_name.c_str()) != nullptr;
}

bool EndsWithRoot(const std::filesystem::path &path) {
  return path.extension() == ".root";
}

std::vector<std::string> FindInputFiles(const std::string &input) {
  std::vector<std::string> out;
  std::filesystem::path path(input);
  if (std::filesystem::is_regular_file(path)) {
    if (EndsWithRoot(path) && HasTree(input, "Truth_Spill") &&
        HasTree(input, "Truth_Info")) {
      out.push_back(input);
    }
  } else if (std::filesystem::is_directory(path)) {
    for (const auto &entry : std::filesystem::recursive_directory_iterator(path)) {
      if (!entry.is_regular_file()) continue;
      if (!EndsWithRoot(entry.path())) continue;
      const std::string filename = entry.path().string();
      if (HasTree(filename, "Truth_Spill") && HasTree(filename, "Truth_Info")) {
        out.push_back(filename);
      }
    }
  } else {
    std::cerr << "Input is neither a ROOT file nor a directory: " << input << std::endl;
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::string DefaultOutputDir() {
  const char *user = std::getenv("USER");
  std::string base = "/exp/dune/data/users/";
  base += user ? user : "unknown";
  base += "/dune-tms/Validation/RecoEffIsolated";
  return base;
}

void EnsureDir(const std::string &path) {
  std::filesystem::create_directories(path);
}

bool InKeBounds(double value) {
  return value >= muon_ke_bins[0] && value < muon_ke_bins[n_muon_ke_bins];
}

bool DefaultLike(double value) {
  return value <= kDefaultLike;
}

long long Quantize(float value) {
  return static_cast<long long>(std::llround(value * 1000.0));
}

std::string RecoTruthFingerprint(const Truth_Info &truth, int track) {
  std::ostringstream key;
  key << truth.RecoTrackPrimaryParticlePDG[track];
  for (int dim = 0; dim < 4; ++dim)
    key << ":" << Quantize(truth.RecoTrackPrimaryParticleTrueMomentum[track][dim]);
  for (int dim = 0; dim < 4; ++dim)
    key << ":" << Quantize(truth.RecoTrackPrimaryParticleTruePositionStart[track][dim]);
  for (int dim = 0; dim < 4; ++dim)
    key << ":" << Quantize(truth.RecoTrackPrimaryParticleTruePositionEnd[track][dim]);
  return key.str();
}

TH1D *MakeKeHist(const std::string &name, const std::string &title) {
  TH1D *hist = new TH1D(name.c_str(), title.c_str(), n_muon_ke_bins, muon_ke_bins);
  hist->GetXaxis()->SetTitle("True muon KE entering TMS (GeV)");
  hist->GetYaxis()->SetTitle("N muons / GeV");
  hist->Sumw2();
  return hist;
}

TH1D *MakeBoolHist(const std::string &name, const std::string &title) {
  TH1D *hist = new TH1D(name.c_str(), title.c_str(), 2, -0.5, 1.5);
  hist->GetXaxis()->SetBinLabel(1, "false");
  hist->GetXaxis()->SetBinLabel(2, "true");
  return hist;
}

TH1D *MakeCutflow(const std::string &name, const std::string &title,
                  const std::vector<std::string> &labels) {
  TH1D *hist = new TH1D(name.c_str(), title.c_str(), labels.size(), 0.5,
                        labels.size() + 0.5);
  for (size_t i = 0; i < labels.size(); ++i) {
    hist->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  hist->LabelsOption("v", "X");
  return hist;
}

void FillCut(TH1D *hist, const std::string &label) {
  hist->Fill(label.c_str(), 1.0);
}

void FillStrictReason(TH1D *hist, bool lar_start, bool tms_end) {
  if (lar_start && tms_end) {
    FillCut(hist, "pass");
  } else if (!lar_start && tms_end) {
    FillCut(hist, "no LAr start");
  } else if (lar_start && !tms_end) {
    FillCut(hist, "no TMS end");
  } else {
    FillCut(hist, "no LAr start or TMS end");
  }
}

void SaveCanvas(TCanvas &canvas, const std::string &path) {
  canvas.Print(path.c_str());
}

void DrawAndSave(TH1 *hist, TCanvas &canvas, const std::string &dir,
                 const std::string &name, const std::string &option = "") {
  canvas.Clear();
  hist->Draw(option.c_str());
  SaveCanvas(canvas, dir + "/" + name + ".png");
}

TH1D *MakeEfficiency(const std::string &name, const std::string &title,
                     TH1D *num, TH1D *den) {
  TH1D *eff = dynamic_cast<TH1D *>(num->Clone(name.c_str()));
  eff->SetTitle(title.c_str());
  eff->Divide(num, den, 1.0, 1.0, "B");
  eff->GetYaxis()->SetTitle("Reconstruction efficiency");
  eff->GetYaxis()->SetRangeUser(0.0, 1.2);
  return eff;
}

void DrawEfficiency(TH1D *eff, TCanvas &canvas, const std::string &dir,
                    const std::string &name) {
  canvas.Clear();
  eff->Draw();
  TLine line(eff->GetXaxis()->GetXmin(), 1.0, eff->GetXaxis()->GetXmax(), 1.0);
  line.SetLineStyle(2);
  line.Draw();
  SaveCanvas(canvas, dir + "/" + name + ".png");
}

void AddExample(std::vector<std::string> &examples, const std::string &text) {
  if (examples.size() < 50) examples.push_back(text);
}

void WriteSummary(const std::string &filename, const Counters &counts,
                  const std::vector<std::string> &examples,
                  const std::vector<std::string> &input_files) {
  std::ofstream out(filename);
  out << "RecoEffIsolated diagnostics\n";
  out << "Input files: " << input_files.size() << "\n";
  for (const auto &file : input_files) out << "  " << file << "\n";
  out << "\nEntries\n";
  out << "  Truth_Spill entries: " << counts.truth_spill_entries << "\n";
  out << "  Truth_Info entries:  " << counts.truth_info_entries << "\n";
  out << "  Reco_Tree entries:   " << counts.reco_entries << "\n";
  out << "\nDenominator from Truth_Spill\n";
  out << "  True particles:      " << counts.truth_spill_particles << "\n";
  out << "  Muons:               " << counts.denom_muons << "\n";
  out << "  Muons TMS-touching:  " << counts.denom_tms_touch << "\n";
  out << "  TMS-touching with LAr start: " << counts.denom_tms_touch_lar_start << "\n";
  out << "  TMS-touching with TMS end:   " << counts.denom_tms_touch_tms_end << "\n";
  out << "  TMS-touching LAr only:       " << counts.denom_tms_touch_lar_only << "\n";
  out << "  TMS-touching TMS only:       " << counts.denom_tms_touch_tms_only << "\n";
  out << "  TMS-touching neither flag:   " << counts.denom_tms_touch_neither_strict_flag << "\n";
  out << "  Strict muons:        " << counts.denom_strict << "\n";
  out << "  Negative KE:         " << counts.denom_negative_ke << "\n";
  out << "  Default-like KE:     " << counts.denom_default_ke << "\n";
  out << "  Out-of-bounds KE:    " << counts.denom_out_of_bounds_ke << "\n";
  out << "\nNumerator from Truth_Info reco-track summaries\n";
  out << "  Truth_Info reco tracks: " << counts.truth_info_reco_tracks << "\n";
  out << "  Reco_Tree tracks:       " << counts.reco_tree_tracks << "\n";
  out << "  RecoTrackN mismatches:  " << counts.reco_track_count_mismatch << "\n";
  out << "  Reco primary muons:     " << counts.numerator_muons << "\n";
  out << "  TMS-touching muons:     " << counts.numerator_tms_touch << "\n";
  out << "  Unique all numerators:  " << counts.numerator_all_unique << "\n";
  out << "  Unique with LAr start:  " << counts.numerator_unique_lar_start << "\n";
  out << "  Unique with TMS end:    " << counts.numerator_unique_tms_end << "\n";
  out << "  Unique LAr only:        " << counts.numerator_unique_lar_only << "\n";
  out << "  Unique TMS only:        " << counts.numerator_unique_tms_only << "\n";
  out << "  Unique neither flag:    " << counts.numerator_unique_neither_strict_flag << "\n";
  out << "  Unique strict numerators: " << counts.numerator_strict_unique << "\n";
  out << "  Duplicate reco muons:   " << counts.numerator_duplicates << "\n";
  out << "  Negative KE:            " << counts.numerator_negative_ke << "\n";
  out << "  Default-like KE:        " << counts.numerator_default_ke << "\n";
  out << "  Out-of-bounds KE:       " << counts.numerator_out_of_bounds_ke << "\n";
  out << "  RecoTrackN over max:    " << counts.numerator_reco_trackn_over_max << "\n";
  out << "\nExamples\n";
  for (const auto &example : examples) out << "  " << example << "\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <input.root|input_directory> [output_dir]\n";
    return 1;
  }

  const std::string input = argv[1];
  const std::string output_dir = argc >= 3 ? argv[2] : DefaultOutputDir();
  const std::string raw_dir = output_dir + "/raw_hists";
  const std::string diag_dir = output_dir + "/diagnostics";
  EnsureDir(output_dir);
  EnsureDir(raw_dir);
  EnsureDir(diag_dir);

  std::vector<std::string> input_files = FindInputFiles(input);
  if (input_files.empty()) {
    std::cerr << "No ROOT files with Truth_Spill and Truth_Info found in " << input
              << std::endl;
    return 2;
  }

  TChain truth_spill_chain("Truth_Spill");
  TChain truth_info_chain("Truth_Info");
  TChain reco_chain("Reco_Tree");
  for (const auto &file : input_files) {
    truth_spill_chain.Add(file.c_str());
    truth_info_chain.Add(file.c_str());
    if (HasTree(file, "Reco_Tree")) reco_chain.Add(file.c_str());
  }

  auto truth = std::make_unique<Truth_Info>(&truth_info_chain);
  auto reco = std::make_unique<Reco_Tree>(&reco_chain);

  Counters counts;
  std::vector<std::string> examples;
  std::map<std::string, int> reconstructed_fingerprints;

  TH1D *h_all_den = MakeKeHist("all_muon_ke_tms_enter_denominator",
                               "All muons denominator");
  TH1D *h_all_num = MakeKeHist("all_muon_ke_tms_enter_numerator",
                               "All muons numerator");
  TH1D *h_strict_den = MakeKeHist("muon_ke_tms_enter_denominator",
                                  "LAr-start/TMS-end muons denominator");
  TH1D *h_strict_num = MakeKeHist("muon_ke_tms_enter_numerator",
                                  "LAr-start/TMS-end muons numerator");

  TH1D *h_den_visible = new TH1D("den_muon_true_visible_energy",
                                 "Denominator muon true visible energy;True visible energy (MeV);N",
                                 100, -50, 2000);
  TH1D *h_num_visible = new TH1D("num_muon_true_visible_energy",
                                 "Numerator muon true visible energy;True visible energy (MeV);N",
                                 100, -50, 2000);
  TH1D *h_den_ke = new TH1D("den_muon_ke_tms_enter_diagnostic",
                            "Denominator muon KE entering TMS;KE (GeV);N", 120,
                            -1, 11);
  TH1D *h_num_ke = new TH1D("num_muon_ke_tms_enter_diagnostic",
                            "Numerator muon KE entering TMS;KE (GeV);N", 120,
                            -1, 11);
  TH1D *h_den_lar = MakeBoolHist("den_lar_fiducial_start",
                                 "Denominator muon LAr fiducial start");
  TH1D *h_den_tms = MakeBoolHist("den_tms_fiducial_end",
                                 "Denominator muon TMS fiducial end");
  TH1D *h_num_lar = MakeBoolHist("num_lar_fiducial_start",
                                 "Numerator muon LAr fiducial start");
  TH1D *h_num_tms = MakeBoolHist("num_tms_fiducial_end",
                                 "Numerator muon TMS fiducial end");
  TH1D *h_reco_trackn_compare = new TH1D(
      "truth_info_recotrackn_minus_reco_tree_ntracks",
      "Truth_Info RecoTrackN - Reco_Tree nTracks;difference;N entries", 21, -10.5,
      10.5);
  TH1D *h_duplicate = MakeBoolHist("num_duplicate_fingerprint",
                                   "Numerator duplicate fingerprint");
  TH1D *h_den_passfail = MakeCutflow(
      "denominator_cut_pass_fail", "Truth_Spill denominator cut pass/fail",
      {"muon pdg pass", "muon pdg fail", "TMS visible pass", "TMS visible fail",
       "KE bounds pass", "KE bounds fail", "LAr start pass", "LAr start fail",
       "TMS end pass", "TMS end fail"});
  TH1D *h_num_passfail = MakeCutflow(
      "numerator_cut_pass_fail", "Truth_Info numerator cut pass/fail",
      {"muon pdg pass", "muon pdg fail", "TMS visible pass", "TMS visible fail",
       "KE bounds pass", "KE bounds fail", "duplicate pass", "duplicate fail",
       "LAr start pass", "LAr start fail", "TMS end pass", "TMS end fail"});
  TH1D *h_den_strict_reason = MakeCutflow(
      "denominator_strict_cut_reason",
      "Truth_Spill strict denominator cut reason",
      {"pass", "no LAr start", "no TMS end", "no LAr start or TMS end"});
  TH1D *h_num_strict_reason = MakeCutflow(
      "numerator_strict_cut_reason",
      "Truth_Info strict numerator cut reason",
      {"pass", "no LAr start", "no TMS end", "no LAr start or TMS end"});
  TH1D *h_den_cutflow = MakeCutflow(
      "denominator_cutflow", "Truth_Spill denominator cutflow",
      {"all particles", "muon pdg", "TMS visible >= 5 MeV", "KE in hist bounds",
       "LAr fid start", "TMS fid end", "strict denominator"});
  TH1D *h_num_cutflow = MakeCutflow(
      "numerator_cutflow", "Truth_Info numerator cutflow",
      {"all reco summaries", "muon pdg", "TMS visible >= 5 MeV",
       "KE in hist bounds", "not duplicate", "LAr fid start", "TMS fid end",
       "strict numerator"});

  const Long64_t max_truth_spill_particles =
      FindMaxTruthSpillParticles(&truth_spill_chain);
  auto spill = std::make_unique<TruthSpill>(&truth_spill_chain,
                                            max_truth_spill_particles);

  counts.truth_spill_entries = spill->GetEntries();
  for (Long64_t spill_entry = 0; spill_entry < spill->GetEntries(); ++spill_entry) {
    spill->GetEntry(spill_entry);
    if (spill->nTrueParticles > spill->max_particles) {
      AddExample(examples, "Truth_Spill entry " + std::to_string(spill_entry) +
                               " SpillNo=" + std::to_string(spill->SpillNo) +
                               " reports nTrueParticles=" +
                               std::to_string(spill->nTrueParticles) +
                               " beyond allocated max=" +
                               std::to_string(spill->max_particles));
    }
    const int n_particles = static_cast<int>(
        std::max<Long64_t>(0, std::min<Long64_t>(spill->nTrueParticles,
                                                spill->max_particles)));
    counts.truth_spill_particles += n_particles;
    for (int ip = 0; ip < n_particles; ++ip) {
      FillCut(h_den_cutflow, "all particles");
      const bool is_muon = std::abs(spill->PDG[ip]) == 13;
      FillCut(h_den_passfail, is_muon ? "muon pdg pass" : "muon pdg fail");
      if (!is_muon) continue;
      FillCut(h_den_cutflow, "muon pdg");
      counts.denom_muons++;

      const double visible = spill->TrueVisibleEnergy[ip];
      const bool tms_touch = visible >= kMinVisibleEnergy;
      const bool lar_start = spill->LArFiducialStart[ip];
      const bool tms_end = spill->TMSFiducialEnd[ip];
      const double ke_tms_enter = spill->MomentumTMSStartE(ip) * kGeV;
      const bool ke_in_bounds = InKeBounds(ke_tms_enter);

      h_den_visible->Fill(visible);
      h_den_ke->Fill(ke_tms_enter);
      h_den_lar->Fill(lar_start ? 1 : 0);
      h_den_tms->Fill(tms_end ? 1 : 0);
      FillCut(h_den_passfail, tms_touch ? "TMS visible pass" : "TMS visible fail");
      FillCut(h_den_passfail, ke_in_bounds ? "KE bounds pass" : "KE bounds fail");
      FillCut(h_den_passfail, lar_start ? "LAr start pass" : "LAr start fail");
      FillCut(h_den_passfail, tms_end ? "TMS end pass" : "TMS end fail");

      if (ke_tms_enter < 0) {
        counts.denom_negative_ke++;
        AddExample(examples, "Denom negative KE: spill entry=" +
                                 std::to_string(spill_entry) + " particle=" +
                                 std::to_string(ip) + " KE=" +
                                 std::to_string(ke_tms_enter));
      }
      if (DefaultLike(spill->MomentumTMSStartE(ip)))
        counts.denom_default_ke++;
      if (!ke_in_bounds) counts.denom_out_of_bounds_ke++;

      if (!tms_touch) continue;
      FillCut(h_den_cutflow, "TMS visible >= 5 MeV");
      counts.denom_tms_touch++;
      if (ke_in_bounds) FillCut(h_den_cutflow, "KE in hist bounds");
      h_all_den->Fill(ke_tms_enter);

      FillStrictReason(h_den_strict_reason, lar_start, tms_end);
      if (lar_start) counts.denom_tms_touch_lar_start++;
      if (tms_end) counts.denom_tms_touch_tms_end++;
      if (lar_start && !tms_end) counts.denom_tms_touch_lar_only++;
      if (!lar_start && tms_end) counts.denom_tms_touch_tms_only++;
      if (!lar_start && !tms_end) counts.denom_tms_touch_neither_strict_flag++;
      if (!(lar_start && tms_end)) {
        AddExample(examples,
                   "Denom strict fail: spill entry=" + std::to_string(spill_entry) +
                       " SpillNo=" + std::to_string(spill->SpillNo) +
                       " particle=" + std::to_string(ip) +
                       " KE=" + std::to_string(ke_tms_enter) +
                       " visible=" + std::to_string(visible) +
                       " LArStart=" + std::to_string(lar_start) +
                       " TMSEnd=" + std::to_string(tms_end));
      }

      if (lar_start) FillCut(h_den_cutflow, "LAr fid start");
      if (tms_end) FillCut(h_den_cutflow, "TMS fid end");
      if (lar_start && tms_end) {
        FillCut(h_den_cutflow, "strict denominator");
        counts.denom_strict++;
        h_strict_den->Fill(ke_tms_enter);
      }
    }
  }

  counts.truth_info_entries = truth->GetEntriesFast();
  counts.reco_entries = reco->GetEntriesFast();
  int last_truth_tree = -1;
  int last_spill = -999999999;
  for (Long64_t entry = 0; entry < truth->GetEntriesFast(); ++entry) {
    truth->GetEntry(entry);
    if (entry < reco->GetEntriesFast()) reco->GetEntry(entry);
    const int current_truth_tree =
        truth->fChain ? truth->fChain->GetTreeNumber() : -1;
    if (current_truth_tree != last_truth_tree || truth->SpillNo != last_spill) {
      reconstructed_fingerprints.clear();
      last_truth_tree = current_truth_tree;
      last_spill = truth->SpillNo;
    }

    const int reco_ntracks = entry < reco->GetEntriesFast() ? reco->nTracks : -1;
    counts.truth_info_reco_tracks += truth->RecoTrackN;
    if (reco_ntracks >= 0) counts.reco_tree_tracks += reco_ntracks;
    if (reco_ntracks >= 0 && truth->RecoTrackN != reco_ntracks) {
      counts.reco_track_count_mismatch++;
      h_reco_trackn_compare->Fill(truth->RecoTrackN - reco_ntracks);
      AddExample(examples, "RecoTrackN mismatch entry=" + std::to_string(entry) +
                               " truth.RecoTrackN=" +
                               std::to_string(truth->RecoTrackN) +
                               " reco.nTracks=" + std::to_string(reco_ntracks));
    }
    if (truth->RecoTrackN > Truth_Info::kMaxRecoTracks) {
      counts.numerator_reco_trackn_over_max++;
      AddExample(examples, "RecoTrackN over max entry=" + std::to_string(entry) +
                               " RecoTrackN=" +
                               std::to_string(truth->RecoTrackN));
    }

    const int n_tracks =
        std::max(0, std::min(truth->RecoTrackN, Truth_Info::kMaxRecoTracks));
    for (int it = 0; it < n_tracks; ++it) {
      FillCut(h_num_cutflow, "all reco summaries");
      const bool is_muon = std::abs(truth->RecoTrackPrimaryParticlePDG[it]) == 13;
      FillCut(h_num_passfail, is_muon ? "muon pdg pass" : "muon pdg fail");
      if (!is_muon) continue;
      FillCut(h_num_cutflow, "muon pdg");
      counts.numerator_muons++;

      const double visible = truth->RecoTrackPrimaryParticleTrueVisibleEnergy[it];
      const bool tms_touch = visible >= kMinVisibleEnergy;
      const bool lar_start = truth->RecoTrackPrimaryParticleLArFiducialStart[it];
      const bool tms_end = truth->RecoTrackPrimaryParticleTMSFiducialEnd[it];
      const double ke_tms_enter =
          truth->RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3] * kGeV;
      const bool ke_in_bounds = InKeBounds(ke_tms_enter);

      h_num_visible->Fill(visible);
      h_num_ke->Fill(ke_tms_enter);
      h_num_lar->Fill(lar_start ? 1 : 0);
      h_num_tms->Fill(tms_end ? 1 : 0);
      FillCut(h_num_passfail, tms_touch ? "TMS visible pass" : "TMS visible fail");
      FillCut(h_num_passfail, ke_in_bounds ? "KE bounds pass" : "KE bounds fail");
      FillCut(h_num_passfail, lar_start ? "LAr start pass" : "LAr start fail");
      FillCut(h_num_passfail, tms_end ? "TMS end pass" : "TMS end fail");

      if (ke_tms_enter < 0) {
        counts.numerator_negative_ke++;
        AddExample(examples, "Numerator negative KE: truth entry=" +
                                 std::to_string(entry) + " track=" +
                                 std::to_string(it) + " KE=" +
                                 std::to_string(ke_tms_enter));
      }
      if (DefaultLike(truth->RecoTrackPrimaryParticleTrueMomentumEnteringTMS[it][3]))
        counts.numerator_default_ke++;
      if (!ke_in_bounds) counts.numerator_out_of_bounds_ke++;

      if (!tms_touch) continue;
      FillCut(h_num_cutflow, "TMS visible >= 5 MeV");
      counts.numerator_tms_touch++;
      if (ke_in_bounds) FillCut(h_num_cutflow, "KE in hist bounds");

      const std::string fingerprint = RecoTruthFingerprint(*truth, it);
      reconstructed_fingerprints[fingerprint]++;
      const bool not_duplicate = reconstructed_fingerprints[fingerprint] == 1;
      h_duplicate->Fill(not_duplicate ? 0 : 1);
      FillCut(h_num_passfail,
              not_duplicate ? "duplicate pass" : "duplicate fail");
      if (!not_duplicate) {
        counts.numerator_duplicates++;
        continue;
      }
      FillCut(h_num_cutflow, "not duplicate");

      h_all_num->Fill(ke_tms_enter);
      counts.numerator_all_unique++;

      FillStrictReason(h_num_strict_reason, lar_start, tms_end);
      if (lar_start) counts.numerator_unique_lar_start++;
      if (tms_end) counts.numerator_unique_tms_end++;
      if (lar_start && !tms_end) counts.numerator_unique_lar_only++;
      if (!lar_start && tms_end) counts.numerator_unique_tms_only++;
      if (!lar_start && !tms_end) counts.numerator_unique_neither_strict_flag++;
      if (!(lar_start && tms_end)) {
        std::ostringstream msg;
        msg << "Numerator strict fail: truth entry=" << entry
            << " SpillNo=" << truth->SpillNo
            << " track=" << it
            << " KE=" << ke_tms_enter
            << " visible=" << visible
            << " LArStart=" << lar_start
            << " TMSEnd=" << tms_end
            << " start=("
            << truth->RecoTrackPrimaryParticleTruePositionStart[it][0] << ","
            << truth->RecoTrackPrimaryParticleTruePositionStart[it][1] << ","
            << truth->RecoTrackPrimaryParticleTruePositionStart[it][2] << ")"
            << " end=("
            << truth->RecoTrackPrimaryParticleTruePositionEnd[it][0] << ","
            << truth->RecoTrackPrimaryParticleTruePositionEnd[it][1] << ","
            << truth->RecoTrackPrimaryParticleTruePositionEnd[it][2] << ")"
            << " enterTMS=("
            << truth->RecoTrackPrimaryParticleTruePositionEnteringTMS[it][0]
            << ","
            << truth->RecoTrackPrimaryParticleTruePositionEnteringTMS[it][1]
            << ","
            << truth->RecoTrackPrimaryParticleTruePositionEnteringTMS[it][2]
            << ")";
        AddExample(examples, msg.str());
      }

      if (lar_start) FillCut(h_num_cutflow, "LAr fid start");
      if (tms_end) FillCut(h_num_cutflow, "TMS fid end");
      if (lar_start && tms_end) {
        FillCut(h_num_cutflow, "strict numerator");
        h_strict_num->Fill(ke_tms_enter);
        counts.numerator_strict_unique++;
      }
    }
  }

  TH1D *h_all_eff = MakeEfficiency("all_muon_ke_tms_enter",
                                   "Reco efficiency, all TMS-touching muons",
                                   h_all_num, h_all_den);
  TH1D *h_strict_eff = MakeEfficiency(
      "muon_ke_tms_enter",
      "Reco efficiency, LAr-start/TMS-end muons", h_strict_num, h_strict_den);

  TCanvas canvas("canvas", "canvas", 900, 700);
  DrawEfficiency(h_all_eff, canvas, output_dir, "all_muon_ke_tms_enter");
  DrawEfficiency(h_strict_eff, canvas, output_dir, "muon_ke_tms_enter");

  DrawAndSave(h_all_den, canvas, raw_dir, "all_muon_ke_tms_enter_denominator");
  DrawAndSave(h_all_num, canvas, raw_dir, "all_muon_ke_tms_enter_numerator");
  DrawAndSave(h_strict_den, canvas, raw_dir, "muon_ke_tms_enter_denominator");
  DrawAndSave(h_strict_num, canvas, raw_dir, "muon_ke_tms_enter_numerator");

  std::vector<std::pair<TH1 *, std::string>> diagnostics = {
      {h_den_visible, "den_muon_true_visible_energy"},
      {h_num_visible, "num_muon_true_visible_energy"},
      {h_den_ke, "den_muon_ke_tms_enter_diagnostic"},
      {h_num_ke, "num_muon_ke_tms_enter_diagnostic"},
      {h_den_lar, "den_lar_fiducial_start"},
      {h_den_tms, "den_tms_fiducial_end"},
      {h_num_lar, "num_lar_fiducial_start"},
      {h_num_tms, "num_tms_fiducial_end"},
      {h_reco_trackn_compare, "truth_info_recotrackn_minus_reco_tree_ntracks"},
      {h_duplicate, "num_duplicate_fingerprint"},
      {h_den_passfail, "denominator_cut_pass_fail"},
      {h_num_passfail, "numerator_cut_pass_fail"},
      {h_den_strict_reason, "denominator_strict_cut_reason"},
      {h_num_strict_reason, "numerator_strict_cut_reason"},
      {h_den_cutflow, "denominator_cutflow"},
      {h_num_cutflow, "numerator_cutflow"},
  };
  for (auto &item : diagnostics) {
    DrawAndSave(item.first, canvas, diag_dir, item.second);
  }

  const std::string root_output = output_dir + "/RecoEffIsolated.root";
  TFile out(root_output.c_str(), "RECREATE");
  h_all_den->Write();
  h_all_num->Write();
  h_strict_den->Write();
  h_strict_num->Write();
  h_all_eff->Write();
  h_strict_eff->Write();
  for (auto &item : diagnostics) item.first->Write();
  out.Close();

  WriteSummary(output_dir + "/diagnostics.txt", counts, examples, input_files);

  std::cout << "RecoEffIsolated output written to " << output_dir << std::endl;
  std::cout << "  efficiencies: " << output_dir << std::endl;
  std::cout << "  raw hists:    " << raw_dir << std::endl;
  std::cout << "  diagnostics:  " << diag_dir << std::endl;
  std::cout << "  root file:    " << root_output << std::endl;
  return 0;
}
