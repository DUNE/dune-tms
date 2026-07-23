#include <TCanvas.h>
#include <TDirectory.h>
#include <TError.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLegend.h>
#include <TROOT.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Comparison {
  std::string path;
  double chi2 = 0.0;
  int ndf = 0;
  double integral_a = 0.0;
  double integral_b = 0.0;
  double max_abs_bin_delta = 0.0;
  bool only_a = false;
  bool only_b = false;
};

std::string DefaultDataDir() {
  const char *user = std::getenv("USER");
  std::string base = "/exp/dune/data/users/";
  base += user ? user : "unknown";
  base += "/dune-tms/Validation/BranchHists/data";
  return base;
}

std::string DefaultComparisonDir(const std::string &file_a_path,
                                 const std::string &file_b_path) {
  const char *user = std::getenv("USER");
  std::string base = "/exp/dune/data/users/";
  base += user ? user : "unknown";
  base += "/dune-tms/Validation/BranchHists/comparisons";
  const std::string stem_a = std::filesystem::path(file_a_path).stem().string();
  const std::string stem_b = std::filesystem::path(file_b_path).stem().string();
  return (std::filesystem::path(base) / (stem_a + "_vs_" + stem_b)).string();
}

std::string PrepareComparisonDir(const std::string &file_a_path,
                                 const std::string &file_b_path) {
  std::error_code ec;
  std::string output_dir = DefaultComparisonDir(file_a_path, file_b_path);
  std::filesystem::create_directories(output_dir, ec);
  if (!ec) return output_dir;

  const std::string stem_a = std::filesystem::path(file_a_path).stem().string();
  const std::string stem_b = std::filesystem::path(file_b_path).stem().string();
  output_dir =
      (std::filesystem::path(file_a_path).parent_path() /
       ("comparison_" + stem_a + "_vs_" + stem_b)).string();
  ec.clear();
  std::filesystem::create_directories(output_dir, ec);
  if (ec) {
    throw std::runtime_error("Could not create comparison plot directory: " +
                             output_dir + " (" + ec.message() + ")");
  }
  return output_dir;
}

std::string ResolveProjectFile(const std::string &project_or_file) {
  const std::filesystem::path path(project_or_file);
  if (std::filesystem::exists(path)) return path.string();
  if (path.has_parent_path() || path.extension() == ".root") return path.string();
  return (std::filesystem::path(DefaultDataDir()) /
          (project_or_file + ".root")).string();
}

std::string SanitizeFilename(const std::string &input) {
  std::string output;
  output.reserve(input.size());
  for (char c : input) {
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
        (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.') {
      output.push_back(c);
    } else {
      output.push_back('_');
    }
  }
  return output;
}

std::string JoinPath(const std::string &dir, const std::string &name) {
  if (dir.empty()) return name;
  return dir + "/" + name;
}

void CollectHistPaths(TDirectory *dir, const std::string &prefix,
                      std::set<std::string> &paths) {
  TIter next(dir->GetListOfKeys());
  while (auto *obj = next()) {
    auto *key = dynamic_cast<TKey *>(obj);
    if (!key) continue;
    const std::string name = key->GetName();
    const std::string class_name = key->GetClassName();
    const std::string path = JoinPath(prefix, name);

    if (class_name == "TDirectoryFile" || class_name == "TDirectory") {
      auto *subdir = dynamic_cast<TDirectory *>(dir->Get(name.c_str()));
      if (subdir) CollectHistPaths(subdir, path, paths);
      continue;
    }

    std::unique_ptr<TObject> object(key->ReadObj());
    if (object && object->InheritsFrom(TH1::Class())) paths.insert(path);
  }
}

std::unique_ptr<TH1> ReadHist(TFile &file, const std::string &path) {
  auto *object = file.Get(path.c_str());
  if (!object || !object->InheritsFrom(TH1::Class())) return nullptr;
  auto *hist = dynamic_cast<TH1 *>(object);
  if (!hist) return nullptr;
  auto *clone = dynamic_cast<TH1 *>(hist->Clone());
  if (!clone) return nullptr;
  clone->SetDirectory(nullptr);
  return std::unique_ptr<TH1>(clone);
}

std::unique_ptr<TH1D> ProjectToAxis(const TH1 &source, const std::string &name,
                                    int bins, double low, double high) {
  auto out = std::make_unique<TH1D>(name.c_str(), source.GetTitle(), bins, low, high);
  out->Sumw2();
  out->SetDirectory(nullptr);
  out->SetBinContent(0, source.GetBinContent(0));
  out->SetBinError(0, source.GetBinError(0));
  out->SetBinContent(bins + 1, source.GetBinContent(source.GetNbinsX() + 1));
  out->SetBinError(bins + 1, source.GetBinError(source.GetNbinsX() + 1));

  for (int ib = 1; ib <= source.GetNbinsX(); ++ib) {
    const double content = source.GetBinContent(ib);
    const double error = source.GetBinError(ib);
    if (content == 0.0 && error == 0.0) continue;
    const int target_bin = out->FindBin(source.GetXaxis()->GetBinCenter(ib));
    const double old_content = out->GetBinContent(target_bin);
    const double old_error = out->GetBinError(target_bin);
    out->SetBinContent(target_bin, old_content + content);
    out->SetBinError(target_bin, std::hypot(old_error, error));
  }
  return out;
}

Comparison CompareHistograms(const std::string &path, const TH1 &a, const TH1 &b) {
  Comparison result;
  result.path = path;
  result.integral_a = a.Integral(0, a.GetNbinsX() + 1);
  result.integral_b = b.Integral(0, b.GetNbinsX() + 1);

  const double low =
      std::min(a.GetXaxis()->GetXmin(), b.GetXaxis()->GetXmin());
  const double high =
      std::max(a.GetXaxis()->GetXmax(), b.GetXaxis()->GetXmax());
  const int bins = std::max(a.GetNbinsX(), b.GetNbinsX());
  auto aa = ProjectToAxis(a, "comparison_a", bins, low, high);
  auto bb = ProjectToAxis(b, "comparison_b", bins, low, high);

  for (int ib = 0; ib <= bins + 1; ++ib) {
    const double delta = aa->GetBinContent(ib) - bb->GetBinContent(ib);
    const double variance =
        std::pow(aa->GetBinError(ib), 2) + std::pow(bb->GetBinError(ib), 2);
    if (variance > 0.0) {
      result.chi2 += delta * delta / variance;
      result.ndf++;
    } else if (delta != 0.0) {
      result.chi2 += std::abs(delta);
      result.ndf++;
    }
    result.max_abs_bin_delta =
        std::max(result.max_abs_bin_delta, std::abs(delta));
  }
  return result;
}

void StyleHist(TH1 &hist, int color, int marker, const std::string &title) {
  hist.SetLineColor(color);
  hist.SetMarkerColor(color);
  hist.SetMarkerStyle(marker);
  hist.SetLineWidth(2);
  hist.SetTitle(title.c_str());
  hist.GetYaxis()->SetTitle("entries");
}

std::filesystem::path PlotPath(const std::string &output_dir,
                               const std::string &hist_path) {
  const size_t slash = hist_path.find('/');
  std::filesystem::path out(output_dir);
  if (slash != std::string::npos) {
    out /= SanitizeFilename(hist_path.substr(0, slash));
    out /= SanitizeFilename(hist_path.substr(slash + 1)) + ".png";
  } else {
    out /= SanitizeFilename(hist_path) + ".png";
  }
  return out;
}

void DrawComparison(TFile &file_a, TFile &file_b, const Comparison &comparison,
                    const std::string &label_a, const std::string &label_b,
                    const std::string &output_dir) {
  auto hist_a = ReadHist(file_a, comparison.path);
  auto hist_b = ReadHist(file_b, comparison.path);
  if (!hist_a && !hist_b) return;

  std::unique_ptr<TH1> draw_a;
  std::unique_ptr<TH1> draw_b;
  if (hist_a && hist_b) {
    const double low =
        std::min(hist_a->GetXaxis()->GetXmin(), hist_b->GetXaxis()->GetXmin());
    const double high =
        std::max(hist_a->GetXaxis()->GetXmax(), hist_b->GetXaxis()->GetXmax());
    const int bins = std::max(hist_a->GetNbinsX(), hist_b->GetNbinsX());
    draw_a = ProjectToAxis(*hist_a, "plot_a", bins, low, high);
    draw_b = ProjectToAxis(*hist_b, "plot_b", bins, low, high);
  } else if (hist_a) {
    draw_a.reset(dynamic_cast<TH1 *>(hist_a->Clone("plot_a")));
    draw_a->SetDirectory(nullptr);
  } else {
    draw_b.reset(dynamic_cast<TH1 *>(hist_b->Clone("plot_b")));
    draw_b->SetDirectory(nullptr);
  }

  if (draw_a) StyleHist(*draw_a, kBlue + 1, 20, comparison.path);
  if (draw_b) StyleHist(*draw_b, kRed + 1, 24, comparison.path);

  const double max_a = draw_a ? draw_a->GetMaximum() : 0.0;
  const double max_b = draw_b ? draw_b->GetMaximum() : 0.0;
  const double ymax = std::max(max_a, max_b) * 1.25 + 1.0;
  if (draw_a) draw_a->SetMaximum(ymax);
  if (draw_b) draw_b->SetMaximum(ymax);

  TCanvas canvas("branch_compare_canvas", "branch_compare_canvas", 1100, 800);
  canvas.SetTicks();
  if (draw_a) {
    draw_a->Draw("hist e");
    if (draw_b) draw_b->Draw("hist e same");
  } else if (draw_b) {
    draw_b->Draw("hist e");
  }

  TLegend legend(0.62, 0.72, 0.88, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  if (draw_a) legend.AddEntry(draw_a.get(), label_a.c_str(), "l");
  if (draw_b) legend.AddEntry(draw_b.get(), label_b.c_str(), "l");
  legend.Draw();

  const auto plot_path = PlotPath(output_dir, comparison.path);
  std::filesystem::create_directories(plot_path.parent_path());
  canvas.Print(plot_path.string().c_str());
}

bool Changed(const Comparison &comparison, double min_chi2,
             double min_integral_delta) {
  if (comparison.only_a || comparison.only_b) return true;
  if (comparison.chi2 >= min_chi2) return true;
  return std::abs(comparison.integral_b - comparison.integral_a) >=
         min_integral_delta;
}

void PrintUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " <project_a.root|project_a> <project_b.root|project_b> "
               "[min_chi2=1e-9] [min_integral_delta=1e-9]\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    PrintUsage(argv[0]);
    return 1;
  }

  TH1::AddDirectory(false);
  gErrorIgnoreLevel = kWarning;
  const std::string file_a_path = ResolveProjectFile(argv[1]);
  const std::string file_b_path = ResolveProjectFile(argv[2]);
  const double min_chi2 = argc >= 4 ? std::atof(argv[3]) : 1e-9;
  const double min_integral_delta = argc >= 5 ? std::atof(argv[4]) : 1e-9;

  TFile file_a(file_a_path.c_str(), "READ");
  TFile file_b(file_b_path.c_str(), "READ");
  if (file_a.IsZombie()) {
    std::cerr << "Could not open " << file_a_path << std::endl;
    return 2;
  }
  if (file_b.IsZombie()) {
    std::cerr << "Could not open " << file_b_path << std::endl;
    return 3;
  }

  std::set<std::string> paths_a;
  std::set<std::string> paths_b;
  CollectHistPaths(&file_a, "", paths_a);
  CollectHistPaths(&file_b, "", paths_b);

  std::set<std::string> all_paths = paths_a;
  all_paths.insert(paths_b.begin(), paths_b.end());

  std::vector<Comparison> changed;
  for (const auto &path : all_paths) {
    const bool has_a = paths_a.count(path) != 0;
    const bool has_b = paths_b.count(path) != 0;
    Comparison comparison;
    comparison.path = path;
    comparison.only_a = has_a && !has_b;
    comparison.only_b = has_b && !has_a;

    if (has_a && has_b) {
      auto hist_a = ReadHist(file_a, path);
      auto hist_b = ReadHist(file_b, path);
      if (!hist_a || !hist_b) continue;
      comparison = CompareHistograms(path, *hist_a, *hist_b);
    }

    if (Changed(comparison, min_chi2, min_integral_delta))
      changed.push_back(comparison);
  }

  std::sort(changed.begin(), changed.end(),
            [](const Comparison &a, const Comparison &b) {
              if (a.only_a != b.only_a) return a.only_a > b.only_a;
              if (a.only_b != b.only_b) return a.only_b > b.only_b;
              return a.chi2 > b.chi2;
            });

  std::cout << "Comparing branch histograms\n";
  std::cout << "  A: " << file_a_path << "\n";
  std::cout << "  B: " << file_b_path << "\n";
  std::cout << "  Changed histograms: " << changed.size() << "\n";
  const std::string output_dir = PrepareComparisonDir(file_a_path, file_b_path);
  std::cout << "  Plots: " << output_dir << "\n";
  std::cout << std::setprecision(8);
  for (const auto &comparison : changed) {
    if (comparison.only_a) {
      std::cout << "ONLY_A " << comparison.path << "\n";
      continue;
    }
    if (comparison.only_b) {
      std::cout << "ONLY_B " << comparison.path << "\n";
      continue;
    }
    std::cout << "CHANGED " << comparison.path << " chi2=" << comparison.chi2
              << " ndf=" << comparison.ndf
              << " chi2_ndf="
              << (comparison.ndf > 0 ? comparison.chi2 / comparison.ndf : 0.0)
              << " integral_a=" << comparison.integral_a
              << " integral_b=" << comparison.integral_b
              << " integral_delta="
              << (comparison.integral_b - comparison.integral_a)
              << " max_abs_bin_delta=" << comparison.max_abs_bin_delta << "\n";
  }

  const std::string label_a = std::filesystem::path(file_a_path).stem().string();
  const std::string label_b = std::filesystem::path(file_b_path).stem().string();
  for (const auto &comparison : changed) {
    DrawComparison(file_a, file_b, comparison, label_a, label_b, output_dir);
  }

  return 0;
}
