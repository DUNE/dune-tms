#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

const std::vector<std::string> kTreeNames = {"Reco_Tree", "Truth_Info",
                                             "Truth_Spill"};

struct BranchSpec {
  std::string tree;
  std::string branch;
  std::string leaf;
  std::string expression;
  std::string type;
  std::string component;
  std::string hist_name;
  std::string title;
};

struct AxisSpec {
  int bins = 100;
  double low = 0.0;
  double high = 1.0;
  std::string source = "dynamic";
};

std::string Sanitize(const std::string &input) {
  std::string out;
  out.reserve(input.size());
  for (char c : input) {
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
        (c >= '0' && c <= '9')) {
      out.push_back(c);
    } else {
      out.push_back('_');
    }
  }
  return out;
}

bool IsRootFile(const std::filesystem::path &path) {
  return path.extension() == ".root";
}

bool HasTree(const std::string &filename, const std::string &tree_name) {
  TFile file(filename.c_str(), "READ");
  if (file.IsZombie()) return false;
  return file.Get(tree_name.c_str()) != nullptr;
}

std::vector<std::string> FindInputFiles(const std::string &input) {
  std::vector<std::string> files;
  std::error_code ec;
  const std::filesystem::path path(input);

  if (std::filesystem::is_regular_file(path, ec)) {
    if (IsRootFile(path)) {
      for (const auto &tree : kTreeNames) {
        if (HasTree(path.string(), tree)) {
          files.push_back(path.string());
          break;
        }
      }
    }
  } else if (std::filesystem::is_directory(path, ec)) {
    for (const auto &entry : std::filesystem::recursive_directory_iterator(
             path, std::filesystem::directory_options::skip_permission_denied,
             ec)) {
      if (ec) {
        std::cerr << "Warning while scanning " << path << ": " << ec.message()
                  << std::endl;
        ec.clear();
        continue;
      }
      if (!entry.is_regular_file(ec) || !IsRootFile(entry.path())) continue;
      for (const auto &tree : kTreeNames) {
        if (HasTree(entry.path().string(), tree)) {
          files.push_back(entry.path().string());
          break;
        }
      }
    }
  } else {
    std::cerr << "Input is neither a ROOT file nor a directory: " << input
              << std::endl;
  }

  std::sort(files.begin(), files.end());
  files.erase(std::unique(files.begin(), files.end()), files.end());
  return files;
}

std::string DefaultOutputDir() {
  const char *user = std::getenv("USER");
  std::string base = "/exp/dune/data/users/";
  base += user ? user : "unknown";
  base += "/dune-tms/Validation/BranchHists/data";
  return base;
}

std::string OutputPath(const std::string &project) {
  std::filesystem::path path(project);
  if (path.has_parent_path() || path.extension() == ".root") return project;
  return (std::filesystem::path(DefaultOutputDir()) / (project + ".root")).string();
}

bool IsNumericType(const std::string &type) {
  static const std::set<std::string> numeric = {
      "Bool_t",  "Char_t",   "UChar_t",  "Short_t", "UShort_t",
      "Int_t",   "UInt_t",   "Float_t",  "Double_t", "Long64_t",
      "ULong64_t", "Long_t", "ULong_t",  "bool",    "char",
      "unsigned char",       "short",    "unsigned short",
      "int",     "unsigned int",         "float",   "double",
      "long",    "unsigned long",        "long long",
      "unsigned long long"};
  return numeric.count(type) != 0;
}

bool IsBoolType(const std::string &type) {
  return type == "Bool_t" || type == "bool";
}

std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

bool Contains(const std::string &haystack, const std::string &needle) {
  return haystack.find(needle) != std::string::npos;
}

bool EndsWith(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> ParseLeafDimensions(const std::string &title) {
  std::vector<std::string> dimensions;
  size_t pos = 0;
  while ((pos = title.find('[', pos)) != std::string::npos) {
    const size_t end = title.find(']', pos + 1);
    if (end == std::string::npos) break;
    dimensions.push_back(title.substr(pos + 1, end - pos - 1));
    pos = end + 1;
  }
  return dimensions;
}

int NumericDimension(const std::string &dimension) {
  if (dimension.empty()) return -1;
  if (!std::all_of(dimension.begin(), dimension.end(),
                   [](unsigned char c) { return std::isdigit(c); })) {
    return -1;
  }
  return std::atoi(dimension.c_str());
}

std::string ComponentExpression(const std::string &leaf_name, size_t ndims,
                                int component) {
  std::string expression = leaf_name;
  for (size_t idim = 1; idim < ndims; ++idim) expression += "[]";
  expression += "[" + std::to_string(component) + "]";
  return expression;
}

std::vector<std::pair<std::string, std::string>> ComponentExpressions(
    const std::string &leaf_name, const std::vector<std::string> &dimensions) {
  std::vector<std::pair<std::string, std::string>> components;
  if (dimensions.empty()) return components;
  const int last_dim = NumericDimension(dimensions.back());
  if (last_dim != 3 && last_dim != 4) return components;

  static const std::vector<std::string> labels = {"X", "Y", "Z", "T"};
  for (int component = 0; component < last_dim; ++component) {
    components.push_back({labels[component],
                          ComponentExpression(leaf_name, dimensions.size(),
                                              component)});
  }
  return components;
}

std::vector<BranchSpec> CollectBranches(TChain &chain,
                                        const std::string &tree_name,
                                        std::ostream &manifest) {
  std::vector<BranchSpec> specs;
  TObjArray *branches = chain.GetListOfBranches();
  if (!branches) return specs;

  for (int ib = 0; ib < branches->GetEntries(); ++ib) {
    auto *branch = dynamic_cast<TBranch *>(branches->At(ib));
    if (!branch) continue;

    TObjArray *leaves = branch->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() == 0) {
      manifest << "SKIP " << tree_name << "/" << branch->GetName()
               << " no leaves\n";
      continue;
    }

    for (int il = 0; il < leaves->GetEntries(); ++il) {
      auto *leaf = dynamic_cast<TLeaf *>(leaves->At(il));
      if (!leaf) continue;
      const std::string type = leaf->GetTypeName();
      const std::string leaf_name = leaf->GetName();
      const std::string branch_name = branch->GetName();
      if (!IsNumericType(type)) {
        manifest << "SKIP " << tree_name << "/" << branch_name << "/"
                 << leaf_name << " type=" << type << "\n";
        continue;
      }

      BranchSpec spec;
      spec.tree = tree_name;
      spec.branch = branch_name;
      spec.leaf = leaf_name;
      spec.type = type;
      spec.expression = leaf_name;
      if (leaves->GetEntries() > 1 && leaf_name != branch_name)
        spec.expression = branch_name + "." + leaf_name;

      const auto dimensions = ParseLeafDimensions(leaf->GetTitle());
      const auto components = ComponentExpressions(leaf_name, dimensions);
      if (!components.empty() && leaves->GetEntries() == 1) {
        for (const auto &component : components) {
          BranchSpec component_spec = spec;
          component_spec.component = component.first;
          component_spec.expression = component.second;
          component_spec.hist_name = Sanitize(tree_name + "__" + branch_name +
                                              "__" + leaf_name + "__" +
                                              component.first);
          component_spec.title = tree_name + "/" + branch_name + "/" +
                                 leaf_name + "[" + component.first +
                                 "];value;entries";
          specs.push_back(component_spec);
        }
      } else {
        spec.hist_name = Sanitize(tree_name + "__" + branch_name + "__" +
                                  leaf_name);
        spec.title = tree_name + "/" + branch_name + "/" + leaf_name +
                     ";value;entries";
        specs.push_back(spec);
      }
    }
  }
  return specs;
}

double PadLow(double min, double max) {
  if (min == max) return min - (min == 0 ? 0.5 : std::abs(min) * 0.05);
  return min - 0.02 * (max - min);
}

double PadHigh(double min, double max) {
  if (min == max) return max + (max == 0 ? 0.5 : std::abs(max) * 0.05);
  return max + 0.02 * (max - min);
}

double NiceStep(double raw_step) {
  if (!(raw_step > 0.0) || !std::isfinite(raw_step)) return 1.0;
  const double exponent = std::floor(std::log10(raw_step));
  const double scale = std::pow(10.0, exponent);
  const double normalized = raw_step / scale;
  if (normalized <= 1.0) return scale;
  if (normalized <= 2.0) return 2.0 * scale;
  if (normalized <= 5.0) return 5.0 * scale;
  return 10.0 * scale;
}

AxisSpec DynamicAxis(double min, double max) {
  AxisSpec axis;
  axis.source = "dynamic_nice";
  if (min == max) {
    const double width = min == 0.0 ? 1.0 : std::abs(min) * 0.1;
    min -= width;
    max += width;
  }
  const double step = NiceStep((max - min) / 50.0);
  axis.low = std::floor(min / step) * step;
  axis.high = std::ceil(max / step) * step;
  if (!(axis.low < axis.high)) {
    axis.low = PadLow(min, max);
    axis.high = PadHigh(min, max);
  }
  return axis;
}

AxisSpec ChooseAxis(const BranchSpec &spec, double min, double max) {
  AxisSpec axis;
  if (IsBoolType(spec.type)) {
    axis.bins = 2;
    axis.low = -0.5;
    axis.high = 1.5;
    axis.source = "bool";
    return axis;
  }

  const std::string name = Lower(spec.branch + " " + spec.leaf);
  std::string compact_name = name;
  compact_name.erase(std::remove_if(compact_name.begin(), compact_name.end(),
                                    [](unsigned char c) {
                                      return std::isspace(c) || c == '_';
                                    }),
                     compact_name.end());
  std::string component = Lower(spec.component);
  if (component.empty() &&
      (Contains(compact_name, "pos") || Contains(compact_name, "position") ||
       Contains(compact_name, "vertex") || Contains(compact_name, "vtx") ||
       Contains(compact_name, "hit"))) {
    if (EndsWith(compact_name, "x")) component = "x";
    if (EndsWith(compact_name, "y")) component = "y";
    if (EndsWith(compact_name, "z")) component = "z";
    if (EndsWith(compact_name, "t")) component = "t";
  }

  if (Contains(name, "runno") || name == "run") {
    axis.bins = 220;
    axis.low = -0.5;
    axis.high = 1100000000.5;
    axis.source = "fixed_run";
  } else if (Contains(name, "globalid") || Contains(name, "global_id")) {
    axis.bins = 200;
    axis.low = -0.5;
    axis.high = 1100000000000000.5;
    axis.source = "fixed_global_id";
  } else if (Contains(name, "pdg")) {
    axis.bins = 600;
    axis.low = -3000.5;
    axis.high = 3000.5;
    axis.source = "fixed_pdg";
  } else if (Contains(name, "charge")) {
    axis.bins = 8;
    axis.low = -3.5;
    axis.high = 4.5;
    axis.source = "fixed_charge";
  } else if (Contains(name, "ntruehits") || Contains(name, "nhits")) {
    axis.bins = 200;
    axis.low = -0.5;
    axis.high = 20000.5;
    axis.source = "fixed_hit_count";
  } else if (Contains(name, "ntrueparticles")) {
    axis.bins = 200;
    axis.low = -0.5;
    axis.high = 20000.5;
    axis.source = "fixed_particle_count";
  } else if (Contains(name, "ntracks") || Contains(name, "recotrackn")) {
    axis.bins = 101;
    axis.low = -0.5;
    axis.high = 100.5;
    axis.source = "fixed_track_count";
  } else if (Contains(name, "truevtxn") || Contains(name, "nprimaryvertices")) {
    axis.bins = 200;
    axis.low = -0.5;
    axis.high = 5000.5;
    axis.source = "fixed_vertex_count";
  } else if (Contains(name, "trackid") || Contains(name, "parent") ||
             Contains(name, "vertexid") || Contains(name, "vtxid") ||
             Contains(name, "index")) {
    axis.bins = 200;
    axis.low = -100.5;
    axis.high = 1000000.5;
    axis.source = "fixed_local_id";
  } else if (Contains(name, "momentum") || Contains(name, "p4") ||
             Contains(name, "px") || Contains(name, "py") ||
             Contains(name, "pz")) {
    axis.bins = 200;
    if (component == "e" || component == "t" || Contains(name, "energy")) {
      axis.low = -100.0;
      axis.high = 50000.0;
      axis.source = "fixed_momentum_e";
    } else {
      axis.low = -50000.0;
      axis.high = 50000.0;
      axis.source = "fixed_momentum_xyz";
    }
  } else if (Contains(name, "energy") || Contains(name, "evis")) {
    axis.bins = 200;
    axis.low = -100.0;
    axis.high = 50000.0;
    axis.source = "fixed_energy";
  } else if (EndsWith(compact_name, "e") || Contains(compact_name, "pe") ||
             Contains(compact_name, "dedx")) {
    axis.bins = 200;
    axis.low = -100.0;
    axis.high = 50000.0;
    axis.source = "fixed_energy_like";
  } else if (Contains(name, "length")) {
    axis.bins = 200;
    axis.low = -100.0;
    axis.high = 50000.0;
    axis.source = "fixed_length";
  } else if (Contains(name, "chi2")) {
    axis.bins = 200;
    axis.low = 0.0;
    axis.high = 1000.0;
    axis.source = "fixed_chi2";
  } else if (Contains(name, "direction")) {
    axis.bins = 120;
    axis.low = -1.2;
    axis.high = 1.2;
    axis.source = "fixed_direction";
  } else if (Contains(name, "pos") || Contains(name, "position") ||
             Contains(name, "vertex") || Contains(name, "vtx") ||
             Contains(compact_name, "hitx") || Contains(compact_name, "hity") ||
             Contains(compact_name, "hitz") || Contains(compact_name, "hitt")) {
    axis.bins = 200;
    if (component == "x" || component == "y") {
      axis.low = -5000.0;
      axis.high = 5000.0;
      axis.source = "fixed_position_xy";
    } else if (component == "z") {
      axis.low = 0.0;
      axis.high = 25000.0;
      axis.source = "fixed_position_z";
    } else if (component == "t") {
      axis.low = -1000.0;
      axis.high = 50000.0;
      axis.source = "fixed_position_t";
    } else {
      axis = DynamicAxis(min, max);
    }
  } else {
    axis = DynamicAxis(min, max);
  }
  return axis;
}

std::string BranchLabel(const BranchSpec &spec) {
  std::string label = spec.tree + "/" + spec.branch + "/" + spec.leaf;
  if (!spec.component.empty()) label += "[" + spec.component + "]";
  return label;
}

TH1D *MakeHist(TChain &chain, const BranchSpec &spec, std::ostream &manifest) {
  const double min = chain.GetMinimum(spec.expression.c_str());
  const double max = chain.GetMaximum(spec.expression.c_str());
  if (!std::isfinite(min) || !std::isfinite(max)) {
    manifest << "SKIP " << BranchLabel(spec)
             << " could not determine finite range\n";
    return nullptr;
  }

  const AxisSpec axis = ChooseAxis(spec, min, max);
  if (!(axis.low < axis.high)) {
    manifest << "SKIP " << BranchLabel(spec) << " invalid range min=" << min
             << " max=" << max << " axis=[" << axis.low << ","
             << axis.high << "]\n";
    return nullptr;
  }

  auto *hist = new TH1D(spec.hist_name.c_str(), spec.title.c_str(), axis.bins,
                        axis.low, axis.high);
  hist->Sumw2();
  hist->SetDirectory(gDirectory);
  const std::string draw = spec.expression + ">>" + spec.hist_name;
  const Long64_t filled = chain.Draw(draw.c_str(), "", "goff");
  hist->SetDirectory(nullptr);
  if (filled < 0) {
    manifest << "SKIP " << BranchLabel(spec)
             << " draw failed expression=" << spec.expression << "\n";
    delete hist;
    return nullptr;
  }

  hist->SetTitle(spec.title.c_str());
  hist->GetXaxis()->SetTitle("value");
  hist->GetYaxis()->SetTitle("entries");
  hist->SetEntries(filled);
  const double integral = hist->Integral(0, axis.bins + 1);
  const double underflow = hist->GetBinContent(0);
  const double overflow = hist->GetBinContent(axis.bins + 1);
  const double flow_fraction =
      integral > 0.0 ? (underflow + overflow) / integral : 0.0;
  manifest << "HIST " << BranchLabel(spec)
           << " name=" << spec.hist_name << " type=" << spec.type
           << " expression=" << spec.expression << " filled=" << filled
           << " bins=" << axis.bins << " range=[" << std::setprecision(12)
           << axis.low << "," << axis.high << "]"
           << " root_reported_range=[" << min << "," << max << "]"
           << " axis_source=" << axis.source << " integral=" << integral
           << " underflow=" << underflow << " overflow=" << overflow
           << " flow_fraction=" << flow_fraction << "\n";
  if (flow_fraction > 0.20) {
    const std::string warning =
        "WARNING " + BranchLabel(spec) + " has " +
        std::to_string(100.0 * flow_fraction) +
        "% of entries in under/overflow";
    manifest << warning << "\n";
    std::cerr << warning << std::endl;
  }
  return hist;
}

void WriteInputs(TFile &out, const std::vector<std::string> &files,
                 const std::string &project) {
  out.mkdir("metadata");
  out.cd("metadata");
  auto *inputs = new TNamed("input_files", "");
  std::ostringstream body;
  body << "project=" << project << "\n";
  for (const auto &file : files) body << file << "\n";
  inputs->SetTitle(body.str().c_str());
  inputs->Write();
  delete inputs;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0]
              << " <input.root|input_directory> <project_name|output.root>\n";
    return 1;
  }

  TH1::AddDirectory(false);
  const std::string input = argv[1];
  const std::string project = argv[2];
  const std::string output_path = OutputPath(project);
  std::filesystem::create_directories(std::filesystem::path(output_path).parent_path());

  const std::vector<std::string> files = FindInputFiles(input);
  if (files.empty()) {
    std::cerr << "No ROOT files with Reco_Tree, Truth_Info, or Truth_Spill found in "
              << input << std::endl;
    return 2;
  }

  const std::string manifest_path =
      (std::filesystem::path(output_path).replace_extension(".manifest.txt")).string();
  std::ofstream manifest(manifest_path);
  manifest << "BranchHists project=" << project << "\n";
  manifest << "output=" << output_path << "\n";
  manifest << "inputs=" << files.size() << "\n";
  for (const auto &file : files) manifest << "INPUT " << file << "\n";

  TFile out(output_path.c_str(), "RECREATE");
  if (out.IsZombie()) {
    std::cerr << "Could not create " << output_path << std::endl;
    return 3;
  }
  WriteInputs(out, files, project);

  for (const auto &tree_name : kTreeNames) {
    TChain chain(tree_name.c_str());
    int added = 0;
    for (const auto &file : files) {
      if (HasTree(file, tree_name)) {
        chain.Add(file.c_str());
        ++added;
      }
    }
    if (added == 0 || chain.GetEntries() == 0) {
      manifest << "SKIP_TREE " << tree_name << " added_files=" << added
               << " entries=" << chain.GetEntries() << "\n";
      continue;
    }

    manifest << "TREE " << tree_name << " files=" << added
             << " entries=" << chain.GetEntries() << "\n";
    out.mkdir(tree_name.c_str());
    out.cd(tree_name.c_str());

    const auto specs = CollectBranches(chain, tree_name, manifest);
    for (const auto &spec : specs) {
      std::unique_ptr<TH1D> hist(MakeHist(chain, spec, manifest));
      if (!hist) continue;
      out.cd(tree_name.c_str());
      hist->Write();
    }
  }

  out.Close();
  std::cout << "Branch histograms written to " << output_path << std::endl;
  std::cout << "Manifest written to " << manifest_path << std::endl;
  return 0;
}
