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
  std::string hist_name;
  std::string title;
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
      spec.hist_name = Sanitize(tree_name + "__" + branch_name + "__" +
                                leaf_name);
      spec.title = tree_name + "/" + branch_name + "/" + leaf_name +
                   ";value;entries";
      specs.push_back(spec);
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

TH1D *MakeHist(TChain &chain, const BranchSpec &spec, std::ostream &manifest) {
  double min = 0.0;
  double max = 1.0;
  int bins = 100;
  if (IsBoolType(spec.type)) {
    bins = 2;
    min = -0.5;
    max = 1.5;
  } else {
    min = chain.GetMinimum(spec.expression.c_str());
    max = chain.GetMaximum(spec.expression.c_str());
    if (!std::isfinite(min) || !std::isfinite(max)) {
      manifest << "SKIP " << spec.tree << "/" << spec.branch << "/"
               << spec.leaf << " could not determine finite range\n";
      return nullptr;
    }
    const double low = PadLow(min, max);
    const double high = PadHigh(min, max);
    if (!(low < high)) {
      manifest << "SKIP " << spec.tree << "/" << spec.branch << "/"
               << spec.leaf << " invalid range min=" << min << " max=" << max
               << "\n";
      return nullptr;
    }
    min = low;
    max = high;
  }

  auto *hist = new TH1D(spec.hist_name.c_str(), spec.title.c_str(), bins, min, max);
  hist->Sumw2();
  hist->SetDirectory(gDirectory);
  const std::string draw = spec.expression + ">>" + spec.hist_name;
  const Long64_t filled = chain.Draw(draw.c_str(), "", "goff");
  hist->SetDirectory(nullptr);
  if (filled < 0) {
    manifest << "SKIP " << spec.tree << "/" << spec.branch << "/" << spec.leaf
             << " draw failed expression=" << spec.expression << "\n";
    delete hist;
    return nullptr;
  }

  hist->SetTitle(spec.title.c_str());
  hist->GetXaxis()->SetTitle("value");
  hist->GetYaxis()->SetTitle("entries");
  hist->SetEntries(filled);
  manifest << "HIST " << spec.tree << "/" << spec.branch << "/" << spec.leaf
           << " name=" << spec.hist_name << " type=" << spec.type
           << " expression=" << spec.expression << " filled=" << filled
           << " bins=" << bins << " range=[" << std::setprecision(12) << min
           << "," << max << "]"
           << " integral=" << hist->Integral(0, bins + 1) << "\n";
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
