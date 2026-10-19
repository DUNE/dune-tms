// Unit test for the bad channel list (TMS_BadChannels) and the decoder's use of it. Writes small
// lists to the directory given as the first argument (default: the current directory). Exits
// non-zero if any check fails. Needs TMS_DIR set, as the managers do.
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TMS_BadChannels.h"
#include "TMS_Decoder.h"

namespace {

int n_failed = 0;

void Check(bool ok, const char* what) {
  std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
  if (!ok) n_failed++;
}

std::string Write(const std::string& dir, const std::string& name, const std::string& text) {
  const std::string path = dir + "/" + name;
  std::ofstream(path) << text;
  return path;
}

template <typename F> bool Throws(F f) {
  try { f(); } catch (const std::exception&) { return true; }
  return false;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string dir = argc > 1 ? argv[1] : ".";

  Check(TMS_BadChannels("").Size() == 0, "empty file name: no bad channels");
  Check(TMS_BadChannels(Write(dir, "bc_empty.toml", "# nothing\n")).Size() == 0, "file without entries: no bad channels");

  const TMS_BadChannels list(Write(dir, "bc_two.toml",
      "[[BadChannel]]\nPlane = 12\nBar = 40\nView = 3\nStatus = \"dead\"\nReason = \"broken fiber\"\n"
      "[[BadChannel]]\nPlane = 12\nBar = 41\nView = 3\nStatus = \"missing\"\n"));
  Check(list.Size() == 2, "two entries read");
  Check(list.GetStatus({12, 40, 3}) == TMS_BadChannels::Status::Dead, "dead entry");
  Check(list.GetStatus({12, 41, 3}) == TMS_BadChannels::Status::Missing, "missing entry");
  Check(list.GetChannels().at({12, 40, 3}).Reason == "broken fiber", "reason kept");
  Check(list.IsBad({12, 40, 3}) && list.IsBad({12, 41, 3}), "dead and missing are bad");
  Check(!list.IsBad({12, 40, 1}), "same plane and bar, other view: good");
  Check(!list.IsBad({13, 40, 3}), "same bar number, other plane: good");

  Check(Throws([&] { TMS_BadChannels(Write(dir, "bc_state.toml", "[[BadChannel]]\nPlane=1\nBar=1\nView=1\nStatus=\"flaky\"\n")); }),
        "unknown status is an error");
  Check(Throws([&] { TMS_BadChannels(Write(dir, "bc_dup.toml",
        "[[BadChannel]]\nPlane=1\nBar=1\nView=1\nStatus=\"dead\"\n[[BadChannel]]\nPlane=1\nBar=1\nView=1\nStatus=\"missing\"\n")); }),
        "channel listed twice is an error");
  Check(Throws([&] { TMS_BadChannels(Write(dir, "bc_key.toml", "[[BadChannel]]\nPlane=1\nBar=1\nStatus=\"dead\"\n")); }),
        "missing key is an error");
  Check(Throws([&] { TMS_BadChannels(dir + "/does_not_exist.toml"); }), "missing file is an error");

  // The shipped list is empty, so the decoder keeps every readout (the decoder's own use of
  // the list is exercised by the real-spill comparison, as it needs the geometry for TMS_Bar).
  Check(TMS_BadChannels::GetInstance().Size() == 0, "list named in the main TOML is empty by default");

  std::cout << (n_failed ? "FAILED" : "ALL PASSED") << "\n";
  return n_failed ? 1 : 0;
}
