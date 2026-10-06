#ifndef TMS_STAGETIMER_H
#define TMS_STAGETIMER_H

// Wall-clock time per named processing stage, for finding where a run spends
// its time (2026-10-04 timing study). Usage:
//
//   TMS_StageTimer::Clock lap;          // starts timing
//   ...first stage...
//   lap.Lap("det_sim");                 // charges the time since the last lap to "det_sim"
//   ...second stage...
//   lap.Lap("slicing");
//   ...
//   TMS_StageTimer::PrintSummary();     // once, at the end of the run
//
// Stage names may be hierarchical ("c3d/dbscan"); a name is charged however
// many times and from wherever it is lapped, and the summary lists stages in
// the order they were first seen. Not thread-safe: all stages here run on the
// main thread. The cost is two steady_clock reads per lap, negligible next to
// stages that take milliseconds.

#include <chrono>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace TMS_StageTimer {

struct Entry {
  std::string Name;
  double Seconds = 0.0;
  long Calls = 0;
};

// The one registry. Function-local static so it exists whenever it is first used.
inline std::vector<Entry> &Registry() {
  static std::vector<Entry> entries;
  return entries;
}

inline void Add(const std::string &name, double seconds) {
  std::vector<Entry> &entries = Registry();
  for (Entry &entry : entries) {
    if (entry.Name == name) { entry.Seconds += seconds; ++entry.Calls; return; }
  }
  entries.push_back({name, seconds, 1});
}

class Clock {
 public:
  Clock() : fLast(std::chrono::steady_clock::now()) {}
  // Charge the time since construction or the previous Lap to `name`.
  void Lap(const std::string &name) {
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    Add(name, std::chrono::duration<double>(now - fLast).count());
    fLast = now;
  }
  // Restart without charging anyone (for time that belongs to nobody, e.g. a progress printout).
  void Skip() { fLast = std::chrono::steady_clock::now(); }
 private:
  std::chrono::steady_clock::time_point fLast;
};

// Stages sharing a prefix before the first '/' are shown as a group; `total`
// is the wall time of the whole loop, for the percentage column. Lines are
// prefixed so a grid log can be grepped: "STAGETIME".
inline void PrintSummary(double totalSeconds, long nEntries) {
  const std::vector<Entry> &entries = Registry();
  double accounted = 0.0;
  for (const Entry &entry : entries)
    if (entry.Name.find('/') == std::string::npos) accounted += entry.Seconds;
  std::printf("STAGETIME %-28s %10s %7s %12s %10s\n", "stage", "seconds", "%loop", "calls", "ms/call");
  for (const Entry &entry : entries) {
    const bool nested = entry.Name.find('/') != std::string::npos;
    std::printf("STAGETIME %s%-*s %10.2f %6.1f%% %12ld %10.3f\n", nested ? "  " : "", nested ? 26 : 28, entry.Name.c_str(),
                entry.Seconds, totalSeconds > 0 ? 100.0 * entry.Seconds / totalSeconds : 0.0, entry.Calls,
                entry.Calls > 0 ? 1000.0 * entry.Seconds / entry.Calls : 0.0);
  }
  std::printf("STAGETIME top-level stages account for %.2f of %.2f s (%.1f%%), %ld entries, %.4f s/entry\n", accounted,
              totalSeconds, totalSeconds > 0 ? 100.0 * accounted / totalSeconds : 0.0, nEntries,
              nEntries > 0 ? totalSeconds / nEntries : 0.0);
}

}  // namespace TMS_StageTimer

#endif
