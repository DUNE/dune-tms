// Unit test for the pure-math passage/re-segmentation helpers in TMS_Passage.h. No
// geometry or edep-sim input needed. Exits non-zero on the first failed check.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "TMS_Passage.h"

using TMS_Passage::Segment;
using TMS_Passage::OpticalDeposit;

namespace {

int n_failed = 0;

void Check(bool ok, const char* what) {
  std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
  if (!ok) n_failed++;
}

bool Near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }

// Straight segment along z from z0 to z1 in a given bar.
Segment Seg(double z0, double z1, double energy, double t0, double t1, int bar = 1, int traj = 1) {
  return Segment{{0, 0, z0}, {0, 0, z1}, t0, t1, energy, bar, traj};
}

double TotalEnergy(const std::vector<OpticalDeposit>& bins) {
  double e = 0;
  for (const auto& b : bins) e += b.energy;
  return e;
}

// Split one uniform 10 mm / 2 MeV crossing into n equal steps.
std::vector<Segment> UniformSplit(int n) {
  std::vector<Segment> segs;
  for (int i = 0; i < n; ++i) {
    const double f0 = static_cast<double>(i) / n, f1 = static_cast<double>(i + 1) / n;
    segs.push_back(Seg(10 * f0, 10 * f1, 2.0 / n, 10 * f0 / 300.0, 10 * f1 / 300.0));
  }
  return segs;
}

}  // namespace

int main() {
  // 1. A single 10 mm segment cut into 1 mm bins: 10 equal bins, energy conserved.
  {
    auto bins = TMS_Passage::Resegment({Seg(0, 10, 2.0, 0, 1)}, 1.0);
    Check(bins.size() == 10, "single segment -> 10 bins of 1mm");
    Check(Near(TotalEnergy(bins), 2.0), "single segment conserves energy");
    Check(Near(bins[3].energy, 0.2) && Near(bins[3].dx, 1.0), "single segment: uniform 0.2 MeV per 1mm bin");
  }

  // 2. Overlap, not equal division: 4mm/1MeV then 6mm/3MeV into two 5mm bins.
  //    bin0 = all of seg1 (1.0) + 1mm of seg2 (0.5) = 1.5; bin1 = 5mm of seg2 = 2.5.
  {
    auto bins = TMS_Passage::Resegment({Seg(0, 4, 1.0, 0, 1), Seg(4, 10, 3.0, 1, 2, 1, 2)}, 5.0);
    Check(bins.size() == 2, "two segments -> 2 bins of 5mm");
    Check(Near(bins[0].energy, 1.5) && Near(bins[1].energy, 2.5), "energy assigned by geometric overlap");
    Check(Near(TotalEnergy(bins), 4.0), "two segments conserve energy");
    double share_sum = 0;
    for (const auto& s : bins[0].shares) share_sum += s.second;
    Check(Near(share_sum, bins[0].energy) && bins[0].shares.size() == 2, "provenance shares sum to bin energy");
  }

  // 3. The point of the whole exercise: how Geant4 split the same physical deposit must not
  //    change the re-segmented result.
  {
    auto ref = TMS_Passage::Resegment(UniformSplit(1), 1.0);
    bool same = true;
    for (int n : {2, 5, 10, 20, 50}) {
      auto bins = TMS_Passage::Resegment(UniformSplit(n), 1.0);
      if (bins.size() != ref.size()) { same = false; continue; }
      for (size_t i = 0; i < bins.size(); ++i) {
        if (!Near(bins[i].energy, ref[i].energy, 1e-9) || !Near(bins[i].t, ref[i].t, 1e-9)) same = false;
      }
    }
    Check(same, "re-segmented bins identical for 1/2/5/10/20/50-step splits of the same deposit");
  }

  // 4. Non-uniform deposit is preserved (localized high dE/dx is not smeared out).
  {
    // 1 MeV in the first 1 mm, 0.1 MeV in the other 9 mm.
    auto bins = TMS_Passage::Resegment({Seg(0, 1, 1.0, 0, 0.1), Seg(1, 10, 0.1, 0.1, 1.0)}, 1.0);
    Check(Near(bins[0].energy, 1.0) && Near(bins[5].energy, 0.1 / 9.0), "localized energy loss preserved");
  }

  // 5. Zero-length step keeps its energy instead of being dropped.
  {
    auto bins = TMS_Passage::Resegment({Seg(0, 10, 1.0, 0, 1), Seg(5, 5, 0.5, 0.5, 0.5)}, 1.0);
    Check(Near(TotalEnergy(bins), 1.5), "zero-length step energy retained");
  }

  // 6. Midpoint position/time interpolate along the passage.
  {
    auto bins = TMS_Passage::Resegment({Seg(0, 10, 1.0, 0, 10)}, 5.0);
    Check(Near(bins[0].position[2], 2.5) && Near(bins[1].position[2], 7.5), "bin midpoints along the passage");
    Check(Near(bins[0].t, 2.5) && Near(bins[1].t, 7.5), "bin times interpolated at midpoints");
  }

  // 7. Passage building: same bar+trajectory contiguous -> one passage; a gap splits it;
  //    a different trajectory or bar is separate.
  {
    std::vector<Segment> segs = {
      Seg(0, 3, 1, 0, 1), Seg(3, 6, 1, 1, 2),          // contiguous -> one passage
      Seg(20, 23, 1, 5, 6),                            // same bar/traj but 14mm away -> new passage
      Seg(1, 2, 1, 0.5, 0.7, 1, 2),                    // delta ray (different trajectory)
      Seg(0, 3, 1, 0, 1, 2, 1),                        // different bar
    };
    auto passages = TMS_Passage::BuildPassages(segs, 0.5);
    Check(passages.size() == 4, "gap, trajectory, and bar each separate passages (4 total)");
    size_t two_seg = 0;
    for (const auto& p : passages) if (p.size() == 2) two_seg++;
    Check(two_seg == 1, "contiguous same-bar same-trajectory steps stay one passage");
    auto shuffled = segs;
    std::swap(shuffled[0], shuffled[1]);
    Check(TMS_Passage::BuildPassages(shuffled, 0.5).size() == 4, "passage count independent of input order");
  }

  std::cout << (n_failed ? "FAILED: " : "ALL PASSED: ") << n_failed << " failure(s)\n";
  return n_failed ? 1 : 0;
}
