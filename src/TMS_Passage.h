#ifndef _TMS_PASSAGE_H_
#define _TMS_PASSAGE_H_

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

// Plain-data helpers for the detector-response redesign: stitch the raw edep-sim steps
// inside one scintillator bar into physical "passages", then re-segment a passage onto a
// fixed spatial scale so the detector model no longer depends on where Geant4 happened to
// end a step. Deliberately free of ROOT/geometry/TMS_Event types so the arithmetic can be
// unit-tested in isolation (see app/PassageUtilsTest.cpp). Not yet consumed by the
// simulation itself -- the optical/timing phases that use OpticalDeposit come later.
namespace TMS_Passage {

// One raw edep-sim step inside a bar, as plain data.
struct Segment {
  std::array<double, 3> start;  // mm
  std::array<double, 3> stop;   // mm
  double t_start;               // ns
  double t_stop;                // ns
  double energy;                // MeV
  int bar_key;                  // identifies one physical bar (caller collapses plane/bar/view)
  int trajectory_id;            // contributing trajectory; provenance and passage grouping
};

// A fixed-length slice of a passage carrying the energy that geometrically overlaps it.
struct OpticalDeposit {
  std::array<double, 3> position;                 // mm, midpoint of the slice along the passage
  double t;                                       // ns, interpolated at the midpoint
  double dx;                                      // mm, slice length along the passage
  double energy;                                  // MeV
  std::vector<std::pair<int, double>> shares;     // (trajectory_id, MeV); sums to energy
};

double Length(const Segment& s);

// Groups segments into passages: same bar and same trajectory, ordered by start time, and
// split wherever the next segment starts more than max_gap_mm from where the previous one
// stopped (a particle that left the bar and came back is two passages). Returns index
// groups into `segments`, each in time order; groups are ordered by first start time.
std::vector<std::vector<size_t>> BuildPassages(const std::vector<Segment>& segments,
                                               double max_gap_mm);

// Re-segments one passage (segments in order along the trajectory) into ceil(L/max_bin)
// equal-length bins, assigning each original segment's energy to bins by geometric overlap
// (E_ij = E_j * L_ij / L_j), not by dividing the total equally. Total energy is conserved.
std::vector<OpticalDeposit> Resegment(const std::vector<Segment>& passage, double max_bin_length_mm);

}  // namespace TMS_Passage

#endif
