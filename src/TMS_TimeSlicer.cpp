#include <cmath>
#include <functional>
#include <map>
#include <vector>

#include "TMS_TimeSlicer.h"
#include "TMS_Event.h"
#include "TMS_Hit.h"
#include "TMS_Manager.h"

#include <algorithm>
#include <numeric>
#include <utility>

namespace {

// The energy-window slicing of RunTimeSlicer(), on any subset of hits: label
// each DT-wide time unit with a slice index (0 = none). A slice opens when the
// energy in a sliding window of slidingWindowWidth units reaches thresholdStart,
// runs at least minimumSliceWidth units, and closes once the window's energy
// drops below thresholdEnd. Returns the number of slices found.
int WindowSlices(const std::vector<const TMS_Hit *> &hits, double DT, int nUnits, double thresholdStart,
                 double thresholdEnd, int slidingWindowWidth, int minimumSliceWidth, std::vector<int> &labels) {
  std::vector<double> energy(nUnits, 0.0);
  for (const TMS_Hit *hit : hits) {
    const int index = hit->GetT() / DT;
    if (index >= 0 && index < nUnits) energy[index] += hit->GetE();
  }
  labels.assign(nUnits, 0);
  int minimumIndex = 0, sliceIndex = 1;
  bool inSlice = false;
  for (int i = 0; i < nUnits; i++) {
    double inWindow = 0;
    for (int j = 0; i + j + slidingWindowWidth < nUnits && j < slidingWindowWidth; j++) inWindow += energy[i + j];
    if (!inSlice && inWindow >= thresholdStart) {
      inSlice = true;
      minimumIndex = i + minimumSliceWidth;
    }
    if (inSlice && inWindow < thresholdEnd && i > minimumIndex) {
      inSlice = false;
      for (int j = 0; i + j + slidingWindowWidth < nUnits && j < slidingWindowWidth - 1; j++) labels[i + j] = sliceIndex;
      i += slidingWindowWidth - 1;
      sliceIndex += 1;
    }
    if (inSlice) labels[i] = sliceIndex;
  }
  return sliceIndex - 1 + (inSlice ? 1 : 0);
}

}  // namespace

int TMS_TimeSlicer::SimpleTimeSlicer(TMS_Event &event) {
  int nslices = 1;
  // For now do the simplest thing and divide into N chunks
  int nsliceswithoneormorehit = 0;
  double spill_time = 10000; // ns
  int n_slices_target = 1; //52;
  double dt = spill_time / n_slices_target; // ns
  auto hits = event.GetHitsRaw();
  //std::cout<<"Running time slicer with n="<<hits.size()<<std::endl;
  for (int slice_number = 1; slice_number <= n_slices_target; slice_number++) {
    double start_time = dt * (slice_number - 1);
    double end_time = dt * slice_number;
    int nhitsinslice = 0;
    //for (auto hit : hits) {
    //for (std::vector<TMS_Hit>::iterator it = hits.begin(); it != hits.end(); it++) {
    for (size_t i = 0; i < hits.size(); i++) {
      //auto hit = (*it);
      auto hit = hits[i];
      // todo add back?
      //if (hit.GetPedSup()) continue; // Skip ped supped hits
      double hit_time = hit.GetT();
      // If in time with slice, add this hit to slice
      bool hit_is_in_slice = false;
      if (start_time <= hit_time && hit_time < end_time) hit_is_in_slice = true;
      // Special cases for hit times outside of standard range
      // todo, understand why true hit_time < 0 is possible
      if (slice_number == 1 && hit_time < 0) hit_is_in_slice = true;
      if (slice_number == n_slices_target && hit_time >= spill_time) hit_is_in_slice = true;
      if (hit_is_in_slice) {
        if (hit.GetSlice() != 0) { std::cout<<"Trying to change a hit slice from "<<hit.GetSlice()<<" to "<<slice_number<<", t="<<hit_time<<std::endl; exit(0); }
        //std::cout<<"Trying to change a hit slice from "<<hit.GetSlice()<<" to "<<slice_number<<", t="<<hit_time;
        //hit.SetSlice(slice_number);
        auto hit_pointer = &hits[i];
        hit_pointer->SetSlice(slice_number);
        nhitsinslice++;
        //std::cout<<", Checking new slice number: "<<hit.GetSlice()<<", "<<hits[i].GetSlice()<<std::endl;
      }
    }
    if (nhitsinslice > 0) nsliceswithoneormorehit += 1;
    nslices += 1;
  }
  // Need to explicitly change the raw hits in the event since we're not dealing with pointers
  event.SetHitsRaw(hits);
  //std::cout<<"Found "<<nslices<<" slices. "<<nsliceswithoneormorehit<<" have more than one hit."<<std::endl;
  
  
  auto hits2 = event.GetHitsRaw();
  int n_hits_outside_slice0 = 0;
  int n_hits_inside_slice0 = 0;
  for (auto hit : hits2) {
    //if (hit.GetSlice() != 0) std::cout<<"Checking new slice number: "<<hit.GetSlice()<<std::endl;
    if (hit.GetSlice() != 0) n_hits_outside_slice0 += 1;
    if (hit.GetSlice() == 0) n_hits_inside_slice0 += 1;
    if (hit.GetSlice() == 0) std::cout<<"Hit in slice 0, T="<<hit.GetT()<<std::endl;
  }
  //std::cout<<"Found "<<n_hits_outside_slice0<<" hits with slice number != 0, and "<<n_hits_inside_slice0<<" inside slice 0"<<std::endl;
  
  
  event.SetNSlices(nslices);
  return nslices;
}

int TMS_TimeSlicer::RunTimeSlicer(TMS_Event &event) {
  int nslices = 1;
  // Sort by T so slices are easier to find
  event.SortHits(TMS_Hit::SortByT);
  
  
  bool RunTimeSlicer = TMS_Manager::GetInstance().Get_Reco_TIME_RunTimeSlicer();
  bool RunSimpleTimeSlicer = TMS_Manager::GetInstance().Get_Reco_TIME_RunSimpleTimeSlicer();
  const double SPILL_LENGTH = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerMaxTime();
  if (!RunTimeSlicer) {
    event.SetNSlices(nslices);
    event.AddTimeSliceInformation({std::make_pair(0.0, SPILL_LENGTH)});
    return nslices;
  }
  if (RunTimeSlicer && RunSimpleTimeSlicer) nslices = SimpleTimeSlicer(event);
  if (RunTimeSlicer && !RunSimpleTimeSlicer && TMS_Manager::GetInstance().Get_RECO_TIME_PerViewSlicing())
    return PerViewTimeSlicer(event);
  if (RunTimeSlicer && !RunSimpleTimeSlicer) {
    // Here are all the constants
    double threshold1 = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerThresholdStart();
    double threshold2 = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerThresholdEnd();
    int sliding_window_width = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerEnergyWindowInUnits();
    int minimum_slice_width = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerMinimumSliceWidthInUnits();
    const double DT = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerSliceUnit();
    const int NUMBER_OF_SLICES = std::ceil(SPILL_LENGTH / DT);

    // First initialize an array of energy and slice labels
    std::vector<double> energy_slices(NUMBER_OF_SLICES, 0.0);
    std::vector<int> time_slices(NUMBER_OF_SLICES, 0);
    
    // Add all hit energy to array
    auto hits = event.GetHitsRaw();
    for (auto hit : hits) {
      // Only include hits that are not pedestal subtracted
      if (!hit.GetPedSup()) {
        int index = hit.GetT() / DT;
        // Make sure we're within bounds, and add energy
        if (index >= 0 && index < NUMBER_OF_SLICES) energy_slices[index] += hit.GetE();
      }
    }
    
    // Now make a sliding window;
    // After starting a slice, go at least this far
    // This allows for a slice that's at least as wide as minimum_slice_width
    int minimum_index = 0;
    int slice_index = 1;
    bool in_slice = false;
    for (int i = 0; i < NUMBER_OF_SLICES; i++) {
      double energy_in_window = 0;
      for (int j = 0; i + j + sliding_window_width < NUMBER_OF_SLICES && j < sliding_window_width; j++) 
        energy_in_window += energy_slices[i + j];
      
      // Reached threshold to start making slice
      if (!in_slice && energy_in_window >= threshold1) { 
        //std::cout<<"Starting slice at i="<<i<<", energy_in_window="<<energy_in_window<<std::endl;
        in_slice = true;
        minimum_index = i + minimum_slice_width;
      }
      
      // Reached below threshold then stop recording slice
      // But only if we reached minimum_index
      if (in_slice && energy_in_window < threshold2 && i > minimum_index) {
        //std::cout<<"Ending slice at i="<<i<<", energy_in_window="<<energy_in_window<<std::endl;
        in_slice = false;
        // Finish writing that window. 
        // So that time_slices[i:i+sliding_window_width-1] all are in slice_index
        // This will mean that time_slices[i+sliding_window_width] will be first array index not in time slice
        for (int j = 0; i + j + sliding_window_width < NUMBER_OF_SLICES && j < sliding_window_width - 1; j++) 
          time_slices[i + j] = slice_index;
        i += sliding_window_width - 1;
        slice_index += 1;
      }
      
      // If in a slice, record that index.
      if (in_slice) {
        time_slices[i] = slice_index;
      }
    }
    // Write out the remaining slices
    nslices = slice_index;
    if (nslices > 0 && false) {
      std::cout<<"Found "<<nslices<<" slices. ";//<<std::endl;
      double occupancy = 0;
      double energy_outside_slice_0 = 0;
      double energy_total = 0;
      for (int i = 0; i < NUMBER_OF_SLICES; i++) {
        double energy = energy_slices[i];
        if (time_slices[i] != 0) { 
          occupancy += 1; 
          energy_outside_slice_0 += energy;
        }
        energy_total += energy;
      }
      occupancy /= NUMBER_OF_SLICES;
      std::cout<<occupancy<<" of time slice units have different slice labels. "; //<<std::endl;
      double percent = 100 * energy_outside_slice_0 / energy_total;
      std::cout<<percent<<"% of the energy"<<std::endl;
    }
    
    // Make a measurement of each slice location
    std::vector<std::pair<double, double>> slice_bounds;
    // Want the indices to match up, so add pair for slice 0
    slice_bounds.push_back(std::make_pair(0, SPILL_LENGTH));
    int prev_slice = -1;
    bool have_prev_slice = false;
    double slice_start_time = -999;
    double slice_end_time = -999;
    for (int i = 0; i < NUMBER_OF_SLICES; i++) {
      int current_slice = time_slices[i];
      double current_slice_time = i * DT;
      if (current_slice != prev_slice) {
        // Write out the current slice and then start a new slice
        if (have_prev_slice) {
          // Last slice end time = window before this one, so subtract DT
          slice_end_time = current_slice_time - DT;
          slice_bounds.push_back(std::make_pair(slice_start_time, slice_end_time));
          have_prev_slice = false;
        }
        
        // Only start tracking the next slice if not equal to zero
        if (!have_prev_slice && current_slice != 0) {
          prev_slice = current_slice;
          slice_end_time = current_slice_time;
          slice_start_time = current_slice_time;
          have_prev_slice = true;
        }
      }
    }
    
    // Finally assign hits based on slice
    std::vector<TMS_Hit> changed_hits;
    for (auto hit : hits) {
      int index = hit.GetT() / DT;
      int slice = 0; 
      // Make sure we're within bounds
      if (index >= 0 && index < NUMBER_OF_SLICES) slice = time_slices[index];
      hit.SetSlice(slice);
      changed_hits.push_back(hit);
    }
    event.SetHitsRaw(changed_hits);
    event.AddTimeSliceInformation(slice_bounds);
  }
  return nslices;
}

// Per-view slicing ([Recon.Time] PerViewSlicing).
//
// Why (2026-09-27, 15 MiniProdN5 files): a hit's time includes the light's
// travel along its bar to the readout, and the two views' bars run in
// different directions -- the x-measuring hits of a muon are delayed by
// ~6-8 ns per meter of its height, the y-measuring ones barely at all -- so a
// muon's two views can sit 20-40 ns apart. With slices ~55 ns wide, a slice
// boundary between them split 2.8% of all muons by view (all x-measuring hits
// in one slice, all y-measuring in the next); neither slice can pair them into
// space points, and ~90% of those muons were lost (2.1% of ND-LAr muons, a
// quarter of all the missed ones). Within one view the delays are consistent.
//
// So: run the energy-window slicer on each view's hits separately (thresholds
// scaled per view, below), then link each view slice to the other view's slice
// it overlaps most in time (widened by PerViewMatchToleranceNs) among those
// overlapping it in z (within PerViewMatchZMarginMM): the connected groups are
// the slices. Two unrelated interactions overlapping in both time and z can end
// up in one slice -- the trackers still separate them in 3D.
//
// Files 1-15 (best-link matching, automatic thresholds) vs. RunTimeSlicer:
// Cluster3D found 91.1 -> 93.0% of ND-LAr muons, single-view slices 3.1 ->
// 1.6%, but 83 vs. 110 slices per spill, so more slices hold >1 interaction.
int TMS_TimeSlicer::PerViewTimeSlicer(TMS_Event &event) {
  TMS_Manager &manager = TMS_Manager::GetInstance();
  const double SPILL_LENGTH = manager.Get_RECO_TIME_TimeSlicerMaxTime();
  const double DT = manager.Get_RECO_TIME_TimeSlicerSliceUnit();
  const int nUnits = std::ceil(SPILL_LENGTH / DT);
  const double scaleSetting = manager.Get_RECO_TIME_PerViewThresholdScale();
  const bool bestLinkMatching = manager.Get_RECO_TIME_PerViewBestLinkMatching();
  const int windowWidth = manager.Get_RECO_TIME_TimeSlicerEnergyWindowInUnits();
  const int minimumWidth = manager.Get_RECO_TIME_TimeSlicerMinimumSliceWidthInUnits();
  const double tolerance = manager.Get_RECO_TIME_PerViewMatchToleranceNs();
  const double zMargin = manager.Get_RECO_TIME_PerViewMatchZMarginMM();

  std::vector<TMS_Hit> hits = event.GetHitsRaw();
  // View 0: y-measuring hits (X-type bars); view 1: everything else (x-measuring).
  auto viewOf = [](const TMS_Hit &hit) { return hit.GetBar().GetBarType() == TMS_Bar::kXBar ? 0 : 1; };
  std::vector<const TMS_Hit *> viewHits[2];
  for (const TMS_Hit &hit : hits)
    if (!hit.GetPedSup()) viewHits[viewOf(hit)].push_back(&hit);

  // Per-view slices: labels per time unit, and each slice's time and z range.
  struct ViewSlice {
    double t0 = 1e30, t1 = -1e30, z0 = 1e30, z1 = -1e30;
  };
  // Each view's thresholds: the full-detector ones scaled by PerViewThresholdScale,
  // or -- if that is <= 0 -- by the view's own share of this spill's hit
  // energy. The views are not equal: y-measuring hits come from 30 of the 82
  // planes, so a muon leaves ~37% of its energy there, and a common factor of
  // 0.5 under-triggered the y view (2026-09-27: 41% of slices single-view).
  double viewEnergy[2] = {0.0, 0.0};
  for (int v = 0; v < 2; ++v)
    for (const TMS_Hit *hit : viewHits[v]) viewEnergy[v] += hit->GetE();
  std::vector<int> labels[2];
  std::vector<ViewSlice> slices[2];
  for (int v = 0; v < 2; ++v) {
    const double total = viewEnergy[0] + viewEnergy[1];
    const double scale = scaleSetting > 0.0 ? scaleSetting : (total > 0.0 ? viewEnergy[v] / total : 0.5);
    const double thresholdStart = scale * manager.Get_RECO_TIME_TimeSlicerThresholdStart();
    const double thresholdEnd = scale * manager.Get_RECO_TIME_TimeSlicerThresholdEnd();
    const int n = WindowSlices(viewHits[v], DT, nUnits, thresholdStart, thresholdEnd, windowWidth, minimumWidth, labels[v]);
    slices[v].assign(n + 1, ViewSlice());  // index 0 unused
    for (int i = 0; i < nUnits; ++i) {
      const int l = labels[v][i];
      if (l <= 0 || l > n) continue;
      slices[v][l].t0 = std::min(slices[v][l].t0, i * DT);
      slices[v][l].t1 = std::max(slices[v][l].t1, (i + 1) * DT);
    }
    for (const TMS_Hit *hit : viewHits[v]) {
      const int index = hit->GetT() / DT;
      if (index < 0 || index >= nUnits) continue;
      const int l = labels[v][index];
      if (l <= 0 || l > n) continue;
      slices[v][l].z0 = std::min(slices[v][l].z0, hit->GetZ());
      slices[v][l].z1 = std::max(slices[v][l].z1, hit->GetZ());
    }
  }

  // For every view slice, find its best partner in the other view: the one it
  // overlaps most in time (after widening by the tolerance, and requiring z
  // overlap). Linking every overlapping pair instead let pairs chain (x1-y1-x2-
  // y2...) and glued successive interactions into one slice: half as many
  // slices on the first spill tried (46 vs 99). The goal is one slice per
  // interaction.
  const int offset = slices[0].size();
  std::vector<int> parent(offset + slices[1].size());
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int a) { return parent[a] == a ? a : parent[a] = find(parent[a]); };
  auto overlap = [&](const ViewSlice &x, const ViewSlice &y) {
    if (x.z0 > x.z1 || y.z0 > y.z1) return 0.0;  // no hits
    if (!(x.z0 - zMargin < y.z1 && y.z0 < x.z1 + zMargin)) return 0.0;
    return std::max(0.0, std::min(x.t1 + tolerance, y.t1) - std::max(x.t0 - tolerance, y.t0));
  };
  std::vector<int> bestY(slices[1].size(), -1), bestX(slices[0].size(), -1);
  std::vector<double> bestYOverlap(slices[1].size(), 0.0), bestXOverlap(slices[0].size(), 0.0);
  for (std::size_t a = 1; a < slices[1].size(); ++a)
    for (std::size_t b = 1; b < slices[0].size(); ++b) {
      const double o = overlap(slices[1][a], slices[0][b]);
      if (o <= 0.0) continue;
      if (o > bestYOverlap[a]) { bestYOverlap[a] = o; bestY[a] = static_cast<int>(b); }
      if (o > bestXOverlap[b]) { bestXOverlap[b] = o; bestX[b] = static_cast<int>(a); }
    }
  // PerViewBestLinkMatching (default): every view slice adds its own best link,
  // so an interaction whose view was cut into several slices can reassemble;
  // one link per slice chains far less than linking every overlap. Otherwise
  // only mutual best pairs are linked -- tested 2026-09-27 and worse: ~40% of
  // slices were left single-view and 18% of muons had their views split.
  for (std::size_t a = 1; a < slices[1].size(); ++a) {
    const int b = bestY[a];
    if (b > 0 && (bestLinkMatching || bestX[b] == static_cast<int>(a))) parent[find(offset + static_cast<int>(a))] = find(b);
  }
  if (bestLinkMatching)
    for (std::size_t b = 1; b < slices[0].size(); ++b) {
      const int a = bestX[b];
      if (a > 0) parent[find(static_cast<int>(b))] = find(offset + a);
    }

  // Final slices: the groups, numbered 1.. in order of their start time.
  std::map<int, std::pair<double, double>> groupTimes;  // root -> (start, end)
  for (int v = 0; v < 2; ++v)
    for (std::size_t l = 1; l < slices[v].size(); ++l) {
      if (slices[v][l].z0 > slices[v][l].z1) continue;
      const int root = find(v == 0 ? static_cast<int>(l) : offset + static_cast<int>(l));
      auto it = groupTimes.find(root);
      if (it == groupTimes.end()) groupTimes[root] = {slices[v][l].t0, slices[v][l].t1};
      else {
        it->second.first = std::min(it->second.first, slices[v][l].t0);
        it->second.second = std::max(it->second.second, slices[v][l].t1);
      }
    }
  std::vector<std::pair<double, int>> order;
  for (const auto &kv : groupTimes) order.push_back({kv.second.first, kv.first});
  std::sort(order.begin(), order.end());
  std::map<int, int> sliceOfRoot;
  std::vector<std::pair<double, double>> sliceBounds;
  sliceBounds.push_back(std::make_pair(0.0, SPILL_LENGTH));  // slice 0
  for (std::size_t k = 0; k < order.size(); ++k) {
    sliceOfRoot[order[k].second] = static_cast<int>(k) + 1;
    sliceBounds.push_back(groupTimes[order[k].second]);
  }

  // Each hit takes the final slice of its own view's slice at its time
  // (pedestal-suppressed hits too, as RunTimeSlicer labels every hit).
  for (TMS_Hit &hit : hits) {
    const int v = viewOf(hit);
    const int index = hit.GetT() / DT;
    int slice = 0;
    if (index >= 0 && index < nUnits) {
      const int l = labels[v][index];
      if (l > 0 && l < static_cast<int>(slices[v].size()) && slices[v][l].z0 <= slices[v][l].z1)
        slice = sliceOfRoot[find(v == 0 ? l : offset + l)];
    }
    hit.SetSlice(slice);
  }
  event.SetHitsRaw(hits);
  event.AddTimeSliceInformation(sliceBounds);
  const int nslices = static_cast<int>(order.size()) + 1;
  return nslices;
}
