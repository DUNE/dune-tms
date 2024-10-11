#include "TMS_TimeSlicer.h"
#include "TMS_Event.h"
#include "TMS_Hit.h"
#include "TMS_Manager.h"

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
  if (RunTimeSlicer && RunSimpleTimeSlicer) nslices = SimpleTimeSlicer(event);
  if (RunTimeSlicer && !RunSimpleTimeSlicer) {
    // Here are all the constants
    double threshold1 = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerThresholdStart();
    double threshold2 = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerThresholdEnd();
    int sliding_window_width = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerEnergyWindowInUnits();
    int minimum_slice_width = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerMinimumSliceWidthInUnits();
    const double SPILL_LENGTH = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerMaxTime();
    const double DT = TMS_Manager::GetInstance().Get_RECO_TIME_TimeSlicerSliceUnit();
    const int NUMBER_OF_SLICES = std::ceil(SPILL_LENGTH / DT);
    
    // First initialize an array of energy and slice labels
    double energy_slices[NUMBER_OF_SLICES];
    int time_slices[NUMBER_OF_SLICES];
    for (int i = 0; i < NUMBER_OF_SLICES; i++) {
      energy_slices[i] = 0;
      time_slices[i] = 0;
    }
    
    // Add all hit energy to array
    auto hits = event.GetHitsRaw();
    for (auto hit : hits) {
      int index = hit.GetT() * DT;
      // Make sure we're within bounds, and add energy
      if (index >= 0 && index < NUMBER_OF_SLICES) energy_slices[index] += hit.GetE();
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
      double current_slice_time = i / DT;
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
      int index = hit.GetT() * DT;
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
