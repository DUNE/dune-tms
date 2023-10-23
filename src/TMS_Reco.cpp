#include "TMS_Reco.h"

TMS_TrackFinder::TMS_TrackFinder() :
  nIntercept(TMS_Manager::GetInstance().Get_Reco_HOUGH_NInter()),
  nSlope(TMS_Manager::GetInstance().Get_Reco_HOUGH_NSlope()),
  InterceptMin(TMS_Manager::GetInstance().Get_Reco_HOUGH_MinInterp()),
  InterceptMax(TMS_Manager::GetInstance().Get_Reco_HOUGH_MaxInterp()),
  SlopeMin(TMS_Manager::GetInstance().Get_Reco_HOUGH_MinSlope()),
  SlopeMax(TMS_Manager::GetInstance().Get_Reco_HOUGH_MaxSlope()),
  InterceptWidth((InterceptMax-InterceptMin)/nIntercept),
  SlopeWidth((SlopeMax-SlopeMin)/nSlope),
  // Max z for us to do Hough in, here choose transition layer
  zMinHough(TMS_Const::TMS_Thin_Start),
  //zMaxHough(TMS_Const::TMS_Thick_Start),
  zMaxHough(TMS_Const::TMS_Thick_End),
  nMaxHough(TMS_Manager::GetInstance().Get_Reco_HOUGH_MaxHough()),
  nThinCont(10),
  nHits_Tol(TMS_Manager::GetInstance().Get_Reco_HOUGH_HitMult()),
  // Minimum number of hits required to run track finding
  nMinHits(TMS_Manager::GetInstance().Get_Reco_MinHits()),
  // Maximum number of merges for one hit
  nMaxMerges(1),
  MinDistHough(TMS_Manager::GetInstance().Get_Reco_HOUGH_MinDist()), // Minimum distance for Hough in mm
  // Is AStar greedy
  IsGreedy(TMS_Manager::GetInstance().Get_Reco_ASTAR_IsGreedy()),
  // Use DBSCAN clustering after track finding
  UseClustering(TMS_Manager::GetInstance().Get_Reco_Clustering())
{
  // Apply the maximum Hough transform to the zx not zy: all bending happens in zx
  Accumulator = new int*[nSlope];
  for (int i = 0; i < nSlope; ++i) {
    Accumulator[i] = new int[nIntercept];
  }

  // Keep track of the tracking efficiency
  Efficiency = new TH1D("Efficiency", "Efficiency;T_{#mu} (GeV); Efficiency", 30, 0, 6);
  Total = new TH1D("Total", "Total;T_{#mu} (GeV); Total", 30, 0, 6);

  HoughLineOne = new TF1("LinearHough", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineOne->SetLineStyle(kDashed);
  HoughLineOne->SetLineColor(kMagenta-9);

  HoughLineOther = new TF1("LinearHough2", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineOther->SetLineStyle(kDashed);
  HoughLineOther->SetLineColor(kMagenta-8);

  DBSCAN.SetEpsilon(TMS_Manager::GetInstance().Get_Reco_DBSCAN_Epsilon());
  DBSCAN.SetMinPoints(TMS_Manager::GetInstance().Get_Reco_DBSCAN_MinPoints());

  // Set up the tracker algorithm
  std::string trackname = TMS_Manager::GetInstance().Get_Reco_TrackMethod();
  if      (trackname == "Hough")  kTrackMethod = TrackMethod::kHough;
  else if (trackname == "AStar")  kTrackMethod = TrackMethod::kAStar;
  else if (trackname == "DBSCAN") kTrackMethod = TrackMethod::kDBSCAN;
  else {
    std::cerr << "Invalid track reconstruction method provided to TMS reconstruction" << std::endl;
    std::cerr << "You provided: " << trackname << std::endl;
    std::cerr << "Options are: Hough, AStar, DBSCAN" << std::endl;
    kTrackMethod = TrackMethod::kUnknown;
    throw;
  }

  std::string heuristicname = TMS_Manager::GetInstance().Get_Reco_ASTAR_CostMetric();
  if      (heuristicname == "Euclidean") kHeuristic = HeuristicType::kEuclidean;
  else if (heuristicname == "Manhattan") kHeuristic = HeuristicType::kManhattan;
  else if (heuristicname == "DetectorZ") kHeuristic = HeuristicType::kDetectorZ;
  else if (heuristicname == "DetectorNotZ") kHeuristic = HeuristicType::kDetectorNotZ;
  else {
    std::cerr << "Invalid AStar heuristic method provided to TMS reconstruction" << std::endl;
    std::cerr << "You provided: " << heuristicname << std::endl;
    std::cerr << "Options are: Euclidean, Manhattan, DetectorZ, DetectorNotZ" << std::endl;
    kHeuristic = HeuristicType::kUnknown;
    throw;
  }

  std::cout << "Using " << trackname << " for main track finding reconstruction" << std::endl;
  std::cout << "Using AStar to clean up HoughTrack? " << TMS_Manager::GetInstance().Get_Reco_HOUGH_RunAStar() << std::endl;
  std::cout << "Using DBSCAN for clustering after main track finding? " << UseClustering << std::endl;
  if (kTrackMethod == TrackMethod::kAStar) {
    std::cout << "Using " << heuristicname << " for A* heuristic cost calculation" << std::endl;
  }


}

void TMS_TrackFinder::ClearClass() {

  // Check through the Houghlines
  for (auto &i: HoughLinesOne) {
    delete i.second;
  }
  for (auto &i: HoughLinesOther) {
    delete i.second;
  }

  // Reset the candidate vector
  Candidates.clear();
  RawHits.clear();
  TotalCandidates.clear();
  HoughLinesOne.clear();
  HoughLinesOther.clear();
  HoughLinesOne_Upstream.clear();
  HoughLinesOther_Upstream.clear();
  HoughLinesOne_Downstream.clear();
  HoughLinesOther_Downstream.clear();
  HoughCandidatesOne.clear();
  HoughCandidatesOther.clear();
  ClusterCandidatesOne.clear();
  ClusterCandidatesOther.clear();
  TrackLengthOne.clear();
  TrackLengthOther.clear();
  TrackEnergyOne.clear();
  TrackEnergyOther.clear();
}

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
  }
  return nslices;
}

// The generic track finder
void TMS_TrackFinder::FindTracks(TMS_Event &event) {

  ClearClass();

  // Get the raw unmerged and untracked hits
  RawHits = event.GetHits();
  //std::cout<<"Working on raw hits with n="<<RawHits.size()<<std::endl;
  //if (RawHits.size() > 0) std::cout<<"Slice "<<slice<<" has "<<RawHits.size()<<" hits."<<std::endl;
  
  
  double min_time = 1e9;
  double max_time = -1e9;
  int slice = event.GetSliceNumber();
  for (auto hit : RawHits) {
    if (hit.GetPedSup()) std::cout<<"Raw hits, found a ped supped hit in slice"<<std::endl;
    min_time = std::min(min_time, hit.GetT());
    max_time = std::max(max_time, hit.GetT());
  }
  //if (max_time - min_time > (10000.0 / 52)) std::cout<<"In Reco RawHits, found a time range larger than expected: "<<(max_time - min_time)<<", min="<<min_time<<", max="<<max_time<<", slice="<<slice<<std::endl;

  // Clean hits (check duplicate hits, and energy threshold)
  CleanedHits = CleanHits(RawHits);
  // Require N hits after cleaning
  if (CleanedHits.size() < nMinHits) return;
  
  double clean_min_time = 1e9;
  double clean_max_time = -1e9;
  int n_in_slice = 0;
  int n_in_wrong_slice = 0;
  for (auto hit : CleanedHits) {
    if (hit.GetPedSup()) std::cout<<"Cleaned hits, found a ped supped hit in slice"<<std::endl;
    clean_min_time = std::min(clean_min_time, hit.GetT());
    clean_max_time = std::max(clean_max_time, hit.GetT());
    if (hit.GetSlice() == slice) n_in_slice += 1;
    if (hit.GetSlice() != slice) n_in_wrong_slice += 1;
  }
  //if (clean_max_time - clean_min_time > (10000.0 / 52)) std::cout<<"In Reco CleanedHits, found a time range larger than expected: "<<(clean_max_time - clean_min_time)<<", min="<<clean_min_time<<", max="<<clean_max_time<<", slice="<<slice<<std::endl;
  if (n_in_wrong_slice > 0) std::cout<<"Found "<<n_in_wrong_slice<<" in wrong slice and "<<n_in_slice<<" in correct slice for "<<slice<<std::endl;

#ifdef DEBUG
  std::cout << "Raw hits: " << RawHits.size() << std::endl;
  std::cout << "Cleaned hits: " << CleanedHits.size() << std::endl;
#endif

  // separate planes into different groups
  // 3 degree stereo -> tilted into +3 degree in one group and into -3 degree in other group
  // 90 degree rotated -> vertical layers in one group and horizontal layers in other group
  
  // for loop over hits
  for (auto hit : CleanedHits) {
  // get PlaneNumber per hit hit.GetBar().GetPlaneNumber()
  // sorting hits into orientation groups
    if (hit.GetBar().GetPlaneNumber() % LayerOrientation) OneHitGroup.push_back(hit); // add hit to one group
    else if (!(hit.GetBar().GetPlaneNumber() % LayerOrientation)) OtherHitGroup.push_back(hit); // add hit to other group
  }

  // Hough transform
  if (kTrackMethod == TrackMethod::kHough) {
    // Do we first run clustering algorithm to separate hits, then hand off to A*?
    if (TMS_Manager::GetInstance().Get_Reco_HOUGH_FirstCluster()) {
      // Let's run a DBSCAN first to cluster up, then run Hough transform on clusters
      //std::vector<std::vector<TMS_Hit> > DBScanCandidates = FindClusters(CleanedHits);
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesOne = FindClusters(OneHitGroup);
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesOther = FindClusters(OtherHitGroup);
      // Hand over each cluster from DBSCAN to a Hough transform
      for (std::vector<std::vector<TMS_Hit> >::iterator it = DBScanCandidatesOne.begin(); it != DBScanCandidatesOne.end(); ++it) {
        std::vector<TMS_Hit> hits = *it;
        std::vector<std::vector<TMS_Hit> > Lines = HoughTransform(hits, 1);
        for (auto jt = Lines.begin(); jt != Lines.end(); ++jt) {
          HoughCandidatesOne.emplace_back(std::move(*jt));
        }
      }
      for (std::vector<std::vector<TMS_Hit> >::iterator it = DBScanCandidatesOther.begin(); it != DBScanCandidatesOther.end(); ++it) {
        std::vector<TMS_Hit> hits = *it;
	std::vector<std::vector<TMS_Hit> > Lines = HoughTransform(hits, 2);
	for (auto jt = Lines.begin(); jt != Lines.end(); ++jt) {
	  HoughCandidatesOther.emplace_back(std::move(*jt));
	}
      }
    } else {
      //HoughCandidates = HoughTransform(CleanedHits);
      HoughCandidatesOne = HoughTransform(OneHitGroup, 1);
      HoughCandidatesOther = HoughTransform(OtherHitGroup, 2);
    }
  } else if (kTrackMethod == TrackMethod::kAStar) {
    //BestFirstSearch(CleanedHits);
    BestFirstSearch(OneHitGroup, 1);
    BestFirstSearch(OtherHitGroup, 2);
  }

  //std::vector<TMS_Hit> Masked = CleanedHits;
  std::vector<TMS_Hit> MaskedOne = OneHitGroup;
  std::vector<TMS_Hit> MaskedOther = OtherHitGroup;
  // Loop over the Hough candidates
  for (auto Lines: HoughCandidatesOne) {
#ifdef DEBUG
    std::cout << "Masked (one) size bef: " << MaskedOne.size() << std::endl;
#endif
    MaskHits(MaskedOne, Lines);
#ifdef DEBUG
    std::cout << "Masked (one) size aft: " << MaskedOne.size() << std::endl;
#endif
  }

  for (auto Lines: HoughCandidatesOther) {
#ifdef DEBUG
    std::cout << "Masked (other) size bef: " << MaskedOther.size() << std::endl;
#endif
    MaskHits(MaskedOther, Lines);
#ifdef DEBUG
    std::cout << "Masked (other) size aft: " << MaskedOther.size() << std::endl;
#endif
  }
  
#ifdef DEBUG
  std::cout << "Masked (one) hits: " << MaskedOne.size() << std::endl;
  std::cout << "Masked (other) hits: " << MaskedOther.size() << std::endl;
#endif
  
  // Now we've got our tracks, refit the upstream and downstream separately with the Hough transform
  int lineno = 0;
  //std::cout << "Event " << event.GetEventNumber() << std::endl;
  for (auto Lines: HoughCandidatesOne) {
    //std::cout << "line  " << lineno << std::endl;
    std::pair<bool, TF1*> houghline = HoughLinesOne[lineno];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
        fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      //std::cout << "Old slope: " << houghline.second->GetParameter(1) << std::endl;
      //std::cout << "New slope: " << slope << std::endl;
      //std::cout << "Old intercept: " << houghline.second->GetParameter(0) << std::endl;
      //std::cout << "New intercept: " << intercept << std::endl;

      HoughLinesOne[lineno].second->SetParameter(0, intercept);
      HoughLinesOne[lineno].second->SetParameter(1, slope);
    }

    // The number of hits in this track, take 20% and call upstream and dowstream segments
    int nrescanhits = 0.3*Lines.size()+1;
    // If there are only a few hits, use all of them
    if (nrescanhits < 5) nrescanhits = Lines.size();
    std::vector<TMS_Hit> upstream;
    std::vector<TMS_Hit> downstream;
    for (int i = 0; i < nrescanhits; ++i) {
      upstream.push_back(Lines[Lines.size()-1-i]);
      downstream.push_back(Lines[i]);
      //std::cout << upstream.back().GetZ() << " " << downstream.back().GetZ() << std::endl;
    }
    //std::cout << "nhits: " << Lines.size() << std::endl;
    //std::cout << nrescanhits << std::endl;
    //std::cout << upstream.size() << std::endl;
    //std::cout << downstream.size() << std::endl;

    double upstreamslope, upstreamintercept = 0;
    double downstreamslope, downstreamintercept = 0;
    GetHoughLine(upstream, upstreamslope, upstreamintercept);
    GetHoughLine(downstream, downstreamslope, downstreamintercept);

    std::pair<double, double> upstreamline = std::pair<double,double>(upstreamintercept, upstreamslope);
    std::pair<double, double> downstreamline = std::pair<double,double>(downstreamintercept, downstreamslope);

    HoughLinesOne_Upstream.push_back(upstreamline);
    HoughLinesOne_Downstream.push_back(downstreamline);

    lineno++;
  }

  lineno = 0;
  for (auto Lines: HoughCandidatesOther) {
    std::pair<bool, TF1*> houghline = HoughLinesOther[lineno];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
	fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      HoughLinesOther[lineno].second->SetParameter(0, intercept);
      HoughLinesOther[lineno].second->SetParameter(1, slope);
    }
    
    // The number of hits in this track, trake 20% and call upstream and downstream segments
    int nrescanhits = 0.3*Lines.size()+1;
    // If there are only a few hits, use all of them
    if (nrescanhits < 5) nrescanhits = Lines.size();
    std::vector<TMS_Hit> upstream;
    std::vector<TMS_Hit> downstream;
    for (int i = 0; i < nrescanhits; ++i) {
      upstream.push_back(Lines[Lines.size()-1-i]);
      downstream.push_back(Lines[i]);
   }

   double upstreamslope, upstreamintercept = 0;
   double downstreamslope, downstreamintercept = 0;
   GetHoughLine(upstream, upstreamslope, upstreamintercept);
   GetHoughLine(downstream, downstreamslope, downstreamintercept);

   std::pair<double, double> upstreamline = std::pair<double,double>(upstreamintercept, upstreamslope);
   std::pair<double, double> downstreamline = std::pair<double,double>(downstreamintercept, downstreamslope);

   HoughLinesOther_Upstream.push_back(upstreamline);
   HoughLinesOther_Downstream.push_back(downstreamline);

   lineno++;
  }

  // Try finding some clusters after the Hough Transform
  if (UseClustering) {
    ClusterCandidatesOne = FindClusters(MaskedOne);
    ClusterCandidatesOther = FindClusters(MaskedOther);
  }

  // Let's try to find a vertex now, just looking at most upstream point, or if there are multiple tracks let's see where they intersect
  //if (nLines > 0) {
    //if (nLines == 1) Vertex =;
    //else {
    //}
  //} else {
    //Vertex = -999;
  //}

  // Now calculate the track length for each track
  CalculateTrackLengthOne(HoughCandidatesOne);
  CalculateTrackEnergyOne(HoughCandidatesOne);
  CalculateTrackLengthOther(HoughCandidatesOther);
  CalculateTrackEnergyOther(HoughCandidatesOther);

  // For future probably want to move track candidates into the TMS_Event class
  //EvaluateTrackFinding(event);

  // Find if the event may have started outside the TMS
  // Look at the first hits of each of the Hough lines
  // Also check that the hits are continuous

  // Skip the Kalman filter for now
  return;

  // Now have the TotalCandidates filled
  // Start some reconstruction chain
  for (auto &i : TotalCandidates) {
    // Get the xz and yz hits
    std::vector<TMS_Hit> xz_hits = ProjectHits(i, TMS_Bar::kYBar);
    size_t nHits = xz_hits.size();
    if (nHits < 1) continue; 
    KalmanFitter = TMS_Kalman(xz_hits);

    /*
       std::vector<TMS_Hit> yz_hits = ProjectHits(i, TMS_Bar::kXBar);
       std::cout << "yz hits: " << yz_hits.size() << std::endl;
       KalmanFitter = TMS_Kalman(yz_hits);
       */
  }
    

}

// Let's try to evaluate the goodness of the trackfinding
void TMS_TrackFinder::EvaluateTrackFinding(TMS_Event &event) {
  // Have the TMS_Event, which has a vector of TMS_TrueParticles, which has a vector of hits from each particle (so target the muon for instance)

  // Get the true muon
  std::vector<TMS_TrueParticle> TrueParticles = event.GetTrueParticles();
  TMS_TrueParticle muon;
  bool FoundMuon = false;
  for (auto &particle: TrueParticles) {
    if (particle.GetPDG() == 13 && 
        particle.GetParent() == -1) {
      muon = particle;
      FoundMuon = true;
    }
  }

  if (!FoundMuon) return;

  // Now compare this to our tracks and see how many of the muon hits are included
  // Loop over the candidate tracks that have been built
  int nTrueHits = muon.GetPositionPoints().size();
  if (nTrueHits == 0) return;

  int nFoundHits = 0;
  // Loop over all the muon hits
  for (auto &TrueHit: muon.GetPositionPoints()) {
    // Find if the true hit sits inside this hit's bar
    double x_true = TrueHit.X();
    //double y_true = TrueHit.Y();
    double z_true = TrueHit.Z();
    //double t_true = TrueHit.T();
    //std::cout << "Finding " << x_true << " " << z_true << std::endl;
    for (auto &Track: TotalCandidates) {
      // Loop over each hit in each track
      for (auto &Hit: Track) {
        // Get the bar of this hit
        TMS_Bar bar = Hit.GetBar();
        if (bar.Contains(x_true, z_true)) {
          nFoundHits++;
        }
      }
    }
  }

  double efficiency = double(nFoundHits)/nTrueHits;
  double muonmom = muon.GetInitialMomentum().Mag();

  double muonke = sqrt(muonmom*muonmom+TMS_KinConst::mass_mu*TMS_KinConst::mass_mu) - TMS_KinConst::mass_mu;
  muonke /= 1E3;
  if (efficiency > 0) {
    Efficiency->Fill(muonke, efficiency);
    Total->Fill(muonke);
  }
}

// Calculate the total track energy
void TMS_TrackFinder::CalculateTrackEnergyOne() {
  // Look at the reconstructed tracks
  if (HoughCandidatesOne.size() == 0) return;

  // Loop over each Hough Candidate and find the track length
  for (auto it = HoughCandidatesOne.begin(); it != HoughCandidatesOne.end(); ++it) {
    double total = 0;
    // Sort by increasing z
    std::sort((*it).begin(), (*it).end(), TMS_Hit::SortByZInc);
    for (auto hit = (*it).begin(); hit != (*it).end(); ++hit) {
      total += (*hit).GetE();
    }
    TrackEnergyOne.push_back(total);
  }
}

void TMS_TrackFinder::CalculateTrackEnergyOther() {
  //Look at the reconstructed tracks
  if (HoughCandidatesOther.size() == 0) return;

  //Lopp over each Hough Candidate and find the track length
  for (auto it = HoughCandidatesOther.begin(); it != HoughCandidatesOther.end(); ++it) {
    double total = 0;
    // Sort by increasing z
    std::sort((*it).begin(), (*it).end(), TMS_Hit::SortByZInc);
    for (auto hit = (*it).begin(); hit != (*it).end(); ++hit) {
      total += (*hit).GetE();
    }
    TrackEnergyOther.push_back(total);
  }
}

// Calculate the track length for each track
void TMS_TrackFinder::CalculateTrackLengthOne() {
  // Look at the reconstructed tracks
  if (HoughCandidatesOne.size() == 0) return;

  // Loop over each Hough Candidate and find the track length
  for (auto it = HoughCandidatesOne.begin(); it != HoughCandidatesOne.end(); ++it) {
    double total = 0;

    // Sort by increasing z
    std::sort((*it).begin(), (*it).end(), TMS_Hit::SortByZInc);

    for (auto hit = (*it).begin(); hit != (*it).end(); ++hit) {
      auto nexthit = *(hit+1);
      // Use the geometry to calculate the track length between hits
      TVector3 point1((*hit).GetNotZ(), -200, (*hit).GetZ());
      TVector3 point2(nexthit.GetNotZ(), -200, nexthit.GetZ());
      if ((point2-point1).Mag() > 100) continue;
      double tracklength = TMS_Geom::GetInstance().GetTrackLength(point1, point2);
      total += tracklength;
    }
    TrackLengthOne.push_back(total);
  }
}

void TMS_TrackFinder::CalculateTrackLengthOther() {
  // Look at the reconstructed tracks
  if (HoughCandidatesOther.size() == 0) return;
  
  // Loop over each Hough Candidate and finde the track length
  for (auto it = HoughCandidatesOther.begin(); it != HoughCandidatesOther.end(); ++it) {
    double total = 0;

    // Sort by increasing z
    std::sort((*it).begin(), (*it).end(), TMS_Hit::SortByZInc);
    
    for (auto hit = (*it).begin(); hit != (*it).end(); ++hit) {
      auto nexthit = *(hit+1);
      // Ue the geometry to calculate the track length between hits
      TVector3 point1((*hit).GetNotZ(), -200, (*hit).GetZ());
      TVector3 point2(nexthit.GetNotZ(), -200, nexthit.GetZ());
      if ((point2-point1).Mag() > 100) continue;
      double tracklength = TMS_Geom::GetInstance().GetTrackLength(point1, point2);
      total += tracklength;
    }
    TrackLengthOther.push_back(total);
  }
}

//void TMS_TrackFinder::HoughTransform(const std::vector<TMS_Hit> &TMS_Hits) {
std::vector<std::vector<TMS_Hit> > TMS_TrackFinder::HoughTransform(const std::vector<TMS_Hit> &TMS_Hits, const int &hitgroup) {

  // The returned vector of tracks
  std::vector<std::vector<TMS_Hit> > LineCandidates;

  // Check it's not empty
  if (TMS_Hits.empty()) return LineCandidates;

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned = CleanHits(TMS_Hits);
  if (TMS_Hits_Cleaned.empty()) return LineCandidates;

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kYBar);

  // Do a spatial analysis of the hits in y and x around low z to ignore hits that are disconnected from other hits
  // Includes a simple sort in decreasing z (plane number)
  SpatialPrio(TMS_xz);

  int nXZ_Hits_Start = TMS_xz.size();
#ifdef DEBUG
  std::cout << "Starting Hough transform with: " << nXZ_Hits_Start << " hits" << std::endl;
#endif

  // We'll be moving out TMS_xz and TMS_yz and putting them into candidates
  // Keep running successive Hough transforms until we've covered 80% of hits (allow for maximum 4 runs)
  int nRuns = 0;
  //std::cout << "All hits" << std::endl;
  //for (auto &hit: TMS_xz) {
    //std::cout << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << std::endl;
  //}

  while (double(TMS_xz.size()) > nHits_Tol*nXZ_Hits_Start && 
      TMS_xz.size() > nMinHits && 
      nRuns < nMaxHough) {

    // The candidate vectors
    std::vector<TMS_Hit> TMS_xz_cand;
    if (TMS_xz.size() > 0) TMS_xz_cand = RunHough(TMS_xz, hitgroup);
    
    if (TMS_xz_cand.size() == 0) {
      nRuns++;
      if (hitgroup == 1) { // hitgroup 1 means OneHitGroup
        delete HoughLinesOne.back().second;

        HoughLinesOne.pop_back(); // Remove the built Hough line
      } else if (hitgroup == 2) { // hitgroup 2 means OtherHitGroup
        delete HoughLinesOther.back().second;
	
	HoughLinesOther.pop_back();
      }
      break;
    }

    // Move into the candidate vector
    for (auto &i: TMS_xz_cand) Candidates.push_back(std::move(i));

    // Loop over vector and remove used hits
    for (auto jt = TMS_xz_cand.begin(); jt != TMS_xz_cand.end();++jt) {
      for (auto it = TMS_xz.begin(); it!= TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    // Push back the candidates into the total candidates
    LineCandidates.push_back(std::move(Candidates));
    nRuns++;
  }
#ifdef DEBUG
  std::cout << "Ran " << nRuns << " Hough algo" << std::endl;
#endif
  //std::cout << LineCandidates.size() << " hough lines" << std::endl;

  //std::cout << "Hough lines before cleaning: " << std::endl;
  //for (auto &line: LineCandidates) {
    //std::cout << "New line" << std::endl;
    //for (auto &hit : line) {
      //std::cout << hit.GetZ() << ", " << hit.GetNotZ() << "(" << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << ")" << std::endl;
    //}
  //}

  // Can clean up hough hits a little bit
  // We may have tracks that have their last candidate adjacent to the start of another track. If so, we should probably merge tracks

  // Loop over our Hough Candidates and merge 
  // The Hough candidates aren't necessarily ordered after any metric
  // Skip for 0 or single track events
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_MergeTracks() && LineCandidates.size() > 1) {

    // Loop over the Hough candidates
    int lineit = 0;
    for (std::vector<std::vector<TMS_Hit> >::iterator it = LineCandidates.begin(); it != LineCandidates.end(); ++it, ++lineit) {
      if ((*it).size() == 0) continue;

      // Sort each in decreasing z (already done)
      SpatialPrio((*it));

      // Get first and last hits
      TMS_Hit firsthit = (*it).front();
      TMS_Hit lasthit = (*it).back();
      int last_hit_z = lasthit.GetPlaneNumber();
      int last_hit_notz = lasthit.GetBarNumber();

      //std::cout << "Running on track with size: " << (*it).size() << std::endl;

      double HoughOneInter_1 = 0.000;
      double HoughOneSlope_1 = 0.000;
      double HoughOtherInter_1 = 0.000;
      double HoughOtherSlope_1 = 0.000;

      if (hitgroup == 1) {
        HoughOneInter_1 = HoughLinesOne[lineit].second->GetParameter(0);
        HoughOneSlope_1 = HoughLinesOne[lineit].second->GetParameter(1);
      } else if (hitgroup == 2) {
        HoughOtherInter_1 = HoughLinesOther[lineit].second->GetParameter(0);
        HoughOtherSlope_1 = HoughLinesOther[lineit].second->GetParameter(1);
      }

      // Now loop over the remaining hits
      // Need to keep a conventional iterator to remove the HoughLine vector entries
      int lineit_2 = 0;
      for (std::vector<std::vector<TMS_Hit> >::iterator jt = LineCandidates.begin(); jt != LineCandidates.end(); ++jt, ++lineit_2) {

        if ((*jt).size() == 0) continue;
        // Let's get the first and last hits of each track
        // Sort each in decreasing z
        SpatialPrio(*jt);

        // Get the first and last hit for this second track
        TMS_Hit firsthit_2 = (*jt).front();
        TMS_Hit lasthit_2 = (*jt).back();
        // Check that it's not the same hit
        if (firsthit_2 == firsthit && lasthit_2 == lasthit) {
          continue;
        }

        int first_hit_z_2 = firsthit_2.GetPlaneNumber();
        int first_hit_notz_2 = firsthit_2.GetBarNumber();

        //std::cout << last_hit_z << ", " << last_hit_notz << std::endl;
        //std::cout << "against" << std::endl;
        //std::cout << first_hit_z_2 << ", " << first_hit_notz_2 << std::endl;

        // Now check how far away the hits are
        bool mergehits = (abs(first_hit_z_2 - last_hit_z) <= 2 && 
                          abs(first_hit_notz_2 - last_hit_notz) <= 2);

        double HoughOneInter_2 = 0.000;
        double HoughOneSlope_2 = 0.000;
        double HoughOtherInter_2 = 0.000;
        double HoughOtherSlope_2 = 0.000;

        if (hitgroup == 1) {
          HoughOneInter_2 = HoughLinesOne[lineit_2].second->GetParameter(0);
          HoughOneSlope_2 = HoughLinesOne[lineit_2].second->GetParameter(1);
	      } else if (hitgroup == 2) {
  	      HoughOtherInter_2 = HoughLinesOther[lineit_2].second->GetParameter(0);
      	  HoughOtherSlope_2 = HoughLinesOther[lineit_2].second->GetParameter(1);
      	}

        //std::cout << "Hough pars: " << std::endl;
        //std::cout << HoughInter_1 << ", " << HoughInter_2 << std::endl;
        //std::cout << HoughSlope_1 << ", " << HoughSlope_2 << std::endl;
        // Now check how similar the Hough lines are
        bool mergehough;
	if (hitgroup == 1) {
	  mergehough = (fabs(HoughOneInter_2 - HoughOneInter_1) < 100 && 
                           fabs(HoughOneSlope_2 - HoughOneSlope_1) < 0.01);
	} else if (hitgroup == 2) {
	  mergehough = (fabs(HoughOtherInter_2 - HoughOtherInter_1) < 100 &&
			     fabs(HoughOtherSlope_2 - HoughOtherSlope_1) < 0.01);
	}

        // Check if we should merge or not
        if (!mergehits && !mergehough) {
          continue;
        }

        //std::cout << "merged track" << std::endl;

        // Copy over the contents of hits_2 into hits
        for (TMS_Hit movehits: (*jt)) {
          (*it).push_back(movehits);
        }

        // Clear out the original vector
        (*jt).clear();
        //std::cout << "After merge size: " << (*it).size() << std::endl;
        //std::cout << "After merge del vector: " << (*jt).size() << std::endl;
        // After a merge we need to re-sort the original vector and recalculate the first and last hits
        SpatialPrio((*it));
        firsthit = (*it).front();
        lasthit = (*it).back();
        last_hit_z = lasthit.GetPlaneNumber();
        last_hit_notz = lasthit.GetBarNumber();
      }
      //std::cout << "Line " << lineit << " after check" << std::endl;
      //for (std::vector<TMS_Hit>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt) {
        //std::cout << (*jt).GetPlaneNumber() << ", " << (*jt).GetBarNumber() << std::endl;
      //}
    }
  }

  //std::cout << "Hough lines right after track merge: " << std::endl;
  //for (auto &line: LineCandidates) {
    //std::cout << "New line" << std::endl;
    //for (auto &hit : line) {
      //std::cout << hit.GetZ() << ", " << hit.GetNotZ() << "(" << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << ")" << std::endl;
    //}
  //}

  // Loop over and remove the empty tracks
  int linenumber = 0;
  for (std::vector<std::vector<TMS_Hit> >::iterator it = LineCandidates.begin(); it != LineCandidates.end(); ) {
    if ((*it).size() == 0) {
      it = LineCandidates.erase(it);
      if (hitgroup == 1) {
        HoughLinesOne.erase(HoughLinesOne.begin()+linenumber);
      } else if (hitgroup == 2) {
        HoughLinesOther.erase(HoughLinesOther.begin()+linenumber);
      }
    } else {
      ++it;
      ++linenumber;
    }
  }

  //std::cout << "Hough lines after track merge: " << std::endl;
  //for (auto &line: LineCandidates) {
    //std::cout << "New line" << std::endl;
    //for (auto &hit : line) {
      //std::cout << hit.GetZ() << ", " << hit.GetNotZ() << "(" << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << ")" << std::endl;
    //}
  //}

  // Now finally connect the start and end points of each Hough line with Astar to get the most efficient path along the hough hits
  // This helps remove biases in the greedy hough adjacent hit merging
  // For A* after Hough, better to run with Euclidean metric, over-riding the previous
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_RunAStar()) {
    int tracknumber = 0;
    for (std::vector<std::vector<TMS_Hit> >::iterator it = LineCandidates.begin(); it != LineCandidates.end(); ) {
      //int nhoughhits = track.size();
      // Need to sort each line in z
      SpatialPrio(*it);
      std::vector<TMS_Hit> CleanedHough = RunAstar(*it); 
      unsigned int ncleaned = CleanedHough.size();

      // Now replace the old hits with this cleaned version
      // Can remove Hough hits if not big enough
      // Also remember to remove the line
      if (ncleaned < 1) {
        it = LineCandidates.erase(it);
	      if (hitgroup == 1) {
          HoughLinesOne.erase(HoughLinesOne.begin()+tracknumber);
	      } else if (hitgroup == 2) {
	        HoughLinesOther.erase(HoughLinesOther.begin()+tracknumber);
        }
      } else {
        *it = CleanedHough;
        it++;
        tracknumber++;
      }
    }
  }

  //std::cout << "Hough lines after A* slimming: " << std::endl;
  //for (auto &line: LineCandidates) {
    //std::cout << "New line" << std::endl;
    //for (auto &hit : line) {
      //std::cout << hit.GetZ() << ", " << hit.GetNotZ() << "(" << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << ")" << std::endl;
    //}
  //}

  //std::cout << "Hough lines after cleaning: " << std::endl;
  //for (auto &line: LineCandidates) {
    //std::cout << "New line" << std::endl;
    //for (auto &hit : line) {
      //std::cout << hit.GetZ() << ", " << hit.GetNotZ() << "(" << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << ")" << std::endl;
    //}
  //}
  //std::cout << LineCandidates.size() << " hough lines at end" << std::endl;

  return LineCandidates;
};

void TMS_TrackFinder::MaskHits(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask) {

#ifdef DEBUG
  int sizebef = Orig.size();
  std::cout << "Size before masking: " << Orig.size() << std::endl;
  std::cout << "Mask size: " << Mask.size() << std::endl;
#endif

  //std::cout << "Before masking: " << std::endl;
  //std::cout << "Original size: " << Orig.size() << std::endl;
  //std::cout << "Masked size: " << Mask.size() << std::endl;

  // Loop over each hit in the mask
  for (auto &MaskHit: Mask) {
    for (auto it = Orig.begin(); it != Orig.end(); ) {
      if (*it == MaskHit) it = Orig.erase(it);
      else it++;
    }
  }
  //std::cout << "After masking: " << std::endl;
  //std::cout << "Original size: " << Orig.size() << std::endl;
  //std::cout << "Masked size: " << Mask.size() << std::endl;

#ifdef DEBUG
  int sizeaft = Orig.size();
  std::cout << "Removed " << sizebef - sizeaft << " in masking (" << sizebef << "-" << sizeaft << ")" << std::endl;
#endif
}

std::vector<std::vector<TMS_Hit> > TMS_TrackFinder::FindClusters(const std::vector<TMS_Hit> &TMS_Hits) {
  // First clean up double hits in bars
  std::vector<TMS_Hit> MaskedHits = CleanHits(TMS_Hits);
  // The vector of DBSCAN points
  std::vector<TMS_DBScan_Point> DB_Points;
  // For each of the hits make a temporary DBSCAN point
  for (auto &i: MaskedHits) {
    // Use the plane number
    double x = i.GetPlaneNumber();
    // And the height of each bar
    double y = i.GetBarNumber();
#ifdef DEBUG
    double y2 = i.GetNotZ()/i.GetNotZw();
    std::cout << "DBSCAN on " << x << " " << y << " y2: " << y2 << std::endl;
#endif
    // Not using z so set to zero (won't impact distance between points)
    // Maybe update this if we have 3D info somehow
    double z = 0;
    TMS_DBScan_Point point(x, y, z);
    DB_Points.push_back(point);
  }

  // Create the DBSCAN instance
  DBSCAN.SetPoints(DB_Points); //TODO how to get this working with separation? Where is this referring to?
  DBSCAN.RunDBScan();

  //std::vector<TMS_DBScan_Point> NoisePoints = DBSCAN.GetNoise();
  std::vector<std::vector<TMS_DBScan_Point> > ClusterPoints = DBSCAN.GetClusters();

  // Convert the ClusterPoints to Clusters of TMS hits
  std::vector<std::vector<TMS_Hit> > ClusterHits;
  ClusterHits.resize(ClusterPoints.size());
  // For each hit find its cluster
  for (auto &hit: MaskedHits) {
    for (auto &Clusters: ClusterPoints) {
      for (auto &Point: Clusters) {
        // Check the matching hit
        if (Point.x == hit.GetPlaneNumber() && 
            Point.y == hit.GetBarNumber()) {
          ClusterHits[Point.ClusterID-1].push_back(std::move(hit));
        }
      }
    }
  }

#ifdef DEBUG
  DBSCAN.Print();
#endif
  return ClusterHits;
}

// Requires hits to be ordered in z
std::vector<TMS_Hit> TMS_TrackFinder::RunHough(const std::vector<TMS_Hit> &TMS_Hits, const int &hitgroup) {

  // Check if we're in XZ view
  bool IsXZ = ((TMS_Hits[0].GetBar()).GetBarType() == TMS_Bar::kYBar);

  // Recalculate Hough parameters event by event... not fully tested!
  bool VariableHough = false;
  if (VariableHough) {
    // Recalculate the min and max of slope and intercept for every event
    // Hits are ordered in z already
    double maxz = TMS_Const::TMS_Thin_Start;
    double minz = TMS_Const::TMS_Thick_End;
    double min_notz = TMS_Const::TMS_End_Exact[0];
    double max_notz =  TMS_Const::TMS_Start_Exact[0];
    for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it!=TMS_Hits.end(); ++it) {
      double z = (*it).GetZ();
      double not_z = (*it).GetNotZ();
      if (z > maxz) maxz = z;
      if (z < minz) minz = z;
      if (not_z > max_notz) max_notz = not_z;
      if (not_z < min_notz) min_notz = not_z;
    }
    // 1.8 comes from wanting to cover vertices that are maximally 80% between min and max z
    double slope = (max_notz-min_notz)*1.8/(maxz-minz);
    double intercept = -1*slope*(minz+0.8*(maxz-minz));
    //std::cout << "***" << std::endl;
    //std::cout << "z range: " << minz << " " << maxz << std::endl;
    //std::cout << "notz range: " << min_notz << " " << max_notz << std::endl;
    //std::cout << "slope: " << slope << std::endl;
    //std::cout << "intercept: " << intercept << std::endl;

    // now try setting these as the max/min
    SlopeMin = -1*slope;
    SlopeMax = slope;
    InterceptMin = -1*intercept;
    InterceptMax = intercept;
    InterceptWidth = (InterceptMax-InterceptMin)/nIntercept;
    SlopeWidth = (SlopeMax-SlopeMin)/nSlope;
  }

  //std::cout << "All hits in Hough" << std::endl;
  //for (auto &hit: TMS_Hits) {
    //std::cout << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << std::endl;
  //}

  double slope, intercept;
  // Calculate the Hough lines
  GetHoughLine(TMS_Hits, slope, intercept);
  if (hitgroup == 1) {
    HoughLineOne->SetParameter(0, intercept);
    HoughLineOne->SetParameter(1, slope);
  } else if (hitgroup == 2) {
    HoughLineOther->SetParameter(0, intercept);
    HoughLineOther->SetParameter(1, slope);
  }

  // Different fitting regions for XZ and YZ views: 
  // Most of the bending happens in xz, so fit only until the transition region. 
  // For the yz view, fit the entire region
  //if (IsXZ) HoughLine->SetRange(zMinHough, zMaxHough);
  //else HoughLine->SetRange(TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_End);
  
  TF1* HoughCopy = (TF1*)HoughLineOne->Clone();

  if (hitgroup == 1) {
    HoughLineOne->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineOne->Clone();
  } else if (hitgroup == 2) {
    HoughLineOther->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineOther->Clone();
  } else {
#ifdef DEBUG
    std::cout << "Something is going wrong with the assigning of hitgroup" << std::endl;
#endif
    return TMS_Hits;
  }

  std::pair<bool, TF1*> HoughPairs = std::make_pair(IsXZ, HoughCopy);
  if (hitgroup == 1) {
    HoughLinesOne.push_back(std::move(HoughPairs));
  } else if (hitgroup == 2) {
    HoughLinesOther.push_back(std::move(HoughPairs));
  }

  // Then run a clustering on the Hough Transform
  // Hough transform is most likely to pick out straigh features at begining of track candidate, so start looking there

  // Make a hard copy since we don't want to remove the original hits
  // Just make a pool of hits, and tracked hits
  // We'll pull hit out of the HitPool and put them into Candidates
  std::vector<TMS_Hit> HitPool = TMS_Hits;
  std::vector<TMS_Hit> returned;

  // Move hits from the pool into the candidates, and remove the candidates from the pool
  // HitPool is going to shrink or stay the same here
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it != HitPool.end();) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = (*it).GetBar();
    double zhit = hit.GetZ();
    // If z position is above region of interest, ignore hit
    if (IsXZ && zhit > zMaxHough) {
      ++it;
      continue;
    }
    
    double HoughPoint = HoughLineOne->Eval(zhit);
    if (hitgroup == 1) {
      HoughPoint = HoughLineOne->Eval(zhit);
    } else if (hitgroup == 2) {
      HoughPoint = HoughLineOther->Eval(zhit);
    }

    // Hough point is inside bar -> start clustering around bar
    if (( bar.Contains(HoughPoint, zhit) ||
          bar.Contains(HoughPoint-bar.GetNotZw(), zhit) ||
          bar.Contains(HoughPoint+bar.GetNotZw(), zhit) )) {
      // Move into line vector and remove from pool
      returned.push_back(std::move(hit));
      it = HitPool.erase(it);
    } else {
      ++it;
    }
  }

  //std::cout << "New line " << std::endl;
  //for (auto &hits : returned) std::cout << hits.GetPlaneNumber() << " " << hits.GetBarNumber() << std::endl;
//
  //std::cout << "Hit pool before start/end dist" << std::endl;
  //for (auto &hit: HitPool) {
    //std::cout << hit.GetPlaneNumber() << "," << hit.GetBarNumber() << std::endl;
  //}

  if (returned.empty()) return returned;

  /*
  // Clean up tracks by asking that the first and last hit of a track has neighbours
  // Otherwise we might have just included this by the Hough line being drawn far and just happened to hit another active cell
  // Sort each in decreasing z
  int indexstart = -1;
  double dist = 999;
  int vecsize = returned.size();
  while (dist > 3 && indexstart < vecsize) {
    indexstart++;
    int startz = returned[indexstart].GetPlaneNumber();
    int startx = returned[indexstart].GetBarNumber();
    // Check adjacent
    int distz = 0;
    int distx = 0;
    int startz_1 = startz;
    int startx_1 = startx;
    // which index are we checking against
    int index = indexstart;
    std::cout << "Running on " << startz << ", " << startx << std::endl;
    while (abs(distz)+abs(distx) == 0 && index < vecsize) {
      // Look at the next hit
      index++;
      startz_1 = returned[index].GetPlaneNumber();
      startx_1 = returned[index].GetBarNumber();
      distz = startz_1-startz;
      distx = startx_1-startx;
      std::cout << index << " " << startz_1 << "," << startx_1 << std::endl;
    }
    std::cout << "Found good comparison at " << startz_1 << ", " << startx_1 << std::endl;
    dist = sqrt(distx*distx+distz*distz);
    std::cout << dist << std::endl;
  }
  std::cout << "Will remove from " << indexstart << " to " << returned.size() << std::endl;
  for (std::vector<TMS_Hit>::iterator it = returned.begin(); it != returned.begin()+indexstart; it++) {
    HitPool.push_back(std::move(*it));
  }
  returned.erase(returned.begin(), returned.begin()+indexstart);

  if (returned.empty()) return returned;

  SpatialPrio(returned);
  std::cout << "Line before end point" << std::endl;
  for (auto &hit: returned) {
    std::cout << hit.GetPlaneNumber() << "," << hit.GetBarNumber() << std::endl;
  }

  // Do the same for the end point 
  dist = 999;
  vecsize = returned.size();
  indexstart = vecsize;
  std::cout << indexstart << std::endl;
  while (dist > 3 && indexstart > 0) {
    indexstart--;
    int startz = returned[indexstart].GetPlaneNumber();
    int startx = returned[indexstart].GetBarNumber();
    // Check adjacent
    int distz = 0;
    int distx = 0;
    int startz_1 = startz;
    int startx_1 = startx;
    // which index are we checking against
    int index = indexstart;
    std::cout << "Running on " << startz << ", " << startx << std::endl;
    while (abs(distz)+abs(distx) == 0 && index > 0) {
      index--;
      startz_1 = returned[index].GetPlaneNumber();
      startx_1 = returned[index].GetBarNumber();
      std::cout << index << " " << startz_1 << "," << startx_1 << std::endl;
      distz = startz_1-startz;
      distx = startx_1-startx;
    }
    std::cout << "Found good comparison at " << startz_1 << ", " << startx_1 << std::endl;
    dist = sqrt(distx*distx+distz*distz);
    std::cout << dist << std::endl;
  }
  indexstart++;
  std::cout << "Will remove from " << indexstart << " to " << returned.size() << std::endl;
  for (std::vector<TMS_Hit>::iterator it = returned.begin()+indexstart; it != returned.end(); it++) {
    HitPool.push_back(std::move(*it));
  }
  returned.erase(returned.begin()+indexstart, returned.end());

  if (returned.empty()) return returned;
  */

  //std::cout << "Hit pool before walking" << std::endl;
  //for (auto &hit: HitPool) {
    //std::cout << hit.GetPlaneNumber() << "," << hit.GetBarNumber() << std::endl;
  //}

  //std::cout << "Line before walking" << std::endl;
  //for (auto &hit: returned) {
    //std::cout << hit.GetPlaneNumber() << "," << hit.GetBarNumber() << std::endl;
  //}

  // Now walk along the Hough hits and add on adjacent hits
  // Hough is most likely to find hits upstream due to bending being less there
  //std::cout << returned.size() << " line hits before walking downstream" << std::endl;
  WalkDownStream(returned, HitPool);
  //std::cout << returned.size() << " line hits after walking downstream" << std::endl;
  WalkUpStream(returned, HitPool);
  //std::cout << returned.size() << " line hits after walking upstream" << std::endl;

  // The vector of DBSCAN points
  std::vector<TMS_DBScan_Point> DB_Points;
  // For each of the hits make a temporary DBSCAN point
  for (auto &i: returned) {
    // Use the plane number
    double x = i.GetPlaneNumber();
    // And the height of each bar
    double y = i.GetBarNumber();

    // Not using z so set to zero (won't impact distance between points)
    // Maybe update this if we have 3D info somehow
    double z = 0;
    TMS_DBScan_Point point(x, y, z);
    DB_Points.push_back(point);
  }
  // See if DBSCAN can find different clusters
  TMS_DBScan LineDB(3, 5, DB_Points);
  LineDB.RunDBScan();
  std::vector<std::vector<TMS_DBScan_Point> > ClusterPoints = LineDB.GetClusters();
  if (ClusterPoints.empty()) {
    std::vector<TMS_Hit> emptyvec;
    return emptyvec;
  }

  // Convert the ClusterPoints to Clusters of TMS hits
  std::vector<std::vector<TMS_Hit> > ClusterHits;
  ClusterHits.resize(ClusterPoints.size());
  // For each hit find its cluster
  for (auto &hit: returned) {
    for (auto &Clusters: ClusterPoints) {
      for (auto &Point: Clusters) {
        // Check the matching hit
        if (Point.x == hit.GetPlaneNumber() && 
            Point.y == hit.GetBarNumber()) {
          ClusterHits[Point.ClusterID-1].push_back(std::move(hit));
        }
      }
    }
  }

  //std::cout << "Found " << ClusterHits.size() << " DBScan clusters" << std::endl;
  if (ClusterHits.size() > 0) {
    unsigned int largest = 0;
    for (unsigned int i = 0; i < ClusterHits.size(); ++i) if (ClusterHits[i].size() > ClusterHits[largest].size()) largest = i;
    //std::cout << "Cluster " << largest << " is largest" << std::endl;
    //for (unsigned int i = 0; i < ClusterHits.size(); ++i) {
      //std::cout << "Cluster " << i << std::endl;
      //for (auto &hit: ClusterHits[i]) std::cout << hit.GetPlaneNumber() << ", " << hit.GetBarNumber() << std::endl;
    //}
    // Now remove the hits that aren't in the cluster hits
    for (auto it = returned.begin(); it != returned.end(); ) {
      bool match = false;
      for (auto jt = ClusterHits[largest].begin(); jt != ClusterHits[largest].end(); ++jt) {
        if ((*it) == (*jt)) match = true;
      }
      if (!match) {
        // Move back into the hitpool
        HitPool.push_back(std::move(*it));
        it = returned.erase(it);
      }
      else it++;
    }
  }

  // Finally run A* to find the shortest path from start to end
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_RunAStar()) {
    SpatialPrio(returned);
    std::vector<TMS_Hit> vec = RunAstar(returned);
    //std::cout << HitPool.size() << " before set_difference" << std::endl;
    //std::cout << "return before a*" << returned.size() << std::endl;

    // Only overwrite when necessary
    if (vec.size() != returned.size()) {
      // Put the missing hits back into the hit pool
      for (auto it = vec.begin(); it != vec.end(); ++it) {
        for (auto jt = returned.begin(); jt != returned.end(); ) {
          if ((*it) == (*jt)) {
            jt = returned.erase(jt);
          } else {
            ++jt;
          }
        }
      }
      for (auto &hit: returned) HitPool.emplace_back(std::move(hit));
      //std::cout << "return after a*" << std::endl;
      //std::cout << returned.size() << std::endl;
      //std::cout << "a* vec" << std::endl;
      //std::cout << vec.size() << std::endl;
      returned = vec;
    }
    //std::cout << HitPool.size() << " after set_difference" << std::endl;
  }

  // Now finally check if any hough tracks are too short, or have too few hits
  // Makes sense to do this right at the end when all the merging and cleaning has been run
  // Order them in z
  SpatialPrio(returned);
  // Now check that the Hough candidates are long enough
  /*
  double xend = (returned).back().GetNotZ();
  double zend = (returned).back().GetZ();
  double xstart = (returned).front().GetNotZ();
  double zstart = (returned).front().GetZ();
  double xdist = xstart-xend;
  double zdist = zstart-zend;
  double dist = sqrt(xdist*xdist+zdist*zdist);
  */
  double xend = (returned).back().GetBarNumber();
  double zend = (returned).back().GetPlaneNumber();
  double xstart = (returned).front().GetBarNumber();
  double zstart = (returned).front().GetPlaneNumber();
  double xdist = xstart-xend;
  double zdist = zstart-zend;
  double dist = sqrt(xdist*xdist+zdist*zdist);
  unsigned int nhits = (returned).size();
  // Calculate the minimum distance in planes and bars instead of physical distance
  if (dist < MinDistHough || nhits < nMinHits) {
    //std::cout << "Failed distance or min hits cut" << std::endl;
    std::vector<TMS_Hit> emptyvec;
    return emptyvec;
  }

  return returned;
}

// Find the bin for the accumulator
int TMS_TrackFinder::FindBin(double c) {
  // Since we're using uniform binning no need for binary search or similar
  if (c > InterceptMax) c = InterceptMax;
  int bin = (c-InterceptMin)/InterceptWidth;
  return bin;
}


// Convert Accumulator to a TH2D
TH2D *TMS_TrackFinder::AccumulatorToTH2D(bool zy) {

  std::string Name;
  if (zy) {
    Name = "TMS_TrackFinder_Accumulator_zy";
  } else {
    Name = "TMS_TrackFinder_Accumulator_zx";
  }
  TH2D *accumulator = new TH2D(Name.c_str(), (Name+";m (slope);c (intercept) (mm)").c_str(), nSlope, SlopeMin, SlopeMax, nIntercept, InterceptMin, InterceptMax);

  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      accumulator->SetBinContent(i+1, j+1, Accumulator[i][j]);
    }
  }

  int maxx, maxy, maxz;
  accumulator->GetMaximumBin(maxx, maxy, maxz);
  double maxtheta = accumulator->GetXaxis()->GetBinCenter(maxx);
  double maxrho = accumulator->GetYaxis()->GetBinCenter(maxy);
  accumulator->SetTitle(Form("#splitline{%s}{m_{max}=%.2f c_{max}=%.2f}", accumulator->GetTitle(), maxtheta, maxrho));

  // Set the minimum (easier to draw)
  double maxcounts = accumulator->GetMaximum();
  accumulator->SetMinimum(maxcounts*0.8);

  return accumulator;
}

// Implement A* algorithm for track finding, starting with most upstream to most downstream hit
void TMS_TrackFinder::BestFirstSearch(const std::vector<TMS_Hit> &TMS_Hits, const int &hitgroup) {

  // Set the Heuristic cost calculator

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned = CleanHits(TMS_Hits);

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kYBar);
  //std::vector<TMS_Hit> TMS_yz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);

  // Do a spatial analysis of the hits in y and x around low z to ignore hits that are disconnected from other hits
  // Includes a simple sort in decreasing z (plane number)
  SpatialPrio(TMS_xz);
  //SpatialPrio(TMS_yz);

  int nRuns = 0;
  int nXZ_Hits_Start = TMS_xz.size();

  while (double(TMS_xz.size()) > nHits_Tol*nXZ_Hits_Start && 
      TMS_xz.size() > nMinHits && 
      nRuns < nMaxHough) {
    //std::cout << "yz = " << TMS_yz.size() << " before a = " << a << std::endl;
    //std::cout << "xz = " << TMS_xz.size() << " before a = " << a << std::endl;
    //std::vector<TMS_Hit> AStarHits_yz;
    std::vector<TMS_Hit> AStarHits_xz;

    // Run on x-z first since that's where the gaps occur, leading to broken tracks
    // We can save where this happens in z and make the yz reconstruction aware of the gap
    if (TMS_xz.size() > 0) AStarHits_xz = RunAstar(TMS_xz);
    else return;
    //if (TMS_yz.size() > 0) AStarHits_yz = RunAstar(TMS_yz);

    // Then order them in z
    SpatialPrio(AStarHits_xz);
    //SpatialPrio(AStarHits_yz);

    // Copy over to the candidates
    //for (auto &i : AStarHits_xz) Candidates.push_back(std::move(i));
    //for (auto &i : AStarHits_yz) Candidates.push_back(std::move(i));

    // Loop over vector and remove used hits
    /*
       for (auto jt = AStarHits_yz.begin(); jt!=AStarHits_yz.end();++jt) {
       for (auto it = TMS_yz.begin(); it!= TMS_yz.end();) {
       if ((*it) == (*jt)) it = TMS_yz.erase(it);
       else it++;
       }
       }
       */

    // Loop over vector and remove used hits
    for (auto jt = AStarHits_xz.begin(); jt!=AStarHits_xz.end();++jt) {
      for (auto it = TMS_xz.begin(); it!=TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    // Only push back if we have more than one candidate
    if (AStarHits_xz.size() > nMinHits) {
      if (hitgroup == 1) {    
        HoughCandidatesOne.push_back(std::move(AStarHits_xz));
      } else if (hitgroup == 2) {
        HoughCandidatesOther.push_back(std::move(AStarHits_xz));
      }
    }
    nRuns++;
  }
#ifdef DEBUG
  std::cout << "Ran " << nRuns << " Astar algo" << std::endl;
#endif
}

// Remove duplicate hits
// and hits with low energy threshold
std::vector<TMS_Hit> TMS_TrackFinder::CleanHits(const std::vector<TMS_Hit> &TMS_Hits) {

  std::vector<TMS_Hit> TMS_Hits_Cleaned;
  if (TMS_Hits.empty()) return TMS_Hits_Cleaned;
  TMS_Hits_Cleaned = TMS_Hits;
  //  Clean hits has functional overlap with ped sup and merging hits steps, which are now done earlier
  //return TMS_Hits_Cleaned; // TODO (Jeffrey) is this safe?

  // First sort in z so overlapping hits are next to each other in the arrays
  SpatialPrio(TMS_Hits_Cleaned);

  // Loop over the original hits
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits_Cleaned.begin(); 
      it != TMS_Hits_Cleaned.end(); ) {

    // Look ahead to find duplicates
    int nDuplicates = 0;
    for (std::vector<TMS_Hit>::iterator jt = it+1; jt != TMS_Hits_Cleaned.end(); ++jt) {

      // Maybe this hit has already been counted
      double z = (*it).GetZ();
      double y = (*it).GetNotZ();
      //double e = hit.GetE();
      //double t = hit.GetT();

      TMS_Hit hit2 = *(jt);
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      double e2 = hit2.GetE();
      double t2 = hit2.GetT();

      // Merge
      //if (z == z2 && y == y2 && fabs(t2-(*it).GetT()) < TMS_Const::TMS_TimeThreshold)
      if (z == z2 && y == y2) {
        (*it).SetE((*it).GetE()+e2);
        (*it).SetT(((*it).GetT()+t2)/2);
        nDuplicates++;
      }
    }
    // Now remove the duplicates
    if (nDuplicates > 0) {
      it = TMS_Hits_Cleaned.erase(it, it+nDuplicates);
    } else it++;
  }

  // Remove zero entries
  // Strip out hits that are outside the actual volume 
  // This is probably some bug in the geometry that sometimes gives hits in the z=30k mm (i.e. 10m downstream of the end of the TMS)
  // Figure out why these happen?
  for (std::vector<TMS_Hit>::iterator jt = TMS_Hits_Cleaned.begin(); 
      jt != TMS_Hits_Cleaned.end(); ) {

    if ( (*jt).GetZ() > TMS_Const::TMS_End[2] ||  // Sometimes a hit downstream of the end geometry
        (*jt).GetZ() < TMS_Const::TMS_Start[2] ||  // Or upstream of the start...
        (*jt).GetE() < TMS_Const::TMS_EnThres) { // Check energy threshold
      jt = TMS_Hits_Cleaned.erase(jt);
    } else {
      jt++;
    }
  }
  return TMS_Hits_Cleaned;
}

std::vector<TMS_Hit> TMS_TrackFinder::ProjectHits(const std::vector<TMS_Hit> &TMS_Hits, TMS_Bar::BarType bartype) {
  std::vector<TMS_Hit> returned;
  if (TMS_Hits.empty()) return returned;

  for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Hit hit = (*it);
    if (hit.GetBar().GetBarType() == bartype) {
      returned.push_back(std::move(hit));
    }
  }
  return returned;
}

// Needs hits ordered in z
std::vector<TMS_Hit> TMS_TrackFinder::RunAstar(const std::vector<TMS_Hit> &TMS_xz, bool ConnectAll) {

  // Remember which orientation these hits are
  // needed when we potentially skip the air gap in xz (but not in yz!)
  bool IsXZ = ((TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kYBar);
  // Reset remembering where gaps are in xz
  if (IsXZ) PlanesNearGap.clear();

  // Set the first and last hit to calculate the heuristic to
  //aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetNotZ(), TMS_xz.back().GetNotZw());
  aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetBarNumber());

  // Also convert the TMS_Hit to our path-node
  std::vector<aNode> Nodes;

  // Give each node an ID
  int NodeID = 0;
  for (std::vector<TMS_Hit>::const_iterator it = TMS_xz.begin(); it != TMS_xz.end(); ++it, NodeID++) {

    // Use the x position as plane number
    double x = (*it).GetPlaneNumber();
    double y = (*it).GetBarNumber();

    // Make the node
    aNode TempNode(x, y, NodeID);
    // Calculate the Heuristic cost for each of the nodes to the last point
    // This could probably be updated less often...
    TempNode.SetHeuristic(kHeuristic);
    TempNode.SetHeuristicCost(Last);

    Nodes.push_back(std::move(TempNode));
  }

  // Now that we have the Heuristic cost calculated, we can set the first node as the node in the most upstream layer with the smallest Heuristic cost
  int lowest_index = 0;
  double lowest_heur = 999;
  int firstlayer = 999;
  for (size_t i = 0; i < Nodes.size(); ++i) {
    // Check nodes is in first layers
    int layer = Nodes[i].x;
    if (layer > firstlayer) continue;
    // Save the first layer
    firstlayer = layer;
    // Check if the heuristic cost is smaller
    if (Nodes[i].HeuristicCost < lowest_heur && Nodes[i].Neighbours.size() > 0) {
      lowest_heur = Nodes[i].HeuristicCost;
      lowest_index = i;
    }
  }
  //std::cout << "Index of lowest heuristic hit: " << lowest_index << " in layer " << firstlayer << " with heuristic = " << lowest_heur << std::endl;
  std::swap(Nodes[0], Nodes[lowest_index]);

  // Can probably evaluate the last node here: if heuristic is large for the 5 nearest hits it's likely wrong?

  // Now find the neighbours of each node
  for (std::vector<aNode>::iterator it = Nodes.begin(); it != Nodes.end(); ++it) {
    for (std::vector<aNode>::iterator jt = Nodes.begin(); jt != Nodes.end(); ++jt) {

      // First check this node isn't itself!
      if ((*jt) == (*it)) continue;

      // Get the address
      aNode* Pointer = &(*jt);
      // Only connect nearby nodes
      // to make calculation much faster
      // This important can't be 1 though, because that would not connect split hits
      if (!ConnectAll) {
        if (abs((*jt).x - (*it).x) > 3 || 
            abs((*jt).y - (*jt).y) > 3) continue;
      }
      double GroundCost = (*it).CalculateGroundCost(*jt);

      // If a greedy algorithm, have no ground cost (only condition on getting closer to the goal)
      if (IsGreedy) GroundCost = 0;

      // Remember how much it is for this node to connect to this neighbour
      (*it).Neighbours[Pointer] = GroundCost;
    }
  }


  // Now finally loop over all hits and output them
  /*
     for (std::vector<aNode>::iterator it = Nodes.begin(); it != Nodes.end(); ++it) {
     std::cout << "****" << std::endl;
     (*it).Print();
     std::cout << "Neighbours: " << std::endl;
     for (auto i: (*it).Neighbours) {
     std::cout << "  ";
     i.first->Print();
     std::cout << "  move cost: " << i.second << std::endl;
     }
     }
     */

  // Remove hits that only have one neighbour?

  // Look at the first hit and see how many neighbours its neighbour has
  unsigned int nrem = 0;
  unsigned int total = Nodes.size();
  while (total > nrem && Nodes[nrem].Neighbours.size() < 1) {
    nrem++;
  }
  //std::cout << "Removed " << nrem << " nodes from front without neighbours" << std::endl;
  //std::cout << "total: " << total << std::endl;
  if (nrem == total) {
    std::vector<TMS_Hit> emptyvec;
    return emptyvec;
  }

  /*
  // If there's only one neighbour and its only neigbour is the start
  if (Nodes[nrem].Neighbours.size() == 1 &&
  Nodes[nrem+1].Neighbours.find(&Nodes[nrem]) != Nodes[nrem+1].Neighbours.end()) {
  nrem++;
  }
  */

  /*
  // Do the same for the back: make sure end point has neighbours
  unsigned int nrem_front = 0;
  while (total-nrem > nrem_front && Nodes[total-nrem_front].Neighbours.size() < 1) {
  nrem_front++;
  }
  std::cout << "Remove " << nrem_front << " nodes from front without neighbours" << std::endl;
  */

  // Return if there are no nodes left
  /*
     if (nrem >= total || nrem_front >= total) {
     std::vector<TMS_Hit> empt;
     return empt;
     }
     */

  // Recalculate the heuristic
  //if (nrem_front > 1) nrem_front--;
  /*
     for (auto &i: Nodes) {
     i.SetHeuristicCost(Nodes.at(total-nrem_front-1));
     }
     std::cout << "Re-calculating Heuristic cost relative entry " << total-nrem_front-1 << " of " << total << std::endl;
     Nodes.at(total-nrem_front-1).Print();
     */

  // Keep a map of the cost so far for reaching an aNode
  std::unordered_map<int, double> cost_so_far;
  // Keep a map of the where a given aNode was reached from
  std::unordered_map<int, int> came_from;
  // The priority queue of which node to next evaluate
  std::priority_queue<aNode, std::vector<aNode> > pq;

  //if (nrem > 1) nrem--;
  pq.push(Nodes.at(nrem));
  cost_so_far[nrem] = 0;
  came_from[nrem] = 0;
  // Keep track when we hit the last point
  bool LastPoint = false;
  while (!pq.empty()) {
    aNode current = pq.top();
    if (LastPoint) break;
    pq.pop();

    // If this is the last point, need to calculate the final cost
    if (current == Nodes.back()) LastPoint = true;

    // Loop over this node's neighbours
    for (auto &neighbour: current.Neighbours) {
      // Check the cost for getting to this node
      double new_cost = cost_so_far[current.NodeID] + // The current cost for getting to this node (irrespective of neighbours)
        neighbour.first->HeuristicCost + // Add up the Heuristic cost of moving to this neighbour
        neighbour.second; // Add up the ground cost of moving to this neighbour

      // If we've never reached this node, or if this path to the node has less cost
      if (cost_so_far.find(neighbour.first->NodeID) == cost_so_far.end() ||
          new_cost < cost_so_far[neighbour.first->NodeID]) {
        cost_so_far[neighbour.first->NodeID] = new_cost; // Update the cost to get to this node via this path
        pq.push(Nodes.at(neighbour.first->NodeID)); // Add to the priority list
        came_from[neighbour.first->NodeID] = current.NodeID; // Update where this node came from
      }
    }
  }

  NodeID = pq.top().NodeID;
  // Check how big the candidate is of total
  //std::cout << "NodeID of top: " << NodeID << std::endl;

  //std::cout << "Path: " << std::endl;
  //if (NodeID != Last.NodeID) {
  //std::cout << "Did not find end of path" << std::endl;
  //}


  // Now push back the candidates by tracing the path
  std::vector<TMS_Hit> returned;
  //std::cout << "Printing priority queue: " << std::endl;
  while (NodeID != 0) {
    returned.push_back(TMS_xz[NodeID]);
    //std::cout << "Node: " << std::endl;
    //Nodes[NodeID].Print();
    //std::cout << "Came from: " << std::endl;
    NodeID = came_from[NodeID];
    //Nodes[NodeID].Print();
    // If the current node came from itself, we've reached the end
    if (NodeID == came_from[NodeID]) {
      returned.push_back(TMS_xz[NodeID]);
      break;
    }
  }

  // Could walk upstream and downstream from here on

  // Put in some check to see if we connected the nodes from 0 to NodeID_{highest}?

  //if (NodeID != 0) {
  //std::cout << "did not find last point" << std::endl;
  //}

  //std::cout << "found path from " << pq.top().NodeID << std::endl;
  //aNode pqtop = pq.top();
  //pqtop.Print();
  //std::cout << "to " << NodeID << std::endl;
  //Nodes[NodeID].Print();

  //std::cout << " ********* " << std::endl;
  return returned;
}

// Gaps are in xz view at x = -1717 to -1804
// barpos is position of bar, width is width of bar
bool TMS_TrackFinder::NextToGap(double barpos, double width) {
  // Get the center of the adjacent bar
  // If this is inside one of the gaps, this bar is next to a gap
  double pos = barpos+width;
  double neg = barpos-width;

  // Check the top
  if ((pos > TMS_Const::TMS_Dead_Top[0] && pos < TMS_Const::TMS_Dead_Top[1]) ||
      (neg < TMS_Const::TMS_Dead_Top[1] && neg > TMS_Const::TMS_Dead_Top[0])) return true;

  // Check the center
  else if ((pos > TMS_Const::TMS_Dead_Center[0] && pos < TMS_Const::TMS_Dead_Center[1]) ||
      (neg < TMS_Const::TMS_Dead_Center[1] && neg > TMS_Const::TMS_Dead_Center[0])) return true;

  // Check the bottom
  else if ((pos > TMS_Const::TMS_Dead_Bottom[0] && pos < TMS_Const::TMS_Dead_Bottom[1]) ||
      (neg < TMS_Const::TMS_Dead_Bottom[1] && neg > TMS_Const::TMS_Dead_Bottom[0])) return true;

  else return false;

  return false;
}

// Do a quick analysis of the hits that are sorted increasing in z
void TMS_TrackFinder::SpatialPrio(std::vector<TMS_Hit> &TMS_Hits) {
  if (TMS_Hits.empty()) return;

  // First do a normal sort decreasing in z
  std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByZInc);

  /*
  // Get the first 10 hits
  std::vector<double> Hits_First;
  double sum = 0;
  double sum2 = 0;
  int firstplane = TMS_Hits.back().GetPlaneNumber();
  // Read in from the back
  for (std::vector<TMS_Hit>::reverse_iterator it = TMS_Hits.rbegin(); it != TMS_Hits.rend(); ++it) {
  double planeval = (*it).GetPlaneNumber();
  // Look at the first 3 planes in each view
  if (planeval > firstplane+6) break;
  double xval = (*it).GetNotZ();
  Hits_First.push_back(xval);
  sum += xval;
  sum2 += xval*xval;
  }

  // Make the median of the entries
  int nEntries = Hits_First.size();
  if (nEntries == 0) return;
  double avg = sum/nEntries;
  double rms = sqrt(sum2/nEntries - (sum/nEntries)*(sum/nEntries));
  //std::sort(Hits_First.begin(), Hits_First.end());
  //double median = 0;
  //if (nEntries % 2 != 0) median = Hits_First[nEntries/2];
  //else median = 0.5*(Hits_First[(nEntries-1)/2] + Hits_First[(nEntries+1)/2]);
  std::cout << "Looked at " << nEntries << " hits" << std::endl;
  std::cout << "avg: " << avg << std::endl;
  //std::cout << "med: " << median << std::endl;
  std::cout << "rms: " << rms << std::endl;

  // Don't remove hits; only reprioritise them if they're right at the front/back
  // Since we run on back to front we need to do both
  int nTotal = TMS_Hits.size();
  // Figure out how deep in we have to go
  int good = 999;
  for (int i = 0; i < nTotal; ++i) {
  double entry = TMS_Hits[i].GetNotZ();
  // Only allow to switch points in first plane
  if (TMS_Hits[i].GetPlaneNumber() > firstplane) break;
  if (entry < avg+rms && entry > avg-rms) {
  good = i;
  std::cout << "Swapping " << TMS_Hits[i].GetNotZ() << " for " << TMS_Hits[0].GetNotZ() << std::endl;
  TMS_Hits[good] = TMS_Hits[nTotal];
  TMS_Hits[nTotal] = TMS_Hits[i];
  break;
  }
  }

  int lastplane = TMS_Hits.front().GetPlaneNumber();
  std::vector<double> Hits_Last;
  sum = 0;
  sum2 = 0;
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
  double planeval = (*it).GetPlaneNumber();
  // Look at the first 3 planes in each view
  if (planeval > lastplane+6) break;
  double xval = (*it).GetNotZ();
  Hits_Last.push_back(xval);
  sum += xval;
  sum2 += xval*xval;
  }
  nEntries = Hits_Last.size();
  if (nEntries == 0) return;
  avg = sum/nEntries;
  rms = sqrt(sum2/nEntries - (sum/nEntries)*(sum/nEntries));
  std::cout << "Looked at " << nEntries << " hits" << std::endl;
  std::cout << "avg: " << avg << std::endl;
  //std::cout << "med: " << median << std::endl;
  std::cout << "rms: " << rms << std::endl;

  // Figure out how deep in we have to go
  good = 999;
  for (int i = 0; i < nTotal; ++i) {
    double entry = TMS_Hits[i].GetNotZ();
    // Only allow to switch points in first plane
    if (TMS_Hits[i].GetPlaneNumber() < lastplane-4) break;
    if (entry < avg+rms && entry > avg-rms) {
      good = i;
      std::cout << "Swapping " << TMS_Hits[i].GetNotZ() << " for " << TMS_Hits[0].GetNotZ() << std::endl;
      TMS_Hits[good] = TMS_Hits[0];
      TMS_Hits[0] = TMS_Hits[i];
      break;
    }
  }
  */

}


// Walk along the "vec" hits and add in appropriate hits from the "Full"
void TMS_TrackFinder::WalkDownStream(std::vector<TMS_Hit> &vec, std::vector<TMS_Hit> &full) {

  // First mask out the vec hits from the full vector so we have no entries from "vec" in "full"
  MaskHits(full, vec);
  //std::cout << "Line: " << std::endl;
  //for (auto &hits: vec) std::cout << hits.GetZ() << " " << hits.GetNotZ() << "(" << hits.GetPlaneNumber() << ", " << hits.GetBarNumber() << ")" << std::endl;
  //std::cout << "Masked: " << std::endl;
  //for (auto &hits: full) std::cout << hits.GetZ() << " " << hits.GetNotZ() << "(" << hits.GetPlaneNumber() << ", " << hits.GetBarNumber() << ")" << std::endl;

  // Sort in z
  SpatialPrio(vec);
  SpatialPrio(full);

  std::vector<TMS_Hit>::size_type size = vec.size();
  // Now walk along vec (start at second element due to derivative calculation)
  for (std::vector<TMS_Hit>::size_type i = 1; i < size; ++i) {
    // Use to calculate gradient to next point
    //double x = vec[i].GetZ();
    //double y = vec[i].GetNotZ();
    double x = vec[i].GetPlaneNumber();
    double y = vec[i].GetBarNumber();
    int PlaneNumber = vec[i].GetPlaneNumber();
    int BarNumber = vec[i].GetBarNumber();

    // Require continous hits for derivative calculation
    // If we can't find the previous hit at first go (might be two hits in the same plane, so adjacent hits would have an inifite slop), try moving down in the order
    int NeighbourIndex = i-1;
    // Previous plane number
    int PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
    //std::cout << "Comparing hit " << x << "," << y << std::endl;

    bool MoveOn = false;

    // Only look at upstream adjacent planes
    while (PlaneNumber - PlaneNumber_prev == 0) {
      NeighbourIndex--;
      // Check that we aren't going beyond unphysical planes
      if (NeighbourIndex < 0) {
        MoveOn = true;
        break;
      }
      // This loop could probably be improved to just be included in the while?
      // Or do we even need a loop, doesn't look like it...
      PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
      if (PlaneNumber - PlaneNumber_prev > 2) {
        MoveOn = true;
        break;
      }
    }

    // Couldn't find the right hit
    if (MoveOn) continue;

    // Calculate expected gradient from upstream hit
    //double xprev = vec[NeighbourIndex].GetZ();
    //double yprev = vec[NeighbourIndex].GetNotZ();
    double xprev = vec[NeighbourIndex].GetPlaneNumber();
    double yprev = vec[NeighbourIndex].GetBarNumber();

    double grad_exp = (y-yprev)/(x-xprev);

    //std::cout << "Comparing hit " << x << "," << y << " to " << xprev << "," << yprev << " expecing gradient " << grad_exp << std::endl;
    //std::cout << "(" << PlaneNumber << ", " << vec[i].GetBarNumber() << "-" << PlaneNumber_prev << ", " << vec[NeighbourIndex].GetBarNumber() << ")" << std::endl;

    // Save the indices of which particle was best
    // Only allow for one merge

    // Match based on the gradient between previous points
    // Here the "full" vector has been masked to not contain the Hough hits
    for (std::vector<TMS_Hit>::iterator it = full.begin(); it != full.end(); ) {
      TMS_Hit hits = *it;
      int PlaneNumber_cand = hits.GetPlaneNumber();
      int BarNumber_cand = hits.GetBarNumber();

      // If candidate is upstream, continue
      if (PlaneNumber_cand < PlaneNumber) {
        ++it;
        continue;
      }

      if (abs(BarNumber_cand - BarNumber) > 2) {
        ++it;
        continue;
      }

      // Since the candidates plane numbers are ordered in z, once we encounter a higher z, all the remaining ones also won't be mergeable, so save some time by not scanning them
      if (PlaneNumber_cand - PlaneNumber > 2) break;

      // At this point we should have a Hough hit and an adjacent non-Hough hit

      // Allow for broken track in z?
      //double xcand = hits.GetZ();
      //double ycand = hits.GetNotZ();
      double xcand = hits.GetPlaneNumber();
      double ycand = hits.GetBarNumber();

      // Calculate the new gradient between the points
      double grad_new = (ycand-y)/(xcand-x);

      //std::cout << "testing merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;
      //std::cout << "With hit info " << hits.GetZ() << ", " << hits.GetNotZ() << std::endl;
      //std::cout << "and new gradient: " << grad_new << std::endl;

      // If gradient is within 0.6 and the correct sign, accept
      if (fabs(grad_new-grad_exp) <= 2.0) {
        // If the gradient is zero we shouldn't do a sign check
        // But if it's not, check the sign of the gradient doesn't flip
        if (fabs(grad_new) > TMS_Const::TMS_Small_Num && 
            fabs(grad_exp) > TMS_Const::TMS_Small_Num &&
            TMS_Utils::sgn(grad_new) != TMS_Utils::sgn(grad_exp)) {
          ++it;
          continue;
        }

        //std::cout << "Will be merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;

        // Now need to find where we insert this hit since we need to guarantee sorted in z
        // Need to find how many duplicates in z there are; don't necessarily want to insert right after this hit if it is followed by a hit in the same z
        // If this is the largest z, just pop to the back
        //if (xcand >= vec.back().GetZ()) {
        if (xcand >= vec.back().GetPlaneNumber()) {
          vec.push_back(std::move(hits));
          // Need to decrement iterator over line
          // Sometimes we've missed one hit between Hough hits
        } else {
          for (std::vector<TMS_Hit>::iterator jt = vec.begin()+i; jt != vec.end(); ++jt) {
            //double compx = (*jt).GetZ();
            // Look one ahead
            //double compx1 = (*(jt-1)).GetZ();
            double compx = (*jt).GetPlaneNumber();
            double compx1 = (*(jt-1)).GetPlaneNumber();
            if (compx == x) continue;
            if (xcand >= compx) {
              vec.insert(jt+1, std::move(hits));
              break;
            } else if (xcand < compx && xcand >= compx1) {
              // Put in the previous
              vec.insert(jt, std::move(hits));
              break;
            } else {
              std::cout << "Didn't find anywhere to insert " << xcand << std::endl;
            }
          }
        }
        ++size;
        // Remove from the pool
        it = full.erase(it);
      } else {
        it++;
      }
      }
    }
  }

  void TMS_TrackFinder::WalkUpStream(std::vector<TMS_Hit> &vec, std::vector<TMS_Hit> &full) {
    // First mask out the vec hits from the full vector so we have no entries from "vec" in "full"
    MaskHits(full, vec);

    // Sort in *decreasing* z (starting at highest z, i.e. most downstream)
    std::sort(vec.begin(), vec.end(), TMS_Hit::SortByZ);
    std::sort(full.begin(), full.end(), TMS_Hit::SortByZ);

    //std::cout << "Line: " << std::endl;
    //for (auto &hits: vec) std::cout << hits.GetZ() << " " << hits.GetNotZ() << "(" << hits.GetPlaneNumber() << ", " << hits.GetBarNumber() << ")" << std::endl;
    //std::cout << "Masked: " << std::endl;
    //for (auto &hits: full) std::cout << hits.GetZ() << " " << hits.GetNotZ() << "(" << hits.GetPlaneNumber() << ", " << hits.GetBarNumber() << ")" << std::endl;

    // Now walk along vec (start at second element due to derivative calculation
    std::vector<TMS_Hit>::size_type size = vec.size();
    for (std::vector<TMS_Hit>::size_type i = 1; i < size; ++i) {
      // Use to calculate gradient to next point
      //double x = vec[i].GetZ();
      //double y = vec[i].GetNotZ();
      double x = vec[i].GetPlaneNumber();
      double y = vec[i].GetBarNumber();
      int PlaneNumber = vec[i].GetPlaneNumber();
      int BarNumber = vec[i].GetBarNumber();

      // Require continous hits for derivative calculation
      // If we can't find the previous hit at first go (might be two hits in the same plane, so adjacent hits would have an inifite slop), try moving down in the order
      int NeighbourIndex = i-1;
      // Previous plane number
      int PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
      //std::cout << "Comparing hit " << x << "," << y << std::endl;

      bool MoveOn = false;
      // Always look for a hit with different plane number, i.e. we're moving in z
      while (PlaneNumber_prev - PlaneNumber == 0) {
        NeighbourIndex--;
        if (NeighbourIndex < 0) {
          MoveOn = true;
          break;
        }
        PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
        if (PlaneNumber_prev - PlaneNumber > 2) {
          MoveOn = true;
          break;
        }
      }

      // Couldn't find the right hit
      if (MoveOn) continue;

      // Calculate expected derivative from upstream hit
      //double xprev = vec[NeighbourIndex].GetZ();
      //double yprev = vec[NeighbourIndex].GetNotZ();
      double xprev = vec[NeighbourIndex].GetPlaneNumber();
      double yprev = vec[NeighbourIndex].GetBarNumber();

      double grad_exp = (y-yprev)/(x-xprev);

      //std::cout << "Comparing hit " << x << "," << y << " to " << xprev << "," << yprev << " expecing gradient " << grad_exp << std::endl;
      //std::cout << "(" << PlaneNumber << ", " << vec[i].GetBarNumber() << "-" << PlaneNumber_prev << ", " << vec[NeighbourIndex].GetBarNumber() << ")" << std::endl;

      // Matching 
      for (std::vector<TMS_Hit>::iterator it = full.begin(); it != full.end(); ) {
        TMS_Hit hits = *it;
        int PlaneNumber_cand = hits.GetPlaneNumber();
        int BarNumber_cand = hits.GetBarNumber();
        // Always look upstream
        if (PlaneNumber < PlaneNumber_cand) {
          ++it;
          continue;
        }

        if (abs(BarNumber_cand - BarNumber) > 2) {
          ++it;
          continue;
        }

        // Finally also check the diagonal distance
        //double disty = BarNumber_cand-BarNumber;
        //double distx = PlaneNumber_cand-PlaneNumber;
        //if (sqrt(distx*distx + disty*disty) > 2.5) {
        //++it;
        //continue;
        //}

        // Since the candidates plane numbers are ordered in z, once we encounter a lower z, all the remaining ones also won't be mergeable, so save some time by not scanning them
        if (PlaneNumber - PlaneNumber_cand > 2) break;

        // Allow for broken track in z?
        //double xcand = hits.GetZ();
        //double ycand = hits.GetNotZ();
        double xcand = hits.GetPlaneNumber();
        double ycand = hits.GetBarNumber();

        double grad_new = (ycand-y)/(xcand-x);

        //std::cout << "testing merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;
        //std::cout << "With hit info " << hits.GetZ() << ", " << hits.GetNotZ() << std::endl;
        //std::cout << "and new gradient: " << grad_new << std::endl;

        // If gradient is within 1 and the correct sign, accept
        if (fabs(grad_new-grad_exp) <= 1.5) {
          // If the gradient is zero we shouldn't do a sign check
          // But if it's not, check the sign of the gradient doesn't flip
          if (fabs(grad_new) > TMS_Const::TMS_Small_Num && 
              fabs(grad_exp) > TMS_Const::TMS_Small_Num &&
              TMS_Utils::sgn(grad_new) != TMS_Utils::sgn(grad_exp)) {
            ++it;
            continue;
          }
          //std::cout << "Will be merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;
          // Now need to find position to insert
          // Guaranteed sorted in z
          // Need to find how many duplicates in z there are; don't necessarily want to insert right after this hit if it is followed by a hit in the same z
          // If this is the largest z, just pop to the back
          if (xcand <= vec.back().GetPlaneNumber()) {
            vec.push_back(std::move(hits));
            // Need to decrement iterator over line
            // Sometimes we've missed one hit between Hough hits
          } else {
            for (std::vector<TMS_Hit>::iterator jt = vec.begin()+i; jt != vec.end(); ++jt) {
              //double compx = (*jt).GetZ();
              // Look one behind
              //double compx1 = (*(jt-1)).GetZ();
              double compx = (*jt).GetPlaneNumber();
              // Look one behind
              double compx1 = (*(jt-1)).GetPlaneNumber();
              if (compx == x) continue;
              //std::cout << "compared " << compx << " to candidate x: " << xcand << std::endl;
              if (xcand <= compx) {
                vec.insert(jt+1, std::move(hits));
                break;
              } else if (xcand > compx && xcand <= compx1) {
                // Put in the previous
                vec.insert(jt, std::move(hits));
                break;
              } else {
                std::cout << "Didn't find anywhere to insert " << xcand << std::endl;
              }
            }
          }
          ++size;
          // Remove off the end
          it = full.erase(it);
        } else {
          it++;
        }
      }
    }

  }

