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
  zMaxHough(TMS_Const::TMS_Double_End),
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

  HoughLineU = new TF1("LinearHough", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineU->SetLineStyle(kDashed);
  HoughLineU->SetLineColor(kMagenta-9);

  HoughLineV = new TF1("LinearHough2", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineV->SetLineStyle(kDashed);
  HoughLineV->SetLineColor(kMagenta-8);

  HoughLineX = new TF1("LinearHough3", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineX->SetLineStyle(kDashed);
  HoughLineX->SetLineColor(kMagenta-7);

  HoughLineY = new TF1("LinearHough4", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLineY->SetLineStyle(kDashed);
  HoughLineY->SetLineColor(kMagenta-6);

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
  for (auto &i: HoughLinesU) {
    delete i.second;
  }
  for (auto &i: HoughLinesV) {
    delete i.second;
  }

  for (auto &i: HoughLinesX) {
    delete i.second;
  }

  for (auto &i: HoughLinesY) {
    delete i.second;
  }

  // Reset the candidate vector
  CandidatesU.clear();
  CandidatesV.clear();
  CandidatesX.clear();
  CandidatesY.clear();
  RawHits.clear();
  TotalCandidatesU.clear();
  TotalCandidatesV.clear();
  TotalCandidatesX.clear();
  TotalCandidatesY.clear();
  HoughLinesU.clear();
  HoughLinesV.clear();
  HoughLinesX.clear();
  HoughLinesY.clear();
  HoughLinesU_Upstream.clear();
  HoughLinesV_Upstream.clear();
  HoughLinesX_Upstream.clear();
  HoughLinesY_Upstream.clear();
  HoughLinesU_Downstream.clear();
  HoughLinesV_Downstream.clear();
  HoughLinesX_Downstream.clear();
  HoughLinesY_Downstream.clear();
  HoughCandidatesU.clear();
  HoughCandidatesV.clear();
  HoughCandidatesX.clear();
  HoughCandidatesY.clear();
  SortedHoughCandidatesU.clear();
  SortedHoughCandidatesV.clear();
  SortedHoughCandidatesX.clear();
  SortedHoughCandidatesY.clear();
  ClusterCandidatesU.clear();
  ClusterCandidatesV.clear();
  ClusterCandidatesX.clear();
  ClusterCandidatesY.clear();
  TrackLengthU.clear();
  TrackLengthV.clear();
  TrackLengthX.clear();
  TrackLengthY.clear();
  TrackEnergyU.clear();
  TrackEnergyV.clear();
  TrackEnergyX.clear();
  TrackEnergyY.clear();

  UHitGroup.clear();
  VHitGroup.clear();
  XHitGroup.clear();
  YHitGroup.clear();

  HoughTracks3D.clear();
}

// The generic track finder
void TMS_TrackFinder::FindTracks(TMS_Event &event) {

  ClearClass();

  // Get the raw unmerged and untracked hits
  RawHits = event.GetHits(); //GetHits needs variables transferred. Try out GetHitsRaw that doesn't need those
  //std::cout<<"Working on raw hits with n="<<RawHits.size()<<std::endl;
  //if (RawHits.size() > 0) std::cout<<"Slice "<<slice<<" has "<<RawHits.size()<<" hits."<<std::endl;

  double min_time = 1e9;
  double max_time = -1e9;
  int slice = event.GetSliceNumber();
  for (auto hit : RawHits) {
#ifdef DEBUG
    if (hit.GetPedSup()) std::cout<<"Raw hits, found a ped supped hit in slice"<<std::endl;
#endif
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
#ifdef DEBUG
    if (hit.GetPedSup()) std::cout << "Cleaned hits, found a ped supped hit in slice" << std::endl;
#endif
    clean_min_time = std::min(clean_min_time, hit.GetT());
    clean_max_time = std::max(clean_max_time, hit.GetT());
    if (hit.GetSlice() == slice) n_in_slice += 1;
    if (hit.GetSlice() != slice) n_in_wrong_slice += 1;
  }
  //if (clean_max_time - clean_min_time > (10000.0 / 52)) std::cout<<"In Reco CleanedHits, found a time range larger than expected: "<<(clean_max_time - clean_min_time)<<", min="<<clean_min_time<<", max="<<clean_max_time<<", slice="<<slice<<std::endl;
  if (n_in_wrong_slice > 0) std::cout << "Found " << n_in_wrong_slice << " in wrong slice and " << n_in_slice << " in correct slice for " << slice << std::endl;

#ifdef DEBUG
  std::cout << "Raw hits: " << RawHits.size() << std::endl;
  std::cout << "Cleaned hits: " << CleanedHits.size() << std::endl;
#endif

  // Separate planes into different groups
  // 3 degree stereo -> tilted into +3 degree in U group and into -3 degree in V group
  // 90 degree rotated -> horizontal layers in X group

  for (auto hit : CleanedHits) {
    // Sorting hits into orientation groups  
    if (hit.GetBar().GetBarType() == TMS_Bar::kUBar) {
      UHitGroup.push_back(hit);
//      std::cout << "U hit: " << hit.GetNotZ() << "," << hit.GetZ() << "," << hit.GetT() << std::endl;
    }
    else if (hit.GetBar().GetBarType() == TMS_Bar::kVBar) {
      VHitGroup.push_back(hit); 
//      std::cout << "V hit: " << hit.GetNotZ() << "," << hit.GetZ() << "," << hit.GetT() << std::endl;
    }
    else if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) {
      XHitGroup.push_back(hit);
//      std::cout << "X hit: " << hit.GetNotZ() << " | " << hit.GetZ() << "," << hit.GetT() << std::endl;
    }
    else if (hit.GetBar().GetBarType() == TMS_Bar::kYBar) {
      YHitGroup.push_back(hit);
//      std::cout << "Y hit: " << hit.GetNotZ() << " | " << hit.GetZ() << std::endl;
    }
    else {
      hit.GetBar().Print();
    }
  }
 
  if ( (UHitGroup.size() + VHitGroup.size() + XHitGroup.size() + YHitGroup.size()) != CleanedHits.size() ) {
    std::cout << "Not all hits in separated hit groups!" << std::endl;
    return;
  }
   
  // Hough transform
  if (kTrackMethod == TrackMethod::kHough) {
    // Do we first run clustering algorithm to separate hits, then hand off to A*?
    if (TMS_Manager::GetInstance().Get_Reco_HOUGH_FirstCluster()) {
      // Let's run a DBSCAN first to cluster up, then run Hough transform on clusters
      // Change to different DBSCAN parameters for this pre-clustering
      DBSCAN.SetEpsilon(TMS_Manager::GetInstance().Get_Reco_DBSCAN_PreDBDistance());
      DBSCAN.SetMinPoints(TMS_Manager::GetInstance().Get_Reco_DBSCAN_PreDBNeighbours());
      
      // actual call to DBScan
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesU = FindClusters(UHitGroup);
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesV = FindClusters(VHitGroup);
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesX = FindClusters(XHitGroup);
      std::vector<std::vector<TMS_Hit> > DBScanCandidatesY = FindClusters(YHitGroup);

      // Prepare for only running track finding (Hough transform and A* algorithm) on track-like structures
      // for this the pre-clusters found by the DBScan need to be filtered out
      std::vector<TMS_Hit> TrackCandidatesU = UHitGroup;
      std::vector<TMS_Hit> TrackCandidatesV = VHitGroup;
      std::vector<TMS_Hit> TrackCandidatesX = XHitGroup;
      std::vector<TMS_Hit> TrackCandidatesY = YHitGroup;
      
      for (auto it : DBScanCandidatesU) {
        MaskHits(TrackCandidatesU, it);
        //for (auto jt : it) std::cout << "U pre-clustered: " << jt.GetNotZ() << "," << jt.GetZ() << "," << jt.GetT() << std::endl;
      }
      //for (auto it : TrackCandidatesU) std::cout << "U track-like: " << it.GetNotZ() << "," << it.GetZ() << "," << it.GetT() << std::endl;
      for (auto it : DBScanCandidatesV) {
        MaskHits(TrackCandidatesV, it);
        //for (auto jt : it) std::cout << "V pre-clustered: " << jt.GetNotZ() << "," << jt.GetZ() << "," << jt.GetT() << std::endl;
      }
      //for (auto it : TrackCandidatesV) std::cout << "V track-like: " << it.GetNotZ() << "," << it.GetZ() << "," << it.GetT() << std::endl;
      for (auto it : DBScanCandidatesX) {
        MaskHits(TrackCandidatesX, it);
      }
      for (auto it : DBScanCandidatesY) {
        MaskHits(TrackCandidatesY, it);
      }

      // Now run the Hough transformation on the track-like structures
      HoughCandidatesU = HoughTransform(TrackCandidatesU, 'U');
      HoughCandidatesV = HoughTransform(TrackCandidatesV, 'V');
      HoughCandidatesX = HoughTransform(TrackCandidatesX, 'X');
      HoughCandidatesY = HoughTransform(TrackCandidatesY, 'Y');
#ifdef DEBUG
      std::cout << "Back extension" << std::endl;
      std::cout << "Found " << HoughCandidatesU.size() << " U simple tracks and " << HoughCandidatesV.size() << " V simple tracks" << std::endl;
#endif
      // Extend potentially back into the pre-clusters now
      for (auto it = HoughCandidatesU.begin(); it != HoughCandidatesU.end(); ++it) {
        for (auto jt : DBScanCandidatesU) {
          std::vector<TMS_Hit> HoughU = (*it);
          std::vector<TMS_Hit> preU = jt;
          
          SpatialPrio(HoughU);

          // Call the extrapolation function
          if (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_Extrapolation()) HoughU = Extrapolation(HoughU, preU);
          // Now add the new hits to the original track
          for (auto HoughIt = HoughU.begin(); HoughIt != HoughU.end(); ) {
            bool match = false;
            for (auto TrackHits = (*it).begin(); TrackHits != (*it).end(); ++TrackHits) {
              if ((*HoughIt) == (*TrackHits)) match = true;
            }
            if (!match) (*it).push_back(std::move(*HoughIt));
            else HoughIt++;
          }
        }
      }
      // Same as for U above
      for (auto it = HoughCandidatesV.begin(); it != HoughCandidatesV.end(); ++it) {
        for (auto jt : DBScanCandidatesV) {
          std::vector<TMS_Hit> HoughV = (*it);
          std::vector<TMS_Hit> preV = jt;

          SpatialPrio(HoughV);

          if (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_Extrapolation()) HoughV = Extrapolation(HoughV, preV);
          for (auto HoughIt = HoughV.begin(); HoughIt != HoughV.end(); ) {
            bool match = false;
            for (auto TrackHits = (*it).begin(); TrackHits != (*it).end(); ++TrackHits) {
              if ((*HoughIt) == (*TrackHits)) match = true;
            }
            if (!match) (*it).push_back(std::move(*HoughIt));
            else HoughIt++;
          }
        }
      }
      // Same as for U above
      for (auto it = HoughCandidatesX.begin(); it != HoughCandidatesX.end(); ++it) {
        for (auto jt : DBScanCandidatesX) { 
          std::vector<TMS_Hit> HoughX = (*it);
          std::vector<TMS_Hit> preX = jt;

          SpatialPrio(HoughX);

          if (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_Extrapolation()) HoughX = Extrapolation(HoughX, preX);
          for (auto HoughIt = HoughX.begin(); HoughIt != HoughX.end(); ) {
            bool match = false;
            for (auto TrackHits = (*it).begin(); TrackHits != (*it).end(); ++TrackHits) {
              if ((*HoughIt) == (*TrackHits)) match = true;
            }
            if (!match) (*it).push_back(std::move(*HoughIt));
            else HoughIt++;
          }
        }
      }
      // Same as for U above
      for (auto it = HoughCandidatesY.begin(); it != HoughCandidatesY.end(); ++it) {
        for (auto jt : DBScanCandidatesY) { 
          std::vector<TMS_Hit> HoughY = (*it);
          std::vector<TMS_Hit> preY = jt;

          SpatialPrio(HoughY);

          if (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_Extrapolation()) HoughY = Extrapolation(HoughY, preY);
          for (auto HoughIt = HoughY.begin(); HoughIt != HoughY.end(); ) {
            bool match = false;
            for (auto TrackHits = (*it).begin(); TrackHits != (*it).end(); ++TrackHits) {
              if ((*HoughIt) == (*TrackHits)) match = true;
            }
            if (!match) (*it).push_back(std::move(*HoughIt));
            else HoughIt++;
          }
        }
      }

      // Restore overwritten DBSCAN parameters to their previous values for final clustering
      DBSCAN.SetEpsilon(TMS_Manager::GetInstance().Get_Reco_DBSCAN_Epsilon());
      DBSCAN.SetMinPoints(TMS_Manager::GetInstance().Get_Reco_DBSCAN_MinPoints());
      // Hand over each cluster from DBSCAN to a Hough transform
      // TODO (Asa / N) no idea whether this is actually necessary for now
      /*for (std::vector<std::vector<TMS_Hit> >::iterator it = DBScanCandidatesU.begin(); it != DBScanCandidatesU.end(); ++it) {
        std::vector<TMS_Hit> hits = *it;
        std::vector<std::vector<TMS_Hit> > LinesU = HoughTransform(hits, 'U');
        for (auto jt = LinesU.begin(); jt != LinesU.end(); ++jt) {
          HoughCandidatesU.emplace_back(std::move(*jt));
        }
      }
      for (std::vector<std::vector<TMS_Hit> >::iterator it = DBScanCandidatesV.begin(); it != DBScanCandidatesV.end(); ++it) {
        std::vector<TMS_Hit> hits = *it;
	      std::vector<std::vector<TMS_Hit> > LinesV = HoughTransform(hits, 'V');
      	for (auto jt = LinesV.begin(); jt != LinesV.end(); ++jt) {
      	  HoughCandidatesV.emplace_back(std::move(*jt));
      	}
      }
      for (std::vector<std::vector<TMS_Hit> >::iterator it = DBScanCandidatesX.begin(); it != DBScanCandidatesX.end(); ++it) {
        std::vector<TMS_Hit> hits = *it;
        std::vector<std::vector<TMS_Hit> > LinesX = HoughTransform(hits, 'X');
        for (auto jt = LinesX.begin(); jt != LinesX.end(); ++jt) {
          HoughCandidatesX.emplace_back(std::move(*jt));
        }
      }*/
    } else {
      HoughCandidatesU = HoughTransform(UHitGroup, 'U');
      HoughCandidatesV = HoughTransform(VHitGroup, 'V');
      HoughCandidatesX = HoughTransform(XHitGroup, 'X');
      HoughCandidatesY = HoughTransform(YHitGroup, 'Y');
    }
  } else if (kTrackMethod == TrackMethod::kAStar) {
    BestFirstSearch(UHitGroup, 'U');
    BestFirstSearch(VHitGroup, 'V');
    BestFirstSearch(XHitGroup, 'X');
    BestFirstSearch(YHitGroup, 'Y');
  }

  std::vector<TMS_Hit> MaskedU = UHitGroup;
  std::vector<TMS_Hit> MaskedV = VHitGroup;
  std::vector<TMS_Hit> MaskedX = XHitGroup;
  std::vector<TMS_Hit> MaskedY = YHitGroup;
  // Loop over the Hough candidates
  for (auto Lines: HoughCandidatesU) {
#ifdef DEBUG
    std::cout << "Masked (u) size bef: " << MaskedU.size() << std::endl;
#endif
    MaskHits(MaskedU, Lines);
#ifdef DEBUG
    std::cout << "Masked (u) size aft: " << MaskedU.size() << std::endl;
#endif
  }

  for (auto Lines: HoughCandidatesV) {
#ifdef DEBUG
    std::cout << "Masked (v) size bef: " << MaskedV.size() << std::endl;
#endif
    MaskHits(MaskedV, Lines);
#ifdef DEBUG
    std::cout << "Masked (v) size aft: " << MaskedV.size() << std::endl;
#endif
  }

  for (auto Lines: HoughCandidatesX) {
#ifdef DEBUG
    std::cout << "Masked (x) size bef: " << MaskedX.size() << std::endl;
#endif
    MaskHits(MaskedX, Lines);
#ifdef DEBUG
    std::cout << "Masked (x) size aft: " << MaskedX.size() << std::endl;
#endif
  }

for (auto Lines: HoughCandidatesY) {
#ifdef DEBUG
  std::cout << "Masked (x) size bef: " << MaskedY.size() << std::endl;
#endif
  MaskHits(MaskedY, Lines);
#ifdef DEBUG
  std::cout << "Masked (x) size aft: " << MaskedY.size() << std::endl;
#endif
}

#ifdef DEBUG
  std::cout << "Masked (u) hits: " << MaskedU.size() << std::endl;
  std::cout << "Masked (v) hits: " << MaskedV.size() << std::endl;
  std::cout << "Masked (x) hits: " << MaskedX.size() << std::endl;
  std::cout << "Masked (x) hits: " << MaskedY.size() << std::endl;
#endif
  
  // Now we've got our tracks, refit the upstream and downstream separately with the Hough transform
  int linenoU = 0;
  for (auto Lines: HoughCandidatesU) {
    std::pair<bool, TF1*> houghline = HoughLinesU[linenoU];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
        fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      //std::cout << "Old slope: " << houghline.second->GetParameter(1) << std::endl;
      //std::cout << "New slope: " << slope << std::endl;
      //std::cout << "Old intercept: " << houghline.second->GetParameter(0) << std::endl;
      //std::cout << "New intercept: " << intercept << std::endl;

      HoughLinesU[linenoU].second->SetParameter(0, intercept);
      HoughLinesU[linenoU].second->SetParameter(1, slope);
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
    }

    double upstreamslope, upstreamintercept = 0;
    double downstreamslope, downstreamintercept = 0;
    GetHoughLine(upstream, upstreamslope, upstreamintercept);
    GetHoughLine(downstream, downstreamslope, downstreamintercept);

    std::pair<double, double> upstreamline = std::pair<double,double>(upstreamintercept, upstreamslope);
    std::pair<double, double> downstreamline = std::pair<double,double>(downstreamintercept, downstreamslope);

    HoughLinesU_Upstream.push_back(upstreamline);
    HoughLinesU_Downstream.push_back(downstreamline);

    linenoU++;
  }

  int linenoV = 0;
  for (auto Lines: HoughCandidatesV) {
    std::pair<bool, TF1*> houghline = HoughLinesV[linenoV];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
	      fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      HoughLinesV[linenoV].second->SetParameter(0, intercept);
      HoughLinesV[linenoV].second->SetParameter(1, slope);
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

   HoughLinesV_Upstream.push_back(upstreamline);
   HoughLinesV_Downstream.push_back(downstreamline);

   linenoV++;
  }

  int linenoX = 0;
  for (auto Lines: HoughCandidatesX) {
    std::pair<bool, TF1*> houghline = HoughLinesX[linenoX];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
	      fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      HoughLinesX[linenoX].second->SetParameter(0, intercept);
      HoughLinesX[linenoX].second->SetParameter(1, slope);
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

   HoughLinesX_Upstream.push_back(upstreamline);
   HoughLinesX_Downstream.push_back(downstreamline);

   linenoX++;
  }

  int linenoY = 0;
  for (auto Lines: HoughCandidatesY) {
    std::pair<bool, TF1*> houghline = HoughLinesY[linenoY];
    double slope, intercept = 0;
    GetHoughLine(Lines, slope, intercept);
    if (fabs(houghline.second->GetParameter(0) - intercept) > 1E2 ||
	      fabs(houghline.second->GetParameter(1) - slope) > 1E-2) {
      HoughLinesY[linenoY].second->SetParameter(0, intercept);
      HoughLinesY[linenoY].second->SetParameter(1, slope);
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

   HoughLinesY_Upstream.push_back(upstreamline);
   HoughLinesY_Downstream.push_back(downstreamline);

   linenoY++;
  }

  // Run a pseudo track finding for X hits, if they exist but not enough for a track to be found
  if (!XHitGroup.empty() && HoughCandidatesX.empty()) {
    FindPseudoXTrack();
  }

  // Try finding some clusters after the Hough Transform
  if (UseClustering) {
    ClusterCandidatesU = FindClusters(MaskedU);
    ClusterCandidatesV = FindClusters(MaskedV);
    ClusterCandidatesX = FindClusters(MaskedX);
    ClusterCandidatesY = FindClusters(MaskedY);
  }

  // Call TrackMatching3D
  if (!UHitGroup.empty()){
      HoughTracks3D = TrackMatching3D();
  }
  else{
      HoughTracks3D = TrackMatching3D_XY();
  }

  // Let's try to find a vertex now, just looking at most upstream point, or if there are multiple tracks let's see where they intersect
  //if (nLines > 0) {
    //if (nLines == 1) Vertex =;
    //else {
    //}
  //} else {
    //Vertex = -999;
  //}

  // Now calculate the track length and energy for each track
  for (auto it = HoughCandidatesU.begin(); it != HoughCandidatesU.end(); ++it) {
    double TrackEnergy = CalculateTrackEnergy((*it));
    double TrackLength = CalculateTrackLength((*it));
    TrackEnergyU.push_back(TrackEnergy);
    TrackLengthU.push_back(TrackLength);
  }

  for (auto it = HoughCandidatesV.begin(); it != HoughCandidatesV.end(); ++it) {
    double TrackEnergy = CalculateTrackEnergy((*it));
    double TrackLength = CalculateTrackLength((*it));
    TrackEnergyV.push_back(TrackEnergy);
    TrackLengthV.push_back(TrackLength);
  }

  for (auto it = HoughCandidatesX.begin(); it != HoughCandidatesX.end(); ++it) {
    double TrackEnergy = CalculateTrackEnergy((*it));
    double TrackLength = CalculateTrackLength((*it));
    TrackEnergyX.push_back(TrackEnergy);
    TrackLengthX.push_back(TrackLength);
  }
  for (auto it = HoughCandidatesY.begin(); it != HoughCandidatesY.end(); ++it) {
    double TrackEnergy = CalculateTrackEnergy((*it));
    double TrackLength = CalculateTrackLength((*it));
    TrackEnergyY.push_back(TrackEnergy);
    TrackLengthY.push_back(TrackLength);
  }

  //EvaluateTrackFinding(event);

  // Run Kalman filter if requested
  if (TMS_Manager::GetInstance().Get_Reco_Kalman_Run()) {
    double kalman_reco_mom, kalman_chi2_plus, kalman_chi2_minus;
    double assumed_charge = TMS_Manager::GetInstance().Get_Reco_Kalman_Assumed_Charge();
    for (auto &trk : HoughTracks3D) {

      KalmanFilter_plus = TMS_Kalman(trk.Hits, assumed_charge);
      KalmanFilter_minus = TMS_Kalman(trk.Hits, -assumed_charge);

      kalman_chi2_plus = KalmanFilter_plus.GetTrackChi2();
      trk.SetChi2_plus(kalman_chi2_plus);
      kalman_chi2_minus = KalmanFilter_minus.GetTrackChi2();
      trk.SetChi2_minus(kalman_chi2_minus);

      if (kalman_chi2_minus < kalman_chi2_plus) {
          trk.Charge_Kalman = -13;}
      else{
          trk.Charge_Kalman = 13;
      }

      trk.Charge_Kalman_curvature=KalmanFilter_plus.GetCharge_curvature();
      kalman_reco_mom = KalmanFilter_plus.GetMomentum();

      bool verbose_kalman = false;
      if (verbose_kalman) std::cout << "Kalman filter momentum: " << kalman_reco_mom << " MeV" << std::endl;
      trk.SetMomentum(kalman_reco_mom); // Fill the momentum of the TMS_Track obj

      if (verbose_kalman) std::cout << "Kalman filter start pos : " << KalmanFilter_plus.Start[0] << ", " << KalmanFilter_plus.Start[1] << ", "  << KalmanFilter_plus.Start[2] << std::endl;
      trk.SetStartPosition(KalmanFilter_plus.Start[0], KalmanFilter_plus.Start[1], KalmanFilter_plus.Start[2]); // Fill the momentum of the TMS_Track obj

      if (verbose_kalman) std::cout << "Kalman filter end pos : " << KalmanFilter_plus.End[0] << ", " << KalmanFilter_plus.End[1] << ", "  << KalmanFilter_plus.End[2] << std::endl;
      trk.SetEndPosition(KalmanFilter_plus.End[0], KalmanFilter_plus.End[1], KalmanFilter_plus.End[2]); // Fill the momentum of the TMS_Track obj

      if (verbose_kalman) std::cout << "Kalman filter start dir : " << KalmanFilter.StartDirection[0] << ", " << KalmanFilter.StartDirection[1] << ", "  << KalmanFilter.StartDirection[2] << std::endl;
      trk.SetStartDirection(KalmanFilter_plus.StartDirection[0], KalmanFilter_plus.StartDirection[1], KalmanFilter_plus.StartDirection[2]); // Fill the momentum of the TMS_Track obj

      if (verbose_kalman) std::cout << "Kalman filter end dir : " << KalmanFilter_plus.EndDirection[0] << ", " << KalmanFilter_plus.EndDirection[1] << ", "  << KalmanFilter_plus.EndDirection[2] << std::endl;
      trk.SetEndDirection(KalmanFilter_plus.EndDirection[0], KalmanFilter_plus.EndDirection[1], KalmanFilter_plus.EndDirection[2]); // Fill the momentum of the TMS_Track obj

      //When Kalman fitting the charge is not involved, so any assumed charge is okay
      trk.KalmanNodes = KalmanFilter_plus.GetKalmanNodes();
      //for chi2 tracking
      trk.KalmanNodes_plus = KalmanFilter_plus.GetKalmanNodes();
      trk.KalmanNodes_minus = KalmanFilter_minus.GetKalmanNodes();
      trk.Length = CalculateTrackLengthKalman(trk);
    }
  } 
  
  return;
}

double TMS_TrackFinder::CalculateTrackLengthKalman(const TMS_Track &Track3D) {
  // Look at the reconstructed tracks
  if (Track3D.nKalmanNodes == 0) return -999999999.;

  double final_total = 0;
  int max_n_nodes_used = 0;
  double total = 0;
  int n_nodes = 0;
  // Loop over each Kalman Node and find the track length
  for (auto it = (Track3D.KalmanNodes).rbegin(); it != (Track3D.KalmanNodes).rend() && (it+1) != (Track3D.KalmanNodes).rend(); ++it) { // turn direction around
    auto nextnode = *(it+1); //+
    // Use the geometry to calculate the track length between hits
    TVector3 point1((*it).RecoX, (*it).RecoY, (*it).z);
    TVector3 point2(nextnode.RecoX, nextnode.RecoY, nextnode.z);
    total += TMS_Geom::GetInstance().GetTrackLength(point1, point2);
    n_nodes += 1;
  }
  if (n_nodes > max_n_nodes_used) {
    final_total = total;
    max_n_nodes_used = n_nodes;
  }

  return final_total;
}

void TMS_TrackFinder::FindPseudoXTrack() {
  // Check if U and V Track exist
  if (HoughCandidatesU.empty() || HoughCandidatesV.empty() || HoughCandidatesU.size() == HoughCandidatesV.size()) return;

  // Check if X hits exist
  if (XHitGroup.empty()) return;

  // Find first and last hit of U/V track to compare X hits to
  int first_z[HoughCandidatesU.size()] = {100000};
  int last_z[HoughCandidatesU.size()] = {0};
  double first_t[HoughCandidatesU.size()] = {20000};
  double last_t[HoughCandidatesU.size()] = {-999};
  int i = 0;
  for (auto UTracks: HoughCandidatesU) {
    // Sort to make it computational less expensive
    SpatialPrio(UTracks);
    // Make sure that the orientation is always the same (inverted!)
    if (UTracks.back().GetZ() > UTracks.front().GetZ()) std::reverse(UTracks.begin(), UTracks.end());
    // Now assign the first and last hit from the U track
    if (UTracks.front().GetZ() > last_z[i]) last_z[i] = UTracks.front().GetZ();
    if (UTracks.back().GetZ() < first_z[i]) first_z[i] = UTracks.back().GetZ();
    if (UTracks.front().GetT() > last_t[i]) last_t[i] = UTracks.front().GetT();
    if (UTracks.back().GetT() < first_t[i]) first_t[i] = UTracks.back().GetT();
    ++i;
  }
  i = 0;
  for (auto VTracks: HoughCandidatesV) {
    // Sort to make it computational less expensive
    SpatialPrio(VTracks);
    // Make sure that the orientation is always the same (inverted!)
    if (VTracks.back().GetZ() > VTracks.front().GetZ()) std::reverse(VTracks.begin(), VTracks.end());
    // Now assign the first and last hit from the V track
    if (VTracks.front().GetZ() > last_z[i]) last_z[i] = VTracks.front().GetZ();
    if (VTracks.back().GetZ() < first_z[i]) first_z[i] = VTracks.back().GetZ();
    if (VTracks.front().GetT() > last_t[i]) last_t[i] = VTracks.front().GetT();
    if (VTracks.back().GetT() < first_t[i]) first_t[i] = VTracks.back().GetT();
    ++i;
  }
  for (long unsigned int i = 0; i < HoughCandidatesU.size(); ++i) {
    // Check for each X hit if between first and last U/V track hits in time plus some margin and plane number
    std::vector<TMS_Hit> CheckedXHits;
    for (auto UncheckedXHit: XHitGroup) {
      if (UncheckedXHit.GetZ() >= (first_z[i] - TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit())
          && UncheckedXHit.GetZ() <= (last_z[i] + TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit())
          && UncheckedXHit.GetT() >= (first_t[i] - TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit())
          && UncheckedXHit.GetT() <= (last_z[i] + TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit())){
        // If so, add the hit to the pseudo X track
        CheckedXHits.push_back(UncheckedXHit);
      }
    }
#ifdef DEBUG
    std::cout << "X hits added to Pseudo X track" << std::endl;
#endif
    // Push pseudo X track into HoughCandidatesX
    if (CheckedXHits.size() > 0) HoughCandidatesX.push_back(CheckedXHits);
  }
  return;
}

double TMS_TrackFinder::CompareY(TMS_Hit &UHit, TMS_Hit &VHit, TMS_Hit &XHit) {
  // Calculate UV reconstruction of y and X hit y position
//  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();
//  Double_t localU[3] = {UHit.GetNotZ(), 0, UHit.GetZ()};
//  Double_t localV[3] = {VHit.GetNotZ(), 0, VHit.GetZ()};
//  Double_t positionU[3];
//  Double_t positionV[3];
//  geom->GetCurrentMatrix()->LocalToMaster(localU, positionU);     // This gives us the position of the bars in x, y and z. Use only the y for now!
//  geom->GetCurrentMatrix()->LocalToMaster(localV, positionV);
//  TMS_Geom::GetInstance().Scale(positionU[0], positionU[1], positionU[2]);
//  TMS_Geom::GetInstance().Scale(positionV[0], positionV[1], positionV[2]);
  double UV_y = 2000;
  bool above = false;
  if ((UHit.GetNotZ() > 0 && VHit.GetNotZ() > 0) || (UHit.GetNotZ() < 0 && VHit.GetNotZ() < 0)) {
    if (std::abs(UHit.GetNotZ()) >= std::abs(VHit.GetNotZ()) && UHit.GetNotZ() < 0) above = true;
    else if (std::abs(UHit.GetNotZ()) >= std::abs(VHit.GetNotZ()) && UHit.GetNotZ() > 0) above = false;
    else if (std::abs(UHit.GetNotZ()) < std::abs(VHit.GetNotZ()) && UHit.GetNotZ() < 0) above = false;
    else if (std::abs(UHit.GetNotZ()) < std::abs(VHit.GetNotZ()) && UHit.GetNotZ() > 0) above = true;
  }
  else if (UHit.GetNotZ() > 0 && VHit.GetNotZ() < 0) above = false;
  else if (UHit.GetNotZ() < 0 && VHit.GetNotZ() > 0) above = true;

  if (above) {
    UV_y = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() * 1000 + 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ());
  } else {
    UV_y = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() * 1000 - 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ());
  }
  //double UV_y = 244 - 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ()); //-1350/-2949
  
  double X_y = XHit.GetNotZ();

//  std::cout<<"UV_y:"<<UV_y<<"X_y:"<<X_y<<"difference:"<<std::abs(UV_y - X_y)<<std::endl;
  double compared_y = -999999;
  // If distance between expectation of UV reconstruction of y and X hit y position is too big
  if (std::abs(UV_y - X_y) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_YDifference()) {
    // Use UV reconstruction of y for U/V hits
  
//      std::cout<<"Use UV_y"<<std::endl;
      compared_y = UV_y;
  } 
  // Otherwise use X hit y position for y
  else {
      compared_y = X_y;
//      std::cout<<"Use X_y"<<std::endl;
  }
//  compared_y = UV_y;

  return compared_y;
}

std::vector<TMS_Track> TMS_TrackFinder::TrackMatching3D() {
//#ifdef DEBUG
  std::cout << "3D matching" << std::endl;
  std::cout << "size Candidates: U: " << HoughCandidatesU.size() << " | V: " << HoughCandidatesV.size() << " | X: " << HoughCandidatesX.size() << std::endl;
//#endif

  // Sorting the candidate tracks descending by their hit number
  if (HoughCandidatesU.size() > 1) {
    std::sort(HoughCandidatesU.begin(), HoughCandidatesU.end(), SortByHitNumber);
    SortedHoughCandidatesU = HoughCandidatesU;
  } else {
    SortedHoughCandidatesU = HoughCandidatesU;
  }
  if (HoughCandidatesV.size() > 1) {
    std::sort(HoughCandidatesV.begin(), HoughCandidatesV.end(), SortByHitNumber);
    SortedHoughCandidatesV = HoughCandidatesV;
  } else {
    SortedHoughCandidatesV = HoughCandidatesV;
  }
  if (HoughCandidatesX.size() > 1) {
    std::sort(HoughCandidatesX.begin(), HoughCandidatesX.end(), SortByHitNumber);
    SortedHoughCandidatesX = HoughCandidatesX;
  } else {
    SortedHoughCandidatesX = HoughCandidatesX;
  }

  std::vector<TMS_Track> returned;

  bool TimeSlicing = TMS_Manager::GetInstance().Get_Reco_TIME_RunTimeSlicer();

  // 3D matching of tracks
  std::vector<std::vector<TMS_Hit> >::iterator Uhelper = SortedHoughCandidatesU.begin();
  std::vector<std::vector<TMS_Hit> >::iterator Vhelper = SortedHoughCandidatesV.begin();
  std::vector<std::vector<TMS_Hit> >::iterator helper = SortedHoughCandidatesX.begin();

  while (Uhelper != SortedHoughCandidatesU.end()) {
    if (SortedHoughCandidatesV.empty()) {
#ifdef DEBUG
      std::cout << "Not enough V tracks" << std::endl;
#endif
      break;
    }

    std::vector<TMS_Hit> UTracks = *Uhelper;
    std::vector<TMS_Hit> VTracks = *Vhelper;

    // Run with matching of X tracks, if X tracks exist. Otherwise match without
    bool Xrun = true;
#ifdef DEBUG
    std::cout << "No X tracks? " << HoughCandidatesX.empty() << std::endl;
#endif

    if (SortedHoughCandidatesX.empty()) Xrun = false;
    std::vector<TMS_Hit> XTracks = *Uhelper;
    if (Xrun) {
      XTracks = *helper;
      if (XTracks.empty()) Xrun = false;
    }
        

    // Run spatial prio just because one last time
    SpatialPrio(UTracks);
    SpatialPrio(VTracks);
    if (Xrun) {
      SpatialPrio(XTracks);
      std::cout << "SpatialPrio X" << std::endl;
    }
#ifdef DEBUG
    std::cout << "U track: " << UTracks.size() << std::endl;
    std::cout << "UTrack back: " << UTracks.front().GetPlaneNumber() << " | " << UTracks.front().GetBarNumber() << " | " << UTracks.front().GetT() << " front: " << UTracks.back().GetPlaneNumber() << " | " << UTracks.back().GetBarNumber() << " | " << UTracks.back().GetT() << std::endl;
    std::cout << "V track: " << VTracks.size() << std::endl;
    std::cout << "VTrack back: " << VTracks.front().GetPlaneNumber() << " | " << VTracks.front().GetBarNumber() << " | " << VTracks.front().GetT() << " front: " << VTracks.back().GetPlaneNumber() << " | " << VTracks.back().GetBarNumber() << " | " << VTracks.back().GetT() << std::endl;
    if (Xrun) std::cout << "X track: " << XTracks.size() << std::endl;
    if (Xrun) std::cout << "XTrack back: " << XTracks.front().GetPlaneNumber() << " | " << XTracks.front().GetBarNumber() << " | " << XTracks.front().GetT() << " front: " << XTracks.back().GetPlaneNumber() << " | " << XTracks.back().GetBarNumber() << " | " << XTracks.back().GetT() << std::endl;
#endif

    // Conditions for close enough tracks: within +/-3 plane numbers, +/-12 bar numbers and in same time slice within 30ns (Asa note: changed to 15 ns for now and just one of time and plane number condition needs to be met)
    bool back_match = false;
    bool front_match = false;
    bool Xback_match = false;
    bool Xfront_match = false;
    int strike = 0;
    // Check for the front and back whether both conditions (plane limit and time limit) are met or not
    if (TimeSlicing) {
      // front of simple tracks: set strike to 1
      if (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit()
          || std::abs(UTracks.front().GetT() - VTracks.front().GetT()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit()) strike = 1;
      // end of simple tracks: add 2 to strike
      if (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit()
          || std::abs(UTracks.back().GetT() - VTracks.back().GetT()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit()) strike += 2;
    }  
    if (Xrun && TimeSlicing) {
      bool bar_front = (std::abs(UTracks.front().GetBarNumber() - VTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      bool bar_back = (std::abs(UTracks.back().GetBarNumber() - VTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      bool slice_front = (UTracks.front().GetSlice() == VTracks.front().GetSlice());
      bool slice_back = (UTracks.back().GetSlice() == VTracks.back().GetSlice());
      bool time_front = (std::abs(UTracks.front().GetT() - VTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      bool time_back = (std::abs(UTracks.back().GetT() - VTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      bool plane_front = (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
      bool plane_back = (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
      if (strike == 0) {
        // front and end match
        back_match = (bar_back && slice_back && time_back && plane_back);
        front_match = (bar_front && slice_front && time_front && plane_front);
      } else if (strike == 1) {
        // front needs exemption, end matches
        back_match = (bar_back && slice_back && time_back && plane_back);
        front_match = (bar_front && slice_front && (time_front || plane_front));
      } else if (strike == 2) {
        // front matches, end needs exemption
        back_match = (bar_back && slice_back && (time_back || plane_back));
        front_match = (bar_front && slice_front && time_front && plane_front);
      } //else if (strike == 3) // both do not match -> ignore
      /*back_match = (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit() 
          && std::abs(UTracks.front().GetBarNumber() - VTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit() 
          && UTracks.front().GetSlice() == VTracks.front().GetSlice()
          && std::abs(UTracks.front().GetT() - VTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      front_match = (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit() 
          && std::abs(UTracks.back().GetBarNumber() - VTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit() 
          && UTracks.back().GetSlice() == VTracks.back().GetSlice()
          && std::abs(UTracks.back().GetT() - VTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());*/
      Xback_match = (std::abs(UTracks.front().GetT() - XTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit()
          || std::abs(VTracks.front().GetT() - XTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit());
      Xfront_match = (std::abs(UTracks.back().GetT() - XTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit()
          || std::abs(VTracks.back().GetT() - XTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_XTimeLimit());
      if (SortedHoughCandidatesU.size() > 1 || SortedHoughCandidatesV.size() > 1) {
        if (Xback_match && Xfront_match) {
          // check if this would actually be the best match for the X track
          // calculate time differences of start and end for X to U and V
          std::vector<float> time_difference_frontU(SortedHoughCandidatesU.size());
          std::vector<float> time_difference_frontV(SortedHoughCandidatesV.size());
          std::vector<float> time_difference_backU(SortedHoughCandidatesU.size());
          std::vector<float> time_difference_backV(SortedHoughCandidatesV.size());
          int time_it = 0;
          for (std::vector<std::vector<TMS_Hit> >::iterator tracks = SortedHoughCandidatesU.begin(); tracks != SortedHoughCandidatesU.end(); ++tracks) {
            SpatialPrio(*tracks);
            time_difference_frontU[time_it] = std::abs((*tracks).front().GetT() - XTracks.front().GetT());
            time_difference_backU[time_it] = std::abs((*tracks).back().GetT() - XTracks.back().GetT());
            ++time_it;
          }
          time_it = 0;
          for (std::vector<std::vector<TMS_Hit> >::iterator tracks = SortedHoughCandidatesV.begin(); tracks != SortedHoughCandidatesV.end(); ++tracks) {
            SpatialPrio(*tracks);
            time_difference_frontV[time_it] = std::abs((*tracks).front().GetT() - XTracks.front().GetT());
            time_difference_backV[time_it] = std::abs((*tracks).back().GetT() - XTracks.back().GetT());
            ++time_it;
          }
          // only allow matching, if time differences are the lowest
          // find lowest element
          float min_frontU = 9999.;
          float min_frontV = 9999.;
          float min_backU = 9999.;
          float min_backV = 9999.;
          float sum_min_U = 9999.;
          float sum_min_V = 9999.;
          for (int it = 0; it < time_difference_frontU.size(); ++it) {
            if (min_frontU > time_difference_frontU[it]) min_frontU = time_difference_frontU[it];
            if (min_backU > time_difference_backU[it]) min_backU = time_difference_backU[it];
            if (sum_min_U > time_difference_frontU[it] + time_difference_backU[it]) sum_min_U = time_difference_frontU[it] + time_difference_backU[it];
          }
          for (int it = 0; it < time_difference_frontV.size(); ++it) {
            if (min_frontV > time_difference_frontV[it]) min_frontV = time_difference_frontV[it];
            if (min_backV > time_difference_backV[it]) min_backV = time_difference_backV[it];
            if (sum_min_V > time_difference_frontV[it] + time_difference_backV[it]) sum_min_V = time_difference_frontV[it] + time_difference_backV[it];
          }
          if ((min_frontU != std::abs(UTracks.front().GetT() - XTracks.front().GetT())
                || min_frontV != std::abs(VTracks.front().GetT() - XTracks.front().GetT()))
              && sum_min_U == (std::abs(UTracks.front().GetT() - XTracks.front().GetT()) + std::abs(UTracks.back().GetT() - XTracks.back().GetT()))
              && sum_min_V == (std::abs(VTracks.front().GetT() - XTracks.front().GetT()) + std::abs(VTracks.back().GetT() - XTracks.back().GetT()))) 
            Xback_match = false;
          if ((min_backU != std::abs(UTracks.front().GetT() - XTracks.front().GetT())
              || min_backV != std::abs(VTracks.back().GetT() - XTracks.back().GetT()))
              && sum_min_U == (std::abs(UTracks.front().GetT() - XTracks.front().GetT()) + std::abs(UTracks.back().GetT() - XTracks.back().GetT()))
              && sum_min_V == (std::abs(VTracks.front().GetT() - XTracks.front().GetT()) + std::abs(VTracks.back().GetT() - XTracks.back().GetT())))
            Xfront_match = false;
        }
      }
    } else if (!Xrun && TimeSlicing) {
      bool bar_front = (std::abs(UTracks.front().GetBarNumber() - VTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      bool bar_back = (std::abs(UTracks.back().GetBarNumber() - VTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      bool slice_front = (UTracks.front().GetSlice() == VTracks.front().GetSlice());
      bool slice_back = (UTracks.back().GetSlice() == VTracks.back().GetSlice());
      bool time_front = (std::abs(UTracks.front().GetT() - VTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      bool time_back = (std::abs(UTracks.back().GetT() - VTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      bool plane_front = (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
      bool plane_back = (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
      if (strike == 0) {
        // front and end match
        back_match = (bar_back && slice_back && time_back && plane_back);
        front_match = (bar_front && slice_front && time_front && plane_front);
      } else if (strike == 1) {
        // front needs exemption, end matches
        back_match = (bar_back && slice_back && time_back && plane_back);
        front_match = (bar_front && slice_front && (time_front || plane_front));
      } else if (strike == 2) {
        // front matches, end needs exemption
        back_match = (bar_back && slice_back && (time_back || plane_back));
        front_match = (bar_front && slice_front && time_front && plane_front);
      } //else if (strike == 3) // both do not match -> ignore
      /*back_match = (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit() 
          && std::abs(UTracks.front().GetBarNumber() - VTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit() 
          && UTracks.front().GetSlice() == VTracks.front().GetSlice()
          && std::abs(UTracks.front().GetT() - VTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
      front_match = (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit() 
          && std::abs(UTracks.back().GetBarNumber() - VTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit() 
          && UTracks.back().GetSlice() == VTracks.back().GetSlice()
          && std::abs(UTracks.back().GetT() - VTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());*/
    } else {
      back_match = (std::abs(UTracks.front().GetPlaneNumber() - VTracks.front().GetPlaneNumber()) < TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit() 
          && std::abs(UTracks.front().GetBarNumber() - VTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      front_match = (std::abs(UTracks.back().GetPlaneNumber() - VTracks.back().GetPlaneNumber()) < TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit()
          && std::abs(UTracks.back().GetBarNumber() - VTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
      Xback_match = true; // matching is pretty bad with time slicing turned off. Enable at least some matching with X 
      Xfront_match = true;
    }
    if (back_match && front_match) { 
      // end condition THIS IS ACTUALLY THE START CONDITION
      TMS_Track aTrack;

      // Make sure that the hits are in the correct order
      if (UTracks.back().GetZ() > UTracks.front().GetZ()) std::reverse(UTracks.begin(), UTracks.end());
      if (VTracks.back().GetZ() > VTracks.front().GetZ()) std::reverse(VTracks.begin(), VTracks.end());
      if (Xrun && (XTracks.back().GetZ() > XTracks.front().GetZ())) std::reverse(XTracks.begin(), XTracks.end());
#ifdef DEBUG          
      std::cout << "UTrack FRONT: " << UTracks.back().GetNotZ() << ", " << UTracks.back().GetZ() << " (" << UTracks.back().GetPlaneNumber() << ") BACK: " << UTracks.front().GetNotZ() << ", " << UTracks.front().GetZ() << " (" << UTracks.front().GetPlaneNumber() << ")" << std::endl;
      std::cout << "VTrack FRONT: " << VTracks.back().GetNotZ() << ", " << VTracks.back().GetZ() << " (" << VTracks.back().GetPlaneNumber() << ") BACK: " << VTracks.front().GetNotZ() << ", " << VTracks.front().GetZ() << " (" << VTracks.front().GetPlaneNumber() << ")" << std::endl;
      if (Xrun) std::cout << "XTrack FRONT: " << XTracks.back().GetNotZ() << ", " << XTracks.back().GetZ() << " (" << XTracks.back().GetPlaneNumber() << ") Back: " << XTracks.front().GetNotZ() << ", " << XTracks.front().GetZ() << " (" << XTracks.front().GetPlaneNumber() << ")" << std::endl;
#endif          
      // If same plane number for start but different for end (U/V)
      if (UTracks.front().GetPlaneNumber() != VTracks.front().GetPlaneNumber()) {
        // If UTrack ends after VTrack
        if (UTracks.front().GetZ() > VTracks.front().GetZ()) {
          bool stereo_view = true;
          if (Xrun && Xback_match && Xfront_match) {
            // Check if X track might be the end of the overall track
            if (UTracks.front().GetZ() < XTracks.front().GetZ() || VTracks.front().GetZ() < XTracks.front().GetZ()) stereo_view = false;
            //if (UTracks.front().GetPlaneNumber() < 20 && std::abs(UTracks.front().GetZ() - XTracks.front().GetZ()) < 65) stereo_view = false;
            //else if (UTracks.front().GetPlaneNumber() >= 20 && std::abs(UTracks.front().GetZ() - XTracks.front().GetZ()) < 90) stereo_view = false;
          }
          if (stereo_view) {
            CalculateRecoY(UTracks.front(), UTracks.front(), VTracks.front());
            CalculateRecoX(UTracks.front(), VTracks.front(), UTracks.front());
            aTrack.End[0] = UTracks.front().GetRecoX();
            aTrack.End[1] = UTracks.front().GetRecoY();
            aTrack.End[2] = UTracks.front().GetZ();
            aTrack.End[3] = UTracks.front().GetT();
#ifdef DEBUG                
            std::cout << "UTrack ends after VTrack" << std::endl;
#endif            
            (aTrack.Hits).push_back(UTracks.front());
          } else {
            if (UTracks.front().GetZ() > XTracks.front().GetZ()) {
              CalculateRecoX(UTracks.front(), VTracks.front(), UTracks.front());
              UTracks.front().SetRecoY(CompareY(UTracks.front(), VTracks.front(), XTracks.front()));
              aTrack.End[0] = 0.5 * (UTracks.front().GetNotZ() + VTracks.front().GetNotZ());
              aTrack.End[1] = CompareY(UTracks.front(), VTracks.front(), XTracks.front());//XTracks.front().GetNotZ();
              aTrack.End[2] = UTracks.front().GetZ();
              aTrack.End[3] = UTracks.front().GetT();
#ifdef DEBUG
              std::cout << "UTrack ends after VTrack and XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(UTracks.front());
            } else if (UTracks.front().GetZ() < XTracks.front().GetZ()) {
              CalculateRecoX(UTracks.front(), VTracks.front(), XTracks.front());
              XTracks.front().SetRecoY(XTracks.front().GetNotZ());
              aTrack.End[0] = 0.5 * (UTracks.front().GetNotZ() + VTracks.front().GetNotZ());
              aTrack.End[1] = XTracks.front().GetNotZ();
              aTrack.End[2] = XTracks.front().GetZ();
              aTrack.End[3] = XTracks.front().GetT();
#ifdef DEBUG
              std::cout << "UTrack ends after VTrack, before XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(XTracks.front());
            } 
          } 
        } // If UTrack ends before VTrack
        else if (UTracks.front().GetZ() < VTracks.front().GetZ()) {
          bool stereo_view = true;
          if (Xrun && Xback_match && Xfront_match) {
            // Check if the X track might be the end of the overall track
            if (UTracks.front().GetZ() < XTracks.front().GetZ() || VTracks.front().GetZ() < XTracks.front().GetZ()) stereo_view = false;
            //if (VTracks.front().GetPlaneNumber() < 20 && std::abs(VTracks.front().GetZ() - XTracks.front().GetZ()) < 65) stereo_view = false;
            //else if (VTracks.front().GetPlaneNumber() >= 20 && std::abs(VTracks.front().GetZ() - XTracks.front().GetZ()) < 90) stereo_view = false;
          }
          if (stereo_view) {
            CalculateRecoY(VTracks.front(), UTracks.front(), VTracks.front());
            CalculateRecoX(UTracks.front(), VTracks.front(), VTracks.front());
            aTrack.End[0] = VTracks.front().GetRecoX();
            aTrack.End[1] = VTracks.front().GetRecoY();
            aTrack.End[2] = VTracks.front().GetZ();
            aTrack.End[3] = VTracks.front().GetT();
#ifdef DEBUG                
            std::cout << "VTracks ends after UTrack, no XTrack" << std::endl;
#endif              
            (aTrack.Hits).push_back(VTracks.front());
          } else {
            if (VTracks.front().GetZ() > XTracks.front().GetZ()) {
              CalculateRecoX(VTracks.front(), UTracks.front(), VTracks.front());
              VTracks.front().SetRecoY(CompareY(UTracks.front(), VTracks.front(), XTracks.front()));
              aTrack.End[0] = 0.5 * (VTracks.front().GetNotZ() + UTracks.front().GetNotZ());
              aTrack.End[1] = CompareY(UTracks.front(), VTracks.front(), XTracks.front());//XTracks.front().GetNotZ();
              aTrack.End[2] = VTracks.front().GetZ();
              aTrack.End[3] = VTracks.front().GetT();
#ifdef DEBUG
              std::cout << "VTracks ends after UTrack and XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(VTracks.front());
            } else if (VTracks.front().GetZ() < XTracks.front().GetZ()) {
              XTracks.front().SetRecoX(0.5 * (UTracks.front().GetNotZ() + VTracks.front().GetNotZ()));
              CalculateRecoX(UTracks.front(), VTracks.front(), XTracks.front());
              XTracks.front().SetRecoY(XTracks.front().GetNotZ());
              aTrack.End[0] = 0.5 * (VTracks.front().GetNotZ() + UTracks.front().GetNotZ());
              aTrack.End[1] = XTracks.front().GetNotZ();
              aTrack.End[2] = XTracks.front().GetZ();
              aTrack.End[3] = XTracks.front().GetT();
#ifdef DEBUG
              std::cout << "VTracks ends after UTrack, before XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(XTracks.front());
            }
          } 
        }
      } 
      // If different plane number for start but same for end (U/V)
      if (UTracks.back().GetPlaneNumber() != VTracks.back().GetPlaneNumber()) {
        // If UTrack starts after VTrack  
        if (UTracks.back().GetZ() > VTracks.back().GetZ()) {
          bool stereo_view = true;
          if (Xrun && Xback_match && Xfront_match) {
            // Check if the X track might be the overall start of the track
            if (UTracks.back().GetZ() > XTracks.back().GetZ() || VTracks.back().GetZ() > XTracks.back().GetZ()) stereo_view = false;
            //if (VTracks.back().GetPlaneNumber() < 20 && std::abs(VTracks.back().GetZ() - XTracks.back().GetZ()) < 65) stereo_view = false;
            //else if (VTracks.back().GetPlaneNumber() >= 20 && std::abs(VTracks.back().GetZ() - XTracks.back().GetZ()) < 90) stereo_view = false;
          }
          if (stereo_view) {
            CalculateRecoY(VTracks.back(), UTracks.back(), VTracks.back());
            CalculateRecoX(UTracks.back(), VTracks.back(), VTracks.back());
            aTrack.Start[0] = VTracks.back().GetRecoX();
            aTrack.Start[1] = VTracks.back().GetRecoY();
            aTrack.Start[2] = VTracks.back().GetZ();
            aTrack.Start[3] = VTracks.back().GetT();
#ifdef DEBUG
            std::cout << "VTrack starts before UTrack, no XTrack" << std::endl;
#endif                
            (aTrack.Hits).push_back(VTracks.back());
          } else {
            if (VTracks.back().GetZ() < XTracks.back().GetZ()) {
              VTracks.back().SetRecoX(0.5 * (VTracks.back().GetNotZ() + UTracks.back().GetNotZ()));
              VTracks.back().SetRecoY(CompareY(UTracks.back(), VTracks.back(), XTracks.back()));
              aTrack.Start[0] = 0.5 * (VTracks.back().GetNotZ() + UTracks.back().GetNotZ());
              aTrack.Start[1] = CompareY(UTracks.back(), VTracks.back(), XTracks.back());//XTracks.back().GetNotZ();
              aTrack.Start[2] = VTracks.back().GetZ();
              aTrack.Start[3] = VTracks.back().GetT();
#ifdef DEBUG
              std::cout << "VTrack starts before UTrack and XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(VTracks.back());
            } else if (VTracks.back().GetZ() > XTracks.back().GetZ()) {
              CalculateRecoX(VTracks.back(), UTracks.back(), XTracks.back());
              XTracks.back().SetRecoY(XTracks.back().GetNotZ());
              aTrack.Start[0] = 0.5 * (VTracks.back().GetNotZ() + UTracks.back().GetNotZ());
              aTrack.Start[1] = XTracks.back().GetNotZ();
              aTrack.Start[2] = XTracks.back().GetZ();
              aTrack.Start[3] = XTracks.back().GetT();
#ifdef DEBUG
              std::cout << "VTrack starts before UTrack, after XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(XTracks.back());
            }
          }
        } 
        else if (UTracks.back().GetZ() < VTracks.back().GetZ()) {
          bool stereo_view = true;
          if (Xrun && Xback_match && Xfront_match) {
            // Check if the X track might be the overall start of the track
            if (UTracks.back().GetZ() > XTracks.back().GetZ() || VTracks.back().GetZ() > XTracks.back().GetZ()) stereo_view = false;
            //if (UTracks.back().GetPlaneNumber() < 20 && std::abs(UTracks.back().GetZ() - XTracks.back().GetZ()) < 65) stereo_view = false;
            //else if (UTracks.back().GetPlaneNumber() >= 20 && std::abs(UTracks.back().GetZ() - XTracks.back().GetZ()) < 90) stereo_view = false;
          }
          if (stereo_view) {
            CalculateRecoY(UTracks.back(), UTracks.back(), VTracks.back());
            CalculateRecoX(UTracks.back(), VTracks.back(), UTracks.back());
            aTrack.Start[0] = UTracks.back().GetRecoX();
            aTrack.Start[1] = UTracks.back().GetRecoY();
            aTrack.Start[2] = UTracks.back().GetZ();
            aTrack.Start[3] = UTracks.back().GetT();
#ifdef DEBUG                
            std::cout << "UTrack starts before VTrack, no XTrack" << std::endl;
#endif                
            (aTrack.Hits).push_back(UTracks.back());
          } else {
            if (UTracks.back().GetZ() < XTracks.back().GetZ()) {
              CalculateRecoX(VTracks.back(), UTracks.back(), UTracks.back());
              UTracks.back().SetRecoY(CompareY(UTracks.back(), VTracks.back(), XTracks.back()));
              aTrack.Start[0] = 0.5 * (UTracks.back().GetNotZ() + VTracks.back().GetNotZ());
              aTrack.Start[1] = CompareY(UTracks.back(), VTracks.back(), XTracks.back());//XTracks.back().GetNotZ();
              aTrack.Start[2] = UTracks.back().GetZ();
              aTrack.Start[3] = UTracks.back().GetT();
#ifdef DEBUG
              std::cout << "UTrack starts before VTrack and XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(UTracks.back());
            } else if (UTracks.back().GetZ() > XTracks.back().GetZ()) {
              CalculateRecoX(VTracks.back(), UTracks.back(), XTracks.back());
              XTracks.back().SetRecoY(XTracks.back().GetNotZ());
              aTrack.Start[0] = 0.5 * (UTracks.back().GetNotZ() + VTracks.back().GetNotZ());
              aTrack.Start[1] = XTracks.back().GetNotZ();
              aTrack.Start[2] = XTracks.back().GetZ();
              aTrack.Start[3] = XTracks.back().GetT();
#ifdef DEBUG
              std::cout << "UTrack starts before VTrack, after XTrack" << std::endl;
#endif
              (aTrack.Hits).push_back(XTracks.back());
            }
          } 
        } 
      }
      if (Xrun && Xback_match && Xfront_match) {
        // If different plane number for start (X)
        if (XTracks.back().GetZ() < UTracks.back().GetZ() && XTracks.back().GetZ() < VTracks.back().GetZ()) {
          CalculateRecoX(VTracks.back(), UTracks.back(), XTracks.back());
          XTracks.back().SetRecoY(XTracks.back().GetNotZ());
          aTrack.Start[0] = 0.5 * (UTracks.back().GetNotZ() + VTracks.back().GetNotZ());
          aTrack.Start[1] = XTracks.back().GetNotZ();
          aTrack.Start[2] = XTracks.back().GetZ();
          aTrack.Start[3] = XTracks.back().GetT();
#ifdef DEBUG
          std::cout << "XTrack starts before all" << std::endl;
#endif
          (aTrack.Hits).push_back(XTracks.back());
        }
        // If different plane number for end (X)
        if (XTracks.front().GetZ() > UTracks.front().GetZ() && XTracks.front().GetZ() > VTracks.front().GetZ()) {
          CalculateRecoX(UTracks.front(), VTracks.front(), XTracks.front());
          XTracks.front().SetRecoY(XTracks.front().GetNotZ());
          aTrack.End[0] = 0.5 * (UTracks.front().GetNotZ() + VTracks.front().GetNotZ());
          aTrack.End[1] = XTracks.front().GetNotZ();
          aTrack.End[2] = XTracks.front().GetZ();
          aTrack.End[3] = XTracks.front().GetT();
#ifdef DEBUG
          std::cout << "XTrack ends after all" << std::endl;
#endif
          (aTrack.Hits).push_back(XTracks.front());
        }
      }

      // Add hits to track
      int itU = UTracks.size() - 1;
      int itV = VTracks.size() - 1;
      int itX = XTracks.size() - 1;

//      std::cout<<"itV:"<<itV<<", itU:"<<itU<<", itX:"<<itX<<std::endl;

      bool sane = (itU > 0 || itV > 0);
      if (Xrun && Xback_match && Xfront_match) sane = (itU > 0 || itV > 0 || itX > 0);

      while (sane) {  // Track seems to be backwards, so run adding of hits backwards
//        std::cout <<"itV:"<<itV<<" ,X: " << VTracks[itV].GetNotZ() <<  ", Z: " << VTracks[itV].GetZ() <<", PlaneNumber: " << VTracks[itV].GetPlaneNumber()<< std::endl;
//        std::cout <<"itU:"<<itU<<" ,X: " << UTracks[itU].GetNotZ() <<  ", Z: " << UTracks[itU].GetZ() <<", PlaneNumber: " << UTracks[itU].GetPlaneNumber()<< std::endl;
//        if (Xrun) std::cout <<"itX:"<<itX<<" ,Y: " << XTracks[itX].GetNotZ() <<  ", Z: " << XTracks[itX].GetZ() <<", PlaneNumber: " << XTracks[itX].GetPlaneNumber()<< std::endl;
        // Sanity check for hits being in the detector volume (x and z)
        bool hit_outside = false;
        if (std::abs(UTracks[itU].GetNotZ()) > 4000.0 || UTracks[itU].GetNotZ() == 0. 
            || UTracks[itU].GetZ() < 11000 || UTracks[itU].GetZ() > 20000) {
          --itU;
          hit_outside = true;
        }
        if (std::abs(VTracks[itV].GetNotZ()) > 4000.0 or VTracks[itV].GetNotZ() == 0.
            || VTracks[itV].GetZ() < 11000 || VTracks[itV].GetZ() > 20000) {
          --itV;
          hit_outside = true;
        }
        if (hit_outside) continue;

        // Stereo check
        bool stereo_view = true;
//        bool stereo_view = false;
        if (Xrun && Xback_match && Xfront_match) {
          // Check if a neighbouring hit is from a X layer
          //if ((UTracks[itU].GetPlaneNumber() < TMS_Const::TMS_nThinPlanes / 2 && std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Thin_gap * 2 + 10))
          //    || (VTracks[itV].GetPlaneNumber() < TMS_Const::TMS_nThinPlanes / 2 && std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Thin_gap * 2 + 10))) {
          //  stereo_view = false;
          //} else if ((UTracks[itU].GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 && std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Thick_gap * 2 + 10))
          //    || (VTracks[itV].GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 && std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Thick_gap * 2 + 10))) {
          //  stereo_view = false;
          //} 
          //if ((UTracks[itU].GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 + TMS_Const::TMS_nThickPlanes / 2 && std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Double_gap * 2 + 10))
          //    || (VTracks[itV].GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 + TMS_Const::TMS_nThickPlanes / 2 && std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ() < TMS_Const::TMS_Double_gap * 2 + 10))) {
          //  stereo_view = false;
          //}
          if (((UTracks[itU].GetZ() <= TMS_Const::TMS_Thick_Start && XTracks[itX].GetZ() <= TMS_Const::TMS_Thick_Start) && (std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Thin_gap * 2 + 10))
              || ((VTracks[itV].GetZ() <= TMS_Const::TMS_Thick_Start) && (std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Thin_gap * 2 + 10))) {
            stereo_view = false;
          } else if (((UTracks[itU].GetZ() <= TMS_Const::TMS_Double_Start && XTracks[itX].GetZ() <= TMS_Const::TMS_Double_Start) && (std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Thick_gap * 2 + 10))
              || ((VTracks[itV].GetZ() <= TMS_Const::TMS_Double_Start) && (std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Thick_gap * 2 + 10))) {
            stereo_view = false;
          } else if (((UTracks[itU].GetZ() >= TMS_Const::TMS_Double_Start && XTracks[itX].GetZ() >= TMS_Const::TMS_Double_Start) && (std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Double_gap * 2 + 10))
              || ((VTracks[itV].GetZ() >= TMS_Const::TMS_Double_Start) && (std::abs(VTracks[itV].GetZ() - XTracks[itX].GetZ()) < TMS_Const::TMS_Double_gap * 2 + 10))) {
            stereo_view = false;
          }
          if (itX > 0 && itU == 0 && itV == 0) stereo_view = false;
        }
#ifdef DEBUG
        std::cout << "itU: " << itU << " | itV: " << itV << std::endl;
        std::cout << "U: " << UTracks[itU].GetNotZ() << " / " << UTracks[itU].GetZ() << " | " << UTracks[itU].GetT() << " | V: " << VTracks[itV].GetNotZ() << " / " << VTracks[itV].GetZ() << " | " << VTracks[itV].GetT() << std::endl;
        if (Xrun) std::cout << "itX: " << itX <<  "| X: " << XTracks[itX].GetNotZ() << " / " << XTracks[itX].GetZ() << " | " << XTracks[itX].GetT() << std::endl;
#endif 
        //std::cout << "Stereo_view:" << stereo_view << std::endl;
        if ((UTracks[itU]).GetPlaneNumber() == (VTracks[itV]).GetPlaneNumber()) {
          // Calculate Y info from bar crossing
          if (stereo_view) {
            CalculateRecoY(UTracks[itU], UTracks[itU], VTracks[itV]);
            CalculateRecoY(VTracks[itV], UTracks[itU], VTracks[itV]);
            CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
            CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
#ifdef DEBUG
            std::cout << "same" << std::endl;
            std::cout << "Hit: " << UTracks[itU].GetRecoX() << " | " << UTracks[itU].GetRecoY() << " | " << UTracks[itU].GetZ() << " than: " << VTracks[itV].GetRecoX() << " | " << VTracks[itV].GetRecoY() << " | " << VTracks[itV].GetZ() << std::endl;
#endif
            // Handling cases of two hits in same plane to be matched
            // By adding a loop into these statements one could also take care of more than 2 hits in the same plane/with the same z coordinate

            if (itU > 0 && UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
              CalculateRecoY(UTracks[itU - 1], UTracks[itU - 1], VTracks[itV]);
              CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
              (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
              if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
            }
            if (itV > 0 && VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
              CalculateRecoY(VTracks[itV - 1], UTracks[itU], VTracks[itV - 1]);
              CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
              (aTrack.Hits).push_back(VTracks[itV]);
              if (itV > 0) --itV;
            }
            (aTrack.Hits).push_back(UTracks[itU]);
            (aTrack.Hits).push_back(VTracks[itV]);
          } else { 
            float max_distance = 0;
            if (UTracks[itU].GetZ() <= TMS_Const::TMS_Thick_Start) max_distance = TMS_Const::TMS_Thin_gap;
            else if (UTracks[itU].GetZ() <= TMS_Const::TMS_Double_Start) max_distance = TMS_Const::TMS_Thick_gap;
            else if (UTracks[itU].GetZ() <= TMS_Const::TMS_Double_End) max_distance = TMS_Const::TMS_Double_gap;
            if ((std::abs(UTracks[itU].GetZ() - XTracks[itX].GetZ()) <= max_distance) || (itX > 0 && itU == 0 && itV == 0)) {
              UTracks[itU].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              VTracks[itV].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
              CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);

              // Handling cases of two hits in same plane
              if (itU > 0 && UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                UTracks[itU - 1].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV], XTracks[itX]));
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]); // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be aded with the next push_back statement
              }
              if (itV > 0 && VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                VTracks[itV - 1].SetRecoY(CompareY(UTracks[itU], VTracks[itV - 1], XTracks[itX]));
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);
                if (itV > 0) --itV;
              }
              if (itX > 0 && XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }

              (aTrack.Hits).push_back(UTracks[itU]);
              (aTrack.Hits).push_back(VTracks[itV]);
              (aTrack.Hits).push_back(XTracks[itX]);
#ifdef DEBUG
              std::cout << "same in all" << std::endl;
              std::cout << "Hit U: " << UTracks[itU].GetRecoX() << " | " << UTracks[itU].GetRecoY() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetRecoX() << " | " << VTracks[itV].GetRecoY() << " | " << VTracks[itV].GetZ() << " / X: " << XTracks[itX].GetRecoX() << " | " << XTracks[itX].GetNotZ() << " | " << XTracks[itX].GetZ() << std::endl;
#endif
              if (itX > 0) --itX;
            } else {// Gap in X handling
                if (UTracks[itU].GetZ() < XTracks[itX].GetZ()) {
                    UTracks[itU].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX - 1]));
                    VTracks[itV].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX - 1]));
                    CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
                    CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
#ifdef DEBUG
                    std::cout << "same in UV, gap in X" << std::endl;
                    std::cout << "Hit U: " << UTracks[itU].GetRecoX() << " | " << UTracks[itU].GetRecoY() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetRecoX() << " | " << VTracks[itV].GetRecoY() << " | " << VTracks[itV].GetZ() << " / X: " << XTracks[itX].GetNotZ() << " | " << XTracks[itX].GetZ() << std::endl;
#endif
                    if (itU > 0 && UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                        CalculateRecoY(UTracks[itU - 1], UTracks[itU - 1], VTracks[itV]);
                        CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                        (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                        if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
                    }
                    if (itV > 0 && VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                        CalculateRecoY(VTracks[itV - 1], UTracks[itU], VTracks[itV - 1]);
                        CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
                        (aTrack.Hits).push_back(VTracks[itV]);
                        if (itV > 0) --itV;
                    }

                    (aTrack.Hits).push_back(UTracks[itU]);
                    (aTrack.Hits).push_back(VTracks[itV]);
                } else if ((VTracks[itV]).GetZ() > (XTracks[itX]).GetZ()) {
                    CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], XTracks[itX]);
                    XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
#ifdef DEBUG
                    std::cout << "gap in UV, X" << std::endl;
                    std::cout << "Hit U: " << UTracks[itU].GetNotZ() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetNotZ() << " | " << VTracks[itV].GetZ() << " / X: " << XTracks[itX].GetNotZ() << " | " << XTracks[itX].GetT() << std::endl;
#endif
                    if (itX > 0 && XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                        XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                        CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], XTracks[itX - 1]);
                        (aTrack.Hits).push_back(XTracks[itX]);  // This adds the original hit
                        if (itX > 0) --itX; // and this allows for the other hit then to be added with the next push_back statement
                    }

                    (aTrack.Hits).push_back(XTracks[itX]);
                    if (itX > 0) --itX;
                }
            }
          }
          if (itU > 0) --itU;
          if (itV > 0) --itV;
        } else if ((UTracks[itU]).GetPlaneNumber() > (VTracks[itV]).GetPlaneNumber()) {
          if (stereo_view) {
#ifdef DEBUG
            std::cout << "Gap in U" << std::endl;
            std::cout << "Hit U: " << UTracks[itU].GetNotZ() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetNotZ() << " | " << VTracks[itV].GetZ() << std::endl;
#endif   
            if (itU > 0 && itV > 0) {
              CalculateRecoY(VTracks[itV], UTracks[itU - 1], VTracks[itV]);
              CalculateRecoX(UTracks[itU - 1], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                CalculateRecoY(VTracks[itV - 1], UTracks[itU - 1], VTracks[itV - 1]);
                CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allows for the other hit then to be added with the next push_back statement
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              if (itV > 0) --itV;
            } else if (itU == 0 && itV > 0) {
              CalculateRecoY(VTracks[itV], UTracks[itU], VTracks[itV]);
              CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                CalculateRecoY(VTracks[itV - 1], UTracks[itU], VTracks[itV - 1]);
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allows for the other hit then to be added with the next push_back statement
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              if (itV > 0) --itV;
            } else if (itU > 0 && itV == 0) {
              CalculateRecoY(UTracks[itU], UTracks[itU], VTracks[itV]);
              CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                CalculateRecoY(UTracks[itU - 1], UTracks[itU - 1], VTracks[itV]);
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }
              (aTrack.Hits).push_back(UTracks[itU]);
              if (itU > 0) --itU;
            }
          } else {
#ifdef DEBUG
            std::cout << "Gap in U, XTrack" << std::endl;
            std::cout << "Hit U: " << UTracks[itU].GetNotZ() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetNotZ() << " | " << VTracks[itV].GetZ() << " / X: " << XTracks[itX].GetNotZ() << " | " << XTracks[itX].GetZ() << std::endl;
#endif
            if (itU > 0 && itV > 0) {
              VTracks[itV].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU - 1], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU - 1], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                VTracks[itV - 1].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV - 1], XTracks[itX]));
                CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allows for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itV > 0) --itV;
            } else if (itU == 0 && itV > 0) {
              VTracks[itV].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                VTracks[itV - 1].SetRecoY(CompareY(UTracks[itU], VTracks[itV - 1], XTracks[itX]));
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);

                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allso for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itV > 0) --itV;
            } else if (itU > 0 && itV == 0) {
              UTracks[itU].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                UTracks[itU - 1].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV], XTracks[itX]));
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(UTracks[itU]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itU > 0) --itU;
            }
            if (itX > 0) --itX;
          }
        } else if ((UTracks[itU]).GetPlaneNumber() < (VTracks[itV]).GetPlaneNumber()) {
          if (stereo_view) {
#ifdef DEBUG
            std::cout << "Gap in V" << std::endl;
            std::cout << "Hit U: " << UTracks[itU].GetNotZ() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetNotZ() << " | " << VTracks[itV].GetZ() << std::endl;
#endif
            if (itV > 0 && itU > 0) {
              CalculateRecoY(UTracks[itU], UTracks[itU], VTracks[itV - 1]);
              CalculateRecoX(UTracks[itU], VTracks[itV - 1], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                CalculateRecoY(UTracks[itU - 1], UTracks[itU - 1], VTracks[itV - 1]);
                CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit                  if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }

              (aTrack.Hits).push_back(UTracks[itU]);
              if (itU > 0) --itU;
            } else if (itV == 0 && itU > 0) {
              CalculateRecoY(UTracks[itU], UTracks[itU], VTracks[itV]);
              CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                CalculateRecoY(UTracks[itU - 1], UTracks[itU - 1], VTracks[itV]);
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }
              (aTrack.Hits).push_back(UTracks[itU]);
              if (itU > 0) --itU;
            } else if (itV > 0 && itU == 0) {
              CalculateRecoY(VTracks[itV], UTracks[itU], VTracks[itV]);
              CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                CalculateRecoY(VTracks[itV - 1], UTracks[itU], VTracks[itV - 1]);
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allows for the other hit then to be added with the next push_back statement
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              if (itV > 0) --itV;
            }

          } else {
#ifdef DEBUG
            std::cout << "Gap in V, XTrack" << std::endl;
            std::cout << "Hit U: " << UTracks[itU].GetNotZ() << " | " << UTracks[itU].GetZ() << " / V: " << VTracks[itV].GetNotZ() << " | " << VTracks[itV].GetZ() << " / X: " << XTracks[itX].GetNotZ() << " | " << XTracks[itX].GetZ() << std::endl;
#endif
            if (itV > 0 && itU > 0) {
              UTracks[itU].SetRecoY(CompareY(UTracks[itU], VTracks[itV - 1], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV - 1], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV - 1], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                UTracks[itU - 1].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV - 1], XTracks[itX]));
                CalculateRecoX(UTracks[itU - 1], VTracks[itV - 1], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(UTracks[itU]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itU > 0) --itU;
            } else if (itV == 0 && itU > 0) {
              UTracks[itU].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV], UTracks[itU]);
              if (UTracks[itU].GetZ() == UTracks[itU - 1].GetZ()) {
                UTracks[itU - 1].SetRecoY(CompareY(UTracks[itU - 1], VTracks[itV], XTracks[itX]));
                CalculateRecoX(UTracks[itU - 1], VTracks[itV], UTracks[itU - 1]);
                (aTrack.Hits).push_back(UTracks[itU]);  // This adds the original hit
                if (itU > 0) --itU; // and this allows for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(UTracks[itU]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itU > 0) --itU;
            } else if (itV > 0 && itU == 0) {
              VTracks[itV].SetRecoY(CompareY(UTracks[itU], VTracks[itV], XTracks[itX]));//XTracks[itX].GetNotZ());
              XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
              CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX]);
              CalculateRecoX(UTracks[itU], VTracks[itV], VTracks[itV]);
              if (VTracks[itV].GetZ() == VTracks[itV - 1].GetZ()) {
                VTracks[itV - 1].SetRecoY(CompareY(UTracks[itU], VTracks[itV - 1], XTracks[itX]));
                CalculateRecoX(UTracks[itU], VTracks[itV - 1], VTracks[itV - 1]);
                (aTrack.Hits).push_back(VTracks[itV]);  // This adds the original hit
                if (itV > 0) --itV; // and this allows for the other hit then to be added with the next push_back statement
              }
              if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                XTracks[itX - 1].SetRecoY(XTracks[itX - 1].GetNotZ());
                CalculateRecoX(UTracks[itU], VTracks[itV], XTracks[itX - 1]);
                (aTrack.Hits).push_back(XTracks[itX]);
                if (itX > 0) --itX;
              }
              (aTrack.Hits).push_back(VTracks[itV]);
              (aTrack.Hits).push_back(XTracks[itX]);
              if (itV > 0) --itV;
            }
            if (itX > 0) --itX;
          }
        }
        sane = (itU > 0 || itV > 0);
        if (Xrun && Xback_match && Xfront_match) sane = (itU > 0 || itV > 0 || itX > 0);
      }
  
      // If same start and end, assign start and end hit in track
      if (UTracks.front().GetPlaneNumber() == VTracks.front().GetPlaneNumber()) {
        bool stereo_view = true;
        if (Xrun && Xback_match && Xfront_match) {
          //if (UTracks.front().GetPlaneNumber() < TMS_Const::TMS_nThinPlanes / 2 && std::abs(UTracks.front().GetZ() - XTracks.front().GetZ() < TMS_Const::TMS_Thin_gap * 4 +10)) stereo_view = false;
          //else if (UTracks.front().GetPlaneNumber() >= TMS_Const::TMS_nThickPlanes / 2 && std::abs(UTracks.front().GetZ() - XTracks.front().GetZ() < TMS_Const::TMS_Thick_gap * 4 + 10)) stereo_view = false;
          //if (UTracks.front().GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 + TMS_Const::TMS_nThickPlanes / 2 && std::abs(UTracks.front().GetZ() - XTracks.front().GetZ() < TMS_Const::TMS_Double_gap * 4 + 10)) stereo_view = false;
          if ((UTracks.front().GetZ() <= TMS_Const::TMS_Thick_Start && XTracks.front().GetZ() <= TMS_Const::TMS_Thick_Start) 
              && (std::abs(UTracks.front().GetZ() - XTracks.front().GetZ()) < TMS_Const::TMS_Thin_gap * 4 + 10)) 
            stereo_view = false;
          else if ((UTracks.front().GetZ() <= TMS_Const::TMS_Double_Start && XTracks.front().GetZ() <= TMS_Const::TMS_Double_Start) 
              && (std::abs(UTracks.front().GetZ() - XTracks.front().GetZ()) < TMS_Const::TMS_Thick_gap * 4 + 10)) 
            stereo_view = false;
          else if ((UTracks.front().GetZ() <= TMS_Const::TMS_Double_End && XTracks.front().GetZ() <= TMS_Const::TMS_Double_End) 
              && (std::abs(UTracks.front().GetZ() - XTracks.front().GetZ()) < TMS_Const::TMS_Double_gap * 4 + 10))
            stereo_view = false;
        }
        if (stereo_view) {
          CalculateRecoY(UTracks.front(), UTracks.front(), VTracks.front());
          CalculateRecoX(UTracks.front(), VTracks.front(), UTracks.front());
          aTrack.End[0] = UTracks.front().GetRecoX();
          aTrack.End[1] = UTracks.front().GetRecoY();
          aTrack.End[2] = UTracks.front().GetZ();
          aTrack.End[3] = UTracks.front().GetT();
#ifdef DEBUG              
          std::cout << "End equal assigned" << std::endl;
#endif              
        } else {
          VTracks.front().SetRecoX(0.5 * (UTracks.front().GetNotZ() + VTracks.front().GetNotZ()));
          VTracks.front().SetRecoY(CompareY(UTracks.front(), VTracks.front(), XTracks.front()));
          aTrack.End[0] = 0.5 * (VTracks.front().GetNotZ() + UTracks.front().GetNotZ());
          aTrack.End[1] = CompareY(UTracks.front(), VTracks.front(), XTracks.front());//XTracks.front().GetNotZ();
          aTrack.End[2] = UTracks.front().GetZ();
          aTrack.End[3] = UTracks.front().GetT();
#ifdef DEBUG
          std::cout << "End equal assigned, XTrack" << std::endl;
#endif
        }
      }
      if (UTracks.back().GetPlaneNumber() == VTracks.back().GetPlaneNumber()) {
        bool stereo_view = true;
        if (Xrun && Xback_match && Xfront_match) {
          //if (UTracks.back().GetPlaneNumber() < TMS_Const::TMS_nThinPlanes / 2 && std::abs(UTracks.back().GetZ() - XTracks.back().GetZ() < TMS_Const::TMS_Thin_gap * 4 + 10)) stereo_view = false;
          //else if (UTracks.back().GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 && std::abs(UTracks.back().GetZ() - XTracks.back().GetZ() < TMS_Const::TMS_Thick_gap * 4 + 10)) stereo_view = false;
          //if (UTracks.back().GetPlaneNumber() >= TMS_Const::TMS_nThinPlanes / 2 + TMS_Const::TMS_nThickPlanes / 2 && std::abs(UTracks.back().GetZ() < TMS_Const::TMS_Double_gap * 4 + 10)) stereo_view = false;
          if ((UTracks.back().GetZ() <= TMS_Const::TMS_Thick_Start && XTracks.back().GetZ() <= TMS_Const::TMS_Thick_Start) 
              && (std::abs(UTracks.back().GetZ() - XTracks.back().GetZ()) < TMS_Const::TMS_Thin_gap * 4 + 10))
            stereo_view = false;
          else if ((UTracks.back().GetZ() <= TMS_Const::TMS_Double_Start && XTracks.back().GetZ() <= TMS_Const::TMS_Double_Start) 
              && (std::abs(UTracks.back().GetZ() - XTracks.back().GetZ()) < TMS_Const::TMS_Thick_gap * 4 + 10)) 
            stereo_view = false;
          else if ((UTracks.back().GetZ() <= TMS_Const::TMS_Double_End && XTracks.back().GetZ() <= TMS_Const::TMS_Double_End)
              && (std::abs(UTracks.back().GetZ() - XTracks.back().GetZ()) < TMS_Const::TMS_Double_gap * 4 + 10))
            stereo_view = false;
        }
        if (stereo_view) {
          CalculateRecoY(VTracks.back(), UTracks.back(), VTracks.back());
          CalculateRecoX(UTracks.back(), VTracks.back(), VTracks.back());
          aTrack.Start[0] = VTracks.back().GetRecoX();
          aTrack.Start[1] = VTracks.back().GetRecoY();
          aTrack.Start[2] = VTracks.back().GetZ();
          aTrack.Start[3] = VTracks.back().GetT();
#ifdef DEBUG              
          std::cout << "Start equal assigned" << std::endl;
#endif              
        } else {
          XTracks.back().SetRecoX(0.5 * (VTracks.back().GetRecoX() + UTracks.back().GetNotZ()));
          XTracks.back().SetRecoY(CompareY(UTracks.back(), VTracks.back(), XTracks.back()));
          aTrack.Start[0] = 0.5 * (VTracks.back().GetRecoX() + UTracks.back().GetNotZ());
          aTrack.Start[1] = CompareY(UTracks.back(), VTracks.back(), XTracks.back());//XTracks.back().GetNotZ();
          aTrack.Start[2] = XTracks.back().GetZ();
          aTrack.Start[3] = XTracks.back().GetT();
#ifdef DEBUG
          std::cout << "Start equal assigned, XTrack" << std::endl;
#endif                
        }
      }
      // Sort track
      SpatialPrio(aTrack.Hits);  
      if (TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UseTrackSmoothing()) {
        aTrack.ApplyTrackSmoothing();
        aTrack.Start[0] = aTrack.Hits.front().GetRecoX();
        aTrack.Start[1] = aTrack.Hits.front().GetRecoY();
        aTrack.Start[2] = aTrack.Hits.front().GetZ();
        aTrack.Start[3] = aTrack.Hits.front().GetT();
        aTrack.End[0] = aTrack.Hits.back().GetRecoX();
        aTrack.End[1] = aTrack.Hits.back().GetRecoY();
        aTrack.End[2] = aTrack.Hits.back().GetZ();
        aTrack.End[3] = aTrack.Hits.back().GetT();
      }

      // Smoothing of start and end of track in case of too much 'flailing around' in the y direction
      // end
      /*bool SameSign = true;
      if ((aTrack.End[1] > 0 && aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY() < 0) || (aTrack.End[1] < 0 && aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY() > 0)) SameSign = false;
      if ((SameSign &&std::abs(aTrack.End[1] - aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) >= 676.6) || (!SameSign && std::abs(aTrack.End[1]) + std::abs(aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) >= 674.6)) {
        aTrack.End[1] = (aTrack.End[1] + aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) / 2;
        if (aTrack.End[1] > 244.0) aTrack.End[1] = 244.0;
        else if (aTrack.End[1] < -2949.0) aTrack.End[1] = -2949.0;
        aTrack.Hits[aTrack.Hits.size() - 1].SetRecoY(aTrack.End[1]);
      }
      // start
      SameSign = true;
      if ((aTrack.Start[1] > 0 && aTrack.Hits[2].GetRecoY() < 0) || (aTrack.Start[1] < 0 && aTrack.Hits[2].GetRecoY() > 0)) SameSign = false;
      if ((SameSign && std::abs(aTrack.Start[1] - aTrack.Hits[2].GetRecoY()) >= 676.6) || (!SameSign && std::abs(aTrack.Start[1]) + std::abs(aTrack.Hits[2].GetRecoY()) >= 674.6)) {
        aTrack.Start[1] = (aTrack.Start[1] + aTrack.Hits[2].GetRecoY()) / 2;
        if (aTrack.Start[1] > 244.0) aTrack.Start[1] = 244.0;
        else if (aTrack.Start[1] < -2949.0) aTrack.Start[1] = -2949.0;
        aTrack.Hits[0].SetRecoY(aTrack.Start[1]);
      }*/
#ifdef DEBUG
      for (auto hits: aTrack.Hits) {
        std::cout << "Match: " << hits.GetRecoX() << "," << hits.GetRecoY() << "," << hits.GetZ() << "(" << hits.GetBar().GetBarType() << ")," << hits.GetT() << std::endl;
      }
#endif
      // Charge ID
      aTrack.Charge = ChargeID.ID_Track_Charge(aTrack.Hits);

      // Track Length
      aTrack.Length = CalculateTrackLength3D(aTrack);
#ifdef DEBUG
      std::cout << "Added TrackLength: " << aTrack.Length << std::endl;
#endif          
      // Track Energy
      aTrack.EnergyDeposit = CalculateTrackEnergy3D(aTrack);
#ifdef DEBUG
      std::cout << "Added TrackEnergyDeposit: " << aTrack.EnergyDeposit << std::endl;
#endif
      // Track Direction
      if (TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance() >= aTrack.Hits.size()) {
        double direction_x = aTrack.Start[0] - aTrack.End[0];
        double direction_y = aTrack.Start[1] - aTrack.End[1];
        double direction_z = aTrack.Start[2] - aTrack.End[2];
        aTrack.SetStartDirection(direction_x, direction_y, direction_z);
      } else {
        double direction_x = aTrack.Start[0] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetRecoX();
        double direction_y = aTrack.Start[1] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetRecoY();
        double direction_z = aTrack.Start[2] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetZ();
        aTrack.SetStartDirection(direction_x, direction_y, direction_z);
      }
//#ifdef DEBUG          
      std::cout << "Start: " << aTrack.Start[0] << "mm(x) | " << aTrack.Start[1] << "mm(y) | " << aTrack.Start[2] << "mm(z) | " << aTrack.Start[3] << "s" << std::endl;
      std::cout << "End: " << aTrack.End[0] << "mm(x) | " << aTrack.End[1] << "mm(y) | " << aTrack.End[2] << "mm(z) | " << aTrack.End[3] << "s" << std::endl;
      std::cout << "Added Direction: " << aTrack.StartDirection[0] << " | " << aTrack.StartDirection[1] << " | " << aTrack.StartDirection[2] << std::endl;
//#endif       
  
      returned.push_back(aTrack);

      if (HoughCandidatesU.size() == 1) break;
      if (HoughCandidatesV.size() == 1) break;

      if (Xrun && !XTracks.empty()) {
        if (Xback_match && Xfront_match) {
          // If match was made, remove the candidate (simple) track from candidate list
          SortedHoughCandidatesV.erase(Vhelper);
          SortedHoughCandidatesX.erase(helper);
          Vhelper = SortedHoughCandidatesV.begin();
          helper = SortedHoughCandidatesX.begin();
          // Set iterator U tracks to next track
          ++Uhelper;
        } else {
          SortedHoughCandidatesV.erase(Vhelper);
          Vhelper = SortedHoughCandidatesV.begin();
          ++Uhelper;
        }
      } else {
        // If match was made, remove the candidate (simple) track from candidate list
        SortedHoughCandidatesV.erase(Vhelper);
        if (SortedHoughCandidatesV.size() > 1) Vhelper = SortedHoughCandidatesV.begin();
        // Set iterator for U tracks to next track
        ++Uhelper;
      }
    } else {
      if (Xrun && !XTracks.empty()) {
        if (helper == SortedHoughCandidatesX.end()-1) {
          helper = SortedHoughCandidatesX.begin();
          if (Vhelper == SortedHoughCandidatesV.end()-1) {
            ++Uhelper;
            Vhelper = SortedHoughCandidatesV.begin();
          } else {
            ++Vhelper;
          }
        } else ++helper;
      } else {
        if (Vhelper == SortedHoughCandidatesV.end()-1) {
          ++Uhelper;
          Vhelper = SortedHoughCandidatesV.begin();
        } else ++Vhelper;
        if (SortedHoughCandidatesV.size() < 2) break;
      }
    }
  }
  return returned;
}

std::vector<TMS_Track> TMS_TrackFinder::TrackMatching3D_XY() {
#ifdef DEBUG
  std::cout << "3D matching" << std::endl;
  std::cout << "size Candidates: U: " << HoughCandidatesU.size() << " | V: " << HoughCandidatesV.size() << std::endl;//" | X: " << HoughCandidatesX.size() << std::endl;
#endif

  // Sorting the candidate tracks descending by their hit numbers
  if (HoughCandidatesY.size() > 1) {
    std::sort(HoughCandidatesY.begin(), HoughCandidatesY.end(), SortByHitNumber);
    SortedHoughCandidatesY = HoughCandidatesY;
  } else SortedHoughCandidatesY = HoughCandidatesY;
  if (HoughCandidatesX.size() > 1) {
    std::sort(HoughCandidatesX.begin(), HoughCandidatesX.end(), SortByHitNumber);
    SortedHoughCandidatesX = HoughCandidatesX;
  } else SortedHoughCandidatesX = HoughCandidatesX;

  std::vector<TMS_Track> returned;

  bool TimeSlicing = TMS_Manager::GetInstance().Get_Reco_TIME_RunTimeSlicer();

  // 3D matching of tracks
  std::vector<std::vector<TMS_Hit> >::iterator Yhelper = SortedHoughCandidatesY.begin();
  std::vector<std::vector<TMS_Hit> >::iterator Xhelper = SortedHoughCandidatesX.begin();

  while (Xhelper != SortedHoughCandidatesX.end()) {
    if (SortedHoughCandidatesY.empty()) {
      std::cout << "Not enough Y tracks" << std::endl;
      break;
    }
    std::vector<TMS_Hit> YTracks = *Yhelper;
    std::vector<TMS_Hit> XTracks = *Xhelper;
        
    SpatialPrio(XTracks);
    SpatialPrio(YTracks);
#ifdef DEBUG
        std::cout << "XTrack back: " << XTracks.front().GetPlaneNumber() << " | " << XTracks.front().GetBarNumber() << " | " << XTracks.front().GetT() << " front: " << XTracks.back().GetPlaneNumber() << " | " << XTracks.back().GetBarNumber() << " | " << XTracks.back().GetT() << std::endl;
        std::cout << "YTrack back: " << YTracks.front().GetPlaneNumber() << " | " << YTracks.front().GetBarNumber() << " | " << YTracks.front().GetT() << " front: " << YTracks.back().GetPlaneNumber() << " | " << YTracks.back().GetBarNumber() << " | " << YTracks.back().GetT() << std::endl;
#endif
 
        // Conditions for close enough tracks: within +/-3 plane numbers, +/-12 bar numbers and in same time slice within 30ns (Asa note: changed to 15 ns for now and just one of time and plane number condition needs to be met)
        bool back_match = false;
        bool front_match = false;
        int strike = 0;
        // Check for the front and back whether both conditions (plane limit and time limit) are met or not
        if (TimeSlicing) {
          // front of simple tracks: set strike to 1
          if (std::abs(XTracks.front().GetPlaneNumber() - YTracks.front().GetPlaneNumber()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit()
              || std::abs(XTracks.front().GetT() - YTracks.front().GetT()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit()) strike = 1;
          // end of simple tracks: add 2 to strike
          if (std::abs(XTracks.back().GetPlaneNumber() - YTracks.back().GetPlaneNumber()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit()
              || std::abs(XTracks.back().GetT() - YTracks.back().GetT()) > TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit()) strike += 2;
        }  
        if (TimeSlicing) {
            //unlike UV, it is hard to compare bar number directly in XY, ignore it.
//          bool bar_front = (std::abs(XTracks.front().GetBarNumber() - YTracks.front().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
//          bool bar_back = (std::abs(XTracks.back().GetBarNumber() - YTracks.back().GetBarNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_BarLimit());
          bool slice_front = (XTracks.front().GetSlice() == YTracks.front().GetSlice());
          bool slice_back = (XTracks.back().GetSlice() == YTracks.back().GetSlice());
          //ignore time matching for XY for now.
//          bool time_front = (std::abs(XTracks.front().GetT() - YTracks.front().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
//          bool time_back = (std::abs(XTracks.back().GetT() - YTracks.back().GetT()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TimeLimit());
          bool time_front = true;
          bool time_back = true;
          bool plane_front = (std::abs(XTracks.front().GetPlaneNumber() - YTracks.front().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
          bool plane_back = (std::abs(XTracks.back().GetPlaneNumber() - YTracks.back().GetPlaneNumber()) <= TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
          if (strike == 0) {
            // front and end match
            back_match =  (slice_back && time_back && plane_back);
            front_match = (slice_front && time_front && plane_front);
          } else if (strike == 1) {
            // front needs exemption, end matches
            back_match =  (slice_back && time_back && plane_back);
            front_match = (slice_front && (time_front || plane_front));
          } else if (strike == 2) {
            // front matches, end needs exemption
            back_match =  (slice_back && (time_back || plane_back));
            front_match = (slice_front && time_front && plane_front);
          }
        }
        else {
          back_match = (std::abs(XTracks.front().GetPlaneNumber() - YTracks.front().GetPlaneNumber()) < TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
          front_match = (std::abs(XTracks.back().GetPlaneNumber() - YTracks.back().GetPlaneNumber()) < TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_PlaneLimit());
        } 

        if (back_match&&front_match) { 
            TMS_Track aTrack;

            // Make sure that the hits are in the correct order
            if (XTracks.back().GetZ() > XTracks.front().GetZ()) std::reverse(XTracks.begin(), XTracks.end());
            if (YTracks.back().GetZ() > YTracks.front().GetZ()) std::reverse(YTracks.begin(), YTracks.end());
//            std::cout << "XTrack FRONT: " << XTracks.back().GetPlaneNumber() << " BACK: " << XTracks.front().GetPlaneNumber() << std::endl;
//            std::cout << "YTrack FRONT: " << YTracks.back().GetPlaneNumber() << " BACK: " << YTracks.front().GetPlaneNumber() << std::endl;
//            std::cout<<"back_Match:"<<back_match<<"front_match:"<<front_match<<std::endl;
#ifdef DEBUG          
            std::cout << "XTrack FRONT: " << XTracks.back().GetPlaneNumber() << " BACK: " << XTracks.front().GetPlaneNumber() << std::endl;
            std::cout << "YTrack FRONT: " << YTracks.back().GetPlaneNumber() << " BACK: " << YTracks.front().GetPlaneNumber() << std::endl;
#endif          
            // If same plane number for start but different for end (U/V)
            if (XTracks.front().GetPlaneNumber() != YTracks.front().GetPlaneNumber()) {
              // If XTrack ends after YTrack
              if (XTracks.front().GetPlaneNumber() > YTracks.front().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                  XTracks.front().SetRecoX(YTracks.front().GetNotZ());
                  XTracks.front().SetRecoY(XTracks.front().GetNotZ());
                  aTrack.End[0] = YTracks.front().GetNotZ();
                  aTrack.End[1] = XTracks.front().GetNotZ();
                  aTrack.End[2] = XTracks.front().GetZ();
                  aTrack.End[3] = XTracks.front().GetT();
                  (aTrack.Hits).push_back(XTracks.front());
                }

              } // If UTrack ends before VTrack
                else if (XTracks.front().GetPlaneNumber() < YTracks.front().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                  YTracks.front().SetRecoX(YTracks.front().GetNotZ());
                  YTracks.front().SetRecoY(XTracks.front().GetNotZ());
                  aTrack.End[0] = YTracks.front().GetNotZ();
                  aTrack.End[1] = XTracks.front().GetNotZ();
                  aTrack.End[2] = YTracks.front().GetZ();
                  aTrack.End[3] = YTracks.front().GetT();
                  (aTrack.Hits).push_back(YTracks.front());
                } 
              }
            } 
            // If different plane number for start but same for end (U/V)
            if (XTracks.back().GetPlaneNumber() != YTracks.back().GetPlaneNumber()) {
              // If XTrack starts after YTrack  
              if (XTracks.back().GetPlaneNumber() > YTracks.back().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                  aTrack.Start[0] = YTracks.back().GetNotZ();
                  aTrack.Start[1] = XTracks.back().GetNotZ();
                  aTrack.Start[2] = YTracks.back().GetZ();
                  aTrack.Start[3] = YTracks.back().GetT();
#ifdef DEBUG
                  std::cout << "YTrack starts before XTrack, no XTrack" << std::endl;
#endif                
                  YTracks.back().SetRecoX(YTracks.back().GetNotZ());
                  YTracks.back().SetRecoY(XTracks.back().GetNotZ());
                  (aTrack.Hits).push_back(YTracks.back());
                } 
              } 
                else if (XTracks.back().GetPlaneNumber() < YTracks.back().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                  aTrack.Start[0] = YTracks.back().GetNotZ();
                  aTrack.Start[1] = XTracks.back().GetNotZ();
                  aTrack.Start[2] = XTracks.back().GetZ();
                  aTrack.Start[3] = XTracks.back().GetT();
#ifdef DEBUG                
                  std::cout << "XTrack starts before YTrack, no XTrack" << std::endl;
#endif                
                  XTracks.back().SetRecoX(YTracks.back().GetNotZ());
                  XTracks.back().SetRecoY(XTracks.back().GetNotZ());
                  (aTrack.Hits).push_back(XTracks.back());
                }
              } 
            }
            // If same start and end, assign start and end hit in track
            if (XTracks.front().GetPlaneNumber() == YTracks.front().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                    XTracks.front().SetRecoX(YTracks.front().GetNotZ());
                    XTracks.front().SetRecoY(XTracks.front().GetNotZ());
                    aTrack.End[0] = YTracks.front().GetNotZ();
                    aTrack.End[1] = XTracks.front().GetNotZ();
                    aTrack.End[2] = YTracks.front().GetZ();
                    aTrack.End[3] = YTracks.front().GetT();
                    (aTrack.Hits).push_back(XTracks.front());
                }
            }
            if (XTracks.back().GetPlaneNumber() == YTracks.back().GetPlaneNumber()) {
                bool stereo_view = true;
                if (stereo_view) {
                    XTracks.back().SetRecoX(YTracks.back().GetNotZ());
                    XTracks.back().SetRecoY(XTracks.back().GetNotZ());
                    aTrack.Start[0] = YTracks.back().GetNotZ();
                    aTrack.Start[1] = XTracks.back().GetNotZ();
                    aTrack.Start[2] = XTracks.back().GetZ();
                    aTrack.Start[3] = XTracks.back().GetT();
                    (aTrack.Hits).push_back(XTracks.back());
                } 
            }

            // Add hits to track
            int itY = YTracks.size() - 1;
            int itX = XTracks.size() - 1;
//            std::cout<<"itY:"<<itY<<", itX:"<<itX<<std::endl;

            bool sane = (itY > 0 || itX > 0);

            while (sane) {  // Track seems to be backwards, so run adding of hits backwards
//                std::cout <<"itY:"<<itY<<" ,X: " << YTracks[itY].GetNotZ() <<  ", Z: " << YTracks[itY].GetZ() <<", PlaneNumber: " << YTracks[itY].GetPlaneNumber()<< std::endl;
//                std::cout <<"itX:"<<itX<<" ,Y: " << XTracks[itX].GetNotZ() <<  ", Z: " << XTracks[itX].GetZ() <<", PlaneNumber: " << XTracks[itX].GetPlaneNumber()<< std::endl;
                // Sanity check for hits being in the detector volume (x and z)
                bool hit_outside = false;
                if (XTracks[itX].GetNotZ() > 500.0 || XTracks[itX].GetNotZ() == 0. || XTracks[itX].GetNotZ() < -4000.0
                        || XTracks[itX].GetZ() < 11000 || XTracks[itX].GetZ() > 20000) {
                    --itX;
                    hit_outside = true;
                }
                if (std::abs(YTracks[itY].GetNotZ()) > 4000.0 or YTracks[itY].GetNotZ() == 0.
                        || YTracks[itY].GetZ() < 11000 || YTracks[itY].GetZ() > 20000) {
                    --itY;
                    hit_outside = true;
                }
                if (hit_outside) continue;

                // Stereo check
                bool stereo_view = true;
                // Check if a neighbouring hit is from a X layer
                if ((XTracks[itX]).GetPlaneNumber() == (YTracks[itY]).GetPlaneNumber()) {
                    if (stereo_view) {
                        XTracks[itX].SetRecoX(YTracks[itY].GetNotZ());
                        XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
                        YTracks[itY].SetRecoX(YTracks[itY].GetNotZ());
                        YTracks[itY].SetRecoY(XTracks[itX].GetNotZ());
//                        std::cout << "same" << std::endl;
//                        std::cout << "YTracks X:: " << YTracks[itY].GetRecoX() << " Y: " << YTracks[itY].GetRecoY() << " Z: " << YTracks[itY].GetZ() << " than: " << XTracks[itX].GetRecoX() << " | " << XTracks[itX].GetRecoY() << " | " << XTracks[itX].GetZ() << std::endl;
#ifdef DEBUG
                        std::cout << "same" << std::endl;
                        std::cout << "Hit: " << YTracks[itY].GetRecoX() << " | " << YTracks[itY].GetRecoY() << " | " << YTracks[itY].GetZ() << " than: " << XTracks[itX].GetRecoX() << " | " << XTracks[itX].GetRecoY() << " | " << XTracks[itX].GetZ() << std::endl;
#endif
                        // Handling cases of two hits in same plane to be matched
                        // By adding a loop into these statements one could also take care of more than 2 hits in the same plane/with the same z coordinate

                        
                        if (itX > 0 && XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                            XTracks[itX-1].SetRecoX(YTracks[itY].GetNotZ());
                            XTracks[itX-1].SetRecoY(XTracks[itX-1].GetNotZ());
                            (aTrack.Hits).push_back(XTracks[itX]);  // This adds the original hit
                            if (itX > 0) --itX; // and this allows for the other hit then to be added with the next push_back statement
                        }
                        if (itY > 0 && YTracks[itY].GetZ() == YTracks[itY - 1].GetZ()) {
                            YTracks[itY-1].SetRecoX(YTracks[itY-1].GetNotZ());
                            YTracks[itY-1].SetRecoY(XTracks[itX].GetNotZ());
                            (aTrack.Hits).push_back(YTracks[itY]);
                            if (itY > 0) --itY;
                        }
                        (aTrack.Hits).push_back(XTracks[itX]);
                        (aTrack.Hits).push_back(YTracks[itY]);
                    } 
                    if (itX > 0) --itX;
                    if (itY > 0) --itY;
                }
                //just change the following index from original code
                else if ((XTracks[itX]).GetPlaneNumber() > (YTracks[itY]).GetPlaneNumber()) {
                        if (itX > 0 && itY > 0) {
                            YTracks[itY].SetRecoX(YTracks[itY].GetNotZ());
                            YTracks[itY].SetRecoY(XTracks[itX-1].GetNotZ());
                            if (YTracks[itY].GetZ() == YTracks[itY - 1].GetZ()) {
                                YTracks[itY-1].SetRecoX(YTracks[itY-1].GetNotZ());
                                YTracks[itY-1].SetRecoY(XTracks[itX-1].GetNotZ());
                                (aTrack.Hits).push_back(YTracks[itY]);  // This adds the original hit
                                if (itY > 0) --itY; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(YTracks[itY]);
                            if (itY > 0) --itY;
                        } else if (itX == 0 && itY > 0) {
                            YTracks[itY].SetRecoX(YTracks[itY].GetNotZ());
                            YTracks[itY].SetRecoY(XTracks[itX].GetNotZ());
                            if (YTracks[itY].GetZ() == YTracks[itY - 1].GetZ()) {
                                YTracks[itY-1].SetRecoX(YTracks[itY-1].GetNotZ());
                                YTracks[itY-1].SetRecoY(XTracks[itX].GetNotZ());
                                (aTrack.Hits).push_back(YTracks[itY]);  // This adds the original hit
                                if (itY > 0) --itY; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(YTracks[itY]);
                            if (itY > 0) --itY;
                        } else if (itX > 0 && itY == 0) {
                            XTracks[itX].SetRecoX(YTracks[itY].GetNotZ());
                            XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
                            if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                                XTracks[itX-1].SetRecoX(YTracks[itY].GetNotZ());
                                XTracks[itX-1].SetRecoY(XTracks[itX-1].GetNotZ());
                                (aTrack.Hits).push_back(XTracks[itX]);  // This adds the original hit
                                if (itX > 0) --itX; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(XTracks[itX]);
                            if (itX > 0) --itX;
                        }
                } else if ((XTracks[itX]).GetPlaneNumber() < (YTracks[itY]).GetPlaneNumber()) {
                    if (stereo_view) {
                        if (itY > 0 && itX > 0) {
                            XTracks[itX].SetRecoX(YTracks[itY-1].GetNotZ());
                            XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
                            if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                                XTracks[itX-1].SetRecoX(YTracks[itY-1].GetNotZ());
                                XTracks[itX-1].SetRecoY(XTracks[itX-1].GetNotZ());
                                (aTrack.Hits).push_back(XTracks[itX]);  // This adds the original hit
                                if (itX > 0) --itX; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(XTracks[itX]);
                            if (itX > 0) --itX;
                        } else if (itY == 0 && itX > 0) {
                            XTracks[itX].SetRecoX(YTracks[itY].GetNotZ());
                            XTracks[itX].SetRecoY(XTracks[itX].GetNotZ());
                            if (XTracks[itX].GetZ() == XTracks[itX - 1].GetZ()) {
                                XTracks[itX-1].SetRecoX(YTracks[itY].GetNotZ());
                                XTracks[itX-1].SetRecoY(XTracks[itX-1].GetNotZ());
                                (aTrack.Hits).push_back(XTracks[itX]);  // This adds the original hit
                                if (itX > 0) --itX; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(XTracks[itX]);
                            if (itX > 0) --itX;
                        } else if (itY > 0 && itX == 0) {
                            YTracks[itY].SetRecoX(YTracks[itY].GetNotZ());
                            YTracks[itY].SetRecoY(XTracks[itX].GetNotZ());
                            if (YTracks[itY].GetZ() == YTracks[itY - 1].GetZ()) {
                                YTracks[itY-1].SetRecoX(YTracks[itY-1].GetNotZ());
                                YTracks[itY-1].SetRecoY(XTracks[itX].GetNotZ());
                                (aTrack.Hits).push_back(YTracks[itY]);  // This adds the original hit
                                if (itY > 0) --itY; // and this allows for the other hit then to be added with the next push_back statement
                            }
                            (aTrack.Hits).push_back(YTracks[itY]);
                            if (itY > 0) --itY;
                        }
                    }
                }
                sane = (itX > 0 || itY > 0);
            }//sane

            // Sort track
            SpatialPrio(aTrack.Hits);


            ///at least no stmoothing for XY
//            if (TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UseTrackSmoothing()) aTrack.ApplyTrackSmoothing();

            // Smoothing of start and end of track in case of too much 'flailing around' in the y direction
            // end
            //            bool SameSign = true;
            //            if ((aTrack.End[1] > 0 && aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY() < 0) || (aTrack.End[1] < 0 && aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY() > 0)) SameSign = false;
            //            if ((SameSign &&std::abs(aTrack.End[1] - aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) >= 676.6) || (!SameSign && std::abs(aTrack.End[1]) + std::abs(aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) >= 674.6)) {
            //              aTrack.End[1] = (aTrack.End[1] + aTrack.Hits[aTrack.Hits.size() - 3].GetRecoY()) / 2;
            //              if (aTrack.End[1] > 244.0) aTrack.End[1] = 244.0;
            //              else if (aTrack.End[1] < -2949.0) aTrack.End[1] = -2949.0;
            //              aTrack.Hits[aTrack.Hits.size() - 1].SetRecoY(aTrack.End[1]);
            //            }
            //            // start
            //             // std::cout<<"3"<<std::endl;
            //            SameSign = true;
            //            if ((aTrack.Start[1] > 0 && aTrack.Hits[2].GetRecoY() < 0) || (aTrack.Start[1] < 0 && aTrack.Hits[2].GetRecoY() > 0)) SameSign = false;
            //            if ((SameSign && std::abs(aTrack.Start[1] - aTrack.Hits[2].GetRecoY()) >= 676.6) || (!SameSign && std::abs(aTrack.Start[1]) + std::abs(aTrack.Hits[2].GetRecoY()) >= 674.6)) {
            //              aTrack.Start[1] = (aTrack.Start[1] + aTrack.Hits[2].GetRecoY()) / 2;
            //              if (aTrack.Start[1] > 244.0) aTrack.Start[1] = 244.0;
            //              else if (aTrack.Start[1] < -2949.0) aTrack.Start[1] = -2949.0;
            //              aTrack.Hits[0].SetRecoY(aTrack.Start[1]);
            //            }
//            for (auto hits: aTrack.Hits) {
//                std::cout << "Match: " << hits.GetRecoX() << "," << hits.GetRecoY() << "," << hits.GetZ() << std::endl;
//            }
#ifdef DEBUG
            for (auto hits: aTrack.Hits) {
                std::cout << "Match: " << hits.GetRecoX() << "," << hits.GetRecoY() << "," << hits.GetZ() << "," << hits.GetT() << std::endl;
            }
#endif
            // Charge ID
            aTrack.Charge = ChargeID.ID_Track_Charge(aTrack.Hits);

            // Track Length
            aTrack.Length = CalculateTrackLength3D(aTrack);
#ifdef DEBUG
            std::cout << "Added TrackLength: " << aTrack.Length << std::endl;
#endif          
            // Track Energy
            aTrack.EnergyDeposit = CalculateTrackEnergy3D(aTrack);
#ifdef DEBUG
            std::cout << "Added TrackEnergyDeposit: " << aTrack.EnergyDeposit << std::endl;
#endif
            // Track Direction
            if (TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance() >= aTrack.Hits.size()) {
              double direction_x = aTrack.Start[0] - aTrack.End[0];
              double direction_y = aTrack.Start[1] - aTrack.End[1];
              double direction_z = aTrack.Start[2] - aTrack.End[2];
              aTrack.SetStartDirection(direction_x, direction_y, direction_z);
            } else {
              double direction_x = aTrack.Start[0] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetRecoX();
              double direction_y = aTrack.Start[1] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetRecoY();
              double direction_z = aTrack.Start[2] - aTrack.Hits[TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_DirectionDistance()].GetZ();
              aTrack.SetStartDirection(direction_x, direction_y, direction_z);
            }
#ifdef DEBUG          
            std::cout << "Start: " << aTrack.Start[0] << "mm(x) | " << aTrack.Start[1] << "mm(y) | " << aTrack.Start[2] << "mm(z) | " << aTrack.Start[3] << "s" << std::endl;
            std::cout << "End: " << aTrack.End[0] << "mm(x) | " << aTrack.End[1] << "mm(y) | " << aTrack.End[2] << "mm(z) | " << aTrack.End[3] << "s" << std::endl;
            //std::cout << "Added Direction: " << aTrack.Direction[0] << " | " << aTrack.Direction[1] << " | " << aTrack.Direction[2] << std::endl;
            std::cout << "Added StartDirection: " << aTrack.StartDirection[0] << " | " << aTrack.StartDirection[1] << " | " << aTrack.StartDirection[2] << std::endl;
#endif          

            returned.push_back(aTrack);

            if (HoughCandidatesX.size() == 1) break;
            if (HoughCandidatesY.size() == 1) break;

            // If match was made, remove the candidate (simple) track from candidate list
            SortedHoughCandidatesY.erase(Yhelper);
            if (SortedHoughCandidatesY.size() > 1) Yhelper = SortedHoughCandidatesY.begin();
            // Set iterator for X tracks to next track
            ++Xhelper;
        } else {
            if (Yhelper == SortedHoughCandidatesY.end()-1) {
                ++Xhelper;
                Yhelper = SortedHoughCandidatesY.begin();
            } else ++Yhelper;
            if (SortedHoughCandidatesY.size() < 2) break;
        }   
  }
  return returned;
}

void TMS_TrackFinder::CalculateRecoY(TMS_Hit &OneHit, TMS_Hit &UHit, TMS_Hit &VHit) {
//  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();
//  Double_t localOne[3] = {OneHit.GetNotZ(), 0, OneHit.GetZ()};
//  Double_t localOther[3] = {OtherHit.GetNotZ(), 0, OtherHit.GetZ()};
//  Double_t positionOne[3];
//  Double_t positionOther[3];
//  geom->GetCurrentMatrix()->LocalToMaster(localOne, positionOne);     // This gives us the position of the bars in x, y and z. Use only the y for now!
//  geom->GetCurrentMatrix()->LocalToMaster(localOther, positionOther);
//  TMS_Geom::GetInstance().Scale(positionOne[0], positionOne[1], positionOne[2]);
//  TMS_Geom::GetInstance().Scale(positionOther[0], positionOther[1], positionOther[2]);

  bool above = false;
  if ((UHit.GetNotZ() > 0 && VHit.GetNotZ() > 0) || (UHit.GetNotZ() < 0 && VHit.GetNotZ() < 0)) {
    if (std::abs(UHit.GetNotZ()) >= std::abs(VHit.GetNotZ()) && UHit.GetNotZ() < 0) above = true;
    else if (std::abs(UHit.GetNotZ()) >=std::abs( VHit.GetNotZ()) && UHit.GetNotZ() > 0) above = false;
    else if (std::abs(UHit.GetNotZ()) < std::abs(VHit.GetNotZ()) && UHit.GetNotZ() < 0) above = false;
    else if (std::abs(UHit.GetNotZ()) < std::abs(VHit.GetNotZ()) && UHit.GetNotZ() > 0) above = true;
  }
  else if (UHit.GetNotZ() > 0 && VHit.GetNotZ() < 0) above = true;
  else if (UHit.GetNotZ() < 0 && VHit.GetNotZ() > 0) above = false;

  if (above) {
    if (std::abs(UHit.GetNotZ() - VHit.GetNotZ()) > 167.0) {
      OneHit.SetRecoY(244.0);
    } else {
      OneHit.SetRecoY(TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() * 1000 + 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ()));
    }
  } else {
    if (std::abs(UHit.GetNotZ() - VHit.GetNotZ()) > 167.0) {
      OneHit.SetRecoY(-2949.0);
    } else {
      OneHit.SetRecoY(TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() * 1000 - 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ()));
    }
  }

  //OneHit.SetRecoY(244 - 0.5 * TMS_Manager::GetInstance().Get_Reco_TRACKMATCH_TiltAngle() * std::abs(UHit.GetNotZ() - VHit.GetNotZ()));  //positionOne[1], positionOther[1] //0.5 * (OneHit.GetY() + OtherHit.GetY())//-1350
  // USING THIS -2949/1350 ALSO IN CompareY FUNCTION. CHANGE THERE AS WELL IF CHANGED HERE!!!!!
  return;
}

void TMS_TrackFinder::CalculateRecoX(TMS_Hit &UHit, TMS_Hit &VHit, TMS_Hit &XHit) {
  XHit.SetRecoX(0.5 * (UHit.GetNotZ() + VHit.GetNotZ()));
  return;
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
    for (auto &Track: TotalCandidatesU) {
      // Loop over each hit in each track
      for (auto &Hit: Track) {
        // Get the bar of this hit
        TMS_Bar bar = Hit.GetBar();
        if (bar.Contains(x_true, z_true)) {
          nFoundHits++;
        }
      }
    }
    for (auto &Track: TotalCandidatesV) {
      // Loop over each hit in each track
      for (auto &Hit: Track) {
        TMS_Bar bar = Hit.GetBar();
        if (bar.Contains(x_true, z_true)) {
          nFoundHits++;
        }
      }
    }

    for (auto &Track: TotalCandidatesX) {
      for (auto &Hit: Track) {
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
double TMS_TrackFinder::CalculateTrackEnergy(const std::vector<TMS_Hit> &Candidate) {
  // Look at the reconstructred tracks
  if (Candidate.size() == 0) return -999.;

  // Sort by increasing z
  //std::sort((*Candidate).begin(), (*Candidate).end(), TMS_Hit::SortByZInc);

  double total = 0;
  // Find the track energy
  for (auto hit = Candidate.begin(); hit != Candidate.end(); ++hit) {
    total += (*hit).GetE();
  }
  return total;
}

double TMS_TrackFinder::CalculateTrackEnergy3D(const TMS_Track &Track3D) {
  // Look at the reconstructed tracks
  if ((Track3D.Hits).size() == 0) return -999.;
  
  double total = 0;
  // Loop over each Hough Candidate and find the track energy
  for (auto it = (Track3D.Hits).begin(); it != (Track3D.Hits).end(); ++it) {
    total += (*it).GetE();
  }  
   return total;
}

double TMS_TrackFinder::CalculateTrackKEByRange(const TMS_Track &Track3D) {
  if (Track3D.Length <= 0.0) return -999.0;

  return 82. + 1.75*Track3D.Length; // Magic Clarence number
}

// Calculate the track length for each track
double TMS_TrackFinder::CalculateTrackLength(const std::vector<TMS_Hit> &Candidate) {
  // Look at the reconstructed track
  if (Candidate.size() == 0) return -999.;

  double final_total = 0;
  int max_n_nodes_used = 0;

  // Loop through all bar types so we're only calculating along the same plane rotation type
  // Otherwise it would overestimate the track length by zig-zagging between u and v
  TMS_Bar::BarType bartypes[] = {TMS_Bar::kXBar, TMS_Bar::kYBar, TMS_Bar::kUBar, TMS_Bar::kVBar, TMS_Bar::kError};
  for (auto itbartype = std::begin(bartypes); itbartype != std::end(bartypes); ++itbartype) {
    // Get only the nodes with the current bar type
    auto track_hits_only_matching_bar = ProjectHits(Candidate, (*itbartype));
    double total = 0;
    int n_nodes = 0;
    for (auto hit = track_hits_only_matching_bar.begin(); hit != track_hits_only_matching_bar.end()
        && (hit+1) != track_hits_only_matching_bar.end(); ++hit) {
      auto nexthit = *(hit+1);
      // Use the geometry to calculate the track length between hits
      TVector3 point1((*hit).GetNotZ(), -200, (*hit).GetZ());
      TVector3 point2(nexthit.GetNotZ(), -200, nexthit.GetZ());
      // handle X bars correctly
      if ((*itbartype) == TMS_Bar::kXBar) {
        point1.SetXYZ(0, (*hit).GetNotZ(), (*hit).GetZ());
        point2.SetXYZ(0, nexthit.GetNotZ(), nexthit.GetZ());
      }
      double tracklength = TMS_Geom::GetInstance().GetTrackLength(point1, point2);
      total += tracklength;
      n_nodes += 1;
    }
    if (n_nodes > max_n_nodes_used) {
      final_total = total;
      max_n_nodes_used = n_nodes;
    }
  }
  return final_total;
}

double TMS_TrackFinder::CalculateTrackLength3D(const TMS_Track &Track3D) {
  //std::cout << "CalculateTrackLength3D" << std::endl;

  // Look at the reconstructed tracks
  if ((Track3D.Hits).size() == 0) return -999.;
  
  double final_total = 0;
  int max_n_nodes_used = 0;
  double total = 0;
  int n_nodes = 0;
  // Loop over each Hough Candidate and find the track length
  for (auto it = (Track3D.Hits).rbegin(); it != (Track3D.Hits).rend()
      && (it+1) != (Track3D.Hits).rend(); ++it) { // turn direction around
    auto nexthit = *(it+1); //+
    if ((*it).GetBar().GetBarType() != TMS_Bar::kXBar) {
      // Use the geometry to calculate the track length between hits
      TVector3 point1((*it).GetRecoX(), (*it).GetRecoY(), (*it).GetZ());
      TVector3 point2(nexthit.GetRecoX(), nexthit.GetRecoY(), nexthit.GetZ());
      if (nexthit.GetBar().GetBarType() == TMS_Bar::kXBar) point2.SetY(nexthit.GetNotZ());
      total += TMS_Geom::GetInstance().GetTrackLength(point2, point1);
      //std::cout << "new total (1): " << total << std::endl;
      n_nodes += 1;
    } else if ((*it).GetBar().GetBarType() == TMS_Bar::kXBar) {
      // Use the geometry to calculate the track length between hits
      TVector3 point1((*it).GetRecoX(), (*it).GetNotZ(), (*it).GetZ());
      TVector3 point2(nexthit.GetRecoX(), nexthit.GetNotZ(), nexthit.GetZ());
      if (nexthit.GetBar().GetBarType() != TMS_Bar::kXBar) point2.SetY(nexthit.GetRecoY());
      double tracklength = TMS_Geom::GetInstance().GetTrackLength(point2, point1);
      total += tracklength;
      //std::cout << "new total (2): " << total << std::endl;
      n_nodes += 1;
    }
  }
  if (n_nodes > max_n_nodes_used) {
    final_total = total;
    max_n_nodes_used = n_nodes;
  }

  return final_total;
}

std::vector<std::vector<TMS_Hit> > TMS_TrackFinder::HoughTransform(const std::vector<TMS_Hit> &TMS_Hits, const char &hitgroup) {
  // The returned vector of tracks
  std::vector<std::vector<TMS_Hit> > LineCandidates;

  // Check it's not empty
  if (TMS_Hits.empty()) return LineCandidates;

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned = CleanHits(TMS_Hits);
  if (TMS_Hits_Cleaned.empty()) return LineCandidates;

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kUBar);
  if (hitgroup == 'V') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kVBar);
  } else if (hitgroup == 'X') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);
  } else if (hitgroup == 'Y') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kYBar);
  }

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

  while (double(TMS_xz.size()) > nHits_Tol*nXZ_Hits_Start && 
      TMS_xz.size() > nMinHits && 
      nRuns < nMaxHough) {

    // The candidate vectors
    std::vector<TMS_Hit> TMS_xz_cand;
    if (TMS_xz.size() > 0) TMS_xz_cand = RunHough(TMS_xz, hitgroup);

    if (TMS_xz_cand.size() == 0) {
      nRuns++;
      if (hitgroup == 'U') {
        delete HoughLinesU.back().second;

        HoughLinesU.pop_back(); // Remove the built Hough line
      } else if (hitgroup == 'V') {
        delete HoughLinesV.back().second;
	
	      HoughLinesV.pop_back();
      } else if (hitgroup == 'X') {
        delete HoughLinesX.back().second;

        HoughLinesX.pop_back();
      } else if (hitgroup == 'Y') {
        delete HoughLinesY.back().second;

          HoughLinesY.pop_back();
      } else {
          std::cout << "Removing built Hough lines goes wrong for hitgroups: hitgroup = " << hitgroup << std::endl;
          break;
        }      
      break;
      }

    if (hitgroup == 'U') {
      // Move into the candidate vector
      for (auto &i: TMS_xz_cand) CandidatesU.push_back(std::move(i));
    } else if (hitgroup == 'V') {
      for (auto &i: TMS_xz_cand) CandidatesV.push_back(std::move(i));
    } else if (hitgroup == 'X') {
      for (auto &i: TMS_xz_cand) CandidatesX.push_back(std::move(i));
    } else if (hitgroup == 'Y') {
      for (auto &i: TMS_xz_cand) CandidatesY.push_back(std::move(i));
    }

    // Loop over vector and remove used hits
    for (auto jt = TMS_xz_cand.begin(); jt != TMS_xz_cand.end();++jt) {
      for (auto it = TMS_xz.begin(); it!= TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    
    if (hitgroup == 'U') {
      // Push back the candidates into the total candidates
      LineCandidates.push_back(std::move(CandidatesU));
    } else if (hitgroup == 'V') {
      LineCandidates.push_back(std::move(CandidatesV));
    } else if (hitgroup == 'X') {
      LineCandidates.push_back(std::move(CandidatesX));
    } else if (hitgroup == 'Y') {
      LineCandidates.push_back(std::move(CandidatesY));
    }
      
    nRuns++;
  }
#ifdef DEBUG
  std::cout << "Ran " << nRuns << " Hough algo" << std::endl;
#endif

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

      double HoughUInter_1 = 0.000;
      double HoughUSlope_1 = 0.000;
      double HoughVInter_1 = 0.000;
      double HoughVSlope_1 = 0.000;
      double HoughXInter_1 = 0.000;
      double HoughXSlope_1 = 0.000;
      double HoughYInter_1 = 0.000;
      double HoughYSlope_1 = 0.000;

      if (hitgroup == 'U') {
        HoughUInter_1 = HoughLinesU[lineit].second->GetParameter(0);
        HoughUSlope_1 = HoughLinesU[lineit].second->GetParameter(1);
      } else if (hitgroup == 'V') {
        HoughVInter_1 = HoughLinesV[lineit].second->GetParameter(0);
        HoughVSlope_1 = HoughLinesV[lineit].second->GetParameter(1);
      } else if (hitgroup == 'X') {
        HoughXInter_1 = HoughLinesX[lineit].second->GetParameter(0);
        HoughXSlope_1 = HoughLinesX[lineit].second->GetParameter(1);
      } else if (hitgroup == 'Y') {
        HoughYInter_1 = HoughLinesY[lineit].second->GetParameter(0);
        HoughYSlope_1 = HoughLinesY[lineit].second->GetParameter(1);
      }

      // Now loop over the remaining hits
      // Need to keep a conventional iterator to remove the HoughLine vector entries
      int lineit_2 = 0;
      for (std::vector<std::vector<TMS_Hit> >::iterator jt = LineCandidates.begin(); jt != LineCandidates.end(); ++jt, ++lineit_2) {
        
        // Make sure that we don't try to merge a track with itself
        if ((*jt) == (*it)) continue;

        if ((*jt).size() == 0) continue;

        // Make sure to only merge the smaller track into the bigger and not the other way around
        if ((*jt).size() > (*it).size()) continue;

        // Let's get the first and last hits of each track
        // Sort each in decreasing z
        SpatialPrio(*jt);

        // Get the first and last hit for this second track
        TMS_Hit firsthit_2 = (*jt).front();
        TMS_Hit lasthit_2 = (*jt).back();
        // Check that it's not the same hit (this doesn't seem to work. Added check for same track above)
        if (firsthit_2 == firsthit && lasthit_2 == lasthit) continue;

        int first_hit_z_2 = firsthit_2.GetPlaneNumber();
        int first_hit_notz_2 = firsthit_2.GetBarNumber();

        //std::cout << "Last hit (running track): " << last_hit_z << ", " << last_hit_notz << std::endl;
        //std::cout << "against" << std::endl;
        //std::cout << "First hit (t.b.del. track): " << first_hit_z_2 << ", " << first_hit_notz_2 << std::endl;

        // Now check how far away the hits are
        bool mergehits = (abs(first_hit_z_2 - last_hit_z) <= 2 && 
                          abs(first_hit_notz_2 - last_hit_notz) <= 2);

        double HoughUInter_2 = 0.000;
        double HoughUSlope_2 = 0.000;
        double HoughVInter_2 = 0.000;
        double HoughVSlope_2 = 0.000;
        double HoughXInter_2 = 0.000;
        double HoughXSlope_2 = 0.000;
        double HoughYInter_2 = 0.000;
        double HoughYSlope_2 = 0.000;

        if (hitgroup == 'U') {
          HoughUInter_2 = HoughLinesU[lineit_2].second->GetParameter(0);
          HoughUSlope_2 = HoughLinesU[lineit_2].second->GetParameter(1);
	      } else if (hitgroup == 'V') {
  	      HoughVInter_2 = HoughLinesV[lineit_2].second->GetParameter(0);
      	  HoughVSlope_2 = HoughLinesV[lineit_2].second->GetParameter(1);
      	} else if (hitgroup == 'X') {
          HoughXInter_2 = HoughLinesX[lineit_2].second->GetParameter(0);
          HoughXSlope_2 = HoughLinesX[lineit_2].second->GetParameter(1);
      	} else if (hitgroup == 'Y') {
          HoughYInter_2 = HoughLinesY[lineit_2].second->GetParameter(0);
          HoughYSlope_2 = HoughLinesY[lineit_2].second->GetParameter(1);
        }

        // Now check how similar the Hough lines are
        bool mergehough = false;
      	if (hitgroup == 'U') {
      	  mergehough = (fabs(HoughUInter_2 - HoughUInter_1) < 1000 &&   // 100 -> 1000
                       fabs(HoughUSlope_2 - HoughUSlope_1) < 0.1);  // 0.01 -> 0.1
      	} else if (hitgroup == 'V') {
      	  mergehough = (fabs(HoughVInter_2 - HoughVInter_1) < 1000 && // 100 -> 1000
			                 fabs(HoughVSlope_2 - HoughVSlope_1) < 0.1);    // 0.01 -> 0.1
      	} else if (hitgroup == 'X') {
          mergehough = (fabs(HoughXInter_2 - HoughXInter_1) < 1000 &&
                       fabs(HoughXSlope_2 - HoughXSlope_1) < 0.1);
      	} else if (hitgroup == 'Y') {
          mergehough = (fabs(HoughYInter_2 - HoughYInter_1) < 1000 &&
                       fabs(HoughYSlope_2 - HoughYSlope_1) < 0.1);
        }

        // Check if we should merge or not
        if (!mergehits && !mergehough) {
          continue;
        }

        // Copy over the contents of hits_2 into hits
        for (TMS_Hit movehits: (*jt)) {
          (*it).push_back(movehits);
        }

        // Clear out the original vector
        (*jt).clear();
        // After a merge we need to re-sort the original vector and recalculate the first and last hits
        SpatialPrio((*it));
        firsthit = (*it).front();
        lasthit = (*it).back();
        last_hit_z = lasthit.GetPlaneNumber();
        last_hit_notz = lasthit.GetBarNumber();
        // TODO recalculate slope and intercept of merged track?
      }
    }
  }

  // Loop over and remove the empty tracks
  int linenumber = 0;
  for (std::vector<std::vector<TMS_Hit> >::iterator it = LineCandidates.begin(); it != LineCandidates.end(); ) {
    if ((*it).size() == 0) {
      it = LineCandidates.erase(it);
      if (hitgroup == 'U') {
        HoughLinesU.erase(HoughLinesU.begin()+linenumber);
      } else if (hitgroup == 'V') {
        HoughLinesV.erase(HoughLinesV.begin()+linenumber);
      } else if (hitgroup == 'X') {
        HoughLinesX.erase(HoughLinesX.begin()+linenumber);
      } else if (hitgroup == 'Y') {
        HoughLinesY.erase(HoughLinesY.begin()+linenumber);
      }
    } else {
      ++it;
      ++linenumber;
    }
  }

  // Now finally connect the start and end points of each Hough line with Astar to get the most efficient path along the hough hits
  // This helps remove biases in the greedy hough adjacent hit merging
  // For A* after Hough, better to run with Euclidean metric, over-riding the previous
  // HoughTransform
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_RunAStar()) {
    int tracknumber = 0;
    for (std::vector<std::vector<TMS_Hit> >::iterator it = LineCandidates.begin(); it != LineCandidates.end(); ) {
      // Need to sort each line in z
      SpatialPrio(*it);
      std::vector<TMS_Hit> CleanedHough = RunAstar(*it);
      unsigned int ncleaned = CleanedHough.size();

      // Now replace the old hits with this cleaned version
      // Can remove Hough hits if not big enough
      // Also remember to remove the line
      if (ncleaned < 1) {
        it = LineCandidates.erase(it);
	      if (hitgroup == 'U') {
          HoughLinesU.erase(HoughLinesU.begin()+tracknumber);
	      } else if (hitgroup == 'V') {
	        HoughLinesV.erase(HoughLinesV.begin()+tracknumber);
        } else if (hitgroup == 'X') {
          HoughLinesX.erase(HoughLinesX.begin()+tracknumber);
        } else if (hitgroup == 'Y') {
          HoughLinesY.erase(HoughLinesY.begin()+tracknumber);
        }
      } else {
        *it = CleanedHough;
        it++;
        tracknumber++;
      }
    }
  }
//  std::cout<<"hit group:"<<hitgroup<<"# of Hough Line:"<<LineCandidates.size()<<std::endl;

  return LineCandidates;
};


void TMS_TrackFinder::MaskHits(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask) {

#ifdef DEBUG
  int sizebef = Orig.size();
  std::cout << "Size before masking: " << Orig.size() << std::endl;
  std::cout << "Mask size: " << Mask.size() << std::endl;
#endif

  // Loop over each hit in the mask
  for (auto &MaskHit: Mask) {
    for (auto it = Orig.begin(); it != Orig.end(); ) {
      if (*it == MaskHit) it = Orig.erase(it);
      else it++;
    }
  }

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
  DBSCAN.SetPoints(DB_Points);
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
std::vector<TMS_Hit> TMS_TrackFinder::RunHough(const std::vector<TMS_Hit> &TMS_Hits, const char &hitgroup) {
//TODO new orientation
  // Check if we're in XZ view
  bool IsXZ = (TMS_Hits.front().GetBar().GetBarType() == TMS_Bar::kUBar || TMS_Hits.front().GetBar().GetBarType() == TMS_Bar::kVBar || TMS_Hits.front().GetBar().GetBarType() == TMS_Bar::kXBar|| TMS_Hits.front().GetBar().GetBarType() == TMS_Bar::kYBar);

  // Recalculate Hough parameters event by event... not fully tested!
  bool VariableHough = false;
  if (VariableHough) {
    // Recalculate the min and max of slope and intercept for every event
    // Hits are ordered in z already
    double maxz = TMS_Const::TMS_Thin_Start;
    double minz = TMS_Const::TMS_Double_End;
    double min_notz = TMS_Const::TMS_End_Exact[0];
    double max_notz =  TMS_Const::TMS_Start_Exact[0];
    bool print = true;
    for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
      double z = (*it).GetZ();
      double not_z = (*it).GetNotZ();
      if (print) {
        print = false;
      }
      if (z > maxz) maxz = z;
      if (z < minz) minz = z;
      if (not_z > max_notz) max_notz = not_z;
      if (not_z < min_notz) min_notz = not_z;
    }
    // 1.8 comes from wanting to cover vertices that are maximally 80% between min and max z
    double slope = (max_notz - min_notz) * 1.8 / (maxz - minz);
    double intercept = -1 * slope * (minz + 0.8 * (maxz - minz));
#ifdef DEBUG
    std::cout << "***" << std::endl;
    std::cout << "z range: " << minz << " " << maxz << std::endl;
    std::cout << "notz range: " << min_notz << " " << max_notz << std::endl;
    std::cout << "slope: " << slope << std::endl;
    std::cout << "intercept: " << intercept << std::endl;
#endif
    // now try setting these as the max/min
    SlopeMin = -1 * slope;
    SlopeMax = slope;
    InterceptMin = -1 * intercept;
    InterceptMax = intercept;
    InterceptWidth = (InterceptMax - InterceptMin) / nIntercept;
    SlopeWidth = (SlopeMax - SlopeMin) / nSlope;
  }

  double slope, intercept;
  // Calculate the Hough lines
  GetHoughLine(TMS_Hits, slope, intercept);
  if (hitgroup == 'U') {
    HoughLineU->SetParameter(0, intercept);
    HoughLineU->SetParameter(1, slope);
  } else if (hitgroup == 'V') {
    HoughLineV->SetParameter(0, intercept);
    HoughLineV->SetParameter(1, slope);
  } else if (hitgroup == 'X') {
    HoughLineX->SetParameter(0, intercept);
    HoughLineX->SetParameter(1, slope);
  } else if (hitgroup == 'Y') {
    HoughLineY->SetParameter(0, intercept);
    HoughLineY->SetParameter(1, slope);
  }

  // Different fitting regions for XZ and YZ views: 
  // Most of the bending happens in xz, so fit only until the transition region. 
  // For the yz view, fit the entire region
  //if (IsXZ) HoughLine->SetRange(zMinHough, zMaxHough);
  //else HoughLine->SetRange(TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Double_End);
  
  TF1* HoughCopy = (TF1*)HoughLineU->Clone();

  if (hitgroup == 'U') {
    HoughLineU->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineU->Clone();
  } else if (hitgroup == 'V') {
    HoughLineV->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineV->Clone();
  } else if (hitgroup == 'X') {
    HoughLineX->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineX->Clone();
  } else if (hitgroup == 'Y') {
    HoughLineY->SetRange(zMinHough, zMaxHough);
    HoughCopy = (TF1*)HoughLineY->Clone();
  } else {
#ifdef DEBUG
    std::cout << "Something is going wrong with the assigning of hitgroup" << std::endl;
#endif
    return TMS_Hits;
  }

  std::pair<bool, TF1*> HoughPairs = std::make_pair(IsXZ, HoughCopy);
  if (hitgroup == 'U') {
    HoughLinesU.push_back(std::move(HoughPairs));
  } else if (hitgroup == 'V') {
    HoughLinesV.push_back(std::move(HoughPairs));
  } else if (hitgroup == 'X') {
    HoughLinesX.push_back(std::move(HoughPairs));
  } else if (hitgroup == 'Y') {
    HoughLinesY.push_back(std::move(HoughPairs));
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
    
    // Calculate 'x'-point of hit with Hough line
    double HoughPoint = -9999999999; //HoughLineOne->Eval(zhit);
    if (hitgroup == 'U') {
      HoughPoint = HoughLineU->Eval(zhit);
    } else if (hitgroup == 'V') {
      HoughPoint = HoughLineV->Eval(zhit);
    } else if (hitgroup == 'X') {
      HoughPoint = HoughLineX->Eval(zhit);   
    } else if (hitgroup == 'Y') {
      HoughPoint = HoughLineY->Eval(zhit);   
    }
    // Hough point is inside bar -> start clustering around bar
    // (Check if 'x'-point is inside hit bar)
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
  WalkDownStream(returned, HitPool);
  WalkUpStream(returned, HitPool);

  // (Asa's understanding)
  // Run DB scan on returned (Hough track hits) and remove all hits that are not
  // found with this scan from returned and move them back to the
  // HitPool (hits not in Hough track)
  // Afterwards the A* algorithm runs on returned hits and connects them before
  // finally runnning the extrapolation function on returned and HitPool
  //
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
          ClusterHits[Point.ClusterID - 1].push_back(std::move(hit));
        }
      }
    }
  }

  // After running DB scan on returned, removing hits from returned not 'connected' and add back to HitPool
  if (ClusterHits.size() > 0) {
    unsigned int largest = 0;
    for (unsigned int i = 0; i < ClusterHits.size(); ++i) if (ClusterHits[i].size() > ClusterHits[largest].size()) largest = i;
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
  // RunHough
  if (TMS_Manager::GetInstance().Get_Reco_HOUGH_RunAStar()) {
    SpatialPrio(returned);
    std::vector<TMS_Hit> vec = RunAstar(returned);

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
      returned = vec;
    }
  }

  // Now finally check if any hough tracks are too short, or have too few hits
  // Makes sense to do this right at the end when all the merging and cleaning has been run
  // Order them in z
  SpatialPrio(returned);

  // Extrapolation of tracks to catch missing hits at end/start of tracks
  if (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_Extrapolation()) {
    returned = Extrapolation(returned, HitPool);
  }

  if (returned.empty()) {
    std::vector<TMS_Hit> emptyvec;
    return emptyvec;
  }

  // Now check that the Hough candidates are long enough
  double xend = (returned).back().GetBarNumber();     //GetNotZ();
  double zend = (returned).back().GetPlaneNumber();   //GetZ();
  double xstart = (returned).front().GetBarNumber();  //GetNotZ();
  double zstart = (returned).front().GetPlaneNumber();//GetZ();
  double xdist = xstart - xend;
  double zdist = zstart - zend;
  double dist = sqrt(xdist * xdist + zdist * zdist);
  unsigned int nhits = (returned).size();
  // Calculate the minimum distance in planes and bars instead of physical distance
  if (dist < MinDistHough || nhits < nMinHits) {
    std::vector<TMS_Hit> emptyvec;
    return emptyvec;
  }

  return returned;
}

std::vector<TMS_Hit> TMS_TrackFinder::Extrapolation(const std::vector<TMS_Hit> &TrackHits, const std::vector<TMS_Hit> &Hitpool) {
  std::vector<TMS_Hit> returned = TrackHits;
  SpatialPrio(returned);
  // Get direction for correct extrapolation from first/last three hits
  struct {
    double slope;
    double intercept;
  } front, end;

  // Get first three hits
  std::vector<TMS_Hit> front_three;
  std::vector<TMS_Hit>::const_iterator it = TrackHits.begin();
  int loop_iterator = 0;
  while (loop_iterator < 3 && it != TrackHits.end()) {
    if (front_three.size() >= 1 && (*it).GetZ() == front_three[loop_iterator-1].GetZ()) {
      ++it;
      continue;
    }
    front_three.push_back((*it));
    ++loop_iterator;
    ++it;
  };
  
  // Not enough hits to do the extrapolation
  if (front_three.size() < 3) return returned;

  // TODO Shorten this with an external function that uses a flag for last/three hits
  // Get last three hits
  std::vector<TMS_Hit> last_three;
  loop_iterator = 0;
  std::vector<TMS_Hit>::const_reverse_iterator ir = TrackHits.rbegin();
  while (loop_iterator < 3 && ir != TrackHits.rend()) {
    if (last_three.size() >= 1 && (*ir).GetZ() == last_three[loop_iterator-1].GetZ()) {
      ++ir;
      continue;
    }
    last_three.push_back((*ir));
    ++loop_iterator;
    ++ir;
  }
  
  // Not enough hits to do the extrapolation
  if (last_three.size() < 3) return returned;

  // Calculate slopes and intercepts of connecting lines of last/first hits
  double slopes_front[2], slopes_end[2], intercepts_front[2], intercepts_end[2];
  for (int i = 0; i < 2; ++i) {
    slopes_front[i] = (front_three[i+1].GetNotZ() - front_three[i].GetNotZ()) / (front_three[i+1].GetZ() - front_three[i].GetZ());
    slopes_end[i] = (last_three[i+1].GetNotZ() - last_three[i].GetNotZ()) / (last_three[i+1].GetZ() - last_three[i].GetZ());
    intercepts_front[i] = front_three[i].GetNotZ() - front_three[i].GetZ() * slopes_front[i];
    intercepts_end[i] = last_three[i].GetNotZ() - last_three[i].GetZ() * slopes_end[i];
  }

  // Take average of connecting lines as direction
  front.slope = (slopes_front[0] + slopes_front[1]) / 2;
  front.intercept = (intercepts_front[0] + intercepts_front[1]) / 2;

  end.slope = (slopes_end[0] + slopes_end[1]) / 2;
  end.intercept = (intercepts_end[0] + intercepts_end[1]) / 2;

  // Calculate new candidate hits that are at most ExtrapolateDist + ExtrapolateLimit from end of track away (heuristic cost)
  // and with a higher z at most +/- 2 bar widths away from the direction line
  std::vector<TMS_Hit> end_extrapolation_cand;
  for (std::vector<TMS_Hit>::const_iterator it = Hitpool.begin(); it != Hitpool.end(); ++it) {
    // Check if hit is after the end of the track
    if ((*it).GetZ() > TrackHits.back().GetZ()) {
      // Check if within 4 bar widths above or below the direction line
      bool CloseBars = ((*it).GetNotZ() <= ((*it).GetZ() * end.slope + end.intercept + TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsEnd() * (*it).GetBar().GetNotZw()) &&
            (*it).GetNotZ() >= ((*it).GetZ() * end.slope + end.intercept - TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsEnd() * (*it).GetBar().GetNotZw()));
      if ((*it).GetBar().GetBarType() == TMS_Bar::kXBar) { // increase Distance limit by 2 to reflect the difference in BarNumber for X layers
        CloseBars = ((*it).GetNotZ() <= ((*it).GetZ() * end.slope + end.intercept + 2 * TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsEnd() * (*it).GetBar().GetNotZw()) &&
            (*it).GetNotZ() >= ((*it).GetZ() * end.slope + end.intercept - 2 * TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsEnd() * (*it).GetBar().GetNotZw()));
      }
      if (CloseBars) {
        // Calculate temporary node to check for distance
        aNode candidate((*it).GetPlaneNumber(), (*it).GetBarNumber());
        aNode track_end(TrackHits.back().GetPlaneNumber(), TrackHits.back().GetBarNumber());
        candidate.SetHeuristic(kHeuristic);
        candidate.SetHeuristicCost(track_end);
#ifdef DEBUG
        std::cout << "Track end heuristic Cost: " << candidate.HeuristicCost << " "  << TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()
          << " + " << TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateLimit() << std::endl;
#endif
        // Check if node is within ExtrapolateDist + ExtrapolateLimit from end of track
        bool MaxDistance = (candidate.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist() + 
              TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateLimit());
        if ((*it).GetBar().GetBarType() == TMS_Bar::kXBar) {
          MaxDistance = (candidate.HeuristicCost <= 2 * (TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist() + 
              TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateLimit()));
        }
        if (MaxDistance) {
          // Move hit now into candidate hits
#ifdef DEBUG
          std::cout << "Added to candidates" << std::endl;
#endif
          end_extrapolation_cand.push_back((*it));
        }
      }
    }
  }

  // If more than 2 candidate hits, run A* algorithm to connect the correct ones
  if (end_extrapolation_cand.size() > 2) {
#ifdef DEBUG
    std::cout << "more than 2 candidates: " << end_extrapolation_cand.size() << std::endl;
    std::cout << "track initial size: " << TrackHits.size() << std::endl;
#endif
    // Make sure the hits are ordered
    SpatialPrio(end_extrapolation_cand);

    // Make sure the start is at most ExtrapolateDist away from end of track
    aNode test(end_extrapolation_cand.front().GetPlaneNumber(), end_extrapolation_cand.front().GetBarNumber());
    aNode track_end(returned.back().GetPlaneNumber(), returned.back().GetBarNumber());
    test.SetHeuristic(kHeuristic);
    test.SetHeuristicCost(track_end);
    if (test.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()) {
      std::vector<TMS_Hit> vec = RunAstar(end_extrapolation_cand);

      // Now add the connected hits into the existing track
      for (auto hit = vec.rbegin(); hit != vec.rend(); ++hit) {
        returned.push_back((*hit));
      }
    }
#ifdef DEBUG
    std::cout << "Hits added. Size now: " << returned.size() << std::endl;
#endif
    // Now order the hits in the existing track
    SpatialPrio(returned);
  } else if (!end_extrapolation_cand.empty()) { // If less than 2 candidate hits, but there are some, just add them
#ifdef DEBUG
    std::cout << "less than 2 candidates: " << end_extrapolation_cand.size() << std::endl;
    std::cout << "track initial size: " << TrackHits.size() << std::endl;
#endif

    // Make sure the hits are still ordered
    SpatialPrio(end_extrapolation_cand);

    // Make sure the start is at most ExtrapolateDist away from end of track
    aNode test(end_extrapolation_cand.front().GetPlaneNumber(), end_extrapolation_cand.front().GetBarNumber());
    aNode track_end(returned.back().GetPlaneNumber(), returned.back().GetBarNumber());
    test.SetHeuristic(kHeuristic);
    test.SetHeuristicCost(track_end);
    if (test.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()) {
      // Add them
      for (auto hit = end_extrapolation_cand.rbegin(); hit != end_extrapolation_cand.rend(); ++hit) {
        returned.push_back((*hit));
      }
#ifdef DEBUG
      std::cout << "Hits added. Size now: " << returned.size() << std::endl;
#endif

      // Now order the hits in the existing track
      SpatialPrio(returned);
    }
  }

    // Do the same as for the end of the track just the other direction for the start
    // Calculate new candidate hits that are at most ExtrapolateDist + ExtrapolateLimit from start of track away (Heuristic cost)
    // and with a smaller z at most +/- 2 bar widths away from the direction line
    std::vector<TMS_Hit> front_extrapolation_cand;
    for (std::vector<TMS_Hit>::const_iterator it = Hitpool.begin(); it != Hitpool.end(); ++it) {
      if (returned.empty()) {
        break;
      }
      // Check if hit is before the start of the track
      if ((*it).GetZ() < returned.front().GetZ()) {
        // Check if within 2 bar widths above or below the direction line
        bool CloseBars = ((*it).GetNotZ() <= ((*it).GetZ() * front.slope + front.intercept + TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsStart() * (*it).GetBar().GetNotZw()) &&
              (*it).GetNotZ() >= ((*it).GetZ() * front.slope + front.intercept - TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsStart() * (*it).GetBar().GetNotZw()));
        if ((*it).GetBar().GetBarType() == TMS_Bar::kXBar) {
          CloseBars = ((*it).GetNotZ() <= ((*it).GetZ() * front.slope + front.intercept + 2 * TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsStart() * (*it).GetBar().GetNotZw()) &&
              (*it).GetNotZ() >= ((*it).GetZ() * front.slope + front.intercept - 2 * TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_NumBarsStart() * (*it).GetBar().GetNotZw()));
        }
        if (CloseBars) {
          // Calculate temporary node to check for distance
          aNode candidate((*it).GetPlaneNumber(), (*it).GetBarNumber());
          aNode track_start((returned).front().GetPlaneNumber(), (returned).front().GetBarNumber());
          candidate.SetHeuristic(kHeuristic);
          candidate.SetHeuristicCost(track_start);
#ifdef DEBUG
          std::cout << "Heuristic Cost: " << candidate.HeuristicCost << " " << TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()
            << " + " << TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateLimit() << " (" << (*it).GetNotZ() << "|" << (*it).GetZ() << ")" << std::endl;
#endif
          // Check if node is within ExtrapolateDist + ExtrapolateLimit from start of track
          if (candidate.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist() +
              TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateLimit()) {
            // Move hit now into candidate hits
#ifdef DEBUG
              std::cout << "Added to candidates" << std::endl;
#endif
            front_extrapolation_cand.push_back((*it));
          }
        }
      }
    }

    // If more than 2 candidate hits, run A* algorithm to connect the correct ones
    if (front_extrapolation_cand.size() > 2 ) {
#ifdef DEBUG
      std::cout << "more than 2 candidates: " << front_extrapolation_cand.size() << std::endl;
      std::cout << "track initial size: " << returned.size() << std::endl;
#endif
  
      // Make sure the hits are ordered
      SpatialPrio(front_extrapolation_cand); 

      // Make sure the end is at most ExtrapolateDist away from start of track
      aNode test(front_extrapolation_cand.back().GetPlaneNumber(), front_extrapolation_cand.back().GetBarNumber());
      aNode track_start(returned.front().GetPlaneNumber(), returned.front().GetBarNumber());
      test.SetHeuristic(kHeuristic);
      test.SetHeuristicCost(track_start);
      if (test.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()) {
        std::vector<TMS_Hit> vec = RunAstar(front_extrapolation_cand);
  
        // Now add the connected hits into the existing track
        for (auto hit = vec.begin(); hit != vec.end(); ++hit)  {
          returned.push_back((*hit));
        }
#ifdef DEBUG
        std::cout << "Hits added. Size now: " << returned.size() << std::endl;
#endif

        // Now order the hits in the existing track
        SpatialPrio(returned);
      }
    } else if (!front_extrapolation_cand.empty()) { // If less than 2 candidate hits, but there are some, just add them
#ifdef DEBUG
      std::cout << "less than 2 candidates: " << front_extrapolation_cand.size() << std::endl;
      std::cout << "track initial size: " << returned.size() << std::endl;
#endif
      // Make sure the hits are still ordered
      SpatialPrio(front_extrapolation_cand);
    
      // Make sure the end is at most ExtrapolateDist away from start of track
      aNode test(front_extrapolation_cand.back().GetPlaneNumber(), front_extrapolation_cand.back().GetBarNumber());
      aNode track_start(returned.front().GetPlaneNumber(), returned.front().GetBarNumber());
      test.SetHeuristic(kHeuristic);
      test.SetHeuristicCost(track_start);
      if (test.HeuristicCost <= TMS_Manager::GetInstance().Get_Reco_EXTRAPOLATION_ExtrapolateDist()) {
        // Add them
        for (auto hit = front_extrapolation_cand.begin(); hit != front_extrapolation_cand.end(); ++hit) {
          returned.push_back((*hit));
        } 
#ifdef DEBUG
        std::cout << "Hits added. Size now: " << returned.size() << std::endl;
#endif

        // Now order the hits in the existing track
        SpatialPrio(returned);
    }
  }
  return returned;
}


// Find the bin for the accumulator
int TMS_TrackFinder::FindBin(double c) {
  // Since we're using uniform binning no need for binary search or similar
  if (c > InterceptMax) c = InterceptMax;
  int bin = (c - InterceptMin) / InterceptWidth;
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
      accumulator->SetBinContent(i + 1, j + 1, Accumulator[i][j]);
    }
  }

  int maxx, maxy, maxz;
  accumulator->GetMaximumBin(maxx, maxy, maxz);
  double maxtheta = accumulator->GetXaxis()->GetBinCenter(maxx);
  double maxrho = accumulator->GetYaxis()->GetBinCenter(maxy);
  accumulator->SetTitle(Form("#splitline{%s}{m_{max}=%.2f c_{max}=%.2f}", accumulator->GetTitle(), maxtheta, maxrho));

  // Set the minimum (easier to draw)
  double maxcounts = accumulator->GetMaximum();
  accumulator->SetMinimum(maxcounts * 0.8);

  return accumulator;
}

// Implement A* algorithm for track finding, starting with most upstream to most downstream hit
void TMS_TrackFinder::BestFirstSearch(const std::vector<TMS_Hit> &TMS_Hits, const char &hitgroup) {

  // Set the Heuristic cost calculator

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned = CleanHits(TMS_Hits);

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kUBar);
  if (hitgroup == 'V') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kVBar);
  } else if (hitgroup == 'X') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);
  } else if (hitgroup == 'Y') {
    TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kYBar);
  }
  //std::vector<TMS_Hit> TMS_yz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);

  // Do a spatial analysis of the hits in y and x around low z to ignore hits that are disconnected from other hits
  // Includes a simple sort in decreasing z (plane number)
  SpatialPrio(TMS_xz);
  //SpatialPrio(TMS_yz);

  int nRuns = 0;
  int nXZ_Hits_Start = TMS_xz.size();

  while (double(TMS_xz.size()) > nHits_Tol * nXZ_Hits_Start && 
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
    for (auto jt = AStarHits_xz.begin(); jt != AStarHits_xz.end(); ++jt) {
      for (auto it = TMS_xz.begin(); it != TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    // Only push back if we have more than one candidate
    if (AStarHits_xz.size() > nMinHits) {
      if (hitgroup == 'U') {    
        HoughCandidatesU.push_back(std::move(AStarHits_xz));
      } else if (hitgroup == 'V') {
        HoughCandidatesV.push_back(std::move(AStarHits_xz));
      } else if (hitgroup == 'X') {
        HoughCandidatesX.push_back(std::move(AStarHits_xz));
      } else if (hitgroup == 'Y') {
        HoughCandidatesY.push_back(std::move(AStarHits_xz));
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
        (*it).SetE((*it).GetE() + e2);
        (*it).SetT(((*it).GetT() + t2) / 2);
        nDuplicates++;
      }
    }
    // Now remove the duplicates
    if (nDuplicates > 0) {
      it = TMS_Hits_Cleaned.erase(it, it + nDuplicates);
    } else it++;
  }

  // Remove zero entries
  // Strip out hits that are outside the actual volume 
  // This is probably some bug in the geometry that sometimes gives hits in the z=30k mm (i.e. 10m downstream of the end of the TMS)
  // Figure out why these happen?
  for (std::vector<TMS_Hit>::iterator jt = TMS_Hits_Cleaned.begin(); 
      jt != TMS_Hits_Cleaned.end(); ) {

    if ( (*jt).GetZ() > TMS_Const::TMS_End[2] ||  // Sometimes a hit downstream of the end geometry
        (*jt).GetZ() < TMS_Const::TMS_Start[2]) { // Or upstream of the start...
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
  bool IsXZ = ((TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kUBar || (TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kVBar ||(TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kYBar);
  bool IsX = (TMS_xz.front().GetBar().GetBarType() == TMS_Bar::kXBar);
  // Reset remembering where gaps are in xz
  if (IsXZ) PlanesNearGap.clear();
  if (IsX) PlanesNearGap.clear();

  // Set the first and last hit to calculate the heuristic to
  //aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetNotZ(), TMS_xz.back().GetNotZw());
  aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetBarNumber());
//  if (IsXZ) std::cout<<"Last:"<<TMS_xz.back().GetPlaneNumber()<<std::endl;

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
#ifdef DEBUG
    std::cout << "Node creation" << std::endl;
    std::cout << "coordinates: x, y = " << (*it).GetZ() << " " << (*it).GetNotZ() << std::endl;
    TempNode.Print();
#endif

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
//    if (IsXZ) std::cout<<"layer:"<<layer<<"first layer"<<firstlayer<<std::endl;
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
            abs((*jt).y - (*it).y) > 3) continue; // (*jt).y - (*jt).y doesn't really make sense...
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
//        cost_so_far[neighbour.first->NodeID] = new_cost; // Update the cost to get to this node via this path
        cost_so_far[neighbour.first->NodeID] = cost_so_far[current.NodeID] + neighbour.second; // Update the cost to get to this node via this path
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
    //std::cout << "Node: " << NodeID << std::endl;
    //Nodes[NodeID].Print();
    //std::cout << "Came from: " << std::endl;
    NodeID = came_from[NodeID];
    //Nodes[NodeID].Print();
    // If the current node came from itself, we've reached the end
    if (NodeID == came_from[NodeID]) {
      //std::cout << "broke" << std::endl;
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


// This is only usable for U or V layers!!!!
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
    int NeighbourIndex = i - 1;
    // Previous plane number
    int PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();

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

    double grad_exp = (y - yprev) / (x - xprev);

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
      double grad_new = (ycand - y) / (xcand - x);

      //std::cout << "testing merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;
      //std::cout << "With hit info " << hits.GetZ() << ", " << hits.GetNotZ() << std::endl;
      //std::cout << "and new gradient: " << grad_new << std::endl;

      // If gradient is within 0.6 and the correct sign, accept
      if (fabs(grad_new - grad_exp) <= 2.0) { 
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

      double grad_exp = (y - yprev) / (x - xprev);

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

        double grad_new = (ycand - y) / (xcand - x);

        //std::cout << "testing merging " << PlaneNumber_cand << "," << hits.GetBarNumber() << std::endl;
        //std::cout << "With hit info " << hits.GetZ() << ", " << hits.GetNotZ() << std::endl;
        //std::cout << "and new gradient: " << grad_new << std::endl;

        // If gradient is within 1 and the correct sign, accept
        if (fabs(grad_new - grad_exp) <= 1.5) {
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

void TMS_TrackFinder::GetHoughLine(const std::vector<TMS_Hit> &TMS_Hits, double &slope, double &intercept) {
  // Reset the accumulator
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      Accumulator[i][j] = 0;
    }
  }

  // First run a simple Hough Transform
  for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Hit hit = (*it);
    double xhit = hit.GetNotZ();
    double zhit = hit.GetZ();

    // If z position is above region of interest, ignore hit
    //if (IsXZ && zhit > zMaxHough) continue;
    if (zhit > zMaxHough) continue;

    Accumulate(xhit, zhit);
  }

  // Find the maximum of the accumulator and which m,c bin the maximum occurs in
  double max_zy = 0;
  int max_zy_slope_bin = 0;
  int max_zy_inter_bin = 0;
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      if (Accumulator[i][j] > max_zy) {
        max_zy = Accumulator[i][j];
        max_zy_slope_bin = i;
        max_zy_inter_bin = j;
      }
    }
  }
  intercept = InterceptMin + max_zy_inter_bin * (InterceptMax - InterceptMin) / nIntercept;
  slope = SlopeMin + max_zy_slope_bin * (SlopeMax - SlopeMin) / nSlope;
}

void TMS_TrackFinder::Accumulate(double xhit, double zhit) {
  // Could probably multi-thread this operation
  // Now do the Hough
  for (int i = 0; i < nSlope; ++i) {
    double m = SlopeMin + i * SlopeWidth;
    if (m > SlopeMax) m = SlopeMax;
    
    // Now calculate rho
    double c = xhit - m  * zhit;
    if (c > InterceptMax) c = InterceptMax;

    // Find which rho bin this corresponds to
    int c_bin = FindBin(c);
      /*
      if (i > nSlope || c_bin > nIntercept) {
      std::cout << "c: " << c << std::endl;
      std::cout << "m: " << m << std::endl;
      std::cout << "i: " << i << std::endl;
      std::cout << "cbin: " << c_bin << std::endl;
      }
      */
    // Fill the accumulator, but only within bounds
    // We don't care about lines outside of bounds
    if (i >= 0 && c_bin >= 0 && i < nSlope && c_bin < nIntercept) Accumulator[i][c_bin]++;
  }
}

