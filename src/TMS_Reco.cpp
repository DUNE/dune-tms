#include "TMS_Reco.h"

TMS_TrackFinder::TMS_TrackFinder() :
  nIntercept(2E3),
  nSlope(2E3),
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
  IsGreedy(TMS_Manager::GetInstance().Get_Reco_ASTAR_IsGreedy())
{
  // Apply the maximum Hough transform to the zx not zy: all bending happens in zx
  Accumulator = new int*[nSlope];
  for (int i = 0; i < nSlope; ++i) {
    Accumulator[i] = new int[nIntercept];
  }

  // Keep track of the tracking efficiency
  Efficiency = new TH1D("Efficiency", "Efficiency;T_{#mu} (GeV); Efficiency", 30, 0, 6);
  Total = new TH1D("Total", "Total;T_{#mu} (GeV); Total", 30, 0, 6);

  HoughLine = new TF1("LinearHough", "[0]+[1]*x", TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_Start);
  HoughLine->SetLineStyle(kDashed);
  HoughLine->SetLineColor(kMagenta-9);

  DBSCAN.SetEpsilon(TMS_Manager::GetInstance().Get_Reco_DBSCAN_Epsilon());
  DBSCAN.SetMinPoints(TMS_Manager::GetInstance().Get_Reco_DBSCAN_MinPoints());

}

// The generic track finder
void TMS_TrackFinder::FindTracks(TMS_Event &event) {

  // Check through the Houghlines
  for (auto &i: HoughLines) {
    delete i.second;
  }

  // Reset the candidate vector
  Candidates.clear();
  RawHits.clear();
  TotalCandidates.clear();
  HoughLines.clear();
  HoughCandidates.clear();
  ClusterCandidates.clear();

  // Get the raw unmerged and untracked hits
  RawHits = event.GetHits();

  // Clean hits (check duplicate hits, and energy threshold)
  CleanedHits = CleanHits(RawHits);
  // Require N hits after cleaning
  if (RawHits.size() < nMinHits) return;

#ifdef DEBUG
  std::cout << "Raw hits: " << RawHits.size() << std::endl;
  std::cout << "Cleaned hits: " << CleanedHits.size() << std::endl;
#endif

  // Set the number of merges in the Hough transform to be zero
  // Hough
  //HoughTransform(CleanedHits);
  BestFirstSearch(CleanedHits);

  std::vector<TMS_Hit> Masked = CleanedHits;
  // Loop over the Hough candidates
  for (auto Lines: HoughCandidates) {
#ifdef DEBUG
    std::cout << "Masked size bef: " << Masked.size() << std::endl;
#endif
    MaskHits(Masked, Lines);
#ifdef DEBUG
    std::cout << "Masked size aft: " << Masked.size() << std::endl;
#endif
  }

#ifdef DEBUG
  std::cout << "Masked hits: " << Masked.size() << std::endl;
#endif

  // Try finding some clusters after the Hough Transform
  ClusterCandidates = FindClusters(Masked);

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
  // Also
  
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
  //
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
  double muonmom = muon.GetInitialTMSMomentum().Mag();

  double muonke = sqrt(muonmom*muonmom+TMS_KinConst::mass_mu*TMS_KinConst::mass_mu) - TMS_KinConst::mass_mu;
  muonke /= 1E3;
  if (efficiency > 0) {
    Efficiency->Fill(muonke, efficiency);
    Total->Fill(muonke);
  }
}

void TMS_TrackFinder::HoughTransform(const std::vector<TMS_Hit> &TMS_Hits) {

  // Check it's not empty
  if (TMS_Hits.empty()) return;

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned = CleanHits(TMS_Hits);
  if (TMS_Hits_Cleaned.empty()) return;

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kYBar);
  //std::vector<TMS_Hit> TMS_yz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);

  // Do a spatial analysis of the hits in y and x around low z to ignore hits that are disconnected from other hits
  // Includes a simple sort in decreasing z (plane number)
  SpatialPrio(TMS_xz);
  //SpatialPrio(TMS_yz);

  int nXZ_Hits_Start = TMS_xz.size();
#ifdef DEBUG
  std::cout << "Starting Hough transform with: " << nXZ_Hits_Start << " hits" << std::endl;
#endif
  //int nYZ_Hits_Start = TMS_yz.size();

  // First find where to do the Hough transform
  // Check how many continous hits there are in the first layers from the first hit
  /*
  int nCont = 0;
  int HitNumber = 0;
  for (auto &hitref: TMS_xz) {
    // The previous layer
    int PrevLayer = hitref.GetPlaneNumber();
    nCont = 0;
    for (auto &hit: TMS_xz) {
      double planeno = hit.GetPlaneNumber();
      double z = hit.GetZ();
      // Run only until we hit the thick target
      // Could also check planenumber
      if (z > TMS_Const::TMS_Thick_Start || planeno > TMS_Const::TMS_nThinPlanes) break;
      // Check continuity
      if (planeno != PrevLayer+1) continue;

      // Increment the number of continuous hits
      nCont++;
      // Update the previous layer
      PrevLayer = planeno;
    }
    // If we've found a sequence of hits that are more than needed
    if (nCont > nThinCont) {
      break;
    }
    HitNumber++;
  }

  // If we have continuos hits in the thin region, run the Hough transform only in that region
  // if we don't , run over whole region
  zMinHough = TMS_xz[HitNumber].GetZ();
  if (nCont < nThinCont) {
    zMaxHough = TMS_Const::TMS_Thick_End;
  } else {
    zMaxHough = TMS_Const::TMS_Thick_Start;
  }
  */

  // We'll be moving out TMS_xz and TMS_yz and putting them into candidates
  // Keep running successive Hough transforms until we've covered 80% of hits (allow for maximum 4 runs)
  int nRuns = 0;

  while (double(TMS_xz.size()) > nHits_Tol*nXZ_Hits_Start && 
         TMS_xz.size() > nMinHits && 
         nRuns < nMaxHough) {

    // The candidate vectors
    std::vector<TMS_Hit> TMS_xz_cand;
    if (TMS_xz.size() > 0) TMS_xz_cand = RunHough(TMS_xz);

    // If we're running out of hits (Hough transform doesn't have enough hits)
    if (TMS_xz_cand.size() < nMinHits) {
      nRuns++;
      delete HoughLines.back().second;
      HoughLines.pop_back(); // Remove the built Hough line
      break;
    }

    // Check if there are any candidates
    // If not, they might all be downstream, i.e. this is not a LAr event, it's a TMS event

    //std::vector<TMS_Hit> TMS_yz_cand;
    //if (TMS_yz.size() > 0) TMS_yz_cand = RunHough(TMS_yz, nRuns);

    // Then order them in z
    SpatialPrio(TMS_xz_cand);
    //SpatialPrio(TMS_yz_cand);

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
    HoughCandidates.push_back(std::move(Candidates));
    nRuns++;
  }
#ifdef DEBUG
  std::cout << "Ran " << nRuns << " Hough algo" << std::endl;
#endif
};

void TMS_TrackFinder::MaskHits(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask) {

#ifdef DEBUG
  int sizebef = Orig.size();
  std::cout << "size before masking: " << sizebef << std::endl;
  std::cout << "Mask size: " << Mask.size() << std::endl;
#endif
  // Loop over the
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
    // Get the point from the hit

    // Only call it "x" since we're plotting z on the x-axis
    //double x = i.GetZ()/i.GetZw();

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

  std::vector<TMS_DBScan_Point> NoisePoints = DBSCAN.GetNoise();
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

std::vector<TMS_Hit> TMS_TrackFinder::RunHough(const std::vector<TMS_Hit> &TMS_Hits) {

  // Check if we're in XZ view
  bool IsXZ = ((TMS_Hits[0].GetBar()).GetBarType() == TMS_Bar::kYBar);

  // Reset the accumulator
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      Accumulator[i][j] = 0;
    }
  }

  // First run a simple Hough Transform
  for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it!=TMS_Hits.end(); ++it) {
    TMS_Hit hit = (*it);
    double xhit = hit.GetNotZ();
    double zhit = hit.GetZ();

    // If z position is above region of interest, ignore hit
    if (IsXZ && zhit > zMaxHough) continue;

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

  double InterceptOpt_zy = InterceptMin+max_zy_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
  double SlopeOpt_zy = SlopeMin+max_zy_slope_bin*(SlopeMax-SlopeMin)/nSlope;
  HoughLine->SetParameter(0, InterceptOpt_zy);
  HoughLine->SetParameter(1, SlopeOpt_zy);

  // Different fitting regions for XZ and YZ views: 
  // Most of the bending happens in xz, so fit only until the transition region. 
  // For the yz view, fit the entire region
  //if (IsXZ) HoughLine->SetRange(zMinHough, zMaxHough);
  //else HoughLine->SetRange(TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_End);

  HoughLine->SetRange(zMinHough, zMaxHough);
  TF1 *HoughCopy = (TF1*)HoughLine->Clone();

  std::pair<bool, TF1*> HoughPairs = std::make_pair(IsXZ, HoughCopy);
  HoughLines.push_back(std::move(HoughPairs));

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

    double HoughPoint = HoughLine->Eval(zhit);

    // Hough point is inside bar -> start clustering around bar
    if (( bar.Contains(HoughPoint, zhit) ||
        bar.Contains(HoughPoint-bar.GetNotZw(), zhit) ||
        bar.Contains(HoughPoint+bar.GetNotZw(), zhit) )) {
      bool IsGood = true;
      if (returned.size() > 5 && 
          abs(returned.back().GetPlaneNumber() - hit.GetPlaneNumber()) > 5) {
        IsGood = false;
      }

      if (IsGood) {
        returned.push_back(std::move(hit));
        // Remove from pool of hits
        it = HitPool.erase(it);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }

  // Now walk along the Hough hits and add on adjacent hits
  // Hough is most likely to find hits upstream due to bending being less there
  WalkDownStream(returned, HitPool);
  WalkUpStream(returned, HitPool);

  /*
  // Loop over the candidates, and add new adjacent candidates to the end
  size_t CandSize = returned.size();
  for (size_t i = 0; i < CandSize; ++i) {
    TMS_Hit Candidate = returned[i];
    TMS_Bar CandidateBar = Candidate.GetBar();
    int CandidatePlaneNumber = Candidate.GetPlaneNumber();

    // Count the number of times a candidate merges adjacent hits
    unsigned int nMerges = 0;
    // Is this hit close to an airgap?
    //bool InGap = Candidate.NextToGap();

    // Now loop over each hit
    for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {

      // If we've already merged more than we're allowed to
      if (nMerges >= nMaxMerges) break;
      TMS_Hit PoolHit = (*jt);
      //TMS_Bar PoolBar = PoolHit.GetBar();
      int PoolPlaneNumber = PoolHit.GetPlaneNumber();

      // Now check the distance in x or y depending on bar
      double PoolPos = PoolHit.GetNotZ();
      //double PoolPosWidth = PoolHit.GetNotZw();

      // Ensure adjacent or matching in z
      // Make an exception for the airgaps; allow any range in z
      //if (!InGap && abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
      //if (!InGap && abs(CandidatePlaneNumber-PoolPlaneNumber) > 1) {
      if (abs(CandidatePlaneNumber-PoolPlaneNumber) > 1) {
        ++jt;
        continue;
      }

      // Should we merge this hit
      bool Merge = false;
      // If this hasn't been merged already, attempt to do so
      if (nMerges < nMaxMerges) {
        // If the Candidate bar contains the pool bar + it's width, they're adjacent
        // Or even better, the adjacent hits have maxing not z

        // Preferentially merge in z (same NotZ position)
        if (CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
          Merge = true;
        // Then check merge in NotZ
        } 
        else if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || 
            CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ())) {
          Merge = true;
          // Then check two bars away for the xz view
        } else if (IsXZ && 
          (CandidateBar.Contains(PoolPos + 2*PoolPosWidth, CandidateBar.GetZ()) || 
          CandidateBar.Contains(PoolPos - 2*PoolPosWidth, CandidateBar.GetZ()))) {
          Merge = true;
          } else if (CandidateBarType == TMS_Bar::kYBar && 
          (CandidateBar.Contains(PoolPos + 3*PoolPosWidth, CandidateBar.GetZ()) || 
          CandidateBar.Contains(PoolPos - 3*PoolPosWidth, CandidateBar.GetZ()))) {
          Merge = true;
          }

        // Make a special arrangement for if we're next to the gap
           else if (InGap && ( 
           PoolBar.Contains(TMS_Const::TMS_Dead_Top[1]+PoolPosWidth, PoolBar.GetZ()) ||
           PoolBar.Contains(TMS_Const::TMS_Dead_Top[0]-PoolPosWidth, PoolBar.GetZ()) ||
           PoolBar.Contains(TMS_Const::TMS_Dead_Center[1]+PoolPosWidth, PoolBar.GetZ()) ||
           PoolBar.Contains(TMS_Const::TMS_Dead_Center[0]-PoolPosWidth, PoolBar.GetZ()) ||
           PoolBar.Contains(TMS_Const::TMS_Dead_Bottom[1]+PoolPosWidth, PoolBar.GetZ()) ||
           PoolBar.Contains(TMS_Const::TMS_Dead_Bottom[0]-PoolPosWidth, PoolBar.GetZ()) )) {
           Merge = true;
           }
      }

      if (Merge) {
        returned.push_back(std::move(PoolHit));
        // Increment the number of candidates in the vector
        CandSize++;
        jt = HitPool.erase(jt);
        // Increment the merges for this candidate
        nMerges++;
      } else {
        ++jt;
      }
    }
    }
    */

    return returned;

    // Now see if the yz and xz views roughly agree on stop and start positions of the main track
    // If not, it implies there may be two tracks and the yz view has selected one whereas xz selected theother, or a broken track in one of the views
    /*
       int StartPlane_xz = 99999999;
       int EndPlane_xz = -99999999;

       for (std::vector<TMS_Hit>::iterator it = Candidates.begin(); it != Candidates.end(); ++it) {
       TMS_Bar Bar = (*it).GetBar();
    // Compare the plane number
    int PlaneNumber = Bar.GetPlaneNumber();
    if (PlaneNumber < StartPlane_yz) StartPlane_yz = PlaneNumber;
    if (PlaneNumber > EndPlane_yz) EndPlane_yz = PlaneNumber;
    }

    // Propagate the xz/yz start and endpoints here
    if (abs(StartPlane_yz - StartPlane_xz) != 1 || abs(EndPlane_yz - EndPlane_xz) != 1) {
    //Matched_yz_xz = false;
    std::cout << "Not matching end/start plane" << std::endl;
    std::cout << "StartPlane yz: " << StartPlane_yz << ", EndPlane yz: " << EndPlane_yz << std::endl;
    std::cout << "StartPlane xz: " << StartPlane_xz << ", EndPlane xz: " << EndPlane_xz << std::endl;
    }
    */

    /*
    // Do a Hough transform on the remaining hits
    // If this Hough transform matches within some number the original Hough transform parameters, merge the tracks
    // Reset the accumulator (or make a new one?)
    for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
    Accumulator_zy[i][j] = 0;
    Accumulator_zx[i][j] = 0;
    }
    }

    // First run a simple Hough Transform
    for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end(); ++it) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    //double EnergyDeposit = hit.GetE();
    //double Time = hit.GetT();
    double xhit = bar.GetX();
    double yhit = bar.GetY();
    double zhit = bar.GetZ();

    TMS_Bar::BarType Type = bar.GetBarType();
    // Scan over the entire range this time
    // Or maybe the z region in which the other view's track has been found?

    Accumulate(xhit, yhit, zhit, Type);
    }

    max_zy = 0;
    max_zx = 0;
    max_zy_slope_bin = 0;
    max_zy_inter_bin = 0;
    max_zx_slope_bin = 0;
    max_zx_inter_bin = 0;
    for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
    if (Accumulator_zy[i][j] > max_zy) {
    max_zy = Accumulator_zy[i][j];
    max_zy_slope_bin = i;
    max_zy_inter_bin = j;
    }
    if (Accumulator_zx[i][j] > max_zx) {
    max_zx = Accumulator_zx[i][j];
    max_zx_slope_bin = i;
    max_zx_inter_bin = j;
    }
    }
    }

    double InterceptOpt_zx_2 = InterceptMin+max_zx_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
    double SlopeOpt_zx_2 = SlopeMin+max_zx_slope_bin*(SlopeMax-SlopeMin)/nSlope;

    double InterceptOpt_zy_2 = InterceptMin+max_zy_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
    double SlopeOpt_zy_2 = SlopeMin+max_zy_slope_bin*(SlopeMax-SlopeMin)/nSlope;

    TF1 *HoughLine_zy_2 = (TF1*)HoughLine_zy->Clone("HZY_2");
    HoughLine_zy_2->SetLineColor(kWhite);
    HoughLine_zy_2->SetParameter(0, InterceptOpt_zy_2);
    HoughLine_zy_2->SetParameter(1, SlopeOpt_zy_2);

    TF1 *HoughLine_zx_2 = (TF1*)HoughLine_zx->Clone("HZX_2");
    HoughLine_zx_2->SetLineColor(kWhite);
    HoughLine_zx_2->SetParameter(0, InterceptOpt_zx_2);
    HoughLine_zx_2->SetParameter(1, SlopeOpt_zx_2);

    // Get the hits that intersect the second hough line
    // First make a second candidate vector
    std::vector<TMS_Hit> Candidates2;
    for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end(); ) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    double zhit = bar.GetZ();
    TMS_Bar::BarType Type = bar.GetBarType();
    if (Type == TMS_Bar::kXBar) {
      HoughLine = HoughLine_zy_2;
    } else if (Type == TMS_Bar::kYBar) {
      HoughLine = HoughLine_zx_2;
    }

    double HoughPoint = HoughLine->Eval(zhit);

    // Hough point is inside bar -> start clustering around bar
    if (bar.Contains(HoughPoint, zhit)) {
      Candidates2.push_back(std::move(hit));
      // Remove from pool of hits
      it = HitPool.erase(it);
    } else {
      ++it;
    }
}

// Loop over the candidates, and add new adjacent candidates to the end
CandSize = Candidates2.size();
for (size_t i = 0; i < CandSize; ++i) {
  TMS_Hit Candidate = Candidates2[i];
  TMS_Bar CandidateBar = Candidate.GetBar();
  int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
  TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

  // Count the number of times a candidate merges adjacent hits
  unsigned int nMerges = 0;
  // Now loop over each hit
  for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {
    TMS_Hit PoolHit = (*jt);
    TMS_Bar PoolBar = PoolHit.GetBar();
    int PoolPlaneNumber = PoolBar.GetPlaneNumber();
    TMS_Bar::BarType PoolBarType = PoolBar.GetBarType();

    // Only match the same type of bars (x bars with x bars, y bars with ybars)
    if (PoolBarType != CandidateBarType) {
      ++jt;
      continue;
    }

    // Ensure adjacent or matching in z
    // Need some exception here for the airgaps! FIGURE IT OUT
    if (abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
      ++jt;
      continue;
    }

    // Now check the distance in x or y depending on bar
    double PoolPos = PoolBar.GetNotZ();
    double PoolPosWidth = PoolBar.GetNotZw();
    // If the Candidate bar contains the pool bar + it's width, they're adjacent
    // Or even better, the adjacent hits have maxing not z
    if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
      Candidates2.push_back(std::move(PoolHit));
      // Increment the number of candidates in the vector
      CandSize++;
      jt = HitPool.erase(jt);
      // Increment the merges for this candidate
      nMerges++;
    } else {
      ++jt;
    }
  }
}

// Now we *should* have two candidates?
std::cout << "Total hits: " << TMS_Hits.size() << std::endl;
std::cout << "Candidates (Hough+Cluster): " << Candidates.size() << std::endl;
std::cout << "Candidates2 (Leftover+Hough+cluster): " << Candidates2.size() << std::endl;
std::cout << "Hitpool: " << HitPool.size() << std::endl;

// Intersection in z (x-axis) of the two Hough lines
double intersection_zx = (InterceptOpt_zx-InterceptOpt_zx_2)/(SlopeOpt_zx_2-SlopeOpt_zx);
double intersection_zy = (InterceptOpt_zy-InterceptOpt_zy_2)/(SlopeOpt_zy_2-SlopeOpt_zy);
std::cout << "Intersection zx: " << intersection_zx << std::endl;
std::cout << "Intersection zy: " << intersection_zy << std::endl;

TotalCandidates.push_back(Candidates);
//TotalCandidates.push_back(Candidates2);

HoughLines_zy.push_back(HoughLine_zy);
//HoughLines_zy.push_back(HoughLine_zy_2);

HoughLines_zx.push_back(HoughLine_zx);
//HoughLines_zx.push_back(HoughLine_zx_2);

// If this Hough transform matches the start and end position in z of the other view, choose the longest track as the muon and make 2 candidate tracks



// Last scan tries looks at hits that occur around the gap regions
// Then use this to inform the yz reconstruction (since broken tracks will be due to holes in x, not y)


// Run a Hough transform on the remaining HitPool to find second tracks
*/


/*
   std::cout << "After Hough+adjacent checks: " << std::endl;
   std::cout << "Of " << TMS_Hits.size()  << " hits " << Candidates.size() << " are candidates and " << HitPool.size() << " are left" << std::endl;
   */


}

/*
   void TMS_TrackFinder::MergeAdjacent(std::vector<TMS_Hit> Candidates, std::vector<TMS_Hit> Pool) {

// Loop over the candidates, and add new adjacent candidates to the end
size_t CandSize = Candidates.size();
for (size_t i = 0; i < CandSize; ++i) {
TMS_Hit Candidate = Candidates[i];
TMS_Bar CandidateBar = Candidate.GetBar();
int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

// Count the number of times a candidate merges adjacent hits
int nMerges = 0;
// Now loop over each hit
for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {
TMS_Hit PoolHit = (*jt);
TMS_Bar PoolBar = PoolHit.GetBar();
int PoolPlaneNumber = PoolBar.GetPlaneNumber();
TMS_Bar::BarType PoolBarType = PoolBar.GetBarType();

// Only match the same type of bars (x bars with x bars, y bars with ybars)
if (PoolBarType != CandidateBarType) {
++jt;
continue;
}

// Ensure adjacent or matching in z
// Need some exception here for the airgaps! FIGURE IT OUT
if (abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
++jt;
continue;
}

// Now check the distance in x or y depending on bar
double PoolPos = PoolBar.GetNotZ();
double PoolPosWidth = PoolBar.GetNotZw();
// If the Candidate bar contains the pool bar + it's width, they're adjacent
// Or even better, the adjacent hits have maxing not z
if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
Candidates.push_back(std::move(PoolHit));
// Increment the number of candidates in the vector
CandSize++;
jt = HitPool.erase(jt);
// Increment the merges for this candidate
nMerges++;
} else {
++jt;
}
}
}
}
*/

// Find the bin for the accumulator
int TMS_TrackFinder::FindBin(double c) {
  // Since we're using uniform binning no need for binary search or similar
  if (c > InterceptMax) c = InterceptMax;
  int bin = (c-InterceptMin)/InterceptWidth;
  return bin;
}

// xvalue is x-axis, y value is y-axis
void TMS_TrackFinder::Accumulate(double xhit, double zhit) {

  // Could probably multi-thread this operation
  // Now do the Hough
  for (int i = 0; i < nSlope; ++i) {
    double m = SlopeMin+i*SlopeWidth;
    if (m > SlopeMax) m = SlopeMax;

    // Now calculate rho
    double c = xhit-m*zhit;
    if (c > InterceptMax) c = InterceptMax;

    // Find which rho bin this corresponds to
    int c_bin = FindBin(c);

    /*
    if (i > nSlope || c_bin > nIntercept) {
      std::cout << "c: " << c << std::endl;
      std::cout << "m: " << m << std::endl;
      std::cout << "i: " <<  i << std::endl;
      std::cout << "cbin: " << c_bin << std::endl;
    }
    */

    // Fill the accumulator
    Accumulator[i][c_bin]++;
  }
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
void TMS_TrackFinder::BestFirstSearch(const std::vector<TMS_Hit> &TMS_Hits) {

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
    if (AStarHits_xz.size() > nMinHits) HoughCandidates.push_back(std::move(AStarHits_xz));
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
      if (z == z2 && y == y2 && fabs(t2-(*it).GetT()) < TMS_Const::TMS_TimeThreshold) {
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
std::vector<TMS_Hit> TMS_TrackFinder::RunAstar(const std::vector<TMS_Hit> &TMS_xz) {

  // Remember which orientation these hits are
  // needed when we potentially skip the air gap in xz (but not in yz!)
  bool IsXZ = ((TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kYBar);
  // Reset remembering where gaps are in xz
  if (IsXZ) PlanesNearGap.clear();

  // Set the first and last hit to calculate the heuristic to
  //aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetNotZ(), TMS_xz.back().GetNotZw());
  aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetBarNumber());

  //std::cout << "Last hit we calculate heuristic cost to" << std::endl;
  //Last.Print();

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
      // Only connect adjacent nodes to make calculation much faster
      if (abs((*jt).x - (*it).x) > 2 || abs((*jt).y - (*jt).y) > 2) continue;
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

// Do a quick analysis of the hits that are sorted *DESCENDING* in z
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
  SpatialPrio(full);
  SpatialPrio(vec);

  // Now walk along vec (start at second element due to derivative calculation
  std::vector<TMS_Hit>::size_type size = vec.size();
  for (std::vector<TMS_Hit>::size_type i = 1; i < size; ++i) {
    // Use to calculate gradient to next point
    double x = vec[i].GetZ();
    double y = vec[i].GetNotZ();
    int PlaneNumber = vec[i].GetPlaneNumber();

    // Require continous hits for derivative calculation
    // If we can't find the previous hit at first go (might be two hits in the same plane, so adjacent hits would have an inifite slop), try moving down in the order
    int NeighbourIndex = i-1;
    // Previous plane number
    int PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();

    bool MoveOn = false;
    while (PlaneNumber - PlaneNumber_prev != 1) {
      NeighbourIndex--;
      if (NeighbourIndex < 0) {
        MoveOn = true;
        break;
      }
      PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
      if (PlaneNumber - PlaneNumber_prev > 1) {
        MoveOn = true;
        break;
      }
    }

    // Couldn't find the right hit
    if (MoveOn) continue;

    // Calculate expected derivative from upstream hit
    double xprev = vec[NeighbourIndex].GetZ();
    double yprev = vec[NeighbourIndex].GetNotZ();

    double grad_exp = (y-yprev)/(x-xprev);

    // Save the indices of which particle was best
    // Only allow for one merge

    // Matching 
    for (std::vector<TMS_Hit>::iterator it = full.begin(); it != full.end(); ) {
      TMS_Hit hits = *it;
      int PlaneNumber_cand = hits.GetPlaneNumber();
      // Always look downstream
      if (PlaneNumber_cand - PlaneNumber != 1) {
        if (PlaneNumber_cand <= PlaneNumber) {
          ++it;
          continue;
          // Since the candidates plane numbers are ordered in z, once we encounter a higher z, all the remaining ones also won't be mergeable, so save some time by not scanning them
        } else {
          break;
        }
      }
      // Allow for broken track in z?
      double xcand = hits.GetZ();
      double ycand = hits.GetNotZ();

      double grad_new = (ycand-y)/(xcand-x);

      // If gradient is within 1 and the correct sign, accept
      if (fabs(grad_new-grad_exp) <= 0.6) {
        // If the gradient is zero we shouldn't do a sign check
        // But if it's not, check the sign of the gradient doesn't flip
        if (fabs(grad_new) > TMS_Const::TMS_Small_Num && 
            fabs(grad_exp) > TMS_Const::TMS_Small_Num &&
            TMS_Utils::sgn(grad_new) != TMS_Utils::sgn(grad_exp)) {
          ++it;
          continue;
        }
        // Now need to find position to insert
        // Guaranteed sorted in z
        // Need to find how many duplicates in z there are; don't necessarily want to insert right after this hit if it is followed by a hit in the same z
        // If this is the largest z, just pop to the back
        if (xcand >= vec.back().GetZ()) {
          vec.push_back(std::move(hits));
          // Need to decrement iterator over line
          // Sometimes we've missed one hit between Hough hits
        } else {
          for (std::vector<TMS_Hit>::iterator jt = vec.begin()+i; jt != vec.end(); ++jt) {
            double compx = (*jt).GetZ();
            // Look one ahead
            double compx1 = (*(jt-1)).GetZ();
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
        // Remove off the end
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

  // Sort in *decreasing* z (starting at highest z)
  std::sort(vec.begin(), vec.end(), TMS_Hit::SortByZ);
  std::sort(full.begin(), full.end(), TMS_Hit::SortByZ);

  // Now walk along vec (start at second element due to derivative calculation
  std::vector<TMS_Hit>::size_type size = vec.size();
  for (std::vector<TMS_Hit>::size_type i = 1; i < size; ++i) {
    // Use to calculate gradient to next point
    double x = vec[i].GetZ();
    double y = vec[i].GetNotZ();
    int PlaneNumber = vec[i].GetPlaneNumber();

    // Require continous hits for derivative calculation
    // If we can't find the previous hit at first go (might be two hits in the same plane, so adjacent hits would have an inifite slop), try moving down in the order
    int NeighbourIndex = i-1;
    // Previous plane number
    int PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();

    bool MoveOn = false;
    while (PlaneNumber_prev - PlaneNumber != 1) {
      NeighbourIndex--;
      if (NeighbourIndex < 0) {
        MoveOn = true;
        break;
      }
      PlaneNumber_prev = vec[NeighbourIndex].GetPlaneNumber();
      if (PlaneNumber_prev - PlaneNumber > 1) {
        MoveOn = true;
        break;
      }
    }

    // Couldn't find the right hit
    if (MoveOn) continue;

    // Calculate expected derivative from upstream hit
    double xprev = vec[NeighbourIndex].GetZ();
    double yprev = vec[NeighbourIndex].GetNotZ();

    double grad_exp = (y-yprev)/(x-xprev);
    // Matching 
    for (std::vector<TMS_Hit>::iterator it = full.begin(); it != full.end(); ) {
      TMS_Hit hits = *it;
      int PlaneNumber_cand = hits.GetPlaneNumber();
      // Always look downstream
      if (PlaneNumber - PlaneNumber_cand != 1) {
        if (PlaneNumber_cand >= PlaneNumber) {
          ++it;
          continue;
          // Since the candidates plane numbers are ordered in z, once we encounter a higher z, all the remaining ones also won't be mergeable, so save some time by not scanning them
        } else {
          break;
        }
      }
      // Allow for broken track in z?
      double xcand = hits.GetZ();
      double ycand = hits.GetNotZ();

      double grad_new = (ycand-y)/(xcand-x);

      // If gradient is within 1 and the correct sign, accept
      if (fabs(grad_new-grad_exp) <= 0.6) {
        // If the gradient is zero we shouldn't do a sign check
        // But if it's not, check the sign of the gradient doesn't flip
        if (fabs(grad_new) > TMS_Const::TMS_Small_Num && 
            fabs(grad_exp) > TMS_Const::TMS_Small_Num &&
            TMS_Utils::sgn(grad_new) != TMS_Utils::sgn(grad_exp)) {
          ++it;
          continue;
        }
        // Now need to find position to insert
        // Guaranteed sorted in z
        // Need to find how many duplicates in z there are; don't necessarily want to insert right after this hit if it is followed by a hit in the same z
        // If this is the largest z, just pop to the back
        if (xcand <= vec.back().GetZ()) {
          vec.push_back(std::move(hits));
          // Need to decrement iterator over line
          // Sometimes we've missed one hit between Hough hits
        } else {
          for (std::vector<TMS_Hit>::iterator jt = vec.begin()+i; jt != vec.end(); ++jt) {
            double compx = (*jt).GetZ();
            // Look one behind
            double compx1 = (*(jt-1)).GetZ();
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

