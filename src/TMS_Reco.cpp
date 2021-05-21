#include "TMS_Reco.h"

TMS_TrackFinder::TMS_TrackFinder() :

  nIntercept(2E3),
  nSlope(2E3),
  InterceptMin(-40E3),
  InterceptMax(+40E3),
  SlopeMin(-2.0),
  SlopeMax(2.0),
  InterceptWidth((InterceptMax-InterceptMin)/nIntercept),
  SlopeWidth((SlopeMax-SlopeMin)/nSlope),
  // Max z for us to do Hough in, here choose transition layer
  zMinHough(TMS_Const::TMS_Thin_Start),
  zMaxHough(TMS_Const::TMS_Thick_Start+1000),
  nMaxHough(2),
  nThinCont(10),
  nHits_Tol(0.5),
  // Minimum number of hits required to run track finding
  nMinHits(10),
  // Maximum number of merges for one hit
  nMaxMerges(1),
  IsGreedy(false),
  // Initialise Highest cost to be very large
  HighestCost(999)
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
}

// The generic track finder
void TMS_TrackFinder::FindTracks(TMS_Event &event) {

  // Check through the Houghlines
  for (auto i : HoughLines) {
    delete i.second;
  }

  // Reset the candidate vector
  Candidates.clear();
  RawHits.clear();
  TotalCandidates.clear();
  HoughLines.clear();

  // Get the raw unmerged and untracked hits
  std::vector<TMS_Hit> TMS_Hits = event.GetHits();

  // Require 10 hits
  if (TMS_Hits.size() < nMinHits) return;

  // A star
  //BestFirstSearch(TMS_Hits);

  // Hough
  HoughTransform(TMS_Hits);

  // For future probably want to move track candidates into the TMS_Event class
  //EvaluateTrackFinding(event);

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
  std::vector<TMS_TrueParticle> TrueParticles = event.GetParticles();
  TMS_TrueParticle muon;
  bool FoundMuon = false;
  for (auto particle: TrueParticles) {
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
  for (auto TrueHit : muon.GetPositionPoints()) {
    // Find if the true hit sits inside this hit's bar
    double x_true = TrueHit.X();
    //double y_true = TrueHit.Y();
    double z_true = TrueHit.Z();
    //double t_true = TrueHit.T();
    //std::cout << "Finding " << x_true << " " << z_true << std::endl;
    for (auto Track : TotalCandidates) {
      // Loop over each hit in each track
      for (auto Hit : Track) {
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

  //int nYZ_Hits_Start = TMS_yz.size();
  int nXZ_Hits_Start = TMS_xz.size();

  // Check how many continous hits there are in the first layers
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
  std::cout << "Found " << nCont << " continuous hits" << std::endl;

  // If we have continuos hits in the thin region, run the Hough transform only in that region
  // if we don't , run over whole region
  zMinHough = TMS_xz[HitNumber].GetZ();
  if (nCont < nThinCont) {
    zMaxHough = TMS_Const::TMS_Thick_End;
  } else {
    zMaxHough = TMS_Const::TMS_Thick_Start+1000;
  }

  // We'll be moving out TMS_xz and TMS_yz and putting them into candidates
  // Keep running successive Hough transforms until we've covered 80% of hits (allow for maximum 4 runs)
  int nRuns = 0;

  while (double(TMS_xz.size()) > nHits_Tol*nXZ_Hits_Start && nRuns < nMaxHough) {

    // The candidate vectors
    std::vector<TMS_Hit> TMS_xz_cand;
    //std::vector<TMS_Hit> TMS_yz_cand;
    if (TMS_xz.size() > 0) TMS_xz_cand = RunHough(TMS_xz);
    //if (TMS_yz.size() > 0) TMS_yz_cand = RunHough(TMS_yz, nRuns);

    // Check if there are any candidates
    // If not, they might all be downstream, i.e. this is not a LAr event, it's a TMS event

    // Then order them in z
    SpatialPrio(TMS_xz_cand);
    //SpatialPrio(TMS_yz_cand);

    // Move into the candidate vector
    for (auto &i : TMS_xz_cand) Candidates.push_back(std::move(i));
    //for (auto &i : TMS_yz_cand) Candidates.push_back(i);
    // Loop over vector and remove used hits
    /*
       for (auto jt = TMS_yz_cand.begin(); jt != TMS_yz_cand.end();++jt) {
       for (auto it = TMS_yz.begin(); it!= TMS_yz.end();) {
       if ((*it) == (*jt)) it = TMS_yz.erase(it);
       else it++;
       }
       }
       */

    // Loop over vector and remove used hits
    for (auto jt = TMS_xz_cand.begin(); jt != TMS_xz_cand.end();++jt) {
      for (auto it = TMS_xz.begin(); it!= TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    // Push back the candidates into the total candidates
    TotalCandidates.push_back(std::move(Candidates));
    nRuns++;
  }
  std::cout << "Ran " << nRuns << " Hough transforms" << std::endl;
}


std::vector<TMS_Hit> TMS_TrackFinder::RunHough(const std::vector<TMS_Hit> &TMS_Hits) {

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
  if (IsXZ) HoughLine->SetRange(zMinHough, zMaxHough);
  else HoughLine->SetRange(TMS_Const::TMS_Thin_Start, TMS_Const::TMS_Thick_End);

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
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end();) {
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
    if (bar.Contains(HoughPoint, zhit)) {
      returned.push_back(std::move(hit));
      // Remove from pool of hits
      it = HitPool.erase(it);
    } else {
      ++it;
    }
  }

  // Loop over the candidates, and add new adjacent candidates to the end
  size_t CandSize = returned.size();
  for (size_t i = 0; i < CandSize; ++i) {
    TMS_Hit Candidate = returned[i];
    TMS_Bar CandidateBar = Candidate.GetBar();
    int CandidatePlaneNumber = Candidate.GetPlaneNumber();

    // Count the number of times a candidate merges adjacent hits
    unsigned int nMerges = 0;
    // Is this hit close to an airgap?
    bool InGap = Candidate.NextToGap();

    // Now loop over each hit
    for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {

      // If we've already merged more than we're allowed to
      if (nMerges >= nMaxMerges) break;
      TMS_Hit PoolHit = (*jt);
      TMS_Bar PoolBar = PoolHit.GetBar();
      int PoolPlaneNumber = PoolHit.GetPlaneNumber();

      // Now check the distance in x or y depending on bar
      double PoolPos = PoolHit.GetNotZ();
      double PoolPosWidth = PoolHit.GetNotZw();

      // Ensure adjacent or matching in z
      // Make an exception for the airgaps; allow any range in z
      //if (!InGap && abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
      if (!InGap && abs(CandidatePlaneNumber-PoolPlaneNumber) > 1) {
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
        } else if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || 
            CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ())) {
          Merge = true;
          // Then check two bars away for the xz view
        } else if (IsXZ && 
            (CandidateBar.Contains(PoolPos + 2*PoolPosWidth, CandidateBar.GetZ()) || 
             CandidateBar.Contains(PoolPos - 2*PoolPosWidth, CandidateBar.GetZ()))) {
          Merge = true;
        } /*else if (CandidateBarType == TMS_Bar::kYBar && 
            (CandidateBar.Contains(PoolPos + 3*PoolPosWidth, CandidateBar.GetZ()) || 
            CandidateBar.Contains(PoolPos - 3*PoolPosWidth, CandidateBar.GetZ()))) {
            Merge = true;
            }
            */

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

    if (i > nSlope || c_bin > nIntercept) {
      std::cout << "c: " << c << std::endl;
      std::cout << "m: " << m << std::endl;
      std::cout << "i: " <<  i << std::endl;
      std::cout << "cbin: " << c_bin << std::endl;
    }

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
  std::vector<TMS_Hit> TMS_yz = ProjectHits(TMS_Hits_Cleaned, TMS_Bar::kXBar);

  // Do a spatial analysis of the hits in y and x around low z to ignore hits that are disconnected from other hits
  // Includes a simple sort in decreasing z (plane number)
  SpatialPrio(TMS_xz);
  SpatialPrio(TMS_yz);

  int nRuns = 0;
  int nXZ_Hits_Start = TMS_xz.size();

  while (TMS_xz.size() > 0.2*nXZ_Hits_Start && nRuns < 4) {
    //std::cout << "yz = " << TMS_yz.size() << " before a = " << a << std::endl;
    //std::cout << "xz = " << TMS_xz.size() << " before a = " << a << std::endl;
    std::vector<TMS_Hit> AStarHits_yz;
    std::vector<TMS_Hit> AStarHits_xz;

    // Run on x-z first since that's where the gaps occur, leading to broken tracks
    // We can save where this happens in z and make the yz reconstruction aware of the gap
    if (TMS_xz.size() > 0) AStarHits_xz = RunAstar(TMS_xz);
    if (TMS_yz.size() > 0) AStarHits_yz = RunAstar(TMS_yz);

    // Then order them in z
    SpatialPrio(AStarHits_xz);
    //SpatialPrio(AStarHits_yz);

    // Copy over to the candidates
    for (auto &i : AStarHits_yz) Candidates.push_back(std::move(i));
    for (auto &i : AStarHits_xz) Candidates.push_back(std::move(i));

    // Loop over vector and remove used hits
    for (auto jt = AStarHits_yz.begin(); jt!=AStarHits_yz.end();++jt) {
      for (auto it = TMS_yz.begin(); it!= TMS_yz.end();) {
        if ((*it) == (*jt)) it = TMS_yz.erase(it);
        else it++;
      }
    }

    // Loop over vector and remove used hits
    for (auto jt = AStarHits_xz.begin(); jt!=AStarHits_xz.end();++jt) {
      for (auto it = TMS_xz.begin(); it!= TMS_xz.end();) {
        if ((*it) == (*jt)) it = TMS_xz.erase(it);
        else it++;
      }
    }
    TotalCandidates.push_back(std::move(Candidates));
    nRuns++;
  }
  std::cout << "Ran " << nRuns << " Astar algo" << std::endl;
}

// Remove duplicate hits
std::vector<TMS_Hit> TMS_TrackFinder::CleanHits(const std::vector<TMS_Hit> &TMS_Hits) {

  std::vector<TMS_Hit> TMS_Hits_Cleaned;
  if (TMS_Hits.empty()) return TMS_Hits_Cleaned;

  for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Hit hit = *it;
    TMS_Bar bar = hit.GetBar();
    // Maybe this hit has already been counted
    double z = bar.GetZ();
    double y = bar.GetNotZ();
    std::pair<double, double> temp = std::make_pair(z,y);
    bool Duplicate = false;

    // Skim out the duplicates
    for (std::vector<TMS_Hit>::const_iterator jt = TMS_Hits_Cleaned.begin(); jt != TMS_Hits_Cleaned.end(); ++jt) {
      TMS_Hit hit2 = *jt;
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      std::pair<double, double> temp2 = std::make_pair(z2,y2);
      if (temp == temp2) {
        Duplicate = true;
        break;
      }
    }

    if (!Duplicate) {
      TMS_Hits_Cleaned.push_back(std::move(hit));
    }
  }

  // Remove zero entries
  // Strip out hits that are outside the actual volume 
  // This is probably some bug in the geometry that sometimes gives hits in the z=30k mm (i.e. 10m downstream of the end of the TMS)
  // Figure out why these happen?
  for (std::vector<TMS_Hit>::iterator jt = TMS_Hits_Cleaned.begin(); jt != TMS_Hits_Cleaned.end(); ) {
    if ((*jt).GetZ() == 0 || (*jt).GetNotZ() == 0 || (*jt).GetZ() > TMS_Const::TMS_End_z) {
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

std::vector<TMS_Hit> TMS_TrackFinder::RunAstar(const std::vector<TMS_Hit> &TMS_xz) {

  // Remember which orientation these hits are
  // needed when we potentially skip the air gap in xz (but not in yz!)
  bool IsXZ = ((TMS_xz[0].GetBar()).GetBarType() == TMS_Bar::kYBar);
  // Reset remembering where gaps are in xz
  if (IsXZ) PlanesNearGap.clear();

  // Set the first and last hit to calculate the heuristic to
  aNode Last(TMS_xz.back().GetPlaneNumber(), TMS_xz.back().GetNotZ(), TMS_xz.back().GetNotZw());

  // Also convert the TMS_Hit to our path-node
  std::vector<aNode> Nodes;

  // Give each node an ID
  int NodeID = 0;
  for (std::vector<TMS_Hit>::const_iterator it = TMS_xz.begin(); it != TMS_xz.end(); ++it, NodeID++) {

    // Use the x position as plane number
    double x = (*it).GetPlaneNumber();
    double y = (*it).GetNotZ();
    double yw = (*it).GetNotZw();

    // Make the node
    aNode TempNode(x, y, yw, NodeID);
    // Calculate the Heuristic cost for each of the nodes to the last point
    TempNode.SetHeuristicCost(Last);

    Nodes.push_back(std::move(TempNode));
  }

  // Can probably evaluate the last node here: if heuristic is large for the 5 nearest hits it's likely wrong?


  // Now find the neighbours of each node
  for (std::vector<aNode>::iterator it = Nodes.begin(); it != Nodes.end(); ++it) {

    double x = (*it).x;
    double y = (*it).y;

    // Check if we're close to the gap if in xz view
    bool IsNextToGap = false;
    if (IsXZ) {
      IsNextToGap = NextToGap(y, (*it).yw);
      if (IsNextToGap) {
        PlanesNearGap.push_back(x);
      }
      // If we're in yz view, see if this z could be close to a gap
    } else {
      IsNextToGap = (std::find(PlanesNearGap.begin(), PlanesNearGap.end(), x-1) != PlanesNearGap.end() ||
          std::find(PlanesNearGap.begin(), PlanesNearGap.end(), x+1) != PlanesNearGap.end() );
    }

    // If a hit is close to a gap and connects to another hit close to a gap, ONLY allow for it to connect to the closest hit
    bool FoundNeighbourNearGap = false;

    for (std::vector<aNode>::iterator jt = Nodes.begin(); jt != Nodes.end(); ++jt) {

      // First check this node isn't itself!
      if ((*jt) == (*it)) continue;

      // Get the address
      aNode* Pointer = &(*jt);

      double xcan = Pointer->x;
      int DeltaPlane = xcan-x;

      double y_n = (*jt).y;
      double DeltaY = fabs(y_n-y);

      // Get the number of vertical bars (NotZw gives bar width)
      int DeltaY_Trunc = int(DeltaY/(*jt).yw);

      // Allow for three bars in y at maximum
      if (abs(DeltaY_Trunc) > 3) continue;

      // Find the other hits that are next to the gaps
      bool IsNextToGap_cand = false;
      // If in xz view we just check the x coordinate and see if it's next to the gap
      // Only allow it tracking downstream
      if (IsXZ) {
        //IsNextToGap_cand = ( DeltaPlane < 2 && NextToGap(y_n, (*jt).yw) );
        IsNextToGap_cand = ( DeltaPlane < 1 && NextToGap(y_n, (*jt).yw) );
      } else {
        // Only allow the gap exception if really needed in z
        // Only allow it tracking downstream
        //IsNextToGap_cand = ( DeltaPlane < 2 && (
        IsNextToGap_cand = ( DeltaPlane < 1 && (
              std::find(PlanesNearGap.begin(), PlanesNearGap.end(), xcan-1) != PlanesNearGap.end() ||
              std::find(PlanesNearGap.begin(), PlanesNearGap.end(), xcan+1) != PlanesNearGap.end() ) );
      }

      // Only look at adjacent bars in z (remember planes are alternating so adjacent xz bars are every two planes)
      //if (!(IsNextToGap_cand && IsNextToGap && !FoundNeighbourNearGap) && abs(DeltaPlane) > 2) continue;
      if (!(IsNextToGap_cand && IsNextToGap && !FoundNeighbourNearGap) && abs(DeltaPlane) > 1) continue;

      // Remember that we've found a neighbour near the gap
      if (IsNextToGap_cand && IsNextToGap) FoundNeighbourNearGap = true;

      // And the ground cost will update as we go along?
      // The ground cost depends on if a diagonal or not is connected
      double GroundCost = HighestCost;
      // Make a special exception when to candidates are at the gap; make the gap jump cheap
      //if (IsNextToGap && IsNextToGap_cand) GroundCost = abs(DeltaPlane)*5;
      if (IsNextToGap && IsNextToGap_cand) GroundCost = abs(DeltaPlane)*10;
      else {
        // left/right
        //if (abs(DeltaY_Trunc) == 0 && abs(DeltaPlane) == 2) GroundCost = 10;
        if (abs(DeltaY_Trunc) == 0 && abs(DeltaPlane) == 1) GroundCost = 10;
        // up/down
        else if (abs(DeltaY_Trunc) == 1 && abs(DeltaPlane) == 0) GroundCost = 10;
        // diagonal->prefer left/right followed by up/down
        //else if (abs(DeltaY_Trunc) == 1 && abs(DeltaPlane) == 2) GroundCost = 40;
        else if (abs(DeltaY_Trunc) == 1 && abs(DeltaPlane) == 1) GroundCost = 40;
        // up/down by 2 (missing one bar) -- make less preferntial than two single moves
        // necessary for missing deposits
        else if (abs(DeltaY_Trunc) == 2 && abs(DeltaPlane) == 0) GroundCost = 40;
        // diagonal and missing 1 y-bar
        // necessary for missing deposits
        //else if (abs(DeltaY_Trunc) == 2 && abs(DeltaPlane) == 2) GroundCost = 60;
        else if (abs(DeltaY_Trunc) == 2 && abs(DeltaPlane) == 1) GroundCost = 60;
        // diagonal and missing 2 bar
        // necessary for missing deposits
        //else if (abs(DeltaY_Trunc) == 3 && abs(DeltaPlane) == 2) GroundCost = 120;
        else if (abs(DeltaY_Trunc) == 3 && abs(DeltaPlane) == 1) GroundCost = 120;
        // up/down by 3 (missing three bars)
        // necessary for missing deposits
        else if (abs(DeltaY_Trunc) == 3 && abs(DeltaPlane) == 0) GroundCost = 90;
        // If we've truncated
        else if (abs(DeltaY_Trunc) == 0 && abs(DeltaPlane) == 0) GroundCost = 0;
      }
      // Catch when we're on an airgap; incur standard penalty
      //else GroundCost = 1000; // Put a large penalty in the case where we need to violate the above to cross a gap

      // This should never happen but is a very very useful check!
      if (GroundCost == HighestCost) {
        std::cout << "Missed exceptioN!" << std::endl;
        std::cout << "DeltaY_Trunc: " << DeltaY_Trunc << std::endl;
        std::cout << "DeltaY: " << DeltaY << std::endl;
        std::cout << "DeltaPlane: " << DeltaPlane << std::endl;
        std::cout << "getnotzW: " << (*jt).yw << std::endl;
        throw;
      }

      // If a greedy algorithm, have no ground cost (only condition on getting closer to the goal)
      if (IsGreedy) GroundCost = 0;

      // Remember how much it is for this node to connect to this neighbour
      (*it).Neighbours[Pointer] = GroundCost;
    }
  }

  // Remove hits that only have one neighbour?
  // Look at the first hit and see how many neighbours its neighbour has
  unsigned int nrem = 0;
  unsigned int total = Nodes.size();
  while (total > nrem && Nodes[nrem].Neighbours.size() < 1) {
    nrem++;
  }

  // If there's only one neighbour and its only neigbour is the start
  /*
     if (Nodes[nrem].Neighbours.size() == 1 &&
     Nodes[nrem+1].Neighbours.find(&Nodes[nrem]) != Nodes[nrem+1].Neighbours.end()) {
     nrem++;
     }
     */

  // Do the same for the front
  unsigned int nrem_front = 0;
  while (total > nrem_front && Nodes[total-nrem_front].Neighbours.size() < 1) {
    nrem_front++;
  }

  // Return if there are no nodes left
  if (nrem >= total || nrem_front >= total) {
    std::vector<TMS_Hit> empt;
    return empt;
  }
  // Recalculate the heuristic
  //if (nrem_front > 1) nrem_front--;
  for (auto i : Nodes) {
    i.SetHeuristicCost(Nodes.at(total-nrem_front-1));
  }

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
    for (auto& neighbour : current.Neighbours) {
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


  // Now push back the candidates!
  std::vector<TMS_Hit> returned;
  while (NodeID != 0) {
    returned.push_back(TMS_xz[NodeID]);
    NodeID = came_from[NodeID];
    if (NodeID == came_from[NodeID]) {
      returned.push_back(TMS_xz[NodeID]);
      break;
    }
  }
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
