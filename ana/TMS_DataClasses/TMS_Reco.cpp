#include "TMS_Reco.h"

TMS_TrackFinder::TMS_TrackFinder() :
  nTheta(2E3), // 1 degree accuracy
  nRho(2E3), // Somewhat arbitrary: comes from checking overlays with events
  RhoMin(-20E3),
  RhoMax(20E3),
  ThetaMin(0),
  ThetaMax(2*M_PI),
  ThetaWidth((ThetaMax-ThetaMin)/nTheta),
  RhoWidth((RhoMax-RhoMin)/nRho),
  zMaxHough((TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10) // Max z for us to do Hough in, here choose transition layer
{
  // Apply the maximum Hough transform to the zx not zy: all bending happens in zx

  Accumulator_zy = new int*[nTheta];
  Accumulator_zx = new int*[nTheta];
  for (int i = 0; i < nTheta; ++i) {
    Accumulator_zy[i] = new int[nRho];
    Accumulator_zx[i] = new int[nRho];
  }

  // Make the TF1
  HoughLine_zy = new TF1("LinearHough_zy", "[0]*x+[2]/[1]", (735+411)*10, (1450+411)*10);
  HoughLine_zx = new TF1("LinearHough_zx", "[0]*x+[2]/[1]", (735+411)*10, zMaxHough);
  //HoughLine_zy = new TF1("LinearHough_zy", "[0]+[1]*x", (735+411)*10, (1450+411)*10);
  //HoughLine_zx = new TF1("LinearHough_zx", "[0]+[1]*x", (735+411)*10, zMaxHough);
  HoughLine_zy->SetLineStyle(kDashed);
  HoughLine_zx->SetLineStyle(kDashed);
  HoughLine_zy->SetLineColor(kWhite);
  HoughLine_zx->SetLineColor(kWhite);
}

void TMS_TrackFinder::FindTracks(TMS_Event &event) {

  // Reset the candidate vector
  Candidates.clear();

  // Get the hits
  std::vector<TMS_Hit> TMS_Hits = event.GetHits();

  // Maybe split up into zx and zy hits here

  // Multi-thread this?
  for (int i = 0; i < nTheta; ++i) {
    for (int j = 0; j < nRho; ++j) {
      Accumulator_zy[i][j] = 0;
      Accumulator_zx[i][j] = 0;
    }
  }

  // First run a simple Hough Transform
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it!=TMS_Hits.end(); ++it) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    //double EnergyDeposit = hit.GetE();
    //double Time = hit.GetT();
    double xhit = bar.GetX();
    double yhit = bar.GetY();
    double zhit = bar.GetZ();

    TMS_Bar::BarType Type = bar.GetBarType();
    // If z position is above region of interest, ignore hit
    if (Type == TMS_Bar::kYBar && zhit > zMaxHough) continue;

    Accumulate(xhit, yhit, zhit, Type);
  }

  // Get the maximum in the Accumulator
  double max_zy = 0;
  double max_zx = 0;
  int max_zy_theta_bin = 0;
  int max_zy_rho_bin = 0;
  int max_zx_theta_bin = 0;
  int max_zx_rho_bin = 0;
  for (int i = 0; i < nTheta; ++i) {
    for (int j = 0; j < nRho; ++j) {
      if (Accumulator_zy[i][j] > max_zy) {
        max_zy = Accumulator_zy[i][j];
        max_zy_theta_bin = i;
        max_zy_rho_bin = j;
      }
      if (Accumulator_zx[i][j] > max_zx) {
        max_zx = Accumulator_zx[i][j];
        max_zx_theta_bin = i;
        max_zx_rho_bin = j;
      }
    }
  }

  // Transform back to a y=mx+c equivalent
  double ThetaOpt_zx = ThetaMin+max_zx_theta_bin*2*M_PI/nTheta;
  double RhoOpt_zx = RhoMin+max_zx_rho_bin*(RhoMax-RhoMin)/nRho;

  double ThetaOpt_zy = ThetaMin+max_zy_theta_bin*2*M_PI/nTheta;
  double RhoOpt_zy = RhoMin+max_zy_rho_bin*(RhoMax-RhoMin)/nRho;

  HoughLine_zy->SetParameter(0, tan(ThetaOpt_zy));
  HoughLine_zy->SetParameter(1, cos(ThetaOpt_zy));
  HoughLine_zy->SetParameter(2, RhoOpt_zy);

  HoughLine_zx->SetParameter(0, tan(ThetaOpt_zx));
  HoughLine_zx->SetParameter(1, cos(ThetaOpt_zx));
  HoughLine_zx->SetParameter(2, RhoOpt_zx);

  // Then run a clustering on the Hough Transform
  // Hough transform is most likely to pick out straigh features at begining of track candidate, so start looking there
  TF1 *HoughLine = NULL;

  // Make a hard copy since we don't want to remove the original hits
  // Just make a pool of hits, and tracked hits
  // We'll pull hit out of the HitPool and put them into Candidates
  std::vector<TMS_Hit> HitPool = TMS_Hits;

  // Move hits from the pool into the candidates, and remove the candidates from the pool
  // HitPool is going to shrink or stay the same here
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end();) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    double zhit = bar.GetZ();
    TMS_Bar::BarType Type = bar.GetBarType();
    if (Type == TMS_Bar::kXBar) {
      HoughLine = HoughLine_zy;
    } else if (Type == TMS_Bar::kYBar) {
      HoughLine = HoughLine_zx;
    }

    double HoughPoint = HoughLine->Eval(zhit);

    // Hough point is inside bar -> start clustering around bar
    if (bar.Contains(HoughPoint, zhit)) {
      Candidates.push_back(std::move(hit));
      // Remove from pool of hits
      it = HitPool.erase(it);
    } else {
      ++it;
    }
  }
  /*
  std::cout << "Event " << event.GetEventNumber() << std::endl;
  std::cout << "After Hough checks: " << std::endl;
  std::cout << "Of " << TMS_Hits.size()  << " hits " << Candidates.size() << " are candidates and " << HitPool.size() << " are left" << std::endl;
  */
  // Loop over the candidates, and add new adjacent candidates to the end
  size_t CandSize = Candidates.size();
  for (size_t i = 0; i < CandSize; ++i) {
    TMS_Hit Candidate = Candidates[i];
    TMS_Bar CandidateBar = Candidate.GetBar();
    int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
    TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

    // Get the y or x coordinate
    double CandidatePos = CandidateBar.GetNotZ();
    double CandidatePosWidth = CandidateBar.GetNotZw();
    //std::cout << "Matching " << CandidatePos << "+/-" << CandidatePosWidth << " in plane " << CandidatePlaneNumber << " (" << CandidateBar.GetZ() << "+/-" << CandidateBar.GetZw() << ")" << std::endl;

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
        //std::cout << "  Found adj " << PoolPos << "+/-" << PoolPosWidth << " in plane " << PoolPlaneNumber << " (" << PoolBar.GetZ() << "+/-" << PoolBar.GetZw() << ")" << std::endl;
        // Increment the number of candidates in the vector
        CandSize++;
        jt = HitPool.erase(jt);
      } else {
        //std::cout << "  Did not find adj " << PoolPos << "+/-" << PoolPosWidth << " in plane " << PoolPlaneNumber << " (" << PoolBar.GetZ() << "+/-" << PoolBar.GetZw() << ")" << std::endl;
        ++jt;
      }
    }
  }

  // Last scan tries looks at hits that occur around the gap regions

  /*
  std::cout << "After Hough+adjacent checks: " << std::endl;
  std::cout << "Of " << TMS_Hits.size()  << " hits " << Candidates.size() << " are candidates and " << HitPool.size() << " are left" << std::endl;
  */


}

// Find the bin for the accumulator
int TMS_TrackFinder::FindBin(double rho) {
  // Since we're using uniform binning no need for binary search or similar
  int bin = (rho-RhoMin)/RhoWidth;
  return bin;
}

// xvalue is x-axis, y value is y-axis
void TMS_TrackFinder::Accumulate(double xhit, double yhit, double zhit, TMS_Bar::BarType Type) {

  float yvar = -9999999999;
  int **Accumulator = NULL;
  if (Type == TMS_Bar::kXBar) {
    yvar = yhit;
    Accumulator = Accumulator_zy;
  } else if (Type == TMS_Bar::kYBar) {
    yvar = xhit;
    Accumulator = Accumulator_zx;
  }

  // Could probably multi-thread this operation
  // Now do the Hough
  for (int i = 0; i < nTheta; ++i) {
    float theta = ThetaMin+i*ThetaWidth;

    // Now calculate rho
    float rho = (yvar-zhit*tan(theta))*cos(theta);

    // Find which rho bin this corresponds to
    int rho_bin = FindBin(rho);

    // Fill the accumulator
    Accumulator[i][rho_bin]++;
  }
}

// Convert Accumulator to a TH2D
TH2D *TMS_TrackFinder::AccumulatorToTH2D(bool zy) {

  std::string Name;
  int **Accumulator;
  if (zy) {
    Name = "TMS_TrackFinder_Accumulator_zy";
    Accumulator = Accumulator_zy;
  } else {
    Name = "TMS_TrackFinder_Accumulator_zx";
    Accumulator = Accumulator_zx;
  }
  TH2D *accumulator = new TH2D(Name.c_str(), (Name+";#theta (rad.);#rho (mm)").c_str(), nTheta, ThetaMin, ThetaMax, nRho, RhoMin, RhoMax);

  for (int i = 0; i < nTheta; ++i) {
    for (int j = 0; j < nRho; ++j) {
      accumulator->SetBinContent(i+1, j+1, Accumulator[i][j]);
    }
  }
  int maxx, maxy, maxz;
  accumulator->GetMaximumBin(maxx, maxy, maxz);
  double maxtheta = accumulator->GetXaxis()->GetBinCenter(maxx);
  double maxrho = accumulator->GetYaxis()->GetBinCenter(maxy);
  accumulator->SetTitle(Form("#splitline{%s}{#theta_{max}=%.2f #rho_{max}=%.2f}", accumulator->GetTitle(), maxtheta, maxrho));

  // Set the minimum (easier to draw)
  double maxcounts = accumulator->GetMaximum();
  accumulator->SetMinimum(maxcounts*0.8);

  return accumulator;
}
