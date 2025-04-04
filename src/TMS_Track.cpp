#include "TMS_Track.h"

#include "TMS_Bar.h"

void TMS_Track::Print()
{
  //0x90; // TODO: add a function here
};


// Set the start direction of the track object, normalised so magnitude == 1
void TMS_Track::SetStartDirection(double ax, double ay, double az) {
  double mag = sqrt(ax*ax + ay*ay + az*az);
 
  StartDirection[0] = ax/mag;
  StartDirection[1] = ay/mag;
  StartDirection[2] = az/mag;
};

// Set the end direction of the track object, normalised so magnitude == 1
void TMS_Track::SetEndDirection(double ax, double ay, double az) {
  double mag = sqrt(ax*ax + ay*ay + az*az);

  EndDirection[0] = ax/mag;
  EndDirection[1] = ay/mag;
  EndDirection[2] = az/mag;
};

std::vector<size_t> TMS_Track::findYTransitionPoints() {
  // Finds and corrects the reco y positions for all points
  // where U and V transition to another y value

  // Loop over hits but don't update positions right away
  // so we don't influence downstream hits
  std::map<size_t, std::pair<double, double>> new_y_positions;
  const double uncertainty_transition_point = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UncertaintyForUVTransitionPoints();
  const double uncertainty_good = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UncertaintyGoodDirection();
  for (size_t i = 0; i + 1 < Hits.size(); i++) {
    auto a = Hits[i];
    auto b = Hits[i+1];
    // A transition point is a point where a.RecoY != b.RecoY. Check for 0.001 due to floating point imprecision
    if (abs(a.GetRecoY() - b.GetRecoY()) > 0.001) {
      // This is a hit transition
      bool foundU = false;
      bool foundV = false;
      bool foundX = false;
      bool foundY = false;
      if (a.GetBar().GetBarType() == TMS_Bar::kUBar || b.GetBar().GetBarType() == TMS_Bar::kUBar) foundU = true;
      if (a.GetBar().GetBarType() == TMS_Bar::kVBar || b.GetBar().GetBarType() == TMS_Bar::kVBar) foundV = true;
      if (a.GetBar().GetBarType() == TMS_Bar::kXBar || b.GetBar().GetBarType() == TMS_Bar::kXBar) foundX = true;
      if (a.GetBar().GetBarType() == TMS_Bar::kYBar || b.GetBar().GetBarType() == TMS_Bar::kYBar) foundY = true;
      if (foundU && foundV) {
        // Found a target transition
        // Take a simple avg for now
        // In the future, we may consider adding in the best estimate for slope * dz
        double shared_y = (a.GetRecoY() + b.GetRecoY()) * 0.5;
        double shared_z = (a.GetZ() + b.GetZ()) * 0.5;
        double slope_estimate = 0;
        double ya = slope_estimate * (a.GetZ() - shared_z) + shared_y;
        double yb = slope_estimate * (b.GetZ() - shared_z) + shared_y;
        double ya_uncertainty = uncertainty_transition_point;
        double yb_uncertainty = uncertainty_transition_point;
        // todo, if the positions already exist in map, do something clever like weighted avg
        new_y_positions[i] = std::make_pair(ya, ya_uncertainty);
        new_y_positions[i+1] = std::make_pair(yb, yb_uncertainty);
      }
      if (foundX && (foundU || foundV || foundY)) {
        // Found a target transition
        // In this case, X is treated as correct for y position
        // We're not updating x positions for now
        size_t index_to_use = i;
        auto hit_that_needs_y_info = a;
        auto hit_that_has_y_info = b;
        if (b.GetBar().GetBarType() != TMS_Bar::kXBar) {
         index_to_use = i+1;
         hit_that_needs_y_info = b;
         hit_that_has_y_info = a;
        }
        new_y_positions[index_to_use] = std::make_pair(hit_that_has_y_info.GetRecoY(), uncertainty_good);
      }
    }
  }
  // Now update the positions and uncertainties
  std::vector<size_t> out;
  for (auto update : new_y_positions) {
    auto index = update.first;
    auto position_and_uncertainty = update.second;
    double y = position_and_uncertainty.first;
    double y_uncertainty = position_and_uncertainty.second;
    Hits.at(index).SetRecoY(y);
    Hits.at(index).SetRecoYUncertainty(y_uncertainty);
    out.push_back(index);
  }
  return out;
}

double TMS_Track::getAvgYSlopeBetween(size_t ia, size_t ib) const {
  // Returns the avg slope between hit ia and hit ib
  double out = 0;
  auto a = Hits.at(ia);
  auto b = Hits.at(ib);
  double ya = a.GetRecoY();
  double yb = b.GetRecoY();
  double za = a.GetZ();
  double zb = b.GetZ();
  if (std::abs(za - zb) > 0.0001) out = (ya - yb) / (za - zb);
  return out;
}

double TMS_Track::getMaxAllowedSlope(size_t ia, size_t ib) const {
  // For a pure UV detector, the max allowed slope is 30cm / dz assuming everything stays within same section of y hits
  double dz = Hits.at(ia).GetZ() - Hits.at(ib).GetZ();
  // This shouldn't happen but if it does, return a default 10 which is very steep
  if (std::abs(dz) < 0.00001) return 10;
  double sign = 1;
  if (Hits.at(ia).GetRecoY() - Hits.at(ib).GetRecoY() > 0) sign = -1;
  const double max_y_distance = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_MaxYDistanceBetweenUVTransitionPoints();
  return sign * max_y_distance / dz;
}

void TMS_Track::simpleTrackSmoothing() {
  // Simple strategy
  // Find all points where we know y
  // Those are points where y jumps
  // Second, do linear interpolation between all known y positions
  // For front and back of track, assume avg slope behavior
  // (continuing local slope can be too sensitive to local fluctuations)

  // First find all the points where there's a transition in UV
  // So points where Y "jumps"
  auto points = findYTransitionPoints();

  // Now calculate the average slope
  double avg_slope = getAvgYSlopeBetween(0, Hits.size() - 1);
  // If we have multiple transition points in UV, then we can use them to get a more accurate slope estimate
  // But it's only accurate if the two points are from either side of the track, and not right next to each other
  if (points.size() > 1 && points.back() > points.front() + 2) avg_slope = getAvgYSlopeBetween(points.front(), points.back());

  // Assuming we have nonzero transition points, we can calculate y positions more accurately.
  // For front and back, assume that we start with one accurate y position, and then lerp using the avg slope.
  // But for very long tracks, the avg slope can't be so long that we'd have seen another transition.
  // So check for a max allowed slope
  // Yes, we may want to extend the slope as it was at the last transition, but
  // I found this was too sensitive to small changes
  // It may be better to take a local avg of some sort
  if (points.size() > 0) {
    // Fix beginning of track
    double avg_slope_to_use_front = avg_slope;
    double max_allowed_slope_front = getMaxAllowedSlope(0, points.front());
    if (std::abs(avg_slope_to_use_front) > std::abs(max_allowed_slope_front)) 
      // Rescale to equal the max slope, but retain the sign
      avg_slope_to_use_front *= (std::abs(max_allowed_slope_front) / std::abs(avg_slope_to_use_front));
    // Use these points as anchor and lerp from there
    double yf = Hits[points.front()].GetRecoY();
    double zf = Hits[points.front()].GetZ();
    for (size_t i = 0; i < points.front() && i < Hits.size(); i++) {
      auto a = Hits[i];
      double ya = avg_slope_to_use_front * (a.GetZ() - zf) + yf;
      Hits[i].SetRecoY(ya);
    }

    // Fix end of track
    double avg_slope_to_use_back = avg_slope;
    double max_allowed_slope_back = getMaxAllowedSlope(points.back(), Hits.size() - 1);
    if (std::abs(avg_slope_to_use_back) > std::abs(max_allowed_slope_back)) 
      // Rescale to equal the max slope, but retain the sign
      avg_slope_to_use_back *= (std::abs(max_allowed_slope_back) / std::abs(avg_slope_to_use_back));
    // Use these points as anchor and lerp from there
    double yb = Hits[points.back()].GetRecoY();
    double zb = Hits[points.back()].GetZ();
    for (size_t i = points.back() + 1; i < Hits.size(); i++) {
      auto a = Hits[i];
      double ya = avg_slope_to_use_back * (a.GetZ() - zb) + yb;
      Hits[i].SetRecoY(ya);
    }
  }

  // Presumably, the best estimate is a straight line between all known y positions
  // This code goes through each pair of known y positions, creates the avg slope,
  // and applies the corrected y positions
  if (points.size() > 1) {
    // Loop through each intermediate pair
    for (size_t i = 0; i+1 < points.size(); i++) {
      // Get avg slope between pair
      double avg_slope_to_use = getAvgYSlopeBetween(points.at(i), points.at(i+1));
      // Use this as starting point anchor
      double y = Hits[points.at(i)].GetRecoY();
      double z = Hits[points.at(i)].GetZ();
      // Loop through all hits between i and i+1, and set y based on avg slope
      for (size_t j = points.at(i)+1; j < points.at(i+1); j++) {
        auto& a = Hits[j];
        double ya = avg_slope_to_use * (a.GetZ() - z) + y;
        a.SetRecoY(ya);
      }
    }
  }
}

double TMS_Track::CalculateTrackSmoothnessY() {
  double out = 0;
  for (size_t i = 0; i+1 < Hits.size(); i++) {
    double dy = Hits.at(i).GetRecoY() - Hits.at(i+1).GetRecoY();
    out += std::abs(dy);
  }
  return out / Hits.size();
}

void TMS_Track::LookForHitsOutsideTMS() {
  for (auto& hit : Hits) {
    TVector3 position(hit.GetRecoX(), hit.GetRecoY(), hit.GetZ());
    if (!TMS_Geom::StaticIsInsideTMS(position)) {
      std::cout<<"Found point outside TMS with x,y,z="<<position.X()<<","<<position.Y()<<","<<position.Z()<<std::endl;
    }
  }
}

void TMS_Track::setDefaultUncertainty() {
  for (auto& hit : Hits) {
    // If this is a Y, V, or U hit, the uncertainty is ~30cm
    // If it's X, the uncertainty in y is ~5cm
    // And vice versa
    auto bar_type = hit.GetBar().GetBarType();
    double uncertainty_y = -999.0;
    const double uncertainty_good = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UncertaintyGoodDirection();
    const double uncertainty_bad = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_UncertaintyBadDirection();
    if (bar_type == TMS_Bar::kXBar) uncertainty_y = uncertainty_good; // mm
    if (bar_type == TMS_Bar::kUBar) uncertainty_y = uncertainty_bad; // mm
    if (bar_type == TMS_Bar::kVBar) uncertainty_y = uncertainty_bad; // mm
    if (bar_type == TMS_Bar::kYBar) uncertainty_y = uncertainty_bad; // mm
    if (uncertainty_y < 0) throw std::runtime_error("This shouldn't happen. Didn't find uncertainty");
    double uncertainty_x = -999.0;
    if (bar_type == TMS_Bar::kXBar) uncertainty_x = uncertainty_bad; // mm
    if (bar_type == TMS_Bar::kUBar) uncertainty_x = uncertainty_good; // mm
    if (bar_type == TMS_Bar::kVBar) uncertainty_x = uncertainty_good; // mm
    if (bar_type == TMS_Bar::kYBar) uncertainty_x = uncertainty_good; // mm
    if (uncertainty_x < 0) throw std::runtime_error("This shouldn't happen. Didn't find uncertainty");
    hit.SetRecoYUncertainty(uncertainty_y);
    hit.SetRecoXUncertainty(uncertainty_x);
  }
}

void TMS_Track::ApplyTrackSmoothing() {
  //LookForHitsOutsideTMS();
  //double initial_track_smoothness = CalculateTrackSmoothnessY();
  // Ideally this would be done in reco somewhere but idk where
  setDefaultUncertainty();
  std::string strategy = TMS_Manager::GetInstance().Get_Reco_TRACKSMOOTHING_TrackSmoothingStrategy();
  if (strategy == "simple") simpleTrackSmoothing();
  // The next level would be to do a minimization that minimizes curvature + chi2, 
  // where chi2 takes into account the uncertainty of each point. Basically almost kalman filter
  //double final_track_smoothness = CalculateTrackSmoothnessY();
  /*std::cout<<"Track smoothness initial: "<<initial_track_smoothness;
  std::cout<<",\ttrack smoothness final: "<<final_track_smoothness;
  std::cout<<",\tn hits: "<<Hits.size()<<std::endl;*/
}


