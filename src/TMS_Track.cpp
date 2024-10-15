#include "TMS_Track.h"

#include "TMS_Bar.h"

void TMS_Track::Print()
{
  //0x90; // TODO: add a function here
};

std::vector<size_t> TMS_Track::findYTransitionPoints() {
  // Finds and corrects the reco y positions for all points
  // where U and V transition to another y value

  // Loop over hits but don't update positions right away
  // so we don't influence downstream hits
  std::map<size_t, std::pair<double, double>> new_y_positions;
  for (size_t i = 0; i + 1 < Hits.size(); i++) {
    auto a = Hits[i];
    auto b = Hits[i+1];
    //std::cout<<"a.GetRecoY(): "<<a.GetRecoY()<<", b.GetRecoY(): "<<b.GetRecoY()<<std::endl;
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
        double ya_uncertainty = 0;
        double yb_uncertainty = 0;
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
          auto hit_that_has_y_info = a;
        }
        new_y_positions[index_to_use] = std::make_pair(hit_that_has_y_info.GetRecoY(), 1);
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
    if (index > Hits.size()) throw std::runtime_error("Got index outside track vector size");
    Hits[index].SetRecoY(y);
    Hits[index].SetRecoYUncertainty(y_uncertainty);
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
  if (za != zb) out = (ya - yb) / (za - zb);
  return out;
}

double TMS_Track::getMaxAllowedSlope(size_t ia, size_t ib) const {
  // For a pure UV detector, the max allowed slope is 30cm / dz assuming everything stays within same section of y hits
  double dz = Hits.at(ia).GetZ() - Hits.at(ib).GetZ();
  double sign = 1;
  if (Hits.at(ia).GetRecoY() - Hits.at(ib).GetRecoY() > 0) sign = -1;
  return sign * 300 / dz;
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
  if (points.size() > 0) {
    // Fix beginning of track
    double avg_slope_to_use_front = avg_slope;
    double max_allowed_slope_front = getMaxAllowedSlope(0, points.front());
    if (avg_slope_to_use_front > max_allowed_slope_front) avg_slope_to_use_front = max_allowed_slope_front;
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
    if (avg_slope_to_use_back > max_allowed_slope_back) avg_slope_to_use_back = max_allowed_slope_back;
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
        auto a = Hits[j];
        double ya = avg_slope_to_use * (a.GetZ() - z) + y;
        a.SetRecoY(ya);
      }
    }
  }
}

void TMS_Track::ApplyTrackSmoothing() {
  std::string strategy = "simple";
  if (strategy == "simple") simpleTrackSmoothing();
  // The next level would be to do a minimization that minimizes curvature + chi2, 
  // where chi2 takes into account the uncertainty of each point. Basically almost kalman filter
}













