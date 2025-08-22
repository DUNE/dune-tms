#include "TMS_ChargeID.h"
#include "TMS_Constants.h"

TMS_ChargeID::TMS_ChargeID() {

}

// Check if hit is in the far negative region of x (detector)
bool TMS_ChargeID::region1(const TMS_Hit &Hit) {
  bool is_region1 = true;
  if (Hit.GetRecoX() < TMS_Const::TMS_Magnetic_region_1_outer_edge || Hit.GetRecoX() > TMS_Const::TMS_Magnetic_region_1_and_2_border) is_region1 = false; // these numbers are from TMS_Constants.h 
  return is_region1;
}

// Check if hit is in the middle region of x (detector)
bool TMS_ChargeID::region2(const TMS_Hit &Hit) {
  bool is_region2 = true;
  if (Hit.GetRecoX() < TMS_Const::TMS_Magnetic_region_1_and_2_border || Hit.GetRecoX() > TMS_Const::TMS_Magnetic_region_2_and_3_border) is_region2 = false;  // these numbers are from TMS_Constants.h 
  return is_region2;
}

// Check if hit is in the far positive region of x (detector)
bool TMS_ChargeID::region3(const TMS_Hit &Hit) {
  bool is_region3 = true;
  if (Hit.GetRecoX() < TMS_Const::TMS_Magnetic_region_2_and_3_border || Hit.GetRecoX() > TMS_Const::TMS_Magnetic_region_3_outer_edge) is_region3 = false;   // these numbers are from TMS_Constants.h 
  return is_region3;
}

// This takes the hit positions of a track as output and outputs a number representing the charge
// Output number 13 -> muon, -13 -> antimuon, -999999999 -> particle cannot be identified by this method
int TMS_ChargeID::ID_Track_Charge(const std::vector<TMS_Hit> &Track) {
  // Initialize some local variables and empty arrays/vectors to be used
  int chargeID = -999999999;
  std::vector<TMS_Hit> muon_hits = Track;
  std::vector<TMS_Hit> region1_hits;
  std::vector<TMS_Hit> region2_hits;
  std::vector<TMS_Hit> region3_hits;
  std::vector<TMS_Hit>::iterator n_changeregion = muon_hits.begin();
  int n_plus = 0;
  int n_minus = 0;

  // Check if there are any hits in the track
  if (muon_hits.size() < 3) return chargeID;

  // Check the muon starting region
  // below are the algorithms that cut the muon hits that are in different regions,
  // and store them separately to do the calculations afterwards
  // If muon starts in region 1
  if (region1(muon_hits.front())) {
    for (std::vector<TMS_Hit>::iterator hit = muon_hits.begin(); hit != muon_hits.end(); ++hit) {
      // Check it the muon is still in the same region as the starting region.
      // Mark down the point where the muon crosses into a different region,
      // or if there is an invalid entry.
      if (!region1(*hit) || (*hit).GetRecoX() == -999) {
        n_changeregion = hit;
        break;
      }
      // Now we have only the hits in region 1
      region1_hits.push_back(*hit);
    }
    if (region2(*n_changeregion)) {
       // Fill region 2 hits if there are any
       for (std::vector<TMS_Hit>::iterator it = n_changeregion; it != muon_hits.end(); ++it) {
        if (!region2(*it) || (*it).GetRecoX() == -999) break;
        region2_hits.push_back(*it);
      }
    }
  }
  // If muon starts in region 3, similar to the above algorithm
  if (region3(muon_hits.front())) {
    for (std::vector<TMS_Hit>::iterator hit = muon_hits.begin(); hit != muon_hits.end(); ++hit) {
      if (!region3(*hit) || (*hit).GetRecoX() == -999) {
        n_changeregion = hit;
        break;
      }
      region3_hits.push_back(*hit);
    }
    if (region2(*n_changeregion)) {
      for (std::vector<TMS_Hit>::iterator it = n_changeregion; it != muon_hits.end(); ++it) {
        if (!region2(*it) || (*it).GetRecoX() == -999) break;
        region2_hits.push_back(*it);
      }
    }
  }
  // If muon starts in region 2, similar to the above algorithm
  if (region2(muon_hits.front())) {
    for (std::vector<TMS_Hit>::iterator hit = muon_hits.begin(); hit != muon_hits.end(); ++hit) {
      if (!region2(*hit) || (*hit).GetRecoX() == -999) {
        n_changeregion = hit;
        break;
      }
      region2_hits.push_back(*hit);
    }
    if (region1(*n_changeregion)) {
      for (std::vector<TMS_Hit>::iterator it = n_changeregion; it != muon_hits.end(); ++it) {
        if (!region1(*it) || (*it).GetRecoX() == -999) break;
        region1_hits.push_back(*it);
      }
    }
    if (region3(*n_changeregion)) {
      for (std::vector<TMS_Hit>::iterator it = n_changeregion; it != muon_hits.end(); ++it) {
        if (!region3(*it) || (*it).GetRecoX() == -999) break;
        region3_hits.push_back(*it);
      }
    }
  }

  // Check how many hits there are, three entries means one hit TODO????
  int total_hit_region1 = region1_hits.size();
  int total_hit_region2 = region2_hits.size();
  int total_hit_region3 = region3_hits.size();

  // Now the muon hits are collected in three different regions, do the calculation
  // The calculations use the following method: draw a line from the first to the last hit
  // and then count if there are more hits to the right of the line
  // or more to the left of the line, to decide which direction the muon is bending

  // Region 1 calculation
  if (total_hit_region1 > 2 && region1_hits.front().GetZ() < region1_hits.back().GetZ()) {
    double m = (region1_hits.front().GetRecoX() - region1_hits.back().GetRecoX()) / (region1_hits.front().GetZ() - region1_hits.back().GetZ());
    for (int i = 1; i != (total_hit_region1 - 1); ++i) {
      double x_interpolation = m * (region1_hits[i].GetZ() - region1_hits.front().GetZ()) + region1_hits.front().GetRecoX();
      double signed_dist = region1_hits[i].GetRecoX() - x_interpolation;
      if (signed_dist > 0) n_plus += 1;
      if (signed_dist < 0) n_minus += 1;
    }
  }

  // Region 2 calculation
  if (total_hit_region2 > 2 && region2_hits.front().GetZ() < region2_hits.back().GetZ()) {
    double m = (region2_hits.front().GetRecoX() - region2_hits.back().GetRecoX()) / (region2_hits.front().GetZ() - region2_hits.back().GetZ());
    for (int i = 1; i != (total_hit_region2 - 1); ++i) {
      double x_interpolation = m * (region2_hits[i].GetZ() - region2_hits.front().GetZ()) + region2_hits.front().GetRecoX();
      double signed_dist = region2_hits[i].GetRecoX() - x_interpolation;
      if (signed_dist < 0) n_plus += 1;
      if (signed_dist > 0) n_minus += 1;
    }
  }

  // Region 3 calculation
  if (total_hit_region3 > 2 && region3_hits.front().GetZ() < region3_hits.back().GetZ()) {
    double m = (region3_hits.front().GetRecoX() - region3_hits.back().GetRecoX()) / (region3_hits.front().GetZ() - region3_hits.back().GetZ());
    for (int i = 1; i != (total_hit_region3 - 1); ++i) {
      double x_interpolation = m * (region3_hits[i].GetZ() - region3_hits.front().GetZ()) + region3_hits.front().GetRecoX();
      double signed_dist = region3_hits[i].GetRecoX() - x_interpolation;
      if (signed_dist > 0) n_plus += 1;
      if (signed_dist < 0) n_minus += 1;
    }
  }

  // Now the calculation is done, identify whether this particle is a muon or an antimuon.
  // After the identificationon you can assign the chargeID to be 13 or -13
  // depending whether it is a muon or antimuon
  // if no calculation can be done then the chargeID here is defaulted to be -999999999,
  // which means this particle cannot be identified with this method
  
  // Now judge whether the particle is a muon or an antimuon
  if (n_plus < n_minus) chargeID = 13;
  if (n_plus > n_minus) chargeID = -13;

  return chargeID;
}


