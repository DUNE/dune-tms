{
  // Reco loses PE relative to true, so need different bounds
  REGISTER_AXIS(pe_reco, std::make_tuple("PE_{Reco}", 50, 0, 100)); 
  REGISTER_AXIS(pe_true, std::make_tuple("PE_{True}", 50, 0, 200));
  
  REGISTER_AXIS(plane, std::make_tuple("Plane", 61, 0, 60)); // Goes up to 50 for new geom
  REGISTER_AXIS(bar, std::make_tuple("Bar", 251, 0, 250)); // Goes up to ~200, but leave room for "dune simulation" tag in 2d
  
  for (int hi = 0; hi < truth.NTrueHits; hi++) {
    if (!truth.TrueRecoHitIsPedSupped[hi]) {
      GetHist("basic__hit_information__pe_reco", 
                      "Reco PE", "pe_reco")->Fill(truth.TrueRecoHitPE[hi]);
      GetHist("basic__hit_information__pe_true", 
                      "True PE", "pe_true")->Fill(truth.TrueHitPE[hi]);
      GetHist("basic__hit_information__pe_comparison", "True PE", "pe_true",
              "pe_reco")->Fill(truth.TrueHitPE[hi], truth.TrueRecoHitPE[hi]);
              
      GetHist("basic__hit_information__occupancy_plane", "Hit Plane Number",
              "plane")->Fill(truth.TrueHitPlane[hi]);
      GetHist("basic__hit_information__occupancy_bar", "Hit Bar Number",
              "bar")->Fill(truth.TrueHitBar[hi]);
      GetHist("basic__hit_information__occupancy_plane_bar", "Hit Plane vs Bar Number",
              "plane", "bar")->Fill(truth.TrueHitPlane[hi], truth.TrueHitBar[hi]);
              
    } // end if (!truth.TrueRecoHitIsPedSupped[hi])
  } // end for (int hi = 0; hi < truth.NTrueHits; hi++)
}
