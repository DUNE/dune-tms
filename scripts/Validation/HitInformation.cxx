{
  // Reco loses PE relative to true, so need different bounds
  REGISTER_AXIS(pe_reco, std::make_tuple("Hit PE_{Reco}", 50, 0, 100)); 
  REGISTER_AXIS(pe_true, std::make_tuple("Hit PE_{True}", 50, 0, 200));
  
  REGISTER_AXIS(e_reco, std::make_tuple("Hit E_{Reco} (MeV)", 50, 0, 10)); 
  REGISTER_AXIS(e_true, std::make_tuple("Hit E_{True} (MeV)", 50, 0, 10));
  
  REGISTER_AXIS(plane, std::make_tuple("Hit Plane", 61, 0, 60)); // Goes up to 50 for new geom
  REGISTER_AXIS(bar, std::make_tuple("Hit Bar", 251, 0, 250)); // Goes up to ~200, but leave room for "dune simulation" tag in 2d
  
  
  REGISTER_AXIS(hit_delta_t, std::make_tuple("Hit #Delta T (ns, reco - true)", 50, -100, 100));
  
  for (int hi = 0; hi < truth.NTrueHits; hi++) {
    if (!truth.TrueRecoHitIsPedSupped[hi]) {
      GetHist("basic__hit_information__pe_reco", 
                      "Reco Hit PE", "pe_reco")->Fill(truth.TrueRecoHitPE[hi]);
      GetHist("basic__hit_information__pe_true", 
                      "True Hit PE", "pe_true")->Fill(truth.TrueHitPE[hi]);
      GetHist("basic__hit_information__pe_comparison", "Reco vs True Hit PE", "pe_true",
              "pe_reco")->Fill(truth.TrueHitPE[hi], truth.TrueRecoHitPE[hi]);
              
      GetHist("basic__hit_information__e_reco", 
                      "Reco Hit E", "e_reco")->Fill(truth.TrueRecoHitE[hi]);
      GetHist("basic__hit_information__e_true", 
                      "True Hit E", "e_true")->Fill(truth.TrueHitE[hi]);
      GetHist("basic__hit_information__e_comparison", "Reco vs True Hit E", "e_true",
              "e_reco")->Fill(truth.TrueHitE[hi], truth.TrueRecoHitE[hi]);
              
      GetHist("basic__hit_information__occupancy_plane", "Hit Plane Number",
              "plane")->Fill(truth.TrueHitPlane[hi]);
      GetHist("basic__hit_information__occupancy_bar", "Hit Bar Number",
              "bar")->Fill(truth.TrueHitBar[hi]);
      GetHist("basic__hit_information__occupancy_plane_bar", "Hit Plane vs Bar Number",
              "plane", "bar")->Fill(truth.TrueHitPlane[hi], truth.TrueHitBar[hi]);
              
      GetHist("basic__hit_information__delta_t", "Hit #Delta T (Reco - True), All Hits",
              "hit_delta_t")->Fill(truth.TrueRecoHitT[hi] - truth.TrueHitT[hi]);
      if (truth.TrueHitView[hi] == 0)
          GetHist("basic__hit_information__delta_t_xbar_only", "Hit #Delta T (Reco - True), X Hits Only",
                  "hit_delta_t")->Fill(truth.TrueRecoHitT[hi] - truth.TrueHitT[hi]);
      if (truth.TrueHitView[hi] == 1)
          GetHist("basic__hit_information__delta_t_ybar_only", "Hit #Delta T (Reco - True), Y Hits Only",
                  "hit_delta_t")->Fill(truth.TrueRecoHitT[hi] - truth.TrueHitT[hi]);
      if (truth.TrueHitView[hi] == 2)
          GetHist("basic__hit_information__delta_t_ubar_only", "Hit #Delta T (Reco - True), U Hits Only",
                  "hit_delta_t")->Fill(truth.TrueRecoHitT[hi] - truth.TrueHitT[hi]);
      if (truth.TrueHitView[hi] == 3)
          GetHist("basic__hit_information__delta_t_vbar_only", "Hit #Delta T (Reco - True), V Hits Only",
                  "hit_delta_t")->Fill(truth.TrueRecoHitT[hi] - truth.TrueHitT[hi]);
              
    } // end if (!truth.TrueRecoHitIsPedSupped[hi])
  } // end for (int hi = 0; hi < truth.NTrueHits; hi++)
}
