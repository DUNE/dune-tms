// Add scope to avoid cross talk with other scripts
{
// 1.0 doesn't end up in the histogram
#define bound_within_1(x) (x > 0.999) ? 0.999 : x
  REGISTER_AXIS(
      slice_completeness,
      std::make_tuple("Primary Int. E_{Vis}, Slice / Spill", 20, 0.0, 1.0));
  REGISTER_AXIS(
      slice_cleanliness,
      std::make_tuple("Slice E_{Vis}, Primary / All Ints.", 20, 0.0, 1.0));
  REGISTER_AXIS(slice_cleanliness_times_completeness,
                std::make_tuple("Cleanliness x Completeness", 20, 0.0, 1.0));
  GetHist("time_slicing__slice_completeness",
          "Slice Cleanliness, ie. How much of E_{Vis} in Slice is from Primary "
          "Int?",
          "slice_cleanliness", "#Number of Slices")
      ->Fill(bound_within_1(truth.PrimaryVertexVisibleEnergyFraction));
  GetHist(
      "time_slicing__slice_cleanliness",
      "Slice Completeness, ie. How much of Primary Int.'s E_{Vis} is in Slice?",
      "slice_completeness", "#Number of Slices")
      ->Fill(bound_within_1(truth.VertexVisibleEnergyFractionInSlice));
  GetHist("time_slicing__slice_cleanliness_times_completeness",
          "Slice Completeness Times Completeness",
          "slice_cleanliness_times_completeness", "#Number of Slices")
      ->Fill(bound_within_1(truth.VertexVisibleEnergyFractionInSlice *
                            truth.PrimaryVertexVisibleEnergyFraction));
}
