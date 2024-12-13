// Add scope to avoid cross talk with other scripts
{
  // 1.0 doesn't end up in the histogram
  #define bound_within_1(x) (x > 0.999) ? 0.999 : x
  REGISTER_AXIS(slice_completeness, std::make_tuple("Slice Completeness (Primary Int Vis E, in Slice / in Spill)", 20, 0.0, 1.0));
  REGISTER_AXIS(slice_cleanliness, std::make_tuple("Slice Cleanliness  (Slice Vis E, Primary Int / All Ints)", 20, 0.0, 1.0));
  GetHist("time_slicing__slice_completeness", "Slice Cleanliness",
    "slice_cleanliness")->Fill(bound_within_1(truth.PrimaryVertexVisibleEnergyFraction));
  GetHist("time_slicing__slice_cleanliness", "Slice Completeness",
    "slice_completeness")->Fill(bound_within_1(truth.VertexVisibleEnergyFractionInSlice));
}
