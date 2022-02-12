void TMS_CAF_converter::TMS_CAF_converter() {
}

void TMS_CAF_converter::ConvertEvent(TMS_event &event) {
  std::vector<std::vector<TMS_Hit> > LineCandidates = TMS_TrackFinder::GetFinder().GetHoughCandidates();

  // CAF object for TMS
  caf::SRTMS caf;
  caf.ntracks = LineCandidates.size();
  // Loop over tracks and convert them to caf::SRTrack format
  std::vector<caf::SRTrack> tracks = caf.tracks;
  tracks.reserve(caf.ntracks);
  // Fill the number of tracks from the track finder
  for (auto it = LineCandidates.begin(); it != LineCandidates.end(); ++it) {
    // Get the track
    auto track = *it;
    // caf tracks are pretty simple right now; 
    // only start point, end point, start, direction and end direction
    TMS_Hit firsthit = track.first();
    TMS_Hit lasthit = track.back();
    // Convert the TMS_Hit hit information to SRVector3D
    SRVector3 firsthit_conv(firsthit.GetNotZ(), 200, firsthit.GetZ()); // Have no y info in the TMS
    SRVector3 lasthit_conv(lasthit.GetNotZ(), 200, lasthit.GetZ()); // Have no y info in the TMS

    tracks.start = firsthit_conv;
    tracks.end = lasthit_conv;
  }
}
