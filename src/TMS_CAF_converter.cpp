#include "TMS_CAF_converter.h"

//caf::SRTMS TMS_CAF_converter::ConvertEvent(TMS_Event &event) {
caf::SRTMS TMS_CAF_converter::ConvertEvent() {
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
    TMS_Hit firsthit = track.front();
    TMS_Hit lasthit = track.back();
    // Convert the TMS_Hit hit information to SRVector3D
    caf::SRVector3D firsthit_conv(firsthit.GetNotZ(), 200, firsthit.GetZ()); // Have no y info in the TMS
    caf::SRVector3D lasthit_conv(lasthit.GetNotZ(), 200, lasthit.GetZ()); // Have no y info in the TMS

    caf::SRTrack srtrack;
    srtrack.start = firsthit_conv;
    srtrack.end = lasthit_conv;
    //srtrack.push_back(track);
  }

  std::vector<std::pair<bool, TF1*> > HoughLines = TMS_TrackFinder::GetFinder().GetHoughLines();
  int nit = 0;
  for (auto it = HoughLines.begin(); it != HoughLines.end(); ++it, ++nit) {
    //double intercept = (*it).second->GetParameter(0);
    //double slope = (*it).second->GetParameter(1);

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double zlow = TMS_Const::TMS_Thin_Start;
    double zhi = TMS_Const::TMS_Thick_Start;
    double xlow = (*it).second->Eval(zlow);
    double xhi = (*it).second->Eval(zhi);

    double zlen = zhi-zlow;
    double xlen = xhi-xlow;
    double len = sqrt(xlen*xlen+zlen*zlen);
    zlen = zlen/len;
    xlen = xlen/len;

    caf::SRVector3D dir(xlen, 0, zlen); // Make the converted direction unit vector
    tracks[nit].dir = dir;
  }

  return caf;
}
